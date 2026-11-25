// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unified Stage-1 Timing engine:
//   Step 0: resolve_timing_roi()           -- uniform Acquire / Locked ROI map
//   Step 1: accumulate_cir_segment()       -- Overlap-Core + K=2 progressive
//                                             normalized CIR accumulation
//   Step 2: localize_leading_edge_boxcar() -- adaptive-threshold O(M) boxcar
//                                             leading-edge localization
//   Step 3: refine_anchors_robust()        -- <=16 K=2 super-anchors + 2-pass
//                                             outlier-rejected linear fit
#include <openuwb/phy/hrp/rx/demod_timing.h>
#include <openuwb/phy/hrp/rx/demod_sync.h>
#include <openuwb/phy/hrp/rx/rx_matched_filter.h>
#include <openuwb/phy/hrp/rx/rx_sync_opt.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace openuwb::demod::core {
namespace {

constexpr int64_t kGuardBandSamples = 40;
constexpr size_t kBlockGroups = 32; // 64 symbols at K_g = 2 (SFO-safe <= 1.6 samples @ 50 ppm)
constexpr size_t kMaxAnchorCount = 16;
constexpr int64_t kAnchorSearchRadius = 8;
constexpr size_t kMaxGridPoints = 1024;

// 32-point bit-reversed order so early batches sample 0%, 50%, 25%, 75%, ...
constexpr std::array<size_t, kBlockGroups> kBitRev32 = {
    0,  16, 8,  24, 4,  20, 12, 28, 2,  18, 10, 26, 6,  22, 14, 30,
    1,  17, 9,  25, 5,  21, 13, 29, 3,  19, 11, 27, 7,  23, 15, 31
};

struct PreparedTimingFilter {
    rxopt::SyncMatchedFilter filter;
    float Et = 0.0f;
    rxopt::SyncMatchedFilter step1_filter;
    float step1_Et = 0.0f;
    const rxopt::SparseCirTemplate* sparse = nullptr;
    const rxopt::ShapedSparseTemplate* shaped = nullptr;
    const rxopt::RealSyncTemplate* real = nullptr;
};

struct TimingRoi {
    bool ok = false;
    bool is_locked = false;
    bool fast_fail = false;
    int64_t lo = 0;
    int64_t hi = 0;
    int64_t seed = -1;
};

struct TimingSegment {
    int64_t seg_lo = 0;
    int64_t seg_hi = 0;
    int64_t core_base = 0;
    int64_t phys_base = 0;
    size_t scan_len = 0;
    size_t bg_scan_len = 0;
    size_t legal_offset = 0;
    bool narrow = false;
};

struct CirDetectResult {
    bool ok = false;
    int64_t anchor_pos = -1;
    int64_t seg_lo = 0;
    int64_t seg_hi = 0;
    int64_t core_base = 0;
    size_t tau_star = 0;
    float s_peak = 0.0f;
    float bg_mean = 0.0f;
    float par = 0.0f;
    bool narrow = false;
};

PreparedTimingFilter prepare_timing_filter(
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    DemodScratch& scratch)
{
    const size_t L = template_wf.size();
    const auto* sparse_tmpl =
        detail::prepare_sparse_sync_template(profile, L, scratch);
    const rxopt::SparseCirTemplate* step1_sparse = sparse_tmpl;
    const rxopt::ShapedSparseTemplate* shaped_tmpl = nullptr;
    if (profile.rx_timing_shaped_sparse && L == 2 * kQm35ChipsPerSymbol) {
        if (!step1_sparse) {
            if (scratch.cached_seq.code_index != profile.code_index ||
                scratch.cached_seq.pcode.empty())
                scratch.cached_seq.update(profile.code_index, profile.sfd_mode);
            const auto& pcode = scratch.cached_seq.pcode;
            rxopt::ensure_sparse_cir_template(pcode.data(), pcode.size(),
                                              scratch.cir_template);
            if (scratch.cir_template.code_energy > 0.0f)
                step1_sparse = &scratch.cir_template;
        }
        if (step1_sparse) {
            shaped_tmpl = detail::ensure_shaped_sparse_template_cached(
                template_wf, L, *step1_sparse, scratch);
        }
    }
    float dense_Et = 0.0f;
    if (!sparse_tmpl && !shaped_tmpl) {
        scratch.work.assign(template_wf.begin(), template_wf.end());
        openuwb::core::uwb_l2_normalize(scratch.work);
        dense_Et = openuwb::core::uwb_template_energy(scratch.work);
        rxopt::ensure_real_sync_template(scratch.work,
                                         scratch.real_sync_template);
    }
    const float Et = shaped_tmpl
        ? shaped_tmpl->energy
        : (sparse_tmpl ? sparse_tmpl->code_energy : dense_Et);
    const auto& real = scratch.real_sync_template;
    const rxopt::SyncMatchedFilter filter(
        scratch.work.data(), L, real.is_real ? real.taps_dup.data() : nullptr,
        sparse_tmpl, true, shaped_tmpl);
    // Step 1 uses ternary sparse filter when available (~5x faster than
    // shaped_sparse with identical lag geometry and 1/P noise floor).
    const rxopt::SyncMatchedFilter s1_filter = step1_sparse
        ? rxopt::SyncMatchedFilter(nullptr, L, nullptr, step1_sparse, false,
                                   nullptr)
        : filter;
    const float s1_Et = step1_sparse ? step1_sparse->code_energy : Et;
    return PreparedTimingFilter{ filter, Et, s1_filter, s1_Et, sparse_tmpl,
                                 shaped_tmpl, &real };
}

TimingRoi resolve_timing_roi(size_t n,
                             size_t L,
                             const Qm35825Profile& profile,
                             int64_t seed_start,
                             TimingResult& out)
{
    TimingRoi roi;
    roi.seed = seed_start;
    const int64_t max_start =
        static_cast<int64_t>(n) - static_cast<int64_t>(L);
    if (max_start < 0)
        return roi;

    if (profile.sync_strategy == SyncStrategy::Locked) {
        roi.is_locked = true;
        out.bounded_attempts = 1;
        const int64_t r = std::max<int64_t>(
            1, static_cast<int64_t>(profile.locked_roi_samples));
        out.bounded_roi_samples = static_cast<size_t>(r);
        if (seed_start < 0) {
            out.bounded_failures = 1;
            return roi;
        }
        roi.lo = std::clamp<int64_t>(seed_start - r, 0, max_start);
        roi.hi = std::clamp<int64_t>(seed_start + r, 0, max_start);
        if (roi.lo > roi.hi) {
            out.bounded_failures = 1;
            return roi;
        }
        roi.ok = true;
        return roi;
    }

    roi.fast_fail =
        profile.rx_wide_interval_probe && !profile.rx_wide_probe_fallback;
    if (seed_start >= 0) {
        const int64_t margin = std::max<int64_t>(
            256, static_cast<int64_t>(profile.timing_search_margin));
        roi.lo = std::clamp<int64_t>(seed_start - margin, 0, max_start);
        roi.hi = std::clamp<int64_t>(seed_start + margin, 0, max_start);
    } else {
        roi.lo = 0;
        roi.hi = max_start;
    }
    roi.ok = (roi.lo <= roi.hi);
    return roi;
}

struct ParStats {
    size_t best_lag = 0; // relative to [legal_offset, legal_offset + scan_len)
    float s_peak = 0.0f;
    float bg_mean = 0.0f;
    float par = 0.0f;
};

ParStats evaluate_cir_par(const float* acc,
                          size_t count_groups,
                          size_t bg_scan_len,
                          size_t legal_offset,
                          size_t scan_len)
{
    ParStats st;
    if (count_groups == 0 || scan_len == 0 ||
        legal_offset + scan_len > bg_scan_len)
        return st;

    const float inv_g = 1.0f / static_cast<float>(count_groups);
    size_t best_idx = legal_offset;
    float best_sum = acc[legal_offset];
    for (size_t tau = legal_offset + 1; tau < legal_offset + scan_len; ++tau) {
        if (acc[tau] > best_sum) {
            best_sum = acc[tau];
            best_idx = tau;
        }
    }

    const bool cyclic = (bg_scan_len == kQm35SamplesPerSymbol);
    const int64_t p_int = static_cast<int64_t>(bg_scan_len);
    const int64_t peak_i = static_cast<int64_t>(best_idx);
    double bg_sum = 0.0;
    size_t bg_cnt = 0;
    for (size_t i = 0; i < bg_scan_len; ++i) {
        int64_t d = std::llabs(static_cast<int64_t>(i) - peak_i);
        if (cyclic)
            d = std::min(d, p_int - d);
        if (d > kGuardBandSamples) {
            bg_sum += static_cast<double>(acc[i]);
            ++bg_cnt;
        }
    }
    if (bg_cnt == 0) {
        for (size_t i = 0; i < bg_scan_len; ++i) {
            if (i != best_idx) {
                bg_sum += static_cast<double>(acc[i]);
                ++bg_cnt;
            }
        }
    }
    st.best_lag = best_idx - legal_offset;
    st.s_peak = best_sum * inv_g;
    st.bg_mean = bg_cnt > 0
        ? static_cast<float>((bg_sum / static_cast<double>(bg_cnt)) * inv_g)
        : 0.0f;
    st.par = st.s_peak / std::max(st.bg_mean, 1e-12f);
    return st;
}

bool batch_gate_passes(size_t groups_acc, const ParStats& st)
{
    if (groups_acc == 0)
        return false;
    if (groups_acc <= 2)
        return st.par >= 16.0f && st.s_peak >= 0.025f;
    if (groups_acc <= 4)
        return st.par >= 12.5f && st.s_peak >= 0.018f;
    if (groups_acc <= 8)
        return st.par >= 9.5f && st.s_peak >= 0.014f;
    if (groups_acc <= 16)
        return st.par >= 7.0f && st.s_peak >= 0.011f;
    if (groups_acc <= 32)
        return st.par >= 5.5f && st.s_peak >= 0.009f;
    return st.par >= 4.5f && st.s_peak >= 0.007f;
}

void accumulate_group_range(const std::complex<float>* rx,
                            const PreparedTimingFilter& pf,
                            size_t L,
                            int64_t phys_base,
                            size_t bg_scan_len,
                            size_t k_g,
                            const size_t* sym_starts,
                            size_t count,
                            float* acc,
                            DemodScratch& scratch)
{
    constexpr size_t P = kQm35SamplesPerSymbol;
    const size_t buf_len = bg_scan_len + L;
    if (k_g == 2)
        scratch.fold_buf.resize(buf_len);
    scratch.timing_fold_lag_metrics.resize(bg_scan_len);
    for (size_t i = 0; i < count; ++i) {
        const int64_t base0 =
            phys_base + static_cast<int64_t>(sym_starts[i] * P);
        const std::complex<float>* scan_ptr = rx + base0;
        if (k_g == 2) {
            const int64_t base1 = base0 + static_cast<int64_t>(P);
            rxopt::vector_add3_complex(scratch.fold_buf.data(), rx + base0,
                                       rx + base1, buf_len);
            scan_ptr = scratch.fold_buf.data();
        }
        const float pw0 = window_power_8acc(scan_ptr, L);
        rxopt::scan_normalized_sync(pf.step1_filter, scan_ptr, bg_scan_len,
                                    pw0, pf.step1_Et,
                                    scratch.timing_fold_lag_metrics.data());
        for (size_t tau = 0; tau < bg_scan_len; ++tau)
            acc[tau] += scratch.timing_fold_lag_metrics[tau];
    }
}

struct BlockRange {
    size_t group_start = 0;
    size_t group_count = 0;
};

size_t collect_block_groups(const BlockRange& blk,
                            size_t k_g,
                            size_t* out_sym_starts)
{
    size_t n_out = 0;
    for (size_t idx : kBitRev32) {
        if (idx < blk.group_count)
            out_sym_starts[n_out++] = (blk.group_start + idx) * k_g;
    }
    return n_out;
}

void align_and_merge_blocks(float* primary,
                            const float* secondary,
                            size_t bg_scan_len,
                            bool cyclic)
{
    int best_shift = 0;
    double best_dot = -1.0;
    constexpr int kMaxShift = 8;
    for (int sh = -kMaxShift; sh <= kMaxShift; ++sh) {
        double d = 0.0;
        for (size_t i = 0; i < bg_scan_len; ++i) {
            const int64_t j = static_cast<int64_t>(i) + sh;
            size_t jj = 0;
            if (cyclic) {
                const int64_t m = static_cast<int64_t>(bg_scan_len);
                jj = static_cast<size_t>((j % m + m) % m);
            } else {
                if (j < 0 || j >= static_cast<int64_t>(bg_scan_len))
                    continue;
                jj = static_cast<size_t>(j);
            }
            d += static_cast<double>(primary[i]) *
                 static_cast<double>(secondary[jj]);
        }
        if (d > best_dot) {
            best_dot = d;
            best_shift = sh;
        }
    }
    for (size_t i = 0; i < bg_scan_len; ++i) {
        const int64_t j = static_cast<int64_t>(i) + best_shift;
        if (cyclic) {
            const int64_t m = static_cast<int64_t>(bg_scan_len);
            const size_t jj = static_cast<size_t>((j % m + m) % m);
            primary[i] += secondary[jj];
        } else if (j >= 0 && j < static_cast<int64_t>(bg_scan_len)) {
            primary[i] += secondary[static_cast<size_t>(j)];
        }
    }
}

CirDetectResult accumulate_cir_segment(const std::complex<float>* rx,
                                       size_t n,
                                       size_t L,
                                       size_t n_rep,
                                       const PreparedTimingFilter& pf,
                                       const TimingSegment& seg,
                                       bool fast_fail,
                                       DemodScratch& scratch)
{
    CirDetectResult res;
    res.seg_lo = seg.seg_lo;
    res.seg_hi = seg.seg_hi;
    res.core_base = seg.core_base;
    res.narrow = seg.narrow;

    constexpr int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t max_tail_need =
        static_cast<int64_t>(seg.bg_scan_len + L);
    if (seg.phys_base < 0 ||
        seg.phys_base + max_tail_need > static_cast<int64_t>(n))
        return res;

    // Max symbols starting at core_base that fit [phys_base + bg_scan_len + L <= n].
    const int64_t rem =
        static_cast<int64_t>(n) - seg.phys_base - max_tail_need;
    const size_t max_syms =
        std::min(n_rep, static_cast<size_t>(rem / P) + 1);
    if (max_syms == 0)
        return res;

    const size_t k_g = (max_syms >= 4) ? 2 : 1;
    const size_t total_groups = max_syms / k_g;
    if (total_groups == 0)
        return res;

    // Select up to 3 SFO-safe 32-group (64-symbol) blocks, center block first.
    std::array<BlockRange, 3> blocks{};
    size_t num_blocks = 0;
    if (total_groups <= kBlockGroups) {
        blocks[num_blocks++] = { 0, total_groups };
    } else if (total_groups <= 2 * kBlockGroups) {
        const size_t second = total_groups - kBlockGroups;
        blocks[num_blocks++] = { second, kBlockGroups };
        blocks[num_blocks++] = { 0, second };
    } else {
        const size_t mid = (total_groups - kBlockGroups) / 2;
        blocks[num_blocks++] = { mid, kBlockGroups };
        blocks[num_blocks++] = { 0, kBlockGroups };
        blocks[num_blocks++] = { total_groups - kBlockGroups, kBlockGroups };
    }

    const size_t bg_len = seg.bg_scan_len;
    scratch.timing_cir_acc.assign(bg_len, 0.0f);
    float* acc = scratch.timing_cir_acc.data();

    size_t active_blk_idx = 0;
    auto fill_result = [&](const ParStats& st, size_t blk_idx) {
        const size_t mid_grp =
            blocks[blk_idx].group_start + blocks[blk_idx].group_count / 2;
        const int64_t mid_sym_off = static_cast<int64_t>(mid_grp * k_g) * P;
        res.ok = true;
        res.tau_star = st.best_lag;
        res.anchor_pos =
            seg.core_base + mid_sym_off + static_cast<int64_t>(st.best_lag);
        res.s_peak = st.s_peak;
        res.bg_mean = st.bg_mean;
        res.par = st.par;
        return res;
    };

    std::array<size_t, kBlockGroups> order0{};
    size_t len0 = collect_block_groups(blocks[0], k_g, order0.data());
    size_t done0 = 0;

    // Batch 0a: 1st K=2 group ultra-high-SNR instant exit (PAR>=28 & S_peak>=0.25).
    accumulate_group_range(rx, pf, L, seg.phys_base, bg_len, k_g,
                           order0.data(), 1, acc, scratch);
    done0 = 1;
    ParStats st = evaluate_cir_par(acc, done0, bg_len, seg.legal_offset,
                                   seg.scan_len);
    if (st.par >= 28.0f && st.s_peak >= 0.25f)
        return fill_result(st, active_blk_idx);

    // Batch 0b: complete 2 groups (4 symbols) fast exit.
    const size_t b0 = std::min<size_t>(2, len0);
    if (b0 > done0) {
        accumulate_group_range(rx, pf, L, seg.phys_base, bg_len, k_g,
                               order0.data() + done0, b0 - done0, acc, scratch);
        done0 = b0;
        st = evaluate_cir_par(acc, done0, bg_len, seg.legal_offset,
                              seg.scan_len);
    }
    if (batch_gate_passes(done0, st))
        return fill_result(st, active_blk_idx);
    // Truncated tiny buffer fallback (< 4 groups available).
    if (len0 < 4 && done0 == len0 && st.par >= 14.0f && st.s_peak >= 0.040f)
        return fill_result(st, active_blk_idx);

    // Fast-fail for blind/SIC windows (rx_wide_probe_fallback=false): bail
    // right after 2 groups if there is no clear peak, keeping empty-window
    // rejection <= ~70 us across all segments.
    if (fast_fail && (st.par < 5.5f || st.s_peak < 0.010f))
        return res;

    // Batch 1: up to 4 groups (8 symbols).
    const size_t target1 = std::min<size_t>(4, len0);
    if (target1 > done0) {
        accumulate_group_range(rx, pf, L, seg.phys_base, bg_len, k_g,
                               order0.data() + done0, target1 - done0, acc,
                               scratch);
        done0 = target1;
        st = evaluate_cir_par(acc, done0, bg_len, seg.legal_offset,
                              seg.scan_len);
        if (batch_gate_passes(done0, st))
            return fill_result(st, active_blk_idx);
    }
    if (fast_fail)
        return res;

    // If primary block is weak (e.g. covered by burst interference) and other
    // blocks exist, 4-group probe them and pick the cleanest block to deepen.
    if (num_blocks > 1 && st.par < 7.0f) {
        float best_score = st.par * std::max(st.s_peak, 1e-4f);
        for (size_t bi = 1; bi < num_blocks; ++bi) {
            std::array<size_t, kBlockGroups> ord_b{};
            const size_t len_b =
                collect_block_groups(blocks[bi], k_g, ord_b.data());
            const size_t pb = std::min<size_t>(4, len_b);
            if (pb == 0)
                continue;
            scratch.timing_cir_acc_alt.assign(bg_len, 0.0f);
            float* alt = scratch.timing_cir_acc_alt.data();
            accumulate_group_range(rx, pf, L, seg.phys_base, bg_len, k_g,
                                   ord_b.data(), pb, alt, scratch);
            const ParStats st_b = evaluate_cir_par(
                alt, pb, bg_len, seg.legal_offset, seg.scan_len);
            if (batch_gate_passes(pb, st_b)) {
                std::copy_n(alt, bg_len, acc);
                return fill_result(st_b, bi);
            }
            const float score_b = st_b.par * std::max(st_b.s_peak, 1e-4f);
            if (score_b > best_score) {
                best_score = score_b;
                active_blk_idx = bi;
                std::copy_n(alt, bg_len, acc);
                order0 = ord_b;
                len0 = len_b;
                done0 = pb;
                st = st_b;
            }
        }
    }

    // Fast-fail for blind SIC empty windows: if even the cleanest 4-group block
    // looks purely like noise, bail before deep accumulation.
    if (fast_fail && st.par < 4.0f && st.s_peak < 0.008f)
        return res;

    // Batches 2, 3, 4: 8 -> 16 -> 32 groups within the chosen SFO-safe block.
    const size_t targets[] = { 8, 16, 32 };
    for (size_t tgt : targets) {
        const size_t lim = std::min(tgt, len0);
        if (lim <= done0)
            continue;
        accumulate_group_range(rx, pf, L, seg.phys_base, bg_len, k_g,
                               order0.data() + done0, lim - done0, acc,
                               scratch);
        done0 = lim;
        st = evaluate_cir_par(acc, done0, bg_len, seg.legal_offset,
                              seg.scan_len);
        if (batch_gate_passes(done0, st))
            return fill_result(st, active_blk_idx);
    }

    // Batch 5: cross-block +/-8-sample aligned merge up to 64 groups (128 syms).
    if (num_blocks >= 2 && k_g == 2) {
        const size_t second_bi = (active_blk_idx == 0) ? 1 : 0;
        std::array<size_t, kBlockGroups> ord2{};
        const size_t len2 =
            collect_block_groups(blocks[second_bi], k_g, ord2.data());
        if (len2 >= 16) {
            scratch.timing_cir_acc_alt.assign(bg_len, 0.0f);
            float* alt = scratch.timing_cir_acc_alt.data();
            accumulate_group_range(rx, pf, L, seg.phys_base, bg_len, k_g,
                                   ord2.data(), len2, alt, scratch);
            const bool cyclic = (bg_len == kQm35SamplesPerSymbol);
            align_and_merge_blocks(acc, alt, bg_len, cyclic);
            st = evaluate_cir_par(acc, done0 + len2, bg_len, seg.legal_offset,
                                  seg.scan_len);
            if (batch_gate_passes(done0 + len2, st))
                return fill_result(st, active_blk_idx);
        }
    }

    return res;
}

struct SymSample {
    float rho = 0.0f;
    float pwr = 0.0f;
};

SymSample single_symbol_sample(const std::complex<float>* rx,
                               const PreparedTimingFilter& pf,
                               size_t L,
                               int64_t pos)
{
    const float pw = window_power_8acc(rx + pos, L);
    return SymSample{
        detail::normalized_correlation_power(pf.step1_filter.dot(rx + pos), pw,
                                             pf.step1_Et),
        pw
    };
}

float single_symbol_metric(const std::complex<float>* rx,
                           const PreparedTimingFilter& pf,
                           size_t L,
                           int64_t pos)
{
    return single_symbol_sample(rx, pf, L, pos).rho;
}

int64_t sfo_adjust_local_pos(const std::complex<float>* rx,
                             size_t n,
                             size_t L,
                             size_t n_rep,
                             const PreparedTimingFilter& pf,
                             int64_t pos,
                             int dir,
                             int64_t max_radius,
                             float min_avg_metric)
{
    constexpr int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t max_start =
        static_cast<int64_t>(n) - static_cast<int64_t>(L);
    const size_t avg_cnt = std::min<size_t>(4, std::max<size_t>(1, n_rep / 8));

    float best_m = -1.0f;
    int64_t best_d = 0;
    for (int64_t d = -max_radius; d <= max_radius; ++d) {
        float sum_m = 0.0f;
        size_t cnt = 0;
        for (size_t j = 0; j < avg_cnt; ++j) {
            const int64_t p =
                pos + d + static_cast<int64_t>(dir) * static_cast<int64_t>(j) * P;
            if (p < 0 || p > max_start)
                continue;
            sum_m += single_symbol_sample(rx, pf, L, p).rho;
            ++cnt;
        }
        if (cnt > 0) {
            const float avg_m = sum_m / static_cast<float>(cnt);
            if (avg_m > best_m) {
                best_m = avg_m;
                best_d = d;
            }
        }
    }
    if (best_d != 0 && best_m >= min_avg_metric)
        return pos + best_d;
    return pos;
}

int64_t localize_leading_edge_boxcar(const std::complex<float>* rx,
                                     size_t n,
                                     size_t L,
                                     size_t n_rep,
                                     const PreparedTimingFilter& pf,
                                     const TimingRoi& roi,
                                     const CirDetectResult& det,
                                     DemodScratch& scratch)
{
    if (det.narrow)
        return det.seg_lo + static_cast<int64_t>(det.tau_star);

    constexpr int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t max_start =
        static_cast<int64_t>(n) - static_cast<int64_t>(L);

    const int64_t left_bound =
        std::max<int64_t>(0, roi.lo - 8 * P);
    const int64_t right_bound = std::min<int64_t>(
        max_start, roi.hi + static_cast<int64_t>(n_rep + 12) * P);

    // Step outward from anchor_pos in both directions, adjusting +/-2 samples
    // every 32 symbols to follow up to +/-50 ppm SFO without losing the peak.
    const float sfo_track_floor = std::max(4.0f * det.bg_mean, 0.008f);
    // Re-center anchor locally (+/-4 samples over 8 symbols) so block-averaged
    // anchor_pos starts dead-center on the pulse even under large SFO.
    const int64_t aligned_anchor = sfo_adjust_local_pos(
        rx, n, L, n_rep, pf, det.anchor_pos, +1, 4, sfo_track_floor);

    const bool high_snr = (det.s_peak >= 0.25f && det.par >= 25.0f);
    const size_t sym_stride = (n_rep >= 64) ? (high_snr ? 4 : 2) : 1;
    const int64_t step_p = static_cast<int64_t>(sym_stride) * P;
    const size_t sfo_period_steps = std::max<size_t>(1, 32 / sym_stride);

    std::array<int64_t, kMaxGridPoints> left_pos{};
    std::array<float, kMaxGridPoints> left_rho{};
    std::array<float, kMaxGridPoints> left_pwr{};
    size_t n_left = 0;
    {
        int64_t cur = aligned_anchor;
        size_t steps = 0;
        while (n_left < kMaxGridPoints / 2 + 128) {
            const int64_t p = cur - step_p;
            if (p < left_bound)
                break;
            cur = p;
            ++steps;
            if ((steps % sfo_period_steps) == 0u)
                cur = sfo_adjust_local_pos(rx, n, L, n_rep, pf, cur, -1,
                                           static_cast<int64_t>(sym_stride),
                                           sfo_track_floor);
            const auto ss = single_symbol_sample(rx, pf, L, cur);
            left_pos[n_left] = cur;
            left_rho[n_left] = ss.rho;
            left_pwr[n_left] = ss.pwr;
            ++n_left;
        }
    }

    std::array<int64_t, kMaxGridPoints> right_pos{};
    std::array<float, kMaxGridPoints> right_rho{};
    std::array<float, kMaxGridPoints> right_pwr{};
    size_t n_right = 0;
    {
        int64_t cur = aligned_anchor;
        bool in_sig = false;
        if (cur >= left_bound && cur <= right_bound) {
            const auto ss = single_symbol_sample(rx, pf, L, cur);
            in_sig = (ss.rho >= sfo_track_floor);
            right_pos[n_right] = cur;
            right_rho[n_right] = ss.rho;
            right_pwr[n_right] = ss.pwr;
            ++n_right;
        }
        size_t steps = 0;
        while (n_left + n_right < kMaxGridPoints) {
            const int64_t p = cur + step_p;
            if (p > right_bound)
                break;
            cur = p;
            ++steps;
            auto ss = single_symbol_sample(rx, pf, L, cur);
            if (!in_sig && ss.rho >= sfo_track_floor) {
                in_sig = true;
                const int64_t adj = sfo_adjust_local_pos(
                    rx, n, L, n_rep, pf, cur, +1, 4, sfo_track_floor);
                if (adj != cur) {
                    cur = adj;
                    ss = single_symbol_sample(rx, pf, L, cur);
                }
                steps = 0;
            } else if ((steps % sfo_period_steps) == 0u) {
                const int64_t adj = sfo_adjust_local_pos(
                    rx, n, L, n_rep, pf, cur, +1,
                    static_cast<int64_t>(sym_stride), sfo_track_floor);
                if (adj != cur) {
                    cur = adj;
                    ss = single_symbol_sample(rx, pf, L, cur);
                }
            }
            right_pos[n_right] = cur;
            right_rho[n_right] = ss.rho;
            right_pwr[n_right] = ss.pwr;
            ++n_right;
        }
    }

    const size_t M = n_left + n_right;
    if (M == 0)
        return det.anchor_pos;

    std::array<int64_t, kMaxGridPoints> grid_pos{};
    std::array<float, kMaxGridPoints> grid_pwr{};
    scratch.timing_fold_cand_metrics.resize(M);
    float* rho = scratch.timing_fold_cand_metrics.data();
    for (size_t i = 0; i < n_left; ++i) {
        const size_t src = n_left - 1 - i;
        grid_pos[i] = left_pos[src];
        rho[i] = left_rho[src];
        grid_pwr[i] = left_pwr[src];
    }
    for (size_t i = 0; i < n_right; ++i) {
        grid_pos[n_left + i] = right_pos[i];
        rho[n_left + i] = right_rho[i];
        grid_pwr[n_left + i] = right_pwr[i];
    }

    // 75th percentile of single-symbol correlation inside the known core.
    std::array<float, 64> core_samples{};
    size_t n_core = 0;
    const int64_t core_end =
        det.core_base + static_cast<int64_t>(n_rep) * P;
    for (size_t i = 0; i < M && n_core < core_samples.size(); ++i) {
        if (grid_pos[i] >= det.core_base && grid_pos[i] < core_end)
            core_samples[n_core++] = rho[i];
    }
    if (n_core < 4) {
        n_core = std::min(M, core_samples.size());
        std::copy_n(rho, n_core, core_samples.begin());
    }
    std::sort(core_samples.begin(), core_samples.begin() + n_core);
    const float mu_sig_p75 =
        core_samples[std::min(n_core - 1, (3 * n_core) / 4)];

    const float theta = std::max(6.0f * det.bg_mean, 0.20f * mu_sig_p75);

    std::array<float, kMaxGridPoints + 1> pref{};
    pref[0] = 0.0f;
    for (size_t i = 0; i < M; ++i) {
        const float s = (rho[i] >= theta) ? 1.0f : -1.0f;
        pref[i + 1] = pref[i] + s;
    }

    const int64_t cand_lo = roi.lo;
    const int64_t cand_hi = roi.hi;
    float j_max = -1e9f;
    size_t best_m = 0;
    bool any_legal = false;
    std::array<float, kMaxGridPoints> j_scores{};
    std::array<uint8_t, kMaxGridPoints> is_legal{};
    const size_t box_syms = (n_rep + sym_stride - 1) / sym_stride;
    for (size_t m = 0; m < M; ++m) {
        if (grid_pos[m] < cand_lo || grid_pos[m] > cand_hi)
            continue;
        is_legal[m] = 1;
        const size_t end_m = std::min(m + box_syms, M);
        const float j = pref[end_m] - pref[m];
        j_scores[m] = j;
        if (!any_legal || j > j_max) {
            j_max = j;
            best_m = m;
            any_legal = true;
        }
    }
    if (!any_legal)
        return std::clamp(det.anchor_pos, roi.lo, roi.hi);

    // Disambiguate plateau [m_left, m_right]: skip leading pure-noise slots
    // (shifted left due to tail truncation/interference), while preserving
    // high-power burst-jammed leading preamble symbols.
    size_t m_left = best_m;
    size_t m_right = best_m;
    for (size_t m = 0; m < M; ++m) {
        if (is_legal[m] && j_scores[m] >= j_max - 0.5f) {
            m_left = m;
            break;
        }
    }
    for (size_t m = M; m > m_left; --m) {
        if (is_legal[m - 1] && j_scores[m - 1] >= j_max - 0.5f) {
            m_right = m - 1;
            break;
        }
    }
    float pwr_floor = 0.0f;
    if (rho[m_left] < theta && m_right > m_left) {
        std::array<float, 64> pwr_sample{};
        const size_t p_cnt = std::min(M, pwr_sample.size());
        for (size_t i = 0; i < p_cnt; ++i)
            pwr_sample[i] = grid_pwr[(i * M) / p_cnt];
        std::sort(pwr_sample.begin(), pwr_sample.begin() + p_cnt);
        pwr_floor = pwr_sample[p_cnt / 8];
        const float jam_pwr_thr = std::max(2.5f * pwr_floor, 1e-12f);
        size_t cur_m = m_left;
        while (cur_m < m_right && rho[cur_m] < theta &&
               grid_pwr[cur_m] <= jam_pwr_thr) {
            ++cur_m;
        }
        m_left = cur_m;
    }

    int64_t coarse_start = grid_pos[m_left];
    if (sym_stride == 1)
        return coarse_start;

    // Fine 1-symbol edge check within [-sym_stride, +sym_stride] around coarse_start
    // so stride=2 achieves exact 0-symbol leading-edge resolution.
    if (pwr_floor <= 0.0f) {
        std::array<float, 64> pwr_sample{};
        const size_t p_cnt = std::min(M, pwr_sample.size());
        for (size_t i = 0; i < p_cnt; ++i)
            pwr_sample[i] = grid_pwr[(i * M) / p_cnt];
        std::sort(pwr_sample.begin(), pwr_sample.begin() + p_cnt);
        pwr_floor = pwr_sample[p_cnt / 8];
    }
    const float jam_pwr_thr = std::max(2.5f * pwr_floor, 1e-12f);
    // Check one symbol earlier (m_left - 1 in 1-symbol units): if it is part
    // of the preamble (signal or high-power burst), step back by 1 symbol.
    for (int64_t d = 1; d < static_cast<int64_t>(sym_stride); ++d) {
        const int64_t cand = coarse_start - P;
        if (cand < cand_lo)
            break;
        const auto ss = single_symbol_sample(rx, pf, L, cand);
        if (ss.rho >= theta || ss.pwr > jam_pwr_thr)
            coarse_start = cand;
        else
            break;
    }
    return coarse_start;
}

void build_anchor_indices(size_t n_rep,
                          size_t& out_count,
                          std::array<size_t, kMaxAnchorCount>& out_idx)
{
    if (n_rep <= 2) {
        out_count = 1;
        out_idx[0] = 0;
        return;
    }
    if (n_rep >= 64) {
        // Align anchors to 4-symbol multiples so adjacent spacing mixes 16 and
        // 20 symbols (GCD = 4), matching fit_cfo()'s sparse-anchor unwrap grid.
        const size_t groups = n_rep / 4;
        const size_t max_g = groups - 1;
        const size_t na = std::min<size_t>(kMaxAnchorCount, groups);
        out_count = na;
        for (size_t i = 0; i < na; ++i)
            out_idx[i] = ((i * max_g) / (na - 1)) * 4;
        return;
    }
    const size_t max_k = n_rep - 2; // room for K_a=2 pairs
    const size_t na = std::min<size_t>(kMaxAnchorCount, max_k + 1);
    out_count = na;
    if (na == 1) {
        out_idx[0] = 0;
        return;
    }
    for (size_t i = 0; i < na; ++i)
        out_idx[i] = (i * max_k) / (na - 1);
}

bool refine_anchors_robust(const std::complex<float>* rx,
                           size_t n,
                           size_t L,
                           size_t n_rep,
                           const PreparedTimingFilter& pf,
                           const TimingRoi& roi,
                           const CirDetectResult& det,
                           int64_t start_cand,
                           const Qm35825Profile& profile,
                           TimingResult& out,
                           DemodScratch& scratch)
{
    constexpr double kNomPeriod = static_cast<double>(kQm35SamplesPerSymbol);
    const int64_t max_start =
        static_cast<int64_t>(n) - static_cast<int64_t>(L);

    std::array<size_t, kMaxAnchorCount> anchor_k{};
    size_t na = 0;
    build_anchor_indices(n_rep, na, anchor_k);

    const float anchor_floor = std::max(4.0f * det.bg_mean, 0.015f);
    scratch.fold_sup_buf.resize(L + static_cast<size_t>(2 * kAnchorSearchRadius + 2));
    scratch.timing_fold_lag_metrics.resize(
        static_cast<size_t>(2 * kAnchorSearchRadius + 1));

    struct AnchorMeas {
        size_t k = 0;
        int64_t sym_end = 0;
        float metric = 0.0f;
        std::complex<float> corr{};
    };

    std::array<AnchorMeas, kMaxAnchorCount> meas{};
    size_t n_meas = 0;
    size_t k_a_used = 1;

    auto measure_pass = [&](int64_t s_hint, double p_hint, bool fast_k1) {
        n_meas = 0;
        k_a_used = 1;
        float max_met = 0.0f;
        std::array<AnchorMeas, kMaxAnchorCount> raw_m{};
        size_t n_raw = 0;
        const int64_t rad = fast_k1 ? 4 : kAnchorSearchRadius;

        for (size_t i = 0; i < na; ++i) {
            const size_t k = anchor_k[i];
            const int64_t pred =
                s_hint + static_cast<int64_t>(std::llround(k * p_hint));
            const int64_t w0 = std::max<int64_t>(0, pred - rad);
            const size_t p_int = static_cast<size_t>(std::llround(p_hint));
            const bool can_pair =
                !fast_k1 && (k + 1 < n_rep) &&
                (pred + rad + static_cast<int64_t>(p_int) <= max_start);
            const int64_t max_w1 = can_pair
                ? (max_start - static_cast<int64_t>(p_int))
                : max_start;
            const int64_t w1 = std::min<int64_t>(max_w1, pred + rad);
            if (w1 < w0)
                continue;

            const size_t lag_cnt = static_cast<size_t>(w1 - w0 + 1);
            const std::complex<float>* scan_ptr = rx + w0;
            if (can_pair) {
                const size_t buf_len = lag_cnt + L;
                rxopt::vector_add3_complex(scratch.fold_sup_buf.data(), rx + w0,
                                           rx + w0 + p_int, buf_len);
                scan_ptr = scratch.fold_sup_buf.data();
                k_a_used = 2;
            }
            const float pw0 = window_power_8acc(scan_ptr, L);
            const auto pk = rxopt::scan_normalized_sync(
                pf.step1_filter, scan_ptr, lag_cnt, pw0, pf.step1_Et,
                scratch.timing_fold_lag_metrics.data());
            const int64_t sym_start = w0 + static_cast<int64_t>(pk.lag);
            raw_m[n_raw++] = AnchorMeas{
                k, sym_start + static_cast<int64_t>(L) - 1, pk.metric, pk.corr
            };
            if (pk.metric > max_met)
                max_met = pk.metric;
        }

        const float gate = std::max(anchor_floor, 0.15f * max_met);
        for (size_t i = 0; i < n_raw; ++i) {
            if (raw_m[i].metric >= gate)
                meas[n_meas++] = raw_m[i];
        }
    };

    auto ols_fit = [](const double* xs,
                      const double* ys,
                      size_t count,
                      double& slope_out,
                      double& intercept_out) -> bool {
        if (count < 2)
            return false;
        detail::LinearFit lf;
        for (size_t i = 0; i < count; ++i)
            lf.add(xs[i], ys[i]);
        const double det = lf.determinant();
        if (std::fabs(det) <= 1e-12)
            return false;
        slope_out = lf.slope(det);
        intercept_out = lf.intercept(slope_out);
        return std::isfinite(slope_out) && std::isfinite(intercept_out);
    };

    auto fit_robust = [&](double& period_out,
                          double& intercept_end_out,
                          std::array<AnchorMeas, kMaxAnchorCount>& kept,
                          size_t& n_kept) -> bool {
        if (n_meas < 2)
            return false;
        std::array<double, kMaxAnchorCount> xs{};
        std::array<double, kMaxAnchorCount> ys{};
        for (size_t i = 0; i < n_meas; ++i) {
            xs[i] = static_cast<double>(meas[i].k);
            ys[i] = static_cast<double>(meas[i].sym_end);
        }
        double slope0 = 0.0;
        double intercept0 = 0.0;
        if (!ols_fit(xs.data(), ys.data(), n_meas, slope0, intercept0))
            return false;

        n_kept = 0;
        if (n_meas >= 5) {
            for (size_t i = 0; i < n_meas; ++i) {
                const double pred = intercept0 + slope0 * xs[i];
                if (std::fabs(ys[i] - pred) <= 3.5) {
                    xs[n_kept] = xs[i];
                    ys[n_kept] = ys[i];
                    kept[n_kept++] = meas[i];
                }
            }
        }
        if (n_kept < 3) {
            n_kept = n_meas;
            for (size_t i = 0; i < n_meas; ++i) {
                xs[i] = static_cast<double>(meas[i].k);
                ys[i] = static_cast<double>(meas[i].sym_end);
                kept[i] = meas[i];
            }
        }
        return ols_fit(xs.data(), ys.data(), n_kept, period_out,
                       intercept_end_out);
    };

    const bool high_snr = (det.s_peak >= 0.25f && det.par >= 25.0f);
    measure_pass(start_cand, kNomPeriod, high_snr);

    double fit_period = kNomPeriod;
    double fit_end_intercept =
        static_cast<double>(start_cand + static_cast<int64_t>(L) - 1);
    std::array<AnchorMeas, kMaxAnchorCount> kept{};
    size_t n_kept = 0;
    bool fit_ok = fit_robust(fit_period, fit_end_intercept, kept, n_kept);

    const bool period_plausible =
        fit_ok && fit_period >= kNomPeriod * 0.995 &&
        fit_period <= kNomPeriod * 1.005;

    // Pass 2: if large SFO (>4 samples across preamble) is detected, re-center
    // all 16 anchor search windows with the fitted slope/intercept.
    if (period_plausible && n_kept >= 4 &&
        std::fabs(fit_period - kNomPeriod) * static_cast<double>(n_rep) > 4.0) {
        const int64_t s1 = static_cast<int64_t>(
            std::llround(fit_end_intercept - static_cast<double>(L - 1)));
        measure_pass(s1, fit_period, false);
        double p2 = fit_period;
        double e2 = fit_end_intercept;
        std::array<AnchorMeas, kMaxAnchorCount> kept2{};
        size_t n_kept2 = 0;
        if (fit_robust(p2, e2, kept2, n_kept2) && p2 >= kNomPeriod * 0.995 &&
            p2 <= kNomPeriod * 1.005) {
            fit_period = p2;
            fit_end_intercept = e2;
            kept = kept2;
            n_kept = n_kept2;
        }
    }

    const size_t min_anchors =
        std::min<size_t>(3, std::max<size_t>(1, profile.min_valid_peaks));
    int64_t final_start = start_cand;
    double final_period = kNomPeriod;

    if (fit_ok && n_kept >= min_anchors && fit_period >= kNomPeriod * 0.995 &&
        fit_period <= kNomPeriod * 1.005) {
        final_start = static_cast<int64_t>(
            std::llround(fit_end_intercept - static_cast<double>(L - 1)));
        final_period = fit_period;
    } else if (n_meas >= 1) {
        n_kept = n_meas;
        for (size_t i = 0; i < n_meas; ++i)
            kept[i] = meas[i];
        final_start = start_cand;
        final_period = kNomPeriod;
    } else {
        // Guaranteed non-empty anchor from Step 1/2 peak.
        n_kept = 1;
        const float pw = window_power_8acc(rx + start_cand, L);
        const auto c = pf.filter.dot(rx + start_cand);
        const float m = detail::normalized_correlation_power(c, pw, pf.Et);
        kept[0] = AnchorMeas{ 0, start_cand + static_cast<int64_t>(L) - 1, m,
                              c };
        final_start = start_cand;
        final_period = kNomPeriod;
    }

    final_start = std::clamp<int64_t>(final_start, 0, max_start);
    if (roi.is_locked) {
        const int64_t max_drift =
            static_cast<int64_t>(out.bounded_roi_samples) + 16;
        if (std::llabs(final_start - roi.seed) > max_drift) {
            out.bounded_failures = 1;
            return false;
        }
    }

    out.peak_samples.clear();
    out.peak_metrics.clear();
    out.peak_corr.clear();
    out.peak_samples.reserve(n_kept);
    out.peak_metrics.reserve(n_kept);
    out.peak_corr.reserve(n_kept);

    std::array<float, kMaxAnchorCount> met_sort{};
    for (size_t i = 0; i < n_kept; ++i) {
        const int64_t sym_start =
            kept[i].sym_end - static_cast<int64_t>(L) + 1;
        // Re-evaluate final kept anchors (<=16 single-point dots) with pf.filter
        // so out.metric/peak_corr keep exact shaped/dense golden calibration.
        float met = kept[i].metric;
        std::complex<float> corr = kept[i].corr;
        if (pf.shaped && met >= 0.30f && sym_start >= 0 &&
            sym_start <= max_start) {
            const float pw = window_power_8acc(rx + sym_start, L);
            corr = pf.filter.dot(rx + sym_start);
            met = detail::normalized_correlation_power(corr, pw, pf.Et);
        }
        out.peak_samples.push_back(kept[i].sym_end);
        out.peak_metrics.push_back(met);
        out.peak_corr.push_back(corr);
        met_sort[i] = met;
    }
    std::nth_element(met_sort.begin(), met_sort.begin() + n_kept / 2,
                     met_sort.begin() + n_kept);

    std::complex<float> lag1_sum(0.0f, 0.0f);
    const size_t phase_cnt = std::min<size_t>(n_kept, 8);
    for (size_t i = 0; i < phase_cnt; ++i) {
        const int64_t s0 = kept[i].sym_end - static_cast<int64_t>(L) + 1;
        const int64_t s1 = s0 + static_cast<int64_t>(kQm35SamplesPerSymbol);
        if (s0 >= 0 && s1 <= max_start) {
            const auto d0 = pf.step1_filter.dot(rx + s0);
            const auto d1 = pf.step1_filter.dot(rx + s1);
            lag1_sum += std::conj(d0) * d1;
        }
    }
    out.fold_phase_per_symbol =
        (std::norm(lag1_sum) > 0.0f) ? static_cast<double>(std::arg(lag1_sum)) : 0.0;

    out.ok = true;
    out.preamble_start_sample = final_start;
    out.measured_period = final_period;
    out.metric = met_sort[n_kept / 2];
    out.detected_peaks = out.peak_samples.size();
    out.expected_peaks = n_rep;
    out.cir_grid_repetitions = n_rep;
    out.sync_peaks_folded = (k_a_used > 1);
    out.bounded_success = roi.is_locked ? 1u : 0u;
    out.bounded_legacy_scans = 0u;
    return true;
}

} // namespace

bool stage_timing(const std::complex<float>* rx,
                  size_t n,
                  const Qm35825Profile& profile,
                  const std::vector<std::complex<float>>& template_wf,
                  int64_t seed_start,
                  TimingResult& out,
                  DemodScratch& scratch)
{
    out = TimingResult{};
    if (!rx || template_wf.empty() || n < template_wf.size())
        return false;

    const size_t L = template_wf.size();
    const size_t n_rep = std::max<size_t>(1, profile.preamble_repetitions);
    const PreparedTimingFilter pf =
        prepare_timing_filter(profile, template_wf, scratch);

    const TimingRoi roi = resolve_timing_roi(n, L, profile, seed_start, out);
    if (!roi.ok)
        return false;

    const auto t_step1_begin = std::chrono::steady_clock::now();
    constexpr int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t max_start =
        static_cast<int64_t>(n) - static_cast<int64_t>(L);
    const int64_t span = roi.hi - roi.lo;

    std::array<TimingSegment, 16> segments{};
    size_t n_seg = 0;
    if (span < P) {
        TimingSegment& seg = segments[n_seg++];
        seg.narrow = true;
        seg.seg_lo = roi.lo;
        seg.seg_hi = roi.hi;
        seg.core_base = roi.lo;
        seg.scan_len = static_cast<size_t>(span + 1);
        seg.bg_scan_len = std::max<size_t>(
            seg.scan_len, std::min<size_t>(static_cast<size_t>(P), 192));
        const int64_t pad = static_cast<int64_t>(
            (seg.bg_scan_len - seg.scan_len) / 2);
        const int64_t max_phys =
            max_start - static_cast<int64_t>(seg.bg_scan_len);
        seg.phys_base = std::clamp<int64_t>(
            seg.core_base - pad, 0, std::max<int64_t>(0, max_phys));
        if (seg.core_base < seg.phys_base)
            seg.phys_base = seg.core_base;
        if (seg.phys_base + static_cast<int64_t>(seg.bg_scan_len) > max_start)
            seg.bg_scan_len = std::max(
                seg.scan_len,
                static_cast<size_t>(
                    std::max<int64_t>(0, max_start - seg.phys_base)));
        seg.legal_offset = static_cast<size_t>(seg.core_base - seg.phys_base);
    } else {
        const int64_t w_seg = std::max<int64_t>(
            P, static_cast<int64_t>(n_rep / 2) * P);
        const size_t min_core_syms =
            std::max<size_t>(1, std::min<size_t>(8, n_rep / 4));
        const int64_t max_core = std::max<int64_t>(
            roi.lo, max_start - static_cast<int64_t>(min_core_syms) * P);

        int64_t cur = roi.lo;
        while (cur <= roi.hi && n_seg < segments.size()) {
            int64_t end = std::min(roi.hi, cur + w_seg);
            if (end > max_core && cur <= max_core)
                end = max_core;
            TimingSegment& seg = segments[n_seg++];
            seg.narrow = false;
            seg.seg_lo = cur;
            seg.seg_hi = end;
            if (cur > max_core) {
                seg.core_base = cur;
                seg.phys_base = cur;
                seg.scan_len = static_cast<size_t>(
                    std::min<int64_t>(P, max_start - cur + 1));
                seg.bg_scan_len = seg.scan_len;
            } else {
                seg.core_base = end;
                seg.phys_base = end;
                seg.scan_len = static_cast<size_t>(P);
                seg.bg_scan_len = static_cast<size_t>(P);
            }
            seg.legal_offset = 0;
            if (end >= roi.hi || cur > max_core)
                break;
            cur = end + 1;
        }
        if (roi.seed >= 0 && n_seg > 1) {
            std::stable_sort(
                segments.begin(), segments.begin() + n_seg,
                [&](const TimingSegment& a, const TimingSegment& b) {
                    const int64_t da = std::max<int64_t>(
                        0, std::max(a.seg_lo - roi.seed, roi.seed - a.seg_hi));
                    const int64_t db = std::max<int64_t>(
                        0, std::max(b.seg_lo - roi.seed, roi.seed - b.seg_hi));
                    return da < db;
                });
        }
    }

    CirDetectResult det;
    for (size_t si = 0; si < n_seg; ++si) {
        det = accumulate_cir_segment(rx, n, L, n_rep, pf, segments[si],
                                     roi.fast_fail, scratch);
        if (det.ok)
            break;
    }

    const auto t_step1_end = std::chrono::steady_clock::now();
    out.coarse_search_us = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            t_step1_end - t_step1_begin)
            .count());

    if (!det.ok) {
        if (roi.is_locked)
            out.bounded_failures = 1;
        return false;
    }

    const auto t_step23_begin = std::chrono::steady_clock::now();
    const int64_t start_cand = localize_leading_edge_boxcar(
        rx, n, L, n_rep, pf, roi, det, scratch);

    const bool ok = refine_anchors_robust(rx, n, L, n_rep, pf, roi, det,
                                          start_cand, profile, out, scratch);
    const auto t_step23_end = std::chrono::steady_clock::now();
    out.fine_track_us = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            t_step23_end - t_step23_begin)
            .count());
    return ok;
}

} // namespace openuwb::demod::core
