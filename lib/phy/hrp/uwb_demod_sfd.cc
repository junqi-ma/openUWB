// SPDX-License-Identifier: GPL-3.0-or-later
//
// Stage-3 SFD engine (R1): template preparation, the full-rate window search
// and the low-SNR discrete integer-symbol grid.  Split out of
// uwb_demod_cfo_sfd.cc on 2026-11-24 as a mechanical move of the previous
// stage_sfd() body plus three planned edits:
//   S1  one coherent weighted soft polarity check for both call sites,
//   S2  the 0.008 discrete-grid peak floor also opens the low-floor window
//       path,
//   S3  zero-allocation discrete grid (flat scratch slices) and a
//       batch-derived sliding window energy.
#include <openuwb/phy/hrp/rx/demod_sfd.h>
#include <openuwb/phy/hrp/rx/demod_grid.h>
#include <openuwb/phy/hrp/rx/demod_sync.h>
#include <openuwb/phy/hrp/rx/rx_cfo_opt.h>
#include <openuwb/phy/hrp/rx/rx_sync_opt.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <vector>

namespace openuwb::demod::core {
namespace {

// ---------------------------------------------------------------------------
// S1 — coherent weighted SFD polarity check.
//
// The SFD symbols are +-1 (0 = empty), so a genuine SFD position produces
// per-SYMBOL correlations whose phase follows the signed sequence.  Projecting
// each correlation onto the sequence-weighted coherent reference and comparing
// the signed/absolute projection energies distinguishes a real SFD from a
// shifted-SFD sidelobe or a noise coincidence.
//
// At most one weak symbol may flip (< 17.5% of the projection energy may sit
// anti-phase); the 4-non-zero-symbol IEEE legacy SFD keeps the strict all-same
// rule because it has no redundancy to spend.
// ---------------------------------------------------------------------------
bool check_sfd_polarity_soft(const std::complex<float>* sym_corr,
                             const std::vector<int8_t>& seq,
                             size_t n_check)
{
    if (sym_corr == nullptr || n_check == 0)
        return false;
    std::complex<double> ref(0.0, 0.0);
    for (size_t i = 0; i < n_check; ++i)
        ref += static_cast<double>(seq[i]) * std::complex<double>(sym_corr[i]);
    if (!(std::norm(ref) > 1e-24))
        return false;
    size_t n_nz = 0, pos_cnt = 0;
    double signed_sum = 0.0, abs_sum = 0.0;
    for (size_t i = 0; i < n_check; ++i) {
        const int8_t s = seq[i];
        if (s == 0)
            continue;
        ++n_nz;
        const double proj =
            (static_cast<double>(sym_corr[i].real()) * ref.real() +
             static_cast<double>(sym_corr[i].imag()) * ref.imag()) *
            static_cast<double>(s);
        if (proj > 0.0)
            ++pos_cnt;
        signed_sum += proj;
        abs_sum += std::abs(proj);
    }
    if (n_nz == 0)
        return false;
    if (n_nz <= 4)
        return pos_cnt == n_nz;      // IEEE legacy 4-symbol SFD
    if (pos_cnt < n_nz - 1)
        return false;                // allow at most one weak flip
    return signed_sum >= 0.65 * abs_sum; // <= 17.5% anti-phase energy
}

// Dependency-free |x|^2 series.  The sliding-window energy of scan_window()
// needs the leaving/entering sample powers for every lag; computing them in one
// linear pass keeps the correlation loop free of serial scalar updates.
inline void sample_power_series(const std::complex<float>* x, size_t count,
                                float* out)
{
    const float* f = reinterpret_cast<const float*>(x);
    for (size_t i = 0; i < count; ++i)
        out[i] = f[2 * i] * f[2 * i] + f[2 * i + 1] * f[2 * i + 1];
}

// Everything stage_sfd() needs from the SFD template build: the correlation
// tap set, its length and the fast-path normalization.
struct SfdFilterPrep {
    const std::complex<float>* sfd_tmpl = nullptr;
    size_t sfd_len = 0;
    float sfd_scale = 1.0f;
    bool use_real_sfd_fast = false;
    const rxopt::SparseCirTemplate* sparse_tmpl = nullptr;
    const float* real_taps = nullptr;
    size_t len = 0;
};

// Build (or reuse) the SFD template: kron(sfd_sequence, preamble_waveform),
// then L2-normalize — same as MATLAB refineTimingWithNsSfd.  The fast path
// caches the normalized template in scratch across packets; the legacy path
// keeps constructing it in scratch.corr exactly as before.  Returns false on
// the same inputs the original inline block rejected.
bool prepare_sfd_filter(const Qm35825Profile& profile, size_t sym,
                        const std::vector<std::complex<float>>& template_wf,
                        const std::vector<int8_t>& sfd_sequence,
                        DemodScratch& scratch, SfdFilterPrep& out)
{
    out = SfdFilterPrep{};
    const size_t len = sfd_sequence.size();
    out.len = len;
    const std::complex<float>* sfd_tmpl = nullptr;
    size_t sfd_len = len * sym;
    float sfd_scale = 1.0f;

    const auto* sparse_tmpl =
        detail::prepare_sparse_sync_template(profile, sym, scratch);
    if (!sparse_tmpl && template_wf.size() < sym)
        return false;
    // Prime the shared Stage-1/2/3 shaped-template cache.  The shaped template
    // is NOT used for the SFD dot: the SFD must keep the proven sparse/dense
    // numerics.  Priming here only saves the next stage a re-hash.
    if (sparse_tmpl) {
        (void)detail::ensure_shaped_sparse_template_cached(
            template_wf, sym, *sparse_tmpl, scratch);
    }
    if (!sparse_tmpl && profile.rx_fast_sfd &&
        (!scratch.real_sync_template.is_real ||
         scratch.real_sync_template.length != sym)) {
        rxopt::build_real_sync_template(template_wf,
                                        scratch.real_sync_template);
    }
    const bool use_real_sfd_fast =
        sparse_tmpl || (profile.rx_fast_sfd && scratch.real_sync_template.is_real);
    if (sparse_tmpl || profile.rx_fast_sfd) {
        float sfd_norm_sq = 0.0f;
        for (int8_t s : sfd_sequence)
            sfd_norm_sq += static_cast<float>(s * s);
        const float energy = sparse_tmpl
            ? sparse_tmpl->code_energy
            : scratch.real_sync_template.energy;
        if (sfd_norm_sq > 0.0f && energy > 0.0f)
            sfd_scale = 1.0f / std::sqrt(sfd_norm_sq * energy);
    }

    if (!use_real_sfd_fast) {
        if (profile.rx_fast_sfd) {
            const uint64_t key =
                rxopt::sfd_template_key(template_wf, sfd_sequence);
            if (key != scratch.sfd_template_key ||
                scratch.sfd_template.size() != len * sym) {
                rxopt::build_sfd_template(template_wf, sfd_sequence,
                                          scratch.sfd_template);
                scratch.sfd_template_key =
                    (scratch.sfd_template.size() == len * sym) ? key : 0;
            }
            sfd_tmpl = scratch.sfd_template.data();
            sfd_len = scratch.sfd_template.size();
        } else {
            rxopt::build_sfd_template(template_wf, sfd_sequence, scratch.corr);
            sfd_tmpl = scratch.corr.data();
            sfd_len = scratch.corr.size();
        }
        if (sfd_len == 0)
            return false;
    }

    out.sfd_tmpl = sfd_tmpl;
    out.sfd_len = sfd_len;
    out.sfd_scale = sfd_scale;
    out.use_real_sfd_fast = use_real_sfd_fast;
    out.sparse_tmpl = sparse_tmpl;
    out.real_taps = use_real_sfd_fast && !sparse_tmpl
        ? scratch.real_sync_template.taps_dup.data()
        : nullptr;
    return true;
}

} // namespace

bool stage_sfd(const std::complex<float>* rx,
               size_t n,
               const Qm35825Profile& profile,
               const TimingResult& timing,
               const std::vector<int8_t>& sfd_sequence,
               const std::vector<std::complex<float>>& template_wf,
               SfdResult& out,
               DemodScratch& scratch,
               const std::function<void(size_t, size_t)>& ensure_derotated)
{
    out = SfdResult{};
    if (sfd_sequence.empty() || timing.peak_samples.empty() ||
        template_wf.empty())
        return false;
    const size_t sym = kQm35SamplesPerSymbol;

    SfdFilterPrep prep;
    if (!prepare_sfd_filter(profile, sym, template_wf, sfd_sequence, scratch,
                            prep))
        return false;
    const size_t len = prep.len;
    const std::complex<float>* const sfd_tmpl = prep.sfd_tmpl;
    const size_t sfd_len = prep.sfd_len;
    const float sfd_scale = prep.sfd_scale;
    const bool use_real_sfd_fast = prep.use_real_sfd_fast;
    const rxopt::SparseCirTemplate* const sparse_tmpl = prep.sparse_tmpl;
    const float* const real_taps = prep.real_taps;

    // The fallback uses the full complex SFD template. The fast search
    // correlates a folded window against a single SYNC symbol instead.
    auto sfd_dot = [&](size_t j) {
        return profile.rx_fast_sfd
            ? rxopt::conjugate_dot(rx + j, sfd_tmpl, sfd_len)
            : detail::conjugate_dot_scalar(rx + j, sfd_tmpl, sfd_len);
    };

    // MATLAB strategy (refineTimingWithNsSfd): SFD start = preamble start +
    // preamble_repetitions × measured_period.  This is physically correct even
    // when the trailing preamble SYNC(s) score below the acquisition threshold
    // and are not tracked (partial train).
    // Live search: ±64-sample narrow window at that center, then the same
    // narrow window at ±1..±4 symbol offsets if the nominal center misses.
    const double period = timing.measured_period > 0.0
        ? timing.measured_period : static_cast<double>(sym);
    const int64_t expected = timing.preamble_start_sample +
        std::llround(static_cast<double>(profile.preamble_repetitions) * period);
    out.expected_start_sample = expected;
    // Score every sample in each candidate window. The fast path folds
    // the signed SFD symbols before the dot, while power covers the raw SFD.
    constexpr auto kSfdNarrowHalf = rxopt::kSfdNarrowHalf;
    constexpr auto kSfdSymbolOffsets = rxopt::kSfdSymbolOffsets;
    float best = -1.0f;
    size_t best_j = 0;
    const bool par_only = profile.rx_sfd_par_only;

    const rxopt::SyncMatchedFilter sync_filter(
        template_wf.data(), sym, real_taps, sparse_tmpl, /*volk_complex=*/false);

    // Coherent sign check over the first min(len,16) SFD symbols.  This is the
    // hard lock that keeps a low-floor / discrete-grid hit from being a
    // shifted-SFD sidelobe or a noise coincidence: 8 random symbols all match
    // a fixed +-1 pattern with probability 2^-8.  Shared by the existing
    // low-floor window path and the discrete grid.
    auto polarity_ok = [&](int64_t pos) -> bool {
        const size_t n_check = std::min<size_t>(len, 16);
        if (pos < 0 ||
            pos + static_cast<int64_t>(n_check * sym) >
                static_cast<int64_t>(n))
            return false;
        std::complex<float> sym_corr[16];
        for (size_t si = 0; si < n_check; ++si) {
            const size_t p = static_cast<size_t>(pos) + si * sym;
            sym_corr[si] = sync_filter.dot(rx + p);
        }
        return check_sfd_polarity_soft(sym_corr, sfd_sequence, n_check);
    };

    auto scan_window = [&](int64_t center, int64_t half_width,
                           int64_t backtrack_symbols, double raw_par_min = 25.0) {
        if (center < 0 || n < sfd_len)
            return false;
        const int64_t lo = std::max<int64_t>(0, center - half_width);
        const int64_t hi = std::min<int64_t>(
            static_cast<int64_t>(n - sfd_len), center + half_width);
        if (hi < lo)
            return false;
        ++out.search_windows;
        const size_t j_lo = static_cast<size_t>(lo);
        const size_t j_end = static_cast<size_t>(hi);
        float window_best = -1.0f;
        size_t window_best_j = j_lo;
        float pwr = par_only ? 0.0f : window_power_8acc(rx + j_lo, sfd_len);
        const size_t lag_count = j_end - j_lo + 1;
        // Leaving/entering sample powers for the sliding window energy, built
        // in one batch so the score loop has no serial |x|^2 dependency.
        float* off_pwr = nullptr;
        float* on_pwr = nullptr;
        if (!par_only && lag_count > 1) {
            const size_t nb = lag_count - 1;
            scratch.energy.resize(2 * nb);
            off_pwr = scratch.energy.data();
            on_pwr = off_pwr + nb;
            sample_power_series(rx + j_lo, nb, off_pwr);
            sample_power_series(rx + j_lo + sfd_len, nb, on_pwr);
        }
        scratch.metric.resize(lag_count);
        auto score_lag = [&](size_t lag, std::complex<float> corr, bool coarse) {
            const size_t j = j_lo + lag;
            if (!par_only && lag > 0)
                pwr += on_pwr[lag - 1] - off_pwr[lag - 1];
            if (coarse)
                ++out.coarse_correlations;
            else
                ++out.fine_correlations;
            const float metric = par_only ? std::norm(corr)
                : detail::normalized_correlation_power(corr, pwr);
            scratch.metric[lag] = metric;
            if (metric > window_best) {
                window_best = metric;
                window_best_j = j;
            }
        };
        if (use_real_sfd_fast) {
            const size_t fold_len = lag_count - 1 + sym;
            if (scratch.sfd_folded.size() < fold_len)
                scratch.sfd_folded.resize(fold_len);
            rxopt::fold_sfd_window(rx + j_lo, sym, sfd_sequence.data(), len,
                                   fold_len, scratch.sfd_folded.data());
            sync_filter.for_each_lag(
                scratch.sfd_folded.data(), lag_count,
                [&](size_t lag, std::complex<float> corr, bool batch) {
                    score_lag(lag, par_only ? corr : corr * sfd_scale, batch);
                });
        } else {
            for (size_t lag = 0; lag < lag_count; ++lag)
                score_lag(lag, sfd_dot(j_lo + lag), (lag & 1u) == 0u);
        }
        if (out.fine_correlations == 0)
            ++out.fine_correlations;
        bool window_passed = !par_only && (window_best >= profile.sfd_detection_threshold);
        if (par_only) {
            const rxopt::MatchedFilterPeak peak{window_best_j - j_lo, window_best, {}};
            window_passed = rxopt::raw_lag_par_passes(
                scratch.metric.data(), lag_count, peak, 12, 16, raw_par_min,
                profile.rx_sfd_par_three_samples) &&
                polarity_ok(static_cast<int64_t>(window_best_j));
        } else if (!window_passed && window_best >= kSfdGridPeakFloor &&
                   scratch.metric.size() == (j_end - j_lo + 1)) {
            const size_t best_off = window_best_j - j_lo;
            double bg_sum = 0.0;
            size_t bg_cnt = 0;
            for (size_t idx = 0; idx < scratch.metric.size(); ++idx) {
                if (idx + 12 < best_off || idx > best_off + 12) {
                    bg_sum += static_cast<double>(scratch.metric[idx]);
                    ++bg_cnt;
                }
            }
            if (bg_cnt >= 16) {
                const double bg_mean = bg_sum / static_cast<double>(bg_cnt);
                const double par =
                    static_cast<double>(window_best) / (bg_mean + 1e-12);
                if (par >= 25.0) {
                    if (window_best >= 0.08f) {
                        window_passed = true;
                    } else if (polarity_ok(static_cast<int64_t>(window_best_j))) {
                        // Low-floor fallback (kSfdGridPeakFloor <= pk < 0.08):
                        // verify that all SFD symbols coherently match
                        // sfd_sequence signs to reject shifted-SFD sidelobes
                        // (e.g. 0.0625).
                        window_passed = true;
                    }
                }
            }
        }
        if (window_passed) {
            best = window_best;
            best_j = window_best_j;
        }
        if (out.first_threshold_backtrack_symbols < 0 && window_passed) {
            out.first_threshold_backtrack_symbols = backtrack_symbols;
        }
        return window_passed;
    };

    bool sfd_passed = scan_window(expected, kSfdNarrowHalf, 0);
    // A weak centre window says nothing about an offset window. Retain full
    // sample resolution here for errors with both symbol and sample offsets.
    for (int k = 1; k <= kSfdSymbolOffsets && !sfd_passed; ++k) {
        const int64_t step = std::llround(static_cast<double>(k) * period);
        sfd_passed = scan_window(expected + step, kSfdNarrowHalf, k);
        if (sfd_passed)
            break;
        sfd_passed = scan_window(expected - step, kSfdNarrowHalf, k);
    }

    // ---- Low-SNR discrete integer-symbol grid fallback --------------------
    // Reached ONLY when the fast path above failed, so high and medium SNR
    // (and multipath) keep the existing search bit-for-bit.  At -10..-14 dB
    // in-band the SFD peak is 0.008..0.05, below the 0.3 detection threshold,
    // and the origin can be several SYMBOLS off (the low-SNR timing fold
    // resolves the grid phase but not always the absolute SYNC index).  Scan
    // the integer symbol grid +-36 symbols at delta = 0, refine +-1 sample
    // around the best grid point, then gate on: peak >= 0.008, a weak
    // discrete-grid PAR >= 1.5, and the shared 8-symbol coherent polarity
    // check.
    //
    // The PAR gate is deliberately weak (the plan suggested 5.5).  Measured on
    // this fixture the true SFD peak is only ~2x the discrete-grid background,
    // because the background grid points sit on SYNC symbols and the 4z2
    // sequence (sum = -4) correlates with a plain SYNC at 1/4 power, not at
    // zero.  A 5.5 gate would therefore reject the -14 dB peak it is meant to
    // accept.  The real false-alarm guard is the 8-symbol polarity lock
    // (2^-8 per candidate) plus the 0.008 floor, which is 15x the pure-noise
    // maximum (0.00053) measured over 200 windows.
    if (!sfd_passed) {
        constexpr auto kGridHalf = rxopt::kSfdGridHalf;
        const int64_t max_pos =
            static_cast<int64_t>(n) - static_cast<int64_t>(sfd_len);
        // Per-SYMBOL decomposition of the SFD correlation.  SFD =
        // kron(sfd_seq, sync), so a window at `pos` is the sfd_seq-weighted sum
        // of `len` single-SYMBOL correlations, and its power is the sum of the
        // same `len` single-SYMBOL powers.
        //
        // The grid spacing is `period`, which need not equal `sym`: with SFO
        // pos = expected + round(k*period), so the symbol offset relative to
        // the nominal k*sym grid differs by delta_k = round(k*period)-k*sym.
        // delta_k takes only a handful of integer values (<= 2 for +27 ppm
        // over +-36 symbols), so one symbol table per delta covers the whole
        // grid.  The old code disabled this optimisation whenever
        // round(kGridHalf*(period-sym)) != 0, i.e. exactly on an untrimmed
        // clock; the per-delta cache keeps it enabled with identical numerics.
        auto sfd_metric_at = [&](int64_t pos) -> float {
            if (pos < 0 || pos > max_pos)
                return -1.0f;
            if (use_real_sfd_fast) {
                float pwr = 0.0f;
                std::complex<float> acc(0.0f, 0.0f);
                for (size_t s = 0; s < len; ++s) {
                    const int64_t p = pos + static_cast<int64_t>(s) *
                                                static_cast<int64_t>(sym);
                    if (!par_only)
                        pwr += window_power_8acc(rx + p, sym);
                    const int8_t sg = sfd_sequence[s];
                    if (sg == 0)
                        continue;
                    const auto d = sync_filter.dot(rx + p);
                    if (sg > 0)
                        acc += d;
                    else
                        acc -= d;
                }
                if (par_only)
                    return std::norm(acc);
                acc *= sfd_scale;
                return detail::normalized_correlation_power(acc, pwr);
            }
            const std::complex<float> acc = sfd_dot(static_cast<size_t>(pos));
            if (par_only)
                return std::norm(acc);
            const float pwr = window_power_8acc(rx + pos, sfd_len);
            return detail::normalized_correlation_power(acc, pwr);
        };
        struct GridCand {
            int64_t pos;
            float metric;
        };
        constexpr size_t kMaxGridCands = 512;
        auto run_discrete_grid = [&](int64_t k_lo, int64_t k_hi) {
            if (k_hi < k_lo || max_pos < 0)
                return;
            const size_t span = static_cast<size_t>(k_hi - k_lo);
            if (span + 1 > kMaxGridCands)
                return;
            const size_t n_sym_slots = span + len;
            // Only the samples this pass actually reads need a derotation.
            // The narrow windows may have rotated a subset already; the caller
            // fills just the uncovered gaps.
            if (ensure_derotated) {
                const int64_t margin = rxopt::kSfdNarrowHalf + 4;
                const int64_t req_lo = expected +
                    std::llround(static_cast<double>(k_lo) * period) - margin;
                const int64_t req_hi = expected +
                    std::llround(static_cast<double>(k_hi) * period) +
                    static_cast<int64_t>(sfd_len) + margin;
                const size_t clo =
                    static_cast<size_t>(std::max<int64_t>(0, req_lo));
                const size_t chi = static_cast<size_t>(std::min<int64_t>(
                    static_cast<int64_t>(n), req_hi));
                if (chi > clo)
                    ensure_derotated(clo, chi);
            }
            // Phase-0 flat scratch: up to 8 per-delta slices of n_sym_slots,
            // sized once per pass so a warm packet allocates nothing.
            scratch.sfd_grid_sym_corr.resize(8 * n_sym_slots);
            scratch.sfd_grid_sym_valid.resize(8 * n_sym_slots);
            if (!par_only)
                scratch.sfd_grid_sym_pwr.resize(8 * n_sym_slots);
            std::complex<float>* const grid_corr =
                scratch.sfd_grid_sym_corr.data();
            float* const grid_pwr = scratch.sfd_grid_sym_pwr.data();
            uint8_t* const grid_valid = scratch.sfd_grid_sym_valid.data();
            std::array<uint8_t, 8> sym_built{};
            std::array<int64_t, 8> delta_value{};
            size_t n_delta = 0;
            auto delta_index = [&](int64_t k) -> int {
                const int64_t d =
                    std::llround(static_cast<double>(k) * period) -
                    k * static_cast<int64_t>(sym);
                for (size_t i = 0; i < n_delta; ++i)
                    if (delta_value[i] == d)
                        return static_cast<int>(i);
                if (n_delta >= delta_value.size())
                    return -1;
                delta_value[n_delta] = d;
                return static_cast<int>(n_delta++);
            };
            auto ensure_sym_cache = [&](size_t d) {
                if (sym_built[d])
                    return;
                sym_built[d] = 1;
                const size_t d_base = d * n_sym_slots;
                std::fill_n(grid_corr + d_base, n_sym_slots,
                            std::complex<float>(0.0f, 0.0f));
                if (!par_only)
                    std::fill_n(grid_pwr + d_base, n_sym_slots, 0.0f);
                std::fill_n(grid_valid + d_base, n_sym_slots, uint8_t{ 0 });
                for (size_t i = 0; i < n_sym_slots; ++i) {
                    const int64_t m = static_cast<int64_t>(i) + k_lo;
                    const int64_t spos = expected + delta_value[d] +
                                         m * static_cast<int64_t>(sym);
                    if (spos < 0 ||
                        spos + static_cast<int64_t>(sym) >
                            static_cast<int64_t>(n))
                        continue;
                    grid_corr[d_base + i] = sync_filter.dot(rx + spos);
                    if (!par_only)
                        grid_pwr[d_base + i] = window_power_8acc(rx + spos, sym);
                    grid_valid[d_base + i] = 1;
                }
            };
            std::array<GridCand, kMaxGridCands> cands{};
            std::array<int8_t, kMaxGridCands> cand_delta{};
            cand_delta.fill(-1);
            size_t nc = 0;
            for (int64_t k = k_lo; k <= k_hi; ++k) {
                const int64_t pos =
                    expected + std::llround(static_cast<double>(k) * period);
                if (pos < 0 || pos > max_pos)
                    continue;
                int d = -1;
                if (use_real_sfd_fast) {
                    d = delta_index(k);
                }
                if (d >= 0) {
                    ensure_sym_cache(static_cast<size_t>(d));
                    const size_t d_base = static_cast<size_t>(d) * n_sym_slots;
                    const size_t base = static_cast<size_t>(k - k_lo);
                    bool ok = true;
                    float pwr = 0.0f;
                    std::complex<float> acc_unscaled(0.0f, 0.0f);
                    for (size_t s = 0; s < len; ++s) {
                        const size_t idx = base + s;
                        if (!grid_valid[d_base + idx]) {
                            ok = false;
                            break;
                        }
                        if (!par_only)
                            pwr += grid_pwr[d_base + idx];
                        const int8_t sg = sfd_sequence[s];
                        if (sg == 0)
                            continue;
                        if (sg > 0)
                            acc_unscaled += grid_corr[d_base + idx];
                        else
                            acc_unscaled -= grid_corr[d_base + idx];
                    }
                    if (!ok)
                        continue;
                    const float metric = par_only
                        ? std::norm(acc_unscaled)
                        : detail::normalized_correlation_power(
                              acc_unscaled * sfd_scale, pwr);
                    cand_delta[nc] = static_cast<int8_t>(d);
                    cands[nc++] = { pos, metric };
                } else {
                    cands[nc++] = { pos, sfd_metric_at(pos) };
                }
            }
            // Reuse the per-symbol cache for the polarity check: when the
            // winner is exactly the best grid point, its `len` symbol
            // correlations are already in the grid slice, so the 8-16 dot
            // products of polarity_ok() are free.
            auto cache_complete_for = [&](int d, int64_t k) -> bool {
                if (d < 0)
                    return false;
                const size_t n_check = std::min<size_t>(len, 16);
                const int64_t base = k - k_lo;
                if (base < 0 ||
                    static_cast<size_t>(base) + n_check > n_sym_slots)
                    return false;
                const size_t d_base = static_cast<size_t>(d) * n_sym_slots;
                for (size_t si = 0; si < n_check; ++si)
                    if (!grid_valid[d_base + static_cast<size_t>(base) + si])
                        return false;
                return true;
            };
            auto polarity_ok_cached = [&](int d, int64_t k) -> bool {
                const size_t n_check = std::min<size_t>(len, 16);
                const size_t base = static_cast<size_t>(k - k_lo);
                const size_t d_base = static_cast<size_t>(d) * n_sym_slots;
                std::complex<float> local[16];
                for (size_t si = 0; si < n_check; ++si)
                    local[si] = grid_corr[d_base + base + si];
                return check_sfd_polarity_soft(local, sfd_sequence, n_check);
            };
            if (nc >= 8) {
                ++out.search_windows;
                size_t bi = 0;
                for (size_t i = 1; i < nc; ++i)
                    if (cands[i].metric > cands[bi].metric)
                        bi = i;
                int64_t gpos = cands[bi].pos;
                // Refine around the best grid point.  The grid already
                // evaluated delta = 0, so only the two neighbours are scored;
                // refining only the winner keeps the noise maximum at 75
                // evaluations.
                constexpr int64_t refine_rad = 1;
                int64_t best_pos = gpos;
                float best_metric = cands[bi].metric;
                for (int64_t d = -refine_rad; d <= refine_rad; ++d) {
                    if (d == 0)
                        continue;
                    const int64_t pos = gpos + d;
                    if (pos < 0 || pos > max_pos)
                        continue;
                    const float m = sfd_metric_at(pos);
                    if (m > best_metric) {
                        best_metric = m;
                        best_pos = pos;
                    }
                }
                // Background mean over the discrete grid, excluding the best
                // grid point and its +-2 symbol neighbours.
                double bg_sum = 0.0;
                size_t bg_cnt = 0;
                for (size_t i = 0; i < nc; ++i) {
                    if (std::llabs(cands[i].pos - gpos) <=
                        std::llround(2.0 * period))
                        continue;
                    bg_sum += static_cast<double>(cands[i].metric);
                    ++bg_cnt;
                }
                if (bg_cnt >= 16) {
                    const bool passes_gate = par_only
                        ? rxopt::raw_par_passes(best_metric, bg_sum, bg_cnt,
                                                kSfdGridParFloor)
                        : best_metric >= kSfdGridPeakFloor &&
                              static_cast<double>(best_metric) /
                                      std::max(bg_sum /
                                                   static_cast<double>(bg_cnt),
                                               1e-12) >=
                                  kSfdGridParFloor;
                    if (passes_gate && par_only && profile.rx_par_confirm) {
                        // On the symbol grid an ordinary SYNC correlates with
                        // 4z2, so its PAR threshold must stay low. Confirm on
                        // off-peak sample lags to distinguish SFD from noise.
                        sfd_passed =
                            scan_window(best_pos, kSfdNarrowHalf, -1, 20.0);
                    } else if (passes_gate) {
                        const int64_t g_k = std::llround(
                            static_cast<double>(gpos - expected) / period);
                        const int g_delta = cand_delta[bi];
                        const bool reuse = best_pos == gpos && g_delta >= 0 &&
                                           cache_complete_for(g_delta, g_k);
                        const bool pol_ok = reuse
                            ? polarity_ok_cached(g_delta, g_k)
                            : polarity_ok(best_pos);
                        if (pol_ok) {
                            best = best_metric;
                            best_j = static_cast<size_t>(best_pos);
                            sfd_passed = true;
                        }
                    }
                }
            }
        };
        run_discrete_grid(-kGridHalf, kGridHalf);
        // Late-acquisition rescue: when timing started >36 symbols late (e.g.
        // front of the preamble buried in burst interference), the true SFD
        // lies earlier than expected, between timing.preamble_start_sample and
        // expected - 37 * period.  Only reached on a discrete-grid miss.
        if (!sfd_passed &&
            profile.sync_strategy != SyncStrategy::LockedBounded &&
            profile.preamble_repetitions > static_cast<size_t>(kGridHalf + 1)) {
            const int64_t k_lo =
                -static_cast<int64_t>(profile.preamble_repetitions) + 1;
            const int64_t k_hi = -kGridHalf - 1;
            run_discrete_grid(k_lo, k_hi);
        }
    }
    if (!sfd_passed)
        return false;

    out.ok = true;
    out.sfd_mode = profile.sfd_mode;
    out.sfd_start_sample = static_cast<int64_t>(best_j);
    out.sfd_end_sample = static_cast<int64_t>(best_j) +
                         static_cast<int64_t>(sfd_len) - 1;
    // Keep the reported metric in normalized power units.
    if (par_only) {
        if (use_real_sfd_fast)
            best *= sfd_scale * sfd_scale;
        const float power = window_power_8acc(rx + best_j, sfd_len);
        best /= power + openuwb::core::kUwbEpsilon;
    }
    out.metric = best;
    out.polarity = 1;
    return true;
}

bool stage_sfd(const std::complex<float>* rx,
               size_t n,
               const Qm35825Profile& profile,
               const TimingResult& timing,
               const std::vector<int8_t>& sfd_sequence,
               const std::vector<std::complex<float>>& template_wf,
               SfdResult& out,
               DemodScratch& scratch,
               const std::function<void()>& prepare_grid)
{
    // Thin forwarder: the legacy callback ran once, lazily, right before the
    // first discrete-grid pass.  Bridge it to the range-aware callback by
    // firing it on that first request.
    if (!prepare_grid)
        return stage_sfd(rx, n, profile, timing, sfd_sequence, template_wf, out,
                         scratch, std::function<void(size_t, size_t)>{});
    bool prepared = false;
    const std::function<void(size_t, size_t)> ensure_derotated =
        [&](size_t, size_t) {
            if (prepared)
                return;
            prepared = true;
            prepare_grid();
        };
    return stage_sfd(rx, n, profile, timing, sfd_sequence, template_wf, out,
                     scratch, ensure_derotated);
}

} // namespace openuwb::demod::core
