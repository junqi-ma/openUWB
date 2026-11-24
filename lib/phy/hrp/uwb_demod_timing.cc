// SPDX-License-Identifier: GPL-3.0-or-later
//
// Out-of-line definitions for the HRP timing stage.  The declarations stay in
// the public headers (demod_timing.h / demod_sync.h / rx_sync_opt.h /
// rx_matched_filter.h); this translation unit hosts the coarse stage entry
// point and the matched-filter scan loops.
//
// stage_timing() is deliberately only the orchestrator: filter preparation,
// the two acquisition strategies, the locked grid tracker and the legacy
// per-SYMBOL train live in file-local helpers below so each policy (bounded
// tracking, blind acquisition, low-SNR fold rescue) can be read on its own.
// Every helper keeps the original data flow -- the shared Result lives in
// `TimingEnv::out`, the per-search acquisition state in `AcquisitionState` --
// so the decomposition is purely structural: no field, ordering or threshold
// changes.
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

// ---------------------------------------------------------------------------
// Prepared matched-filter state
// ---------------------------------------------------------------------------
// Everything the timing stage derives from the packet's SYNC template in one
// immutable bundle: the (dense / sparse / shaped / real-dup) matched filter,
// its normalizer Et, and the template views the fold helpers need.  Building it
// once keeps the dense-template copy / L2 normalize / shaped-cache probe off the
// helper signatures.
struct PreparedTimingFilter {
    rxopt::SyncMatchedFilter filter; // dense or sparse/shaped, chosen per profile
    float Et = 0.0f;                 // matched-filter normalizer
    const rxopt::SparseCirTemplate* sparse = nullptr;
    const rxopt::ShapedSparseTemplate* shaped = nullptr;
    const rxopt::RealSyncTemplate* real = nullptr;
};

// Immutable per-call environment shared by every helper.
struct TimingEnv {
    const std::complex<float>* rx;
    size_t n;
    size_t L;
    const Qm35825Profile& profile;
    const std::vector<std::complex<float>>& template_wf;
    DemodScratch& scratch;
    TimingResult& out;
    PreparedTimingFilter pf;
    // Seeded search window (absolute coordinates).  Defaults to the whole
    // buffer; narrowed around seed_start by stage_timing().
    int64_t roi_start = 0;
    int64_t roi_end = 0;
};

// The acquisition outcome that is shared between the two strategy branches and
// the common backtrack / tracking tail.  `first_*` is the best confirmed SYNC
// start so far; `fold_used` records that it came from the K-symbol fold;
// `acq_seed` is the non-bounded seed hoisted for the post-train fold retry.
struct AcquisitionState {
    bool raw_par_hit = false;
    int64_t first_start = -1;
    size_t first_hit_k = 0;
    float first_metric = 0.0f;
    std::complex<float> first_corr{ 0.0f, 0.0f };
    bool fold_used = false;
    int64_t acq_seed = -1;
};

// Build the matched-filter state exactly once per packet.
//
// The sparse path never reads the dense template, so building it would be pure
// waste (copy + L2 normalize + duplicate real taps on every call).  Only the
// dense path materializes the normalized template, and the real duplicate is
// reused when its content key matches.
PreparedTimingFilter prepare_timing_filter(
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    DemodScratch& scratch)
{
    const size_t L = template_wf.size();
    const auto* sparse_tmpl =
        detail::prepare_sparse_sync_template(profile, L, scratch);
    // MF-4: factorize the dense template into (r * code) * p.  This replaces
    // the dense matched filter (same metric scale, same lag convention, same
    // pulse-matched gain) at a fraction of its cost; it does not change the
    // sparse path's peak-sample convention.  The builder rejects any template
    // that is not code*p, so a custom waveform silently keeps its old path.
    // The code template is needed even when the sparse timing path is off.
    const rxopt::ShapedSparseTemplate* shaped_tmpl = nullptr;
    if (profile.rx_timing_shaped_sparse && L == 2 * kQm35ChipsPerSymbol) {
        const rxopt::SparseCirTemplate* code_tmpl = sparse_tmpl;
        if (!code_tmpl) {
            if (scratch.cached_seq.code_index != profile.code_index ||
                scratch.cached_seq.pcode.empty())
                scratch.cached_seq.update(profile.code_index, profile.sfd_mode);
            const auto& pcode = scratch.cached_seq.pcode;
            rxopt::ensure_sparse_cir_template(pcode.data(), pcode.size(),
                                              scratch.cir_template);
            if (scratch.cir_template.code_energy > 0.0f)
                code_tmpl = &scratch.cir_template;
        }
        if (code_tmpl)
            shaped_tmpl = detail::ensure_shaped_sparse_template_cached(
                template_wf, L, *code_tmpl, scratch);
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
    return PreparedTimingFilter{ filter, Et, sparse_tmpl, shaped_tmpl, &real };
}

// Local normalized-correlation helper over a small window around `want`:
// finds the SYNC START that best matches the template, returns position,
// metric, and the complex matched-filter value at the best position (the
// CFO stage uses its phase).
bool local_best(TimingEnv& env,
                AcquisitionState& st,
                int64_t want,
                int64_t radius,
                float min_metric,
                int64_t* pos,
                float* met,
                std::complex<float>* corr_out,
                bool allow_par_fallback = false,
                float* raw_max_out = nullptr)
{
    const auto& profile = env.profile;
    const auto& filter = env.pf.filter;
    const float Et = env.pf.Et;
    const size_t L = env.L;
    const auto* rx = env.rx;
    const size_t n = env.n;
    const int64_t w0 = std::max<int64_t>(0, want - radius);
    const int64_t w1 =
        std::min<int64_t>(static_cast<int64_t>(n) - static_cast<int64_t>(L),
                          want + radius);
    if (w1 <= w0)
        return false;
    const int64_t count = w1 - w0;
    const bool track_par =
        allow_par_fallback && count >= 64 && count <= 1024;
    constexpr float kSyncParFloor = 0.06f;
    constexpr float kSyncParMin = 20.0f;
    constexpr int64_t kSyncParGuard = 12;
    float lag_metrics[1024];
    if (track_par && profile.rx_sync_par_only) {
        const auto peak = rxopt::scan_raw_sync(
            filter, rx + w0, static_cast<size_t>(count), lag_metrics);
        const bool accepted = rxopt::raw_lag_par_passes(
            lag_metrics, static_cast<size_t>(count), peak, kSyncParGuard, 1,
            kSyncParMin, profile.rx_sync_par_three_samples);
        const int64_t bp = w0 + static_cast<int64_t>(peak.lag);
        // Keep the normalized metric for downstream tracking and the
        // locked probe-width heuristic. Only score the winner.
        float normalized = 0.0f;
        if (accepted || raw_max_out)
            normalized = detail::normalized_correlation_power(
                peak.corr, window_power_8acc(rx + bp, L), Et);
        if (raw_max_out)
            *raw_max_out = normalized;
        if (!accepted)
            return false;
        st.raw_par_hit = true;
        *pos = bp;
        *met = normalized;
        *corr_out = peak.corr;
        return true;
    }
    const auto peak = rxopt::scan_normalized_sync(
        filter, rx + w0, static_cast<size_t>(count),
        window_power_8acc(rx + w0, L), Et, track_par ? lag_metrics : nullptr);
    const float bm = peak.metric;
    const int64_t bp = w0 + static_cast<int64_t>(peak.lag);
    // Report the window's raw maximum before the acceptance gate: the
    // low-SNR early-stop decisions below need to know how far the best
    // candidate was from the threshold, not just whether it passed.
    if (raw_max_out)
        *raw_max_out = bm;
    if (bm < min_metric) {
        bool par_accept = false;
        if (track_par && bm >= kSyncParFloor) {
            double bg_sum = 0.0;
            int64_t bg_cnt = 0;
            const int64_t bp_idx = bp - w0;
            for (int64_t idx = 0; idx < count; ++idx) {
                if (std::abs(idx - bp_idx) > kSyncParGuard) {
                    bg_sum += static_cast<double>(lag_metrics[idx]);
                    ++bg_cnt;
                }
            }
            if (bg_cnt > 0) {
                const double bg_mean =
                    std::max(bg_sum / static_cast<double>(bg_cnt), 1e-9);
                par_accept =
                    (static_cast<double>(bm) / bg_mean) >= kSyncParMin;
            }
        }
        if (!par_accept)
            return false;
    }
    *pos = bp;
    *met = bm;
    *corr_out = peak.corr;
    return true;
}

// Anti-false-alarm floor for the local K=4 fast fold (see try_fold_acquire's
// fast_only gate).  peak_metrics in fold_anchor_fit are AMPLITUDE-normalized
// (|corr|/sqrt(E*Et)); 16-group K=4 anchors at -12/-14 dB sit ~0.3-0.5 while
// 4-symbol-folded noise sits ~0.045, so 0.08 is a strict anti-false-alarm floor.
constexpr float kFoldFastAnchorMedianFloor = 0.08f;

// ---- K-symbol periodic-fold acquisition -----------------------------------
// Reached ONLY after every single-symbol path failed.  The fold costs
// ~1 ms (fold + full-window correlation), so it must never run on a packet
// the normal path already acquired.  Both strategies share this helper;
// only the search window differs (bounded ROI vs a capped acquisition
// margin), and the fold keeps the absolute position because it folds over
// the window rather than into one symbol period.
//
// `fast_only` is the low-SNR local K=4 fast path (optimization A): it never
// probes the CFO phase or escalates to K=64, and it additionally requires a
// confirmed leading edge plus >= 4 anchor groups with a plausible anchor
// median.  A fast-only rejection clears everything it wrote and leaves the
// caller to fall through to the ORIGINAL wide/TM-1/coarse/grid/full-fold path
// unchanged.  `fast_only == false` is the original, gate-free behaviour.
bool try_fold_acquire(TimingEnv& env,
                      AcquisitionState& st,
                      int64_t centre,
                      int64_t radius,
                      bool fast_only)
{
    const auto& profile = env.profile;
    auto& out = env.out;
    auto& scratch = env.scratch;
    const auto& pf = env.pf;
    const size_t L = env.L;
    const float Et = pf.Et;
    const auto* sparse_tmpl = pf.sparse;
    const auto* shaped_tmpl = pf.shaped;
    const auto& real_tmpl = *pf.real;
    // A previous packet's shadow CFO grid must never survive into this call.
    scratch.timing_fold_cfo_samples.clear();
    scratch.timing_fold_cfo_corr.clear();
    auto acquire_fold = [&](size_t max_syms, double phase) {
        return detail::fold_acquire(
            env.rx, env.n, centre, centre - radius, centre + radius,
            scratch.work.data(), L, Et, real_tmpl,
            max_syms, profile.preamble_repetitions,
            detail::kFoldThreshold, phase,
            sparse_tmpl, &scratch.fold_buf, shaped_tmpl,
            &scratch.timing_fold_lag_metrics);
    };
    // Fast K=4 coherent fold first: immune to typical CFO (< 5 deg over
    // 4 us) and avoids probe_fold_phase overhead; falls back to K=64
    // only when K=4 + PAR misses.
    constexpr size_t kFoldFastSymbols = 4;
    double fold_phase = 0.0;
    auto fold = acquire_fold(kFoldFastSymbols, 0.0);
    if (!fold.ok && !fast_only) {
        const auto phase_probe = detail::probe_fold_phase(
            env.rx, env.n, centre, profile.preamble_repetitions);
        const double coarse_phase = phase_probe.phase();
        const bool could_alias = std::abs(coarse_phase) >=
            M_PI / detail::kFoldAnchorGroupSymbols;
        fold_phase = could_alias && phase_probe.coherent()
            ? coarse_phase : 0.0;
        fold = acquire_fold(detail::kFoldMaxSymbols, fold_phase);
        if (!fold.ok && !could_alias &&
            std::abs(coarse_phase) >= M_PI / (2 * detail::kFoldMaxSymbols) &&
            phase_probe.coherent()) {
            fold_phase = coarse_phase;
            fold = acquire_fold(detail::kFoldMaxSymbols, fold_phase);
        }
    }
    if (!fold.ok)
        return false;
    // ---- Absolute-symbol disambiguation ------------------------------
    // The preamble is 256 VERBATIM copies of one SYNC symbol, so the fold
    // metric is exactly P-periodic: `fold.start` is only the grid phase
    // nearest the seed, and the live frontend seed is accurate to a few
    // SYMBOLS (measured -6.5..+4.4 SYNC periods at -14 dB in-band), not a
    // few samples.  The one non-periodic feature available before the SFD
    // is the preamble's leading edge: immediately before the true start
    // there is noise (every acquisition window carries a pre-guard), while
    // every other candidate has a SYNC before it.  Score each candidate by
    // metric(c) - metric(c-P) and take the best.  At -14 dB in-band the
    // single-SYNC metric is ~0.038 against ~0.001 in noise.
    int64_t fold_start = fold.start;
    const auto edge = detail::disambiguate_fold_leading_edge(
        env.rx, env.n, pf.filter, L, Et, fold, radius,
        scratch.timing_fold_cand_metrics);
    if (fast_only) {
        // The local fast attempt only accepts an unambiguous leading edge;
        // `crossed` is equivalent to max_here >= 0.003f (a noise-only window
        // never produces a slot above the 0.003 floor).
        //
        // Scan-boundary guard: fold_acquire clamps the window to +/-P/2 around
        // the seed, so a seed that is DEEP INSIDE the preamble yields a strong
        // periodic fold hit at a symbol-shifted origin, and
        // disambiguate_fold_leading_edge only scans max_j = radius/P + 2
        // candidates.  In that late-seed case the very FIRST scanned candidate
        // already clears the edge threshold, so edge.start lands exactly on the
        // scan-window boundary fold.start - max_j*P -- not on the real packet
        // start.  Accepting it would skip the coarse scan + backtrack and lose
        // the packet, so a crossing at the first candidate means the true
        // leading edge was not observed inside this window: reject the fast
        // attempt and let the original wide/coarse/grid/full-K=64 path (whose
        // 40-symbol backtrack does handle late seeds) run.  This costs only
        // speed for large late seeds, never recall.
        const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
        const int64_t max_j = radius / P + 2;
        const int64_t boundary = fold.start - max_j * P;
        if (!edge.crossed || edge.start <= boundary) {
            out.peak_samples.clear();
            out.peak_metrics.clear();
            out.peak_corr.clear();
            return false;
        }
        fold_start = edge.start;
    } else if (edge.max_here >= 0.003f) {
        fold_start = edge.start;
    }
    // Period/origin from the 16-SYMBOL super-anchor fit.  A degenerate fit
    // still leaves a usable origin (the refined fold argmax) with the
    // nominal period, so the packet is not thrown away for a fit failure.
    double period_fit = static_cast<double>(kQm35SamplesPerSymbol);
    int64_t origin_fit = fold_start;
    const size_t groups = std::min<size_t>(
        64,
        profile.preamble_repetitions / detail::kFoldAnchorGroupSymbols);
    bool anchor_fit_ok = false;
    if (groups >= 2) {
        anchor_fit_ok = detail::fold_anchor_fit(
            env.rx, env.n, fold_start, scratch.work.data(), L, Et,
            real_tmpl, groups, &period_fit,
            &origin_fit, out.peak_samples, out.peak_metrics,
            out.peak_corr, fold_phase,
            sparse_tmpl, &scratch.fold_sup_buf, shaped_tmpl,
            &scratch.timing_fold_cfo_samples,
            &scratch.timing_fold_cfo_corr);
        if (!anchor_fit_ok) {
            out.peak_samples.clear();
            out.peak_metrics.clear();
            out.peak_corr.clear();
        }
    }
    if (fast_only) {
        // Strict anti-false-alarm gate (see kFoldFastAnchorMedianFloor):
        // a valid period fit, at least four measured anchors, and a median
        // anchor amplitude above the 4-symbol-folded noise floor.  Any
        // failure clears what fold_anchor_fit wrote and falls through to the
        // original acquisition path.
        if (!anchor_fit_ok || out.peak_samples.size() < 4) {
            out.peak_samples.clear();
            out.peak_metrics.clear();
            out.peak_corr.clear();
            scratch.timing_fold_cfo_samples.clear();
            scratch.timing_fold_cfo_corr.clear();
            return false;
        }
        const size_t n_peaks = out.peak_metrics.size();
        std::array<float, detail::kFoldAnchorMaxGroups> median_buf{};
        const size_t cnt = std::min(n_peaks, median_buf.size());
        std::copy_n(out.peak_metrics.begin(), cnt, median_buf.begin());
        std::nth_element(median_buf.begin(),
                         median_buf.begin() + static_cast<ptrdiff_t>(cnt / 2),
                         median_buf.begin() + static_cast<ptrdiff_t>(cnt));
        if (!(median_buf[cnt / 2] >= kFoldFastAnchorMedianFloor)) {
            out.peak_samples.clear();
            out.peak_metrics.clear();
            out.peak_corr.clear();
            scratch.timing_fold_cfo_samples.clear();
            scratch.timing_fold_cfo_corr.clear();
            return false;
        }
    } else if (out.peak_samples.empty()) {
        out.peak_samples.push_back(fold_start +
                                   static_cast<int64_t>(L) - 1);
        out.peak_metrics.push_back(fold.metric);
        out.peak_corr.push_back(fold.corr);
    }
    out.measured_period = period_fit;
    out.fold_phase_per_symbol = fold_phase;
    out.preamble_start_sample = origin_fit;
    out.metric = fold.metric;
    st.first_start = fold_start;
    st.first_metric = fold.metric;
    st.first_corr = fold.corr;
    st.first_hit_k = 0;
    st.fold_used = true;
    return true;
}

// ---------------------------------------------------------------------------
// SyncStrategy::LockedBounded acquisition
// ---------------------------------------------------------------------------
// Locked tracking contract: the target lock predicted this origin, so only a
// bounded ROI around the seed is legal.  A miss is a tracking failure -- never
// a wide scan and never a long backtrack.
bool acquire_locked_bounded(TimingEnv& env,
                            AcquisitionState& st,
                            int64_t seed_start)
{
    const auto& profile = env.profile;
    auto& out = env.out;
    const size_t L = env.L;
    const float thr = profile.radar_verification_threshold;
    ++out.bounded_attempts;
    if (seed_start < 0) {
        ++out.bounded_failures;
        return false;
    }
    const int64_t roi = std::max<int64_t>(
        1, static_cast<int64_t>(profile.locked_roi_samples));
    out.bounded_roi_samples = static_cast<size_t>(roi);
    // Probe the predicted SYNC first, then a bounded number of inner SYNCs
    // on the known symbol grid.  A weak startup SYNC must not cost the
    // packet, but the probe set stays finite: this is NOT a scan.
    const size_t inner = profile.locked_inner_symbols;
    int64_t pos = -1;
    float met = 0.0f;
    std::complex<float> corr(0.0f, 0.0f);
    size_t used_inner = 0;
    bool found = false;
    // Low-SNR early stop.  m = 0 keeps the full ROI (it is the predicted
    // SYNC, so a large residual must still be catchable).  If it misses
    // AND its raw window maximum is far below the PAR floor, the whole
    // leading region is deep-noise / deep-low-SNR, not merely mis-seeded:
    // the locked residual is a few samples, so the inner SYNCs can then be
    // probed at +/-32 instead of +/-roi (4 x 513 -> 4 x 65 correlations).
    // Anything this skips is still covered by the tail probe and, failing
    // that, by the fold -- both of which keep the full ROI.
    float m0_raw_max = 0.0f;
    for (size_t m = 0; m <= inner; ++m) {
        const int64_t want =
            static_cast<int64_t>(seed_start) +
            static_cast<int64_t>(m) * static_cast<int64_t>(kQm35SamplesPerSymbol);
        if (want + static_cast<int64_t>(L) >
            static_cast<int64_t>(env.n))
            break;
        const bool narrow_inner = (m != 0) && (m0_raw_max < kLockedInnerNarrowRawGate);
        const int64_t inner_rad =
            narrow_inner ? std::min<int64_t>(roi, 32) : roi;
        float raw_max = 0.0f;
        bool inner_hit = false;
        if (m == 0 && roi > 32) {
            // Steady-state lock residual is a few samples.  Probe +/-32
            // first (65 lags, the PAR background minimum) and only pay for
            // the full +/-roi scan when that misses.  A miss keeps the old
            // behaviour exactly, including m0_raw_max for the early-stop
            // heuristic below.
            inner_hit = local_best(env, st, want, 32, thr, &pos, &met, &corr,
                                   /*allow_par_fallback=*/true, &raw_max);
            if (!inner_hit) {
                float wide_raw = 0.0f;
                inner_hit = local_best(env, st, want, roi, thr, &pos, &met,
                                       &corr, /*allow_par_fallback=*/true,
                                       &wide_raw);
                raw_max = std::max(raw_max, wide_raw);
            }
        } else {
            inner_hit = local_best(env, st, want, inner_rad, thr, &pos, &met,
                                   &corr, /*allow_par_fallback=*/true,
                                   &raw_max);
        }
        if (inner_hit) {
            used_inner = m;
            found = true;
            break;
        }
        if (m == 0) {
            m0_raw_max = raw_max;
            out.inner_probe_raw_max = raw_max;
        }
    }
    // Failure path only.  The leading SYNCs are all below threshold --
    // typically a burst interferer parked on the first symbols while the
    // tail of the preamble is clean (measured on the 5.2 m capture:
    // first_5_max 0.029..0.053, med_tail_32_255 0.14..0.20).  The old code
    // declared a tracking failure without ever looking past
    // `locked_inner_symbols`.  Keep the bounded ROI (the grid is known, so
    // only the residual is unknown) and probe sparse tail anchors on the
    // symbol grid.
    //
    // The threshold is the TRACKING threshold, not the acquisition one: a
    // low-SNR tail sits at 0.14..0.20, well under
    // radar_verification_threshold (0.25) but far above the per-SYNC
    // tracking floor.  The neighbour gate is what keeps a +/-roi false
    // alarm from being accepted at this lower threshold.
    const float tail_thr =
        std::min(profile.radar_verification_threshold, 0.08f);
    const size_t reps = profile.preamble_repetitions;
    const size_t step =
        std::max<size_t>(1, profile.locked_tail_probe_step);
    const size_t need_confirm =
        std::max<size_t>(1, profile.locked_tail_probe_confirm);
    // Two probe widths.  The locked residual is a few samples, so the
    // narrow (+/-32) sweep finds a clean tail cheaply.  The wide sweep
    // (+/-roi) only matters when the residual is larger than the narrow
    // window -- and it costs 32 x 513 = 16,416 single-SYMBOL correlations,
    // so it is NOT run before the fold: the fold covers the identical
    // window with a 64-SYMBOL coherent gain (~18 dB more sensitive) and is
    // the correct rescue in that regime.  The wide sweep is therefore the
    // last resort, reached only when the fold itself could not rescue the
    // packet.  Measured at -12 dB in-band (200 packets): the narrow sweep
    // never produced a confirmable hit and the fold rescued every packet,
    // so this ordering removes ~1.2 ms/packet of doomed work without
    // changing any outcome.
    auto tail_sweep = [&](int64_t trad) {
        for (size_t m = inner + 1; m < reps; m += step) {
            const int64_t want = static_cast<int64_t>(seed_start) +
                static_cast<int64_t>(m) *
                    static_cast<int64_t>(kQm35SamplesPerSymbol);
            if (want + static_cast<int64_t>(L) >
                static_cast<int64_t>(env.n))
                break;
            ++out.tail_probe_attempts;
            int64_t cpos = -1;
            float cmet = 0.0f;
            std::complex<float> ccorr(0.0f, 0.0f);
            float raw_max = 0.0f;
            const bool hit =
                local_best(env, st, want, trad, tail_thr, &cpos, &cmet,
                           &ccorr, /*allow_par_fallback=*/true, &raw_max);
            if (raw_max > out.tail_probe_raw_max)
                out.tail_probe_raw_max = raw_max;
            if (!hit)
                continue;
            if (cmet > out.tail_probe_best_metric)
                out.tail_probe_best_metric = cmet;
            // Consistency gate: a real grid hit is confirmed by its
            // neighbouring SYMBOLS.  A false alarm (wrong lag inside the
            // ROI) is not, so probing continues instead of locking onto a
            // bogus origin.
            size_t confirm = 0;
            const int64_t nb_delta[3] = {1, 2, -1};
            for (int64_t d : nb_delta) {
                const int64_t nb = static_cast<int64_t>(m) + d;
                if (nb < 0 || nb >= static_cast<int64_t>(reps))
                    continue;
                const int64_t nb_want =
                    static_cast<int64_t>(seed_start) +
                    nb * static_cast<int64_t>(kQm35SamplesPerSymbol);
                if (nb_want + static_cast<int64_t>(L) >
                    static_cast<int64_t>(env.n))
                    continue;
                int64_t np = -1;
                float nm = 0.0f;
                std::complex<float> nc(0.0f, 0.0f);
                float nb_raw = 0.0f;
                const bool nb_hit =
                    local_best(env, st, nb_want, trad, tail_thr, &np, &nm,
                               &nc, /*allow_par_fallback=*/false, &nb_raw);
                if (nb_raw > out.tail_probe_raw_max)
                    out.tail_probe_raw_max = nb_raw;
                if (nb_hit)
                    ++confirm;
            }
            if (confirm < need_confirm)
                continue;
            ++out.tail_probe_hits;
            out.tail_probe_k = static_cast<int64_t>(m);
            pos = cpos;
            met = cmet;
            corr = ccorr;
            used_inner = m;
            found = true;
            return;
        }
    };
    if (!found && profile.locked_tail_probe)
        tail_sweep(std::min<int64_t>(roi, 32));
    if (!found) {
        // Low-SNR fallback: every single-symbol probe inside the bounded
        // ROI missed.  The fold keeps the same bounded window -- it is not
        // an escalation to a wide scan.
        //
        // Gate: LockedBounded's failure-path rescue is a single opt-out
        // tier.  `locked_tail_probe == false` is how a caller asks for the
        // original hard-fail tracking contract (no rescue at all), so the
        // fold honours the same switch instead of silently replacing the
        // tail probe with a different rescue.  qa_uwb_demod_core's
        // "tail probe rescues leading interference" case pins that
        // contract: with rescue disabled the packet MUST stay lost.
        const bool fold_ok =
            profile.locked_tail_probe &&
            try_fold_acquire(env, st, static_cast<int64_t>(seed_start), roi,
                             /*fast_only=*/false);
        if (!fold_ok && profile.locked_tail_probe) {
            // Last resort (see the note above): the fold could not rescue
            // this packet either, so the wide single-SYMBOL sweep is the
            // only remaining option.  Paying its cost here cannot slow down
            // any packet the fold already handles.
            tail_sweep(roi);
        }
        if (!fold_ok && !found) {
            ++out.bounded_failures;
            return false;
        }
    }
    if (!st.fold_used) {
        // Convert the inner-SYNC hit back to the packet start on the grid.
        st.first_start = pos - static_cast<int64_t>(used_inner) *
                                static_cast<int64_t>(kQm35SamplesPerSymbol);
        st.first_hit_k = used_inner;
        st.first_metric = met;
        st.first_corr = corr;
        if (used_inner != 0)
            ++out.bounded_inner_used;
    }
    ++out.bounded_success;
    return true;
}

// ---------------------------------------------------------------------------
// SyncStrategy::AcquireLegacy acquisition
// ---------------------------------------------------------------------------
bool acquire_unbounded_seed(TimingEnv& env,
                            AcquisitionState& st,
                            int64_t seed_start)
{
    const auto& profile = env.profile;
    auto& out = env.out;
    auto& scratch = env.scratch;
    const auto& pf = env.pf;
    const size_t L = env.L;
    const float Et = pf.Et;
    const float thr = profile.radar_verification_threshold;
    const int64_t seed =
        (seed_start >= 0) ? seed_start
                          : (int64_t)(env.roi_start +
                                      (env.roi_end - env.roi_start) / 2);
    st.acq_seed = seed;
    const int64_t margin =
        static_cast<int64_t>(profile.timing_search_margin);
    const int64_t stride = std::max<int64_t>(
        1, static_cast<int64_t>(profile.timing_coarse_stride));
    const int64_t j_lo = std::max<int64_t>(0, seed - margin);
    const int64_t j_hi =
        std::min<int64_t>(static_cast<int64_t>(env.n) - static_cast<int64_t>(L),
                          seed + margin);

    auto consider = [&](int64_t pos, float met, std::complex<float> corr) {
        if (met >= thr && met > st.first_metric) {
            st.first_metric = met;
            st.first_start = pos;
            st.first_corr = corr;
        }
    };

    // MATLAB detectSeededPreamble: prefer a narrow FULL-RATE search around
    // the seed (max(256, ceil(symbolLength/32)) samples each side).  A wide
    // stride-S coarse scan can miss a narrow symbol-0 peak (sub-symbol seed
    // offset) and then prefer a later STRONGER SYNC -- exactly what made
    // the DW1000 preamble_start land ~6 SYNCs late.  The seed-nearby
    // result is authoritative; the wide scan is only a fallback.
    bool seed_acquired = false;
    const int64_t legacy_seed_radius =
        std::max<int64_t>(256, stride);
    if (profile.timing_seed_refine_radius != 0) {
        const int64_t narrow = rxopt::seed_refine_radius(
            profile.timing_seed_refine_radius, stride);
        if (narrow < legacy_seed_radius) {
            const auto narrow_begin = std::chrono::steady_clock::now();
            int64_t pos = -1;
            float met = 0.0f;
            std::complex<float> corr(0.0f, 0.0f);
            // PAR fallback is enabled on the two NARROW local searches
            // only.  A 256-SYNC preamble at -4..-8 dB in-band SNR has
            // every SYNC peak above 0.08 with PAR 232..253, but only
            // 11-24 of 256 clear the 0.25 absolute threshold, so the
            // absolute gate alone loses the packet.  The wide stride scan
            // and the whole-SYNC grid probes stay absolute-only: at
            // thousands of lags a pure-noise false alarm would clear any
            // PAR gate.
            const bool narrow_ok =
                local_best(env, st, seed, narrow, thr, &pos, &met, &corr,
                           /*allow_par_fallback=*/true);
            out.seed_search_us = static_cast<uint64_t>(
                std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::steady_clock::now() - narrow_begin).count());
            ++out.seed_refine_narrow_attempts;
            if (narrow_ok) {
                st.first_start = pos;
                st.first_metric = met;
                st.first_corr = corr;
                seed_acquired = true;
                ++out.seed_refine_narrow_success;
            } else {
                ++out.seed_refine_fallbacks;
            }
        }
    }
    if (!seed_acquired) {
        // Legacy first seeded refine.  Also the fallback after a failed
        // narrow search: a recovered packet must not be lost just because
        // the local peak was outside the narrow radius.
        int64_t pos = -1;
        float met = 0.0f;
        std::complex<float> corr(0.0f, 0.0f);
        if (local_best(env, st, seed, legacy_seed_radius, thr, &pos, &met,
                       &corr, /*allow_par_fallback=*/true)) {
            st.first_start = pos;
            st.first_metric = met;
            st.first_corr = corr;
            seed_acquired = true;
        }
    }

    // ---- Low-SNR local K=4 fold fast path (optimization A) --------------
    // The single-symbol paths above miss at -12/-14 dB in-band, but the live
    // frontend seed is still good to a few SYMBOLS.  A local K=4 coherent fold
    // over seed +/- min(margin, 8192) buys +6 dB and costs one 1016-lag
    // correlation (~15-20 us); when it is confirmed by the leading-edge
    // disambiguation and the sparse anchor fit, the entire TM-1 / +/-131072
    // coarse scan / whole-SYNC grid / K=64 fold tail is skipped for this
    // packet.  Any rejection falls through to that original path untouched.
    // No seed (blind search) or the wide interval probe (SIC fast-fail) keep
    // their existing behaviour and never enter here.
    if (!seed_acquired &&
        profile.sync_strategy == SyncStrategy::AcquireLegacy &&
        seed_start >= 0 && !profile.rx_wide_interval_probe) {
        const int64_t fast_radius = std::min<int64_t>(margin, 8192);
        if (try_fold_acquire(env, st, seed_start, fast_radius,
                             /*fast_only=*/true))
            return true;
    }

    // ---- Wide-window interval probe (K-symbol local fold) --------------
    // Reached only for AcquireLegacy when the blind search spans many
    // SYMBOL periods.  Each probe coherently folds K SYNC periods
    // (+3 dB at K=2, +6 dB at K=4) and correlates over one period.
    // The coarse probe grid only skips *start* positions, never a whole
    // preamble: a probe
    // anywhere inside the preamble still sees the same P-periodic code
    // (its peak lag resolves back to a true symbol boundary), so the
    // normal backtrack / SYNC-train / SFD path recovers the origin.
    //
    // On an all-miss, the historical stride-14 coarse scan, whole-SYNC
    // grid and fold rescue are skipped so an SIC candidate that contains
    // no target preamble fails fast to L0_Bypass instead of paying the
    // full wide acquisition.  rx_wide_probe_fallback restores the old
    // behaviour for callers that prefer recall over latency.
    if (!seed_acquired && profile.rx_wide_interval_probe) {
        const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
        if (P > 0 && j_hi - j_lo >= 16 * P) {
            const size_t fold_k = std::max<size_t>(
                2, std::min<size_t>(profile.rx_wide_probe_fold_symbols, 4));
            size_t step_symbols = profile.rx_wide_probe_step_symbols
                ? profile.rx_wide_probe_step_symbols : 42;
            const size_t m_cap = std::max<size_t>(
                8, profile.preamble_repetitions / 3);
            step_symbols = std::min(step_symbols, m_cap);
            if (profile.preamble_repetitions > fold_k)
                step_symbols = std::min(
                    step_symbols, profile.preamble_repetitions - fold_k);
            else
                step_symbols = 0;
            // The probe may land up to `step_symbols` into the preamble;
            // the backtrack must be able to walk that far back.
            if (profile.timing_max_backtrack_symbols > 0)
                step_symbols = std::min(step_symbols,
                                        profile.timing_max_backtrack_symbols);
            const int64_t last_probe =
                std::min(j_hi,
                         static_cast<int64_t>(env.n) -
                             static_cast<int64_t>(fold_k) * P -
                             static_cast<int64_t>(L) + 1);
            if (step_symbols > 0 && last_probe >= j_lo) {
                const int64_t step =
                    static_cast<int64_t>(step_symbols) * P;
                int64_t best_pos = -1;
                float best_met = -1.0f;
                std::complex<float> best_corr(0.0f, 0.0f);
                int64_t last_scanned = -1;
                auto probe_one = [&](int64_t p) {
                    if (p == last_scanned)
                        return;
                    last_scanned = p;
                    const auto r = detail::fold_probe_acquire(
                        env.rx, env.n, p, scratch.work.data(), L, Et,
                        *pf.real, fold_k, thr,
                        pf.sparse, &scratch.fold_buf, pf.shaped);
                    if (r.ok && r.start >= j_lo && r.start <= j_hi &&
                        r.metric > best_met) {
                        best_met = r.metric;
                        best_pos = r.start;
                        best_corr = r.corr;
                    }
                };
                for (int64_t p = j_lo; p <= last_probe; p += step)
                    probe_one(p);
                probe_one(last_probe);
                if (best_pos >= 0) {
                    st.first_start = best_pos;
                    st.first_metric = best_met;
                    st.first_corr = best_corr;
                    seed_acquired = true;
                } else if (!profile.rx_wide_probe_fallback) {
                    return false;
                }
            }
        }
    }

    // ---- Phase probe + symbol-grid acquisition (TM-1) ------------------
    // The preamble is N verbatim copies of one P=1016-sample SYNC symbol.
    // A window of P consecutive lags therefore contains one code-phase
    // match, so the exact sample phase can be recovered from a handful of
    // probe windows instead of correlating all 2*margin lags.  A strict
    // threshold keeps this zero-risk: if the probe does not see a solid
    // SYNC (low SNR, packet outside the probes), everything below falls
    // through to the original coarse scan unchanged.
    //
    // Cheap pre-gate: the frontend seed is accurate to a few SYMBOLS, so
    // the code phase at the seed is preserved and one single-lag dot tells
    // us whether this window is deep noise.  The 0.2 floor is close to the
    // tracking threshold on purpose: at -12/-14 dB the single-SYNC metric
    // is ~0.04-0.1, so the three P-length probes (which would all fail)
    // are skipped and the fold path is not slowed down by them.
    float probe_gate = -1.0f;
    {
        const int64_t gp = std::min(std::max(seed, j_lo), j_hi);
        if (gp >= 0 && gp + static_cast<int64_t>(L) <=
                           static_cast<int64_t>(env.n)) {
            const auto acc = pf.filter.dot(env.rx + gp);
            probe_gate = detail::normalized_correlation_power(
                acc, window_power_8acc(env.rx + gp, L), Et);
        }
    }
    if (probe_gate >= 0.2f && !seed_acquired && j_hi > j_lo &&
        j_hi - j_lo >= 2 * static_cast<int64_t>(kQm35SamplesPerSymbol)) {
        const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
        int64_t probes[3];
        size_t nprobe = 0;
        auto push_probe = [&](int64_t p) {
            if (p < j_lo || p > j_hi)
                return;
            if (p + static_cast<int64_t>(L) + 64 > static_cast<int64_t>(env.n))
                return;
            if (nprobe < 3)
                probes[nprobe++] = p;
        };
        push_probe(std::min(std::max(seed, j_lo), j_hi));
        if (margin >= 32 * P) {
            push_probe(seed - 32 * P);
            push_probe(seed + 32 * P);
        }
        float best_probe = -1.0f;
        int64_t anchor = -1;
        for (size_t i = 0; i < nprobe; ++i) {
            const int64_t p = probes[i];
            const int64_t avail =
                static_cast<int64_t>(env.n) - static_cast<int64_t>(L) - p + 1;
            const size_t lags = static_cast<size_t>(
                std::min<int64_t>(P, std::max<int64_t>(0, avail)));
            if (lags < 64)
                continue;
            const auto pk = rxopt::scan_normalized_sync(
                pf.filter, env.rx + p, lags,
                window_power_8acc(env.rx + p, L), Et);
            if (pk.metric > best_probe) {
                best_probe = pk.metric;
                anchor = p + static_cast<int64_t>(pk.lag);
            }
        }
        if (anchor >= 0 && best_probe >= thr) {
            // Same phase modulo P for every SYNC, so only the symbol grid
            // remains: ~2*margin/P grid points at +/-2 samples each.
            const int64_t rem = ((anchor - j_lo) % P + P) % P;
            for (int64_t want = j_lo + rem; want <= j_hi; want += P) {
                int64_t pos = -1;
                float met = 0.0f;
                std::complex<float> corr(0.0f, 0.0f);
                if (local_best(env, st, want, 2, thr, &pos, &met, &corr))
                    consider(pos, met, corr);
            }
            if (st.first_start >= 0)
                seed_acquired = true;
        }
    }

    // (a) Coarse scan (fallback only) -- batched SIMD when stride is 1 or 2,
    // sliding power conjugate dots at stride otherwise.
    const auto coarse_begin = std::chrono::steady_clock::now();
    if (!seed_acquired && j_hi > j_lo) {
        float best_coarse_m = -1.0f;
        int64_t best_coarse_j = j_lo;
        std::complex<float> best_coarse_c(0.0f, 0.0f);
        if (stride == 1) {
            const size_t lag_count = static_cast<size_t>(j_hi - j_lo + 1);
            const auto pk = rxopt::scan_normalized_sync(
                pf.filter, env.rx + j_lo, lag_count,
                window_power_8acc(env.rx + j_lo, L), Et);
            best_coarse_j = j_lo + static_cast<int64_t>(pk.lag);
            best_coarse_m = pk.metric;
            best_coarse_c = pk.corr;
        } else if (stride == 2 && pf.sparse && !pf.shaped && (L % 2 == 0)) {
            const size_t half_lags =
                static_cast<size_t>((j_hi - j_lo) / 2 + 1);
            const size_t half_len = L / 2;
            const size_t even_count = (half_lags - 1) + half_len;
            scratch.corr.resize(even_count);
            rxopt::pack_even_complex(env.rx + j_lo, even_count,
                                     scratch.corr.data());
            rxopt::ensure_half_rate_sparse_template(
                *pf.sparse, scratch.half_cir_template);
            const rxopt::SyncMatchedFilter half_filter(
                nullptr, half_len, nullptr, &scratch.half_cir_template);
            const auto pk = rxopt::scan_normalized_sync(
                half_filter, scratch.corr.data(), half_lags,
                window_power_8acc(scratch.corr.data(), half_len), Et);
            best_coarse_j = j_lo + static_cast<int64_t>(pk.lag) * 2;
            best_coarse_c = pk.corr;
            best_coarse_m = detail::normalized_correlation_power(
                pk.corr, window_power_8acc(env.rx + best_coarse_j, L), Et);
        } else {
            float pw = window_power_8acc(env.rx + j_lo, L);
            // Slide the window power forward by `cnt` samples with eight
            // independent accumulators, so the two |x|^2 streams vectorise and
            // only ONE serial float add per stride block remains.  The
            // per-sample form of this update was the single most expensive loop
            // in low-SNR acquisition: 262,144 iterations (the live search
            // margin is +-131,072 samples) at ~8 cycles each = ~2.2 ms/packet,
            // 85% of the whole acquisition search.  The sum is mathematically
            // identical; only the float association order changes.
            auto advance_power = [&](int64_t j, int64_t cnt) {
                const float* a =
                    reinterpret_cast<const float*>(env.rx + j + static_cast<int64_t>(L));
                const float* b = reinterpret_cast<const float*>(env.rx + j);
                const size_t m = static_cast<size_t>(2 * cnt);
                float s0 = 0.f, s1 = 0.f, s2 = 0.f, s3 = 0.f;
                float s4 = 0.f, s5 = 0.f, s6 = 0.f, s7 = 0.f;
                size_t i = 0;
                for (; i + 8 <= m; i += 8) {
                    s0 += a[i + 0] * a[i + 0] - b[i + 0] * b[i + 0];
                    s1 += a[i + 1] * a[i + 1] - b[i + 1] * b[i + 1];
                    s2 += a[i + 2] * a[i + 2] - b[i + 2] * b[i + 2];
                    s3 += a[i + 3] * a[i + 3] - b[i + 3] * b[i + 3];
                    s4 += a[i + 4] * a[i + 4] - b[i + 4] * b[i + 4];
                    s5 += a[i + 5] * a[i + 5] - b[i + 5] * b[i + 5];
                    s6 += a[i + 6] * a[i + 6] - b[i + 6] * b[i + 6];
                    s7 += a[i + 7] * a[i + 7] - b[i + 7] * b[i + 7];
                }
                float s = ((s0 + s1) + (s2 + s3)) + ((s4 + s5) + (s6 + s7));
                for (; i < m; ++i)
                    s += a[i] * a[i] - b[i] * b[i];
                pw += s;
            };
            for (int64_t j = j_lo; j <= j_hi; j += stride) {
                const auto acc = pf.filter.dot(env.rx + j);
                const float m = detail::normalized_correlation_power(acc, pw, Et);
                if (m > best_coarse_m) {
                    best_coarse_m = m;
                    best_coarse_j = j;
                    best_coarse_c = acc;
                }
                const int64_t cnt = std::min<int64_t>(stride, j_hi + 1 - j);
                // Only slide when another search lag remains: the final
                // advance would read rx[j + L] one complex sample past the
                // n-sample window (and its value is never consumed).
                if (j + cnt <= j_hi)
                    advance_power(j, cnt);
            }
        }
        out.coarse_search_us = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - coarse_begin).count());
        // (b) Fine refine around the coarse peak.
        int64_t pos = -1;
        float met = 0.0f;
        std::complex<float> corr(0.0f, 0.0f);
        if (local_best(env, st, best_coarse_j, stride, thr, &pos, &met, &corr))
            consider(pos, met, corr);
        else if (best_coarse_m >= thr)
            consider(best_coarse_j, best_coarse_m, best_coarse_c);
    }

    // (c) Whole-SYNC grid probes (fallback only; covers a seed that lands
    // mid-preamble and a sub-symbol confirmation that must not replace an
    // SFD-refined anchor with a stronger neighbouring SYNC).
    if (!seed_acquired) {
        int64_t pos = -1;
        float met = 0.0f;
        std::complex<float> corr(0.0f, 0.0f);
        const int64_t grid_half = margin / kQm35SamplesPerSymbol;
        for (int64_t n_sym = -grid_half; n_sym <= grid_half; ++n_sym) {
            if (local_best(env, st, seed + n_sym * static_cast<int64_t>(
                                           kQm35SamplesPerSymbol),
                           8, thr, &pos, &met, &corr))
                consider(pos, met, corr);
        }
    }
    if (st.first_start < 0) {
        // Low-SNR fallback: narrow seed refine, legacy seed refine, stride
        // coarse scan and whole-SYNC grid probes all missed.  Fold over a
        // capped acquisition window (never the full 40k margin: the fold
        // cost grows with the window and the seed is already close).
        if (!try_fold_acquire(env, st, seed, std::min<int64_t>(margin, 8192),
                              /*fast_only=*/false))
            return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// LockedBounded grid tracking: sparse anchors + single point
// ---------------------------------------------------------------------------
// The locked contract already knows the symbol grid (the target lock
// predicted this origin), so the only unknowns are the residual origin and
// the residual period.  Measure a sparse set of SYNC anchors with the
// normal +/-track_radius search, fit period+origin from the anchors, then
// evaluate every other SYNC with ONE correlation at start+round(k*period).
// The peak train keeps the same meaning (CFO needs peak_corr phases, CIR
// needs the fitted period), at ~1/17 of the per-symbol search cost.
// measured_period MUST come from the anchor fit: it feeds
// integer_chip_grid()/spc_ppm and the fused-CMF admission, so a period
// that degenerates to the nominal 1016 would silently disable the
// fractional-grid fallback for an untrimmed clock.
bool track_locked_grid(TimingEnv& env,
                       AcquisitionState& st,
                       int64_t start,
                       int64_t period,
                       double* grid_period,
                       int64_t* grid_start)
{
    const auto& profile = env.profile;
    auto& out = env.out;
    auto& scratch = env.scratch;
    const size_t L = env.L;
    const float Et = env.pf.Et;
    const size_t reps = profile.preamble_repetitions;
    const size_t first_hit_k = st.first_hit_k;
    const float first_metric = st.first_metric;
    const std::complex<float> first_corr = st.first_corr;
    if (reps < 8)
        return false;
    const float gthr =
        std::min(profile.radar_verification_threshold, 0.08f);
    const int64_t arad = std::max<int64_t>(
        8, static_cast<int64_t>(profile.timing_track_radius));
    const double nominal = static_cast<double>(period);

    // Sparse anchor set on the symbol grid: the confirmed first SYNC, the
    // last SYNC of the train, and evenly spaced points in between.
    std::array<size_t, 32> ak{};
    size_t na = 0;
    ak[na++] = first_hit_k;
    const size_t last_k = reps - 1;
    if (last_k > first_hit_k) {
        const size_t want_anchors = std::min<size_t>(
            std::max<size_t>(4, profile.locked_grid_anchor_count), ak.size());
        const size_t intervals = (want_anchors > 1) ? (want_anchors - 1) : 1;
        const double span = static_cast<double>(last_k - first_hit_k);
        for (size_t i = 1; i <= intervals && na < ak.size(); ++i) {
            const size_t k = first_hit_k + static_cast<size_t>(std::llround(
                span * static_cast<double>(i) /
                static_cast<double>(intervals)));
            if (k != ak[na - 1])
                ak[na++] = k;
        }
    }

    // One least-squares pass over the anchors: peak(k) = start + k*period,
    // in the same symbol-END convention as the legacy fit.
    auto fit_anchors = [&](double period_hint, int64_t start_hint,
                           int64_t radius, double* period_out,
                           int64_t* start_out, float* median_out,
                           size_t* hits_out) -> bool {
        detail::LinearFit fit;
        std::array<float, 32> mets{};
        size_t hits = 0;
        size_t kmin = 0, kmax = 0;
        for (size_t i = 0; i < na; ++i) {
            const size_t k = ak[i];
            const int64_t want = start_hint + static_cast<int64_t>(
                std::llround(static_cast<double>(k) * period_hint));
            if (want < 0 ||
                want + static_cast<int64_t>(L) > static_cast<int64_t>(env.n))
                continue;
            ++out.grid_anchor_attempts;
            int64_t pos = -1;
            float met = 0.0f;
            std::complex<float> c(0.0f, 0.0f);
            if (!local_best(env, st, want, radius, gthr, &pos, &met, &c))
                continue;
            const double x = static_cast<double>(k);
            const double y = static_cast<double>(pos + (int64_t)(L - 1));
            fit.add(x, y);
            if (hits == 0) {
                kmin = k;
                kmax = k;
            } else {
                kmin = std::min(kmin, k);
                kmax = std::max(kmax, k);
            }
            mets[hits++] = met;
        }
        *hits_out = hits;
        if (hits < std::max<size_t>(2, profile.locked_grid_min_anchors))
            return false;
        if (kmax - kmin < 8)
            return false; // anchors too bunched to measure a slope
        const double denom = fit.determinant();
        if (!(std::abs(denom) > 1e-9))
            return false;
        const double p = fit.slope(denom);
        if (!std::isfinite(p) || !(p > 0.0))
            return false;
        const double intercept_end = fit.intercept(p);
        *period_out = p;
        *start_out = std::max<int64_t>(
            0, static_cast<int64_t>(std::llround(intercept_end)) -
                   static_cast<int64_t>(L - 1));
        std::sort(mets.begin(), mets.begin() + hits);
        *median_out = mets[hits / 2];
        return true;
    };

    double p1 = nominal;
    int64_t s1 = start;
    float med1 = 0.0f;
    size_t hits1 = 0;
    if (!fit_anchors(nominal, start, arad, &p1, &s1, &med1, &hits1))
        return false;
    // A fitted period outside 0.5% of the nominal work-grid SYNC period
    // means the anchors did not land on the expected symbol indices.
    if (std::abs(p1 - nominal) > 0.005 * nominal)
        return false;
    double p2 = p1;
    int64_t s2 = s1;
    float med2 = med1;
    size_t hits2 = hits1;
    // Second pass only when the nominal grid would put the last anchor
    // more than 4 samples off (an untrimmed clock).  With a trimmed
    // crystal the first pass is already unbiased and this costs nothing.
    if (std::abs(p1 - nominal) * static_cast<double>(last_k) > 4.0) {
        if (!fit_anchors(p1, s1, std::max<int64_t>(arad, 16), &p2, &s2,
                         &med2, &hits2))
            return false;
        if (std::abs(p2 - nominal) > 0.005 * nominal)
            return false;
    }
    out.grid_anchor_hits = static_cast<uint32_t>(hits2);

    // Single-point stage: evaluate eight consecutive lags around the fitted
    // grid position and keep the best.  The 8-lag kernel has independent
    // accumulators and is far cheaper per lag than the single-lag kernel
    // (which is one long FMA dependency chain), so eight lags cost barely
    // more than one -- and taking the argmax keeps the peak positions and
    // the CFO phases identical to the legacy per-SYNC search instead of
    // pinning them to the fitted grid.
    scratch.metric.assign(reps, 0.0f);
    size_t valid = 0;
    // Least-squares accumulators over every refined peak (see the final
    // fit note after the loop).
    detail::LinearFit final_fit;
    size_t fn = 0;
    for (size_t k = 0; k < reps; ++k) {
        if (k == first_hit_k) {
            // The bounded probe already measured this SYNC with a real
            // +/-roi search; keep that measurement instead of
            // re-projecting it onto the fitted grid.
            scratch.metric[k] = first_metric;
            ++out.grid_single_points;
            if (first_metric >= gthr) {
                ++valid;
                out.peak_samples.push_back(start + static_cast<int64_t>(L - 1));
                out.peak_metrics.push_back(first_metric);
                out.peak_corr.push_back(first_corr);
                const double x = static_cast<double>(k);
                const double y = static_cast<double>(start + (int64_t)(L - 1));
                final_fit.add(x, y);
                ++fn;
            }
            continue;
        }
        const int64_t centre = s2 + static_cast<int64_t>(
            std::llround(static_cast<double>(k) * p2));
        int64_t base = centre - 3;
        int nl = 8;
        if (base + 7 + static_cast<int64_t>(L) > static_cast<int64_t>(env.n)) {
            // Window edge: fall back to the single lag at the prediction
            // so the last SYNCs are still measured.
            if (centre + static_cast<int64_t>(L) >
                static_cast<int64_t>(env.n))
                continue;
            base = centre;
            nl = 1;
        }
        if (base < 0)
            continue;
        float pw = window_power_8acc(env.rx + base, L);
        rxopt::MatchedFilterPeak peak;
        env.pf.filter.for_each_lag(env.rx + base, static_cast<size_t>(nl),
                                   [&](size_t lag, std::complex<float> corr, bool) {
            if (lag != 0)
                pw += std::norm(env.rx[base + lag - 1 + L]) -
                      std::norm(env.rx[base + lag - 1]);
            const float metric = detail::normalized_correlation_power(corr, pw, Et);
            if (lag == 0 || metric > peak.metric)
                peak = { lag, metric, corr };
        });
        const int64_t want = base + static_cast<int64_t>(peak.lag);
        const float m = peak.metric;
        scratch.metric[k] = m;
        ++out.grid_single_points;
        if (m < gthr)
            continue;
        ++valid;
        out.peak_samples.push_back(want + static_cast<int64_t>(L - 1));
        out.peak_metrics.push_back(m);
        out.peak_corr.push_back(peak.corr);
        const double x = static_cast<double>(k);
        const double y = static_cast<double>(want + (int64_t)(L - 1));
        final_fit.add(x, y);
        ++fn;
    }
    out.grid_single_valid = static_cast<uint32_t>(valid);
    if (valid < 4)
        return false;
    // Grid sanity: if the fitted grid does not land on real peaks, the
    // single-point metrics collapse relative to the anchor peaks.  Use the
    // 90th percentile so a weak (interfered) part of the preamble does not
    // by itself reject a correct grid.
    if (reps <= 512) {
        const size_t idx = static_cast<size_t>(
            0.9 * static_cast<double>(reps - 1));
        std::nth_element(scratch.metric.begin(),
                         scratch.metric.begin() + idx,
                         scratch.metric.end());
        if (!(scratch.metric[idx] >= 0.5f * med2))
            return false;
    }
    // Final period/origin: least squares over EVERY refined peak, which is
    // exactly the estimator the legacy 256-point fit uses.  The anchors are
    // the bootstrap that puts the single-point windows in the right place
    // (and the fallback if this fit is degenerate); they are not the final
    // estimate, because a 16-anchor fit is ~4x noisier than the legacy
    // 256-point fit and that noise measurably flips marginal payloads
    // (measured: 1 packet / 2000 at 2.5 m).  The refined positions are real
    // argmax measurements, NOT the quantised grid, so measured_period stays
    // a genuine measurement and integer_chip_grid()/spc_ppm keeps working.
    if (fn >= 8) {
        const double denom = final_fit.determinant();
        if (std::abs(denom) > 1e-9) {
            const double pr = final_fit.slope(denom);
            const double intercept_end = final_fit.intercept(pr);
            if (std::isfinite(pr) && pr > 0.0 &&
                std::abs(pr - p2) <= 0.01) { // 10 ppm sanity vs anchors
                p2 = pr;
                s2 = std::max<int64_t>(
                    0, static_cast<int64_t>(std::llround(intercept_end)) -
                           static_cast<int64_t>(L - 1));
            }
        }
    }
    *grid_period = p2;
    *grid_start = s2;
    return true;
}

// ---------------------------------------------------------------------------
// Legacy per-SYMBOL SYNC train
// ---------------------------------------------------------------------------
// Track all SYNC starts: expected = start + k*period, then convert each to
// symbol-END convention (start + L - 1) to match the detector / MATLAB.
// When bounded inner-SYNC probing hit symbol `first_hit_k > 0`, anchor the
// first peak at its true symbol coordinate and continue from `first_hit_k + 1`
// rather than restarting at k=1 (which already failed).  Also allow up to 2
// consecutive missed SYNCs so a brief 1-2 symbol interference dip does not
// prematurely truncate the 256-symbol train.
void track_sync_train(TimingEnv& env,
                      AcquisitionState& st,
                      int64_t start,
                      int64_t period,
                      size_t backtracked_symbols)
{
    const auto& profile = env.profile;
    auto& out = env.out;
    auto& scratch = env.scratch;
    const size_t L = env.L;
    const bool bounded = profile.sync_strategy == SyncStrategy::LockedBounded;
    std::vector<size_t>& peak_indices = scratch.timing_peak_indices;
    constexpr size_t kMaxConsecutiveSyncMisses = 2;
    // One tracking pass at a given floor.  `weak` selects the 0.08 floor
    // LockedBounded already uses for its tail probe; the legacy absolute
    // threshold is used otherwise.
    auto track_train = [&](bool weak) {
        out.peak_samples.clear();
        out.peak_metrics.clear();
        out.peak_corr.clear();
        peak_indices.clear();
        const int64_t first_hit_start =
            start + static_cast<int64_t>(st.first_hit_k) * period;
        out.peak_samples.push_back(first_hit_start +
                                   static_cast<int64_t>(L - 1));
        out.peak_metrics.push_back(st.first_metric);
        out.peak_corr.push_back(st.first_corr);
        peak_indices.push_back(st.first_hit_k);
        size_t consecutive_misses = 0;
        for (size_t k = st.first_hit_k + 1;
             k < profile.preamble_repetitions; ++k) {
            const int64_t want = start + static_cast<int64_t>(k * period);
            int64_t pk = -1;
            float pm = 0.0f;
            std::complex<float> pk_corr(0.0f, 0.0f);
            const float tracking_threshold =
                (bounded || weak)
                    ? std::min(profile.radar_verification_threshold,
                               bounded ? 0.08f : 0.06f)
                    : ((k < backtracked_symbols)
                           ? 0.7f * profile.radar_verification_threshold
                           : profile.radar_verification_threshold);
            const int64_t track_rad =
                (bounded && consecutive_misses > 0)
                    ? std::max<int64_t>(
                          12, static_cast<int64_t>(
                                  profile.timing_track_radius))
                    : static_cast<int64_t>(profile.timing_track_radius);
            if (!local_best(env, st, want, track_rad, tracking_threshold, &pk,
                            &pm, &pk_corr)) {
                ++consecutive_misses;
                if (!bounded &&
                    consecutive_misses > kMaxConsecutiveSyncMisses)
                    break; // lost the train
                continue;
            }
            consecutive_misses = 0;
            out.peak_samples.push_back(pk + static_cast<int64_t>(L - 1));
            out.peak_metrics.push_back(pm);
            out.peak_corr.push_back(pk_corr);
            peak_indices.push_back(k);
        }
        return out.peak_samples.size();
    };
    // First pass at the legacy threshold: high-SNR packets are complete
    // here and pay nothing.  Only a train that the absolute threshold
    // truncated (a low-SNR acquisition whose SYNCs all sit at the same
    // sub-threshold level) is retried at the 0.08 floor.  Without this the
    // packet is lost one stage later anyway (fit_cfo needs >= 4 peaks and
    // stage_cir_softchips needs more than cir_skip_initial_repetitions).
    const size_t min_peaks =
        std::max<size_t>(4, profile.min_valid_peaks);
    size_t tracked = track_train(false);
    if (!bounded && tracked < min_peaks)
        tracked = track_train(true);
    (void)tracked;
}

} // namespace

bool stage_timing(const std::complex<float>* rx,
                  size_t n,
                  const Qm35825Profile& profile,
                  const std::vector<std::complex<float>>& template_wf,
                  int64_t seed_start, // absolute predicted start, or -1
                  TimingResult& out,
                  DemodScratch& scratch)
{
    out = TimingResult{};
    const auto timing_begin = std::chrono::steady_clock::now();
    auto finish_timing_profile = [&]() {
        const uint64_t total = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - timing_begin).count());
        out.fine_track_us =
            (total >= out.coarse_search_us) ? total - out.coarse_search_us : 0;
    };
    const size_t L = template_wf.size();
    if (L == 0 || n < L) {
        finish_timing_profile();
        return false;
    }

    TimingEnv env{ rx, n, L, profile, template_wf, scratch, out,
                   prepare_timing_filter(profile, template_wf, scratch) };

    // Determine the ROI around the seed (absolute coordinates).
    env.roi_start = 0;
    env.roi_end = static_cast<int64_t>(n);
    if (seed_start >= 0) {
        const int64_t seed = seed_start;
        const int64_t margin = static_cast<int64_t>(profile.timing_search_margin);
        const int64_t lo = std::max<int64_t>(0, seed - margin);
        const int64_t hi =
            std::min<int64_t>(static_cast<int64_t>(n), seed + margin + (int64_t)L);
        env.roi_start = lo;
        env.roi_end = hi;
    }

    // Locate the first SYNC start near the seed.
    AcquisitionState st;
    st.acq_seed = seed_start;
    const bool bounded = profile.sync_strategy == SyncStrategy::LockedBounded;
    if (bounded) {
        if (!acquire_locked_bounded(env, st, seed_start)) {
            finish_timing_profile();
            return false;
        }
    } else {
        if (!acquire_unbounded_seed(env, st, seed_start)) {
            finish_timing_profile();
            return false;
        }
    }

    // A real QM35825 has a repeatable startup transient: the first few SYNCs
    // can score below the normal verification threshold, so the detector seeds
    // the first strong SYNC.  Once locked, walk backwards on the known grid
    // with a relaxed threshold, stopping at the first missing predecessor.
    // LockedBounded disables this by default: walking backwards is exactly how
    // a lock silently drifts onto SYNC #2/#3.
    //
    // The fold path disables it unconditionally: its origin came from a
    // 64-symbol coherent fold, which is far more reliable than a single-SYNC
    // relaxed-threshold probe, and a single SYNC at -12..-14 dB cannot clear
    // 0.7*0.25 anyway.
    size_t backtracked_symbols = 0;
    const size_t backtrack_limit =
        (st.fold_used || (bounded && !profile.locked_allow_backtrack))
            ? 0
            : profile.timing_max_backtrack_symbols;
    for (size_t back = 0; back < backtrack_limit; ++back) {
        int64_t previous_start = -1;
        float previous_metric = 0.0f;
        std::complex<float> previous_corr(0.0f, 0.0f);
        const float relaxed = 0.7f * profile.radar_verification_threshold;
        if (local_best(env, st, st.first_start -
                                      static_cast<int64_t>(kQm35SamplesPerSymbol),
                       8, relaxed, &previous_start, &previous_metric,
                       &previous_corr)) {
            st.first_start = previous_start;
            st.first_metric = previous_metric;
            st.first_corr = previous_corr;
            ++backtracked_symbols;
        } else
            break;
    }
    out.backtrack_symbols = static_cast<uint32_t>(backtracked_symbols);
    const int64_t start = st.first_start;
    if (start < 0) {
        finish_timing_profile();
        return false;
    }
    const int64_t period = static_cast<int64_t>(kQm35SamplesPerSymbol);

    bool grid_applied = false;
    double grid_period = static_cast<double>(period);
    int64_t grid_start = start;
    if (bounded && profile.locked_grid_timing && !st.fold_used) {
        if (!track_locked_grid(env, st, start, period, &grid_period,
                               &grid_start)) {
            ++out.grid_fallbacks;
            out.peak_samples.clear();
            out.peak_metrics.clear();
            out.peak_corr.clear();
        } else {
            grid_applied = true;
        }
    }

    // Hot-path zero allocation: the tracked-peak repetition indices reuse the
    // worker scratch buffer (cleared here exactly as the old local vector was)
    // instead of constructing a per-packet std::vector.
    std::vector<size_t>& peak_indices = scratch.timing_peak_indices;
    peak_indices.clear();

    // LockedBounded's per-SYMBOL tracking floor.  A first peak below it cannot
    // be the head of a usable train: the tail probe already swept 32 tail
    // symbols at exactly this floor and found nothing confirmable, so running
    // the 256-SYMBOL train here is provably doomed (it would truncate to a
    // single peak and fall through to the post-train fold retry anyway).
    // Skipping it saves 256 x 27 single-SYMBOL correlations; the fold below is
    // the path that actually rescues these packets.  A first peak ABOVE the
    // floor (e.g. the 5.2 m tail-probe rescue at 0.14..0.20) still runs the
    // full train, because there the clean tail really does track.
    const float bounded_track_floor =
        std::min(profile.radar_verification_threshold, 0.08f);
    const bool skip_doomed_bounded_train =
        bounded && profile.locked_tail_probe && !grid_applied && !st.fold_used &&
        (st.first_metric < bounded_track_floor);
    if (!grid_applied && !st.fold_used && !skip_doomed_bounded_train)
        track_sync_train(env, st, start, period, backtracked_symbols);

    // ---- Post-train fold retry -------------------------------------------
    // "All single-symbol paths failed" is also observable one stage later: at
    // -4..-8 dB in-band SNR the narrow PAR search DOES acquire SYNC #0 (Task 1)
    // but the per-SYNC tracking threshold (radar_verification_threshold = 0.25)
    // still truncates the train to a single peak.  That result is unusable --
    // fit_cfo needs >= 4 peaks and stage_cir_softchips needs more than
    // cir_skip_initial_repetitions -- so the packet was already lost, just at
    // CfoFailed instead of TimingFailed.  Retry the fold, which measures a
    // usable train directly.  profile.min_valid_peaks ("min SYNC peaks to
    // accept timing") is the existing contract for "this train is too short to
    // be a timing result"; it was declared but never enforced, and this is
    // where it belongs.
    const bool allow_post_train_fold = !bounded || profile.locked_tail_probe;
    size_t min_train_peaks = std::max<size_t>(4, profile.min_valid_peaks);
    if (st.raw_par_hit && profile.rx_par_confirm) {
        // A weak single-SYNC PAR hit may start tracking with too few phase
        // points for a stable CFO estimate. Give the existing fold a chance
        // to build its coherent anchors in that case.
        min_train_peaks = std::max(
            min_train_peaks,
            std::min(profile.cfo_min_fit_repetitions, profile.preamble_repetitions));
    }
    if (allow_post_train_fold && !st.fold_used && !grid_applied &&
        out.peak_samples.size() < min_train_peaks) {
        const int64_t fold_centre =
            bounded ? static_cast<int64_t>(seed_start) : st.acq_seed;
        const int64_t fold_radius = bounded
            ? std::max<int64_t>(1, static_cast<int64_t>(profile.locked_roi_samples))
            : std::min<int64_t>(static_cast<int64_t>(profile.timing_search_margin),
                                8192);
        if (fold_centre >= 0) {
            std::vector<int64_t> saved_ps = std::move(out.peak_samples);
            std::vector<float> saved_pm = std::move(out.peak_metrics);
            std::vector<std::complex<float>> saved_pc = std::move(out.peak_corr);
            out.peak_samples.clear();
            out.peak_metrics.clear();
            out.peak_corr.clear();
            if (!try_fold_acquire(env, st, fold_centre, fold_radius,
                                  /*fast_only=*/false)) {
                out.peak_samples = std::move(saved_ps);
                out.peak_metrics = std::move(saved_pm);
                out.peak_corr = std::move(saved_pc);
            }
        }
    }

    out.detected_peaks = out.peak_samples.size();
    out.sync_peaks_folded = st.fold_used;
    if (st.fold_used)
        out.cir_grid_repetitions = profile.preamble_repetitions;
    out.expected_peaks = profile.preamble_repetitions;
    out.metric = st.first_metric;

    // Linear fit of peak positions vs repetition index -> measured period and
    // the preamble start.  peak(k) = start + k*period  =>  least squares
    // slope = period, intercept = the SYNC-start at index 0.  This matches
    // MATLAB trackCandidatePreamble: startSample = round(polyfit(...,1).coef
    // (2)), i.e. the start is derived from the WHOLE tracked train rather than
    // the first peak, so a slightly-off first peak does not bias the origin.
    const size_t np = out.peak_samples.size();
    if (st.fold_used) {
        // measured_period / preamble_start_sample were already set by the
        // fold path's 16-super-anchor fit (see try_fold_acquire).  The legacy
        // per-SYNC train never ran, so there is nothing to refit here.
    } else if (grid_applied) {
        // LockedBounded grid: period and origin come from the anchor
        // least-squares fit, NOT from the single-point train (whose positions
        // are `period` by construction and would make integer_chip_grid()
        // blind to a real clock offset).
        out.measured_period = grid_period;
        out.preamble_start_sample = grid_start;
    } else if (np >= 2) {
        detail::LinearFit fit;
        for (size_t k = 0; k < np; ++k) {
            const double x = static_cast<double>(peak_indices[k]);
            const double y = static_cast<double>(out.peak_samples[k]);
            fit.add(x, y);
        }
        const double denom = fit.determinant();
        if (std::abs(denom) > 1e-12) {
            out.measured_period = fit.slope(denom);
            // peak_samples use the symbol-END convention (start + L - 1), so
            // the intercept is also symbol-END; convert back to the SYNC start.
            const double intercept_end = fit.intercept(out.measured_period);
            out.preamble_start_sample = std::max<int64_t>(
                0, static_cast<int64_t>(std::llround(intercept_end)) -
                       (static_cast<int64_t>(L) - 1));
        } else {
            out.measured_period = static_cast<double>(period);
            out.preamble_start_sample = start;
        }
    } else {
        out.measured_period = static_cast<double>(period);
        out.preamble_start_sample = start;
    }

    out.ok = true;
    finish_timing_profile();
    return true;
}

namespace detail {

size_t
normalized_mf(const std::complex<float>* rx,
              size_t n,
              size_t roi_start,
              size_t roi_end,
              const std::complex<float>* tmpl,
              size_t L,
              float tmpl_energy,
              std::vector<float>& metric, // caller-owned, >= n
              float* best_metric)
{
    const size_t jmax = (roi_end >= L) ? roi_end - L + 1 : roi_start;
    size_t best = roi_start;
    float bm = -1.0f;
    if (jmax <= roi_start)
        return best;
    // window power of rx over the same L-length window (running sum)
    float wsum = sample_power_sum(rx + roi_start, L);
    for (size_t j = roi_start; j < jmax; ++j) {
        const auto corr = conjugate_dot_scalar(rx + j, tmpl, L);
        const float m = normalized_correlation_power(corr, wsum, tmpl_energy);
        metric[j] = m;
        if (m > bm) {
            bm = m;
            best = j;
        }
        if (j + 1 < jmax) {
            wsum -= std::norm(rx[j]);
            wsum += std::norm(rx[j + L]);
        }
    }
    if (best_metric)
        *best_metric = bm;
    return best;
}

} // namespace detail

} // namespace openuwb::demod::core
