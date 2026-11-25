// SPDX-License-Identifier: GPL-3.0-or-later
#include <openuwb/phy/hrp/rx/demod_core.h>
#include <openuwb/phy/hrp/rx/softchip_fir_diag.h>
#include <openuwb/phy/hrp/rx/demod_execution.h>
#include <openuwb/phy/hrp/rx/demod_finalize.h>
namespace openuwb::demod::core {
namespace {

struct StageClock {
    std::chrono::steady_clock::time_point t_prev = std::chrono::steady_clock::now();

    void lap(uint64_t& dst, uint64_t& total_us) {
        const auto now = std::chrono::steady_clock::now();
        dst = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(now - t_prev)
                .count());
        total_us += dst;
        t_prev = now;
    }

    void add_elapsed(uint64_t& dst, uint64_t& total_us) {
        const auto now = std::chrono::steady_clock::now();
        const uint64_t dt = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(now - t_prev)
                .count());
        dst += dt;
        total_us += dt;
        t_prev = now;
    }
};

struct RotationStepResult {
    uint64_t rotate_us = 0;
    uint32_t kernel = 0;
};

struct SyncBootstrapState {
    size_t sfd_len = 0;
    CfoResult stage2_cfo;
    DerotatedSpan predicted_span;
    DerotatedSpan derotated_span;
};

struct SfdAnchorConfig {
    size_t full_repetitions = 0;
    int64_t bootstrap_sfd_start = 0;
    double period_hint = 0.0;
    size_t template_len = 0;
};

// The derotation kernels (the vector block form and the legacy recurrence)
// restart their phase from the FIRST sample of the requested range, so two
// overlapping rotations started at different origins differ by up to ~2e-5.
// Rounding every rotation outward to the kernel's 1024-sample absolute block
// grid makes each output sample a pure function of its index, which is what
// lets the SFD stage extend its coverage incrementally without perturbing the
// bits it has already rotated.
constexpr size_t kRotationAlign = 1024;

SampleRange align_rotation_range(const SampleRange& range, size_t n)
{
    if (range.empty())
        return range;
    SampleRange aligned = range;
    aligned.lo = (aligned.lo / kRotationAlign) * kRotationAlign;
    const size_t hi_up =
        ((aligned.hi + kRotationAlign - 1) / kRotationAlign) * kRotationAlign;
    aligned.hi = std::min(n, hi_up);
    return aligned;
}

RotationStepResult rotate_window_range(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    double cfo_hz,
    const SampleRange& range,
    uint32_t kernel_flag,
    DemodScratch& scratch,
    DemodResult& res)
{
    RotationStepResult step;
    const SampleRange aligned = align_rotation_range(range, n);
    derotate_cfo_range(rx, n, profile, cfo_hz, aligned.lo, aligned.hi,
                       scratch, step.rotate_us, step.kernel);
    res.stage_cfo_rotate_us += step.rotate_us;
    res.stage_cfo_rotation_kernels |= (step.kernel ? kernel_flag : 0u);
    res.stage_cfo_rotated_samples += aligned.size();
    ++res.stage_cfo_rotation_passes;
    return step;
}

uint64_t ensure_derotated_range(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    double cfo_hz,
    int64_t req_lo,
    int64_t req_hi,
    SyncBootstrapState& sync_state,
    DemodScratch& scratch,
    DemodResult& res)
{
    // `span.range` is the contract "every sample in this interval is already
    // derotated with `cfo_hz`/origin".  It must therefore stay a single,
    // CONTIGUOUS interval: growing it across an unrotated gap would let a
    // later caller read stale samples.  The incremental merge below is only
    // valid when `req` overlaps or adjoins the span; a disjoint request is
    // handled as a fresh span instead of a union.
    const SampleRange req =
        align_rotation_range(clamp_sample_range(req_lo, req_hi, n), n);
    if (req.empty())
        return 0;
    DerotatedSpan& span = sync_state.derotated_span;
    const bool compatible =
        span.valid && span.cfo_hz == cfo_hz && span.phase_origin_sample == 0;
    const bool disjoint =
        !compatible || req.hi < span.range.lo || req.lo > span.range.hi;
    uint64_t total_us = 0;
    if (disjoint) {
        // Rotate the request on its own and replace the span.  Re-anchoring
        // (rather than bridging the gap) keeps the invariant exact and avoids
        // rotating samples no stage is going to read.
        const RotationStepResult step = rotate_window_range(
            rx, n, profile, cfo_hz, req, 2u, scratch, res);
        total_us += step.rotate_us;
        span = DerotatedSpan{req, cfo_hz, 0, true};
        return total_us;
    }
    if (span.range.covers(req))
        return 0;
    if (req.lo < span.range.lo) {
        const SampleRange gap{req.lo, std::min(req.hi, span.range.lo)};
        if (!gap.empty())
            total_us += rotate_window_range(
                rx, n, profile, cfo_hz, gap, 2u, scratch, res).rotate_us;
    }
    if (req.hi > span.range.hi) {
        const SampleRange gap{std::max(req.lo, span.range.hi), req.hi};
        if (!gap.empty())
            total_us += rotate_window_range(
                rx, n, profile, cfo_hz, gap, 2u, scratch, res).rotate_us;
    }
    span.range.lo = std::min(span.range.lo, req.lo);
    span.range.hi = std::max(span.range.hi, req.hi);
    return total_us;
}

bool search_sfd_with_wide_rescue(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<int8_t>& sfd_seq,
    const std::vector<std::complex<float>>& template_wf,
    size_t sfd_len,
    SyncBootstrapState& sync_state,
    DemodScratch& scratch,
    DemodResult& res)
{
    (void)sfd_len;
    return stage_sfd(
        scratch.derotated.data(), n, profile, res.timing, sfd_seq, template_wf,
        res.sfd, scratch,
        [&](size_t lo, size_t hi) {
            (void)ensure_derotated_range(
                rx, n, profile, res.cfo.cfo_hz,
                static_cast<int64_t>(lo), static_cast<int64_t>(hi),
                sync_state, scratch, res);
        });
}

// ---------------------------------------------------------------------------
// ML CFO helper (see docs/OpenCode执行方案_CFO估计架构重构_2026-11-24.md §2).
//
// Runs the complex-domain periodogram estimator on an SFD/timing-anchored
// preamble window and publishes it through the public CfoResult so the rest of
// the pipeline (rotation, refit contract, telemetry) is unchanged.
// `peak_ref` only supplies the diagnostic first/last sample coordinates.
// ---------------------------------------------------------------------------
bool estimate_stage_cfo_ml(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    int64_t origin_sample,
    size_t symbol_begin,
    size_t symbol_end,
    bool has_cfo_hint,
    double cfo_hint_hz,
    const TimingResult& peak_ref,
    DemodScratch& scratch,
    CfoResult& out)
{
    if (rx == nullptr || origin_sample < 0 || symbol_end <= symbol_begin)
        return false;
    if (symbol_end - symbol_begin < 8)
        return false;

    CfoMlConfig cfg;
    cfg.symbol_begin = symbol_begin;
    cfg.symbol_end = symbol_end;
    cfg.origin_sample = origin_sample;
    cfg.nominal_period = static_cast<double>(kQm35SamplesPerSymbol);
    cfg.period_hint = peak_ref.measured_period;
    cfg.cfo_hint_hz = cfo_hint_hz;
    cfg.has_cfo_hint = has_cfo_hint;
    cfg.sfo_from_cfo = profile.rx_sfo_from_cfo;
    cfg.carrier_freq_hz = profile.carrier_freq_hz;
    cfg.sample_rate_hz = profile.sample_rate;
    cfg.max_abs_cfo_hz = 0.0;

    const auto t0 = std::chrono::steady_clock::now();
    CfoMlResult ml;
    if (!estimate_cfo_ml(rx, n, profile, template_wf, cfg, scratch, ml) || !ml.ok)
        return false;
    const uint64_t fit_us = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - t0).count());

    out = CfoResult{};
    out.ok = true;
    out.cfo_hz = ml.cfo_hz;
    out.peaks_used = ml.symbols_used;
    out.skipped_peaks = symbol_begin;
    if (peak_ref.peak_samples.size() > symbol_begin)
        out.fit_first_peak_sample = peak_ref.peak_samples[symbol_begin];
    if (!peak_ref.peak_samples.empty())
        out.fit_last_peak_sample = peak_ref.peak_samples.back();
    out.fit_us = fit_us;
    return true;
}

bool run_initial_cfo_and_bootstrap_sfd(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    StageClock& clock,
    DemodScratch& scratch,
    DemodResult& res,
    SyncBootstrapState& sync_state)
{
    // Stage 2: ML CFO fit + bootstrap SFD search window rotation (Task A1).
    // The safe interior window keeps the startup transient and any leading
    // acquisition error away from the estimator; a short/abnormal frame falls
    // back to the legacy fit so every edge case stays covered.
    //
    // In the default legacy refit mode Stage 3 always runs the SFD-anchored
    // 232-SYMBOL full-preamble estimate, so Stage 2 only has to bootstrap the
    // SFO and the first rotation: a single 64-SYMBOL centred ML pass (48 on a
    // 64-SYNC fixture) buys <50 Hz accuracy for ~1/3 of the work.  Modes 1/2
    // can skip Stage 3, so they keep the full interior window to satisfy
    // cfo_refit_min_peaks and the auto/legacy rotation contract.
    const size_t n_rep = profile.preamble_repetitions;
    const size_t margin = (n_rep >= 16) ? std::min<size_t>(32, n_rep / 8) : 0;
    size_t ml_begin = margin;
    size_t ml_end = n_rep - margin;
    if (profile.cfo_refit_mode == 0 && n_rep >= 64 && ml_end > ml_begin) {
        const size_t interior = ml_end - ml_begin;
        const size_t window = std::min<size_t>(64, interior);
        ml_begin += (interior - window) / 2;
        ml_end = ml_begin + window;
    }
    bool cfo_ok = false;
    if (res.timing.preamble_start_sample >= 0 && ml_end > ml_begin) {
        cfo_ok = estimate_stage_cfo_ml(
            rx, n, profile, template_wf, res.timing.preamble_start_sample,
            ml_begin, ml_end, /*has_cfo_hint=*/false, 0.0, res.timing,
            scratch, res.cfo);
    }
    if (!cfo_ok && !fit_cfo(rx, n, profile, res.timing, res.cfo, scratch)) {
        clock.lap(res.stage_cfo_us, res.stage_total_us);
        res.stage_cfo_fit_us += res.cfo.fit_us;
        res.status = DemodStatus::CfoFailed;
        return false;
    }
    if (profile.rx_sfo_from_cfo &&
        profile.carrier_freq_hz > 0.0 &&
        std::isfinite(profile.carrier_freq_hz) &&
        std::isfinite(res.cfo.cfo_hz)) {
        const double nominal_period = static_cast<double>(kQm35SamplesPerSymbol);
        const double denom = 1.0 + res.cfo.cfo_hz / profile.carrier_freq_hz;
        if (denom > 0.0 && std::isfinite(denom)) {
            res.timing.measured_period = nominal_period / denom;
        }
    }
    scratch.cached_seq.update(profile.code_index, profile.sfd_mode);
    const std::vector<int8_t>& sfd_seq = scratch.cached_seq.sfd_seq;
    sync_state.sfd_len = sfd_seq.size() * kQm35SamplesPerSymbol;
    const SampleRange sfd1_range = compute_sfd_window_range(
        n, res.timing, profile.preamble_repetitions, sync_state.sfd_len);
    // Stage-2 CFO, kept so the single-pass rotation mode can restore it after
    // the SFD-anchored refit has been evaluated.
    sync_state.stage2_cfo = res.cfo;
    res.cfo_initial_hz = sync_state.stage2_cfo.cfo_hz;
    res.cfo_initial_peaks =
        static_cast<uint32_t>(sync_state.stage2_cfo.peaks_used);
    // Predicted single-pass rotation range.  `scratch.derotated` is the only
    // CFO-compensated buffer that the downstream stages read, so one rotation
    // with the Stage-2 estimate can replace the SFD-anchored second pass only
    // if it already covers everything those stages will read.  The SFD anchor
    // is not known yet, so the range is predicted from the Stage-1 timing with
    // an explicit symbol/backtrack slack and re-validated against the actual
    // [sec_range.lo, sec_range.hi) after the bootstrap SFD (a miss falls back
    // to the legacy second rotation).
    const SampleRange predicted_range = compute_predicted_single_pass_range(
        n, profile, res.timing, sfd1_range);
    sync_state.predicted_span = DerotatedSpan{
        predicted_range, sync_state.stage2_cfo.cfo_hz, 0,
        !predicted_range.empty()};
    const SampleRange rot1_range = sync_state.predicted_span.valid
        ? sync_state.predicted_span.range
        : sfd1_range;
    const RotationStepResult rot1 = rotate_window_range(
        rx, n, profile, res.cfo.cfo_hz, rot1_range, 1u, scratch, res);
    sync_state.derotated_span = DerotatedSpan{
        rot1_range, res.cfo.cfo_hz, 0, !rot1_range.empty()};
    res.cfo.rotate_us = rot1.rotate_us;
    res.cfo.rotation_kernel = rot1.kernel;
    clock.lap(res.stage_cfo_us, res.stage_total_us);
    res.stage_cfo_fit_us += res.cfo.fit_us;

    // Stage 3: bootstrap SFD timing.  Timing acquisition can start several
    // SYNCs late when a weak leading repetition falls below its tracking
    // threshold.  Its 64-point phase train would then spill into the SFD and
    // corrupt a CFO fit, even though the SFD search can recover the true frame
    // origin by testing earlier hypotheses.  Use this first SFD only to anchor
    // the preamble; the final reported CFO is re-fitted below from that anchor.
    scratch.initial_timing = res.timing;
    const auto sfd_bootstrap_t0 = std::chrono::steady_clock::now();
    const bool bootstrap_sfd_ok = search_sfd_with_wide_rescue(
        rx, n, profile, sfd_seq, template_wf, sync_state.sfd_len, sync_state,
        scratch, res);
    res.stage_sfd_bootstrap_us = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - sfd_bootstrap_t0).count());
    res.sfd_bootstrap_windows = res.sfd.search_windows;
    res.sfd_bootstrap_coarse_correlations = res.sfd.coarse_correlations;
    res.sfd_bootstrap_fine_correlations = res.sfd.fine_correlations;
    res.sfd_initial_predicted_sample = res.sfd.expected_start_sample;
    res.sfd_bootstrap_detected_sample = res.sfd.sfd_start_sample;
    res.sfd_bootstrap_first_threshold_backtrack_symbols =
        res.sfd.first_threshold_backtrack_symbols;
    if (!bootstrap_sfd_ok) {
        clock.lap(res.stage_sfd_us, res.stage_total_us);
        res.status = DemodStatus::SfdFailed;
        return false;
    }
    return true;
}

void build_sfd_anchored_timing(
    const SfdAnchorConfig& cfg,
    const TimingResult& prev_timing,
    const TimingResult& initial_timing,
    TimingResult& cfo_timing)
{
    // The SFD anchor gives the exact end of the last SYNC, so the complete
    // preamble grid is reconstructed from it.  No tail re-tracking runs any
    // more: the ML CFO estimator samples this grid directly (full preamble).
    const double measured_period = cfg.period_hint;
    const int64_t full_start = cfg.bootstrap_sfd_start -
        static_cast<int64_t>(std::llround(
            static_cast<double>(cfg.full_repetitions) * measured_period));
    // The SFD-anchored grid replaces the seeded search result, but the
    // synchronization-strategy diagnostics describe the FIRST seeded search
    // (the only one the strategy controls), so carry them over.
    cfo_timing.bounded_attempts = prev_timing.bounded_attempts;
    cfo_timing.bounded_success = prev_timing.bounded_success;
    cfo_timing.bounded_failures = prev_timing.bounded_failures;
    cfo_timing.bounded_legacy_scans = prev_timing.bounded_legacy_scans;
    cfo_timing.bounded_inner_used = prev_timing.bounded_inner_used;
    cfo_timing.bounded_roi_samples = prev_timing.bounded_roi_samples;
    cfo_timing.backtrack_symbols = prev_timing.backtrack_symbols;
    // Diagnostic-only carry-over (no interval math): the locked-tracking
    // diagnostics describe the FIRST seeded search, which is the only one
    // the strategy controls, and res.timing is replaced below.
    cfo_timing.grid_anchor_attempts = prev_timing.grid_anchor_attempts;
    cfo_timing.grid_anchor_hits = prev_timing.grid_anchor_hits;
    cfo_timing.grid_single_points = prev_timing.grid_single_points;
    cfo_timing.grid_single_valid = prev_timing.grid_single_valid;
    cfo_timing.grid_fallbacks = prev_timing.grid_fallbacks;
    cfo_timing.tail_probe_attempts = prev_timing.tail_probe_attempts;
    cfo_timing.tail_probe_hits = prev_timing.tail_probe_hits;
    cfo_timing.tail_probe_k = prev_timing.tail_probe_k;
    cfo_timing.tail_probe_best_metric = prev_timing.tail_probe_best_metric;
    cfo_timing.tail_probe_raw_max = prev_timing.tail_probe_raw_max;
    cfo_timing.inner_probe_raw_max = prev_timing.inner_probe_raw_max;
    cfo_timing.initial_coarse_search_us = prev_timing.coarse_search_us;
    cfo_timing.initial_fine_track_us = prev_timing.fine_track_us;
    cfo_timing.ok = true;
    cfo_timing.preamble_start_sample = full_start;
    cfo_timing.preamble_start_uncropped = full_start;
    cfo_timing.measured_period = measured_period;
    cfo_timing.metric = prev_timing.metric;
    cfo_timing.detected_peaks = initial_timing.detected_peaks;
    cfo_timing.sync_peaks_folded = initial_timing.sync_peaks_folded;
    cfo_timing.cir_grid_repetitions = cfg.full_repetitions;
    cfo_timing.expected_peaks = cfg.full_repetitions;
    cfo_timing.coarse_search_us = prev_timing.coarse_search_us;
    cfo_timing.fine_track_us = prev_timing.fine_track_us;
    cfo_timing.seed_search_us = prev_timing.seed_search_us;
    cfo_timing.seed_refine_narrow_attempts =
        prev_timing.seed_refine_narrow_attempts;
    cfo_timing.seed_refine_narrow_success =
        prev_timing.seed_refine_narrow_success;
    cfo_timing.seed_refine_fallbacks = prev_timing.seed_refine_fallbacks;
    cfo_timing.peak_samples.resize(cfg.full_repetitions);
    cfo_timing.peak_metrics.assign(cfg.full_repetitions, 0.0f);
    cfo_timing.peak_corr.assign(
        cfg.full_repetitions, std::complex<float>(0.0f, 0.0f));
    for (size_t k = 0; k < cfg.full_repetitions; ++k) {
        cfo_timing.peak_samples[k] = full_start +
            static_cast<int64_t>(std::llround(
                static_cast<double>(k) * measured_period)) +
            static_cast<int64_t>(cfg.template_len - 1);
    }
}

bool refit_tail_cfo_and_rotate_anchored_prefix(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    bool integer_grid_ok,
    const SfdAnchorConfig& cfg,
    DemodScratch& scratch,
    DemodResult& res,
    SyncBootstrapState& sync_state,
    TimingResult& cfo_timing)
{
    const auto sfd_cfo_t0 = std::chrono::steady_clock::now();
    // SFD-anchored full-preamble ML refit.  The SFD pins the end of the last
    // SYNC, so the whole preamble grid (minus the startup exclusion) is a
    // measured interval: no 40-SYMBOL tail retrack, no zero placeholders, and
    // a coherent baseline ~sqrt(232/40) longer than the old tail fit.
    // TODO(cfo-perf): If Stage-3 latency (~50 us, ~40 us in 232-symbol SIMD
    // sampling) becomes a bottleneck, consider stride-2 subsampling across the
    // same 232-symbol span (116 symbols + 512-pt FFT, saving ~20-24 us while
    // keeping +/-250 kHz unambiguous range; verify low-SNR CRLB / Std impact).
    const size_t n_rep = profile.preamble_repetitions;
    const size_t k_min = std::max<size_t>(8, profile.cfo_min_fit_repetitions);
    const size_t avail_skip = (n_rep > k_min) ? (n_rep - k_min) : 0;
    const size_t k_skip = std::min(profile.cfo_skip_initial_repetitions,
                                   avail_skip);
    const double p_ref = (cfg.period_hint > 0.0)
        ? cfg.period_hint
        : static_cast<double>(kQm35SamplesPerSymbol);
    const int64_t origin = cfg.bootstrap_sfd_start -
        static_cast<int64_t>(std::llround(
            static_cast<double>(n_rep) * p_ref));
    if (k_skip < n_rep && n_rep - k_skip >= 8) {
        CfoResult ml_cfo;
        if (estimate_stage_cfo_ml(
                rx, n, profile, template_wf, origin, k_skip, n_rep,
                /*has_cfo_hint=*/sync_state.stage2_cfo.ok,
                sync_state.stage2_cfo.cfo_hz, cfo_timing, scratch, ml_cfo)) {
            res.cfo = ml_cfo;
            res.cfo_refit_hz = ml_cfo.cfo_hz;
            res.cfo_refit_peaks = static_cast<uint32_t>(ml_cfo.peaks_used);
            res.cfo_refit_fit_us = ml_cfo.fit_us;
            res.stage_cfo_fit_us += ml_cfo.fit_us;
            // Publish the measured correlations on the anchored grid for
            // downstream diagnostics (the grid is no longer a placeholder
            // beyond the skipped startup symbols).
            const size_t used = std::min<size_t>(
                ml_cfo.peaks_used, scratch.cfo_sym_z.size());
            for (size_t i = 0; i < used; ++i) {
                const size_t dst = k_skip + i;
                if (dst < cfo_timing.peak_corr.size())
                    cfo_timing.peak_corr[dst] = scratch.cfo_sym_z[i];
            }
        }
    }
    res.stage_sfd_cfo_refit_us += static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - sfd_cfo_t0).count());
    if (!res.cfo.ok)
        return false;
    if (res.cfo_refit_hz == 0.0 && res.cfo_refit_peaks == 0) {
        // ML refit unavailable (extreme window): keep the Stage-2 estimate and
        // let the refit-mode state machine below decide about the rotation.
        res.cfo_refit_hz = res.cfo.cfo_hz;
        res.cfo_refit_peaks = static_cast<uint32_t>(res.cfo.peaks_used);
        res.cfo_refit_fit_us = res.cfo.fit_us;
    }

    // Second rotation (Task A1): covers final SFD, CIR (54 SYNCs starting
    // from 10th), phase estimation, and initial soft-chips / PHR.
    const SampleRange final_sfd_range = compute_sfd_window_range(
        n, cfo_timing, cfg.full_repetitions, sync_state.sfd_len);
    const SampleRange sec_range = compute_anchored_prefix_range(
        n, profile, integer_grid_ok, cfo_timing, final_sfd_range);

    // ---- Single-pass rotation decision ---------------------------------
    // The second rotation exists only because the SFD-anchored refit moves
    // the estimate.  When the Stage-2 estimate is good enough, the pass-1
    // rotation (already widened to the predicted union range) is the only
    // derotation the packet needs.  Every condition below is a hard
    // requirement; any miss keeps the legacy refit + second rotation.
    const uint64_t refit_span = static_cast<uint64_t>(sec_range.size());
    const double refit_delta_hz = res.cfo.cfo_hz - sync_state.stage2_cfo.cfo_hz;
    const double refit_drift_rad = std::abs(
        2.0 * M_PI * refit_delta_hz *
        static_cast<double>(refit_span) / profile.sample_rate);
    res.cfo_refit_span_samples = refit_span;
    res.cfo_refit_phase_drift_rad = refit_drift_rad;
    bool want_skip = false;
    uint32_t skip_action = 0;
    if (profile.cfo_refit_mode == 2) {
        want_skip = true;
        skip_action = 2;
    } else if (profile.cfo_refit_mode == 1) {
        want_skip =
            sync_state.stage2_cfo.peaks_used >= profile.cfo_refit_min_peaks &&
            refit_drift_rad <= profile.cfo_refit_skip_rad;
        skip_action = 1;
    }
    if (want_skip &&
        sync_state.predicted_span.covers_with_cfo(
            sec_range, sync_state.stage2_cfo.cfo_hz, 0)) {
        // Keep the Stage-2 estimate: pass 1 already applied it to the whole
        // predicted span, which provably contains [sec_range.lo, sec_range.hi).
        res.cfo = sync_state.stage2_cfo;
        res.cfo_refit_action = skip_action;
        sync_state.derotated_span = sync_state.predicted_span;
        return true;
    }
    if (profile.cfo_refit_mode != 0)
        res.cfo_refit_action = 3; // requested skip, fell back (fail-closed)
    const RotationStepResult rot2 = rotate_window_range(
        rx, n, profile, res.cfo.cfo_hz, sec_range, 2u, scratch, res);
    res.cfo.rotate_us = rot2.rotate_us;
    res.cfo.rotation_kernel = rot2.kernel;
    res.cfo_refit_rotate_us = rot2.rotate_us;
    res.stage_sfd_cfo_refit_us += rot2.rotate_us;
    sync_state.derotated_span = DerotatedSpan{
        sec_range, res.cfo.cfo_hz, 0, !sec_range.empty()};
    return true;
}

bool anchor_and_verify_final_sfd(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    bool integer_grid_ok,
    SyncBootstrapState& sync_state,
    StageClock& clock,
    DemodScratch& scratch,
    DemodResult& res)
{
    // Use the SFD as the frame anchor.  The first SYNCs can be weak or can be
    // missed by acquisition, so reconstruct the complete SYNC grid from the
    // SFD coordinate instead of re-tracking the preamble from its first
    // repetition.  The ML CFO estimator samples that full-preamble grid.
    const size_t full_repetitions = profile.preamble_repetitions;
    if (full_repetitions < 2 || res.timing.measured_period <= 0.0) {
        clock.lap(res.stage_sfd_us, res.stage_total_us);
        res.status = DemodStatus::CfoFailed;
        return false;
    }

    SfdAnchorConfig cfg;
    cfg.full_repetitions = full_repetitions;
    cfg.bootstrap_sfd_start = res.sfd.sfd_start_sample;
    cfg.period_hint = res.timing.measured_period;
    cfg.template_len = template_wf.size();

    const std::vector<int8_t>& sfd_seq = scratch.cached_seq.sfd_seq;
    auto run_final_sfd = [&]() -> bool {
        const auto sfd_final_t0 = std::chrono::steady_clock::now();
        const bool ok = search_sfd_with_wide_rescue(
            rx, n, profile, sfd_seq, template_wf, sync_state.sfd_len,
            sync_state, scratch, res);
        res.stage_sfd_final_us += static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - sfd_final_t0).count());
        return ok;
    };

    TimingResult cfo_timing;
    build_sfd_anchored_timing(cfg, res.timing, scratch.initial_timing,
                              cfo_timing);
    if (!refit_tail_cfo_and_rotate_anchored_prefix(
            rx, n, profile, template_wf, integer_grid_ok, cfg, scratch, res,
            sync_state, cfo_timing)) {
        clock.lap(res.stage_sfd_us, res.stage_total_us);
        res.status = DemodStatus::CfoFailed;
        return false;
    }
    res.timing = std::move(cfo_timing);

    // Re-run SFD on the final MATLAB-aligned CFO compensation.  This removes
    // any bootstrap-CFO bias from the SFD/CIR path while retaining the same
    // full-rate refinement contract as MATLAB refineTimingWithNsSfd.
    //
    // When the refit reused the Stage-2 rotation (actions 1/2), or the refit
    // CFO moved the 8-symbol SFD phase by less than 0.05 rad, the bootstrap SFD
    // result is already the final one: the derotated bits did not change (or
    // changed too little to move the peak by a sample).  Skipping the complete
    // +/-64 + offset-window + grid search is a pure win on those packets.
    bool final_sfd_ok = false;
    const bool rotation_reused =
        res.cfo_refit_action == 1 || res.cfo_refit_action == 2;
    bool tiny_refit_delta = false;
    if (!rotation_reused && res.cfo.ok) {
        const double drift_rad = std::abs(
            2.0 * M_PI * (res.cfo.cfo_hz - sync_state.stage2_cfo.cfo_hz) *
            (8.0 * static_cast<double>(kQm35SamplesPerSymbol)) /
            profile.sample_rate);
        tiny_refit_delta = drift_rad < 0.05;
    }
    if (rotation_reused || tiny_refit_delta) {
        final_sfd_ok = res.sfd.ok;
    } else {
        final_sfd_ok = run_final_sfd();
    }
    // No legacy-radius retry on a failed final SFD.  The old retry re-ran the
    // identical CFO refit on the same raw frame and then the identical final
    // SFD on identically derotated data, so it could never change
    // `final_sfd_ok`.  The retry path (and its `seed_radius` plumbing) is gone.
    res.sfd_final_windows = res.sfd.search_windows;
    res.sfd_final_coarse_correlations = res.sfd.coarse_correlations;
    res.sfd_final_fine_correlations = res.sfd.fine_correlations;
    if (!final_sfd_ok) {
        clock.lap(res.stage_sfd_us, res.stage_total_us);
        res.status = DemodStatus::SfdFailed;
        return false;
    }
    // Refresh the SYNC period from the final full-preamble CFO.  With a
    // coupled crystal SFO = -CFO/fc, so this period is exact to the CFO
    // accuracy and is far better conditioned than the Stage-1 OLS period
    // (which the old 40-symbol tail fit could not improve).
    if (profile.rx_sfo_from_cfo && profile.carrier_freq_hz > 0.0 &&
        std::isfinite(profile.carrier_freq_hz) && res.cfo.ok &&
        std::isfinite(res.cfo.cfo_hz)) {
        const double denom = 1.0 + res.cfo.cfo_hz / profile.carrier_freq_hz;
        if (denom > 0.0 && std::isfinite(denom)) {
            res.timing.measured_period =
                static_cast<double>(kQm35SamplesPerSymbol) / denom;
        }
    }
    res.timing.preamble_start_sample =
        res.sfd.sfd_start_sample - static_cast<int64_t>(std::llround(
            static_cast<double>(profile.preamble_repetitions) *
            res.timing.measured_period));
    clock.lap(res.stage_sfd_us, res.stage_total_us);
    // Diagnostic remainder of the SFD parent stage.  The four child sub-timers
    // above are strictly children of `stage_sfd_us` — the parent already
    // includes them, so this value is diagnostics only and is never added to
    // `stage_sfd_us` or `stage_total_us`.  Clamp at 0 so an independent-clock
    // overshoot can never produce a negative (wrapped) value.
    {
        const uint64_t sfd_children_us =
            res.stage_sfd_bootstrap_us + res.stage_sfd_tail_retrack_us +
            res.stage_sfd_cfo_refit_us + res.stage_sfd_final_us;
        res.stage_sfd_anchor_other_us =
            (res.stage_sfd_us >= sfd_children_us)
                ? res.stage_sfd_us - sfd_children_us
                : 0;
    }
    return true;
}

// Diagnostic-only segmented FIR replay.  It is intentionally outside the
// lap/add_elapsed parents and only runs when the explicit diagnostic
// profile flag is set; formal timing runs leave it off.
void maybe_diagnose_fir(
    const Qm35825Profile& profile,
    size_t n,
    DemodScratch& scratch,
    DemodResult& res)
{
    if (!scratch.diagnose_fir_segments)
        return;
    if (!(res.cir.ok && res.ns_sfd.ok && res.phr.ok && res.payload.ok))
        return;
    const size_t tap_count = res.cir.cir_filter_taps.size();
    if (tap_count == 0 || tap_count > 64 || res.cir.soft_chip_count == 0)
        return;
    const int64_t chip_start =
        res.timing.preamble_start_sample + static_cast<int64_t>(tap_count) -
        static_cast<int64_t>(profile.cir_pre_samples) - 1;
    if (chip_start < 0)
        return;
    size_t payload_symbols = 0;
    if (res.phr.psdu_length > 0) {
        const size_t psdu_bits =
            static_cast<size_t>(res.phr.psdu_length) * 8;
        const size_t num_blocks = (psdu_bits + 329) / 330;
        payload_symbols = psdu_bits + 48 * num_blocks;
    }
    detail::SoftFirSegments segments;
    if (!detail::diagnose_softchip_fir_segments(
            scratch.derotated.data(), n, res.cir.cir_filter_taps.data(),
            tap_count, chip_start, res.cir.soft_chip_count,
            res.cir.samples_per_chip, res.ns_sfd, res.phr,
            payload_symbols, scratch.corr, segments))
        return;
    res.diag_fir_valid = 1;
    res.diag_fir_preamble_us = segments.preamble_us;
    res.diag_fir_sfd_us = segments.sfd_us;
    res.diag_fir_phr_us = segments.phr_us;
    res.diag_fir_payload_us = segments.payload_us;
    res.diag_fir_extra_us = segments.extra_us;
    res.diag_fir_preamble_chips = static_cast<uint32_t>(
        segments.preamble_chips);
    res.diag_fir_sfd_chips = static_cast<uint32_t>(segments.sfd_chips);
    res.diag_fir_phr_chips = static_cast<uint32_t>(segments.phr_chips);
    res.diag_fir_payload_chips = static_cast<uint32_t>(
        segments.payload_chips);
    res.diag_fir_extra_chips = static_cast<uint32_t>(segments.extra_chips);
}

struct FullSoftchipDecodeOutputs {
    CirResult& cir;
    NsSfdResult& ns_sfd;
    PhrResult& phr;
    PayloadResult& payload;
};

bool run_full_softchip_decode_stages(
    size_t n,
    const Qm35825Profile& cir_profile,
    const Qm35825Profile& decode_profile,
    const TimingResult& timing,
    const std::vector<int8_t>& pcode,
    const std::vector<int8_t>& nsfd_seq,
    StageClock* clock,
    PayloadSubTimers* payload_sub,
    DemodScratch& scratch,
    FullSoftchipDecodeOutputs out,
    DemodResult& res,
    DemodStatus& status_out)
{
    const auto lap_if = [&](uint64_t& dst) {
        if (clock != nullptr)
            clock->lap(dst, res.stage_total_us);
    };

    if (out.cir.ok) {
        if (!stage_softchips_from_cir(scratch.derotated.data(), n, cir_profile,
                                      timing, out.cir, scratch)) {
            if (clock != nullptr)
                clock->add_elapsed(res.stage_cir_us, res.stage_total_us);
            status_out = DemodStatus::CirFailed;
            return false;
        }
        if (clock != nullptr)
            clock->add_elapsed(res.stage_cir_us, res.stage_total_us);
    } else if (!stage_cir_softchips(scratch.derotated.data(), n, cir_profile,
                                    timing, pcode, out.cir, scratch)) {
        lap_if(res.stage_cir_us);
        status_out = DemodStatus::CirFailed;
        return false;
    } else {
        lap_if(res.stage_cir_us);
    }

    if (!stage_ns_sfd(scratch.soft_chips, decode_profile, nsfd_seq,
                      kQm35ChipsPerSymbol, out.ns_sfd, scratch)) {
        lap_if(res.stage_ns_sfd_us);
        status_out = DemodStatus::SfdFailed;
        return false;
    }
    lap_if(res.stage_ns_sfd_us);

    if (!stage_phr(scratch.soft_chips, decode_profile, out.ns_sfd,
                   kQm35ChipsPerSymbol, out.phr, scratch)) {
        lap_if(res.stage_phr_us);
        status_out = DemodStatus::PhrFailed;
        return false;
    }
    lap_if(res.stage_phr_us);

    if (!stage_payload_fcs(scratch.soft_chips, decode_profile, out.phr,
                           out.ns_sfd, kQm35ChipsPerSymbol, out.payload,
                           scratch, payload_sub)) {
        lap_if(res.stage_payload_us);
        status_out = DemodStatus::PayloadFailed;
        return false;
    }
    lap_if(res.stage_payload_us);

    status_out = (out.payload.ok && out.payload.fcs_pass)
                     ? DemodStatus::Success
                     : DemodStatus::FcsFailed;
    return true;
}

void run_primary_full_stream_decode(
    size_t n,
    const Qm35825Profile& profile,
    const Qm35825Profile& eff,
    StageClock& clock,
    PayloadSubTimers& payload_sub,
    DemodScratch& scratch,
    DemodResult& res)
{
    FullSoftchipDecodeOutputs out{res.cir, res.ns_sfd, res.phr, res.payload};
    if (run_full_softchip_decode_stages(
            n, eff, profile, res.timing,
            scratch.cached_seq.pcode, scratch.cached_seq.sfd_seq,
            &clock, &payload_sub, scratch, out, res, res.status)) {
        maybe_diagnose_fir(profile, n, scratch, res);
    }
}

// `run_full_stream_fallback` retired: fused / on-demand decode fails directly
// without re-rotating the entire stream or re-running full-stream softchips.

uint64_t rotate_uncovered_tail_range(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const SampleRange& pay_range,
    DerotatedSpan& derotated_span,
    DemodScratch& scratch,
    DemodResult& res)
{
    uint64_t rot_pay_us = 0;
    uint32_t rk_pay = 0;
    derotate_cfo_range(rx, n, profile, res.cfo.cfo_hz,
                       pay_range.lo, pay_range.hi, scratch,
                       rot_pay_us, rk_pay);
    res.stage_cfo_rotate_us += rot_pay_us;
    res.stage_cfo_rotated_samples += pay_range.size();
    ++res.stage_cfo_rotation_passes;
    if (!derotated_span.valid ||
        derotated_span.cfo_hz != res.cfo.cfo_hz) {
        derotated_span = DerotatedSpan{
            pay_range, res.cfo.cfo_hz, 0, true};
    } else {
        derotated_span.range.lo = std::min(
            derotated_span.range.lo, pay_range.lo);
        derotated_span.range.hi = std::max(
            derotated_span.range.hi, pay_range.hi);
    }
    return rot_pay_us;
}

struct PrefixPhrAttempt {
    size_t prefix_chips = 0;
    bool phr_fused_used = false;
};

bool run_prefix_cir_ns_sfd_and_phr(
    size_t n,
    const Qm35825Profile& profile,
    const Qm35825Profile& eff,
    bool integer_grid_ok,
    StageClock& clock,
    DemodScratch& scratch,
    DemodResult& res,
    PrefixPhrAttempt& out)
{
    const std::vector<int8_t>& pcode = scratch.cached_seq.pcode;
    const std::vector<int8_t>& nsfd_seq = scratch.cached_seq.sfd_seq;
    const bool phr_fusion_requested = profile.rx_fused_phr_cmf;
    const bool phr_fusion_cfg_supported =
        profile.rx_fused_payload_cmf &&
        profile.cir_rake_top_k == 0 &&
        profile.cir_soft_chip_mode != CirSoftChipMode::Rake &&
        profile.cir_soft_chip_mode != CirSoftChipMode::Bypass;
    const size_t prefix_budget =
        phr_fusion_requested ? rxopt::fused_phr_prefix_budget(profile)
                             : rxopt::initial_chip_budget(profile);
    bool ok = stage_softchips_from_cir(scratch.derotated.data(), n, eff,
                                       res.timing, res.cir, scratch,
                                       prefix_budget, 0);
    clock.add_elapsed(res.stage_cir_us, res.stage_total_us);
    out.prefix_chips = res.cir.soft_chip_count;
    if (ok) {
        ok = stage_ns_sfd(scratch.soft_chips, profile, nsfd_seq,
                          kQm35ChipsPerSymbol, res.ns_sfd, scratch);
        clock.lap(res.stage_ns_sfd_us, res.stage_total_us);
    }
    out.phr_fused_used = false;
    if (ok) {
        if (phr_fusion_requested && integer_grid_ok) {
            const bool supported =
                phr_fusion_cfg_supported &&
                !res.cir.cir_filter_taps.empty() &&
                (profile.rx_fused_sfo_track ||
                 integer_chip_grid_at(res.cir.samples_per_chip,
                                      rxfusion::kSamplesPerChip,
                                      profile.fused_gate_ppm * 1e-6)) &&
                scratch.softchip.prepared;
            if (supported) {
                ++res.stage_phr_fused_attempts;
                uint64_t phr_template_us = 0;
                uint64_t phr_correlate_us = 0;
                uint64_t phr_decode_us = 0;
                uint64_t phr_decode_ns = 0;
                size_t phr_symbols = 0;
                const bool fused_ok = stage_phr_fused_cmf(
                    scratch.derotated.data(), n, profile, res.timing,
                    res.cir, scratch.soft_chips, res.ns_sfd, res.phr,
                    scratch, phr_template_us, phr_correlate_us,
                    phr_decode_us, phr_decode_ns, phr_symbols);
                res.stage_phr_fused_template_us += phr_template_us;
                res.stage_phr_fused_correlate_us += phr_correlate_us;
                res.stage_phr_fused_decode_us += phr_decode_us;
                res.stage_phr_fused_decode_ns += phr_decode_ns;
                res.stage_phr_fused_symbols =
                    static_cast<uint32_t>(phr_symbols);
                clock.lap(res.stage_phr_us, res.stage_total_us);
                // A SECDED correction needs the full-stream tail-energy
                // guard.  The reduced prefix cannot reproduce it, so keep
                // the existing conservative full-stream retry.
                if (fused_ok && res.phr.ok && !res.phr.secded_corrected) {
                    out.phr_fused_used = true;
                } else {
                    ++res.stage_phr_fused_fallbacks;
                    ok = false;
                }
            } else {
                // The explicit option only supports the full integer-grid
                // CMF.  Do not mix a fused PHR request with the legacy PHR
                // kernel; let the single full-stream retry below decide.
                ok = false;
            }
        } else {
            ok = stage_phr(scratch.soft_chips, profile, res.ns_sfd,
                           kQm35ChipsPerSymbol, res.phr, scratch);
            clock.lap(res.stage_phr_us, res.stage_total_us);
        }
    }
    if (ok && !out.phr_fused_used && res.phr.secded_corrected) {
        // The legacy PHR path can inspect the full-range tail energy when
        // SECDED corrected a longer length.  The prefix cannot, so retry
        // the complete stream rather than silently changing that decision.
        ok = false;
    }
    return ok;
}

enum class FusedPayloadOutcome {
    NotRun,
    Succeeded,
    FailedToFallback,
};

FusedPayloadOutcome try_fused_payload_decode(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    bool integer_grid_ok,
    bool phr_fused_used,
    DerotatedSpan& derotated_span,
    StageClock& clock,
    PayloadSubTimers& payload_sub,
    DemodScratch& scratch,
    DemodResult& res)
{
    if (!profile.rx_fused_payload_cmf || !integer_grid_ok)
        return FusedPayloadOutcome::NotRun;

    const bool fused_supported =
        !res.cir.cir_filter_taps.empty() &&
        (profile.rx_fused_sfo_track ||
         integer_chip_grid_at(res.cir.samples_per_chip,
                              rxfusion::kSamplesPerChip,
                              profile.fused_gate_ppm * 1e-6)) &&
        profile.cir_rake_top_k == 0 &&
        profile.cir_soft_chip_mode != CirSoftChipMode::Rake &&
        profile.cir_soft_chip_mode != CirSoftChipMode::Bypass &&
        scratch.softchip.prepared;
    if (!fused_supported)
        return FusedPayloadOutcome::NotRun;

    // Payload-region derotation is fused exactly into the CMF template and
    // per-symbol phase below, so this path deliberately does not rotate the
    // payload interval into `scratch.derotated`.  `pay_derot` therefore stays
    // at 0 and the cost is visible where it is actually spent (`pay_soft`).
    // Fail-safe: the full-stream fallback and the non-fused on-demand path
    // still rotate their own ranges before reading `scratch.derotated`.
    (void)derotated_span;

    ++res.stage_fused_attempts;
    uint64_t fused_template_us = 0;
    uint64_t fused_correlate_us = 0;
    uint64_t fused_decode_us = 0;
    uint64_t fused_decode_ns = 0;
    size_t fused_symbols = 0;
    const std::vector<int8_t>* cached_g0 =
        (phr_fused_used && scratch.fused_phr.decisions_valid)
            ? &scratch.fused_phr.g0
            : nullptr;
    const std::vector<int8_t>* cached_g1 =
        (phr_fused_used && scratch.fused_phr.decisions_valid)
            ? &scratch.fused_phr.g1
            : nullptr;
    // Same cache, soft form: the fused PHR LLRs feed the payload
    // soft retry without re-demodulating the PHR field.
    const std::vector<float>* cached_llr0 =
        (phr_fused_used && scratch.fused_phr.decisions_valid)
            ? &scratch.fused_phr.llr0
            : nullptr;
    const std::vector<float>* cached_llr1 =
        (phr_fused_used && scratch.fused_phr.decisions_valid)
            ? &scratch.fused_phr.llr1
            : nullptr;
    const double cfo_w_rad =
        2.0 * M_PI * res.cfo.cfo_hz / profile.sample_rate;
    const bool fused_ok = stage_payload_fused_cmf(
        rx, n, profile, res.timing, res.cir,
        scratch.soft_chips, res.ns_sfd, res.phr, res.payload,
        scratch, fused_template_us, fused_correlate_us,
        fused_decode_us, fused_decode_ns, fused_symbols,
        cached_g0, cached_g1, &payload_sub, cached_llr0,
        cached_llr1, cfo_w_rad);
    res.stage_fused_template_us += fused_template_us;
    res.stage_fused_correlate_us += fused_correlate_us;
    res.stage_fused_decode_us += fused_decode_us;
    res.stage_fused_decode_ns += fused_decode_ns;
    res.stage_fused_symbols = static_cast<uint32_t>(fused_symbols);
    clock.lap(res.stage_payload_us, res.stage_total_us);
    if (fused_ok && res.payload.fcs_pass) {
        res.status = DemodStatus::Success;
        maybe_diagnose_fir(profile, n, scratch, res);
        return FusedPayloadOutcome::Succeeded;
    }
    ++res.stage_fused_fallbacks;
    return FusedPayloadOutcome::FailedToFallback;
}

bool try_ondemand_extension_and_payload_decode(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const Qm35825Profile& eff,
    size_t prefix_chips,
    DerotatedSpan& derotated_span,
    StageClock& clock,
    PayloadSubTimers& payload_sub,
    DemodScratch& scratch,
    DemodResult& res)
{
    const size_t final_budget = rxopt::final_chip_budget(
        profile, res.ns_sfd.sfd_end_chip, res.phr.psdu_length);
    if (final_budget > prefix_chips) {
        // Compensate extension interval with refit CFO (Task A1)
        const size_t covered_hi =
            derotated_span.covers_with_cfo(
                derotated_span.range, res.cfo.cfo_hz, 0)
                ? derotated_span.range.hi
                : 0;
        const SampleRange pay_range =
            compute_ondemand_extension_rotation_range(
                n, profile, res.timing, res.cir.cir_filter_taps.size(),
                final_budget, covered_hi);
        if (!pay_range.empty()) {
            (void)rotate_uncovered_tail_range(
                rx, n, profile, pay_range, derotated_span, scratch, res);
        }

        CirResult ext = res.cir;
        const bool ext_ok = stage_softchips_from_cir(
            scratch.derotated.data(), n, eff, res.timing,
            ext, scratch, final_budget, prefix_chips);
        clock.add_elapsed(res.stage_cir_us, res.stage_total_us);
        if (!ext_ok)
            return false;

        res.cir.soft_extend_us =
            ext.soft_fir_us + ext.postprocess_us;
        res.cir.soft_fir_us += ext.soft_fir_us;
        res.cir.postprocess_us += ext.postprocess_us;
        res.cir.soft_chip_count = ext.soft_chip_count;
        res.cir.extended_chip_count = ext.extended_chip_count;
        res.cir.on_demand_extended = ext.on_demand_extended;
        res.cir.samples_per_chip = ext.samples_per_chip;
    }

    const bool ok = stage_payload_fcs(
        scratch.soft_chips, profile, res.phr, res.ns_sfd,
        kQm35ChipsPerSymbol, res.payload, scratch, &payload_sub);
    clock.lap(res.stage_payload_us, res.stage_total_us);
    if (!ok || !res.payload.fcs_pass)
        return false;

    res.status = DemodStatus::Success;
    maybe_diagnose_fir(profile, n, scratch, res);
    return true;
}

void run_on_demand_with_full_stream_fallback(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const Qm35825Profile& eff,
    bool integer_grid_ok,
    DerotatedSpan& derotated_span,
    StageClock& clock,
    PayloadSubTimers& payload_sub,
    DemodScratch& scratch,
    DemodResult& res)
{
    PrefixPhrAttempt prefix_phr;
    bool ok = run_prefix_cir_ns_sfd_and_phr(
        n, profile, eff, integer_grid_ok, clock, scratch, res, prefix_phr);
    if (ok) {
        const FusedPayloadOutcome fused_outcome = try_fused_payload_decode(
            rx, n, profile, integer_grid_ok, prefix_phr.phr_fused_used,
            derotated_span, clock, payload_sub, scratch, res);
        if (fused_outcome == FusedPayloadOutcome::Succeeded)
            return;
        if (fused_outcome == FusedPayloadOutcome::FailedToFallback) {
            ok = false;
        } else {
            ok = try_ondemand_extension_and_payload_decode(
                rx, n, profile, eff, prefix_phr.prefix_chips,
                derotated_span, clock, payload_sub, scratch, res);
        }
    }
    if (!ok) {
        if (!res.cir.ok)
            res.status = DemodStatus::CirFailed;
        else if (!res.ns_sfd.ok)
            res.status = DemodStatus::SfdFailed;
        else if (!res.phr.ok)
            res.status = DemodStatus::PhrFailed;
        else if (!res.payload.ok)
            res.status = DemodStatus::PayloadFailed;
        else
            res.status = DemodStatus::FcsFailed;
    }
}

} // namespace

DemodResult demodulate_one(const std::complex<float>* rx,
                                  size_t n,
                                  const Qm35825Profile& profile,
                                  uint64_t packet_id,
                                  int64_t predicted_start,
                                  int64_t window_start,
                                  const std::vector<std::complex<float>>& template_wf,
                                  DemodScratch& scratch)
{
    DemodResult res;
    res.packet_id = packet_id;
    res.predicted_start_sample = predicted_start;
    res.window_start_sample = window_start;
    res.sync_strategy = profile.sync_strategy;
    // The fused PHR decisions are packet-scoped; never let the joint payload
    // Viterbi consume a previous packet's cache.
    scratch.begin_packet();

    // predicted_start is ABSOLUTE.  The core stages index rx as a 0-based
    // window buffer, so seed with the window-relative position.  All stage
    // outputs are window-relative; `rebase` re-bases the absolute sample
    // fields (timing / SFD sample positions / CIR first path) by window_start
    // before each return, so result coordinates match the schema (absolute).
    const int64_t rel_seed =
        (predicted_start >= 0) ? predicted_start - window_start : -1;
    // Diagnostic payload sub-stage timers.  They are children of
    // `stage_payload_us` (never added to `stage_total_us`) and are published
    // from `rebase()` so every exit path reports them without touching any
    // parent lap interval.
    PayloadSubTimers payload_sub;
    auto rebase = [&]() { finalize_demod_result(res, payload_sub, window_start); };

    // Per-stage wall-clock timing.  `tick` advances a lap clock; each lap
    // records the elapsed µs into the stage field AND accumulates the running
    // total, so early-failure paths still attribute the failed stage's time.
    StageClock clock;

    // Stage 1: timing (seed = rel_seed if inside this window, else -1)
    const bool t_ok = stage_timing(rx, n, profile, template_wf, rel_seed,
                                   res.timing, scratch);
    clock.lap(res.stage_timing_us, res.stage_total_us);
    res.stage_timing_coarse_us = res.timing.coarse_search_us;
    res.stage_timing_fine_track_us = res.timing.fine_track_us;
    if (!t_ok) {
        res.status = DemodStatus::TimingFailed;
        rebase();
        return res;
    }
    SyncBootstrapState sync_state;
    if (!run_initial_cfo_and_bootstrap_sfd(
            rx, n, profile, template_wf, clock, scratch, res, sync_state)) {
        rebase();
        return res;
    }

    // Integer chip-grid admission: use the period estimate (derived from Stage-2
    // CFO when rx_sfo_from_cfo is enabled, or Stage-1 timing OLS otherwise).
    // The decision must be made here and not later, because it controls three
    // coupled things:
    //   * the second CFO rotation range (a narrowed range is only valid when
    //     the on-demand prefix is actually used),
    //   * whether the on-demand prefix is generated at all,
    //   * whether the fused PHR/payload CMF may run.
    // Mixing them -- e.g. narrowing the rotation and shortening the prefix but
    // then refusing the fused CMF -- leaves the payload decoder without PHR
    // soft chips, which fails the FCS and triggers a full-stream retry
    // (~16 ms/packet instead of ~2.6 ms).  Measured on the X410 capture with a
    // 1.5 ppm gate: 16% of packets hit exactly that path.
    const auto execution = plan_packet_execution(profile, res.timing);
    const bool integer_grid_ok = execution.integer_grid_ok;
    const auto& eff = execution.effective_profile;
    if (!integer_grid_ok)
        ++res.stage_grid_rejects;

    if (!anchor_and_verify_final_sfd(
            rx, n, eff, template_wf, integer_grid_ok,
            sync_state, clock, scratch, res)) {
        rebase();
        return res;
    }

    // Early termination if caller only requested preamble (SYNC + SFD) detection.
    if (profile.rx_stop_after_sfd) {
        res.status = DemodStatus::Success;
        rebase();
        res.stage_total_us = res.stage_timing_us + res.stage_cfo_us + res.stage_sfd_us;
        return res;
    }

    // Standalone CIR estimation after SFD (independent of CMF / despreading).
    // The estimator folds the CFO rotation into the SYNC average, so it reads
    // the raw frame and the rotation ranges above no longer need to cover the
    // preamble for it.
    if (!estimate_cir(rx, n, eff, res.timing,
                      scratch.cached_seq.pcode, res.cir, scratch,
                      res.cfo.cfo_hz)) {
        clock.lap(res.stage_cir_us, res.stage_total_us);
        res.status = DemodStatus::CirFailed;
        rebase();
        return res;
    }
    clock.lap(res.stage_cir_us, res.stage_total_us);

    if (profile.rx_stop_after_cir) {
        res.status = DemodStatus::Success;
        rebase();
        return res;
    }

    // Stage 4: soft-chip CMF + PHR + payload/FCS on the CFO-compensated frame.
    const bool on_demand_path = integer_grid_ok &&
                                (eff.rx_softchip_on_demand ||
                                 eff.rx_fused_payload_cmf ||
                                 eff.rx_fused_phr_cmf);
    if (on_demand_path) {
        run_on_demand_with_full_stream_fallback(
            rx, n, profile, eff, integer_grid_ok,
            sync_state.derotated_span, clock, payload_sub, scratch, res);
    } else {
        run_primary_full_stream_decode(
            n, profile, eff, clock, payload_sub, scratch, res);
    }
    rebase();
    return res;
}


} // namespace openuwb::demod::core
