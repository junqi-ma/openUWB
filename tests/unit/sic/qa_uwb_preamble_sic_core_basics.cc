// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 1/3 of the preamble-SIC core QA: integer/fractional/CFO and
// S1 prefix/support isolation cases.
#include "preamble_sic_test_helpers.h"

// ---------------------------------------------------------------------------
// 1. integer single path
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_integer_single_path)
{
    Fixture f =
        make_fixture(10, 32, static_cast<int64_t>(kPre) + 8, 64, 64, 4, 8, 8, 100);
    Signal s = build_signal(
        f.model, f.n, f.comm_start, { { 0.0, cplx(1.f, 0.f) } }, 0.0, kQaNoise, 1);
    SicResult r = run(f, s);
    std::printf("[sic-qa] single-path detected_start=%lld (true %lld) n_paths=%zu "
                "coarse=%llu fine=%llu corr=%llu scratch_peak_bytes=%zu\n",
                static_cast<long long>(r.detected_start),
                static_cast<long long>(f.comm_start),
                r.n_paths,
                static_cast<unsigned long long>(r.coarse_candidates),
                static_cast<unsigned long long>(r.fine_candidates),
                static_cast<unsigned long long>(r.correlations),
                f.scratch.peak_bytes);
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK(std::abs(r.detected_start - f.comm_start) <= 1);
    // S0 optimization iter-2 detection-equivalence guard: for this fixed
    // synthetic window the detection phase output is fully deterministic, so
    // pin the exact argmax and the coarse-search geometry.  An inner-loop
    // change that shifts the peak by one sample or alters the coarse search
    // fails here; the golden cross-check below independently re-checks
    // detected_start against the NumPy reference for every golden case.
    BOOST_CHECK_MESSAGE(r.detected_start == f.comm_start,
                        "detected_start " << r.detected_start << " != exact "
                                          << f.comm_start);
    BOOST_CHECK_EQUAL(r.coarse_candidates, uint64_t(2332)); // ceil(32640/14)
    const double e = nmse(s.i, f.ihat, 0, f.n);
    std::printf("[sic-qa] single-path NMSE=%.3e suppression=%.2f dB\n",
                e,
                10.0 * std::log10(1.0 / std::max(e, 1e-30)));
    BOOST_CHECK_LT(e, 1e-3);
    BOOST_CHECK(!r.uses_truth);
    BOOST_CHECK_GT(suppression_db(s.i, f.ihat, 0, f.n), 30.0);
    // S1: the support is now the CONFIGURED NOMINAL prediction (origin +
    // sync_reps*sps + post_tail), not the confirmed-repetition support, so it is
    // 32 symbols long here -- not 16 (fit+validation).
    BOOST_CHECK_EQUAL(std::string(r.support_source ? r.support_source : ""),
                      std::string("configured_nominal_prediction"));
    const int64_t nominal_len =
        static_cast<int64_t>(f.cfg.sync_reps * kSps + f.cfg.post_tail);
    BOOST_CHECK_EQUAL(r.support_hi - r.support_lo, nominal_len);
    BOOST_CHECK_GT(r.support_hi - r.support_lo,
                   static_cast<int64_t>((f.cfg.fit_reps + f.cfg.validation_reps) * kSps));
    // Outside the permitted cancellation mask out == y byte-for-byte.  With no
    // sub_roi configured the mask is support ∩ search ROI (the legacy mask
    // SHAPE, now over the nominal support).
    const int64_t mask_lo =
        std::max<int64_t>(r.support_lo, std::max<int64_t>(f.cfg.roi_lo, 0));
    const int64_t mask_hi = std::min<int64_t>(
        r.support_hi, std::min<int64_t>(static_cast<int64_t>(f.n), f.cfg.roi_hi));
    std::printf("[sic-qa] single-path support=[%lld,%lld) mask=[%lld,%lld) "
                "coverage=%.1f%%\n",
                static_cast<long long>(r.support_lo),
                static_cast<long long>(r.support_hi),
                static_cast<long long>(mask_lo),
                static_cast<long long>(mask_hi),
                100.0 * static_cast<double>(mask_hi - mask_lo) /
                    static_cast<double>(nominal_len));
    for (size_t k = 0; k < f.n; ++k) {
        if (static_cast<int64_t>(k) < mask_lo || static_cast<int64_t>(k) >= mask_hi)
            BOOST_CHECK_EQUAL(f.out[k], s.y[k]);
    }
}

// ---------------------------------------------------------------------------
// 2. integer four paths
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_integer_four_path)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = {
        { 0.0, cplx(1.f, 0.f) },
        { 7.0, cplx(0.f, 0.5f) },
        { 23.0, cplx(0.276f, -0.116f) },
        { 51.0, cplx(0.2f, 0.f) },
    };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 2);
    SicResult r = run(f, s);
    const double e = nmse(s.i, f.ihat, 0, f.n);
    std::printf("[sic-qa] four-path(state=%s) n_paths=%zu cand=%zu cond=%.2f "
                "NMSE=%.3e suppression=%.2f dB\n",
                state_name(r.state),
                r.n_paths,
                r.n_path_candidates,
                r.gain_condition,
                e,
                10.0 * std::log10(1.0 / std::max(e, 1e-30)));
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_GE(r.n_paths, size_t(4));
    BOOST_CHECK_LT(e, 0.02);
}

// ---------------------------------------------------------------------------
// 3. fractional delays (cross-model placement)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_fractional)
{
    Fixture f =
        make_fixture(10, 32, static_cast<int64_t>(kPre) + 1, 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = {
        { 0.4, cplx(1.f, 0.f) },
        { 7.4, cplx(0.f, 0.5f) },
        { 23.0, cplx(0.276f, -0.116f) },
        { 51.2, cplx(0.2f, 0.f) },
    };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 3);
    SicResult r = run(f, s);
    const double e = nmse(s.i, f.ihat, 0, f.n);
    std::printf("[sic-qa] fractional(state=%s) n_paths=%zu NMSE=%.3e "
                "suppression=%.2f dB frac_us=%llu\n",
                state_name(r.state),
                r.n_paths,
                e,
                10.0 * std::log10(1.0 / std::max(e, 1e-30)),
                static_cast<unsigned long long>(r.timings.fractional_refine_us));
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    // Cross-model (Blackman generator vs Hamming core): declared looser bound.
    BOOST_CHECK_MESSAGE(e < 0.05, "fractional cross-model NMSE too high: " << e);
    // S0 optimization iter-3 fractional-refine guard: the optional-window stage
    // now scores candidate offsets through the float32 VOLK projection path (and
    // a branch-free shifted-template build) instead of rebuilding a double
    // column per offset.  This fixed synthetic window is deterministic and the
    // double reference measured NMSE=1.309e-06, so a regression that shifts a
    // selected fractional delay by one 0.05-sample grid step fails here far
    // below the declared 0.05 cross-model bound.  detected_start / n_paths are
    // pinned separately by the golden cross-check.
    BOOST_CHECK_MESSAGE(e < 1e-4, "fractional refine NMSE regressed: " << e);
    BOOST_CHECK_EQUAL(r.n_paths, size_t(4));
}

// ---------------------------------------------------------------------------
// 4. CFO (reconstruction-basis rotation)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_cfo)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = {
        { 0.0, cplx(1.f, 0.f) },
        { 7.0, cplx(0.f, 0.5f) },
        { 23.0, cplx(0.276f, -0.116f) },
        { 51.0, cplx(0.2f, 0.f) },
    };
    const double cfo = 50000.0;
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, cfo, kQaNoise, 4);
    SicResult r = run(f, s);
    const double e = nmse(s.i, f.ihat, 0, f.n);
    std::printf("[sic-qa] cfo(state=%s) cfo_hat=%.1f Hz (true %.1f) NMSE=%.3e "
                "suppression=%.2f dB\n",
                state_name(r.state),
                r.cfo_hz,
                cfo,
                e,
                10.0 * std::log10(1.0 / std::max(e, 1e-30)));
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_LT(std::abs(r.cfo_hz - cfo), 3000.0);
    BOOST_CHECK_LT(e, 0.05);
}

// ---------------------------------------------------------------------------
// 5. K sweep
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_roi_narrow_vs_wide)
{
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) },
                                      { 11.0, cplx(0.3f, 0.1f) } };
    Fixture fn =
        make_fixture(10, 32, static_cast<int64_t>(kPre), 32, 32, 4, 8, 8, 100, false);
    Fixture fw =
        make_fixture(10, 32, static_cast<int64_t>(kPre), 32, 32, 4, 8, 8, 100, true);
    Signal sn = build_signal(fn.model, fn.n, fn.comm_start, paths, 0.0, kQaNoise, 21);
    Signal sw = build_signal(fw.model, fw.n, fw.comm_start, paths, 0.0, kQaNoise, 21);
    SicResult rn = run(fn, sn);
    SicResult rw = run(fw, sw);
    std::printf("[sic-qa] ROI narrow roi_samples=%zu coarse=%llu | wide "
                "roi_samples=%zu coarse=%llu\n",
                rn.roi_samples,
                static_cast<unsigned long long>(rn.coarse_candidates),
                rw.roi_samples,
                static_cast<unsigned long long>(rw.coarse_candidates));
    BOOST_CHECK_EQUAL(static_cast<int>(rn.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_EQUAL(static_cast<int>(rw.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_GT(rw.roi_samples, rn.roi_samples);
    BOOST_CHECK_GT(rw.coarse_candidates, rn.coarse_candidates);
}

// ---------------------------------------------------------------------------
// 7b. real comm present but OUTSIDE the sensing processing ROI -> bypass
// ---------------------------------------------------------------------------
// Regression for the R10 no_overlap false commit.  A REAL, template-consistent
// comm preamble is injected entirely OUTSIDE the sensing processing ROI, and the
// SEARCH ROI stays wide (it includes the comm) so the core still detects it.
// The subtraction/commit ROI is separate: the comm support does not overlap it,
// so the core must fail closed (Bypassed, out == y byte-for-byte).  The positive
// control moves the sensing ROI over the comm and must then commit.
BOOST_AUTO_TEST_CASE(test_sic_core_comm_outside_sensing_roi_bypasses)
{
    const size_t comm_reps = 16;
    const size_t sensing_reps = 32;
    Fixture f;
    f.model = build_model(10, comm_reps);
    const int64_t sensing_lo = static_cast<int64_t>(kPre);
    const int64_t sensing_hi = sensing_lo + static_cast<int64_t>(sensing_reps) * kSps;
    f.comm_start = sensing_hi + 400; // fully outside the sensing ROI
    f.n =
        static_cast<size_t>(f.comm_start + static_cast<int64_t>(comm_reps) * kSps + 2000);
    f.cfg.mode = SicMode::Estimate;
    f.cfg.comm_code = 10;
    f.cfg.sync_reps = comm_reps;
    f.cfg.samples_per_symbol = kSps;
    f.cfg.sample_rate = kFs;
    f.cfg.k_paths = 4;
    f.cfg.fit_reps = 8;
    f.cfg.validation_reps = 8;
    f.cfg.roi_lo = 0; // WIDE search ROI (sees the comm)
    f.cfg.roi_hi = static_cast<int64_t>(f.n);
    f.cfg.coarse_stride = 14;
    f.cfg.path_min_spacing = 0;
    f.cfg.post_tail = 100;
    f.cfg.enable_fractional = true;
    f.cfg.enable_cfo = true;
    f.cfg.detect_threshold = 0.5;
    f.cfg.validate_threshold = 0.25;
    f.cfg.sub_roi_lo = sensing_lo; // SEPARATE sensing processing ROI
    f.cfg.sub_roi_hi = sensing_hi;
    f.cfg.sub_min_overlap = 0;
    f.out.assign(f.n, cplx(0.f, 0.f));
    f.ihat.assign(f.n, cplx(0.f, 0.f));
    std::string err;
    BOOST_REQUIRE_MESSAGE(sic::prepare_sic_core(f.cfg,
                                                f.model.sync.data(),
                                                f.model.sync.size(),
                                                f.model.preamble.data(),
                                                f.model.preamble.size(),
                                                f.scratch,
                                                &err),
                          "prepare_sic_core failed: " << err);

    Signal s = build_signal(
        f.model, f.n, f.comm_start, { { 0.0, cplx(1.f, 0.f) } }, 0.0, kQaNoise, 77);
    SicResult r = run(f, s);
    std::printf("[sic-qa] outside-sensing-roi state=%s reason=%s detected=%d "
                "start=%lld sub_overlap=%lld\n",
                state_name(r.state),
                r.bypass_reason ? r.bypass_reason : "-",
                static_cast<int>(r.detected),
                static_cast<long long>(r.detected_start),
                static_cast<long long>(r.sub_overlap));
    // The comm is real: the WIDE search ROI still detects it...
    BOOST_CHECK_MESSAGE(r.detected,
                        "real comm outside the sensing ROI must still be "
                        "detected by the wide search ROI");
    BOOST_CHECK(std::abs(r.detected_start - f.comm_start) <= 1);
    // ...but it is NOT actionable: no overlap with the sensing ROI -> bypass.
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Bypassed));
    BOOST_CHECK(!r.committed);
    BOOST_CHECK_EQUAL(r.sub_overlap, int64_t(0));
    BOOST_CHECK_EQUAL(std::string(r.bypass_reason ? r.bypass_reason : ""),
                      std::string("outside_sensing_roi"));
    BOOST_CHECK_EQUAL_COLLECTIONS(f.out.begin(), f.out.end(), s.y.begin(), s.y.end());

    // Positive control: the SAME window, sensing ROI moved to cover the comm,
    // must commit and change the output.  (sub_roi is a runtime field not used
    // by prepare, so the prepared scratch is reused.)
    f.cfg.sub_roi_lo = f.comm_start - 64;
    f.cfg.sub_roi_hi = f.comm_start + static_cast<int64_t>(comm_reps) * kSps + 1000;
    SicResult r2 = run(f, s);
    std::printf("[sic-qa] overlap-control state=%s sub_overlap=%lld\n",
                state_name(r2.state),
                static_cast<long long>(r2.sub_overlap));
    BOOST_CHECK_EQUAL(static_cast<int>(r2.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK(r2.committed);
    BOOST_CHECK_GT(r2.sub_overlap, int64_t(0));
    bool changed = false;
    for (size_t k = 0; k < f.n; ++k)
        changed = changed || (f.out[k] != s.y[k]);
    BOOST_CHECK_MESSAGE(changed, "committed control must change the output");
}

// ---------------------------------------------------------------------------
// 8. wrong comm code must not commit on a code-9 signal
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_wrong_code_no_commit)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    // Signal uses code 9 but the templates are code 10.
    Model c9 = build_model(9, 32);
    Signal s;
    s.i.assign(f.n, cplx(0.f, 0.f));
    add_path(s.i, c9.preamble, f.comm_start, 0.0, cplx(1.f, 0.f), 0.0);
    s.y.assign(f.n, cplx(0.f, 0.f));
    for (size_t k = 0; k < f.n; ++k)
        s.y[k] = s.i[k];
    SicResult r = run(f, s);
    std::printf("[sic-qa] wrong-code state=%s reason=%s detected=%d\n",
                state_name(r.state),
                r.bypass_reason ? r.bypass_reason : "-",
                static_cast<int>(r.detected));
    BOOST_CHECK_NE(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK(!r.committed);
    // out must be bit-identical to y
    BOOST_CHECK_EQUAL_COLLECTIONS(f.out.begin(), f.out.end(), s.y.begin(), s.y.end());
    BOOST_CHECK(!r.uses_truth);
}

// ---------------------------------------------------------------------------
// 9. bypass bit-exactness (Off mode) + estimate never uses truth
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_bypass_bit_exact)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    Signal s = build_signal(
        f.model, f.n, f.comm_start, { { 0.0, cplx(1.f, 0.f) } }, 0.0, 1e-3, 31);
    f.cfg.mode = SicMode::Off;
    SicResult r = run(f, s);
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Bypassed));
    BOOST_CHECK(!r.committed);
    BOOST_CHECK(!r.uses_truth);
    BOOST_CHECK_EQUAL_COLLECTIONS(f.out.begin(), f.out.end(), s.y.begin(), s.y.end());
    for (size_t k = 0; k < f.n; ++k)
        BOOST_CHECK_EQUAL(f.ihat[k], cplx(0.f, 0.f));
}

// ---------------------------------------------------------------------------
// 10. oracle vs estimate truth isolation
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_oracle_truth_isolation)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) },
                                      { 7.4, cplx(0.0f, 0.5f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, 0.0, 41);

    // Estimate: never truth.
    SicResult re = run(f, s, nullptr);
    BOOST_CHECK(!re.uses_truth);
    BOOST_CHECK_EQUAL(static_cast<int>(re.state), static_cast<int>(SicState::Committed));

    // Oracle: consumes the explicit truth struct, flags uses_truth.
    SicTruth truth;
    truth.valid = true;
    truth.start_sample = static_cast<double>(f.comm_start);
    truth.cfo_hz = 0.0;
    truth.n_paths = paths.size();
    for (size_t k = 0; k < paths.size(); ++k) {
        truth.tau[k] = paths[k].tau;
        truth.gain[k] = paths[k].g;
    }
    f.cfg.mode = SicMode::Oracle;
    SicResult ro = run(f, s, &truth);
    const double eo = nmse(s.i, f.ihat, 0, f.n);
    std::printf("[sic-qa] oracle uses_truth=%d state=%s NMSE=%.3e\n",
                static_cast<int>(ro.uses_truth),
                state_name(ro.state),
                eo);
    BOOST_CHECK(ro.uses_truth);
    BOOST_CHECK_EQUAL(static_cast<int>(ro.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_LT(eo, 0.1);

    // Oracle with no truth struct must fail closed and not commit.
    SicResult rn = run(f, s, nullptr);
    BOOST_CHECK(!rn.uses_truth);
    BOOST_CHECK_NE(static_cast<int>(rn.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_EQUAL_COLLECTIONS(f.out.begin(), f.out.end(), s.y.begin(), s.y.end());
}

// ---------------------------------------------------------------------------
// 11. timings: children are not double-counted against core_total
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_timings_no_double_count)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) },
                                      { 7.4, cplx(0.0f, 0.5f) },
                                      { 23.0, cplx(0.2f, -0.1f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 25000.0, 1e-4, 51);
    SicResult r = run(f, s);
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    const sic::SicTimings& t = r.timings;
    const uint64_t children = t.detect_coarse_us + t.detect_fine_us + t.cfo_us +
                              t.cir_us + t.path_select_us + t.fractional_refine_us +
                              t.matrix_build_us + t.gain_fit_us + t.reconstruct_us +
                              t.validate_commit_us;
    std::printf("[sic-qa] timings children_sum=%llu core_total=%llu "
                "(coarse=%llu fine=%llu cfo=%llu cir=%llu sel=%llu frac=%llu "
                "mat=%llu fit=%llu rec=%llu val=%llu) | S1 estimate_total=%llu "
                "model_prepare=%llu apply=%llu backfill=%llu forward=%llu\n",
                static_cast<unsigned long long>(children),
                static_cast<unsigned long long>(t.core_total_us),
                static_cast<unsigned long long>(t.detect_coarse_us),
                static_cast<unsigned long long>(t.detect_fine_us),
                static_cast<unsigned long long>(t.cfo_us),
                static_cast<unsigned long long>(t.cir_us),
                static_cast<unsigned long long>(t.path_select_us),
                static_cast<unsigned long long>(t.fractional_refine_us),
                static_cast<unsigned long long>(t.matrix_build_us),
                static_cast<unsigned long long>(t.gain_fit_us),
                static_cast<unsigned long long>(t.reconstruct_us),
                static_cast<unsigned long long>(t.validate_commit_us),
                static_cast<unsigned long long>(t.estimate_total_us),
                static_cast<unsigned long long>(t.model_prepare_us),
                static_cast<unsigned long long>(t.apply_us),
                static_cast<unsigned long long>(t.backfill_us),
                static_cast<unsigned long long>(t.apply_forward_us));
    BOOST_CHECK_GT(t.core_total_us, uint64_t(0));
    for (uint64_t c : { t.detect_coarse_us,
                        t.detect_fine_us,
                        t.cfo_us,
                        t.cir_us,
                        t.path_select_us,
                        t.fractional_refine_us,
                        t.matrix_build_us,
                        t.gain_fit_us,
                        t.reconstruct_us,
                        t.validate_commit_us })
        BOOST_CHECK_LE(c, t.core_total_us);
    // Sequential children: their sum cannot exceed the outer total (allow a
    // couple of microseconds of per-measurement truncation rounding).
    BOOST_CHECK_LE(children, t.core_total_us + 12);
    // S1 decomposition: the two apply passes partition the window and the model
    // preparation is inside the estimation total; none of them is added to the
    // stage children.
    BOOST_CHECK_GT(t.estimate_total_us, uint64_t(0));
    BOOST_CHECK_LE(t.estimate_total_us, t.core_total_us + 12);
    BOOST_CHECK_LE(t.model_prepare_us, t.estimate_total_us + 12);
    BOOST_CHECK_LE(t.backfill_us + t.apply_forward_us, t.core_total_us + 12);
    BOOST_CHECK_LE(t.apply_us, t.backfill_us + t.apply_forward_us + 12);
    BOOST_CHECK_GT(t.reconstruct_us, uint64_t(0)); // apply reconstructs both passes
}

// ---------------------------------------------------------------------------
// 12. unexplained/negative paths: low SNR should still either cancel or
//     bypass cleanly (never a partial write)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_gain_condition_semantics)
{
    const std::vector<Path> paths = {
        { 0.0, cplx(1.f, 0.f) },
        { 7.0, cplx(0.f, 0.5f) },
        { 23.0, cplx(0.276f, -0.116f) },
        { 51.0, cplx(0.2f, 0.f) },
    };
    for (size_t K : { 1u, 4u }) {
        Fixture f =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, K, 8, 8, 100);
        Signal s =
            build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 900 + K);
        SicResult r = run(f, s);
        const double e = nmse(s.i, f.ihat, 0, f.n);
        std::printf("[sic-qa] gain-cond K=%zu state=%s n_paths=%zu cond=%.6f "
                    "fit_res=%.4f residual_corr=%.4f NMSE=%.3e\n",
                    K,
                    state_name(r.state),
                    r.n_paths,
                    r.gain_condition,
                    r.fit_residual,
                    r.residual_corr,
                    e);
        BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                          static_cast<int>(SicState::Committed));
        BOOST_CHECK(std::isfinite(r.gain_condition));
        BOOST_CHECK_GE(r.gain_condition, 1.0 - 1e-9);
        BOOST_CHECK_LT(r.gain_condition, 1e6);
        // S1 default path: no holdout wait, so residual_corr is NaN BY DESIGN and
        // the model-quality gate is the training-region fit_residual.
        BOOST_CHECK(std::isnan(r.residual_corr));
        BOOST_CHECK(std::isfinite(r.fit_residual));
        BOOST_CHECK_LE(r.fit_residual, f.cfg.fit_residual_max);
        if (K == 1) {
            BOOST_CHECK_SMALL(r.gain_condition - 1.0, 1e-9);
            BOOST_CHECK_EQUAL(r.n_paths, size_t(1));
        } else {
            BOOST_CHECK_GE(r.n_paths, size_t(4));
            BOOST_CHECK_LT(e, 0.02);
        }
    }
}

// ---------------------------------------------------------------------------
// 14. fractional selection: direct taus differential (on vs off)
// ---------------------------------------------------------------------------
// The fractional-score path was only ever validated through NMSE bounds.  This
// is a DIRECT differential: the same window is solved with enable_fractional
// true and false, and the exposed absolute fit delays (SicResult::taus) must
// agree one-for-one within the +/-0.5-sample refine bracket, with identical
// detection and path count.  A wrong fractional selection therefore shows up as
// a delay mismatch, not merely as a slightly worse reconstruction.
BOOST_AUTO_TEST_CASE(test_sic_core_fractional_selection_differential)
{
    const std::vector<Path> paths = {
        { 0.4, cplx(1.f, 0.f) },
        { 7.4, cplx(0.f, 0.5f) },
        { 23.0, cplx(0.276f, -0.116f) },
        { 51.2, cplx(0.2f, 0.f) },
    };
    Fixture f =
        make_fixture(10, 32, static_cast<int64_t>(kPre) + 1, 64, 64, 4, 8, 8, 100);
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 77);

    SicResult rf = run(f, s); // fractional on (fixture default)
    BOOST_REQUIRE_EQUAL(static_cast<int>(rf.state),
                        static_cast<int>(SicState::Committed));
    BOOST_REQUIRE_EQUAL(rf.n_taus, rf.n_paths);

    f.cfg.enable_fractional = false;
    SicResult ri = run(f, s); // integer-only on the SAME window
    f.cfg.enable_fractional = true;

    const double ef = nmse(s.i, f.ihat, 0, f.n); // f.ihat now holds the integer run
    // Re-run fractional to leave f.ihat in the fractional state for the print.
    SicResult rf2 = run(f, s);
    const double ef2 = nmse(s.i, f.ihat, 0, f.n);
    BOOST_CHECK_EQUAL(rf2.detected_start, ri.detected_start);
    BOOST_CHECK_EQUAL(rf2.n_paths, ri.n_paths);
    std::printf("[sic-qa] frac-diff detected_start=%lld n_paths=%zu (int NMSE=%.3e "
                "frac NMSE=%.3e)\n",
                static_cast<long long>(rf2.detected_start),
                rf2.n_paths,
                ef,
                ef2);
    double max_delta = 0.0;
    for (size_t k = 0; k < rf2.n_taus && k < ri.n_taus; ++k) {
        const double d = std::abs(rf2.taus[k] - ri.taus[k]);
        max_delta = std::max(max_delta, d);
        std::printf("[sic-qa] frac-diff   path %zu: int=%.4f frac=%.4f |d|=%.4f\n",
                    k,
                    ri.taus[k],
                    rf2.taus[k],
                    d);
        // The refine search is bounded to +/-0.5 sample around the integer delay.
        BOOST_CHECK_MESSAGE(d <= 0.5 + 1e-9,
                            "fractional delay " << rf2.taus[k]
                                                << " left the +/-0.5 "
                                                   "bracket of "
                                                << ri.taus[k]);
    }
    // The fractional stage must actually move at least one delay for this
    // deliberately fractional-delay signal (otherwise it silently did nothing).
    BOOST_CHECK_GT(max_delta, 0.0);
    // Integer-only must still close, and fractional must not make it worse.
    BOOST_CHECK_LT(ef, 0.05);
    BOOST_CHECK_LE(ef2, ef + 1e-9);
}

// ---------------------------------------------------------------------------
// 15. gain-fit condition number: cond2(A) is NOT the R diagonal ratio (R05)
// ---------------------------------------------------------------------------
// The commit gate is `cond < 1e6`, so `cond` needs a DEFINED meaning.  This case
// pins it to cond2(A) -- the ratio of the largest to smallest singular value,
// i.e. the NumPy `linalg.cond` quantity -- and drives detail::qr_cond_solve
// directly.  It reproduces the R05 counterexample
//
//     A = [[1, 1], [0, 1.5e-6]]
//     old max|R_kk|/min|R_kk| = 666666.7  (would PASS the gate)
//     true cond2(A)          = 1333333.3 (> 1e6, must be REJECTED)
//
// and additionally checks, for column scales / near-collinear columns / column
// permutation / both sides of the 1e6 gate / an ill-conditioned design, that
// (a) the reported cond matches an independent oracle and (b) the least-squares
// solution has a machine-level residual.
//
// The oracle is a self-contained long-double eigenvalue decomposition of the
// Gram (classical cyclic Jacobi), a DIFFERENT algorithm and a wider precision
// from the core's double one-sided complex Jacobi on the QR factor R.  It is
// validated below against designs built from a KNOWN singular-value spectrum.

// Independent reference cond2(A) = sqrt(lambda_max/lambda_min) of the Hermitian
// Gram G = A^H A, eigenvalues by the classical cyclic Jacobi method in long
// double.  (Different algorithm from the core; used only as a cross-check.)
static long double ref_cond2(const std::vector<dcp>& A, size_t M, size_t K)
{
    using lc = std::complex<long double>;
    std::vector<lc> G(K * K, lc(0, 0));
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = i; j < K; ++j) {
            long double rr = 0, ii = 0;
            for (size_t m = 0; m < M; ++m) {
                const long double ar = A[i * M + m].real();
                const long double ai = A[i * M + m].imag();
                const long double br = A[j * M + m].real();
                const long double bi = A[j * M + m].imag();
                rr += ar * br + ai * bi;
                ii += ar * bi - ai * br;
            }
            G[i * K + j] = lc(rr, ii);
            G[j * K + i] = lc(rr, -ii);
        }
    }
    for (int sweep = 0; sweep < 100; ++sweep) {
        long double off = 0;
        for (size_t p = 0; p + 1 < K; ++p) {
            for (size_t q = p + 1; q < K; ++q) {
                const lc gpq = G[p * K + q];
                const long double b = std::abs(gpq);
                const long double app = G[p * K + p].real();
                const long double aqq = G[q * K + q].real();
                off = std::max(off, b);
                if (b <= 1e-30L * std::max(1.0L, std::max(std::abs(app), std::abs(aqq))))
                    continue;
                const long double theta = 0.5L * (aqq - app);
                const long double r = std::hypot(theta, b);
                const long double np = std::hypot(b, theta + r);
                const long double nm = std::hypot(b, theta - r);
                const lc c0 = gpq / np, c2 = lc((theta + r) / np, 0);
                const lc c1 = gpq / nm, c3 = lc((theta - r) / nm, 0);
                for (size_t i = 0; i < K; ++i) { // G <- G J
                    const lc gp = G[i * K + p], gq = G[i * K + q];
                    G[i * K + p] = gp * c0 + gq * c2;
                    G[i * K + q] = gp * c1 + gq * c3;
                }
                for (size_t j = 0; j < K; ++j) { // G <- J^H G
                    const lc gp = G[p * K + j], gq = G[q * K + j];
                    G[p * K + j] = std::conj(c0) * gp + std::conj(c2) * gq;
                    G[q * K + j] = std::conj(c1) * gp + std::conj(c3) * gq;
                }
            }
        }
        if (off < 1e-25L)
            break;
    }
    long double lmin = 1e300L, lmax = 0;
    for (size_t i = 0; i < K; ++i) {
        long double l = G[i * K + i].real();
        if (l < 0)
            l = 0; // clamp a tiny negative rounding residue
        lmin = std::min(lmin, l);
        lmax = std::max(lmax, l);
    }
    return (lmin > 0) ? std::sqrt(lmax / lmin) : 1e300L;
}

// Build a column-major M x K design with EXACTLY the prescribed singular values:
// A = U diag(sigma) V^H, U (M x K) and V (K x K) with orthonormal columns.
static std::vector<dcp>
make_design(size_t M, size_t K, const std::vector<double>& sigma, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    auto orthonormalize = [&](std::vector<dcp>& B, size_t rows, size_t cols) {
        for (size_t k = 0; k < cols; ++k) {
            for (size_t j = 0; j < k; ++j) {
                dcp d(0, 0);
                for (size_t m = 0; m < rows; ++m)
                    d += std::conj(B[j * rows + m]) * B[k * rows + m];
                for (size_t m = 0; m < rows; ++m)
                    B[k * rows + m] -= d * B[j * rows + m];
            }
            double n = 0;
            for (size_t m = 0; m < rows; ++m)
                n += std::norm(B[k * rows + m]);
            n = std::sqrt(n);
            for (size_t m = 0; m < rows; ++m)
                B[k * rows + m] /= n;
        }
    };
    std::vector<dcp> U(M * K), V(K * K);
    for (auto& z : U)
        z = dcp(nd(rng), nd(rng));
    for (auto& z : V)
        z = dcp(nd(rng), nd(rng));
    orthonormalize(U, M, K);
    orthonormalize(V, K, K);
    std::vector<dcp> A(M * K, dcp(0, 0));
    for (size_t k = 0; k < K; ++k)
        for (size_t j = 0; j < K; ++j) {
            const dcp vj = std::conj(V[k * K + j]); // (V^H)[j][k]
            for (size_t m = 0; m < M; ++m)
                A[k * M + m] += U[j * M + m] * sigma[j] * vj;
        }
    return A;
}

static std::vector<dcp> make_h_true(size_t K, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<dcp> h(K);
    for (auto& z : h)
        z = dcp(nd(rng), nd(rng));
    return h;
}

struct SolveOut {
    bool ok = false;
    double cond = 0.0;
    double old_pivot_ratio = 0.0;
    double ref_cond = 0.0;
    double h_err = 0.0;
    double resid = 0.0;
};

static SolveOut run_solve_case(const std::vector<dcp>& A,
                               size_t M,
                               size_t K,
                               const std::vector<dcp>& h_true)
{
    SolveOut o;
    std::vector<dcp> Ac = A;
    std::vector<dcp> y(M), y_orig(M);
    for (size_t m = 0; m < M; ++m) {
        dcp s(0, 0);
        for (size_t k = 0; k < K; ++k)
            s += A[k * M + m] * h_true[k];
        y[m] = s;
        y_orig[m] = s;
    }
    std::vector<dcp> h(K, dcp(0, 0)), work(M, dcp(0, 0));
    o.ok = sic::detail::qr_cond_solve(
        Ac.data(), M, K, y.data(), h.data(), work.data(), &o.cond);
    o.ref_cond = static_cast<double>(ref_cond2(A, M, K));
    if (!o.ok)
        return o;
    // Old (wrong) convention for contrast: the diagonal ratio of the factor R.
    double dmax = 0.0, dmin = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < K; ++k) {
        const double d = std::abs(Ac[k * M + k]);
        dmax = std::max(dmax, d);
        dmin = std::min(dmin, d);
    }
    o.old_pivot_ratio =
        (dmin > 0.0) ? dmax / dmin : std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < K; ++k)
        o.h_err = std::max(o.h_err, std::abs(h[k] - h_true[k]));
    double rn = 0.0, yn = 0.0;
    for (size_t m = 0; m < M; ++m) {
        dcp s(0, 0);
        for (size_t k = 0; k < K; ++k)
            s += A[k * M + m] * h[k];
        rn += std::norm(s - y_orig[m]);
        yn += std::norm(y_orig[m]);
    }
    o.resid = (yn > 0.0) ? std::sqrt(rn / yn) : 0.0;
    return o;
}

BOOST_AUTO_TEST_CASE(test_sic_core_gain_fit_ill_conditioned)
{
    // (a) exact R05 counterexample.
    {
        const size_t M = 2, K = 2;
        std::vector<dcp> A(M * K, dcp(0, 0));
        A[0 * M + 0] = dcp(1.0, 0.0);
        A[1 * M + 0] = dcp(1.0, 0.0);
        A[1 * M + 1] = dcp(1.5e-6, 0.0);
        const std::vector<dcp> ht = { dcp(0.7, -0.2), dcp(0.3, 0.45) };
        SolveOut o = run_solve_case(A, M, K, ht);
        std::printf("[sic-qa] R05 counterexample A=[[1,1],[0,1.5e-6]]: "
                    "old_ratio=%.6g (PASSES <1e6), cond2(A)=%.6g (REJECTS >=1e6), "
                    "ref=%.6g, resid=%.2e\n",
                    o.old_pivot_ratio,
                    o.cond,
                    o.ref_cond,
                    o.resid);
        BOOST_REQUIRE(o.ok);
        // The rejected old-vs-true gap: the diagonal ratio would let this
        // through the gate, the true cond(A) must not.
        BOOST_CHECK_LT(o.old_pivot_ratio, 1e6);
        BOOST_CHECK_GE(o.cond, 1e6);
        BOOST_CHECK_GT(o.cond / o.old_pivot_ratio, 1.5);
        BOOST_CHECK_CLOSE(o.cond, 1333333.333, 1.0);
        BOOST_CHECK_CLOSE(o.ref_cond, 1333333.333, 1.0);
        BOOST_CHECK_SMALL(o.resid, 1e-9);
    }

    // (b) the same counterexample embedded in a realistic-size design: two
    //     columns u and u + eps*w with w orthogonal to u (M = 4096).
    {
        const size_t M = 4096, K = 2;
        const double eps = 1.5e-6;
        std::mt19937_64 rng(0xC0FFEEu);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<dcp> u(M), w(M);
        for (auto& z : u)
            z = dcp(nd(rng), nd(rng));
        for (auto& z : w)
            z = dcp(nd(rng), nd(rng));
        auto dotc = [&](const std::vector<dcp>& a, const std::vector<dcp>& b) {
            dcp s(0, 0);
            for (size_t i = 0; i < M; ++i)
                s += std::conj(a[i]) * b[i];
            return s;
        };
        const double nu = std::sqrt(dotc(u, u).real());
        for (auto& z : u)
            z /= nu;
        const dcp uw = dotc(u, w);
        for (size_t i = 0; i < M; ++i)
            w[i] -= uw * u[i];
        const double nw = std::sqrt(dotc(w, w).real());
        for (auto& z : w)
            z /= nw;
        std::vector<dcp> A(M * K, dcp(0, 0));
        for (size_t i = 0; i < M; ++i) {
            A[0 * M + i] = u[i];
            A[1 * M + i] = u[i] + eps * w[i];
        }
        const std::vector<dcp> ht = { dcp(0.7, -0.2), dcp(0.3, 0.45) };
        SolveOut o = run_solve_case(A, M, K, ht);
        std::printf("[sic-qa] R05 embedded (M=4096, near-collinear): "
                    "old_ratio=%.6g cond2=%.6g ref=%.6g gate=%s\n",
                    o.old_pivot_ratio,
                    o.cond,
                    o.ref_cond,
                    (o.cond < 1e6) ? "accept" : "reject");
        BOOST_REQUIRE(o.ok);
        BOOST_CHECK_LT(o.old_pivot_ratio, 1e6);
        BOOST_CHECK_GE(o.cond, 1e6);
        BOOST_CHECK_CLOSE(o.cond, o.ref_cond, 5.0);
        BOOST_CHECK_SMALL(o.resid, 1e-8);
    }

    // (c) known-spectrum designs: both sides of the 1e6 gate + a clearly
    //     ill-conditioned one.  cond_target is exact by construction.
    struct Case {
        const char* name;
        std::vector<double> sigma;
        double cond_target;
        bool expect_commit;
        double resid_tol;
    };
    const std::vector<Case> cases = {
        { "well-conditioned", { 1.0, 0.5, 0.25, 0.125 }, 8.0, true, 1e-13 },
        { "gate just under (5e5)", { 1.0, 1.0, 1.0, 2e-6 }, 5e5, true, 1e-8 },
        { "gate just over (2e6)", { 1.0, 1.0, 1.0, 5e-7 }, 2e6, false, 1e-7 },
        { "ill-conditioned 1e8", { 1.0, 1.0, 1.0, 1e-8 }, 1e8, false, 1e-5 },
    };
    for (const Case& c : cases) {
        const size_t M = 2048, K = c.sigma.size();
        const std::vector<dcp> A = make_design(M, K, c.sigma, 0xABCDEFull + K);
        const std::vector<dcp> ht = make_h_true(K, 0x1234567ull + K);
        SolveOut o = run_solve_case(A, M, K, ht);
        double hmax = 0.0;
        for (const dcp& z : ht)
            hmax = std::max(hmax, std::abs(z));
        std::printf("[sic-qa] cond-case '%s': cond=%.6g (target %.3g, rel %.1e) "
                    "ref=%.6g gate=%s h_err=%.2e resid=%.2e\n",
                    c.name,
                    o.cond,
                    c.cond_target,
                    std::abs(o.cond - c.cond_target) / c.cond_target,
                    o.ref_cond,
                    (o.cond < 1e6) ? "accept" : "reject",
                    o.h_err,
                    o.resid);
        BOOST_REQUIRE_MESSAGE(o.ok, c.name << ": qr_cond_solve failed");
        BOOST_CHECK_GE(o.cond, 1.0 - 1e-12);
        BOOST_CHECK_CLOSE(o.cond, c.cond_target, 1.0);
        BOOST_CHECK_CLOSE(o.cond, o.ref_cond, 1.0);
        BOOST_CHECK_EQUAL((o.cond < 1e6), c.expect_commit);
        BOOST_CHECK_SMALL(o.resid, c.resid_tol);
        // Forward error of the solved gains must obey the standard
        // O(cond * eps * ||h||) QR least-squares bound (a grossly wrong
        // solution, or a cond that did not track the design, fails here).
        BOOST_CHECK_LE(o.h_err, 1000.0 * o.cond * 2.3e-16 * hmax);
    }

    // (d) column permutation must not change cond2(A) (the SVD is invariant);
    //     the old diagonal ratio is NOT invariant, which is the point.
    {
        const size_t M = 2048, K = 4;
        const std::vector<double> sigma = { 1.0, 1.0, 1.0, 1e-8 };
        const std::vector<dcp> A = make_design(M, K, sigma, 0x55AAull);
        std::vector<dcp> P(M * K, dcp(0, 0));
        for (size_t k = 0; k < K; ++k)
            for (size_t m = 0; m < M; ++m)
                P[k * M + m] = A[(K - 1 - k) * M + m]; // reverse column order
        const std::vector<dcp> ht = make_h_true(K, 77);
        SolveOut oa = run_solve_case(A, M, K, ht);
        SolveOut op = run_solve_case(P, M, K, ht);
        std::printf("[sic-qa] cond permutation: cond=%.6g -> %.6g (rel %.1e), "
                    "old_ratio=%.6g -> %.6g\n",
                    oa.cond,
                    op.cond,
                    std::abs(oa.cond - op.cond) / oa.cond,
                    oa.old_pivot_ratio,
                    op.old_pivot_ratio);
        BOOST_REQUIRE(oa.ok && op.ok);
        BOOST_CHECK_CLOSE(oa.cond, op.cond, 1e-6);
        BOOST_CHECK_CLOSE(oa.cond, oa.ref_cond, 1.0);
        BOOST_CHECK_CLOSE(op.cond, op.ref_cond, 1.0);
    }

    // (e) different column scales: scaling columns by s_k changes cond2(A);
    //     the reported value must still track the independent oracle.
    {
        const size_t M = 2048, K = 4;
        const std::vector<double> sigma = { 1.0, 0.5, 0.25, 0.125 };
        const std::vector<dcp> A = make_design(M, K, sigma, 0x9900ull);
        const double sc[4] = { 1e3, 1.0, 1e-2, 10.0 };
        std::vector<dcp> S(M * K, dcp(0, 0));
        for (size_t k = 0; k < K; ++k)
            for (size_t m = 0; m < M; ++m)
                S[k * M + m] = A[k * M + m] * sc[k];
        const std::vector<dcp> ht = make_h_true(K, 91);
        SolveOut oa = run_solve_case(A, M, K, ht);
        SolveOut os = run_solve_case(S, M, K, ht);
        std::printf("[sic-qa] cond column-scales: cond=%.6g -> %.6g, ref=%.6g "
                    "(rel %.1e)\n",
                    oa.cond,
                    os.cond,
                    os.ref_cond,
                    std::abs(os.cond - os.ref_cond) / os.ref_cond);
        BOOST_REQUIRE(oa.ok && os.ok);
        BOOST_CHECK_CLOSE(os.cond, os.ref_cond, 1.0);
        BOOST_CHECK_NE(oa.cond, os.cond); // scaling really changed cond
    }
}

// ---------------------------------------------------------------------------
// 16. adjacent-path selection with path_min_spacing = 1 and post_tail >= sps
// ---------------------------------------------------------------------------
// Coverage for the largest-post_tail / minimum-spacing geometry: the periodic
// CIR peak one SYNC period after the origin is selectable, so the fit runs with
// closely spaced (conditioning-relevant) columns.  The core must commit, keep
// the gate consistent with its own reported cond, and expose taus == n_paths.
BOOST_AUTO_TEST_CASE(test_sic_core_adjacent_path_conditioning)
{
    const size_t pt = kSps + 84;
    const std::vector<Path> paths = {
        { 0.0, cplx(1.f, 0.f) },
        { static_cast<double>(kSps), cplx(1.f, 0.f) },
    };
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, pt);
    f.cfg.path_min_spacing = 1;
    std::string err;
    BOOST_REQUIRE_MESSAGE(sic::prepare_sic_core(f.cfg,
                                                f.model.sync.data(),
                                                f.model.sync.size(),
                                                f.model.preamble.data(),
                                                f.model.preamble.size(),
                                                f.scratch,
                                                &err),
                          "adjacent-path prepare failed: " << err);
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 313);
    SicResult r = run(f, s);
    std::printf("[sic-qa] adjacent-path state=%s n_paths=%zu n_cand=%zu cond=%.6g "
                "fit_res=%.4f res_corr=%.4f\n",
                state_name(r.state),
                r.n_paths,
                r.n_path_candidates,
                r.gain_condition,
                r.fit_residual,
                r.residual_corr);
    for (size_t k = 0; k < r.n_taus; ++k)
        std::printf("[sic-qa] adjacent-path   tau[%zu]=%.4f\n", k, r.taus[k]);
    BOOST_CHECK_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    BOOST_CHECK_EQUAL(r.n_taus, r.n_paths);
    BOOST_CHECK(std::isfinite(r.gain_condition));
    BOOST_CHECK_GE(r.gain_condition, 1.0 - 1e-9);
    // Committed implies the default S1 gate passed on the core's OWN reported
    // cond and training fit residual.
    BOOST_CHECK_LT(r.gain_condition, 1e6);
    BOOST_CHECK_LE(r.fit_residual, f.cfg.fit_residual_max);
    BOOST_CHECK(std::isnan(r.residual_corr)); // no holdout in the default path
}

// ---------------------------------------------------------------------------
// 17. reconstruction window low edge for a negative relative delay (F7)
// ---------------------------------------------------------------------------
// path_select returns paths magnitude-sorted, so tau_abs[0] is not necessarily
// the smallest delay; the reconstruction window must start at the earliest
// SELECTED delay (union_lo = min_k floor(tau_k)), matching the reference, which
// builds i_hat over the full window.  Two paths a few samples apart with
// opposite signs (and path_min_spacing = 1) make the earlier arrival a distinct
// CIR peak, so a selected path has a NEGATIVE delay relative to the detected
// origin and i_hat must cover it.
//
// Honest scope note: the pre-fix `union_lo = floor(tau_abs[0])` could not be
// made to differ from the fixed value on any geometry reachable through this
// API -- the mainlobe-spaced CIR peak selection makes the strongest arrival the
// earliest selected peak (so tau_abs[0] == min tau), and the later
// `min(union_lo, sub_lo)` clamp masks the rest.  This is therefore a low-edge
// coverage/regression guard, not a reproduction of the old line's failure.
BOOST_AUTO_TEST_CASE(test_sic_core_union_low_edge_negative_delay)
{
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) },
                                      { 3.0, cplx(-0.6f, 0.f) } };
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    f.cfg.path_min_spacing = 1;
    std::string err;
    BOOST_REQUIRE_MESSAGE(sic::prepare_sic_core(f.cfg,
                                                f.model.sync.data(),
                                                f.model.sync.size(),
                                                f.model.preamble.data(),
                                                f.model.preamble.size(),
                                                f.scratch,
                                                &err),
                          "union-low-edge prepare failed: " << err);
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 61);
    SicResult r = run(f, s);
    double tmin = r.taus[0];
    for (size_t k = 1; k < r.n_taus; ++k)
        tmin = std::min(tmin, r.taus[k]);
    std::printf("[sic-qa] union-low-edge state=%s ds=%lld np=%zu tmin=%.4f "
                "ihat[floor(tmin)]=%s\n",
                state_name(r.state),
                static_cast<long long>(r.detected_start),
                r.n_taus,
                tmin,
                f.ihat[static_cast<size_t>(std::floor(tmin))] == cplx(0.f, 0.f)
                    ? "zero"
                    : "nonzero");
    for (size_t k = 0; k < r.n_taus; ++k)
        std::printf("[sic-qa] union-low-edge   tau[%zu]=%.4f\n", k, r.taus[k]);
    BOOST_REQUIRE_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    BOOST_REQUIRE_GT(r.n_taus, size_t(1));
    // A selected path has a negative delay relative to the detected origin.
    BOOST_REQUIRE_LT(tmin, static_cast<double>(r.detected_start));
    // The reconstruction buffer must cover it (F7: union_lo = min_k floor(tau_k)).
    BOOST_CHECK(f.ihat[static_cast<size_t>(std::floor(tmin))] != cplx(0.f, 0.f));
}

// ---------------------------------------------------------------------------
// 18. energy_prefix fallback + empty-prefix prepare (F3)
// ---------------------------------------------------------------------------
// prepare_sic_core sizes the |y|^2 prefix from the ROI it is given.  A runtime
// ROI that reaches past it must fall back to the per-lag VOLK self-dot and still
// produce the SAME detection as a run whose prefix covers the whole search; and
// preparing with roi 0,0 (prefix disabled) with a full-window runtime ROI must
// agree with it too.
BOOST_AUTO_TEST_CASE(test_sic_s1_prefix_isolation)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) }, { 7.4, cplx(0.f, 0.5f) } };
    Signal a = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 401);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;

    SicResult ra;
    SicModel ma = estimate_model_checked(f, a, oe, &ra);
    BOOST_REQUIRE_EQUAL(static_cast<int>(ra.state),
                        static_cast<int>(SicState::Committed));
    BOOST_REQUIRE_EQUAL(ma.estimation_ready_sample, oe);
    BOOST_REQUIRE_EQUAL(ma.training_lo, f.comm_start);
    BOOST_REQUIRE_EQUAL(ma.training_hi, oe);

    // (a) large additive mutation of the whole tail
    Signal b = a;
    for (size_t k = static_cast<size_t>(oe); k < b.y.size(); ++k)
        b.y[k] = b.y[k] + cplx(7.0f, -3.0f);
    SicModel mb = estimate_model_checked(f, b, oe, nullptr);
    check_models_identical(ma, mb, "prefix-additive-mutation");

    // (b) zeroed tail
    Signal c = a;
    for (size_t k = static_cast<size_t>(oe); k < c.y.size(); ++k)
        c.y[k] = cplx(0.f, 0.f);
    SicModel mc = estimate_model_checked(f, c, oe, nullptr);
    check_models_identical(ma, mc, "prefix-zeroed-tail");

    // (c) a wholly different tail signal (noise only)
    Signal d = a;
    std::mt19937_64 rng(402);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (size_t k = static_cast<size_t>(oe); k < d.y.size(); ++k)
        d.y[k] = cplx(static_cast<float>(nd(rng)), static_cast<float>(nd(rng)));
    SicModel md = estimate_model_checked(f, d, oe, nullptr);
    check_models_identical(ma, md, "prefix-unrelated-tail");

    std::printf("[sic-qa] s1-prefix-isolation: origin=%lld ready=%lld "
                "training=[%lld,%lld) model_id=%016llx 3 tail variants identical\n",
                static_cast<long long>(ma.origin),
                static_cast<long long>(ma.estimation_ready_sample),
                static_cast<long long>(ma.training_lo),
                static_cast<long long>(ma.training_hi),
                static_cast<unsigned long long>(ma.model_id));
}

// ---------------------------------------------------------------------------
// S1-2. CFO isolation: a post-prefix phase/CFO discontinuity leaves the training
// estimate untouched, while the independent-region effect degrades realistically.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_s1_cfo_isolation)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
    Signal a = build_signal(f.model, f.n, f.comm_start, paths, 25000.0, kQaNoise, 403);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;

    SicModel ma = estimate_model_checked(f, a, oe, nullptr);

    // Tail gets an extra +80 kHz rotation on top of the +25 kHz carrier: the
    // prefix (used for the estimate) is byte-identical, the tail is not.  The
    // modified interferer and the SAME noise realisation are used to rebuild y,
    // so the prefix bytes are exactly a.y while i_true stays consistent.
    Signal b = a;
    for (size_t k = static_cast<size_t>(oe); k < b.y.size(); ++k) {
        const double ph =
            2.0 * kPi * 80000.0 * static_cast<double>(k - static_cast<size_t>(oe)) / kFs;
        const cplx rot(static_cast<float>(std::cos(ph)),
                       static_cast<float>(std::sin(ph)));
        const cplx nse = a.y[k] - a.i[k];
        b.i[k] = a.i[k] * rot;
        b.y[k] = b.i[k] + nse;
    }
    SicModel mb = estimate_model_checked(f, b, oe, nullptr);
    check_models_identical(ma, mb, "cfo-isolation");

    // Effect over the independent region (reps 8..64, past the pulse tail).
    const size_t ind_lo = static_cast<size_t>(oe) + 200;
    const size_t ind_hi = static_cast<size_t>(f.comm_start) + 64 * kSps;
    SicResult r = run(f, a);
    BOOST_REQUIRE_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    const double supp_a = actual_suppression_db(a.i, a.y, f.out, ind_lo, ind_hi);
    SicResult rb = run(f, b);
    BOOST_REQUIRE_EQUAL(static_cast<int>(rb.state),
                        static_cast<int>(SicState::Committed));
    const double supp_b = actual_suppression_db(b.i, b.y, f.out, ind_lo, ind_hi);
    std::printf("[sic-qa] s1-cfo-isolation: models identical; independent-region "
                "actual suppression clean=%.2f dB tail-phase-jump=%.2f dB\n",
                supp_a,
                supp_b);
    BOOST_CHECK_GT(supp_a, 20.0);
    BOOST_CHECK_LT(supp_b, supp_a - 10.0);
}

// ---------------------------------------------------------------------------
// S1-3 / S1-4. Full-nominal support + tail actual subtraction.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_s1_full_support_and_tail)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 404);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;

    SicResult r = run(f, s);
    BOOST_REQUIRE_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    const int64_t nominal_len =
        static_cast<int64_t>(f.cfg.sync_reps * kSps + f.cfg.post_tail);
    BOOST_CHECK_EQUAL(std::string(r.support_source ? r.support_source : ""),
                      std::string("configured_nominal_prediction"));
    BOOST_CHECK_EQUAL(r.support_hi - r.support_lo, nominal_len);
    BOOST_CHECK_GT(r.support_hi - r.support_lo, 60 * static_cast<int64_t>(kSps));
    // Coverage of the TRUE interferer support (its extent is primer tail included).
    const int64_t true_hi =
        std::min<int64_t>(static_cast<int64_t>(f.n),
                          f.comm_start + static_cast<int64_t>(f.model.preamble.size()));
    const int64_t cov_lo = std::max(r.support_lo, f.comm_start);
    const int64_t cov_hi = std::min(r.support_hi, true_hi);
    const double coverage = static_cast<double>(cov_hi - cov_lo) /
                            static_cast<double>(true_hi - f.comm_start);

    const size_t ind_lo = static_cast<size_t>(oe);
    const size_t ind_hi = static_cast<size_t>(f.comm_start) + 64 * kSps;
    const size_t t4_lo = static_cast<size_t>(f.comm_start) + 60 * kSps;
    const size_t t4_hi = static_cast<size_t>(f.comm_start) + 64 * kSps;
    const double supp_ind = actual_suppression_db(s.i, s.y, f.out, ind_lo, ind_hi);
    const double supp_all = actual_suppression_db(s.i, s.y, f.out, 0, f.n);
    const double supp_tail = actual_suppression_db(s.i, s.y, f.out, t4_lo, t4_hi);

    // Differential against the legacy confirmed-support mode: it can only cover
    // confirmed (fit+online=16) symbols, so the last 4 SYNC stay untouched.
    f.cfg.support_mode = SicSupportMode::ConfirmedLegacy;
    f.cfg.online_validation_reps = f.cfg.validation_reps;
    run(f, s);
    const double supp_tail_legacy = actual_suppression_db(s.i, s.y, f.out, t4_lo, t4_hi);
    const double ei_tail = region_energy(s.i, t4_lo, t4_hi);

    std::printf("[sic-qa] s1-full-support: support=[%lld,%lld) nominal_len=%lld "
                "coverage=%.4f | actual suppression all=%.2f dB indep=%.2f dB "
                "tail4=%.2f dB | legacy tail4=%.2f dB (raw tail E_i=%.6g)\n",
                static_cast<long long>(r.support_lo),
                static_cast<long long>(r.support_hi),
                static_cast<long long>(nominal_len),
                coverage,
                supp_all,
                supp_ind,
                supp_tail,
                supp_tail_legacy,
                ei_tail);
    BOOST_CHECK_GE(coverage, 0.99);
    BOOST_CHECK_GT(supp_ind, 20.0);
    BOOST_CHECK_GT(supp_tail, 20.0);
    std::printf("[sic-qa] s1-full-support timing: core_total=%llu "
                "estimate_total=%llu detect=%llu/%llu cfo=%llu cir=%llu "
                "sel=%llu frac=%llu mat=%llu fit=%llu model_prepare=%llu "
                "apply=%llu backfill=%llu forward=%llu recon=%llu "
                "scratch_peak=%zu bytes\n",
                static_cast<unsigned long long>(r.timings.core_total_us),
                static_cast<unsigned long long>(r.timings.estimate_total_us),
                static_cast<unsigned long long>(r.timings.detect_coarse_us),
                static_cast<unsigned long long>(r.timings.detect_fine_us),
                static_cast<unsigned long long>(r.timings.cfo_us),
                static_cast<unsigned long long>(r.timings.cir_us),
                static_cast<unsigned long long>(r.timings.path_select_us),
                static_cast<unsigned long long>(r.timings.fractional_refine_us),
                static_cast<unsigned long long>(r.timings.matrix_build_us),
                static_cast<unsigned long long>(r.timings.gain_fit_us),
                static_cast<unsigned long long>(r.timings.model_prepare_us),
                static_cast<unsigned long long>(r.timings.apply_us),
                static_cast<unsigned long long>(r.timings.backfill_us),
                static_cast<unsigned long long>(r.timings.apply_forward_us),
                static_cast<unsigned long long>(r.timings.reconstruct_us),
                f.scratch.peak_bytes);
    // The legacy differential PROVES the tail was previously not subtracted.
    BOOST_CHECK_LT(supp_tail_legacy, 3.0);
    BOOST_CHECK_GT(supp_tail - supp_tail_legacy, 15.0);
    // Raw arrays for an independent recompute: out == y - i_hat within the mask.
    {
        int64_t tested = 0;
        double worst = 0.0;
        for (int64_t k = t4_lo; k < static_cast<int64_t>(t4_hi); ++k) {
            const cplx rhs = s.y[static_cast<size_t>(k)] - f.ihat[static_cast<size_t>(k)];
            worst = std::max(
                worst,
                static_cast<double>(std::abs(f.out[static_cast<size_t>(k)] - rhs)));
            ++tested;
        }
        std::printf("[sic-qa] s1-full-support: out==y-i_hat over %lld tail samples "
                    "(max|diff|=%.3e)\n",
                    static_cast<long long>(tested),
                    worst);
        BOOST_CHECK_LT(worst, 2e-6);
    }

    // Nominal-load reality check (MEASURED, not asserted as a target): with a
    // co-channel sensing signal at ISR=10 dB and SNR=30 dB, the 8-symbol CFO
    // estimate is perturbed and the far-tail depth becomes CFO-limited.  It is
    // reported next to the clean-environment depth so the SUPPORT mechanism
    // (proven above) is not confused with the CFO-limited DEPTH.
    {
        Fixture g =
            make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
        Signal s2 = build_signal(g.model, g.n, g.comm_start, paths, 0.0, 0.0, 505);
        const double ei = region_energy(s2.i, 0, g.n);
        const double e_tgt = ei / 10.0; // total target: ISR 10 dB
        std::mt19937_64 rng(506);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<cplx> tgt(g.n);
        double e0 = 0.0;
        for (size_t k = 0; k < g.n; ++k) {
            tgt[k] = cplx(static_cast<float>(nd(rng)), static_cast<float>(nd(rng)));
            e0 += std::norm(tgt[k]);
        }
        const double sc = std::sqrt(e_tgt / std::max(e0, 1e-30));
        // SNR = 30 dB: total noise energy = e_tgt/1000, split over 2 components.
        const double sigma = std::sqrt(e_tgt / (2000.0 * static_cast<double>(g.n)));
        for (size_t k = 0; k < g.n; ++k) {
            tgt[k] = cplx(static_cast<float>(tgt[k].real() * sc),
                          static_cast<float>(tgt[k].imag() * sc));
            s2.y[k] = s2.i[k] + tgt[k] +
                      cplx(static_cast<float>(nd(rng) * sigma),
                           static_cast<float>(nd(rng) * sigma));
        }
        const int64_t oe2 = g.comm_start + 8 * static_cast<int64_t>(kSps);
        g.cfg.observed_end_sample = oe2;
        SicResult rn = run(g, s2);
        const double tail_n = actual_suppression_db(s2.i, s2.y, g.out, t4_lo, t4_hi);
        const double ind_n = actual_suppression_db(s2.i, s2.y, g.out, ind_lo, ind_hi);
        std::printf("[sic-qa] s1-full-support nominal-load: state=%s reason=%s "
                    "cfo_hat=%.0f Hz "
                    "fit_res=%.4f indep=%.2f dB tail4=%.2f dB "
                    "(8-symbol CFO limited; clean env above)\n",
                    state_name(rn.state),
                    rn.bypass_reason ? rn.bypass_reason : "-",
                    rn.cfo_hz,
                    rn.fit_residual,
                    ind_n,
                    tail_n);
        BOOST_CHECK_EQUAL(static_cast<int>(rn.state),
                          static_cast<int>(SicState::Committed));
        BOOST_CHECK_GT(tail_n, 0.0); // actually subtracted, not ~0 dB
    }
}

// ---------------------------------------------------------------------------
// S1-5. Chunk invariance: whole-window vs fixed chunks vs random chunks produce
// bit-identical out and i_hat (absolute phase anchoring, no cross-chunk state).
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_s1_chunk_invariance)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) },
                                      { 11.4, cplx(0.4f, -0.2f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 30000.0, kQaNoise, 405);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;
    SicModel m = estimate_model_checked(f, s, oe, nullptr);
    BOOST_REQUIRE(m.valid);

    std::vector<cplx> out_ref, ih_ref;
    apply_chunks(f, m, s, { { 0, s.y.size() } }, out_ref, ih_ref);

    // sic_one (estimate + two-pass backfill/forward) must agree bit-for-bit.
    SicResult r = run(f, s);
    BOOST_REQUIRE_EQUAL(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
    size_t diff_one = 0;
    for (size_t k = 0; k < s.y.size(); ++k)
        diff_one += (f.out[k] != out_ref[k]) ? 1 : 0;

    // Fixed chunks aligned to the training boundary + odd sizes.
    std::vector<std::pair<int64_t, size_t>> fixed;
    for (int64_t p = 0; p < static_cast<int64_t>(s.y.size());) {
        const size_t len = static_cast<size_t>(std::min<int64_t>(
            1500 + (p % 3) * 137, static_cast<int64_t>(s.y.size()) - p));
        fixed.push_back({ p, len });
        p += static_cast<int64_t>(len);
    }
    std::vector<cplx> out_fix, ih_fix;
    apply_chunks(f, m, s, fixed, out_fix, ih_fix);

    // Random chunks.
    std::vector<std::pair<int64_t, size_t>> rnd;
    std::mt19937_64 rng(406);
    for (int64_t p = 0; p < static_cast<int64_t>(s.y.size());) {
        const size_t len =
            static_cast<size_t>(std::min<int64_t>(1 + static_cast<int64_t>(rng() % 4000),
                                                  static_cast<int64_t>(s.y.size()) - p));
        rnd.push_back({ p, len });
        p += static_cast<int64_t>(len);
    }
    std::vector<cplx> out_rnd, ih_rnd;
    apply_chunks(f, m, s, rnd, out_rnd, ih_rnd);

    size_t diff_f = 0, diff_r = 0, diff_ih_f = 0, diff_ih_r = 0;
    for (size_t k = 0; k < s.y.size(); ++k) {
        diff_f += (out_fix[k] != out_ref[k]) ? 1 : 0;
        diff_r += (out_rnd[k] != out_ref[k]) ? 1 : 0;
        diff_ih_f += (ih_fix[k] != ih_ref[k]) ? 1 : 0;
        diff_ih_r += (ih_rnd[k] != ih_ref[k]) ? 1 : 0;
    }
    std::printf("[sic-qa] s1-chunk-invariance: n=%zu chunks fixed=%zu random=%zu | "
                "out diffs: whole-vs-sic_one=%zu fixed=%zu random=%zu | "
                "i_hat diffs: fixed=%zu random=%zu\n",
                s.y.size(),
                fixed.size(),
                rnd.size(),
                diff_one,
                diff_f,
                diff_r,
                diff_ih_f,
                diff_ih_r);
    // Documented tolerance: bit-identical (0 differing samples).  If a platform
    // ever fails this, the printed counts are the measured deviation.
    BOOST_CHECK_EQUAL(diff_one, size_t(0));
    BOOST_CHECK_EQUAL(diff_f, size_t(0));
    BOOST_CHECK_EQUAL(diff_r, size_t(0));
    BOOST_CHECK_EQUAL(diff_ih_f, size_t(0));
    BOOST_CHECK_EQUAL(diff_ih_r, size_t(0));
}

// ---------------------------------------------------------------------------
// S1-6. Preamble boundary: start +/- offsets, fractional delay, finite template
// tail and window clipping; and the nominal-vs-measured support distinction.
// ---------------------------------------------------------------------------
