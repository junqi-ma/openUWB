// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 1/4 of the split realtime-demodulator core QA: R1 timing/CFO/SFD.
// Golden reference: testdata/realtime_demod_golden/ (see the helpers header).
#include "demod_core_test_helpers.h"

// ---------------------------------------------------------------------------
// Integer chip-grid admission (fused CMF / AVX2 soft-chip kernel).
//
// The grid tolerance is derived from the worst-case payload drift, not from an
// arbitrary epsilon.  These cases pin the two ends that matter on real
// hardware: a 27 ppm clock must still be rejected, a sub-ppm clock accepted,
// and the boundary must sit at the documented drift budget.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(integer_chip_grid_tolerance_matches_the_drift_budget)
{
    using namespace openuwb::demod::core;

    // Default admission is 1.5 ppm = 0.23 samples of drift over a 127 B
    // payload.  The hard ceiling is 1 sample = 6.468 ppm.
    BOOST_CHECK_CLOSE(kIntegerGridMaxClockPpm, 1.5, 1e-9);
    BOOST_CHECK_CLOSE(kIntegerGridCeilingPpm, 6.468, 0.1);
    BOOST_CHECK_CLOSE(kIntegerGridRelTolerance, 1.5e-6, 1e-9);
    BOOST_CHECK_CLOSE(kIntegerGridMaxDriftSamples, 0.232, 0.5);

    // Exactly 2 samples/chip is always accepted.
    BOOST_CHECK(integer_chip_grid(2.0));
    BOOST_CHECK(integer_chip_grid_at(2.0, 2));

    // 0.66 ppm (the best clock measured on hardware) -> 0.10 samples of
    // drift -> accepted with 2.3x margin.
    const double spc_066 = 2.0 * (1.0 - 0.66e-6);
    BOOST_CHECK(integer_chip_grid(spc_066));
    BOOST_CHECK(integer_chip_grid_at(spc_066, 2));

    // 1.24 ppm (an earlier trim step) -> 0.19 samples -> still accepted.
    BOOST_CHECK(integer_chip_grid_at(2.0 * (1.0 - 1.24e-6), 2));

    // Just inside / just outside the default gate.
    BOOST_CHECK(integer_chip_grid_at(2.0 * (1.0 - 1.4e-6), 2));
    BOOST_CHECK(!integer_chip_grid_at(2.0 * (1.0 - 1.6e-6), 2));

    // 27 ppm (the untrimmed clock) -> 4.2 samples of drift -> rejected at
    // every tolerance up to the physical ceiling.
    BOOST_CHECK(!integer_chip_grid_at(2.0 * (1.0 - 27.0e-6), 2));
    BOOST_CHECK(!integer_chip_grid(2.0 * (1.0 - 27.0e-6)));
    BOOST_CHECK(
        !integer_chip_grid_at(2.0 * (1.0 - 27.0e-6), 2, kIntegerGridCeilingPpm * 1e-6));

    // An explicit tolerance is honoured: the physical ceiling admits a 6 ppm
    // clock that the default gate refuses.
    const double spc_6 = 2.0 * (1.0 - 6.0e-6);
    BOOST_CHECK(!integer_chip_grid(spc_6));
    BOOST_CHECK(integer_chip_grid(spc_6, kIntegerGridCeilingPpm * 1e-6));

    // Degenerate inputs fail closed.
    BOOST_CHECK(!integer_chip_grid(0.0));
    BOOST_CHECK(!integer_chip_grid(-2.0));
    BOOST_CHECK(!integer_chip_grid_at(2.0, 0));
    BOOST_CHECK(!integer_chip_grid_at(2.0, 2, 0.0));
    BOOST_CHECK(!integer_chip_grid_at(2.0, 2, -1e-6));
    BOOST_CHECK(!integer_chip_grid_at(2.0, 2, std::numeric_limits<double>::quiet_NaN()));

    // A different integer grid is not silently accepted as 2 samples/chip.
    BOOST_CHECK(!integer_chip_grid_at(1.0, 2));
    BOOST_CHECK(!integer_chip_grid_at(4.0, 2));
    BOOST_CHECK(integer_chip_grid(1.0)); // its own nearest integer, though
}


BOOST_AUTO_TEST_CASE(test_demod_core_r1_timing_matches_golden)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    BOOST_REQUIRE_EQUAL(tmpl.size(), size_t(1016));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();

    TimingResult tr;
    // seed = first SYNC start (9984) as the detector would supply.
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    // Absolute coordinate start (golden: 9984).
    BOOST_CHECK_CLOSE(static_cast<double>(tr.preamble_start_sample), 9984.0, 0.05);
    // Period (golden: 1016).
    BOOST_CHECK_CLOSE(tr.measured_period, 1016.0, 0.1);
    // Metric (golden: 1.0).
    BOOST_CHECK_GE(tr.metric, 0.95f);
    // Peaks (golden: 64-SYNC grid, sparse anchor fit <= 16 measured anchors).
    BOOST_CHECK_EQUAL(tr.cir_grid_repetitions, size_t(64));
    BOOST_CHECK_GE(tr.detected_peaks, size_t(8));
    // First peak is symbol-END = start + L - 1 (golden: 10999).
    if (!tr.peak_samples.empty())
        BOOST_CHECK_CLOSE(static_cast<double>(tr.peak_samples.front()), 10999.0, 0.1);
}

BOOST_AUTO_TEST_CASE(test_demod_core_tail_seed_refine_radius_fallback)
{
    // SFD-anchored tail profile: the first tail SYNC prediction is within a
    // few samples, so the experiment tries +/-8 before the legacy wide search.
    // A miss at +12 samples must fall back and still recover all 40 peaks.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    BOOST_REQUIRE_EQUAL(tmpl.size(), size_t(1016));

    auto prof = Qm35825Profile::Default();
    prof.preamble_repetitions = 40;
    prof.timing_search_margin = 0;
    prof.timing_coarse_stride = 1;
    prof.timing_max_backtrack_symbols = 0;
    prof.timing_track_radius = 8;
    prof.radar_verification_threshold = 0.20f;
    prof.timing_seed_refine_radius = 8;
    const int64_t tail_start = 9984 + 24 * static_cast<int64_t>(kQm35SamplesPerSymbol);

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    TimingResult narrow;
    BOOST_REQUIRE(core::stage_timing(
        iq.data(), iq.size(), prof, tmpl, tail_start, narrow, scratch));
    BOOST_CHECK_EQUAL(narrow.cir_grid_repetitions, size_t(40));
    BOOST_CHECK_GE(narrow.detected_peaks, size_t(8));
    BOOST_CHECK_SMALL(
        std::abs(static_cast<double>(narrow.preamble_start_sample) -
                 static_cast<double>(tail_start)),
        2.0);
    BOOST_CHECK_SMALL(std::abs(narrow.measured_period - 1016.0), 0.1);

    TimingResult fallback;
    BOOST_REQUIRE(core::stage_timing(
        iq.data(), iq.size(), prof, tmpl, tail_start + 12, fallback, scratch));
    BOOST_CHECK_EQUAL(fallback.cir_grid_repetitions, size_t(40));
    BOOST_CHECK_GE(fallback.detected_peaks, size_t(8));
    BOOST_CHECK_SMALL(
        std::abs(static_cast<double>(fallback.preamble_start_sample) -
                 static_cast<double>(tail_start)),
        2.0);

    prof.timing_seed_refine_radius = 0;
    TimingResult legacy;
    BOOST_REQUIRE(core::stage_timing(
        iq.data(), iq.size(), prof, tmpl, tail_start, legacy, scratch));
    BOOST_CHECK_EQUAL(legacy.cir_grid_repetitions, size_t(40));
    BOOST_CHECK_GE(legacy.detected_peaks, size_t(8));
}

BOOST_AUTO_TEST_CASE(test_locked_grid_timing_matches_legacy_train)
{
    // LockedBounded grid tracking replaces the per-SYNC +/-8 search with a
    // sparse anchor fit plus one correlation per SYNC.  On a clean train it
    // must reproduce the legacy peak train and the fitted period, because CFO
    // (peak_corr phases) and CIR (measured_period) are downstream consumers.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    auto prof = Qm35825Profile::Default();
    prof.sync_strategy = SyncStrategy::Locked;
    prof.locked_roi_samples = 256;
    prof.locked_inner_symbols = 4;
    prof.locked_allow_backtrack = false;

    core::DemodScratch scratch_a;
    scratch_a.reserve(iq.size());
    prof.sync_strategy = SyncStrategy::Acquire;
    TimingResult legacy;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, legacy, scratch_a));

    core::DemodScratch scratch_b;
    scratch_b.reserve(iq.size());
    prof.sync_strategy = SyncStrategy::Locked;
    TimingResult grid;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, grid, scratch_b));
    BOOST_CHECK_EQUAL(grid.cir_grid_repetitions, prof.preamble_repetitions);
    BOOST_CHECK_EQUAL(legacy.cir_grid_repetitions, prof.preamble_repetitions);
    BOOST_CHECK_GE(grid.detected_peaks, size_t(8));
    BOOST_CHECK_EQUAL(grid.detected_peaks, legacy.detected_peaks);
    BOOST_CHECK_SMALL(std::abs(grid.measured_period - legacy.measured_period), 0.02);
    BOOST_CHECK_SMALL(std::abs(static_cast<double>(grid.preamble_start_sample) -
                               static_cast<double>(legacy.preamble_start_sample)),
                      2.0);
    // Locked contract still holds: no backtrack, no legacy scan.
    BOOST_CHECK_EQUAL(grid.backtrack_symbols, 0u);
    BOOST_CHECK_EQUAL(grid.bounded_legacy_scans, 0u);
    BOOST_CHECK_EQUAL(grid.bounded_success, 1u);

    // End to end: the grid path must still decode the golden packet, including
    // the SFD-anchored tail re-track (which runs stage_timing a second time on
    // the 40-SYNC tail with the same bounded profile).
    auto p_legacy = prof;
    p_legacy.locked_grid_timing = false;
    p_legacy.sfd_mode = "ieee";
    auto p_grid = prof;
    p_grid.locked_grid_timing = true;
    p_grid.sfd_mode = "ieee";
    core::DemodScratch s_legacy, s_grid;
    s_legacy.reserve(iq.size());
    s_grid.reserve(iq.size());
    auto r_legacy =
        core::demodulate_one(iq.data(), iq.size(), p_legacy, 1, 9984, 0, tmpl, s_legacy);
    auto r_grid =
        core::demodulate_one(iq.data(), iq.size(), p_grid, 1, 9984, 0, tmpl, s_grid);
    BOOST_CHECK(r_grid.status == DemodStatus::Success);
    BOOST_CHECK(r_legacy.status == r_grid.status);
    BOOST_CHECK(r_grid.payload.bytes == r_legacy.payload.bytes);
    BOOST_CHECK(r_grid.payload.fcs_pass == r_legacy.payload.fcs_pass);
    BOOST_CHECK_EQUAL(r_grid.timing.grid_fallbacks, 0u);
    BOOST_CHECK_SMALL(
        std::abs(r_grid.timing.measured_period - r_legacy.timing.measured_period), 0.02);
    BOOST_CHECK_SMALL(std::abs(r_grid.cfo.cfo_hz - r_legacy.cfo.cfo_hz), 15.0);
    BOOST_CHECK_SMALL(std::abs(static_cast<double>(r_grid.sfd.sfd_start_sample) -
                               static_cast<double>(r_legacy.sfd.sfd_start_sample)),
                      2.0);
}

BOOST_AUTO_TEST_CASE(test_locked_grid_timing_period_is_measured_not_nominal)
{
    // The grid path must MEASURE the period (anchors bootstrap the grid, then a
    // least-squares fit over every refined peak), never report the nominal 1016
    // grid: measured_period feeds integer_chip_grid()/spc_ppm and the fused-CMF
    // admission, so a degenerate constant would silently disable the
    // Full->Bypass (grid_reject) protection for an untrimmed clock.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    // Compress the sample clock by 50 ppm: out[n] = in[n*rate] maps a signal
    // feature at input index m to output index m/rate, so the SYNC period
    // becomes 1016/rate = 1015.949 samples (drift over 63 symbols is 3.2
    // samples, inside the +/-8 anchor radius) and the integer chip grid must
    // be rejected.
    const double ppm = 50.0;
    const double rate = 1.0 + ppm * 1e-6;
    std::vector<gr_complex> stretched;
    stretched.resize(
        static_cast<size_t>(std::floor(static_cast<double>(iq.size() - 1) / rate)));
    for (size_t n = 0; n < stretched.size(); ++n) {
        const double x = static_cast<double>(n) * rate;
        const size_t i0 = static_cast<size_t>(x);
        const double f = x - static_cast<double>(i0);
        const size_t i1 = std::min(i0 + 1, iq.size() - 1);
        stretched[n] =
            iq[i0] * static_cast<float>(1.0 - f) + iq[i1] * static_cast<float>(f);
    }

    auto prof = Qm35825Profile::Default();
    prof.sync_strategy = SyncStrategy::Locked;
    prof.locked_roi_samples = 256;
    prof.locked_inner_symbols = 4;
    prof.locked_grid_timing = true;

    core::DemodScratch scratch;
    scratch.reserve(stretched.size());
    TimingResult tr;
    BOOST_REQUIRE(core::stage_timing(
        stretched.data(), stretched.size(), prof, tmpl, 9984, tr, scratch));
    BOOST_CHECK_GE(tr.detected_peaks, size_t(8));
    BOOST_CHECK_EQUAL(tr.bounded_legacy_scans, 0u);
    const double expected = 1016.0 / rate;
    BOOST_CHECK_SMALL(std::abs(tr.measured_period - expected), 0.05);
    // Not the nominal grid.
    BOOST_CHECK_GT(std::abs(tr.measured_period - 1016.0), 0.02);
    // 50 ppm is far outside the integer chip-grid drift budget, so the fused
    // path must still be refused (this is the Full->Bypass protection).
    BOOST_CHECK(
        !core::integer_chip_grid(tr.measured_period / double(kQm35ChipsPerSymbol)));
}

BOOST_AUTO_TEST_CASE(test_locked_tail_probe_rescues_leading_interference)
{
    // A burst interferer can cover only the first few SYNCs while the tail is
    // clean (5.2 m capture: first_5_max 0.029..0.053, med_tail 0.14..0.20).
    // The bounded probe must then keep probing sparse tail anchors instead of
    // declaring a tracking failure, and the rescued origin must still be the
    // true first SYNC start.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    constexpr size_t kStart = 9984;
    constexpr size_t kSym = kQm35SamplesPerSymbol;
    constexpr size_t kAtten = 20;

    std::vector<gr_complex> hit = iq;
    // Replace the leading SYNCs with uncorrelated noise of similar power: the
    // normalised metric is scale-invariant, so attenuating them would NOT
    // lower it.  Real burst interference decorrelates the SYNC, which is what
    // the 5.2 m capture shows (first_5_max 0.029..0.053).
    std::mt19937 rng(20261116u);
    std::normal_distribution<float> nd(0.0f, 0.35f);
    float rms = 0.0f;
    for (size_t i = 0; i < kSym; ++i)
        rms += std::norm(iq[kStart + i]);
    rms = std::sqrt(rms / static_cast<float>(kSym));
    for (size_t k = 0; k < kAtten; ++k) {
        const size_t base = kStart + k * kSym;
        BOOST_REQUIRE_LT(base + kSym, hit.size());
        for (size_t i = 0; i < kSym; ++i)
            hit[base + i] = gr_complex(nd(rng), nd(rng)) * rms;
    }

    auto prof = Qm35825Profile::Default();
    prof.sync_strategy = SyncStrategy::Locked;
    prof.locked_roi_samples = 256;
    prof.locked_inner_symbols = 4;
    prof.locked_grid_timing = true;

    // Unified Stage 1 CIR accumulation rescues leading-edge interference without
    // relying on a separate tail_probe fallback toggle.
    core::DemodScratch scratch_b;
    scratch_b.reserve(hit.size());
    TimingResult res;
    BOOST_REQUIRE(
        core::stage_timing(hit.data(), hit.size(), prof, tmpl, kStart, res, scratch_b));
    BOOST_CHECK_EQUAL(res.bounded_failures, 0u);
    BOOST_CHECK_EQUAL(res.cir_grid_repetitions, prof.preamble_repetitions);
    BOOST_CHECK_GE(res.detected_peaks, size_t(8));
    // The anchor fit extrapolates back to the true first SYNC start.
    BOOST_CHECK_SMALL(std::abs(static_cast<double>(res.preamble_start_sample) -
                               static_cast<double>(kStart)),
                      4.0);
    BOOST_CHECK_SMALL(std::abs(res.measured_period - 1016.0), 0.05);
    BOOST_CHECK_EQUAL(res.backtrack_symbols, 0u);
    BOOST_CHECK_EQUAL(res.bounded_legacy_scans, 0u);
}

BOOST_AUTO_TEST_CASE(test_demod_core_rxopt_on_demand_rejects_fractional_grid)
{
    // The on-demand experiment must not round a fractional SFO chip grid or
    // change RAKE/bypass semantics; it declines and leaves the full legacy
    // path available.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    auto prof = Qm35825Profile::Default();
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    TimingResult fractional = tr;
    fractional.measured_period += 0.5;
    const int8_t* pc = GetPreambleCode(prof.code_index);
    std::vector<int8_t> pcode(pc, pc + kQm35CodeLength);
    prof.rx_softchip_on_demand = true;
    CirResult prefix;
    BOOST_CHECK(!core::stage_cir_softchips(iq.data(),
                                           iq.size(),
                                           prof,
                                           fractional,
                                           pcode,
                                           prefix,
                                           scratch,
                                           rxopt::initial_chip_budget(prof),
                                           0));
    BOOST_CHECK(!scratch.softchip.prepared);

    prof.rx_softchip_on_demand = false;
    CirResult full;
    BOOST_REQUIRE(core::stage_cir_softchips(
        iq.data(), iq.size(), prof, fractional, pcode, full, scratch));
    BOOST_CHECK_GT(full.soft_chip_count, 0u);
    BOOST_CHECK_GT(full.samples_per_chip, 0.0);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_backtracks_weak_qm35_first_sync)
{
    std::vector<gr_complex> tmpl;
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    BOOST_REQUIRE_EQUAL(tmpl.size(), size_t(1016));

    auto prof = Qm35825Profile::Default();
    constexpr int64_t start = 2048;
    const size_t sym = tmpl.size();
    const auto sfd = GetSfdSequence("4z2");
    std::vector<gr_complex> iq(static_cast<size_t>(start) +
                                   (prof.preamble_repetitions + sfd.size()) * sym,
                               gr_complex(0.0f, 0.0f));

    // The real QM35 startup transient can leave the first few SYNCs below the
    // normal 0.3 threshold.  Make SYNC #1/#2 metric 0.25 and seed at strong
    // SYNC #3; timing must lock strong first, then recover both predecessors.
    std::vector<gr_complex> orth(sym);
    for (size_t i = 0; i < sym; ++i)
        orth[i] = gr_complex(std::sin(0.173f * static_cast<float>(i)),
                             std::cos(0.119f * static_cast<float>(i)));
    gr_complex projection(0.0f, 0.0f);
    float et = 0.0f;
    for (size_t i = 0; i < sym; ++i) {
        projection += std::conj(tmpl[i]) * orth[i];
        et += std::norm(tmpl[i]);
    }
    for (size_t i = 0; i < sym; ++i)
        orth[i] -= (projection / et) * tmpl[i];
    float eo = 0.0f;
    for (const auto& v : orth)
        eo += std::norm(v);
    const float orth_scale = std::sqrt(et / eo);
    for (size_t rep = 0; rep < 2; ++rep)
        for (size_t i = 0; i < sym; ++i)
            iq[static_cast<size_t>(start) + rep * sym + i] =
                0.5f * tmpl[i] + std::sqrt(0.75f) * orth_scale * orth[i];
    for (size_t rep = 2; rep < prof.preamble_repetitions; ++rep)
        std::copy(tmpl.begin(),
                  tmpl.end(),
                  iq.begin() + start + static_cast<int64_t>(rep * sym));
    for (size_t k = 0; k < sfd.size(); ++k)
        for (size_t i = 0; i < sym; ++i)
            iq[static_cast<size_t>(start) + (prof.preamble_repetitions + k) * sym + i] =
                static_cast<float>(sfd[k]) * tmpl[i];

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    TimingResult tr;
    BOOST_REQUIRE(core::stage_timing(iq.data(),
                                     iq.size(),
                                     prof,
                                     tmpl,
                                     start + static_cast<int64_t>(2 * sym),
                                     tr,
                                     scratch));
    BOOST_CHECK_GE(tr.preamble_start_sample, start);
    BOOST_CHECK_LE(tr.preamble_start_sample,
                   start + static_cast<int64_t>(2 * sym));
    BOOST_CHECK_EQUAL(tr.cir_grid_repetitions, prof.preamble_repetitions);
    BOOST_CHECK_GE(tr.detected_peaks, size_t(8));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_cfo_matches_golden)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();

    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    CfoResult cf;
    BOOST_REQUIRE(core::stage_cfo(iq.data(), iq.size(), prof, tr, cf, scratch));
    // Golden: 0 Hz.  Allow tight tolerance.
    BOOST_CHECK_LT(std::abs(cf.cfo_hz), 50.0);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_cfo_skips_24_startup_syncs)
{
    // Mirror MATLAB compensateCarrierOffset: the first 24 matched-filter
    // phases can be a front-end transient; retain the later 40 clean peaks.
    auto prof = Qm35825Profile::Default();
    constexpr size_t kPeaks = 64;
    constexpr size_t kTransient = 24;
    constexpr double kCfoHz = 1900.0;

    TimingResult tr;
    tr.peak_samples.reserve(kPeaks);
    tr.peak_corr.reserve(kPeaks);
    for (size_t k = 0; k < kPeaks; ++k) {
        const int64_t sample = 4096 + static_cast<int64_t>(k * kQm35SamplesPerSymbol);
        tr.peak_samples.push_back(sample);
        // Deliberately flat/invalid startup phase.  The clean segment starts
        // at zero phase so phase unwrap across the discarded boundary is safe.
        const double phase =
            k < kTransient
                ? 0.0
                : 2.0 * M_PI * kCfoHz *
                      static_cast<double>(sample - tr.peak_samples[kTransient]) /
                      prof.sample_rate;
        tr.peak_corr.emplace_back(static_cast<float>(std::cos(phase)),
                                  static_cast<float>(std::sin(phase)));
    }

    std::vector<gr_complex> iq(static_cast<size_t>(tr.peak_samples.back()) + 1,
                               gr_complex(0.0f, 0.0f));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    CfoResult cf;
    BOOST_REQUIRE(core::stage_cfo(iq.data(), iq.size(), prof, tr, cf, scratch));
    BOOST_CHECK_EQUAL(cf.peaks_used, kPeaks - kTransient);
    BOOST_CHECK_EQUAL(cf.skipped_peaks, kTransient);
    BOOST_CHECK_EQUAL(cf.fit_first_peak_sample, tr.peak_samples[kTransient]);
    BOOST_CHECK_EQUAL(cf.fit_last_peak_sample, tr.peak_samples.back());
    BOOST_CHECK_CLOSE(cf.cfo_hz, kCfoHz, 0.1);
}

BOOST_AUTO_TEST_CASE(test_demod_core_low_snr_sparse_anchor_cfo)
{
    // The low-SNR fold path reports 16 super-anchors spaced 16 SYNC periods
    // (16 * 1016 samples).  The legacy unwrap gate only accumulated pairs with
    // dt <= 1.5 symbols, so it collected ZERO pairs and the slope silently
    // stayed 0; the LS fit then ran on wrapped phases and happened to work
    // only because |CFO| * 256 symbols stays inside +-2*pi for small offsets.
    // The sparse gate (dt <= 20 symbols) must recover the true slope.
    auto prof = Qm35825Profile::Default();
    constexpr size_t kAnchors = 16;
    constexpr size_t kGroup = 16;
    const int64_t start = 4096;

    for (double cfo : { 1900.0, 4000.0, -2500.0 }) {
        TimingResult tr;
        for (size_t g = 0; g < kAnchors; ++g) {
            const int64_t sample =
                start + static_cast<int64_t>(g * kGroup) * kQm35SamplesPerSymbol;
            tr.peak_samples.push_back(sample);
            const double phase =
                2.0 * M_PI * cfo * static_cast<double>(sample - start) / prof.sample_rate;
            tr.peak_corr.emplace_back(static_cast<float>(std::cos(phase)),
                                      static_cast<float>(std::sin(phase)));
            tr.peak_metrics.push_back(0.56f);
        }
        tr.detected_peaks = prof.preamble_repetitions;
        tr.measured_period = static_cast<double>(kQm35SamplesPerSymbol);

        std::vector<gr_complex> iq(static_cast<size_t>(tr.peak_samples.back()) + 1,
                                   gr_complex(0.0f, 0.0f));
        core::DemodScratch scratch;
        scratch.reserve(iq.size());
        CfoResult cf;
        BOOST_REQUIRE(core::fit_cfo(iq.data(), iq.size(), prof, tr, cf, scratch));
        BOOST_CHECK_EQUAL(cf.peaks_used, kAnchors);
        BOOST_CHECK_CLOSE(cf.cfo_hz, cfo, 0.1);
    }

    // With one missing super-anchor (32-symbol gap) at 18 kHz CFO, the gap
    // rotates by > pi (~3.68 rad), so zero-slope unwrap fails while the
    // dt <= 20*1016 slope estimator recovers the exact 18 kHz frequency.
    {
        constexpr double kGapCfo = 18000.0;
        TimingResult tr;
        for (size_t g = 0; g < kAnchors; ++g) {
            if (g == 5)
                continue;
            const int64_t sample =
                start + static_cast<int64_t>(g * kGroup) * kQm35SamplesPerSymbol;
            tr.peak_samples.push_back(sample);
            const double phase = 2.0 * M_PI * kGapCfo *
                                 static_cast<double>(sample - start) / prof.sample_rate;
            tr.peak_corr.emplace_back(static_cast<float>(std::cos(phase)),
                                      static_cast<float>(std::sin(phase)));
            tr.peak_metrics.push_back(0.56f);
        }
        tr.detected_peaks = prof.preamble_repetitions;
        tr.measured_period = static_cast<double>(kQm35SamplesPerSymbol);
        std::vector<gr_complex> iq(static_cast<size_t>(tr.peak_samples.back()) + 1,
                                   gr_complex(0.0f, 0.0f));
        core::DemodScratch scratch;
        scratch.reserve(iq.size());
        CfoResult cf;
        BOOST_REQUIRE(core::fit_cfo(iq.data(), iq.size(), prof, tr, cf, scratch));
        BOOST_CHECK_EQUAL(cf.peaks_used, kAnchors - 1);
        BOOST_CHECK_CLOSE(cf.cfo_hz, kGapCfo, 0.1);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_matches_golden)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();

    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));

    // Golden uses IEEE legacy SFD (the lrwpan generator's default SFD).
    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    // SFD begins right after the 64 SYNC symbols: start + 64*1016 = 75008.
    BOOST_CHECK_CLOSE(static_cast<double>(sr.expected_start_sample), 75008.0, 0.1);
    BOOST_CHECK_CLOSE(static_cast<double>(sr.sfd_start_sample), 75008.0, 0.1);
    BOOST_CHECK_GE(sr.metric, 0.95f);
    BOOST_CHECK_EQUAL(sr.search_windows, 1u);
    BOOST_CHECK_GT(sr.coarse_correlations, 0u);
    BOOST_CHECK_GT(sr.fine_correlations, 0u);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_symmetric_one_symbol_forward)
{
    // MATLAB-style symmetric +-1-symbol search: an SFD one symbol LATER than
    // the nominal position (preamble start one SYNC early) is recovered inside
    // the expected +- samples_per_symbol window.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();

    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    for (auto& peak : tr.peak_samples)
        peak -= static_cast<int64_t>(kQm35SamplesPerSymbol);
    tr.preamble_start_sample -= static_cast<int64_t>(kQm35SamplesPerSymbol);

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
    // One symbol later: ±64 at nominal misses, then +1 symbol hits (2 windows).
    BOOST_CHECK_EQUAL(sr.search_windows, 2u);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_symmetric_partial_train_backward)
{
    // An SFD one symbol EARLIER than the nominal position (preamble start one
    // SYNC late) is recovered by the same symmetric window, even from a
    // 63/64 tracked train (np no longer narrows the half-width).
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    for (auto& peak : tr.peak_samples)
        peak += static_cast<int64_t>(kQm35SamplesPerSymbol);
    tr.preamble_start_sample += static_cast<int64_t>(kQm35SamplesPerSymbol);
    BOOST_REQUIRE_GE(tr.peak_samples.size(), size_t(8));
    BOOST_REQUIRE_EQUAL(tr.cir_grid_repetitions, size_t(64));

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
    // One symbol earlier: nominal miss, +1 miss, then -1 hits (3 windows).
    BOOST_CHECK_EQUAL(sr.search_windows, 3u);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_recovers_four_symbol_forward)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    for (auto& peak : tr.peak_samples)
        peak -= static_cast<int64_t>(2 * kQm35SamplesPerSymbol);
    tr.preamble_start_sample -= static_cast<int64_t>(2 * kQm35SamplesPerSymbol);

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_recovers_four_symbol_backward)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    for (auto& peak : tr.peak_samples)
        peak += static_cast<int64_t>(2 * kQm35SamplesPerSymbol);
    tr.preamble_start_sample += static_cast<int64_t>(2 * kQm35SamplesPerSymbol);

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_grid_rescues_five_symbol_forward)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    for (auto& peak : tr.peak_samples)
        peak -= static_cast<int64_t>(5 * kQm35SamplesPerSymbol);
    tr.preamble_start_sample -= static_cast<int64_t>(5 * kQm35SamplesPerSymbol);

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    // The +/-4-symbol direct windows miss, but the +/-36-symbol discrete grid
    // is the documented low-SNR rescue and must find the true SFD.  (The old
    // `par_only &&` guard in polarity_ok made every ieee/decawave SFD fail
    // that check, which is why this test used to expect a miss.)
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_sfd_grid_rescues_five_symbol_backward)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    for (auto& peak : tr.peak_samples)
        peak += static_cast<int64_t>(5 * kQm35SamplesPerSymbol);
    tr.preamble_start_sample += static_cast<int64_t>(5 * kQm35SamplesPerSymbol);

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_sfd_sfo_grid_cache_equivalence)
{
    // The per-delta symbol cache must reproduce the direct per-symbol grid on
    // an untrimmed clock (+27 ppm), where round(k*period) - k*sym is not zero
    // and the old code disabled the decomposition entirely.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    tr.measured_period += 0.027; // +27 ppm on the 1016-sample SYNC period
    const auto sfd_seq = GetSfdSequence("4z2");

    Qm35825Profile dense = prof;
    dense.rx_fast_sfd = false;
    dense.rx_timing_sparse = false;
    SfdResult direct;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), dense, tr, sfd_seq, tmpl, direct, scratch));

    Qm35825Profile cached = prof;
    cached.rx_fast_sfd = true;
    cached.rx_timing_sparse = false;
    SfdResult fast;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), cached, tr, sfd_seq, tmpl, fast, scratch));

    BOOST_CHECK_EQUAL(fast.sfd_start_sample, direct.sfd_start_sample);
    BOOST_CHECK_SMALL(std::abs(fast.metric - direct.metric), 1e-4f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r1_bad_input_fails_cleanly)
{
    std::vector<gr_complex> noise(4000); // too short for a full SYNC
    for (size_t i = 0; i < noise.size(); ++i)
        noise[i] = gr_complex((float)(i % 17) / 17.0f, 0.0f);
    std::vector<gr_complex> tmpl(1016, gr_complex(0.0f, 0.0f));
    core::DemodScratch scratch;
    scratch.reserve(4096);
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_CHECK(
        !core::stage_timing(noise.data(), noise.size(), prof, tmpl, 0, tr, scratch));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r2_cir_matches_golden)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iq, scratch, prof, tr, cf, cir));

    // Golden CIR: 38 normalized taps, pre=8, post=30, first path = 9986.
    std::vector<float> gcir;
    BOOST_REQUIRE(load_f32(golden_dir() + "/stage_cir.f32", gcir));
    BOOST_REQUIRE_EQUAL(gcir.size(), size_t(38));
    BOOST_REQUIRE_EQUAL(cir.cir_values.size(), gcir.size());
    BOOST_CHECK_EQUAL(cir.pre_samples, size_t(8));
    BOOST_CHECK_EQUAL(cir.post_samples, size_t(30));
    BOOST_CHECK_EQUAL(cir.first_path_sample, size_t(9986));

    size_t pk = 0;
    float pm = -1.0f;
    for (size_t i = 0; i < cir.cir_values.size(); ++i) {
        const float a = std::abs(cir.cir_values[i]);
        if (a > pm) {
            pm = a;
            pk = i;
        }
    }
    BOOST_CHECK_EQUAL(pk, size_t(10)); // offset +2
    BOOST_CHECK_CLOSE(cir.cir_peak_metric, 0.732856f, 1.0);

    float l2 = 0.0f;
    for (size_t i = 0; i < gcir.size(); ++i) {
        const float d = cir.cir_values[i] - gcir[i];
        l2 += d * d;
    }
    BOOST_CHECK_SMALL(std::sqrt(l2), 0.01f); // actual ~6e-8
}

BOOST_AUTO_TEST_CASE(test_demod_core_cir_fused_rotation_matches_derotated)
{
    // CIR-1: folding the CFO rotation into the SYNC average (raw rx + per-SYNC
    // scalar + one per-sample ramp) must reproduce the legacy "derotate the
    // whole preamble, then average" result, and the tail average must give the
    // same sync_phase_gain statistic without re-reading 32 SYNC windows.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const double cfo_hz = 25000.0;
    const double w = 2.0 * M_PI * cfo_hz / kQm35SampleRate;
    std::vector<gr_complex> iqc(iq.size());
    for (size_t i = 0; i < iq.size(); ++i) {
        const double ph = w * static_cast<double>(i);
        iqc[i] = iq[i] * gr_complex(static_cast<float>(std::cos(ph)),
                                    static_cast<float>(std::sin(ph)));
    }

    core::DemodScratch scratch_a;
    core::DemodScratch scratch_b;
    scratch_a.reserve(iqc.size());
    scratch_b.reserve(iqc.size());
    auto prof = Qm35825Profile::Default();
    std::vector<gr_complex> tmpl;
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    TimingResult tr;
    CfoResult cf;
    BOOST_REQUIRE(
        core::stage_timing(iqc.data(), iqc.size(), prof, tmpl, 9984, tr, scratch_a));
    BOOST_REQUIRE(core::stage_cfo(iqc.data(), iqc.size(), prof, tr, cf, scratch_a));
    BOOST_CHECK_CLOSE(cf.cfo_hz, cfo_hz, 10.0);
    std::vector<int8_t> code9(kPreambleCode9.begin(), kPreambleCode9.end());

    CirResult legacy;
    CirResult fused;
    BOOST_REQUIRE(core::estimate_cir(
        scratch_a.derotated.data(), iqc.size(), prof, tr, code9, legacy, scratch_a));
    BOOST_REQUIRE(core::estimate_cir(
        iqc.data(), iqc.size(), prof, tr, code9, fused, scratch_b, cf.cfo_hz));
    BOOST_REQUIRE_EQUAL(fused.cir_values.size(), legacy.cir_values.size());
    BOOST_CHECK_EQUAL(fused.first_path_sample, legacy.first_path_sample);
    BOOST_CHECK_SMALL(std::abs(fused.cir_peak_metric - legacy.cir_peak_metric), 1e-4f);
    double l2 = 0.0;
    for (size_t i = 0; i < legacy.cir_values.size(); ++i) {
        const double d = static_cast<double>(fused.cir_values[i] - legacy.cir_values[i]);
        l2 += d * d;
    }
    BOOST_CHECK_SMALL(std::sqrt(l2), 1e-4);

    BOOST_REQUIRE(scratch_b.cir_avg_tail_valid);
    const auto spread = BuildSampledCode(code9.data(), code9.size());
    const size_t taps = fused.cir_filter_taps.size();
    const int64_t chip_start = tr.preamble_start_sample + static_cast<int64_t>(taps) -
                               static_cast<int64_t>(prof.cir_pre_samples) - 1;
    const size_t phase_reps = std::min<size_t>(32, prof.preamble_repetitions);
    const size_t phase_first =
        (prof.preamble_repetitions - phase_reps) * kQm35ChipsPerSymbol;
    const int64_t base = chip_start - static_cast<int64_t>(taps - 1) +
                         static_cast<int64_t>(2 * phase_first);
    std::vector<gr_complex> sum;
    gr_complex gain_legacy{};
    gr_complex gain_tail{};
    BOOST_REQUIRE(rxopt::sync_phase_gain(scratch_a.derotated.data(),
                                         iqc.size(),
                                         base,
                                         fused.cir_filter_taps.data(),
                                         taps,
                                         spread,
                                         phase_reps,
                                         sum,
                                         gain_legacy));
    BOOST_REQUIRE(rxopt::sync_phase_gain_from_tail(scratch_b.cir_avg_tail.data(),
                                                   fused.cir_filter_taps.data(),
                                                   taps,
                                                   spread,
                                                   gain_tail));
    BOOST_CHECK_SMALL(std::abs(gain_tail - gain_legacy) /
                          std::max(1e-6f, std::abs(gain_legacy)),
                      1e-3f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r2_softchips_match_golden)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iq, scratch, prof, tr, cf, cir));

    // Golden soft-chip stream: 154578 chips, samples_per_chip = 2, first chip
    // at absolute sample 10013 (= start + post - 1).
    std::vector<float> gsoft;
    BOOST_REQUIRE(load_f32(golden_dir() + "/stage_softchips.f32", gsoft));
    BOOST_REQUIRE_EQUAL(gsoft.size(), size_t(154578));
    BOOST_REQUIRE_EQUAL(scratch.soft_chips.size(), gsoft.size());
    BOOST_CHECK_EQUAL(cir.soft_chip_count, size_t(154578));
    BOOST_CHECK_CLOSE(cir.samples_per_chip, 2.0, 0.5);
    BOOST_CHECK_EQUAL(tr.preamble_start_sample + (int64_t)cir.post_samples - 1,
                      int64_t(10013));

    float mx = 0.0f;
    for (size_t i = 0; i < gsoft.size(); ++i)
        mx = std::max(mx, std::abs(scratch.soft_chips[i] - gsoft[i]));
    BOOST_CHECK_LT(mx, 1e-3f); // actual ~6e-7
    // First soft chip is near the normalized max (preamble present).
    BOOST_CHECK_GT(std::abs(scratch.soft_chips[0]), 0.9f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r2_cfo_compensated_matches_golden)
{
    // Inject +1 kHz CFO into the golden window; after stage_cfo derotation the
    // CIR and soft chips must still match the (CFO-free) golden reference.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    auto prof = Qm35825Profile::Default();
    const double w = 2.0 * M_PI * 1000.0 / prof.sample_rate;
    std::vector<gr_complex> iqc(iq.size());
    for (size_t i = 0; i < iq.size(); ++i) {
        const double ph = w * static_cast<double>(i);
        iqc[i] = iq[i] * gr_complex(static_cast<float>(std::cos(ph)),
                                    static_cast<float>(std::sin(ph)));
    }

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iqc, scratch, prof, tr, cf, cir));

    // Estimated CFO should be near the injected +1 kHz.
    BOOST_CHECK_CLOSE(cf.cfo_hz, 1000.0, 10.0);

    std::vector<float> gcir;
    BOOST_REQUIRE(load_f32(golden_dir() + "/stage_cir.f32", gcir));
    float l2 = 0.0f;
    for (size_t i = 0; i < gcir.size(); ++i) {
        const float d = cir.cir_values[i] - gcir[i];
        l2 += d * d;
    }
    BOOST_CHECK_SMALL(std::sqrt(l2), 0.05f);

    std::vector<float> gsoft;
    BOOST_REQUIRE(load_f32(golden_dir() + "/stage_softchips.f32", gsoft));
    float mx = 0.0f;
    for (size_t i = 0; i < gsoft.size(); ++i)
        mx = std::max(mx, std::abs(scratch.soft_chips[i] - gsoft[i]));
    BOOST_CHECK_LT(mx, 0.05f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r2_bad_input_fails_cleanly)
{
    std::vector<gr_complex> noise(4000);
    for (size_t i = 0; i < noise.size(); ++i)
        noise[i] = gr_complex(static_cast<float>(i % 17) / 17.0f, 0.0f);
    core::DemodScratch scratch;
    scratch.reserve(4096);
    auto prof = Qm35825Profile::Default();
    // A too-short / silent buffer must fail cleanly (no CIR / no soft chips).
    std::vector<int8_t> code9(kPreambleCode9.begin(), kPreambleCode9.end());
    TimingResult tr;
    tr.preamble_start_sample = 0;
    tr.measured_period = 1016.0;
    tr.detected_peaks = 1; // too few repetitions
    CirResult cir;
    BOOST_CHECK(!core::stage_cir_softchips(
        noise.data(), noise.size(), prof, tr, code9, cir, scratch));
}

// ---------------------------------------------------------------------------
// R3: NS-SFD location + BPRF PHR (locateNsSfd / helperUWBBPRFDemod /
// helperUWBConvDec / helperUWBPHRDecode).
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_demod_core_r3_ns_sfd_matches_golden)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iq, scratch, prof, tr, cf, cir));

    // Golden window uses the IEEE legacy SFD.
    const auto sfd_seq = GetSfdSequence("ieee");
    NsSfdResult ns;
    BOOST_REQUIRE(core::stage_ns_sfd(
        scratch.soft_chips, prof, sfd_seq, kQm35ChipsPerSymbol, ns, scratch));
    // Golden: start_chip=32512, end_chip=36575, polarity=1, correlation=1.0.
    BOOST_CHECK_EQUAL(ns.sfd_start_chip, int64_t(32512));
    BOOST_CHECK_EQUAL(ns.sfd_end_chip, int64_t(36575));
    BOOST_CHECK_EQUAL(ns.polarity, 1);
    BOOST_CHECK_GE(ns.metric, 0.95f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r3_phr_matches_golden)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iq, scratch, prof, tr, cf, cir));

    const auto sfd_seq = GetSfdSequence("ieee");
    NsSfdResult ns;
    BOOST_REQUIRE(core::stage_ns_sfd(
        scratch.soft_chips, prof, sfd_seq, kQm35ChipsPerSymbol, ns, scratch));
    PhrResult phr;
    BOOST_REQUIRE(
        core::stage_phr(scratch.soft_chips, prof, ns, kQm35ChipsPerSymbol, phr, scratch));
    // Golden: psdu_length=127, secded_pass=1, data rate 6.81.
    BOOST_CHECK_EQUAL(phr.psdu_length, uint32_t(127));
    BOOST_CHECK(!phr.secded_uncorrectable);
    BOOST_CHECK_CLOSE(phr.data_rate_mbps, 6.81f, 0.5);
    // Decoded 19 PHR bits match golden (stage_phr_decoded.bin).
    std::vector<int8_t> gdec;
    BOOST_REQUIRE(load_i8(golden_dir() + "/stage_phr_decoded.bin", gdec));
    BOOST_REQUIRE_EQUAL(gdec.size(), size_t(19));
    BOOST_REQUIRE_EQUAL(phr.phr_bits.size(), gdec.size());
    for (size_t i = 0; i < 19; ++i)
        BOOST_CHECK_EQUAL(phr.phr_bits[i], static_cast<uint8_t>(gdec[i]));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r3_full_pipeline_reaches_phr)
{
    // Run the whole demodulate_one chain on the golden window with the SFD
    // switched to the IEEE legacy SFD the window was generated with.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    auto res =
        core::demodulate_one(iq.data(), iq.size(), prof, 1, 9984, 0, tmpl, scratch);
    BOOST_CHECK(res.status == DemodStatus::Success);
    // The complete path (including SFD-anchored CFO retracking) must keep
    // MATLAB's 24-SYNC startup exclusion and 40-point LS fit.
    BOOST_CHECK_EQUAL(res.cfo.skipped_peaks, size_t(24));
    BOOST_CHECK_EQUAL(res.cfo.peaks_used, size_t(40));
    BOOST_CHECK_EQUAL(res.cfo.fit_first_peak_sample, res.timing.peak_samples[24]);
    BOOST_CHECK_EQUAL(res.cfo.fit_last_peak_sample, res.timing.peak_samples.back());
    BOOST_CHECK_EQUAL(res.ns_sfd.sfd_start_chip, int64_t(32512));
    BOOST_CHECK_EQUAL(res.phr.psdu_length, uint32_t(127));
    BOOST_CHECK(!res.phr.secded_uncorrectable);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r3_sfd_child_subtimers_are_consistent)
{
    // S0-A: the SFD parent stage publishes diagnostic child sub-timers.  They
    // must partition `stage_sfd_us` and must NOT be double-counted into
    // `stage_total_us`.  Decoding behaviour itself must be unchanged.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    const auto res = core::demodulate_one(
        iq.data(), iq.size(), prof, 1, 9984, 0, load_reference_template(), scratch);

    // (c) the decode verdict is unchanged (same fixture as the r3 pipeline
    // test above).
    BOOST_REQUIRE(res.status == DemodStatus::Success);
    BOOST_REQUIRE(res.payload.ok);
    BOOST_CHECK(res.payload.fcs_pass);

    // (a) children partition the parent.  `stage_sfd_anchor_other_us` is the
    // clamped remainder, so the sum telescopes back to `stage_sfd_us`; an
    // independent-clock overshoot would clamp to 0 and break the equality.
    const uint64_t children_us = res.stage_sfd_bootstrap_us +
                                 res.stage_sfd_tail_retrack_us +
                                 res.stage_sfd_cfo_refit_us + res.stage_sfd_final_us;
    BOOST_TEST_MESSAGE("stage_sfd_us="
                       << res.stage_sfd_us << " bootstrap=" << res.stage_sfd_bootstrap_us
                       << " tail_retrack=" << res.stage_sfd_tail_retrack_us
                       << " cfo_refit=" << res.stage_sfd_cfo_refit_us
                       << " final=" << res.stage_sfd_final_us
                       << " anchor_other=" << res.stage_sfd_anchor_other_us);
    BOOST_CHECK_LE(children_us, res.stage_sfd_us);
    BOOST_CHECK_EQUAL(children_us + res.stage_sfd_anchor_other_us, res.stage_sfd_us);
    BOOST_CHECK_GT(res.stage_sfd_bootstrap_us, uint64_t(0));
    // The final SFD is skipped when the refit reused the Stage-2 rotation or
    // moved the 8-symbol phase by < 0.05 rad; otherwise it must have run.
    const double drift8_rad =
        std::abs(2.0 * M_PI * (res.cfo.cfo_hz - res.cfo_initial_hz) *
                 (8.0 * static_cast<double>(kQm35SamplesPerSymbol)) / prof.sample_rate);
    BOOST_CHECK(res.stage_sfd_final_us > 0 || res.cfo_refit_action != 0 ||
                drift8_rad < 0.05);

    // (b) `stage_total_us` is the sum of exactly the fields written by the
    // `lap()` lambda (timing, cfo, sfd, cir, ns_sfd, phr, payload) plus the
    // separately timed soft-chip fallback bucket.  The SFD children are NOT
    // separately accumulated, so this equality would fail if they leaked into
    // the total.
    const uint64_t parent_sum = res.stage_timing_us + res.stage_cfo_us +
                                res.stage_sfd_us + res.stage_cir_us +
                                res.stage_ns_sfd_us + res.stage_phr_us +
                                res.stage_payload_us + res.stage_softchip_fallback_us;
    BOOST_CHECK_EQUAL(res.stage_total_us, parent_sum);
}

BOOST_AUTO_TEST_CASE(test_demod_core_sfd_late_acquisition_backward_rescue)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));

    // Shift the detected preamble start 45 symbols late (beyond the +/-36
    // discrete grid) as happens when burst interference buries the leading
    // SYNCs during acquire.
    tr.preamble_start_sample += static_cast<int64_t>(45 * kQm35SamplesPerSymbol);
    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    BOOST_REQUIRE(
        core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);

    // LockedBounded keeps the +/-36 bound and does not scan further backward.
    auto bounded = prof;
    bounded.sync_strategy = SyncStrategy::Locked;
    SfdResult sr_bounded;
    BOOST_CHECK(!core::stage_sfd(
        iq.data(), iq.size(), bounded, tr, sfd_seq, tmpl, sr_bounded, scratch));
}


// ---------------------------------------------------------------------------
// SFD-anchored CFO refit policy (single-pass rotation).
