// SPDX-License-Identifier: GPL-3.0-or-later
//
// Stage-3 SFD refactor QA (2026-11-24):
//   * S1 coherent-weighted soft polarity rescues a single flipped SFD symbol;
//   * S1 rejects shifted-SFD sidelobes and pure-noise false alarms;
//   * S2 the 0.008 low floor hits the centre window without the offset scan;
//   * S4 incremental derotation reproduces the full-range rotation to 1e-6;
//   * late-acquisition (>36 symbol) backward rescue still works through the
//     new range-aware entry point;
//   * S3 the discrete grid is zero-allocation once warm.
//
// Fixtures are the 1016-sample reference SYNC tiled on the 998.4 MS/s WORK
// grid plus white AWGN scaled to an in-band SNR over the 500 MHz UWB channel
// inside the 998.4 MHz work band, the same convention as
// qa_uwb_demod_core_heavy_grid.cc.
#include "demod_core_test_helpers.h"

#include <openuwb/phy/hrp/rx/demod_cfo.h>
#include <openuwb/phy/hrp/rx/demod_sfd.h>
#include <openuwb/phy/hrp/rx/rx_cfo_opt.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {

using namespace openuwb::demod;

constexpr size_t kSym = kQm35SamplesPerSymbol;
constexpr double kWorkBandHz = 998.4e6;
constexpr double kInBandHz = 500.0e6;

struct SfdFixture {
    std::vector<gr_complex> sync;
    std::vector<int8_t> sfd_seq;
    std::vector<gr_complex> sfd_wf;
    size_t sym = kSym;
    size_t reps = 64;
    int64_t pre = 4096;
    int64_t true_start = 4096;
    int64_t true_sfd = 0;
    size_t n = 0;
    double p_sig = 0.0;

    explicit SfdFixture(const char* mode, size_t reps_ = 64)
        : reps(reps_)
    {
        BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", sync));
        BOOST_REQUIRE_EQUAL(sync.size(), sym);
        sfd_seq = GetSfdSequence(mode);
        BOOST_REQUIRE_GT(sfd_seq.size(), 0u);
        sfd_wf.assign(sfd_seq.size() * sym, gr_complex(0.0f, 0.0f));
        for (size_t i = 0; i < sfd_seq.size(); ++i)
            for (size_t k = 0; k < sym; ++k)
                sfd_wf[i * sym + k] =
                    static_cast<float>(sfd_seq[i]) * sync[k];
        true_sfd = true_start + static_cast<int64_t>(reps * sym);
        n = static_cast<size_t>(true_sfd) + sfd_wf.size() + 4096;
        for (size_t i = 0; i < sym; ++i)
            p_sig += std::norm(sync[i]);
        p_sig /= static_cast<double>(sym);
    }

    // Build a frame with the preamble + SFD.  `flip_symbol >= 0` inverts that
    // (non-zero) SFD symbol, i.e. one anti-phase symbol of the coherent train.
    std::vector<gr_complex> make_frame(double snr_db, uint32_t seed,
                                       int flip_symbol = -1) const
    {
        std::vector<gr_complex> rx(n, gr_complex(0.0f, 0.0f));
        for (size_t i = 0; i < reps; ++i)
            std::copy(sync.begin(), sync.end(),
                      rx.begin() + true_start + static_cast<int64_t>(i * sym));
        std::copy(sfd_wf.begin(), sfd_wf.end(), rx.begin() + true_sfd);
        if (flip_symbol >= 0) {
            size_t nz = 0;
            size_t target = 0;
            for (size_t i = 0; i < sfd_seq.size(); ++i) {
                if (sfd_seq[i] == 0)
                    continue;
                if (static_cast<int>(nz) == flip_symbol) {
                    target = i;
                    break;
                }
                ++nz;
            }
            for (size_t k = 0; k < sym; ++k)
                rx[true_sfd + static_cast<int64_t>(target * sym + k)] =
                    -rx[true_sfd + static_cast<int64_t>(target * sym + k)];
        }
        if (std::isfinite(snr_db)) {
            const double snr_lin = std::pow(10.0, snr_db / 10.0);
            const double p_noise =
                p_sig * (kWorkBandHz / kInBandHz) / snr_lin;
            const float sigma = static_cast<float>(std::sqrt(p_noise / 2.0));
            std::mt19937 rng(seed);
            std::normal_distribution<float> g(0.0f, 1.0f);
            for (auto& v : rx)
                v += gr_complex(sigma * g(rng), sigma * g(rng));
        }
        return rx;
    }

    TimingResult make_tr(int64_t start) const
    {
        TimingResult tr;
        tr.preamble_start_sample = start;
        tr.measured_period = static_cast<double>(sym);
        tr.detected_peaks = reps;
        tr.peak_samples.push_back(start + static_cast<int64_t>(sym - 1));
        tr.peak_metrics.push_back(0.5f);
        tr.peak_corr.push_back(gr_complex(1.0f, 0.0f));
        tr.ok = true;
        return tr;
    }
};

Qm35825Profile make_profile(const SfdFixture& fix, const char* mode)
{
    auto prof = Qm35825Profile::Default();
    prof.preamble_repetitions = fix.reps;
    prof.sfd_mode = mode;
    return prof;
}

struct CapSnapshot {
    std::vector<size_t> caps;
    bool operator==(const CapSnapshot& o) const { return caps == o.caps; }
};

CapSnapshot capacities(const core::DemodScratch& s)
{
    CapSnapshot snap;
    const auto add = [&snap](size_t c) { snap.caps.push_back(c); };
    add(s.work.capacity());
    add(s.derotated.capacity());
    add(s.energy.capacity());
    add(s.metric.capacity());
    add(s.corr.capacity());
    add(s.soft_chips.capacity());
    add(s.cir_taps.capacity());
    add(s.peaks.capacity());
    add(s.peak_values.capacity());
    add(s.softchip_corr.capacity());
    add(s.sfd_template.capacity());
    add(s.fold_buf.capacity());
    add(s.fold_sup_buf.capacity());
    add(s.timing_fold_cand_metrics.capacity());
    add(s.timing_fold_lag_metrics.capacity());
    add(s.timing_peak_indices.capacity());
    add(s.timing_fold_cfo_samples.capacity());
    add(s.timing_fold_cfo_corr.capacity());
    add(s.cir_avg_tail.capacity());
    add(s.cir_phase_rot.capacity());
    add(s.cfo_sym_z.capacity());
    add(s.cfo_coarse_z.capacity());
    add(s.cfo_fft_buf.capacity());
    add(s.cfo_mag_buf.capacity());
    add(s.sfd_folded.capacity());
    add(s.sfd_grid_sym_corr.capacity());
    add(s.sfd_grid_sym_pwr.capacity());
    add(s.sfd_grid_sym_valid.capacity());
    add(s.cached_seq.pcode.capacity());
    add(s.cached_seq.sfd_seq.capacity());
    add(s.cached_seq.sampled_pcode.capacity());
    add(s.cached_seq.ns_sfd_sparse.capacity());
    return snap;
}

} // namespace

// ---------------------------------------------------------------------------
// S1: one flipped (anti-phase) SFD symbol must no longer veto the whole lock.
// At -12 dB the centre-window metric sits in the 0.008..0.08 soft-polarity
// band, so the rescue really exercises check_sfd_polarity_soft().
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sfd_soft_polarity_rescues_single_flipped_symbol)
{
    const SfdFixture fix("4z2", 64);
    const Qm35825Profile prof = make_profile(fix, "4z2");
    const auto sfd_seq = fix.sfd_seq;

    const auto clean = fix.make_frame(-12.0, 1234u);
    core::DemodScratch scratch_clean;
    scratch_clean.reserve(clean.size());
    SfdResult sr_clean;
    BOOST_REQUIRE(core::stage_sfd(clean.data(), clean.size(), prof,
                                  fix.make_tr(fix.true_start), sfd_seq,
                                  fix.sync, sr_clean, scratch_clean));
    BOOST_CHECK_LE(std::llabs(sr_clean.sfd_start_sample - fix.true_sfd), 1);

    // Flip one non-zero 4z2 symbol (there are eight).  A hard all-symbol vote
    // rejects this; the coherent-weighted rule accepts 7/8 at 0.75 >= 0.65.
    const auto flipped = fix.make_frame(-12.0, 1234u, /*flip_symbol=*/0);
    core::DemodScratch scratch_flip;
    scratch_flip.reserve(flipped.size());
    SfdResult sr_flip;
    BOOST_REQUIRE(core::stage_sfd(flipped.data(), flipped.size(), prof,
                                  fix.make_tr(fix.true_start), sfd_seq,
                                  fix.sync, sr_flip, scratch_flip));
    BOOST_CHECK_LE(std::llabs(sr_flip.sfd_start_sample - fix.true_sfd), 1);
    // The centre-window metric is below the 0.08 direct floor, so the decision
    // came through the soft polarity check rather than the >=0.08 bypass.
    BOOST_CHECK_LT(sr_flip.metric, 0.08f);
    BOOST_CHECK_GE(sr_flip.metric, core::kSfdGridPeakFloor);
}

// ---------------------------------------------------------------------------
// S1: shifted-SFD hypotheses must never lock to a wrong sample, and pure noise
// must stay at 0 false alarms over 200 windows.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sfd_soft_polarity_rejects_shifted_sidelobes_and_noise)
{
    const SfdFixture fix("4z2", 32);
    const Qm35825Profile prof = make_profile(fix, "4z2");
    const auto sfd_seq = fix.sfd_seq;
    const auto frame = fix.make_frame(-8.0, 77u);

    // Timing hypotheses shifted by +-1..+-7 symbols.  The search may recover
    // the true SFD (offset windows / grid) or fail, but it must never report a
    // sidelobe position.
    for (int shift = -7; shift <= 7; ++shift) {
        if (shift == 0)
            continue;
        const int64_t start =
            fix.true_start + static_cast<int64_t>(shift) *
                                  static_cast<int64_t>(kSym);
        core::DemodScratch scratch;
        scratch.reserve(frame.size());
        SfdResult sr;
        const bool ok = core::stage_sfd(
            frame.data(), frame.size(), prof, fix.make_tr(start), sfd_seq,
            fix.sync, sr, scratch);
        if (ok) {
            const int64_t err = std::llabs(sr.sfd_start_sample - fix.true_sfd);
            BOOST_TEST_CONTEXT("shift=" << shift << " err=" << err)
            {
                BOOST_CHECK_LE(err, 2);
            }
        }
    }

    // 200 pure-noise windows: zero false alarms.
    int false_alarms = 0;
    for (int w = 0; w < 200; ++w) {
        std::vector<gr_complex> noise(fix.n, gr_complex(0.0f, 0.0f));
        std::mt19937 rng(90000u + static_cast<uint32_t>(w));
        std::normal_distribution<float> g(0.0f, 1.0f);
        for (auto& v : noise)
            v = gr_complex(g(rng), g(rng));
        core::DemodScratch scratch;
        scratch.reserve(noise.size());
        SfdResult sr;
        if (core::stage_sfd(noise.data(), noise.size(), prof,
                            fix.make_tr(fix.true_start), sfd_seq, fix.sync,
                            sr, scratch))
            ++false_alarms;
    }
    BOOST_TEST_MESSAGE("sfd refactor noise false alarms=" << false_alarms);
    BOOST_CHECK_EQUAL(false_alarms, 0);
}

// ---------------------------------------------------------------------------
// S2: a weak (-14 dB) packet whose centre-window metric is below the legacy
// 0.02 floor but above 0.008 must pass in the FIRST window -- no +-1..+-4
// symbol rescan, no discrete grid.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sfd_low_floor_hits_center_window_without_offset_scan)
{
    const SfdFixture fix("4z2", 64);
    const Qm35825Profile prof = make_profile(fix, "4z2");
    const auto sfd_seq = fix.sfd_seq;
    const auto frame = fix.make_frame(-14.0, 4242u);

    core::DemodScratch scratch;
    scratch.reserve(frame.size());
    SfdResult sr;
    BOOST_REQUIRE(core::stage_sfd(frame.data(), frame.size(), prof,
                                  fix.make_tr(fix.true_start), sfd_seq,
                                  fix.sync, sr, scratch));
    BOOST_CHECK_LE(std::llabs(sr.sfd_start_sample - fix.true_sfd), 2);
    BOOST_CHECK_EQUAL(sr.search_windows, 1u);
    // Not a >= 0.3 direct pass: this is the low-floor + polarity branch.
    BOOST_CHECK_LT(sr.metric, 0.3f);
}

// ---------------------------------------------------------------------------
// S4: rotating disjoint gaps in non-monotonic order (the incremental contract)
// must reproduce one full-range rotation sample by sample.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sfd_incremental_derotate_matches_full_rotation)
{
    Qm35825Profile prof = Qm35825Profile::Default();
    constexpr size_t n = 20000;
    std::vector<gr_complex> rx(n);
    std::mt19937 rng(5u);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (auto& v : rx)
        v = gr_complex(g(rng), g(rng));

    const double cfo_hz = 50000.0;
    core::DemodScratch full;
    full.reserve(n);
    uint64_t us = 0;
    uint32_t kernel = 0;
    core::derotate_cfo_range(rx.data(), n, prof, cfo_hz, 0, n, full, us, kernel);

    // Same total interval, applied as middle -> left gap -> right gap.
    core::DemodScratch inc;
    inc.reserve(n);
    const size_t lo = 4096, hi = 12288;
    core::derotate_cfo_range(rx.data(), n, prof, cfo_hz, lo, hi, inc, us, kernel);
    core::derotate_cfo_range(rx.data(), n, prof, cfo_hz, 0, lo, inc, us, kernel);
    core::derotate_cfo_range(rx.data(), n, prof, cfo_hz, hi, n, inc, us, kernel);

    BOOST_REQUIRE_EQUAL(full.derotated.size(), n);
    BOOST_REQUIRE_EQUAL(inc.derotated.size(), n);
    float max_err = 0.0f;
    for (size_t i = 0; i < n; ++i)
        max_err = std::max(max_err, std::abs(full.derotated[i] - inc.derotated[i]));
    BOOST_TEST_MESSAGE("incremental derotate max error=" << max_err);
    BOOST_CHECK_LE(max_err, 1e-6f);
}

// ---------------------------------------------------------------------------
// Late acquisition: timing started >36 symbols late, so the true SFD is earlier
// than expected; the (-preamble_repetitions+1 .. -37) backward grid must still
// find it through the new range-aware entry point (no-op callback here).
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sfd_late_acquisition_rescue_still_works)
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
    tr.preamble_start_sample += static_cast<int64_t>(45 * kQm35SamplesPerSymbol);

    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    const auto no_rotate = [](size_t, size_t) {};
    BOOST_REQUIRE(core::stage_sfd(iq.data(), iq.size(), prof, tr, sfd_seq, tmpl,
                                  sr, scratch, no_rotate));
    BOOST_CHECK_EQUAL(sr.sfd_start_sample, int64_t(75008));
    BOOST_CHECK_GE(sr.metric, 0.95f);
}

// ---------------------------------------------------------------------------
// S3: forcing the discrete grid twice must leave every DemodScratch capacity
// unchanged (warm path zero allocation).
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sfd_discrete_grid_warm_call_zero_allocation)
{
    const SfdFixture fix("4z2", 32);
    const Qm35825Profile prof = make_profile(fix, "4z2");
    const auto sfd_seq = fix.sfd_seq;
    // A >4-symbol timing error forces the +/-36 discrete grid (the +/-4 windows
    // miss), so the flat grid scratch is exercised.
    TimingResult tr = fix.make_tr(
        fix.true_start + static_cast<int64_t>(6 * kSym));
    const auto frame = fix.make_frame(-12.0, 909u);

    core::DemodScratch scratch;
    scratch.reserve(frame.size());
    SfdResult sr;
    const bool first_ok = core::stage_sfd(frame.data(), frame.size(), prof, tr,
                                          sfd_seq, fix.sync, sr, scratch);
    BOOST_REQUIRE(first_ok);
    BOOST_CHECK_LE(std::llabs(sr.sfd_start_sample - fix.true_sfd), 2);
    // A 6-symbol seed cannot be reached by the +-4-symbol windows, so any
    // count above 1 means the discrete-grid pass actually ran.
    BOOST_CHECK_GE(sr.search_windows, 2u);

    const CapSnapshot before = capacities(scratch);

    SfdResult sr2;
    const bool second_ok = core::stage_sfd(frame.data(), frame.size(), prof, tr,
                                           sfd_seq, fix.sync, sr2, scratch);
    BOOST_REQUIRE(second_ok);
    BOOST_CHECK_EQUAL(sr2.sfd_start_sample, sr.sfd_start_sample);

    const CapSnapshot after = capacities(scratch);
    BOOST_CHECK(before == after);
}
