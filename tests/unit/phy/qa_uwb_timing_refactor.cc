// SPDX-License-Identifier: GPL-3.0-or-later
//
// Stage-1 Timing refactor QA: the low-SNR local K=4 fold fast path, its strict
// acceptance gates, the seamless K=64 fallback, and the warm-path zero-alloc
// contract.
//
// Fixtures are 256-SYNC preambles tiled from testdata/reference_preamble.bin on
// the 998.4 MS/s WORK grid (one SYNC = 1016 samples) plus white AWGN scaled to
// an IN-BAND SNR over the 500 MHz UWB channel inside the 998.4 MHz work band --
// the same convention as qa_uwb_radar_timing_core.cc, so the numbers are
// directly comparable.
#include "demod_core_test_helpers.h"

#include <openuwb/phy/hrp/rx/demod_cfo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace {

using namespace openuwb::demod;
using openuwb::demod::core::DemodScratch;
using openuwb::demod::core::stage_timing;

constexpr size_t kReps = 256;
constexpr double kWorkBandHz = 998.4e6;
constexpr double kInBandHz = 500.0e6;
constexpr size_t kSym = kQm35SamplesPerSymbol;

std::vector<gr_complex> load_normalized_sync()
{
    std::vector<gr_complex> raw;
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", raw));
    BOOST_REQUIRE_EQUAL(raw.size(), kSym);
    auto tmpl = raw;
    openuwb::core::uwb_l2_normalize(tmpl);
    return tmpl;
}

struct LowSnrFrame {
    std::vector<gr_complex> rx;
    int64_t start = 0;
};

// Build a 256-SYNC preamble with optional SFO (symbol spacing scaled by
// 1 + ppm*1e-6) and white AWGN at `snr_db` in-band.
LowSnrFrame make_low_snr_frame(const std::vector<gr_complex>& sync,
                               int64_t pre, int64_t tail, double snr_db,
                               uint32_t seed, double sfo_ppm = 0.0)
{
    LowSnrFrame f;
    f.start = pre;
    const double scale = 1.0 + sfo_ppm * 1e-6;
    const int64_t span = static_cast<int64_t>(
        std::llround(static_cast<double>((kReps - 1) * kSym) * scale)) +
        static_cast<int64_t>(kSym);
    f.rx.assign(static_cast<size_t>(pre + span + tail),
                gr_complex(0.0f, 0.0f));
    for (size_t i = 0; i < kReps; ++i) {
        const int64_t at = pre + static_cast<int64_t>(std::llround(
            static_cast<double>(i * kSym) * scale));
        for (size_t k = 0; k < kSym; ++k)
            f.rx[static_cast<size_t>(at) + k] = sync[k];
    }
    double p_sig = 0.0;
    for (size_t i = 0; i < kSym; ++i)
        p_sig += std::norm(sync[i]);
    p_sig /= static_cast<double>(kSym);
    const double snr_lin = std::pow(10.0, snr_db / 10.0);
    const double p_noise = p_sig * (kWorkBandHz / kInBandHz) / snr_lin;
    const float sigma = static_cast<float>(std::sqrt(p_noise / 2.0));
    std::mt19937 rng(seed);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (auto& v : f.rx)
        v += gr_complex(sigma * g(rng), sigma * g(rng));
    return f;
}

void apply_cfo(std::vector<gr_complex>& rx, double cfo_hz)
{
    const double w = 2.0 * M_PI * cfo_hz / kWorkBandHz;
    for (size_t i = 0; i < rx.size(); ++i) {
        const double ph = w * static_cast<double>(i);
        rx[i] *= gr_complex(static_cast<float>(std::cos(ph)),
                            static_cast<float>(std::sin(ph)));
    }
}

Qm35825Profile low_snr_profile(size_t margin, bool sparse)
{
    auto p = Qm35825Profile::Default();
    p.preamble_repetitions = kReps;
    p.sync_strategy = SyncStrategy::AcquireLegacy;
    p.timing_search_margin = margin;
    p.timing_seed_refine_radius = 0;
    p.rx_timing_sparse = sparse;
    return p;
}

float median_peak_metric(const TimingResult& tr)
{
    if (tr.peak_metrics.empty())
        return 0.0f;
    std::array<float, 32> buf{};
    const size_t n = std::min(tr.peak_metrics.size(), buf.size());
    std::copy_n(tr.peak_metrics.begin(), n, buf.begin());
    std::nth_element(buf.begin(), buf.begin() + n / 2, buf.begin() + n);
    return buf[n / 2];
}

} // namespace

BOOST_AUTO_TEST_CASE(test_timing_k4_fast_path_hits_low_snr_accurate_seed)
{
    const auto sync = load_normalized_sync();
    const double snrs[] = { -12.0, -14.0 };
    for (double snr : snrs) {
        for (bool sparse : { false, true }) {
            const auto prof = low_snr_profile(/*margin=*/8192, sparse);
            auto frame = make_low_snr_frame(sync, 8192, 2048, snr, 771u);
            DemodScratch scratch;
            scratch.reserve(frame.rx.size());
            TimingResult tr;
            BOOST_REQUIRE_MESSAGE(
                stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                             frame.start, tr, scratch),
                "fast path failed at " << snr << " dB (sparse=" << sparse << ")");
            std::printf("k4 fast: snr=%.0f sparse=%d start_err=%lld metric=%.4f "
                        "peaks=%zu median=%.3f coarse=%llu\n",
                        snr, int(sparse),
                        static_cast<long long>(tr.preamble_start_sample -
                                               frame.start),
                        static_cast<double>(tr.metric), tr.detected_peaks,
                        static_cast<double>(median_peak_metric(tr)),
                        static_cast<unsigned long long>(tr.coarse_search_us));
            BOOST_CHECK_EQUAL(tr.coarse_search_us, uint64_t(0));
            BOOST_CHECK(tr.sync_peaks_folded);
            BOOST_CHECK_EQUAL(tr.detected_peaks, tr.peak_samples.size());
            BOOST_CHECK_LE(tr.detected_peaks, size_t(16));
            BOOST_CHECK_GE(tr.detected_peaks, size_t(4));
            BOOST_CHECK_LE(std::llabs(tr.preamble_start_sample - frame.start),
                           int64_t(2));
            BOOST_CHECK_GE(median_peak_metric(tr), 0.08f);
            BOOST_CHECK_EQUAL(tr.expected_peaks, kReps);
            BOOST_CHECK_EQUAL(tr.cir_grid_repetitions, kReps);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_timing_k4_fast_path_warm_call_does_not_grow_scratch)
{
    const auto sync = load_normalized_sync();
    const auto prof = low_snr_profile(/*margin=*/8192, /*sparse=*/true);
    auto frame = make_low_snr_frame(sync, 8192, 2048, -14.0, 991u);
    DemodScratch scratch;
    scratch.reserve(frame.rx.size());
    TimingResult tr;
    BOOST_REQUIRE(stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                               frame.start, tr, scratch));

    // Snapshot every scratch buffer that participates in Stage 1.  The second
    // (warm) call must reuse the exact same capacities: the fast path is not
    // allowed to allocate per packet.
    const auto cap = [](const auto& v) { return v.capacity(); };
    const std::array<size_t, 16> before = {
        cap(scratch.work),           cap(scratch.derotated),
        cap(scratch.metric),         cap(scratch.corr),
        cap(scratch.fold_buf),       cap(scratch.fold_sup_buf),
        cap(scratch.timing_fold_cand_metrics),
        cap(scratch.timing_fold_lag_metrics),
        cap(scratch.timing_peak_indices),
        cap(scratch.peaks),          cap(scratch.peak_values),
        cap(scratch.cir_template.entries),
        cap(scratch.shaped_sync_template.code.entries),
        cap(scratch.cfo_sym_z),
        cap(scratch.timing_fold_cfo_samples),
        cap(scratch.timing_fold_cfo_corr),
    };
    TimingResult tr2;
    BOOST_REQUIRE(stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                               frame.start, tr2, scratch));
    const std::array<size_t, 16> after = {
        cap(scratch.work),           cap(scratch.derotated),
        cap(scratch.metric),         cap(scratch.corr),
        cap(scratch.fold_buf),       cap(scratch.fold_sup_buf),
        cap(scratch.timing_fold_cand_metrics),
        cap(scratch.timing_fold_lag_metrics),
        cap(scratch.timing_peak_indices),
        cap(scratch.peaks),          cap(scratch.peak_values),
        cap(scratch.cir_template.entries),
        cap(scratch.shaped_sync_template.code.entries),
        cap(scratch.cfo_sym_z),
        cap(scratch.timing_fold_cfo_samples),
        cap(scratch.timing_fold_cfo_corr),
    };
    BOOST_CHECK(before == after);
}

BOOST_AUTO_TEST_CASE(test_timing_k4_fast_path_rejects_late_seed_boundary)
{
    // A seed deep inside the preamble (well beyond the fast path's +/-8192
    // sample local window) still folds strongly because the preamble is
    // P-periodic.  The leading-edge disambiguation then crosses at its very
    // first candidate, i.e. at the scan-window boundary, which is NOT the true
    // packet start.  The fast path MUST reject that result and fall through to
    // the wide path: if it reports coarse_search_us == 0, the origin must be
    // correct (guard failed), never a symbol-shifted origin.
    const auto sync = load_normalized_sync();
    const auto prof = low_snr_profile(/*margin=*/40960, /*sparse=*/true);
    auto frame = make_low_snr_frame(sync, 8192, 4096, -12.0, 5150u);
    DemodScratch scratch;
    scratch.reserve(frame.rx.size());
    TimingResult tr;
    const int64_t late_seed = frame.start + 20 * static_cast<int64_t>(kSym);
    const bool ok = stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                                 late_seed, tr, scratch);
    if (ok) {
        std::printf("k4 late-seed: ok=%d coarse=%llu start_err=%lld peaks=%zu\n",
                    int(ok),
                    static_cast<unsigned long long>(tr.coarse_search_us),
                    static_cast<long long>(tr.preamble_start_sample -
                                           frame.start),
                    tr.detected_peaks);
        if (tr.coarse_search_us == 0) {
            // Fast path accepted: it must have found the real leading edge.
            BOOST_CHECK_LE(
                std::llabs(tr.preamble_start_sample - frame.start), int64_t(2));
        }
        BOOST_CHECK_LT(
            std::llabs(tr.preamble_start_sample - frame.start), int64_t(40 * kSym));
    }
}

BOOST_AUTO_TEST_CASE(test_timing_k4_fast_path_pure_noise_no_false_alarm)
{
    const auto sync = load_normalized_sync();
    const auto prof = low_snr_profile(/*margin=*/8192, /*sparse=*/true);
    constexpr int kWindows = 50;
    int false_alarms = 0;
    for (int w = 0; w < kWindows; ++w) {
        std::vector<gr_complex> noise(200000);
        std::mt19937 rng(424200u + static_cast<uint32_t>(w));
        std::normal_distribution<float> g(0.0f, 1.0f);
        for (auto& v : noise)
            v = gr_complex(g(rng), g(rng));
        DemodScratch scratch;
        scratch.reserve(noise.size());
        TimingResult tr;
        if (stage_timing(noise.data(), noise.size(), prof, sync, 8192, tr,
                         scratch))
            ++false_alarms;
    }
    std::printf("k4 fast noise: %d/%d false alarms\n", false_alarms, kWindows);
    BOOST_CHECK_EQUAL(false_alarms, 0);
}

BOOST_AUTO_TEST_CASE(test_timing_fallback_survives_cfo_and_sfo)
{
    // A large CFO/SFO must not lose the packet: either K=4 still passes or the
    // full path (probe_fold_phase + K=64) rescues it.  The origin must stay
    // exact in both cases.
    const auto sync = load_normalized_sync();
    const double cfos[] = { -40000.0, 40000.0 };
    for (double cfo : cfos) {
        auto frame = make_low_snr_frame(sync, 8192, 4096, -12.0, 90210u);
        apply_cfo(frame.rx, cfo);
        const auto prof = low_snr_profile(/*margin=*/8192, /*sparse=*/true);
        DemodScratch scratch;
        scratch.reserve(frame.rx.size());
        TimingResult tr;
        BOOST_REQUIRE_MESSAGE(
            stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                         frame.start, tr, scratch),
            "CFO " << cfo << " Hz lost the packet");
        BOOST_CHECK_LE(std::llabs(tr.preamble_start_sample - frame.start),
                       int64_t(2));
    }
    for (double ppm : { -27.0, 27.0 }) {
        auto frame = make_low_snr_frame(sync, 8192, 4096, -12.0, 7331u, ppm);
        const auto prof = low_snr_profile(/*margin=*/8192, /*sparse=*/true);
        DemodScratch scratch;
        scratch.reserve(frame.rx.size());
        TimingResult tr;
        BOOST_REQUIRE_MESSAGE(
            stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                         frame.start, tr, scratch),
            "SFO " << ppm << " ppm lost the packet");
        BOOST_CHECK_LE(std::llabs(tr.preamble_start_sample - frame.start),
                       int64_t(4));
    }
}

BOOST_AUTO_TEST_CASE(test_timing_wide_interval_probe_keeps_k64_fallback)
{
    // SIC's blind profile sets rx_wide_interval_probe.  The K=4 fast path is
    // disabled there, and the historical wide-probe / K=64 fallback must still
    // acquire the same low-SNR packet when the caller opts into the fallback.
    // (The all-miss fast-fail to L0_Bypass is pinned separately by
    // qa_uwb_live_packet_processor_sic, which sets rx_wide_probe_fallback=false
    // and asserts the target is rejected without the full wide acquisition.)
    const auto sync = load_normalized_sync();
    auto prof = low_snr_profile(/*margin=*/8192, /*sparse=*/true);
    prof.rx_wide_interval_probe = true;
    prof.rx_wide_probe_step_symbols = 42;
    prof.rx_wide_probe_fold_symbols = 4;
    prof.rx_wide_probe_fallback = true;
    auto frame = make_low_snr_frame(sync, 8192, 4096, -12.0, 6006u);
    DemodScratch scratch;
    scratch.reserve(frame.rx.size());
    TimingResult tr;
    BOOST_REQUIRE(stage_timing(frame.rx.data(), frame.rx.size(), prof, sync,
                               frame.start, tr, scratch));
    BOOST_CHECK_LE(std::llabs(tr.preamble_start_sample - frame.start),
                   int64_t(2));
    BOOST_CHECK(tr.sync_peaks_folded);
}
