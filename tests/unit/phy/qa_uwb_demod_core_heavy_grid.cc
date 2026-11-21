// SPDX-License-Identifier: GPL-3.0-or-later
// Heavy low-SNR SFD discrete-grid scan, isolated from the rx_sic fast loop.
// Split into per-SNR and noise cases so CTest can run all four in parallel on
// one compiled binary while keeping the original seeds and assertions identical.
#include "demod_core_test_helpers.h"

namespace {

struct HeavyGridFixture {
    std::vector<gr_complex> sync;
    std::vector<int8_t> sfd_seq;
    std::vector<gr_complex> sfd_wf;
    size_t sym = kQm35SamplesPerSymbol;
    size_t reps = 256;
    int64_t pre = 8192;
    int64_t true_start = 8192;
    int64_t true_sfd = 0;
    size_t n = 0;
    double p_sig = 0.0;

    HeavyGridFixture()
    {
        BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", sync));
        BOOST_REQUIRE_EQUAL(sync.size(), kQm35SamplesPerSymbol);
        sfd_seq = GetSfdSequence("4z2");
        sfd_wf.resize(sfd_seq.size() * sym);
        for (size_t i = 0; i < sfd_seq.size(); ++i)
            for (size_t k = 0; k < sym; ++k)
                sfd_wf[i * sym + k] = static_cast<float>(sfd_seq[i]) * sync[k];
        true_sfd = true_start + static_cast<int64_t>(reps * sym);
        n = static_cast<size_t>(true_sfd) + sfd_wf.size() + 4096;
        for (size_t i = 0; i < sym; ++i)
            p_sig += std::norm(sync[i]);
        p_sig /= static_cast<double>(sym);
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

void run_snr_sweep(double snr, bool check_shifted_backtrack)
{
    const HeavyGridFixture fix;
    auto prof = Qm35825Profile::Default();
    prof.preamble_repetitions = fix.reps;
    prof.sfd_mode = "4z2";

    constexpr int kFrames = 50;
    const double snr_lin = std::pow(10.0, snr / 10.0);
    const double p_noise = fix.p_sig * (998.4 / 500.0) / snr_lin;
    const float sigma = static_cast<float>(std::sqrt(p_noise / 2.0));

    for (bool fast_sfd : { false, true }) {
        prof.rx_fast_sfd = fast_sfd;
        int passed = 0;
        int within2 = 0;
        int64_t worst = 0;
        for (int f = 0; f < kFrames; ++f) {
            std::vector<gr_complex> rx(fix.n, gr_complex(0.0f, 0.0f));
            for (size_t i = 0; i < fix.reps; ++i)
                std::copy(fix.sync.begin(),
                          fix.sync.end(),
                          rx.begin() + fix.true_start +
                              static_cast<int64_t>(i * fix.sym));
            std::copy(fix.sfd_wf.begin(),
                      fix.sfd_wf.end(),
                      rx.begin() + fix.true_sfd);
            std::mt19937 rng(9000u + static_cast<uint32_t>(f));
            std::normal_distribution<float> g(0.0f, 1.0f);
            for (auto& v : rx)
                v += gr_complex(sigma * g(rng), sigma * g(rng));
            if (check_shifted_backtrack && f == 0) {
                auto shifted = fix.make_tr(
                    fix.true_start - 4 * static_cast<int64_t>(fix.sym) - 20);
                shifted.metric = 0.04f;
                core::DemodScratch shifted_scratch;
                SfdResult shifted_sfd;
                BOOST_REQUIRE(core::stage_sfd(rx.data(),
                                              rx.size(),
                                              prof,
                                              shifted,
                                              fix.sfd_seq,
                                              fix.sync,
                                              shifted_sfd,
                                              shifted_scratch));
                BOOST_CHECK_EQUAL(shifted_sfd.sfd_start_sample, fix.true_sfd);
                BOOST_CHECK_EQUAL(shifted_sfd.first_threshold_backtrack_symbols, 4);
            }
            const int64_t off =
                static_cast<int64_t>(static_cast<int>(f % 21) - 10) *
                static_cast<int64_t>(fix.sym);
            TimingResult tr = fix.make_tr(fix.true_start + off);
            core::DemodScratch scratch;
            scratch.reserve(rx.size());
            SfdResult sr;
            if (!core::stage_sfd(rx.data(),
                                 rx.size(),
                                 prof,
                                 tr,
                                 fix.sfd_seq,
                                 fix.sync,
                                 sr,
                                 scratch))
                continue;
            ++passed;
            const int64_t err = std::llabs(sr.sfd_start_sample - fix.true_sfd);
            worst = std::max(worst, err);
            if (err <= 2)
                ++within2;
        }
        std::printf("sfd low-snr: fast=%d snr=%.0f passed=%d/%d within2=%d "
                    "worst=%lld\n",
                    fast_sfd ? 1 : 0,
                    snr,
                    passed,
                    kFrames,
                    within2,
                    static_cast<long long>(worst));
        BOOST_CHECK_GE(passed, 48);
        BOOST_CHECK_GE(within2, 48);
    }
}

void run_noise_sweep(bool fast_sfd)
{
    const HeavyGridFixture fix;
    auto prof = Qm35825Profile::Default();
    prof.preamble_repetitions = fix.reps;
    prof.sfd_mode = "4z2";
    prof.rx_fast_sfd = fast_sfd;

    int false_alarms = 0;
    for (int w = 0; w < 200; ++w) {
        std::vector<gr_complex> rx(fix.n, gr_complex(0.0f, 0.0f));
        std::mt19937 rng(70000u + static_cast<uint32_t>(w));
        std::normal_distribution<float> g(0.0f, 1.0f);
        for (auto& v : rx)
            v = gr_complex(g(rng), g(rng));
        TimingResult tr = fix.make_tr(fix.true_start);
        core::DemodScratch scratch;
        scratch.reserve(rx.size());
        SfdResult sr;
        if (core::stage_sfd(rx.data(),
                            rx.size(),
                            prof,
                            tr,
                            fix.sfd_seq,
                            fix.sync,
                            sr,
                            scratch))
            ++false_alarms;
    }
    std::printf("sfd low-snr noise: fast=%d %d/200 false alarms\n",
                fast_sfd ? 1 : 0,
                false_alarms);
    BOOST_CHECK_EQUAL(false_alarms, 0);
}

} // namespace

BOOST_AUTO_TEST_CASE(test_demod_core_low_snr_sfd_discrete_grid_m10db)
{
    run_snr_sweep(-10.0, false);
}

BOOST_AUTO_TEST_CASE(test_demod_core_low_snr_sfd_discrete_grid_m12db)
{
    run_snr_sweep(-12.0, true);
}

BOOST_AUTO_TEST_CASE(test_demod_core_low_snr_sfd_discrete_grid_m14db)
{
    run_snr_sweep(-14.0, false);
}

BOOST_AUTO_TEST_CASE(test_demod_core_low_snr_sfd_discrete_grid_noise_scalar)
{
    run_noise_sweep(false);
}

BOOST_AUTO_TEST_CASE(test_demod_core_low_snr_sfd_discrete_grid_noise_fast)
{
    run_noise_sweep(true);
}
