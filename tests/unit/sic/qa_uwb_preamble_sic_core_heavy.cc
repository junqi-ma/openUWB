// SPDX-License-Identifier: GPL-3.0-or-later
// Heavy S1 preamble-boundary scan, isolated from the rx_sic fast loop.
// Keeps the original assertions; tagged `heavy` (not `rx_sic`).
#include "preamble_sic_test_helpers.h"

BOOST_AUTO_TEST_CASE(test_sic_s1_preamble_boundary)
{
    for (int64_t off : { -3, -1, 0, 1, 3 }) {
        Fixture f =
            make_fixture(10, 64, static_cast<int64_t>(kPre) + off, 64, 64, 4, 8, 8, 100);
        const std::vector<Path> paths = { { 0.4, cplx(1.f, 0.f) },
                                          { 11.4, cplx(0.f, 0.5f) } };
        Signal s = build_signal(f.model,
                                f.n,
                                f.comm_start,
                                paths,
                                0.0,
                                kQaNoise,
                                410 + static_cast<uint64_t>(off + 3));
        const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
        f.cfg.observed_end_sample = oe;
        SicResult r = run(f, s);
        const double supp = actual_suppression_db(s.i, s.y, f.out, 0, f.n);
        std::printf("[sic-qa] s1-boundary off=%lld state=%s ds=%lld np=%zu "
                    "support=[%lld,%lld) actual=%.2f dB fit_res=%.4f\n",
                    static_cast<long long>(off),
                    state_name(r.state),
                    static_cast<long long>(r.detected_start),
                    r.n_paths,
                    static_cast<long long>(r.support_lo),
                    static_cast<long long>(r.support_hi),
                    supp,
                    r.fit_residual);
        BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                          static_cast<int>(SicState::Committed));
        BOOST_CHECK_LE(std::abs(r.detected_start - f.comm_start), 1);
        BOOST_CHECK_GT(supp, 15.0);
        BOOST_CHECK_EQUAL(std::string(r.support_source ? r.support_source : ""),
                          std::string("configured_nominal_prediction"));
    }

    // Nominal vs measured support: the profile claims 64 SYNC but the actual
    // interferer is only 20 SYNC long.  The NOMINAL support must still be the
    // full 64-symbol prediction (it is a prediction, not a measurement), which
    // is exactly the over-reach the extrapolation test quantifies.
    {
        Fixture f =
            make_fixture(10, 64, static_cast<int64_t>(kPre) + 2, 64, 64, 4, 8, 8, 100);
        Model gen20 = build_model(10, 20);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(gen20, f.n, f.comm_start, paths, 0.0, kQaNoise, 429);
        const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
        f.cfg.observed_end_sample = oe;
        SicModel m = estimate_model_checked(f, s, oe, nullptr);
        BOOST_REQUIRE(m.valid);
        std::printf(
            "[sic-qa] s1-boundary nominal-vs-measured: nominal support "
            "len=%lld (measured interferer len=%lld) over-reach=%lld\n",
            static_cast<long long>(m.predicted_support_hi - m.predicted_support_lo),
            static_cast<long long>(gen20.preamble.size()),
            static_cast<long long>((m.predicted_support_hi - m.predicted_support_lo) -
                                   static_cast<int64_t>(gen20.preamble.size())));
        BOOST_CHECK_EQUAL(m.predicted_support_hi - m.predicted_support_lo,
                          64 * static_cast<int64_t>(kSps) + 100);
        BOOST_CHECK_GT((m.predicted_support_hi - m.predicted_support_lo),
                       20 * static_cast<int64_t>(kSps));
    }

    // Window clipping: a window shorter than the nominal preamble clips the
    // prediction to [0,n) but still commits and rebuilds the visible part.
    {
        Fixture f =
            make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 430);
        const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
        f.cfg.observed_end_sample = oe;
        SicModel m = estimate_model_checked(f, s, oe, nullptr);
        BOOST_REQUIRE(m.valid);
        // A deliberately short apply chunk must stay inside the chunk.
        std::vector<cplx> out(1000, cplx(0.f, 0.f));
        std::vector<cplx> ih(1000, cplx(0.f, 0.f));
        const int64_t off = f.comm_start - 100;
        BOOST_REQUIRE(sic::apply_model(m,
                                       s.y.data() + off,
                                       1000,
                                       off,
                                       f.templates(),
                                       f.cfg,
                                       f.scratch,
                                       out.data(),
                                       ih.data(),
                                       nullptr));
        size_t changed = 0;
        for (size_t k = 0; k < 1000; ++k)
            changed += (out[k] != s.y[static_cast<size_t>(off) + k]) ? 1 : 0;
        std::printf("[sic-qa] s1-boundary clipped-chunk: 1000-sample chunk at "
                    "offset %lld changed %zu samples\n",
                    static_cast<long long>(off),
                    changed);
        BOOST_CHECK_GT(changed, size_t(0));
        // Outside the mask the chunk is byte-identical.
        for (size_t k = 0; k < 1000; ++k) {
            const int64_t abs = off + static_cast<int64_t>(k);
            if (abs < m.sub_lo || abs >= m.sub_hi)
                BOOST_CHECK_EQUAL(out[k], s.y[static_cast<size_t>(abs)]);
        }
    }
}

// ---------------------------------------------------------------------------
// S1-7. Extrapolation mismatch: the model must NOT just "commit"; the
// over-subtraction and the target damage must be MEASURED.
// ---------------------------------------------------------------------------
