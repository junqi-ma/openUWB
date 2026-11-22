/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * QA for the discrete equivalent CIR full-packet reconstruction.
 *
 * The oracle here is an INDEPENDENT double-precision direct convolution with
 * an absolute-phase CFO loop written from the header math, not a golden file
 * and not the estimator.  Packets are built locally from the public modulator
 * (or, where the modulator would hide the arithmetic, from a deterministic
 * synthetic waveform).  No sidecar and no injected truth is used.
 */
#include <boost/test/unit_test.hpp>
#include <gnuradio/uwb/uwb_discrete_cir_sic.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace openuwb;
using cf = std::complex<float>;
using cd = std::complex<double>;

constexpr double kPi = 3.14159265358979323846;
constexpr double kWorkRateHz = 998.4e6;
constexpr size_t kSps = 1016;
constexpr size_t kSyncSamples = 1016; // one shaped SYNC on the work grid
constexpr float kRefScale = 0.5f;     // exactly representable -> exact template

// Independent reconstruction: y_hat[j] = exp(j*w*(j-r)) * sum_l h[l]*x[j-q-l]
// with the single waveform scale applied to x in double.  It walks the DENSE
// hull (inactive taps are zero), so it also checks that the active set matches
// the hull the model declares.
std::vector<cf> independent_reconstruction(const std::vector<cf>& wave,
                                           const sic::DiscreteCirModel& m,
                                           double scale,
                                           size_t n)
{
    std::vector<cf> ref(n, cf(0.f, 0.f));
    const int64_t lo = m.q + m.lag_lo;
    const int64_t hi = lo + int64_t(m.tap_count) + int64_t(wave.size());
    const double omega = 2.0 * kPi * m.cfo_hz / kWorkRateHz;
    for (int64_t j = lo; j < hi; ++j) {
        if (j < 0 || j >= int64_t(n))
            continue;
        cd acc(0.0, 0.0);
        for (size_t i = 0; i < m.tap_count; ++i) {
            const int64_t x = j - m.q - (m.lag_lo + int64_t(i));
            if (x < 0 || x >= int64_t(wave.size()))
                continue;
            acc += cd(m.taps[i]) * (cd(wave[static_cast<size_t>(x)]) * scale);
        }
        const cd rot =
            std::exp(cd(0.0, omega * static_cast<double>(j - m.phase_reference_sample)));
        ref[static_cast<size_t>(j)] = cf(acc * rot);
    }
    return ref;
}

sic::DiscreteCirModel make_model(int64_t q, int64_t lag_lo, size_t taps)
{
    sic::DiscreteCirModel m;
    m.valid = true;
    m.status = "committed";
    m.model_id = 0x5eed;
    m.q = q;
    m.phase_reference_sample = q;
    m.observed_end_sample = q + int64_t(8 * kSps);
    m.lag_lo = lag_lo;
    m.tap_count = taps;
    m.fit_residual = 0.01;
    m.gain_condition = 12.0;
    m.effective_rank = taps;
    m.support_source = "coarse_cir_peak_groups";
    return m;
}

// Active set in increasing lag order; every lag must lie inside the dense hull.
void set_active(sic::DiscreteCirModel& m,
                std::initializer_list<std::pair<int64_t, cf>> entries)
{
    m.active_count = 0;
    for (const auto& e : entries) {
        const int64_t idx = e.first - m.lag_lo;
        BOOST_REQUIRE_GE(idx, 0);
        BOOST_REQUIRE_LT(idx, int64_t(m.tap_count));
        BOOST_REQUIRE_LT(m.active_count, sic::kDiscreteCirMaxActive);
        m.active_lags[m.active_count++] = static_cast<int16_t>(e.first);
        m.taps[static_cast<size_t>(idx)] = e.second;
    }
}

size_t count_nonzero(const std::vector<cf>& v, int64_t lo, int64_t hi)
{
    size_t bad = 0;
    for (int64_t j = 0; j < int64_t(v.size()); ++j)
        if ((j < lo || j >= hi) && v[static_cast<size_t>(j)] != cf(0.f, 0.f))
            ++bad;
    return bad;
}

// A locally modulated packet.  `reference` is the public one-SYNC reference
// scaled by kRefScale, so the exact waveform scale is kRefScale.
struct Packet {
    mod::HrpModConfig profile;
    std::vector<uint8_t> psdu; // includes the FCS appended by the QA
    std::vector<uint8_t> phr;
    std::vector<cf> wave;

    Packet(size_t wire_len, uint32_t seed)
    {
        profile.code_index = 10;
        profile.sync_repetitions = 64;
        profile.sfd_mode = "4z2";
        profile.peak_amplitude = 0.8f;
        BOOST_REQUIRE_GT(wire_len, 2);
        psdu = mod::make_random_psdu(wire_len - 2, seed, true);
        BOOST_REQUIRE_EQUAL(psdu.size(), wire_len);
        std::vector<int8_t> bits;
        mod::detail::encode_phr19(psdu.size(), profile.sync_repetitions, false, bits);
        phr.assign(bits.begin(), bits.end());
        mod::HrpModScratch ms;
        const size_t cap =
            mod::packet_samples_998p4(64, 8, psdu.size(), false, &profile.pulse) + 64;
        ms.reserve(127, cap);
        wave.resize(cap);
        size_t produced = 0;
        BOOST_REQUIRE(mod::modulate_one(psdu.data(), psdu.size(), profile, ms,
                                        wave.data(), cap, produced));
        BOOST_REQUIRE_GT(produced, kSyncSamples);
        wave.resize(produced);
    }

    std::vector<cf> reference(float k = kRefScale) const
    {
        std::vector<cf> r(wave.begin(), wave.begin() + kSyncSamples);
        for (auto& v : r)
            v *= k;
        return r;
    }
};

// Apply-level rig with a synthetic (modulator-free) observation window.
struct ApplyRig {
    size_t n;
    int64_t q;
    sic::DiscreteCirModel model;
    std::vector<cf> wave, ref, y, out, ihat;
    sic::DiscreteCirScratch scratch;

    ApplyRig(size_t wave_n, int64_t origin)
        : n(wave_n + 4096), q(origin), model(make_model(origin, -1, 3))
    {
        wave.assign(wave_n, cf(0.8f, 0.f));
        ref.assign(wave.begin(), wave.begin() + kSyncSamples);
        for (auto& v : ref)
            v *= kRefScale;
        y.assign(n, cf(0.01f, 0.f));
        out.assign(n, cf(0.f, 0.f));
        ihat.assign(n, cf(0.f, 0.f));
        set_active(model, { { -1, cf(1.f, 0.f) },
                            { 0, cf(0.5f, 0.f) },
                            { 1, cf(0.25f, 0.f) } });
        reset_scratch();
    }

    void reset_scratch()
    {
        scratch.prepared = true;
        scratch.sps = kSps;
        scratch.max_window_samples = n;
        scratch.accumulator.assign(n, cf(0.f, 0.f));
        scratch.fir_tmp.assign(n, cf(0.f, 0.f));
        scratch.max_packet_samples = wave.size() + 64;
        scratch.scaled_waveform.assign(scratch.max_packet_samples, cf(0.f, 0.f));
    }

    bool run(sic::DiscreteCirResult& r)
    {
        return sic::apply_discrete_cir(y.data(), n, wave.data(), wave.size(),
                                       ref.data(), ref.size(), model, 0,
                                       int64_t(n), scratch, out.data(),
                                       ihat.data(), &r);
    }

    void expect_bypass() const
    {
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }
};

// Transactional rig: the observation is the remodulated packet through one
// integer-lag tap plus a CFO, so a committed run must actually cancel it.
struct TransactionRig {
    Packet pkt;
    size_t n;
    int64_t q;
    sic::DiscreteCirModel model;
    std::vector<cf> ref, y, out, ihat;
    sic::DiscreteCirScratch scratch;

    TransactionRig(size_t wire_len, uint32_t seed)
        : pkt(wire_len, seed),
          n(pkt.wave.size() + 4096),
          q(2048),
          model(make_model(2048, 3, 1))
    {
        set_active(model, { { 3, cf(0.5f, 0.25f) } });
        model.cfo_hz = 1500.0;
        ref = pkt.reference();
        y.assign(n, cf(0.002f, -0.001f));
        const double omega = 2.0 * kPi * model.cfo_hz / kWorkRateHz;
        for (size_t i = 0; i < pkt.wave.size(); ++i) {
            // CFO phase is anchored at r == q, so the observed packet sample at
            // output index q+3+i carries exp(j*omega*(3+i)).
            const double rel = 3.0 + double(i);
            const cf rot = cf(static_cast<float>(std::cos(omega * rel)),
                              static_cast<float>(std::sin(omega * rel)));
            // ihat = (kRefScale * wave) * h, h = 0.5 + 0.25j
            y[static_cast<size_t>(q + 3 + int64_t(i))] +=
                pkt.wave[i] * cf(kRefScale * 0.5f, kRefScale * 0.25f) * rot;
        }
        out.assign(n, cf(0.f, 0.f));
        ihat.assign(n, cf(0.f, 0.f));
        scratch.prepared = true;
        scratch.sps = kSps;
        scratch.max_window_samples = n;
        scratch.accumulator.assign(n, cf(0.f, 0.f));
        scratch.fir_tmp.assign(n, cf(0.f, 0.f));
        scratch.max_packet_samples = pkt.wave.size() + 64;
        scratch.waveform.assign(scratch.max_packet_samples, cf(0.f, 0.f));
        scratch.scaled_waveform.assign(scratch.max_packet_samples, cf(0.f, 0.f));
        scratch.expected_phr.reserve(19);
        scratch.modulator.reserve(127, scratch.max_packet_samples);
    }

    sic::DiscreteCirResult cancel(bool fcs = true,
                                  const mod::HrpModConfig* prof = nullptr,
                                  const std::vector<uint8_t>* phr = nullptr,
                                  const std::vector<uint8_t>* psdu = nullptr)
    {
        return sic::fullpacket_cancel_discrete(
            y.data(), n, phr ? *phr : pkt.phr, psdu ? *psdu : pkt.psdu, fcs,
            prof ? *prof : pkt.profile, model, ref.data(), ref.size(), 0,
            int64_t(n), scratch, out.data(), ihat.data());
    }

    void expect_bypass() const
    {
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }
};

// Sparse-pulse rig: a unit-impulse packet (External pulse {1, 0}) plus the
// Ipatov scratch the sparse reconstruction consumes.  `events` is a stable
// copy because modulate_pulse_list() returns a view into the modulator
// scratch.
struct SparseRig {
    std::vector<uint8_t> wire;
    mod::HrpModConfig unit;
    mod::HrpModScratch ms;
    mod::PulseList pl;
    std::vector<mod::PulseEvent> events;
    std::vector<cf> unit_wave; // the same packet as a real unit-impulse waveform
    size_t n = 0;
    int64_t q = 2048;
    sic::DiscreteCirScratch scratch;

    explicit SparseRig(size_t wire_len = 30, uint32_t seed = 909)
    {
        wire = mod::make_random_psdu(wire_len, seed, true);
        unit.code_index = 10;
        unit.sync_repetitions = 64;
        unit.sfd_mode = "4z2";
        unit.peak_amplitude = 1.f;
        unit.pulse.shape = mod::PulseShape::External;
        unit.pulse.external_n_taps = 2;
        unit.external_taps = { 1.f, 0.f };
        const size_t cap = 300000;
        ms.reserve(127, cap);
        BOOST_REQUIRE(
            mod::modulate_pulse_list(wire.data(), wire.size(), unit, ms, cap, pl));
        events.assign(pl.events, pl.events + pl.count);
        pl.events = events.data();
        n = pl.packet_samples + 4096;
        unit_wave.assign(pl.packet_samples, cf{});
        for (const auto& e : events)
            unit_wave[e.sample_offset] = cf(float(e.sign), 0.f);
        reset_scratch();
    }

    void reset_scratch()
    {
        scratch.reserve(n);
        scratch.prepared = true;
        scratch.sps = kSps;
        scratch.ipatov_correlation = true;
        scratch.sparse_pulse_reconstruct = true;
        scratch.impulse_sync.assign(kSps, cf{});
        const int8_t* code = demod::GetPreambleCode(10);
        for (size_t c = 0; c < demod::kQm35CodeLength; ++c)
            scratch.impulse_sync[8 * c] = cf(float(code[c]), 0.f);
        scratch.max_packet_samples = pl.packet_samples + 64;
    }

    sic::DiscreteCirModel model(int64_t origin, int64_t lag_lo,
                                size_t taps) const
    {
        sic::DiscreteCirModel m = make_model(origin, lag_lo, taps);
        m.impulse_template = true;
        return m;
    }
};

} // namespace

BOOST_AUTO_TEST_CASE(discrete_apply_matches_independent_direct_convolution)
{
    Packet pkt(32, 11);
    const std::vector<cf> ref = pkt.reference();
    const int64_t q = 2048;
    const size_t n = pkt.wave.size() + 4096;

    for (double cfo : { 0.0, -50000.0, 50000.0 }) {
        sic::DiscreteCirModel m = make_model(q, -2, 10);
        m.cfo_hz = cfo;
        set_active(m, { { -2, cf(0.9f, -0.2f) },  // lowest lag, negative
                        { 0, cf(0.4f, 0.25f) },   // packet origin
                        { 7, cf(-0.15f, 0.1f) } }); // inside the hull, far right

        sic::DiscreteCirScratch scratch;
        scratch.prepared = true;
        scratch.sps = kSps;
        scratch.max_window_samples = n;
        scratch.accumulator.assign(n, cf(0.f, 0.f));
        scratch.fir_tmp.assign(n, cf(0.f, 0.f));
        scratch.max_packet_samples = pkt.wave.size() + 64;
        scratch.scaled_waveform.assign(scratch.max_packet_samples, cf(0.f, 0.f));

        std::vector<cf> y(n, cf(0.001f, -0.002f)), out(n), ihat(n);
        const double omega = 2.0 * kPi * cfo / kWorkRateHz;
        for (size_t i = 0; i < pkt.wave.size(); ++i) {
            const cf rot = cf(static_cast<float>(std::cos(omega * double(i))),
                              static_cast<float>(std::sin(omega * double(i))));
            y[static_cast<size_t>(q + int64_t(i))] += pkt.wave[i] * 0.5f * rot;
        }

        sic::DiscreteCirResult r;
        BOOST_REQUIRE_MESSAGE(
            sic::apply_discrete_cir(y.data(), n, pkt.wave.data(), pkt.wave.size(),
                                    ref.data(), ref.size(), m, 0, int64_t(n),
                                    scratch, out.data(), ihat.data(), &r),
            r.reason);
        BOOST_CHECK(r.committed);
        BOOST_CHECK_EQUAL(std::string(r.reason), "committed");
        BOOST_CHECK_CLOSE(r.waveform_scale, double(kRefScale), 1e-6);
        BOOST_CHECK_EQUAL(r.packet_samples, pkt.wave.size());
        BOOST_CHECK_EQUAL(r.model_id, m.model_id);
        BOOST_CHECK_EQUAL(r.q, q);
        BOOST_CHECK_EQUAL(r.tap_count, m.tap_count);
        BOOST_CHECK(r.scale_us >= 0.0);
        BOOST_CHECK(r.convolution_us >= 0.0);
        BOOST_CHECK(r.cfo_us >= 0.0);
        BOOST_CHECK(r.subtract_us >= 0.0);
        BOOST_CHECK_GE(r.channel_us, r.scale_us + r.convolution_us + r.cfo_us);
        BOOST_CHECK_GE(r.apply_total_us, r.channel_us + r.subtract_us);

        const std::vector<cf> want =
            independent_reconstruction(pkt.wave, m, double(kRefScale), n);
        const int64_t lo = q + m.lag_lo;
        const int64_t hi = lo + int64_t(m.tap_count) + int64_t(pkt.wave.size());
        BOOST_CHECK_EQUAL(r.support_lo, lo);
        BOOST_CHECK_EQUAL(r.support_hi, hi);
        BOOST_CHECK_EQUAL(r.sub_lo, lo);
        BOOST_CHECK_EQUAL(r.sub_hi, hi);

        float worst = 0.f;
        size_t leaks = 0, out_bad = 0, sub_bad = 0;
        for (int64_t j = 0; j < int64_t(n); ++j) {
            const size_t s = static_cast<size_t>(j);
            worst = std::max(worst, std::abs(ihat[s] - want[s]));
            if ((j < lo || j >= hi) && ihat[s] != cf(0.f, 0.f))
                ++leaks;
            if (j < r.sub_lo || j >= r.sub_hi) {
                if (!(out[s] == y[s]))
                    ++out_bad;
            } else if (std::abs(out[s] - (y[s] - ihat[s])) > 2e-6f) {
                ++sub_bad;
            }
        }
        // float FIR + float absolute-grid rotation vs a double oracle.  The
        // measured worst case over the three CFO cells is ~9e-8, so 1e-5 leaves
        // more than 100x margin for compiler/vectorization differences.
        BOOST_CHECK_LT(worst, 1e-5f);
        BOOST_CHECK_EQUAL(leaks, 0);
        BOOST_CHECK_EQUAL(out_bad, 0);
        BOOST_CHECK_EQUAL(sub_bad, 0);
    }
}

BOOST_AUTO_TEST_CASE(discrete_apply_tail_and_unaligned_offsets)
{
    // 1019 = 254 four-wide vectors + one two-wide + one scalar remainder.
    // Lags -1/0/1 shift src/dst off the 16-byte vector boundary, so this
    // checks the fused AXPY remainder path and unaligned loads/stores against
    // the independent double convolution.
    ApplyRig rig(1019, 2048);
    sic::DiscreteCirResult r;
    BOOST_REQUIRE_MESSAGE(rig.run(r), r.reason);
    BOOST_CHECK(r.committed);
    const std::vector<cf> want =
        independent_reconstruction(rig.wave, rig.model, double(kRefScale), rig.n);
    float worst = 0.f;
    for (size_t j = 0; j < rig.n; ++j)
        worst = std::max(worst, std::abs(rig.ihat[j] - want[j]));
    BOOST_CHECK_LT(worst, 1e-6f);
}

BOOST_AUTO_TEST_CASE(discrete_apply_truncation_roi_and_packet_tail)
{
    Packet pkt(24, 13);
    const std::vector<cf> ref = pkt.reference();
    const size_t n = pkt.wave.size() + 1024;
    std::vector<cf> y(n, cf(0.01f, 0.02f)), out(n), ihat(n);

    sic::DiscreteCirScratch scratch;
    scratch.prepared = true;
    scratch.sps = kSps;
    scratch.max_window_samples = n;
    scratch.accumulator.assign(n, cf(0.f, 0.f));
    scratch.fir_tmp.assign(n, cf(0.f, 0.f));
    scratch.max_packet_samples = pkt.wave.size() + 64;
    scratch.scaled_waveform.assign(scratch.max_packet_samples, cf(0.f, 0.f));

    // q near the end: a complete observed packet is mandatory, never clipped.
    {
        sic::DiscreteCirModel m = make_model(int64_t(n) - 512, -2, 10);
        set_active(m, { { 0, cf(0.5f, 0.f) } });
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                             pkt.wave.size(), ref.data(),
                                             ref.size(), m, 0, int64_t(n),
                                             scratch, out.data(), ihat.data(), &r));
        BOOST_CHECK(!r.committed);
        BOOST_CHECK_EQUAL(std::string(r.reason), "context_truncated");
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }

    // Negative lowest lag with a packet origin smaller than |lag_lo| would put
    // the support before sample 0: also a truncation, not a silent clip.
    {
        sic::DiscreteCirModel m = make_model(1, -8, 10);
        set_active(m, { { 0, cf(0.5f, 0.f) } });
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                             pkt.wave.size(), ref.data(),
                                             ref.size(), m, 0, int64_t(n),
                                             scratch, out.data(), ihat.data(), &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "context_truncated");
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }

    // Consumer ROI that ends at the support start has no overlap.
    {
        sic::DiscreteCirModel m = make_model(64, -2, 10);
        set_active(m, { { 0, cf(0.5f, 0.f) } });
        const int64_t lo = m.q + m.lag_lo;
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                             pkt.wave.size(), ref.data(),
                                             ref.size(), m, 0, lo,
                                             scratch, out.data(), ihat.data(), &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "no_overlap");
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }

    // A waveform shorter than the reference SYNC cannot be scaled and must not
    // read past its end.
    {
        sic::DiscreteCirModel m = make_model(64, -2, 10);
        set_active(m, { { 0, cf(0.5f, 0.f) } });
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                             kSyncSamples - 1, ref.data(),
                                             ref.size(), m, 0, int64_t(n),
                                             scratch, out.data(), ihat.data(), &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "invalid_buffers");
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }

    // Packet tail exactly at the window end commits, and a narrower ROI still
    // subtracts only inside its own interval.
    {
        const int64_t wave_n = int64_t(pkt.wave.size());
        // hi = q + lag_lo + tap_count + wave_n = q - 2 + 10 + wave_n.
        const int64_t q = int64_t(n) - 8 - wave_n; // hi == n
        sic::DiscreteCirModel m = make_model(q, -2, 10);
        set_active(m, { { 0, cf(0.5f, 0.f) } });
        sic::DiscreteCirResult r;
        const int64_t roi_hi = int64_t(n) - 500;
        BOOST_REQUIRE_MESSAGE(sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                                      pkt.wave.size(), ref.data(),
                                                      ref.size(), m, 0, roi_hi,
                                                      scratch, out.data(), ihat.data(),
                                                      &r),
                              r.reason);
        BOOST_CHECK_EQUAL(r.support_hi, int64_t(n));
        BOOST_CHECK_EQUAL(r.sub_hi, roi_hi);
        BOOST_CHECK_LT(r.sub_hi, r.support_hi);
        size_t out_bad = 0;
        for (int64_t j = 0; j < int64_t(n); ++j)
            if ((j < r.sub_lo || j >= r.sub_hi) &&
                !(out[static_cast<size_t>(j)] == y[static_cast<size_t>(j)]))
                ++out_bad;
        BOOST_CHECK_EQUAL(out_bad, 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, r.support_lo, r.support_hi), 0);
    }
}

BOOST_AUTO_TEST_CASE(discrete_apply_scale_once_and_wrong_reference_bypass)
{
    Packet pkt(32, 17);
    const std::vector<cf> ref = pkt.reference();
    const size_t n = pkt.wave.size() + 2048;

    sic::DiscreteCirModel m = make_model(1024, -1, 5);
    set_active(m, { { -1, cf(0.7f, 0.1f) },
                    { 0, cf(0.2f, -0.05f) },
                    { 3, cf(0.05f, 0.f) } });
    m.cfo_hz = 12000.0;

    sic::DiscreteCirScratch scratch;
    scratch.prepared = true;
    scratch.sps = kSps;
    scratch.max_window_samples = n;
    scratch.accumulator.assign(n, cf(0.f, 0.f));
    scratch.fir_tmp.assign(n, cf(0.f, 0.f));
    scratch.max_packet_samples = pkt.wave.size() + 64;
    scratch.scaled_waveform.assign(scratch.max_packet_samples, cf(0.f, 0.f));

    std::vector<cf> y(n, cf(0.0f, 0.0f));
    std::vector<cf> out1(n), ihat1(n), out3(n), ihat3(n);

    sic::DiscreteCirResult r1;
    BOOST_REQUIRE_MESSAGE(sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                                  pkt.wave.size(), ref.data(),
                                                  ref.size(), m, 0, int64_t(n),
                                                  scratch, out1.data(), ihat1.data(),
                                                  &r1),
                          r1.reason);

    // Same packet observed three times louder: h is fitted against the public
    // template, so the recovered scale must drop by 3 and the reconstruction
    // must be unchanged (taps are not scaled a second time).
    std::vector<cf> wave3 = pkt.wave;
    for (auto& v : wave3)
        v *= 3.f;
    sic::DiscreteCirResult r3;
    BOOST_REQUIRE_MESSAGE(sic::apply_discrete_cir(y.data(), n, wave3.data(),
                                                  wave3.size(), ref.data(),
                                                  ref.size(), m, 0, int64_t(n),
                                                  scratch, out3.data(), ihat3.data(),
                                                  &r3),
                          r3.reason);
    BOOST_CHECK_CLOSE(r1.waveform_scale, 3.0 * r3.waveform_scale, 1e-3);

    float worst = 0.f, worst_out = 0.f;
    for (size_t j = 0; j < n; ++j) {
        worst = std::max(worst, std::abs(ihat1[j] - ihat3[j]));
        worst_out = std::max(worst_out, std::abs(out1[j] - out3[j]));
    }
    BOOST_CHECK_LT(worst, 1e-6f);
    // Same reconstruction, but the pre-scaled float waveform differs by its
    // last bits, so out matches to float rounding rather than bit-for-bit.
    BOOST_CHECK_LT(worst_out, 2e-6f);

    // A perturbed public reference is a template mismatch, never a silent
    // rescale of a different waveform.
    std::vector<cf> bad_ref = ref;
    bad_ref[80] += 0.25f;
    std::vector<cf> outb(n, cf(0.f, 0.f)), ihatb(n, cf(0.f, 0.f));
    sic::DiscreteCirResult rb;
    BOOST_CHECK(!sic::apply_discrete_cir(y.data(), n, pkt.wave.data(),
                                         pkt.wave.size(), bad_ref.data(),
                                         bad_ref.size(), m, 0, int64_t(n),
                                         scratch, outb.data(), ihatb.data(), &rb));
    BOOST_CHECK(!rb.committed);
    BOOST_CHECK_EQUAL(std::string(rb.reason), "template_scale_mismatch");
    BOOST_CHECK_EQUAL(std::memcmp(y.data(), outb.data(), n * sizeof(cf)), 0);
    BOOST_CHECK_EQUAL(count_nonzero(ihatb, 0, 0), 0);
}

BOOST_AUTO_TEST_CASE(discrete_transaction_validates_fcs_phr_profile_and_length)
{
    TransactionRig rig(32, 19);
    auto r = rig.cancel();
    BOOST_REQUIRE_MESSAGE(r.committed, r.reason);
    BOOST_CHECK_EQUAL(std::string(r.reason), "committed");
    BOOST_CHECK_CLOSE(r.waveform_scale, double(kRefScale), 1e-6);
    BOOST_CHECK_EQUAL(r.packet_samples, rig.pkt.wave.size());
    BOOST_CHECK_GT(r.validate_us, 0.0);
    BOOST_CHECK_GT(r.tx_us, 0.0);
    BOOST_CHECK_GE(r.apply_total_us, r.channel_us + r.subtract_us);

    // The committed transaction must actually cancel the matching observation
    // inside its ROI and leave the floor untouched elsewhere.
    float worst_residual = 0.f;
    size_t out_bad = 0;
    for (int64_t j = 0; j < int64_t(rig.n); ++j) {
        const size_t s = static_cast<size_t>(j);
        if (j < r.sub_lo || j >= r.sub_hi) {
            if (!(rig.out[s] == rig.y[s]))
                ++out_bad;
        } else if (j >= rig.q + 3 &&
                   j < rig.q + 3 + int64_t(rig.pkt.wave.size())) {
            worst_residual =
                std::max(worst_residual, std::abs(rig.out[s] - cf(0.002f, -0.001f)));
        }
    }
    BOOST_CHECK_EQUAL(out_bad, 0);
    // Measured ~3.4e-8: the observation term and the reconstruction differ only
    // by float rounding of the same product, so 1e-6 is a ~30x margin.
    BOOST_CHECK_LT(worst_residual, 1e-6f);
    BOOST_CHECK_EQUAL(count_nonzero(rig.ihat, r.support_lo, r.support_hi), 0);

    // The wrapper's remodulation capacity is the real packet buffer length; a
    // stale planning bound must not reject a packet that fits.
    {
        const size_t saved = rig.scratch.max_packet_samples;
        rig.scratch.max_packet_samples = 1024;
        auto again = rig.cancel();
        BOOST_REQUIRE_MESSAGE(again.committed, again.reason);
        BOOST_CHECK_EQUAL(again.packet_samples, rig.pkt.wave.size());
        rig.scratch.max_packet_samples = saved;
    }

    // FCS flag.
    BOOST_CHECK_EQUAL(std::string(rig.cancel(false).reason), "fcs_failed");
    rig.expect_bypass();

    // PHR bit flip and non-binary bit.
    {
        auto phr = rig.pkt.phr;
        phr[9] ^= 1;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, nullptr, &phr).reason),
                          "unsupported_phr");
        rig.expect_bypass();
        phr[9] ^= 1;
        phr[9] = 2;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, nullptr, &phr).reason),
                          "unsupported_phr");
        rig.expect_bypass();
    }

    // Unsupported profiles (same gate as the frozen full-packet core).
    {
        mod::HrpModConfig p = rig.pkt.profile;
        p.sfd_mode = "ieee";
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, &p).reason),
                          "unsupported_profile");
        rig.expect_bypass();
        p = rig.pkt.profile;
        p.sync_repetitions = 16;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, &p).reason),
                          "unsupported_profile");
        rig.expect_bypass();
        p = rig.pkt.profile;
        p.code_index = 9;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, &p).reason),
                          "unsupported_profile");
        rig.expect_bypass();
        p = rig.pkt.profile;
        p.insert_sts = true;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, &p).reason),
                          "unsupported_profile");
        rig.expect_bypass();
        p = rig.pkt.profile;
        p.ranging = true;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, &p).reason),
                          "unsupported_profile");
        rig.expect_bypass();
        p = rig.pkt.profile;
        p.pulse.shape = mod::PulseShape::Gaussian;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, &p).reason),
                          "unsupported_profile");
        rig.expect_bypass();
    }

    // Wire length outside [2, 127].
    {
        std::vector<uint8_t> one(1, 0);
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, nullptr, nullptr, &one).reason),
                          "invalid_frame_length");
        rig.expect_bypass();
        std::vector<uint8_t> big(128, 0);
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, nullptr, nullptr, &big).reason),
                          "invalid_frame_length");
        rig.expect_bypass();
        // Correct length but a zeroed PHR: the encoder-derived PHR cannot match.
        std::vector<uint8_t> zero_phr(19, 0);
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, nullptr, &zero_phr).reason),
                          "unsupported_phr");
        rig.expect_bypass();
    }

    // CRC mismatch inside a valid length.
    {
        auto psdu = rig.pkt.psdu;
        psdu[0] ^= 1;
        BOOST_CHECK_EQUAL(std::string(rig.cancel(true, nullptr, nullptr, &psdu).reason),
                          "fcs_failed");
        rig.expect_bypass();
    }

    // Model gate: anchor, active set and finiteness.
    {
        rig.model.phase_reference_sample = rig.model.q + 1;
        BOOST_CHECK_EQUAL(std::string(rig.cancel().reason), "invalid_model");
        rig.expect_bypass();
        rig.model.phase_reference_sample = rig.model.q;

        rig.model.active_count = 0;
        BOOST_CHECK_EQUAL(std::string(rig.cancel().reason), "invalid_model");
        rig.expect_bypass();
        rig.model.active_count = 1;

        const int16_t saved_lag = rig.model.active_lags[0];
        rig.model.active_lags[0] = static_cast<int16_t>(rig.model.lag_lo + 1);
        BOOST_CHECK_EQUAL(std::string(rig.cancel().reason), "invalid_model");
        rig.expect_bypass();
        rig.model.active_lags[0] = saved_lag;

        rig.model.taps[0] = cf(std::numeric_limits<float>::quiet_NaN(), 0.f);
        BOOST_CHECK_EQUAL(std::string(rig.cancel().reason), "invalid_model");
        rig.expect_bypass();
        rig.model.taps[0] = cf(0.5f, 0.25f);

        rig.model.valid = false;
        BOOST_CHECK_EQUAL(std::string(rig.cancel().reason), "invalid_model");
        rig.expect_bypass();
        rig.model.valid = true;
    }

    // ROI and template gates.
    {
        BOOST_CHECK_EQUAL(std::string(sic::fullpacket_cancel_discrete(
                                          rig.y.data(), rig.n, rig.pkt.phr, rig.pkt.psdu,
                                          true, rig.pkt.profile, rig.model, rig.ref.data(),
                                          rig.ref.size(), 10, 10, rig.scratch,
                                          rig.out.data(), rig.ihat.data())
                                          .reason),
                          "invalid_roi");
        rig.expect_bypass();
        BOOST_CHECK_EQUAL(std::string(sic::fullpacket_cancel_discrete(
                                          rig.y.data(), rig.n, rig.pkt.phr, rig.pkt.psdu,
                                          true, rig.pkt.profile, rig.model, rig.ref.data(),
                                          rig.ref.size() - 1, 0, int64_t(rig.n),
                                          rig.scratch, rig.out.data(), rig.ihat.data())
                                          .reason),
                          "invalid_template");
        rig.expect_bypass();
    }
}

BOOST_AUTO_TEST_CASE(discrete_transaction_length_edges_16_and_127)
{
    // Shortest and longest accepted wire PSDU.  The FCS is appended by the QA
    // and the PHR derives from the same length, so a 16 B PHR must not validate
    // the 127 B frame.
    Packet small(16, 23);
    TransactionRig short_rig(16, 23);
    auto r16 = short_rig.cancel();
    BOOST_REQUIRE_MESSAGE(r16.committed, r16.reason);
    BOOST_CHECK_EQUAL(r16.packet_samples, short_rig.pkt.wave.size());

    TransactionRig long_rig(127, 29);
    auto r127 = long_rig.cancel();
    BOOST_REQUIRE_MESSAGE(r127.committed, r127.reason);
    BOOST_CHECK_EQUAL(r127.packet_samples, long_rig.pkt.wave.size());
    BOOST_CHECK_GT(r127.packet_samples, r16.packet_samples);

    // Same wire length, but the PHR of the other frame length.
    auto cross = long_rig.cancel(true, nullptr, &small.phr);
    BOOST_CHECK(!cross.committed);
    BOOST_CHECK_EQUAL(std::string(cross.reason), "unsupported_phr");
    long_rig.expect_bypass();
}

BOOST_AUTO_TEST_CASE(discrete_numeric_safety_and_scratch_capacity_bypass)
{
    ApplyRig rig(4096, 1024);
    sic::DiscreteCirResult r;

    // Sanity: the unmodified synthetic rig commits, so the bypass checks below
    // are not vacuous.
    BOOST_REQUIRE_MESSAGE(rig.run(r), r.reason);
    BOOST_CHECK_EQUAL(r.reason, "committed");

    // max_packet_samples is the estimator's planning bound only: a zero or
    // stale value must not reject a waveform that fits the real buffers.
    {
        const std::vector<cf> reference_out = rig.out;
        const size_t saved = rig.scratch.max_packet_samples;
        rig.scratch.max_packet_samples = 0;
        BOOST_REQUIRE_MESSAGE(rig.run(r), r.reason);
        BOOST_CHECK_EQUAL(std::string(r.reason), "committed");
        BOOST_CHECK_EQUAL(std::memcmp(rig.out.data(), reference_out.data(),
                                      rig.n * sizeof(cf)),
                          0);
        rig.scratch.max_packet_samples = 1024; // far below the packet length
        BOOST_REQUIRE_MESSAGE(rig.run(r), r.reason);
        BOOST_CHECK_EQUAL(std::memcmp(rig.out.data(), reference_out.data(),
                                      rig.n * sizeof(cf)),
                          0);
        rig.scratch.max_packet_samples = saved;
    }
    // Committed reference image; the failure cases below intentionally reset
    // out to y, so keep a copy instead of re-reading rig.out.
    const std::vector<cf> committed_out = rig.out;

    // NaN in the observation.
    {
        std::vector<cf> saved = rig.y;
        rig.y[3] = cf(std::numeric_limits<float>::quiet_NaN(), 0.f);
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "nonfinite_input");
        rig.expect_bypass();
        rig.y = saved;
    }

    // Non-finite waveform outside the reference prefix (the scale pass only
    // sees the first SYNC, so this exercises the full pre-scale check).
    {
        std::vector<cf> saved = rig.wave;
        rig.wave.back() = cf(std::numeric_limits<float>::infinity(), 0.f);
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "nonfinite_input");
        rig.expect_bypass();
        rig.wave = saved;
    }

    // Non-finite taps are a model failure, not a silent truncation.
    {
        rig.model.taps[2] = cf(std::numeric_limits<float>::quiet_NaN(), 0.f);
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "invalid_model");
        rig.expect_bypass();
        rig.model.taps[2] = cf(0.25f, 0.f);
    }

    // Finite taps whose accumulated reconstruction overflows float.
    {
        const float huge = std::numeric_limits<float>::max();
        rig.model.taps[0] = cf(huge, 0.f);
        rig.model.taps[1] = cf(huge, 0.f);
        rig.model.taps[2] = cf(huge, 0.f);
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "nonfinite_reconstruction");
        rig.expect_bypass();
        rig.model.taps[0] = cf(1.f, 0.f);
        rig.model.taps[1] = cf(0.5f, 0.f);
        rig.model.taps[2] = cf(0.25f, 0.f);
    }

    // Optional ihat: the same reconstruction must still be produced, out must
    // still equal y outside the subtraction interval, and a failure must still
    // restore out exactly.
    {
        std::vector<cf> out_only(rig.n, cf(0.f, 0.f));
        sic::DiscreteCirResult rn;
        BOOST_REQUIRE(sic::apply_discrete_cir(rig.y.data(), rig.n, rig.wave.data(),
                                              rig.wave.size(), rig.ref.data(),
                                              rig.ref.size(), rig.model, 0,
                                              int64_t(rig.n), rig.scratch,
                                              out_only.data(), nullptr, &rn));
        BOOST_CHECK(rn.committed);
        BOOST_CHECK_EQUAL(std::memcmp(out_only.data(), committed_out.data(),
                                      rig.n * sizeof(cf)),
                          0);
        // With ihat == nullptr the trial reconstruction lives in `out`, so the
        // flanks of the support outside the ROI must be restored to y.
        size_t flank_bad = 0;
        for (int64_t j = rn.support_lo; j < rn.sub_lo; ++j)
            if (!(out_only[static_cast<size_t>(j)] == rig.y[static_cast<size_t>(j)]))
                ++flank_bad;
        for (int64_t j = rn.sub_hi; j < rn.support_hi; ++j)
            if (!(out_only[static_cast<size_t>(j)] == rig.y[static_cast<size_t>(j)]))
                ++flank_bad;
        BOOST_CHECK_EQUAL(flank_bad, 0);

        rig.model.taps[0] = cf(std::numeric_limits<float>::max(), 0.f);
        rig.model.taps[1] = cf(std::numeric_limits<float>::max(), 0.f);
        rig.model.taps[2] = cf(std::numeric_limits<float>::max(), 0.f);
        std::vector<cf> out_fail(rig.n, cf(0.f, 0.f));
        BOOST_CHECK(!sic::apply_discrete_cir(rig.y.data(), rig.n, rig.wave.data(),
                                             rig.wave.size(), rig.ref.data(),
                                             rig.ref.size(), rig.model, 0,
                                             int64_t(rig.n), rig.scratch,
                                             out_fail.data(), nullptr, &rn));
        BOOST_CHECK_EQUAL(std::string(rn.reason), "nonfinite_reconstruction");
        BOOST_CHECK_EQUAL(std::memcmp(rig.y.data(), out_fail.data(),
                                      rig.n * sizeof(cf)),
                          0);
        rig.model.taps[0] = cf(1.f, 0.f);
        rig.model.taps[1] = cf(0.5f, 0.f);
        rig.model.taps[2] = cf(0.25f, 0.f);
    }

    // Scratch capacity contract.
    {
        rig.scratch.prepared = false;
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "scratch_unprepared");
        rig.expect_bypass();
        rig.scratch.prepared = true;

        const size_t saved_window = rig.scratch.max_window_samples;
        rig.scratch.max_window_samples = rig.n - 1;
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "scratch_capacity");
        rig.expect_bypass();
        rig.scratch.max_window_samples = saved_window;

        rig.scratch.scaled_waveform.resize(1024);
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "scratch_capacity");
        rig.expect_bypass();
        rig.scratch.scaled_waveform.resize(rig.scratch.max_packet_samples);

        // The tap-major FIR accumulates in scratch.accumulator, which must
        // cover the whole window.
        rig.scratch.accumulator.resize(rig.n - 1);
        BOOST_CHECK(!rig.run(r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "scratch_capacity");
        rig.expect_bypass();
        rig.scratch.accumulator.assign(rig.n, cf(0.f, 0.f));
    }

    // Template gate still holds for the synthetic rig.
    {
        std::vector<cf> bad = rig.ref;
        bad.back() += 0.5f;
        sic::DiscreteCirResult rt;
        BOOST_CHECK(!sic::apply_discrete_cir(rig.y.data(), rig.n, rig.wave.data(),
                                             rig.wave.size(), bad.data(), bad.size(),
                                             rig.model, 0, int64_t(rig.n),
                                             rig.scratch, rig.out.data(),
                                             rig.ihat.data(), &rt));
        BOOST_CHECK_EQUAL(std::string(rt.reason), "template_scale_mismatch");
        rig.expect_bypass();
    }
}

BOOST_AUTO_TEST_CASE(b1_avoid_out_double_write_and_ihat_semantics)
{
    ApplyRig rig(4096, 1024);
    const int64_t roi_lo = 1500;
    const int64_t roi_hi = 3000;
    sic::DiscreteCirResult r;
    BOOST_REQUIRE(sic::apply_discrete_cir(
        rig.y.data(), rig.n, rig.wave.data(), rig.wave.size(),
        rig.ref.data(), rig.ref.size(), rig.model, roi_lo, roi_hi,
        rig.scratch, rig.out.data(), rig.ihat.data(), &r));
    BOOST_CHECK(r.committed);
    BOOST_CHECK_EQUAL(r.sub_lo, roi_lo);
    BOOST_CHECK_EQUAL(r.sub_hi, roi_hi);

    // 1. Outside [sub_lo, sub_hi), out must be byte-identical to y.
    BOOST_CHECK_EQUAL(std::memcmp(rig.out.data(), rig.y.data(), roi_lo * sizeof(cf)), 0);
    BOOST_CHECK_EQUAL(std::memcmp(rig.out.data() + roi_hi, rig.y.data() + roi_hi,
                                  (rig.n - roi_hi) * sizeof(cf)), 0);

    // 2. Inside [sub_lo, sub_hi), out == y - ihat.
    for (int64_t j = roi_lo; j < roi_hi; ++j) {
        BOOST_CHECK_EQUAL(rig.out[j], rig.y[j] - rig.ihat[j]);
    }

    // 3. ihat must be strictly zero outside [support_lo, support_hi).
    for (int64_t j = 0; j < r.support_lo; ++j) {
        BOOST_CHECK_EQUAL(rig.ihat[j], cf(0.f, 0.f));
    }
    for (int64_t j = r.support_hi; j < int64_t(rig.n); ++j) {
        BOOST_CHECK_EQUAL(rig.ihat[j], cf(0.f, 0.f));
    }

    // 4. In-place y == out test:
    std::vector<cf> in_place = rig.y;
    sic::DiscreteCirResult r_ip;
    BOOST_REQUIRE(sic::apply_discrete_cir(
        in_place.data(), rig.n, rig.wave.data(), rig.wave.size(),
        rig.ref.data(), rig.ref.size(), rig.model, roi_lo, roi_hi,
        rig.scratch, in_place.data(), nullptr, &r_ip));
    BOOST_CHECK(r_ip.committed);
    BOOST_CHECK_EQUAL(std::memcmp(in_place.data(), rig.out.data(), rig.n * sizeof(cf)), 0);

    // 5. Failure restore test:
    std::vector<cf> out_fail(rig.n, cf(999.f, -999.f));
    std::vector<cf> ihat_fail(rig.n, cf(999.f, -999.f));
    rig.model.cfo_hz = 998.4e6; // beyond Nyquist -> invalid_model
    BOOST_CHECK(!sic::apply_discrete_cir(
        rig.y.data(), rig.n, rig.wave.data(), rig.wave.size(),
        rig.ref.data(), rig.ref.size(), rig.model, roi_lo, roi_hi,
        rig.scratch, out_fail.data(), ihat_fail.data(), &r));
    BOOST_CHECK(!r.committed);
    BOOST_CHECK_EQUAL(std::memcmp(out_fail.data(), rig.y.data(), rig.n * sizeof(cf)), 0);
    BOOST_CHECK_EQUAL(count_nonzero(ihat_fail, 0, 0), 0);
}

BOOST_AUTO_TEST_CASE(sparse_pulse_apply_matches_exact_cfo_reference)
{
    SparseRig rig;
    const int64_t q = rig.q;
    const size_t n = rig.n;
    const std::vector<mod::PulseEvent>& events = rig.events;

    for (double cfo : { 0.0, -50000.0, 50000.0, 100000.0 }) {
        sic::DiscreteCirModel m = rig.model(q, -2, 10);
        m.cfo_hz = cfo;
        set_active(m, { { -2, cf(0.9f, -0.2f) }, { 0, cf(0.4f, 0.25f) },
                        { 7, cf(-0.15f, 0.1f) } });
        const int64_t lo = q + m.lag_lo;
        const int64_t hi =
            lo + int64_t(m.tap_count) + int64_t(rig.pl.packet_samples);

        std::vector<cf> y(n, cf(0.01f, -0.02f)), out(n), ihat(n);
        sic::DiscreteCirResult r;
        BOOST_REQUIRE_MESSAGE(
            sic::apply_discrete_cir_pulses(
                y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(n), rig.scratch,
                out.data(), ihat.data(), &r),
            r.reason);
        BOOST_CHECK(r.committed);
        BOOST_CHECK_EQUAL(std::string(r.reason), "committed");
        BOOST_CHECK_EQUAL(r.cfo_us, 0.0);
        BOOST_CHECK_EQUAL(r.pulse_event_count, rig.pl.count);
        BOOST_CHECK_EQUAL(r.pulse_packet_samples, rig.pl.packet_samples);

        // Independent double oracle of the EXACT discrete model.  Both the
        // event phase exp(j*w*m_i) and the per-lag phase exp(j*w*l) must be
        // present; the old approximate form only had exp(j*w*m_i).
        std::vector<cd> want(n, cd{});
        const double omega = 2.0 * kPi * cfo / kWorkRateHz;
        for (const auto& e : events) {
            const double mph = omega * double(int64_t(e.sample_offset));
            for (size_t d = 0; d < m.tap_count; ++d) {
                const int64_t lag = m.lag_lo + int64_t(d);
                const int64_t j = q + int64_t(e.sample_offset) + lag;
                if (j >= 0 && j < int64_t(n))
                    want[static_cast<size_t>(j)] +=
                        double(e.sign) *
                        std::polar(1.0, mph + omega * double(lag)) *
                        cd(m.taps[d]);
            }
        }
        float worst = 0.f;
        for (size_t j = 0; j < n; ++j)
            worst = std::max(worst, std::abs(ihat[j] - cf(want[j])));
        BOOST_CHECK_LT(worst, 2e-5f);

        for (int64_t j = 0; j < int64_t(n); ++j) {
            const size_t u = static_cast<size_t>(j);
            if (j < lo || j >= hi) {
                BOOST_CHECK_EQUAL(ihat[u], cf{});
                BOOST_CHECK_EQUAL(out[u], y[u]);
            } else {
                BOOST_CHECK_LT(std::abs(out[u] - (y[u] - ihat[u])), 1e-6f);
            }
        }

        // The L1D tile path must reproduce the full-buffer ihat path
        // BIT-FOR-BIT (same events, same kernels, same addition order).
        std::vector<cf> out_only(n);
        sic::DiscreteCirResult rn;
        BOOST_REQUIRE_MESSAGE(
            sic::apply_discrete_cir_pulses(
                y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(n), rig.scratch,
                out_only.data(), nullptr, &rn),
            rn.reason);
        BOOST_CHECK(rn.committed);
        BOOST_CHECK_EQUAL(std::memcmp(out_only.data(), out.data(),
                                      n * sizeof(cf)),
                          0);

        if (cfo == 100000.0) {
            // 100 kHz high-CFO regression against the full-waveform path,
            // which rotates exp(j*w*(j-q)) exactly.  The sparse model must
            // match to float rounding only: relative L2 residual < 1e-4.
            std::vector<cf> full_out(n), full_ihat(n);
            sic::DiscreteCirResult rf;
            BOOST_REQUIRE_MESSAGE(
                sic::apply_discrete_cir(
                    y.data(), n, rig.unit_wave.data(), rig.unit_wave.size(),
                    rig.scratch.impulse_sync.data(),
                    rig.scratch.impulse_sync.size(), m, 0, int64_t(n),
                    rig.scratch, full_out.data(), full_ihat.data(), &rf),
                rf.reason);
            BOOST_CHECK(rf.committed);
            double num = 0.0, den = 0.0;
            for (int64_t j = lo; j < hi; ++j) {
                const size_t u = static_cast<size_t>(j);
                num += std::norm(cd(ihat[u]) - cd(full_ihat[u]));
                den += std::norm(cd(full_ihat[u]));
            }
            BOOST_REQUIRE_GT(den, 0.0);
            BOOST_CHECK_LT(std::sqrt(num / den), 1e-4);
        }
    }

    // Representation failures must bypass transactionally.
    {
        sic::DiscreteCirModel m = make_model(q, -2, 10);
        m.impulse_template = true;
        set_active(m, { { 0, cf(0.5f, 0.0f) } });
        std::vector<cf> y(n, cf(0.02f, 0.01f)), out(n), ihat(n);
        auto check_fail = [&](const mod::PulseList& bad, const char* reason) {
            sic::DiscreteCirResult r;
            BOOST_CHECK(!sic::apply_discrete_cir_pulses(
                y.data(), n, bad, rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(n),
                rig.scratch, out.data(), ihat.data(), &r));
            BOOST_CHECK_EQUAL(std::string(r.reason), reason);
            BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)),
                              0);
            BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
        };
        auto duplicate = events;
        if (duplicate.size() > 1)
            duplicate[1].sample_offset = duplicate[0].sample_offset;
        mod::PulseList bad{ duplicate.data(), duplicate.size(),
                            rig.pl.packet_samples };
        check_fail(bad, "invalid_pulse_list");

        auto wrong_sync = events;
        wrong_sync[0].sign = static_cast<int8_t>(-wrong_sync[0].sign);
        mod::PulseList bad_sync{ wrong_sync.data(), wrong_sync.size(),
                                 rig.pl.packet_samples };
        check_fail(bad_sync, "pulse_template_mismatch");

        auto far = m;
        far.q = int64_t(n) - 1;
        far.phase_reference_sample = far.q;
        sic::DiscreteCirResult rf;
        BOOST_CHECK(!sic::apply_discrete_cir_pulses(
            y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), far, 0, int64_t(n), rig.scratch,
            out.data(), ihat.data(), &rf));
        BOOST_CHECK_EQUAL(std::string(rf.reason), "context_truncated");
        BOOST_CHECK_EQUAL(std::memcmp(y.data(), out.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
    }
}

BOOST_AUTO_TEST_CASE(b2_b3_sparse_pulse_avx2_fma_20tap_and_fallback_parity)
{
    SparseRig rig;
    const int64_t q = rig.q;
    const size_t n = rig.n;
    std::vector<cf> y(n, cf(0.01f, -0.02f));

    // 20 taps keeps the fully unrolled register-resident fast path; every
    // other width must use the generic AVX2+FMA kernels instead of silently
    // falling back to the noinline scalar fused_axpy.
    for (size_t tap_count : { 1, 7, 16, 19, 20, 24, 32, 48 }) {
        for (double cfo : { 0.0, 35000.0 }) {
            sic::DiscreteCirModel m = rig.model(q, -1, tap_count);
            m.cfo_hz = cfo;
            for (size_t t = 0; t < tap_count; ++t)
                m.taps[t] = cf(0.1f * float(t + 1), -0.05f * float(t));
            m.active_count = tap_count;
            for (size_t t = 0; t < tap_count; ++t)
                m.active_lags[t] =
                    static_cast<int16_t>(m.lag_lo + int64_t(t));

            // Run with B3 SIMD enabled.
            sic::set_b3_simd_disabled(0);
            std::vector<cf> out_opt(n), ihat_opt(n);
            sic::DiscreteCirResult r_opt;
            BOOST_REQUIRE(sic::apply_discrete_cir_pulses(
                y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(n), rig.scratch,
                out_opt.data(), ihat_opt.data(), &r_opt));
            BOOST_CHECK(r_opt.committed);

            // Run with B3 SIMD disabled (scalar reference).
            sic::set_b3_simd_disabled(1);
            std::vector<cf> out_fallback(n), ihat_fallback(n);
            sic::DiscreteCirResult r_fallback;
            BOOST_REQUIRE(sic::apply_discrete_cir_pulses(
                y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(n), rig.scratch,
                out_fallback.data(), ihat_fallback.data(), &r_fallback));
            BOOST_CHECK(r_fallback.committed);
            sic::set_b3_simd_disabled(-1);

            if (cfo == 0.0) {
                // Zero CFO must be bit-for-bit identical between the SIMD
                // kernels and the scalar reference.
                BOOST_CHECK_EQUAL(std::memcmp(out_opt.data(),
                                              out_fallback.data(),
                                              n * sizeof(cf)),
                                  0);
                BOOST_CHECK_EQUAL(std::memcmp(ihat_opt.data(),
                                              ihat_fallback.data(),
                                              n * sizeof(cf)),
                                  0);
            } else {
                // Non-zero CFO uses FMA vs non-FMA; the accepted error scales
                // with the tap count (the 20-tap bar stays at 1e-5).
                float max_err = 0.f;
                for (size_t j = 0; j < n; ++j)
                    max_err =
                        std::max(max_err, std::abs(out_opt[j] - out_fallback[j]));
                const float limit =
                    1e-5f * std::max(1.f, float(tap_count) / 20.f);
                BOOST_CHECK_LT(max_err, limit);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(tile_path_matches_full_ihat_path_bit_exact)
{
    SparseRig rig;
    std::vector<cf> y(rig.n, cf(0.01f, -0.02f));

    for (size_t tap_count : { size_t(7), size_t(20) }) {
        for (double cfo : { 0.0, 35000.0 }) {
            sic::DiscreteCirModel m = rig.model(rig.q, -1, tap_count);
            m.cfo_hz = cfo;
            for (size_t t = 0; t < tap_count; ++t)
                m.taps[t] = cf(0.1f * float(t + 1), -0.05f * float(t));
            m.active_count = tap_count;
            for (size_t t = 0; t < tap_count; ++t)
                m.active_lags[t] =
                    static_cast<int16_t>(m.lag_lo + int64_t(t));

            const int64_t roi_lo = 1000;
            const int64_t roi_hi = int64_t(rig.n) - 1000;
            sic::DiscreteCirApplyOptions opts;
            opts.compute_roi_energy = true;

            // Reference: full-materialized ihat path.
            std::vector<cf> out_ref(rig.n), ihat_ref(rig.n);
            sic::DiscreteCirResult r_ref;
            BOOST_REQUIRE_MESSAGE(
                sic::apply_discrete_cir_pulses_ex(
                    y.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
                    rig.scratch.impulse_sync.size(), m, roi_lo, roi_hi, opts,
                    rig.scratch, out_ref.data(), ihat_ref.data(), &r_ref),
                r_ref.reason);
            BOOST_CHECK(r_ref.committed);
            BOOST_CHECK(r_ref.roi_energy_valid);

            // L1D tile path, out != y.
            std::vector<cf> out_tile(rig.n, cf(0.f, 0.f));
            sic::DiscreteCirResult r_tile;
            BOOST_REQUIRE_MESSAGE(
                sic::apply_discrete_cir_pulses_ex(
                    y.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
                    rig.scratch.impulse_sync.size(), m, roi_lo, roi_hi, opts,
                    rig.scratch, out_tile.data(), nullptr, &r_tile),
                r_tile.reason);
            BOOST_CHECK(r_tile.committed);
            BOOST_CHECK(r_tile.roi_energy_valid);
            BOOST_CHECK_EQUAL(std::memcmp(out_ref.data(), out_tile.data(),
                                          rig.n * sizeof(cf)),
                              0);

            // L1D tile path, out == y in place.
            std::vector<cf> y_inplace = y;
            sic::DiscreteCirResult r_inplace;
            BOOST_REQUIRE_MESSAGE(
                sic::apply_discrete_cir_pulses_ex(
                    y_inplace.data(), rig.n, rig.pl,
                    rig.scratch.impulse_sync.data(),
                    rig.scratch.impulse_sync.size(), m, roi_lo, roi_hi, opts,
                    rig.scratch, y_inplace.data(), nullptr, &r_inplace),
                r_inplace.reason);
            BOOST_CHECK(r_inplace.committed);
            BOOST_CHECK_EQUAL(std::memcmp(y_inplace.data(), out_tile.data(),
                                          rig.n * sizeof(cf)),
                              0);

            // ROI energies against an independent double accumulation.  The
            // contract (matching compute_roi_energies) is float squares
            // widened to double, accumulated over support ∩ ROI.
            BOOST_REQUIRE_EQUAL(r_ref.sub_lo, r_tile.sub_lo);
            BOOST_REQUIRE_EQUAL(r_ref.sub_hi, r_tile.sub_hi);
            double e_y = 0.0, e_out = 0.0;
            for (int64_t j = r_tile.sub_lo; j < r_tile.sub_hi; ++j) {
                const cf yv = y[static_cast<size_t>(j)];
                const cf ov = out_tile[static_cast<size_t>(j)];
                e_y += double(yv.real() * yv.real()) +
                       double(yv.imag() * yv.imag());
                e_out += double(ov.real() * ov.real()) +
                         double(ov.imag() * ov.imag());
            }
            auto rel = [](double a, double b) {
                return std::fabs(a - b) / std::max(1e-300, std::fabs(b));
            };
            BOOST_CHECK_LT(rel(r_tile.roi_e_y, e_y), 1e-12);
            BOOST_CHECK_LT(rel(r_tile.roi_e_out, e_out), 1e-12);
            BOOST_CHECK_LT(rel(r_ref.roi_e_y, e_y), 1e-12);
            BOOST_CHECK_LT(rel(r_ref.roi_e_out, e_out), 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE(partial_window_clip_left_and_right_edges)
{
    SparseRig rig;
    // q = 2 with lag_lo = -4 puts the support start at -2 (left clip) while a
    // window ending `n = packet_samples - 1000` cuts the packet tail (right
    // clip): hi = q + lag_lo + tap_count + packet_samples > n.
    const int64_t q = 2;
    const size_t n = rig.pl.packet_samples - 1000;
    BOOST_REQUIRE_GT(n, size_t(kSyncSamples));
    const size_t tap_count = 7;
    const double cfo = 20000.0;

    sic::DiscreteCirModel m = rig.model(q, -4, tap_count);
    m.cfo_hz = cfo;
    for (size_t t = 0; t < tap_count; ++t)
        m.taps[t] = cf(0.1f * float(t + 1), -0.05f * float(t));
    m.active_count = tap_count;
    for (size_t t = 0; t < tap_count; ++t)
        m.active_lags[t] = static_cast<int16_t>(m.lag_lo + int64_t(t));
    const int64_t lo = q + m.lag_lo;
    BOOST_REQUIRE_LT(lo, 0);
    BOOST_REQUIRE_GT(lo + int64_t(tap_count) + int64_t(rig.pl.packet_samples),
                     int64_t(n));

    std::vector<cf> y(n, cf(0.003f, -0.004f));
    const cf sentinel(777.f, -777.f);
    const size_t guard = 8;

    // Default: an incomplete observed packet is a hard bypass.
    {
        std::vector<cf> out(n + guard, sentinel), ihat(n, sentinel);
        std::vector<cf> out_guard(out);
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir_pulses(
            y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(n), rig.scratch,
            out.data(), ihat.data(), &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "context_truncated");
        BOOST_CHECK_EQUAL(std::memcmp(out.data(), y.data(), n * sizeof(cf)), 0);
        BOOST_CHECK_EQUAL(count_nonzero(ihat, 0, 0), 0);
        BOOST_CHECK_EQUAL(std::memcmp(out.data() + n, out_guard.data() + n,
                                      guard * sizeof(cf)),
                          0);
    }

    sic::DiscreteCirApplyOptions opts;
    opts.allow_partial_clip = true;
    opts.compute_roi_energy = true;

    // Opt-in: clipped ihat path.
    std::vector<cf> out_ref(n + guard, sentinel), ihat_ref(n, sentinel);
    std::vector<cf> out_ref_guard(out_ref);
    sic::DiscreteCirResult r_ref;
    BOOST_REQUIRE_MESSAGE(
        sic::apply_discrete_cir_pulses_ex(
            y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(n), opts,
            rig.scratch, out_ref.data(), ihat_ref.data(), &r_ref),
        r_ref.reason);
    BOOST_CHECK(r_ref.committed);
    BOOST_CHECK_EQUAL(r_ref.support_lo, lo);
    BOOST_CHECK_EQUAL(r_ref.sub_lo, 0);
    BOOST_CHECK_EQUAL(r_ref.sub_hi, int64_t(n));
    BOOST_CHECK(r_ref.roi_energy_valid);
    BOOST_CHECK_EQUAL(std::memcmp(out_ref.data() + n, out_ref_guard.data() + n,
                                  guard * sizeof(cf)),
                      0);

    // Opt-in: clipped L1D tile path (ihat == nullptr).
    std::vector<cf> out_tile(n + guard, sentinel);
    std::vector<cf> out_tile_guard(out_tile);
    sic::DiscreteCirResult r_tile;
    BOOST_REQUIRE_MESSAGE(
        sic::apply_discrete_cir_pulses_ex(
            y.data(), n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(n), opts,
            rig.scratch, out_tile.data(), nullptr, &r_tile),
        r_tile.reason);
    BOOST_CHECK(r_tile.committed);
    BOOST_CHECK(r_tile.roi_energy_valid);
    BOOST_CHECK_EQUAL(std::memcmp(out_tile.data() + n,
                                  out_tile_guard.data() + n,
                                  guard * sizeof(cf)),
                      0);

    // Independent double oracle: only tap fragments that land inside [0, n)
    // may contribute, and neither boundary may be written out of bounds.
    std::vector<cd> want(n, cd{});
    const double omega = 2.0 * kPi * cfo / kWorkRateHz;
    for (const auto& e : rig.events) {
        const double mph = omega * double(int64_t(e.sample_offset));
        for (size_t d = 0; d < tap_count; ++d) {
            const int64_t lag = m.lag_lo + int64_t(d);
            const int64_t j = q + int64_t(e.sample_offset) + lag;
            if (j < 0 || j >= int64_t(n))
                continue;
            want[static_cast<size_t>(j)] +=
                double(e.sign) * std::polar(1.0, mph + omega * double(lag)) *
                cd(m.taps[d]);
        }
    }
    float worst_ref = 0.f, worst_tile = 0.f;
    for (size_t j = 0; j < n; ++j) {
        const cf expected = cf(cd(y[j]) - want[j]);
        worst_ref = std::max(worst_ref, std::abs(out_ref[j] - expected));
        worst_tile = std::max(worst_tile, std::abs(out_tile[j] - expected));
    }
    BOOST_CHECK_LT(worst_ref, 3e-5f);
    BOOST_CHECK_LT(worst_tile, 3e-5f);

    // Fused energies cover exactly the support ∩ window under a double
    // accumulation of the same float squares the kernel forms.
    double e_y = 0.0, e_out = 0.0;
    for (int64_t j = r_tile.sub_lo; j < r_tile.sub_hi; ++j) {
        const cf yv = y[static_cast<size_t>(j)];
        const cf ov = out_tile[static_cast<size_t>(j)];
        e_y += double(yv.real() * yv.real()) + double(yv.imag() * yv.imag());
        e_out += double(ov.real() * ov.real()) + double(ov.imag() * ov.imag());
    }
    auto rel = [](double a, double b) {
        return std::fabs(a - b) / std::max(1e-300, std::fabs(b));
    };
    BOOST_CHECK_LT(rel(r_tile.roi_e_y, e_y), 1e-10);
    BOOST_CHECK_LT(rel(r_tile.roi_e_out, e_out), 1e-10);
}

BOOST_AUTO_TEST_CASE(copy_out_on_fail_false_skips_buffer_write)
{
    SparseRig rig;
    const size_t tap_count = 7;
    sic::DiscreteCirModel m = rig.model(rig.q, -1, tap_count);
    for (size_t t = 0; t < tap_count; ++t) {
        m.taps[t] = cf(0.1f * float(t + 1), -0.05f * float(t));
        m.active_lags[t] = static_cast<int16_t>(m.lag_lo + int64_t(t));
    }
    m.active_count = tap_count;
    // Support cannot fit at the very end of the window -> context_truncated.
    m.q = int64_t(rig.n) - 1;
    m.phase_reference_sample = m.q;

    std::vector<cf> y(rig.n, cf(0.02f, 0.01f));
    const cf sentinel(999.f, -999.f);

    {
        std::vector<cf> out(rig.n, sentinel);
        sic::DiscreteCirApplyOptions opts;
        opts.copy_out_on_fail = false;
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir_pulses_ex(
            y.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(rig.n), opts,
            rig.scratch, out.data(), nullptr, &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "context_truncated");
        // The opted-out failure must not touch the output buffer at all.
        bool untouched = true;
        for (size_t j = 0; j < rig.n; ++j)
            if (!(out[j] == sentinel))
                untouched = false;
        BOOST_CHECK(untouched);
    }

    {
        // Default copy_out_on_fail == true restores the observation.
        std::vector<cf> out(rig.n, sentinel);
        sic::DiscreteCirApplyOptions opts;
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir_pulses_ex(
            y.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(rig.n), opts,
            rig.scratch, out.data(), nullptr, &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "context_truncated");
        BOOST_CHECK_EQUAL(std::memcmp(out.data(), y.data(),
                                      rig.n * sizeof(cf)),
                          0);
    }

    // The transactional wrapper honours the same option on its own failures.
    {
        TransactionRig tr(32, 19);
        const std::vector<cf> before = tr.out;
        sic::DiscreteCirApplyOptions opts;
        opts.copy_out_on_fail = false;
        auto r = sic::fullpacket_cancel_discrete_ex(
            tr.y.data(), tr.n, tr.pkt.phr, tr.pkt.psdu, false, tr.pkt.profile,
            tr.model, tr.ref.data(), tr.ref.size(), 0, int64_t(tr.n), opts,
            tr.scratch, tr.out.data(), tr.ihat.data());
        BOOST_CHECK_EQUAL(std::string(r.reason), "fcs_failed");
        BOOST_CHECK_EQUAL(std::memcmp(tr.out.data(), before.data(),
                                      tr.n * sizeof(cf)),
                          0);
    }
}

BOOST_AUTO_TEST_CASE(inplace_overflow_injection_is_fail_closed)
{
    SparseRig rig;
    const size_t tap_count = 7;
    std::vector<cf> y(rig.n, cf(0.01f, -0.02f));

    // Sanity: the same in-place call with normal taps commits and modifies the
    // observation inside the ROI.
    {
        sic::DiscreteCirModel m = rig.model(rig.q, -1, tap_count);
        for (size_t t = 0; t < tap_count; ++t) {
            m.taps[t] = cf(0.1f * float(t + 1), -0.05f * float(t));
            m.active_lags[t] = static_cast<int16_t>(m.lag_lo + int64_t(t));
        }
        m.active_count = tap_count;
        std::vector<cf> y_inplace = y;
        sic::DiscreteCirResult r;
        BOOST_REQUIRE_MESSAGE(
            sic::apply_discrete_cir_pulses(
                y_inplace.data(), rig.n, rig.pl,
                rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(rig.n),
                rig.scratch, y_inplace.data(), nullptr, &r),
            r.reason);
        BOOST_CHECK(r.committed);
        BOOST_CHECK(std::memcmp(y_inplace.data(), y.data(),
                                rig.n * sizeof(cf)) != 0);
    }

    // Inject finite taps far above the safe in-place bound.  The two-phase
    // scratch.accumulator fallback must validate before it writes, so the
    // in-place window stays byte-identical on failure.
    sic::DiscreteCirModel m = rig.model(rig.q, -1, tap_count);
    const float huge = std::numeric_limits<float>::max();
    for (size_t t = 0; t < tap_count; ++t) {
        m.taps[t] = cf(huge, 0.f);
        m.active_lags[t] = static_cast<int16_t>(m.lag_lo + int64_t(t));
    }
    m.active_count = tap_count;

    {
        std::vector<cf> y_inplace = y;
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir_pulses(
            y_inplace.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(rig.n),
            rig.scratch, y_inplace.data(), nullptr, &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "nonfinite_reconstruction");
        BOOST_CHECK(!r.committed);
        BOOST_CHECK_EQUAL(std::memcmp(y_inplace.data(), y.data(),
                                      rig.n * sizeof(cf)),
                          0);
    }

    // Out-of-place failure also restores out == y (copy_out_on_fail default).
    {
        std::vector<cf> out(rig.n, cf(5.f, -5.f));
        sic::DiscreteCirResult r;
        BOOST_CHECK(!sic::apply_discrete_cir_pulses(
            y.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
            rig.scratch.impulse_sync.size(), m, 0, int64_t(rig.n),
            rig.scratch, out.data(), nullptr, &r));
        BOOST_CHECK_EQUAL(std::string(r.reason), "nonfinite_reconstruction");
        BOOST_CHECK_EQUAL(std::memcmp(out.data(), y.data(),
                                      rig.n * sizeof(cf)),
                          0);
    }
}
