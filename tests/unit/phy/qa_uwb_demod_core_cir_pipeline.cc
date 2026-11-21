// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 2/4 of the split realtime-demodulator core QA: R2-R4
// CIR/softchips/payload plus the SFD-anchored CFO-refit policy.
#include "demod_core_test_helpers.h"

// ---------------------------------------------------------------------------
// SFD-anchored CFO refit policy (single-pass rotation).
//
// `demodulate_one` fits the carrier offset twice and derotates the downstream
// span twice.  The second pass is the largest single cost inside `stage_sfd`.
// `Qm35825Profile::cfo_refit_mode` can replace it with one rotation that uses
// the Stage-2 estimate.  The contract that must hold:
//   * legacy (default) still runs both fits and both rotations,
//   * skip keeps the Stage-2 estimate and drops one rotation pass,
//   * the decode verdict (status / PSDU / FCS) is unchanged,
//   * auto falls back to the legacy path when the drift budget is exceeded.
// ---------------------------------------------------------------------------
namespace {

struct CfoRefitRun {
    DemodResult res;
    std::vector<uint8_t> psdu;
};

CfoRefitRun run_cfo_refit_case(const std::vector<gr_complex>& iq,
                               const std::vector<gr_complex>& tmpl,
                               int mode,
                               double skip_rad,
                               size_t min_peaks)
{
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    prof.cfo_refit_mode = mode;
    prof.cfo_refit_skip_rad = skip_rad;
    prof.cfo_refit_min_peaks = min_peaks;
    CfoRefitRun out;
    out.res = core::demodulate_one(iq.data(), iq.size(), prof, 1, 9984, 0, tmpl, scratch);
    out.psdu = out.res.payload.bytes;
    return out;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_demod_core_cfo_refit_skip_matches_legacy_decode)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    BOOST_REQUIRE(!tmpl.empty());

    // Legacy: both fits, both rotations.  The Stage-2 fit is published, the
    // refit replaces it, and the second rotation covers the whole window.
    const auto legacy = run_cfo_refit_case(iq, tmpl, 0, 0.10, 64);
    BOOST_REQUIRE(legacy.res.status == DemodStatus::Success);
    BOOST_CHECK(legacy.res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(legacy.res.cfo_refit_action, uint32_t(0));
    BOOST_CHECK_EQUAL(legacy.res.stage_cfo_rotation_passes, uint32_t(2));
    BOOST_CHECK_GT(legacy.res.cfo_refit_span_samples, uint64_t(0));
    BOOST_CHECK_GT(legacy.res.cfo_initial_peaks, uint32_t(0));
    BOOST_CHECK_GT(legacy.res.cfo_refit_peaks, uint32_t(0));
    BOOST_CHECK_EQUAL(legacy.res.cfo.cfo_hz, legacy.res.cfo_refit_hz);
    BOOST_CHECK(std::abs(legacy.res.cfo_initial_hz - legacy.res.cfo_refit_hz) < 1e4);

    // Forced skip: the Stage-2 estimate survives, one rotation pass is gone,
    // and the payload/FCS verdict is identical.
    const auto skipped = run_cfo_refit_case(iq, tmpl, 2, 0.10, 64);
    BOOST_REQUIRE(skipped.res.status == DemodStatus::Success);
    BOOST_CHECK(skipped.res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(skipped.res.cfo_refit_action, uint32_t(2));
    BOOST_CHECK_EQUAL(skipped.res.stage_cfo_rotation_passes, uint32_t(1));
    BOOST_CHECK_EQUAL(skipped.res.cfo.cfo_hz, skipped.res.cfo_initial_hz);
    BOOST_CHECK_EQUAL(skipped.res.cfo_refit_rotate_us, uint64_t(0));
    // The refit is still measured (it is what the decision is based on) but it
    // must not be applied.
    BOOST_CHECK_GT(skipped.res.cfo_refit_peaks, uint32_t(0));
    BOOST_CHECK_GT(skipped.res.cfo_refit_span_samples, uint64_t(0));
    // Same PSDU and same SFD/PHR coordinates: the rotation only changes the
    // carrier phase reference, not the frame geometry.
    BOOST_CHECK(skipped.psdu == legacy.psdu);
    BOOST_CHECK_EQUAL(skipped.res.ns_sfd.sfd_start_chip,
                      legacy.res.ns_sfd.sfd_start_chip);
    BOOST_CHECK_EQUAL(skipped.res.phr.psdu_length, legacy.res.phr.psdu_length);
    BOOST_CHECK_EQUAL(skipped.res.timing.preamble_start_sample,
                      legacy.res.timing.preamble_start_sample);
    BOOST_CHECK_EQUAL(skipped.res.sfd.sfd_start_sample, legacy.res.sfd.sfd_start_sample);
}

BOOST_AUTO_TEST_CASE(test_demod_core_cfo_refit_auto_falls_back_on_drift)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    BOOST_REQUIRE(!tmpl.empty());

    // A zero drift budget can never be met, so auto must keep the legacy path
    // and report the fallback instead of silently skipping.
    const auto strict = run_cfo_refit_case(iq, tmpl, 1, 1e-12, 64);
    BOOST_REQUIRE(strict.res.status == DemodStatus::Success);
    BOOST_CHECK(strict.res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(strict.res.cfo_refit_action, uint32_t(3));
    BOOST_CHECK_EQUAL(strict.res.stage_cfo_rotation_passes, uint32_t(2));
    BOOST_CHECK_EQUAL(strict.res.cfo.cfo_hz, strict.res.cfo_refit_hz);

    // A peak-count requirement the Stage-2 fit cannot satisfy must also fall
    // back (the fit only sees preamble_repetitions - 24 peaks).
    const auto too_few = run_cfo_refit_case(iq, tmpl, 1, 1e9, 100000);
    BOOST_REQUIRE(too_few.res.status == DemodStatus::Success);
    BOOST_CHECK_EQUAL(too_few.res.cfo_refit_action, uint32_t(3));
    BOOST_CHECK_EQUAL(too_few.res.stage_cfo_rotation_passes, uint32_t(2));

    // A generous budget and a reachable peak count must skip.
    const auto loose = run_cfo_refit_case(iq, tmpl, 1, 1e9, 1);
    BOOST_REQUIRE(loose.res.status == DemodStatus::Success);
    BOOST_CHECK(loose.res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(loose.res.cfo_refit_action, uint32_t(1));
    BOOST_CHECK_EQUAL(loose.res.stage_cfo_rotation_passes, uint32_t(1));
    BOOST_CHECK_EQUAL(loose.res.cfo.cfo_hz, loose.res.cfo_initial_hz);
    BOOST_CHECK(loose.psdu == strict.psdu);
}

BOOST_AUTO_TEST_CASE(test_demod_core_cfo_refit_legacy_default_is_unchanged)
{
    // The default profile must keep the pre-existing two-fit / two-rotation
    // behaviour so every historical measurement stays comparable.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    BOOST_REQUIRE(!tmpl.empty());

    const auto dflt = run_cfo_refit_case(iq, tmpl, 0, 0.10, 64);
    BOOST_REQUIRE(dflt.res.status == DemodStatus::Success);
    BOOST_CHECK_EQUAL(dflt.res.cfo.skipped_peaks, size_t(24));
    BOOST_CHECK_EQUAL(dflt.res.cfo.peaks_used, size_t(40));
    BOOST_CHECK_EQUAL(dflt.res.cfo_refit_action, uint32_t(0));
    BOOST_CHECK_EQUAL(dflt.res.stage_cfo_rotation_passes, uint32_t(2));
    BOOST_CHECK(dflt.res.cfo_refit_rotate_us <= dflt.res.stage_sfd_cfo_refit_us);
}

BOOST_AUTO_TEST_CASE(test_demod_core_second_rotation_stops_at_payload_start)
{
    // The payload pass derotates from min(sec_hi, base0) onwards.  The second
    // rotation used to run to the legacy 40-PHR-head budget, which sits PAST
    // the payload start, so ~35k samples were derotated twice with the same
    // CFO.  The bound must now stop at the real requirement (soft-chip prefix
    // and fused-PHR reads) and the decode must be unchanged.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    BOOST_REQUIRE(!tmpl.empty());

    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    prof.rx_fast_sfd = true;
    prof.rx_fast_cfo_rotation = true;
    prof.rx_softchip_on_demand = true;
    prof.rx_fused_payload_cmf = true;
    prof.rx_fused_phr_cmf = true;

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    const auto res =
        core::demodulate_one(iq.data(), iq.size(), prof, 1, 9984, 0, tmpl, scratch);
    BOOST_REQUIRE(res.status == DemodStatus::Success);
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));

    // The legacy budget would need 2 * initial_chip_budget samples; the real
    // requirement is smaller because the fused PHR reads raw IQ.
    const size_t legacy_span = 2 * rxopt::initial_chip_budget(prof);
    BOOST_CHECK_LT(res.cfo_refit_span_samples, legacy_span);
    BOOST_CHECK_GT(res.cfo_refit_span_samples, uint64_t(0));
}

BOOST_AUTO_TEST_CASE(test_demod_core_single_pass_rotation_reduces_rotated_samples)
{
    // Single-pass mode must derotate strictly fewer samples than the legacy
    // two-pass chain while producing the same PSDU.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    BOOST_REQUIRE(!tmpl.empty());

    auto make = [&](int mode) {
        auto prof = Qm35825Profile::Default();
        prof.sfd_mode = "ieee";
        prof.rx_fast_sfd = true;
        prof.rx_fast_cfo_rotation = true;
        prof.rx_softchip_on_demand = true;
        prof.rx_fused_payload_cmf = true;
        prof.rx_fused_phr_cmf = true;
        prof.cfo_refit_mode = mode;
        return prof;
    };
    auto run = [&](int mode) {
        core::DemodScratch scratch;
        scratch.reserve(iq.size());
        return core::demodulate_one(
            iq.data(), iq.size(), make(mode), 1, 9984, 0, tmpl, scratch);
    };

    const auto legacy = run(0);
    const auto single = run(2);
    BOOST_REQUIRE(legacy.status == DemodStatus::Success);
    BOOST_REQUIRE(single.status == DemodStatus::Success);
    BOOST_CHECK(legacy.payload.fcs_pass);
    BOOST_CHECK(single.payload.fcs_pass);
    BOOST_CHECK(single.payload.bytes == legacy.payload.bytes);
    BOOST_CHECK_EQUAL(single.cfo_refit_action, uint32_t(2));
    BOOST_CHECK_EQUAL(single.stage_cfo_rotation_passes + 1,
                      legacy.stage_cfo_rotation_passes);
    BOOST_CHECK_LT(single.stage_cfo_rotated_samples, legacy.stage_cfo_rotated_samples);
    BOOST_CHECK_EQUAL(single.cfo.cfo_hz, single.cfo_initial_hz);
    // Same frame geometry: the derotation only changes the phase reference.
    BOOST_CHECK_EQUAL(single.ns_sfd.sfd_start_chip, legacy.ns_sfd.sfd_start_chip);
    BOOST_CHECK_EQUAL(single.phr.psdu_length, legacy.phr.psdu_length);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r3_bad_input_fails_cleanly)
{
    std::vector<float> short_soft(5000, 0.0f); // no SFD / no PHR
    auto prof = Qm35825Profile::Default();
    core::DemodScratch scratch;
    scratch.reserve(8192);
    const auto sfd_seq = GetSfdSequence("ieee");
    NsSfdResult ns;
    BOOST_CHECK(
        !core::stage_ns_sfd(short_soft, prof, sfd_seq, kQm35ChipsPerSymbol, ns, scratch));
    // A fabricated NS-SFD beyond the buffer must fail the PHR stage.
    NsSfdResult bad_ns;
    bad_ns.ok = true;
    bad_ns.sfd_end_chip = 100000;
    bad_ns.polarity = 1;
    PhrResult phr;
    BOOST_CHECK(
        !core::stage_phr(short_soft, prof, bad_ns, kQm35ChipsPerSymbol, phr, scratch));
}

// ---------------------------------------------------------------------------
// R4: payload BPM-BPSK + RS + FCS (helperUWBPayloadDecode / hrpRS / CRC16).
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_demod_core_r4_crc16_matches_golden)
{
    // IEEE 802.15.4 CRC-16 (reflected 0x8408): golden FCS of the 125 data
    // bytes = 0x584b.
    std::vector<uint8_t> data;
    BOOST_REQUIRE(load_u8(golden_dir() + "/stage_payload_bytes.bin", data));
    BOOST_REQUIRE_EQUAL(data.size(), size_t(127));
    const uint16_t crc = core::detail::crc16_802154(data.data(), 125);
    BOOST_CHECK_EQUAL(crc, uint16_t(0x584b));
    // The received FCS in the last two bytes (little-endian) must match.
    const uint16_t recv =
        static_cast<uint16_t>(data[125]) | (static_cast<uint16_t>(data[126]) << 8);
    BOOST_CHECK_EQUAL(recv, uint16_t(0x584b));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r4_payload_fcs_matches_golden)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    const auto res = core::demodulate_one(
        iq.data(), iq.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
    BOOST_REQUIRE(res.status == DemodStatus::Success);
    BOOST_REQUIRE(res.phr.ok);
    BOOST_REQUIRE(res.payload.ok);
    // Golden: 127 payload bytes, FCS pass, received == calculated == 0x584b.
    std::vector<uint8_t> gbytes;
    BOOST_REQUIRE(load_u8(golden_dir() + "/stage_payload_bytes.bin", gbytes));
    BOOST_REQUIRE_EQUAL(res.payload.bytes.size(), gbytes.size());
    for (size_t i = 0; i < gbytes.size(); ++i)
        BOOST_CHECK_EQUAL(res.payload.bytes[i], gbytes[i]);
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
    BOOST_CHECK_EQUAL(res.payload.calculated_fcs, uint16_t(0x584b));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r4_bad_input_fails_cleanly)
{
    std::vector<float> short_soft(5000, 0.0f);
    auto prof = Qm35825Profile::Default();
    core::DemodScratch scratch;
    scratch.reserve(8192);
    PhrResult phr;
    phr.ok = true;
    phr.psdu_length = 127;
    NsSfdResult ns;
    ns.ok = true;
    ns.sfd_end_chip = 100000; // PHR end beyond the buffer
    PayloadResult pay;
    BOOST_CHECK(!core::stage_payload_fcs(
        short_soft, prof, phr, ns, kQm35ChipsPerSymbol, pay, scratch));
}

// ---------------------------------------------------------------------------
// P0: CFO sweep — estimate + full demod across 0/±1/±5/±10/±25/±50 kHz.
// The full chain (timing → cfo → sfd(derotated) → cir → ns_sfd → phr →
// payload) must still decode the golden payload at every CFO.
// ---------------------------------------------------------------------------


BOOST_AUTO_TEST_CASE(test_demod_core_p0_cfo_sweep_clean)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const double kfreqs[] = { -50000.0, -25000.0, -10000.0, -5000.0, -1000.0, 0.0,
                              1000.0,   5000.0,   10000.0,  25000.0, 50000.0 };
    for (double f : kfreqs) {
        auto res = demod_window_with_cfo(iq, f);
        BOOST_TEST_MESSAGE("CFO sweep f=" << f << " est=" << res.cfo.cfo_hz
                                          << " status=" << (int)res.status
                                          << " sfd_metric=" << res.sfd.metric
                                          << " fcs=" << res.payload.fcs_pass);
        // CFO estimate within 5% + 50 Hz of the injected value.
        if (f == 0.0)
            BOOST_CHECK_LT(std::abs(res.cfo.cfo_hz), 50.0);
        else
            BOOST_CHECK_CLOSE(res.cfo.cfo_hz, f, 5.0);
        // The full chain must still decode the golden payload.
        BOOST_CHECK(res.status == DemodStatus::Success);
        BOOST_CHECK(res.payload.fcs_pass);
        BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
        // SFD position must be stable across the sweep (golden: 75008).
        BOOST_CHECK_CLOSE(static_cast<double>(res.sfd.sfd_start_sample), 75008.0, 0.2);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_p0_cfo_sweep_awgn)
{
    // Same sweep on a 20 dB-SNR AWGN version of the golden window: CFO
    // estimate must stay accurate and the payload must still decode.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const double snr_db = 20.0;
    double sp = 0.0;
    for (size_t i = 0; i < iq.size(); ++i)
        sp += static_cast<double>(std::norm(iq[i]));
    sp /= static_cast<double>(iq.size());
    const double sigma = std::sqrt(sp / std::pow(10.0, snr_db / 10.0) / 2.0);
    std::mt19937 rng(12345);
    std::normal_distribution<float> g(0.0f, static_cast<float>(sigma));
    std::vector<gr_complex> iqn(iq.size());
    for (size_t i = 0; i < iq.size(); ++i)
        iqn[i] = iq[i] + gr_complex(g(rng), g(rng));

    const double kfreqs[] = { -10000.0, 0.0, 10000.0 };
    for (double f : kfreqs) {
        auto res = demod_window_with_cfo(iqn, f);
        BOOST_TEST_MESSAGE("CFO sweep AWGN f=" << f << " est=" << res.cfo.cfo_hz
                                               << " status=" << (int)res.status
                                               << " fcs=" << res.payload.fcs_pass);
        // Wider tolerance under noise: 10% + 200 Hz.
        if (f == 0.0)
            BOOST_CHECK_LT(std::abs(res.cfo.cfo_hz), 200.0);
        else
            BOOST_CHECK_CLOSE(res.cfo.cfo_hz, f, 10.0);
        BOOST_CHECK(res.status == DemodStatus::Success);
        BOOST_CHECK(res.payload.fcs_pass);
        BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_sfo_from_cfo_default_and_override)
{
    const auto def = Qm35825Profile::Default();
    BOOST_CHECK(def.rx_sfo_from_cfo);
    BOOST_CHECK_CLOSE(def.carrier_freq_hz, 6489.6e6, 1e-12);

    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();

    // Coupled 27 ppm crystal offset: CFO = +175 kHz at 6489.6 MHz, SFO = +27 ppm.
    constexpr double kCfoHz = 175000.0;
    auto coupled =
        apply_cfo_sfo_constraint(iq, kCfoHz, def.sample_rate, def.carrier_freq_hz);
    auto p_on = Qm35825Profile::Default();
    p_on.sfd_mode = "ieee";
    core::DemodScratch scratch;
    scratch.reserve(coupled.size());
    const auto res_on = core::demodulate_one(
        coupled.data(), coupled.size(), p_on, 1, 9984, 0, tmpl, scratch);
    BOOST_CHECK(res_on.status == DemodStatus::Success);
    BOOST_CHECK(res_on.payload.fcs_pass);
    const double expected_period = static_cast<double>(kQm35SamplesPerSymbol) /
                                   (1.0 + res_on.cfo.cfo_hz / p_on.carrier_freq_hz);
    BOOST_CHECK_CLOSE(res_on.timing.measured_period, expected_period, 1e-9);

    // Explicit opt-out: decoupled CFO-only signal (0 ppm SFO) decodes cleanly
    // when rx_sfo_from_cfo = false, using Stage-1 timing OLS period (~1016).
    auto decoupled = apply_cfo_sfo_constraint(iq, kCfoHz, def.sample_rate, 0.0);
    auto p_off = p_on;
    p_off.rx_sfo_from_cfo = false;
    const auto res_off = core::demodulate_one(
        decoupled.data(), decoupled.size(), p_off, 2, 9984, 0, tmpl, scratch);
    BOOST_CHECK(res_off.status == DemodStatus::Success);
    BOOST_CHECK(res_off.payload.fcs_pass);
    BOOST_CHECK_SMALL(std::abs(res_off.timing.measured_period -
                               static_cast<double>(kQm35SamplesPerSymbol)),
                      0.05);
}

// ---------------------------------------------------------------------------
// P3: SFD ROI compression robustness — the search center adapts to the last
// tracked SYNC + measured_period with a narrow ±half-width window, so SFD must
// still be located under timing offset, SFO, AWGN and a multipath echo.
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_demod_core_p3_sfd_offset)
{
    // Shift the whole golden window right by +7 samples; through the full chain
    // the SFD must be found at 75008 + 7.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    std::vector<gr_complex> iq2(iq.size() + 7, gr_complex(0.0f, 0.0f));
    std::copy(iq.begin(), iq.end(), iq2.begin() + 7);
    core::DemodScratch scratch;
    scratch.reserve(iq2.size());
    auto tmpl = load_reference_template();
    auto res =
        core::demodulate_one(iq2.data(), iq2.size(), prof, 1, 9984 + 7, 0, tmpl, scratch);
    BOOST_CHECK(res.status == DemodStatus::Success);
    BOOST_CHECK_CLOSE(static_cast<double>(res.sfd.sfd_start_sample), 75015.0, 0.1);
    BOOST_CHECK_EQUAL(res.timing.preamble_start_sample,
                      res.sfd.sfd_start_sample -
                          static_cast<int64_t>(std::llround(
                              static_cast<double>(prof.preamble_repetitions) *
                              res.timing.measured_period)));
}

BOOST_AUTO_TEST_CASE(test_demod_core_p3_sfd_sfo_adaptation)
{
    // Apply the full preamble drift before rounding the SFD search centre.
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    BOOST_REQUIRE(
        core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch));
    tr.measured_period = 1016.0 * (1.0 + 200e-6); // ~13-sample drift over 64 SYNCs
    const auto sfd_seq = GetSfdSequence("ieee");
    SfdResult sr;
    const std::complex<float>* sfd_rx =
        scratch.derotated.empty() ? iq.data() : scratch.derotated.data();
    BOOST_REQUIRE(
        core::stage_sfd(sfd_rx, iq.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_EQUAL(sr.expected_start_sample,
                      tr.preamble_start_sample +
                          std::llround(prof.preamble_repetitions * tr.measured_period));
    BOOST_CHECK_CLOSE(static_cast<double>(sr.sfd_start_sample), 75008.0, 0.5);
    BOOST_CHECK_GE(sr.metric, 0.9f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_p3_sfd_awgn)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const double snr_db = 20.0;
    double sp = 0.0;
    for (size_t i = 0; i < iq.size(); ++i)
        sp += static_cast<double>(std::norm(iq[i]));
    sp /= static_cast<double>(iq.size());
    const double sigma = std::sqrt(sp / std::pow(10.0, snr_db / 10.0) / 2.0);
    std::mt19937 rng(7);
    std::normal_distribution<float> g(0.0f, static_cast<float>(sigma));
    std::vector<gr_complex> iqn(iq.size());
    for (size_t i = 0; i < iq.size(); ++i)
        iqn[i] = iq[i] + gr_complex(g(rng), g(rng));

    core::DemodScratch scratch;
    scratch.reserve(iqn.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    SfdResult sr;
    std::vector<gr_complex> tmpl;
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    BOOST_REQUIRE(
        core::stage_timing(iqn.data(), iqn.size(), prof, tmpl, 9984, tr, scratch));
    const auto sfd_seq = GetSfdSequence("ieee");
    const std::complex<float>* sfd_rx =
        scratch.derotated.empty() ? iqn.data() : scratch.derotated.data();
    BOOST_REQUIRE(
        core::stage_sfd(sfd_rx, iqn.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_CLOSE(static_cast<double>(sr.sfd_start_sample), 75008.0, 2.0);
    BOOST_CHECK_GE(sr.metric, 0.5f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_p3_sfd_multipath)
{
    // First path + a 0.4x echo delayed by 150 samples (≈1.5 chips): the SFD
    // must still lock onto the first-path peak.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    std::vector<gr_complex> m(iq.size() + 150, gr_complex(0.0f, 0.0f));
    std::copy(iq.begin(), iq.end(), m.begin());
    for (size_t i = 0; i + 150 < m.size() && i < iq.size(); ++i)
        m[i + 150] += 0.4f * iq[i];

    core::DemodScratch scratch;
    scratch.reserve(m.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    SfdResult sr;
    std::vector<gr_complex> tmpl;
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));
    BOOST_REQUIRE(core::stage_timing(m.data(), m.size(), prof, tmpl, 9984, tr, scratch));
    const auto sfd_seq = GetSfdSequence("ieee");
    const std::complex<float>* sfd_rx =
        scratch.derotated.empty() ? m.data() : scratch.derotated.data();
    BOOST_REQUIRE(
        core::stage_sfd(sfd_rx, m.size(), prof, tr, sfd_seq, tmpl, sr, scratch));
    BOOST_CHECK_CLOSE(static_cast<double>(sr.sfd_start_sample), 75008.0, 2.0);
    BOOST_CHECK_GE(sr.metric, 0.5f);
}

// ---------------------------------------------------------------------------
// P2: CIR soft-chip FIR SIMD — kernel consistency vs multi-acc reference,
// golden soft-chip / CIR tolerances, and a small microbench of the three
// candidates (a=MultiAcc8, b=VOLK, c=AVX2 fixed-38).
// ---------------------------------------------------------------------------
