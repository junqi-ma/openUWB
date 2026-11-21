/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * End-to-end TX->RX demodulator loopback QA split out of qa_uwb_hrp_mod_core.cc
 * so UWB_DEV_SCOPE=tx can build and test the HRP modulator without compiling
 * the receive demodulation pipeline.
 */

#include "hrp_mod_test_helpers.h"

#include <boost/test/unit_test.hpp>

#include <gnuradio/uwb/uwb_demod_core.h>
#include <gnuradio/uwb/uwb_gr_complex.h>
#include <gnuradio/uwb/uwb_hrp_mod_core.h>
#include <gnuradio/uwb/uwb_phy_profile.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace {

using openuwb::demod::Qm35825Profile;
using openuwb::demod::core::demodulate_one;
using openuwb::hrp_mod_test::golden_dir;
using openuwb::hrp_mod_test::load_cf32;
using openuwb::hrp_mod_test::load_f32;
using openuwb::hrp_mod_test::load_template;
using openuwb::hrp_mod_test::load_u8;
using openuwb::hrp_mod_test::kTestdata;
using openuwb::mod::HrpModConfig;
using openuwb::mod::HrpModScratch;
using openuwb::mod::modulate_one;
using openuwb::mod::packet_samples_998p4;
using openuwb::mod::PulseShape;

void run_external_live_template_demod(const std::string& rel_path,
                                      size_t n_taps, size_t center,
                                      size_t tail_extra, uint32_t seed)
{
    std::vector<float> taps;
    BOOST_REQUIRE(load_f32(std::string(kTestdata) + rel_path, taps));
    BOOST_REQUIRE_EQUAL(taps.size(), n_taps);

    auto bytes = openuwb::mod::make_random_psdu(8, seed, true);
    HrpModConfig cfg;
    cfg.sfd_mode = "4z2";
    cfg.pulse.shape = PulseShape::External;
    cfg.pulse.external_n_taps = taps.size();
    cfg.pulse.external_center = center;
    cfg.external_taps = taps;
    HrpModScratch scratch;
    const size_t want = packet_samples_998p4(64, 8, bytes.size(), false,
                                             &cfg.pulse);
    BOOST_CHECK_EQUAL(want,
                      packet_samples_998p4(64, 8, bytes.size()) + tail_extra);
    scratch.reserve(bytes.size(), want);
    std::vector<gr_complex> tx(want);
    size_t n = 0;
    BOOST_REQUIRE(modulate_one(bytes.data(), bytes.size(), cfg, scratch,
                               tx.data(), tx.size(), n));

    std::vector<gr_complex> tmpl(tx.begin(), tx.begin() + 1016);
    double e = 0.0;
    for (const auto& z : tmpl)
        e += static_cast<double>(z.real()) * z.real() +
             static_cast<double>(z.imag()) * z.imag();
    const float inv = (e > 0.0) ? static_cast<float>(1.0 / std::sqrt(e)) : 1.f;
    for (auto& z : tmpl)
        z *= inv;

    const size_t t0 = 9984;
    std::vector<gr_complex> iq(t0 + n + 8192, gr_complex(0.f, 0.f));
    std::copy(tx.begin(), tx.end(), iq.begin() + static_cast<std::ptrdiff_t>(t0));
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "4z2";
    openuwb::demod::core::DemodScratch dscratch;
    dscratch.reserve(iq.size());
    const auto res = demodulate_one(iq.data(), iq.size(), prof, 1,
                                    static_cast<int64_t>(t0), 0, tmpl, dscratch);
    BOOST_REQUIRE_MESSAGE(res.status == openuwb::demod::DemodStatus::Success,
                          "external live-template demod status=" +
                              std::to_string(static_cast<int>(res.status)));
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_REQUIRE_EQUAL(res.payload.bytes.size(), bytes.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(res.payload.bytes.begin(),
                                  res.payload.bytes.end(),
                                  bytes.begin(), bytes.end());
}

} // namespace

BOOST_AUTO_TEST_CASE(test_hrp_mod_gaussian_live_template_demod)
{
    auto bytes = openuwb::mod::make_random_psdu(8, 4242u, true);
    HrpModConfig cfg;
    cfg.sfd_mode = "4z2";
    cfg.pulse.shape = PulseShape::Gaussian;
    cfg.pulse.gaussian_sigma_ns = 2.5f;
    HrpModScratch scratch;
    const size_t want = packet_samples_998p4(64, 8, bytes.size(), false,
                                             &cfg.pulse);
    scratch.reserve(bytes.size(), want);
    std::vector<gr_complex> tx(want);
    size_t n = 0;
    BOOST_REQUIRE(modulate_one(bytes.data(), bytes.size(), cfg, scratch,
                               tx.data(), tx.size(), n));

    std::vector<gr_complex> tmpl(tx.begin(), tx.begin() + 1016);
    double e = 0.0;
    for (const auto& z : tmpl)
        e += static_cast<double>(z.real()) * z.real() +
             static_cast<double>(z.imag()) * z.imag();
    const float inv = (e > 0.0) ? static_cast<float>(1.0 / std::sqrt(e)) : 1.f;
    for (auto& z : tmpl)
        z *= inv;

    const size_t t0 = 9984;
    std::vector<gr_complex> iq(t0 + n + 8192, gr_complex(0.f, 0.f));
    std::copy(tx.begin(), tx.end(), iq.begin() + static_cast<std::ptrdiff_t>(t0));
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "4z2";
    openuwb::demod::core::DemodScratch dscratch;
    dscratch.reserve(iq.size());
    const auto res = demodulate_one(iq.data(), iq.size(), prof, 1,
                                    static_cast<int64_t>(t0), 0, tmpl, dscratch);
    BOOST_REQUIRE_MESSAGE(res.status == openuwb::demod::DemodStatus::Success,
                          "gaussian live-template demod status=" +
                              std::to_string(static_cast<int>(res.status)));
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_REQUIRE_EQUAL(res.payload.bytes.size(), bytes.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(res.payload.bytes.begin(),
                                  res.payload.bytes.end(),
                                  bytes.begin(), bytes.end());
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_external_live_template_demod)
{
    run_external_live_template_demod(
        "/uwb_hrp_tx/pulse_minphase_rc100_215.f32", 257, 4, 255, 777u);
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_minphase_live_template_demod)
{
    run_external_live_template_demod(
        "/uwb_hrp_tx/pulse_minphase_rc160_240.f32", 1025, 4, 1023, 778u);
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_external_default_live_template_demod)
{
    run_external_live_template_demod(
        "/uwb_hrp_tx/pulse_trunc_linear_rc183_240.f32", 1025, 32, 1023,
        779u);
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_rs_roundtrip_demod)
{
    std::vector<uint8_t> bytes;
    BOOST_REQUIRE(load_u8(golden_dir() + "/stage_payload_bytes.bin", bytes));
    BOOST_REQUIRE_EQUAL(bytes.size(), size_t(127));

    HrpModScratch scratch;
    scratch.reserve(127, 249280);
    openuwb::mod::detail::bytes_to_lsb_bits(bytes.data(), bytes.size(),
                                            scratch.psdu_bits);
    BOOST_REQUIRE(openuwb::mod::detail::rs_encode_stream(
        scratch.psdu_bits.data(), scratch.psdu_bits.size(), scratch.rs_cw));

    std::vector<int8_t> decoded;
    BOOST_REQUIRE(openuwb::demod::core::detail::rs_decode_stream(
        scratch.rs_cw, scratch.psdu_bits.size(), decoded));
    BOOST_CHECK_EQUAL_COLLECTIONS(decoded.begin(), decoded.end(),
                                  scratch.psdu_bits.begin(),
                                  scratch.psdu_bits.end());
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_demodulate_one_fcs)
{
    std::vector<uint8_t> bytes;
    BOOST_REQUIRE(load_u8(golden_dir() + "/stage_payload_bytes.bin", bytes));

    HrpModConfig cfg;
    HrpModScratch scratch;
    scratch.reserve(127, 249280);
    std::vector<gr_complex> tx(249280);
    size_t n = 0;
    BOOST_REQUIRE(modulate_one(bytes.data(), bytes.size(), cfg, scratch,
                               tx.data(), tx.size(), n));

    const size_t t0 = 9984;
    std::vector<gr_complex> iq(t0 + n + 4096, gr_complex(0.f, 0.f));
    std::copy(tx.begin(), tx.end(), iq.begin() + static_cast<std::ptrdiff_t>(t0));

    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    openuwb::demod::core::DemodScratch dscratch;
    dscratch.reserve(iq.size());
    const auto res = demodulate_one(iq.data(), iq.size(), prof, 1,
                                    static_cast<int64_t>(t0), 0,
                                    load_template(), dscratch);
    BOOST_REQUIRE(res.status == openuwb::demod::DemodStatus::Success);
    BOOST_REQUIRE(res.payload.ok);
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_REQUIRE_EQUAL(res.payload.bytes.size(), bytes.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(res.payload.bytes.begin(),
                                  res.payload.bytes.end(),
                                  bytes.begin(), bytes.end());
    BOOST_CHECK_EQUAL(res.phr.psdu_length, 127u);
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_4z2_no_sts_demod_fcs)
{
    auto bytes = openuwb::mod::make_random_psdu(8, 20260904u, true);
    BOOST_REQUIRE_EQUAL(bytes.size(), size_t(10));
    HrpModConfig cfg;
    cfg.sfd_mode = "4z2";
    HrpModScratch scratch;
    const size_t want = packet_samples_998p4(64, 8, bytes.size());
    scratch.reserve(bytes.size(), want);
    std::vector<gr_complex> tx(want);
    size_t n = 0;
    BOOST_REQUIRE(modulate_one(bytes.data(), bytes.size(), cfg, scratch,
                               tx.data(), tx.size(), n));
    const size_t t0 = 9984;
    std::vector<gr_complex> iq(t0 + n + 8192, gr_complex(0.f, 0.f));
    std::copy(tx.begin(), tx.end(), iq.begin() + static_cast<std::ptrdiff_t>(t0));
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "4z2";
    openuwb::demod::core::DemodScratch dscratch;
    dscratch.reserve(iq.size());
    const auto res = demodulate_one(iq.data(), iq.size(), prof, 1,
                                    static_cast<int64_t>(t0), 0,
                                    load_template(), dscratch);
    BOOST_REQUIRE_MESSAGE(res.status == openuwb::demod::DemodStatus::Success,
                          "4z2 demod status=" +
                              std::to_string(static_cast<int>(res.status)));
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_REQUIRE_EQUAL(res.payload.bytes.size(), bytes.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(res.payload.bytes.begin(),
                                  res.payload.bytes.end(),
                                  bytes.begin(), bytes.end());
}

BOOST_AUTO_TEST_CASE(test_hrp_mod_code10_ieee_fcs)
{
    auto bytes = openuwb::mod::make_random_psdu(8, 7u, true);
    HrpModConfig cfg;
    cfg.code_index = 10;
    HrpModScratch scratch;
    const size_t want = packet_samples_998p4(64, 8, bytes.size());
    scratch.reserve(bytes.size(), want);
    std::vector<gr_complex> tx(want);
    size_t n = 0;
    BOOST_REQUIRE(modulate_one(bytes.data(), bytes.size(), cfg, scratch,
                               tx.data(), tx.size(), n));
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    prof.code_index = 10;
    std::vector<gr_complex> tmpl;
    BOOST_REQUIRE(load_cf32(std::string(kTestdata) +
                                "/reference_preamble.bin",
                            tmpl));
    // code-10 needs its own SYNC template; first generated symbol, L2-norm.
    tmpl.assign(tx.begin(), tx.begin() + 1016);
    double e = 0.0;
    for (const auto& z : tmpl)
        e += static_cast<double>(z.real()) * z.real() +
             static_cast<double>(z.imag()) * z.imag();
    const float inv = (e > 0.0) ? static_cast<float>(1.0 / std::sqrt(e)) : 1.f;
    for (auto& z : tmpl)
        z *= inv;
    const size_t t0 = 9984;
    std::vector<gr_complex> iq(t0 + n + 8192, gr_complex(0.f, 0.f));
    std::copy(tx.begin(), tx.end(), iq.begin() + static_cast<std::ptrdiff_t>(t0));
    openuwb::demod::core::DemodScratch dscratch;
    dscratch.reserve(iq.size());
    const auto res = demodulate_one(iq.data(), iq.size(), prof, 1,
                                    static_cast<int64_t>(t0), 0, tmpl, dscratch);
    BOOST_REQUIRE_MESSAGE(res.status == openuwb::demod::DemodStatus::Success,
                          "code10 demod status=" +
                              std::to_string(static_cast<int>(res.status)));
    BOOST_CHECK(res.payload.fcs_pass);
}
