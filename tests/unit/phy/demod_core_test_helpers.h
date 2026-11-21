// SPDX-License-Identifier: GPL-3.0-or-later
//
// Shared fixtures for the split realtime-demodulator core QA shards.
//
// The original single-TU qa_uwb_demod_core.cc (5,254 lines / 90 cases) was
// split into four executables so the 20-core QA run no longer serialises on
// one 32 s compile and one 7 s test process:
//
//   qa_uwb_demod_core_timing_cfo_sfd.cc   R1 timing/CFO/SFD
//   qa_uwb_demod_core_cir_pipeline.cc     R2-R4 CIR/softchips/payload + refit
//   qa_uwb_demod_core_kernel_grid.cc      P0/P3 sweeps + SIMD kernel equivalence
//   qa_uwb_demod_core_fec_joint_ml.cc     FEC + Wave-1 joint-ML Viterbi
//
// All common loaders and the CFO/SFO waveform helpers live here; every test
// case and assertion was moved 1:1, nothing was dropped.

#ifndef OPENUWB_QA_DEMOD_CORE_TEST_HELPERS_H
#define OPENUWB_QA_DEMOD_CORE_TEST_HELPERS_H

#include <boost/test/unit_test.hpp>

#include <gnuradio/uwb/uwb_demod_core.h>
#include <gnuradio/uwb/uwb_gr_complex.h>
#include <gnuradio/uwb/uwb_phy_profile.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

using namespace openuwb::demod;

// Load a CF32 (interleaved float I/Q) file into gr_complex.
bool load_cf32(const std::string& path, std::vector<gr_complex>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const size_t bytes = static_cast<size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    if (bytes == 0 || bytes % 8 != 0)
        return false;
    out.resize(bytes / 8);
    f.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(bytes));
    return true;
}

// Path to the golden directory (repo-root/testdata/realtime_demod_golden).
std::string golden_dir()
{
    // QA runs with CWD = gr-uwb/build; golden lives at repo-root/testdata/...
    return "../../../testdata/realtime_demod_golden";
}

// Load a raw float32 (f32) file.
bool load_f32(const std::string& path, std::vector<float>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const size_t bytes = static_cast<size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    if (bytes == 0 || bytes % 4 != 0)
        return false;
    out.resize(bytes / 4);
    f.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(bytes));
    return true;
}

// Load a raw int8 file.
bool load_i8(const std::string& path, std::vector<int8_t>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const size_t bytes = static_cast<size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    if (bytes == 0)
        return false;
    out.resize(bytes);
    f.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(bytes));
    return true;
}

// Load a raw uint8 file.
bool load_u8(const std::string& path, std::vector<uint8_t>& out)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    f.seekg(0, std::ios::end);
    const size_t bytes = static_cast<size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    if (bytes == 0)
        return false;
    out.resize(bytes);
    f.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(bytes));
    return true;
}

// The 1016-sample code-9 SYNC waveform template (reference_preamble.bin).
std::vector<gr_complex> load_reference_template()
{
    std::vector<gr_complex> tmpl;
    load_cf32("../../../testdata/reference_preamble.bin", tmpl);
    return tmpl;
}

} // namespace

// ---------------------------------------------------------------------------
// Integer chip-grid admission (fused CMF / AVX2 soft-chip kernel).
//
// The grid tolerance is derived from the worst-case payload drift, not from an
// arbitrary epsilon.  These cases pin the two ends that matter on real
// hardware: a 27 ppm clock must still be rejected, a sub-ppm clock accepted,
// and the boundary must sit at the documented drift budget.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// P0 CFO/SFO waveform helpers (originally defined just before the P0 cases).
// Shared by the cir_pipeline and kernel_grid shards.
// ---------------------------------------------------------------------------
namespace {

// Apply co-sourced crystal CFO and matching SFO (rate = 1 + f_hz / carrier_freq_hz)
// to a baseband signal at sample_rate. Preserves buffer length (zero-padded).
std::vector<gr_complex> apply_cfo_sfo_constraint(const std::vector<gr_complex>& iq,
                                                 double f_hz,
                                                 double sample_rate,
                                                 double carrier_freq_hz)
{
    std::vector<gr_complex> out(iq.size(), gr_complex(0.0f, 0.0f));
    if (iq.empty())
        return out;
    const double rate = (carrier_freq_hz > 0.0) ? (1.0 + f_hz / carrier_freq_hz) : 1.0;
    const double w = 2.0 * M_PI * f_hz / sample_rate;
    for (size_t i = 0; i < iq.size(); ++i) {
        gr_complex s(0.0f, 0.0f);
        if (rate == 1.0) {
            s = iq[i];
        } else {
            const double x = static_cast<double>(i) * rate;
            const size_t i0 = static_cast<size_t>(x);
            if (i0 + 1 < iq.size()) {
                const double frac = x - static_cast<double>(i0);
                s = iq[i0] * static_cast<float>(1.0 - frac) +
                    iq[i0 + 1] * static_cast<float>(frac);
            } else if (i0 < iq.size()) {
                s = iq[i0];
            }
        }
        if (f_hz != 0.0) {
            const double ph = w * static_cast<double>(i);
            s *= gr_complex(static_cast<float>(std::cos(ph)),
                            static_cast<float>(std::sin(ph)));
        }
        out[i] = s;
    }
    return out;
}

// Inject a carrier rotation of +f Hz (e^{+j2π f n/fs}) together with its
// co-sourced crystal SFO (1 + f/fc) into the golden window and run the full
// demod chain (IEEE SFD, seeded at the golden preamble).
DemodResult demod_window_with_cfo(const std::vector<gr_complex>& iq, double f_hz)
{
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    auto iqc = apply_cfo_sfo_constraint(iq, f_hz, prof.sample_rate, prof.carrier_freq_hz);
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto tmpl = load_reference_template();
    return core::demodulate_one(iqc.data(), iqc.size(), prof, 1, 9984, 0, tmpl, scratch);
}

} // namespace

// ---------------------------------------------------------------------------
// Shared R2 CIR/soft-chip driver (used by the timing and kernel-grid shards).
// ---------------------------------------------------------------------------
// Run timing + CFO + CIR/soft chips on `iq` and return the CIR + soft chips.
static bool run_cir_softchips(const std::vector<gr_complex>& iq,
                              core::DemodScratch& scratch,
                              Qm35825Profile& prof,
                              TimingResult& tr,
                              CfoResult& cf,
                              CirResult& cir)
{
    std::vector<gr_complex> tmpl;
    if (!load_cf32("../../../testdata/reference_preamble.bin", tmpl))
        return false;
    if (!core::stage_timing(iq.data(), iq.size(), prof, tmpl, 9984, tr, scratch))
        return false;
    if (!core::stage_cfo(iq.data(), iq.size(), prof, tr, cf, scratch))
        return false;
    std::vector<int8_t> code9(kPreambleCode9.begin(), kPreambleCode9.end());
    return core::stage_cir_softchips(
        scratch.derotated.data(), iq.size(), prof, tr, code9, cir, scratch);
}

#endif // OPENUWB_QA_DEMOD_CORE_TEST_HELPERS_H
