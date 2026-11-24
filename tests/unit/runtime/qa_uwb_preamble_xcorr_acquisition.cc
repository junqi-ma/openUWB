/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * QA for the folded-preamble cross-correlation acquisition frontend
 * (openuwb/runtime/preamble_xcorr_acquisition.h).
 *
 * Acceptance items from
 * docs/OpenCode执行方案_Acquire前端折叠前导码互相关设计与验收_2026-11-24.md §4/§6
 * (Subagent-3):
 *   1. pure AWGN and strong CW tones (-232.96 MHz, -124.4 MHz, multi-tone)
 *      across >= 50 ms of 737.28 MS/s samples emit no seed;
 *   2. 256-SYNC code-10 preamble at CFO = 0/±50/±100 kHz and in-band
 *      SNR = -12/-14 dB is detected 100% of the time with exactly one seed
 *      per packet (no duplicate);
 *   3. leading-edge alignment: the seed stays within 16 * 751 samples of the
 *      true preamble start for the four 128-symbol hop phases;
 *   4. a strong bandlimited code-9 preamble scanned with code_index = 10 is
 *      rejected by the Stage-2 relative code-identity gate
 *      (code_par_ratio_min; the optional absolute code_rho_min gate is off by
 *      default) at -8/0/+10 dB in-band SNR; code_rejected counts the
 *      rejections, so the pass is not an "energy too low" pass; the same 0 dB
 *      interferer is detected with code_index = 9 (matched template), so the
 *      rejection is code-specific;
 *   5. ring-wrap fold, post-seed holdoff and reset() state clearing;
 *   6. the AVX2 and scalar kernels agree on par / peak_power / code_rho /
 *      peak_tau.
 *
 * Synthetic preamble: testdata/reference_preamble_code10_998p4.cf32 (one
 * 998.4 MS/s code-10 SYNC symbol) is resampled to the native 737.28 MS/s grid
 * with a periodic Kaiser-windowed sinc (48/65, half-length 10, beta 5, the same
 * parameters as testdata/generate_qm35_reference_737p28.py) and repeated 256
 * times on the exact period grid llround(n * 750.276923...).  The wrong-code
 * interferer uses the same path on testdata/reference_preamble.bin (code-9).
 *
 * SNR is in-band, as in qa_uwb_cfo_ml: mean signal power over the preamble
 * span versus complex AWGN power (noise added over the whole stream).
 */
#define BOOST_TEST_MODULE qa_uwb_preamble_xcorr_acquisition
#include <boost/test/unit_test.hpp>

#include <openuwb/runtime/preamble_xcorr_acquisition.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace {

using cf = std::complex<float>;
using sc16 = std::complex<int16_t>;
using openuwb::live::PreambleXcorrAcquisition;
using openuwb::live::PreambleXcorrConfig;
using openuwb::live::PreambleXcorrDiagnostics;
using openuwb::live::PreambleXcorrStats;

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 737280000.0;                 // native X410 rate
constexpr double kPeriod = 750.276923076923;        // 1016 * 48 / 65
constexpr size_t kSymLen = 750;
constexpr size_t kNativeSymLen = 751;               // one native SYNC symbol
constexpr size_t kPreambleReps = 256;
constexpr size_t kHopSymbols = 128;
constexpr int64_t kHopSamples = 96035;              // llround(128 * kPeriod)
constexpr int64_t kSeedTolerance = 16 * 751;        // design §4.2.5
constexpr size_t kRingSignal = 1u << 19;            // signal tests (2 MB)
constexpr size_t kRingNoise = 1u << 18;             // streaming noise test
constexpr size_t kRingWrap = 1u << 16;              // ring-wrap test (>= 48767)
constexpr size_t kBlock = 8192;

#ifndef UWB_TESTDATA_DIR
#define UWB_TESTDATA_DIR "../../../testdata"
#endif

// ---------------------------------------------------------------------------
// Fixtures / signal construction
// ---------------------------------------------------------------------------
std::string testdata_path(const std::string& rel)
{
    const std::string from_def = std::string(UWB_TESTDATA_DIR) + "/" + rel;
    {
        std::ifstream f(from_def, std::ios::binary);
        if (f)
            return from_def;
    }
    const char* prefixes[] = { "../../../testdata/", "../../testdata/",
                               "../testdata/", "testdata/" };
    for (const char* p : prefixes) {
        std::ifstream f(std::string(p) + rel, std::ios::binary);
        if (f)
            return std::string(p) + rel;
    }
    return from_def;
}

std::vector<cf> load_cf32(const std::string& rel)
{
    const std::string path = testdata_path(rel);
    std::ifstream f(path, std::ios::binary);
    BOOST_REQUIRE_MESSAGE(f, "cannot open " << path);
    f.seekg(0, std::ios::end);
    const size_t bytes = static_cast<size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    BOOST_REQUIRE_MESSAGE(bytes > 0 && bytes % 8 == 0,
                          "bad cf32 size in " << path);
    std::vector<cf> out(bytes / 8);
    f.read(reinterpret_cast<char*>(out.data()),
           static_cast<std::streamsize>(bytes));
    BOOST_REQUIRE_MESSAGE(static_cast<bool>(f), "short read on " << path);
    return out;
}

double sinc(double d)
{
    if (std::fabs(d) < 1e-12)
        return 1.0;
    const double a = kPi * d;
    return std::sin(a) / a;
}

// Periodic Kaiser-windowed sinc interpolation of the 998.4 MS/s SYNC template
// onto the native 737.28 MS/s grid (positions n * 65/48).
std::vector<cf> resample_symbol_65_48(const std::vector<cf>& ref,
                                      size_t n_out)
{
    constexpr int kHalf = 10;
    constexpr double kBeta = 5.0;
    const double i0_beta = std::cyl_bessel_i(0.0, kBeta);
    const int64_t n = static_cast<int64_t>(ref.size());
    std::vector<cf> out(n_out);
    for (size_t i = 0; i < n_out; ++i) {
        const double u = static_cast<double>(i) * 65.0 / 48.0;
        const int64_t i0 = static_cast<int64_t>(std::floor(u));
        double re = 0.0;
        double im = 0.0;
        double wsum = 0.0;
        for (int k = -kHalf + 1; k <= kHalf; ++k) {
            const double d = u - static_cast<double>(i0 + k);
            const double r = d / kHalf;
            const double win = std::cyl_bessel_i(
                0.0, kBeta * std::sqrt(std::max(0.0, 1.0 - r * r))) /
                i0_beta;
            const double w = sinc(d) * win;
            const int64_t idx = ((i0 + k) % n + n) % n;
            re += static_cast<double>(ref[static_cast<size_t>(idx)].real()) * w;
            im += static_cast<double>(ref[static_cast<size_t>(idx)].imag()) * w;
            wsum += w;
        }
        out[i] = cf(static_cast<float>(re / wsum), static_cast<float>(im / wsum));
    }
    return out;
}

double symbol_peak(const std::vector<cf>& x)
{
    double p = 0.0;
    for (const cf& v : x)
        p = std::max(p, static_cast<double>(std::abs(v)));
    return p;
}

// Unit-energy SYNC symbol: the preamble mean power over one period becomes
// 1 / 750.28, which is what the in-band SNR below refers to.
std::vector<cf> normalize_energy(std::vector<cf> sym)
{
    double e = 0.0;
    for (const cf& v : sym)
        e += std::norm(v);
    const double s = 1.0 / std::sqrt(e);
    for (cf& v : sym)
        v *= static_cast<float>(s);
    return sym;
}

std::vector<cf> native_symbol(const std::string& rel)
{
    const std::vector<cf> ref = load_cf32(rel);
    BOOST_REQUIRE_MESSAGE(ref.size() == 1016,
                          "unexpected reference length " << ref.size());
    return normalize_energy(resample_symbol_65_48(ref, kNativeSymLen));
}

// Overlap-add the SYNC symbol on the exact period grid.
std::vector<cf> build_preamble(const std::vector<cf>& sym, size_t reps)
{
    const size_t len = static_cast<size_t>(
                           std::llround(static_cast<double>(reps - 1) *
                                        kPeriod)) +
        sym.size();
    std::vector<cf> out(len, cf(0.0f, 0.0f));
    for (size_t r = 0; r < reps; ++r) {
        const size_t st =
            static_cast<size_t>(std::llround(static_cast<double>(r) * kPeriod));
        for (size_t i = 0; i < sym.size(); ++i)
            out[st + i] += sym[i];
    }
    return out;
}

std::vector<cf> embed(const std::vector<cf>& pre, size_t lead, size_t tail)
{
    std::vector<cf> x(lead + pre.size() + tail, cf(0.0f, 0.0f));
    std::copy(pre.begin(), pre.end(), x.begin() + lead);
    return x;
}

// In-band SNR (repo convention, see qa_uwb_cfo_ml): mean signal power over
// [begin, begin + count) sets the complex AWGN power; noise is added over the
// whole stream.
void add_awgn_inband(std::vector<cf>& x, size_t begin, size_t count,
                     double snr_db, uint32_t seed)
{
    double power = 0.0;
    for (size_t i = 0; i < count; ++i)
        power += std::norm(x[begin + i]);
    const double p_sig = power / static_cast<double>(count);
    const double sigma =
        std::sqrt(p_sig / std::pow(10.0, snr_db / 10.0) / 2.0);
    std::mt19937 rng(seed);
    std::normal_distribution<float> g(0.0f, static_cast<float>(sigma));
    for (cf& v : x)
        v += cf(g(rng), g(rng));
}

void mix_cfo(std::vector<cf>& x, double cfo_hz)
{
    if (cfo_hz == 0.0)
        return;
    const double w = 2.0 * kPi * cfo_hz / kFs;
    for (size_t n = 0; n < x.size(); ++n)
        x[n] *= cf(static_cast<float>(std::cos(w * static_cast<double>(n))),
                   static_cast<float>(std::sin(w * static_cast<double>(n))));
}

inline int16_t sat_sc16(float v)
{
    const long r = std::lround(v);
    if (r > 32767)
        return 32767;
    if (r < -32767)
        return -32767;
    return static_cast<int16_t>(r);
}

std::vector<sc16> quantize(const std::vector<cf>& x, double gain)
{
    std::vector<sc16> out(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        out[i] = sc16(sat_sc16(static_cast<float>(x[i].real() * gain)),
                      sat_sc16(static_cast<float>(x[i].imag() * gain)));
    return out;
}

// Copy the whole (short) test stream into a power-of-two ring so that the
// diagnostics entry can be called at any absolute window position.
std::vector<sc16> make_ring(const std::vector<sc16>& data, size_t ring_size)
{
    std::vector<sc16> ring(ring_size, sc16(0, 0));
    const size_t mask = ring_size - 1;
    BOOST_REQUIRE(data.size() <= ring_size);
    for (size_t i = 0; i < data.size(); ++i)
        ring[i & mask] = data[i];
    return ring;
}

struct WindowDiag {
    double par = 0.0;
    double code_rho = 0.0;
    double code_par_ratio = 0.0;
};

WindowDiag window_diag(PreambleXcorrAcquisition& acq,
                       const std::vector<sc16>& ring, int64_t win_start)
{
    PreambleXcorrDiagnostics diag;
    BOOST_REQUIRE(acq.evaluate_window(ring, 0, win_start, false, &diag));
    return WindowDiag{ diag.par, diag.code_rho, diag.code_par_ratio };
}

// ---------------------------------------------------------------------------
// Streaming drive loop, mirroring the openuwb_rx session contract:
// feed a block, advance end, valid_begin = max(0, end - ring.size()), then
// scan while [probe, probe + required) is inside the valid data.
// ---------------------------------------------------------------------------
struct ScanOutcome {
    std::vector<int64_t> seeds;
    std::vector<double> metrics;
    PreambleXcorrStats stats{};
};

ScanOutcome stream_scan(PreambleXcorrAcquisition& acq,
                        const std::vector<sc16>& data, size_t ring_size,
                        size_t block = kBlock)
{
    BOOST_REQUIRE_MESSAGE((ring_size & (ring_size - 1)) == 0,
                          "ring must be a power of two");
    ScanOutcome out;
    std::vector<sc16> ring(ring_size, sc16(0, 0));
    const size_t mask = ring_size - 1;
    const int64_t req = acq.required_window_samples();
    const int64_t total = static_cast<int64_t>(data.size());
    int64_t end = 0;
    while (end < total) {
        const size_t n = static_cast<size_t>(
            std::min<int64_t>(static_cast<int64_t>(block), total - end));
        for (size_t i = 0; i < n; ++i)
            ring[static_cast<size_t>(end + static_cast<int64_t>(i)) & mask] =
                data[static_cast<size_t>(end) + i];
        end += static_cast<int64_t>(n);
        const int64_t valid_begin =
            std::max<int64_t>(0, end - static_cast<int64_t>(ring_size));
        while (acq.probe() + req <= end) {
            const int64_t before = acq.probe();
            double metric = 0.0;
            const int64_t seed = acq.scan(ring, valid_begin, &metric);
            if (seed >= 0) {
                out.seeds.push_back(seed);
                out.metrics.push_back(metric);
            }
            if (acq.probe() <= before)
                break; // defensive: scan must always advance
        }
    }
    out.stats = acq.stats();
    return out;
}

// Ring wrap geometry: with a 65536-sample ring and probe = 150000, fold symbol
// 62 starts at 150000 + round(62 * 750.276923) = 196517, i.e. it crosses
// absolute sample 196608 (3 * 65536) inside its 750-sample span.
constexpr int64_t kWrapProbe = 150000;
constexpr int64_t kWrapFoldStart = 196517;
constexpr int64_t kWrapBoundary = 3 * 65536;

BOOST_AUTO_TEST_CASE(test_noise_and_tones_no_false_alarm)
{
    // >= 50 ms at 737.28 MS/s.
    const size_t total = static_cast<size_t>(0.05 * kFs);
    BOOST_TEST_MESSAGE("noise scan: " << total << " samples ("
                                      << (total / static_cast<size_t>(kHopSamples))
                                      << " ~hop windows)");

    struct Scenario {
        const char* name;
        std::vector<double> tones_hz;
        double tone_amp;
    };
    const std::vector<Scenario> cases = {
        { "awgn", {}, 0.0 },
        { "-232.96 MHz CW", { -232.96e6 }, 8000.0 },
        { "-124.4 MHz CW", { -124.4e6 }, 8000.0 },
        { "multi-tone", { -232.96e6, -124.4e6, 40e6 }, 5000.0 },
    };

    const float sigma = 1000.0f; // per component, 1 sigma
    for (const Scenario& sc : cases) {
        PreambleXcorrAcquisition acq;
        acq.reset(0);
        const int64_t req = acq.required_window_samples();

        std::vector<sc16> ring(kRingNoise, sc16(0, 0));
        const size_t mask = kRingNoise - 1;
        std::mt19937 rng(4242u + static_cast<uint32_t>(sc.tones_hz.size()));
        std::normal_distribution<float> gauss(0.0f, sigma);

        std::vector<cf> rot(sc.tones_hz.size(), cf(1.0f, 0.0f));
        std::vector<cf> step(sc.tones_hz.size(), cf(1.0f, 0.0f));
        for (size_t j = 0; j < sc.tones_hz.size(); ++j) {
            const double w = 2.0 * kPi * sc.tones_hz[j] / kFs;
            step[j] = cf(static_cast<float>(std::cos(w)),
                         static_cast<float>(std::sin(w)));
        }

        int64_t end = 0;
        while (end < static_cast<int64_t>(total)) {
            const size_t n = static_cast<size_t>(std::min<int64_t>(
                static_cast<int64_t>(kBlock),
                static_cast<int64_t>(total) - end));
            for (size_t i = 0; i < n; ++i) {
                float re = gauss(rng);
                float im = gauss(rng);
                for (size_t j = 0; j < rot.size(); ++j) {
                    re += static_cast<float>(sc.tone_amp * rot[j].real());
                    im += static_cast<float>(sc.tone_amp * rot[j].imag());
                    rot[j] *= step[j];
                }
                ring[static_cast<size_t>(end + static_cast<int64_t>(i)) & mask] =
                    sc16(sat_sc16(re), sat_sc16(im));
            }
            end += static_cast<int64_t>(n);
            const int64_t valid_begin =
                std::max<int64_t>(0, end - static_cast<int64_t>(kRingNoise));
            while (acq.probe() + req <= end) {
                const int64_t before = acq.probe();
                const int64_t seed = acq.scan(ring, valid_begin);
                BOOST_CHECK_MESSAGE(seed < 0, "false seed " << seed << " on "
                                                            << sc.name);
                if (acq.probe() <= before)
                    break;
            }
        }

        BOOST_TEST_MESSAGE(sc.name << ": windows=" << acq.stats().windows_scanned
                                   << " stage2=" << acq.stats().stage2_promoted
                                   << " seeds=" << acq.stats().seeds_emitted);
        BOOST_TEST(acq.stats().seeds_emitted == 0u);
        BOOST_TEST(acq.stats().windows_scanned > 300u);
    }
}

BOOST_AUTO_TEST_CASE(test_cfo_sensitivity)
{
    const std::vector<cf> sym = native_symbol(
        "reference_preamble_code10_998p4.cf32");
    const std::vector<cf> pre = build_preamble(sym, kPreambleReps);
    const size_t lead = 20000;
    const size_t tail = 80000;
    const std::vector<cf> base = embed(pre, lead, tail);
    const double gain = 4000.0 / symbol_peak(sym);

    const std::array<double, 5> cfos = { 0.0, 50000.0, -50000.0, 100000.0,
                                         -100000.0 };
    const std::array<double, 2> snrs = { -12.0, -14.0 };
    constexpr int kTrials = 10;

    for (double snr : snrs) {
        for (double cfo : cfos) {
            double par_min = 1e30;
            double par_sum = 0.0;
            int64_t err_max = 0;
            int detected = 0;
            std::array<double, kTrials> ratio{};
            for (int tr = 0; tr < kTrials; ++tr) {
                std::vector<cf> rx = base;
                add_awgn_inband(rx, lead, pre.size(), snr,
                                1000u + static_cast<uint32_t>(tr));
                mix_cfo(rx, cfo);
                const std::vector<sc16> data = quantize(rx, gain);

                PreambleXcorrAcquisition acq;
                acq.reset(0);
                const ScanOutcome out = stream_scan(acq, data, kRingSignal);

                BOOST_REQUIRE_MESSAGE(
                    out.seeds.size() == 1,
                    "snr=" << snr << " dB cfo=" << cfo / 1e3
                           << " kHz trial=" << tr
                           << " seeds=" << out.seeds.size());
                BOOST_CHECK_EQUAL(out.stats.seeds_emitted, 1u);
                BOOST_CHECK_GE(static_cast<int64_t>(out.stats.stage2_promoted),
                               1);
                BOOST_CHECK_EQUAL(out.stats.code_rejected, 0u);
                par_min = std::min(par_min, out.metrics[0]);
                par_sum += out.metrics[0];
                const int64_t err =
                    out.seeds[0] - static_cast<int64_t>(lead);
                err_max = std::max<int64_t>(
                    err_max, err < 0 ? -err : err);

                // The scan commits at the first probe inside the preamble.  The
                // default Stage-2 gate is the relative code PAR ratio
                // (code_par_ratio_min = 2.0); the absolute code_rho_min gate is
                // off by default and only reported.
                const std::vector<sc16> ring = make_ring(data, kRingSignal);
                const WindowDiag diag = window_diag(acq, ring, kHopSamples);
                BOOST_CHECK_GE(diag.code_par_ratio, 2.0);
                ratio[static_cast<size_t>(tr)] = diag.code_par_ratio;
                ++detected;
            }
            std::sort(ratio.begin(), ratio.end());
            BOOST_TEST_MESSAGE("snr=" << snr << " dB cfo=" << cfo / 1e3
                                      << " kHz: detected=" << detected << "/"
                                      << kTrials << " par[med/min]="
                                      << par_sum / kTrials << "/" << par_min
                                      << " code_par_ratio[med/min]="
                                      << ratio[kTrials / 2] << "/" << ratio[0]
                                      << " max|seed-true|=" << err_max);
            BOOST_TEST(detected == kTrials);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_front_edge_alignment)
{
    const std::vector<cf> sym = native_symbol(
        "reference_preamble_code10_998p4.cf32");
    const std::vector<cf> pre = build_preamble(sym, kPreambleReps);
    const double gain = 4000.0 / symbol_peak(sym);
    const double snr = -6.0;
    const size_t base_lead = 20000;
    // 0/25/50/75% of one 128-symbol hop.
    const std::array<int64_t, 4> phases = { 0, 24009, 48018, 72026 };

    for (int64_t off : phases) {
        const size_t lead = base_lead + static_cast<size_t>(off);
        const std::vector<cf> base = embed(pre, lead, 80000);
        int64_t err = 0;
        for (int tr = 0; tr < 2; ++tr) {
            std::vector<cf> rx = base;
            add_awgn_inband(rx, lead, pre.size(), snr,
                            2000u + static_cast<uint32_t>(tr));
            const std::vector<sc16> data = quantize(rx, gain);

            PreambleXcorrAcquisition acq;
            acq.reset(0);
            const ScanOutcome out = stream_scan(acq, data, kRingSignal);
            BOOST_REQUIRE_MESSAGE(out.seeds.size() == 1,
                                  "phase offset " << off << " trial " << tr
                                                  << " seeds="
                                                  << out.seeds.size());
            err = out.seeds[0] - static_cast<int64_t>(lead);
            BOOST_TEST_MESSAGE("hop phase +" << off << " samples: seed="
                                             << out.seeds[0] << " true=" << lead
                                             << " err=" << err);
            BOOST_TEST(std::llabs(err) <= kSeedTolerance);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_wrong_code_rejected)
{
    // Bandlimited reference-style code-9 preamble (same 48/65 resample path as
    // the target) with a 60 kHz CFO and a half-sample fractional delay (a
    // distinct transmitter never lands on the receiver's sample phase).  It is
    // injected at -8/0/+10 dB in-band SNR: PAR alone cannot reject it (an ideal
    // foreign Ipatov code reaches PAR > 6.5 around 0 dB), so the default Stage-2
    // relative code-identity gate (code_par_ratio_min = 2.0; the optional
    // absolute code_rho_min gate is 0/off) must carry the rejection and count
    // it in stats().code_rejected.
    std::vector<cf> sym = native_symbol("reference_preamble.bin");
    {
        // Half-sample delay via the same periodic windowed-sinc kernel.
        std::vector<cf> d(sym.size());
        constexpr int kHalf = 10;
        for (size_t i = 0; i < sym.size(); ++i) {
            const double u = static_cast<double>(i) + 0.5;
            const int64_t i0 = static_cast<int64_t>(std::floor(u));
            double re = 0.0;
            double im = 0.0;
            double wsum = 0.0;
            for (int k = -kHalf + 1; k <= kHalf; ++k) {
                const double dd = u - static_cast<double>(i0 + k);
                const double w = sinc(dd);
                const size_t idx = static_cast<size_t>(
                    ((i0 + k) % static_cast<int64_t>(sym.size()) +
                     static_cast<int64_t>(sym.size())) %
                    static_cast<int64_t>(sym.size()));
                re += static_cast<double>(sym[idx].real()) * w;
                im += static_cast<double>(sym[idx].imag()) * w;
                wsum += w;
            }
            d[i] = cf(static_cast<float>(re / wsum),
                      static_cast<float>(im / wsum));
        }
        sym = normalize_energy(d);
    }

    const std::vector<cf> pre = build_preamble(sym, kPreambleReps);
    const size_t lead = 20000;
    const std::vector<cf> base = embed(pre, lead, 80000);
    const double gain = 4000.0 / symbol_peak(sym);
    std::vector<sc16> data_0db;

    for (double snr : { -8.0, 0.0, 10.0 }) {
        std::vector<cf> rx = base;
        add_awgn_inband(rx, lead, pre.size(), snr, 777u);
        mix_cfo(rx, 60000.0);
        const std::vector<sc16> data = quantize(rx, gain);

        PreambleXcorrAcquisition acq; // default code_index = 10
        acq.reset(0);
        const ScanOutcome wrong = stream_scan(acq, data, kRingSignal);

        // Frontend's own full-window view of the interferer, three windows
        // spread across the preamble (fixed absolute positions).
        const std::vector<sc16> ring = make_ring(data, kRingSignal);
        double par_max = 0.0;
        double rho_max = 0.0;
        for (size_t k : { size_t(80), size_t(100), size_t(120) }) {
            const int64_t win = static_cast<int64_t>(lead) +
                static_cast<int64_t>(
                    std::llround(static_cast<double>(k) * kPeriod));
            const WindowDiag diag = window_diag(acq, ring, win);
            par_max = std::max(par_max, diag.par);
            rho_max = std::max(rho_max, diag.code_rho);
        }
        BOOST_TEST_MESSAGE("code-9 interferer at "
                           << snr << " dB: seeds=" << wrong.stats.seeds_emitted
                           << " stage2=" << wrong.stats.stage2_promoted
                           << " code_rejected=" << wrong.stats.code_rejected
                           << " window PAR max=" << par_max
                           << " window rho max=" << rho_max);
        BOOST_TEST(wrong.seeds.empty());
        BOOST_TEST(wrong.stats.seeds_emitted == 0u);
        // The foreign template is incoherent with the code-10 taps.
        BOOST_TEST(rho_max < 0.5);
        if (snr >= 0.0) {
            // PAR alone would have committed a seed here; the relative
            // code-identity gate is what rejected it.
            BOOST_TEST(wrong.stats.code_rejected >= 1u);
            if (snr == 0.0)
                data_0db = data;
        }

        // Positive control at the matching code on the -8 dB stream.
        if (snr == -8.0) {
            PreambleXcorrConfig cfg9;
            cfg9.code_index = 9;
            PreambleXcorrAcquisition acq9;
            acq9.configure(cfg9);
            acq9.reset(0);
            const ScanOutcome match = stream_scan(acq9, data, kRingSignal);
            BOOST_REQUIRE_MESSAGE(match.seeds.size() == 1,
                                  "code-9 control seeds=" << match.seeds.size());
            BOOST_TEST(match.metrics[0] > 6.5);
        }
    }

    // Positive control on the same 0 dB stream: matched code is detected and
    // the matched-template rho is high, so the rejection above is
    // code-specific, not energy-limited.
    BOOST_REQUIRE(!data_0db.empty());
    PreambleXcorrConfig cfg9;
    cfg9.code_index = 9;
    PreambleXcorrAcquisition acq9;
    acq9.configure(cfg9);
    acq9.reset(0);
    const ScanOutcome match = stream_scan(acq9, data_0db, kRingSignal);
    BOOST_REQUIRE_MESSAGE(match.seeds.size() == 1,
                          "code-9 0 dB control seeds=" << match.seeds.size());
    const std::vector<sc16> ring0 = make_ring(data_0db, kRingSignal);
    const WindowDiag match_diag = window_diag(acq9, ring0, kHopSamples);
    BOOST_TEST_MESSAGE("code-9 0 dB control: seed="
                       << match.seeds[0] << " true=" << lead
                       << " par=" << match.metrics[0]
                       << " code_rho=" << match_diag.code_rho);
    BOOST_TEST(match.metrics[0] > 6.5);
    BOOST_TEST(match_diag.code_rho > 0.8);
}

BOOST_AUTO_TEST_CASE(test_ring_wrap_holdoff_reset)
{
    const std::vector<cf> sym = native_symbol(
        "reference_preamble_code10_998p4.cf32");
    const std::vector<cf> pre = build_preamble(sym, kPreambleReps);
    const double gain = 4000.0 / symbol_peak(sym);

    // --- (a) fold crossing the 65536-sample ring boundary -------------------
    BOOST_REQUIRE(kWrapFoldStart < kWrapBoundary + 750);
    BOOST_REQUIRE((kWrapFoldStart & (kWrapBoundary - 1)) > 65536 - 750);
    {
        const size_t lead = 145000; // preamble start; probe 150000 is inside
        std::vector<cf> rx = embed(pre, lead, 20000);
        add_awgn_inband(rx, lead, pre.size(), 0.0, 11u);
        const std::vector<sc16> data = quantize(rx, gain);

        std::vector<sc16> ring(kRingWrap, sc16(0, 0));
        const size_t mask = kRingWrap - 1;
        const int64_t end = 210000;
        for (int64_t i = 0; i < end; ++i)
            ring[static_cast<size_t>(i) & mask] = data[static_cast<size_t>(i)];

        PreambleXcorrAcquisition acq;
        acq.reset(kWrapProbe);
        const int64_t valid_begin = end - static_cast<int64_t>(kRingWrap);
        double metric = 0.0;
        const int64_t seed = acq.scan(ring, valid_begin, &metric);
        BOOST_TEST_MESSAGE("wrap: seed=" << seed << " par=" << metric
                                         << " valid_begin=" << valid_begin);
        BOOST_TEST(seed >= 0);
        BOOST_CHECK_EQUAL(acq.stats().seeds_emitted, 1u);
    }

    // --- (b) holdoff suppression + (c) reset() clears the state -------------
    {
        const size_t lead_a = 20000;
        const size_t lead_b = 250000; // inside A's holdoff window
        std::vector<cf> rx(lead_b + pre.size() + 20000, cf(0.0f, 0.0f));
        for (size_t i = 0; i < pre.size(); ++i) {
            rx[lead_a + i] += pre[i];
            rx[lead_b + i] += pre[i];
        }
        add_awgn_inband(rx, lead_a, pre.size(), -6.0, 12u);
        const std::vector<sc16> data = quantize(rx, gain);

        PreambleXcorrAcquisition acq;
        acq.reset(0);
        const ScanOutcome out = stream_scan(acq, data, kRingSignal);
        BOOST_REQUIRE_MESSAGE(out.seeds.size() == 1,
                              "holdoff: seeds=" << out.seeds.size());
        const int64_t seed_a = out.seeds[0];
        BOOST_TEST_MESSAGE("holdoff: seed=" << seed_a
                                            << " probe=" << acq.probe()
                                            << " (second preamble at "
                                            << lead_b << " suppressed)");
        BOOST_CHECK_EQUAL(out.stats.seeds_emitted, 1u);
        BOOST_CHECK(acq.probe() >= seed_a + 589824);

        // (c) reset(start): stats/probe/holdoff cleared, reacquire from B.
        acq.reset(static_cast<int64_t>(lead_b));
        BOOST_CHECK_EQUAL(acq.stats().seeds_emitted, 0u);
        BOOST_CHECK_EQUAL(acq.stats().windows_scanned, 0u);
        BOOST_CHECK_EQUAL(acq.probe(), static_cast<int64_t>(lead_b));

        std::vector<sc16> ring(kRingSignal, sc16(0, 0));
        const size_t mask = kRingSignal - 1;
        for (size_t i = 0; i < data.size(); ++i)
            ring[i & mask] = data[i];
        const int64_t valid_begin =
            std::max<int64_t>(0, static_cast<int64_t>(data.size()) -
                                     static_cast<int64_t>(kRingSignal));
        BOOST_REQUIRE(acq.probe() + acq.required_window_samples() <=
                      static_cast<int64_t>(data.size()));
        double metric = 0.0;
        const int64_t seed_b = acq.scan(ring, valid_begin, &metric);
        BOOST_TEST_MESSAGE("reset: seed=" << seed_b << " true=" << lead_b
                                          << " err=" << (seed_b - static_cast<int64_t>(lead_b))
                                          << " par=" << metric);
        BOOST_TEST(seed_b >= 0);
        BOOST_TEST(std::llabs(seed_b - static_cast<int64_t>(lead_b)) <=
                   kSeedTolerance);
        BOOST_CHECK_EQUAL(acq.stats().seeds_emitted, 1u);
    }
}

BOOST_AUTO_TEST_CASE(test_avx2_scalar_consistency)
{
    const std::vector<cf> sym = native_symbol(
        "reference_preamble_code10_998p4.cf32");
    const std::vector<cf> pre = build_preamble(sym, kPreambleReps);
    const size_t lead = 20000;
    std::vector<cf> rx = embed(pre, lead, 40000);
    add_awgn_inband(rx, lead, pre.size(), -8.0, 999u);
    const double gain = 4000.0 / symbol_peak(sym);
    const std::vector<sc16> data = quantize(rx, gain);

    const size_t ring_size = 1u << 18;
    std::vector<sc16> ring(ring_size, sc16(0, 0));
    const size_t mask = ring_size - 1;
    const int64_t end = static_cast<int64_t>(data.size());
    for (int64_t i = 0; i < end; ++i)
        ring[static_cast<size_t>(i) & mask] = data[static_cast<size_t>(i)];
    const int64_t valid_begin = 0;
    const int64_t win = 40000; // inside the preamble, full history available

    PreambleXcorrAcquisition acq;
    PreambleXcorrDiagnostics fast;
    PreambleXcorrDiagnostics slow;
    BOOST_REQUIRE(acq.evaluate_window(ring, valid_begin, win, false, &fast));
    BOOST_REQUIRE(acq.evaluate_window(ring, valid_begin, win, true, &slow));

    auto rel_err = [](double a, double b) {
        const double d = std::max(std::fabs(a), std::fabs(b));
        return d > 0.0 ? std::fabs(a - b) / d : 0.0;
    };

    BOOST_TEST_MESSAGE("evaluate_window: par=" << fast.par
                                               << " peak=" << fast.peak_power
                                               << " code_rho=" << fast.code_rho
                                               << " tau=" << fast.peak_tau
                                               << " avx2="
                                               << PreambleXcorrAcquisition::avx2_available());
    BOOST_TEST(fast.peak_tau == slow.peak_tau);
    BOOST_TEST(rel_err(fast.par, slow.par) < 1e-5);
    BOOST_TEST(rel_err(fast.peak_power, slow.peak_power) < 1e-5);
    BOOST_TEST(rel_err(fast.code_rho, slow.code_rho) < 1e-5);

    // Both paths must also still see the preamble at this window, and the
    // matched-template rho must clear the default gate.
    BOOST_TEST(fast.par > 6.5);
    BOOST_TEST(slow.par > 6.5);
    BOOST_TEST(fast.code_rho > 0.5);
    BOOST_TEST(slow.code_rho > 0.5);

    if (!PreambleXcorrAcquisition::avx2_available()) {
        BOOST_TEST_MESSAGE(
            "AVX2 kernel not compiled in (scalar-debug): scalar self-check only");
        PreambleXcorrDiagnostics again;
        BOOST_REQUIRE(acq.evaluate_window(ring, valid_begin, win, true, &again));
        BOOST_TEST(again.par == slow.par);
        BOOST_TEST(again.code_rho == slow.code_rho);
        BOOST_TEST(again.peak_tau == slow.peak_tau);
    }
}

} // namespace
