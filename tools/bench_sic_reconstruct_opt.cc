/* -*- c++ -*- */
/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * SIC reconstruct/subtract benchmark (WP4).
 *
 * This tool is deliberately self-contained and deterministic: it builds every
 * fixture from the public HRP modulator and from fixed seeds, so it can be
 * compiled unmodified against both the pre-optimization baseline library and
 * the optimized library.  It only uses the legacy public entry points:
 *
 *   sic::prepare_sic_unified
 *   sic::sic_process_unified
 *   sic::reconstruct_and_subtract_unified
 *   sic::compute_roi_energies
 *   sic::estimate_discrete_cir / apply_discrete_cir / apply_discrete_cir_pulses
 *   mod::modulate_one / modulate_pulse_list / make_random_psdu
 *
 * The new *_ex / DiscreteCirApplyOptions entry points are intentionally NOT
 * used, so the same translation unit links against the baseline library.
 *
 * Scenarios:
 *   live-hit   n=506000, code 11, 128 SYNC, "decawave", 32-byte PSDU,
 *              20 kHz CFO, sparse + shaped-sparse, out-of-place + in-place.
 *   short-l1   n=135000, code 10, 64 SYNC, "4z2", PreambleOnly (L1).
 *   short-l2   n=135000, code 10, 64 SYNC, "4z2", FullPacket (L2).
 *   cfo        apply-layer sparse-pulse vs non-sparse discrete CIR for
 *              CFO in {0, 20k, 50k, 100k} Hz with a deterministic multi-lag
 *              model; reports suppression in dB for both paths.
 *   taps       apply-layer sparse-pulse path for tap_count in {16,20,24,32}.
 *   edge       model-level L1 reconstruction with q=2 (left lo=-2) and
 *              support_hi = n + 2000 (right edge); clip off/on.
 *
 * Every measured call is:
 *   - warmed up at least `--warmup` (default 10) iterations,
 *   - followed by `--reps` (default 30) measured iterations,
 *   - validated for finite output samples (hard failure otherwise).
 *
 * One machine-readable JSON line is printed per configuration, plus an
 * optional per-iteration raw JSONL file via `--raw <path>`.
 */

#include <openuwb/sic/sic_unified.h>
#include <openuwb/sic/discrete_cir_sic.h>
#include <openuwb/phy/hrp/tx/hrp_mod_core.h>
#include <openuwb/phy/hrp/common/profile.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace openuwb;
using cf = std::complex<float>;
using cd = std::complex<double>;
using Clock = std::chrono::steady_clock;

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 998.4e6;
constexpr size_t kSps = 1016;
constexpr double kFloor = 1e-4;   // deterministic DC floor for packet windows
constexpr double kCfoFloor = 1e-6; // near-zero floor for suppression headroom

struct Args {
    std::string scenario;
    bool inplace_only = false;
    double cfo_hz = -1.0;
    int taps = -1;
    int reps = 30;
    int warmup = 10;
    bool json = false;
    std::string raw_path;
};

[[noreturn]] void die(const std::string& msg)
{
    std::fprintf(stderr, "bench_sic_reconstruct_opt: %s\n", msg.c_str());
    std::exit(2);
}

double us_since(Clock::time_point t0)
{
    return std::chrono::duration<double, std::micro>(Clock::now() - t0).count();
}

cf phase_of(double omega, double rel)
{
    return cf(static_cast<float>(std::cos(omega * rel)),
              static_cast<float>(std::sin(omega * rel)));
}

// Nearest-rank percentile; p in [0,1].  Empty -> 0.
double percentile(std::vector<double> v, double p)
{
    if (v.empty())
        return 0.0;
    std::sort(v.begin(), v.end());
    if (p <= 0.0)
        return v.front();
    if (p >= 1.0)
        return v.back();
    const size_t idx = static_cast<size_t>(std::ceil(p * double(v.size()))) - 1;
    return v[std::min(idx, v.size() - 1)];
}

struct Series {
    std::vector<double> v;
    void add(double x) { v.push_back(x); }
    double p50() const { return percentile(v, 0.50); }
    double p95() const { return percentile(v, 0.95); }
    double last() const { return v.empty() ? 0.0 : v.back(); }
};

struct Measured {
    Series wall, core, recon, sub, roi_extra, back_half;
    Series tx_remod, zero, scale, conv, cfo_rot;
    double suppression_db = 0.0;
    bool committed = false;
    bool fused_energy = false;
    const char* level = "L0_Bypass";
    const char* reason = "not_run";
    const char* fallback = "";
    double cfo_hz = 0.0;
    size_t pulse_events = 0;
    double estimate_us = 0.0, comm_rx_us = 0.0;
    int64_t n = 0;
    int64_t q = -1;
    int64_t support_lo = 0;
    int64_t support_hi = 0;
    int64_t sub_lo = 0;
    int64_t sub_hi = 0;
    size_t wave_samples = 0;
    size_t tap_count = 0;
    int64_t lag_lo = 0;
    bool inplace = false;
    const char* roi_extra_kind = "none";
    std::string scenario;
    std::string mode;
    bool finite_ok = true;
    bool consistent = true;
};

const char* level_name(sic::SicCancellationLevel l)
{
    switch (l) {
    case sic::SicCancellationLevel::L1_Preamble:
        return "L1_Preamble";
    case sic::SicCancellationLevel::L2_FullPacket:
        return "L2_FullPacket";
    default:
        return "L0_Bypass";
    }
}

void emit_json(const Measured& m, int reps, int warmup, std::FILE* out)
{
    std::fprintf(
        out,
        "{\"scenario\":\"%s\",\"mode\":\"%s\",\"reps\":%d,\"warmup\":%d,"
        "\"n\":%lld,\"q\":%lld,\"wave_samples\":%zu,\"tap_count\":%zu,"
        "\"lag_lo\":%lld,\"span\":%lld,\"support_span\":%lld,"
        "\"cfo_hz\":%.3f,\"pulse_events\":%zu,"
        "\"committed\":%s,\"level\":\"%s\",\"reason\":\"%s\","
        "\"fallback_reason\":\"%s\",\"inplace\":%s,\"fused_energy\":%s,"
        "\"finite_ok\":%s,\"consistent\":%s,"
        "\"core_total_p50_us\":%.6f,\"core_total_p95_us\":%.6f,"
        "\"wall_p50_us\":%.6f,\"wall_p95_us\":%.6f,"
        "\"reconstruct_p50_us\":%.6f,\"reconstruct_p95_us\":%.6f,"
        "\"subtract_p50_us\":%.6f,\"subtract_p95_us\":%.6f,"
        "\"roi_energy_extra_us\":%.6f,\"roi_extra_kind\":\"%s\","
        "\"back_half_p50_us\":%.6f,\"back_half_p95_us\":%.6f,"
        "\"tx_remod_us\":%.6f,\"recon_zero_us\":%.6f,"
        "\"recon_scale_us\":%.6f,\"recon_conv_us\":%.6f,"
        "\"recon_cfo_rot_us\":%.6f,\"suppression_db\":%.6f,"
        "\"estimate_us\":%.6f,\"comm_rx_us\":%.6f}\n",
        m.scenario.c_str(), m.mode.c_str(), reps, warmup,
        (long long)m.n, (long long)m.q, m.wave_samples, m.tap_count,
        (long long)m.lag_lo, (long long)(m.sub_hi - m.sub_lo),
        (long long)(m.support_hi - m.support_lo),
        m.cfo_hz, m.pulse_events,
        m.committed ? "true" : "false", m.level, m.reason, m.fallback,
        m.inplace ? "true" : "false", m.fused_energy ? "true" : "false",
        m.finite_ok ? "true" : "false", m.consistent ? "true" : "false",
        m.core.p50(), m.core.p95(), m.wall.p50(), m.wall.p95(),
        m.recon.p50(), m.recon.p95(), m.sub.p50(), m.sub.p95(),
        m.roi_extra.last(), m.roi_extra_kind,
        m.back_half.p50(), m.back_half.p95(),
        m.tx_remod.last(), m.zero.last(), m.scale.last(), m.conv.last(),
        m.cfo_rot.last(), m.suppression_db, m.estimate_us, m.comm_rx_us);
    std::fflush(out);
}

void emit_text(const Measured& m, int reps, int warmup)
{
    std::printf(
        "scenario=%s mode=%s reps=%d warmup=%d n=%lld q=%lld wave=%zu "
        "taps=%zu lag_lo=%lld span=%lld committed=%d level=%s reason=%s "
        "core_p50=%.3f core_p95=%.3f recon=%.3f sub=%.3f roi_extra=%.3f "
        "back_half=%.3f zero=%.3f scale=%.3f conv=%.3f cfo=%.3f "
        "suppression=%.3f inplace=%d fused=%d\n",
        m.scenario.c_str(), m.mode.c_str(), reps, warmup,
        (long long)m.n, (long long)m.q, m.wave_samples, m.tap_count,
        (long long)m.lag_lo, (long long)(m.sub_hi - m.sub_lo),
        m.committed ? 1 : 0, m.level, m.reason, m.core.p50(), m.core.p95(),
        m.recon.p50(), m.sub.p50(), m.roi_extra.last(), m.back_half.p50(),
        m.zero.last(), m.scale.last(), m.conv.last(), m.cfo_rot.last(),
        m.suppression_db, m.inplace ? 1 : 0, m.fused_energy ? 1 : 0);
    std::fflush(stdout);
}

// ---------------------------------------------------------------------------
// Deterministic packet fixtures (public modulator only).
// ---------------------------------------------------------------------------

struct Packet {
    mod::HrpModConfig profile;
    mod::HrpModScratch scratch;
    std::vector<uint8_t> wire;
    std::vector<cf> wave;
    std::vector<cf> sync;

    Packet(size_t code, size_t reps, const char* sfd, size_t wire_bytes,
           uint32_t seed, size_t cap)
    {
        profile.code_index = code;
        profile.sync_repetitions = reps;
        profile.sfd_mode = sfd;
        profile.peak_amplitude = 0.8f;
        wire = mod::make_random_psdu(wire_bytes - 2, seed, true);
        scratch.reserve(127, cap);
        wave.resize(cap);
        size_t produced = 0;
        if (!mod::modulate_one(wire.data(), wire.size(), profile, scratch,
                               wave.data(), cap, produced))
            die("modulate_one failed");
        wave.resize(produced);
        if (wave.size() < kSps)
            die("modulated packet shorter than one SYNC");
        sync.assign(wave.begin(), wave.begin() + kSps);
    }

    sic::SicUnifiedTemplates templates() const
    {
        sic::SicUnifiedTemplates t;
        t.sync = sync.data();
        t.sync_len = sync.size();
        t.preamble = wave.data();
        t.preamble_len = wave.size();
        return t;
    }
};

std::vector<cf> make_packet_window(size_t n, int64_t q,
                                   const std::vector<cf>& wave, double cfo_hz,
                                   double floor_amp, uint32_t seed)
{
    std::vector<cf> y(n, cf(static_cast<float>(floor_amp),
                            static_cast<float>(-floor_amp)));
    // A tiny deterministic dither keeps the floor from being a pure DC term.
    uint32_t s = seed;
    auto next = [&s, floor_amp]() {
        s = s * 1664525u + 1013904223u;
        return (static_cast<float>((s >> 8) & 0xffffu) / 65536.f - 0.5f) *
               static_cast<float>(floor_amp * 0.1);
    };
    for (auto& v : y)
        v += cf(next(), next());
    const double omega = 2.0 * kPi * cfo_hz / kFs;
    const int64_t end = std::min<int64_t>(q + int64_t(wave.size()),
                                          int64_t(n));
    for (int64_t j = std::max<int64_t>(q, 0); j < end; ++j) {
        const double rel = double(j - q);
        y[static_cast<size_t>(j)] +=
            wave[static_cast<size_t>(j - q)] * phase_of(omega, rel);
    }
    return y;
}

bool window_is_finite(const cf* p, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (!std::isfinite(p[i].real()) || !std::isfinite(p[i].imag()))
            return false;
    }
    return true;
}

double suppression_from(const std::vector<cf>& y, const std::vector<cf>& out,
                        int64_t lo, int64_t hi)
{
    if (hi <= lo || lo < 0 || hi > int64_t(y.size()))
        return 0.0;
    double e_y = 0.0, e_out = 0.0;
    sic::compute_roi_energies(y.data(), out.data(), lo, hi, e_y, e_out);
    if (e_y > 0.0 && e_out > 0.0)
        return 10.0 * std::log10(e_y / e_out);
    return 0.0;
}

// ---------------------------------------------------------------------------
// Unified-pipeline measurement (live-hit / short / edge).
// ---------------------------------------------------------------------------

struct PipelineCase {
    std::string scenario;
    std::string mode;
    size_t n = 0;
    size_t wave_samples = 0;
    size_t tap_count = 0;
    int64_t lag_lo = 0;
    bool inplace = false;
    std::vector<cf> y; // pristine input
    std::vector<cf> out;
    cf* in_ptr = nullptr; // per-iteration input base (out-of-place)
};

Measured run_pipeline_case(const Args& args, sic::SicUnifiedConfig cfg,
                           const sic::SicUnifiedTemplates& t,
                           sic::SicUnifiedScratch& scratch, PipelineCase& c)
{
    std::string err;
    if (!sic::prepare_sic_unified(cfg, t, scratch, &err))
        die("prepare_sic_unified failed for " + c.mode + ": " + err);

    Measured m;
    m.scenario = c.scenario;
    m.mode = c.mode;
    m.n = int64_t(c.n);
    m.q = cfg.roi_lo + (cfg.roi_hi - cfg.roi_lo) / 2;
    m.wave_samples = c.wave_samples;
    m.tap_count = c.tap_count;
    m.lag_lo = c.lag_lo;
    m.inplace = c.inplace;

    std::vector<cf> work;
    if (c.inplace)
        work.assign(c.n, cf{});

    const int warmup = std::max(0, args.warmup);
    const int reps = std::max(1, args.reps);
    bool first = true;

    for (int it = 0; it < warmup + reps; ++it) {
        const bool measure = (it >= warmup);
        cf* in = nullptr;
        cf* out = nullptr;
        if (c.inplace) {
            std::memcpy(work.data(), c.y.data(), c.n * sizeof(cf));
            in = work.data();
            out = work.data();
        } else {
            in = c.y.data();
            std::fill(c.out.begin(), c.out.end(), cf{});
            out = c.out.data();
        }

        sic::SicUnifiedResult r;
        r.out = out;
        r.i_hat = nullptr;
        r.copy_out_on_bypass = !c.inplace;

        const auto t0 = Clock::now();
        const bool ok =
            sic::sic_process_unified(in, c.n, t, cfg, scratch, r);
        const double wall = us_since(t0);
        if (!ok)
            die("sic_process_unified returned false for " + c.mode);

        const cf* outp = c.inplace ? work.data() : c.out.data();
        const bool finite = window_is_finite(outp, c.n);
        if (!finite) {
            std::fprintf(stderr,
                         "FATAL: non-finite output sample in %s/%s iter=%d\n",
                         c.scenario.c_str(), c.mode.c_str(), it);
            std::exit(4);
        }

        if (first) {
            m.committed = r.committed;
            m.level = level_name(r.level);
            m.reason = r.reason ? r.reason : "null";
            m.fallback = r.fallback_reason ? r.fallback_reason : "";
            m.support_lo = r.support_lo;
            m.support_hi = r.support_hi;
            m.sub_lo = r.sub_lo;
            m.sub_hi = r.sub_hi;
            m.cfo_hz = r.channel.cfo_hz;
            m.estimate_us = r.timings.estimate_us;
            m.comm_rx_us = r.timings.comm_rx_us;
            m.suppression_db = r.suppression_db;
            m.fused_energy = r.suppression_energy_valid;
            if (r.committed && r.sub_hi > r.sub_lo) {
                m.roi_extra_kind = r.suppression_energy_valid ? "fused_zero"
                                      : (c.inplace ? "post_read_estimate"
                                                   : "standalone");
            }
            first = false;
        } else {
            m.consistent = m.consistent &&
                           (r.committed == m.committed) &&
                           (std::strcmp(level_name(r.level), m.level) == 0);
        }

        if (!measure)
            continue;

        double roi_extra = 0.0;
        if (!r.committed || r.sub_hi <= r.sub_lo) {
            m.suppression_db = r.suppression_db;
            roi_extra = 0.0;
        } else if (r.suppression_energy_valid) {
            m.suppression_db = r.suppression_db;
            roi_extra = 0.0;
        } else if (c.inplace) {
            // Baseline in-place: the pipeline re-reads y/out after the write,
            // so its own metric is 0 dB.  Time the same two-read pass to keep
            // the back-half accounting honest; flag it in the report.
            const auto e0 = Clock::now();
            double e_y = 0.0, e_out = 0.0;
            sic::compute_roi_energies(work.data(), work.data(), r.sub_lo,
                                      r.sub_hi, e_y, e_out);
            roi_extra = us_since(e0);
            m.suppression_db = r.suppression_db;
        } else {
            const auto e0 = Clock::now();
            double e_y = 0.0, e_out = 0.0;
            sic::compute_roi_energies(c.y.data(), c.out.data(), r.sub_lo,
                                      r.sub_hi, e_y, e_out);
            roi_extra = us_since(e0);
            m.suppression_db = r.suppression_db;
        }

        m.wall.add(wall);
        m.core.add(r.timings.core_total_us);
        m.recon.add(r.timings.reconstruct_us);
        m.sub.add(r.timings.subtract_us);
        m.roi_extra.add(roi_extra);
        m.back_half.add(r.timings.reconstruct_us + r.timings.subtract_us +
                        roi_extra);
        m.tx_remod.add(r.timings.tx_remod_us);
        m.zero.add(r.timings.recon_zero_us);
        m.scale.add(r.timings.recon_scale_us);
        m.conv.add(r.timings.recon_conv_us);
        m.cfo_rot.add(r.timings.recon_cfo_rot_us);
        m.finite_ok = m.finite_ok && finite;

        if (!args.raw_path.empty()) {
            std::ofstream raw(args.raw_path, std::ios::app);
            raw << "{\"scenario\":\"" << c.scenario << "\",\"mode\":\""
                << c.mode << "\",\"iter\":" << it
                << ",\"inplace\":" << (c.inplace ? "true" : "false")
                << ",\"committed\":" << (r.committed ? "true" : "false")
                << ",\"wall_us\":" << wall
                << ",\"core_total_us\":" << r.timings.core_total_us
                << ",\"reconstruct_us\":" << r.timings.reconstruct_us
                << ",\"subtract_us\":" << r.timings.subtract_us
                << ",\"recon_zero_us\":" << r.timings.recon_zero_us
                << ",\"recon_scale_us\":" << r.timings.recon_scale_us
                << ",\"recon_conv_us\":" << r.timings.recon_conv_us
                << ",\"recon_cfo_rot_us\":" << r.timings.recon_cfo_rot_us
                << ",\"roi_extra_us\":" << roi_extra
                << ",\"back_half_us\":"
                << (r.timings.reconstruct_us + r.timings.subtract_us +
                    roi_extra)
                << ",\"suppression_db\":" << m.suppression_db << "}\n";
        }
    }
    return m;
}

int run_live_hit(const Args& args)
{
    const size_t n = 506000;
    const size_t q = 1024;
    const double cfo = 20000.0;
    Packet pkt(11, 128, "decawave", 32, 0xC0FFEEu, 400000);
    const auto t = pkt.templates();
    const std::vector<cf> y0 =
        make_packet_window(n, int64_t(q), pkt.wave, cfo, kFloor, 0x1234u);

    struct Mode {
        const char* name;
        bool shaped;
        bool inplace;
    };
    const Mode modes[] = {
        { "sparse_outplace", false, false },
        { "sparse_inplace", false, true },
        { "shaped_sparse_outplace", true, false },
        { "shaped_sparse_inplace", true, true },
    };

    sic::SicUnifiedScratch scratch;
    for (const auto& mode : modes) {
        if (args.inplace_only && !mode.inplace)
            continue;
        sic::SicUnifiedConfig cfg;
        cfg.policy = sic::SicPolicy::AutoFallback;
        cfg.backend = sic::SicEstimatorBackend::Ipatov;
        cfg.sparse_pulse_reconstruct = true;
        cfg.comm_code = 11;
        cfg.sync_reps = 128;
        cfg.samples_per_symbol = kSps;
        cfg.sample_rate = kFs;
        cfg.sfd_mode = "decawave";
        cfg.roi_lo = int64_t(q) - 32;
        cfg.roi_hi = int64_t(q) + 32;
        cfg.sub_roi_lo = 0;
        cfg.sub_roi_hi = int64_t(n);
        cfg.max_window_samples = n;
        cfg.max_packet_samples = n;
        cfg.rx_timing_shaped_sparse = mode.shaped;

        PipelineCase c;
        c.scenario = "live-hit";
        c.mode = mode.name;
        c.n = n;
        c.wave_samples = pkt.wave.size();
        c.tap_count = cfg.tap_count;
        c.lag_lo = cfg.lag_lo;
        c.inplace = mode.inplace;
        c.y = y0;
        c.out.assign(n, cf{});
        Measured m = run_pipeline_case(args, cfg, t, scratch, c);
        if (args.json)
            emit_json(m, args.reps, args.warmup, stdout);
        else
            emit_text(m, args.reps, args.warmup);
    }
    return 0;
}

int run_short(const Args& args, bool l2)
{
    const size_t n = 135000;
    const size_t q = 512;
    Packet pkt(10, 64, "4z2", 32, 0x5EEDu, 250000);
    const auto t = pkt.templates();
    const std::vector<cf> y0 =
        make_packet_window(n, int64_t(q), pkt.wave, 0.0, kFloor, 0x77u);

    sic::SicUnifiedConfig cfg;
    cfg.policy = l2 ? sic::SicPolicy::FullPacket
                    : sic::SicPolicy::PreambleOnly;
    cfg.backend = sic::SicEstimatorBackend::Ipatov;
    cfg.sparse_pulse_reconstruct = true;
    cfg.comm_code = 10;
    cfg.sync_reps = 64;
    cfg.samples_per_symbol = kSps;
    cfg.sample_rate = kFs;
    cfg.sfd_mode = "4z2";
    cfg.roi_lo = int64_t(q) - 32;
    cfg.roi_hi = int64_t(q) + 32;
    cfg.sub_roi_lo = 0;
    cfg.sub_roi_hi = int64_t(n);
    cfg.max_window_samples = n;
    cfg.max_packet_samples = n;

    PipelineCase c;
    c.scenario = l2 ? "short-l2" : "short-l1";
    c.mode = l2 ? "l2_fullpacket" : "l1_preamble";
    c.n = n;
    c.wave_samples = pkt.wave.size();
    c.tap_count = cfg.tap_count;
    c.lag_lo = cfg.lag_lo;
    c.inplace = false;
    c.y = y0;
    c.out.assign(n, cf{});

    sic::SicUnifiedScratch scratch;
    Measured m = run_pipeline_case(args, cfg, t, scratch, c);
    if (args.json)
        emit_json(m, args.reps, args.warmup, stdout);
    else
        emit_text(m, args.reps, args.warmup);
    return 0;
}

// ---------------------------------------------------------------------------
// Apply-layer rig (sparse pulses vs non-sparse full waveform).
// ---------------------------------------------------------------------------

struct SparseRig {
    size_t n = 0;
    int64_t q = 2048;
    std::vector<uint8_t> wire;
    mod::HrpModConfig unit;
    mod::HrpModScratch ms;
    mod::PulseList pl{};
    std::vector<mod::PulseEvent> events;
    std::vector<cf> unit_wave;
    std::vector<cf> impulse_sync;
    sic::DiscreteCirScratch scratch;

    explicit SparseRig(size_t wire_len, uint32_t seed, size_t extra)
    {
        wire = mod::make_random_psdu(wire_len, seed, true);
        unit.code_index = 10;
        unit.sync_repetitions = 64;
        unit.sfd_mode = "4z2";
        unit.peak_amplitude = 1.f;
        unit.pulse.shape = mod::PulseShape::External;
        unit.pulse.external_n_taps = 2;
        unit.external_taps = { 1.f, 0.f };
        const size_t cap = 400000;
        ms.reserve(127, cap);
        if (!mod::modulate_pulse_list(wire.data(), wire.size(), unit, ms, cap, pl))
            die("modulate_pulse_list failed");
        events.assign(pl.events, pl.events + pl.count);
        pl.events = events.data();
        n = pl.packet_samples + extra;
        unit_wave.assign(pl.packet_samples, cf{});
        for (const auto& e : events)
            unit_wave[e.sample_offset] = cf(float(e.sign), 0.f);

        scratch.reserve(n);
        scratch.prepared = true;
        scratch.sps = kSps;
        scratch.ipatov_correlation = true;
        scratch.sparse_pulse_reconstruct = true;
        scratch.impulse_sync.assign(kSps, cf{});
        const int8_t* code = demod::GetPreambleCode(10);
        if (!code)
            die("GetPreambleCode(10) failed");
        for (size_t c = 0; c < demod::kQm35CodeLength; ++c)
            scratch.impulse_sync[8 * c] = cf(float(code[c]), 0.f);
        scratch.max_packet_samples = pl.packet_samples + 64;
    }

    sic::DiscreteCirModel model(int64_t origin, int64_t lag_lo, size_t taps,
                                double cfo_hz) const
    {
        sic::DiscreteCirModel m;
        m.valid = true;
        m.status = "model_ready";
        m.model_id = 0xB0B0B0B0ull;
        m.q = origin;
        m.phase_reference_sample = origin;
        m.observed_end_sample = int64_t(n);
        m.cfo_hz = cfo_hz;
        m.lag_lo = lag_lo;
        m.tap_count = taps;
        m.fit_residual = 0.01;
        m.gain_condition = 12.0;
        m.effective_rank = taps;
        m.support_source = "bench_synthetic";
        m.impulse_template = true;
        return m;
    }

    static void set_active(sic::DiscreteCirModel& m,
                           const std::vector<std::pair<int64_t, cf>>& entries)
    {
        m.active_count = entries.size();
        for (size_t k = 0; k < entries.size(); ++k) {
            const int64_t lag = entries[k].first;
            const int64_t idx = lag - m.lag_lo;
            if (idx < 0 || idx >= int64_t(m.tap_count))
                die("active lag outside hull");
            m.active_lags[k] = static_cast<int16_t>(lag);
            m.taps[static_cast<size_t>(idx)] = entries[k].second;
        }
    }

    // Exact observation: y[q + m + l] += sign * h[l] with the CFO phase
    // exp(j*w*(m+l)) anchored at q, i.e. the discrete model exactly.
    std::vector<cf> synthesize(const sic::DiscreteCirModel& m,
                               double floor_amp, uint32_t seed) const
    {
        std::vector<cf> y(n, cf(static_cast<float>(floor_amp),
                                static_cast<float>(-floor_amp)));
        uint32_t s = seed;
        auto next = [&s, floor_amp]() {
            s = s * 1664525u + 1013904223u;
            return (static_cast<float>((s >> 8) & 0xffffu) / 65536.f - 0.5f) *
                   static_cast<float>(floor_amp * 0.1);
        };
        for (auto& v : y)
            v += cf(next(), next());
        const double omega = 2.0 * kPi * m.cfo_hz / kFs;
        for (const auto& e : events) {
            const double m_rel = double(e.sample_offset);
            for (size_t d = 0; d < m.tap_count; ++d) {
                const cf h = m.taps[d];
                if (h == cf{})
                    continue;
                const int64_t lag = m.lag_lo + int64_t(d);
                const int64_t j = m.q + int64_t(e.sample_offset) + lag;
                if (j < 0 || j >= int64_t(n))
                    continue;
                y[static_cast<size_t>(j)] +=
                    phase_of(omega, m_rel + double(lag)) *
                    (float(e.sign) * h);
            }
        }
        return y;
    }
};

Measured run_cfo_case(const Args& args, SparseRig& rig, double cfo_hz)
{
    sic::DiscreteCirModel m = rig.model(rig.q, -4, 20, cfo_hz);
    // Deterministic multi-lag CIR with the energy at large lags, so the
    // legacy intra-tap CFO approximation (phase from the event origin only)
    // is visibly wrong at 100 kHz while the exact per-lag path is not.
    SparseRig::set_active(m,
                          { { -4, cf(0.30f, 0.0f) },
                            { 0, cf(0.10f, 0.0f) },
                            { 10, cf(0.30f, 0.0f) },
                            { 15, cf(0.90f, 0.0f) } });
    const int64_t lo = rig.q + m.lag_lo;
    const int64_t hi = lo + int64_t(m.tap_count) + int64_t(rig.pl.packet_samples);
    const std::vector<cf> y = rig.synthesize(m, kCfoFloor, 0x9911u);

    Measured out_m, ns_m;
    for (int pass = 0; pass < 2; ++pass) {
        const bool sparse = (pass == 0);
        Measured& mm = sparse ? out_m : ns_m;
        mm.scenario = "cfo";
        mm.mode = sparse ? "sparse_pulse" : "nonsparse_discrete";
        mm.n = int64_t(rig.n);
        mm.q = rig.q;
        mm.cfo_hz = cfo_hz;
        mm.pulse_events = rig.pl.count;
        mm.wave_samples = rig.pl.packet_samples;
        mm.tap_count = m.tap_count;
        mm.lag_lo = m.lag_lo;
        mm.support_lo = lo;
        mm.support_hi = hi;
        mm.sub_lo = lo;
        mm.sub_hi = hi;
        std::vector<cf> o(rig.n), ihat(rig.n);

        bool first = true;
        for (int it = 0; it < args.warmup + args.reps; ++it) {
            const bool measure = (it >= args.warmup);
            std::fill(o.begin(), o.end(), cf{});
            std::fill(ihat.begin(), ihat.end(), cf{});
            sic::DiscreteCirResult r;
            const auto t0 = Clock::now();
            const bool ok =
                sparse
                    ? sic::apply_discrete_cir_pulses(
                          y.data(), rig.n, rig.pl,
                          rig.scratch.impulse_sync.data(),
                          rig.scratch.impulse_sync.size(), m, 0,
                          int64_t(rig.n), rig.scratch, o.data(), ihat.data(),
                          &r)
                    : sic::apply_discrete_cir(
                          y.data(), rig.n, rig.unit_wave.data(),
                          rig.unit_wave.size(), rig.scratch.impulse_sync.data(),
                          rig.scratch.impulse_sync.size(), m, 0,
                          int64_t(rig.n), rig.scratch, o.data(), ihat.data(),
                          &r);
            const double wall = us_since(t0);
            if (!ok || !r.committed)
                die(std::string("apply failed: ") +
                    (sparse ? "sparse_pulse" : "nonsparse_discrete") +
                    " reason=" + (r.reason ? r.reason : "null"));
            if (!window_is_finite(o.data(), rig.n))
                die("non-finite apply output");

            const auto e0 = Clock::now();
            double e_y = 0.0, e_out = 0.0;
            sic::compute_roi_energies(y.data(), o.data(), lo, hi, e_y, e_out);
            const double roi_us = us_since(e0);
            const double supp =
                (e_y > 0.0 && e_out > 0.0)
                    ? 10.0 * std::log10(e_y / e_out)
                    : 0.0;

            if (first) {
                mm.committed = true;
                mm.level = "apply_layer";
                mm.reason = r.reason ? r.reason : "committed";
                mm.suppression_db = supp;
                mm.roi_extra_kind = "standalone";
                first = false;
            }
            if (!measure)
                continue;
            mm.wall.add(wall);
            mm.core.add(r.apply_total_us);
            mm.recon.add(r.channel_us);
            mm.sub.add(r.subtract_us);
            mm.roi_extra.add(roi_us);
            mm.back_half.add(r.channel_us + r.subtract_us + roi_us);
            mm.zero.add(r.zero_us);
            mm.scale.add(r.scale_us);
            mm.conv.add(r.convolution_us);
            mm.cfo_rot.add(r.cfo_us + r.phase_us);
            mm.suppression_db = supp;
            if (!args.raw_path.empty()) {
                std::ofstream raw(args.raw_path, std::ios::app);
                raw << "{\"scenario\":\"cfo\",\"mode\":\""
                    << (sparse ? "sparse_pulse" : "nonsparse_discrete")
                    << "\",\"cfo_hz\":" << cfo_hz << ",\"iter\":" << it
                    << ",\"wall_us\":" << wall
                    << ",\"apply_total_us\":" << r.apply_total_us
                    << ",\"channel_us\":" << r.channel_us
                    << ",\"subtract_us\":" << r.subtract_us
                    << ",\"roi_extra_us\":" << roi_us
                    << ",\"suppression_db\":" << supp << "}\n";
            }
        }
    }

    // Emit one line per path, with the CFO encoded in the mode label and the
    // suppression value from the stable measured pass.
    for (Measured* p : { &out_m, &ns_m }) {
        std::string mode = p->mode + "@" + std::to_string(int(std::lround(cfo_hz))) + "Hz";
        p->mode = mode;
        if (args.json)
            emit_json(*p, args.reps, args.warmup, stdout);
        else
            emit_text(*p, args.reps, args.warmup);
    }
    return out_m; // sparse result (for the tap-scaling helper style reuse)
}

int run_cfo(const Args& args)
{
    SparseRig rig(30, 909u, 4096);
    std::vector<double> cfos = { 0.0, 20000.0, 50000.0, 100000.0 };
    if (args.cfo_hz >= 0.0)
        cfos = { args.cfo_hz };
    for (double c : cfos)
        run_cfo_case(args, rig, c);
    return 0;
}

int run_taps(const Args& args)
{
    SparseRig rig(30, 909u, 4096);
    std::vector<size_t> counts = { 16, 20, 24, 32 };
    if (args.taps > 0)
        counts = { static_cast<size_t>(args.taps) };
    for (size_t taps : counts) {
        if (taps < 8 || taps >= kSps)
            die("tap_count outside [8,1016)");
        sic::DiscreteCirModel m = rig.model(rig.q, -4, taps, 0.0);
        // Active lags that fit every hull.
        std::vector<std::pair<int64_t, cf>> active = {
            { -4, cf(0.35f, -0.10f) },
            { 0, cf(0.90f, 0.20f) },
            { 4, cf(-0.25f, 0.15f) },
            { 8, cf(0.15f, 0.05f) },
        };
        if (taps > 16)
            active.push_back({ int64_t(taps) - 5, cf(0.10f, -0.05f) });
        SparseRig::set_active(m, active);
        const std::vector<cf> y = rig.synthesize(m, kCfoFloor, 0x9911u);
        const int64_t lo = rig.q + m.lag_lo;
        const int64_t hi =
            lo + int64_t(m.tap_count) + int64_t(rig.pl.packet_samples);

        std::vector<cf> o(rig.n), ihat(rig.n);
        Measured mm;
        mm.scenario = "taps";
        mm.mode = "sparse_pulse_t" + std::to_string(taps);
        mm.n = int64_t(rig.n);
        mm.q = rig.q;
        mm.pulse_events = rig.pl.count;
        mm.wave_samples = rig.pl.packet_samples;
        mm.tap_count = taps;
        mm.lag_lo = m.lag_lo;
        mm.support_lo = lo;
        mm.support_hi = hi;
        mm.sub_lo = lo;
        mm.sub_hi = hi;
        mm.roi_extra_kind = "standalone";
        bool first = true;
        for (int it = 0; it < args.warmup + args.reps; ++it) {
            const bool measure = (it >= args.warmup);
            std::fill(o.begin(), o.end(), cf{});
            // i_hat == nullptr exercises the live streaming (L1D tiled) path.
            sic::DiscreteCirResult r;
            const auto t0 = Clock::now();
            const bool ok = sic::apply_discrete_cir_pulses(
                y.data(), rig.n, rig.pl, rig.scratch.impulse_sync.data(),
                rig.scratch.impulse_sync.size(), m, 0, int64_t(rig.n),
                rig.scratch, o.data(), nullptr, &r);
            const double wall = us_since(t0);
            if (!ok || !r.committed)
                die("taps apply failed for " + std::to_string(taps));
            if (!window_is_finite(o.data(), rig.n))
                die("non-finite taps output");
            const auto e0 = Clock::now();
            double e_y = 0.0, e_out = 0.0;
            sic::compute_roi_energies(y.data(), o.data(), lo, hi, e_y, e_out);
            const double roi_us = us_since(e0);
            if (first) {
                mm.committed = true;
                mm.level = "apply_layer";
                mm.reason = r.reason ? r.reason : "committed";
                mm.suppression_db =
                    (e_y > 0.0 && e_out > 0.0)
                        ? 10.0 * std::log10(e_y / e_out)
                        : 0.0;
                first = false;
            }
            if (!measure)
                continue;
            mm.wall.add(wall);
            mm.core.add(r.apply_total_us);
            mm.recon.add(r.channel_us);
            mm.sub.add(r.subtract_us);
            mm.roi_extra.add(roi_us);
            mm.back_half.add(r.channel_us + r.subtract_us + roi_us);
            mm.zero.add(r.zero_us);
            mm.scale.add(r.scale_us);
            mm.conv.add(r.convolution_us);
            mm.cfo_rot.add(r.cfo_us + r.phase_us);
            if (!args.raw_path.empty()) {
                std::ofstream raw(args.raw_path, std::ios::app);
                raw << "{\"scenario\":\"taps\",\"mode\":\"" << mm.mode
                    << "\",\"taps\":" << taps << ",\"iter\":" << it
                    << ",\"wall_us\":" << wall
                    << ",\"apply_total_us\":" << r.apply_total_us
                    << ",\"channel_us\":" << r.channel_us
                    << ",\"subtract_us\":" << r.subtract_us
                    << ",\"roi_extra_us\":" << roi_us
                    << ",\"suppression_db\":" << mm.suppression_db << "}\n";
            }
        }
        if (args.json)
            emit_json(mm, args.reps, args.warmup, stdout);
        else
            emit_text(mm, args.reps, args.warmup);
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Edge clipping: model-level L1 reconstruction through the public unified
// reconstruct entry (no demodulation dependency).
// ---------------------------------------------------------------------------

Measured run_edge_case(const Args& args, bool right_edge, bool clip,
                       sic::SicUnifiedScratch& scratch,
                       const sic::SicUnifiedTemplates& t)
{
    // Prepare the Ipatov pulse caches for code 10 / 64 SYNC / 4z2.
    sic::SicUnifiedConfig cfg;
    cfg.policy = sic::SicPolicy::PreambleOnly;
    cfg.backend = sic::SicEstimatorBackend::Ipatov;
    cfg.sparse_pulse_reconstruct = true;
    cfg.comm_code = 10;
    cfg.sync_reps = 64;
    cfg.samples_per_symbol = kSps;
    cfg.sample_rate = kFs;
    cfg.sfd_mode = "4z2";
    cfg.roi_lo = 0;
    cfg.roi_hi = 64;
    cfg.sub_roi_lo = 0;
    cfg.max_window_samples = 200000;
    cfg.max_packet_samples = 200000;
    cfg.allow_partial_window_clip = clip;

    std::string err;
    if (!sic::prepare_sic_unified(cfg, t, scratch, &err))
        die("edge prepare failed: " + err);
    if (scratch.cached_preamble_events.empty())
        die("edge prepare did not cache Ipatov pulse events");

    const size_t cached = scratch.cached_preamble_samples;
    const int64_t q = right_edge ? 512 : 2;
    const int64_t lag_lo = -4;
    const size_t taps = 20;
    const int64_t support_lo = q + lag_lo;
    const int64_t support_hi =
        support_lo + int64_t(taps) + int64_t(cached);
    const size_t n =
        right_edge ? static_cast<size_t>(support_hi - 2000)
                   : static_cast<size_t>(support_hi + 4096);
    cfg.max_window_samples = n;
    cfg.max_packet_samples = n;
    cfg.sub_roi_hi = int64_t(n);
    // Re-prepare with the final window size (no hot-path allocation later).
    if (!sic::prepare_sic_unified(cfg, t, scratch, &err))
        die("edge re-prepare failed: " + err);

    sic::SicUnifiedChannel ch;
    ch.valid = true;
    ch.backend = sic::SicEstimatorBackend::Ipatov;
    ch.status = "model_ready";
    ch.q = q;
    ch.cfo_hz = 0.0;
    ch.phase_reference_sample = q;
    ch.observed_end_sample = int64_t(n);
    ch.lag_lo = lag_lo;
    ch.tap_count = taps;
    ch.active_count = 4;
    ch.impulse_template = true;
    ch.model_id = 0xED6Eull;
    ch.fit_residual = 0.01;
    ch.gain_condition = 12.0;
    const std::pair<int64_t, cf> active[] = {
        { -4, cf(0.40f, -0.10f) },
        { 0, cf(0.90f, 0.20f) },
        { 7, cf(-0.25f, 0.15f) },
        { 15, cf(0.15f, 0.05f) },
    };
    for (size_t k = 0; k < 4; ++k) {
        ch.active_lags[k] = static_cast<int16_t>(active[k].first);
        ch.taps[static_cast<size_t>(active[k].first - lag_lo)] = active[k].second;
    }

    // Build y from the exact same event/tap model, truncated to the window.
    std::vector<cf> y(n, cf(float(kCfoFloor), float(-kCfoFloor)));
    for (const auto& e : scratch.cached_preamble_events) {
        for (size_t d = 0; d < taps; ++d) {
            const cf h = ch.taps[d];
            if (h == cf{})
                continue;
            const int64_t lag = lag_lo + int64_t(d);
            const int64_t j = q + int64_t(e.sample_offset) + lag;
            if (j < 0 || j >= int64_t(n))
                continue;
            y[static_cast<size_t>(j)] += float(e.sign) * h;
        }
    }

    Measured mm;
    mm.scenario = "edge";
    mm.mode = std::string(right_edge ? "right_n_plus_2000" : "left_lo_minus_2") +
              (clip ? "_clip_on" : "_clip_off");
    mm.n = int64_t(n);
    mm.q = q;
    mm.wave_samples = cached;
    mm.tap_count = taps;
    mm.lag_lo = lag_lo;
    mm.pulse_events = scratch.cached_preamble_events.size();
    mm.inplace = false;
    mm.roi_extra_kind = "none";
    std::vector<cf> out(n, cf{});
    bool first = true;

    for (int it = 0; it < args.warmup + args.reps; ++it) {
        const bool measure = (it >= args.warmup);
        std::fill(out.begin(), out.end(), cf{});
        sic::SicUnifiedResult r;
        r.out = out.data();
        r.i_hat = nullptr;
        const auto t0 = Clock::now();
        const bool ok = sic::reconstruct_and_subtract_unified(
            y.data(), n, ch, sic::SicCancellationLevel::L1_Preamble,
            sic::SicKnownFields::PreambleSfd, t.preamble, t.preamble_len, t,
            cfg, scratch, out.data(), nullptr, r);
        const double wall = us_since(t0);
        // A fail-closed transaction is a legitimate edge result (e.g. the
        // baseline returns false with reason "context_truncated").  It is not
        // a harness error.
        (void)ok;
        if (!window_is_finite(out.data(), n))
            die("non-finite edge output");

        if (first) {
            mm.committed = r.committed;
            mm.level = level_name(r.level);
            mm.reason = r.reason ? r.reason : "null";
            mm.support_lo = r.support_lo;
            mm.support_hi = r.support_hi;
            mm.sub_lo = r.sub_lo;
            mm.sub_hi = r.sub_hi;
            mm.fused_energy = r.suppression_energy_valid;
            mm.suppression_db = r.suppression_db;
            if (r.committed && r.sub_hi > r.sub_lo)
                mm.roi_extra_kind = r.suppression_energy_valid ? "fused_zero"
                                                               : "standalone";
            first = false;
        } else {
            mm.consistent = mm.consistent &&
                            (r.committed == mm.committed) &&
                            (std::strcmp(level_name(r.level), mm.level) == 0);
        }

        if (!measure)
            continue;

        double roi_extra = 0.0;
        if (r.committed && r.sub_hi > r.sub_lo &&
            !r.suppression_energy_valid) {
            const auto e0 = Clock::now();
            double e_y = 0.0, e_out = 0.0;
            sic::compute_roi_energies(y.data(), out.data(), r.sub_lo,
                                      r.sub_hi, e_y, e_out);
            roi_extra = us_since(e0);
        }
        mm.wall.add(wall);
        // reconstruct_and_subtract_unified does not populate the top-level
        // core_total timer (that belongs to sic_process_unified); use the
        // caller-side steady_clock for the edge row.
        mm.core.add(wall);
        mm.recon.add(r.timings.reconstruct_us);
        mm.sub.add(r.timings.subtract_us);
        mm.roi_extra.add(roi_extra);
        mm.back_half.add(r.timings.reconstruct_us + r.timings.subtract_us +
                         roi_extra);
        mm.zero.add(r.timings.recon_zero_us);
        mm.scale.add(r.timings.recon_scale_us);
        mm.conv.add(r.timings.recon_conv_us);
        mm.cfo_rot.add(r.timings.recon_cfo_rot_us);
        if (!args.raw_path.empty()) {
            std::ofstream raw(args.raw_path, std::ios::app);
            raw << "{\"scenario\":\"edge\",\"mode\":\"" << mm.mode
                << "\",\"iter\":" << it << ",\"committed\":"
                << (r.committed ? "true" : "false")
                << ",\"wall_us\":" << wall
                << ",\"core_total_us\":" << wall
                << ",\"reconstruct_us\":" << r.timings.reconstruct_us
                << ",\"subtract_us\":" << r.timings.subtract_us
                << ",\"roi_extra_us\":" << roi_extra
                << ",\"suppression_db\":" << mm.suppression_db << "}\n";
        }
    }
    return mm;
}

int run_edge(const Args& args)
{
    Packet pkt(10, 64, "4z2", 16, 0xED6Eu, 120000);
    const auto t = pkt.templates();
    sic::SicUnifiedScratch scratch;
    for (bool right : { false, true }) {
        for (bool clip : { false, true }) {
            Measured m = run_edge_case(args, right, clip, scratch, t);
            if (args.json)
                emit_json(m, args.reps, args.warmup, stdout);
            else
                emit_text(m, args.reps, args.warmup);
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------

void usage()
{
    std::fprintf(
        stderr,
        "usage: bench_sic_reconstruct_opt --scenario "
        "live-hit|short-l1|short-l2|cfo|taps|edge [--inplace] "
        "[--cfo-hz HZ] [--taps N] [--reps N] [--warmup N] [--json] "
        "[--raw PATH]\n");
}

Args parse_args(int argc, char** argv)
{
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc)
                die(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (k == "--scenario")
            a.scenario = need("--scenario");
        else if (k == "--inplace")
            a.inplace_only = true;
        else if (k == "--cfo-hz")
            a.cfo_hz = std::atof(need("--cfo-hz"));
        else if (k == "--taps")
            a.taps = std::atoi(need("--taps"));
        else if (k == "--reps")
            a.reps = std::atoi(need("--reps"));
        else if (k == "--warmup")
            a.warmup = std::atoi(need("--warmup"));
        else if (k == "--json")
            a.json = true;
        else if (k == "--raw")
            a.raw_path = need("--raw");
        else if (k == "--help" || k == "-h") {
            usage();
            std::exit(0);
        } else {
            usage();
            die("unknown argument: " + k);
        }
    }
    if (a.scenario.empty()) {
        usage();
        die("--scenario is required");
    }
    if (a.reps < 1)
        die("--reps must be >= 1");
    if (a.warmup < 0)
        die("--warmup must be >= 0");
    return a;
}

} // namespace

int main(int argc, char** argv)
{
    const Args args = parse_args(argc, argv);
    if (args.scenario == "live-hit")
        return run_live_hit(args);
    if (args.scenario == "short-l1")
        return run_short(args, false);
    if (args.scenario == "short-l2")
        return run_short(args, true);
    if (args.scenario == "cfo")
        return run_cfo(args);
    if (args.scenario == "taps")
        return run_taps(args);
    if (args.scenario == "edge")
        return run_edge(args);
    usage();
    die("unknown scenario: " + args.scenario);
}
