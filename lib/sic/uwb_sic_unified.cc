/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Unified UWB SIC Implementation (Stage M2 Contract).
 * Namespace openuwb::sic.
 */

#include <openuwb/sic/sic_unified.h>
#include <openuwb/support/timing.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

namespace openuwb {
namespace sic {

namespace {

using Clock = openuwb::timing::Clock;

double us(Clock::time_point t)
{
    return std::chrono::duration<double, std::micro>(Clock::now() - t).count();
}

constexpr uint64_t kFnvOffset = 14695981039346656037ull;
constexpr uint64_t kFnvPrime  = 1099511628211ull;

inline void hash_u64(uint64_t& h, uint64_t v)
{
    h = (h ^ v) * kFnvPrime;
}

inline void hash_bytes(uint64_t& h, const void* data, size_t n)
{
    const auto* p = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < n; ++i) {
        h = (h ^ static_cast<uint64_t>(p[i])) * kFnvPrime;
    }
}

inline void hash_words64(uint64_t& h, const void* data, size_t n)
{
    const auto* p64 = static_cast<const uint64_t*>(data);
    const size_t n64 = n / 8;
    for (size_t i = 0; i < n64; ++i) {
        h = (h ^ p64[i]) * kFnvPrime;
    }
    const auto* p8 = static_cast<const uint8_t*>(data) + n64 * 8;
    for (size_t i = 0; i < (n & 7); ++i) {
        h = (h ^ static_cast<uint64_t>(p8[i])) * kFnvPrime;
    }
}

} // namespace

// ---------------------------------------------------------------------------
// SicUnifiedChannel Conversions
// ---------------------------------------------------------------------------

void SicUnifiedChannel::from_discrete_cir(const DiscreteCirModel& m)
{
    valid = m.valid;
    backend = m.impulse_template ? SicEstimatorBackend::Ipatov : SicEstimatorBackend::DiscreteCir;
    status = m.status;
    reason = m.reason;
    model_id = m.model_id;
    q = m.q;
    cfo_hz = m.cfo_hz;
    phase_reference_sample = m.phase_reference_sample;
    observed_end_sample = m.observed_end_sample;
    training_lo = m.training_lo;
    training_hi = m.training_hi;
    fit_residual = m.fit_residual;
    gain_condition = m.gain_condition;
    repeats_confirmed = m.repeats_confirmed;
    lag_lo = m.lag_lo;
    tap_count = m.tap_count;
    active_count = m.active_count;
    active_lags = m.active_lags;
    taps = m.taps;
    impulse_template = m.impulse_template;
    cir_repetitions = m.cir_repetitions;
    n_paths = 0;
}

void SicUnifiedChannel::from_parametric(const SicModel& m)
{
    valid = m.valid;
    backend = SicEstimatorBackend::PathFractional;
    status = m.status;
    reason = m.reason;
    model_id = m.model_id;
    q = m.origin;
    cfo_hz = m.cfo_hz;
    phase_reference_sample = m.phase_reference_sample;
    observed_end_sample = m.observed_end_sample;
    training_lo = m.training_lo;
    training_hi = m.training_hi;
    fit_residual = m.fit_residual;
    gain_condition = m.gain_condition;
    repeats_confirmed = m.training_reps_confirmed;
    lag_lo = 0;
    tap_count = 0;
    active_count = 0;
    impulse_template = false;
    cir_repetitions = 0;
    n_paths = m.n_paths;
    for (size_t k = 0; k < kSicMaxPaths; ++k) {
        tau_abs[k] = m.tau_abs[k];
        gain[k] = m.gain[k];
    }
}

void SicUnifiedChannel::to_discrete_cir(DiscreteCirModel& m) const
{
    m.valid = valid;
    m.status = status;
    m.reason = reason;
    m.model_id = model_id;
    m.q = q;
    m.cfo_hz = cfo_hz;
    m.phase_reference_sample = phase_reference_sample;
    m.observed_end_sample = observed_end_sample;
    m.training_lo = training_lo;
    m.training_hi = training_hi;
    m.fit_residual = fit_residual;
    m.gain_condition = gain_condition;
    m.repeats_confirmed = repeats_confirmed;
    m.lag_lo = lag_lo;
    m.tap_count = tap_count;
    m.active_count = active_count;
    m.active_lags = active_lags;
    m.taps = taps;
    m.impulse_template = impulse_template;
    m.cir_repetitions = cir_repetitions;
}

void SicUnifiedChannel::to_parametric(SicModel& m) const
{
    m.valid = valid;
    m.state = valid ? SicState::Committed : SicState::Failed;
    m.status = status;
    m.reason = reason;
    m.model_id = model_id;
    m.origin = q;
    m.cfo_hz = cfo_hz;
    m.phase_reference_sample = phase_reference_sample;
    m.observed_end_sample = observed_end_sample;
    m.training_lo = training_lo;
    m.training_hi = training_hi;
    m.fit_residual = fit_residual;
    m.gain_condition = gain_condition;
    m.training_reps_confirmed = repeats_confirmed;
    m.n_paths = n_paths;
    for (size_t k = 0; k < kSicMaxPaths; ++k) {
        m.tau_abs[k] = tau_abs[k];
        m.gain[k] = gain[k];
    }
}

// ---------------------------------------------------------------------------
// SicUnifiedConfig Conversions
// ---------------------------------------------------------------------------

SicConfig SicUnifiedConfig::to_sic_config() const
{
    SicConfig c;
    c.mode = SicMode::Estimate;
    c.comm_code = comm_code;
    c.sync_reps = sync_reps;
    c.samples_per_symbol = samples_per_symbol;
    c.sample_rate = sample_rate;
    c.k_paths = k_paths;
    c.fit_reps = fit_reps;
    c.roi_lo = roi_lo;
    c.roi_hi = roi_hi;
    c.sub_roi_lo = sub_roi_lo;
    c.sub_roi_hi = sub_roi_hi;
    c.sub_min_overlap = sub_min_overlap;
    c.coarse_stride = coarse_stride;
    c.post_tail = post_tail;
    c.enable_fractional = enable_fractional;
    c.enable_cfo = enable_cfo;
    c.detect_threshold = detect_threshold;
    // Forward the S1 holdout gate.  Without it the PathFractional backend
    // rejects every window with "uncalibrated_validate_threshold".
    c.validate_threshold = validate_threshold;
    c.observed_end_sample = observed_end_sample;
    c.fit_residual_max = fit_residual_max;
    return c;
}

DiscreteCirConfig SicUnifiedConfig::to_discrete_cir_config() const
{
    DiscreteCirConfig c;
    c.ipatov_correlation = (backend == SicEstimatorBackend::Ipatov);
    c.allow_configured_profile = true;
    c.sparse_pulse_reconstruct = (backend == SicEstimatorBackend::Ipatov && sparse_pulse_reconstruct);
    c.comm_code = comm_code;
    c.sync_reps = sync_reps;
    c.samples_per_symbol = samples_per_symbol;
    c.sfd_mode = sfd_mode;
    c.sample_rate = sample_rate;
    c.fit_reps = fit_reps;
    c.fit_skip_reps = fit_skip_reps;
    c.roi_lo = roi_lo;
    c.roi_hi = roi_hi;
    c.lag_lo = lag_lo;
    c.tap_count = tap_count;
    c.peak_floor = peak_floor;
    c.peak_min_spacing = peak_min_spacing;
    c.active_radius = active_radius;
    c.lambda = lambda;
    c.coarse_stride = coarse_stride;
    c.detect_threshold = detect_threshold;
    c.fit_residual_max = fit_residual_max;
    c.gain_condition_max = gain_condition_max;
    c.max_window_samples = max_window_samples;
    c.max_packet_samples = max_packet_samples;
    return c;
}

mod::HrpModConfig SicUnifiedConfig::to_hrp_mod_config() const
{
    mod::HrpModConfig c;
    c.code_index = comm_code;
    c.sync_repetitions = sync_reps;
    c.sfd_mode = sfd_mode;
    c.peak_amplitude = 0.8f;
    return c;
}

// ---------------------------------------------------------------------------
// SicUnifiedTemplates Conversions
// ---------------------------------------------------------------------------

SicTemplates SicUnifiedTemplates::to_sic_templates() const
{
    SicTemplates t;
    t.sync = sync;
    t.sync_len = sync_len;
    t.preamble = preamble;
    t.preamble_len = preamble_len;
    return t;
}

DiscreteCirTemplates SicUnifiedTemplates::to_discrete_cir_templates() const
{
    DiscreteCirTemplates t;
    t.sync = sync;
    t.sync_len = sync_len;
    t.train = preamble;
    t.train_len = preamble_len;
    return t;
}

// ---------------------------------------------------------------------------
// Hash Identity Computations
// ---------------------------------------------------------------------------

uint64_t compute_unified_config_hash(const SicUnifiedConfig& cfg)
{
    uint64_t h = kFnvOffset;
    hash_u64(h, 0x534943554e494631ull); // "SICUNIF1"
    hash_u64(h, static_cast<uint64_t>(cfg.comm_code));
    hash_u64(h, static_cast<uint64_t>(cfg.sync_reps));
    hash_u64(h, static_cast<uint64_t>(cfg.samples_per_symbol));
    uint64_t fs_bits = 0;
    std::memcpy(&fs_bits, &cfg.sample_rate, sizeof(double));
    hash_u64(h, fs_bits);
    if (cfg.sfd_mode) {
        hash_bytes(h, cfg.sfd_mode, std::strlen(cfg.sfd_mode));
    }
    hash_u64(h, static_cast<uint64_t>(cfg.backend));
    hash_u64(h, cfg.sparse_pulse_reconstruct ? 1ull : 0ull);
    hash_u64(h, static_cast<uint64_t>(cfg.fit_reps));
    hash_u64(h, static_cast<uint64_t>(cfg.fit_skip_reps));
    hash_u64(h, cfg.allow_partial_window_clip ? 1ull : 0ull);
    hash_u64(h, cfg.enable_empirical_phase_tracking ? 1ull : 0ull);
    hash_u64(h, static_cast<uint64_t>(cfg.lag_lo));
    hash_u64(h, static_cast<uint64_t>(cfg.tap_count));
    hash_u64(h, static_cast<uint64_t>(cfg.k_paths));
    hash_u64(h, static_cast<uint64_t>(cfg.post_tail));
    return h;
}

uint64_t compute_unified_template_hash(const SicUnifiedTemplates& tmpl)
{
    uint64_t h = kFnvOffset;
    hash_u64(h, 0x54454d504c415431ull); // "TEMPLAT1"
    hash_u64(h, static_cast<uint64_t>(tmpl.sync_len));
    if (tmpl.sync && tmpl.sync_len > 0) {
        hash_words64(h, tmpl.sync, tmpl.sync_len * sizeof(std::complex<float>));
    }
    hash_u64(h, static_cast<uint64_t>(tmpl.preamble_len));
    hash_u64(h, reinterpret_cast<uintptr_t>(tmpl.preamble));
    hash_u64(h, static_cast<uint64_t>(tmpl.impulse_sync_len));
    if (tmpl.impulse_sync && tmpl.impulse_sync_len > 0) {
        hash_words64(h, tmpl.impulse_sync, tmpl.impulse_sync_len * sizeof(std::complex<float>));
    }
    return h;
}

// ---------------------------------------------------------------------------
// SicUnifiedScratch Reservation & Cache Management
// ---------------------------------------------------------------------------

void SicUnifiedScratch::reserve(size_t max_window_samples, size_t max_packet_samples)
{
    if (max_packet_samples == 0) {
        max_packet_samples = max_window_samples;
    }
    if (max_window_samples > capacity_window_samples) {
        capacity_window_samples = max_window_samples;
    }
    if (max_packet_samples > capacity_packet_samples) {
        capacity_packet_samples = max_packet_samples;
    }

    discrete.reserve(capacity_window_samples);
    fullpacket.reserve(capacity_window_samples);
    modulator.reserve(127, capacity_packet_samples);
    demod.reserve(capacity_window_samples);

    if (work_waveform.capacity() < capacity_packet_samples) {
        work_waveform.reserve(capacity_packet_samples);
    }
    work_waveform.resize(capacity_packet_samples);

    if (work_accum.capacity() < capacity_window_samples) {
        work_accum.reserve(capacity_window_samples);
    }
    work_accum.resize(capacity_window_samples);

    if (work_rot.capacity() < capacity_window_samples) {
        work_rot.reserve(capacity_window_samples);
    }
    work_rot.resize(capacity_window_samples);

    if (work_phr_bits.capacity() < 19) {
        work_phr_bits.reserve(19);
    }
    work_phr_bits.resize(19);

    if (work_wire_psdu.capacity() < radar_meta::kMaxPsduBytes) {
        work_wire_psdu.reserve(radar_meta::kMaxPsduBytes);
    }
    work_wire_psdu.resize(radar_meta::kMaxPsduBytes);
}

bool SicUnifiedScratch::cache_matches(const SicUnifiedConfig& cfg, const SicUnifiedTemplates& tmpl) const
{
    if (!prepared || !preamble_cache_valid) {
        return false;
    }
    const uint64_t ch = compute_unified_config_hash(cfg);
    const uint64_t th = compute_unified_template_hash(tmpl);
    return (config_hash == ch && template_hash == th &&
            preamble_cache_config_hash == ch &&
            preamble_cache_template_hash == th);
}

void SicUnifiedScratch::invalidate_cache()
{
    preamble_cache_valid = false;
    preamble_cache_config_hash = 0;
    preamble_cache_template_hash = 0;
    cached_preamble_samples = 0;
    cached_sync_vec.clear();
    cached_preamble_waveform.clear();
    cached_preamble_events.clear();
    cached_preamble_unit_waveform.clear();
    cached_sfd_signs.fill(0);
    cached_sfd_len = 0;
    config_hash = 0;
    template_hash = 0;
    prepared = false;
}

// ---------------------------------------------------------------------------
// prepare_sic_unified
// ---------------------------------------------------------------------------

bool prepare_sic_unified(
    const SicUnifiedConfig& cfg,
    const SicUnifiedTemplates& templates,
    SicUnifiedScratch& scratch,
    std::string* err)
{
    scratch.invalidate_cache();

    if (!templates.sync || templates.sync_len != cfg.samples_per_symbol) {
        if (err) *err = "invalid_sync_template";
        return false;
    }
    if (!templates.preamble || templates.preamble_len == 0) {
        if (err) *err = "invalid_preamble_template";
        return false;
    }

    scratch.reserve(cfg.max_window_samples, cfg.max_packet_samples);
    // The unified layer validates the caller's public profile; downstream
    // full-packet transactions may therefore use the configured code/SYNC/SFD
    // instead of the historical code10/64/4z2 gate.
    scratch.fullpacket.allow_configured_profile = true;

    if (cfg.backend == SicEstimatorBackend::PathFractional) {
        SicConfig sc = cfg.to_sic_config();
        if (!prepare_sic_core(sc, templates.sync, templates.sync_len,
                              templates.preamble, templates.preamble_len,
                              scratch.fractional, err)) {
            return false;
        }
    } else {
        DiscreteCirConfig dc = cfg.to_discrete_cir_config();
        DiscreteCirTemplates dt = templates.to_discrete_cir_templates();
        if (!prepare_discrete_cir(dc, dt, scratch.discrete, err)) {
            return false;
        }
    }

    // 1. Cached sync_vec for receiver demod
    scratch.cached_sync_vec.assign(templates.sync, templates.sync + templates.sync_len);

    // 2. Materialized prefix complex waveform (SYNC + SFD)
    mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
    if (scratch.cached_preamble_waveform.capacity() < scratch.capacity_packet_samples) {
        scratch.cached_preamble_waveform.reserve(scratch.capacity_packet_samples);
    }
    scratch.cached_preamble_waveform.resize(scratch.capacity_packet_samples);
    size_t n_prefix = 0;
    if (mod::modulate_prefix_waveform(modcfg, scratch.modulator,
                                      scratch.cached_preamble_waveform.data(),
                                      scratch.cached_preamble_waveform.size(),
                                      n_prefix) && n_prefix > 0) {
        scratch.cached_preamble_waveform.resize(n_prefix);
    } else {
        scratch.cached_preamble_waveform.clear();
    }

    // 3. Dedicated cached preamble pulse events for sparse pulse reconstruction
    mod::HrpModConfig imp_cfg = modcfg;
    imp_cfg.peak_amplitude = 1.f;
    imp_cfg.pulse.shape = mod::PulseShape::External;
    imp_cfg.external_taps = {1.f, 0.f};
    imp_cfg.pulse.external_n_taps = 2;

    const mod::HrpModConfig& prof = (cfg.backend == SicEstimatorBackend::Ipatov &&
                                     scratch.discrete.impulse_profile.code_index != 0)
        ? scratch.discrete.impulse_profile
        : imp_cfg;
    auto& mod_scratch = (cfg.backend == SicEstimatorBackend::Ipatov)
        ? scratch.discrete.modulator
        : scratch.modulator;

    mod::PulseList pulses;
    if (mod::modulate_prefix_pulse_list(prof, mod_scratch,
                                        scratch.capacity_packet_samples, pulses) && pulses.count > 0) {
        scratch.cached_preamble_events.assign(pulses.events, pulses.events + pulses.count);
        scratch.cached_preamble_samples = pulses.packet_samples;
    } else {
        scratch.cached_preamble_events.clear();
        scratch.cached_preamble_samples = 0;
    }

    // 4. Ipatov unit-impulse prefix waveform (SYNC + SFD)
    // The Ipatov model is fitted against the unit impulse reference and its
    // taps already contain the pulse shaping, so the L1 reconstruction input
    // is the unit impulse prefix, not the shaped one.  Materializing it here
    // keeps the hot path allocation-free and lets the non-sparse apply path
    // run the same transaction as the shaped discrete backend.
    scratch.cached_preamble_unit_waveform.clear();
    if (cfg.backend == SicEstimatorBackend::Ipatov &&
        scratch.cached_preamble_samples > 0 &&
        !scratch.cached_preamble_events.empty()) {
        scratch.cached_preamble_unit_waveform.assign(
            scratch.cached_preamble_samples, std::complex<float>(0.f, 0.f));
        bool ok = true;
        for (const mod::PulseEvent& e : scratch.cached_preamble_events) {
            if (e.sample_offset >= scratch.cached_preamble_samples) {
                ok = false;
                break;
            }
            scratch.cached_preamble_unit_waveform[e.sample_offset] =
                std::complex<float>(static_cast<float>(e.sign), 0.f);
        }
        if (!ok)
            scratch.cached_preamble_unit_waveform.clear();
    }

    const auto sfd_seq = demod::GetSfdSequence(cfg.sfd_mode ? cfg.sfd_mode : "4z2");
    scratch.cached_sfd_len = std::min(sfd_seq.size(), scratch.cached_sfd_signs.size());
    for (size_t i = 0; i < scratch.cached_sfd_len; ++i) {
        scratch.cached_sfd_signs[i] = sfd_seq[i];
    }

    scratch.config_hash = compute_unified_config_hash(cfg);
    scratch.template_hash = compute_unified_template_hash(templates);
    scratch.preamble_cache_config_hash = scratch.config_hash;
    scratch.preamble_cache_template_hash = scratch.template_hash;
    scratch.preamble_cache_valid = true;
    scratch.prepared = true;
    return true;
}

// ---------------------------------------------------------------------------
// estimate_channel_unified
// ---------------------------------------------------------------------------

bool estimate_channel_unified(
    const std::complex<float>* y,
    size_t n,
    int64_t observed_end,
    const SicUnifiedTemplates& templates,
    const SicUnifiedConfig& cfg,
    SicUnifiedScratch& scratch,
    SicUnifiedChannel& channel,
    SicUnifiedTimings* timings)
{
    channel.valid = false;
    channel.backend = cfg.backend;
    if (!y || n == 0 || !scratch.prepared) {
        channel.status = "invalid_input";
        return false;
    }
    const int64_t obs = (observed_end <= 0) ? int64_t(n) : std::min(int64_t(n), observed_end);

    const auto t_start = Clock::now();

    if (cfg.backend == SicEstimatorBackend::PathFractional) {
        SicModel model;
        SicTimings st;
        SicConfig sc = cfg.to_sic_config();
        sc.observed_end_sample = obs;
        SicTemplates stemp = templates.to_sic_templates();
        bool ok = estimate_model(y, n, obs, stemp, sc, scratch.fractional, model, &st);
        channel.from_parametric(model);
        if (timings) {
            timings->estimate_us = us(t_start);
            timings->est_detect_us = static_cast<double>(st.detect_coarse_us + st.detect_fine_us);
            timings->est_cfo_us = static_cast<double>(st.cfo_us);
            timings->est_cir_us = static_cast<double>(st.cir_us + st.path_select_us);
            timings->est_solve_us = static_cast<double>(st.matrix_build_us + st.gain_fit_us);
            timings->est_fit_check_us = static_cast<double>(st.model_prepare_us);
        }
        return ok;
    } else {
        DiscreteCirModel dm;
        DiscreteCirTimings dt;
        DiscreteCirConfig dc = cfg.to_discrete_cir_config();
        DiscreteCirTemplates dtemp = templates.to_discrete_cir_templates();
        bool ok = estimate_discrete_cir(y, n, obs, dtemp, dc, scratch.discrete, dm, &dt);
        channel.from_discrete_cir(dm);
        if (timings) {
            timings->estimate_us = us(t_start);
            timings->est_detect_us = dt.detect_us;
            timings->est_cfo_us = dt.cfo_us;
            timings->est_cir_us = dt.support_us;
            timings->est_solve_us = dt.projection_us + dt.solve_us;
            timings->est_fit_check_us = dt.fit_check_us;
        }
        return ok;
    }
}

// ---------------------------------------------------------------------------
// reconstruct_and_subtract_unified
// ---------------------------------------------------------------------------

bool reconstruct_and_subtract_unified(
    const std::complex<float>* y,
    size_t n,
    const SicUnifiedChannel& channel,
    SicCancellationLevel level,
    SicKnownFields known_fields,
    const std::complex<float>* s_known,
    size_t s_known_len,
    const SicUnifiedTemplates& templates,
    const SicUnifiedConfig& cfg,
    SicUnifiedScratch& scratch,
    std::complex<float>* out,
    std::complex<float>* i_hat,
    SicUnifiedResult& result)
{
    // True once some backend has already seeded `out` from `y` for this
    // transaction (apply_model's own fail-closed copy).  Keeps the unified
    // fail_closed lambda to AT MOST one bypass copy.
    bool bypass_copy_done = false;
    auto fail_closed = [&](const char* reason) {
        if (out && y && n > 0 && out != y && result.copy_out_on_bypass &&
            !bypass_copy_done) {
            std::copy(y, y + n, out);
        }
        if (i_hat && i_hat != y && i_hat != out) {
            std::fill(i_hat, i_hat + n, std::complex<float>(0.f, 0.f));
        }
        result.committed = false;
        result.reason = reason;
        result.level = SicCancellationLevel::L0_Bypass;
        result.cancelled_fields = SicKnownFields::None;
        return false;
    };

    if (!y || !out || n == 0 || !channel.valid) {
        return fail_closed("invalid_input_or_channel");
    }

    if (channel.backend == SicEstimatorBackend::PathFractional) {
        SicModel model;
        channel.to_parametric(model);
        SicConfig sc = cfg.to_sic_config();
        SicTemplates stemp = templates.to_sic_templates();
        SicTimings st;

        if (level == SicCancellationLevel::L1_Preamble) {
            // apply_model() owns a fail-closed y -> out copy (and i_hat zero
            // fill) unless the caller preinitializes both.  Seed them exactly
            // once here so fail_closed never repeats the copy; skipping the
            // seed when out == y keeps the in-place mode memcpy-free.
            if (out != y) {
                std::copy(y, y + n, out);
            }
            if (i_hat && i_hat != y && i_hat != out) {
                std::fill(i_hat, i_hat + n, std::complex<float>(0.f, 0.f));
            }
            bypass_copy_done = true;
            bool ok = apply_model(model, y, n, 0, stemp, sc, scratch.fractional,
                                  out, i_hat, &st, /*out_preinitialized=*/true);
            // apply_model() reports `reconstruct_us` over the whole per-chunk
            // synthesis INCLUDING the masked subtraction, and `apply_us` as
            // the parent total.  Attribute the remainder to subtraction
            // instead of timing the caller's wall clock, which used to fold
            // the pre-seed copy into subtract_us.
            result.timings.reconstruct_us = static_cast<double>(st.reconstruct_us);
            result.timings.subtract_us =
                (st.apply_us > st.reconstruct_us)
                    ? static_cast<double>(st.apply_us - st.reconstruct_us)
                    : 0.0;
            if (!ok) {
                return fail_closed("apply_model_failed");
            }
            result.committed = true;
            result.reason = "committed_preamble_fractional";
            result.level = SicCancellationLevel::L1_Preamble;
            result.cancelled_fields = SicKnownFields::Preamble;
            result.support_lo = model.predicted_support_lo;
            result.support_hi = model.predicted_support_hi;
            result.sub_lo = model.sub_lo;
            result.sub_hi = model.sub_hi;
            result.sub_overlap = (result.sub_hi > result.sub_lo) ? (result.sub_hi - result.sub_lo) : 0;
            return true;
        } else {
            // FullPacket with PathFractional
            mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
            FullPacketResult f = fullpacket_cancel(
                y, n, scratch.work_phr_bits, scratch.work_wire_psdu,
                true, modcfg, model, templates.sync, templates.sync_len,
                cfg.sub_roi_lo, cfg.sub_roi_hi, scratch.fullpacket, out, i_hat);
            result.timings.reconstruct_us = f.channel_us;
            result.timings.subtract_us = f.subtract_us;
            result.waveform_scale = f.waveform_scale;
            if (!f.committed) {
                return fail_closed(f.reason);
            }
            result.committed = true;
            result.reason = "committed_full_packet_fractional";
            result.level = SicCancellationLevel::L2_FullPacket;
            result.cancelled_fields = SicKnownFields::FullPacket;
            result.sub_lo = f.sub_lo;
            result.sub_hi = f.sub_hi;
            result.sub_overlap = (result.sub_hi > result.sub_lo) ? (result.sub_hi - result.sub_lo) : 0;
            return true;
        }
    } else {
        // Discrete CIR / Ipatov backend
        DiscreteCirModel dm;
        channel.to_discrete_cir(dm);
        DiscreteCirResult dres;

        // A committed *_ex transaction with compute_roi_energy = true has
        // already accumulated ||y||^2 and ||out||^2 over the subtraction ROI
        // BEFORE writing `out`, so the suppression metric never re-reads y and
        // stays correct for in-place (out == y) cancellation.
        auto accept_fused_energy = [&](const DiscreteCirResult& dr) {
            if (!dr.roi_energy_valid) {
                return;
            }
            result.suppression_energy_valid = true;
            if (std::isfinite(dr.roi_e_y) && std::isfinite(dr.roi_e_out) &&
                dr.roi_e_y > 0.0 && dr.roi_e_out > 0.0) {
                result.suppression_db = 10.0 * std::log10(dr.roi_e_y / dr.roi_e_out);
            }
        };

        DiscreteCirApplyOptions opts;
        // The unified fail_closed lambda owns the single bypass copy, so the
        // low-level transaction must never copy y -> out on failure itself.
        opts.copy_out_on_fail = false;
        opts.allow_partial_clip = cfg.allow_partial_window_clip;
        opts.compute_roi_energy = true;

        if (level == SicCancellationLevel::L1_Preamble) {
            // Check if sparse pulse reconstruction is enabled for Ipatov backend
            if (dm.impulse_template && scratch.discrete.sparse_pulse_reconstruct) {
                mod::PulseList pulses;
                const bool cache_hit =
                    scratch.preamble_cache_valid && !scratch.cached_preamble_events.empty();
                if (cache_hit) {
                    pulses.events = scratch.cached_preamble_events.data();
                    pulses.count = scratch.cached_preamble_events.size();
                    pulses.packet_samples = scratch.cached_preamble_samples;
                } else {
                    const bool tx_ok = mod::modulate_prefix_pulse_list(
                        scratch.discrete.impulse_profile, scratch.discrete.modulator,
                        scratch.discrete.waveform.size(), pulses);
                    if (!tx_ok || pulses.count == 0) {
                        pulses = mod::PulseList{};
                    }
                }
                if (pulses.count > 0) {
                    // Only the cached event list is known prepared (order,
                    // signs and the estimator's own y-finite scan); a freshly
                    // modulated pulse list keeps the public rescan.
                    opts.transaction_prepared = cache_hit;
                    bool ok = apply_discrete_cir_pulses_ex(
                        y, n, pulses, templates.sync, templates.sync_len, dm,
                        cfg.sub_roi_lo, cfg.sub_roi_hi, opts, scratch.discrete, out, i_hat,
                        &dres);
                    result.timings.reconstruct_us = dres.channel_us;
                    result.timings.subtract_us = dres.subtract_us;
                    result.timings.recon_zero_us = dres.zero_us;
                    result.timings.recon_scale_us = dres.scale_us;
                    result.timings.recon_conv_us = dres.convolution_us;
                    result.timings.recon_cfo_rot_us = dres.phase_us + dres.cfo_us;
                    result.waveform_scale = dres.waveform_scale;
                    if (!ok || !dres.committed) {
                        return fail_closed(dres.reason ? dres.reason : "apply_discrete_pulses_failed");
                    }
                    result.committed = true;
                    result.reason = "committed_preamble_sparse_pulse";
                    result.level = SicCancellationLevel::L1_Preamble;
                    result.cancelled_fields = known_fields;
                    result.support_lo = dres.support_lo;
                    result.support_hi = dres.support_hi;
                    result.sub_lo = dres.sub_lo;
                    result.sub_hi = dres.sub_hi;
                    result.sub_overlap = (result.sub_hi > result.sub_lo) ? (result.sub_hi - result.sub_lo) : 0;
                    accept_fused_energy(dres);
                    return true;
                }
            }

            // Ipatov: the reconstruction input must be the UNIT impulse prefix
            // (the taps already contain the shaping).  Using the shaped prefix
            // here would fail the waveform-scale consistency check against the
            // unit reference SYNC, so it is never a valid fallback.
            const std::complex<float>* wf = s_known ? s_known : templates.preamble;
            size_t wf_n = s_known ? s_known_len : templates.preamble_len;
            const char* commit_reason = "committed_preamble_discrete";
            if (dm.impulse_template) {
                if (scratch.cached_preamble_unit_waveform.empty()) {
                    return fail_closed("impulse_prefix_unavailable");
                }
                wf = scratch.cached_preamble_unit_waveform.data();
                wf_n = scratch.cached_preamble_unit_waveform.size();
                commit_reason = "committed_preamble_ipatov";
            }

            bool ok = apply_discrete_cir(
                y, n, wf, wf_n, templates.sync, templates.sync_len, dm,
                cfg.sub_roi_lo, cfg.sub_roi_hi, scratch.discrete, out, i_hat, &dres);
            result.timings.reconstruct_us = dres.channel_us;
            result.timings.subtract_us = dres.subtract_us;
            result.timings.recon_zero_us = dres.zero_us;
            result.timings.recon_scale_us = dres.scale_us;
            result.timings.recon_conv_us = dres.convolution_us;
            result.timings.recon_cfo_rot_us = dres.phase_us + dres.cfo_us;
            result.waveform_scale = dres.waveform_scale;
            if (!ok || !dres.committed) {
                return fail_closed(dres.reason ? dres.reason : "apply_discrete_failed");
            }
            result.committed = true;
            result.reason = commit_reason;
            result.level = SicCancellationLevel::L1_Preamble;
            result.cancelled_fields = known_fields;
            result.support_lo = dres.support_lo;
            result.support_hi = dres.support_hi;
            result.sub_lo = dres.sub_lo;
            result.sub_hi = dres.sub_hi;
            result.sub_overlap = (result.sub_hi > result.sub_lo) ? (result.sub_hi - result.sub_lo) : 0;
            accept_fused_energy(dres);
            return true;
        } else {
            // FullPacket with DiscreteCir / Ipatov.  Sparse experiments still
            // validate the remodulated pulse list themselves; the cached
            // prefix events stay untouched.
            mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
            dres = fullpacket_cancel_discrete_ex(
                y, n, scratch.work_phr_bits, scratch.work_wire_psdu,
                true, modcfg, dm, templates.sync, templates.sync_len,
                cfg.sub_roi_lo, cfg.sub_roi_hi, opts, scratch.discrete, out, i_hat);
            result.timings.tx_remod_us = dres.tx_us;
            result.timings.reconstruct_us = dres.channel_us;
            result.timings.subtract_us = dres.subtract_us;
            result.timings.recon_zero_us = dres.zero_us;
            result.timings.recon_scale_us = dres.scale_us;
            result.timings.recon_conv_us = dres.convolution_us;
            result.timings.recon_cfo_rot_us = dres.phase_us + dres.cfo_us;
            result.waveform_scale = dres.waveform_scale;
            if (!dres.committed) {
                return fail_closed(dres.reason);
            }
            result.committed = true;
            result.reason = "committed_full_packet_discrete";
            result.level = SicCancellationLevel::L2_FullPacket;
            result.cancelled_fields = SicKnownFields::FullPacket;
            result.support_lo = dres.support_lo;
            result.support_hi = dres.support_hi;
            result.sub_lo = dres.sub_lo;
            result.sub_hi = dres.sub_hi;
            result.sub_overlap = (result.sub_hi > result.sub_lo) ? (result.sub_hi - result.sub_lo) : 0;
            accept_fused_energy(dres);
            return true;
        }
    }
}

// ---------------------------------------------------------------------------
// compute_roi_energies (SIMD float load + double accumulation)
// ---------------------------------------------------------------------------

void compute_roi_energies(
    const std::complex<float>* y,
    const std::complex<float>* out,
    int64_t start,
    int64_t end,
    double& e_y,
    double& e_out)
{
    e_y = 0.0;
    e_out = 0.0;
    if (!y || !out || end <= start || start < 0) {
        return;
    }
    const size_t count = static_cast<size_t>(end - start);

#if defined(__AVX2__)
    const float* py = reinterpret_cast<const float*>(y + start);
    const float* pout = reinterpret_cast<const float*>(out + start);

    __m256d y_acc0 = _mm256_setzero_pd();
    __m256d y_acc1 = _mm256_setzero_pd();
    __m256d y_acc2 = _mm256_setzero_pd();
    __m256d y_acc3 = _mm256_setzero_pd();

    __m256d out_acc0 = _mm256_setzero_pd();
    __m256d out_acc1 = _mm256_setzero_pd();
    __m256d out_acc2 = _mm256_setzero_pd();
    __m256d out_acc3 = _mm256_setzero_pd();

    size_t i = 0;
    const size_t unroll8 = count & ~size_t(7);
    for (; i < unroll8; i += 8) {
        __m256 vy0 = _mm256_loadu_ps(py + 2 * i);
        __m256 vy1 = _mm256_loadu_ps(py + 2 * i + 8);
        __m256 vout0 = _mm256_loadu_ps(pout + 2 * i);
        __m256 vout1 = _mm256_loadu_ps(pout + 2 * i + 8);

        __m256 vy0_sq = _mm256_mul_ps(vy0, vy0);
        __m256 vy1_sq = _mm256_mul_ps(vy1, vy1);
        __m256 vout0_sq = _mm256_mul_ps(vout0, vout0);
        __m256 vout1_sq = _mm256_mul_ps(vout1, vout1);

        y_acc0 = _mm256_add_pd(y_acc0, _mm256_cvtps_pd(_mm256_castps256_ps128(vy0_sq)));
        y_acc1 = _mm256_add_pd(y_acc1, _mm256_cvtps_pd(_mm256_extractf128_ps(vy0_sq, 1)));
        y_acc2 = _mm256_add_pd(y_acc2, _mm256_cvtps_pd(_mm256_castps256_ps128(vy1_sq)));
        y_acc3 = _mm256_add_pd(y_acc3, _mm256_cvtps_pd(_mm256_extractf128_ps(vy1_sq, 1)));

        out_acc0 = _mm256_add_pd(out_acc0, _mm256_cvtps_pd(_mm256_castps256_ps128(vout0_sq)));
        out_acc1 = _mm256_add_pd(out_acc1, _mm256_cvtps_pd(_mm256_extractf128_ps(vout0_sq, 1)));
        out_acc2 = _mm256_add_pd(out_acc2, _mm256_cvtps_pd(_mm256_castps256_ps128(vout1_sq)));
        out_acc3 = _mm256_add_pd(out_acc3, _mm256_cvtps_pd(_mm256_extractf128_ps(vout1_sq, 1)));
    }

    if (i + 4 <= count) {
        __m256 vy0 = _mm256_loadu_ps(py + 2 * i);
        __m256 vout0 = _mm256_loadu_ps(pout + 2 * i);

        __m256 vy0_sq = _mm256_mul_ps(vy0, vy0);
        __m256 vout0_sq = _mm256_mul_ps(vout0, vout0);

        y_acc0 = _mm256_add_pd(y_acc0, _mm256_cvtps_pd(_mm256_castps256_ps128(vy0_sq)));
        y_acc1 = _mm256_add_pd(y_acc1, _mm256_cvtps_pd(_mm256_extractf128_ps(vy0_sq, 1)));

        out_acc0 = _mm256_add_pd(out_acc0, _mm256_cvtps_pd(_mm256_castps256_ps128(vout0_sq)));
        out_acc1 = _mm256_add_pd(out_acc1, _mm256_cvtps_pd(_mm256_extractf128_ps(vout0_sq, 1)));
        i += 4;
    }

    __m256d y_sum = _mm256_add_pd(_mm256_add_pd(y_acc0, y_acc1),
                                  _mm256_add_pd(y_acc2, y_acc3));
    alignas(32) double y_buf[4];
    _mm256_storeu_pd(y_buf, y_sum);
    double sum_y = y_buf[0] + y_buf[1] + y_buf[2] + y_buf[3];

    __m256d out_sum = _mm256_add_pd(_mm256_add_pd(out_acc0, out_acc1),
                                    _mm256_add_pd(out_acc2, out_acc3));
    alignas(32) double out_buf[4];
    _mm256_storeu_pd(out_buf, out_sum);
    double sum_out = out_buf[0] + out_buf[1] + out_buf[2] + out_buf[3];

    for (; i < count; ++i) {
        sum_y += std::norm(y[start + i]);
        sum_out += std::norm(out[start + i]);
    }
    e_y = sum_y;
    e_out = sum_out;
#else
    double sum_y = 0.0;
    double sum_out = 0.0;
    for (size_t i = 0; i < count; ++i) {
        sum_y += std::norm(y[start + i]);
        sum_out += std::norm(out[start + i]);
    }
    e_y = sum_y;
    e_out = sum_out;
#endif
}

namespace {

/** Verify that the visible complex prefix contains only finite float32 samples. */
bool scan_prefix_finite(const std::complex<float>* y, size_t n_cplx)
{
    if (!y || n_cplx == 0) {
        return true;
    }
    const float* data = reinterpret_cast<const float*>(y);
    const size_t count = 2 * n_cplx;
    size_t i = 0;
#if defined(__AVX2__)
    const __m256i exp_mask = _mm256_set1_epi32(0x7F800000);
    __m256i or_acc = _mm256_setzero_si256();
    const size_t limit32 = (count >= 32) ? (count - 31) : 0;
    for (; i < limit32; i += 32) {
        const __m256i v0 = _mm256_castps_si256(_mm256_loadu_ps(data + i));
        const __m256i v1 = _mm256_castps_si256(_mm256_loadu_ps(data + i + 8));
        const __m256i v2 = _mm256_castps_si256(_mm256_loadu_ps(data + i + 16));
        const __m256i v3 = _mm256_castps_si256(_mm256_loadu_ps(data + i + 24));
        const __m256i m0 = _mm256_cmpeq_epi32(_mm256_and_si256(v0, exp_mask), exp_mask);
        const __m256i m1 = _mm256_cmpeq_epi32(_mm256_and_si256(v1, exp_mask), exp_mask);
        const __m256i m2 = _mm256_cmpeq_epi32(_mm256_and_si256(v2, exp_mask), exp_mask);
        const __m256i m3 = _mm256_cmpeq_epi32(_mm256_and_si256(v3, exp_mask), exp_mask);
        or_acc = _mm256_or_si256(
            or_acc, _mm256_or_si256(_mm256_or_si256(m0, m1), _mm256_or_si256(m2, m3)));
    }
    const size_t limit8 = (count >= 8) ? (count - 7) : 0;
    for (; i < limit8; i += 8) {
        const __m256i v = _mm256_castps_si256(_mm256_loadu_ps(data + i));
        or_acc = _mm256_or_si256(
            or_acc, _mm256_cmpeq_epi32(_mm256_and_si256(v, exp_mask), exp_mask));
    }
    if (!_mm256_testz_si256(or_acc, or_acc)) {
        return false;
    }
#endif
    for (; i < count; ++i) {
        if (!std::isfinite(data[i])) {
            return false;
        }
    }
    return true;
}

#if defined(__AVX2__) && defined(__FMA__)
inline double hsum256_ps(__m256 v)
{
    const __m128 lo = _mm256_castps256_ps128(v);
    const __m128 hi = _mm256_extractf128_ps(v, 1);
    const __m128 s4 = _mm_add_ps(lo, hi);
    const __m128 s2 = _mm_add_ps(s4, _mm_movehl_ps(s4, s4));
    const __m128 s1 = _mm_add_ss(s2, _mm_shuffle_ps(s2, s2, 0x55));
    return static_cast<double>(_mm_cvtss_f32(s1));
}
#endif

/** Accumulate dst[k] += (ar + j*ai) * pk[k] for interleaved complex floats. */
inline void accum_rotated_symbol(
    const float* __restrict pk,
    float ar,
    float ai,
    float* __restrict dst,
    size_t n_cplx)
{
    size_t i = 0;
#if defined(__AVX2__) && defined(__FMA__)
    const __m256 var = _mm256_set1_ps(ar);
    const __m256 vai = _mm256_setr_ps(-ai, ai, -ai, ai, -ai, ai, -ai, ai);
    for (; i + 8 <= n_cplx; i += 8) {
        const __m256 x0 = _mm256_loadu_ps(pk + 2 * i);
        const __m256 x1 = _mm256_loadu_ps(pk + 2 * i + 8);
        const __m256 d0 = _mm256_loadu_ps(dst + 2 * i);
        const __m256 d1 = _mm256_loadu_ps(dst + 2 * i + 8);
        const __m256 x0s = _mm256_permute_ps(x0, 0xB1);
        const __m256 x1s = _mm256_permute_ps(x1, 0xB1);
        const __m256 r0 = _mm256_fmadd_ps(x0s, vai, _mm256_mul_ps(x0, var));
        const __m256 r1 = _mm256_fmadd_ps(x1s, vai, _mm256_mul_ps(x1, var));
        _mm256_storeu_ps(dst + 2 * i, _mm256_add_ps(d0, r0));
        _mm256_storeu_ps(dst + 2 * i + 8, _mm256_add_ps(d1, r1));
    }
    if (i + 4 <= n_cplx) {
        const __m256 x0 = _mm256_loadu_ps(pk + 2 * i);
        const __m256 d0 = _mm256_loadu_ps(dst + 2 * i);
        const __m256 x0s = _mm256_permute_ps(x0, 0xB1);
        const __m256 r0 = _mm256_fmadd_ps(x0s, vai, _mm256_mul_ps(x0, var));
        _mm256_storeu_ps(dst + 2 * i, _mm256_add_ps(d0, r0));
        i += 4;
    }
#endif
    for (; i < n_cplx; ++i) {
        const float xr = pk[2 * i];
        const float xi = pk[2 * i + 1];
        dst[2 * i] += xr * ar - xi * ai;
        dst[2 * i + 1] += xi * ar + xr * ai;
    }
}

/** Compute sum_i a[i] * conj(b[i]) over interleaved complex floats. */
inline void complex_dot_conj(
    const float* __restrict a,
    const float* __restrict b,
    size_t n_cplx,
    double& dre_out,
    double& dim_out)
{
    size_t i = 0;
    double dre = 0.0;
    double dim = 0.0;
#if defined(__AVX2__) && defined(__FMA__)
    const __m256 sign_even = _mm256_setr_ps(-0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f);
    __m256 re0 = _mm256_setzero_ps();
    __m256 re1 = _mm256_setzero_ps();
    __m256 im0 = _mm256_setzero_ps();
    __m256 im1 = _mm256_setzero_ps();
    for (; i + 8 <= n_cplx; i += 8) {
        const __m256 va0 = _mm256_loadu_ps(a + 2 * i);
        const __m256 vb0 = _mm256_loadu_ps(b + 2 * i);
        const __m256 va1 = _mm256_loadu_ps(a + 2 * i + 8);
        const __m256 vb1 = _mm256_loadu_ps(b + 2 * i + 8);
        const __m256 vb0r = _mm256_xor_ps(_mm256_permute_ps(vb0, 0xB1), sign_even);
        const __m256 vb1r = _mm256_xor_ps(_mm256_permute_ps(vb1, 0xB1), sign_even);
        re0 = _mm256_fmadd_ps(va0, vb0, re0);
        re1 = _mm256_fmadd_ps(va1, vb1, re1);
        im0 = _mm256_fmadd_ps(va0, vb0r, im0);
        im1 = _mm256_fmadd_ps(va1, vb1r, im1);
    }
    if (i + 4 <= n_cplx) {
        const __m256 va0 = _mm256_loadu_ps(a + 2 * i);
        const __m256 vb0 = _mm256_loadu_ps(b + 2 * i);
        const __m256 vb0r = _mm256_xor_ps(_mm256_permute_ps(vb0, 0xB1), sign_even);
        re0 = _mm256_fmadd_ps(va0, vb0, re0);
        im0 = _mm256_fmadd_ps(va0, vb0r, im0);
        i += 4;
    }
    dre = hsum256_ps(_mm256_add_ps(re0, re1));
    dim = hsum256_ps(_mm256_add_ps(im0, im1));
#endif
    for (; i < n_cplx; ++i) {
        const double ar = static_cast<double>(a[2 * i]);
        const double ai = static_cast<double>(a[2 * i + 1]);
        const double br = static_cast<double>(b[2 * i]);
        const double bi = static_cast<double>(b[2 * i + 1]);
        dre += ar * br + ai * bi;
        dim += ai * br - ar * bi;
    }
    dre_out = dre;
    dim_out = dim;
}

/** Measure ||y||^2 and sum y * conj(s_ref) on a symbol segment in one AVX2+FMA pass. */
inline void measure_symbol_energy_and_dot(
    const float* pyf,
    const float* __restrict psf,
    size_t len,
    double& ey_out,
    double& dot_re_out,
    double& dot_im_out)
{
    size_t j = 0;
    double ey = 0.0;
    double dre = 0.0;
    double dim = 0.0;
#if defined(__AVX2__) && defined(__FMA__)
    const __m256 sign_even = _mm256_setr_ps(-0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f);
    __m256 ey0 = _mm256_setzero_ps();
    __m256 ey1 = _mm256_setzero_ps();
    __m256 re0 = _mm256_setzero_ps();
    __m256 re1 = _mm256_setzero_ps();
    __m256 im0 = _mm256_setzero_ps();
    __m256 im1 = _mm256_setzero_ps();
    for (; j + 8 <= len; j += 8) {
        const __m256 vy0 = _mm256_loadu_ps(pyf + 2 * j);
        const __m256 vs0 = _mm256_loadu_ps(psf + 2 * j);
        const __m256 vy1 = _mm256_loadu_ps(pyf + 2 * j + 8);
        const __m256 vs1 = _mm256_loadu_ps(psf + 2 * j + 8);
        const __m256 vs0r = _mm256_xor_ps(_mm256_permute_ps(vs0, 0xB1), sign_even);
        const __m256 vs1r = _mm256_xor_ps(_mm256_permute_ps(vs1, 0xB1), sign_even);
        ey0 = _mm256_fmadd_ps(vy0, vy0, ey0);
        ey1 = _mm256_fmadd_ps(vy1, vy1, ey1);
        re0 = _mm256_fmadd_ps(vy0, vs0, re0);
        re1 = _mm256_fmadd_ps(vy1, vs1, re1);
        im0 = _mm256_fmadd_ps(vy0, vs0r, im0);
        im1 = _mm256_fmadd_ps(vy1, vs1r, im1);
    }
    if (j + 4 <= len) {
        const __m256 vy0 = _mm256_loadu_ps(pyf + 2 * j);
        const __m256 vs0 = _mm256_loadu_ps(psf + 2 * j);
        const __m256 vs0r = _mm256_xor_ps(_mm256_permute_ps(vs0, 0xB1), sign_even);
        ey0 = _mm256_fmadd_ps(vy0, vy0, ey0);
        re0 = _mm256_fmadd_ps(vy0, vs0, re0);
        im0 = _mm256_fmadd_ps(vy0, vs0r, im0);
        j += 4;
    }
    ey = hsum256_ps(_mm256_add_ps(ey0, ey1));
    dre = hsum256_ps(_mm256_add_ps(re0, re1));
    dim = hsum256_ps(_mm256_add_ps(im0, im1));
#endif
    for (; j < len; ++j) {
        const double yr = static_cast<double>(pyf[2 * j]);
        const double yi = static_cast<double>(pyf[2 * j + 1]);
        const double sr = static_cast<double>(psf[2 * j]);
        const double si = static_cast<double>(psf[2 * j + 1]);
        ey += yr * yr + yi * yi;
        dre += yr * sr + yi * si;
        dim += yi * sr - yr * si;
    }
    ey_out = ey;
    dot_re_out = dre;
    dot_im_out = dim;
}

/** Subtract (kr + j*ki) * s_ref from y, optionally emit i_hat, and measure ||out||^2. */
template <bool WriteIHat>
inline double subtract_scaled_symbol_impl(
    const float* pyf,
    const float* __restrict psf,
    float kr,
    float ki,
    size_t len,
    float* poutf,
    float* __restrict pihatf)
{
    size_t j = 0;
    double eout = 0.0;
#if defined(__AVX2__) && defined(__FMA__)
    const __m256 vkr = _mm256_set1_ps(kr);
    const __m256 vki = _mm256_setr_ps(-ki, ki, -ki, ki, -ki, ki, -ki, ki);
    __m256 eo0 = _mm256_setzero_ps();
    __m256 eo1 = _mm256_setzero_ps();
    for (; j + 8 <= len; j += 8) {
        const __m256 vs0 = _mm256_loadu_ps(psf + 2 * j);
        const __m256 vs1 = _mm256_loadu_ps(psf + 2 * j + 8);
        const __m256 vy0 = _mm256_loadu_ps(pyf + 2 * j);
        const __m256 vy1 = _mm256_loadu_ps(pyf + 2 * j + 8);
        const __m256 vs0s = _mm256_permute_ps(vs0, 0xB1);
        const __m256 vs1s = _mm256_permute_ps(vs1, 0xB1);
        const __m256 vi0 = _mm256_fmadd_ps(vs0s, vki, _mm256_mul_ps(vs0, vkr));
        const __m256 vi1 = _mm256_fmadd_ps(vs1s, vki, _mm256_mul_ps(vs1, vkr));
        const __m256 vo0 = _mm256_sub_ps(vy0, vi0);
        const __m256 vo1 = _mm256_sub_ps(vy1, vi1);
        eo0 = _mm256_fmadd_ps(vo0, vo0, eo0);
        eo1 = _mm256_fmadd_ps(vo1, vo1, eo1);
        if constexpr (WriteIHat) {
            _mm256_storeu_ps(pihatf + 2 * j, vi0);
            _mm256_storeu_ps(pihatf + 2 * j + 8, vi1);
        }
        _mm256_storeu_ps(poutf + 2 * j, vo0);
        _mm256_storeu_ps(poutf + 2 * j + 8, vo1);
    }
    if (j + 4 <= len) {
        const __m256 vs0 = _mm256_loadu_ps(psf + 2 * j);
        const __m256 vy0 = _mm256_loadu_ps(pyf + 2 * j);
        const __m256 vs0s = _mm256_permute_ps(vs0, 0xB1);
        const __m256 vi0 = _mm256_fmadd_ps(vs0s, vki, _mm256_mul_ps(vs0, vkr));
        const __m256 vo0 = _mm256_sub_ps(vy0, vi0);
        eo0 = _mm256_fmadd_ps(vo0, vo0, eo0);
        if constexpr (WriteIHat) {
            _mm256_storeu_ps(pihatf + 2 * j, vi0);
        }
        _mm256_storeu_ps(poutf + 2 * j, vo0);
        j += 4;
    }
    eout = hsum256_ps(_mm256_add_ps(eo0, eo1));
#endif
    for (; j < len; ++j) {
        const float sr = psf[2 * j];
        const float si = psf[2 * j + 1];
        const float ir = kr * sr - ki * si;
        const float ii = kr * si + ki * sr;
        const float or_ = pyf[2 * j] - ir;
        const float oi_ = pyf[2 * j + 1] - ii;
        if constexpr (WriteIHat) {
            pihatf[2 * j] = ir;
            pihatf[2 * j + 1] = ii;
        }
        poutf[2 * j] = or_;
        poutf[2 * j + 1] = oi_;
        eout += static_cast<double>(or_) * static_cast<double>(or_) +
                static_cast<double>(oi_) * static_cast<double>(oi_);
    }
    return eout;
}

// Empirical steady-state SYNC template + long-baseline CFO + per-symbol phase
// tracking for L1 Preamble+SFD cancellation.  Bypasses pulse-shaping model
// mismatch (~1.3 dB) and suppresses both startup TX turn-on phase/amplitude
// transients (k < fit_skip_reps) and symbol-to-symbol phase noise while
// preserving the zero-allocation, fail-closed, in-place contract.
bool reconstruct_and_subtract_empirical_l1(
    const std::complex<float>* y,
    size_t n,
    int64_t obs,
    const SicUnifiedChannel& channel,
    const SicUnifiedConfig& cfg,
    const SicUnifiedScratch& scratch,
    std::complex<float>* out,
    std::complex<float>* i_hat,
    SicUnifiedResult& result)
{
    using cf = std::complex<float>;
    using cd = std::complex<double>;
    constexpr size_t kSps = 1016;
    constexpr double kTwoPi = 6.283185307179586476925286766559;

    if (!y || !out || n == 0 || !channel.valid || cfg.samples_per_symbol != kSps ||
        channel.q < 0 || cfg.sample_rate <= 0.0 || cfg.sync_reps == 0) {
        return false;
    }

    const int64_t q = channel.q;
    const int64_t vis_end = (obs > 0) ? std::min<int64_t>(int64_t(n), obs) : int64_t(n);
    if (vis_end <= q) {
        return false;
    }
    const size_t avail_sync = std::min<size_t>(
        cfg.sync_reps, static_cast<size_t>((vis_end - q) / static_cast<int64_t>(kSps)));
    const size_t k_skip = std::max<size_t>(1, cfg.fit_skip_reps);
    const size_t k_train_end = std::min(avail_sync, cfg.fit_reps);
    if (k_train_end < k_skip + 4) {
        return false;
    }

    const int64_t roi_lo = std::max<int64_t>(0, cfg.sub_roi_lo);
    const int64_t roi_hi = (cfg.sub_roi_hi <= 0)
        ? static_cast<int64_t>(n)
        : std::min<int64_t>(static_cast<int64_t>(n), cfg.sub_roi_hi);
    const size_t sfd_len = scratch.cached_sfd_len;
    const size_t total_syms = cfg.sync_reps + sfd_len;
    const int64_t supp_lo = q;
    const int64_t supp_hi = q + static_cast<int64_t>(total_syms * kSps);
    const int64_t sub_lo = std::max<int64_t>(roi_lo, std::max<int64_t>(0, supp_lo));
    const int64_t sub_hi = std::min<int64_t>(roi_hi, std::min<int64_t>(static_cast<int64_t>(n), supp_hi));
    if (sub_hi <= sub_lo) {
        return false;
    }
    if (!cfg.allow_partial_window_clip && (supp_lo < 0 || supp_hi > static_cast<int64_t>(n))) {
        return false;
    }

    const auto t_recon0 = Clock::now();

    // 1. Steady-state residual CFO refinement using two 16-symbol coherent blocks
    // separated across the visible steady-state SYNC span.  Averaging 16 symbols
    // at each end suppresses per-symbol phase jitter by 4x without endpoint
    // telescoping and costs only 32 symbol accumulations.
    double cfo_hz = channel.cfo_hz;
    if (avail_sync >= k_skip + 32) {
        constexpr size_t kBlk = 16;
        const size_t blk0 = k_skip;
        const size_t blk1 = avail_sync - kBlk;
        const double omega0 = kTwoPi * channel.cfo_hz / cfg.sample_rate;
        alignas(32) std::array<float, 2 * kSps> b0{}, b1{};
        auto accum_block = [&](size_t start_k, float* __restrict dst) {
            for (size_t k = start_k; k < start_k + kBlk; ++k) {
                const double ph = -omega0 * static_cast<double>(k * kSps);
                const float ar = static_cast<float>(std::cos(ph));
                const float ai = static_cast<float>(std::sin(ph));
                const float* __restrict pk = reinterpret_cast<const float*>(y + q + k * kSps);
                accum_rotated_symbol(pk, ar, ai, dst, kSps);
            }
        };
        accum_block(blk0, b0.data());
        accum_block(blk1, b1.data());
        double dre = 0.0, dim = 0.0;
        complex_dot_conj(b1.data(), b0.data(), kSps, dre, dim);
        const cd cross(dre, dim);
        if (std::norm(cross) > 1e-24) {
            const double dt_s = static_cast<double>((blk1 - blk0) * kSps) / cfg.sample_rate;
            const double delta_cfo = std::arg(cross) / (kTwoPi * dt_s);
            if (std::isfinite(delta_cfo) && std::abs(delta_cfo) < 5000.0) {
                cfo_hz += delta_cfo;
            }
        }
    }

    const double omega = kTwoPi * cfo_hz / cfg.sample_rate;

    // 2. Build 1016-sample steady-state template s_bar[0..1015] * exp(+j*w*i).
    alignas(32) std::array<float, 2 * kSps> s_acc{};
    for (size_t k = k_skip; k < k_train_end; ++k) {
        const double ph = -omega * static_cast<double>(k * kSps);
        const float ar = static_cast<float>(std::cos(ph));
        const float ai = static_cast<float>(std::sin(ph));
        const float* __restrict pk = reinterpret_cast<const float*>(y + q + k * kSps);
        accum_rotated_symbol(pk, ar, ai, s_acc.data(), kSps);
    }
    const float inv_cnt_f = static_cast<float>(1.0 / static_cast<double>(k_train_end - k_skip));
    double e_ref = 0.0;
    alignas(32) std::array<cf, kSps> s_ref_cf{};
    for (size_t i = 0; i < kSps; ++i) {
        const float vr = s_acc[2 * i] * inv_cnt_f;
        const float vi = s_acc[2 * i + 1] * inv_cnt_f;
        s_ref_cf[i] = cf(vr, vi);
        e_ref += static_cast<double>(vr) * vr + static_cast<double>(vi) * vi;
    }
    if (!std::isfinite(e_ref) || e_ref <= 0.0) {
        return false;
    }

    const auto t_sub0 = Clock::now();

    if (out != y) {
        if (sub_lo > 0) {
            std::copy(y, y + sub_lo, out);
        }
        if (sub_hi < static_cast<int64_t>(n)) {
            std::copy(y + sub_hi, y + n, out + sub_hi);
        }
    }
    if (i_hat && i_hat != y && i_hat != out) {
        std::fill(i_hat, i_hat + n, cf(0.f, 0.f));
    }

    double e_y_tot = 0.0;
    double e_out_tot = 0.0;
    bool any_subtracted = false;

    const double csym = std::cos(omega * static_cast<double>(kSps));
    const double ssym = std::sin(omega * static_cast<double>(kSps));
    double sym_r = 1.0;
    double sym_i = 0.0;

    for (size_t k = 0; k < total_syms; ++k) {
        const cd sym_rot_base(sym_r, sym_i);
        const double nr = sym_r * csym - sym_i * ssym;
        const double ni = sym_r * ssym + sym_i * csym;
        sym_r = nr;
        sym_i = ni;

        const int8_t sgn = (k < cfg.sync_reps) ? int8_t(1) : scratch.cached_sfd_signs[k - cfg.sync_reps];
        const int64_t sym_lo = q + static_cast<int64_t>(k * kSps);
        const int64_t sym_hi = sym_lo + static_cast<int64_t>(kSps);
        const int64_t a = std::max(sym_lo, sub_lo);
        const int64_t b = std::min(sym_hi, sub_hi);
        if (b <= a) {
            continue;
        }
        const size_t off0 = static_cast<size_t>(a - sym_lo);
        const size_t len = static_cast<size_t>(b - a);

        if (sgn == 0) {
            double ey = 0.0;
            for (size_t j = 0; j < len; ++j) {
                const cf yv = y[a + j];
                if (out != y) {
                    out[a + j] = yv;
                }
                ey += std::norm(yv);
            }
            e_y_tot += ey;
            e_out_tot += ey;
            continue;
        }

        const cd sym_rot = sym_rot_base * static_cast<double>(sgn);

        double ey_sym = 0.0;
        double dot_re = 0.0;
        double dot_im = 0.0;
        const float* pyf = reinterpret_cast<const float*>(y + a);
        const float* __restrict psf = reinterpret_cast<const float*>(s_ref_cf.data() + off0);
        measure_symbol_energy_and_dot(pyf, psf, len, ey_sym, dot_re, dot_im);

        const double e_seg = (len == kSps) ? e_ref : [&]() {
            double es = 0.0;
            for (size_t j = 0; j < len; ++j) {
                es += std::norm(s_ref_cf[off0 + j]);
            }
            return es;
        }();

        const cd dot_raw(dot_re, dot_im);
        const cd dot = dot_raw * std::conj(sym_rot);
        cd beta(1.0, 0.0);
        if (len >= 256 && e_seg > 0.0) {
            if (k < k_skip) {
                beta = dot / e_seg;
                const double mag = std::abs(beta);
                if (mag > 1.5) {
                    beta *= (1.5 / mag);
                }
            } else {
                const double mag = std::abs(dot);
                if (mag > 0.0) {
                    beta = dot / mag;
                }
            }
        }

        const double delta_e = -2.0 * (beta.real() * dot.real() + beta.imag() * dot.imag()) +
                               std::norm(beta) * e_seg;
        const cd scale_cd = beta * sym_rot;
        if (!std::isfinite(delta_e) || !std::isfinite(scale_cd.real()) ||
            !std::isfinite(scale_cd.imag()) || delta_e >= 0.0) {
            if (out != y) {
                std::copy(y + a, y + a + len, out + a);
            }
            e_y_tot += ey_sym;
            e_out_tot += ey_sym;
            continue;
        }

        const float kr = static_cast<float>(scale_cd.real());
        const float ki = static_cast<float>(scale_cd.imag());
        float* poutf = reinterpret_cast<float*>(out + a);
        const double eout_sym = i_hat
            ? subtract_scaled_symbol_impl<true>(
                  pyf, psf, kr, ki, len, poutf, reinterpret_cast<float*>(i_hat + a))
            : subtract_scaled_symbol_impl<false>(
                  pyf, psf, kr, ki, len, poutf, nullptr);
        e_y_tot += ey_sym;
        e_out_tot += eout_sym;
        any_subtracted = true;
    }

    if (!any_subtracted || !std::isfinite(e_y_tot) || !std::isfinite(e_out_tot) ||
        !(e_y_tot > 0.0) || !(e_out_tot > 0.0)) {
        return false;
    }

    result.committed = true;
    result.reason = "committed_preamble_empirical_phase";
    result.level = SicCancellationLevel::L1_Preamble;
    result.cancelled_fields = SicKnownFields::PreambleSfd;
    result.support_lo = supp_lo;
    result.support_hi = supp_hi;
    result.sub_lo = sub_lo;
    result.sub_hi = sub_hi;
    result.sub_overlap = sub_hi - sub_lo;
    result.waveform_scale = 1.0;
    result.suppression_energy_valid = true;
    result.suppression_db = 10.0 * std::log10(e_y_tot / e_out_tot);
    result.channel.cfo_hz = cfo_hz;
    result.timings.reconstruct_us = std::chrono::duration<double, std::micro>(t_sub0 - t_recon0).count();
    result.timings.subtract_us = us(t_sub0);
    return true;
}

} // namespace

// ---------------------------------------------------------------------------
// sic_process_unified
// ---------------------------------------------------------------------------

bool sic_process_unified(
    const std::complex<float>* y,
    size_t n,
    const SicUnifiedTemplates& templates,
    const SicUnifiedConfig& cfg,
    SicUnifiedScratch& scratch,
    SicUnifiedResult& result)
{
    const auto t_core_start = Clock::now();
    result.policy = cfg.policy;
    result.committed = false;
    result.level = SicCancellationLevel::L0_Bypass;
    result.cancelled_fields = SicKnownFields::None;
    result.suppression_db = 0.0;
    result.suppression_energy_valid = false;
    result.fallback_reason = nullptr;

    if (!y || !result.out || n == 0 || !scratch.prepared ||
        !scratch.cache_matches(cfg, templates)) {
        result.reason = (!y || !result.out || n == 0 || !scratch.prepared)
                            ? "invalid_buffers_or_scratch"
                            : "config_or_template_mismatch";
        if (result.copy_out_on_bypass && result.out && y && n > 0 &&
            result.out != y) {
            std::copy(y, y + n, result.out);
        }
        if (result.i_hat) {
            std::fill(result.i_hat, result.i_hat + n, std::complex<float>(0.f, 0.f));
        }
        return false;
    }

    auto fail_closed_before_recon = [&](const char* reason) {
        if (result.copy_out_on_bypass && result.out != y) {
            std::copy(y, y + n, result.out);
        }
        if (result.i_hat) {
            std::fill(result.i_hat, result.i_hat + n, std::complex<float>(0.f, 0.f));
        }
        result.reason = reason;
        result.level = SicCancellationLevel::L0_Bypass;
        result.cancelled_fields = SicKnownFields::None;
        result.timings.core_total_us = us(t_core_start);
        return true;
    };

    const int64_t obs = (cfg.observed_end_sample <= 0)
                            ? static_cast<int64_t>(n)
                            : std::min<int64_t>(static_cast<int64_t>(n), cfg.observed_end_sample);

    SicCancellationLevel target_level = SicCancellationLevel::L0_Bypass;
    SicKnownFields target_fields = SicKnownFields::None;
    const std::complex<float>* s_known_ptr = templates.preamble;
    size_t s_known_len = templates.preamble_len;
    const std::vector<std::complex<float>>& sync_vec = scratch.cached_sync_vec;

    auto demod_prof = demod::Qm35825Profile::Default();
    demod_prof.code_index = cfg.comm_code;
    demod_prof.preamble_repetitions = cfg.sync_reps;
    if (cfg.sync_reps > demod_prof.cir_skip_initial_repetitions)
        demod_prof.cir_repetitions = cfg.sync_reps - demod_prof.cir_skip_initial_repetitions;
    demod_prof.sfd_mode = cfg.sfd_mode;
    demod_prof.rx_fast_sfd = cfg.rx_fast_sfd;
    demod_prof.rx_fast_cfo_rotation = cfg.rx_fast_cfo;
    demod_prof.rx_seed_refine_tail = cfg.rx_seed_refine_tail;
    demod_prof.rx_softchip_on_demand = cfg.rx_softchip_on_demand;
    demod_prof.rx_fused_payload_cmf = cfg.rx_fused_payload_cmf;
    demod_prof.rx_fused_phr_cmf = cfg.rx_fused_phr_cmf;
    demod_prof.rx_timing_sparse = cfg.rx_timing_sparse;
    demod_prof.rx_timing_shaped_sparse = cfg.rx_timing_shaped_sparse;
    demod_prof.cfo_refit_mode = cfg.cfo_refit_mode;
    demod_prof.rx_wide_interval_probe = cfg.rx_wide_interval_probe;
    demod_prof.rx_wide_probe_step_symbols = cfg.rx_wide_probe_step_symbols;
    demod_prof.rx_wide_probe_fold_symbols = cfg.rx_wide_probe_fold_symbols;
    demod_prof.rx_wide_probe_fallback = cfg.rx_wide_probe_fallback;

    SicUnifiedChannel channel;

    bool empirical_done = false;
    if (cfg.backend == SicEstimatorBackend::Ipatov) {
        // Unified Ipatov path: run standard demodulate_one FIRST for blind
        // timing (AVX-512 16-lag), CFO, SFD and optional PHR/payload, then
        // run only the strict-prefix Ipatov CIR + physical residual gate.
        int64_t predicted_start = -1;
        if (cfg.roi_hi > cfg.roi_lo) {
            predicted_start = cfg.roi_lo + (cfg.roi_hi - cfg.roi_lo) / 2;
            const size_t half_roi = static_cast<size_t>((cfg.roi_hi - cfg.roi_lo + 1) / 2);
            demod_prof.timing_search_margin = std::max<size_t>(256, half_roi + 64);
        }

        const bool preamble_only = (cfg.policy == SicPolicy::PreambleOnly);
        demod_prof.rx_stop_after_sfd = preamble_only;
        // A blind / whole-window SIC search can acquire a later SYNC of a long
        // interferer preamble (the stride-14 coarse scan is intentionally
        // cheap).  The default 40-symbol backtrack cannot reach the true first
        // SYNC of a 128/256-SYNC DW1000-style preamble, which would then place
        // the SFD 128 symbols too late and fail detection.  The SIC knows the
        // configured repetition count, so allow the walk-back to cover it.
        demod_prof.timing_max_backtrack_symbols = std::max(
            demod_prof.timing_max_backtrack_symbols,
            static_cast<size_t>(cfg.sync_reps));
        const size_t demod_len = preamble_only ? static_cast<size_t>(obs) : n;

        const auto t_rx = Clock::now();
        demod::DemodResult decoded = demod::core::demodulate_one(
            y, demod_len, demod_prof, 0, predicted_start, 0, sync_vec, scratch.demod);
        result.timings.comm_rx_us = us(t_rx);
        result.timings.rx_timing_us = static_cast<double>(decoded.stage_timing_us);
        result.timings.rx_cfo_us = static_cast<double>(decoded.stage_cfo_us);
        result.timings.rx_sfd_us = static_cast<double>(decoded.stage_sfd_us);
        if (!preamble_only) {
            result.timings.rx_cir_us = static_cast<double>(decoded.stage_cir_us);
            result.timings.rx_phr_us = static_cast<double>(decoded.stage_phr_us);
            result.timings.rx_payload_us = static_cast<double>(decoded.stage_payload_us);
        }

        result.comm_detected = (decoded.status != demod::DemodStatus::TimingFailed &&
                                decoded.status != demod::DemodStatus::InvalidInput);
        const bool sfd_ok = decoded.sfd.ok && result.comm_detected;
        const int64_t q_det = decoded.timing.preamble_start_sample;
        const size_t sfd_syms = (scratch.cached_sfd_len > 0) ? scratch.cached_sfd_len : 8;
        const size_t grid_reps = std::max(decoded.timing.detected_peaks,
                                          decoded.timing.cir_grid_repetitions);
        const bool clipped_sfd_l1_ok =
            cfg.allow_partial_window_clip && result.comm_detected &&
            decoded.timing.ok && !sfd_ok && q_det >= 0 && decoded.cfo.ok &&
            decoded.timing.metric >= cfg.detect_threshold &&
            grid_reps >= std::min<size_t>(16, cfg.sync_reps / 4) &&
            (q_det + static_cast<int64_t>((cfg.sync_reps + sfd_syms) * cfg.samples_per_symbol) >
             static_cast<int64_t>(demod_len));

        if (preamble_only) {
            if (!((decoded.status == demod::DemodStatus::Success && decoded.sfd.ok) ||
                  clipped_sfd_l1_ok)) {
                return fail_closed_before_recon("sfd_detect_failed");
            }
            target_level = SicCancellationLevel::L1_Preamble;
            target_fields = SicKnownFields::PreambleSfd;
        } else {
            result.comm_phr_pass = (decoded.phr.psdu_length > 0 && !decoded.phr.secded_uncorrectable);
            result.comm_fcs_pass = decoded.payload.fcs_pass;
            result.phr_secded_corrected = decoded.phr.secded_corrected;
            result.phr_secded_uncorrectable = decoded.phr.secded_uncorrectable;
            result.comm_phr_length = decoded.phr.psdu_length;
            result.comm_phr_rate_mbps = decoded.phr.data_rate_mbps;

            if (decoded.payload.bytes.size() <= 127) {
                result.wire_psdu_len = decoded.payload.bytes.size();
                std::copy(decoded.payload.bytes.begin(), decoded.payload.bytes.end(), result.wire_psdu);
                scratch.work_wire_psdu = decoded.payload.bytes;
            }
            scratch.work_phr_bits = decoded.phr.phr_bits;

            const bool demod_success = (decoded.status == demod::DemodStatus::Success &&
                                        decoded.payload.fcs_pass &&
                                        decoded.phr.phr_bits.size() == 19);
            if (demod_success) {
                target_level = SicCancellationLevel::L2_FullPacket;
                target_fields = SicKnownFields::FullPacket;
            } else {
                const char* fb = !result.comm_detected ? "comm_detect_failed"
                               : (!result.comm_phr_pass ? "comm_phr_failed"
                               : (!result.comm_fcs_pass ? "comm_fcs_failed"
                               : "comm_demod_failed"));
                if (cfg.policy == SicPolicy::FullPacket || (!sfd_ok && !clipped_sfd_l1_ok)) {
                    return fail_closed_before_recon(fb);
                }
                result.fallback_reason = fb;
                target_level = SicCancellationLevel::L1_Preamble;
                target_fields = SicKnownFields::PreambleSfd;
            }
        }

        const bool try_empirical_l1 =
            cfg.enable_empirical_phase_tracking &&
            target_level == SicCancellationLevel::L1_Preamble &&
            cfg.samples_per_symbol == 1016;

        // Empirical L1 cancels directly from the steady-state per-symbol SYNC
        // stack without building a sparse discrete CIR or theoretical fit_check.
        // Fall through to full discrete CIR estimation if preconditions or
        // empirical subtraction fail.
        if (try_empirical_l1 && q_det >= 0 && obs > q_det && std::isfinite(decoded.cfo.cfo_hz)) {
            const auto t_est = Clock::now();
            if (!scan_prefix_finite(y, static_cast<size_t>(obs))) {
                result.timings.estimate_us = us(t_est);
                return fail_closed_before_recon("nonfinite_input");
            }
            const size_t avail_sync = std::min<size_t>(
                cfg.sync_reps, static_cast<size_t>((obs - q_det) / 1016));
            const size_t k_skip = std::max<size_t>(1, cfg.fit_skip_reps);
            const size_t k_train_end = std::min(avail_sync, cfg.fit_reps);
            channel.valid = true;
            channel.backend = cfg.backend;
            channel.status = "estimated";
            channel.q = q_det;
            channel.cfo_hz = decoded.cfo.cfo_hz;
            channel.phase_reference_sample = q_det;
            channel.observed_end_sample = obs;
            channel.impulse_template = true;
            channel.repeats_confirmed = avail_sync;
            if (k_train_end > k_skip) {
                channel.training_lo = q_det + static_cast<int64_t>(k_skip * 1016);
                channel.training_hi = q_det + static_cast<int64_t>(k_train_end * 1016);
                channel.cir_repetitions = k_train_end - k_skip;
            }
            result.channel = channel;
            result.timings.estimate_us = us(t_est);

            empirical_done = reconstruct_and_subtract_empirical_l1(
                y, n, obs, channel, cfg, scratch, result.out, result.i_hat, result);
        }

        if (!empirical_done) {
            if (target_level == SicCancellationLevel::L1_Preamble) {
                if (scratch.preamble_cache_valid && !scratch.cached_preamble_waveform.empty()) {
                    s_known_ptr = scratch.cached_preamble_waveform.data();
                    s_known_len = scratch.cached_preamble_waveform.size();
                } else {
                    mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
                    size_t n_prefix = 0;
                    if (mod::modulate_prefix_waveform(modcfg, scratch.modulator,
                                                      scratch.work_waveform.data(),
                                                      scratch.work_waveform.size(),
                                                      n_prefix) && n_prefix > 0) {
                        s_known_ptr = scratch.work_waveform.data();
                        s_known_len = n_prefix;
                    }
                }
            }

            const auto t_est = Clock::now();
            DiscreteCirModel dm;
            DiscreteCirTimings dt;
            DiscreteCirConfig dc = cfg.to_discrete_cir_config();
            dc.known_q = decoded.timing.preamble_start_sample;
            dc.known_cfo_hz = decoded.cfo.cfo_hz;
            dc.has_known_cfo = true;
            DiscreteCirTemplates dtemp = templates.to_discrete_cir_templates();
            const bool est_ok = estimate_discrete_cir(
                y, n, obs, dtemp, dc, scratch.discrete, dm, &dt);
            channel.from_discrete_cir(dm);
            result.channel = channel;
            result.timings.estimate_us = us(t_est);
            result.timings.est_detect_us = dt.detect_us;
            result.timings.est_cfo_us = dt.cfo_us;
            result.timings.est_cir_us = dt.support_us;
            result.timings.est_solve_us = dt.projection_us + dt.solve_us;
            result.timings.est_fit_check_us = dt.fit_check_us;

            if (!est_ok || !channel.valid) {
                return fail_closed_before_recon(channel.reason ? channel.reason : channel.status);
            }
        }
    } else {
        // 1. Strict-prefix Channel Estimation (DiscreteCir / PathFractional)
        const bool est_ok = estimate_channel_unified(
            y, n, cfg.observed_end_sample, templates, cfg, scratch,
            channel, &result.timings);
        result.channel = channel;

        if (!est_ok || !channel.valid) {
            return fail_closed_before_recon(channel.reason ? channel.reason : channel.status);
        }

        // 2. Cancellation Policy & Comm RX Demodulation
        demod_prof.timing_search_margin = 256;
        demod_prof.timing_coarse_stride = 1;
        demod_prof.timing_seed_refine_radius = 16;

        if (cfg.policy == SicPolicy::PreambleOnly) {
            demod_prof.rx_stop_after_sfd = true;
            const auto t_rx = Clock::now();
            demod::DemodResult decoded = demod::core::demodulate_one(
                y, n, demod_prof, 0, channel.q, 0, sync_vec, scratch.demod);
            result.timings.comm_rx_us = us(t_rx);
            result.timings.rx_timing_us = static_cast<double>(decoded.stage_timing_us);
            result.timings.rx_cfo_us = static_cast<double>(decoded.stage_cfo_us);
            result.timings.rx_sfd_us = static_cast<double>(decoded.stage_sfd_us);

            result.comm_detected = (decoded.status != demod::DemodStatus::TimingFailed &&
                                    decoded.status != demod::DemodStatus::InvalidInput);

            const bool sfd_ok = (decoded.status == demod::DemodStatus::Success && decoded.sfd.ok);
            if (sfd_ok) {
                target_level = SicCancellationLevel::L1_Preamble;
                target_fields = SicKnownFields::PreambleSfd;

                if (scratch.preamble_cache_valid && !scratch.cached_preamble_waveform.empty()) {
                    s_known_ptr = scratch.cached_preamble_waveform.data();
                    s_known_len = scratch.cached_preamble_waveform.size();
                } else {
                    mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
                    size_t n_prefix = 0;
                    if (mod::modulate_prefix_waveform(modcfg, scratch.modulator,
                                                      scratch.work_waveform.data(),
                                                      scratch.work_waveform.size(),
                                                      n_prefix) && n_prefix > 0) {
                        s_known_ptr = scratch.work_waveform.data();
                        s_known_len = n_prefix;
                    }
                }
            } else {
                return fail_closed_before_recon("sfd_detect_failed");
            }
        } else {
            demod_prof.rx_stop_after_sfd = false;
            const auto t_rx = Clock::now();
            demod::DemodResult decoded = demod::core::demodulate_one(
                y, n, demod_prof, 0, channel.q, 0, sync_vec, scratch.demod);
            result.timings.comm_rx_us = us(t_rx);
            result.timings.rx_timing_us = static_cast<double>(decoded.stage_timing_us);
            result.timings.rx_cfo_us = static_cast<double>(decoded.stage_cfo_us);
            result.timings.rx_sfd_us = static_cast<double>(decoded.stage_sfd_us);
            result.timings.rx_cir_us = static_cast<double>(decoded.stage_cir_us);
            result.timings.rx_phr_us = static_cast<double>(decoded.stage_phr_us);
            result.timings.rx_payload_us = static_cast<double>(decoded.stage_payload_us);

            result.comm_detected = (decoded.status != demod::DemodStatus::TimingFailed &&
                                    decoded.status != demod::DemodStatus::InvalidInput);
            result.comm_phr_pass = (decoded.phr.psdu_length > 0 && !decoded.phr.secded_uncorrectable);
            result.comm_fcs_pass = decoded.payload.fcs_pass;
            result.phr_secded_corrected = decoded.phr.secded_corrected;
            result.phr_secded_uncorrectable = decoded.phr.secded_uncorrectable;
            result.comm_phr_length = decoded.phr.psdu_length;
            result.comm_phr_rate_mbps = decoded.phr.data_rate_mbps;

            if (decoded.payload.bytes.size() <= 127) {
                result.wire_psdu_len = decoded.payload.bytes.size();
                std::copy(decoded.payload.bytes.begin(), decoded.payload.bytes.end(), result.wire_psdu);
                scratch.work_wire_psdu = decoded.payload.bytes;
            }
            scratch.work_phr_bits = decoded.phr.phr_bits;

            const bool demod_success = (decoded.status == demod::DemodStatus::Success &&
                                        decoded.payload.fcs_pass &&
                                        decoded.phr.phr_bits.size() == 19);

            if (demod_success) {
                target_level = SicCancellationLevel::L2_FullPacket;
                target_fields = SicKnownFields::FullPacket;
            } else {
                const char* fb = !result.comm_detected ? "comm_detect_failed"
                               : (!result.comm_phr_pass ? "comm_phr_failed"
                               : (!result.comm_fcs_pass ? "comm_fcs_failed"
                               : "comm_demod_failed"));

                if (cfg.policy == SicPolicy::FullPacket) {
                    return fail_closed_before_recon(fb);
                } else {
                    result.fallback_reason = fb;
                    target_level = SicCancellationLevel::L1_Preamble;
                    target_fields = SicKnownFields::PreambleSfd;

                    if (scratch.preamble_cache_valid && !scratch.cached_preamble_waveform.empty()) {
                        s_known_ptr = scratch.cached_preamble_waveform.data();
                        s_known_len = scratch.cached_preamble_waveform.size();
                    } else {
                        mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
                        size_t n_prefix = 0;
                        if (mod::modulate_prefix_waveform(modcfg, scratch.modulator,
                                                          scratch.work_waveform.data(),
                                                          scratch.work_waveform.size(),
                                                          n_prefix) && n_prefix > 0) {
                            s_known_ptr = scratch.work_waveform.data();
                            s_known_len = n_prefix;
                        }
                    }
                }
            }
        }
    }

    // 3. Reconstruction and Subtraction
    bool ok = empirical_done;
    if (!ok && cfg.enable_empirical_phase_tracking &&
        cfg.backend == SicEstimatorBackend::Ipatov &&
        target_level == SicCancellationLevel::L1_Preamble &&
        cfg.samples_per_symbol == 1016) {
        ok = reconstruct_and_subtract_empirical_l1(
            y, n, obs, channel, cfg, scratch, result.out, result.i_hat, result);
    }
    if (!ok) {
        ok = reconstruct_and_subtract_unified(
            y, n, channel, target_level, target_fields,
            s_known_ptr, s_known_len,
            templates, cfg, scratch, result.out, result.i_hat, result);
    }

    // 3b. AutoFallback: a full-packet transaction whose support leaves the
    // window (edge interferer / pulse tail past the last sample) fails closed
    // with "context_truncated" before anything is written.  The preamble is
    // shorter, so retry at L1 with the same known prefix instead of losing the
    // whole window.  Other failures keep the strict L0 contract.
    if (!ok && cfg.policy == SicPolicy::AutoFallback &&
        target_level == SicCancellationLevel::L2_FullPacket && result.reason &&
        std::strcmp(result.reason, "context_truncated") == 0) {
        result.fallback_reason = "l2_context_truncated";
        s_known_ptr = templates.preamble;
        s_known_len = templates.preamble_len;
        if (scratch.preamble_cache_valid && !scratch.cached_preamble_waveform.empty()) {
            s_known_ptr = scratch.cached_preamble_waveform.data();
            s_known_len = scratch.cached_preamble_waveform.size();
        } else {
            mod::HrpModConfig modcfg = cfg.to_hrp_mod_config();
            size_t n_prefix = 0;
            if (mod::modulate_prefix_waveform(modcfg, scratch.modulator,
                                              scratch.work_waveform.data(),
                                              scratch.work_waveform.size(),
                                              n_prefix) && n_prefix > 0) {
                s_known_ptr = scratch.work_waveform.data();
                s_known_len = n_prefix;
            }
        }
        ok = reconstruct_and_subtract_unified(
            y, n, channel, SicCancellationLevel::L1_Preamble,
            SicKnownFields::PreambleSfd, s_known_ptr, s_known_len, templates,
            cfg, scratch, result.out, result.i_hat, result);
    }

    if (!ok) {
        result.timings.core_total_us = us(t_core_start);
        return true;
    }

    // 4. Calculate actual suppression dB in subtraction ROI: 10 * log10(||y||^2 / ||out||^2)
    // When the fused subtraction pass already measured both ROI energies the
    // metric is final: re-reading y here would be wrong for in-place output
    // (out == y) and would add a second full-ROI pass on the hot path.  The
    // post-hoc fallback is therefore restricted to backends that did not fuse
    // the energies (PathFractional / non-sparse L1) AND kept y separate.
    if (result.committed && !result.suppression_energy_valid &&
        result.out != y && result.sub_hi > result.sub_lo &&
        result.sub_lo >= 0 && result.sub_hi <= static_cast<int64_t>(n)) {
        double e_y = 0.0;
        double e_out = 0.0;
        compute_roi_energies(y, result.out, result.sub_lo, result.sub_hi, e_y, e_out);
        if (std::isfinite(e_y) && std::isfinite(e_out) && e_y > 0.0 && e_out > 0.0) {
            result.suppression_db = 10.0 * std::log10(e_y / e_out);
        }
    }

    result.timings.core_total_us = us(t_core_start);
    return true;
}

} // namespace sic
} // namespace openuwb
