/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Folded-preamble cross-correlation acquisition for the native X410 live
 * receiver (openuwb_rx continuous mode).
 *
 * Design contract:
 * docs/OpenCode执行方案_Acquire前端折叠前导码互相关设计与验收_2026-11-24.md §2/§4
 *
 *   - Template: the 64 non-zero chips of the HRP Ipatov code (9..12) map to
 *     the native 750-sample cycle offset
 *         off = llround(8 * c_i * 48 / 65) % 750
 *     so the sparse ternary correlation is add/sub only, zero multiplies.
 *   - Coherent fold of K_coh = 4 SYNC symbols (sinc^2 loss -2.5 dB at
 *     +-100 kHz CFO), non-coherent power accumulation over G_full = 16 groups
 *     (64 symbols, ~65 us).  Stage 1 (G_fast = 6) early-exits silent windows
 *     before the full cost is paid.
 *   - stride = 1, all 750 lags.  lag_stride = 2 is forbidden: the sparse code
 *     has no pulse shaping, so the odd lags carry the main lobe.
 *   - Every fold start is re-aligned with the exact period
 *         st = win_start + llround((g * K_coh + k) * 750.276923...)
 *     so the 17.7-sample fractional drift over 64 symbols never defocuses the
 *     750-point cyclic peak.
 *   - Leading-edge backtrack on a confirmed hit so the committed seed lands on
 *     the first SYNC of the preamble rather than somewhere inside it.
 *   - Stage-2 identity gate: PAR alone is not code-selective (a strong foreign
 *     Ipatov code cross-correlates at up to ~20/64 of the template norm), so a
 *     confirmed window must also win a relative code test: its last K_coh fold
 *     is correlated against all four HRP templates and the target template's
 *     single-group PAR must beat the best foreign template by
 *     code_par_ratio_min (default 2.0).  That test is a ratio of two PARs taken
 *     from the same fold, so it is insensitive to SNR and CFO; the optional
 *     absolute code_rho_min gate is off by default.
 *
 * Hot path (scan) performs no allocation: every working buffer is sized once in
 * configure().  The kernel is selected at compile time (UWB_SIMD); the AVX2
 * path is used only when the compiler flags provide it, and force_scalar (QA /
 * diagnostics) switches to the logically equivalent scalar path at runtime.
 */
#pragma once
#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>
// Same compile-time ISA convention as live_acquisition.h: the build selects the
// backend through UWB_SIMD (never -march=native); without AVX2+FMA the same
// kernels fall back to the portable scalar path, so -DUWB_SIMD=SCALAR stays
// buildable.
#if defined(__AVX2__) && defined(__FMA__)
#define UWB_PREAMBLE_XCORR_HAVE_AVX2 1
#else
#define UWB_PREAMBLE_XCORR_HAVE_AVX2 0
#endif

namespace openuwb {
namespace live {

// ---------------------------------------------------------------------------
// Configuration.  Defaults are the calibrated ±100 kHz CFO operating point of
// the design note (§2/§4); CLI validation lives in apps/live_rx_config.
// ---------------------------------------------------------------------------
struct PreambleXcorrConfig {
    size_t code_index = 10;          // 9..12 (IEEE 802.15.4a/z Ipatov 127)
    size_t preamble_reps = 256;      // preamble length; caps the backtrack
    size_t k_coh = 4;                // coherent fold symbols (4 -> -2.5 dB @100 kHz)
    size_t g_fast = 6;               // Stage-1 early-stop groups
    size_t g_full = 16;              // Stage-2 confirmation groups
    size_t hop_symbols = 128;        // scan step in symbols (0 = auto)
    double par_fast_thr = 4.5;       // Stage-1 PAR gate
    double par_thr = 6.5;            // Stage-2 PAR gate
    // Relative code-identity gate (primary): the last coherent fold is scored
    // against all four HRP templates (9..12) and must beat the best foreign
    // template by this factor.  Both PARs come from the same fold, so the test
    // is insensitive to SNR and CFO: at high SNR the on-code peak (64 coherent
    // taps) beats a foreign Ipatov template (~20 taps) by 10-20x, while a
    // pure-noise fold floors the foreign PAR near 8, i.e. ~2-3x at -14 dB.
    // 0 disables the relative gate.
    double code_par_ratio_min = 2.0;
    // Optional absolute gate: normalized correlation of the last fold with the
    // target template.  Off by default; the relative gate above is the primary
    // check.  0 disables.
    double code_rho_min = 0.0;
    int64_t holdoff_len = 589824;    // blind zone after a seed (~0.80 ms)
    bool force_scalar = false;       // QA/diagnostics: force the scalar kernel
};

struct PreambleXcorrStats {
    uint64_t windows_scanned = 0;
    uint64_t stage2_promoted = 0;
    uint64_t seeds_emitted = 0;
    uint64_t code_rejected = 0;      // Stage-2 PAR pass rejected by a code gate
};

// QA/diagnostics snapshot of one full G_full window evaluation.
struct PreambleXcorrDiagnostics {
    double peak_power = 0.0;
    double mean_power = 0.0;
    double par = 0.0;
    // Code-identity read-out from the last K_coh fold: the target template's
    // single-group PAR, the strongest foreign template's PAR, and their ratio.
    double target_group_par = 0.0;
    double other_group_par_max = 0.0;
    double code_par_ratio = 0.0;
    // Optional absolute rho (0 = no energy).
    double code_rho = 0.0;
    int64_t peak_tau = -1;
    size_t groups_used = 0;
};

// ---------------------------------------------------------------------------
// Streaming folded-preamble cross-correlation acquirer.
//
// The caller owns the power-of-two ring of native-rate SC16 samples and keeps
// [probe(), probe() + required_window_samples()) valid before each scan(), as
// with openuwb::live::Acquisition.  `valid_begin` is the oldest valid absolute
// sample coordinate in the ring; it bounds the leading-edge backtrack.  A gap
// invalidates all history: call reset(start).
//
// scan() returns the absolute seed (>= 0) of a confirmed preamble start, or -1.
// Inside its own holdoff the frontend neither scans nor commits.
// ---------------------------------------------------------------------------
class PreambleXcorrAcquisition {
public:
    static constexpr size_t kSymLen = 750;
    static constexpr double kExactPeriod = 1016.0 * 48.0 / 65.0;
    static constexpr size_t kMaxTaps = 64;   // non-zero Ipatov chips

    static bool avx2_available() noexcept;   // compiled AVX2+FMA kernel?

    PreambleXcorrAcquisition();

    // Fail-closed: throws std::invalid_argument on an invalid configuration.
    void configure(const PreambleXcorrConfig& cfg);
    void reset(int64_t start_sample);

    // Fold span of G_full*K_coh symbols plus one symbol of tail margin.
    int64_t required_window_samples() const noexcept;
    int64_t probe() const noexcept { return probe_; }
    const PreambleXcorrStats& stats() const noexcept { return stats_; }

    int64_t scan(const std::vector<std::complex<int16_t>>& ring,
                 int64_t valid_begin, double* out_metric = nullptr);

    // QA/diagnostics entry: full G_full evaluation at win_start without
    // touching probe_/holdoff_/stats_.  Returns false when the window is not
    // addressable (win_start before valid_begin, negative, or a ring too small
    // to hold one evaluation window).
    bool evaluate_window(const std::vector<std::complex<int16_t>>& ring,
                         int64_t valid_begin, int64_t win_start,
                         bool force_scalar, PreambleXcorrDiagnostics* out);

private:
    static constexpr size_t kNumCodes = 4;   // HRP code 9..12

    // Leading-edge probe: K_coh x G_probe fold at cand_pos, PAR of the
    // peak_tau triple against 32 spread non-peak lags.
    double probe_par(const std::complex<int16_t>* ring, size_t ring_size,
                     size_t mask, int64_t cand_pos, int64_t peak_tau,
                     bool use_avx);

    // Single-group PAR of fold_re_/fold_im_ for each HRP template (index =
    // code - 9).  Uses the preallocated scratch_power_ buffer.
    void code_group_pars(double* out_par, bool use_avx);

    PreambleXcorrConfig cfg_{};
    PreambleXcorrStats stats_{};

    // One ternary tap table per HRP code 9..12 (index = code - 9).
    std::array<std::array<uint16_t, kMaxTaps>, kNumCodes> pos_taps_{};
    std::array<std::array<uint16_t, kMaxTaps>, kNumCodes> neg_taps_{};
    std::array<size_t, kNumCodes> n_pos_{};
    std::array<size_t, kNumCodes> n_neg_{};
    size_t target_code_ = 0;   // cfg_.code_index - 9

    // Linear double-width fold buffers [0, 1500): the second half mirrors the
    // first so the 750-lag correlation needs no modulo.
    std::vector<float> fold_re_;
    std::vector<float> fold_im_;
    std::vector<float> acc_power_;    // 750 non-coherent accumulators
    std::vector<float> probe_re_;     // backtrack fold (K_coh x G_probe)
    std::vector<float> probe_im_;
    std::vector<float> scratch_power_;  // 750, Stage-2 per-code PAR scoring

    int64_t probe_ = 0;
    int64_t holdoff_ = 0;
    int64_t hop_samples_ = 0;
    int64_t required_samples_ = 0;
    bool use_avx_ = false;
};

} // namespace live
} // namespace openuwb
