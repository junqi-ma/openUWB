// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <openuwb/phy/hrp/rx/demod_sync.h>
#include <openuwb/phy/hrp/rx/rx_cfo_opt.h>

namespace openuwb::demod::core {
// ---------------------------------------------------------------------------
// Stage 2 — CFO estimation + compensation (R1).
// Uses the stable SYNC peaks after the startup transient to fit phase vs time
// linearly; the slope / 2pi is the CFO in Hz.  MATLAB skips up to its first
// 24 peaks while retaining at least 32, controlled here by the separate CFO
// profile fields (not the CIR skip policy). The frame is then derotated by
// exp(-j2pi*f*n/fs) and phase-resolved against the preamble waveform. Golden
// reference: CFO = 0 Hz for the clean baseband signal.
// ---------------------------------------------------------------------------
// Decoupled CFO estimation (Task A1).
// Fits the CFO frequency from SYNC peaks without derotating samples.

bool fit_cfo(const std::complex<float>* rx,
             size_t n,
             const Qm35825Profile& profile,
             const TimingResult& timing,
             CfoResult& out,
             DemodScratch& scratch);

// Decoupled CFO range rotation (Task A1).
// Compensates [lo, hi) where phase origin is strictly sample 0.
void derotate_cfo_range(const std::complex<float>* rx,
                        size_t n,
                        const Qm35825Profile& profile,
                        double cfo_hz,
                        size_t lo,
                        size_t hi,
                        DemodScratch& scratch,
                        uint64_t& rotate_us,
                        uint32_t& rotation_kernel);

bool stage_cfo(const std::complex<float>* rx,
               size_t n,
               const Qm35825Profile& profile,
               const TimingResult& timing,
               CfoResult& out,
               DemodScratch& scratch);

// ---------------------------------------------------------------------------
// ML CFO estimator (complex-domain periodogram, no phase unwrapping).
//
// Replaces the angle-unwrap + unweighted OLS path for the pipeline stages while
// `fit_cfo()` / `stage_cfo()` above keep their public signature and behaviour
// for the single-stage QA entry points.
//
// Sampling is decoupled from the timing argmax: the SYNC symbols are sampled on
// a uniform (optionally fractional) grid anchored at `origin_sample`, so
// multipath tap hopping cannot inject phase steps into the estimator.
// ---------------------------------------------------------------------------
struct CfoMlConfig {
    size_t symbol_begin = 0;      // first SYNC symbol index (inclusive)
    size_t symbol_end = 0;        // one past the last SYNC symbol index
    int64_t origin_sample = 0;    // sample of SYNC symbol #0 relative to rx[0]
    double nominal_period = static_cast<double>(kQm35SamplesPerSymbol);
    double period_hint = 0.0;     // optional Stage 1 period (decoupled SFO)
    double cfo_hint_hz = 0.0;     // optional Stage 2 CFO (skips the coarse step)
    bool has_cfo_hint = false;
    bool sfo_from_cfo = true;     // mirrors profile.rx_sfo_from_cfo
    double carrier_freq_hz = 6489.6e6;
    double sample_rate_hz = 998.4e6;
    double max_abs_cfo_hz = 0.0;  // <= 0 searches the full Nyquist band
    double frac_drift_threshold_samples = 0.35;
    double huber_median_scale = 2.5; // <= 0 disables soft clipping
    size_t coarse_window_symbols = 64;
};

struct CfoMlResult {
    bool ok = false;
    double cfo_hz = 0.0;
    double effective_period = 0.0;
    double coherence = 0.0;       // [0, 1]
    double center_phase_rad = 0.0;
    size_t symbols_used = 0;
    bool used_fractional_grid = false;
    bool used_coherence_fallback = false;
};

// Collects [symbol_begin, symbol_end) complex SYNC correlations into `out_z`
// (capacity >= symbol_end - symbol_begin) on a uniform grid of `period_samples`
// anchored at `origin_sample`.  Uses the existing SIMD matched-filter kernels
// only.  Returns the number of contiguously sampled symbols; the index of the
// first sampled symbol is written to `*out_first_symbol` when non-null.
size_t sample_sync_symbols_simd(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    int64_t origin_sample,
    size_t symbol_begin,
    size_t symbol_end,
    double period_samples,
    bool use_fractional_interp,
    DemodScratch& scratch,
    std::complex<float>* out_z,
    size_t* out_first_symbol = nullptr);

// ML frequency estimate from an evenly spaced complex symbol sequence
// (zero heap allocation, no phase unwrapping).
bool estimate_cfo_from_symbols_ml(
    const std::complex<float>* z,
    size_t count,
    double symbol_period_samples,
    double sample_rate_hz,
    double max_abs_cfo_hz,
    double huber_median_scale,
    DemodScratch& scratch,
    CfoMlResult& out);

// End-to-end two-step adaptive ML CFO estimate: optional coarse centre window
// for the SFO bootstrap, drift-gated integer/fractional grid sampling, and the
// coherence guard for a misleading `period_hint`.
bool estimate_cfo_ml(
    const std::complex<float>* rx,
    size_t n,
    const Qm35825Profile& profile,
    const std::vector<std::complex<float>>& template_wf,
    const CfoMlConfig& cfg,
    DemodScratch& scratch,
    CfoMlResult& out);

} // namespace openuwb::demod::core
