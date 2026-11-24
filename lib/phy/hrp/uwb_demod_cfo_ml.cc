// SPDX-License-Identifier: GPL-3.0-or-later
//
// Complex-domain ML CFO estimator + large-SFO orchestration.
//
// This translation unit implements the two entry points declared in
// include/openuwb/phy/hrp/rx/demod_cfo.h:
//   * estimate_cfo_from_symbols_ml() - zero-unwrap periodogram frequency
//     estimator over an evenly spaced complex symbol sequence.
//   * estimate_cfo_ml()              - two-step SFO bootstrap / drift-gated
//     sampling / coherence-guard orchestration around the estimator.
//
// Design notes (see docs/OpenCode执行方案_CFO估计架构重构_2026-11-24.md §2.2/§2.3):
//   * No phase unwrapping: the ML periodogram is immune to ±π jumps.
//   * Zero heap allocation on the normal path: every intermediate buffer lives
//     in DemodScratch (cfo_sym_z / cfo_coarse_z / cfo_fft_buf / cfo_mag_buf).
//     Only counts that exceed the reserved capacities grow them once.
//   * No I/O, logging or timing here; the pipeline owns instrumentation.

#include <openuwb/phy/hrp/rx/demod_cfo.h>
#include <openuwb/phy/hrp/rx/demod_scratch.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace openuwb::demod::core {
namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;

// Smallest power of two >= v (v == 0 -> 1).  Local to avoid pulling in the dsp
// header's helper (and its exception-throwing fft).
size_t next_pow2_size(size_t v)
{
    size_t p = 1;
    while (p < v)
        p <<= 1;
    return p;
}

// In-place, allocation-free radix-2 forward FFT (e^{-j2πkn/N} convention).
// Mirrors openuwb::core::fft_radix2 but reports failure instead of throwing and
// never allocates.  Returns false when `a` is null or `n` is not a power of two.
bool fft_radix2_inplace(std::complex<double>* a, size_t n)
{
    if (a == nullptr || n == 0 || (n & (n - 1)) != 0)
        return false;

    // Decimation-in-time bit-reversal permutation.
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }

    // Iterative butterflies.  Angles use double-precision cos/sin.
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = -kTwoPi / static_cast<double>(len);
        const std::complex<double> wlen(std::cos(ang), std::sin(ang));
        const size_t half = len >> 1;
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (size_t j = 0; j < half; ++j) {
                const std::complex<double> u = a[i + j];
                const std::complex<double> v = a[i + j + half] * w;
                a[i + j] = u + v;
                a[i + j + half] = u - v;
                w *= wlen;
            }
        }
    }
    return true;
}

// Map an FFT bin index onto the angular frequency in [-π, π) (rad/symbol).
double bin_angular_freq(size_t k, size_t nfft)
{
    const double f = kTwoPi * static_cast<double>(k) / static_cast<double>(nfft);
    return (k >= nfft / 2) ? (f - kTwoPi) : f;
}

// Wrap an angle to [-π, π).
double wrap_angle(double w)
{
    while (w >= kPi)
        w -= kTwoPi;
    while (w < -kPi)
        w += kTwoPi;
    return w;
}

} // namespace

// ---------------------------------------------------------------------------
// estimate_cfo_from_symbols_ml
// ---------------------------------------------------------------------------
// `z[k]` are complex SYNC correlations sampled at a uniform symbol period.
// Returns the maximum-likelihood single-tone frequency via:
//   1. Huber/median soft clipping (burst-interference robustness, no angle ops)
//   2. zero-padded radix-2 FFT global coarse search (no ±2π ambiguity)
//   3. Newton-Raphson refinement of the periodogram objective with backtracking
// ---------------------------------------------------------------------------
bool estimate_cfo_from_symbols_ml(const std::complex<float>* z,
                                  size_t count,
                                  double symbol_period_samples,
                                  double sample_rate_hz,
                                  double max_abs_cfo_hz,
                                  double huber_median_scale,
                                  DemodScratch& scratch,
                                  CfoMlResult& out)
{
    out = CfoMlResult{};

    // ---- defensive validation -------------------------------------------
    if (z == nullptr || count < 8)
        return false;
    if (!(symbol_period_samples > 0.0) || !std::isfinite(symbol_period_samples))
        return false;
    if (!(sample_rate_hz > 0.0) || !std::isfinite(sample_rate_hz))
        return false;
    for (size_t k = 0; k < count; ++k) {
        if (!std::isfinite(z[k].real()) || !std::isfinite(z[k].imag()))
            return false;
    }

    // ---- Step 1: Huber / median soft clipping ---------------------------
    // Magnitudes go into the scratch buffer purely to compute the median; the
    // per-element clip scale is re-derived from z[k] while filling the FFT
    // buffer so the nth_element permutation cannot corrupt index association.
    auto& mags = scratch.cfo_mag_buf;
    mags.resize(count); // reserved to 512; grows only for >512 symbols
    for (size_t k = 0; k < count; ++k)
        mags[k] = std::abs(z[k]);

    const size_t mid = count / 2;
    std::nth_element(mags.begin(), mags.begin() + mid, mags.end());
    double med = static_cast<double>(mags[mid]);
    if (!(med > 1e-12)) {
        // Degenerate median (mostly-zero sequence): fall back to the mean.
        double sum = 0.0;
        for (size_t k = 0; k < count; ++k)
            sum += static_cast<double>(mags[k]);
        med = sum / static_cast<double>(count);
        if (!(med > 0.0))
            return false; // all-zero / all-non-positive magnitude input
    }

    const bool do_clip = (huber_median_scale > 0.0) && std::isfinite(huber_median_scale);
    const double clip = do_clip ? huber_median_scale * med : 0.0;

    // Soft-clipped complex symbol, scaled in place when later re-filling.
    auto make_zt = [&](size_t k) -> std::complex<double> {
        double scale = 1.0;
        if (do_clip) {
            const double m = static_cast<double>(std::abs(z[k]));
            if (m > 0.0) {
                const double lim = clip / m;
                if (lim < 1.0)
                    scale = lim;
            }
        }
        return std::complex<double>(static_cast<double>(z[k].real()) * scale,
                                    static_cast<double>(z[k].imag()) * scale);
    };

    // ---- Step 2: zero-padded radix-2 FFT coarse search ------------------
    const size_t nfft = (count > 512) ? next_pow2_size(2 * count) : size_t{ 1024 };
    auto& buf = scratch.cfo_fft_buf;
    buf.assign(nfft, std::complex<double>(0.0, 0.0)); // reserved 2048 normally

    double sum_mag = 0.0;
    for (size_t k = 0; k < count; ++k) {
        const std::complex<double> zt = make_zt(k);
        buf[k] = zt;
        sum_mag += std::abs(zt);
    }
    if (!(sum_mag > 0.0))
        return false;

    if (!fft_radix2_inplace(buf.data(), nfft))
        return false;

    // Legitimate search band in rad/symbol.
    double max_w = kPi;
    if (max_abs_cfo_hz > 0.0) {
        const double mw = kTwoPi * max_abs_cfo_hz * symbol_period_samples /
                          sample_rate_hz;
        if (std::isfinite(mw) && mw < kPi)
            max_w = std::max(0.0, mw);
    }

    size_t k_best = 0;
    double best_pow = -1.0;
    for (size_t k = 0; k < nfft; ++k) {
        const double w = bin_angular_freq(k, nfft);
        if (std::abs(w) > max_w)
            continue;
        const double p = std::norm(buf[k]);
        if (p > best_pow) {
            best_pow = p;
            k_best = k;
        }
    }
    // Bin 0 (w = 0) always satisfies |w| <= max_w, so k_best is valid here.

    // Three-point parabolic interpolation on log-magnitudes (cyclic neighbours).
    const size_t mask = nfft - 1; // nfft is a power of two
    const size_t km1 = (k_best + nfft - 1) & mask;
    const size_t kp1 = (k_best + 1) & mask;
    const double a = std::log(std::abs(buf[km1]) + 1e-300);
    const double b = std::log(std::abs(buf[k_best]) + 1e-300);
    const double c = std::log(std::abs(buf[kp1]) + 1e-300);
    const double denom = a - 2.0 * b + c;
    double delta = 0.0;
    if (std::abs(denom) > 1e-300)
        delta = 0.5 * (a - c) / denom;
    if (!(delta > -0.5)) // also rejects NaN
        delta = -0.5;
    if (delta > 0.5)
        delta = 0.5;

    double w = wrap_angle(bin_angular_freq(k_best, nfft) +
                          delta * kTwoPi / static_cast<double>(nfft));

    // The FFT destroyed the time-domain samples; overwrite the front of the
    // buffer with the soft-clipped sequence for the Newton refinement.
    for (size_t k = 0; k < count; ++k)
        buf[k] = make_zt(k);

    // ---- Step 3: Newton-Raphson periodogram refinement ------------------
    // Zero-centred index for conditioning; phase origin is the window centre.
    const double u0 = -0.5 * static_cast<double>(count - 1);
    auto eval_spectrum = [&](double omega, std::complex<double>& s0,
                             std::complex<double>& s1,
                             std::complex<double>& s2) {
        // Rotate by e^{-j ω u_k}; step is e^{-j ω} per unit u.
        const std::complex<double> step(std::cos(-omega), std::sin(-omega));
        double ph_arg = -omega * u0;
        std::complex<double> ph(std::cos(ph_arg), std::sin(ph_arg));
        s0 = std::complex<double>(0.0, 0.0);
        s1 = std::complex<double>(0.0, 0.0);
        s2 = std::complex<double>(0.0, 0.0);
        for (size_t k = 0; k < count; ++k) {
            if ((k & 31u) == 0) {
                // Re-anchor every 32 points to avoid rotation drift.
                ph_arg = -omega * (u0 + static_cast<double>(k));
                ph = std::complex<double>(std::cos(ph_arg), std::sin(ph_arg));
            }
            const double uk = u0 + static_cast<double>(k);
            const std::complex<double> v = buf[k] * ph;
            s0 += v;
            s1 += uk * v;
            s2 += (uk * uk) * v;
            ph *= step;
        }
    };

    std::complex<double> s0, s1, s2;
    const double max_step = kPi / static_cast<double>(count);
    for (int iter = 0; iter < 5; ++iter) {
        eval_spectrum(w, s0, s1, s2);
        const double J = std::norm(s0);
        const double J1 = 2.0 * std::imag(std::conj(s0) * s1);
        const double J2 = 2.0 * (std::norm(s1) - std::real(std::conj(s0) * s2));
        if (!(J2 < -1e-18))
            break; // not concave / flat: keep the FFT estimate

        double step = -J1 / J2;
        if (!std::isfinite(step))
            break;
        if (step > max_step)
            step = max_step;
        if (step < -max_step)
            step = -max_step;

        // Backtracking line search: accept only a step that does not lower J.
        double accepted = 0.0;
        for (int bt = 0; bt < 3; ++bt) {
            if (!(std::abs(step) > 1e-14))
                break;
            std::complex<double> t0, t1, t2;
            eval_spectrum(w + step, t0, t1, t2);
            if (std::norm(t0) >= J) {
                accepted = step;
                break;
            }
            step *= 0.5;
        }
        if (accepted == 0.0)
            break;
        w += accepted;
        if (std::abs(accepted) < 1e-8)
            break; // converged
    }

    // Final metrics at the refined frequency.
    eval_spectrum(w, s0, s1, s2);
    const double s0_abs = std::abs(s0);
    const double cfo_hz = w * sample_rate_hz / (kTwoPi * symbol_period_samples);

    if (!std::isfinite(cfo_hz) || !std::isfinite(s0_abs))
        return false;

    out.ok = true;
    out.cfo_hz = cfo_hz;
    out.effective_period = symbol_period_samples;
    out.coherence = s0_abs / (sum_mag + 1e-12); // [0, 1] by construction
    out.center_phase_rad = std::arg(s0);
    out.symbols_used = count;
    out.used_fractional_grid = false;
    out.used_coherence_fallback = false;
    return true;
}

// ---------------------------------------------------------------------------
// estimate_cfo_ml
// ---------------------------------------------------------------------------
// Two-step large-SFO orchestration (§2.3 flowchart):
//   * pick the effective symbol period P_eff from an optional Stage-2 CFO hint
//     (coupled SFO), or bootstrap it from a centred coarse window,
//   * sample the whole preamble on a drift-gated integer/fractional grid,
//   * in decoupled SFO mode, fall back to the coarse result when the full
//     window's coherence collapses (misleading period_hint).
// ---------------------------------------------------------------------------
bool estimate_cfo_ml(const std::complex<float>* rx,
                     size_t n,
                     const Qm35825Profile& profile,
                     const std::vector<std::complex<float>>& template_wf,
                     const CfoMlConfig& cfg,
                     DemodScratch& scratch,
                     CfoMlResult& out)
{
    out = CfoMlResult{};

    if (rx == nullptr || cfg.symbol_end <= cfg.symbol_begin)
        return false;
    const size_t total = cfg.symbol_end - cfg.symbol_begin;
    if (total < 8)
        return false;
    if (!(cfg.nominal_period > 0.0) || !std::isfinite(cfg.nominal_period))
        return false;
    if (!(cfg.sample_rate_hz > 0.0) || !std::isfinite(cfg.sample_rate_hz))
        return false;

    const double nominal = cfg.nominal_period;
    double p_eff = nominal;
    CfoMlResult coarse_res; // ok=false unless Step 1 actually ran and succeeded

    const bool hint_path = cfg.has_cfo_hint && cfg.sfo_from_cfo &&
                           cfg.carrier_freq_hz > 0.0;

    if (hint_path) {
        // Coupled SFO: the Stage-2 CFO already pins the true sample period.
        const double denom = 1.0 + cfg.cfo_hint_hz / cfg.carrier_freq_hz;
        if (std::isfinite(denom) && denom > 1e-9)
            p_eff = nominal / denom;
    } else {
        // ---- Step 1: centred coarse window ------------------------------
        size_t W = std::min(cfg.coarse_window_symbols, total);
        if (W < 8)
            W = total; // total >= 8 guaranteed above

        // Decoupled extreme SFO: shrink the window so it still stays inside the
        // pulse main lobe over its whole span.
        if (!cfg.sfo_from_cfo && cfg.period_hint > 0.0) {
            const double diff = std::abs(cfg.period_hint - nominal);
            if (diff > 1e-6) {
                const double raw = std::floor(1.5 / diff);
                size_t wd = (raw >= static_cast<double>(W))
                                ? W
                                : static_cast<size_t>(raw);
                if (wd < 16)
                    wd = std::min<size_t>(16, W);
                W = wd;
            }
        }

        const size_t c_begin = cfg.symbol_begin + (total - W) / 2;
        const size_t c_end = c_begin + W;

        scratch.cfo_coarse_z.resize(W); // reserved 128; coarse W <= 64 normally
        size_t first_coarse = 0;
        size_t used_c = sample_sync_symbols_simd(
            rx, n, profile, template_wf, cfg.origin_sample, c_begin, c_end,
            nominal, /*use_fractional_interp=*/false, scratch,
            scratch.cfo_coarse_z.data(), &first_coarse);
        if (used_c > W)
            used_c = W; // defensive clamp (sampler must never over-report)

        if (used_c >= 8) {
            estimate_cfo_from_symbols_ml(
                scratch.cfo_coarse_z.data(), used_c, nominal,
                cfg.sample_rate_hz, cfg.max_abs_cfo_hz, cfg.huber_median_scale,
                scratch, coarse_res);
        }

        if (cfg.sfo_from_cfo && cfg.carrier_freq_hz > 0.0 && coarse_res.ok) {
            const double denom = 1.0 + coarse_res.cfo_hz / cfg.carrier_freq_hz;
            if (std::isfinite(denom) && denom > 1e-9)
                p_eff = nominal / denom;
        } else if (cfg.period_hint > 0.0 &&
                   std::abs(cfg.period_hint - nominal) <= 500e-6 * nominal) {
            p_eff = cfg.period_hint;
        } else {
            p_eff = nominal;
        }
    }

    if (!std::isfinite(p_eff) || !(p_eff > 0.0))
        p_eff = nominal;

    // ---- drift gate: integer fast path vs fractional (4-lag Catmull-Rom) --
    const double delta_max =
        0.5 * static_cast<double>(total - 1) * std::abs(p_eff - nominal);
    const bool use_frac = delta_max >= cfg.frac_drift_threshold_samples;

    // ---- Step 2: full-window sampling + ML estimate ---------------------
    scratch.cfo_sym_z.resize(total); // reserved 512
    size_t first_full = 0;
    size_t used = sample_sync_symbols_simd(
        rx, n, profile, template_wf, cfg.origin_sample, cfg.symbol_begin,
        cfg.symbol_end, p_eff, use_frac, scratch, scratch.cfo_sym_z.data(),
        &first_full);
    if (used > total)
        used = total; // never read past the reserved window

    CfoMlResult full_res;
    if (used >= 8) {
        estimate_cfo_from_symbols_ml(
            scratch.cfo_sym_z.data(), used, p_eff, cfg.sample_rate_hz,
            cfg.max_abs_cfo_hz, cfg.huber_median_scale, scratch, full_res);
    }

    // ---- Coherence guard (decoupled mode only) --------------------------
    // If the full-length estimate is unavailable or its coherence collapses
    // relative to the centred coarse window, the period_hint was misleading:
    // return the nominal-period coarse result instead.
    if (!cfg.sfo_from_cfo && coarse_res.ok &&
        (!full_res.ok || full_res.coherence < 0.70 * coarse_res.coherence)) {
        out = coarse_res;
        out.used_coherence_fallback = true;
        // The fallback was sampled on the nominal integer grid, so report the
        // period that was actually used rather than p_eff.
        out.effective_period = nominal;
        out.used_fractional_grid = false;
        return out.ok;
    }

    if (!full_res.ok)
        return false;

    out = full_res;
    out.used_fractional_grid = use_frac;
    out.effective_period = p_eff;
    return out.ok;
}

} // namespace openuwb::demod::core
