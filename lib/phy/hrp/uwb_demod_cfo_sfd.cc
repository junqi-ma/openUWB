// SPDX-License-Identifier: GPL-3.0-or-later
//
// Out-of-line definitions for the HRP CFO/SFD stages.  The declarations stay in
// the public headers (demod_cfo.h / demod_sfd.h / rx_cfo_opt.h /
// rx_phase_sfo.h); this translation unit hosts the CFO fit, the derotation
// kernels' batch entry points, the SFD refinement stage and the phase/SFO
// helpers.
#include <openuwb/phy/hrp/rx/demod_cfo.h>
#include <openuwb/phy/hrp/rx/demod_sfd.h>
#include <openuwb/phy/hrp/rx/rx_cfo_opt.h>
#include <openuwb/phy/hrp/rx/rx_phase_sfo.h>

#include <array>
#include <cmath>
#include <limits>
#include <numeric>

namespace openuwb::demod::core {

bool fit_cfo(const std::complex<float>* rx,
             size_t n,
             const Qm35825Profile& profile,
             const TimingResult& timing,
             CfoResult& out,
             DemodScratch& scratch)
{
    out = CfoResult{};
    // The fold path publishes a dense 4-SYMBOL-group CFO grid in scratch
    // (up to 64 phases, alias-free to +/-122.8 kHz, 24-symbol skip -> 40
    // points).  It restores the information the old 64-anchor train gave the
    // legacy Stage-2 cross-check while TimingResult keeps its sparse <= 16
    // detected_peaks contract.  Non-fold / synthetic sparse trains fall back
    // to timing.peak_samples.
    const bool use_fold_cfo =
        timing.sync_peaks_folded &&
        scratch.timing_fold_cfo_samples.size() >= 4 &&
        scratch.timing_fold_cfo_corr.size() ==
            scratch.timing_fold_cfo_samples.size();
    const std::vector<int64_t>& ps =
        use_fold_cfo ? scratch.timing_fold_cfo_samples : timing.peak_samples;
    const std::vector<std::complex<float>>& pc =
        use_fold_cfo ? scratch.timing_fold_cfo_corr : timing.peak_corr;
    const size_t np = ps.size();
    if (np < 4)
        return false;

    // Match MATLAB compensateCarrierOffset: discard up to the first 24 SYNC
    // correlation phases, but retain at least 32 for the least-squares fit.
    // For shorter preambles, keep all available peaks rather than failing.
    const size_t min_fit = std::max<size_t>(2, profile.cfo_min_fit_repetitions);
    const size_t available_to_skip = np > min_fit ? np - min_fit : 0;
    const size_t skip =
        std::min(profile.cfo_skip_initial_repetitions, available_to_skip);
    const size_t nfit = np - skip;
    if (nfit < 2)
        return false;
    out.skipped_peaks = skip;
    out.fit_first_peak_sample = ps[skip];
    out.fit_last_peak_sample = ps[np - 1];

    // Gather per-peak carrier phases from the SYNC matched-filter values (the
    // noise-averaged correlation phase, NOT the raw rx[peak] sample).  Falls
    // back to the raw sample if the timing stage did not populate peak_corr.
    const auto cfo_fit_begin = std::chrono::steady_clock::now();
    const bool have_corr = pc.size() >= np;
    scratch.metric.resize(nfit);
    for (size_t k = 0; k < nfit; ++k) {
        const size_t idx = skip + k;
        scratch.metric[k] =
            have_corr ? std::arg(pc[idx])
                      : std::arg(rx[ps[idx]]);
    }
    // Unwrap: adjacent SYNC phases must stay within ±π.  When LockedBounded
    // bridges a multi-symbol interference gap, estimate the per-sample phase
    // slope from adjacent 1-symbol pairs so multi-symbol gaps unwrap cleanly.
    double unwrap_slope_rad_per_sample = 0.0;
    {
        auto collect_adj = [&](double max_dt, double& mean_dt) -> std::complex<double> {
            std::complex<double> s(0.0, 0.0);
            double dt_sum = 0.0;
            size_t cnt = 0;
            for (size_t k = 1; k < nfit; ++k) {
                const size_t idx0 = skip + k - 1;
                const size_t idx1 = skip + k;
                const double dt = static_cast<double>(
                    ps[idx1] - ps[idx0]);
                if (dt > 0.0 && dt <= max_dt) {
                    if (have_corr) {
                        const std::complex<double> c0(pc[idx0]);
                        const std::complex<double> c1(pc[idx1]);
                        s += c1 * std::conj(c0);
                    } else {
                        const double dphi = static_cast<double>(scratch.metric[k] - scratch.metric[k - 1]);
                        s += std::complex<double>(std::cos(dphi), std::sin(dphi));
                    }
                    dt_sum += dt;
                    ++cnt;
                }
            }
            mean_dt = (cnt > 0) ? (dt_sum / static_cast<double>(cnt)) : 0.0;
            return s;
        };
        double mean_dt = 0.0;
        std::complex<double> adj_sum = collect_adj(
            1.5 * static_cast<double>(kQm35SamplesPerSymbol), mean_dt);
        if (std::norm(adj_sum) > 1e-24) {
            unwrap_slope_rad_per_sample =
                std::arg(adj_sum) / static_cast<double>(kQm35SamplesPerSymbol);
        } else {
            // A dense 4-SYMBOL-group train (the fold path's shadow CFO grid,
            // or the legacy 64-anchor train) has adjacent 4P pairs: the 20P
            // average is alias-free to +/-122.8 kHz and needs no search.
            double min_gap = std::numeric_limits<double>::infinity();
            for (size_t k = 1; k < nfit; ++k)
                min_gap = std::min(
                    min_gap, static_cast<double>(ps[skip + k] -
                                                 ps[skip + k - 1]));
            if (min_gap <= 4.5 * static_cast<double>(kQm35SamplesPerSymbol)) {
                adj_sum = collect_adj(
                    20.0 * static_cast<double>(kQm35SamplesPerSymbol),
                    mean_dt);
                if (std::norm(adj_sum) > 1e-24 && mean_dt > 0.0)
                    unwrap_slope_rad_per_sample =
                        std::arg(adj_sum) / mean_dt;
            } else {
            // Sparse train (the fold path reports up to 16 super-anchors over
            // the full preamble).  There is no 1-symbol pair to unwrap from.
            // Every anchor sits on a 4-SYMBOL group boundary, so the phase
            // train has a hard alias at +/-fs/(4*P) ~= +/-122.8 kHz; the mixed
            // 16/20-symbol gaps shape the peak but do not remove that alias.
            // A timing-side CFO hint is not reliable under the real capture's
            // phase noise, so run a small periodogram over the unambiguous
            // +/-120 kHz range (the Stage-2 spec is +/-40 kHz): the true slope
            // aligns all 16 anchor phasors, an alias does not.  The final OLS
            // below then refines the chosen branch.
            constexpr double kCoarseStepHz = 1000.0;
            constexpr double kFineStepHz = 50.0;
            std::array<double, 512> raw_phase{};
            std::array<std::complex<double>, 512> anchor{};
            std::array<int, 512> group_pos{};
            if (nfit <= raw_phase.size()) {
                const double group_samples =
                    4.0 * static_cast<double>(kQm35SamplesPerSymbol);
                for (size_t k = 0; k < nfit; ++k) {
                    const double ph =
                        static_cast<double>(scratch.metric[skip + k]);
                    raw_phase[k] = ph;
                    anchor[k] = std::polar(1.0, ph);
                    const double dt = static_cast<double>(
                        ps[skip + k] -
                        ps[skip]);
                    group_pos[k] = static_cast<int>(
                        std::llround(dt / group_samples));
                }
                // The unambiguous range is set by the GCD of the anchor gaps:
                // a uniform 16-symbol train (4-group gaps) repeats every
                // fs/(4*4P) ~= 61 kHz, while the fold path's mixed 4/5-group
                // (16/20-symbol) gaps push the first alias to fs/(4P) ~=
                // 245.7 kHz.  Search only +/- half that spacing.
                int gap_gcd = 0;
                for (size_t k = 1; k < nfit; ++k)
                    gap_gcd = std::gcd(gap_gcd,
                                       group_pos[k] - group_pos[k - 1]);
                if (gap_gcd <= 0)
                    gap_gcd = 1;
                const double max_cfo_hz = std::min(
                    120000.0,
                    profile.sample_rate /
                        (2.0 * static_cast<double>(gap_gcd) * group_samples));
                // |sum_k e^{j(phi_k - u*g_k)}| maximised over the CFO grid.
                // g_k counts 4-SYMBOL groups, so one complex multiply per
                // group step evaluates the whole train; no per-(k,f)
                // transcendental is needed.
                auto coherence = [&](double f_hz) -> double {
                    const double u = 2.0 * M_PI * f_hz * group_samples /
                                     profile.sample_rate;
                    const std::complex<double> w = std::polar(1.0, -u);
                    std::complex<double> rot(1.0, 0.0);
                    std::complex<double> acc(0.0, 0.0);
                    int cur = 0;
                    for (size_t k = 0; k < nfit; ++k) {
                        while (cur < group_pos[k]) {
                            rot *= w;
                            ++cur;
                        }
                        acc += anchor[k] * rot;
                    }
                    return std::abs(acc);
                };
                double best_f = 0.0;
                double best_coh = -1.0;
                for (double f = -max_cfo_hz; f <= max_cfo_hz;
                     f += kCoarseStepHz) {
                    const double c = coherence(f);
                    if (c > best_coh) {
                        best_coh = c;
                        best_f = f;
                    }
                }
                const double f_lo = std::max(-max_cfo_hz,
                                             best_f - kCoarseStepHz);
                const double f_hi = std::min(max_cfo_hz,
                                             best_f + kCoarseStepHz);
                for (double f = f_lo; f <= f_hi; f += kFineStepHz) {
                    const double c = coherence(f);
                    if (c > best_coh) {
                        best_coh = c;
                        best_f = f;
                    }
                }
                unwrap_slope_rad_per_sample =
                    2.0 * M_PI * best_f / profile.sample_rate;
            }
            }
        }
    }
    for (size_t k = 1; k < nfit; ++k) {
        const size_t idx0 = skip + k - 1;
        const size_t idx1 = skip + k;
        const double dt = static_cast<double>(
            ps[idx1] - ps[idx0]);
        const float pred = scratch.metric[k - 1] +
                           static_cast<float>(unwrap_slope_rad_per_sample * dt);
        while (scratch.metric[k] - pred > (float)M_PI)
            scratch.metric[k] -= 2.0f * (float)M_PI;
        while (scratch.metric[k] - pred < -(float)M_PI)
            scratch.metric[k] += 2.0f * (float)M_PI;
    }

    // Linear fit of phase vs RELATIVE time (the slope is offset-invariant, so
    // subtracting t0 only improves conditioning).
    const double t0 = (double)ps[skip] / profile.sample_rate;
    detail::LinearFit fit;
    for (size_t k = 0; k < nfit; ++k) {
        const size_t idx = skip + k;
        const double t =
            (double)ps[idx] / profile.sample_rate - t0;
        const double ph = (double)scratch.metric[k];
        fit.add(t, ph);
    }
    const double denom = fit.determinant();
    if (std::abs(denom) < 1e-30)
        return false;
    const double slope = fit.slope(denom);
    out.cfo_hz = slope / (2.0 * M_PI);
    out.peaks_used = nfit;
    out.ok = true;
    out.fit_us = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - cfo_fit_begin).count());
    return true;
}

void derotate_cfo_range(const std::complex<float>* rx,
                        size_t n,
                        const Qm35825Profile& profile,
                        double cfo_hz,
                        size_t lo,
                        size_t hi,
                        DemodScratch& scratch,
                        uint64_t& rotate_us,
                        uint32_t& rotation_kernel)
{
    const auto cfo_rotate_begin = std::chrono::steady_clock::now();
    if (scratch.derotated.size() < n)
        scratch.derotated.resize(n);
    hi = std::min(hi, n);
    const double w = 2.0 * M_PI * cfo_hz / profile.sample_rate;
    if (lo < hi) {
        if (profile.rx_fast_cfo_rotation) {
            rxopt::derotate_range_block(scratch.derotated.data(), rx, n, lo, hi, w);
            rotation_kernel = 1;
        } else {
            rxopt::derotate_range_recurrence(scratch.derotated.data(), rx, n, lo, hi, w);
            rotation_kernel = 0;
        }
    } else {
        rotation_kernel = profile.rx_fast_cfo_rotation ? 1 : 0;
    }
    rotate_us = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - cfo_rotate_begin).count());
}

bool stage_cfo(const std::complex<float>* rx,
               size_t n,
               const Qm35825Profile& profile,
               const TimingResult& timing,
               CfoResult& out,
               DemodScratch& scratch)
{
    if (!fit_cfo(rx, n, profile, timing, out, scratch))
        return false;
    uint32_t rk = 0;
    derotate_cfo_range(rx, n, profile, out.cfo_hz, 0, n, scratch, out.rotate_us, rk);
    out.rotation_kernel = rk;
    return true;
}

} // namespace openuwb::demod::core

namespace openuwb::demod::rxopt {

namespace detail {

void derotate_block_anchored(std::complex<float>* dst,
                             const std::complex<float>* src,
                             const std::complex<float>* table,
                             std::size_t len,
                             float ar,
                             float ai)
{
    std::size_t k = 0;
#if UWB_RX_CFO_HAVE_AVX2
    // Complex multiply, interleaved [re, im, re, im, ...] lanes:
    //   permute(x, 0xB1) swaps re/im, addsub(a, b) = (a_even - b_even,
    //   a_odd + b_odd), which yields re = xr*yr - xi*yi, im = xi*yr + xr*yi.
    const __m256 are = _mm256_set1_ps(ar);
    const __m256 aim = _mm256_set1_ps(ai);
    for (; k + 4 <= len; k += 4) {
        const __m256 t = _mm256_loadu_ps(reinterpret_cast<const float*>(table + k));
        const __m256 ts = _mm256_permute_ps(t, 0xB1);
        const __m256 a =
            _mm256_addsub_ps(_mm256_mul_ps(t, are), _mm256_mul_ps(ts, aim));
        const __m256 s = _mm256_loadu_ps(reinterpret_cast<const float*>(src + k));
        const __m256 ss = _mm256_permute_ps(s, 0xB1);
        const __m256 r =
            _mm256_addsub_ps(_mm256_mul_ps(s, _mm256_moveldup_ps(a)),
                             _mm256_mul_ps(ss, _mm256_movehdup_ps(a)));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst + k), r);
    }
#endif
    for (; k < len; ++k) {
        const std::complex<float> t = table[k];
        const float tr = t.real();
        const float ti = t.imag();
        const float yr = tr * ar - ti * ai;
        const float yi = ti * ar + tr * ai;
        const float sr = src[k].real();
        const float si = src[k].imag();
        dst[k] = std::complex<float>(sr * yr - si * yi, si * yr + sr * yi);
    }
}

} // namespace detail

void derotate_range_block(std::complex<float>* dst,
                          const std::complex<float>* src,
                          std::size_t n,
                          std::size_t lo,
                          std::size_t hi,
                          double w)
{
    if (n == 0 || src == nullptr || dst == nullptr)
        return;
    hi = std::min(hi, n);
    if (lo >= hi)
        return;

    std::array<std::complex<float>, kCfoDerotateBlock> table;
    detail::build_cfo_table(table, w);

    const std::size_t b_start = (lo / kCfoDerotateBlock) * kCfoDerotateBlock;
    for (std::size_t b0 = b_start; b0 < hi; b0 += kCfoDerotateBlock) {
        const std::size_t k_start = (b0 < lo) ? (lo - b0) : 0;
        const std::size_t k_end =
            std::min(hi - b0, std::min(kCfoDerotateBlock, n - b0));
        if (k_start >= k_end)
            continue;
        const std::size_t len = k_end - k_start;

        // Exact double block anchor; b0 == 0 gives (1, 0) and the fold below is
        // then bit-exactly the identity, matching the legacy "no second pass
        // for the first block" behaviour.
        const double phase = -w * static_cast<double>(b0);
        const float ar = static_cast<float>(std::cos(phase));
        const float ai = static_cast<float>(std::sin(phase));

        detail::derotate_block_anchored(dst + b0 + k_start,
                                        src + b0 + k_start,
                                        table.data() + k_start,
                                        len,
                                        ar,
                                        ai);
    }
}

void derotate_range_recurrence(std::complex<float>* dst,
                               const std::complex<float>* src,
                               std::size_t n,
                               std::size_t lo,
                               std::size_t hi,
                               double w)
{
    if (n == 0 || src == nullptr || dst == nullptr)
        return;
    hi = std::min(hi, n);
    if (lo >= hi)
        return;

    constexpr std::size_t kB = 1024;
    const double cw = w;
    const double s4_phase = -4.0 * cw;
    const float s4r = static_cast<float>(std::cos(s4_phase));
    const float s4i = static_cast<float>(std::sin(s4_phase));

    for (std::size_t b0 = lo; b0 < hi; b0 += kB) {
        const std::size_t b1 = std::min(hi, b0 + kB);
        const double ph0 = -cw * static_cast<double>(b0);
        float cr[4], ci[4];
        for (int c = 0; c < 4; ++c) {
            const double ph = ph0 - cw * static_cast<double>(c);
            cr[c] = static_cast<float>(std::cos(ph));
            ci[c] = static_cast<float>(std::sin(ph));
        }
        std::size_t i = b0;
        for (; i + 4 <= b1; i += 4) {
            for (int c = 0; c < 4; ++c) {
                const float sr = src[i + c].real();
                const float si = src[i + c].imag();
                dst[i + c] = std::complex<float>(sr * cr[c] - si * ci[c],
                                                 sr * ci[c] + si * cr[c]);
                const float nr = cr[c] * s4r - ci[c] * s4i;
                const float ni = cr[c] * s4i + ci[c] * s4r;
                cr[c] = nr;
                ci[c] = ni;
            }
        }
        for (; i < b1; ++i) {
            const double ph = -cw * static_cast<double>(i);
            const float c = static_cast<float>(std::cos(ph));
            const float s = static_cast<float>(std::sin(ph));
            const float sr = src[i].real();
            const float si = src[i].imag();
            dst[i] = std::complex<float>(sr * c - si * s, sr * s + si * c);
        }
    }
}

void derotate_range_scalar(std::complex<float>* dst,
                           const std::complex<float>* src,
                           std::size_t n,
                           std::size_t lo,
                           std::size_t hi,
                           double w)
{
    if (n == 0 || src == nullptr || dst == nullptr)
        return;
    hi = std::min(hi, n);
    if (lo >= hi)
        return;
    for (std::size_t i = lo; i < hi; ++i) {
        const double phase = -w * static_cast<double>(i);
        const float c = static_cast<float>(std::cos(phase));
        const float s = static_cast<float>(std::sin(phase));
        const float re = src[i].real();
        const float im = src[i].imag();
        dst[i] = std::complex<float>(re * c - im * s, re * s + im * c);
    }
}

void derotate_block(std::vector<std::complex<float>>& dst,
                    const std::complex<float>* src,
                    std::size_t n,
                    double w)
{
    dst.resize(n);
    derotate_range_block(dst.data(), src, n, 0, n, w);
}

void derotate_reference(std::complex<float>* dst,
                        const std::complex<float>* src,
                        std::size_t n,
                        double w)
{
    if (n == 0 || src == nullptr || dst == nullptr)
        return;
    for (std::size_t i = 0; i < n; ++i) {
        const double phase = -w * static_cast<double>(i);
        const double c = std::cos(phase);
        const double s = std::sin(phase);
        const double re = static_cast<double>(src[i].real());
        const double im = static_cast<double>(src[i].imag());
        dst[i] = std::complex<float>(static_cast<float>(re * c - im * s),
                                     static_cast<float>(re * s + im * c));
    }
}

void derotate_range_reference(std::complex<float>* dst,
                              const std::complex<float>* src,
                              std::size_t n,
                              std::size_t lo,
                              std::size_t hi,
                              double w)
{
    if (n == 0 || src == nullptr || dst == nullptr)
        return;
    hi = std::min(hi, n);
    if (lo >= hi)
        return;
    for (std::size_t i = lo; i < hi; ++i) {
        const double phase = -w * static_cast<double>(i);
        const double c = std::cos(phase);
        const double s = std::sin(phase);
        const double re = static_cast<double>(src[i].real());
        const double im = static_cast<double>(src[i].imag());
        dst[i] = std::complex<float>(static_cast<float>(re * c - im * s),
                                     static_cast<float>(re * s + im * c));
    }
}

PhrPhaseRefineResult refine_phase_cfo_from_phr(
    const std::complex<float>* phr_z_win,
    size_t nsym,
    const int8_t* phr_g1,
    int8_t sfd_polarity,
    std::complex<float> rot_init,
    double phr_symbol_samples,
    double phr_last_to_payload_samples,
    float min_coherence,
    double max_delta_w_rad)
{
    PhrPhaseRefineResult out;

    // --- Argument / value guards -------------------------------------------
    const bool rot_finite =
        std::isfinite(rot_init.real()) && std::isfinite(rot_init.imag());
    const double rot_mag = rot_finite ? static_cast<double>(std::abs(rot_init))
                                      : 0.0;
    const bool rot_usable = rot_finite && std::isfinite(rot_mag) &&
                            rot_mag > 1e-12;
    if (rot_usable)
        out.rot_at_payload_start = rot_init;
    else
        out.rot_at_payload_start = std::complex<float>(1.0f, 0.0f);
    // Every early return must leave the mean-phase rotation at the same value
    // as `rot_at_payload_start` (i.e. rot_init, or unit when rot_init is
    // unusable); Pass 2 overrides it with the centre-symbol mean phase.
    out.rot_mean_phase = out.rot_at_payload_start;

    if (phr_z_win == nullptr || phr_g1 == nullptr)
        return out;
    if (nsym < 4)
        return out;
    if (!std::isfinite(phr_symbol_samples) || phr_symbol_samples <= 0.0)
        return out;
    if (!std::isfinite(phr_last_to_payload_samples))
        return out;
    if (!std::isfinite(min_coherence))
        return out;
    if (!std::isfinite(max_delta_w_rad) || max_delta_w_rad <= 0.0)
        return out;
    if (!rot_usable)
        return out;

    const double rr = static_cast<double>(rot_init.real()) / rot_mag;
    const double ri = static_cast<double>(rot_init.imag()) / rot_mag;

    // tilde[k] = (1-2*g1)*pol*z, expressed in the rot_init frame:
    //   u[k] = rot_init_normalised * tilde[k].
    // Recomputing it per pass avoids any allocation in this per-packet helper.
    auto tilde_at = [&](size_t k, std::complex<double>& u) -> bool {
        const double zr = static_cast<double>(phr_z_win[k].real());
        const double zi = static_cast<double>(phr_z_win[k].imag());
        if (!std::isfinite(zr) || !std::isfinite(zi))
            return false;
        const int sign_bit = 1 - 2 * static_cast<int>(phr_g1[k]);
        const double sgn = static_cast<double>(sign_bit) *
                           static_cast<double>(sfd_polarity);
        const double tr = sgn * zr;
        const double ti = sgn * zi;
        if (!std::isfinite(tr) || !std::isfinite(ti))
            return false;
        u = std::complex<double>(rr * tr - ri * ti, ri * tr + rr * ti);
        return std::isfinite(u.real()) && std::isfinite(u.imag());
    };

    // --- Pass 1: coherence and lag-1 residual CFO --------------------------
    constexpr double kEps = 1e-30;
    std::complex<double> sum(0.0, 0.0);
    std::complex<double> lag(0.0, 0.0);
    double abs_sum = 0.0;
    std::complex<double> prev(0.0, 0.0);
    for (size_t k = 0; k < nsym; ++k) {
        std::complex<double> u;
        if (!tilde_at(k, u))
            return out;
        // |u| == |tilde| because the frame rotation has unit modulus.
        sum += u;
        abs_sum += std::abs(u);
        if (k > 0)
            lag += u * std::conj(prev);
        prev = u;
    }
    const double coherence = std::abs(sum) / (abs_sum + kEps);
    out.coherence = static_cast<float>(coherence);
    if (coherence < static_cast<double>(min_coherence))
        return out; // applied = false, rot = rot_init, delta = 0

    const double dtheta = std::atan2(lag.imag(), lag.real());
    double delta = dtheta / phr_symbol_samples;
    delta = std::max(-max_delta_w_rad, std::min(max_delta_w_rad, delta));
    if (!std::isfinite(delta))
        return out;

    // --- Pass 2: mean phase at the centre symbol, extrapolated to payload ---
    const double kc = 0.5 * static_cast<double>(nsym - 1);
    std::complex<double> acc(0.0, 0.0);
    for (size_t k = 0; k < nsym; ++k) {
        std::complex<double> u;
        if (!tilde_at(k, u))
            return out;
        const double ang = -delta *
                           (static_cast<double>(k) - kc) *
                           phr_symbol_samples;
        acc += u * std::complex<double>(std::cos(ang), std::sin(ang));
    }
    const double theta_c = std::atan2(acc.imag(), acc.real());
    const double d =
        phr_last_to_payload_samples +
        (static_cast<double>(nsym - 1) - kc) * phr_symbol_samples;
    const double theta_payload = theta_c + delta * d;
    if (!std::isfinite(theta_c) || !std::isfinite(theta_payload))
        return out;

    // rot = rot_init * exp(-j*theta_payload) in the unit-modulus frame.
    const double cr = std::cos(-theta_payload);
    const double ci = std::sin(-theta_payload);
    out.rot_at_payload_start = std::complex<float>(
        static_cast<float>(rr * cr - ri * ci),
        static_cast<float>(ri * cr + rr * ci));
    // Mean-phase-only companion: the same unit-modulus construction with
    // theta_payload replaced by theta_c, i.e. no delta*d extrapolation.
    const double mr = std::cos(-theta_c);
    const double mi = std::sin(-theta_c);
    out.rot_mean_phase = std::complex<float>(
        static_cast<float>(rr * mr - ri * mi),
        static_cast<float>(ri * mr + rr * mi));
    out.delta_cfo_w_rad = delta;
    out.applied = true;
    return out;
}

bool build_fractional_cir_bank(const std::complex<float>* fir_taps,
                               size_t tap_count,
                               size_t phase_bins,
                               std::complex<float>* out_bank)
{
    if (fir_taps == nullptr || out_bank == nullptr)
        return false;
    if (tap_count == 0 || phase_bins == 0)
        return false;
    // Sanity bound: keeps tap_count*phase_bins far away from size_t overflow.
    if (phase_bins > (static_cast<size_t>(1) << 16))
        return false;

    // bin 0: bit-exact integer grid.
    std::copy(fir_taps, fir_taps + tap_count, out_bank);
    if (phase_bins == 1)
        return true;

    double ref_norm_sq = 0.0;
    for (size_t m = 0; m < tap_count; ++m) {
        const double re = static_cast<double>(fir_taps[m].real());
        const double im = static_cast<double>(fir_taps[m].imag());
        ref_norm_sq += re * re + im * im;
    }
    const double ref_norm = std::sqrt(ref_norm_sq);

    auto tap_value = [&](long long idx) -> std::complex<double> {
        if (idx < 0 || idx >= static_cast<long long>(tap_count))
            return std::complex<double>(0.0, 0.0);
        const std::complex<float>& v = fir_taps[static_cast<size_t>(idx)];
        return std::complex<double>(static_cast<double>(v.real()),
                                    static_cast<double>(v.imag()));
    };

    for (size_t b = 1; b < phase_bins; ++b) {
        const double delta =
            static_cast<double>(b) / static_cast<double>(phase_bins);
        std::complex<float>* dst = out_bank + b * tap_count;
        double norm_sq = 0.0;
        for (size_t m = 0; m < tap_count; ++m) {
            // Catmull-Rom interpolation of h at x = m - delta: the template
            // peak therefore moves by +delta samples.
            const double x = static_cast<double>(m) - delta;
            const double fl = std::floor(x);
            const long long i = static_cast<long long>(fl);
            const double t = x - fl;
            const double t2 = t * t;
            const double t3 = t2 * t;
            const double w0 = -0.5 * t3 + t2 - 0.5 * t;
            const double w1 = 1.5 * t3 - 2.5 * t2 + 1.0;
            const double w2 = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
            const double w3 = 0.5 * t3 - 0.5 * t2;
            const std::complex<double> v =
                w0 * tap_value(i - 1) + w1 * tap_value(i) +
                w2 * tap_value(i + 1) + w3 * tap_value(i + 2);
            dst[m] = std::complex<float>(static_cast<float>(v.real()),
                                         static_cast<float>(v.imag()));
            const double re = static_cast<double>(dst[m].real());
            const double im = static_cast<double>(dst[m].imag());
            norm_sq += re * re + im * im;
        }
        const double norm = std::sqrt(norm_sq);
        if (ref_norm > 0.0 && norm > 0.0) {
            const float scale = static_cast<float>(ref_norm / norm);
            for (size_t m = 0; m < tap_count; ++m)
                dst[m] *= scale;
        }
    }
    return true;
}

} // namespace openuwb::demod::rxopt
