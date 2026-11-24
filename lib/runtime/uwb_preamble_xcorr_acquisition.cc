/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Out-of-line implementation of the folded-preamble cross-correlation
 * acquisition frontend (declarations in
 * include/openuwb/runtime/preamble_xcorr_acquisition.h).
 *
 * The translation unit contains the four heavy pieces the public header keeps
 * out: the exact-period coherent fold with ring-wrap handling, the zero-multiply
 * 64-tap sparse cyclic correlation (AVX2+FMA 8-wide plus a logically identical
 * scalar path), the Stage-1/Stage-2 PAR scan with two-level early stop, and the
 * leading-edge backtrack that snaps a confirmed window onto the first SYNC
 * symbol.  All scratch storage lives in the object and is sized in configure().
 */

#include <openuwb/runtime/preamble_xcorr_acquisition.h>

#include <openuwb/phy/hrp/common/profile.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>

#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>
#endif

namespace openuwb {
namespace live {
namespace {

constexpr size_t kN = PreambleXcorrAcquisition::kSymLen;      // 750
constexpr size_t kN2 = 2 * kN;                                // 1500
constexpr size_t kProbeGroups = 4;                            // G_probe
constexpr double kProbeParThr = 3.5;                          // backtrack gate
constexpr size_t kBackStepSymbols = 16;
constexpr size_t kBgLags = 32;
constexpr size_t kBgGuard = 4;                                // skip peak_tau +- 4

// llround(symbols * P_nat), the same rounding used for every fold start so the
// group grid never drifts from the sample grid.
inline int64_t sym_span(double symbols)
{
    return static_cast<int64_t>(std::llround(symbols * PreambleXcorrAcquisition::kExactPeriod));
}

// ---------------------------------------------------------------------------
// Fold: add one 750-sample symbol at native index `src` into the float
// accumulators.  int16 -> float keeps the metric scale-free (PAR is a ratio);
// the four coherent symbols sum to at most 4*32767 < 2^24, so the float
// accumulation is exact and the AVX2 and scalar folds agree bit for bit.
// ---------------------------------------------------------------------------
void fold_run_scalar(const std::complex<int16_t>* src, size_t n,
                     float* fre, float* fim)
{
    for (size_t i = 0; i < n; ++i) {
        fre[i] += static_cast<float>(src[i].real());
        fim[i] += static_cast<float>(src[i].imag());
    }
}

#if UWB_PREAMBLE_XCORR_HAVE_AVX2
// One 256-bit load holds 8 complex samples; vcvtepi16_epi32 splits the two
// 128-bit halves (I/Q interleaved), then two shuffle/permute pairs deinterleave
// into re/im 8-vectors.
void fold_run_avx2(const std::complex<int16_t>* src, size_t n,
                   float* fre, float* fim)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        const __m256i iq =
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + i));
        const __m256 lo = _mm256_cvtepi32_ps(
            _mm256_cvtepi16_epi32(_mm256_castsi256_si128(iq)));
        const __m256 hi = _mm256_cvtepi32_ps(
            _mm256_cvtepi16_epi32(_mm256_extracti128_si256(iq, 1)));
        const __m256 re = _mm256_castpd_ps(_mm256_permute4x64_pd(
            _mm256_castps_pd(_mm256_shuffle_ps(lo, hi, _MM_SHUFFLE(2, 0, 2, 0))),
            0xD8));
        const __m256 im = _mm256_castpd_ps(_mm256_permute4x64_pd(
            _mm256_castps_pd(_mm256_shuffle_ps(lo, hi, _MM_SHUFFLE(3, 1, 3, 1))),
            0xD8));
        _mm256_storeu_ps(fre + i, _mm256_add_ps(_mm256_loadu_ps(fre + i), re));
        _mm256_storeu_ps(fim + i, _mm256_add_ps(_mm256_loadu_ps(fim + i), im));
    }
    fold_run_scalar(src + i, n - i, fre + i, fim + i);
}
#endif

inline void fold_run(const std::complex<int16_t>* src, size_t n,
                     float* fre, float* fim, bool use_avx)
{
#if UWB_PREAMBLE_XCORR_HAVE_AVX2
    if (use_avx) {
        fold_run_avx2(src, n, fre, fim);
        return;
    }
#else
    (void)use_avx;
#endif
    fold_run_scalar(src, n, fre, fim);
}

// Fold K_coh symbols of group `g` starting at win_start into fre/fim[0,750),
// then mirror to the upper half (memcpy, no modulo anywhere downstream).  The
// lower half is cleared first: one group is one coherent fold, and the caller
// accumulates the squared magnitude non-coherently across groups.
void fold_group(const std::complex<int16_t>* ring, size_t ring_size, size_t mask,
                int64_t win_start, size_t g, size_t k_coh, bool use_avx,
                float* fre, float* fim)
{
    std::memset(fre, 0, kN * sizeof(float));
    std::memset(fim, 0, kN * sizeof(float));
    for (size_t k = 0; k < k_coh; ++k) {
        const int64_t st = win_start + sym_span(static_cast<double>(g * k_coh + k));
        const size_t p0 = static_cast<size_t>(st) & mask;
        if (ring_size >= kN) {
            size_t n1 = ring_size - p0;
            if (n1 > kN) n1 = kN;
            fold_run(ring + p0, n1, fre, fim, use_avx);
            if (n1 < kN) fold_run(ring, kN - n1, fre + n1, fim + n1, use_avx);
        } else {
            for (size_t i = 0; i < kN; ++i) {
                const auto v = ring[(p0 + i) & mask];
                fre[i] += static_cast<float>(v.real());
                fim[i] += static_cast<float>(v.imag());
            }
        }
    }
    std::memcpy(fre + kN, fre, kN * sizeof(float));
    std::memcpy(fim + kN, fim, kN * sizeof(float));
}

// ---------------------------------------------------------------------------
// Sparse cyclic correlation over all 750 lags, add/sub only:
//   acc[t] += | sum_{+off} f[t+off] - sum_{-off} f[t+off] |^2
// The doubled fold buffer makes t+off < 1500 always addressable.
//
// The 64 taps are split over four independent accumulators per component.  The
// sum is the same 64 terms in the same order for the AVX2 and the scalar path
// (offset i feeds accumulator i&3, combined as (a0+a1)+(a2+a3)), so both
// kernels agree bit for bit up to the final power add; four chains instead of
// one keep the loop load-bound rather than add-latency-bound.
// ---------------------------------------------------------------------------
inline void sparse_sum(const float* fre, const float* fim, size_t t,
                       const uint16_t* pos, size_t npos,
                       const uint16_t* neg, size_t nneg,
                       float* out_re, float* out_im)
{
    float r0 = 0.0f, r1 = 0.0f, r2 = 0.0f, r3 = 0.0f;
    float i0 = 0.0f, i1 = 0.0f, i2 = 0.0f, i3 = 0.0f;
    size_t i = 0;
    for (; i + 4 <= npos; i += 4) {
        const size_t o0 = pos[i];
        const size_t o1 = pos[i + 1];
        const size_t o2 = pos[i + 2];
        const size_t o3 = pos[i + 3];
        r0 += fre[t + o0]; i0 += fim[t + o0];
        r1 += fre[t + o1]; i1 += fim[t + o1];
        r2 += fre[t + o2]; i2 += fim[t + o2];
        r3 += fre[t + o3]; i3 += fim[t + o3];
    }
    for (; i < npos; ++i) {
        r0 += fre[t + pos[i]];
        i0 += fim[t + pos[i]];
    }
    i = 0;
    for (; i + 4 <= nneg; i += 4) {
        const size_t o0 = neg[i];
        const size_t o1 = neg[i + 1];
        const size_t o2 = neg[i + 2];
        const size_t o3 = neg[i + 3];
        r0 -= fre[t + o0]; i0 -= fim[t + o0];
        r1 -= fre[t + o1]; i1 -= fim[t + o1];
        r2 -= fre[t + o2]; i2 -= fim[t + o2];
        r3 -= fre[t + o3]; i3 -= fim[t + o3];
    }
    for (; i < nneg; ++i) {
        r0 -= fre[t + neg[i]];
        i0 -= fim[t + neg[i]];
    }
    // Same combination order as the AVX2 path.
    out_re[0] = r0;
    out_re[1] = r1;
    out_re[2] = r2;
    out_re[3] = r3;
    out_im[0] = i0;
    out_im[1] = i1;
    out_im[2] = i2;
    out_im[3] = i3;
}

inline float combine4(const float* a)
{
    return (a[0] + a[1]) + (a[2] + a[3]);
}

void correlate_scalar(const float* fre, const float* fim,
                      const uint16_t* pos, size_t npos,
                      const uint16_t* neg, size_t nneg, float* acc)
{
    float sr[4];
    float si[4];
    for (size_t t = 0; t < kN; ++t) {
        sparse_sum(fre, fim, t, pos, npos, neg, nneg, sr, si);
        const float re = combine4(sr);
        const float im = combine4(si);
        acc[t] += re * re + im * im;
    }
}

#if UWB_PREAMBLE_XCORR_HAVE_AVX2
// One 8-lag block.  Four independent partial sums per component keep the loop
// load-bound instead of add-latency-bound; the accumulation order is the same
// as the scalar path ((a0+a1)+(a2+a3)), so the two kernels stay bit-identical.
inline void correlate_avx2_block8(const float* fre, const float* fim,
                                  const uint16_t* pos, size_t npos,
                                  const uint16_t* neg, size_t nneg,
                                  float* acc, size_t t)
{
    __m256 r0 = _mm256_setzero_ps();
    __m256 r1 = _mm256_setzero_ps();
    __m256 r2 = _mm256_setzero_ps();
    __m256 r3 = _mm256_setzero_ps();
    __m256 i0 = _mm256_setzero_ps();
    __m256 i1 = _mm256_setzero_ps();
    __m256 i2 = _mm256_setzero_ps();
    __m256 i3 = _mm256_setzero_ps();
    size_t i = 0;
    for (; i + 4 <= npos; i += 4) {
        const size_t o0 = pos[i];
        const size_t o1 = pos[i + 1];
        const size_t o2 = pos[i + 2];
        const size_t o3 = pos[i + 3];
        r0 = _mm256_add_ps(r0, _mm256_loadu_ps(fre + t + o0));
        i0 = _mm256_add_ps(i0, _mm256_loadu_ps(fim + t + o0));
        r1 = _mm256_add_ps(r1, _mm256_loadu_ps(fre + t + o1));
        i1 = _mm256_add_ps(i1, _mm256_loadu_ps(fim + t + o1));
        r2 = _mm256_add_ps(r2, _mm256_loadu_ps(fre + t + o2));
        i2 = _mm256_add_ps(i2, _mm256_loadu_ps(fim + t + o2));
        r3 = _mm256_add_ps(r3, _mm256_loadu_ps(fre + t + o3));
        i3 = _mm256_add_ps(i3, _mm256_loadu_ps(fim + t + o3));
    }
    for (; i < npos; ++i) {
        r0 = _mm256_add_ps(r0, _mm256_loadu_ps(fre + t + pos[i]));
        i0 = _mm256_add_ps(i0, _mm256_loadu_ps(fim + t + pos[i]));
    }
    i = 0;
    for (; i + 4 <= nneg; i += 4) {
        const size_t o0 = neg[i];
        const size_t o1 = neg[i + 1];
        const size_t o2 = neg[i + 2];
        const size_t o3 = neg[i + 3];
        r0 = _mm256_sub_ps(r0, _mm256_loadu_ps(fre + t + o0));
        i0 = _mm256_sub_ps(i0, _mm256_loadu_ps(fim + t + o0));
        r1 = _mm256_sub_ps(r1, _mm256_loadu_ps(fre + t + o1));
        i1 = _mm256_sub_ps(i1, _mm256_loadu_ps(fim + t + o1));
        r2 = _mm256_sub_ps(r2, _mm256_loadu_ps(fre + t + o2));
        i2 = _mm256_sub_ps(i2, _mm256_loadu_ps(fim + t + o2));
        r3 = _mm256_sub_ps(r3, _mm256_loadu_ps(fre + t + o3));
        i3 = _mm256_sub_ps(i3, _mm256_loadu_ps(fim + t + o3));
    }
    for (; i < nneg; ++i) {
        r0 = _mm256_sub_ps(r0, _mm256_loadu_ps(fre + t + neg[i]));
        i0 = _mm256_sub_ps(i0, _mm256_loadu_ps(fim + t + neg[i]));
    }
    const __m256 re = _mm256_add_ps(_mm256_add_ps(r0, r1), _mm256_add_ps(r2, r3));
    const __m256 im = _mm256_add_ps(_mm256_add_ps(i0, i1), _mm256_add_ps(i2, i3));
    __m256 p = _mm256_loadu_ps(acc + t);
    p = _mm256_fmadd_ps(im, im, p);
    p = _mm256_fmadd_ps(re, re, p);
    _mm256_storeu_ps(acc + t, p);
}

void correlate_avx2(const float* fre, const float* fim,
                    const uint16_t* pos, size_t npos,
                    const uint16_t* neg, size_t nneg, float* acc)
{
    constexpr size_t kLane = 8;
    constexpr size_t kLane2 = 2 * kLane;
    // 2 x 8 lags per iteration: the two blocks share the tap-offset loads and
    // halve the outer-loop/address bookkeeping.  Each lane keeps the identical
    // 4-partial accumulation order, so this is numerically the same kernel.
    size_t t = 0;
    for (; t + kLane2 <= kN; t += kLane2) {
        const float* fre0 = fre + t;
        const float* fim0 = fim + t;
        const float* fre1 = fre0 + kLane;
        const float* fim1 = fim0 + kLane;
        __m256 r0 = _mm256_setzero_ps();
        __m256 r1 = _mm256_setzero_ps();
        __m256 r2 = _mm256_setzero_ps();
        __m256 r3 = _mm256_setzero_ps();
        __m256 i0 = _mm256_setzero_ps();
        __m256 i1 = _mm256_setzero_ps();
        __m256 i2 = _mm256_setzero_ps();
        __m256 i3 = _mm256_setzero_ps();
        __m256 R0 = _mm256_setzero_ps();
        __m256 R1 = _mm256_setzero_ps();
        __m256 R2 = _mm256_setzero_ps();
        __m256 R3 = _mm256_setzero_ps();
        __m256 I0 = _mm256_setzero_ps();
        __m256 I1 = _mm256_setzero_ps();
        __m256 I2 = _mm256_setzero_ps();
        __m256 I3 = _mm256_setzero_ps();
        size_t i = 0;
        for (; i + 4 <= npos; i += 4) {
            const size_t o0 = pos[i];
            const size_t o1 = pos[i + 1];
            const size_t o2 = pos[i + 2];
            const size_t o3 = pos[i + 3];
            r0 = _mm256_add_ps(r0, _mm256_loadu_ps(fre0 + o0));
            i0 = _mm256_add_ps(i0, _mm256_loadu_ps(fim0 + o0));
            R0 = _mm256_add_ps(R0, _mm256_loadu_ps(fre1 + o0));
            I0 = _mm256_add_ps(I0, _mm256_loadu_ps(fim1 + o0));
            r1 = _mm256_add_ps(r1, _mm256_loadu_ps(fre0 + o1));
            i1 = _mm256_add_ps(i1, _mm256_loadu_ps(fim0 + o1));
            R1 = _mm256_add_ps(R1, _mm256_loadu_ps(fre1 + o1));
            I1 = _mm256_add_ps(I1, _mm256_loadu_ps(fim1 + o1));
            r2 = _mm256_add_ps(r2, _mm256_loadu_ps(fre0 + o2));
            i2 = _mm256_add_ps(i2, _mm256_loadu_ps(fim0 + o2));
            R2 = _mm256_add_ps(R2, _mm256_loadu_ps(fre1 + o2));
            I2 = _mm256_add_ps(I2, _mm256_loadu_ps(fim1 + o2));
            r3 = _mm256_add_ps(r3, _mm256_loadu_ps(fre0 + o3));
            i3 = _mm256_add_ps(i3, _mm256_loadu_ps(fim0 + o3));
            R3 = _mm256_add_ps(R3, _mm256_loadu_ps(fre1 + o3));
            I3 = _mm256_add_ps(I3, _mm256_loadu_ps(fim1 + o3));
        }
        for (; i < npos; ++i) {
            const size_t o = pos[i];
            r0 = _mm256_add_ps(r0, _mm256_loadu_ps(fre0 + o));
            i0 = _mm256_add_ps(i0, _mm256_loadu_ps(fim0 + o));
            R0 = _mm256_add_ps(R0, _mm256_loadu_ps(fre1 + o));
            I0 = _mm256_add_ps(I0, _mm256_loadu_ps(fim1 + o));
        }
        i = 0;
        for (; i + 4 <= nneg; i += 4) {
            const size_t o0 = neg[i];
            const size_t o1 = neg[i + 1];
            const size_t o2 = neg[i + 2];
            const size_t o3 = neg[i + 3];
            r0 = _mm256_sub_ps(r0, _mm256_loadu_ps(fre0 + o0));
            i0 = _mm256_sub_ps(i0, _mm256_loadu_ps(fim0 + o0));
            R0 = _mm256_sub_ps(R0, _mm256_loadu_ps(fre1 + o0));
            I0 = _mm256_sub_ps(I0, _mm256_loadu_ps(fim1 + o0));
            r1 = _mm256_sub_ps(r1, _mm256_loadu_ps(fre0 + o1));
            i1 = _mm256_sub_ps(i1, _mm256_loadu_ps(fim0 + o1));
            R1 = _mm256_sub_ps(R1, _mm256_loadu_ps(fre1 + o1));
            I1 = _mm256_sub_ps(I1, _mm256_loadu_ps(fim1 + o1));
            r2 = _mm256_sub_ps(r2, _mm256_loadu_ps(fre0 + o2));
            i2 = _mm256_sub_ps(i2, _mm256_loadu_ps(fim0 + o2));
            R2 = _mm256_sub_ps(R2, _mm256_loadu_ps(fre1 + o2));
            I2 = _mm256_sub_ps(I2, _mm256_loadu_ps(fim1 + o2));
            r3 = _mm256_sub_ps(r3, _mm256_loadu_ps(fre0 + o3));
            i3 = _mm256_sub_ps(i3, _mm256_loadu_ps(fim0 + o3));
            R3 = _mm256_sub_ps(R3, _mm256_loadu_ps(fre1 + o3));
            I3 = _mm256_sub_ps(I3, _mm256_loadu_ps(fim1 + o3));
        }
        for (; i < nneg; ++i) {
            const size_t o = neg[i];
            r0 = _mm256_sub_ps(r0, _mm256_loadu_ps(fre0 + o));
            i0 = _mm256_sub_ps(i0, _mm256_loadu_ps(fim0 + o));
            R0 = _mm256_sub_ps(R0, _mm256_loadu_ps(fre1 + o));
            I0 = _mm256_sub_ps(I0, _mm256_loadu_ps(fim1 + o));
        }
        __m256 re = _mm256_add_ps(_mm256_add_ps(r0, r1), _mm256_add_ps(r2, r3));
        __m256 im = _mm256_add_ps(_mm256_add_ps(i0, i1), _mm256_add_ps(i2, i3));
        __m256 p = _mm256_loadu_ps(acc + t);
        p = _mm256_fmadd_ps(im, im, p);
        p = _mm256_fmadd_ps(re, re, p);
        _mm256_storeu_ps(acc + t, p);
        re = _mm256_add_ps(_mm256_add_ps(R0, R1), _mm256_add_ps(R2, R3));
        im = _mm256_add_ps(_mm256_add_ps(I0, I1), _mm256_add_ps(I2, I3));
        p = _mm256_loadu_ps(acc + t + kLane);
        p = _mm256_fmadd_ps(im, im, p);
        p = _mm256_fmadd_ps(re, re, p);
        _mm256_storeu_ps(acc + t + kLane, p);
    }
    for (; t + kLane <= kN; t += kLane)
        correlate_avx2_block8(fre, fim, pos, npos, neg, nneg, acc, t);
    // 750 = 46*16 + 8 + 6: same arithmetic as the vector body.
    float sr[4];
    float si[4];
    for (; t < kN; ++t) {
        sparse_sum(fre, fim, t, pos, npos, neg, nneg, sr, si);
        const float re = combine4(sr);
        const float im = combine4(si);
        acc[t] += re * re + im * im;
    }
}
#endif

inline void correlate_sparse(const float* fre, const float* fim,
                             const uint16_t* pos, size_t npos,
                             const uint16_t* neg, size_t nneg, bool use_avx,
                             float* acc)
{
#if UWB_PREAMBLE_XCORR_HAVE_AVX2
    if (use_avx) {
        correlate_avx2(fre, fim, pos, npos, neg, nneg, acc);
        return;
    }
#else
    (void)use_avx;
#endif
    correlate_scalar(fre, fim, pos, npos, neg, nneg, acc);
}

// Single-lag correlation (backtrack probe: a few dozen lags per candidate).
inline void xcorr_lag(const float* fre, const float* fim,
                      const uint16_t* pos, size_t npos,
                      const uint16_t* neg, size_t nneg, size_t t,
                      float& re, float& im)
{
    re = 0.0f;
    im = 0.0f;
    for (size_t i = 0; i < npos; ++i) {
        re += fre[t + pos[i]];
        im += fim[t + pos[i]];
    }
    for (size_t i = 0; i < nneg; ++i) {
        re -= fre[t + neg[i]];
        im -= fim[t + neg[i]];
    }
}

void power_stats(const float* acc, double& peak, double& mean, int64_t& tau)
{
    double sum = 0.0;
    double mx = 0.0;
    int64_t arg = 0;
    for (size_t t = 0; t < kN; ++t) {
        const double v = static_cast<double>(acc[t]);
        sum += v;
        if (v > mx) {
            mx = v;
            arg = static_cast<int64_t>(t);
        }
    }
    peak = mx;
    mean = sum / static_cast<double>(kN);
    tau = arg;
}

inline double par_of(double peak, double mean)
{
    return peak / std::max(mean, 1e-12);
}

// Code-identity metric on the last coherent fold at the cyclic peak:
//   rho = |sum_{+off} F[t+off] - sum_{-off} F[t+off]|
//         / ( ||F at the taps|| * sqrt(n_pos + n_neg) )
// A correct code is phase-coherent with the template (rho -> 1); a foreign
// Ipatov code only reaches |cross-correlation| / n_taps (0.172 for 11/64).
// Scalar double on the doubled fold buffer, so the AVX2 and the scalar build
// report bit-identical values (t + off < 1500 by construction).
double code_identity(const float* fre, const float* fim,
                     const uint16_t* pos, size_t npos,
                     const uint16_t* neg, size_t nneg, int64_t tau,
                     double* out_num_mag = nullptr,
                     double* out_energy = nullptr)
{
    const size_t t = static_cast<size_t>(tau);
    double re = 0.0;
    double im = 0.0;
    double energy = 0.0;
    for (size_t i = 0; i < npos; ++i) {
        const double fr = static_cast<double>(fre[t + pos[i]]);
        const double fi = static_cast<double>(fim[t + pos[i]]);
        re += fr;
        im += fi;
        energy += fr * fr + fi * fi;
    }
    for (size_t i = 0; i < nneg; ++i) {
        const double fr = static_cast<double>(fre[t + neg[i]]);
        const double fi = static_cast<double>(fim[t + neg[i]]);
        re -= fr;
        im -= fi;
        energy += fr * fr + fi * fi;
    }
    const double num_mag = std::sqrt(re * re + im * im);
    const double n_taps = static_cast<double>(npos + nneg);
    if (out_num_mag) *out_num_mag = num_mag;
    if (out_energy) *out_energy = energy;
    if (energy <= 0.0 || n_taps <= 0.0) return 0.0;
    return num_mag / (std::sqrt(energy) * std::sqrt(n_taps));
}

} // namespace

// ---------------------------------------------------------------------------
// Construction / configuration
// ---------------------------------------------------------------------------
PreambleXcorrAcquisition::PreambleXcorrAcquisition()
{
    configure(PreambleXcorrConfig{});
}

bool PreambleXcorrAcquisition::avx2_available() noexcept
{
#if UWB_PREAMBLE_XCORR_HAVE_AVX2
    return true;
#else
    return false;
#endif
}

void PreambleXcorrAcquisition::configure(const PreambleXcorrConfig& cfg)
{
    if (cfg.code_index < 9 || cfg.code_index > 12) {
        throw std::invalid_argument(
            "preamble xcorr: code_index must be 9..12");
    }
    if (cfg.k_coh < 1) {
        throw std::invalid_argument("preamble xcorr: k_coh must be >= 1");
    }
    if (cfg.g_fast < 1 || cfg.g_full < cfg.g_fast) {
        throw std::invalid_argument(
            "preamble xcorr: require 1 <= g_fast <= g_full");
    }
    if (cfg.k_coh > 256 || cfg.g_full > 1024 ||
        cfg.g_full * cfg.k_coh > 4096) {
        throw std::invalid_argument(
            "preamble xcorr: K_coh/G_full fold span out of range");
    }
    if (cfg.hop_symbols > 256) {
        throw std::invalid_argument("preamble xcorr: hop_symbols must be <= 256");
    }
    if (!(cfg.par_fast_thr > 0.0) || !std::isfinite(cfg.par_fast_thr) ||
        !(cfg.par_thr > 0.0) || !std::isfinite(cfg.par_thr)) {
        throw std::invalid_argument(
            "preamble xcorr: PAR thresholds must be finite and > 0");
    }
    if (!std::isfinite(cfg.code_rho_min) || cfg.code_rho_min < 0.0 ||
        cfg.code_rho_min > 1.0) {
        throw std::invalid_argument(
            "preamble xcorr: code_rho_min must be finite and in [0, 1]");
    }
    if (!std::isfinite(cfg.code_par_ratio_min) || cfg.code_par_ratio_min < 0.0) {
        throw std::invalid_argument(
            "preamble xcorr: code_par_ratio_min must be finite and >= 0");
    }
    if (cfg.holdoff_len < 0) {
        throw std::invalid_argument("preamble xcorr: holdoff_len must be >= 0");
    }

    cfg_ = cfg;
    if (cfg_.hop_symbols == 0) {
        // Automatic step: half the preamble, bounded by the fold length and
        // the 128-symbol scan grid (§4.1).
        size_t hop = cfg_.preamble_reps / 2;
        if (hop > 128) hop = 128;
        if (hop < cfg_.k_coh) hop = cfg_.k_coh;
        cfg_.hop_symbols = hop;
    }

    // Ternary templates for all four HRP codes: 64 non-zero chips -> native
    // cyclic offsets.  Split by sign so the hot loop is pure add/sub and sorted
    // so both kernels walk the fold buffer in the same order (identical
    // rounding).  The target table is the one scan() correlates with; the other
    // three are only touched by the Stage-2 relative code gate.
    target_code_ = cfg_.code_index - 9;
    for (size_t ci = 0; ci < kNumCodes; ++ci) {
        const int8_t* code = demod::GetPreambleCode(9 + ci);
        size_t& np = n_pos_[ci];
        size_t& nn = n_neg_[ci];
        np = 0;
        nn = 0;
        for (size_t c = 0; c < demod::kQm35CodeLength; ++c) {
            const int8_t chip = code[c];
            if (chip == 0) continue;
            if (np + nn >= kMaxTaps) {
                throw std::runtime_error(
                    "preamble xcorr: Ipatov template exceeds 64 non-zero chips");
            }
            const int64_t off =
                static_cast<int64_t>(std::llround(8.0 * static_cast<double>(c) *
                                                 48.0 / 65.0)) %
                static_cast<int64_t>(kN);
            if (chip > 0) {
                pos_taps_[ci][np++] = static_cast<uint16_t>(off);
            } else {
                neg_taps_[ci][nn++] = static_cast<uint16_t>(off);
            }
        }
        if (np + nn == 0) {
            throw std::runtime_error("preamble xcorr: empty Ipatov template");
        }
        std::sort(pos_taps_[ci].begin(), pos_taps_[ci].begin() + np);
        std::sort(neg_taps_[ci].begin(), neg_taps_[ci].begin() + nn);
    }

    // One-time working storage.  Nothing below allocates again.
    fold_re_.assign(kN2, 0.0f);
    fold_im_.assign(kN2, 0.0f);
    acc_power_.assign(kN, 0.0f);
    probe_re_.assign(kN2, 0.0f);
    probe_im_.assign(kN2, 0.0f);
    scratch_power_.assign(kN, 0.0f);

    hop_samples_ = sym_span(static_cast<double>(cfg_.hop_symbols));
    required_samples_ =
        sym_span(static_cast<double>(cfg_.g_full * cfg_.k_coh - 1)) +
        static_cast<int64_t>(kN) + static_cast<int64_t>(kN);
    use_avx_ = !cfg_.force_scalar && avx2_available();

    reset(0);
}

void PreambleXcorrAcquisition::reset(int64_t start_sample)
{
    stats_ = PreambleXcorrStats{};
    probe_ = start_sample;
    holdoff_ = 0;
    std::fill(acc_power_.begin(), acc_power_.end(), 0.0f);
    std::fill(fold_re_.begin(), fold_re_.end(), 0.0f);
    std::fill(fold_im_.begin(), fold_im_.end(), 0.0f);
    std::fill(probe_re_.begin(), probe_re_.end(), 0.0f);
    std::fill(probe_im_.begin(), probe_im_.end(), 0.0f);
    std::fill(scratch_power_.begin(), scratch_power_.end(), 0.0f);
}

int64_t PreambleXcorrAcquisition::required_window_samples() const noexcept
{
    return required_samples_;
}

// ---------------------------------------------------------------------------
// Leading-edge probe: fold K_coh x G_probe symbols at cand_pos, compare the
// peak_tau triple against 32 spread non-peak lags.
// ---------------------------------------------------------------------------
double PreambleXcorrAcquisition::probe_par(const std::complex<int16_t>* ring,
                                           size_t ring_size, size_t mask,
                                           int64_t cand_pos, int64_t peak_tau,
                                           bool use_avx)
{
    std::array<size_t, 3> sig{};
    sig[0] = static_cast<size_t>((peak_tau + static_cast<int64_t>(kN) - 1) %
                                static_cast<int64_t>(kN));
    sig[1] = static_cast<size_t>(peak_tau);
    sig[2] = static_cast<size_t>((peak_tau + 1) % static_cast<int64_t>(kN));

    // Uniformly spread background lags, skipping peak_tau +- 4.
    std::array<size_t, kBgLags> bg{};
    size_t n_bg = 0;
    for (size_t i = 0; i < kBgLags; ++i) {
        const size_t lag = (i * kN) / kBgLags;
        size_t d = (lag + kN - static_cast<size_t>(peak_tau)) % kN;
        if (d > kN / 2) d = kN - d;
        if (d <= kBgGuard) continue;
        bg[n_bg++] = lag;
    }

    std::array<float, 3> psig{};
    std::array<float, kBgLags> pbg{};
    for (size_t g = 0; g < kProbeGroups; ++g) {
        fold_group(ring, ring_size, mask, cand_pos, g, cfg_.k_coh, use_avx,
                   probe_re_.data(), probe_im_.data());
        for (size_t j = 0; j < sig.size(); ++j) {
            float re = 0.0f;
            float im = 0.0f;
            xcorr_lag(probe_re_.data(), probe_im_.data(),
                      pos_taps_[target_code_].data(), n_pos_[target_code_],
                      neg_taps_[target_code_].data(), n_neg_[target_code_],
                      sig[j], re, im);
            psig[j] += re * re + im * im;
        }
        for (size_t j = 0; j < n_bg; ++j) {
            float re = 0.0f;
            float im = 0.0f;
            xcorr_lag(probe_re_.data(), probe_im_.data(),
                      pos_taps_[target_code_].data(), n_pos_[target_code_],
                      neg_taps_[target_code_].data(), n_neg_[target_code_],
                      bg[j], re, im);
            pbg[j] += re * re + im * im;
        }
    }

    const float p_sig = std::max(psig[0], std::max(psig[1], psig[2]));
    double bg_sum = 0.0;
    for (size_t j = 0; j < n_bg; ++j) bg_sum += static_cast<double>(pbg[j]);
    const double p_bg = n_bg ? bg_sum / static_cast<double>(n_bg) : 0.0;
    return static_cast<double>(p_sig) / std::max(p_bg, 1e-12);
}

// ---------------------------------------------------------------------------
// Single-group PAR of the current fold (fold_re_/fold_im_) for every HRP
// template.  Cost is one 750-lag sparse correlation per code and only runs on
// Stage-2 confirmations; scratch_power_ is preallocated in configure().
// ---------------------------------------------------------------------------
void PreambleXcorrAcquisition::code_group_pars(double* out_par, bool use_avx)
{
    for (size_t c = 0; c < kNumCodes; ++c) {
        std::fill(scratch_power_.begin(), scratch_power_.end(), 0.0f);
        correlate_sparse(fold_re_.data(), fold_im_.data(),
                         pos_taps_[c].data(), n_pos_[c], neg_taps_[c].data(),
                         n_neg_[c], use_avx, scratch_power_.data());
        double peak = 0.0;
        double mean = 0.0;
        int64_t tau = 0;
        power_stats(scratch_power_.data(), peak, mean, tau);
        out_par[c] = par_of(peak, mean);
    }
}

// ---------------------------------------------------------------------------
// Streaming scan: Stage 1 (G_fast groups) -> Stage 2 (remaining groups) ->
// code-identity gate -> leading-edge backtrack -> holdoff.
// ---------------------------------------------------------------------------
int64_t PreambleXcorrAcquisition::scan(
    const std::vector<std::complex<int16_t>>& ring, int64_t valid_begin,
    double* out_metric)
{
    if (ring.empty() || (ring.size() & (ring.size() - 1)) != 0) return -1;

    if (probe_ < valid_begin) probe_ = valid_begin;
    if (probe_ < holdoff_) {
        // Inside the blind zone: never commit, and skip the whole zone so the
        // caller's drive loop makes progress.
        probe_ = holdoff_;
        return -1;
    }

    const std::complex<int16_t>* rp = ring.data();
    const size_t ring_size = ring.size();
    const size_t mask = ring_size - 1;
    const int64_t win_start = probe_;

    std::fill(acc_power_.begin(), acc_power_.end(), 0.0f);
    ++stats_.windows_scanned;
    for (size_t g = 0; g < cfg_.g_fast; ++g) {
        fold_group(rp, ring_size, mask, win_start, g, cfg_.k_coh, use_avx_,
                   fold_re_.data(), fold_im_.data());
        correlate_sparse(fold_re_.data(), fold_im_.data(),
                         pos_taps_[target_code_].data(), n_pos_[target_code_],
                         neg_taps_[target_code_].data(), n_neg_[target_code_],
                         use_avx_, acc_power_.data());
    }

    double peak = 0.0;
    double mean = 0.0;
    int64_t tau = 0;
    power_stats(acc_power_.data(), peak, mean, tau);
    if (par_of(peak, mean) < cfg_.par_fast_thr) {
        probe_ += hop_samples_;
        return -1;
    }

    ++stats_.stage2_promoted;
    for (size_t g = cfg_.g_fast; g < cfg_.g_full; ++g) {
        fold_group(rp, ring_size, mask, win_start, g, cfg_.k_coh, use_avx_,
                   fold_re_.data(), fold_im_.data());
        correlate_sparse(fold_re_.data(), fold_im_.data(),
                         pos_taps_[target_code_].data(), n_pos_[target_code_],
                         neg_taps_[target_code_].data(), n_neg_[target_code_],
                         use_avx_, acc_power_.data());
    }
    power_stats(acc_power_.data(), peak, mean, tau);
    const double par_full = par_of(peak, mean);
    if (par_full < cfg_.par_thr) {
        probe_ += hop_samples_;
        return -1;
    }

    // Code-identity gate.  fold_re_/fold_im_ still hold the last group
    // (g_full-1) fold, i.e. exactly the K_coh-symbol template match the PAR
    // test above used.  Score that fold against all four HRP templates and
    // require the target template to win by code_par_ratio_min.  Both PARs
    // come from the same fold, so the test is a ratio of two quantities that
    // scale identically with SNR and CFO: a strong foreign Ipatov code can
    // reach PAR by coincidence, but never beats its own template.
    if (cfg_.code_par_ratio_min > 0.0) {
        double group_par[kNumCodes] = {0.0, 0.0, 0.0, 0.0};
        code_group_pars(group_par, use_avx_);
        const double target = group_par[target_code_];
        double other_max = 0.0;
        for (size_t c = 0; c < kNumCodes; ++c) {
            if (c != target_code_ && group_par[c] > other_max) {
                other_max = group_par[c];
            }
        }
        const bool pass = (other_max <= 0.0)
                              ? (target > 0.0)
                              : (target >=
                                 cfg_.code_par_ratio_min * other_max);
        if (!pass) {
            ++stats_.code_rejected;
            probe_ += hop_samples_;
            return -1;
        }
    }

    // Optional absolute gate (off by default): normalized code correlation of
    // the same fold with the target template.
    if (cfg_.code_rho_min > 0.0) {
        double rho_energy = 0.0;
        const double rho = code_identity(
            fold_re_.data(), fold_im_.data(),
            pos_taps_[target_code_].data(), n_pos_[target_code_],
            neg_taps_[target_code_].data(), n_neg_[target_code_], tau, nullptr,
            &rho_energy);
        if (rho_energy <= 0.0 || rho < cfg_.code_rho_min) {
            ++stats_.code_rejected;
            probe_ += hop_samples_;
            return -1;
        }
    }

    // Leading-edge backtrack (§4.2.5).  m = 0 is win_start itself, which Stage 2
    // has just confirmed, so it is the initial edge; the loop only looks for
    // earlier candidates still inside the preamble.
    const int64_t lower = std::max<int64_t>(valid_begin, holdoff_);
    int64_t edge_pos = win_start;
    for (int64_t back = static_cast<int64_t>(kBackStepSymbols);
         back <= static_cast<int64_t>(cfg_.preamble_reps);
         back += static_cast<int64_t>(kBackStepSymbols)) {
        const int64_t cand = win_start - sym_span(static_cast<double>(back));
        if (cand < lower) break;
        if (probe_par(rp, ring_size, mask, cand, tau, use_avx_) < kProbeParThr) {
            break;
        }
        edge_pos = cand;
    }

    const int64_t seed = std::max<int64_t>(0, edge_pos + tau);
    holdoff_ = seed + cfg_.holdoff_len;
    probe_ = std::max(probe_ + hop_samples_, holdoff_);
    ++stats_.seeds_emitted;
    if (out_metric) *out_metric = par_full;
    return seed;
}

// ---------------------------------------------------------------------------
// QA / diagnostics: complete G_full evaluation of one window, no state change
// beyond the scratch accumulators.
// ---------------------------------------------------------------------------
bool PreambleXcorrAcquisition::evaluate_window(
    const std::vector<std::complex<int16_t>>& ring, int64_t valid_begin,
    int64_t win_start, bool force_scalar, PreambleXcorrDiagnostics* out)
{
    if (out) *out = PreambleXcorrDiagnostics{};
    if (ring.empty() || (ring.size() & (ring.size() - 1)) != 0) return false;
    // Out of bounds: the window starts before the valid data, or the ring
    // cannot hold one full evaluation window.
    if (win_start < 0 || win_start < valid_begin) return false;
    if (static_cast<uint64_t>(required_samples_) > ring.size()) return false;

    const std::complex<int16_t>* rp = ring.data();
    const size_t ring_size = ring.size();
    const size_t mask = ring_size - 1;
    const bool use_avx = !force_scalar && avx2_available();

    std::fill(acc_power_.begin(), acc_power_.end(), 0.0f);
    for (size_t g = 0; g < cfg_.g_full; ++g) {
        fold_group(rp, ring_size, mask, win_start, g, cfg_.k_coh, use_avx,
                   fold_re_.data(), fold_im_.data());
        correlate_sparse(fold_re_.data(), fold_im_.data(),
                         pos_taps_[target_code_].data(), n_pos_[target_code_],
                         neg_taps_[target_code_].data(), n_neg_[target_code_],
                         use_avx, acc_power_.data());
    }

    double peak = 0.0;
    double mean = 0.0;
    int64_t tau = 0;
    power_stats(acc_power_.data(), peak, mean, tau);
    // Same last-group fold the scan() gate uses; diagnostics only, the code
    // gates never make evaluate_window() fail.
    double group_par[kNumCodes] = {0.0, 0.0, 0.0, 0.0};
    code_group_pars(group_par, use_avx);
    double other_max = 0.0;
    for (size_t c = 0; c < kNumCodes; ++c) {
        if (c != target_code_ && group_par[c] > other_max) {
            other_max = group_par[c];
        }
    }
    const double target_par = group_par[target_code_];
    const double rho = code_identity(fold_re_.data(), fold_im_.data(),
                                     pos_taps_[target_code_].data(),
                                     n_pos_[target_code_],
                                     neg_taps_[target_code_].data(),
                                     n_neg_[target_code_], tau);
    if (out) {
        out->peak_power = peak;
        out->mean_power = mean;
        out->par = par_of(peak, mean);
        out->target_group_par = target_par;
        out->other_group_par_max = other_max;
        out->code_par_ratio = target_par / std::max(other_max, 1e-12);
        out->code_rho = rho;
        out->peak_tau = tau;
        out->groups_used = cfg_.g_full;
    }
    return true;
}

} // namespace live
} // namespace openuwb
