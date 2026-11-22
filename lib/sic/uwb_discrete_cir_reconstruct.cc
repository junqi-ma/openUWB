/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Discrete equivalent CIR full-packet reconstruction.
 *
 * Implements the reconstruction half of uwb_discrete_cir_sic.h:
 *
 *   ihat[j] = exp(j*w*(j-r)) * sum_{l in active} h[l] * (s*x)[j-q-l]
 *
 * r == q, x is the COMPLETE packet waveform produced by the public modulator,
 * and h is the frozen finite FIR from estimate_discrete_cir (h already
 * contains the packet amplitude because it was fitted against the public
 * shaped SYNC train).  Only the ACTIVE lags of the dense hull are multiplied;
 * the inactive hull taps are exactly zero by contract.
 *
 * The ONE waveform scale s is recomputed here from the public reference SYNC
 * and applied to the complete waveform in a single pass before the FIR; the
 * taps are never scaled again.  The direct O(N * active_count) FIR accumulates
 * in std::complex<float> and the CFO rotation uses a 1024-sample block-anchored
 * ABSOLUTE phase grid (double block base, float rotation), so no float phase
 * drift can accumulate.
 *
 * Sparse-pulse path (model.impulse_template + scratch.sparse_pulse_reconstruct)
 *
 * The exact discrete model is
 *
 *   y_hat[q + m_i + l] += sign_i * exp(j*w*(m_i + l)) * h[l]
 *
 * so the per-event phase exp(j*w*m_i) is combined with the per-lag phase
 * exp(j*w*l) folded into `taps_eff`.  That is mathematically identical to the
 * full-waveform rotation exp(j*w*(j-q)) at j = q + m_i + l: there is no
 * intra-tap CFO approximation.  With cfo == 0 the taps are used verbatim.
 *
 * L1D tiled reconstruction
 *
 * When the diagnostic ihat buffer is not requested, the reconstruction is
 * never materialized for the whole window.  A 4096-sample, 64-byte-aligned
 * tile (plus the 256-tap overlap tail) walks the clipped support and each tile
 * is immediately subtracted from out with the fused subtract/energy kernel.
 * Events are consumed in strictly increasing order and the carry tail is
 * copied forward, so every output sample receives the same float additions in
 * the same order as the full-buffer path: the two are BIT-IDENTICAL.
 *
 * Fail-closed: every failure leaves out == y byte-for-byte (unless the caller
 * explicitly disabled that copy) and ihat == 0.  Extreme amplitudes / unsafe
 * in-place transactions fall back to a two-phase scratch.accumulator path that
 * validates finite(recon) and finite(y - recon) before it writes out, so an
 * in-place window is never half-written.
 *
 * Nothing here reads sidecar/truth data and there is no fallback to the
 * fractional-path backend.
 */
#include <openuwb/sic/discrete_cir_sic.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <openuwb/phy/hrp/common/hrp_codec.h>
#include <openuwb/support/kernel_backend.h>
#include <openuwb/support/simd_config.h>
#include <openuwb/support/timing.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <volk/volk.h>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

namespace openuwb {
namespace sic {

static std::atomic<int> s_force_disable_b3(-1);

void set_b3_simd_disabled(int disabled)
{
    s_force_disable_b3.store(disabled, std::memory_order_relaxed);
}

bool is_b3_simd_disabled()
{
    const int forced = s_force_disable_b3.load(std::memory_order_relaxed);
    if (forced >= 0)
        return forced != 0;
    const char* env = std::getenv("UWB_SIC_DISABLE_B3");
    return (env != nullptr && *env != '\0' && *env != '0');
}

namespace {

using cf = std::complex<float>;
using cd = std::complex<double>;
using Clock = openuwb::timing::Clock;

constexpr double kPi = 3.14159265358979323846;
// Full-packet coordinates are frozen on the 998.4 MS/s work grid.  The frozen
// model carries a CFO in Hz but no rate, and the scratch carries only sps, so
// this public constant is the only admissible Hz -> rad/sample conversion.
// A model CFO at or beyond Nyquist is rejected instead of silently aliasing.
constexpr double kWorkRateHz = radar_meta::kWorkRateHz;
constexpr size_t kCfoBlock = 1024;
constexpr double kScaleImagMaxRel = 1e-6;
constexpr double kTemplateErrorMaxRel = 1e-10;
// L1D-resident reconstruction tile.  4096 complex samples = 32 KiB, plus the
// 256-tap overlap tail, which keeps the whole working set inside one 48 KiB
// L1D cache.
constexpr size_t kTileSamples = 4096;
// In-place amplitude bound.  A sample is touched by at most tap_count <= 256
// events, so |recon| <= 256 * max_tap_abs; 1e20 keeps ~17 decades of float
// head-room while still catching injected/overflowing inputs BEFORE an
// in-place window can be half-written.
constexpr float kTileAmplitudeLimit = 1e20f;

double us(Clock::time_point t)
{
    return std::chrono::duration<double, std::micro>(Clock::now() - t).count();
}

bool finite(cf v)
{
    return std::isfinite(v.real()) && std::isfinite(v.imag());
}

int64_t floordiv(int64_t a, int64_t b)
{
    return a / b - ((a % b) < 0 ? 1 : 0);
}

// Fused complex AXPY for one tap: dst[i] += h * src[i].
//
// This replaces multiply-into-tmp followed by add (two full passes over the
// destination) with a single read/modify/write pass, removing one vector
// write and one vector read per tap.  The scalar float form is deliberate:
// GCC lowers std::complex operator* through a NaN-preserving __mulsc3 call,
// which defeats vectorization, while this component form compiles to packed
// SSE with the project's baseline flags.  The helper stays out of line
// because inlining it into the large apply transaction suppressed the loop
// vectorizer; one call per active tap is negligible against O(N) samples.
// `dst` and `src` live in distinct preallocated scratch buffers, hence
// __restrict; per-sample tap order is unchanged because the caller still
// runs one tap at a time in ascending active order.
void __attribute__((noinline))
fused_axpy(cf* __restrict dst, const cf* __restrict src, cf h, size_t count)
{
    float* __restrict d = reinterpret_cast<float*>(dst);
    const float* __restrict s = reinterpret_cast<const float*>(src);
    const float hr = h.real();
    const float hi = h.imag();
    // No loop-carried dependence: each sample updates only its own two floats.
#pragma GCC ivdep
    for (size_t i = 0; i < count; ++i) {
        const float sr = s[2 * i];
        const float si = s[2 * i + 1];
        d[2 * i] += hr * sr - hi * si;
        d[2 * i + 1] += hi * sr + hr * si;
    }
}

// Zero-CFO dedicated accumulation.  A unit sign does not need a complex
// multiply; this preserves the old unit-real path's per-tap accumulation
// order while removing the two zero multiplies in fused_axpy.
void __attribute__((noinline))
fused_addsub(cf* __restrict dst, const cf* __restrict src, int sign,
             size_t count)
{
    float* __restrict d = reinterpret_cast<float*>(dst);
    const float* __restrict s = reinterpret_cast<const float*>(src);
    if (sign > 0) {
#pragma GCC ivdep
        for (size_t i = 0; i < count; ++i) {
            d[2 * i] += s[2 * i];
            d[2 * i + 1] += s[2 * i + 1];
        }
    } else {
#pragma GCC ivdep
        for (size_t i = 0; i < count; ++i) {
            d[2 * i] -= s[2 * i];
            d[2 * i + 1] -= s[2 * i + 1];
        }
    }
}

// Add one boundary-crossing pulse only where its taps land inside
// [win_lo, win_hi).  `dst[0]` corresponds to absolute sample `dst_origin`.
// This is the cold path used by the left/right clipped edges (and by the
// full-buffer diagnostic path when clipping is active); the hot tile path
// consumes whole events through the SIMD kernels so it stays bit-identical
// to the full-buffer path.
void add_pulse_fragment(cf* __restrict dst, int64_t dst_origin, int64_t pos,
                        const cf* __restrict taps, size_t tap_count, cf c,
                        int64_t win_lo, int64_t win_hi)
{
    const int64_t d0 = std::max<int64_t>(0, win_lo - pos);
    const int64_t d1 =
        std::min<int64_t>(int64_t(tap_count), win_hi - pos);
    const float cr = c.real();
    const float ci = c.imag();
    for (int64_t d = d0; d < d1; ++d) {
        const cf h = taps[static_cast<size_t>(d)];
        dst[pos + d - dst_origin] +=
            cf(cr * h.real() - ci * h.imag(),
               cr * h.imag() + ci * h.real());
    }
}

// Double-precision block phase shared by every CFO path.  All complex products
// are written out as four scalar multiplies: under the project's no-fast-math
// flags, std::complex operator* lowers to the out-of-line __muldc3/__mulsc3
// runtime (NaN/Inf handling), and 1024 such calls dominated the sparse-pulse
// setup time.  Adjacent blocks advance by the precomputed step; skips (and
// every 16th block) re-anchor with one sin/cos pair so no drift can build up.
struct CfoBlockState {
    // `kNone` cannot collide with a real block index: the non-sparse path
    // starts at q + lag_lo, so its first p can be negative and the first block
    // can legitimately be -1 (the old -1 sentinel silently skipped that
    // anchor and dropped a whole exp(j*w*1024) block).
    static constexpr int64_t kNone = std::numeric_limits<int64_t>::min();
    int64_t block = kNone;
    int since_anchor = 0;
    double base_r = 1.0, base_i = 0.0;
    double step_r = 1.0, step_i = 0.0;

    void anchor(int64_t b, double omega)
    {
        const double ph = omega * double(b * int64_t(kCfoBlock));
        base_r = std::cos(ph);
        base_i = std::sin(ph);
        since_anchor = 0;
        block = b;
    }

    void next_block(int64_t b, double omega)
    {
        if (block != kNone && b == block + 1 && ++since_anchor < 16) {
            const double r = base_r * step_r - base_i * step_i;
            const double i = base_r * step_i + base_i * step_r;
            base_r = r;
            base_i = i;
        } else {
            anchor(b, omega);
        }
        block = b;
    }
};

inline void rotor_mul(double& r, double& i, double sr, double si)
{
    const double nr = r * sr - i * si;
    const double ni = r * si + i * sr;
    r = nr;
    i = ni;
}

// rotation[k] = exp(j*omega*k) for k < kCfoBlock.  Four independent double
// recurrences advance by step^4, so the 1024-entry table has no long serial
// dependency chain and never calls std::complex/cexp.  One sin/cos pair plus
// 256 recurrence steps keeps the double error far below the float storage
// rounding, even for a 100 kHz offset.
inline void build_rotation_table(double omega,
                                 std::array<cf, kCfoBlock>& rotation)
{
    static_assert(kCfoBlock % 4 == 0, "rotation table is 4-way");
    const double c1 = std::cos(omega);
    const double s1 = std::sin(omega);
    const double c2 = c1 * c1 - s1 * s1;
    const double s2 = 2.0 * c1 * s1;
    const double c3 = c2 * c1 - s2 * s1;
    const double s3 = c2 * s1 + s2 * c1;
    const double c4 = c2 * c2 - s2 * s2;
    const double s4 = 2.0 * c2 * s2;

    double r0 = 1.0, i0 = 0.0;
    double r1 = c1, i1 = s1;
    double r2 = c2, i2 = s2;
    double r3 = c3, i3 = s3;
    size_t k = 0;
    for (; k + 4 <= kCfoBlock; k += 4) {
        rotation[k + 0] = cf(static_cast<float>(r0), static_cast<float>(i0));
        rotation[k + 1] = cf(static_cast<float>(r1), static_cast<float>(i1));
        rotation[k + 2] = cf(static_cast<float>(r2), static_cast<float>(i2));
        rotation[k + 3] = cf(static_cast<float>(r3), static_cast<float>(i3));
        rotor_mul(r0, i0, c4, s4);
        rotor_mul(r1, i1, c4, s4);
        rotor_mul(r2, i2, c4, s4);
        rotor_mul(r3, i3, c4, s4);
    }
    for (; k < kCfoBlock; ++k) {
        rotation[k] = cf(static_cast<float>(r0), static_cast<float>(i0));
        rotor_mul(r0, i0, c1, s1);
    }
}

#if defined(__x86_64__) || defined(_M_X64)
#if defined(__GNUC__) || defined(__clang__)

inline bool has_avx2_fma()
{
    static const bool s_has =
        __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    return s_has;
}

// ---------------------------------------------------------------------------
// Sparse-pulse accumulation kernels.
//
// Every kernel takes the destination window base `dst` and `delta` such that
// an event contributes at dst[delta + events[i].sample_offset + d], d in
// [0, tap_count).  The full-buffer path uses delta = lo (dst == ihat) or the
// clipped-window origin; the tile path uses delta = lo - tile_lo.  Because the
// arithmetic is element-wise and the events stay in ascending order, every
// kernel is bit-identical to the scalar reference for zero CFO.
// ---------------------------------------------------------------------------

// 20-tap zero-CFO fast path: five resident YMM registers, one add/sub per
// register and one store per event.
__attribute__((target("avx2")))
void accumulate_zero_cfo_20tap_avx2(
    cf* __restrict dst, const cf* __restrict taps,
    const mod::PulseEvent* __restrict events, size_t count, int64_t delta)
{
    const float* t_flt = reinterpret_cast<const float*>(taps);
    const __m256 t0 = _mm256_loadu_ps(t_flt);
    const __m256 t1 = _mm256_loadu_ps(t_flt + 8);
    const __m256 t2 = _mm256_loadu_ps(t_flt + 16);
    const __m256 t3 = _mm256_loadu_ps(t_flt + 24);
    const __m256 t4 = _mm256_loadu_ps(t_flt + 32);

    for (size_t i = 0; i < count; ++i) {
        const mod::PulseEvent& e = events[i];
        float* d = reinterpret_cast<float*>(
            dst + delta + int64_t(e.sample_offset));
        if (e.sign > 0) {
            _mm256_storeu_ps(d, _mm256_add_ps(_mm256_loadu_ps(d), t0));
            _mm256_storeu_ps(d + 8, _mm256_add_ps(_mm256_loadu_ps(d + 8), t1));
            _mm256_storeu_ps(d + 16,
                             _mm256_add_ps(_mm256_loadu_ps(d + 16), t2));
            _mm256_storeu_ps(d + 24,
                             _mm256_add_ps(_mm256_loadu_ps(d + 24), t3));
            _mm256_storeu_ps(d + 32,
                             _mm256_add_ps(_mm256_loadu_ps(d + 32), t4));
        } else {
            _mm256_storeu_ps(d, _mm256_sub_ps(_mm256_loadu_ps(d), t0));
            _mm256_storeu_ps(d + 8, _mm256_sub_ps(_mm256_loadu_ps(d + 8), t1));
            _mm256_storeu_ps(d + 16,
                             _mm256_sub_ps(_mm256_loadu_ps(d + 16), t2));
            _mm256_storeu_ps(d + 24,
                             _mm256_sub_ps(_mm256_loadu_ps(d + 24), t3));
            _mm256_storeu_ps(d + 32,
                             _mm256_sub_ps(_mm256_loadu_ps(d + 32), t4));
        }
    }
}

// Generic zero-CFO path for tap_count 1..256: four complex (one YMM) per
// iteration plus a scalar tail.  Removes the old silent fall-back to the
// noinline scalar fused_axpy for every non-20 tap width.
__attribute__((target("avx2")))
void accumulate_zero_cfo_generic_avx2(
    cf* __restrict dst, const cf* __restrict taps,
    const mod::PulseEvent* __restrict events, size_t count, int64_t delta,
    size_t tap_count)
{
    const float* __restrict t_flt = reinterpret_cast<const float*>(taps);
    for (size_t i = 0; i < count; ++i) {
        const mod::PulseEvent& e = events[i];
        float* __restrict d = reinterpret_cast<float*>(
            dst + delta + int64_t(e.sample_offset));
        size_t t = 0;
        if (e.sign > 0) {
            for (; t + 4 <= tap_count; t += 4)
                _mm256_storeu_ps(d + 2 * t,
                                 _mm256_add_ps(_mm256_loadu_ps(d + 2 * t),
                                               _mm256_loadu_ps(t_flt + 2 * t)));
            for (; t < tap_count; ++t) {
                d[2 * t] += t_flt[2 * t];
                d[2 * t + 1] += t_flt[2 * t + 1];
            }
        } else {
            for (; t + 4 <= tap_count; t += 4)
                _mm256_storeu_ps(d + 2 * t,
                                 _mm256_sub_ps(_mm256_loadu_ps(d + 2 * t),
                                               _mm256_loadu_ps(t_flt + 2 * t)));
            for (; t < tap_count; ++t) {
                d[2 * t] -= t_flt[2 * t];
                d[2 * t + 1] -= t_flt[2 * t + 1];
            }
        }
    }
}

// 20-tap CFO fast path.  The complex multiply is written by hand (never
// std::complex operator*, which GCC lowers to __mulsc3), the sign mask is
// pre-folded into the pre-permuted taps, and the block phase is a recurrence.
// `block_step` (= exp(j*w*1024)) is precomputed once per apply so the tile
// loop never pays a transcendental per tile.
__attribute__((target("avx2,fma")))
void accumulate_cfo_20tap_avx2(
    cf* __restrict dst, const cf* __restrict taps,
    const mod::PulseEvent* __restrict events, size_t count, int64_t delta,
    const std::array<cf, kCfoBlock>& rotation, double omega,
    const cd& block_step)
{
    static_assert(kCfoBlock == 1024, "block math assumes 1024 samples");
    const float* t_flt = reinterpret_cast<const float*>(taps);
    const __m256 h0 = _mm256_loadu_ps(t_flt);
    const __m256 h1 = _mm256_loadu_ps(t_flt + 8);
    const __m256 h2 = _mm256_loadu_ps(t_flt + 16);
    const __m256 h3 = _mm256_loadu_ps(t_flt + 24);
    const __m256 h4 = _mm256_loadu_ps(t_flt + 32);

    const __m256 sign_mask =
        _mm256_setr_ps(-1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f);
    // hs_k = permute(h_k) * sign_mask = (-imag, +real) per complex sample, so
    // the loop only needs the two FMAs.
    const __m256 hs0 = _mm256_mul_ps(_mm256_permute_ps(h0, 0xB1), sign_mask);
    const __m256 hs1 = _mm256_mul_ps(_mm256_permute_ps(h1, 0xB1), sign_mask);
    const __m256 hs2 = _mm256_mul_ps(_mm256_permute_ps(h2, 0xB1), sign_mask);
    const __m256 hs3 = _mm256_mul_ps(_mm256_permute_ps(h3, 0xB1), sign_mask);
    const __m256 hs4 = _mm256_mul_ps(_mm256_permute_ps(h4, 0xB1), sign_mask);

    CfoBlockState phase;
    phase.step_r = block_step.real();
    phase.step_i = block_step.imag();

    for (size_t i = 0; i < count; ++i) {
        const mod::PulseEvent& e = events[i];
        const uint32_t off = e.sample_offset;
        const int64_t b = int64_t(off >> 10);
        if (b != phase.block)
            phase.next_block(b, omega);
        const cf rot = rotation[off & 1023u];
        const float br = static_cast<float>(phase.base_r);
        const float bi = static_cast<float>(phase.base_i);
        float pr = br * rot.real() - bi * rot.imag();
        float pi = br * rot.imag() + bi * rot.real();
        if (e.sign < 0) {
            pr = -pr;
            pi = -pi;
        }
        const __m256 v_cr = _mm256_set1_ps(pr);
        const __m256 v_ci = _mm256_set1_ps(pi);
        float* d = reinterpret_cast<float*>(dst + delta + int64_t(off));

        __m256 d0 = _mm256_loadu_ps(d);
        d0 = _mm256_fmadd_ps(v_cr, h0, d0);
        d0 = _mm256_fmadd_ps(v_ci, hs0, d0);
        _mm256_storeu_ps(d, d0);

        __m256 d1 = _mm256_loadu_ps(d + 8);
        d1 = _mm256_fmadd_ps(v_cr, h1, d1);
        d1 = _mm256_fmadd_ps(v_ci, hs1, d1);
        _mm256_storeu_ps(d + 8, d1);

        __m256 d2 = _mm256_loadu_ps(d + 16);
        d2 = _mm256_fmadd_ps(v_cr, h2, d2);
        d2 = _mm256_fmadd_ps(v_ci, hs2, d2);
        _mm256_storeu_ps(d + 16, d2);

        __m256 d3 = _mm256_loadu_ps(d + 24);
        d3 = _mm256_fmadd_ps(v_cr, h3, d3);
        d3 = _mm256_fmadd_ps(v_ci, hs3, d3);
        _mm256_storeu_ps(d + 24, d3);

        __m256 d4 = _mm256_loadu_ps(d + 32);
        d4 = _mm256_fmadd_ps(v_cr, h4, d4);
        d4 = _mm256_fmadd_ps(v_ci, hs4, d4);
        _mm256_storeu_ps(d + 32, d4);
    }
}

// Generic CFO path for tap_count 1..256 (4 complex per YMM + scalar tail).
__attribute__((target("avx2,fma")))
void accumulate_cfo_generic_avx2(
    cf* __restrict dst, const cf* __restrict taps,
    const mod::PulseEvent* __restrict events, size_t count, int64_t delta,
    size_t tap_count, const std::array<cf, kCfoBlock>& rotation, double omega,
    const cd& block_step)
{
    static_assert(kCfoBlock == 1024, "block math assumes 1024 samples");
    const float* t_flt = reinterpret_cast<const float*>(taps);
    const __m256 sign_mask =
        _mm256_setr_ps(-1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f);

    CfoBlockState phase;
    phase.step_r = block_step.real();
    phase.step_i = block_step.imag();

    for (size_t i = 0; i < count; ++i) {
        const mod::PulseEvent& e = events[i];
        const uint32_t off = e.sample_offset;
        const int64_t b = int64_t(off >> 10);
        if (b != phase.block)
            phase.next_block(b, omega);
        const cf rot = rotation[off & 1023u];
        const float br = static_cast<float>(phase.base_r);
        const float bi = static_cast<float>(phase.base_i);
        float pr = br * rot.real() - bi * rot.imag();
        float pi = br * rot.imag() + bi * rot.real();
        if (e.sign < 0) {
            pr = -pr;
            pi = -pi;
        }
        const __m256 v_cr = _mm256_set1_ps(pr);
        const __m256 v_ci = _mm256_set1_ps(pi);
        float* d = reinterpret_cast<float*>(dst + delta + int64_t(off));

        size_t t = 0;
        for (; t + 4 <= tap_count; t += 4) {
            const __m256 h = _mm256_loadu_ps(t_flt + 2 * t);
            const __m256 hs =
                _mm256_mul_ps(_mm256_permute_ps(h, 0xB1), sign_mask);
            __m256 dd = _mm256_loadu_ps(d + 2 * t);
            dd = _mm256_fmadd_ps(v_cr, h, dd);
            dd = _mm256_fmadd_ps(v_ci, hs, dd);
            _mm256_storeu_ps(d + 2 * t, dd);
        }
        for (; t < tap_count; ++t) {
            const float hr = t_flt[2 * t];
            const float hi = t_flt[2 * t + 1];
            float* dp = d + 2 * t;
            const float dr = dp[0];
            const float di = dp[1];
            dp[0] = dr + (pr * hr - pi * hi);
            dp[1] = di + (pr * hi + pi * hr);
        }
    }
}

#if UWB_SIMD_HAVE_AVX512
// Production AVX-512 20-tap path: 40 floats = two 512-bit registers plus one
// 256-bit register, same per-element operation order as the AVX2 kernel.
__attribute__((target("avx512f,avx512dq,avx512vl")))
void accumulate_zero_cfo_20tap_avx512(
    cf* __restrict dst, const cf* __restrict taps,
    const mod::PulseEvent* __restrict events, size_t count, int64_t delta)
{
    const float* t_flt = reinterpret_cast<const float*>(taps);
    const __m512 t0 = _mm512_loadu_ps(t_flt);
    const __m512 t1 = _mm512_loadu_ps(t_flt + 16);
    const __m256 t2 = _mm256_loadu_ps(t_flt + 32);

    for (size_t i = 0; i < count; ++i) {
        const mod::PulseEvent& e = events[i];
        float* d = reinterpret_cast<float*>(
            dst + delta + int64_t(e.sample_offset));
        if (e.sign > 0) {
            _mm512_storeu_ps(d, _mm512_add_ps(_mm512_loadu_ps(d), t0));
            _mm512_storeu_ps(d + 16,
                             _mm512_add_ps(_mm512_loadu_ps(d + 16), t1));
            _mm256_storeu_ps(d + 32,
                             _mm256_add_ps(_mm256_loadu_ps(d + 32), t2));
        } else {
            _mm512_storeu_ps(d, _mm512_sub_ps(_mm512_loadu_ps(d), t0));
            _mm512_storeu_ps(d + 16,
                             _mm512_sub_ps(_mm512_loadu_ps(d + 16), t1));
            _mm256_storeu_ps(d + 32,
                             _mm256_sub_ps(_mm256_loadu_ps(d + 32), t2));
        }
    }
}

__attribute__((target("avx512f,avx512dq,avx512vl")))
void accumulate_cfo_20tap_avx512(
    cf* __restrict dst, const cf* __restrict taps,
    const mod::PulseEvent* __restrict events, size_t count, int64_t delta,
    const std::array<cf, kCfoBlock>& rotation, double omega,
    const cd& block_step)
{
    static_assert(kCfoBlock == 1024, "block math assumes 1024 samples");
    const float* t_flt = reinterpret_cast<const float*>(taps);
    const __m512 h0 = _mm512_loadu_ps(t_flt);
    const __m512 h1 = _mm512_loadu_ps(t_flt + 16);
    const __m256 h2 = _mm256_loadu_ps(t_flt + 32);

    const __m512 sign_mask512 = _mm512_setr_ps(
        -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f,
        -1.f, 1.f, -1.f, 1.f);
    const __m256 sign_mask256 =
        _mm256_setr_ps(-1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f);
    const __m512 hs0 =
        _mm512_mul_ps(_mm512_permute_ps(h0, 0xB1), sign_mask512);
    const __m512 hs1 =
        _mm512_mul_ps(_mm512_permute_ps(h1, 0xB1), sign_mask512);
    const __m256 hs2 =
        _mm256_mul_ps(_mm256_permute_ps(h2, 0xB1), sign_mask256);

    CfoBlockState phase;
    phase.step_r = block_step.real();
    phase.step_i = block_step.imag();

    for (size_t i = 0; i < count; ++i) {
        const mod::PulseEvent& e = events[i];
        const uint32_t off = e.sample_offset;
        const int64_t b = int64_t(off >> 10);
        if (b != phase.block)
            phase.next_block(b, omega);
        const cf rot = rotation[off & 1023u];
        const float br = static_cast<float>(phase.base_r);
        const float bi = static_cast<float>(phase.base_i);
        float pr = br * rot.real() - bi * rot.imag();
        float pi = br * rot.imag() + bi * rot.real();
        if (e.sign < 0) {
            pr = -pr;
            pi = -pi;
        }
        const __m512 v_cr = _mm512_set1_ps(pr);
        const __m512 v_ci = _mm512_set1_ps(pi);
        const __m256 v_cr256 = _mm256_set1_ps(pr);
        const __m256 v_ci256 = _mm256_set1_ps(pi);
        float* d = reinterpret_cast<float*>(dst + delta + int64_t(off));

        __m512 d0 = _mm512_loadu_ps(d);
        d0 = _mm512_fmadd_ps(v_cr, h0, d0);
        d0 = _mm512_fmadd_ps(v_ci, hs0, d0);
        _mm512_storeu_ps(d, d0);

        __m512 d1 = _mm512_loadu_ps(d + 16);
        d1 = _mm512_fmadd_ps(v_cr, h1, d1);
        d1 = _mm512_fmadd_ps(v_ci, hs1, d1);
        _mm512_storeu_ps(d + 16, d1);

        __m256 d2 = _mm256_loadu_ps(d + 32);
        d2 = _mm256_fmadd_ps(v_cr256, h2, d2);
        d2 = _mm256_fmadd_ps(v_ci256, hs2, d2);
        _mm256_storeu_ps(d + 32, d2);
    }
}
#endif // UWB_SIMD_HAVE_AVX512

// Scale + finite check in one AVX2 pass.  Bit-identical to
// `dst[i] = src[i] * scale` (element-wise float multiply) and safe when src
// and dst are the same buffer.
__attribute__((target("avx2")))
bool scale_copy_checked_avx2(const cf* src, cf* dst, size_t count, float scale)
{
    const float* s = reinterpret_cast<const float*>(src);
    float* d = reinterpret_cast<float*>(dst);
    const size_t floats = 2 * count;
    const __m256 vs = _mm256_set1_ps(scale);
    const __m256 sign = _mm256_set1_ps(-0.0f);
    const __m256 inf = _mm256_set1_ps(std::numeric_limits<float>::infinity());
    size_t i = 0;
    for (; i + 8 <= floats; i += 8) {
        const __m256 v = _mm256_loadu_ps(s + i);
        const __m256 mag = _mm256_andnot_ps(sign, v);
        const __m256 fin = _mm256_cmp_ps(mag, inf, _CMP_LT_OQ);
        if (_mm256_movemask_ps(fin) != 0xFF)
            return false;
        _mm256_storeu_ps(d + i, _mm256_mul_ps(v, vs));
    }
    for (; i < floats; ++i) {
        const float v = s[i];
        if (!std::isfinite(v))
            return false;
        d[i] = v * scale;
    }
    return true;
}

#else
inline bool has_avx2_fma() { return false; }
#endif
#else
inline bool has_avx2_fma() { return false; }
#endif

// Dispatch the sparse-pulse accumulation for `count` events into a window
// whose base is `dst`; event i lands at dst[delta + events[i].sample_offset].
// `taps` is already the exact-CFO `taps_eff` array (or model.taps for zero
// CFO).  set_b3_simd_disabled(1) forces the scalar reference so the QA parity
// comparison stays real.
void accumulate_pulses_range(cf* dst, int64_t delta,
                             const mod::PulseEvent* events, size_t count,
                             const cf* taps, size_t tap_count, double cfo_hz,
                             double omega,
                             const std::array<cf, kCfoBlock>& rotation,
                             const cd& block_step, bool disable_simd)
{
    if (count == 0)
        return;
    const bool zero_cfo = (cfo_hz == 0.0);
#if (defined(__x86_64__) || defined(_M_X64)) && \
    (defined(__GNUC__) || defined(__clang__))
    if (!disable_simd && has_avx2_fma()) {
#if UWB_SIMD_HAVE_AVX512
        if (tap_count == 20 &&
            demod::rxopt::use_avx512(core::simd::KernelGroup::Sync)) {
            if (zero_cfo)
                accumulate_zero_cfo_20tap_avx512(dst, taps, events, count,
                                                 delta);
            else
                accumulate_cfo_20tap_avx512(dst, taps, events, count, delta,
                                            rotation, omega, block_step);
            return;
        }
#endif
        if (tap_count == 20) {
            if (zero_cfo)
                accumulate_zero_cfo_20tap_avx2(dst, taps, events, count,
                                               delta);
            else
                accumulate_cfo_20tap_avx2(dst, taps, events, count, delta,
                                          rotation, omega, block_step);
        } else if (zero_cfo) {
            accumulate_zero_cfo_generic_avx2(dst, taps, events, count, delta,
                                             tap_count);
        } else {
            accumulate_cfo_generic_avx2(dst, taps, events, count, delta,
                                        tap_count, rotation, omega,
                                        block_step);
        }
        return;
    }
#endif
    if (zero_cfo) {
        // Zero CFO: sign_i is real; the dedicated add/sub keeps the old
        // unit-real accumulation order and is bit-comparable to it.
        for (size_t i = 0; i < count; ++i) {
            const mod::PulseEvent& e = events[i];
            fused_addsub(dst + delta + int64_t(e.sample_offset), taps,
                         int(e.sign), tap_count);
        }
    } else {
        CfoBlockState phase;
        phase.step_r = block_step.real();
        phase.step_i = block_step.imag();
        for (size_t i = 0; i < count; ++i) {
            const mod::PulseEvent& e = events[i];
            const uint32_t off = e.sample_offset;
            const int64_t b = int64_t(off >> 10);
            if (b != phase.block)
                phase.next_block(b, omega);
            const cf rot = rotation[off & 1023u];
            const float br = static_cast<float>(phase.base_r);
            const float bi = static_cast<float>(phase.base_i);
            float pr = br * rot.real() - bi * rot.imag();
            float pi = br * rot.imag() + bi * rot.real();
            if (e.sign < 0) {
                pr = -pr;
                pi = -pi;
            }
            fused_axpy(dst + delta + int64_t(off), taps, cf(pr, pi),
                       tap_count);
        }
    }
}

// Fail-closed model gate shared by apply and the transactional wrapper.  A
// frozen model is usable only with an integer origin, the matching phase
// anchor, a finite hull strictly shorter than one SYNC period, a finite CFO
// inside Nyquist, a non-empty ACTIVE set that is sorted, inside the hull and
// backed by finite dense taps.  This reads no estimator state.
const char* model_failure(const DiscreteCirModel& m, size_t sps, size_t n)
{
    if (!m.valid || sps < 2)
        return "invalid_model";
    if (m.q < 0 || m.q >= int64_t(n))
        return "invalid_model";
    if (m.phase_reference_sample != m.q)
        return "invalid_model";
    if (m.tap_count == 0 || m.tap_count > kDiscreteCirMaxTaps)
        return "invalid_model";
    if (m.tap_count >= sps)
        return "invalid_model";
    if (m.lag_lo < -int64_t(sps) || m.lag_lo > int64_t(sps))
        return "invalid_model";
    if (m.lag_lo + int64_t(m.tap_count) > int64_t(sps))
        return "invalid_model";
    if (!std::isfinite(m.cfo_hz) || !(std::fabs(m.cfo_hz) < 0.5 * kWorkRateHz))
        return "invalid_model";
    // The dense hull is zero outside the active groups, so every hull entry
    // must be finite even though only the active ones are multiplied.
    for (size_t i = 0; i < m.tap_count; ++i)
        if (!finite(m.taps[i]))
            return "invalid_model";
    if (m.active_count == 0 || m.active_count > kDiscreteCirMaxActive ||
        m.active_count > m.tap_count)
        return "invalid_model";
    for (size_t a = 0; a < m.active_count; ++a) {
        const int64_t lag = int64_t(m.active_lags[a]);
        if (lag < m.lag_lo || lag >= m.lag_lo + int64_t(m.tap_count))
            return "invalid_model";
        if (a && lag <= int64_t(m.active_lags[a - 1]))
            return "invalid_model"; // must be strictly increasing
        if (!finite(m.taps[static_cast<size_t>(lag - m.lag_lo)]))
            return "invalid_model";
    }
    return nullptr;
}

// Capacity contract.  prepare_discrete_cir/reserve own every allocation on
// this path, so an undersized buffer is a hard failure, never a hidden resize.
// `max_packet_samples` is only the estimator's planning bound: the hot call is
// sized by the ACTUAL waveform_n/packet_n, so a stale or zero bound cannot
// reject a waveform that fits the real buffers.
const char* apply_scratch_failure(const DiscreteCirScratch& s,
                                  size_t n,
                                  size_t waveform_n)
{
    if (!s.prepared)
        return "scratch_unprepared";
    if (s.sps == 0 || s.max_window_samples < n)
        return "scratch_capacity";
    if (s.scaled_waveform.size() < waveform_n ||
        s.accumulator.size() < n || s.fir_tmp.size() < n)
        return "scratch_capacity";
    return nullptr;
}

const char* transaction_scratch_failure(const DiscreteCirScratch& s, size_t n)
{
    if (!s.prepared)
        return "scratch_unprepared";
    if (s.sps == 0 || s.max_window_samples < n)
        return "scratch_capacity";
    // modulate_one() is sized by the real buffer length; an empty buffer can
    // never hold a packet.
    if (s.waveform.empty())
        return "scratch_capacity";
    // encode_phr19() resizes its output vector to 19 entries.  Requiring the
    // capacity up front keeps the transactional path allocation-free.
    if (s.expected_phr.capacity() < 19)
        return "scratch_capacity";
    return nullptr;
}

// Sparse event-list gate.  It validates the representation itself (strictly
// increasing unit signs, complete packet length) and checks the first SYNC
// against the prepared unit template.  It never regenerates the full packet
// waveform.  The linear scan is over events, not samples.
// For internal just-generated pulse list transactions (transaction_prepared),
// the full event order/sign scan is skipped while prefix SYNC match and bounds
// are verified.
const char* pulse_list_failure(const mod::PulseList& p,
                               size_t n,
                               const DiscreteCirScratch& s,
                               size_t sync_n,
                               bool transaction_prepared = false)
{
    if (!s.prepared || s.sps == 0 ||
        sync_n != s.sps || s.impulse_sync.size() != sync_n ||
        s.max_window_samples < n)
        return "scratch_capacity";
    if (!s.sparse_pulse_reconstruct || !s.ipatov_correlation)
        return "sparse_pulse_unprepared";
    if (!p.events || p.count == 0 || p.packet_samples == 0 ||
        p.packet_samples > mod::kMaxHrpTxSamples)
        return "invalid_pulse_list";

    if (transaction_prepared) {
        // Trusted internal transaction: skip full event order/sign scan.
        // Still verify the prefix SYNC match and basic packet bounds.
        size_t first_sync_events = 0;
        for (size_t i = 0; i < p.count; ++i) {
            const mod::PulseEvent& e = p.events[i];
            if (e.sample_offset >= sync_n)
                break;
            const cf expected = s.impulse_sync[e.sample_offset];
            if (expected.imag() != 0.f ||
                expected.real() != static_cast<float>(e.sign))
                return "pulse_template_mismatch";
            ++first_sync_events;
        }
        size_t first_sync_nonzero = 0;
        for (size_t i = 0; i < sync_n; ++i)
            if (s.impulse_sync[i] != cf(0.f, 0.f))
                ++first_sync_nonzero;
        if (first_sync_events != first_sync_nonzero)
            return "pulse_template_mismatch";
        return nullptr;
    }

    uint32_t last = 0;
    size_t first_sync_events = 0;
    for (size_t i = 0; i < p.count; ++i) {
        const mod::PulseEvent& e = p.events[i];
        if (e.sign != 1 && e.sign != -1)
            return "invalid_pulse_list";
        if (e.sample_offset >= p.packet_samples)
            return "invalid_pulse_list";
        if (i && e.sample_offset <= last)
            return "invalid_pulse_list";
        last = e.sample_offset;
        if (e.sample_offset < sync_n) {
            const cf expected = s.impulse_sync[e.sample_offset];
            if (expected.imag() != 0.f ||
                expected.real() != static_cast<float>(e.sign))
                return "pulse_template_mismatch";
            ++first_sync_events;
        }
    }
    size_t first_sync_nonzero = 0;
    for (size_t i = 0; i < sync_n; ++i)
        if (s.impulse_sync[i] != cf(0.f, 0.f))
            ++first_sync_nonzero;
    if (first_sync_events != first_sync_nonzero)
        return "pulse_template_mismatch";
    return nullptr;
}

} // namespace

static bool apply_discrete_cir_impl(
    const std::complex<float>* y,
    size_t n,
    const std::complex<float>* waveform,
    size_t waveform_n,
    const std::complex<float>* reference_sync,
    size_t sync_n,
    const DiscreteCirModel& model,
    int64_t roi_lo,
    int64_t roi_hi,
    const DiscreteCirApplyOptions& opts,
    DiscreteCirScratch& scratch,
    std::complex<float>* out,
    std::complex<float>* ihat,
    DiscreteCirResult* result,
    const float* real_impulses)
{
    const auto t_apply = Clock::now();
    DiscreteCirResult local;
    DiscreteCirResult& r = result ? *result : local;
    r = DiscreteCirResult();
    auto fail = [&](const char* reason) {
        if (opts.copy_out_on_fail && out && y && n > 0 && out != y)
            std::copy(y, y + n, out);
        if (ihat && ihat != y && ihat != out)
            std::fill(ihat, ihat + n, cf(0.f, 0.f));
        r.committed = false;
        r.reason = reason;
        r.apply_total_us = us(t_apply);
        return false;
    };

    if (!y || !out || n == 0 || (ihat && (ihat == y || ihat == out)))
        return fail("invalid_buffers");
    if ((!waveform && !real_impulses) || !reference_sync || waveform_n == 0 ||
        sync_n == 0 || waveform_n < sync_n)
        return fail("invalid_buffers");

    const char* sf = apply_scratch_failure(scratch, n, waveform_n);
    if (sf)
        return fail(sf);
    const char* mf = model_failure(model, scratch.sps, n);
    if (mf)
        return fail(mf);
    if (roi_lo < 0 || roi_hi > int64_t(n) || roi_hi <= roi_lo)
        return fail("invalid_roi");
    if (sync_n != scratch.sps)
        return fail("invalid_template");
    if (!opts.transaction_prepared)
        for (size_t i = 0; i < n; ++i)
            if (!finite(y[i]))
                return fail("nonfinite_input");

    if (model.impulse_template != scratch.ipatov_correlation ||
        (real_impulses && !model.impulse_template))
        return fail("template_kind_mismatch");
    if (model.impulse_template) {
        if (scratch.impulse_sync.size() != sync_n)
            return fail("impulse_template_unprepared");
        reference_sync = scratch.impulse_sync.data();
    }

    r.model_id = model.model_id;
    r.q = model.q;
    r.cfo_hz = model.cfo_hz;
    r.tap_count = model.tap_count;
    r.fit_residual = model.fit_residual;
    r.gain_condition = model.gain_condition;

    const int64_t q = model.q;
    const int64_t lo = q + model.lag_lo;
    const int64_t hi = lo + int64_t(model.tap_count) + int64_t(waveform_n);
    r.support_lo = lo;
    r.support_hi = hi;

    const auto t_channel = Clock::now();

    // --- the ONE waveform scale, against the public reference SYNC ----------
    const auto t_scale = Clock::now();
    cf* const scaled = scratch.scaled_waveform.data();
    if (real_impulses) {
        // Private input from modulate_unit_impulses: validated unit pulse
        // profile + integer BPRF symbols guarantee finite {-1,0,+1} samples.
        // Match the public unit reference without constructing/scaling a
        // complete complex waveform. Absolute amplitude remains in the CIR.
        for (size_t i=0; i<sync_n; ++i)
            if (real_impulses[i] != reference_sync[i].real() ||
                reference_sync[i].imag() != 0.f)
                return fail("template_scale_mismatch");
        r.waveform_scale = 1.0;
    } else {
        cd cross(0.0, 0.0);
        double energy = 0.0, ref_energy = 0.0;
        for (size_t i = 0; i < sync_n; ++i) {
            if (!finite(waveform[i]) || !finite(reference_sync[i]))
                return fail("nonfinite_input");
            // conj(w) * ref, written out to avoid __muldc3 (1024 calls here
            // used to dominate the non-sparse scale stage).
            const double wr = waveform[i].real();
            const double wi = waveform[i].imag();
            const double rr = reference_sync[i].real();
            const double ri = reference_sync[i].imag();
            cross += cd(wr * rr + wi * ri, wr * ri - wi * rr);
            energy += wr * wr + wi * wi;
            ref_energy += rr * rr + ri * ri;
        }
        const cd scale = energy > 0.0 ? cross / energy : cd(0.0, 0.0);
        const double scale_real = scale.real();
        if (!(scale_real > 0.0) || !std::isfinite(scale_real) ||
            !std::isfinite(scale.imag()) ||
            std::fabs(scale.imag()) > kScaleImagMaxRel * scale_real ||
            !(ref_energy > 0.0))
            return fail("template_scale_mismatch");
        double template_error = 0.0;
        const double sr = scale.real();
        const double si = scale.imag();
        for (size_t i = 0; i < sync_n; ++i) {
            const double wr = waveform[i].real();
            const double wi = waveform[i].imag();
            const double pr = wr * sr - wi * si;
            const double pi = wr * si + wi * sr;
            const double dr = pr - reference_sync[i].real();
            const double di = pi - reference_sync[i].imag();
            template_error += dr * dr + di * di;
        }
        if (!std::isfinite(template_error) ||
            template_error > kTemplateErrorMaxRel * ref_energy)
            return fail("template_scale_mismatch");
        // Pre-scale the COMPLETE waveform exactly once.  h already contains the
        // packet amplitude, so the taps are never scaled again.  The AVX2 pass
        // fuses the finite check with the float multiply; it reads src before
        // storing dst, so it also stays correct if the caller handed in the
        // same scratch buffer.
        const float scale_f = static_cast<float>(scale_real);
        bool finite_input = true;
#if (defined(__x86_64__) || defined(_M_X64)) && \
    (defined(__GNUC__) || defined(__clang__))
        if (has_avx2_fma()) {
            finite_input =
                scale_copy_checked_avx2(waveform, scaled, waveform_n, scale_f);
        } else
#endif
        {
            for (size_t i = 0; i < waveform_n; ++i) {
                const cf v = waveform[i];
                if (!finite(v)) {
                    finite_input = false;
                    break;
                }
                scaled[i] = v * scale_f;
            }
        }
        if (!finite_input)
            return fail("nonfinite_input");
        r.waveform_scale = scale_real;
    }
    r.scale_us = us(t_scale);

    // --- support, optional clipping and consumer ROI ------------------------
    // By default a complete observed packet is mandatory; truncation never
    // clips.  allow_partial_clip opts into reconstructing only the in-window
    // portion (edge interferers in a sliding capture window).
    int64_t clip_lo = lo, clip_hi = hi;
    if (!opts.allow_partial_clip) {
        if (lo < 0 || hi > int64_t(n)) {
            r.channel_us = us(t_channel);
            return fail("context_truncated");
        }
    } else {
        clip_lo = std::max<int64_t>(0, lo);
        clip_hi = std::min<int64_t>(int64_t(n), hi);
        if (clip_hi <= clip_lo) {
            r.channel_us = us(t_channel);
            return fail("no_overlap");
        }
    }
    const int64_t sub_lo = std::max(clip_lo, roi_lo);
    const int64_t sub_hi = std::min(clip_hi, roi_hi);
    if (sub_hi <= sub_lo) {
        r.channel_us = us(t_channel);
        return fail("no_overlap");
    }
    r.sub_lo = sub_lo;
    r.sub_hi = sub_hi;
    const bool clipped = (clip_lo != lo || clip_hi != hi);

    if (ihat) {
        // Diagnostics: zero the mathematical support flanks (everything the
        // rotation below does not write).
        if (clipped) {
            if (clip_lo > 0)
                std::fill(ihat, ihat + clip_lo, cf(0.f, 0.f));
            if (clip_hi < int64_t(n))
                std::fill(ihat + clip_hi, ihat + n, cf(0.f, 0.f));
        } else {
            if (lo > 0)
                std::fill(ihat, ihat + lo, cf(0.f, 0.f));
            if (hi < int64_t(n))
                std::fill(ihat + hi, ihat + n, cf(0.f, 0.f));
        }
    }

    // --- direct O(N * active_count) FIR, tap-major contiguous passes ------
    // Each active tap is one contiguous run over j.  The inner loop is a
    // vectorized fused complex multiply-add with no per-sample bounds branch;
    // the accumulator lives in preallocated scratch, so the hot path never
    // allocates.  Only the active lags contribute (inactive hull taps are 0).
    // The accumulator is indexed relative to clip_lo so a clipped support that
    // extends before sample 0 cannot underflow it.
    const auto t_conv = Clock::now();
    cf* const recon = ihat ? ihat : scratch.fir_tmp.data();
    const int64_t x_n = int64_t(waveform_n);
    const int64_t win_n = clip_hi - clip_lo;
    cf* const acc = scratch.accumulator.data();
    std::fill(acc, acc + win_n, cf(0.f, 0.f));
    const size_t active = model.active_count;
    if (model.impulse_template) {
        // The unit impulse packet is sparse. Add a contiguous response per
        // nonzero impulse, avoiding N * response_width zero multiplications.
        for (int64_t i=0; i<x_n; ++i) {
            const cf impulse = real_impulses ? cf(real_impulses[i],0.f) : scaled[i];
            if (impulse==cf{}) continue;
            const int64_t base=q+i+model.lag_lo;
            if (clipped)
                add_pulse_fragment(acc, clip_lo, base, model.taps.data(),
                                   model.tap_count, impulse, clip_lo, clip_hi);
            else
                fused_axpy(acc+(base-clip_lo),model.taps.data(),impulse,
                           model.tap_count);
        }
    } else for (size_t a = 0; a < active; ++a) {
        const int64_t lag = int64_t(model.active_lags[a]);
        const cf h =
            model.taps[static_cast<size_t>(lag - model.lag_lo)];
        const int64_t j0 = std::max<int64_t>(clip_lo, q + lag);
        const int64_t j1 = std::min<int64_t>(clip_hi, q + lag + x_n);
        if (j1 <= j0)
            continue;
        const cf* const src = scaled + (j0 - q - lag);
        cf* const dst = acc + (j0 - clip_lo);
        fused_axpy(dst, src, h, static_cast<size_t>(j1 - j0));
    }
    r.convolution_us = us(t_conv);

    // --- block-anchored absolute CFO grid and finite checks -----------------
    const auto t_cfo = Clock::now();
    alignas(32) std::array<cf, kCfoBlock> tmp_block;
    cf* const tmp = tmp_block.data(); // shared scalar-multiply temp
    const double omega = 2.0 * kPi * model.cfo_hz / kWorkRateHz;
    std::array<cf, kCfoBlock> rotation;
    build_rotation_table(omega, rotation);
    const int64_t anchor = model.phase_reference_sample; // == q by the gate
    CfoBlockState bstate;
    bstate.step_r = std::cos(omega * double(kCfoBlock));
    bstate.step_i = std::sin(omega * double(kCfoBlock));
    for (int64_t j = clip_lo; j < clip_hi;) {
        const int64_t p = j - anchor;
        const int64_t b = floordiv(p, int64_t(kCfoBlock));
        const size_t j0 = static_cast<size_t>(p - b * int64_t(kCfoBlock));
        const size_t count =
            std::min<size_t>(kCfoBlock - j0, static_cast<size_t>(clip_hi - j));
        if (b != bstate.block)
            bstate.next_block(b, omega);
        const cf block(static_cast<float>(bstate.base_r),
                       static_cast<float>(bstate.base_i));
        const cf* const src = acc + (j - clip_lo);
        const unsigned nrot = static_cast<unsigned>(count);
        volk_32fc_s32fc_multiply_32fc(tmp, rotation.data() + j0, block, nrot);
        volk_32fc_x2_multiply_32fc(recon + j, src, tmp, nrot);
        j += static_cast<int64_t>(count);
    }
    for (int64_t j = clip_lo; j < clip_hi; ++j) {
        if (!finite(recon[j]) || !finite(y[j] - recon[j])) {
            r.cfo_us = us(t_cfo);
            r.channel_us = us(t_channel);
            return fail("nonfinite_reconstruction");
        }
    }
    r.cfo_us = us(t_cfo);
    r.channel_us = us(t_channel);

    // --- subtract only the consumer ROI ------------------------------------
    const auto t_subtract = Clock::now();
    if (out != y) {
        if (sub_lo > 0)
            std::copy(y, y + sub_lo, out);
        if (sub_hi < int64_t(n))
            std::copy(y + sub_hi, y + n, out + sub_hi);
    }
    // Element-wise out = y - recon.  Routed through the shared D1 element-wise
    // kernels (bit-identical; AVX-512 when the Sync group selects it) instead
    // of a bare loop the compiler only widens to 256 bits.  When ROI energies
    // are requested the fused kernel computes them before out is written, so
    // the result is also correct for out == y (in-place cancellation).
    const size_t sub_count = static_cast<size_t>(sub_hi - sub_lo);
    if (opts.compute_roi_energy) {
        double e_y = 0.0, e_out = 0.0;
        demod::rxopt::vector_sub3_energy_complex(out + sub_lo, y + sub_lo,
                                                 recon + sub_lo, sub_count,
                                                 e_y, e_out);
        r.roi_energy_valid = true;
        r.roi_e_y = e_y;
        r.roi_e_out = e_out;
    } else {
        demod::rxopt::vector_sub3_complex(out + sub_lo, y + sub_lo,
                                          recon + sub_lo, sub_count);
    }
    r.subtract_us = us(t_subtract);

    r.packet_samples = waveform_n;
    r.committed = true;
    r.reason = "committed";
    r.apply_total_us = us(t_apply);
    return true;
}


bool apply_discrete_cir(const std::complex<float>* y, size_t n,
                        const std::complex<float>* waveform, size_t waveform_n,
                        const std::complex<float>* reference_sync, size_t sync_n,
                        const DiscreteCirModel& model, int64_t roi_lo, int64_t roi_hi,
                        DiscreteCirScratch& scratch, std::complex<float>* out,
                        std::complex<float>* ihat, DiscreteCirResult* result)
{
    const DiscreteCirApplyOptions opts;
    return apply_discrete_cir_impl(y,n,waveform,waveform_n,reference_sync,sync_n,
                                  model,roi_lo,roi_hi,opts,scratch,out,ihat,result,
                                  nullptr);
}

static bool apply_discrete_cir_pulses_impl(
    const std::complex<float>* y,
    size_t n,
    const mod::PulseList& pulses,
    const std::complex<float>* reference_sync,
    size_t sync_n,
    const DiscreteCirModel& model,
    int64_t roi_lo,
    int64_t roi_hi,
    const DiscreteCirApplyOptions& opts,
    DiscreteCirScratch& scratch,
    std::complex<float>* out,
    std::complex<float>* ihat,
    DiscreteCirResult* result)
{
    const auto t_apply = Clock::now();
    DiscreteCirResult local;
    DiscreteCirResult& r = result ? *result : local;
    r = DiscreteCirResult();
    auto fail = [&](const char* reason) {
        if (opts.copy_out_on_fail && out && y && n > 0 && out != y)
            std::copy(y, y + n, out);
        if (ihat && ihat != y && ihat != out)
            std::fill(ihat, ihat + n, cf(0.f, 0.f));
        r.committed = false;
        r.reason = reason;
        r.apply_total_us = us(t_apply);
        return false;
    };

    if (!y || !out || n == 0 || (ihat && (ihat == y || ihat == out)))
        return fail("invalid_buffers");

    const char* pf = pulse_list_failure(pulses, n, scratch, sync_n,
                                        opts.transaction_prepared);
    if (pf)
        return fail(pf);
    // The two-phase fallback accumulates in scratch.accumulator; keeping the
    // capacity contract here also keeps the whole path allocation-free.
    if (!ihat && scratch.accumulator.size() < n)
        return fail("scratch_capacity");
    const char* mf = model_failure(model, scratch.sps, n);
    if (mf)
        return fail(mf);
    if (!reference_sync || sync_n != scratch.sps ||
        scratch.impulse_sync.size() != sync_n)
        return fail("invalid_template");
    if (!model.impulse_template || !scratch.ipatov_correlation)
        return fail("template_kind_mismatch");
    if (roi_lo < 0 || roi_hi > int64_t(n) || roi_hi <= roi_lo)
        return fail("invalid_roi");

    // Scan y once (also records the amplitude bound that gates the in-place
    // tile path).  A transaction-prepared caller already guaranteed finiteness.
    float max_abs_y = 0.f;
    if (!opts.transaction_prepared) {
        for (size_t i = 0; i < n; ++i) {
            if (!finite(y[i]))
                return fail("nonfinite_input");
            max_abs_y = std::max(max_abs_y, std::abs(y[i]));
        }
    }

    r.model_id = model.model_id;
    r.q = model.q;
    r.cfo_hz = model.cfo_hz;
    r.tap_count = model.tap_count;
    r.fit_residual = model.fit_residual;
    r.gain_condition = model.gain_condition;
    r.pulse_event_count = pulses.count;
    r.pulse_packet_samples = pulses.packet_samples;

    const int64_t q = model.q;
    const size_t tap_count = model.tap_count;
    const int64_t lo = q + model.lag_lo;
    const int64_t hi = lo + int64_t(tap_count) +
                       int64_t(pulses.packet_samples);
    r.support_lo = lo;
    r.support_hi = hi;

    const auto t_channel = Clock::now();

    // --- unit template identity and support/ROI window ---------------------
    const auto t_scale = Clock::now();
    // The pulse list's first SYNC was checked against impulse_sync above; the
    // public reference is intentionally bypassed here, exactly like the unit
    // real path.  Absolute amplitude lives in the frozen CIR.
    r.waveform_scale = 1.0;
    r.scale_us = us(t_scale);

    int64_t clip_lo = lo, clip_hi = hi;
    if (!opts.allow_partial_clip) {
        if (lo < 0 || hi > int64_t(n)) {
            r.channel_us = us(t_channel);
            return fail("context_truncated");
        }
    } else {
        clip_lo = std::max<int64_t>(0, lo);
        clip_hi = std::min<int64_t>(int64_t(n), hi);
        if (clip_hi <= clip_lo) {
            r.channel_us = us(t_channel);
            return fail("no_overlap");
        }
    }
    const int64_t sub_lo = std::max(clip_lo, roi_lo);
    const int64_t sub_hi = std::min(clip_hi, roi_hi);
    if (sub_hi <= sub_lo) {
        r.channel_us = us(t_channel);
        return fail("no_overlap");
    }
    r.sub_lo = sub_lo;
    r.sub_hi = sub_hi;
    const bool clipped = (clip_lo != lo || clip_hi != hi);

    const double omega = 2.0 * kPi * model.cfo_hz / kWorkRateHz;

    // --- exact per-lag CFO pre-rotation ------------------------------------
    // Source contribution at j = q + m_i + l is
    //   sign_i * exp(j*w*(m_i + l)) * h[l]
    // = (sign_i * exp(j*w*m_i)) * taps_eff[l - lag_lo].
    // With cfo == 0 the taps are used verbatim (bit-identical legacy).
    // Resolve the debug override once, outside the timed phase setup.
    const bool disable_simd = is_b3_simd_disabled();
    const auto t_phase = Clock::now();
    alignas(32) std::array<cf, kDiscreteCirMaxTaps> taps_eff;
    const cf* taps = model.taps.data();
    if (model.cfo_hz != 0.0) {
        // exp(j*w*lag_lo) then a handwritten double recurrence by
        // exp(j*w): two sin/cos pairs replace tap_count cexp() calls.
        const double cw = std::cos(omega);
        const double sw = std::sin(omega);
        const double lag0 = omega * double(model.lag_lo);
        double rr = std::cos(lag0);
        double ri = std::sin(lag0);
        for (size_t d = 0; d < tap_count; ++d) {
            const double tr = model.taps[d].real();
            const double ti = model.taps[d].imag();
            taps_eff[d] = cf(static_cast<float>(tr * rr - ti * ri),
                             static_cast<float>(tr * ri + ti * rr));
            rotor_mul(rr, ri, cw, sw);
        }
        taps = taps_eff.data();
    }
    std::array<cf, kCfoBlock> rotation;
    cd block_step(1.0, 0.0);
    if (model.cfo_hz != 0.0) {
        block_step = cd(std::cos(omega * double(kCfoBlock)),
                        std::sin(omega * double(kCfoBlock)));
        build_rotation_table(omega, rotation);
    }
    r.phase_us = us(t_phase);

    float max_tap_abs = 0.f;
    for (size_t d = 0; d < tap_count; ++d)
        max_tap_abs = std::max(max_tap_abs, std::abs(model.taps[d]));

    // The tile path is safe for in-place use when every sample stays far
    // inside the float range: at most tap_count events touch a sample, so
    // |recon| <= 256 * max_tap_abs.  transaction_prepared callers own the
    // finiteness guarantee for y (the unprepared path just measured it).
    const bool tile_mode = (ihat == nullptr) &&
                           (max_tap_abs <= kTileAmplitudeLimit) &&
                           (opts.transaction_prepared ||
                            max_abs_y <= kTileAmplitudeLimit);

    // Event phase for the cold fragment path.
    auto event_phase = [&](const mod::PulseEvent& e) {
        if (model.cfo_hz == 0.0)
            return cf(static_cast<float>(e.sign), 0.f);
        const double ph = omega * static_cast<double>(e.sample_offset);
        return cf(static_cast<float>(e.sign) * static_cast<float>(std::cos(ph)),
                  static_cast<float>(e.sign) * static_cast<float>(std::sin(ph)));
    };

    double total_ey = 0.0;
    double total_eout = 0.0;

    if (ihat) {
        // --- full-buffer diagnostics: unchanged ihat semantics -------------
        const auto t_zero = Clock::now();
        std::fill(ihat, ihat + n, cf(0.f, 0.f));
        r.zero_us = us(t_zero);

        const auto t_conv = Clock::now();
        if (clipped) {
            for (size_t i = 0; i < pulses.count; ++i) {
                const mod::PulseEvent& e = pulses.events[i];
                const int64_t pos = lo + int64_t(e.sample_offset);
                if (pos >= clip_hi)
                    break; // events are strictly increasing
                if (pos + int64_t(tap_count) <= clip_lo)
                    continue;
                add_pulse_fragment(ihat, 0, pos, taps, tap_count,
                                   event_phase(e), clip_lo, clip_hi);
            }
        } else {
            accumulate_pulses_range(ihat, lo, pulses.events, pulses.count, taps,
                                    tap_count, model.cfo_hz, omega, rotation,
                                    block_step, disable_simd);
        }
        r.convolution_us = us(t_conv);
        r.cfo_us = 0.0; // CFO was folded into each pulse; no packet rotation.

        const auto t_check = Clock::now();
        for (int64_t j = clip_lo; j < clip_hi; ++j) {
            if (!finite(ihat[j]) || !finite(y[j] - ihat[j])) {
                r.check_us = us(t_check);
                r.channel_us = us(t_channel);
                return fail("nonfinite_reconstruction");
            }
        }
        r.check_us = us(t_check);
        r.channel_us = us(t_channel);

        const auto t_subtract = Clock::now();
        if (out != y) {
            if (sub_lo > 0)
                std::copy(y, y + sub_lo, out);
            if (sub_hi < int64_t(n))
                std::copy(y + sub_hi, y + n, out + sub_hi);
        }
        const size_t sub_count = static_cast<size_t>(sub_hi - sub_lo);
        if (opts.compute_roi_energy) {
            demod::rxopt::vector_sub3_energy_complex(
                out + sub_lo, y + sub_lo, ihat + sub_lo, sub_count, total_ey,
                total_eout);
            r.roi_energy_valid = true;
            r.roi_e_y = total_ey;
            r.roi_e_out = total_eout;
        } else {
            demod::rxopt::vector_sub3_complex(out + sub_lo, y + sub_lo,
                                              ihat + sub_lo, sub_count);
        }
        r.subtract_us = us(t_subtract);
    } else if (tile_mode) {
        // --- L1D tiled reconstruction + fused subtract/energy -------------
        // The buffer holds [tile_lo, tile_lo + tile_len) plus the overlap tail
        // that the events of this tile spill into the next one.  The carry
        // copy means every sample receives the same additions in the same
        // (strictly increasing event) order as the full-buffer path, so the
        // result is bit-identical.
        alignas(64) cf tile_buf[kTileSamples + kDiscreteCirMaxTaps];
        const int64_t rel_clip_lo = clip_lo - lo;

        const auto t_zero = Clock::now();
        size_t ev = 0;
        // Skip events whose complete response ends at/before clip_lo.
        while (ev < pulses.count &&
               int64_t(pulses.events[ev].sample_offset) + int64_t(tap_count) <=
                   rel_clip_lo)
            ++ev;
        {
            const size_t first_len = static_cast<size_t>(std::min<int64_t>(
                int64_t(kTileSamples), clip_hi - clip_lo));
            std::fill(tile_buf, tile_buf + first_len + tap_count - 1,
                      cf(0.f, 0.f));
        }
        r.zero_us = us(t_zero);

        const auto t_conv = Clock::now();
        // Left-clip carry: events starting before clip_lo but whose taps reach
        // into the first tile contribute only their in-window fragments.
        for (; ev < pulses.count; ++ev) {
            const int64_t off = int64_t(pulses.events[ev].sample_offset);
            if (off >= rel_clip_lo)
                break;
            add_pulse_fragment(tile_buf, clip_lo, lo + off, taps, tap_count,
                               event_phase(pulses.events[ev]), clip_lo,
                               clip_hi);
        }

        size_t prev_len = 0;
        bool first_tile = true;
        for (int64_t tile_lo = clip_lo; tile_lo < clip_hi;
             tile_lo += int64_t(kTileSamples)) {
            const size_t tile_len = static_cast<size_t>(std::min<int64_t>(
                int64_t(kTileSamples), clip_hi - tile_lo));
            if (!first_tile) {
                const size_t carry = tap_count - 1;
                if (carry)
                    std::copy(tile_buf + prev_len,
                              tile_buf + prev_len + carry, tile_buf);
                std::fill(tile_buf + carry, tile_buf + tile_len + carry,
                          cf(0.f, 0.f));
            }
            first_tile = false;

            const int64_t rel_end = (tile_lo + int64_t(tile_len)) - lo;
            const size_t ev_begin = ev;
            while (ev < pulses.count &&
                   int64_t(pulses.events[ev].sample_offset) < rel_end)
                ++ev;
            if (ev > ev_begin)
                accumulate_pulses_range(tile_buf, lo - tile_lo,
                                        pulses.events + ev_begin,
                                        ev - ev_begin, taps, tap_count,
                                        model.cfo_hz, omega, rotation,
                                        block_step, disable_simd);

            const int64_t s0 = std::max(tile_lo, sub_lo);
            const int64_t s1 = std::min(tile_lo + int64_t(tile_len), sub_hi);
            if (s1 > s0) {
                double tile_ey = 0.0, tile_eout = 0.0;
                demod::rxopt::vector_sub3_energy_complex(
                    out + s0, y + s0, tile_buf + (s0 - tile_lo),
                    static_cast<size_t>(s1 - s0), tile_ey, tile_eout);
                total_ey += tile_ey;
                total_eout += tile_eout;
            }
            prev_len = tile_len;
        }
        r.convolution_us = us(t_conv);
        r.cfo_us = 0.0;
        r.check_us = 0.0;

        // The non-finite flag is the free finite check of the fused pass.
        if (!std::isfinite(total_eout)) {
            r.channel_us = us(t_channel);
            return fail("nonfinite_reconstruction");
        }
        r.channel_us = us(t_channel);

        const auto t_subtract = Clock::now();
        if (out != y) {
            if (sub_lo > 0)
                std::copy(y, y + sub_lo, out);
            if (sub_hi < int64_t(n))
                std::copy(y + sub_hi, y + n, out + sub_hi);
        }
        r.subtract_us = us(t_subtract);
        if (opts.compute_roi_energy) {
            r.roi_energy_valid = true;
            r.roi_e_y = total_ey;
            r.roi_e_out = total_eout;
        }
    } else {
        // --- two-phase scratch.accumulator fallback -----------------------
        // Extreme amplitudes (or an aliasing/capacity corner) never touch the
        // in-place window until every sample has been validated.
        cf* const recon = scratch.accumulator.data();
        const int64_t win_len = clip_hi - clip_lo;
        const auto t_zero = Clock::now();
        std::fill(recon, recon + win_len, cf(0.f, 0.f));
        r.zero_us = us(t_zero);

        const auto t_conv = Clock::now();
        if (clipped) {
            for (size_t i = 0; i < pulses.count; ++i) {
                const mod::PulseEvent& e = pulses.events[i];
                const int64_t pos = lo + int64_t(e.sample_offset);
                if (pos >= clip_hi)
                    break;
                if (pos + int64_t(tap_count) <= clip_lo)
                    continue;
                add_pulse_fragment(recon, clip_lo, pos, taps, tap_count,
                                   event_phase(e), clip_lo, clip_hi);
            }
        } else {
            accumulate_pulses_range(recon, 0, pulses.events, pulses.count, taps,
                                    tap_count, model.cfo_hz, omega, rotation,
                                    block_step, disable_simd);
        }
        r.convolution_us = us(t_conv);
        r.cfo_us = 0.0;

        const auto t_check = Clock::now();
        for (int64_t j = clip_lo; j < clip_hi; ++j) {
            const size_t k = static_cast<size_t>(j - clip_lo);
            if (!finite(recon[k]) || !finite(y[j] - recon[k])) {
                r.check_us = us(t_check);
                r.channel_us = us(t_channel);
                return fail("nonfinite_reconstruction");
            }
        }
        r.check_us = us(t_check);
        r.channel_us = us(t_channel);

        const auto t_subtract = Clock::now();
        if (out != y) {
            if (sub_lo > 0)
                std::copy(y, y + sub_lo, out);
            if (sub_hi < int64_t(n))
                std::copy(y + sub_hi, y + n, out + sub_hi);
        }
        const cf* const sub_recon = recon + (sub_lo - clip_lo);
        const size_t sub_count = static_cast<size_t>(sub_hi - sub_lo);
        if (opts.compute_roi_energy) {
            demod::rxopt::vector_sub3_energy_complex(
                out + sub_lo, y + sub_lo, sub_recon, sub_count, total_ey,
                total_eout);
            r.roi_energy_valid = true;
            r.roi_e_y = total_ey;
            r.roi_e_out = total_eout;
        } else {
            demod::rxopt::vector_sub3_complex(out + sub_lo, y + sub_lo,
                                              sub_recon, sub_count);
        }
        r.subtract_us = us(t_subtract);
    }

    r.packet_samples = pulses.packet_samples;
    r.committed = true;
    r.reason = "committed";
    r.apply_total_us = us(t_apply);
    return true;
}

// Public standalone entry retains its own bypass preparation and IQ scan.
bool apply_discrete_cir_pulses(
    const std::complex<float>* y, size_t n, const mod::PulseList& pulses,
    const std::complex<float>* reference_sync, size_t sync_n,
    const DiscreteCirModel& model, int64_t roi_lo, int64_t roi_hi,
    DiscreteCirScratch& scratch, std::complex<float>* out,
    std::complex<float>* ihat, DiscreteCirResult* result)
{
    const DiscreteCirApplyOptions opts;
    return apply_discrete_cir_pulses_ex(y, n, pulses, reference_sync, sync_n,
        model, roi_lo, roi_hi, opts, scratch, out, ihat, result);
}

bool apply_discrete_cir_pulses_ex(
    const std::complex<float>* y, size_t n, const mod::PulseList& pulses,
    const std::complex<float>* reference_sync, size_t sync_n,
    const DiscreteCirModel& model, int64_t roi_lo, int64_t roi_hi,
    const DiscreteCirApplyOptions& opts, DiscreteCirScratch& scratch,
    std::complex<float>* out, std::complex<float>* ihat,
    DiscreteCirResult* result)
{
    return apply_discrete_cir_pulses_impl(y, n, pulses, reference_sync, sync_n,
        model, roi_lo, roi_hi, opts, scratch, out, ihat, result);
}

DiscreteCirResult fullpacket_cancel_discrete_ex(
    const std::complex<float>* y,
    size_t n,
    const std::vector<uint8_t>& phr_bits,
    const std::vector<uint8_t>& wire_psdu,
    bool fcs_pass,
    const mod::HrpModConfig& profile,
    const DiscreteCirModel& model,
    const std::complex<float>* reference_sync,
    size_t sync_n,
    int64_t roi_lo,
    int64_t roi_hi,
    const DiscreteCirApplyOptions& opts,
    DiscreteCirScratch& scratch,
    std::complex<float>* out,
    std::complex<float>* ihat)
{
    DiscreteCirResult r;
    const auto t_validate = Clock::now();
    auto fail = [&](const char* reason) {
        if (opts.copy_out_on_fail && out && y && n > 0 && out != y)
            std::copy(y, y + n, out);
        if (ihat && ihat != y && ihat != out)
            std::fill(ihat, ihat + n, cf(0.f, 0.f));
        r.committed = false;
        r.reason = reason;
        r.validate_us = us(t_validate);
        return r;
    };

    if (!y || !out || n == 0 || (ihat && (ihat == y || ihat == out)))
        return fail("invalid_buffers");

    const char* sf = transaction_scratch_failure(scratch, n);
    if (sf)
        return fail(sf);
    if (!reference_sync || sync_n != scratch.sps)
        return fail("invalid_template");
    if (!fcs_pass)
        return fail("fcs_failed");
    // The Legacy-pulse gate only protects the shaped-waveform path, which
    // needs profile.pulse to remodulate.  The Ipatov path uses the private
    // unit-impulse profile and must not reject a non-Legacy caller profile.
    if (profile.insert_sts || profile.ranging || !profile.sfd_mode ||
        !radar_meta::code_index_supported(profile.code_index) ||
        !radar_meta::sync_reps_supported(profile.sync_repetitions) ||
        demod::GetSfdSequence(profile.sfd_mode).empty() ||
        (!model.impulse_template &&
         profile.pulse.shape != mod::PulseShape::Legacy))
        return fail("unsupported_profile");
    // Historical contract for direct callers: only code10 / 64 SYNC / 4z2 is
    // exercised by the frozen full-packet core.  Unified SIC opts in through
    // prepare_discrete_cir(allow_configured_profile) and may use its own
    // validated public profile.
    if (!scratch.allow_configured_profile &&
        (profile.code_index != 10 || profile.sync_repetitions != 64 ||
         std::strcmp(profile.sfd_mode, "4z2") != 0))
        return fail("unsupported_profile");
    if (phr_bits.size() != 19 || wire_psdu.size() < 2 ||
        wire_psdu.size() > radar_meta::kMaxPsduBytes)
        return fail("invalid_frame_length");
    const size_t crc_len = wire_psdu.size() - 2;
    const uint16_t crc_expected =
        static_cast<uint16_t>(wire_psdu[crc_len]) |
        static_cast<uint16_t>(static_cast<uint16_t>(wire_psdu[crc_len + 1]) << 8);
    if (codec::crc16_802154(wire_psdu.data(), crc_len) != crc_expected)
        return fail("fcs_failed");
    mod::detail::encode_phr19(wire_psdu.size(),
                              profile.sync_repetitions,
                              profile.ranging,
                              scratch.expected_phr);
    for (size_t i = 0; i < 19; ++i)
        if (phr_bits[i] > 1 ||
            phr_bits[i] != static_cast<uint8_t>(scratch.expected_phr[i]))
            return fail("unsupported_phr");

    const char* mf = model_failure(model, scratch.sps, n);
    if (mf)
        return fail(mf);
    if (roi_lo < 0 || roi_hi > int64_t(n) || roi_hi <= roi_lo)
        return fail("invalid_roi");
    // Sparse SIC experiment trusts the receiver IQ: do not scan the entire
    // window for NaN/Inf. Keep frame/model/bounds validation and the standalone
    // public apply input check. Legacy waveform transactions retain the scan
    // unless the caller already prepared the transaction.
    if (!opts.transaction_prepared &&
        !(model.impulse_template && scratch.sparse_pulse_reconstruct))
        for (size_t i = 0; i < n; ++i)
            if (!finite(y[i]))
                return fail("nonfinite_input");
    r.validate_us = us(t_validate);

    // Remodulate the decoded wire PSDU once into the scratch packet buffer.
    // The capacity is the REAL buffer length, not the estimator's planning
    // bound, so a stale max_packet_samples cannot reject a fitting packet.
    const auto t_tx = Clock::now();
    if (model.impulse_template != scratch.ipatov_correlation) {
        return fail("template_kind_mismatch");
    }
    const bool real_only = model.impulse_template;
    if (real_only && scratch.sparse_pulse_reconstruct) {
        // Sparse experiment: the decoded packet is generated as positions and
        // signs only.  No real/complex waveform is materialized.
        mod::PulseList pulses;
        const bool tx_ok = mod::modulate_pulse_list(
            wire_psdu.data(), wire_psdu.size(), scratch.impulse_profile,
            scratch.modulator, scratch.waveform.size(), pulses);
        if (!tx_ok || pulses.count == 0) {
            r.tx_us = us(t_tx);
            return fail("tx_failed");
        }
        const double tx_us = us(t_tx);

        DiscreteCirApplyOptions prepared = opts;
        // The list was just generated from the validated wire PSDU.
        prepared.transaction_prepared = true;
        DiscreteCirResult applied;
        apply_discrete_cir_pulses_impl(y, n, pulses, reference_sync, sync_n, model,
                                  roi_lo, roi_hi, prepared, scratch, out, ihat,
                                  &applied);
        applied.validate_us = r.validate_us;
        applied.tx_us = tx_us;
        return applied;
    }

    size_t packet_n = 0;
    const bool tx_ok = real_only
        ? mod::modulate_unit_impulses(wire_psdu.data(), wire_psdu.size(),
                                    scratch.impulse_profile, scratch.modulator,
                                    scratch.waveform.size(), packet_n)
        : mod::modulate_one(wire_psdu.data(), wire_psdu.size(), profile,
                            scratch.modulator, scratch.waveform.data(),
                            scratch.waveform.size(), packet_n);
    if (!tx_ok || packet_n == 0 || packet_n > scratch.waveform.size()) {
        r.tx_us = us(t_tx);
        return fail("tx_failed");
    }
    const double tx_us = us(t_tx);

    DiscreteCirResult applied;
    apply_discrete_cir_impl(y,
                       n,
                       real_only ? nullptr : scratch.waveform.data(),
                       packet_n,
                       reference_sync,
                       sync_n,
                       model,
                       roi_lo,
                       roi_hi,
                       opts,
                       scratch,
                       out,
                       ihat,
                       &applied,
                       real_only ? scratch.modulator.iq_re.data() : nullptr);
    // apply() owns channel_us/subtract_us/apply_total_us; the wrapper owns the
    // validation and remodulation stages (never double-counted as children).
    applied.validate_us = r.validate_us;
    applied.tx_us = tx_us;
    return applied;
}

DiscreteCirResult fullpacket_cancel_discrete(
    const std::complex<float>* y,
    size_t n,
    const std::vector<uint8_t>& phr_bits,
    const std::vector<uint8_t>& wire_psdu,
    bool fcs_pass,
    const mod::HrpModConfig& profile,
    const DiscreteCirModel& model,
    const std::complex<float>* reference_sync,
    size_t sync_n,
    int64_t roi_lo,
    int64_t roi_hi,
    DiscreteCirScratch& scratch,
    std::complex<float>* out,
    std::complex<float>* ihat)
{
    const DiscreteCirApplyOptions opts;
    return fullpacket_cancel_discrete_ex(y, n, phr_bits, wire_psdu, fcs_pass,
                                         profile, model, reference_sync, sync_n,
                                         roi_lo, roi_hi, opts, scratch, out,
                                         ihat);
}

} // namespace sic
} // namespace openuwb
