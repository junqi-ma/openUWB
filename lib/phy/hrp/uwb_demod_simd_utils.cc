// SPDX-License-Identifier: GPL-3.0-or-later
//
// rxopt SIMD vector operators, rotation ramps, sparse/shaped template builders
// and CIR correlation primitives.  These used to live at the tail of
// uwb_demod_timing.cc but are shared by the whole RX pipeline, so they get
// their own translation unit.  ISA dispatch (AVX2 / AVX-512 / scalar) and the
// compile-time ISA selection are unchanged bit-for-bit.
#include <openuwb/phy/hrp/rx/demod_timing.h>
#include <openuwb/phy/hrp/rx/demod_sync.h>
#include <openuwb/phy/hrp/rx/rx_matched_filter.h>
#include <openuwb/phy/hrp/rx/rx_sync_opt.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace openuwb::demod::rxopt {

void build_real_sync_template(const std::vector<std::complex<float>>& template_wf,
                              RealSyncTemplate& out)
{
    out.reset();
    if (template_wf.empty())
        return;
    const size_t L = template_wf.size();
    bool real_only = true;
    for (size_t k = 0; k < L; ++k) {
        if (std::abs(template_wf[k].imag()) > 1e-7f) {
            real_only = false;
            break;
        }
    }
    out.is_real = real_only;
    out.length = L;
    if (real_only) {
        out.taps_dup.resize(2 * L);
        float e = 0.0f;
        for (size_t k = 0; k < L; ++k) {
            const float r = template_wf[k].real();
            out.taps_dup[2 * k] = r;
            out.taps_dup[2 * k + 1] = r;
            e += r * r;
        }
        out.energy = e;
    }
}

uint64_t
sfd_template_key(const std::vector<std::complex<float>>& template_wf,
                 const std::vector<int8_t>& sfd_sequence)
{
    static_assert(sizeof(float) == 4, "SFD cache key assumes 32-bit float");
    uint64_t h = detail::kFnvOffsetBasis;
    h = detail::fnv1a_u64(h, static_cast<uint64_t>(template_wf.size()));
    for (const std::complex<float>& v : template_wf) {
        uint32_t re_bits = 0;
        uint32_t im_bits = 0;
        const float re = v.real();
        const float im = v.imag();
        std::memcpy(&re_bits, &re, sizeof(re_bits));
        std::memcpy(&im_bits, &im, sizeof(im_bits));
        h = detail::fnv1a_u32(h, re_bits);
        h = detail::fnv1a_u32(h, im_bits);
    }
    h = detail::fnv1a_u64(h, static_cast<uint64_t>(sfd_sequence.size()));
    if (!sfd_sequence.empty())
        h = detail::fnv1a_bytes(h, sfd_sequence.data(), sfd_sequence.size());
    return h;
}

void
build_sfd_template(const std::vector<std::complex<float>>& template_wf,
                   const std::vector<int8_t>& sfd_sequence,
                   std::vector<std::complex<float>>& out)
{
    out.clear();
    const size_t sym = kQm35SamplesPerSymbol;
    if (sfd_sequence.empty() || template_wf.size() < sym)
        return; // the legacy loop would read past the template: fail closed
    if (sfd_sequence.size() > std::numeric_limits<size_t>::max() / sym)
        return;
    out.resize(sfd_sequence.size() * sym);
    for (size_t i = 0; i < sfd_sequence.size(); ++i) {
        const int8_t s = sfd_sequence[i];
        for (size_t k = 0; k < sym; ++k)
            out[i * sym + k] =
                std::complex<float>(static_cast<float>(s) * template_wf[k].real(),
                                    static_cast<float>(s) * template_wf[k].imag());
    }
    openuwb::core::uwb_l2_normalize(out);
}

namespace {

#if UWB_RXOPT_HAVE_AVX2
// Local AVX2 three-operand add.  Sixteen complex samples per iteration
// (4 x 256-bit) plus one 256-bit register and a scalar tail.  Both sources are
// loaded before dst is stored, so dst == a or dst == b stays safe.
__attribute__((target("avx2")))
void vector_add3_complex_avx2(std::complex<float>* dst,
                              const std::complex<float>* a,
                              const std::complex<float>* b, size_t count)
{
    float* d = reinterpret_cast<float*>(dst);
    const float* pa = reinterpret_cast<const float*>(a);
    const float* pb = reinterpret_cast<const float*>(b);
    size_t i = 0;
    const size_t limit16 = count & ~size_t(15);
    for (; i < limit16; i += 16) {
        const __m256 a0 = _mm256_loadu_ps(pa + 2 * i);
        const __m256 a1 = _mm256_loadu_ps(pa + 2 * i + 8);
        const __m256 a2 = _mm256_loadu_ps(pa + 2 * i + 16);
        const __m256 a3 = _mm256_loadu_ps(pa + 2 * i + 24);
        const __m256 b0 = _mm256_loadu_ps(pb + 2 * i);
        const __m256 b1 = _mm256_loadu_ps(pb + 2 * i + 8);
        const __m256 b2 = _mm256_loadu_ps(pb + 2 * i + 16);
        const __m256 b3 = _mm256_loadu_ps(pb + 2 * i + 24);
        _mm256_storeu_ps(d + 2 * i, _mm256_add_ps(a0, b0));
        _mm256_storeu_ps(d + 2 * i + 8, _mm256_add_ps(a1, b1));
        _mm256_storeu_ps(d + 2 * i + 16, _mm256_add_ps(a2, b2));
        _mm256_storeu_ps(d + 2 * i + 24, _mm256_add_ps(a3, b3));
    }
    for (; i + 4 <= count; i += 4) {
        const __m256 va = _mm256_loadu_ps(pa + 2 * i);
        const __m256 vb = _mm256_loadu_ps(pb + 2 * i);
        _mm256_storeu_ps(d + 2 * i, _mm256_add_ps(va, vb));
    }
    for (; i < count; ++i)
        dst[i] = a[i] + b[i];
}

__attribute__((target("avx2")))
void vector_accumulate_complex_dual_avx2(std::complex<float>* dst0,
                                         std::complex<float>* dst1,
                                         const std::complex<float>* src,
                                         size_t count)
{
    float* d0 = reinterpret_cast<float*>(dst0);
    float* d1 = reinterpret_cast<float*>(dst1);
    const float* s = reinterpret_cast<const float*>(src);
    const size_t num_floats = 2 * count;
    size_t i = 0;
    const size_t limit = num_floats & ~size_t(7);
    for (; i < limit; i += 8) {
        const __m256 vs = _mm256_loadu_ps(s + i);
        const __m256 vd0 = _mm256_loadu_ps(d0 + i);
        const __m256 vd1 = _mm256_loadu_ps(d1 + i);
        _mm256_storeu_ps(d0 + i, _mm256_add_ps(vd0, vs));
        _mm256_storeu_ps(d1 + i, _mm256_add_ps(vd1, vs));
    }
    for (; i < num_floats; ++i) {
        const float sv = s[i];
        d0[i] += sv;
        d1[i] += sv;
    }
}

__attribute__((target("avx2")))
void vector_accumulate_scaled_complex_dual_avx2(std::complex<float>* dst0,
                                                std::complex<float>* dst1,
                                                const std::complex<float>* src,
                                                float c, float s, size_t count)
{
    const __m256 vc = _mm256_set1_ps(c);
    const __m256 vs = _mm256_set1_ps(s);
    size_t j = 0;
    const size_t limit8 = count & ~size_t(7);
    for (; j < limit8; j += 8) {
        const __m256 v0 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(src + j));
        const __m256 v1 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(src + j + 4));
        const __m256 rot0 = _mm256_addsub_ps(
            _mm256_mul_ps(v0, vc),
            _mm256_mul_ps(_mm256_permute_ps(v0, 0xB1), vs));
        const __m256 rot1 = _mm256_addsub_ps(
            _mm256_mul_ps(v1, vc),
            _mm256_mul_ps(_mm256_permute_ps(v1, 0xB1), vs));
        const __m256 d00 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(dst0 + j));
        const __m256 d01 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(dst0 + j + 4));
        const __m256 d10 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(dst1 + j));
        const __m256 d11 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(dst1 + j + 4));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst0 + j),
                         _mm256_add_ps(d00, rot0));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst0 + j + 4),
                         _mm256_add_ps(d01, rot1));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst1 + j),
                         _mm256_add_ps(d10, rot0));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst1 + j + 4),
                         _mm256_add_ps(d11, rot1));
    }
    for (; j + 4 <= count; j += 4) {
        const __m256 v =
            _mm256_loadu_ps(reinterpret_cast<const float*>(src + j));
        const __m256 rot = _mm256_addsub_ps(
            _mm256_mul_ps(v, vc),
            _mm256_mul_ps(_mm256_permute_ps(v, 0xB1), vs));
        const __m256 d0 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(dst0 + j));
        const __m256 d1 =
            _mm256_loadu_ps(reinterpret_cast<const float*>(dst1 + j));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst0 + j),
                         _mm256_add_ps(d0, rot));
        _mm256_storeu_ps(reinterpret_cast<float*>(dst1 + j),
                         _mm256_add_ps(d1, rot));
    }
    for (; j < count; ++j) {
        const auto v = src[j];
        const std::complex<float> rot(c * v.real() - s * v.imag(),
                                      s * v.real() + c * v.imag());
        dst0[j] += rot;
        dst1[j] += rot;
    }
}

// Local AVX2 three-operand subtract.  Sixteen complex samples per iteration
// (4 x 256-bit) plus one 256-bit register and a scalar tail.  Both sources are
// loaded before dst is stored, so dst == a stays safe.
__attribute__((target("avx2")))
void vector_sub3_complex_avx2(std::complex<float>* dst,
                              const std::complex<float>* a,
                              const std::complex<float>* b, size_t count)
{
    float* d = reinterpret_cast<float*>(dst);
    const float* pa = reinterpret_cast<const float*>(a);
    const float* pb = reinterpret_cast<const float*>(b);
    size_t i = 0;
    const size_t limit16 = count & ~size_t(15);
    for (; i < limit16; i += 16) {
        const __m256 a0 = _mm256_loadu_ps(pa + 2 * i);
        const __m256 a1 = _mm256_loadu_ps(pa + 2 * i + 8);
        const __m256 a2 = _mm256_loadu_ps(pa + 2 * i + 16);
        const __m256 a3 = _mm256_loadu_ps(pa + 2 * i + 24);
        const __m256 b0 = _mm256_loadu_ps(pb + 2 * i);
        const __m256 b1 = _mm256_loadu_ps(pb + 2 * i + 8);
        const __m256 b2 = _mm256_loadu_ps(pb + 2 * i + 16);
        const __m256 b3 = _mm256_loadu_ps(pb + 2 * i + 24);
        _mm256_storeu_ps(d + 2 * i, _mm256_sub_ps(a0, b0));
        _mm256_storeu_ps(d + 2 * i + 8, _mm256_sub_ps(a1, b1));
        _mm256_storeu_ps(d + 2 * i + 16, _mm256_sub_ps(a2, b2));
        _mm256_storeu_ps(d + 2 * i + 24, _mm256_sub_ps(a3, b3));
    }
    for (; i + 4 <= count; i += 4) {
        const __m256 va = _mm256_loadu_ps(pa + 2 * i);
        const __m256 vb = _mm256_loadu_ps(pb + 2 * i);
        _mm256_storeu_ps(d + 2 * i, _mm256_sub_ps(va, vb));
    }
    for (; i < count; ++i)
        dst[i] = a[i] - b[i];
}

// Fused AVX2 dst = a - b + double energies.  Eight complex samples per
// iteration (2 x 256-bit loads = four 128-bit halves per energy), four
// independent __m256d accumulators per energy -- the same structure as
// compute_roi_energies() in lib/sic/uwb_sic_unified.cc.  Squares are formed in
// float and widened to double, matching that reference and the AVX-512 kernel.
__attribute__((target("avx2")))
void vector_sub3_energy_complex_avx2(std::complex<float>* dst,
                                     const std::complex<float>* a,
                                     const std::complex<float>* b,
                                     size_t count, double& e_a,
                                     double& e_dst)
{
    float* d = reinterpret_cast<float*>(dst);
    const float* pa = reinterpret_cast<const float*>(a);
    const float* pb = reinterpret_cast<const float*>(b);

    __m256d ea0 = _mm256_setzero_pd();
    __m256d ea1 = _mm256_setzero_pd();
    __m256d ea2 = _mm256_setzero_pd();
    __m256d ea3 = _mm256_setzero_pd();
    __m256d ed0 = _mm256_setzero_pd();
    __m256d ed1 = _mm256_setzero_pd();
    __m256d ed2 = _mm256_setzero_pd();
    __m256d ed3 = _mm256_setzero_pd();

    size_t i = 0;
    const size_t limit8 = count & ~size_t(7);
    for (; i < limit8; i += 8) {
        // Load a before overwriting dst, so dst == a is safe.
        const __m256 va0 = _mm256_loadu_ps(pa + 2 * i);
        const __m256 va1 = _mm256_loadu_ps(pa + 2 * i + 8);
        const __m256 vb0 = _mm256_loadu_ps(pb + 2 * i);
        const __m256 vb1 = _mm256_loadu_ps(pb + 2 * i + 8);
        const __m256 vd0 = _mm256_sub_ps(va0, vb0);
        const __m256 vd1 = _mm256_sub_ps(va1, vb1);
        _mm256_storeu_ps(d + 2 * i, vd0);
        _mm256_storeu_ps(d + 2 * i + 8, vd1);

        const __m256 sa0 = _mm256_mul_ps(va0, va0);
        const __m256 sa1 = _mm256_mul_ps(va1, va1);
        const __m256 sd0 = _mm256_mul_ps(vd0, vd0);
        const __m256 sd1 = _mm256_mul_ps(vd1, vd1);

        ea0 = _mm256_add_pd(
            ea0, _mm256_cvtps_pd(_mm256_castps256_ps128(sa0)));
        ea1 = _mm256_add_pd(
            ea1, _mm256_cvtps_pd(_mm256_extractf128_ps(sa0, 1)));
        ea2 = _mm256_add_pd(
            ea2, _mm256_cvtps_pd(_mm256_castps256_ps128(sa1)));
        ea3 = _mm256_add_pd(
            ea3, _mm256_cvtps_pd(_mm256_extractf128_ps(sa1, 1)));

        ed0 = _mm256_add_pd(
            ed0, _mm256_cvtps_pd(_mm256_castps256_ps128(sd0)));
        ed1 = _mm256_add_pd(
            ed1, _mm256_cvtps_pd(_mm256_extractf128_ps(sd0, 1)));
        ed2 = _mm256_add_pd(
            ed2, _mm256_cvtps_pd(_mm256_castps256_ps128(sd1)));
        ed3 = _mm256_add_pd(
            ed3, _mm256_cvtps_pd(_mm256_extractf128_ps(sd1, 1)));
    }

    alignas(32) double buf[4];
    _mm256_storeu_pd(buf, _mm256_add_pd(_mm256_add_pd(ea0, ea1),
                                        _mm256_add_pd(ea2, ea3)));
    double sum_a = buf[0] + buf[1] + buf[2] + buf[3];
    _mm256_storeu_pd(buf, _mm256_add_pd(_mm256_add_pd(ed0, ed1),
                                        _mm256_add_pd(ed2, ed3)));
    double sum_dst = buf[0] + buf[1] + buf[2] + buf[3];

    for (; i < count; ++i) {
        const float ar = pa[2 * i];
        const float ai = pa[2 * i + 1];
        const float dr = ar - pb[2 * i];
        const float di = ai - pb[2 * i + 1];
        d[2 * i] = dr;
        d[2 * i + 1] = di;
        sum_a += static_cast<double>(ar * ar) +
                 static_cast<double>(ai * ai);
        sum_dst += static_cast<double>(dr * dr) +
                   static_cast<double>(di * di);
    }

    e_a = sum_a;
    e_dst = sum_dst;
}
#endif

} // namespace

void vector_accumulate_complex(std::complex<float>* dst,
                               const std::complex<float>* src,
                               size_t count)
{
    if (count == 0)
        return;
#if UWB_SIMD_HAVE_AVX512
    if (use_avx512(openuwb::core::simd::KernelGroup::Sync)) {
        openuwb::core::simd::avx512::accumulate_cf32(dst, src, count);
        return;
    }
#endif
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_accumulate_complex_avx2(dst, src, count);
        return;
    }
#endif
    for (size_t m = 0; m < count; ++m)
        dst[m] += src[m];
}

void vector_accumulate_complex_dual(std::complex<float>* dst0,
                                    std::complex<float>* dst1,
                                    const std::complex<float>* src,
                                    size_t count)
{
    if (count == 0)
        return;
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_accumulate_complex_dual_avx2(dst0, dst1, src, count);
        return;
    }
#endif
    for (size_t m = 0; m < count; ++m) {
        const std::complex<float> sv = src[m];
        dst0[m] += sv;
        dst1[m] += sv;
    }
}

void vector_add3_complex(std::complex<float>* dst,
                         const std::complex<float>* a,
                         const std::complex<float>* b, size_t count)
{
    if (count == 0)
        return;
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_add3_complex_avx2(dst, a, b, count);
        return;
    }
#endif
    for (size_t m = 0; m < count; ++m)
        dst[m] = a[m] + b[m];
}

void vector_sub_complex(std::complex<float>* dst,
                        const std::complex<float>* src,
                        size_t count)
{
    if (count == 0)
        return;
#if UWB_SIMD_HAVE_AVX512
    if (use_avx512(openuwb::core::simd::KernelGroup::Sync)) {
        openuwb::core::simd::avx512::sub_cf32(dst, src, count);
        return;
    }
#endif
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_sub_complex_avx2(dst, src, count);
        return;
    }
#endif
    for (size_t m = 0; m < count; ++m)
        dst[m] -= src[m];
}

void vector_sub3_complex(std::complex<float>* dst,
                         const std::complex<float>* a,
                         const std::complex<float>* b,
                         size_t count)
{
    if (count == 0)
        return;
#if UWB_SIMD_HAVE_AVX512
    if (use_avx512(openuwb::core::simd::KernelGroup::Sync)) {
        openuwb::core::simd::avx512::sub3_cf32(dst, a, b, count);
        return;
    }
#endif
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_sub3_complex_avx2(dst, a, b, count);
        return;
    }
#endif
    for (size_t m = 0; m < count; ++m)
        dst[m] = a[m] - b[m];
}

void vector_sub3_energy_complex(std::complex<float>* dst,
                                const std::complex<float>* a,
                                const std::complex<float>* b, size_t count,
                                double& e_a, double& e_dst)
{
    // Overwrite (never accumulate into) the caller's energies on entry.
    e_a = 0.0;
    e_dst = 0.0;
    if (count == 0)
        return;
#if UWB_SIMD_HAVE_AVX512
    if (use_avx512(openuwb::core::simd::KernelGroup::Sync)) {
        openuwb::core::simd::avx512::sub3_energy_cf32(dst, a, b, count, &e_a,
                                                      &e_dst);
        return;
    }
#endif
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_sub3_energy_complex_avx2(dst, a, b, count, e_a, e_dst);
        return;
    }
#endif
    // Scalar double fallback: float squares widened to double, NaN/Inf
    // propagate through the IEEE default so isfinite(e_dst) still fails.
    double sum_a = 0.0;
    double sum_dst = 0.0;
    for (size_t m = 0; m < count; ++m) {
        const float ar = a[m].real();
        const float ai = a[m].imag();
        const float dr = ar - b[m].real();
        const float di = ai - b[m].imag();
        dst[m] = std::complex<float>(dr, di);
        sum_a += static_cast<double>(ar * ar) +
                 static_cast<double>(ai * ai);
        sum_dst += static_cast<double>(dr * dr) +
                   static_cast<double>(di * di);
    }
    e_a = sum_a;
    e_dst = sum_dst;
}

void vector_accumulate_scaled_complex(std::complex<float>* dst,
                                      const std::complex<float>* src,
                                      std::complex<float> alpha,
                                      size_t count)
{
    if (count == 0)
        return;
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_accumulate_scaled_complex_avx2(dst, src, alpha.real(),
                                              alpha.imag(), count);
        return;
    }
#endif
    for (size_t j = 0; j < count; ++j)
        dst[j] += alpha * src[j];
}

void vector_accumulate_scaled_complex_dual(std::complex<float>* dst0,
                                           std::complex<float>* dst1,
                                           const std::complex<float>* src,
                                           std::complex<float> alpha,
                                           size_t count)
{
    if (count == 0)
        return;
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_accumulate_scaled_complex_dual_avx2(
            dst0, dst1, src, alpha.real(), alpha.imag(), count);
        return;
    }
#endif
    for (size_t j = 0; j < count; ++j) {
        const std::complex<float> rot = alpha * src[j];
        dst0[j] += rot;
        dst1[j] += rot;
    }
}

void vector_mul_complex(std::complex<float>* dst,
                        const std::complex<float>* mul,
                        size_t count)
{
    if (count == 0)
        return;
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_mul_complex_avx2(dst, mul, count);
        return;
    }
#endif
    for (size_t j = 0; j < count; ++j)
        dst[j] *= mul[j];
}

void vector_scale_complex(std::complex<float>* dst,
                          float scale,
                          size_t count)
{
    if (count == 0)
        return;
#if UWB_SIMD_HAVE_AVX512
    if (use_avx512(openuwb::core::simd::KernelGroup::Sync)) {
        openuwb::core::simd::avx512::scale_cf32(dst, scale, count);
        return;
    }
#endif
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        vector_scale_complex_avx2(dst, scale, count);
        return;
    }
#endif
    for (size_t m = 0; m < count; ++m)
        dst[m] *= scale;
}

void fold_sfd_window(const std::complex<float>* rx_base,
                     size_t sym,
                     const int8_t* sfd_seq,
                     size_t sfd_len,
                     size_t fold_len,
                     std::complex<float>* out)
{
    if (!rx_base || !sfd_seq || !out || sym == 0 || sfd_len == 0 || fold_len == 0)
        return;
#if UWB_SIMD_HAVE_AVX512
    if (use_avx512(openuwb::core::simd::KernelGroup::Sync)) {
        // Batched: the whole symbol loop runs inside the wide translation
        // unit, so the fold pays one call instead of one per symbol.
        openuwb::core::simd::avx512::fold_sfd_window(rx_base, sym, sfd_seq,
                                                     sfd_len, fold_len, out);
        return;
    }
#endif
#if UWB_RXOPT_HAVE_AVX2
    if (has_avx2_fma()) {
        fold_sfd_window_avx2(rx_base, sym, sfd_seq, sfd_len, fold_len, out);
        return;
    }
#endif
    size_t idx[8];
    int8_t sg[8];
    bool have_out = false;
    for (size_t chunk = 0; chunk < sfd_len; chunk += 8) {
        const size_t chunk_end = std::min(chunk + 8, sfd_len);
        size_t cnt = 0;
        for (size_t i = chunk; i < chunk_end; ++i) {
            if (sfd_seq[i] == 0)
                continue;
            idx[cnt] = i;
            sg[cnt] = sfd_seq[i];
            ++cnt;
        }
        if (cnt == 0)
            continue;
        const bool first_chunk = !have_out;
        for (size_t k = 0; k < fold_len; ++k) {
            size_t j = 0;
            std::complex<float> acc;
            if (first_chunk) {
                const std::complex<float> v = rx_base[idx[0] * sym + k];
                acc = (sg[0] == 1) ? v : std::complex<float>(0.0f, 0.0f) - v;
                j = 1;
            } else {
                acc = out[k];
            }
            for (; j < cnt; ++j) {
                const std::complex<float> v = rx_base[idx[j] * sym + k];
                acc = (sg[j] == 1) ? acc + v : acc - v;
            }
            out[k] = acc;
        }
        have_out = true;
    }
    if (!have_out)
        std::fill(out, out + fold_len, std::complex<float>(0.0f, 0.0f));
}

void build_rotation_ramp(double w, size_t count,
                         std::complex<float>* table)
{
    if (count == 0)
        return;
    constexpr size_t kB = 1024;
    const float s4r = static_cast<float>(std::cos(-4.0 * w));
    const float s4i = static_cast<float>(std::sin(-4.0 * w));
    for (size_t b0 = 0; b0 < count; b0 += kB) {
        const size_t b1 = std::min(count, b0 + kB);
        float cr[4], ci[4];
        for (int c = 0; c < 4; ++c) {
            const double ph = -w * static_cast<double>(b0 + c);
            cr[c] = static_cast<float>(std::cos(ph));
            ci[c] = static_cast<float>(std::sin(ph));
        }
        size_t i = b0;
        for (; i + 4 <= b1; i += 4) {
            for (int c = 0; c < 4; ++c) {
                table[i + c] = std::complex<float>(cr[c], ci[c]);
                const float nr = cr[c] * s4r - ci[c] * s4i;
                const float ni = cr[c] * s4i + ci[c] * s4r;
                cr[c] = nr;
                ci[c] = ni;
            }
        }
        for (; i < b1; ++i) {
            const double ph = -w * static_cast<double>(i);
            table[i] = std::complex<float>(static_cast<float>(std::cos(ph)),
                                           static_cast<float>(std::sin(ph)));
        }
    }
}

void build_sparse_cir_template(const int8_t* code,
                               size_t code_len,
                               SparseCirTemplate& tmpl)
{
    tmpl.reset();
    tmpl.key = cir_template_key(code, code_len);
    const std::vector<int8_t> spread = BuildSampledCode(code, code_len);
    tmpl.entries.reserve(code_len);
    for (size_t j = 0; j < spread.size(); ++j) {
        if (spread[j] != 0) {
            const size_t pos = 2 * j;
            tmpl.entries.push_back({ pos, spread[j] });
            if (spread[j] > 0) {
                tmpl.pos_offsets.push_back(2 * pos);
                if (tmpl.pos_count < kMaxSparseCirOffsets)
                    tmpl.pos_off[tmpl.pos_count++] =
                        static_cast<uint16_t>(2 * pos);
            } else {
                tmpl.neg_offsets.push_back(2 * pos);
                if (tmpl.neg_count < kMaxSparseCirOffsets)
                    tmpl.neg_off[tmpl.neg_count++] =
                        static_cast<uint16_t>(2 * pos);
            }
        }
    }
    tmpl.code_energy = static_cast<float>(tmpl.entries.size());
}

void build_half_rate_sparse_template(const SparseCirTemplate& full,
                                     SparseCirTemplate& half)
{
    half.reset();
    half.key = full.key ^ 0x9e3779b97f4a7c15ULL;
    half.code_energy = full.code_energy;
    half.entries.reserve(full.entries.size());
    for (const auto& e : full.entries)
        half.entries.push_back({ e.pos / 2, e.sign });
    half.pos_offsets.reserve(full.pos_offsets.size());
    for (const size_t off : full.pos_offsets)
        half.pos_offsets.push_back(off / 2);
    half.neg_offsets.reserve(full.neg_offsets.size());
    for (const size_t off : full.neg_offsets)
        half.neg_offsets.push_back(off / 2);
    half.pos_count = full.pos_count;
    for (size_t i = 0; i < full.pos_count; ++i)
        half.pos_off[i] = static_cast<uint16_t>(full.pos_off[i] / 2);
    half.neg_count = full.neg_count;
    for (size_t i = 0; i < full.neg_count; ++i)
        half.neg_off[i] = static_cast<uint16_t>(full.neg_off[i] / 2);
}

bool build_shaped_sparse_template(const std::complex<float>* tmpl,
                                  size_t len,
                                  const SparseCirTemplate& code,
                                  ShapedSparseTemplate& out)
{
    out.reset();
    if (!tmpl || len < 32 || code.entries.empty() ||
        code.entries.size() > kMaxSparseCirOffsets)
        return false;
    const SparseCirEntry& first = code.entries.front();
    const size_t base = first.pos;
    const size_t last = code.entries.back().pos;
    if (first.sign == 0 || base >= len || last >= len)
        return false;
    const size_t nt = std::min(ShapedSparseTemplate::kMaxPulse, len - last);
    if (nt < 3)
        return false;

    double energy = 0.0;
    double imag_energy = 0.0;
    for (size_t k = 0; k < len; ++k) {
        energy += static_cast<double>(std::norm(tmpl[k]));
        imag_energy += static_cast<double>(tmpl[k].imag()) *
                       static_cast<double>(tmpl[k].imag());
    }
    if (energy <= 0.0 || imag_energy > 1e-12 * energy)
        return false;

    // First pulse is uncontaminated up to the next chip (8 samples).
    const float sign0 = static_cast<float>(first.sign);
    for (size_t k = 0; k < nt; ++k)
        out.pulse[k] = tmpl[base + k].real() / sign0;

    // Reconstruction gate. The residual is the omitted pulse tail, which is
    // ~0.3% of the energy for kPulse48; anything above 5% is not code*p.
    double resid = 0.0;
    for (const auto& e : code.entries) {
        for (size_t k = 0; k < nt; ++k) {
            const size_t idx = e.pos + k;
            if (idx >= len)
                break;
            const double rec =
                static_cast<double>(e.sign) * out.pulse[k];
            const double d = rec - tmpl[idx].real();
            resid += d * d;
        }
    }
    if (!(resid <= 0.05 * energy))
        return false;

    out.code = code;
    out.nt = static_cast<uint8_t>(nt);
    out.energy = static_cast<float>(energy);
    out.key = shaped_template_key(tmpl, len) ^ code.key;
    out.valid = true;
    return true;
}

const ShapedSparseTemplate* ensure_shaped_sparse_template(
    const std::complex<float>* tmpl, size_t len,
    const SparseCirTemplate& code, ShapedSparseTemplate& out)
{
    const uint64_t key = shaped_template_key(tmpl, len) ^ code.key;
    if (out.valid && out.key == key)
        return &out;
    if (!build_shaped_sparse_template(tmpl, len, code, out))
        return nullptr;
    return &out;
}

MatchedFilterPeak scan_normalized_sync(const SyncMatchedFilter& filter,
                                       const std::complex<float>* rx,
                                       size_t lag_count,
                                       float power,
                                       float template_energy,
                                       float* metrics)
{
    MatchedFilterPeak peak;
    if (lag_count == 0)
        return peak;
    const size_t L = filter.length();
    constexpr size_t kChunk = 256; // multiple of 32; 2 KB stack, zero heap alloc
    float off_pwr[kChunk];
    float on_pwr[kChunk];
    for (size_t base = 0; base < lag_count; base += kChunk) {
        const size_t chunk_len = std::min<size_t>(kChunk, lag_count - base);
        const size_t n_diff =
            (base + chunk_len < lag_count) ? chunk_len : (chunk_len - 1);
        if (n_diff > 0) {
            const float* f_off = reinterpret_cast<const float*>(rx + base);
            const float* f_on = reinterpret_cast<const float*>(rx + base + L);
            for (size_t i = 0; i < n_diff; ++i) {
                const float r0 = f_off[2 * i];
                const float i0 = f_off[2 * i + 1];
                const float r1 = f_on[2 * i];
                const float i1 = f_on[2 * i + 1];
                off_pwr[i] = r0 * r0 + i0 * i0;
                on_pwr[i] = r1 * r1 + i1 * i1;
            }
        }
        filter.for_each_lag(
            rx + base, chunk_len,
            [&](size_t local_lag, std::complex<float> corr, bool) {
                const size_t lag = base + local_lag;
                const float metric =
                    openuwb::core::normalized_correlation_power(
                        corr, power, template_energy);
                if (metrics)
                    metrics[lag] = metric;
                if (metric > peak.metric)
                    peak = { lag, metric, corr };
                if (local_lag < n_diff)
                    power += on_pwr[local_lag] - off_pwr[local_lag];
            });
    }
    return peak;
}

MatchedFilterPeak scan_raw_sync(const SyncMatchedFilter& filter,
                                const std::complex<float>* rx,
                                size_t lag_count,
                                float* powers)
{
    MatchedFilterPeak peak;
    filter.for_each_lag(rx, lag_count, [&](size_t lag, std::complex<float> corr, bool) {
        const float power = std::norm(corr);
        powers[lag] = power;
        if (power > peak.metric)
            peak = { lag, power, corr };
    });
    return peak;
}

void correlate_cir_sparse_raw(const std::complex<float>* corr,
                              const SparseCirTemplate& tmpl,
                              size_t tap_count,
                              std::complex<float>* values)
{
    if (tap_count == 0)
        return;
    const SyncMatchedFilter filter(nullptr, 0, nullptr, &tmpl);
    filter.for_each_lag(corr, tap_count, [&](size_t lag, std::complex<float> acc, bool) {
        values[lag] = acc;
    });
}

void correlate_cir_sparse(const std::complex<float>* corr,
                          const SparseCirTemplate& tmpl,
                          size_t tap_count,
                          std::complex<float>* values)
{
    if (tap_count == 0 || tmpl.code_energy <= 0.0f)
        return;
    correlate_cir_sparse_raw(corr, tmpl, tap_count, values);
    for (size_t i = 0; i < tap_count; ++i)
        values[i] /= tmpl.code_energy;
}

void correlate_cir_sparse(const std::complex<float>* corr,
                          const SparseCirTemplate& tmpl,
                          size_t tap_count,
                          std::vector<std::complex<float>>& values)
{
    values.resize(tap_count);
    if (tmpl.code_energy <= 0.0f || tap_count == 0)
        return;
    correlate_cir_sparse(corr, tmpl, tap_count, values.data());
}

} // namespace openuwb::demod::rxopt
