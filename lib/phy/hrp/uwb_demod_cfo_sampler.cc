// SPDX-License-Identifier: GPL-3.0-or-later
//
// SIMD SYNC-symbol sampler for the ML CFO estimator.
//
// Samples the complex SYNC correlation on the uniform (optionally fractional)
// grid anchored at `origin_sample`:
//     t_k = origin_sample + k * period_samples,  k in [symbol_begin, symbol_end)
// so the CFO estimator never consumes the Stage-1 per-symbol argmax (whose
// multipath tap hopping injects non-physical phase steps).
//
// Every lane goes through an existing SIMD primitive -- no scalar dot product
// is written here:
//   * ShapedSparseTemplate : shaped_ternary_dot_4lag() (one 16-lag batch) or
//                            shaped_ternary_dot_single() near the frame edge;
//   * SparseCirTemplate    : sparse_ternary_dot_4lag() / SyncMatchedFilter;
//   * RealSyncTemplate     : real_template_dot_4lag() / ..._single();
//   * any other waveform   : openuwb::core::conjugate_dot() (VOLK).
//
// Zero per-packet heap allocation: the templates and the code tables live in
// DemodScratch and are only rebuilt when their content key changes.
#include <openuwb/phy/hrp/rx/demod_cfo.h>
#include <openuwb/phy/hrp/rx/rx_matched_filter.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>

namespace openuwb::demod::core {

namespace {

// Template dispatch order, identical to rxopt::SyncMatchedFilter::dot().
enum class SyncPath { Shaped, Sparse, Real, Volk };

// 4-tap Catmull-Rom weights for a fractional delay `d` in [0, 1); the same
// cubic as rxopt::build_fractional_cir_bank():
//   w(-1) = -0.5d + d^2 - 0.5d^3
//   w( 0) =  1   - 2.5d^2 + 1.5d^3
//   w( 1) =  0.5d + 2d^2 - 1.5d^3
//   w( 2) = -0.5d^2 + 0.5d^3
// They sum to 1 and reduce to (0, 1, 0, 0) / (0, 0, 1, 0) at d = 0 / d = 1,
// so the fractional grid matches the integer grid exactly at both ends.
struct CatmullRom4 {
    float wm1 = 0.0f;
    float w0 = 1.0f;
    float w1 = 0.0f;
    float w2 = 0.0f;
};

inline CatmullRom4 catmull_rom_4(float d)
{
    const float d2 = d * d;
    const float d3 = d2 * d;
    CatmullRom4 w;
    w.wm1 = -0.5f * d + d2 - 0.5f * d3;
    w.w0 = 1.0f - 2.5f * d2 + 1.5f * d3;
    w.w1 = 0.5f * d + 2.0f * d2 - 1.5f * d3;
    w.w2 = -0.5f * d2 + 0.5f * d3;
    return w;
}

// double -> int64 conversions that cannot overflow (a wild origin/period would
// otherwise be UB, and UB here is an out-of-bounds read).
inline bool to_int64_floor(double v, int64_t& out)
{
    if (!std::isfinite(v) || v < -9.0e15 || v > 9.0e15)
        return false;
    out = static_cast<int64_t>(std::floor(v));
    return true;
}

inline bool to_int64_round(double v, int64_t& out)
{
    if (!std::isfinite(v) || v < -9.0e15 || v > 9.0e15)
        return false;
    out = static_cast<int64_t>(std::llround(v));
    return true;
}

// `base + extent <= n` without promoting either addend (no size_t overflow).
inline bool fits(size_t base, size_t extent, size_t n)
{
    return base <= n && extent <= n - base;
}

// Largest complex-sample offset of any non-zero chip.  `pos_off` / `neg_off`
// are FLOAT offsets (twice the complex-sample index), matching the wide
// kernels' indexing, so the complex offset is off / 2.
inline size_t max_code_sample_offset(const rxopt::SparseCirTemplate& code)
{
    size_t max_off = 0;
    for (uint16_t i = 0; i < code.pos_count; ++i)
        max_off = std::max<size_t>(max_off, code.pos_off[i] / 2);
    for (uint16_t i = 0; i < code.neg_count; ++i)
        max_off = std::max<size_t>(max_off, code.neg_off[i] / 2);
    return max_off;
}

// One SYNC symbol is kQm35SamplesPerSymbol (1016) complex samples.  Callers in
// this tree pass exactly that one-symbol waveform, but a doubled buffer (the
// template bank form) still describes one symbol, so clamp instead of failing.
inline size_t symbol_length(const std::vector<std::complex<float>>& template_wf)
{
    const size_t sym = kQm35SamplesPerSymbol;
    const size_t len = template_wf.size();
    return len > sym ? sym : len;
}

} // namespace

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
    size_t* out_first_symbol)
{
    // The output is always the contiguous suffix of the requested range that
    // starts at `symbol_begin`: leading out-of-window symbols are NOT skipped
    // (that would change the phase origin of the estimator), the first
    // unsamplable symbol simply ends the run.
    if (out_first_symbol)
        *out_first_symbol = symbol_begin;

    if (!rx || !out_z || n == 0 || template_wf.empty())
        return 0;
    if (symbol_end <= symbol_begin || !std::isfinite(period_samples) ||
        period_samples <= 0.0)
        return 0;

    const size_t L = symbol_length(template_wf);
    if (L == 0 || L > n)
        return 0;

    // --- Template caches (never allocate once warm) -------------------------
    // Full-symbol sparse/shaped templates are only meaningful for the
    // 2 * kQm35ChipsPerSymbol (1016) HRP waveform.
    const bool full_symbol = (L == 2 * kQm35ChipsPerSymbol);
    const rxopt::SparseCirTemplate* sparse =
        detail::prepare_sparse_sync_template(profile, L, scratch);

    const rxopt::ShapedSparseTemplate* shaped = nullptr;
    if (full_symbol) {
        // The sparse code table is needed even when the timing sparse path is
        // off: the shaped factorization is gated by a reconstruction check, so
        // a template that is not code * p simply keeps its old path.
        const rxopt::SparseCirTemplate* code = sparse;
        const bool known_code = profile.code_index >= 9 && profile.code_index <= 12;
        if (!code && known_code) {
            if (scratch.cached_seq.code_index != profile.code_index ||
                scratch.cached_seq.pcode.empty())
                scratch.cached_seq.update(profile.code_index, profile.sfd_mode);
            const auto& pcode = scratch.cached_seq.pcode;
            if (!pcode.empty()) {
                rxopt::ensure_sparse_cir_template(pcode.data(), pcode.size(),
                                                  scratch.cir_template);
                if (scratch.cir_template.code_energy > 0.0f)
                    code = &scratch.cir_template;
            }
        }
        if (code) {
            shaped = detail::ensure_shaped_sparse_template_cached(
                template_wf, L, *code, scratch);
        }
    }

    const rxopt::RealSyncTemplate& real_tmpl = scratch.real_sync_template;
    SyncPath path = SyncPath::Volk;
    if (shaped) {
        path = SyncPath::Shaped;
    } else if (sparse) {
        path = SyncPath::Sparse;
    } else {
        // Only the dense paths need the real-tap duplicate; on the sparse /
        // shaped paths it would be a wasted L2 normalize + copy.
        rxopt::ensure_real_sync_template(template_wf, scratch.real_sync_template);
        if (real_tmpl.is_real && real_tmpl.taps_dup.size() >= 2 * L)
            path = SyncPath::Real;
    }

    const rxopt::SyncMatchedFilter mf(
        template_wf.data(), L,
        path == SyncPath::Real ? real_tmpl.taps_dup.data() : nullptr,
        path == SyncPath::Sparse ? sparse : nullptr,
        /*volk_complex=*/true,
        path == SyncPath::Shaped ? shaped : nullptr);

    // The 16-lag batch reads lag 15 of every sparse offset; a window that ends
    // closer than that to the frame edge falls back to four single-lag calls.
    const size_t shaped_max_off =
        path == SyncPath::Shaped ? max_code_sample_offset(shaped->code) : 0;

    // Exact read extent of one SyncMatchedFilter::dot() at lag 0:
    //   * shaped_ternary_dot_single() runs the 8-lag sparse batch, so it reads
    //     max_code_sample_offset + 8 samples.  A valid template only guarantees
    //     last_chip + nt <= L, so that extent can exceed L by up to 4 samples:
    //     the guard must use the kernel's real width, not L.
    //   * every other kernel reads exactly L.
    const size_t dot_extent =
        path == SyncPath::Shaped ? shaped_max_off + 8 : L;

    // Correlate consecutive integer lags q-1 .. q+2.  Returns false when the
    // window does not fit, which makes the caller stop the sequence instead of
    // reading out of bounds.
    auto gather_4lag = [&](size_t base, std::complex<float> z4[4]) -> bool {
        switch (path) {
        case SyncPath::Shaped:
            if (fits(base, shaped_max_off + 16, n)) {
                rxopt::shaped_ternary_dot_4lag(rx + base, *shaped, z4);
                return true;
            }
            // Safe fallback: lag 3 + pulse tap 7 reaches lag 10 of the sparse
            // batch, i.e. max_code_sample_offset + 11 samples from `base`.
            if (fits(base, shaped_max_off + 11, n)) {
                for (size_t j = 0; j < 4; ++j)
                    z4[j] = rxopt::shaped_ternary_dot_single(rx + base + j,
                                                             *shaped);
                return true;
            }
            return false;
        case SyncPath::Sparse:
            if (fits(base, L + 3, n)) {
                rxopt::sparse_ternary_dot_4lag(rx + base, *sparse, z4);
                return true;
            }
            return false;
        case SyncPath::Real:
            if (fits(base, L + 3, n)) {
                rxopt::real_template_dot_4lag(rx + base,
                                              real_tmpl.taps_dup.data(), L, z4);
                return true;
            }
            return false;
        case SyncPath::Volk:
        default:
            if (fits(base, L + 3, n)) {
                for (size_t j = 0; j < 4; ++j)
                    z4[j] = rxopt::conjugate_dot(rx + base + j,
                                                 template_wf.data(), L);
                return true;
            }
            return false;
        }
    };

    const double origin = static_cast<double>(origin_sample);
    size_t produced = 0;
    for (size_t k = symbol_begin; k < symbol_end; ++k) {
        const double t = origin + static_cast<double>(k) * period_samples;
        std::complex<float> z(0.0f, 0.0f);
        bool ok = false;

        if (!use_fractional_interp) {
            int64_t pos = 0;
            if (to_int64_round(t, pos) && pos >= 0 &&
                fits(static_cast<size_t>(pos), dot_extent, n)) {
                z = mf.dot(rx + pos);
                ok = true;
            }
        } else {
            int64_t q = 0;
            if (to_int64_floor(t, q) && q >= 0 &&
                static_cast<uint64_t>(q) < static_cast<uint64_t>(n)) {
                const float d = static_cast<float>(t - std::floor(t));
                if (d < 1e-4f) {
                    // On the integer grid: no interpolation, no extra lags.
                    if (fits(static_cast<size_t>(q), dot_extent, n)) {
                        z = mf.dot(rx + q);
                        ok = true;
                    }
                } else if (q >= 1) {
                    // z(t) at lag q-1 .. q+2, then 4-tap Catmull-Rom.
                    std::complex<float> z4[4];
                    if (gather_4lag(static_cast<size_t>(q - 1), z4)) {
                        const CatmullRom4 w = catmull_rom_4(d);
                        z = w.wm1 * z4[0] + w.w0 * z4[1] + w.w1 * z4[2] +
                            w.w2 * z4[3];
                        ok = true;
                    }
                }
            }
        }

        // Contiguity contract: the first unsamplable symbol ends the run, so
        // out_z never contains a hole (the ML estimator assumes a uniform
        // grid across every sample it receives).
        if (!ok)
            break;
        out_z[produced++] = z;
    }
    return produced;
}

} // namespace openuwb::demod::core
