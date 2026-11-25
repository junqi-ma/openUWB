// SPDX-License-Identifier: GPL-3.0-or-later
//
// K-symbol periodic-fold acquisition engine for the HRP timing stage.
// Split out of uwb_demod_timing.cc so the fold path (low-SNR fallback) can be
// changed and tested without rebuilding the Stage-1 orchestration unit.
// Every buffer used on the hot path is caller-owned (DemodScratch) so a warm
// packet performs no heap allocation.
#include <openuwb/phy/hrp/rx/demod_timing.h>
#include <openuwb/phy/hrp/rx/demod_sync.h>
#include <openuwb/phy/hrp/rx/rx_matched_filter.h>
#include <openuwb/phy/hrp/rx/rx_sync_opt.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace openuwb::demod::core {
namespace detail {

FoldPhaseProbe probe_fold_phase(const std::complex<float>* rx,
                                size_t n,
                                int64_t centre,
                                size_t repetitions)
{
    const size_t sym = kQm35SamplesPerSymbol;
    if (repetitions < 32)
        return {};
    const int64_t lo = centre + static_cast<int64_t>((repetitions / 4) * sym);
    if (lo < 0 || static_cast<size_t>(lo) >= n ||
        n - static_cast<size_t>(lo) <= sym)
        return {};
    const size_t count = std::min(std::min<size_t>(repetitions / 2, 128) * sym,
                                  n - static_cast<size_t>(lo) - sym);
    if (count < 16 * sym)
        return {};
    const auto a = rx + lo;
    FoldPhaseProbe probe;
    probe.samples = a;
    probe.count = count;
    probe.corr = timing_conjugate_dot(a + sym, a, count);
    probe.ea = timing_conjugate_dot(a, a, count).real();
    float head = 0.0f;
    float tail = 0.0f;
    for (size_t i = 0; i < sym; ++i) {
        head += std::norm(a[i]);
        tail += std::norm(a[count + i]);
    }
    probe.eb = probe.ea - head + tail;
    return probe;
}

FoldAcquireResult fold_acquire(const std::complex<float>* rx,
                               size_t n,
                               int64_t centre,
                               int64_t w0_in,
                               int64_t w1_in,
                               const std::complex<float>* tmpl,
                               size_t L,
                               float Et,
                               const rxopt::RealSyncTemplate& real_tmpl,
                               size_t max_symbols,
                               size_t preamble_repetitions,
                               float threshold,
                               double phase_per_symbol,
                               const rxopt::SparseCirTemplate* sparse_tmpl,
                               std::vector<std::complex<float>>* fold_buf,
                               const rxopt::ShapedSparseTemplate* shaped_tmpl,
                               std::vector<float>* lag_metrics_buf)
{
    FoldAcquireResult res;
    const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t Li = static_cast<int64_t>(L);
    const bool have_sparse = sparse_tmpl && sparse_tmpl->code_energy > 0.0f;
    const bool have_shaped = shaped_tmpl && shaped_tmpl->valid;
    if (!rx || L == 0 || n < L || (!tmpl && !have_sparse && !have_shaped))
        return res;

    int64_t w0 = std::max<int64_t>(0, w0_in);
    int64_t w1 = std::min<int64_t>(w1_in, static_cast<int64_t>(n) - Li);
    if (w1 < w0)
        return res;

    // The deepest sample the fold touches is w1 + L - 1 + (K-1)*P, so the
    // buffer must hold K whole symbol periods starting at w0.
    const int64_t room = static_cast<int64_t>(n) - w0 - Li;
    if (room < 0)
        return res;
    int64_t k_eff = std::min<int64_t>(static_cast<int64_t>(max_symbols),
                                      room / P + 1);
    // Folding past the end of the preamble only adds noise, never signal.
    k_eff = std::min<int64_t>(k_eff,
                              static_cast<int64_t>(preamble_repetitions));
    const int64_t min_syms = std::min<int64_t>(
        static_cast<int64_t>(kFoldMinSymbols),
        static_cast<int64_t>(max_symbols));
    if (k_eff < min_syms)
        return res;
    w1 = std::min<int64_t>(w1, static_cast<int64_t>(n) - Li - (k_eff - 1) * P);
    if (w1 < w0)
        return res;

    // Only lags inside +/-P/2 of the seed can be selected (see the note below),
    // so the fold and the correlation only have to cover that sub-window.
    // Clamping here instead of inside the correlation loop shrinks the fold
    // accumulation by the same factor: under AcquireLegacy the search radius is
    // 8192 samples, i.e. span 16384 -> 1016 lags (-93.8%), and every discarded
    // lag was unreachable by the argmax anyway.  `room` / `k_eff` above are
    // still derived from the ORIGINAL w0 so the fold depth (and therefore the
    // folded vector itself) is unchanged.
    const int64_t half_period = P / 2;
    const int64_t c0 = std::max<int64_t>(w0, centre - half_period);
    const int64_t c1 = std::min<int64_t>(w1, centre + half_period);
    if (c1 < c0)
        return res; // the seed lies entirely outside the legal window
    w0 = c0;
    w1 = c1;

    const int64_t span = w1 - w0;
    const size_t fold_len = static_cast<size_t>(span + Li);
    // Reuse the worker buffer when the caller provides one: the fold path runs
    // at most once per packet, but it is also reached on retries, and a
    // per-packet allocation is pure allocator traffic on the low-SNR path.
    std::vector<std::complex<float>> local_folded;
    std::vector<std::complex<float>>& folded =
        fold_buf ? *fold_buf : local_folded;
    if (folded.size() < fold_len)
        folded.resize(fold_len);
    std::fill(folded.begin(), folded.begin() + fold_len,
              std::complex<float>(0.0f, 0.0f));
    for (int64_t k = 0; k < k_eff; ++k)
        fold_accumulate(folded.data(), rx + w0 + k * P, fold_len,
                        -k * phase_per_symbol);

    const float* real_taps = real_tmpl.is_real && real_tmpl.length == L
        ? real_tmpl.taps_dup.data() : nullptr;
    const rxopt::SyncMatchedFilter filter(tmpl, L, real_taps, sparse_tmpl,
                                          true, shaped_tmpl);
    if (shaped_tmpl && shaped_tmpl->valid)
        Et = shaped_tmpl->energy;
    else if (sparse_tmpl && sparse_tmpl->code_energy > 0.0f)
        Et = sparse_tmpl->code_energy;
    const size_t lags = static_cast<size_t>(span + 1);
    std::vector<float> local_lag_metrics;
    std::vector<float>& lag_metrics =
        lag_metrics_buf ? *lag_metrics_buf : local_lag_metrics;
    if (lag_metrics.size() < lags)
        lag_metrics.resize(lags);
    const auto peak = rxopt::scan_normalized_sync(
        filter, folded.data(), lags, window_power_8acc(folded.data(), L), Et,
        lag_metrics.data());

    res.symbols_folded = static_cast<size_t>(k_eff);
    res.window_lags = lags;
    bool accepted = (peak.metric >= threshold);
    if (!accepted && peak.metric >= 0.015f) {
        constexpr int64_t kFoldParGuard = 40;
        constexpr double kFoldParMin = 17.7827941; // 12.5 dB
        double bg_sum = 0.0;
        size_t bg_cnt = 0;
        const int64_t pk_i = static_cast<int64_t>(peak.lag);
        for (size_t i = 0; i < lags; ++i) {
            if (std::abs(static_cast<int64_t>(i) - pk_i) > kFoldParGuard) {
                bg_sum += static_cast<double>(lag_metrics[i]);
                ++bg_cnt;
            }
        }
        if (bg_cnt >= 64) {
            const double bg_mean =
                std::max(bg_sum / static_cast<double>(bg_cnt), 1e-12);
            accepted = (static_cast<double>(peak.metric) / bg_mean) >= kFoldParMin;
        }
    }
    if (!accepted)
        return res;
    res.ok = true;
    res.start = w0 + static_cast<int64_t>(peak.lag);
    res.metric = peak.metric;
    res.corr = peak.corr;
    return res;
}

FoldAcquireResult fold_probe_acquire(
    const std::complex<float>* rx,
    size_t n,
    int64_t probe,
    const std::complex<float>* tmpl,
    size_t L,
    float Et,
    const rxopt::RealSyncTemplate& real_tmpl,
    size_t fold_symbols,
    float threshold,
    const rxopt::SparseCirTemplate* sparse_tmpl,
    std::vector<std::complex<float>>* fold_buf,
    const rxopt::ShapedSparseTemplate* shaped_tmpl)
{
    FoldAcquireResult res;
    const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t Li = static_cast<int64_t>(L);
    const bool have_sparse = sparse_tmpl && sparse_tmpl->code_energy > 0.0f;
    const bool have_shaped = shaped_tmpl && shaped_tmpl->valid;
    if (!rx || L == 0 || n < L || fold_symbols == 0 || probe < 0 ||
        (!tmpl && !have_sparse && !have_shaped))
        return res;
    if (fold_symbols > static_cast<size_t>(P))
        fold_symbols = static_cast<size_t>(P);
    const int64_t K = static_cast<int64_t>(fold_symbols);
    // Buffer covers one period of possible start phases plus the template.
    const int64_t fold_len = P + Li - 1;
    // The deepest sample is probe + (K-1)*P + fold_len - 1.
    if (probe + K * P + Li - 1 > static_cast<int64_t>(n))
        return res;

    std::vector<std::complex<float>> local_folded;
    std::vector<std::complex<float>>& folded =
        fold_buf ? *fold_buf : local_folded;
    const size_t flen = static_cast<size_t>(fold_len);
    if (folded.size() < flen)
        folded.resize(flen);
    std::fill(folded.begin(), folded.begin() + flen,
              std::complex<float>(0.0f, 0.0f));
    for (int64_t k = 0; k < K; ++k)
        fold_accumulate(folded.data(), rx + probe + k * P, flen, 0.0);

    const float* real_taps = real_tmpl.is_real && real_tmpl.length == L
        ? real_tmpl.taps_dup.data() : nullptr;
    const rxopt::SyncMatchedFilter filter(tmpl, L, real_taps, sparse_tmpl,
                                          true, shaped_tmpl);
    if (shaped_tmpl && shaped_tmpl->valid)
        Et = shaped_tmpl->energy;
    else if (sparse_tmpl && sparse_tmpl->code_energy > 0.0f)
        Et = sparse_tmpl->code_energy;
    const auto peak = rxopt::scan_normalized_sync(
        filter, folded.data(), static_cast<size_t>(P),
        window_power_8acc(folded.data(), L), Et);

    res.symbols_folded = fold_symbols;
    res.window_lags = static_cast<size_t>(P);
    if (!(peak.metric >= threshold))
        return res;
    res.ok = true;
    res.start = probe + static_cast<int64_t>(peak.lag);
    res.metric = peak.metric;
    res.corr = peak.corr;
    return res;
}

FoldLeadingEdge disambiguate_fold_leading_edge(
    const std::complex<float>* rx,
    size_t n,
    const rxopt::SyncMatchedFilter& filter,
    size_t L,
    float Et,
    const FoldAcquireResult& fold,
    int64_t radius,
    std::vector<float>& cand_metric)
{
    FoldLeadingEdge res;
    res.start = fold.start;
    if (!rx || fold.start < 0 || L == 0)
        return res;
    const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t max_j = radius / P + 2;
    // One pass over the candidate grid: the same metrics are needed twice
    // (once for the strongest slot, once for the first crossing of the
    // leading-edge threshold), and each is a full single-SYMBOL correlation.
    // Caching them halves this cost.  The buffer is caller-owned scratch, so
    // the candidate count stays bounded (radius <= 8192 -> 21 entries).
    const size_t ncand = static_cast<size_t>(2 * max_j + 1);
    if (cand_metric.size() < ncand)
        cand_metric.resize(ncand);
    auto ss_metric = [&](int64_t pos) -> float {
        if (pos < 0 ||
            pos + static_cast<int64_t>(L) > static_cast<int64_t>(n))
            return -1.0f;
        const float pw = window_power_8acc(rx + pos, L);
        const auto acc = filter.dot(rx + pos);
        return normalized_correlation_power(acc, pw, Et);
    };
    for (size_t i = 0; i < ncand; ++i) {
        const int64_t j = static_cast<int64_t>(i) -
                          static_cast<int64_t>(max_j);
        const float m = ss_metric(fold.start + j * P);
        cand_metric[i] = m;
        if (m > 0.0f)
            res.max_here = std::max(res.max_here, m);
    }
    // Leading-edge detection: a slot before the true start is noise (~0.001 at
    // -14 dB) and every slot after it is a SYNC (~0.038), so the first crossing
    // IS the packet start.  The threshold is relative to the strongest slot
    // seen (a noise-only window would otherwise never cross) and floored at
    // 0.003, which is ~3x the largest of 21 exponential noise metrics.
    const float edge_thr = std::max(0.003f, 0.2f * res.max_here);
    for (size_t i = 0; i < ncand; ++i) {
        const float here = cand_metric[i];
        if (here < 0.0f)
            continue;
        if (here >= edge_thr) {
            const int64_t j = static_cast<int64_t>(i) -
                              static_cast<int64_t>(max_j);
            res.start = fold.start + j * P;
            res.crossed = true;
            break;
        }
    }
    return res;
}

bool fold_anchor_fit(const std::complex<float>* rx,
                     size_t n,
                     int64_t start,
                     const std::complex<float>* tmpl,
                     size_t L,
                     float Et,
                     const rxopt::RealSyncTemplate& real_tmpl,
                     size_t groups,
                     double* period_out,
                     int64_t* preamble_start_out,
                     std::vector<int64_t>& peak_samples,
                     std::vector<float>& peak_metrics,
                     std::vector<std::complex<float>>& peak_corr,
                     double phase_per_symbol,
                     const rxopt::SparseCirTemplate* sparse_tmpl,
                     std::vector<std::complex<float>>* sup_buf,
                     const rxopt::ShapedSparseTemplate* shaped_tmpl)
{
    peak_samples.clear();
    peak_metrics.clear();
    peak_corr.clear();
    const int64_t P = static_cast<int64_t>(kQm35SamplesPerSymbol);
    const int64_t Li = static_cast<int64_t>(L);
    const int64_t refine = kFoldAnchorRefine;
    const size_t sup_len = L + 2 * static_cast<size_t>(refine);
    if (groups < 2 || L == 0 || !rx ||
        (!tmpl && !(sparse_tmpl && sparse_tmpl->code_energy > 0.0f) &&
         !(shaped_tmpl && shaped_tmpl->valid)))
        return false;

    // Sparse uniform SELECTION for the public peak outputs (<= 16) and for
    // the dense CFO shadow grid.  The period/origin fit itself still measures
    // ALL `groups` blocks: the lock feedback and the SFD-anchored re-timing
    // are demonstrably sensitive to the intercept noise of a 16-point fit on
    // the loose-lock replay, while the 64-point first pass reproduces the
    // historical estimator exactly.  The saving therefore comes from the
    // second pass (16 full +/-refine searches instead of 64) and the cheap
    // one-dot shadow grid, not from dropping first-pass anchors.
    const size_t n_anchors = std::min<size_t>(groups, kFoldAnchorMaxGroups);
    std::array<size_t, kFoldAnchorMaxGroups> g_sel{};
    size_t n_sel = 0;
    const double g_step = static_cast<double>(groups - 1) /
                          static_cast<double>(n_anchors - 1);
    for (size_t i = 0; i < n_anchors; ++i) {
        const size_t g = static_cast<size_t>(
            std::llround(static_cast<double>(i) * g_step));
        // Consecutive duplicates are impossible for n_anchors <= groups, but
        // dropping them keeps the fit well-conditioned if that ever changes.
        if (n_sel == 0 || g != g_sel[n_sel - 1])
            g_sel[n_sel++] = g;
    }
    if (n_sel < 2)
        return false;

    std::vector<std::complex<float>> local_sup;
    std::vector<std::complex<float>>& sup = sup_buf ? *sup_buf : local_sup;
    if (sup.size() < sup_len)
        sup.resize(sup_len);
    const float* real_taps = real_tmpl.is_real && real_tmpl.length == L
        ? real_tmpl.taps_dup.data() : nullptr;
    const rxopt::SyncMatchedFilter filter(tmpl, L, real_taps, sparse_tmpl,
                                          true, shaped_tmpl);
    if (shaped_tmpl && shaped_tmpl->valid)
        Et = shaped_tmpl->energy;
    else if (sparse_tmpl && sparse_tmpl->code_energy > 0.0f)
        Et = sparse_tmpl->code_energy;

    LinearFit fit;
    size_t used = 0;
    for (size_t g = 0; g < groups; ++g) {
        const int64_t base =
            start + static_cast<int64_t>(g) *
                        static_cast<int64_t>(kFoldAnchorGroupSymbols) * P -
            refine;
        if (base < 0)
            continue;
        if (base + static_cast<int64_t>(sup_len) +
                static_cast<int64_t>(kFoldAnchorGroupSymbols - 1) * P >
            static_cast<int64_t>(n))
            continue;
        std::fill(sup.begin(), sup.end(), std::complex<float>(0.0f, 0.0f));
        for (int64_t j = 0; j < static_cast<int64_t>(kFoldAnchorGroupSymbols);
             ++j)
            fold_accumulate(sup.data(), rx + base + j * P, sup_len,
                            -j * phase_per_symbol);
        const auto peak = rxopt::scan_normalized_sync(
            filter, sup.data(), sup_len - L + 1,
            window_power_8acc(sup.data(), L), Et);
        const int64_t pos = base + static_cast<int64_t>(peak.lag);
        const double x = static_cast<double>(
            static_cast<int64_t>(g) *
            static_cast<int64_t>(kFoldAnchorGroupSymbols));
        const double y = static_cast<double>(pos + Li - 1);
        fit.add(x, y);
        ++used;
    }
    if (used < 2)
        return false;
    const double denom = fit.determinant();
    if (!(std::abs(denom) > 1e-9))
        return false;
    const double p = fit.slope(denom);
    if (!std::isfinite(p) || !(p > 0.0))
        return false;
    // Same 0.5% sanity gate the LockedBounded anchor fit uses: a period that
    // far off the nominal grid means the anchors did not land on the expected
    // symbol indices at all.
    if (std::abs(p - static_cast<double>(P)) > 0.005 * static_cast<double>(P))
        return false;
    const double intercept_end = fit.intercept(p);

    // Sample complex correlation on the fitted linear timing grid rather than
    // at each selected block's independent argmax lag, so multipath tap hopping
    // within +-refine does not inject phase jumps into fit_cfo().  When the
    // fitted drift across the preamble is < 1.5 samples, snap to the exact
    // integer period P around the midpoint so llround() cannot step 1 sample
    // across multipath taps mid-preamble.  The span is the TOTAL preamble
    // baseline (groups), not the sparse anchor count, so the snap decision and
    // the intercept conversion are identical to the old all-groups fit.
    const double full_span_syms = static_cast<double>(
        (groups - 1) * kFoldAnchorGroupSymbols);
    const double mid_x = 0.5 * full_span_syms;
    const bool lock_integer_grid =
        std::abs(p - static_cast<double>(P)) * full_span_syms < 1.5;
    const int64_t grid_origin_end = std::llround(
        intercept_end + (p - static_cast<double>(P)) * mid_x);
    // Single refold pass for the public sparse peaks and the shadow dense CFO
    // grid.  Previously the 16 selected groups were refolded twice (once for
    // the public peaks, once for the grid), i.e. up to `n_sel + groups` refolds.
    // Walking every group once and emitting the correlation to the grid (and,
    // when the group is one of the selected anchors, to the public peak_*
    // outputs as well) keeps both outputs bit-identical while capping the work
    // at `groups` refolds.
    //
    // Shadow dense CFO grid rationale: every 4-SYMBOL group on the same fitted
    // timing grid.  Stage-2 fit_cfo() consumes up to 64 phases at 4-symbol
    // spacing (alias-free to +/-122.8 kHz, 24-symbol startup skip -> 40
    // points), which is exactly the baseline 64-anchor estimator's
    // information.  The public peak_samples/peak_corr stay sparse (<= 16) for
    // the detected_peaks contract.  This pass does one refold + one dot per
    // group (no +/-refine search), so it is a fraction of the period-fit cost.
    //
    // Refold one group on the fitted timing grid and correlate it.  Shared by
    // the dense-grid and sparse-only loops so both walk the identical fill /
    // fold_accumulate / dot sequence (bit-identical outputs).
    auto refold_and_dot = [&](size_t g, int64_t& fitted_end,
                              std::complex<float>& corr) -> bool {
        const int64_t sym_idx = static_cast<int64_t>(g) *
                                static_cast<int64_t>(
                                    kFoldAnchorGroupSymbols);
        const double x = static_cast<double>(sym_idx);
        fitted_end = lock_integer_grid
            ? (grid_origin_end + sym_idx * P)
            : std::llround(intercept_end + p * x);
        const int64_t fitted_start = fitted_end - (Li - 1);
        if (fitted_start < 0)
            return false;
        if (fitted_start + Li +
                static_cast<int64_t>(kFoldAnchorGroupSymbols - 1) * P >
            static_cast<int64_t>(n))
            return false;
        std::fill(sup.begin(), sup.begin() + static_cast<ptrdiff_t>(L),
                  std::complex<float>(0.0f, 0.0f));
        for (int64_t j = 0;
             j < static_cast<int64_t>(kFoldAnchorGroupSymbols); ++j)
            fold_accumulate(sup.data(), rx + fitted_start + j * P, L,
                            -j * phase_per_symbol);
        corr = filter.dot(sup.data());
        return true;
    };
    auto corr_metric = [&](const std::complex<float>& corr) -> float {
        // 8-accumulator windowed power (same primitive the rest of the fold
        // engine uses) instead of a scalar std::norm loop.
        const float rx_energy = window_power_8acc(sup.data(), L);
        const float denom_e = std::sqrt(std::max(1e-12f, rx_energy * Et));
        return denom_e > 0.0f ? std::abs(corr) / denom_e : 0.0f;
    };
    for (size_t i = 0; i < n_sel; ++i) {
        int64_t fitted_end = 0;
        std::complex<float> corr;
        if (!refold_and_dot(g_sel[i], fitted_end, corr))
            continue;
        peak_samples.push_back(fitted_end);
        peak_metrics.push_back(corr_metric(corr));
        peak_corr.push_back(corr);
    }
    if (peak_samples.size() < 2)
        return false;
    *period_out = p;
    *preamble_start_out = std::max<int64_t>(
        0, static_cast<int64_t>(std::llround(intercept_end)) - (Li - 1));
    return true;
}

} // namespace detail

} // namespace openuwb::demod::core
