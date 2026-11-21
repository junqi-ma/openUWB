// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 3/4 of the split realtime-demodulator core QA: P0/P3 CFO/SFD sweeps
// and the SIMD kernel-equivalence/microbench cases.
#include "demod_core_test_helpers.h"

#include <gnuradio/uwb/uwb_cir_fir_simd.h>
#include <chrono>
#include <iostream>

namespace {

// Serial scalar reference: sum_q conj(values[q]) * rx[win+q]  (forward form).
std::complex<float> cir_fir_serial_ref(const std::complex<float>* rx_win,
                                       const std::complex<float>* values,
                                       size_t tap_count)
{
    std::complex<float> acc(0.0f, 0.0f);
    for (size_t q = 0; q < tap_count; ++q)
        acc += std::conj(values[q]) * rx_win[q];
    return acc;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_demod_core_timing_volk_dot_agrees_with_scalar)
{
    constexpr size_t N = 1016;
    std::mt19937 rng(0x71A11u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    std::vector<std::complex<float>> input(N + 3), taps(N);
    for (auto& x : input)
        x = { dist(rng), dist(rng) };
    for (auto& x : taps)
        x = { dist(rng), dist(rng) };

    for (size_t off = 0; off < 3; ++off) {
        std::complex<float> ref(0.0f, 0.0f);
        for (size_t k = 0; k < N; ++k)
            ref += input[off + k] * std::conj(taps[k]);
        const auto got =
            core::detail::timing_conjugate_dot(input.data() + off, taps.data(), N);
        BOOST_CHECK_LT(std::abs(got - ref), 2e-3f);
    }
}

BOOST_AUTO_TEST_CASE(test_sync_matched_filter_exact_window_edges)
{
    const std::vector<std::complex<float>> taps{ { 1.0f, 0.0f }, { -1.0f, 0.0f } };
    rxopt::RealSyncTemplate real;
    rxopt::build_real_sync_template(taps, real);
    const rxopt::SyncMatchedFilter filter(
        taps.data(), taps.size(), real.taps_dup.data(), nullptr);
    const std::vector<std::complex<float>> single{ { 1.0f, 0.0f }, { -1.0f, 0.0f } };
    filter.for_each_lag(
        single.data(), 1, [](size_t lag, std::complex<float> value, bool batch) {
            BOOST_CHECK_EQUAL(lag, 0u);
            BOOST_CHECK_EQUAL(value, std::complex<float>(2.0f, 0.0f));
            BOOST_CHECK(!batch);
        });

    // A complete batch followed by the last legal lag in an exact-size buffer.
    std::vector<std::complex<float>> samples(10);
    samples[8] = single[0];
    samples[9] = single[1];
    const auto peak =
        rxopt::scan_normalized_sync(filter, samples.data(), 9, 0.0f, real.energy);
    BOOST_CHECK_EQUAL(peak.lag, 8u);
    BOOST_CHECK_EQUAL(peak.metric, 1.0f);
    BOOST_CHECK_EQUAL(peak.corr, std::complex<float>(2.0f, 0.0f));
}

BOOST_AUTO_TEST_CASE(test_demod_core_topk_rake_kernel_agrees_with_scalar)
{
    constexpr size_t N = 48;
    constexpr size_t K = 8;
    std::mt19937 rng(0x70A4u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    std::vector<std::complex<float>> input(N), weights(K);
    const uint8_t indices[K] = { 1, 4, 9, 15, 23, 31, 36, 37 };
    for (auto& x : input)
        x = { dist(rng), dist(rng) };
    for (auto& x : weights)
        x = { dist(rng), dist(rng) };

    std::complex<float> ref(0.0f, 0.0f);
    for (size_t i = 0; i < K; ++i)
        ref += weights[i] * input[indices[i]];
    const auto got = cir_fir::dot_topk(input.data(), indices, weights.data(), K);
    BOOST_CHECK_LT(std::abs(got - ref), 1e-5f);

#if UWB_CIR_FIR_HAVE_AVX2
    // Four outputs at decimation 2, including deliberately unaligned input.
    for (size_t k : { size_t(4), size_t(8) }) {
        std::complex<float> out[4];
        if (k == 4)
            cir_fir::dot_topk_x4_avx2<4>(input.data() + 1, indices, weights.data(), out);
        else
            cir_fir::dot_topk_x4_avx2<8>(input.data() + 1, indices, weights.data(), out);
        for (size_t o = 0; o < 4; ++o) {
            const auto scalar =
                cir_fir::dot_topk(input.data() + 1 + 2 * o, indices, weights.data(), k);
            BOOST_CHECK_LT(std::abs(out[o] - scalar), 1e-5f);
        }
    }
#endif
}

BOOST_AUTO_TEST_CASE(test_demod_core_topk_rake_golden_fcs)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    for (size_t k : { size_t(4), size_t(8) }) {
        auto prof = Qm35825Profile::Default();
        prof.sfd_mode = "ieee";
        prof.cir_rake_top_k = k;
        core::DemodScratch scratch;
        scratch.reserve(iq.size());
        const auto res =
            core::demodulate_one(iq.data(), iq.size(), prof, 1, 9984, 0, tmpl, scratch);
        BOOST_TEST_CONTEXT("Top-K=" << k)
        {
            BOOST_REQUIRE(res.status == DemodStatus::Success);
            BOOST_CHECK(res.payload.fcs_pass);
            BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_bypass_filter_golden_fcs)
{
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const auto tmpl = load_reference_template();
    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    prof.cir_soft_chip_mode = CirSoftChipMode::Bypass;
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    const auto res =
        core::demodulate_one(iq.data(), iq.size(), prof, 1, 9984, 0, tmpl, scratch);
    BOOST_REQUIRE(res.status == DemodStatus::Success);
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
}

BOOST_AUTO_TEST_CASE(test_demod_core_p2_cir_fir_kernels_agree)
{
    // Synthetic random taps + window: all three kernels must match the serial
    // forward-form reference within a tight abs tolerance (FP reassociation).
    constexpr size_t T = 38;
    std::mt19937 rng(0xC1F12u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    std::vector<std::complex<float>> values(T), rx(T + 8);
    for (size_t i = 0; i < T; ++i)
        values[i] = { dist(rng), dist(rng) };
    for (size_t i = 0; i < rx.size(); ++i)
        rx[i] = { dist(rng), dist(rng) };

    alignas(32) std::complex<float> h_conj[64];
    alignas(32) std::complex<float> values_nat[64];
    for (size_t q = 0; q < T; ++q) {
        values_nat[q] = values[q];
        h_conj[q] = std::conj(values[q]);
    }

    // Probe several base offsets so unaligned AVX loads are exercised.
    for (size_t off = 0; off < 4; ++off) {
        const std::complex<float>* win = rx.data() + off;
        const auto ref = cir_fir_serial_ref(win, values.data(), T);

        const auto a = cir_fir::dot_multi_acc8(win, h_conj, T);
        const auto b = cir_fir::dot_volk(win, values_nat, T);
        const auto c = cir_fir::dot_avx2_fixed38(win, h_conj, T);

        auto check_close = [&](std::complex<float> x, const char* name) {
            const float d = std::abs(x - ref);
            BOOST_CHECK_MESSAGE(
                d < 1e-4f, name << " vs serial ref abs-diff=" << d << " off=" << off);
        };
        check_close(a, "MultiAcc8");
        check_close(b, "VOLK");
        check_close(c, "Avx2Fixed38");
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_p2_cir_softchips_still_golden)
{
    // Full CIR stage must keep the existing golden tolerances after the FIR
    // rewrite (L2 CIR < 0.01, soft max-diff < 1e-3, first_path / peak / count).
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE_EQUAL(iq.size(), size_t(319168));

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iq, scratch, prof, tr, cf, cir));

    BOOST_CHECK_EQUAL(cir.first_path_sample, size_t(9986));
    BOOST_CHECK_CLOSE(cir.cir_peak_metric, 0.732856f, 1.0);
    BOOST_CHECK_EQUAL(cir.soft_chip_count, size_t(154578));
    BOOST_CHECK_CLOSE(cir.samples_per_chip, 2.0, 0.5);
    BOOST_CHECK_EQUAL(tr.preamble_start_sample + (int64_t)cir.post_samples - 1,
                      int64_t(10013));

    std::vector<float> gcir;
    BOOST_REQUIRE(load_f32(golden_dir() + "/stage_cir.f32", gcir));
    float l2 = 0.0f;
    for (size_t i = 0; i < gcir.size(); ++i) {
        const float d = cir.cir_values[i] - gcir[i];
        l2 += d * d;
    }
    BOOST_CHECK_SMALL(std::sqrt(l2), 0.01f);

    std::vector<float> gsoft;
    BOOST_REQUIRE(load_f32(golden_dir() + "/stage_softchips.f32", gsoft));
    BOOST_REQUIRE_EQUAL(scratch.soft_chips.size(), gsoft.size());
    float mx = 0.0f;
    for (size_t i = 0; i < gsoft.size(); ++i)
        mx = std::max(mx, std::abs(scratch.soft_chips[i] - gsoft[i]));
    BOOST_CHECK_LT(mx, 1e-3f);

    // Full pipeline still yields PHR/payload FCS 0x584b.
    auto prof2 = Qm35825Profile::Default();
    prof2.sfd_mode = "ieee";
    core::DemodScratch scratch2;
    scratch2.reserve(iq.size());
    const auto res = core::demodulate_one(
        iq.data(), iq.size(), prof2, 1, 9984, 0, load_reference_template(), scratch2);
    BOOST_REQUIRE(res.status == DemodStatus::Success);
    BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
    BOOST_CHECK_EQUAL(res.payload.calculated_fcs, uint16_t(0x584b));
    BOOST_CHECK(res.payload.fcs_pass);
}

BOOST_AUTO_TEST_CASE(test_demod_core_p2_cir_fir_microbench)
{
    // Microbench the three FIR candidates over the golden chip grid using the
    // real CIR taps.  Prints wall µs so the P2 report can compare a/b/c
    // without touching benchmark_detector.cc.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    auto prof = Qm35825Profile::Default();
    TimingResult tr;
    CfoResult cf;
    CirResult cir;
    BOOST_REQUIRE(run_cir_softchips(iq, scratch, prof, tr, cf, cir));
    BOOST_REQUIRE_EQUAL(cir.cir_values.size(), size_t(38));

    // Rebuild complex CIR taps by re-running the estimate portion is heavy;
    // soft chips already passed golden, so for the microbench re-derive h
    // from the real-valued CIR taps (Q≈0 golden channel) as complex.
    constexpr size_t T = 38;
    alignas(32) std::complex<float> values_nat[64];
    alignas(32) std::complex<float> h_conj[64];
    for (size_t q = 0; q < T; ++q) {
        values_nat[q] = std::complex<float>(cir.cir_values[q], 0.0f);
        h_conj[q] = std::conj(values_nat[q]);
    }

    // Use the CFO-compensated buffer the stage already filled.
    const std::complex<float>* rx = scratch.derotated.data();
    const size_t n = iq.size();
    const int64_t chip_start =
        tr.preamble_start_sample + static_cast<int64_t>(cir.post_samples) - 1;
    const int64_t spc_i = 2;
    const int64_t fwd = static_cast<int64_t>(T) - 1;
    const size_t num_chips = cir.soft_chip_count;
    BOOST_REQUIRE(chip_start >= fwd);
    BOOST_REQUIRE(static_cast<size_t>(chip_start + (int64_t)(num_chips - 1) * spc_i) < n);

    auto time_kernel = [&](cir_fir::Kernel k, const char* name) {
        std::vector<std::complex<float>> out(num_chips);
        // Warmup
        for (size_t i = 0; i < num_chips; i += 1024) {
            const size_t win = static_cast<size_t>(chip_start + (int64_t)i * spc_i - fwd);
            out[i] = cir_fir::dot(k, rx + win, h_conj, values_nat, T);
        }
        auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < num_chips; ++i) {
            const size_t win = static_cast<size_t>(chip_start + (int64_t)i * spc_i - fwd);
            out[i] = cir_fir::dot(k, rx + win, h_conj, values_nat, T);
        }
        auto t1 = std::chrono::steady_clock::now();
        const double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        std::cout << "[P2 FIR microbench] " << name << " chips=" << num_chips
                  << " fir_us=" << us << "  (" << (1e-6 * num_chips / (us * 1e-6))
                  << " Mchips/s)\n";

        // Spot-check vs serial on a few chips (reassociation ok within 1e-3).
        for (size_t i : { size_t(0), size_t(1000), num_chips / 2, num_chips - 1 }) {
            const size_t win = static_cast<size_t>(chip_start + (int64_t)i * spc_i - fwd);
            const auto ref = cir_fir_serial_ref(rx + win, values_nat, T);
            BOOST_CHECK_LT(std::abs(out[i] - ref), 1e-3f);
        }
        return us;
    };

    const double ua = time_kernel(cir_fir::Kernel::MultiAcc8, "a_MultiAcc8");
    const double ub = time_kernel(cir_fir::Kernel::Volk, "b_VOLK");
    const double uc = time_kernel(cir_fir::Kernel::Avx2Fixed, "c_Avx2Fixed38");
    std::cout << "[P2 FIR microbench] summary a=" << ua << " b=" << ub << " c=" << uc
              << " us (default kernel=" << static_cast<int>(cir_fir::kDefaultKernel)
              << ")\n";
    BOOST_CHECK_GT(ua, 0.0);
    BOOST_CHECK_GT(ub, 0.0);
    BOOST_CHECK_GT(uc, 0.0);

    // 203776-sample scheduling window: re-run full CIR stage on a truncated
    // view of the golden buffer and report estimate/fir/post (O(num_chips)).
    {
        constexpr size_t n_sched = 203776;
        BOOST_REQUIRE(iq.size() >= n_sched);
        core::DemodScratch scr2;
        scr2.reserve(n_sched);
        // Reuse timing/cfo from the full window where possible, but CIR stage
        // needs derotated of length n_sched: re-run stages on the prefix.
        TimingResult tr2;
        CfoResult cf2;
        CirResult cir2;
        std::vector<gr_complex> prefix(iq.begin(), iq.begin() + (long)n_sched);
        BOOST_REQUIRE(run_cir_softchips(prefix, scr2, prof, tr2, cf2, cir2));
        std::cout << "[P2 CIR window] n=319168 chips=" << cir.soft_chip_count
                  << " estimate=" << cir.cir_estimate_us << " fir=" << cir.soft_fir_us
                  << " post=" << cir.postprocess_us << " total="
                  << (cir.cir_estimate_us + cir.soft_fir_us + cir.postprocess_us)
                  << " us\n";
        std::cout << "[P2 CIR window] n=" << n_sched << " chips=" << cir2.soft_chip_count
                  << " estimate=" << cir2.cir_estimate_us << " fir=" << cir2.soft_fir_us
                  << " post=" << cir2.postprocess_us << " total="
                  << (cir2.cir_estimate_us + cir2.soft_fir_us + cir2.postprocess_us)
                  << " us\n";
        BOOST_CHECK_LT(cir2.soft_chip_count, cir.soft_chip_count);
    }
}

// ---------------------------------------------------------------------------
// R7: code-10 / DW1000 profile self-consistency + representative robustness.
// Cyclic origin of kPreambleCode10 is locked to the MATLAB-generated
// testdata/uwb_code10_preamble16_payload8.cfile single-symbol template.
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_demod_core_r7_code10_self_consistency)
{
    // code-10: ternary {-1,0,+1}, length 127, 64 non-zeros, energy 64.
    // code-9 path must be UNCHANGED (same pointer content as before R7).
    const int8_t* c9 = GetPreambleCode(9);
    const int8_t* c10 = GetPreambleCode(10);
    BOOST_REQUIRE(c9 != nullptr);
    BOOST_REQUIRE(c10 != nullptr);
    BOOST_CHECK(c9 != c10);

    size_t nnz9 = 0, nnz10 = 0, e9 = 0, e10 = 0;
    for (size_t i = 0; i < kQm35CodeLength; ++i) {
        BOOST_CHECK(c9[i] == -1 || c9[i] == 0 || c9[i] == 1);
        BOOST_CHECK(c10[i] == -1 || c10[i] == 0 || c10[i] == 1);
        if (c9[i] != 0)
            ++nnz9;
        if (c10[i] != 0)
            ++nnz10;
        e9 += static_cast<size_t>(std::abs(c9[i]));
        e10 += static_cast<size_t>(std::abs(c10[i]));
    }
    BOOST_CHECK_EQUAL(nnz9, size_t(64));
    BOOST_CHECK_EQUAL(e9, size_t(64));
    BOOST_CHECK_EQUAL(nnz10, size_t(64));
    BOOST_CHECK_EQUAL(e10, size_t(64));

    // MATLAB R2025b lrwpan.internal.HRPCodes(9:12) lock.  FNV-1a over the
    // exact ternary sequence catches a cyclic origin error as well as a wrong
    // chip value, while keeping this QA compact.
    const std::array<size_t, 4> code_indices = { 9, 10, 11, 12 };
    const std::array<uint64_t, 4> matlab_hashes = { UINT64_C(0xcba009f2c3a1660a),
                                                    UINT64_C(0xb90c5ccecde2762c),
                                                    UINT64_C(0x5ca2e4c081cae608),
                                                    UINT64_C(0xb95c6c21c2aecc56) };
    for (size_t code_i = 0; code_i < code_indices.size(); ++code_i) {
        const int8_t* code = GetPreambleCode(code_indices[code_i]);
        BOOST_REQUIRE(code != nullptr);
        uint64_t hash = UINT64_C(1469598103934665603);
        size_t nnz = 0;
        for (size_t chip = 0; chip < kQm35CodeLength; ++chip) {
            BOOST_CHECK(code[chip] == -1 || code[chip] == 0 || code[chip] == 1);
            if (code[chip] != 0)
                ++nnz;
            hash ^= static_cast<uint8_t>(code[chip] + 1);
            hash *= UINT64_C(1099511628211);
        }
        BOOST_CHECK_EQUAL(nnz, size_t(64));
        BOOST_CHECK_EQUAL(hash, matlab_hashes[code_i]);
    }

    // First 30 output bits of MATLAB R2025b
    // lrwpan.internal.createScrambler(code, 30, 0), mapped 0 -> +1 and
    // 1 -> -1 by bprf_spreading. This validates code-11/12 data de-spreading,
    // not just their preamble sequences.
    const std::array<std::array<int8_t, 30>, 4> matlab_spreading = {
        { { { 1,  -1, 1, 1, 1, 1,  -1, 1, 1,  -1, -1, -1, -1, 1,  -1,
              -1, -1, 1, 1, 1, -1, -1, 1, -1, 1,  1,  1,  -1, -1, 1 } },
          { { 1, 1,  -1, -1, 1, 1,  -1, 1, 1, 1, 1,  -1, -1, -1, -1,
              1, -1, 1,  -1, 1, -1, -1, 1, 1, 1, -1, 1,  1,  1,  -1 } },
          { { -1, -1, -1, -1, 1, 1,  -1, 1, 1,  -1, -1, 1, 1,  -1, -1,
              1,  1,  1,  -1, 1, -1, -1, 1, -1, 1,  -1, 1, -1, 1,  -1 } },
          { { -1, 1, 1,  -1, -1, -1, -1, -1, 1, 1,  -1, -1, -1, 1, 1,
              -1, 1, -1, 1,  1,  1,  1,  -1, 1, -1, 1,  1,  -1, 1, -1 } } }
    };
    for (size_t code_i = 0; code_i < code_indices.size(); ++code_i) {
        std::vector<int8_t> spreading;
        BOOST_REQUIRE(core::detail::bprf_spreading(
            spreading, code_indices[code_i], 0, matlab_spreading[code_i].size()));
        BOOST_CHECK_EQUAL_COLLECTIONS(spreading.begin(),
                                      spreading.end(),
                                      matlab_spreading[code_i].begin(),
                                      matlab_spreading[code_i].end());
    }

    // Sampled code for code-10: non-zeros placed every SpreadingFactor chips.
    auto samp10 = BuildSampledCode(c10, kQm35CodeLength);
    BOOST_REQUIRE_EQUAL(samp10.size(), kQm35ChipsPerSymbol);
    size_t snz = 0, se = 0;
    for (int8_t v : samp10) {
        if (v != 0)
            ++snz;
        se += static_cast<size_t>(std::abs(v));
    }
    BOOST_CHECK_EQUAL(snz, size_t(64));
    BOOST_CHECK_EQUAL(se, size_t(64));

    // code-9 sampled path unchanged: first non-zero is +1 at index 0.
    auto samp9 = BuildSampledCode(c9, kQm35CodeLength);
    BOOST_CHECK_EQUAL(samp9[0], int8_t(1));
    BOOST_CHECK_EQUAL(static_cast<int>(c9[0]), 1);
    BOOST_CHECK_EQUAL(static_cast<int>(c9[7]), -1);

    // Cyclic origin lock: sparse-grid despread of the MATLAB code-10 SYNC
    // template must match kPreambleCode10 (roll-23 origin).  Sample phase 2
    // is the pulse peak of the lrwpan waveform at 998.4 MS/s.
    {
        std::vector<gr_complex> ref;
        BOOST_REQUIRE(
            load_cf32("../../../testdata/uwb_code10_preamble16_payload8.cfile", ref));
        BOOST_REQUIRE_GE(ref.size(), kQm35SamplesPerSymbol);
        // DC-remove first symbol (same as the capture driver).
        gr_complex mean(0.0f, 0.0f);
        for (size_t i = 0; i < kQm35SamplesPerSymbol; ++i)
            mean += ref[i];
        mean /= static_cast<float>(kQm35SamplesPerSymbol);
        std::vector<float> chips(kQm35CodeLength, 0.0f);
        float e_chips = 0.0f;
        for (size_t c = 0; c < kQm35CodeLength; ++c) {
            const size_t idx = 2 + c * 8; // sample phase 2, SF=4, 2 samp/chip
            BOOST_REQUIRE_LT(idx, kQm35SamplesPerSymbol);
            const gr_complex v = ref[idx] - mean;
            chips[c] = v.real(); // waveform is essentially real after align
            e_chips += chips[c] * chips[c];
        }
        // Project to real after a single global phase from the complex grid.
        gr_complex gain(0.0f, 0.0f);
        for (size_t c = 0; c < kQm35CodeLength; ++c) {
            const size_t idx = 2 + c * 8;
            gain += (ref[idx] - mean) * static_cast<float>(c10[c]);
        }
        // Normalized correlation of phase-aligned chips vs code.
        float corr = 0.0f, e_code = 0.0f, e_rx = 0.0f;
        const float ang = std::arg(gain);
        const gr_complex rot(std::cos(-ang), std::sin(-ang));
        for (size_t c = 0; c < kQm35CodeLength; ++c) {
            const size_t idx = 2 + c * 8;
            const gr_complex v = (ref[idx] - mean) * rot;
            const float r = v.real();
            corr += r * static_cast<float>(c10[c]);
            e_rx += r * r;
            e_code += static_cast<float>(c10[c] * c10[c]);
        }
        const float norm_corr = corr / (std::sqrt(e_rx * e_code) + 1e-12f);
        // Wrong cyclic origin (pre-fix draft) scored ~0.64; correct is ~0.999.
        BOOST_CHECK_GT(norm_corr, 0.95f);
        BOOST_CHECK_EQUAL(static_cast<int>(c10[0]), 1);
        BOOST_CHECK_EQUAL(static_cast<int>(c10[1]), 1);
    }

    // Dw1000Profile factory: code_index=10, 64 SYNC, converts to Qm35825 layout.
    auto dw = Dw1000Profile::Default();
    BOOST_CHECK_EQUAL(dw.code_index, size_t(10));
    BOOST_CHECK_EQUAL(dw.preamble_repetitions, size_t(64));
    auto qm = dw.as_qm35825();
    BOOST_CHECK_EQUAL(qm.code_index, size_t(10));
    BOOST_CHECK_EQUAL(qm.preamble_repetitions, size_t(64));
    BOOST_CHECK_CLOSE(qm.sample_rate, kQm35SampleRate, 1e-9);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r7_code10_matlab_waveform_phr)
{
    // Full frame generated by testdata/generate_uwb_comm_sample.m with
    // lrwpanWaveformGenerator: code 10, 16 SYNC, IEEE legacy SFD and an
    // 8-byte PSDU.  This locks the code-index-dependent BPM scrambler, not
    // merely the code-10 preamble cyclic origin.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(
        load_cf32("../../../testdata/uwb_code10_preamble16_payload8.cfile", iq));
    BOOST_REQUIRE_GT(iq.size(), size_t(16 * kQm35SamplesPerSymbol));

    std::vector<gr_complex> tmpl(iq.begin(), iq.begin() + kQm35SamplesPerSymbol);
    // The generator ends exactly on the final payload sample.  Supply the
    // receive guard needed by the CIR matched filter, as a live PDU does.
    iq.resize(iq.size() + 128, gr_complex(0.0f, 0.0f));
    auto prof = Qm35825Profile::Default();
    prof.code_index = 10;
    prof.preamble_repetitions = 16;
    prof.sfd_mode = "ieee";
    prof.cir_skip_initial_repetitions = 0;
    prof.cir_repetitions = 16;

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    const auto result =
        core::demodulate_one(iq.data(), iq.size(), prof, 10, 0, 0, tmpl, scratch);
    BOOST_REQUIRE(result.timing.ok);
    BOOST_REQUIRE(result.sfd.ok);
    BOOST_REQUIRE(result.ns_sfd.ok);
    BOOST_REQUIRE(result.phr.ok);
    BOOST_CHECK(!result.phr.secded_uncorrectable);
    BOOST_CHECK_EQUAL(result.phr.psdu_length, uint32_t(8));
    BOOST_CHECK_CLOSE(result.phr.data_rate_mbps, 6.81f, 0.5);
    BOOST_REQUIRE(result.payload.ok);
    BOOST_REQUIRE_EQUAL(result.payload.bytes.size(), size_t(8));
    for (size_t i = 0; i < result.payload.bytes.size(); ++i)
        BOOST_CHECK_EQUAL(result.payload.bytes[i], static_cast<uint8_t>(i + 1));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r7_long_preamble_cir_honors_explicit_skip)
{
    // MATLAB acceleration estimateCir uses [skip, skip+count) when an
    // explicit cir_skip_initial_repetitions is configured.  For a 256-SYNC
    // DW1000 frame this must not silently change to the final `count` SYNCs.
    auto prof = Qm35825Profile::Default();
    prof.code_index = 10;
    prof.preamble_repetitions = 256;
    prof.cir_skip_initial_repetitions = 10;
    prof.cir_repetitions = 64;
    prof.max_psdu_bytes = 1;

    constexpr int64_t start = 128;
    constexpr size_t late_delay = 12;
    const size_t n = static_cast<size_t>(start) +
                     prof.preamble_repetitions * kQm35SamplesPerSymbol +
                     48 * kQm35SamplesPerSymbol;
    std::vector<gr_complex> iq(n, gr_complex(0.0f, 0.0f));

    const int8_t* code = GetPreambleCode(10);
    const auto spread = BuildSampledCode(code, kQm35CodeLength);
    for (size_t rep = 0; rep < prof.preamble_repetitions; ++rep) {
        const size_t delay = (rep >= 10 && rep < 74) ? late_delay : 0;
        const size_t base =
            static_cast<size_t>(start) + rep * kQm35SamplesPerSymbol + delay;
        for (size_t chip = 0; chip < spread.size(); ++chip) {
            const size_t sample = base + 2 * chip;
            BOOST_REQUIRE_LT(sample, iq.size());
            iq[sample] += gr_complex(static_cast<float>(spread[chip]), 0.0f);
        }
    }

    TimingResult timing;
    timing.ok = true;
    timing.preamble_start_sample = start;
    timing.measured_period = static_cast<double>(kQm35SamplesPerSymbol);
    timing.detected_peaks = prof.preamble_repetitions;

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    CirResult cir;
    std::vector<int8_t> pcode(code, code + kQm35CodeLength);
    BOOST_REQUIRE(core::stage_cir_softchips(
        iq.data(), iq.size(), prof, timing, pcode, cir, scratch));
    BOOST_CHECK_EQUAL(cir.first_path_sample, static_cast<size_t>(start) + late_delay);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r7_fractional_chip_grid_tracks_sfo)
{
    auto prof = Qm35825Profile::Default();
    prof.code_index = 10;
    prof.preamble_repetitions = 64;
    prof.sfd_mode = "decawave";
    prof.max_psdu_bytes = 1;

    constexpr int64_t start = 128;
    constexpr double period = 1016.25;
    const auto sfd = GetSfdSequence("decawave");
    const int8_t* code = GetPreambleCode(10);
    const auto spread = BuildSampledCode(code, kQm35CodeLength);
    const size_t symbol_count = prof.preamble_repetitions + sfd.size();
    const size_t n = static_cast<size_t>(start + std::ceil((symbol_count + 16) * period));
    std::vector<gr_complex> iq(n, gr_complex(0.0f, 0.0f));

    auto add_symbol = [&](size_t symbol, float sign) {
        const double base = static_cast<double>(start) + symbol * period;
        for (size_t chip = 0; chip < spread.size(); ++chip) {
            const float value = sign * static_cast<float>(spread[chip]);
            if (value == 0.0f)
                continue;
            const double pos = base + 2.0 * static_cast<double>(chip);
            const size_t i0 = static_cast<size_t>(std::floor(pos));
            const float frac = static_cast<float>(pos - std::floor(pos));
            BOOST_REQUIRE_LT(i0 + 1, iq.size());
            iq[i0] += gr_complex(value * (1.0f - frac), 0.0f);
            iq[i0 + 1] += gr_complex(value * frac, 0.0f);
        }
    };
    for (size_t rep = 0; rep < prof.preamble_repetitions; ++rep)
        add_symbol(rep, 1.0f);
    for (size_t k = 0; k < sfd.size(); ++k)
        add_symbol(prof.preamble_repetitions + k, static_cast<float>(sfd[k]));

    TimingResult timing;
    timing.ok = true;
    timing.preamble_start_sample = start;
    timing.measured_period = period;
    timing.detected_peaks = prof.preamble_repetitions;

    core::DemodScratch scratch;
    scratch.reserve(iq.size());
    CirResult cir;
    std::vector<int8_t> pcode(code, code + kQm35CodeLength);
    BOOST_REQUIRE(core::stage_cir_softchips(
        iq.data(), iq.size(), prof, timing, pcode, cir, scratch));
    NsSfdResult ns;
    BOOST_REQUIRE(core::stage_ns_sfd(
        scratch.soft_chips, prof, sfd, kQm35ChipsPerSymbol, ns, scratch));
    BOOST_CHECK_GE(ns.metric, 0.75f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_r7_robust_awgn_representative)
{
    // One representative AWGN point (20 dB, N=3 seeds) — full sweep lives in
    // benchmark_detector demod-robust.  Pass = Success + FCS 0x584b.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    double sp = 0.0;
    for (size_t i = 0; i < iq.size(); ++i)
        sp += static_cast<double>(std::norm(iq[i]));
    sp /= static_cast<double>(iq.size());
    const double snr_db = 20.0;
    const double sigma = std::sqrt(sp / std::pow(10.0, snr_db / 10.0) / 2.0);

    int pass = 0;
    for (unsigned seed = 0; seed < 3; ++seed) {
        std::mt19937 rng(1000u + seed);
        std::normal_distribution<float> g(0.0f, static_cast<float>(sigma));
        std::vector<gr_complex> iqn(iq.size());
        for (size_t i = 0; i < iq.size(); ++i)
            iqn[i] = iq[i] + gr_complex(g(rng), g(rng));
        auto res = demod_window_with_cfo(iqn, 0.0);
        if (res.status == DemodStatus::Success && res.payload.fcs_pass &&
            res.payload.received_fcs == uint16_t(0x584b))
            ++pass;
    }
    BOOST_CHECK_GE(pass, 2); // at least 2/3 at 20 dB
}

BOOST_AUTO_TEST_CASE(test_demod_core_rxopt_cfo_block_rotation_boundaries)
{
    // Experimental CFO kernel: 1024-sample blocks, non-multiple lengths and
    // positive/negative CFO.  Check against the independent double reference.
    std::vector<gr_complex> src(4097);
    for (size_t i = 0; i < src.size(); ++i) {
        src[i] = gr_complex(std::sin(0.013f * static_cast<float>(i)) + 0.2f,
                            std::cos(0.021f * static_cast<float>(i)) - 0.1f);
    }
    for (size_t n : { size_t(1000), size_t(1024), size_t(1025), size_t(4097) }) {
        for (double cfo_hz : { 50000.0, -50000.0, 300000.0 }) {
            const double w = 2.0 * M_PI * cfo_hz / kQm35SampleRate;
            std::vector<gr_complex> got, ref(n);
            rxopt::derotate_block(got, src.data(), n, w);
            rxopt::derotate_reference(ref.data(), src.data(), n, w);
            double max_err = 0.0;
            for (size_t i = 0; i < n; ++i)
                max_err =
                    std::max(max_err, static_cast<double>(std::abs(got[i] - ref[i])));
            BOOST_CHECK_LT(max_err, 2e-5);
            BOOST_CHECK_LT(std::abs(got[0] - src[0]), 1e-7f);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_r7_robust_multipath_representative)
{
    // One multipath point: gain 0.4, delay 120 samples — first packet must
    // still decode (FCS 0x584b).
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const float gain = 0.4f;
    const size_t delay = 120;
    std::vector<gr_complex> m(iq.size() + delay, gr_complex(0.0f, 0.0f));
    std::copy(iq.begin(), iq.end(), m.begin());
    for (size_t i = 0; i < iq.size(); ++i)
        m[i + delay] += gain * iq[i];

    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    core::DemodScratch scratch;
    scratch.reserve(m.size());
    auto res = core::demodulate_one(
        m.data(), m.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
    BOOST_CHECK(res.status == DemodStatus::Success);
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
}

BOOST_AUTO_TEST_CASE(test_demod_core_r7_robust_collision_representative)
{
    // Late packet B preamble lands on A's SFD start (offset 0 symbols).
    // Demod of A must not crash; report FCS of A (may pass or fail).
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    const size_t sfd_start = 75008;
    const size_t total = sfd_start + iq.size();
    std::vector<gr_complex> mix(total, gr_complex(0.0f, 0.0f));
    std::copy(iq.begin(), iq.end(), mix.begin());
    for (size_t i = 0; i < iq.size(); ++i)
        mix[sfd_start + i] += iq[i];

    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    core::DemodScratch scratch;
    scratch.reserve(mix.size());
    auto t0 = std::chrono::steady_clock::now();
    auto res = core::demodulate_one(
        mix.data(), mix.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
    auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    // Graceful: no exception (we got here), bounded latency (< 100 ms).
    BOOST_CHECK_LT(ms, 100.0);
    // Status is one of the defined enum values (no crash / no hang).
    BOOST_CHECK(static_cast<int>(res.status) >= 0);
    BOOST_TEST_MESSAGE("collision offset0: status=" << static_cast<int>(res.status)
                                                    << " fcs=" << res.payload.fcs_pass
                                                    << " latency_ms=" << ms);
}
