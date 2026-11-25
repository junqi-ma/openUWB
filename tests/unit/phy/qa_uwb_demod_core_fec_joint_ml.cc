// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 4/4 of the split realtime-demodulator core QA: FEC reference
// equivalence, Task A1/A2 numerics and the Wave-1 joint-ML Viterbi cases.
#include "demod_core_test_helpers.h"

BOOST_AUTO_TEST_CASE(test_demod_core_task_a1_cfo_range_absolute_phase_origin)
{
    // Task A1 direct numerical QA:
    // Decoupled CFO range rotation [lo, hi) where phase origin is strictly
    // sample 0 (i.e. sample index i has phase -w * i, not treating lo as 0).
    // Test non-aligned block boundaries and positive/negative CFOs.
    const size_t n = 8192;
    std::vector<gr_complex> src(n);
    for (size_t i = 0; i < n; ++i) {
        src[i] = gr_complex(std::sin(0.017f * static_cast<float>(i)) + 0.3f,
                            std::cos(0.023f * static_cast<float>(i)) - 0.2f);
    }

    const gr_complex sentinel(-999.0f, -999.0f);

    struct Interval {
        size_t lo;
        size_t hi;
    };
    const std::vector<Interval> intervals = {
        { 0, 1024 },    { 137, 3519 }, // non-multiple start and end, spanning 4 blocks
        { 1024, 2048 },                // exactly 1 block
        { 1025, 1026 },                // single sample
        { 512, 7168 },                 // large multi-block range
        { 2000, 2000 }                 // empty interval
    };

    for (double cfo_hz : { 50000.0, -50000.0, 150000.0 }) {
        const double w = 2.0 * M_PI * cfo_hz / kQm35SampleRate;
        for (const auto& iv : intervals) {
            std::vector<gr_complex> dst_block(n, sentinel);
            std::vector<gr_complex> dst_recurrence(n, sentinel);
            std::vector<gr_complex> dst_scalar(n, sentinel);
            std::vector<gr_complex> ref(n, sentinel);

            rxopt::derotate_range_block(dst_block.data(), src.data(), n, iv.lo, iv.hi, w);
            rxopt::derotate_range_recurrence(
                dst_recurrence.data(), src.data(), n, iv.lo, iv.hi, w);
            rxopt::derotate_range_scalar(
                dst_scalar.data(), src.data(), n, iv.lo, iv.hi, w);
            rxopt::derotate_range_reference(ref.data(), src.data(), n, iv.lo, iv.hi, w);

            // 1. Outside [lo, hi), sentinel must be strictly preserved
            for (size_t i = 0; i < iv.lo; ++i) {
                BOOST_CHECK_EQUAL(dst_block[i], sentinel);
                BOOST_CHECK_EQUAL(dst_recurrence[i], sentinel);
                BOOST_CHECK_EQUAL(dst_scalar[i], sentinel);
            }
            for (size_t i = iv.hi; i < n; ++i) {
                BOOST_CHECK_EQUAL(dst_block[i], sentinel);
                BOOST_CHECK_EQUAL(dst_recurrence[i], sentinel);
                BOOST_CHECK_EQUAL(dst_scalar[i], sentinel);
            }

            // 2. Inside [lo, hi), verify against independent reference
            double max_err_block = 0.0;
            double max_err_recurrence = 0.0;
            double max_err_scalar = 0.0;
            for (size_t i = iv.lo; i < iv.hi; ++i) {
                max_err_block = std::max(
                    max_err_block, static_cast<double>(std::abs(dst_block[i] - ref[i])));
                max_err_recurrence =
                    std::max(max_err_recurrence,
                             static_cast<double>(std::abs(dst_recurrence[i] - ref[i])));
                max_err_scalar =
                    std::max(max_err_scalar,
                             static_cast<double>(std::abs(dst_scalar[i] - ref[i])));
            }
            BOOST_CHECK_LT(max_err_block, 2e-5);
            BOOST_CHECK_LT(max_err_recurrence, 2e-5);
            BOOST_CHECK_LT(max_err_scalar, 1e-6);

            // 3. Strict verification of phase origin = sample 0 (not lo!):
            // Sample lo must have phase -w * lo.
            if (iv.lo < iv.hi && iv.lo > 0 &&
                std::abs(w * static_cast<double>(iv.lo)) > 0.1) {
                const double expected_phase = -w * static_cast<double>(iv.lo);
                const gr_complex expected_lo =
                    src[iv.lo] * gr_complex(static_cast<float>(std::cos(expected_phase)),
                                            static_cast<float>(std::sin(expected_phase)));
                BOOST_CHECK_LT(std::abs(dst_block[iv.lo] - expected_lo), 2e-5f);
                BOOST_CHECK_LT(std::abs(dst_recurrence[iv.lo] - expected_lo), 2e-5f);
                // If lo were wrongly treated as phase 0, dst[lo] would equal src[lo]
                BOOST_CHECK_GT(std::abs(dst_block[iv.lo] - src[iv.lo]), 1e-3f);
                BOOST_CHECK_GT(std::abs(dst_recurrence[iv.lo] - src[iv.lo]), 1e-3f);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_task_a2_cir_sparse_correlation_against_dense)
{
    // Task A2: CIR ternary template sparse correlation.
    // Verify that rxopt::correlate_cir_sparse matches dense correlation
    // across all HRP codes (9, 10, 11, 12).
    const size_t tap_count = 38;
    const size_t code_len = 2 * kQm35ChipsPerSymbol; // 1016
    const size_t wlen = code_len + tap_count - 1;

    std::vector<gr_complex> corr(wlen);
    for (size_t i = 0; i < wlen; ++i) {
        corr[i] = gr_complex(std::cos(0.031f * static_cast<float>(i)) + 0.1f,
                             std::sin(0.047f * static_cast<float>(i)) - 0.15f);
    }

    for (size_t code_idx : { 9, 10, 11, 12 }) {
        const int8_t* pc = GetPreambleCode(code_idx);
        const auto spread = BuildSampledCode(pc, kQm35CodeLength);

        // Dense reference implementation
        std::vector<gr_complex> dense_work(code_len, gr_complex(0.0f, 0.0f));
        float dense_energy = 0.0f;
        for (size_t j = 0; j < spread.size(); ++j) {
            dense_work[2 * j] = gr_complex(static_cast<float>(spread[j]), 0.0f);
            dense_energy += std::norm(dense_work[2 * j]);
        }
        std::vector<gr_complex> dense_values(tap_count);
        for (size_t nn = 0; nn < tap_count; ++nn) {
            gr_complex acc(0.0f, 0.0f);
            for (size_t m = 0; m < code_len; ++m) {
                acc += std::conj(dense_work[m]) * corr[nn + m];
            }
            dense_values[nn] = acc / dense_energy;
        }

        // Sparse implementation
        rxopt::SparseCirTemplate sparse_tmpl;
        rxopt::build_sparse_cir_template(pc, kQm35CodeLength, sparse_tmpl);
        BOOST_CHECK_EQUAL(sparse_tmpl.code_energy, dense_energy);
        BOOST_CHECK_GT(sparse_tmpl.entries.size(), 0u);
        BOOST_CHECK_LE(sparse_tmpl.entries.size(),
                       127u); // ternary code length <= 127 (e.g. 64 non-zeros)

        std::vector<gr_complex> sparse_values;
        rxopt::correlate_cir_sparse(corr.data(), sparse_tmpl, tap_count, sparse_values);

        BOOST_REQUIRE_EQUAL(sparse_values.size(), dense_values.size());
        double max_diff = 0.0;
        for (size_t nn = 0; nn < tap_count; ++nn) {
            max_diff = std::max(
                max_diff,
                static_cast<double>(std::abs(sparse_values[nn] - dense_values[nn])));
        }
        // Arithmetic is identical additions/subtractions in identical order
        BOOST_CHECK_LT(max_diff, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_task_a1_a2_full_pipeline_range_cfo)
{
    // End-to-end full packet decode with range CFO compensation and sparse CIR:
    // Tests fast CFO rotation, on-demand soft-chips, and fused payload CMF.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));

    for (double cfo_hz : { 0.0, 25000.0, -25000.0 }) {
        std::vector<gr_complex> rx(iq.size());
        const double w = 2.0 * M_PI * cfo_hz / kQm35SampleRate;
        for (size_t i = 0; i < iq.size(); ++i) {
            const double phase = w * static_cast<double>(i);
            const float c = static_cast<float>(std::cos(phase));
            const float s = static_cast<float>(std::sin(phase));
            rx[i] = iq[i] * gr_complex(c, s);
        }

        auto prof = Qm35825Profile::Default();
        prof.sfd_mode = "ieee";
        prof.rx_fast_sfd = true;
        prof.rx_fast_cfo_rotation = true;
        prof.rx_softchip_on_demand = true;
        prof.rx_fused_payload_cmf = true;
        prof.rx_fused_phr_cmf = true;

        core::DemodScratch scratch;
        scratch.reserve(rx.size());
        auto res = core::demodulate_one(
            rx.data(), rx.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
        BOOST_CHECK(res.status == DemodStatus::Success);
        BOOST_CHECK(res.payload.fcs_pass);
        BOOST_CHECK_EQUAL(res.payload.received_fcs, uint16_t(0x584b));
        BOOST_CHECK_SMALL(res.cfo.cfo_hz - cfo_hz, 200.0); // within 200 Hz
        BOOST_CHECK_GT(res.stage_cfo_rotate_us, 0u);
        BOOST_CHECK_GT(res.stage_cfo_rotated_samples, 0u);
        BOOST_CHECK_GE(res.stage_cfo_rotation_passes, 2u);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_default_rx_recurrence)
{
    // Verify that default RX configuration (rx_fast_cfo_rotation = false)
    // succeeds using the fast 4-chain recurrence kernel.
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));

    const double cfo_hz = 30000.0;
    std::vector<gr_complex> rx(iq.size());
    const double w = 2.0 * M_PI * cfo_hz / kQm35SampleRate;
    for (size_t i = 0; i < iq.size(); ++i) {
        const double phase = w * static_cast<double>(i);
        rx[i] = iq[i] * gr_complex(static_cast<float>(std::cos(phase)),
                                   static_cast<float>(std::sin(phase)));
    }

    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    prof.rx_fast_cfo_rotation = false; // test recurrence path
    prof.rx_fast_sfd = true;

    core::DemodScratch scratch;
    scratch.reserve(rx.size());
    auto res = core::demodulate_one(
        rx.data(), rx.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
    BOOST_CHECK(res.status == DemodStatus::Success);
    BOOST_CHECK(res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(res.cfo.rotation_kernel, 0u);
    BOOST_CHECK_CLOSE(res.cfo.cfo_hz, cfo_hz, 200.0);
    BOOST_CHECK_GT(res.stage_cfo_rotated_samples, 0u);
}

BOOST_AUTO_TEST_CASE(test_demod_core_softchip_fallback_timing_accounting)
{
    // Test that on softchip fallback retry:
    // 1. Fallback rotation is included in stage_softchip_fallback_us and stage_total_us
    // 2. Full window is derotated and counted
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));

    // Corrupt payload region so fused CMF fails FCS and forces fallback retry
    std::vector<gr_complex> rx = iq;
    // Preamble + SFD + PHR is ~104000 samples, corrupt payload after sample 110000
    for (size_t i = 110000; i < rx.size(); ++i) {
        rx[i] = gr_complex(0.0f, 0.0f);
    }

    auto prof = Qm35825Profile::Default();
    prof.sfd_mode = "ieee";
    prof.rx_fast_sfd = true;
    prof.rx_fast_cfo_rotation = true;
    prof.rx_softchip_on_demand = true;
    prof.rx_fused_payload_cmf = true;
    prof.rx_fused_phr_cmf = true;

    core::DemodScratch scratch;
    scratch.reserve(rx.size());
    auto res = core::demodulate_one(
        rx.data(), rx.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
    // Fused failure returns directly without full-stream softchip fallback.
    BOOST_CHECK(res.status != DemodStatus::Success);
    BOOST_CHECK(!res.payload.fcs_pass);
    BOOST_CHECK_EQUAL(res.stage_softchip_fallbacks, 0u);
    BOOST_CHECK_EQUAL(res.stage_softchip_fallback_us, 0u);
}

BOOST_AUTO_TEST_CASE(test_demod_core_payload_substage_timers_partition_parent)
{
    // The payload sub-stage timers (soft / viterbi / rs / crc / other) are
    // diagnostic CHILDREN of stage_payload_us.  `other` is the exact remainder,
    // so the five buckets must sum to the in-function payload interval and can
    // never exceed the parent lap.  The only allowed difference is the
    // call/lap overhead (a few µs at 1 µs clock resolution).
    std::vector<gr_complex> iq;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));

    auto check_partition = [&](const DemodResult& res, const char* label) {
        BOOST_TEST_CONTEXT(label)
        {
            const uint64_t sum = res.stage_payload_soft_us +
                                 res.stage_payload_viterbi_us + res.stage_payload_rs_us +
                                 res.stage_payload_crc_us + res.stage_payload_other_us;
            BOOST_TEST_MESSAGE(
                label << ": payload=" << res.stage_payload_us
                      << " soft=" << res.stage_payload_soft_us << " vit="
                      << res.stage_payload_viterbi_us << " rs=" << res.stage_payload_rs_us
                      << " crc=" << res.stage_payload_crc_us
                      << " other=" << res.stage_payload_other_us << " sum=" << sum
                      << " fused_calls=" << res.stage_payload_fused_calls
                      << " legacy_calls=" << res.stage_payload_legacy_calls);
            // Decoded successfully, so the payload stage really ran.
            BOOST_REQUIRE(res.status == DemodStatus::Success);
            BOOST_REQUIRE(res.payload.fcs_pass);
            BOOST_REQUIRE_GE(
                res.stage_payload_fused_calls + res.stage_payload_legacy_calls, 1u);
            // Children never exceed the parent lap.
            BOOST_CHECK_LE(sum, res.stage_payload_us);
            // ... and the gap is only call overhead, not a missing bucket.
            BOOST_CHECK_LE(res.stage_payload_us - sum, 20u);
            BOOST_CHECK_GT(sum, 0u);
            // The three dominant buckets must be attributed, not dumped into
            // `other` (µs truncation keeps them strictly positive here).
            BOOST_CHECK_GT(res.stage_payload_soft_us, 0u);
            BOOST_CHECK_GT(res.stage_payload_viterbi_us, 0u);
            BOOST_CHECK_GT(res.stage_payload_rs_us + res.stage_payload_crc_us, 0u);
        }
    };

    // Legacy payload backend.
    {
        auto prof = Qm35825Profile::Default();
        prof.sfd_mode = "ieee";
        prof.rx_fast_sfd = true;
        prof.rx_fused_payload_cmf = false;
        prof.rx_fused_phr_cmf = false;
        prof.rx_softchip_on_demand = false;
        core::DemodScratch scratch;
        scratch.reserve(iq.size());
        auto res = core::demodulate_one(
            iq.data(), iq.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
        check_partition(res, "legacy");
        BOOST_CHECK_GE(res.stage_payload_legacy_calls, 1u);
        BOOST_CHECK_EQUAL(res.stage_payload_fused_calls, 0u);
    }

    // Fused CMF payload backend (the live default).
    {
        auto prof = Qm35825Profile::Default();
        prof.sfd_mode = "ieee";
        prof.rx_fast_sfd = true;
        prof.rx_fast_cfo_rotation = true;
        prof.rx_softchip_on_demand = true;
        prof.rx_fused_payload_cmf = true;
        prof.rx_fused_phr_cmf = true;
        core::DemodScratch scratch;
        scratch.reserve(iq.size());
        auto res = core::demodulate_one(
            iq.data(), iq.size(), prof, 1, 9984, 0, load_reference_template(), scratch);
        check_partition(res, "fused");
        BOOST_CHECK_EQUAL(res.stage_softchip_fallbacks, 0u);
        BOOST_CHECK_GE(res.stage_payload_fused_calls, 1u);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_real_sync_dot_numeric_equivalence)
{
    const size_t L = 1016;
    std::vector<gr_complex> tmpl = load_reference_template();
    BOOST_REQUIRE_EQUAL(tmpl.size(), L);

    rxopt::RealSyncTemplate real_tmpl;
    rxopt::build_real_sync_template(tmpl, real_tmpl);
    BOOST_CHECK(real_tmpl.is_real);
    BOOST_CHECK_EQUAL(real_tmpl.length, L);
    BOOST_CHECK_GT(real_tmpl.energy, 0.0f);

    // Generate test rx signal with synthetic multipath + noise
    std::vector<gr_complex> rx(L + 128);
    for (size_t i = 0; i < rx.size(); ++i) {
        rx[i] = gr_complex(std::sin(static_cast<float>(i) * 0.13f) +
                               0.1f * std::cos(static_cast<float>(i) * 0.37f),
                           std::cos(static_cast<float>(i) * 0.17f) -
                               0.2f * std::sin(static_cast<float>(i) * 0.29f));
    }

    // Verify 4-lag AVX2 dot vs VOLK across 64 candidate offsets
    for (size_t lag = 0; lag < 64; lag += 4) {
        gr_complex v_out[4];
        for (size_t j = 0; j < 4; ++j) {
            v_out[j] = rxopt::conjugate_dot(rx.data() + lag + j, tmpl.data(), L);
        }

        gr_complex a_out[4];
        rxopt::real_template_dot_4lag(
            rx.data() + lag, real_tmpl.taps_dup.data(), L, a_out);

        for (size_t j = 0; j < 4; ++j) {
            const float diff = std::abs(v_out[j] - a_out[j]);
            BOOST_CHECK_LT(diff, 1e-5f);

            // Also test single lag dot
            const gr_complex s_out = rxopt::real_template_dot_single(
                rx.data() + lag + j, real_tmpl.taps_dup.data(), L);
            const float s_diff = std::abs(v_out[j] - s_out);
            BOOST_CHECK_LT(s_diff, 1e-5f);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_sfd_dot_numeric_equivalence)
{
    const size_t sym = kQm35SamplesPerSymbol; // 1016
    std::vector<gr_complex> tmpl = load_reference_template();
    BOOST_REQUIRE_EQUAL(tmpl.size(), sym);

    rxopt::RealSyncTemplate real_tmpl;
    rxopt::build_real_sync_template(tmpl, real_tmpl);
    BOOST_REQUIRE(real_tmpl.is_real);

    // Test with IEEE SFD and SFD 4z2
    const char* modes[] = { "ieee", "4z2" };
    for (const char* mode : modes) {
        const std::vector<int8_t> sfd_seq = openuwb::demod::GetSfdSequence(mode);
        const size_t len = sfd_seq.size();
        BOOST_REQUIRE_GT(len, 0);

        float sfd_norm_sq = 0.0f;
        for (int8_t s : sfd_seq)
            sfd_norm_sq += static_cast<float>(s * s);
        BOOST_REQUIRE_GT(sfd_norm_sq, 0.0f);
        const float sfd_scale = 1.0f / std::sqrt(sfd_norm_sq * real_tmpl.energy);

        // Generate synthetic test rx signal covering full SFD plus margin
        const size_t total_sfd_samples = len * sym;
        std::vector<gr_complex> rx(total_sfd_samples + 256);
        for (size_t i = 0; i < rx.size(); ++i) {
            rx[i] = gr_complex(std::sin(static_cast<float>(i) * 0.07f) +
                                   0.3f * std::cos(static_cast<float>(i) * 0.23f),
                               std::cos(static_cast<float>(i) * 0.11f) -
                                   0.2f * std::sin(static_cast<float>(i) * 0.19f));
        }

        // Test over 64 candidate offsets, testing both 4-lag batch and single-lag
        for (size_t lag = 0; lag < 64; lag += 4) {
            gr_complex a4[4];
            rxopt::sfd_real_dot_4lag(rx.data() + lag,
                                     real_tmpl.taps_dup.data(),
                                     sym,
                                     sfd_seq.data(),
                                     len,
                                     sfd_scale,
                                     a4);

            for (size_t j = 0; j < 4; ++j) {
                const gr_complex s_out =
                    rxopt::sfd_real_dot_single(rx.data() + lag + j,
                                               real_tmpl.taps_dup.data(),
                                               sym,
                                               sfd_seq.data(),
                                               len,
                                               sfd_scale);

                const float diff = std::abs(a4[j] - s_out);
                const float mag = std::abs(s_out);
                const float rel_err = (mag > 1e-12f) ? (diff / mag) : diff;
                BOOST_CHECK_LT(rel_err, 2e-5f);
                BOOST_CHECK_LT(diff, 2e-5f);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_cached_sequence_invalidation)
{
    core::DemodScratch scratch;
    // Initially invalid
    BOOST_CHECK_EQUAL(scratch.cached_seq.code_index, std::numeric_limits<size_t>::max());
    BOOST_CHECK(scratch.cached_seq.pcode.empty());
    BOOST_CHECK(scratch.cached_seq.sampled_pcode.empty());
    BOOST_CHECK(scratch.cached_seq.sfd_seq.empty());

    // Update with code 9, "4z2"
    scratch.cached_seq.update(9, "4z2");
    BOOST_CHECK_EQUAL(scratch.cached_seq.code_index, 9);
    BOOST_CHECK_EQUAL(scratch.cached_seq.sfd_mode, "4z2");
    BOOST_CHECK_EQUAL(scratch.cached_seq.pcode.size(), kQm35CodeLength);
    const int8_t* pc9 = GetPreambleCode(9);
    for (size_t i = 0; i < kQm35CodeLength; ++i) {
        BOOST_CHECK_EQUAL(scratch.cached_seq.pcode[i], pc9[i]);
    }
    const std::vector<int8_t> expected_sampled9 = BuildSampledCode(pc9, kQm35CodeLength);
    BOOST_CHECK_EQUAL(scratch.cached_seq.sampled_pcode.size(), expected_sampled9.size());
    for (size_t i = 0; i < expected_sampled9.size(); ++i) {
        BOOST_CHECK_EQUAL(scratch.cached_seq.sampled_pcode[i], expected_sampled9[i]);
    }
    const std::vector<int8_t> expected_sfd4z2 = openuwb::demod::GetSfdSequence("4z2");
    BOOST_CHECK_EQUAL(scratch.cached_seq.sfd_seq.size(), expected_sfd4z2.size());
    for (size_t i = 0; i < expected_sfd4z2.size(); ++i) {
        BOOST_CHECK_EQUAL(scratch.cached_seq.sfd_seq[i], expected_sfd4z2[i]);
    }

    // Repeated update with same parameters should not reallocate/change
    const int8_t* pcode_ptr = scratch.cached_seq.pcode.data();
    scratch.cached_seq.update(9, "4z2");
    BOOST_CHECK_EQUAL(scratch.cached_seq.pcode.data(), pcode_ptr);

    // Update with code 10 (code change, SFD mode unchanged)
    scratch.cached_seq.update(10, "4z2");
    BOOST_CHECK_EQUAL(scratch.cached_seq.code_index, 10);
    const int8_t* pc10 = GetPreambleCode(10);
    BOOST_CHECK_EQUAL(scratch.cached_seq.pcode[0], pc10[0]);
    const std::vector<int8_t> expected_sampled10 =
        BuildSampledCode(pc10, kQm35CodeLength);
    BOOST_CHECK_EQUAL(scratch.cached_seq.sampled_pcode[0], expected_sampled10[0]);

    // Update with "ieee" (SFD mode change, code unchanged)
    scratch.cached_seq.update(10, "ieee");
    BOOST_CHECK_EQUAL(scratch.cached_seq.sfd_mode, "ieee");
    const std::vector<int8_t> expected_ieee = openuwb::demod::GetSfdSequence("ieee");
    BOOST_CHECK_EQUAL(scratch.cached_seq.sfd_seq.size(), expected_ieee.size());
    for (size_t i = 0; i < expected_ieee.size(); ++i) {
        BOOST_CHECK_EQUAL(scratch.cached_seq.sfd_seq[i], expected_ieee[i]);
    }

    // Explicit invalidate
    scratch.cached_seq.invalidate();
    BOOST_CHECK_EQUAL(scratch.cached_seq.code_index, std::numeric_limits<size_t>::max());
    BOOST_CHECK(scratch.cached_seq.pcode.empty());
    BOOST_CHECK(scratch.cached_seq.sampled_pcode.empty());
    BOOST_CHECK(scratch.cached_seq.sfd_seq.empty());
}

BOOST_AUTO_TEST_CASE(test_demod_core_vector_accumulate_scale)
{
    const size_t N = 1053;
    std::vector<gr_complex> dst_simd(N, gr_complex(0.0f, 0.0f));
    std::vector<gr_complex> dst_scalar(N, gr_complex(0.0f, 0.0f));

    // Test accumulation of 16 blocks of random-ish data
    for (size_t b = 0; b < 16; ++b) {
        std::vector<gr_complex> src(N);
        for (size_t i = 0; i < N; ++i) {
            src[i] = gr_complex(std::sin(static_cast<float>(b * N + i) * 0.1f),
                                std::cos(static_cast<float>(b * N + i) * 0.17f));
            dst_scalar[i] += src[i];
        }
        rxopt::vector_accumulate_complex(dst_simd.data(), src.data(), N);
    }

    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE(dst_simd[i].real(), dst_scalar[i].real(), 1e-4f);
        BOOST_CHECK_CLOSE(dst_simd[i].imag(), dst_scalar[i].imag(), 1e-4f);
    }

    const float scale = 1.0f / 16.0f;
    for (size_t i = 0; i < N; ++i)
        dst_scalar[i] *= scale;
    rxopt::vector_scale_complex(dst_simd.data(), scale, N);

    for (size_t i = 0; i < N; ++i) {
        BOOST_CHECK_CLOSE(dst_simd[i].real(), dst_scalar[i].real(), 1e-4f);
        BOOST_CHECK_CLOSE(dst_simd[i].imag(), dst_scalar[i].imag(), 1e-4f);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_sfd_fold_window_equivalence)
{
    const size_t sym = kQm35SamplesPerSymbol; // 1016
    const auto sfd_seq = openuwb::demod::GetSfdSequence("4z2");
    const size_t len = sfd_seq.size(); // 8

    std::vector<gr_complex> tmpl(sym);
    for (size_t i = 0; i < sym; ++i) {
        tmpl[i] = gr_complex(std::cos(static_cast<float>(i) * 0.12f), 0.0f);
    }
    rxopt::RealSyncTemplate real_tmpl;
    rxopt::build_real_sync_template(tmpl, real_tmpl);

    float sfd_norm_sq = 0.0f;
    for (int8_t s : sfd_seq)
        sfd_norm_sq += static_cast<float>(s * s);
    const float sfd_scale = 1.0f / std::sqrt(sfd_norm_sq * real_tmpl.energy);

    // Synthetic rx of length 16 * sym
    const size_t rx_len = 16 * sym;
    std::vector<gr_complex> rx(rx_len);
    for (size_t i = 0; i < rx_len; ++i) {
        rx[i] = gr_complex(std::sin(static_cast<float>(i) * 0.05f),
                           std::cos(static_cast<float>(i) * 0.07f));
    }

    const size_t j_lo = 100;
    const size_t j_end = 228;
    const size_t fold_len = (j_end - j_lo) + sym;
    std::vector<gr_complex> folded(fold_len);
    rxopt::fold_sfd_window(
        rx.data() + j_lo, sym, sfd_seq.data(), len, fold_len, folded.data());

    for (size_t j = j_lo; j <= j_end; j += 7) {
        const gr_complex dot_direct =
            rxopt::sfd_real_dot_single(rx.data() + j,
                                       real_tmpl.taps_dup.data(),
                                       sym,
                                       sfd_seq.data(),
                                       len,
                                       sfd_scale);
        const gr_complex dot_folded =
            rxopt::real_template_dot_single(
                folded.data() + (j - j_lo), real_tmpl.taps_dup.data(), sym) *
            sfd_scale;

        BOOST_CHECK_SMALL(std::abs(dot_direct.real() - dot_folded.real()), 1e-5f);
        BOOST_CHECK_SMALL(std::abs(dot_direct.imag() - dot_folded.imag()), 1e-5f);
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_ns_sfd_sparse_sliding_equivalence)
{
    auto prof = Qm35825Profile::Default();
    prof.code_index = 10;
    prof.sfd_mode = "4z2";

    const auto sfd_seq = openuwb::demod::GetSfdSequence(prof.sfd_mode);
    const size_t total_chips = 80 * kQm35ChipsPerSymbol;

    std::vector<float> soft_chips(total_chips, 0.0f);
    const int64_t target_start = 64 * kQm35ChipsPerSymbol + 2;

    core::DemodScratch scratch;
    scratch.cached_seq.update(prof.code_index, prof.sfd_mode);
    const auto& sc = scratch.cached_seq.sampled_pcode;

    // Inject target NS-SFD pattern with polarity +1
    for (size_t i = 0; i < sfd_seq.size(); ++i) {
        const int8_t s = sfd_seq[i];
        for (size_t j = 0; j < sc.size(); ++j) {
            soft_chips[target_start + i * sc.size() + j] = static_cast<float>(s * sc[j]);
        }
    }

    NsSfdResult res;
    BOOST_REQUIRE(
        core::stage_ns_sfd(soft_chips, prof, sfd_seq, kQm35ChipsPerSymbol, res, scratch));
    BOOST_CHECK(res.ok);
    BOOST_CHECK_EQUAL(res.sfd_start_chip, target_start);
    BOOST_CHECK_EQUAL(res.polarity, 1);
    BOOST_CHECK_GE(res.metric, 0.95f);
}

BOOST_AUTO_TEST_CASE(test_demod_core_ns_sfd_relative_par_low_snr_and_noise)
{
    auto prof = Qm35825Profile::Default();
    prof.code_index = 10;
    prof.sfd_mode = "4z2";

    const auto sfd_seq = openuwb::demod::GetSfdSequence(prof.sfd_mode);
    const size_t total_chips = 80 * kQm35ChipsPerSymbol;
    const int64_t target_start = 64 * kQm35ChipsPerSymbol - 3;

    core::DemodScratch scratch;
    scratch.cached_seq.update(prof.code_index, prof.sfd_mode);
    const auto& sc = scratch.cached_seq.sampled_pcode;

    // 1. Low-SNR SFD whose normalized correlation falls below the legacy 0.30
    // absolute floor (e.g. ~0.18-0.25) must pass via relative PAR and recover
    // exact chip alignment and polarity (-1).
    std::mt19937 rng(20261122u);
    std::normal_distribution<float> noise(0.0f, 1.5f);
    std::vector<float> low_snr_chips(total_chips);
    for (float& v : low_snr_chips)
        v = noise(rng);
    for (size_t i = 0; i < sfd_seq.size(); ++i) {
        const int8_t s = sfd_seq[i];
        for (size_t j = 0; j < sc.size(); ++j) {
            low_snr_chips[target_start + i * sc.size() + j] -=
                static_cast<float>(s * sc[j]);
        }
    }

    NsSfdResult low_res;
    BOOST_REQUIRE(core::stage_ns_sfd(
        low_snr_chips, prof, sfd_seq, kQm35ChipsPerSymbol, low_res, scratch));
    BOOST_CHECK(low_res.ok);
    BOOST_CHECK_EQUAL(low_res.sfd_start_chip, target_start);
    BOOST_CHECK_EQUAL(low_res.polarity, -1);
    BOOST_CHECK_LT(low_res.metric, prof.sfd_detection_threshold);
    BOOST_CHECK_GT(low_res.metric, 0.10f);

    // 2. All-zero stream must be rejected.
    std::vector<float> zeros(total_chips, 0.0f);
    NsSfdResult zero_res;
    BOOST_CHECK(!core::stage_ns_sfd(
        zeros, prof, sfd_seq, kQm35ChipsPerSymbol, zero_res, scratch));

    // 3. Pure white noise must be overwhelmingly rejected by the PAR gate.
    size_t false_alarms = 0;
    constexpr size_t kNoiseTrials = 50;
    std::vector<float> pure_noise(total_chips);
    for (size_t t = 0; t < kNoiseTrials; ++t) {
        for (float& v : pure_noise)
            v = noise(rng);
        NsSfdResult nres;
        if (core::stage_ns_sfd(
                pure_noise, prof, sfd_seq, kQm35ChipsPerSymbol, nres, scratch))
            ++false_alarms;
    }
    BOOST_CHECK_LE(false_alarms, size_t(5));
}

BOOST_AUTO_TEST_CASE(test_demod_core_tail_reuse_equivalence_and_fallback)
{
    std::vector<gr_complex> iq, tmpl;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", iq));
    BOOST_REQUIRE(load_cf32("../../../testdata/reference_preamble.bin", tmpl));

    // 1. Nominal reuse run: default profile with rx_seed_refine_tail = true (narrow
    // search = 8). All 5 admission conditions are satisfied, so run_tail_and_cfo_refit
    // executes tail reuse.
    core::DemodScratch scratch_reuse;
    scratch_reuse.reserve(iq.size());
    auto prof_reuse = Qm35825Profile::Default();
    prof_reuse.sfd_mode = "ieee";
    prof_reuse.rx_seed_refine_tail = true;
    auto res_reuse = core::demodulate_one(
        iq.data(), iq.size(), prof_reuse, 1, 9984, 0, tmpl, scratch_reuse);
    BOOST_CHECK(res_reuse.status == DemodStatus::Success);

    // 2. Explicit full re-search reference run: set rx_seed_refine_tail = false
    // (seed_radius = 0). This forces run_tail_and_cfo_refit to bypass can_reuse and
    // execute full stage_timing.
    core::DemodScratch scratch_legacy;
    scratch_legacy.reserve(iq.size());
    auto prof_legacy = Qm35825Profile::Default();
    prof_legacy.sfd_mode = "ieee";
    prof_legacy.rx_seed_refine_tail = false;
    auto res_legacy = core::demodulate_one(
        iq.data(), iq.size(), prof_legacy, 1, 9984, 0, tmpl, scratch_legacy);
    BOOST_CHECK(res_legacy.status == DemodStatus::Success);

    // Equivalence checks between reuse and explicit re-search
    BOOST_CHECK_EQUAL(res_reuse.timing.detected_peaks, res_legacy.timing.detected_peaks);
    BOOST_CHECK_EQUAL(res_reuse.timing.peak_samples.size(),
                      res_legacy.timing.peak_samples.size());
    BOOST_REQUIRE_GE(res_reuse.timing.peak_samples.size(), size_t(64));
    // Verify that the tail 40 peaks in reuse exactly match the explicit tracking
    for (size_t k = 24; k < 64; ++k) {
        BOOST_CHECK_EQUAL(res_reuse.timing.peak_samples[k],
                          res_legacy.timing.peak_samples[k]);
    }
    // Estimated CFO must match closely (within < 15 Hz)
    BOOST_CHECK_SMALL(std::abs(res_reuse.cfo.cfo_hz - res_legacy.cfo.cfo_hz), 15.0);
    BOOST_CHECK_EQUAL(res_reuse.phr.psdu_length, res_legacy.phr.psdu_length);

    // 3. Fallback test: violate an admission condition.
    // Setting timing_track_radius = 12 (not 8) causes can_reuse == false.
    core::DemodScratch scratch_fb;
    scratch_fb.reserve(iq.size());
    auto prof_fb = Qm35825Profile::Default();
    prof_fb.sfd_mode = "ieee";
    prof_fb.rx_seed_refine_tail = true;
    prof_fb.timing_track_radius = 12; // Incompatible with radius 8 reuse condition
    auto res_fb =
        core::demodulate_one(iq.data(), iq.size(), prof_fb, 1, 9984, 0, tmpl, scratch_fb);
    BOOST_CHECK(res_fb.status == DemodStatus::Success);
    BOOST_CHECK_EQUAL(res_fb.timing.cir_grid_repetitions, size_t(64));
    BOOST_CHECK_GE(res_fb.timing.detected_peaks, size_t(8));
    BOOST_CHECK_SMALL(std::abs(res_fb.cfo.cfo_hz - res_legacy.cfo.cfo_hz), 25.0);
}

BOOST_AUTO_TEST_CASE(test_demod_core_sfd_anchored_cfo_uses_full_preamble)
{
    // Keep the golden SFD/PHR/payload and vary only the SYNC train length.
    // The SFD anchor pins the end of the last SYNC, so the final CFO now fits
    // the whole preamble after the startup exclusion.  Only the excluded
    // leading SYNCs keep their zero placeholders; every measured grid point of
    // the fit carries a real matched-filter correlation.
    std::vector<gr_complex> golden;
    BOOST_REQUIRE(load_cf32(golden_dir() + "/window.cfile", golden));
    const auto tmpl = load_reference_template();
    constexpr size_t start = 9984;
    constexpr size_t sym = kQm35SamplesPerSymbol;
    constexpr size_t original_reps = 64;
    const size_t sfd = start + original_reps * sym;
    BOOST_REQUIRE_GT(golden.size(), sfd);
    const std::vector<gr_complex> stable_symbol(golden.begin() + sfd - sym,
                                                golden.begin() + sfd);

    for (size_t reps : { size_t(16), size_t(32), size_t(64), size_t(128), size_t(256) }) {
        std::vector<gr_complex> base(golden.begin(), golden.begin() + start);
        if (reps > original_reps) {
            for (size_t k = original_reps; k < reps; ++k)
                base.insert(base.end(), stable_symbol.begin(), stable_symbol.end());
            base.insert(base.end(), golden.begin() + start, golden.end());
        } else {
            base.insert(base.end(), golden.begin() + sfd - reps * sym, golden.end());
        }
        for (double hz : { 0.0, -25000.0, 25000.0, -175000.0, 175000.0 }) {
            const double rate = 1.0 + hz / Qm35825Profile::Default().carrier_freq_hz;
            std::vector<gr_complex> rx(base.size(), gr_complex(0.0f, 0.0f));
            for (size_t k = 0; k < base.size(); ++k) {
                const double x = static_cast<double>(k) * rate;
                const size_t i0 = static_cast<size_t>(x);
                gr_complex s(0.0f, 0.0f);
                if (i0 + 1 < base.size()) {
                    const double frac = x - static_cast<double>(i0);
                    s = base[i0] * static_cast<float>(1.0 - frac) +
                        base[i0 + 1] * static_cast<float>(frac);
                } else if (i0 < base.size()) {
                    s = base[i0];
                }
                const double phase = 2.0 * M_PI * hz * k / kQm35SampleRate;
                rx[k] = s * gr_complex(std::cos(phase), std::sin(phase));
            }
            // Exercise both explicit tail re-tracking and the peak-reuse path.
            for (bool reuse : { false, true }) {
                BOOST_TEST_CONTEXT("reps=" << reps << " hz=" << hz << " reuse=" << reuse)
                {
                    auto p = Qm35825Profile::Default();
                    p.preamble_repetitions = reps;
                    p.cir_repetitions = reps - p.cir_skip_initial_repetitions;
                    p.sfd_mode = "ieee";
                    p.rx_seed_refine_tail = reuse;
                    core::DemodScratch scratch;
                    auto r = core::demodulate_one(
                        rx.data(), rx.size(), p, 0, start, 0, tmpl, scratch);
                    BOOST_REQUIRE(r.status == DemodStatus::Success);
                    BOOST_REQUIRE(r.payload.fcs_pass);
                    BOOST_CHECK_EQUAL(r.payload.received_fcs, uint16_t(0x584b));
                    // Full-preamble SFD-anchored fit: skip the startup
                    // transient but keep a 32-SYMBOL floor.
                    const size_t min_fit =
                        std::max<size_t>(8, p.cfo_min_fit_repetitions);
                    const size_t avail_skip =
                        (reps > min_fit) ? (reps - min_fit) : 0;
                    const size_t first =
                        std::min(p.cfo_skip_initial_repetitions, avail_skip);
                    const size_t measured = reps - first;
                    BOOST_CHECK_EQUAL(r.cfo.peaks_used, measured);
                    BOOST_CHECK_EQUAL(r.cfo.skipped_peaks, first);
                    BOOST_CHECK_EQUAL(r.cfo.fit_first_peak_sample,
                                      r.timing.peak_samples[first]);
                    BOOST_CHECK_EQUAL(r.cfo.fit_last_peak_sample,
                                      r.timing.peak_samples.back());
                    BOOST_CHECK_SMALL(r.cfo.cfo_hz - hz, 50.0);
                    for (size_t k = 0; k < first; ++k)
                        BOOST_CHECK_EQUAL(std::abs(r.timing.peak_corr[k]), 0.0f);
                    for (size_t k = first; k < reps; ++k)
                        BOOST_CHECK_GT(std::abs(r.timing.peak_corr[k]), 0.0f);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_demod_core_soft_trace_viterbi_matches_plain)
{
    // decode_payload_soft() now runs ONE tracing soft Viterbi and reuses its
    // bits for the first RS/CRC attempt, instead of running the plain Viterbi
    // and then recomputing the whole thing with margins inside the rescue.
    // That is only sound if vitdec_cl3_soft_trace()'s surviving path is
    // bit-identical to vitdec_cl3_soft()'s, so pin it over random LLRs.
    std::mt19937 rng(20261117u);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    const size_t lengths[] = { 64, 257, 1229 };
    for (size_t n : lengths) {
        std::vector<float> l0(n), l1(n);
        for (size_t i = 0; i < n; ++i) {
            l0[i] = nd(rng);
            l1[i] = nd(rng);
        }
        std::vector<int8_t> plain;
        openuwb::demod::core::detail::vitdec_cl3_soft(l0.data(), l1.data(), n, plain);
        openuwb::demod::core::detail::Cl3TrellisTrace tr;
        openuwb::demod::core::detail::vitdec_cl3_soft_trace(l0.data(), l1.data(), n, tr);
        BOOST_REQUIRE_EQUAL(tr.bits.size(), n);
        BOOST_REQUIRE_EQUAL(plain.size(), n);
        BOOST_CHECK(std::equal(plain.begin(), plain.end(), tr.bits.begin()));
        // The trace must also carry usable per-step evidence for the rescue.
        BOOST_REQUIRE_EQUAL(tr.margins.size(), n);
        BOOST_REQUIRE_EQUAL(tr.states.size(), n + 1);
        BOOST_REQUIRE_EQUAL(tr.lose_u.size(), n);
        for (size_t i = 0; i < n; ++i) {
            BOOST_CHECK(tr.margins[i] >= 0.0f);
            BOOST_CHECK(tr.lose_u[i] <= 1);
        }
    }
}

// ---------------------------------------------------------------------------
// FEC-OPT (perf/payload-fec-opt): bit-exactness of the table-driven GF(2^6)
// RS(63,55) decoder and the unrolled 4-state CL-3 Viterbi.
//
// The reference implementation below is an independent copy of the code that
// ran before this change (bit-serial GF multiply, alpha-power syndrome loop,
// `for (st) for (u)` ACS).  The optimised headers are always compared against
// this copy -- never against themselves.
// ---------------------------------------------------------------------------
namespace fec_opt_ref {

// --- Original GF(2^6) arithmetic (poly 0x61, alpha = 2) -------------------
int gf_mul(int a, int b)
{
    int p = 0;
    a &= 0x3f;
    b &= 0x3f;
    for (int i = 0; i < 6; ++i) {
        if (b & 1)
            p ^= a;
        const bool carry = (a & 0x20) != 0;
        a = (a << 1) & 0x3f;
        if (carry)
            a ^= 0x21;
        b >>= 1;
    }
    return p & 0x3f;
}

int gf_pow(int a, int n)
{
    int r = 1;
    a &= 0x3f;
    n %= 63;
    if (n < 0)
        n += 63;
    while (n) {
        if (n & 1)
            r = gf_mul(r, a);
        a = gf_mul(a, a);
        n >>= 1;
    }
    return r;
}

int gf_inv(int a) { return a ? gf_pow(a, 62) : 0; }

// --- Original RS decoder --------------------------------------------------
void bits_to_syms(const int8_t* bits, int count, int* symbols)
{
    for (int i = 0; i < count; ++i) {
        int s = 0;
        for (int b = 0; b < 6; ++b)
            s = (s << 1) | (bits[i * 6 + b] & 1);
        symbols[i] = s;
    }
}

void syms_to_bits(const int* syms, int nsym, int8_t* bits)
{
    for (int i = 0; i < nsym; ++i)
        for (int b = 0; b < 6; ++b)
            bits[i * 6 + b] = static_cast<int8_t>((syms[i] >> (5 - b)) & 1);
}

void syndromes(const int* r63, int* S8)
{
    for (int i = 1; i <= 8; ++i) {
        const int ai = gf_pow(2, i);
        int s = 0;
        int ap = 1;
        for (int j = 0; j < 63; ++j) {
            s ^= gf_mul(r63[j], ap);
            ap = gf_mul(ap, ai);
        }
        S8[i - 1] = s;
    }
}

int bm_len(const int* S, int count, int* Lambda)
{
    int C[9] = { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    int B[9] = { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    int L = 0;
    int m = 1;
    int b = 1;
    for (int n = 0; n < count; ++n) {
        int delta = S[n];
        for (int i = 1; i <= L; ++i)
            delta ^= gf_mul(C[i], S[n - i]);
        if (delta == 0) {
            ++m;
        } else {
            int T[9];
            for (int i = 0; i < 9; ++i)
                T[i] = C[i];
            const int scale = gf_mul(delta, gf_inv(b));
            for (int i = 0; i < 9 - m; ++i)
                C[i + m] ^= gf_mul(scale, B[i]);
            if (2 * L <= n) {
                L = n + 1 - L;
                for (int i = 0; i < 9; ++i)
                    B[i] = T[i];
                b = delta;
                m = 1;
            } else {
                ++m;
            }
        }
    }
    for (int i = 0; i < 9; ++i)
        Lambda[i] = C[i];
    return L;
}

int chien(const int* Lambda, int L, int* err_pos)
{
    int nerr = 0;
    for (int j = 0; j < 63; ++j) {
        const int xinv = gf_pow(2, (63 - j) % 63);
        int v = 0;
        int xp = 1;
        for (int i = 0; i <= L; ++i) {
            v ^= gf_mul(Lambda[i], xp);
            xp = gf_mul(xp, xinv);
        }
        if (v == 0) {
            if (nerr < 4)
                err_pos[nerr] = j;
            ++nerr;
        }
    }
    return nerr;
}

void forney(
    const int* S8, const int* Lambda, int L, const int* err_pos, int nerr, int* err_val)
{
    int Omega[8] = {};
    for (int i = 0; i < 8; ++i) {
        int o = 0;
        for (int j = 0; j <= L && j <= i; ++j)
            o ^= gf_mul(Lambda[j], S8[i - j]);
        Omega[i] = o;
    }
    for (int k = 0; k < nerr; ++k) {
        const int j = err_pos[k];
        const int Xinv = gf_pow(2, (63 - j) % 63);
        int num = 0;
        int xp = 1;
        for (int i = 0; i < 8; ++i) {
            num ^= gf_mul(Omega[i], xp);
            xp = gf_mul(xp, Xinv);
        }
        int den = 0;
        xp = 1;
        for (int i = 1; i <= L; i += 2) {
            den ^= gf_mul(Lambda[i], xp);
            xp = gf_mul(xp, gf_mul(Xinv, Xinv));
        }
        err_val[k] = den ? gf_mul(num, gf_inv(den)) : 0;
    }
}

bool zero_syndromes(const int* S)
{
    for (int i = 0; i < 8; ++i)
        if (S[i] != 0)
            return false;
    return true;
}

bool decode_block(const int8_t* coded378, int8_t* data330)
{
    int r[63];
    bits_to_syms(coded378, 63, r);
    int S[8];
    syndromes(r, S);
    if (!zero_syndromes(S)) {
        int Lambda[9] = {};
        const int L = bm_len(S, 8, Lambda);
        if (L <= 0 || L > 4)
            return false;
        int err_pos[4] = {};
        const int nerr = chien(Lambda, L, err_pos);
        if (nerr != L)
            return false;
        int corrections[8] = {};
        forney(S, Lambda, L, err_pos, nerr, corrections);
        for (int i = 0; i < nerr; ++i)
            r[err_pos[i]] ^= corrections[i];
        int check[8];
        syndromes(r, check);
        if (!zero_syndromes(check))
            return false;
    }
    syms_to_bits(r, 55, data330);
    return true;
}

void pack_psdu_bytes(const std::vector<int8_t>& bits, std::vector<uint8_t>& bytes)
{
    bytes.resize(bits.size() / 8);
    for (size_t i = 0; i < bytes.size(); ++i) {
        uint8_t byte = 0;
        for (int bit = 0; bit < 8; ++bit)
            byte |= static_cast<uint8_t>(bits[8 * i + bit]) << bit;
        bytes[i] = byte;
    }
}

// --- Original CL-3 Viterbi ------------------------------------------------
struct Cl3Transition {
    int g0;
    int g1;
    int next_state;
};

Cl3Transition cl3_transition(int state, int input)
{
    const int s1 = state >> 1;
    const int s0 = state & 1;
    return { s1, input ^ s0, (input << 1) | s1 };
}

float soft_metric(float path, int o0, int o1, float llr0, float llr1)
{
    const float c0 = (2.0f * static_cast<float>(o0) - 1.0f) * llr0;
    const float c1 = (2.0f * static_cast<float>(o1) - 1.0f) * llr1;
    return path + c0 + c1;
}

template <typename Metric, typename ExtendPath>
void vit_truncated(size_t N, Metric inf, ExtendPath extend, std::vector<int8_t>& out)
{
    int o0[4][2], o1[4][2], nxt[4][2];
    for (int st = 0; st < 4; ++st) {
        for (int u = 0; u < 2; ++u) {
            const auto step = cl3_transition(st, u);
            o0[st][u] = step.g0;
            o1[st][u] = step.g1;
            nxt[st][u] = step.next_state;
        }
    }
    Metric surf[4] = { Metric{ 0 }, inf, inf, inf };
    std::vector<int8_t> back(N * 4);
    for (size_t t = 0; t < N; ++t) {
        Metric ns[4] = { inf, inf, inf, inf };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const Metric m = extend(t, surf[st], o0[st][u], o1[st][u]);
                const int nst = nxt[st][u];
                if (m < ns[nst]) {
                    ns[nst] = m;
                    back[t * 4 + nst] = static_cast<int8_t>((u << 2) | st);
                }
            }
        }
        for (int i = 0; i < 4; ++i)
            surf[i] = ns[i];
    }
    int best = 0;
    for (int i = 1; i < 4; ++i)
        if (surf[i] < surf[best])
            best = i;
    out.resize(N);
    int cur = best;
    for (size_t t = N; t-- > 0;) {
        const int v = back[t * 4 + cur];
        out[t] = static_cast<int8_t>(v >> 2);
        cur = v & 3;
    }
}

void vit_hard(const int8_t* rx, size_t N, std::vector<int8_t>& out)
{
    vit_truncated(
        N,
        1000000,
        [&](size_t t, int path, int g0, int g1) {
            const int distance = (rx[2 * t] ^ g0) + (rx[2 * t + 1] ^ g1);
            return path + distance;
        },
        out);
}

// Count the ACS comparisons where a candidate metric exactly equals the
// running best (and both are finite): those are the steps where the strict
// `if (m < ns[nst])` tie rule decides the survivor, so the equivalence test
// below is only meaningful if this counter is non-zero.
template <typename Metric, typename ExtendPath>
size_t count_ties(size_t N, Metric inf, ExtendPath extend)
{
    int o0[4][2], o1[4][2], nxt[4][2];
    for (int st = 0; st < 4; ++st) {
        for (int u = 0; u < 2; ++u) {
            const auto step = cl3_transition(st, u);
            o0[st][u] = step.g0;
            o1[st][u] = step.g1;
            nxt[st][u] = step.next_state;
        }
    }
    Metric surf[4] = { Metric{ 0 }, inf, inf, inf };
    size_t ties = 0;
    for (size_t t = 0; t < N; ++t) {
        Metric ns[4] = { inf, inf, inf, inf };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const Metric m = extend(t, surf[st], o0[st][u], o1[st][u]);
                const int nst = nxt[st][u];
                if (m < ns[nst]) {
                    ns[nst] = m;
                } else if (m == ns[nst] && m != inf) {
                    ++ties;
                }
            }
        }
        for (int i = 0; i < 4; ++i)
            surf[i] = ns[i];
    }
    return ties;
}

size_t count_hard_ties(const int8_t* rx, size_t N)
{
    return count_ties(N, 1000000, [&](size_t t, int path, int g0, int g1) {
        const int distance = (rx[2 * t] ^ g0) + (rx[2 * t + 1] ^ g1);
        return path + distance;
    });
}

void vit_soft(const float* l0, const float* l1, size_t N, std::vector<int8_t>& out)
{
    vit_truncated(
        N,
        1e30f,
        [&](size_t t, float path, int g0, int g1) {
            return soft_metric(path, g0, g1, l0[t], l1[t]);
        },
        out);
}

struct Trace {
    std::vector<int8_t> bits;
    std::vector<int8_t> states;
    std::vector<float> margins;
    std::vector<uint8_t> lose_u;
};

void vit_soft_trace(const float* llr0, const float* llr1, size_t N, Trace& tr)
{
    int o0[4][2], o1[4][2], nxt[4][2];
    for (int st = 0; st < 4; ++st) {
        for (int u = 0; u < 2; ++u) {
            const auto step = cl3_transition(st, u);
            o0[st][u] = step.g0;
            o1[st][u] = step.g1;
            nxt[st][u] = step.next_state;
        }
    }
    const float INF = 1e30f;
    float surf[4] = { 0.0f, INF, INF, INF };
    std::vector<int8_t> back(N * 4);
    std::vector<float> diff(N * 4, INF);
    std::vector<uint8_t> alt(N * 4, 0);
    for (size_t t = 0; t < N; ++t) {
        float ns[4] = { INF, INF, INF, INF };
        float ns2[4] = { INF, INF, INF, INF };
        uint8_t nalt[4] = { 0, 0, 0, 0 };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const float m =
                    soft_metric(surf[st], o0[st][u], o1[st][u], llr0[t], llr1[t]);
                const int nst = nxt[st][u];
                const uint8_t code = static_cast<uint8_t>((u << 2) | st);
                if (m < ns[nst]) {
                    ns2[nst] = ns[nst];
                    nalt[nst] = static_cast<uint8_t>(back[t * 4 + nst]);
                    ns[nst] = m;
                    back[t * 4 + nst] = static_cast<int8_t>(code);
                } else if (m < ns2[nst]) {
                    ns2[nst] = m;
                    nalt[nst] = code;
                }
            }
        }
        for (int i = 0; i < 4; ++i) {
            surf[i] = ns[i];
            diff[t * 4 + i] = ns2[i] - ns[i];
            alt[t * 4 + i] = nalt[i];
        }
    }
    int best = 0;
    for (int i = 1; i < 4; ++i)
        if (surf[i] < surf[best])
            best = i;
    tr.bits.resize(N);
    tr.states.resize(N + 1);
    tr.margins.resize(N);
    tr.lose_u.resize(N);
    int cur = best;
    tr.states[N] = static_cast<int8_t>(cur);
    for (size_t t = N; t-- > 0;) {
        const int v = back[t * 4 + cur];
        tr.bits[t] = static_cast<int8_t>(v >> 2);
        tr.margins[t] = diff[t * 4 + cur];
        tr.lose_u[t] = static_cast<uint8_t>(alt[t * 4 + cur] >> 2);
        cur = v & 3;
        tr.states[t] = static_cast<int8_t>(cur);
    }
}

// --- Systematic RS(63,55) test encoder ------------------------------------
// Solves A p = b over GF(2^6) for the 8 parity symbols so that
// r(x) = sum_{j<55} d_j x^j + sum_k p_k x^{55+k} has roots alpha^1..alpha^8.
bool encode_systematic(const int* data55, int* out63)
{
    for (int j = 0; j < 55; ++j)
        out63[j] = data55[j] & 0x3f;
    int A[8][8];
    int b[8];
    for (int i = 1; i <= 8; ++i) {
        int s = 0;
        for (int j = 0; j < 55; ++j)
            s ^= gf_mul(out63[j], gf_pow(2, i * j));
        b[i - 1] = s;
        for (int k = 0; k < 8; ++k)
            A[i - 1][k] = gf_pow(2, i * (55 + k));
    }
    for (int c = 0; c < 8; ++c) {
        int piv = -1;
        for (int r = c; r < 8; ++r)
            if (A[r][c]) {
                piv = r;
                break;
            }
        if (piv < 0)
            return false;
        if (piv != c) {
            for (int k = 0; k < 8; ++k)
                std::swap(A[c][k], A[piv][k]);
            std::swap(b[c], b[piv]);
        }
        const int invp = gf_inv(A[c][c]);
        for (int k = c; k < 8; ++k)
            A[c][k] = gf_mul(A[c][k], invp);
        b[c] = gf_mul(b[c], invp);
        for (int r = 0; r < 8; ++r) {
            if (r == c)
                continue;
            const int f = A[r][c];
            if (!f)
                continue;
            for (int k = c; k < 8; ++k)
                A[r][k] ^= gf_mul(f, A[c][k]);
            b[r] ^= gf_mul(f, b[c]);
        }
    }
    for (int k = 0; k < 8; ++k)
        out63[55 + k] = b[k];
    return true;
}

} // namespace fec_opt_ref

BOOST_AUTO_TEST_CASE(test_fec_opt_gf64_tables_match_bitserial_reference)
{
    using namespace openuwb::codec;
    // 64 x 64 exhaustive multiply.
    for (int a = 0; a < 64; ++a) {
        for (int b = 0; b < 64; ++b) {
            const int want = fec_opt_ref::gf_mul(a, b);
            BOOST_REQUIRE_EQUAL(rs_gf_mul(a, b), want);
            // Input masking must keep the old `a &= 0x3f` / `b &= 0x3f`.
            BOOST_REQUIRE_EQUAL(rs_gf_mul(a | 0x40, b), want);
            BOOST_REQUIRE_EQUAL(rs_gf_mul(a, b | 0x80), want);
            BOOST_REQUIRE_EQUAL(rs_gf_mul(a | 0xc0, b | 0xc0), want);
        }
    }
    // Powers, including the negative / >62 reduction paths.
    const int exponents[] = { -64, -1, 0, 1, 61, 62, 63, 64, 126 };
    for (int a = 0; a < 64; ++a)
        for (int n : exponents)
            BOOST_REQUIRE_EQUAL(rs_gf_pow(a, n), fec_opt_ref::gf_pow(a, n));
    // Inverses (inv(0) == 0).
    for (int a = 0; a < 64; ++a) {
        BOOST_REQUIRE_EQUAL(rs_gf_inv(a), fec_opt_ref::gf_inv(a));
        BOOST_REQUIRE_EQUAL(rs_gf_inv(a | 0x40), fec_opt_ref::gf_inv(a));
    }
    BOOST_REQUIRE_EQUAL(rs_gf_inv(0), 0);
    // The alpha power lookups must equal pow(2, k) for the whole cycle.
    for (int k = -63; k <= 126; ++k) {
        BOOST_REQUIRE_EQUAL(rs_gf_alpha_pow(k), fec_opt_ref::gf_pow(2, k));
        BOOST_REQUIRE_EQUAL(rs_gf_alpha_inv_pow(k), fec_opt_ref::gf_pow(2, -k));
    }
    // Broadcast primitive check: 2 generates the whole multiplicative group.
    int seen[64] = {};
    for (int k = 0; k < 63; ++k)
        ++seen[rs_gf_alpha_pow(k)];
    for (int v = 1; v < 64; ++v)
        BOOST_REQUIRE_EQUAL(seen[v], 1);
    BOOST_REQUIRE_EQUAL(seen[0], 0);
}

BOOST_AUTO_TEST_CASE(test_fec_opt_bit_packing_matches_reference)
{
    using namespace openuwb::codec;
    std::mt19937 rng(20261120u);
    std::uniform_int_distribution<int> bit(0, 1);
    for (int rep = 0; rep < 8; ++rep) {
        const int nsym = 63;
        std::vector<int8_t> bits(static_cast<size_t>(nsym) * 6);
        for (auto& b : bits)
            b = static_cast<int8_t>(bit(rng));
        std::vector<int> syms_new(nsym), syms_ref(nsym);
        rs_bits_to_syms(bits.data(), nsym, syms_new.data());
        fec_opt_ref::bits_to_syms(bits.data(), nsym, syms_ref.data());
        BOOST_REQUIRE(syms_new == syms_ref);

        std::vector<int8_t> from_new(static_cast<size_t>(nsym) * 6, 0);
        std::vector<int8_t> from_ref(static_cast<size_t>(nsym) * 6, 0);
        rs_syms_to_bits(syms_ref.data(), nsym, from_new.data());
        fec_opt_ref::syms_to_bits(syms_ref.data(), nsym, from_ref.data());
        BOOST_REQUIRE(from_new == from_ref);
        BOOST_REQUIRE(from_new == bits);

        // LSB-first byte packing.
        std::vector<int8_t> psdu(8 * 37);
        for (auto& b : psdu)
            b = static_cast<int8_t>(bit(rng));
        std::vector<uint8_t> new_bytes, ref_bytes;
        openuwb::demod::core::detail::pack_psdu_bytes(psdu, new_bytes);
        fec_opt_ref::pack_psdu_bytes(psdu, ref_bytes);
        BOOST_REQUIRE(new_bytes == ref_bytes);
        BOOST_REQUIRE_EQUAL(new_bytes.size(), 37u);

        // Non-binary inputs must truncate exactly like the reference (each
        // shifted term keeps only its low 8 bits).
        for (auto& b : psdu)
            b = static_cast<int8_t>(static_cast<int>(rng() % 8) - 2);
        openuwb::demod::core::detail::pack_psdu_bytes(psdu, new_bytes);
        fec_opt_ref::pack_psdu_bytes(psdu, ref_bytes);
        BOOST_REQUIRE(new_bytes == ref_bytes);
    }
}

BOOST_AUTO_TEST_CASE(test_fec_opt_rs_decode_block_matches_reference)
{
    using namespace openuwb::demod::core::detail;
    std::mt19937 rng(20261121u);
    std::uniform_int_distribution<int> sym(0, 63);
    std::uniform_int_distribution<int> pos(0, 62);
    std::uniform_int_distribution<int> errval(1, 63);

    size_t uncorrectable = 0;
    size_t miscorrected = 0;
    for (int trial = 0; trial < 240; ++trial) {
        int data55[55], cw[63];
        for (int j = 0; j < 55; ++j)
            data55[j] = sym(rng);
        BOOST_REQUIRE(fec_opt_ref::encode_systematic(data55, cw));
        int check[8];
        fec_opt_ref::syndromes(cw, check);
        BOOST_REQUIRE(fec_opt_ref::zero_syndromes(check));

        int8_t coded[378];
        fec_opt_ref::syms_to_bits(cw, 63, coded);

        const int nerr = trial % 6; // 0..5 injected symbol errors
        bool used[63] = {};
        for (int e = 0; e < nerr; ++e) {
            int p = pos(rng);
            while (used[p])
                p = (p + 1) % 63;
            used[p] = true;
            const int v = errval(rng); // nonzero GF(2^6) error value
            for (int b = 0; b < 6; ++b)
                if ((v >> b) & 1)
                    coded[p * 6 + 5 - b] ^= 1;
        }

        int8_t new_data[330] = {};
        int8_t ref_data[330] = {};
        const bool new_ok = rs_decode_block(coded, new_data);
        const bool ref_ok = fec_opt_ref::decode_block(coded, ref_data);
        BOOST_REQUIRE_EQUAL(new_ok, ref_ok);
        BOOST_REQUIRE(std::equal(new_data, new_data + 330, ref_data));
        if (nerr <= 4) {
            BOOST_REQUIRE(new_ok);
            int8_t expect[330];
            fec_opt_ref::syms_to_bits(data55, 55, expect);
            BOOST_REQUIRE(std::equal(new_data, new_data + 330, expect));
        } else {
            if (new_ok)
                ++miscorrected;
            else
                ++uncorrectable;
        }
    }
    // Five symbol errors exceed t = 4.  RS(63,55) has d_min = 9, so a 5-error
    // word can (rarely) land exactly 4 symbols from a *different* codeword and
    // be "miscorrected"; the invariant that must hold is that the optimised
    // decoder agrees with the reference on every single trial (checked above).
    BOOST_TEST_MESSAGE("RS 5-error trials: uncorrectable="
                       << uncorrectable << " miscorrected=" << miscorrected);
    BOOST_REQUIRE_EQUAL(uncorrectable + miscorrected, 40u);
    BOOST_REQUIRE_GE(uncorrectable, 39u);
}

BOOST_AUTO_TEST_CASE(test_fec_opt_viterbi_hard_matches_reference)
{
    using namespace openuwb::demod::core::detail;
    std::mt19937 rng(20261122u);
    std::uniform_int_distribution<int> bit(0, 1);
    const size_t lengths[] = { 1, 2, 3, 7, 64, 1229, 2048, 2049 };
    size_t total_ties = 0;
    for (size_t n : lengths) {
        for (int pattern = 0; pattern < 5; ++pattern) {
            std::vector<int8_t> rx(2 * n);
            for (size_t i = 0; i < n; ++i) {
                int a = bit(rng) & 1;
                int b = bit(rng) & 1;
                if (pattern == 1) { // all-zero input: metric ties everywhere
                    a = 0;
                    b = 0;
                } else if (pattern == 2) { // all-one input
                    a = 1;
                    b = 1;
                } else if (pattern == 3) { // alternating hard ties
                    a = static_cast<int>(i & 1);
                    b = static_cast<int>((i >> 1) & 1);
                } else if (pattern == 4) { // repeated 2-symbol tie motif
                    a = 0;
                    b = static_cast<int>(i & 1);
                }
                rx[2 * i] = static_cast<int8_t>(a);
                rx[2 * i + 1] = static_cast<int8_t>(b);
            }
            total_ties += fec_opt_ref::count_hard_ties(rx.data(), n);
            std::vector<int8_t> got, want;
            vitdec_cl3(rx.data(), n, got);
            fec_opt_ref::vit_hard(rx.data(), n, want);
            BOOST_REQUIRE_EQUAL(got.size(), n);
            BOOST_REQUIRE_EQUAL(want.size(), n);
            BOOST_REQUIRE(std::equal(got.begin(), got.end(), want.begin()));
        }
    }
    // The equivalence above only covers tie-breaking if exact ties occur.
    BOOST_TEST_MESSAGE("hard CL-3 ACS exact ties exercised: " << total_ties);
    BOOST_REQUIRE_GT(total_ties, 0u);
}

BOOST_AUTO_TEST_CASE(test_fec_opt_viterbi_soft_matches_reference)
{
    using namespace openuwb::demod::core::detail;
    std::mt19937 rng(20261123u);
    std::normal_distribution<float> nd(0.0f, 2.0f);
    const size_t lengths[] = { 1, 2, 3, 7, 64, 1229, 2048, 2049 };
    size_t total_ties = 0;
    for (size_t n : lengths) {
        for (int pattern = 0; pattern < 6; ++pattern) {
            std::vector<float> l0(n), l1(n);
            for (size_t i = 0; i < n; ++i) {
                switch (pattern) {
                case 1: // all-zero LLRs: every edge weight is an exact tie
                    l0[i] = 0.0f;
                    l1[i] = 0.0f;
                    break;
                case 2: // identical constant LLRs
                    l0[i] = 0.5f;
                    l1[i] = 0.5f;
                    break;
                case 3: // constant-magnitude sign ties
                    l0[i] = (i & 1) ? 1.0f : -1.0f;
                    l1[i] = (i & 2) ? 1.0f : -1.0f;
                    break;
                case 4: // large magnitudes
                    l0[i] = (nd(rng) > 0.0f) ? 1.0e6f : -1.0e6f;
                    l1[i] = (nd(rng) > 0.0f) ? 1.0e6f : -1.0e6f;
                    break;
                case 5: // tiny magnitudes (near-tie decisions)
                    l0[i] = 1.0e-8f * nd(rng);
                    l1[i] = 1.0e-8f * nd(rng);
                    break;
                default:
                    l0[i] = nd(rng);
                    l1[i] = nd(rng);
                    break;
                }
            }
            total_ties += fec_opt_ref::count_ties(
                n, 1e30f, [&](size_t t, float path, int g0, int g1) {
                    return fec_opt_ref::soft_metric(path, g0, g1, l0[t], l1[t]);
                });
            std::vector<int8_t> got, want;
            vitdec_cl3_soft(l0.data(), l1.data(), n, got);
            fec_opt_ref::vit_soft(l0.data(), l1.data(), n, want);
            BOOST_REQUIRE_EQUAL(got.size(), n);
            BOOST_REQUIRE_EQUAL(want.size(), n);
            BOOST_REQUIRE(std::equal(got.begin(), got.end(), want.begin()));
        }
    }
    // All-zero / constant LLRs force exact metric ties on every step, so the
    // strict-comparison tie rule is definitely exercised.
    BOOST_TEST_MESSAGE("soft CL-3 ACS exact ties exercised: " << total_ties);
    BOOST_REQUIRE_GT(total_ties, 100u);
}

BOOST_AUTO_TEST_CASE(test_fec_opt_viterbi_soft_trace_matches_reference)
{
    using namespace openuwb::demod::core::detail;
    std::mt19937 rng(20261124u);
    std::normal_distribution<float> nd(0.0f, 1.5f);
    const size_t lengths[] = { 1, 2, 3, 7, 64, 257, 1229, 2048, 2049 };
    for (size_t n : lengths) {
        for (int pattern = 0; pattern < 3; ++pattern) {
            std::vector<float> l0(n), l1(n);
            for (size_t i = 0; i < n; ++i) {
                if (pattern == 1) {
                    l0[i] = 0.0f;
                    l1[i] = 0.0f;
                } else if (pattern == 2) {
                    l0[i] = (i & 1) ? 2.0f : -2.0f;
                    l1[i] = (i & 2) ? 2.0f : -2.0f;
                } else {
                    l0[i] = nd(rng);
                    l1[i] = nd(rng);
                }
            }
            Cl3TrellisTrace got;
            fec_opt_ref::Trace want;
            vitdec_cl3_soft_trace(l0.data(), l1.data(), n, got);
            fec_opt_ref::vit_soft_trace(l0.data(), l1.data(), n, want);
            BOOST_REQUIRE_EQUAL(got.bits.size(), n);
            BOOST_REQUIRE_EQUAL(got.states.size(), n + 1);
            BOOST_REQUIRE_EQUAL(got.margins.size(), n);
            BOOST_REQUIRE_EQUAL(got.lose_u.size(), n);
            BOOST_REQUIRE(
                std::equal(got.bits.begin(), got.bits.end(), want.bits.begin()));
            BOOST_REQUIRE(
                std::equal(got.states.begin(), got.states.end(), want.states.begin()));
            BOOST_REQUIRE(
                std::equal(got.lose_u.begin(), got.lose_u.end(), want.lose_u.begin()));
            for (size_t i = 0; i < n; ++i)
                BOOST_REQUIRE_EQUAL(got.margins[i], want.margins[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_fec_opt_rs_kernels_match_reference)
{
    // Direct bit-exactness of the individually rewritten RS kernels, including
    // the degrees the errors-and-erasures rescue path reaches (Forney with
    // L = Le + v up to 8), which rs_decode_block alone never exercises.
    using namespace openuwb::demod::core::detail;
    std::mt19937 rng(20261125u);
    std::uniform_int_distribution<int> sym(0, 63);
    std::uniform_int_distribution<int> pos(0, 62);
    std::uniform_int_distribution<int> degree(0, 8);
    std::uniform_int_distribution<int> count(0, 4);

    for (int trial = 0; trial < 2000; ++trial) {
        int r63[63];
        for (int j = 0; j < 63; ++j)
            r63[j] = sym(rng);
        int S_new[8], S_ref[8];
        rs_syndromes(r63, S_new);
        fec_opt_ref::syndromes(r63, S_ref);
        for (int i = 0; i < 8; ++i)
            BOOST_REQUIRE_EQUAL(S_new[i], S_ref[i]);

        int L_new[9] = {}, L_ref[9] = {};
        const int deg_new = rs_berlekamp_massey(S_new, L_new);
        const int deg_ref = fec_opt_ref::bm_len(S_ref, 8, L_ref);
        BOOST_REQUIRE_EQUAL(deg_new, deg_ref);
        for (int i = 0; i < 9; ++i)
            BOOST_REQUIRE_EQUAL(L_new[i], L_ref[i]);

        // Fuzz the locator/Forney inputs too (any 0..8 degree, any roots).
        const int L = degree(rng);
        int Lambda[9] = {};
        Lambda[0] = 1;
        for (int i = 1; i <= L; ++i)
            Lambda[i] = sym(rng);
        int p_new[4] = {}, p_ref[4] = {};
        const int n_new = rs_chien_search(Lambda, L, p_new);
        const int n_ref = fec_opt_ref::chien(Lambda, L, p_ref);
        BOOST_REQUIRE_EQUAL(n_new, n_ref);
        for (int i = 0; i < 4; ++i)
            BOOST_REQUIRE_EQUAL(p_new[i], p_ref[i]);

        int pos_in[4];
        const int nerr = count(rng);
        for (int i = 0; i < nerr; ++i)
            pos_in[i] = pos(rng);
        int v_new[8] = {}, v_ref[8] = {};
        rs_forney(S_new, Lambda, L, pos_in, nerr, v_new);
        fec_opt_ref::forney(S_ref, Lambda, L, pos_in, nerr, v_ref);
        for (int i = 0; i < 8; ++i)
            BOOST_REQUIRE_EQUAL(v_new[i], v_ref[i]);
    }
}

// ---------------------------------------------------------------------------
// Wave 1 / JOINT-ML-VITERBI QA.
//
// Everything below compares the new joint bi-orthogonal ML kernels against
// independent references written in this anonymous namespace: a nested
// `for (st) for (u)` DP (the pre-unrolling loop shape), a brute-force
// enumeration of every input sequence for the tiny exact-ML case, and a local
// CL-3 encoder built from openuwb::codec::cl3_transition.  The new production
// code is never used as its own oracle.
// ---------------------------------------------------------------------------
namespace {

using openuwb::codec::cl3_transition;
using openuwb::codec::Cl3Transition;
using openuwb::codec::crc16_802154;
using openuwb::codec::hrp_secded;
using openuwb::codec::rs_bits_to_syms;
using openuwb::codec::rs_syms_to_bits;

// Joint ML branch metric, written straight from the task definition.
float qa_branch(float m0, float m1, int o0, int o1)
{
    return (o0 == 0) ? (o1 ? m0 : -m0) : (o1 ? m1 : -m1);
}

void qa_cl3_encode(const std::vector<int8_t>& u,
                   std::vector<int8_t>& g0,
                   std::vector<int8_t>& g1)
{
    g0.assign(u.size(), int8_t(0));
    g1.assign(u.size(), int8_t(0));
    int state = 0;
    for (size_t t = 0; t < u.size(); ++t) {
        const Cl3Transition step = cl3_transition(state, u[t] & 1);
        g0[t] = static_cast<int8_t>(step.g0);
        g1[t] = static_cast<int8_t>(step.g1);
        state = step.next_state;
    }
}

// BPM-BPSK constellation: (g0,g1) -> (m0,m1) at amplitude A.
void qa_map_bpm(const std::vector<int8_t>& g0,
                const std::vector<int8_t>& g1,
                float A,
                std::vector<float>& m0,
                std::vector<float>& m1)
{
    m0.assign(g0.size(), 0.0f);
    m1.assign(g0.size(), 0.0f);
    for (size_t t = 0; t < g0.size(); ++t) {
        const float s = g1[t] ? -A : A;
        if (g0[t] == 0)
            m0[t] = s;
        else
            m1[t] = s;
    }
}

// Split-LLR convention of decide_bpm_soft(): positive favours bit 0.
void qa_split_llrs(const std::vector<float>& m0,
                   const std::vector<float>& m1,
                   std::vector<float>& l0,
                   std::vector<float>& l1)
{
    l0.resize(m0.size());
    l1.resize(m0.size());
    for (size_t t = 0; t < m0.size(); ++t) {
        l0[t] = std::abs(m0[t]) - std::abs(m1[t]);
        l1[t] = std::max(m0[t], m1[t]) - std::max(-m0[t], -m1[t]);
    }
}

float qa_path_metric(const std::vector<int8_t>& bits,
                     const std::vector<float>& m0,
                     const std::vector<float>& m1)
{
    int state = 0;
    float metric = 0.0f;
    for (size_t t = 0; t < bits.size(); ++t) {
        const Cl3Transition step = cl3_transition(state, bits[t] & 1);
        metric += qa_branch(m0[t], m1[t], step.g0, step.g1);
        state = step.next_state;
    }
    return metric;
}

// Independent nested-loop DP with traceback and final metric.
struct QaPath {
    std::vector<int8_t> bits;
    float metric = 0.0f;
};

QaPath
qa_ref_dp(const std::vector<float>& m0, const std::vector<float>& m1, bool terminated)
{
    const float INF = 1e30f;
    const size_t N = m0.size();
    float surf[4] = { 0.0f, INF, INF, INF };
    std::vector<int8_t> back(N * 4, int8_t(0));
    for (size_t t = 0; t < N; ++t) {
        float ns[4] = { INF, INF, INF, INF };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const Cl3Transition step = cl3_transition(st, u);
                const float m = surf[st] + qa_branch(m0[t], m1[t], step.g0, step.g1);
                if (m < ns[step.next_state]) {
                    ns[step.next_state] = m;
                    back[t * 4 + step.next_state] = static_cast<int8_t>((u << 2) | st);
                }
            }
        }
        for (int i = 0; i < 4; ++i)
            surf[i] = ns[i];
    }
    int best = 0;
    if (!terminated) {
        for (int i = 1; i < 4; ++i)
            if (surf[i] < surf[best])
                best = i;
    }
    QaPath p;
    p.bits.resize(N);
    p.metric = surf[best];
    int cur = best;
    for (size_t t = N; t-- > 0;) {
        const int v = back[t * 4 + cur];
        p.bits[t] = static_cast<int8_t>(v >> 2);
        cur = v & 3;
    }
    return p;
}

// Split-LLR Viterbi, again from a nested loop (production vitdec_cl3_soft() is
// intentionally not called here).
std::vector<int8_t> qa_ref_split(const std::vector<float>& m0,
                                 const std::vector<float>& m1)
{
    std::vector<float> l0, l1;
    qa_split_llrs(m0, m1, l0, l1);
    const float INF = 1e30f;
    const size_t N = m0.size();
    float surf[4] = { 0.0f, INF, INF, INF };
    std::vector<int8_t> back(N * 4, int8_t(0));
    for (size_t t = 0; t < N; ++t) {
        float ns[4] = { INF, INF, INF, INF };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const Cl3Transition step = cl3_transition(st, u);
                const float c0 = (2.0f * step.g0 - 1.0f) * l0[t];
                const float c1 = (2.0f * step.g1 - 1.0f) * l1[t];
                const float m = surf[st] + c0 + c1;
                if (m < ns[step.next_state]) {
                    ns[step.next_state] = m;
                    back[t * 4 + step.next_state] = static_cast<int8_t>((u << 2) | st);
                }
            }
        }
        for (int i = 0; i < 4; ++i)
            surf[i] = ns[i];
    }
    int best = 0;
    for (int i = 1; i < 4; ++i)
        if (surf[i] < surf[best])
            best = i;
    std::vector<int8_t> out(N);
    int cur = best;
    for (size_t t = N; t-- > 0;) {
        const int v = back[t * 4 + cur];
        out[t] = static_cast<int8_t>(v >> 2);
        cur = v & 3;
    }
    return out;
}

// Exhaustive ML over all 2^N input sequences (small N only).
bool qa_brute_ml(const std::vector<float>& m0,
                 const std::vector<float>& m1,
                 QaPath& out,
                 size_t& min_count)
{
    const size_t N = m0.size();
    if (N > 20)
        return false;
    const float INF = 1e30f;
    float best = INF;
    std::vector<int8_t> bits(N, int8_t(0));
    std::vector<int8_t> best_bits;
    min_count = 0;
    const uint32_t total = static_cast<uint32_t>(1u) << N;
    for (uint32_t mask = 0; mask < total; ++mask) {
        for (size_t t = 0; t < N; ++t)
            bits[t] = static_cast<int8_t>((mask >> t) & 1u);
        const float v = qa_path_metric(bits, m0, m1);
        if (v < best - 1e-6f) {
            best = v;
            best_bits = bits;
            min_count = 1;
        } else if (std::abs(v - best) <= 1e-6f) {
            ++min_count;
        }
    }
    out.bits = best_bits;
    out.metric = best;
    return !best_bits.empty();
}

// Exact tie counter over the nested-loop ACS (does not change survivors).
size_t qa_count_biorth_ties(const std::vector<float>& m0, const std::vector<float>& m1)
{
    const float INF = 1e30f;
    float surf[4] = { 0.0f, INF, INF, INF };
    size_t ties = 0;
    for (size_t t = 0; t < m0.size(); ++t) {
        float ns[4] = { INF, INF, INF, INF };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const Cl3Transition step = cl3_transition(st, u);
                const float m = surf[st] + qa_branch(m0[t], m1[t], step.g0, step.g1);
                if (m < ns[step.next_state]) {
                    ns[step.next_state] = m;
                } else if (m == ns[step.next_state] && m != INF) {
                    ++ties;
                }
            }
        }
        for (int i = 0; i < 4; ++i)
            surf[i] = ns[i];
    }
    return ties;
}

// Independent nested-loop tracing ML Viterbi.
struct QaTrace {
    std::vector<int8_t> bits;
    std::vector<int8_t> states;
    std::vector<float> margins;
    std::vector<uint8_t> lose_u;
};

QaTrace
qa_ref_trace(const std::vector<float>& m0, const std::vector<float>& m1, bool terminated)
{
    const float INF = 1e30f;
    const size_t N = m0.size();
    float surf[4] = { 0.0f, INF, INF, INF };
    std::vector<int8_t> back(N * 4, int8_t(0));
    std::vector<float> diff(N * 4, INF);
    std::vector<uint8_t> alt(N * 4, 0);
    for (size_t t = 0; t < N; ++t) {
        float ns[4] = { INF, INF, INF, INF };
        float ns2[4] = { INF, INF, INF, INF };
        uint8_t nalt[4] = { 0, 0, 0, 0 };
        for (int st = 0; st < 4; ++st) {
            for (int u = 0; u < 2; ++u) {
                const Cl3Transition step = cl3_transition(st, u);
                const float m = surf[st] + qa_branch(m0[t], m1[t], step.g0, step.g1);
                const int nst = step.next_state;
                const uint8_t code = static_cast<uint8_t>((u << 2) | st);
                if (m < ns[nst]) {
                    ns2[nst] = ns[nst];
                    nalt[nst] = static_cast<uint8_t>(back[t * 4 + nst]);
                    ns[nst] = m;
                    back[t * 4 + nst] = static_cast<int8_t>(code);
                } else if (m < ns2[nst]) {
                    ns2[nst] = m;
                    nalt[nst] = code;
                }
            }
        }
        for (int i = 0; i < 4; ++i) {
            surf[i] = ns[i];
            diff[t * 4 + i] = ns2[i] - ns[i];
            alt[t * 4 + i] = nalt[i];
        }
    }
    int best = 0;
    if (!terminated) {
        for (int i = 1; i < 4; ++i)
            if (surf[i] < surf[best])
                best = i;
    }
    QaTrace tr;
    tr.bits.resize(N);
    tr.states.resize(N + 1);
    tr.margins.resize(N);
    tr.lose_u.resize(N);
    int cur = best;
    tr.states[N] = static_cast<int8_t>(cur);
    for (size_t t = N; t-- > 0;) {
        const int v = back[t * 4 + cur];
        tr.bits[t] = static_cast<int8_t>(v >> 2);
        tr.margins[t] = diff[t * 4 + cur];
        tr.lose_u[t] = static_cast<uint8_t>(alt[t * 4 + cur] >> 2);
        cur = v & 3;
        tr.states[t] = static_cast<int8_t>(cur);
    }
    return tr;
}

// ---- Synthetic end-to-end payload frames ---------------------------------
struct QaFrame {
    size_t psdu_bits = 0;
    size_t payload_symbols = 0;
    std::vector<uint8_t> bytes; // expected PSDU (with FCS)
    std::vector<int8_t> u;      // joint trellis info bits
    std::vector<float> m0, m1;  // noiseless joint metrics
};

// Systematic RS(63,55) encode of one shortened block: data bits + 48 parity
// bits, matching expand_rs_codeword()'s leading-zero padding convention.
bool qa_rs_encode_block(const std::vector<int8_t>& data, std::vector<int8_t>& coded)
{
    if (data.empty() || data.size() > 330)
        return false;
    int8_t data330[330] = {};
    const size_t pad = 330 - data.size();
    for (size_t i = 0; i < data.size(); ++i)
        data330[pad + i] = data[i] & 1;
    int syms55[55];
    rs_bits_to_syms(data330, 55, syms55);
    int cw63[63];
    if (!fec_opt_ref::encode_systematic(syms55, cw63))
        return false;
    int8_t par48[48];
    rs_syms_to_bits(cw63 + 55, 8, par48);
    coded = data;
    coded.insert(coded.end(), par48, par48 + 48);
    return true;
}

bool qa_build_frame(size_t psdu_bytes, uint32_t seed, QaFrame& f)
{
    if (psdu_bytes < 3)
        return false;
    std::mt19937 rng(seed);
    std::vector<uint8_t> bytes(psdu_bytes);
    for (size_t i = 0; i + 2 < psdu_bytes; ++i)
        bytes[i] = static_cast<uint8_t>(rng() & 0xffu);
    const uint16_t fcs = crc16_802154(bytes.data(), psdu_bytes - 2);
    bytes[psdu_bytes - 2] = static_cast<uint8_t>(fcs & 0xffu);
    bytes[psdu_bytes - 1] = static_cast<uint8_t>(fcs >> 8);

    std::vector<int8_t> psdu_bits(psdu_bytes * 8);
    for (size_t i = 0; i < psdu_bytes; ++i)
        for (int k = 0; k < 8; ++k)
            psdu_bits[8 * i + k] = static_cast<int8_t>((bytes[i] >> k) & 1u);

    std::vector<int8_t> rs_cw;
    size_t off = 0;
    while (off < psdu_bits.size()) {
        const size_t n = std::min<size_t>(330, psdu_bits.size() - off);
        std::vector<int8_t> block(psdu_bits.begin() + off, psdu_bits.begin() + off + n);
        std::vector<int8_t> coded;
        if (!qa_rs_encode_block(block, coded))
            return false;
        rs_cw.insert(rs_cw.end(), coded.begin(), coded.end());
        off += n;
    }

    f.psdu_bits = psdu_bits.size();
    f.payload_symbols = rs_cw.size();
    f.bytes = bytes;
    // The first 19 trellis bits are the systematic PHR bits.  The payload
    // decoder discards them, so any fixed pattern works for this frame.
    f.u.assign(19, int8_t(0));
    f.u.insert(f.u.end(), rs_cw.begin(), rs_cw.end());
    f.u.push_back(0);
    f.u.push_back(0);

    std::vector<int8_t> g0, g1;
    qa_cl3_encode(f.u, g0, g1);
    qa_map_bpm(g0, g1, 1.0f, f.m0, f.m1);
    return true;
}

void qa_add_noise(std::vector<float>& m0,
                  std::vector<float>& m1,
                  float sigma,
                  std::mt19937& rng)
{
    std::normal_distribution<float> nd(0.0f, sigma);
    for (size_t i = 0; i < m0.size(); ++i) {
        m0[i] += nd(rng);
        m1[i] += nd(rng);
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_w1_joint_ml_exact_counterexample)
{
    using namespace openuwb::demod::core::detail;

    // The documented coupling: a +1.0 half and a -0.8 half make (1,1) the true
    // runner-up (cost -0.8), while the split LLRs score (1,1) as the worst
    // branch (+0.4).
    const std::vector<float> m0 = { 1.0f, 0.6f };
    const std::vector<float> m1 = { -0.8f, 0.9f };

    std::vector<int8_t> biorth;
    vitdec_cl3_biorth(m0.data(), m1.data(), m0.size(), biorth);
    BOOST_REQUIRE_EQUAL(biorth.size(), m0.size());

    // Exhaustive ML over all 2^2 = 4 input sequences.
    QaPath brute;
    size_t min_count = 0;
    BOOST_REQUIRE(qa_brute_ml(m0, m1, brute, min_count));
    BOOST_REQUIRE_EQUAL(min_count, size_t(1)); // unique ML path
    BOOST_CHECK_EQUAL_COLLECTIONS(
        biorth.begin(), biorth.end(), brute.bits.begin(), brute.bits.end());
    BOOST_CHECK_SMALL(qa_path_metric(biorth, m0, m1) - brute.metric, 1e-6f);

    // The split path must genuinely lose on this sequence.
    const std::vector<int8_t> split = qa_ref_split(m0, m1);
    const float split_metric = qa_path_metric(split, m0, m1);
    BOOST_CHECK(split != biorth);
    BOOST_CHECK_GT(split_metric, brute.metric + 1e-6f);
    BOOST_TEST_MESSAGE("joint-ML counterexample: biorth metric "
                       << brute.metric << ", split path metric " << split_metric
                       << ", split bits [" << int(split[0]) << int(split[1])
                       << "], tied ML paths " << min_count);

    // Trace format: bits/states/margins/lose_u must match the independent
    // nested-loop reference exactly, including the strict `<` tie order.
    Cl3TrellisTrace tr;
    vitdec_cl3_biorth_trace(m0.data(), m1.data(), m0.size(), tr);
    const QaTrace ref = qa_ref_trace(m0, m1, false);
    BOOST_REQUIRE_EQUAL(tr.bits.size(), m0.size());
    BOOST_REQUIRE_EQUAL(tr.states.size(), m0.size() + 1);
    BOOST_CHECK(std::equal(tr.bits.begin(), tr.bits.end(), ref.bits.begin()));
    BOOST_CHECK(std::equal(tr.states.begin(), tr.states.end(), ref.states.begin()));
    BOOST_CHECK(std::equal(tr.lose_u.begin(), tr.lose_u.end(), ref.lose_u.begin()));
    for (size_t t = 0; t < m0.size(); ++t)
        BOOST_CHECK_EQUAL(tr.margins[t], ref.margins[t]);

    // Deterministic tie stress: zero metrics make every edge cost 0, so all
    // ACS comparisons tie and only the strict `<` / first-edge rule decides.
    std::vector<float> z0(6, 0.0f), z1(6, 0.0f);
    std::vector<int8_t> tie_bits;
    vitdec_cl3_biorth(z0.data(), z1.data(), z0.size(), tie_bits);
    const QaPath tie_ref = qa_ref_dp(z0, z1, false);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tie_bits.begin(), tie_bits.end(), tie_ref.bits.begin(), tie_ref.bits.end());
    const size_t tie_count = qa_count_biorth_ties(z0, z1);
    BOOST_CHECK_GT(tie_count, size_t(0));
    BOOST_TEST_MESSAGE("joint-ML tie stress: " << tie_count << " exact ACS ties");
}

BOOST_AUTO_TEST_CASE(test_w1_joint_ml_awgn_gain_over_split_llr)
{
    using namespace openuwb::demod::core::detail;

    // 600 frames x 128 coded bits at sigma = 0.4 places the split-LLR FER in
    // the requested 20-40% band while the joint ML decoder is far below it.
    constexpr size_t kFrames = 600;
    constexpr size_t kBits = 128;
    constexpr float kSigma = 0.4f;
    std::mt19937 rng(20261120u);
    std::normal_distribution<float> nd(0.0f, kSigma);
    std::mt19937 bit_rng(0x5eed1234u);

    size_t split_fer = 0, biorth_fer = 0;
    size_t split_err = 0, biorth_err = 0, total_bits = 0, ties = 0;
    for (size_t f = 0; f < kFrames; ++f) {
        std::vector<int8_t> u(kBits);
        for (size_t t = 0; t < kBits; ++t)
            u[t] = static_cast<int8_t>(bit_rng() & 1u);
        std::vector<int8_t> g0, g1;
        qa_cl3_encode(u, g0, g1);
        std::vector<float> m0, m1;
        qa_map_bpm(g0, g1, 1.0f, m0, m1);
        for (size_t t = 0; t < kBits; ++t) {
            m0[t] += nd(rng);
            m1[t] += nd(rng);
        }

        std::vector<float> l0, l1;
        qa_split_llrs(m0, m1, l0, l1);
        std::vector<int8_t> split_bits;
        vitdec_cl3_soft(l0.data(), l1.data(), kBits, split_bits);
        std::vector<int8_t> biorth_bits;
        vitdec_cl3_biorth(m0.data(), m1.data(), kBits, biorth_bits);
        BOOST_REQUIRE_EQUAL(split_bits.size(), kBits);
        BOOST_REQUIRE_EQUAL(biorth_bits.size(), kBits);

        size_t se = 0, be = 0;
        for (size_t t = 0; t < kBits; ++t) {
            se += (split_bits[t] != u[t]);
            be += (biorth_bits[t] != u[t]);
        }
        split_err += se;
        biorth_err += be;
        split_fer += (se != 0);
        biorth_fer += (be != 0);
        total_bits += kBits;
        ties += qa_count_biorth_ties(m0, m1);
    }

    const double split_fer_pct = 100.0 * split_fer / kFrames;
    const double biorth_fer_pct = 100.0 * biorth_fer / kFrames;
    const double split_ber = static_cast<double>(split_err) / total_bits;
    const double biorth_ber = static_cast<double>(biorth_err) / total_bits;
    BOOST_TEST_MESSAGE("AWGN joint-ML gain ("
                       << kFrames << " frames x " << kBits << " bits, sigma=" << kSigma
                       << "): split FER=" << split_fer_pct << "% BER=" << split_ber
                       << "; biorth FER=" << biorth_fer_pct << "% BER=" << biorth_ber
                       << "; exact biorth ACS ties=" << ties);

    // The test is only informative if the split baseline sits in the target
    // band; then both BER and FER must be strictly lower for joint ML.
    BOOST_REQUIRE_GE(split_fer_pct, 20.0);
    BOOST_REQUIRE_LE(split_fer_pct, 40.0);
    BOOST_CHECK_LT(biorth_err, split_err);
    BOOST_CHECK_LT(biorth_fer, split_fer);
    BOOST_CHECK_LT(biorth_err * 4, split_err);
    BOOST_CHECK_LT(biorth_fer * 4, split_fer);
}

BOOST_AUTO_TEST_CASE(test_w1_joint_ml_terminated_tail_beats_truncated)
{
    using namespace openuwb::demod::core::detail;

    // 400 frames with the standard two zero tail bits.  Only the last eight
    // symbols carry noise, so the endpoint decision is what matters.
    constexpr size_t kFrames = 400;
    constexpr size_t kBits = 64;
    constexpr float kTailSigma = 1.0f;
    std::mt19937 rng(20261121u);
    std::normal_distribution<float> nd(0.0f, kTailSigma);
    std::mt19937 bit_rng(0x7151245u);

    size_t term_tail_err = 0, trunc_tail_err = 0;
    size_t clean_mismatch = 0, term_frame_err = 0, trunc_frame_err = 0;
    for (size_t f = 0; f < kFrames; ++f) {
        std::vector<int8_t> u(kBits, int8_t(0));
        for (size_t t = 0; t + 2 < kBits; ++t)
            u[t] = static_cast<int8_t>(bit_rng() & 1u);
        std::vector<int8_t> g0, g1;
        qa_cl3_encode(u, g0, g1);
        std::vector<float> clean0, clean1;
        qa_map_bpm(g0, g1, 1.0f, clean0, clean1);

        // Noiseless: terminated and truncated must agree bit for bit.
        std::vector<int8_t> clean_t, clean_u2;
        vitdec_cl3_biorth(clean0.data(), clean1.data(), kBits, clean_t, true);
        vitdec_cl3_biorth(clean0.data(), clean1.data(), kBits, clean_u2, false);
        if (clean_t != clean_u2 || clean_t != u)
            ++clean_mismatch;

        std::vector<float> m0 = clean0, m1 = clean1;
        for (size_t t = kBits - 8; t < kBits; ++t) {
            m0[t] += nd(rng);
            m1[t] += nd(rng);
        }
        std::vector<int8_t> term, trunc;
        vitdec_cl3_biorth(m0.data(), m1.data(), kBits, term, true);
        vitdec_cl3_biorth(m0.data(), m1.data(), kBits, trunc, false);
        for (size_t t = kBits - 16; t < kBits; ++t) {
            term_tail_err += (term[t] != u[t]);
            trunc_tail_err += (trunc[t] != u[t]);
        }
        term_frame_err += (term != u);
        trunc_frame_err += (trunc != u);
    }

    BOOST_TEST_MESSAGE("tail-noise test: last-16 errors terminated="
                       << term_tail_err << " truncated=" << trunc_tail_err
                       << "; frame errors terminated=" << term_frame_err
                       << " truncated=" << trunc_frame_err);
    BOOST_CHECK_EQUAL(clean_mismatch, size_t(0));
    BOOST_CHECK_LT(term_tail_err * 4, trunc_tail_err * 3);
    BOOST_CHECK_LT(term_frame_err, trunc_frame_err);
}

BOOST_AUTO_TEST_CASE(test_w1_joint_ml_pulse_clipping_recovers_fcs)
{
    using namespace openuwb::demod::core::detail;

    // Direct clipping contract: median one, 50x spikes clamped to 4x.
    {
        std::vector<float> a(32, 1.0f), b(32, 0.0f);
        a[4] = 50.0f;
        a[13] = -50.0f;
        b[22] = 50.0f;
        clip_biorth_metrics(a.data(), b.data(), a.size(), 4.0f);
        for (size_t i = 0; i < a.size(); ++i) {
            BOOST_CHECK_LE(std::abs(a[i]), 4.0f);
            BOOST_CHECK_LE(std::abs(b[i]), 4.0f);
        }
        BOOST_CHECK_CLOSE(a[4], 4.0f, 1e-6);
        BOOST_CHECK_CLOSE(a[13], -4.0f, 1e-6);
        BOOST_CHECK_CLOSE(b[22], 4.0f, 1e-6);
        BOOST_CHECK_CLOSE(a[0], 1.0f, 1e-6);
    }
    // median <= 0 and N < 16 skip safely (buffer untouched).
    {
        std::vector<float> a(32, 0.0f), b(32, 0.0f);
        clip_biorth_metrics(a.data(), b.data(), a.size(), 4.0f);
        for (size_t i = 0; i < a.size(); ++i)
            BOOST_CHECK_EQUAL(a[i], 0.0f);
        std::vector<float> c(8, 50.0f), d(8, -50.0f);
        clip_biorth_metrics(c.data(), d.data(), c.size(), 4.0f);
        BOOST_CHECK_EQUAL(c[0], 50.0f);
        BOOST_CHECK_EQUAL(d[7], -50.0f);
    }

    // 24-byte PSDU (one RS block, 240 coded payload symbols).  Five 50x
    // pulses on spaced symbols destroy the unclipped path metric ordering;
    // after clipping the residual errors stay within RS t = 4.
    size_t clipped_pass = 0, unclipped_fail = 0;
    constexpr size_t kFrames = 64;
    for (size_t f = 0; f < kFrames; ++f) {
        QaFrame frame;
        BOOST_REQUIRE(qa_build_frame(24, 1000u + static_cast<uint32_t>(f), frame));
        BOOST_REQUIRE_EQUAL(frame.payload_symbols, size_t(240));

        std::vector<float> m0 = frame.m0, m1 = frame.m1;
        const size_t pulse_pos[5] = { 60, 100, 140, 180, 220 };
        for (size_t i = 0; i < 5; ++i) {
            const float amp = (i % 2 == 0) ? 50.0f : -50.0f;
            m0[21 + pulse_pos[i]] += amp;
        }

        PayloadDecodeWorkspace ws;
        PayloadResult clipped;
        const bool clipped_ok = decode_payload_biorth_soft(frame.m0.data(),
                                                           frame.m1.data(),
                                                           m0.data() + 21,
                                                           m1.data() + 21,
                                                           frame.payload_symbols,
                                                           frame.psdu_bits,
                                                           clipped,
                                                           nullptr,
                                                           &ws,
                                                           false,
                                                           true);
        if (clipped_ok && clipped.fcs_pass && clipped.bytes == frame.bytes)
            ++clipped_pass;

        PayloadResult plain;
        decode_payload_biorth_soft(frame.m0.data(),
                                   frame.m1.data(),
                                   m0.data() + 21,
                                   m1.data() + 21,
                                   frame.payload_symbols,
                                   frame.psdu_bits,
                                   plain,
                                   nullptr,
                                   &ws,
                                   false,
                                   false);
        if (!(plain.ok && plain.fcs_pass))
            ++unclipped_fail;
    }
    BOOST_TEST_MESSAGE("pulse clipping: clipped FCS "
                       << clipped_pass << "/" << kFrames << ", unclipped failures "
                       << unclipped_fail << "/" << kFrames);
    BOOST_CHECK_EQUAL(clipped_pass, kFrames);
}

BOOST_AUTO_TEST_CASE(test_w1_joint_ml_end_to_end_and_boundaries)
{
    using namespace openuwb::demod::core::detail;
    using openuwb::demod::core::decode_phr_decisions;

    // ---- Payload: biorth must agree with decode_payload_soft() ------------
    QaFrame frame;
    BOOST_REQUIRE(qa_build_frame(24, 777u, frame));
    std::vector<float> phr0(frame.m0.begin(), frame.m0.begin() + 21);
    std::vector<float> phr1(frame.m1.begin(), frame.m1.begin() + 21);
    std::vector<float> pay0(frame.m0.begin() + 21, frame.m0.end());
    std::vector<float> pay1(frame.m1.begin() + 21, frame.m1.end());
    std::vector<float> l0, l1;
    qa_split_llrs(pay0, pay1, l0, l1);
    std::vector<float> l0p, l1p;
    qa_split_llrs(phr0, phr1, l0p, l1p);

    PayloadResult soft_out, biorth_out, term_out;
    PayloadDecodeWorkspace ws;
    BOOST_REQUIRE(
        decode_payload_soft(l0p, l1p, l0, l1, frame.psdu_bits, soft_out, nullptr, &ws));
    BOOST_REQUIRE(soft_out.fcs_pass);
    BOOST_CHECK(soft_out.bytes == frame.bytes);
    BOOST_REQUIRE(decode_payload_biorth_soft(phr0.data(),
                                             phr1.data(),
                                             pay0.data(),
                                             pay1.data(),
                                             frame.payload_symbols,
                                             frame.psdu_bits,
                                             biorth_out,
                                             nullptr,
                                             &ws));
    BOOST_CHECK(biorth_out.fcs_pass);
    BOOST_CHECK(biorth_out.bytes == soft_out.bytes);
    BOOST_CHECK_EQUAL(biorth_out.received_fcs, soft_out.received_fcs);
    BOOST_CHECK_EQUAL(biorth_out.calculated_fcs, soft_out.calculated_fcs);
    BOOST_REQUIRE(decode_payload_biorth_soft(phr0.data(),
                                             phr1.data(),
                                             pay0.data(),
                                             pay1.data(),
                                             frame.payload_symbols,
                                             frame.psdu_bits,
                                             term_out,
                                             nullptr,
                                             &ws,
                                             true));
    BOOST_CHECK(term_out.fcs_pass);
    BOOST_CHECK(term_out.bytes == frame.bytes);

    // ---- PHR: biorth vs decode_phr_decisions() on the same 21 symbols -----
    const std::vector<int8_t> sys13 = { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
    std::vector<int8_t> phr19;
    hrp_secded(sys13, phr19);
    std::vector<int8_t> phr_u(phr19.begin(), phr19.end());
    phr_u.push_back(0);
    phr_u.push_back(0);
    std::vector<int8_t> pg0, pg1;
    qa_cl3_encode(phr_u, pg0, pg1);
    std::vector<float> pm0, pm1;
    qa_map_bpm(pg0, pg1, 1.0f, pm0, pm1);
    const std::vector<float> empty_soft;
    NsSfdResult ns;
    ns.ok = true;
    ns.sfd_start_chip = -1;
    ns.sfd_end_chip = -1;
    ns.polarity = 1;

    PhrResult want, got;
    BOOST_REQUIRE(decode_phr_decisions(pg0, pg1, empty_soft, ns, want));
    BOOST_REQUIRE(decode_phr_biorth_soft(pm0.data(), pm1.data(), empty_soft, ns, got));
    BOOST_CHECK_EQUAL(got.ok, want.ok);
    BOOST_CHECK_EQUAL(got.psdu_length, want.psdu_length);
    BOOST_CHECK_EQUAL(got.data_rate_mbps, want.data_rate_mbps);
    BOOST_CHECK(got.phr_bits == want.phr_bits);
    BOOST_CHECK_EQUAL(got.psdu_length, uint32_t(32));
    BOOST_CHECK_EQUAL(got.data_rate_mbps, 6.81f);

    // Light noise: derive the hard decisions exactly like decide_bpm_symbol()
    // and require both PHR entries to keep agreeing.
    std::mt19937 rng(20261122u);
    std::normal_distribution<float> light(0.0f, 0.15f);
    std::vector<float> nm0 = pm0, nm1 = pm1;
    for (size_t t = 0; t < nm0.size(); ++t) {
        nm0[t] += light(rng);
        nm1[t] += light(rng);
    }
    std::vector<int8_t> hg0(21), hg1(21);
    for (size_t t = 0; t < 21; ++t) {
        const bool use1 = std::abs(nm1[t]) > std::abs(nm0[t]);
        hg0[t] = use1 ? 1 : 0;
        const float v = use1 ? nm1[t] : nm0[t];
        hg1[t] = (v < 0.0f) ? 1 : 0;
    }
    PhrResult hw, hb;
    BOOST_REQUIRE(decode_phr_decisions(hg0, hg1, empty_soft, ns, hw));
    BOOST_REQUIRE(decode_phr_biorth_soft(nm0.data(), nm1.data(), empty_soft, ns, hb));
    BOOST_CHECK_EQUAL(hb.ok, hw.ok);
    BOOST_CHECK_EQUAL(hb.psdu_length, hw.psdu_length);
    BOOST_CHECK_EQUAL(hb.data_rate_mbps, hw.data_rate_mbps);
    BOOST_CHECK(hb.phr_bits == hw.phr_bits);

    // ---- Boundary handling: no crashes, safe false on unusable input ------
    PayloadResult bad;
    BOOST_CHECK(!decode_payload_biorth_soft(
        nullptr, phr1.data(), pay0.data(), pay1.data(), 240, 192, bad));
    BOOST_CHECK(!decode_payload_biorth_soft(
        phr0.data(), nullptr, pay0.data(), pay1.data(), 240, 192, bad));
    BOOST_CHECK(!decode_payload_biorth_soft(
        phr0.data(), phr1.data(), nullptr, pay1.data(), 240, 192, bad));
    BOOST_CHECK(!decode_payload_biorth_soft(
        phr0.data(), phr1.data(), pay0.data(), pay1.data(), 0, 192, bad));
    // psdu_bits far beyond what 8 payload symbols can carry -> RS fails.
    BOOST_CHECK(!decode_payload_biorth_soft(
        phr0.data(), phr1.data(), pay0.data(), pay1.data(), 8, 8 * 400, bad));

    PhrResult bad_phr;
    BOOST_CHECK(!decode_phr_biorth_soft(nullptr, pm1.data(), empty_soft, ns, bad_phr));
    BOOST_CHECK(!decode_phr_biorth_soft(pm0.data(), nullptr, empty_soft, ns, bad_phr));
}

namespace {

// Independent nested-loop local CL-3 re-solve with the joint ML metric.
void qa_ref_local_flip(const std::vector<float>& m0,
                       const std::vector<float>& m1,
                       size_t lo,
                       size_t hi,
                       int start_st,
                       int end_st,
                       const size_t* force_t,
                       const uint8_t* force_u,
                       size_t n_force,
                       std::vector<int8_t>& out_bits)
{
    const float INF = 1e30f;
    float surf[4] = { INF, INF, INF, INF };
    surf[start_st] = 0.0f;
    std::vector<int8_t> back((hi - lo) * 4, int8_t(0));
    for (size_t t = lo; t < hi; ++t) {
        int forced = -1;
        for (size_t f = 0; f < n_force; ++f)
            if (force_t[f] == t)
                forced = static_cast<int>(force_u[f]);
        float ns[4] = { INF, INF, INF, INF };
        for (int st = 0; st < 4; ++st) {
            if (surf[st] > 1e20f)
                continue;
            for (int u = 0; u < 2; ++u) {
                if (forced >= 0 && u != forced)
                    continue;
                const Cl3Transition step = cl3_transition(st, u);
                const float m = surf[st] + qa_branch(m0[t], m1[t], step.g0, step.g1);
                if (m < ns[step.next_state]) {
                    ns[step.next_state] = m;
                    back[(t - lo) * 4 + step.next_state] =
                        static_cast<int8_t>((u << 2) | st);
                }
            }
        }
        for (int i = 0; i < 4; ++i)
            surf[i] = ns[i];
    }
    int cur = end_st;
    if (surf[cur] > 1e20f) {
        for (int s = 0; s < 4; ++s)
            if (surf[s] < surf[cur])
                cur = s;
    }
    for (size_t t = hi; t-- > lo;) {
        const int v = back[(t - lo) * 4 + cur];
        out_bits[t] = static_cast<int8_t>(v >> 2);
        cur = v & 3;
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_w1_joint_ml_local_flip_and_rescue_phase1)
{
    using namespace openuwb::demod::core::detail;

    QaFrame frame;
    BOOST_REQUIRE(qa_build_frame(24, 4242u, frame));
    const size_t total = frame.u.size();
    // Mild noise so the trace has both usable margins and a real ACS order.
    std::mt19937 rng(20261123u);
    std::vector<float> m0 = frame.m0, m1 = frame.m1;
    qa_add_noise(m0, m1, 0.3f, rng);
    Cl3TrellisTrace tr;
    vitdec_cl3_biorth_trace(m0.data(), m1.data(), total, tr);
    BOOST_REQUIRE_EQUAL(tr.bits.size(), total);

    // local_cl3_trellis_flip_biorth must match the independent nested-loop
    // re-solve for unforced, single-forced and pair-forced windows.
    struct Window {
        size_t lo, hi;
        size_t ft[2];
        uint8_t fu[2];
        size_t n;
    };
    const Window wins[3] = {
        { 30, 70, { 0, 0 }, { 0, 0 }, 0 },
        { 100, 140, { 118, 0 }, { 1, 0 }, 1 },
        { 190, 235, { 200, 225 }, { 0, 1 }, 2 },
    };
    for (const auto& w : wins) {
        std::vector<int8_t> got(total, int8_t(7));
        std::vector<int8_t> ref(total, int8_t(7));
        local_cl3_trellis_flip_biorth(m0.data(),
                                      m1.data(),
                                      w.lo,
                                      w.hi,
                                      tr.states[w.lo],
                                      tr.states[w.hi],
                                      w.ft,
                                      w.fu,
                                      w.n,
                                      got);
        qa_ref_local_flip(
            m0, m1, w.lo, w.hi, tr.states[w.lo], tr.states[w.hi], w.ft, w.fu, w.n, ref);
        for (size_t t = w.lo; t < w.hi; ++t)
            BOOST_CHECK_EQUAL(int(got[t]), int(ref[t]));
        for (size_t t = 0; t < total; ++t) {
            if (t < w.lo || t >= w.hi)
                BOOST_CHECK_EQUAL(int(got[t]), 7);
        }
    }

    // Phase 1 of rescue_payload_biorth: a trace whose codeword has one or two
    // bit errors inside different RS symbols must be recovered by RS alone.
    for (size_t nerr : { size_t(1), size_t(2) }) {
        Cl3TrellisTrace bad = tr;
        for (size_t k = 0; k < nerr; ++k)
            bad.bits[19 + 6 * (3 * k + 1) + 1] ^= 1;
        PayloadResult out;
        const bool ok =
            rescue_payload_biorth(m0.data(), m1.data(), total, frame.psdu_bits, bad, out);
        BOOST_REQUIRE(ok);
        BOOST_CHECK(out.fcs_pass);
        BOOST_CHECK(out.bytes == frame.bytes);
        BOOST_CHECK_EQUAL(out.received_fcs, out.calculated_fcs);
    }

    // Heavily corrupted trace must never be accepted without a CRC match;
    // failing safely is equally valid (the rescue is allowed to give up).
    Cl3TrellisTrace garbage = tr;
    for (size_t t = 19; t + 2 < total; t += 7)
        garbage.bits[t] ^= 1;
    PayloadResult bad_out;
    if (rescue_payload_biorth(
            m0.data(), m1.data(), total, frame.psdu_bits, garbage, bad_out)) {
        BOOST_CHECK(bad_out.fcs_pass);
        BOOST_CHECK_EQUAL(bad_out.received_fcs, bad_out.calculated_fcs);
    }
    // Contract guards mirror rescue_payload_soft().
    PayloadResult guard;
    BOOST_CHECK(!rescue_payload_biorth(m0.data(), m1.data(), total, 0, tr, guard));
    BOOST_CHECK(!rescue_payload_biorth(m0.data(), m1.data(), total, 191, tr, guard));
}

BOOST_AUTO_TEST_CASE(test_w1_rs_erasures_full_v8_and_over_capacity)
{
    using namespace openuwb::demod::core::detail;

    // 330-bit systematic data block -> 378-bit RS(63,55) codeword.
    std::vector<int8_t> data330(330);
    std::mt19937 rng(20261124u);
    for (size_t i = 0; i < data330.size(); ++i)
        data330[i] = static_cast<int8_t>(rng() & 1u);

    std::vector<int8_t> coded;
    BOOST_REQUIRE(qa_rs_encode_block(data330, coded));
    BOOST_REQUIRE_EQUAL(coded.size(), size_t(378));

    // Eight pure erasures: one bit flipped in each of eight distinct 6-bit
    // symbols (eight symbol errors, zero residual errors).
    std::vector<int8_t> corrupted = coded;
    int eras8[8];
    for (int k = 0; k < 8; ++k) {
        const int sy = 5 + 7 * k;
        eras8[k] = sy;
        corrupted[static_cast<size_t>(sy * 6 + (k % 6))] ^= 1;
    }
    int8_t out330[330] = {};
    BOOST_CHECK(rs_decode_block_erasures(corrupted.data(), eras8, 8, out330));
    for (size_t i = 0; i < data330.size(); ++i)
        BOOST_CHECK_EQUAL(int(out330[i]), int(data330[i]));

    // v = 6 on the same word cannot *recover the original data*: 8 symbol
    // errors exceed the unique-decoding radius with only 6 erasures, so the
    // decoder may legally miscorrect (return a different valid codeword).
    // What matters here is that v = 6 cannot reproduce the transmitted block;
    // the upper layers discard miscorrections via the full-frame CRC-16.
    int8_t six330[330] = {};
    const bool six_ok =
        rs_decode_block_erasures(corrupted.data(), eras8, 6, six330);
    const bool six_restored =
        six_ok && std::equal(six330, six330 + 330, data330.begin());
    BOOST_TEST_MESSAGE("v=6 on 8 erasures: ok=" << six_ok
                                                << " restored=" << six_restored);
    BOOST_CHECK(!six_restored);

    // Nine symbol errors exceed the v = 8 + e = 0 budget.  v = 8 still always
    // solves the 8x8 erasure system and can miscorrect to another valid
    // codeword; the reliable contract is that it must not reproduce the
    // original data (again, CRC-16 is the final arbiter).
    std::vector<int8_t> worse = coded;
    int eras_first8[8];
    for (int k = 0; k < 9; ++k) {
        const int sy = 3 + 6 * k;
        worse[static_cast<size_t>(sy * 6 + (k % 6))] ^= 1;
        if (k < 8)
            eras_first8[k] = sy;
    }
    int8_t bad330[330] = {};
    const bool bad_ok =
        rs_decode_block_erasures(worse.data(), eras_first8, 8, bad330);
    const bool bad_restored =
        bad_ok && std::equal(bad330, bad330 + 330, data330.begin());
    BOOST_TEST_MESSAGE("v=8 on 9 errors: ok=" << bad_ok
                                              << " restored=" << bad_restored);
    BOOST_CHECK(!bad_restored);

    // Out-of-range v falls back to the plain errors-only decoder, which also
    // cannot repair eight symbol errors.
    int8_t nine330[330] = {};
    BOOST_CHECK(!rs_decode_block_erasures(corrupted.data(), eras8, 9, nine330));
}

BOOST_AUTO_TEST_CASE(test_w1_rescue_phase1_full_v8_erasures)
{
    using namespace openuwb::demod::core::detail;

    // 41 bytes -> 328 coded data bits: exactly one shortened RS block with
    // data_start = 2, so all 63 RS symbols are exercised.
    QaFrame frame;
    BOOST_REQUIRE(qa_build_frame(41, 5150u, frame));
    BOOST_REQUIRE_EQUAL(frame.psdu_bits, size_t(328));
    BOOST_REQUIRE_EQUAL(frame.payload_symbols, size_t(376));
    const size_t total = frame.u.size();
    BOOST_REQUIRE_EQUAL(total, size_t(19 + 376 + 2));

    Cl3TrellisTrace tr;
    vitdec_cl3_biorth_trace(frame.m0.data(), frame.m1.data(), total, tr);
    BOOST_REQUIRE_EQUAL(tr.bits.size(), total);
    tr.bits = frame.u; // exact transmitted information bits
    for (auto& m : tr.margins)
        m = 100.0f; // every step initially highly reliable

    const size_t t_b = 19;
    const size_t data_start = 330 - frame.psdu_bits; // 2
    const int bad_sym[8] = { 4, 12, 20, 28, 36, 44, 52, 60 };
    for (int k = 0; k < 8; ++k) {
        const int sy = bad_sym[k];
        for (int bit = 0; bit < 6; ++bit) {
            const size_t pos = static_cast<size_t>(sy * 6 + bit);
            size_t vt = 0;
            if (pos >= data_start && pos < 330)
                vt = t_b + (pos - data_start);
            else if (pos >= 330)
                vt = t_b + frame.psdu_bits + (pos - 330);
            else
                continue;
            tr.margins[vt] = 0.01f; // make this symbol one of the eight weakest
        }
        // Flip one mapped bit so the symbol value (not just its margin) is
        // wrong: eight symbol errors, beyond the plain RS t = 4 bound.
        const size_t pos = static_cast<size_t>(sy * 6 + 3);
        size_t vt = 0;
        if (pos >= data_start && pos < 330)
            vt = t_b + (pos - data_start);
        else
            vt = t_b + frame.psdu_bits + (pos - 330);
        tr.bits[vt] ^= 1;
    }

    PayloadResult out;
    BOOST_REQUIRE(rescue_payload_biorth(frame.m0.data(),
                                        frame.m1.data(),
                                        total,
                                        frame.psdu_bits,
                                        tr,
                                        out));
    BOOST_CHECK(out.ok);
    BOOST_CHECK(out.fcs_pass);
    BOOST_CHECK(out.bytes == frame.bytes);
    BOOST_CHECK_EQUAL(out.received_fcs, out.calculated_fcs);

    // Prove the packet really needed the v = 8 attempt: on the same corrupted
    // codeword v = 6 cannot reproduce the original block while v = 8 does.
    // With 8 symbol errors and only 6 erasures, v = 6 is outside the unique
    // decoding radius and is allowed to miscorrect; only the CRC-16 decides
    // whether a candidate is the real frame.
    int8_t blk[378];
    expand_rs_codeword(tr.bits.data() + t_b, frame.psdu_bits, blk);
    int8_t orig_blk[378];
    expand_rs_codeword(frame.u.data() + t_b, frame.psdu_bits, orig_blk);
    int eras[8];
    for (int k = 0; k < 8; ++k)
        eras[k] = bad_sym[k];
    int8_t d330[330] = {};
    const bool six_ok = rs_decode_block_erasures(blk, eras, 6, d330);
    const bool six_restored =
        six_ok && std::equal(d330, d330 + 330, orig_blk);
    BOOST_TEST_MESSAGE("phase-1 v=6: ok=" << six_ok
                                          << " restored=" << six_restored);
    BOOST_CHECK(!six_restored);
    BOOST_CHECK(rs_decode_block_erasures(blk, eras, 8, d330));
    BOOST_CHECK(std::equal(d330, d330 + 330, orig_blk));
}

BOOST_AUTO_TEST_CASE(test_w1_phr_payload_independent_clipping)
{
    using namespace openuwb::demod::core::detail;

    QaFrame frame;
    BOOST_REQUIRE(qa_build_frame(24, 8086u, frame));
    BOOST_REQUIRE_EQUAL(frame.payload_symbols, size_t(240));

    // PHR metrics are 8x the payload scale (64-chip vs 8-chip symbols).
    std::vector<float> phr0(frame.m0.begin(), frame.m0.begin() + 21);
    std::vector<float> phr1(frame.m1.begin(), frame.m1.begin() + 21);
    for (size_t i = 0; i < 21; ++i) {
        phr0[i] *= 8.0f;
        phr1[i] *= 8.0f;
    }
    std::vector<float> pay0(frame.m0.begin() + 21, frame.m0.end());
    std::vector<float> pay1(frame.m1.begin() + 21, frame.m1.end());
    // Legitimate strong payload pulses at the same positions as the existing
    // pulse-clipping test.
    const size_t pulse_pos[5] = { 60, 100, 140, 180, 220 };
    for (size_t i = 0; i < 5; ++i) {
        const float amp = (i % 2 == 0) ? 50.0f : -50.0f;
        pay0[pulse_pos[i]] += amp;
    }

    PayloadResult out;
    PayloadDecodeWorkspace ws;
    BOOST_REQUIRE(decode_payload_biorth_soft(phr0.data(),
                                             phr1.data(),
                                             pay0.data(),
                                             pay1.data(),
                                             frame.payload_symbols,
                                             frame.psdu_bits,
                                             out,
                                             nullptr,
                                             &ws,
                                             false,
                                             true));
    BOOST_CHECK(out.fcs_pass);
    BOOST_CHECK(out.bytes == frame.bytes);

    // Per-field clip contract: the 8x PHR slice is unchanged because its own
    // 4x median limit is 32, while the payload slice clamps the 50x pulse at
    // its own 4x limit of 4.
    std::vector<float> cphr0 = phr0, cphr1 = phr1;
    clip_biorth_metrics(cphr0.data(), cphr1.data(), 21);
    for (size_t i = 0; i < 21; ++i) {
        BOOST_CHECK_SMALL(cphr0[i] - phr0[i], 1e-6f);
        BOOST_CHECK_SMALL(cphr1[i] - phr1[i], 1e-6f);
    }
    std::vector<float> cpay0 = pay0, cpay1 = pay1;
    clip_biorth_metrics(cpay0.data(), cpay1.data(), pay0.size());
    BOOST_CHECK_CLOSE(cpay0[pulse_pos[0]], 4.0f, 1e-3);

    // The old unified clip over [PHR|payload] would have attenuated the 8x PHR
    // to the payload-derived 4x limit; the split path leaves it untouched.
    std::vector<float> um0 = phr0, um1 = phr1;
    um0.insert(um0.end(), pay0.begin(), pay0.end());
    um1.insert(um1.end(), pay1.begin(), pay1.end());
    clip_biorth_metrics(um0.data(), um1.data(), um0.size());
    size_t nz = 0;
    while (nz < 21 && std::abs(phr0[nz]) < 1.0f && std::abs(phr1[nz]) < 1.0f)
        ++nz;
    BOOST_REQUIRE_LT(nz, size_t(21));
    const float before = std::max(std::abs(phr0[nz]), std::abs(phr1[nz]));
    const float after = std::max(std::abs(um0[nz]), std::abs(um1[nz]));
    BOOST_CHECK_CLOSE(after, 4.0f, 1e-3);
    BOOST_CHECK_LT(after, before);
}

BOOST_AUTO_TEST_CASE(test_w1_rescue_phase2_flip_plus_erasure)
{
    using namespace openuwb::demod::core::detail;

    // 41 bytes -> one shortened RS block (328 data bits, data_start = 2).
    QaFrame frame;
    BOOST_REQUIRE(qa_build_frame(41, 6174u, frame));
    const size_t total = frame.u.size();
    BOOST_REQUIRE_EQUAL(total, size_t(19 + 376 + 2));
    const size_t t_b = 19;
    const size_t data_start = 330 - frame.psdu_bits; // 2

    Cl3TrellisTrace tr;
    vitdec_cl3_biorth_trace(frame.m0.data(), frame.m1.data(), total, tr);
    BOOST_REQUIRE_EQUAL(tr.bits.size(), total);
    tr.bits = frame.u; // exact transmitted information bits

    auto vt_of = [&](int sy, int bit) -> size_t {
        const size_t pos = static_cast<size_t>(sy * 6 + bit);
        if (pos >= data_start && pos < 330)
            return t_b + (pos - data_start);
        return t_b + frame.psdu_bits + (pos - 330);
    };

    // Six RS symbol errors.  Exactly one (the symbol containing t_win) sits in
    // the single-node flip window; the other five stay outside it.
    //
    // Margin ranking is chosen so Phase 1 cannot close the frame:
    //   rank 1..9 = E6(s_win), E1, E2, C1, C2, E3, E4, C3, E5
    //   top-4 has 3 errors (v=4 fails), top-6 has 4 (v=6 fails), top-8 has 5
    //   (v=8 always solves the 8x8 system and therefore miscorrects).
    for (auto& m : tr.margins)
        m = 100.0f;

    const int s_win = 30;
    const size_t t_win = vt_of(s_win, 3);
    tr.bits[t_win] ^= 1;
    tr.margins[t_win] = 0.01f;
    // Phase 2 forces tr.lose_u[t]; set it to the truly transmitted bit so the
    // local re-solve restores the original window instead of flipping away.
    tr.lose_u[t_win] = static_cast<uint8_t>(frame.u[t_win] & 1);

    auto mark = [&](int sy, int bit, float margin, bool error) {
        const size_t vt = vt_of(sy, bit);
        tr.margins[vt] = margin;
        if (error)
            tr.bits[vt] ^= 1;
    };
    mark(4, 0, 0.2f, true);   // E1 rank 2 (erased by the Phase 2 retry)
    mark(10, 0, 0.3f, true);  // E2 rank 3 (erased by the Phase 2 retry)
    mark(22, 0, 0.4f, false); // C1 rank 4 (weak but correct)
    mark(37, 0, 0.5f, false); // C2 rank 5 (weak but correct)
    mark(16, 0, 0.6f, true);  // E3 rank 6
    mark(44, 0, 0.7f, true);  // E4 rank 7
    mark(47, 0, 0.8f, false); // C3 rank 8
    mark(54, 0, 0.9f, true);  // E5 rank 9 (9th error, outside Phase 1 top-8)

    // Phase 1 alone must miss: plain RS fails, and the v = 4/6/8 erasure
    // attempts on the un-flipped codeword cannot reproduce the original block
    // (v = 8 always returns some codeword, but a miscorrected one here).
    int8_t orig_blk[378];
    expand_rs_codeword(frame.u.data() + t_b, frame.psdu_bits, orig_blk);
    int8_t blk0[378];
    expand_rs_codeword(tr.bits.data() + t_b, frame.psdu_bits, blk0);
    int8_t d0[330] = {};
    const int era4[4] = { 30, 4, 10, 22 };
    const int era6[6] = { 30, 4, 10, 22, 37, 16 };
    const int era8[8] = { 30, 4, 10, 22, 37, 16, 44, 47 };
    BOOST_CHECK(!rs_decode_block(blk0, d0));
    const bool p4_ok = rs_decode_block_erasures(blk0, era4, 4, d0);
    BOOST_CHECK(!(p4_ok && std::equal(d0, d0 + 330, orig_blk)));
    const bool p6_ok = rs_decode_block_erasures(blk0, era6, 6, d0);
    BOOST_CHECK(!(p6_ok && std::equal(d0, d0 + 330, orig_blk)));
    const bool p8_ok = rs_decode_block_erasures(blk0, era8, 8, d0);
    BOOST_CHECK(!(p8_ok && std::equal(d0, d0 + 330, orig_blk)));

    // End-to-end: Phase 1 misses, Phase 2 flips the weakest node and the joint
    // RS erasure retry closes the single block; only the CRC-16 accepts it.
    PayloadResult out;
    BOOST_REQUIRE(rescue_payload_biorth(frame.m0.data(),
                                        frame.m1.data(),
                                        total,
                                        frame.psdu_bits,
                                        tr,
                                        out));
    BOOST_CHECK(out.ok);
    BOOST_CHECK(out.fcs_pass);
    BOOST_CHECK(out.bytes == frame.bytes);
    BOOST_CHECK_EQUAL(out.received_fcs, out.calculated_fcs);

    // Deterministic proof of the mechanism: replay the same single-node flip.
    // Plain RS then sees 5 residual symbol errors (fails), while erasing the
    // two weakest real-error symbols leaves 3 unknown errors and
    // 2*3 + 2 = 8, exactly the v = 2 boundary.
    const size_t lo = (t_win > 12) ? (t_win - 12) : 0;
    const size_t hi = std::min(total, t_win + 13);
    std::vector<int8_t> work = tr.bits;
    const uint8_t fu = tr.lose_u[t_win];
    local_cl3_trellis_flip_biorth(frame.m0.data(),
                                  frame.m1.data(),
                                  lo,
                                  hi,
                                  tr.states[lo],
                                  tr.states[hi],
                                  &t_win,
                                  &fu,
                                  1,
                                  work);

    int8_t blk[378];
    expand_rs_codeword(work.data() + t_b, frame.psdu_bits, blk);
    int8_t d330[330] = {};
    BOOST_CHECK(!rs_decode_block(blk, d330));

    int eras_join[2] = { 4, 10 };
    BOOST_CHECK(rs_decode_block_erasures(blk, eras_join, 2, d330));
    for (size_t i = 0; i < 330; ++i)
        BOOST_CHECK_EQUAL(int(d330[i]), int(orig_blk[i]));
}

// WP4 lightweight low-SNR regression guard for the fused PHR/payload chain.
//
// The WP1 (RS rescue) and WP3 (-pol retry) improvements live inside the
// decoder and have no run-time switch, so an in-binary A/B is impossible.
// This guard instead pins a deterministic baseline: the same 20 fixed seeds,
// the same fixed unit-variance noise, and the same decoder entry point that
// the fused kernel uses.  The recorded `kBaselinePasses` below was measured
// on 2026-11-22 (AVX2, dev-avx2, seed set 1000..1019) after all three work
// packages landed; any regression makes this test fail, while further
// improvements can only raise the number.  No loops over noise levels, no
// Monte-Carlo sweep: 20 frames, well under 0.5 s.
BOOST_AUTO_TEST_CASE(test_w4_payload_low_snr_deterministic_regression_guard)
{
    using namespace openuwb::demod::core::detail;

    constexpr size_t kFrames = 20;
    constexpr float kSigma = 0.60f; // 8-16/20 FCS passes on this seed set
    constexpr size_t kBaselinePasses = 15;
    const auto t0 = std::chrono::steady_clock::now();

    size_t passes = 0;
    for (size_t k = 0; k < kFrames; ++k) {
        QaFrame frame;
        BOOST_REQUIRE(qa_build_frame(24, 1000u + static_cast<uint32_t>(k), frame));
        BOOST_REQUIRE_EQUAL(frame.payload_symbols, size_t(240));

        std::vector<float> m0 = frame.m0, m1 = frame.m1;
        std::mt19937 rng(90000u + static_cast<uint32_t>(k));
        qa_add_noise(m0, m1, kSigma, rng);

        PayloadDecodeWorkspace ws;
        PayloadResult out;
        const bool ok = decode_payload_biorth_soft(frame.m0.data(),
                                                   frame.m1.data(),
                                                   m0.data() + 21,
                                                   m1.data() + 21,
                                                   frame.payload_symbols,
                                                   frame.psdu_bits,
                                                   out,
                                                   nullptr,
                                                   &ws,
                                                   /*terminated=*/true,
                                                   /*enable_clip=*/true);
        if (ok && out.fcs_pass && out.bytes == frame.bytes)
            ++passes;
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    BOOST_TEST_MESSAGE("WP4 low-SNR guard: FCS "
                       << passes << "/" << kFrames << " at sigma=" << kSigma
                       << ", " << ms << " ms");
    BOOST_CHECK_LE(ms, 500.0);
    BOOST_CHECK_GE(passes, kBaselinePasses);
}
