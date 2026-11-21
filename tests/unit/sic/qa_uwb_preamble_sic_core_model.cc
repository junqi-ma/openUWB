// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 2/3 of the preamble-SIC core QA: the S1 eight-symbol model guard
// cases (extrapolation, bypass, guard page, golden cross-check).
#include "preamble_sic_test_helpers.h"

// ---------------------------------------------------------------------------
// S1-7. Extrapolation mismatch: the model must NOT just "commit"; the
// over-subtraction and the target damage must be MEASURED.
// ---------------------------------------------------------------------------
namespace {
struct ExtrapMetrics {
    double over_sub_db = -300.0; // i_hat energy where i_true is absent
    double damage_db = -300.0;   // ||out - (y - i_true)||^2 / ||target||^2
    double cancel_db = -300.0;   // i_true removed inside the active region
};

ExtrapMetrics measure_extrap(const Signal& s,
                             const std::vector<cplx>& target,
                             const std::vector<cplx>& ihat,
                             const std::vector<cplx>& out,
                             size_t lo,
                             size_t hi)
{
    ExtrapMetrics m;
    double i_mean = region_energy(s.i, lo, hi) / static_cast<double>(hi - lo);
    double e_ih_off = 0.0, e_i_on = 0.0, e_ih_on = 0.0;
    for (size_t k = lo; k < hi; ++k) {
        if (std::norm(s.i[k]) < 1e-6 * i_mean) {
            e_ih_off += std::norm(ihat[k]);
        } else {
            e_i_on += std::norm(s.i[k]);
            e_ih_on += std::norm(ihat[k]);
        }
    }
    double e_tgt = 0.0, e_dmg = 0.0;
    for (size_t k = lo; k < hi; ++k) {
        e_tgt += std::norm(target[k]);
        e_dmg += std::norm(out[k] - (s.y[k] - s.i[k]));
    }
    m.over_sub_db = 10.0 * std::log10((e_ih_off + 1e-30) / (e_i_on + 1e-30));
    m.damage_db = 10.0 * std::log10((e_dmg + 1e-30) / (e_tgt + 1e-30));
    m.cancel_db = 10.0 * std::log10((e_i_on + 1e-30) / (e_dmg + 1e-30));
    (void)e_ih_on;
    return m;
}
} // namespace

BOOST_AUTO_TEST_CASE(test_sic_s1_extrapolation_mismatch)
{
    const size_t reps = 64;
    const int64_t oe_off = 8 * static_cast<int64_t>(kSps);
    struct Case {
        const char* name;
        // 0 = matched, 1 = short preamble, 2 = burst interruption,
        // 3 = phase jump, 4 = time-varying channel gain.
        int kind;
    };
    const Case cases[] = { { "matched", 0 },
                           { "short_preamble", 1 },
                           { "burst", 2 },
                           { "phase_jump", 3 },
                           { "time_varying_gain", 4 } };
    double matched_cancel = 0.0;
    for (const Case& c : cases) {
        Fixture f =
            make_fixture(10, reps, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
        const size_t gen_reps = (c.kind == 1) ? 20 : reps;
        Model gen = build_model(10, gen_reps);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(gen,
                                f.n,
                                f.comm_start,
                                paths,
                                0.0,
                                kQaNoise,
                                500 + static_cast<uint64_t>(c.kind));
        // Deterministic target (sensing stand-in) at ISR = 10 dB.
        std::vector<cplx> target(s.y.size(), cplx(0.f, 0.f));
        {
            std::mt19937_64 rng(900 + static_cast<uint64_t>(c.kind));
            std::normal_distribution<double> nd(0.0, 1.0);
            double pi = 0.0;
            for (size_t k = 0; k < s.y.size(); ++k) {
                cplx v(static_cast<float>(nd(rng)), static_cast<float>(nd(rng)));
                target[k] = v;
                pi += std::norm(v);
            }
            const double rms_i =
                std::sqrt(region_energy(s.i, 0, f.n) / static_cast<double>(f.n));
            const double rms_t = std::sqrt(pi / static_cast<double>(f.n));
            const double scale =
                (rms_t > 0.0) ? rms_i / std::pow(10.0, 10.0 / 20.0) / rms_t : 0.0;
            for (auto& v : target)
                v = cplx(static_cast<float>(v.real() * scale),
                         static_cast<float>(v.imag() * scale));
        }
        // Rebuild y = target + modified interferer + noise.
        std::vector<cplx> i_mod = s.i;
        switch (c.kind) {
        case 2: // burst interruption: the interferer stops after 32 symbols
            for (size_t k = static_cast<size_t>(f.comm_start) + 32 * kSps;
                 k < i_mod.size();
                 ++k)
                i_mod[k] = cplx(0.f, 0.f);
            break;
        case 3: // phase jump of 1.0 rad after 32 symbols
            for (size_t k = static_cast<size_t>(f.comm_start) + 32 * kSps;
                 k < i_mod.size();
                 ++k)
                i_mod[k] = i_mod[k] * cplx(static_cast<float>(std::cos(1.0)),
                                           static_cast<float>(std::sin(1.0)));
            break;
        case 4: // time-varying channel gain: amplitude + phase drift
            for (size_t k = 0; k < i_mod.size(); ++k) {
                const double x = static_cast<double>(k) / static_cast<double>(f.n);
                const double amp = 1.0 + 0.4 * std::sin(2.0 * kPi * 2.0 * x);
                const double ph = 1.2 * kPi * x;
                i_mod[k] = i_mod[k] * cplx(static_cast<float>(amp * std::cos(ph)),
                                           static_cast<float>(amp * std::sin(ph)));
            }
            break;
        default:
            break;
        }
        {
            std::mt19937_64 rng(777 + static_cast<uint64_t>(c.kind));
            std::normal_distribution<double> nd(0.0, kQaNoise);
            for (size_t k = 0; k < f.n; ++k) {
                s.i[k] = i_mod[k];
                s.y[k] = i_mod[k] + target[k] +
                         cplx(static_cast<float>(nd(rng)), static_cast<float>(nd(rng)));
            }
        }
        const int64_t oe = f.comm_start + oe_off;
        f.cfg.observed_end_sample = oe;
        SicResult r = run(f, s);
        const size_t lo = static_cast<size_t>(oe);
        const size_t hi = static_cast<size_t>(f.comm_start) + reps * kSps;
        const ExtrapMetrics mt = measure_extrap(s, target, f.ihat, f.out, lo, hi);
        std::printf("[sic-qa] s1-extrap '%-18s' state=%s fit_res=%.4f | "
                    "cancel=%.2f dB damage=%.2f dB over_sub=%.2f dB\n",
                    c.name,
                    state_name(r.state),
                    r.fit_residual,
                    mt.cancel_db,
                    mt.damage_db,
                    mt.over_sub_db);
        if (c.kind == 0) {
            BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                              static_cast<int>(SicState::Committed));
            BOOST_CHECK_LT(mt.damage_db, -10.0);   // target essentially intact
            BOOST_CHECK_LT(mt.over_sub_db, -10.0); // nothing subtracted where i=0
            BOOST_CHECK_GT(mt.cancel_db, 20.0);
            matched_cancel = mt.cancel_db;
        } else {
            // The point of the case: the mismatch is MEASURED.  Over-subtraction
            // energy and target damage are real numbers, not a "detected" flag.
            // Every mismatch degrades cancellation relative to the matched case
            // and leaves measurable target damage.
            BOOST_CHECK_GT(mt.damage_db, -6.0); // >= 25% of target energy hit
            BOOST_CHECK_LT(mt.cancel_db, matched_cancel - 20.0);
            if (c.kind == 1 || c.kind == 2) // phantom region really exists
                BOOST_CHECK_GT(mt.over_sub_db, -20.0);
        }
    }
}

// ---------------------------------------------------------------------------
// S1-8. Bypass bit-exactness for an UNTRUSTED model and for save-region
// mismatches; and apply_model with a model that was never committed.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_s1_bypass_untrusted_model)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 601);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;

    // (a) A withdrawn/untrusted model (valid=false) must not touch the output.
    {
        SicModel m;
        m.valid = false;
        std::vector<cplx> out(s.y.size(), cplx(0.f, 0.f));
        std::vector<cplx> ih(s.y.size(), cplx(0.f, 0.f));
        const bool ok = sic::apply_model(m,
                                         s.y.data(),
                                         s.y.size(),
                                         0,
                                         f.templates(),
                                         f.cfg,
                                         f.scratch,
                                         out.data(),
                                         ih.data(),
                                         nullptr);
        BOOST_CHECK(ok);
        BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), s.y.begin(), s.y.end());
        for (size_t k = 0; k < ih.size(); ++k)
            BOOST_CHECK_EQUAL(ih[k], cplx(0.f, 0.f));
    }

    // (b) No interference at all: detection must not commit and out == y.
    {
        Signal n0;
        n0.y.assign(f.n, cplx(0.f, 0.f));
        n0.i.assign(f.n, cplx(0.f, 0.f));
        std::mt19937_64 rng(602);
        std::normal_distribution<double> nd(0.0, kQaNoise);
        for (size_t k = 0; k < f.n; ++k)
            n0.y[k] = cplx(static_cast<float>(nd(rng)), static_cast<float>(nd(rng)));
        SicResult r = run(f, n0);
        BOOST_CHECK_NE(static_cast<int>(r.state), static_cast<int>(SicState::Committed));
        BOOST_CHECK_EQUAL_COLLECTIONS(
            f.out.begin(), f.out.end(), n0.y.begin(), n0.y.end());
    }

    // (c) Committed model but a save region entirely outside it: the persistence
    //     of the bypass is guaranteed by the model's own sub region, checked by
    //     applying where the model never overlaps.
    {
        SicModel m = estimate_model_checked(f, s, oe, nullptr);
        BOOST_REQUIRE(m.valid);
        std::vector<cplx> out(512, cplx(0.f, 0.f));
        std::vector<cplx> ih(512, cplx(0.f, 0.f));
        const int64_t off = static_cast<int64_t>(f.n) - 512; // after the mask
        BOOST_REQUIRE(sic::apply_model(m,
                                       s.y.data() + off,
                                       512,
                                       off,
                                       f.templates(),
                                       f.cfg,
                                       f.scratch,
                                       out.data(),
                                       ih.data(),
                                       nullptr));
        for (size_t k = 0; k < 512; ++k)
            BOOST_CHECK_EQUAL(out[k], s.y[static_cast<size_t>(off) + k]);
    }

    // (d) sub_roi disjoint from the predicted support -> bypass, out == y.
    {
        Fixture g =
            make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
        const int64_t sensing_hi = g.comm_start - 2000;
        g.cfg.sub_roi_lo = 0;
        g.cfg.sub_roi_hi = std::max<int64_t>(1, sensing_hi);
        Signal t = build_signal(g.model, g.n, g.comm_start, paths, 0.0, kQaNoise, 603);
        g.cfg.observed_end_sample = g.comm_start + 8 * static_cast<int64_t>(kSps);
        SicResult r = run(g, t);
        std::printf("[sic-qa] s1-bypass: disjoint-sub-roi state=%s reason=%s "
                    "sub_overlap=%lld\n",
                    state_name(r.state),
                    r.bypass_reason ? r.bypass_reason : "-",
                    static_cast<long long>(r.sub_overlap));
        BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                          static_cast<int>(SicState::Bypassed));
        BOOST_CHECK_EQUAL(std::string(r.bypass_reason ? r.bypass_reason : ""),
                          std::string("outside_sensing_roi"));
        BOOST_CHECK_EQUAL_COLLECTIONS(g.out.begin(), g.out.end(), t.y.begin(), t.y.end());
    }
}

// ---------------------------------------------------------------------------
// S1-10. A DEGENERATE / NON-FINITE model must bypass bit-exactly.
//
// apply_model is the only consumer of a frozen model.  A model that is
// `valid` but carries a NaN/Inf gain or delay (a corrupted freeze, or a model
// the estimator never actually committed) must never write a non-finite sample
// into `out`: the whole chunk stays byte-for-byte identical to `y`, exactly
// like the untrusted/invalid-model and sub-ROI bypass paths.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_s1_degenerate_nonfinite_model_bypass)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 611);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;

    SicModel good = estimate_model_checked(f, s, oe, nullptr);
    BOOST_REQUIRE(good.valid);
    BOOST_REQUIRE(good.n_paths > 0);

    const auto apply_expect_bypass = [&](const SicModel& m, const char* what) {
        std::vector<cplx> out(s.y.size(), cplx(9.f, 9.f));
        std::vector<cplx> ih(s.y.size(), cplx(9.f, 9.f));
        const bool ok = sic::apply_model(m,
                                         s.y.data(),
                                         s.y.size(),
                                         0,
                                         f.templates(),
                                         f.cfg,
                                         f.scratch,
                                         out.data(),
                                         ih.data(),
                                         nullptr);
        BOOST_CHECK_MESSAGE(ok, what << ": apply_model returned false");
        for (size_t k = 0; k < out.size(); ++k)
            BOOST_CHECK_EQUAL(out[k], s.y[k]);
        for (size_t k = 0; k < ih.size(); ++k)
            BOOST_CHECK_EQUAL(ih[k], cplx(0.f, 0.f));
    };

    // (a) valid but a NaN complex gain.
    {
        SicModel m = good;
        m.gain[0] = cplx(std::numeric_limits<float>::quiet_NaN(), 0.f);
        apply_expect_bypass(m, "nan_gain");
    }
    // (b) valid but an Inf delay.
    {
        SicModel m = good;
        m.tau_abs[0] = std::numeric_limits<double>::infinity();
        apply_expect_bypass(m, "inf_tau");
    }
    // (c) valid but a NaN CFO.
    {
        SicModel m = good;
        m.cfo_hz = std::numeric_limits<double>::quiet_NaN();
        apply_expect_bypass(m, "nan_cfo");
    }
    // (d) never-committed (valid=false) model.
    {
        SicModel m = good;
        m.valid = false;
        apply_expect_bypass(m, "invalid");
    }

    // Control: the untouched committed model DOES subtract, so the bypass
    // assertions above are not vacuous.
    {
        std::vector<cplx> out(s.y.size(), cplx(0.f, 0.f));
        std::vector<cplx> ih(s.y.size(), cplx(0.f, 0.f));
        BOOST_REQUIRE(sic::apply_model(good,
                                       s.y.data(),
                                       s.y.size(),
                                       0,
                                       f.templates(),
                                       f.cfg,
                                       f.scratch,
                                       out.data(),
                                       ih.data(),
                                       nullptr));
        bool changed = false;
        for (size_t k = 0; k < out.size(); ++k) {
            if (out[k] != s.y[k]) {
                changed = true;
                break;
            }
        }
        BOOST_CHECK_MESSAGE(changed, "committed control model did not subtract");
    }
}

// ---------------------------------------------------------------------------
// S1-9. Strict-prefix memory guard: the estimation must not read at or past the
// declared visibility limit (a guard page catches a vectorised overread beyond a
// 64-sample slack; the mutation test above catches any read within the slack).
// ---------------------------------------------------------------------------
#if defined(__linux__)
BOOST_AUTO_TEST_CASE(test_sic_s1_prefix_guard_page)
{
    Fixture f = make_fixture(10, 64, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 700);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;

    SicModel ref = estimate_model_checked(f, s, oe, nullptr);

    const size_t slack = 64; // readable samples past oe
    const size_t payload = static_cast<size_t>(oe) + slack;
    const size_t page = 4096;
    const size_t bytes = payload * sizeof(cplx);
    const size_t map = ((bytes + page - 1) / page + 1) * page;
    void* p =
        ::mmap(nullptr, map, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    BOOST_REQUIRE(p != MAP_FAILED);
    uint8_t* base = static_cast<uint8_t*>(p);
    const size_t data_off = map - page - bytes;
    cplx* yp = reinterpret_cast<cplx*>(base + data_off);
    for (size_t k = 0; k < payload; ++k)
        yp[k] = s.y[k];
    BOOST_REQUIRE_EQUAL(::mprotect(base + map - page, page, PROT_NONE), 0);

    SicModel m;
    SicResult r;
    const bool ok =
        sic::estimate_model(yp, f.n, oe, f.templates(), f.cfg, f.scratch, m, nullptr, &r);
    BOOST_CHECK(ok);
    check_models_identical(ref, m, "guard-page");
    std::printf("[sic-qa] s1-prefix-guard-page: observed_end=%lld slack=%zu "
                "PROT_NONE after payload -> estimate identical\n",
                static_cast<long long>(oe),
                slack);
    ::munmap(p, map);
}
#endif

// ---------------------------------------------------------------------------
// 18. golden cross-check (optional; skipped when the directory is absent)
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_sic_core_golden_crosscheck)
{
    const std::string root = golden_root();
    if (root.empty() || !dir_exists(root)) {
        std::printf("[sic-qa] SKIP golden cross-check: '%s' not present\n", root.c_str());
        return;
    }
    DIR* d = ::opendir(root.c_str());
    if (!d) {
        std::printf("[sic-qa] SKIP golden cross-check: cannot open '%s'\n", root.c_str());
        return;
    }
    size_t cases = 0;
    while (struct dirent* e = ::readdir(d)) {
        const std::string name = e->d_name;
        if (name == "." || name == "..")
            continue;
        const std::string cdir = root + "/" + name;
        if (!dir_exists(cdir))
            continue;
        const std::string npz_path = cdir + "/io.npz";
        const std::string meta_path = cdir + "/meta.json";
        std::vector<uint8_t> z;
        std::vector<uint8_t> mt;
        if (!read_file(npz_path, z) || !read_file(meta_path, mt)) {
            std::printf("[sic-qa] SKIP golden case '%s': io.npz/meta.json missing\n",
                        name.c_str());
            continue;
        }
        const std::string meta(reinterpret_cast<const char*>(mt.data()), mt.size());
        std::vector<cplx> y, sync, preamble, ihat_exp, out_exp, i_true;
        auto load = [&](const char* key, std::vector<cplx>& dst) -> bool {
            std::vector<uint8_t> b;
            if (!npz_get(z, std::string(key) + ".npy", b))
                return false;
            return npy_complex(b, dst);
        };
        BOOST_REQUIRE_MESSAGE(load("y", y), name << ": cannot load y from io.npz");
        BOOST_REQUIRE_MESSAGE(load("sync", sync), name << ": cannot load sync");
        BOOST_REQUIRE_MESSAGE(load("preamble", preamble),
                              name << ": cannot load preamble");
        const bool have_ihat = load("i_hat_expected", ihat_exp);
        const bool have_out = load("out_expected", out_exp);
        const bool have_i = load("i", i_true);

        double v = 0.0;
        SicConfig cfg;
        cfg.mode = SicMode::Estimate; // golden cases are all Estimate
        if (json_number(meta, "samples_per_symbol", v))
            cfg.samples_per_symbol = static_cast<size_t>(v);
        if (json_number(meta, "sync_reps", v))
            cfg.sync_reps = static_cast<size_t>(v);
        if (json_number(meta, "k_paths", v))
            cfg.k_paths = static_cast<size_t>(v);
        if (json_number(meta, "fit_reps", v))
            cfg.fit_reps = static_cast<size_t>(v);
        if (json_number(meta, "validation_reps", v))
            cfg.validation_reps = static_cast<size_t>(v);
        if (json_number(meta, "roi_lo", v))
            cfg.roi_lo = static_cast<int64_t>(v);
        if (json_number(meta, "roi_hi", v))
            cfg.roi_hi = static_cast<int64_t>(v);
        if (json_number(meta, "coarse_stride", v))
            cfg.coarse_stride = static_cast<size_t>(v);
        if (json_number(meta, "path_min_spacing", v))
            cfg.path_min_spacing = static_cast<size_t>(v);
        if (json_number(meta, "post_tail", v))
            cfg.post_tail = static_cast<size_t>(v);
        if (json_number(meta, "sample_rate", v))
            cfg.sample_rate = v;
        bool b = false;
        if (json_bool(meta, "enable_fractional", b))
            cfg.enable_fractional = b;
        if (json_bool(meta, "enable_cfo", b))
            cfg.enable_cfo = b;
        if (json_number(meta, "detect_threshold", v))
            cfg.detect_threshold = v;
        if (json_number(meta, "validate_threshold", v))
            cfg.validate_threshold = v;

        // ------------------------------------------------------------------
        // S1 legacy replay.
        //
        // The committed S0 vectors were emitted with the S0 contract: support =
        // confirmed repetitions, and the commit gate = disjoint-holdout residual.
        // They carry no support_mode/observed_end keys.  Replaying them through
        // the S1 default (nominal support, no holdout wait) would compare against
        // a different contract, so an S0-era vector is replayed in the EXPLICIT
        // legacy mode it was generated with.  A vector that declares support_mode
        // (future S1 vectors) is honored literally.
        //
        // The one unavoidable difference: S1 estimates the CFO from the TRAINING
        // repetitions only (8 symbols), while S0 averaged over fit+validation.
        // The stored i_hat/out therefore predate the strict 8-symbol CFO
        // estimator; the measured NMSE of that difference is printed and the
        // stored tolerances are widened for the legacy replay only.  Everything
        // else (detection, delays, gains, mask, byte-exactness outside the mask)
        // is still compared with the stored values.
        // ------------------------------------------------------------------
        bool legacy_replay = true;
        std::string smode;
        cfg.support_mode = SicSupportMode::ConfirmedLegacy;
        cfg.online_validation_reps = cfg.validation_reps;
        if (json_string(meta, "support_mode", smode)) {
            if (smode == "configured_nominal" || smode == "ConfiguredNominal") {
                cfg.support_mode = SicSupportMode::ConfiguredNominal;
                cfg.online_validation_reps = 0;
                legacy_replay = false;
            } else if (smode == "confirmed_legacy" || smode == "ConfirmedLegacy") {
                cfg.support_mode = SicSupportMode::ConfirmedLegacy;
                cfg.online_validation_reps = cfg.validation_reps;
                legacy_replay = false;
            }
        }
        if (json_number(meta, "observed_end_sample", v))
            cfg.observed_end_sample = static_cast<int64_t>(v);
        if (json_number(meta, "online_validation_reps", v))
            cfg.online_validation_reps = static_cast<size_t>(v);

        // Tolerances: explicit keys win; fall back to the README guidance.
        double tol_ihat = 1e-3, tol_out = 1e-3, tol_ds = 1.0;
        double estimate_nmse_max = -1.0, suppression_min_db = -1.0;
        double expected_ds = -1.0;
        json_number(meta, "tol_nmse_ihat", tol_ihat);
        json_number(meta, "tol_max_abs_out", tol_out);
        json_number(meta, "detected_start_tol", tol_ds);
        json_number(meta, "estimate_nmse_max", estimate_nmse_max);
        json_number(meta, "suppression_min_db", suppression_min_db);
        json_number(meta, "detected_start", expected_ds);
        const bool require_commit = json_bool(meta, "require_commit", b) ? b : false;
        if (legacy_replay) {
            tol_ihat = std::max(tol_ihat, 5e-3);
            tol_out = std::max(tol_out, 2e-2);
        }

        SicScratch scratch;
        std::string err;
        BOOST_REQUIRE_MESSAGE(sic::prepare_sic_core(cfg,
                                                    sync.data(),
                                                    sync.size(),
                                                    preamble.data(),
                                                    preamble.size(),
                                                    scratch,
                                                    &err),
                              name << ": prepare failed: " << err);
        std::vector<cplx> out(y.size(), cplx(0.f, 0.f));
        std::vector<cplx> ihat(y.size(), cplx(0.f, 0.f));
        SicTemplates t;
        t.sync = sync.data();
        t.sync_len = sync.size();
        t.preamble = preamble.data();
        t.preamble_len = preamble.size();
        SicResult r;
        r.out = out.data();
        r.i_hat = ihat.data();
        const bool ok = sic::sic_one(y.data(), y.size(), t, cfg, scratch, r, nullptr);
        BOOST_REQUIRE_MESSAGE(ok, name << ": sic_one returned false");
        ++cases;
        std::printf("[sic-qa] golden '%s': state=%s detected_start=%lld np=%zu "
                    "cfo=%.1f res=%.4f fit_res=%.4f i_hat=%s out=%s%s\n",
                    name.c_str(),
                    state_name(r.state),
                    static_cast<long long>(r.detected_start),
                    r.n_paths,
                    r.cfo_hz,
                    r.residual_corr,
                    r.fit_residual,
                    have_ihat ? "yes" : "no",
                    have_out ? "yes" : "no",
                    legacy_replay ? "  [S0 legacy replay: strict 8-symbol CFO, "
                                    "tolerances widened]"
                                  : "");
        if (require_commit)
            BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                              static_cast<int>(SicState::Committed));
        if (expected_ds >= 0.0)
            BOOST_CHECK_LE(std::abs(static_cast<double>(r.detected_start) - expected_ds),
                           tol_ds);
        if (have_i && estimate_nmse_max >= 0.0) {
            const double ne = nmse(i_true, ihat, 0, y.size());
            const double sup = 10.0 * std::log10(1.0 / std::max(ne, 1e-30));
            std::printf("[sic-qa] golden '%s': NMSE(i,my_ihat)=%.3e "
                        "suppression=%.2f dB (max %.3e)\n",
                        name.c_str(),
                        ne,
                        sup,
                        estimate_nmse_max);
            BOOST_CHECK_LE(ne, estimate_nmse_max);
            if (suppression_min_db >= 0.0)
                BOOST_CHECK_GE(sup, suppression_min_db);
        }
        if (have_ihat) {
            const double ne = nmse(ihat_exp, ihat, 0, ihat.size());
            std::printf("[sic-qa] golden '%s': i_hat vs reference NMSE=%.3e "
                        "(tol %.3e)\n",
                        name.c_str(),
                        ne,
                        tol_ihat);
            BOOST_CHECK_LE(ne, tol_ihat);
        }
        // F1: direct selection cross-check.  The reference emits its selected
        // ABSOLUTE fractional delays as expected.taus; every reference tau must be
        // matched by a core tau within 0.02 sample (one 0.05 grid step away fails).
        // For a DECLARED stress case (generator and estimator use different
        // band-limited kernels) the core may legitimately select extra paths, so
        // only the reference->core direction is required, with a looser bound.
        {
            std::vector<double> ref_taus;
            bool stress = false;
            json_bool(meta, "declared_stress", stress);
            const double tau_tol = stress ? 0.2 : 0.02;
            if (json_number_array(meta, "taus", ref_taus) && !ref_taus.empty()) {
                double worst = 0.0;
                for (double rt : ref_taus) {
                    double best = 1e30;
                    for (size_t k = 0; k < r.n_taus; ++k)
                        best = std::min(best, std::abs(r.taus[k] - rt));
                    worst = std::max(worst, best);
                    BOOST_CHECK_MESSAGE(best <= tau_tol,
                                        name << ": reference tau " << rt
                                             << " not matched (nearest core tau off by "
                                             << best << ")");
                }
                if (!stress) {
                    BOOST_CHECK_EQUAL(r.n_taus, ref_taus.size());
                    for (size_t k = 0; k < r.n_taus; ++k) {
                        double best = 1e30;
                        for (double rt : ref_taus)
                            best = std::min(best, std::abs(r.taus[k] - rt));
                        worst = std::max(worst, best);
                        BOOST_CHECK_MESSAGE(best <= tau_tol,
                                            name << ": core tau " << r.taus[k]
                                                 << " has no reference match (off by "
                                                 << best << ")");
                    }
                }
                std::printf("[sic-qa] golden '%s': taus vs reference: %zu ref / %zu "
                            "core, max nearest-distance=%.3e (tol %.2g%s)\n",
                            name.c_str(),
                            ref_taus.size(),
                            r.n_taus,
                            worst,
                            tau_tol,
                            stress ? ", declared stress" : "");
            }
        }
        if (have_out) {
            // out = y - i_hat, so comparing out is equivalent to comparing
            // i_hat; use NMSE as the reference README recommends.  Also require
            // out == y byte-for-byte everywhere the reference left it unchanged.
            const double ne = nmse(out_exp, out, 0, out.size());
            size_t outside = 0;
            for (size_t k = 0; k < y.size(); ++k) {
                if (out_exp[k] != y[k]) {
                    continue; // inside the reference mask
                }
                ++outside;
                BOOST_CHECK_EQUAL(out[k], y[k]);
            }
            std::printf("[sic-qa] golden '%s': out vs reference NMSE=%.3e "
                        "(tol %.3e), %zu samples outside mask byte-exact\n",
                        name.c_str(),
                        ne,
                        tol_out,
                        outside);
            BOOST_CHECK_LE(ne, tol_out);
        }
    }
    ::closedir(d);
    std::printf("[sic-qa] golden cross-check: %zu case(s) executed\n", cases);
}

// ---------------------------------------------------------------------------
// Perf E0: column / rotation / reconstruction / self-energy caches.
// Cache-off is the original operation order.  Cache-on must match it
// bit-for-bit, including across model switches and chunk order.
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_sic_core_k_sweep)
{
    const std::vector<Path> paths = {
        { 0.0, cplx(1.f, 0.f) },
        { 7.0, cplx(0.f, 0.5f) },
        { 23.0, cplx(0.276f, -0.116f) },
        { 51.0, cplx(0.2f, 0.f) },
    };
    double supp4 = -1.0;
    for (size_t K : { 1u, 2u, 4u, 8u, 16u }) {
        Fixture f =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, K, 8, 8, 100);
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 10 + K);
        SicResult r = run(f, s);
        const double e = nmse(s.i, f.ihat, 0, f.n);
        const double sup = 10.0 * std::log10(1.0 / std::max(e, 1e-30));
        std::printf("[sic-qa] K=%-2zu state=%s n_paths=%zu cand=%zu cond=%.2f "
                    "suppression=%.2f dB\n",
                    K,
                    state_name(r.state),
                    r.n_paths,
                    r.n_path_candidates,
                    r.gain_condition,
                    sup);
        BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                          static_cast<int>(SicState::Committed));
        if (K == 4)
            supp4 = sup;
        if (K >= 4)
            BOOST_CHECK_GT(sup, 15.0);
    }
    BOOST_CHECK_GT(supp4, 15.0);
}

// ---------------------------------------------------------------------------
// 6. fit_reps sweep (8 / 16 / 32), disjoint validation
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_sic_core_noise_bypass_is_clean)
{
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
    for (uint64_t seed = 0; seed < 8; ++seed) {
        Fixture f =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
        Signal s =
            build_signal(f.model, f.n, f.comm_start, paths, 0.0, 5e-3, 1000 + seed);
        SicResult r = run(f, s);
        if (r.state != SicState::Committed) {
            BOOST_CHECK_EQUAL_COLLECTIONS(
                f.out.begin(), f.out.end(), s.y.begin(), s.y.end());
        } else {
            const double e = nmse(s.i, f.ihat, 0, f.n);
            BOOST_CHECK_LT(e, 0.2);
        }
    }
}

// ---------------------------------------------------------------------------
// 13. gain-fit conditioning indicator (Householder QR; cond2(A) of the design)
// ---------------------------------------------------------------------------
// The final gain fit solves the M x K least squares with a Householder QR and
// reports `gain_condition` = cond2(A), computed from the singular values of the
// K x K factor R (A = QR).  It must stay finite, >= 1 and below the commit
// ceiling; the solve must still produce a committed, model-matching
// reconstruction.  A single fitted column is the sharpest pin: R is 1x1, so
// cond2(A) == 1 exactly.
