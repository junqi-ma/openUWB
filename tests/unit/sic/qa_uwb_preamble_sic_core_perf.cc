// SPDX-License-Identifier: GPL-3.0-or-later
// Shard 3/3 of the preamble-SIC core QA: E0 cache-equivalence and R01
// self-energy geometry performance/consistency cases.
#include "preamble_sic_test_helpers.h"

namespace {

class EnvSet
{
public:
    void set(const char* name, const char* value)
    {
        Item it;
        it.name = name;
        const char* old = std::getenv(name);
        it.had = old != nullptr;
        if (it.had)
            it.old = old;
        if (value)
            ::setenv(name, value, 1);
        else
            ::unsetenv(name);
        items_.push_back(std::move(it));
    }
    ~EnvSet()
    {
        for (auto it = items_.rbegin(); it != items_.rend(); ++it) {
            if (it->had)
                ::setenv(it->name.c_str(), it->old.c_str(), 1);
            else
                ::unsetenv(it->name.c_str());
        }
    }
    EnvSet(const EnvSet&) = delete;
    EnvSet& operator=(const EnvSet&) = delete;
    EnvSet() = default;

private:
    struct Item {
        std::string name;
        bool had = false;
        std::string old;
    };
    std::vector<Item> items_;
};

struct RunSnap {
    SicState state = SicState::Bypassed;
    uint64_t model_id = 0;
    size_t n_paths = 0;
    double cfo_hz = 0.0;
    double taus[8] = {};
    std::vector<cplx> out;
    std::vector<cplx> ihat;
    sic::SicTimings timings;
};

RunSnap run_snap(Fixture& f, const Signal& s)
{
    SicResult r = run(f, s);
    RunSnap o;
    o.state = r.state;
    o.model_id = r.model_id;
    o.n_paths = r.n_paths;
    o.cfo_hz = r.cfo_hz;
    for (size_t k = 0; k < r.n_taus && k < 8; ++k)
        o.taus[k] = r.taus[k];
    o.out = f.out;
    o.ihat = f.ihat;
    o.timings = r.timings;
    return o;
}

size_t n_diff(const std::vector<cplx>& a, const std::vector<cplx>& b)
{
    const size_t n = std::min(a.size(), b.size());
    size_t d = (a.size() == b.size()) ? 0 : (a.size() + b.size());
    for (size_t i = 0; i < n; ++i)
        d += (a[i] != b[i]) ? 1 : 0;
    return d;
}

void expect_same_snap(const RunSnap& a, const RunSnap& b, const char* tag)
{
    BOOST_CHECK_EQUAL(static_cast<int>(a.state), static_cast<int>(b.state));
    BOOST_CHECK_EQUAL(a.model_id, b.model_id);
    BOOST_CHECK_EQUAL(a.n_paths, b.n_paths);
    BOOST_CHECK(std::memcmp(&a.cfo_hz, &b.cfo_hz, sizeof(double)) == 0);
    for (size_t k = 0; k < a.n_paths && k < 8; ++k)
        BOOST_CHECK(std::memcmp(&a.taus[k], &b.taus[k], sizeof(double)) == 0);
    const size_t od = n_diff(a.out, b.out);
    const size_t id = n_diff(a.ihat, b.ihat);
    std::printf("[sic-qa] perf-e0 %s: out_diff=%zu ihat_diff=%zu state=%s "
                "paths=%zu\n",
                tag,
                od,
                id,
                state_name(a.state),
                a.n_paths);
    BOOST_CHECK_EQUAL(od, size_t(0));
    BOOST_CHECK_EQUAL(id, size_t(0));
}

const std::vector<Path> kPerfPaths4 = {
    { 0.0, cplx(1.f, 0.f) },
    { 7.0, cplx(0.f, 0.5f) },
    { 23.0, cplx(0.276f, -0.116f) },
    { 51.0, cplx(0.2f, 0.f) },
};

} // namespace

BOOST_AUTO_TEST_CASE(test_sic_perf_e0_cache_matches_uncached)
{
    const std::vector<Path> paths = { { 0.35, cplx(1.f, 0.f) },
                                      { 14.2, cplx(0.3f, -0.4f) } };
    auto once =
        [&](const char* col, const char* rot, const char* recon, const char* self) {
            EnvSet env;
            env.set("UWB_SIC_COL_CACHE", col);
            env.set("UWB_SIC_ROT_BASE_CACHE", rot);
            env.set("UWB_SIC_RECON_CACHE", recon);
            env.set("UWB_SIC_SELF_ENERGY", self);
            env.set("UWB_SIC_INCR_RESID", "0");
            env.set("UWB_SIC_FUSED_SCORE", "0");
            env.set("UWB_SIC_REFINE_DIAG", "1");
            Fixture f =
                make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
            Signal s =
                build_signal(f.model, f.n, f.comm_start, paths, 25000.0, kQaNoise, 9101);
            f.cfg.observed_end_sample = f.comm_start + 8 * static_cast<int64_t>(kSps);
            return run_snap(f, s);
        };
    const RunSnap off = once("0", "0", "0", "0");
    const RunSnap on = once("1", "1", "1", "1");
    expect_same_snap(off, on, "all-caches");
    BOOST_REQUIRE_EQUAL(static_cast<int>(on.state),
                        static_cast<int>(SicState::Committed));
    const auto& t = on.timings;
    const uint64_t sub = t.refine_init_rot_us + t.refine_init_columns_us +
                         t.refine_init_qr_us + t.refine_resid_us + t.refine_derot_us +
                         t.refine_candidate_us + t.refine_project_us +
                         t.refine_self_energy_us;
    std::printf("[sic-qa] perf-e0 sub_timers=%llu fractional=%llu K=%u M=%u "
                "cand=%u rebuilds=%u hits=%u misses=%u se_hit=%u se_miss=%u "
                "blk_build=%u blk_hit=%u rot_hit=%u rot_miss=%u\n",
                static_cast<unsigned long long>(sub),
                static_cast<unsigned long long>(t.fractional_refine_us),
                t.refine_k,
                t.refine_m,
                t.refine_candidates,
                t.refine_rebuilds,
                t.refine_cache_hits,
                t.refine_cache_misses,
                t.refine_self_energy_hits,
                t.refine_self_energy_misses,
                t.apply_block_builds,
                t.apply_block_hits,
                t.apply_rot_base_hits,
                t.apply_rot_base_misses);
    BOOST_CHECK_LE(sub, t.fractional_refine_us + 16);
    BOOST_CHECK_EQUAL(t.refine_k, static_cast<uint32_t>(on.n_paths));
    BOOST_CHECK_EQUAL(t.refine_candidates, static_cast<uint32_t>(on.n_paths * 21u));
    BOOST_CHECK_EQUAL(t.refine_rebuilds, static_cast<uint32_t>(on.n_paths));
    BOOST_CHECK_EQUAL(t.refine_cache_hits,
                      static_cast<uint32_t>(on.n_paths * on.n_paths));
    BOOST_CHECK_EQUAL(t.refine_cache_misses, uint32_t(0));
    BOOST_CHECK_EQUAL(t.refine_self_energy_hits + t.refine_self_energy_misses,
                      static_cast<uint32_t>(on.n_paths * 21u));
    // The two apply passes share the block that straddles observed_end, or
    // at least the rotation base built while estimating is reused.
    BOOST_CHECK_GT(t.apply_block_builds, uint32_t(0));
    BOOST_CHECK_GT(t.apply_rot_base_hits + t.apply_block_hits, uint32_t(0));

    const RunSnap col_off = once("0", "1", "1", "1");
    const RunSnap rot_off = once("1", "0", "1", "1");
    const RunSnap recon_off = once("1", "1", "0", "1");
    const RunSnap self_off = once("1", "1", "1", "0");
    expect_same_snap(off, col_off, "col-off");
    expect_same_snap(off, rot_off, "rot-off");
    expect_same_snap(off, recon_off, "recon-off");
    expect_same_snap(off, self_off, "self-off");
}

BOOST_AUTO_TEST_CASE(test_sic_perf_e0_k_and_flags)
{
    auto check_k = [](size_t k, const std::vector<Path>& paths, uint64_t seed) {
        Fixture f =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, k, 8, 8, 100);
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, seed);
        RunSnap a = run_snap(f, s);
        BOOST_CHECK_EQUAL(static_cast<int>(a.state),
                          static_cast<int>(SicState::Committed));
        BOOST_CHECK_GE(a.n_paths, size_t(1));
        BOOST_CHECK_LE(a.n_paths, k);
        if (a.n_paths > 0 && a.timings.refine_candidates > 0) {
            BOOST_CHECK_EQUAL(a.timings.refine_rebuilds,
                              static_cast<uint32_t>(a.n_paths));
            BOOST_CHECK_EQUAL(a.timings.refine_cache_hits,
                              static_cast<uint32_t>(a.n_paths * a.n_paths));
        }
        // Same scratch, second window of a different draw must not reuse the
        // first window's output.
        Signal s2 = build_signal(
            f.model, f.n, f.comm_start + 3, paths, 8000.0, kQaNoise, seed + 1);
        RunSnap b = run_snap(f, s2);
        RunSnap b2 = run_snap(f, s2);
        expect_same_snap(b, b2, "repeat-window");
        // Replay the first signal after the second.  A leaked column or
        // reconstruction block would move the output.
        RunSnap a2 = run_snap(f, s);
        expect_same_snap(a, a2, "alternating-model");
    };
    check_k(1, { { 0.0, cplx(1.f, 0.f) } }, 9201);
    check_k(2, { { 0.0, cplx(1.f, 0.f) }, { 18.0, cplx(0.f, 0.45f) } }, 9202);
    check_k(4, kPerfPaths4, 9204);

    // Fractional refine off: no column rebuild, output stable.
    {
        Fixture f =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 2, 8, 8, 100);
        f.cfg.enable_fractional = false;
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) },
                                          { 18.0, cplx(0.2f, 0.2f) } };
        Signal s =
            build_signal(f.model, f.n, f.comm_start, paths, 12000.0, kQaNoise, 9301);
        RunSnap a = run_snap(f, s);
        BOOST_CHECK_EQUAL(a.timings.refine_candidates, uint32_t(0));
        BOOST_CHECK_EQUAL(a.timings.refine_rebuilds, uint32_t(0));
        BOOST_CHECK_EQUAL(a.timings.refine_cache_hits, uint32_t(0));
        RunSnap b = run_snap(f, s);
        expect_same_snap(a, b, "fractional-off");
    }
    // Training span shorter than the prepared fit length: self-energy table
    // must miss, and the window stays deterministic.
    {
        Fixture f =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 2, 8, 8, 100);
        f.cfg.observed_end_sample = f.comm_start + 6 * static_cast<int64_t>(kSps);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 9302);
        RunSnap a = run_snap(f, s);
        BOOST_CHECK_EQUAL(a.timings.refine_self_energy_hits, uint32_t(0));
        RunSnap b = run_snap(f, s);
        expect_same_snap(a, b, "short-prefix");
        if (a.state != SicState::Committed)
            BOOST_CHECK_EQUAL(n_diff(a.out, s.y), size_t(0));
    }
}

BOOST_AUTO_TEST_CASE(test_sic_perf_recon_cache_and_diag)
{
    Fixture f = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 2, 8, 8, 100);
    const std::vector<Path> paths = { { 0.2, cplx(1.f, 0.f) },
                                      { 11.4, cplx(0.25f, 0.1f) } };
    Signal s = build_signal(f.model, f.n, f.comm_start, paths, 15000.0, kQaNoise, 9401);
    const int64_t oe = f.comm_start + 8 * static_cast<int64_t>(kSps);
    f.cfg.observed_end_sample = oe;
    SicModel m = estimate_model_checked(f, s, oe, nullptr);
    BOOST_REQUIRE(m.valid);

    std::vector<cplx> whole, whole_ih;
    apply_chunks(f, m, s, { { 0, s.y.size() } }, whole, whole_ih);

    // Same 128-sample chunk twice: the second call must hit the one block.
    sic::SicTimings tm;
    const size_t nchunk = 128;
    const int64_t chunk_at = f.comm_start;
    std::vector<cplx> c1(nchunk), h1(nchunk), c2(nchunk), h2(nchunk);
    BOOST_REQUIRE(sic::apply_model(m,
                                   s.y.data() + chunk_at,
                                   nchunk,
                                   chunk_at,
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   c1.data(),
                                   h1.data(),
                                   &tm));
    const uint32_t builds = tm.apply_block_builds;
    const uint32_t hits = tm.apply_block_hits;
    BOOST_REQUIRE(sic::apply_model(m,
                                   s.y.data() + chunk_at,
                                   nchunk,
                                   chunk_at,
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   c2.data(),
                                   h2.data(),
                                   &tm));
    BOOST_CHECK_EQUAL(n_diff(c1, c2), size_t(0));
    BOOST_CHECK_EQUAL(n_diff(h1, h2), size_t(0));
    BOOST_CHECK_EQUAL(tm.apply_block_builds, builds);
    BOOST_CHECK_GT(tm.apply_block_hits, hits);

    // A different CFO must miss.  Restoring the CFO rebuilds the original block.
    SicModel mutated = m;
    mutated.cfo_hz += 1000.0;
    std::vector<cplx> cm(nchunk), hm(nchunk);
    const uint32_t builds_before = tm.apply_block_builds;
    BOOST_REQUIRE(sic::apply_model(mutated,
                                   s.y.data() + chunk_at,
                                   nchunk,
                                   chunk_at,
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   cm.data(),
                                   hm.data(),
                                   &tm));
    BOOST_CHECK_GT(tm.apply_block_builds, builds_before);
    BOOST_CHECK(n_diff(c1, cm) > 0 || n_diff(h1, hm) > 0);
    std::vector<cplx> c3(nchunk), h3(nchunk);
    BOOST_REQUIRE(sic::apply_model(m,
                                   s.y.data() + chunk_at,
                                   nchunk,
                                   chunk_at,
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   c3.data(),
                                   h3.data(),
                                   &tm));
    BOOST_CHECK_EQUAL(n_diff(c1, c3), size_t(0));

    // Reverse chunk order and an overlapping pair stay bit-identical.
    std::vector<std::pair<int64_t, size_t>> parts;
    for (int64_t p = 0; p < static_cast<int64_t>(s.y.size());) {
        const size_t len = static_cast<size_t>(
            std::min<int64_t>(3000, static_cast<int64_t>(s.y.size()) - p));
        parts.push_back({ p, len });
        p += static_cast<int64_t>(len);
    }
    std::vector<std::pair<int64_t, size_t>> rev = parts;
    std::reverse(rev.begin(), rev.end());
    std::vector<cplx> out_rev, ih_rev;
    apply_chunks(f, m, s, rev, out_rev, ih_rev);
    BOOST_CHECK_EQUAL(n_diff(whole, out_rev), size_t(0));
    BOOST_CHECK_EQUAL(n_diff(whole_ih, ih_rev), size_t(0));

    std::vector<cplx> overlap(s.y.size()), overlap_ih(s.y.size());
    const size_t mid = std::min<size_t>(5000, s.y.size());
    const size_t back = std::min<size_t>(4000, mid);
    BOOST_REQUIRE(sic::apply_model(m,
                                   s.y.data(),
                                   mid,
                                   0,
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   overlap.data(),
                                   overlap_ih.data()));
    BOOST_REQUIRE(sic::apply_model(m,
                                   s.y.data() + back,
                                   s.y.size() - back,
                                   static_cast<int64_t>(back),
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   overlap.data() + back,
                                   overlap_ih.data() + back));
    BOOST_CHECK_EQUAL(n_diff(whole, overlap), size_t(0));
    BOOST_CHECK_EQUAL(n_diff(whole_ih, overlap_ih), size_t(0));

    // NaN model stays a bypass under the public (copying) entry.
    SicModel bad = m;
    bad.gain[0] = cplx(std::numeric_limits<float>::quiet_NaN(), 0.f);
    std::vector<cplx> bout(nchunk, cplx(9.f, 9.f));
    BOOST_REQUIRE(sic::apply_model(bad,
                                   s.y.data() + chunk_at,
                                   nchunk,
                                   chunk_at,
                                   f.templates(),
                                   f.cfg,
                                   f.scratch,
                                   bout.data(),
                                   nullptr));
    size_t kept = 0;
    for (size_t i = 0; i < nchunk; ++i)
        kept += (bout[i] == s.y[static_cast<size_t>(chunk_at) + i]) ? 1 : 0;
    BOOST_CHECK_EQUAL(kept, nchunk);

    // Diag off: clocks stay zero, counts still fill, IQ matches diag on.
    RunSnap diag_on;
    {
        Fixture fd =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 2, 8, 8, 100);
        Signal sd =
            build_signal(fd.model, fd.n, fd.comm_start, paths, 15000.0, kQaNoise, 9401);
        fd.cfg.observed_end_sample = oe;
        diag_on = run_snap(fd, sd);
    }
    {
        EnvSet env;
        env.set("UWB_SIC_REFINE_DIAG", "0");
        Fixture fd =
            make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 2, 8, 8, 100);
        Signal sd =
            build_signal(fd.model, fd.n, fd.comm_start, paths, 15000.0, kQaNoise, 9401);
        fd.cfg.observed_end_sample = oe;
        RunSnap diag_off = run_snap(fd, sd);
        expect_same_snap(diag_on, diag_off, "diag-off");
        const auto& t = diag_off.timings;
        BOOST_CHECK_EQUAL(t.refine_init_rot_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_init_columns_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_init_qr_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_resid_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_derot_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_candidate_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_project_us, uint64_t(0));
        BOOST_CHECK_EQUAL(t.refine_self_energy_us, uint64_t(0));
        BOOST_CHECK_GT(t.refine_candidates, uint32_t(0));
        BOOST_CHECK_GT(t.fractional_refine_us, uint64_t(0));
    }
}

// PERF-R01.  The self-energy table stores the w0=0 column.  An integer
// relative delay is not enough to reuse it: at a large origin the live
// select_phase float blend can change.  Counts are machine-specific; the
// test requires that a mismatch exists and that a true hit is bitwise.
BOOST_AUTO_TEST_CASE(test_sic_perf_r01_self_energy_geometry)
{
    const auto wrapped = sic::detail::phase_geom(-0.02, 0);
    BOOST_CHECK(wrapped.wrap);
    BOOST_CHECK_EQUAL(wrapped.g, size_t(31));
    const auto plain = sic::detail::phase_geom(-0.05, 0);
    BOOST_CHECK(!plain.wrap);
    BOOST_CHECK(plain.g < size_t(31));
    // Same M/plen.  The wrap bit alone must reject the cached denominator.
    BOOST_CHECK(!sic::detail::self_energy_geometry_matches(
        /*fit_lo=*/0,
        /*rel0=*/0.0,
        /*tau_abs=*/-0.02,
        /*delta=*/0.0,
        /*M=*/64,
        /*plen=*/128));
    // Both sides wrap together: the gate is allowed to hit.
    BOOST_CHECK(sic::detail::self_energy_geometry_matches(
        /*fit_lo=*/0,
        /*rel0=*/-1.0,
        /*tau_abs=*/-1.0,
        /*delta=*/-0.02,
        /*M=*/64,
        /*plen=*/128));
    const auto both = sic::detail::phase_geom(-1.02, 0);
    BOOST_CHECK(both.wrap);

    struct Scan {
        size_t checked = 0;
        size_t geom_hit = 0;
        size_t geom_miss = 0;
        size_t old_rule_wrong = 0;
        size_t hit_bits_bad = 0;
        size_t rel0_miss = 0;
    };
    auto scan_origin = [](SicScratch& sc, int64_t origin) {
        Scan s;
        const int64_t rel_lo = sc.self_energy_rel_min;
        const int64_t rel_hi = rel_lo + static_cast<int64_t>(sc.self_energy_rel_span);
        for (int64_t rel = rel_lo; rel < rel_hi; ++rel) {
            const double tau = static_cast<double>(origin) + static_cast<double>(rel);
            for (int oi = 0; oi <= 20; ++oi) {
                bool geom = false;
                double cached = 0.0;
                double live = 0.0;
                const bool ok = sic::detail::self_energy_denominator_probe(
                    sc, origin, tau, oi, &geom, &cached, &live);
                if (!ok) {
                    s.checked = 0;
                    return s;
                }
                ++s.checked;
                const bool bits = std::memcmp(&cached, &live, sizeof(double)) == 0;
                if (geom) {
                    ++s.geom_hit;
                    if (!bits)
                        ++s.hit_bits_bad;
                } else {
                    ++s.geom_miss;
                    if (rel == 0)
                        ++s.rel0_miss;
                    // The pre-fix gate would have consumed this entry: the
                    // relative delay is an integer inside the table, and the
                    // stored denominator is a positive finite number.
                    if (std::isfinite(cached) && cached > 0.0 && !bits)
                        ++s.old_rule_wrong;
                }
            }
        }
        std::printf("[sic-qa] r01 origin=%lld checked=%zu hit=%zu miss=%zu "
                    "old_wrong=%zu hit_bits_bad=%zu rel0_miss=%zu\n",
                    static_cast<long long>(origin),
                    s.checked,
                    s.geom_hit,
                    s.geom_miss,
                    s.old_rule_wrong,
                    s.hit_bits_bad,
                    s.rel0_miss);
        return s;
    };

    Fixture base = make_fixture(10, 8, static_cast<int64_t>(kPre), 64, 64, 1, 8, 8, 100);
    BOOST_REQUIRE(sic::detail::fill_self_energy_for_test(base.scratch));
    for (int64_t origin : { int64_t(1997), int64_t(122000), int64_t(2000000) }) {
        const Scan s = scan_origin(base.scratch, origin);
        BOOST_REQUIRE_GT(s.checked, size_t(0));
        BOOST_CHECK_EQUAL(s.geom_miss, size_t(0));
        BOOST_CHECK_EQUAL(s.hit_bits_bad, size_t(0));
        BOOST_CHECK_GT(s.geom_hit, size_t(0));
        BOOST_CHECK_EQUAL(s.old_rule_wrong, size_t(0));
    }

    // Random complex template, seed 42.  Not an Ipatov waveform.  The
    // counterexample is "a mismatch can occur", not a frozen diff count.
    {
        std::mt19937 rng(42);
        std::uniform_real_distribution<float> dist(-1.f, 1.f);
        for (cplx& z : base.model.preamble)
            z = cplx(dist(rng), dist(rng));
        base.model.sync.assign(base.model.preamble.begin(),
                               base.model.preamble.begin() +
                                   static_cast<std::ptrdiff_t>(kSps));
        std::string err;
        BOOST_REQUIRE(sic::prepare_sic_core(base.cfg,
                                            base.model.sync.data(),
                                            base.model.sync.size(),
                                            base.model.preamble.data(),
                                            base.model.preamble.size(),
                                            base.scratch,
                                            &err));
        BOOST_REQUIRE(sic::detail::fill_self_energy_for_test(base.scratch));
        for (int64_t origin : { int64_t(1997), int64_t(122000), int64_t(2000000) }) {
            const Scan s = scan_origin(base.scratch, origin);
            BOOST_REQUIRE_GT(s.checked, size_t(0));
            BOOST_CHECK_EQUAL(s.geom_miss, size_t(0));
            BOOST_CHECK_EQUAL(s.hit_bits_bad, size_t(0));
            BOOST_CHECK_GT(s.geom_hit, size_t(0));
        }
        const Scan far = scan_origin(base.scratch, 16777216);
        BOOST_REQUIRE_GT(far.checked, size_t(0));
        BOOST_CHECK_EQUAL(far.hit_bits_bad, size_t(0));
        BOOST_CHECK_GT(far.geom_miss, size_t(0));
        BOOST_CHECK_GT(far.old_rule_wrong, size_t(0));
        BOOST_CHECK_GT(far.rel0_miss, size_t(0));
        BOOST_CHECK_GT(far.geom_hit, size_t(0));
    }

    // Complex-but-still-detectable template: a valid preamble plus a small
    // random imag/real perturbation.  Cache on and off must match, and the
    // ordinary origin must record hits.
    {
        Fixture f = make_fixture(10, 8, static_cast<int64_t>(kPre), 64, 64, 1, 8, 8, 100);
        std::mt19937 rng(42);
        std::uniform_real_distribution<float> dist(-0.05f, 0.05f);
        for (cplx& z : f.model.preamble) {
            z = cplx(z.real() + dist(rng), dist(rng));
        }
        f.model.sync.assign(f.model.preamble.begin(),
                            f.model.preamble.begin() + static_cast<std::ptrdiff_t>(kSps));
        std::string err;
        BOOST_REQUIRE(sic::prepare_sic_core(f.cfg,
                                            f.model.sync.data(),
                                            f.model.sync.size(),
                                            f.model.preamble.data(),
                                            f.model.preamble.size(),
                                            f.scratch,
                                            &err));
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 9501);
        f.cfg.observed_end_sample = f.comm_start + 8 * static_cast<int64_t>(kSps);
        const RunSnap on = run_snap(f, s);
        BOOST_REQUIRE_EQUAL(static_cast<int>(on.state),
                            static_cast<int>(SicState::Committed));
        BOOST_CHECK_GT(on.timings.refine_self_energy_hits, uint32_t(0));
        BOOST_CHECK_EQUAL(on.timings.refine_self_energy_hits +
                              on.timings.refine_self_energy_misses,
                          on.timings.refine_candidates);

        RunSnap off;
        {
            EnvSet env;
            env.set("UWB_SIC_SELF_ENERGY", "0");
            BOOST_REQUIRE(sic::prepare_sic_core(f.cfg,
                                                f.model.sync.data(),
                                                f.model.sync.size(),
                                                f.model.preamble.data(),
                                                f.model.preamble.size(),
                                                f.scratch,
                                                &err));
            off = run_snap(f, s);
        }
        expect_same_snap(on, off, "r01-complex-on-off");
        BOOST_CHECK_EQUAL(off.timings.refine_self_energy_hits, uint32_t(0));
        BOOST_CHECK_EQUAL(off.timings.refine_self_energy_misses,
                          off.timings.refine_candidates);
    }

    // A stale stamp must not use a planted denominator.  Restoring the live
    // VOLK path keeps the committed waveform.
    {
        Fixture f = make_fixture(10, 8, static_cast<int64_t>(kPre), 64, 64, 1, 8, 8, 100);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 9502);
        f.cfg.observed_end_sample = f.comm_start + 8 * static_cast<int64_t>(kSps);
        const RunSnap clean = run_snap(f, s);
        BOOST_REQUIRE_GT(clean.timings.refine_self_energy_hits, uint32_t(0));
        f.scratch.self_energy_template_hash ^= 0x5a5a5a5a5a5a5a5aull;
        for (double& e : f.scratch.self_energy)
            e = 1.0e6;
        const RunSnap poisoned = run_snap(f, s);
        BOOST_CHECK_EQUAL(poisoned.timings.refine_self_energy_hits, uint32_t(0));
        BOOST_CHECK_EQUAL(poisoned.timings.refine_self_energy_misses,
                          poisoned.timings.refine_candidates);
        expect_same_snap(clean, poisoned, "r01-hash-poison");
    }
    {
        Fixture f = make_fixture(10, 8, static_cast<int64_t>(kPre), 64, 64, 1, 8, 8, 100);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 9503);
        f.cfg.observed_end_sample = f.comm_start + 8 * static_cast<int64_t>(kSps);
        const RunSnap clean = run_snap(f, s);
        BOOST_REQUIRE_GT(clean.timings.refine_self_energy_hits, uint32_t(0));
        static const cplx kOtherUp[1] = { cplx(0.f, 0.f) };
        f.scratch.self_energy_up = kOtherUp;
        for (double& e : f.scratch.self_energy)
            e = 1.0e6;
        const RunSnap poisoned = run_snap(f, s);
        BOOST_CHECK_EQUAL(poisoned.timings.refine_self_energy_hits, uint32_t(0));
        expect_same_snap(clean, poisoned, "r01-pointer-poison");
    }

    // Training shorter than the prepared fit length: every candidate misses.
    {
        Fixture f = make_fixture(10, 8, static_cast<int64_t>(kPre), 64, 64, 1, 8, 8, 100);
        f.cfg.observed_end_sample = f.comm_start + 6 * static_cast<int64_t>(kSps);
        const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) } };
        Signal s = build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 9504);
        const RunSnap a = run_snap(f, s);
        BOOST_CHECK_EQUAL(a.timings.refine_self_energy_hits, uint32_t(0));
        if (a.timings.refine_candidates > 0) {
            BOOST_CHECK_EQUAL(a.timings.refine_self_energy_misses,
                              a.timings.refine_candidates);
        }
    }

    // Production search at 2^22.  rel=0 first mixes hits and misses here
    // (4 of 21 grid points).  2^24 stays in the denominator probe above:
    // that probe calls the same self_energy_column_matches the search loop
    // calls, and a 2^24 sic_one would allocate a 16e6-sample window only to
    // re-check that predicate.  The ROI matches the ordinary fixture so the
    // cancellation support is the predicted preamble, not a few dozen samples.
    {
        const int64_t origin = 4194304; // 2^22
        const int64_t train = 8 * static_cast<int64_t>(kSps);
        auto make_far = [&](Fixture& f) {
            f.model = build_model(10, 8);
            f.comm_start = origin;
            f.n = static_cast<size_t>(origin + train + 4000);
            f.cfg = SicConfig{};
            f.cfg.mode = SicMode::Estimate;
            f.cfg.comm_code = 10;
            f.cfg.sync_reps = 8;
            f.cfg.samples_per_symbol = kSps;
            f.cfg.sample_rate = kFs;
            f.cfg.k_paths = 1;
            f.cfg.fit_reps = 8;
            f.cfg.validation_reps = 8;
            f.cfg.roi_lo = origin - 64;
            f.cfg.roi_hi = origin + train + 64;
            f.cfg.coarse_stride = 14;
            f.cfg.path_min_spacing = 0;
            f.cfg.post_tail = 100;
            f.cfg.enable_fractional = true;
            f.cfg.enable_cfo = true;
            f.cfg.detect_threshold = 0.5;
            f.cfg.validate_threshold = 0.25;
            f.cfg.observed_end_sample = origin + train;
            f.out.assign(f.n, cplx(0.f, 0.f));
            f.ihat.assign(f.n, cplx(0.f, 0.f));
            std::string err;
            BOOST_REQUIRE(sic::prepare_sic_core(f.cfg,
                                                f.model.sync.data(),
                                                f.model.sync.size(),
                                                f.model.preamble.data(),
                                                f.model.preamble.size(),
                                                f.scratch,
                                                &err));
        };
        Fixture on_fx;
        make_far(on_fx);
        Signal s;
        s.y.assign(on_fx.n, cplx(0.f, 0.f));
        add_path(s.y, on_fx.model.preamble, origin, 0.0, cplx(1.f, 0.f), 0.0);
        double emax = 0.0;
        double dot_r = 0.0, dot_i = 0.0, e_sig = 0.0, e_tpl = 0.0;
        const std::vector<cplx>& syn = on_fx.model.sync;
        for (size_t k = 0; k < syn.size(); ++k) {
            const cplx yk = s.y[static_cast<size_t>(origin) + k];
            const cplx sk = syn[k];
            if (k < 64)
                emax = std::max(emax, static_cast<double>(std::norm(yk)));
            dot_r += static_cast<double>(yk.real()) * sk.real() +
                     static_cast<double>(yk.imag()) * sk.imag();
            dot_i += static_cast<double>(yk.imag()) * sk.real() -
                     static_cast<double>(yk.real()) * sk.imag();
            e_sig += static_cast<double>(std::norm(yk));
            e_tpl += static_cast<double>(std::norm(sk));
        }
        const double manual_nrm =
            (e_sig > 0.0 && e_tpl > 0.0)
                ? std::sqrt(dot_r * dot_r + dot_i * dot_i) / std::sqrt(e_sig * e_tpl)
                : 0.0;
        const SicResult raw = run(on_fx, s);
        std::printf("[sic-qa] r01 far origin=%lld state=%s reason=%s reject=%s "
                    "coarse=%llu fine=%llu detected=%d start=%lld emax=%.6g "
                    "manual_nrm=%.6g n=%zu hits=%u misses=%u cand=%u\n",
                    static_cast<long long>(origin),
                    state_name(raw.state),
                    raw.bypass_reason ? raw.bypass_reason : "",
                    raw.reject_reason ? raw.reject_reason : "",
                    static_cast<unsigned long long>(raw.coarse_candidates),
                    static_cast<unsigned long long>(raw.fine_candidates),
                    raw.detected ? 1 : 0,
                    static_cast<long long>(raw.detected_start),
                    emax,
                    manual_nrm,
                    s.y.size(),
                    raw.timings.refine_self_energy_hits,
                    raw.timings.refine_self_energy_misses,
                    raw.timings.refine_candidates);
        std::fflush(stdout);
        RunSnap on;
        on.state = raw.state;
        on.model_id = raw.model_id;
        on.n_paths = raw.n_paths;
        on.cfo_hz = raw.cfo_hz;
        for (size_t k = 0; k < raw.n_taus && k < 8; ++k)
            on.taus[k] = raw.taus[k];
        on.out = on_fx.out;
        on.ihat = on_fx.ihat;
        on.timings = raw.timings;
        BOOST_REQUIRE_EQUAL(static_cast<int>(on.state),
                            static_cast<int>(SicState::Committed));
        BOOST_REQUIRE_GT(on.timings.refine_candidates, uint32_t(0));
        BOOST_CHECK_GT(on.timings.refine_self_energy_hits, uint32_t(0));
        BOOST_CHECK_GT(on.timings.refine_self_energy_misses, uint32_t(0));
        BOOST_CHECK_EQUAL(on.timings.refine_self_energy_hits +
                              on.timings.refine_self_energy_misses,
                          on.timings.refine_candidates);

        Fixture off_fx;
        RunSnap off;
        {
            EnvSet env;
            env.set("UWB_SIC_SELF_ENERGY", "0");
            make_far(off_fx);
            off = run_snap(off_fx, s);
        }
        expect_same_snap(on, off, "r01-far-on-off");
        BOOST_REQUIRE_EQUAL(on.out.size(), off.out.size());
        BOOST_REQUIRE_EQUAL(on.ihat.size(), off.ihat.size());
        BOOST_CHECK(std::memcmp(on.out.data(),
                                off.out.data(),
                                on.out.size() * sizeof(cplx)) == 0);
        BOOST_CHECK(std::memcmp(on.ihat.data(),
                                off.ihat.data(),
                                on.ihat.size() * sizeof(cplx)) == 0);
        BOOST_CHECK(std::memcmp(on.out.data(), s.y.data(), sizeof(cplx)) == 0);
        BOOST_CHECK(std::memcmp(&on.out.back(), &s.y.back(), sizeof(cplx)) == 0);
        BOOST_CHECK(std::memcmp(&off.out.back(), &s.y.back(), sizeof(cplx)) == 0);
    }
}

BOOST_AUTO_TEST_CASE(test_sic_core_fit_reps_sweep)
{
    const std::vector<Path> paths = {
        { 0.0, cplx(1.f, 0.f) },
        { 7.0, cplx(0.f, 0.5f) },
    };
    for (size_t fr : { 8u, 16u, 32u }) {
        Fixture f =
            make_fixture(10, 48, static_cast<int64_t>(kPre), 64, 64, 4, fr, 8, 100);
        Signal s =
            build_signal(f.model, f.n, f.comm_start, paths, 0.0, kQaNoise, 100 + fr);
        SicResult r = run(f, s);
        const double e = nmse(s.i, f.ihat, 0, f.n);
        std::printf("[sic-qa] fit_reps=%-2zu state=%s NMSE=%.3e suppression=%.2f "
                    "dB matrix_us=%llu\n",
                    fr,
                    state_name(r.state),
                    e,
                    10.0 * std::log10(1.0 / std::max(e, 1e-30)),
                    static_cast<unsigned long long>(r.timings.matrix_build_us));
        BOOST_CHECK_EQUAL(static_cast<int>(r.state),
                          static_cast<int>(SicState::Committed));
        BOOST_CHECK_LT(e, 0.02);
    }
}

// ---------------------------------------------------------------------------
// 7. narrow vs wide ROI
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_sic_core_energy_prefix_fallback)
{
    const std::vector<Path> paths = { { 0.0, cplx(1.f, 0.f) }, { 7.0, cplx(0.f, 0.5f) } };

    // (a) reference: prepared with the WIDE roi, run with the same wide roi.
    Fixture fw = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    Signal s = build_signal(fw.model, fw.n, fw.comm_start, paths, 0.0, kQaNoise, 55);
    SicResult rw = run(fw, s);

    // (b) prepared with a NARROW roi (tiny prefix), run with the WIDE runtime roi.
    Fixture fn = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    SicConfig narrow = fn.cfg;
    narrow.roi_lo = fn.comm_start - 4;
    narrow.roi_hi = fn.comm_start + 4;
    std::string err;
    BOOST_REQUIRE_MESSAGE(sic::prepare_sic_core(narrow,
                                                fn.model.sync.data(),
                                                fn.model.sync.size(),
                                                fn.model.preamble.data(),
                                                fn.model.preamble.size(),
                                                fn.scratch,
                                                &err),
                          "narrow prepare failed: " << err);
    SicResult rn = run(fn, s); // fn.cfg still holds the WIDE roi

    // (c) prepared with roi 0,0 (EMPTY prefix), run with the full-window roi.
    Fixture fe = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    SicConfig empty = fe.cfg;
    empty.roi_lo = 0;
    empty.roi_hi = 0;
    BOOST_REQUIRE_MESSAGE(sic::prepare_sic_core(empty,
                                                fe.model.sync.data(),
                                                fe.model.sync.size(),
                                                fe.model.preamble.data(),
                                                fe.model.preamble.size(),
                                                fe.scratch,
                                                &err),
                          "empty prepare failed: " << err);
    fe.cfg.roi_lo = 0;
    fe.cfg.roi_hi = 0;
    SicResult re = run(fe, s);

    // (d) prefix sized for the wide roi but a FULL-WINDOW runtime roi.
    Fixture fm = make_fixture(10, 32, static_cast<int64_t>(kPre), 64, 64, 4, 8, 8, 100);
    fm.cfg.roi_lo = 0;
    fm.cfg.roi_hi = 0;
    SicResult rm = run(fm, s);

    std::printf("[sic-qa] prefix-fallback wide: ds=%lld coarse=%llu fine=%llu np=%zu "
                "| narrow-prep wide-run: ds=%lld coarse=%llu fine=%llu np=%zu | "
                "empty-prep full: ds=%lld coarse=%llu np=%zu | wide-prep full: "
                "ds=%lld coarse=%llu np=%zu\n",
                static_cast<long long>(rw.detected_start),
                static_cast<unsigned long long>(rw.coarse_candidates),
                static_cast<unsigned long long>(rw.fine_candidates),
                rw.n_paths,
                static_cast<long long>(rn.detected_start),
                static_cast<unsigned long long>(rn.coarse_candidates),
                static_cast<unsigned long long>(rn.fine_candidates),
                rn.n_paths,
                static_cast<long long>(re.detected_start),
                static_cast<unsigned long long>(re.coarse_candidates),
                re.n_paths,
                static_cast<long long>(rm.detected_start),
                static_cast<unsigned long long>(rm.coarse_candidates),
                rm.n_paths);

    // (b) must reproduce (a) exactly: same detection, same candidate bookkeeping.
    // First prove the fallback is really exercised: the narrow-prepared prefix
    // must not reach the wide runtime ROI's search extent, and the empty prepare
    // must have produced no prefix at all.
    BOOST_REQUIRE(!fn.scratch.energy_prefix.empty());
    BOOST_REQUIRE(!fw.scratch.energy_prefix.empty());
    BOOST_CHECK_LT(fn.scratch.energy_prefix.size(), fw.scratch.energy_prefix.size());
    BOOST_CHECK_LT(fn.scratch.energy_prefix.size() - 1,
                   static_cast<size_t>(fn.cfg.roi_hi - 1) + fn.scratch.sync_len);
    BOOST_CHECK(fe.scratch.energy_prefix.empty());
    BOOST_CHECK_EQUAL(rn.detected_start, rw.detected_start);
    BOOST_CHECK_EQUAL(rn.coarse_candidates, rw.coarse_candidates);
    BOOST_CHECK_EQUAL(rn.fine_candidates, rw.fine_candidates);
    BOOST_CHECK_EQUAL(rn.correlations, rw.correlations);
    BOOST_CHECK_EQUAL(rn.n_paths, rw.n_paths);
    BOOST_CHECK_EQUAL(static_cast<int>(rn.state), static_cast<int>(rw.state));
    // (b) and (c): the two runs that actually consult the fallback (all lags in
    // (c/the tail of (b)) must agree with each other on the full window.
    BOOST_CHECK_EQUAL(re.detected_start, rm.detected_start);
    BOOST_CHECK_EQUAL(re.coarse_candidates, rm.coarse_candidates);
    BOOST_CHECK_EQUAL(re.n_paths, rm.n_paths);
    BOOST_CHECK_EQUAL(static_cast<int>(re.state), static_cast<int>(SicState::Committed));
}


// ---------------------------------------------------------------------------
// S1 QA: strict 8-SYNC estimation, frozen model, full-nominal cancellation.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// S1-1. prefix isolation: mutating/zeroing everything at or after observed_end
// must not change ANY frozen-model field (and ASan/guard-page-wise, no read may
// touch it).
// ---------------------------------------------------------------------------
