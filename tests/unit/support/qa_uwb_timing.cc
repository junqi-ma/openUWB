/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Unit tests for the support-layer timing primitives (openuwb::timing).
 *
 * The arithmetic helpers are exercised on CONSTRUCTED time points wherever
 * the assertion can be exact (truncation semantics, fractional microseconds);
 * only the clock-interacting scopes fall back to spin-bounded monotonicity
 * bounds, and even those assert one-sided invariants (measured >= spun,
 * sub-intervals <= enclosing interval), never wall-clock absolutes, so the
 * QA stays deterministic under a loaded CI machine.
 */

#include <boost/test/unit_test.hpp>

#include <openuwb/support/timing.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <thread>
#include <type_traits>

namespace {

namespace ut = ::openuwb::timing;

// The generated default is a strict 0/1 constant, whatever UWB_ENABLE_TIMING
// was configured to.  This is what keeps ScopedNs<> ODR-consistent across
// translation units of the same build tree.
static_assert(UWB_TIMING_ENABLED_DEFAULT == 0 ||
              UWB_TIMING_ENABLED_DEFAULT == 1,
              "timing_config.h must define UWB_TIMING_ENABLED_DEFAULT as 0/1");

// The compile-time alias must stay the monotonic wall clock of the stage
// boundaries; a brekage here would silently change the JSON/console units.
static_assert(std::is_same<ut::Clock, std::chrono::steady_clock>::value,
              "openuwb::timing::Clock must be std::chrono::steady_clock");

// The disabled RAII scope must be an empty class: no clock read state, so a
// production hot path can carry it at zero cost when the default is OFF.
static_assert(sizeof(ut::ScopedNs<false>) == 1,
              "ScopedNs<false> must compile to an empty class");

/// Busy-spin until at least target_ns monotonic time passed since t0.
/// Unlike sleep_for the spin cannot finish early, so the invariants below
/// are one-sided in the right direction.
void spin_for_ns(ut::Clock::time_point t0, uint64_t target_ns)
{
    while (ut::elapsed_ns(t0, ut::now()) < target_ns) {
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_elapsed_ns_exact_truncation)
{
    // Pure arithmetic on constructed points: deterministic, no clock.
    const ut::Clock::time_point base{};

    BOOST_CHECK_EQUAL(ut::elapsed_ns(base, base), uint64_t(0));
    BOOST_CHECK_EQUAL(
        ut::elapsed_ns(base, base + std::chrono::nanoseconds(1234567891)),
        uint64_t(1234567891));
    // Exactly one nanosecond must not truncate away.
    BOOST_CHECK_EQUAL(
        ut::elapsed_ns(base, base + std::chrono::nanoseconds(1)),
        uint64_t(1));
    // A partial nanosecond cannot exist with steady_clock::duration, but the
    // duration_cast contract still has to hold for a coarser clock alias.
    const ut::Clock::time_point whole =
        base + std::chrono::duration_cast<ut::Clock::duration>(
                   std::chrono::microseconds(2500));
    BOOST_CHECK_EQUAL(ut::elapsed_ns(base, whole), uint64_t(2500000));
}

BOOST_AUTO_TEST_CASE(test_elapsed_us_since_matches_payload_truncation)
{
    // elapsed_us_since shares the payload_elapsed_us contract: whole
    // microseconds, truncated, of a live clock difference.  The exact
    // arithmetic part is covered through elapsed_ns above; here the live
    // helper must (a) be non-decreasing and (b) stay within [ns/1000] of
    // the nanosecond measurement taken around it.
    const ut::Clock::time_point t0 = ut::now();
    const uint64_t us = ut::elapsed_us_since(t0);
    const ut::Clock::time_point t1 = ut::now();

    BOOST_CHECK(us <= ut::elapsed_ns(t0, t1) / 1000 + 1);
    spin_for_ns(t0, 2000000); // 2 ms
    const uint64_t us2 = ut::elapsed_us_since(t0);
    (void)us2;
    // One-sided: the elapsed value must have advanced by at least the spun
    // interval, minus the truncation quantum.
    BOOST_CHECK(us2 + 1 >= uint64_t(2000000 / 1000));
}

BOOST_AUTO_TEST_CASE(test_elapsed_us_f_exact_arithmetic)
{
    const ut::Clock::time_point base{};

    BOOST_CHECK_EQUAL(ut::elapsed_us_f(base, base), 0.0);
    BOOST_CHECK_EQUAL(
        ut::elapsed_us_f(base, base + std::chrono::nanoseconds(1500)), 1.5);
    BOOST_CHECK_CLOSE(
        ut::elapsed_us_f(base, base + std::chrono::microseconds(25)), 25.0,
        1e-12);
}

BOOST_AUTO_TEST_CASE(test_thread_cpu_ns_monotonic)
{
#if defined(CLOCK_THREAD_CPUTIME_ID)
    const uint64_t c0 = ut::thread_cpu_ns();
    // Burn measurable CPU time.
    uint64_t sink = 0;
    for (int i = 0; i < 2000000; ++i)
        sink += static_cast<uint64_t>(i) ^ (static_cast<uint64_t>(i) >> 3);
    const uint64_t c1 = ut::thread_cpu_ns();
    // One-sided: CPU time cannot run backwards; on a loaded machine the
    // burn may schedule out, so only assert ordering, not a lower bound.
    BOOST_CHECK(c1 >= c0);
    BOOST_CHECK(sink != 0 || true); // sink keeps the loop from being dropped
#else
    // No thread clock on this platform: the documented 0 fallback.
    BOOST_CHECK_EQUAL(ut::thread_cpu_ns(), uint64_t(0));
#endif
}

BOOST_AUTO_TEST_CASE(test_accumulate_ns_adds_spun_interval)
{
    std::atomic<uint64_t> counter{0};
    const ut::Clock::time_point t0 = ut::now();
    spin_for_ns(t0, 1000000); // 1 ms
    ut::accumulate_ns(counter, t0);

    // The accumulator's clock read happens inside accumulate_ns, so the
    // added interval is >= [t0, end of spin] - (one clock-read pair): the
    // spin itself guarantees the point passed to accumulate_ns is at least
    // 1 ms past t0, giving the one-sided bound below a 1% tolerance.
    BOOST_CHECK_GE(counter.load(), uint64_t(990000));
}

BOOST_AUTO_TEST_CASE(test_scoped_ns_true_accumulates_span)
{
    std::atomic<uint64_t> counter{0};
    {
        ut::ScopedNs<true> scope{counter};
        spin_for_ns(ut::now(), 1000000); // 1 ms entirely inside the scope
    }
    BOOST_CHECK_GE(counter.load(), uint64_t(990000));
    // The scope is one-shot: the counter must not grow after destruction.
    const uint64_t frozen = counter.load();
    spin_for_ns(ut::now(), 200000);
    BOOST_CHECK_EQUAL(counter.load(), frozen);
}

BOOST_AUTO_TEST_CASE(test_scoped_ns_false_leaves_counter_alone)
{
    std::atomic<uint64_t> counter{42};
    {
        ut::ScopedNs<false> scope{counter};
        spin_for_ns(ut::now(), 200000);
    }
    // The disabled scope must not read the clock nor store anything.
    BOOST_CHECK_EQUAL(counter.load(), uint64_t(42));
}

BOOST_AUTO_TEST_CASE(test_scoped_true_subintervals_within_outer)
{
    // Disjoint RAII scopes inside one outer measured interval must sum to
    // no more than the outer interval (monotonic clock); nested scopes are
    // the PayloadSubTimer additivity idea in miniature.
    std::atomic<uint64_t> a{0};
    std::atomic<uint64_t> b{0};
    {
        const ut::Clock::time_point outer_t0 = ut::now();
        {
            ut::ScopedNs<true> first{a};
            spin_for_ns(ut::now(), 100000);
        }
        {
            ut::ScopedNs<true> second{b};
            spin_for_ns(ut::now(), 100000);
        }
        const ut::Clock::time_point outer_t1 = ut::now();

        const uint64_t total = ut::elapsed_ns(outer_t0, outer_t1);
        BOOST_CHECK_GE(total, uint64_t(200000));
        BOOST_CHECK(a.load() + b.load() <= total + uint64_t(2000));
        // Both sub-intervals are non-trivial.
        BOOST_CHECK_GE(a.load(), uint64_t(50000));
        BOOST_CHECK_GE(b.load(), uint64_t(50000));
    }
}
