/* -*- c++ -*- */
/*
 * Copyright 2026
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * openuwb::timing -- the single home of the wall-clock measurement
 * primitives used by every domain.
 *
 * Layers (see lib/CMakeLists.txt for the domain DAG):
 *   support  THESE primitives: clock aliases, elapsed helpers, thread CPU
 *            time, and the optional RAII accumulation scopes.
 *   runtime  the stage SEMANTICS: the monotonic boundary fields
 *            (resample_end / sic_end / demod_end / compute_end,
 *            rx_push_time / worker_start_time / commit_time) stay where the
 *            results are produced and keep their names, because the JSON
 *            summaries and console stats quote them verbatim.
 *   apps/tools  consumption: ABBA pairing, percentiles, reports.
 *
 * Contract shared with the historical code this header deduplicates:
 *   - Clock is std::chrono::steady_clock everywhere (monotonic, never the
 *     RF plane; see twr_config_util.h HostMonotonic).
 *   - Truncated integer helpers keep the exact duration_cast truncation of
 *     the call sites they replace, so reported microsecond values do not
 *     move (payload_elapsed_us, elapsed_ns).
 *   - accumulate_ns() is a relaxed fetch_add on the caller's atomic, the
 *     lock-free aggregation pattern of scheduled_pipeline.cc and
 *     uwb_fixed_tone_cancel_source.cc: no locks, no allocation, safe from
 *     several DSP workers concurrently.
 *   - ScopedNs<Enabled> is the srsRAN tprof idea in this tree's idiom: the
 *     Enabled=false specialisation is an empty class that never reads the
 *     clock, so production hot paths can carry the scope guard at zero cost
 *     and decide by compile-time default or per-TU <true>.
 *
 * This header is deliberately baseline-ISA (pure <chrono>/<ctime>, no SIMD,
 * no Volk): openuwb_support refuses uwb_simd_apply() and timing must not
 * change that.
 */
#ifndef INCLUDED_OPENUWB_SUPPORT_TIMING_H
#define INCLUDED_OPENUWB_SUPPORT_TIMING_H

#include <openuwb/support/timing_config.h>

#include <atomic>
#include <chrono>
#include <cstdint>

#if !defined(_WIN32)
#include <ctime>
#endif

namespace openuwb {
namespace timing {

/// The wall clock. Monotonic; used for every stage boundary in the tree.
using Clock = std::chrono::steady_clock;

/// Compile-time default of the optional RAII scopes (timing_config.h).
constexpr bool kEnabledByDefault = (UWB_TIMING_ENABLED_DEFAULT == 1);

inline Clock::time_point now() { return Clock::now(); }

/// Nanoseconds between two monotonic points, truncated towards zero.
/// Matches the duration_cast<nanoseconds> of the historical call sites.
inline uint64_t elapsed_ns(Clock::time_point a, Clock::time_point b)
{
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count());
}

/// Whole microseconds elapsed since t0, truncated like the payload sub-timer
/// laps (the payload_elapsed_us contract of demod_fec.h).
inline uint64_t elapsed_us_since(Clock::time_point t0)
{
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(now() - t0)
            .count());
}

/// Fractional microseconds between two points (console/report statistics,
/// the live_rx_session elapsed_us contract).
inline double elapsed_us_f(Clock::time_point a, Clock::time_point b)
{
    return std::chrono::duration<double, std::micro>(b - a).count();
}

/// Add the nanoseconds elapsed since t0 to an atomic counter.  Relaxed is
/// enough: the counters are telemetry, never ordering.  The fetch_add
/// replaces the hand-rolled wait_ns_/copy/fold/probe accumulations.
inline void accumulate_ns(std::atomic<uint64_t>& counter,
                           Clock::time_point t0)
{
    counter.fetch_add(elapsed_ns(t0, now()), std::memory_order_relaxed);
}

#if defined(CLOCK_THREAD_CPUTIME_ID)
/// Thread CPU time in nanoseconds (CLOCK_THREAD_CPUTIME_ID).  Distinguish
/// "actually computing" from "wall clock including waits/scheduling".
inline uint64_t thread_cpu_ns()
{
    ::timespec t{};
    ::clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return static_cast<uint64_t>(t.tv_sec) * 1000000000ull +
           static_cast<uint64_t>(t.tv_nsec);
}
#else
/// No thread-CPU clock on this platform (non-POSIX): report 0 so the
/// console paths keep compiling and the caller sees "not measured".
inline uint64_t thread_cpu_ns() { return 0; }
#endif

/// Optional RAII scope that adds its wall-clock duration to an atomic
/// counter on destruction.  Enabled defaults to UWB_TIMING_ENABLED_DEFAULT;
/// pass ScopedNs<true> explicitly for a diagnostic scope that must survive
/// an -DUWB_ENABLE_TIMING=OFF build.
template <bool Enabled = kEnabledByDefault>
class ScopedNs;

/// Disabled: an empty class.  No clock read, no store, nothing.
template <>
class ScopedNs<false> {
public:
    explicit ScopedNs(std::atomic<uint64_t>&) {}
    ScopedNs(const ScopedNs&) = delete;
    ScopedNs& operator=(const ScopedNs&) = delete;
};

/// Enabled: one clock read at construction, one fetch_add at destruction.
template <>
class ScopedNs<true> {
public:
    explicit ScopedNs(std::atomic<uint64_t>& counter)
        : counter_(counter), t0_(now()) {}
    ~ScopedNs() { accumulate_ns(counter_, t0_); }
    ScopedNs(const ScopedNs&) = delete;
    ScopedNs& operator=(const ScopedNs&) = delete;

private:
    std::atomic<uint64_t>& counter_;
    Clock::time_point t0_;
};

} // namespace timing
} // namespace openuwb

#endif /* INCLUDED_OPENUWB_SUPPORT_TIMING_H */
