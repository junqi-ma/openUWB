// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <openuwb/phy/hrp/rx/demod_sync.h>
#include <functional>

namespace openuwb::demod::core {
// ---------------------------------------------------------------------------
// Stage 3 — SFD timing refinement (R1).
// Selects the best SFD template (from the profile; golden uses IEEE legacy
// SFD #0 = {0,1,0,-1,1,0,0,-1}) and refines the preamble start at full rate.
// For Phase-1 the SFD search is done on the symbol grid: the SFD starts one
// symbol after the last SYNC (i.e. after preamble_repetitions symbols).
// We locate the SFD by correlating the SFD template over a narrow window
// around the expected position, then refine the preamble start.
// ---------------------------------------------------------------------------

// Range-aware variant: before the discrete integer-symbol grid reads a sample
// interval it calls ensure_derotated(lo, hi) so the caller can incrementally
// rotate only the uncovered gaps instead of re-rotating a wide window.  A
// default-constructed callback is a no-op (the caller already rotated every
// interval the search reads), so the narrow ±64 / ±4-symbol windows behave
// exactly as before.
bool stage_sfd(const std::complex<float>* rx,
               size_t n,
               const Qm35825Profile& profile,
               const TimingResult& timing,
               const std::vector<int8_t>& sfd_sequence,
               const std::vector<std::complex<float>>& template_wf,
               SfdResult& out,
               DemodScratch& scratch,
               const std::function<void(size_t, size_t)>& ensure_derotated);

// Backward-compatible 9-parameter form.  The legacy `prepare_grid` callback is
// invoked once, lazily, right before the first discrete-grid pass -- exactly
// where the pre-refactor implementation called it.
bool stage_sfd(const std::complex<float>* rx,
               size_t n,
               const Qm35825Profile& profile,
               const TimingResult& timing,
               const std::vector<int8_t>& sfd_sequence,
               const std::vector<std::complex<float>>& template_wf,
               SfdResult& out,
               DemodScratch& scratch,
               const std::function<void()>& prepare_grid = {});

} // namespace openuwb::demod::core
