// SPDX-License-Identifier: GPL-3.0-or-later
// Run-time selection of the BPM-BPSK algorithm ablation levels (M0..M4).
//
// The switches live in Qm35825Profile so the same binary can produce every
// ablation point without recompiling:
//
//   M0  previous fused baseline (4 candidates, split LLR, truncated)
//   M1  + known-hop 2-candidate correlation
//   M2  + joint bi-orthogonal ML Viterbi (+ terminated tail)
//   M3  + PHR pilot phase refine + symbol DD-PLL
//   M4  + fused fractional-phase SFO tracking
//   (opt-in via UWB_BPM_CIR_DENOISE=1: adaptive CIR-tap noise-floor shrinkage)
//
// Set UWB_BPM_LEVEL=0..4 in the environment of the app.  When unset, the app
// runs the recommended low-SNR payload profile (M4 with the PHR pilot refine
// off); UWB_BPM_LEVEL=0 still restores the previous M0 behaviour for A/B.
#pragma once

#include <gnuradio/uwb/uwb_phy_profile.h>

#include <cstdlib>

namespace openuwb::apps {

// Production default for the live receiver: M4 (known-hop + joint bi-orthogonal
// ML Viterbi + terminated tail + fused SFO tracking).  The PHR pilot refine is
// kept off because it costs ~40 us/packet with no measured FER gain; enable it
// with UWB_BPM_LEVEL=3/4 plus UWB_BPM_PHR_REFINE=1 when running the M3 ablation.
inline constexpr int kDefaultBpmAblationLevel = 4;

inline bool bpm_ablation_level_explicit()
{
    const char* v = std::getenv("UWB_BPM_LEVEL");
    return v && *v;
}

inline int bpm_ablation_level()
{
    if (!bpm_ablation_level_explicit())
        return kDefaultBpmAblationLevel;
    return std::atoi(std::getenv("UWB_BPM_LEVEL"));
}

inline bool env_flag(const char* name, bool dflt)
{
    const char* v = std::getenv(name);
    if (!v || !*v)
        return dflt;
    return std::atoi(v) != 0;
}

// `fused_baseline` forces the previous-round fused path on, so an explicitly
// selected M0 is the 2d8f49b fused behaviour rather than the app's own default.
// The default production profile must respect the app's --fused/--ondemand
// flags, so the forcing only applies when UWB_BPM_LEVEL was set explicitly.
inline void apply_bpm_ablation(demod::Qm35825Profile& p, bool fused_baseline)
{
    const bool level_explicit = bpm_ablation_level_explicit();
    const int level = bpm_ablation_level();
    if (fused_baseline && level_explicit) {
        p.rx_softchip_on_demand = true;
        p.rx_fused_payload_cmf = true;
        p.rx_fused_phr_cmf = true;
    }
    p.rx_bpm_known_hop = level >= 1;
    p.rx_bpm_joint_ml = level >= 2;
    p.rx_bpm_terminated_viterbi = level >= 2;
    // DD-PLL is implemented and QA-verified, but it costs ~40 us/packet in the
    // 1208-symbol kernel and gave no FER gain on the 6x40 fixtures or the
    // -13 dB sweep, so it stays opt-in (UWB_BPM_DDPLL=1).
    p.rx_bpm_ddpll = false;
    // The PHR pilot refine is part of the M3 ablation, but has no measured FER
    // gain.  Keep it off on the default (unset) path; an explicit
    // UWB_BPM_LEVEL=3/4 or UWB_BPM_PHR_REFINE=1 can still turn it on.
    p.rx_phr_phase_refine = level_explicit && level >= 3;
    p.rx_fused_sfo_track = level >= 4;
    // Adaptive CIR-tap noise-floor shrinkage for the fused PHR/payload kernels.
    // TODO(openuwb-cir-denoise): A/B validate UWB_BPM_CIR_DENOISE=1 vs 0 on
    // real-capture and dense-multipath recordings (confirm no regression when
    // multipath energy spans > 50% of the 38-tap window) before enabling by
    // default.
    p.rx_cir_tap_denoise = env_flag("UWB_BPM_CIR_DENOISE", false);
    // Per-switch overrides for switch isolation during the ablation.
    p.rx_bpm_known_hop = env_flag("UWB_BPM_KNOWN_HOP", p.rx_bpm_known_hop);
    p.rx_bpm_joint_ml = env_flag("UWB_BPM_JOINT_ML", p.rx_bpm_joint_ml);
    p.rx_bpm_terminated_viterbi =
        env_flag("UWB_BPM_TERMINATED", p.rx_bpm_terminated_viterbi);
    p.rx_bpm_ddpll = env_flag("UWB_BPM_DDPLL", p.rx_bpm_ddpll);
    p.rx_phr_phase_refine =
        env_flag("UWB_BPM_PHR_REFINE", p.rx_phr_phase_refine);
    p.rx_fused_sfo_track = env_flag("UWB_BPM_SFO", p.rx_fused_sfo_track);
}

} // namespace openuwb::apps
