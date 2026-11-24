# phy/hrp: the HRP demodulation pipeline implementation.
#
# The demodulation stages used to be header-only; Wave 3 (2026-11-21) moved
# the coarse stage/batch entry points into these translation units so that a
# change to an implementation rebuilds one object instead of every consumer.
#
# Scope (2026-11-21): a TX-only build (UWB_DEV_SCOPE=tx) needs only the
# modulator core.  Keeping the receive demodulation TUs out of that archive
# avoids compiling five heavy files that nothing in the TX tree references.
# UWB_ENABLE_RX=ON (the default and every delivery preset) builds the complete
# archive exactly as before.
set(_uwb_phy_hrp_sources
    uwb_hrp_mod_core.cc)
if(UWB_ENABLE_RX)
    list(APPEND _uwb_phy_hrp_sources
        uwb_demod_pipeline.cc
        uwb_demod_simd_utils.cc
        uwb_demod_timing.cc
        uwb_demod_timing_fold.cc
        uwb_demod_cfo_sfd.cc
        uwb_demod_sfd.cc
        uwb_demod_cfo_sampler.cc
        uwb_demod_cfo_ml.cc
        uwb_demod_softchips_payload.cc
        uwb_demod_fec.cc)
endif()

add_library(openuwb_phy_hrp ${UWB_LIBRARY_TYPE} ${_uwb_phy_hrp_sources})
target_link_libraries(openuwb_phy_hrp PUBLIC openuwb_dsp)
set_target_properties(openuwb_phy_hrp PROPERTIES POSITION_INDEPENDENT_CODE ON)
uwb_simd_apply(openuwb_phy_hrp)
