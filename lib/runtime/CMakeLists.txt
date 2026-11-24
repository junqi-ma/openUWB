# runtime: the UHD SC16 scheduled pipeline and its sources.
#
# Scope: the pipeline TU itself carries the Sc16WindowRingPool and
# UhdSc16Pipeline classes the RX QA drives, so it is always compiled.
# `uwb_scheduled_extractor_core.cc` is only referenced by the radar/scheduled
# path, so RX-only builds (UWB_ENABLE_RADAR=OFF) exclude it.  The extractor
# has no sensing dependency either; the conditional keeps the RX archive lean
# without dropping the pipeline tests.
set(_uwb_runtime_sources
    scheduled_pipeline.cc
    uwb_live_packet_processor.cc
    uwb_preamble_xcorr_acquisition.cc
    uwb_target_lock_core.cc
    uwb_timed_rx_plan.cc)
set(_uwb_runtime_deps
    openuwb_phy_hrp openuwb_sic openuwb_radio)
if(UWB_ENABLE_RADAR)
    list(APPEND _uwb_runtime_sources
        uwb_scheduled_extractor_core.cc)
    list(APPEND _uwb_runtime_deps openuwb_sensing)
endif()

add_library(openuwb_runtime ${UWB_LIBRARY_TYPE} ${_uwb_runtime_sources})
target_link_libraries(openuwb_runtime PUBLIC ${_uwb_runtime_deps})
set_target_properties(openuwb_runtime PROPERTIES POSITION_INDEPENDENT_CODE ON)
uwb_simd_apply(openuwb_runtime)
