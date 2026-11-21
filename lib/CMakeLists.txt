# Copyright 2011-2020 Free Software Foundation, Inc.
# SPDX-License-Identifier: GPL-3.0-or-later

########################################################################
# Phase 2: one static library per domain, with a strictly one-way DAG.
#
#   openuwb_support      Threads, Volk            (baseline ISA, no UHD)
#   openuwb_dsp       -> support
#   openuwb_phy_hrp   -> dsp
#   openuwb_sic       -> phy_hrp
#   openuwb_sensing   -> phy_hrp                  (INTERFACE, header-only)
#   openuwb_twr_core  -> phy_hrp
#   openuwb_twr_fake_link -> twr_core             (QA transport, not installed)
#   openuwb_twr_phy   -> twr_core, dsp
#   openuwb_radio     -> dsp                      (+ UHD when available)
#   openuwb_runtime   -> phy_hrp, sic, sensing, radio
#
# Nothing below may depend on UHD; only openuwb_radio does, and only when
# ENABLE_UHD_BACKEND did not turn it off.
#
# uwb_core / uwb_native are compatibility aggregates kept for downstream and
# for the install export.
# TODO(openuwb-compat-removal): drop the aggregates once consumers link the
# domain targets directly.
########################################################################

add_subdirectory(support)
add_subdirectory(dsp)
add_subdirectory(phy/hrp)
if(UWB_ENABLE_SIC)
    add_subdirectory(sic)
endif()
# sensing is header-only (INTERFACE) and is also referenced by the HRP
# modulator QA, so configure it in every scope at no compile cost.
add_subdirectory(sensing)
if(UWB_ENABLE_TWR)
    add_subdirectory(twr)
endif()
if(UWB_ENABLE_RADIO)
    add_subdirectory(radio)
endif()
if(UWB_ENABLE_RADIO)
    add_subdirectory(runtime)
endif()

########################################################################
# Compatibility aggregates.
#
# A static archive needs a source, so both aggregates compile one empty TU: all
# of their content is the PUBLIC link interface, which also carries the SIMD
# compile flags so a header-only consumer instantiates the same kernels.
########################################################################
set(OPENUWB_COMPAT_TU ${CMAKE_CURRENT_BINARY_DIR}/openuwb_compat_empty.cc)
file(WRITE ${OPENUWB_COMPAT_TU}
    "// Phase 2 compatibility aggregate: see lib/CMakeLists.txt.\n")

add_library(uwb_core ${UWB_LIBRARY_TYPE} ${OPENUWB_COMPAT_TU})
set(_uwb_core_deps openuwb_support openuwb_dsp openuwb_phy_hrp)
if(UWB_ENABLE_SIC)
    list(APPEND _uwb_core_deps openuwb_sic)
endif()
if(UWB_ENABLE_RADAR)
    list(APPEND _uwb_core_deps openuwb_sensing)
endif()
target_link_libraries(uwb_core PUBLIC ${_uwb_core_deps})
set_target_properties(uwb_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
uwb_simd_apply(uwb_core)

add_library(uwb_native ${UWB_LIBRARY_TYPE} ${OPENUWB_COMPAT_TU})
set(_uwb_native_deps "")
if(TARGET openuwb_runtime)
    list(APPEND _uwb_native_deps openuwb_runtime)
endif()
if(TARGET openuwb_radio)
    list(APPEND _uwb_native_deps openuwb_radio)
endif()
if(_uwb_native_deps)
    target_link_libraries(uwb_native PUBLIC ${_uwb_native_deps})
endif()
set_target_properties(uwb_native PROPERTIES POSITION_INDEPENDENT_CODE ON)
uwb_simd_apply(uwb_native)

# Stage 2+: shared lightweight STL PCH prototype for all PIC+SIMD libraries
# (openuwb_support uses baseline ISA without uwb_simd_apply and has only two tiny
# TUs; uwb_core/uwb_native compile only openuwb_compat_empty.cc, so skip them).
if(UWB_ENABLE_PCH)
    set(OPENUWB_LIB_PCH_TU ${CMAKE_CURRENT_BINARY_DIR}/openuwb_lib_pch_proto.cc)
    file(WRITE ${OPENUWB_LIB_PCH_TU}
        "// Shared prototype TU for lib/ STL precompiled headers.\n"
        "static_assert(true, \"lib pch proto carries no code\");\n")
    add_library(uwb_lib_pch_simd STATIC ${OPENUWB_LIB_PCH_TU})
    set_target_properties(uwb_lib_pch_simd PROPERTIES POSITION_INDEPENDENT_CODE ON)
    uwb_simd_apply(uwb_lib_pch_simd)
    target_precompile_headers(uwb_lib_pch_simd PUBLIC ${UWB_STL_PCH_HEADERS})
endif()

foreach(_lib IN ITEMS openuwb_dsp openuwb_phy_hrp
        openuwb_sic openuwb_twr_core openuwb_twr_phy openuwb_radio
        openuwb_runtime)
    uwb_apply_stl_pch(${_lib})
endforeach()

########################################################################
# Install and export the native libraries.
#
# The exported names drop the `openuwb_` target prefix, so consumers write
# openuwb::support, openuwb::dsp, openuwb::phy_hrp, openuwb::sic,
# openuwb::twr_core, openuwb::twr_phy, openuwb::radio, openuwb::runtime.
# The two compatibility aggregates keep their own names:
# openuwb::uwb_core and openuwb::uwb_native.
#
# A scoped development build (UWB_DEV_SCOPE != all) sets
# CMAKE_SKIP_INSTALL_RULES, so only the export-name bookkeeping below runs.
########################################################################
set_target_properties(openuwb_support PROPERTIES EXPORT_NAME support)
set_target_properties(openuwb_dsp PROPERTIES EXPORT_NAME dsp)
set_target_properties(openuwb_phy_hrp PROPERTIES EXPORT_NAME phy_hrp)
if(TARGET openuwb_sic)
    set_target_properties(openuwb_sic PROPERTIES EXPORT_NAME sic)
endif()
if(TARGET openuwb_sensing)
    set_target_properties(openuwb_sensing PROPERTIES EXPORT_NAME sensing)
endif()
if(TARGET openuwb_twr_core)
    set_target_properties(openuwb_twr_core PROPERTIES EXPORT_NAME twr_core)
endif()
if(TARGET openuwb_twr_phy)
    set_target_properties(openuwb_twr_phy PROPERTIES EXPORT_NAME twr_phy)
endif()
if(TARGET openuwb_radio)
    set_target_properties(openuwb_radio PROPERTIES EXPORT_NAME radio)
endif()
if(TARGET openuwb_runtime)
    set_target_properties(openuwb_runtime PROPERTIES EXPORT_NAME runtime)
endif()

set(_uwb_install_libs
    openuwb_support openuwb_dsp openuwb_phy_hrp uwb_core uwb_native)
foreach(_optional_lib IN ITEMS
        openuwb_sic openuwb_twr_core openuwb_twr_phy openuwb_radio
        openuwb_runtime)
    if(TARGET ${_optional_lib})
        list(APPEND _uwb_install_libs ${_optional_lib})
    endif()
endforeach()
install(TARGETS ${_uwb_install_libs}
    EXPORT uwb-export
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
# INTERFACE libraries take no destination.
if(TARGET openuwb_sensing)
    install(TARGETS openuwb_sensing EXPORT uwb-export)
endif()
install(EXPORT uwb-export FILE OpenUWBTargets.cmake NAMESPACE openuwb::
    DESTINATION ${OPENUWB_CMAKE_DIR})

########################################################################
# UHD backend dry-run CLI (Radar Step 11)
#
# Pure-contract printer: includes no UHD header and opens no device. Registered
# in CTest in both build configurations (UHD ON and OFF) when the radar scope
# is configured.
########################################################################
if(UWB_ENABLE_RADAR)
add_executable(uwb_uhd_burst_dry_run radio/uwb_uhd_dry_run.cc)
target_compile_definitions(uwb_uhd_burst_dry_run PRIVATE UWB_NO_UHD_LINK=1)
target_link_libraries(uwb_uhd_burst_dry_run PRIVATE uwb_core)

add_test(NAME uwb_uhd_burst_dry_run
    COMMAND uwb_uhd_burst_dry_run
        --device-args addr=192.168.10.2
        --rate 737280000
        --tx-channel 0 --rx-channel 1
        --clock-source internal --time-source internal
        --pri-num 3686400 --pri-den 1
        --pre-guard-ticks 1475
        --tx-samples 140982 --rx-samples 55591
        --max-fragment-size 65536
        --t0-ticks 36864000
        --schedule-index 0)
set_tests_properties(uwb_uhd_burst_dry_run PROPERTIES LABELS "radio_tx")

add_test(NAME uwb_uhd_burst_dry_run_rejects_invalid_config
    COMMAND uwb_uhd_burst_dry_run --pri-num 0)
set_tests_properties(uwb_uhd_burst_dry_run_rejects_invalid_config
    PROPERTIES WILL_FAIL TRUE LABELS "radio_tx")
endif()
