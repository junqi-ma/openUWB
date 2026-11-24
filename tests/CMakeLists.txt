# Copyright 2011-2020 Free Software Foundation, Inc.
# SPDX-License-Identifier: GPL-3.0-or-later

########################################################################
# Phase 1 tests tree.
#
# The unit QA used to sit inside lib/ and the benchmarks/integration harness
# inside apps/, so a "library" directory was half tests and an "apps" directory
# was half measurement rigs.  They live here now (still under gr-uwb/ during
# Phase 1; Phase 3 lifts this subtree to <repo>/tests/).
#
# Registration stays next to the code it exercises: lib/ keeps the UHD dry-run
# CLI case, apps/ keeps the installed-binary console contract and the tools/sic
# script self-checks, because those sources did not move.
#
# Target names are load-bearing.  uwb_register_qa() names each executable (and
# therefore each CTest case) "uwb_<file>.cc"; several uwb_simd_apply() and
# target_compile_definitions() calls address those names, and tools/ scripts
# match them.  Do not "tidy" the .cc suffix away.
########################################################################

set(UWB_UNIT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/unit)
set(UWB_BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
set(UWB_INTEGRATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/integration)
set(UWB_APPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../apps)

# The QA sources resolve fixtures as "../../../testdata" relative to the
# working directory, so keep the run directory exactly three levels below the
# repository root (build/<preset>/qa_work).  qa_uwb_demod_core.cc has no
# UWB_TESTDATA_DIR fallback and is the one case that depends on this depth.
# TODO(openuwb-arch-phase3): replace those relative paths with
# ${UWB_TESTDATA_DIR} once the tests tree moves to the repository root.
set(UWB_QA_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/qa_work")
file(MAKE_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")

########################################################################
# Stage 1: precompiled headers for the QA executables.
#
# 63+ QA executables are all Boost.Test frontends over the same header
# set; PCH removes the repeated frontend parse of <boost/test/unit_test.hpp>,
# the STL headers and the HRP profile header.  Two prototype PCH targets
# carry the shared header set:
#   uwb_qa_pch_plain : bound to the plain QA targets (no explicit ISA
#                      flags; GCC/Clang require the PCH and the consuming
#                      TU to carry the same -D macros and -m flags per the
#                      UwbSimd.cmake compatibility rule).
#   uwb_qa_pch_simd  : bound to the targets that call uwb_simd_apply();
#                      it is compiled with the production ISA flags.
# The REUSE_FROM bindings happen after the uwb_simd_apply() calls (the
# SIMD target set is not known before then), so uwb_register_qa() binds
# the plain PCH first and the SIMD rebinding runs right after.
########################################################################
set(UWB_QA_PCH_PUBLIC_HEADERS
    <boost/test/unit_test.hpp>
    <vector> <complex> <algorithm> <functional> <memory> <string>
    <cstdint> <cstddef> <new> <stdexcept>
    <openuwb/phy/hrp/common/profile.h>)

if(UWB_ENABLE_PCH)
    # The two PCH prototype targets need at least one source file
    # (add_library(STATIC) without sources is rejected); this shared empty
    # TU carries no code, the PCH comes from the header lists below.
    set(OPENUWB_QA_PCH_PROTO_TU
        ${CMAKE_CURRENT_SOURCE_DIR}/unit/openuwb_qa_pch_proto.cc)
    add_library(uwb_qa_pch_plain STATIC ${OPENUWB_QA_PCH_PROTO_TU})
    target_include_directories(uwb_qa_pch_plain PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    # GCC validates a PCH against the consuming TU's macro set and -m flags.
    # The QA executables define the two Boost.Test macros and (through the
    # library interface, and now explicitly) the production ISA flags, so the
    # prototype must carry exactly the same values or every include silently
    # falls back to a plain parse (-Winvalid-pch).
    target_compile_definitions(uwb_qa_pch_plain PUBLIC
        BOOST_TEST_DYN_LINK BOOST_TEST_MAIN)
    target_link_libraries(uwb_qa_pch_plain PRIVATE Boost::unit_test_framework)
    target_precompile_headers(uwb_qa_pch_plain
        PUBLIC ${UWB_QA_PCH_PUBLIC_HEADERS})
    uwb_simd_apply(uwb_qa_pch_plain)

    add_library(uwb_qa_pch_simd STATIC ${OPENUWB_QA_PCH_PROTO_TU})
    target_include_directories(uwb_qa_pch_simd PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    target_compile_definitions(uwb_qa_pch_simd PUBLIC
        BOOST_TEST_DYN_LINK BOOST_TEST_MAIN)
    target_link_libraries(uwb_qa_pch_simd PRIVATE Boost::unit_test_framework)
    target_precompile_headers(uwb_qa_pch_simd
        PUBLIC ${UWB_QA_PCH_PUBLIC_HEADERS} <immintrin.h>)
    uwb_simd_apply(uwb_qa_pch_simd)

    # Filled after the uwb_simd_apply() calls below: the names of the QA
    # targets compiled through the production ISA chain.
    set(UWB_QA_PCH_SIMD_TARGETS "")
endif()

# ---------------------------------------------------------------------------
# Stage 2: per-domain link libraries.
#
# Historically every unit QA linked `uwb_native uwb_twr_phy`, which bound 65
# executables to dsp/phy_hrp/twr and made a one-line DSP edit relink the whole
# tree.  The table below links the smallest archive that owns the symbols the
# file actually uses.  A per-file override covers the few targets that reach
# across the default domain boundary.
# ---------------------------------------------------------------------------
set(UWB_QA_LINK_OVERRIDES
    # support/: uwb_defaults also includes the radar-meta vocabulary.
    qa_uwb_defaults.cc=openuwb_phy_hrp
    # dsp/: the AVX-512 SIC kernel QA includes the PHY matched-filter headers.
    qa_uwb_demod_sic_avx512.cc=openuwb_phy_hrp
    # phy/: the pure HRP modulator QA also exercises the header-only radar SFD core.
    qa_uwb_hrp_mod_core.cc=openuwb_phy_hrp+openuwb_sensing
    # phy/: TX->RX demodulator loopback QA runs whenever RX is enabled.
    qa_uwb_hrp_mod_rx_loopback.cc=openuwb_phy_hrp
    # phy/: the refactor contract walks the live packet processor (runtime).
    qa_uwb_refactor_contracts.cc=openuwb_runtime
    # phy/: the fused CMF/PHR QAs drive the fullpacket pipeline (runtime + sic).
    qa_uwb_rx_fused_cmf.cc=openuwb_runtime
    qa_uwb_rx_fused_phr.cc=openuwb_runtime
    # twr/: only the two native-rate PHY QAs link the PHY loop.
    qa_uwb_twr_native_phy.cc=openuwb_twr_phy
    qa_uwb_twr_phy_impairment.cc=openuwb_twr_phy
    # radio/: the tone-source decorator QA feeds the SC16 pipeline.
    qa_uwb_fixed_tone_cancel_source.cc=openuwb_radio+openuwb_runtime
)

# Map a QA source path to its minimal library set.
function(uwb_qa_link_libs_for qa_source out_var)
    get_filename_component(_name "${qa_source}" NAME)
    foreach(_entry IN LISTS UWB_QA_LINK_OVERRIDES)
        string(FIND "${_entry}" "=" _eq)
        string(SUBSTRING "${_entry}" 0 ${_eq} _key)
        math(EXPR _val_start "${_eq} + 1")
        string(SUBSTRING "${_entry}" ${_val_start} -1 _value)
        if(_key STREQUAL _name)
            string(REPLACE "+" ";" _value "${_value}")
            set(${out_var} ${_value} PARENT_SCOPE)
            return()
        endif()
    endforeach()

    get_filename_component(_dir "${qa_source}" DIRECTORY)
    get_filename_component(_dir "${_dir}" NAME)
    if(_dir STREQUAL "support")
        set(_libs openuwb_support)
    elseif(_dir STREQUAL "dsp")
        set(_libs openuwb_dsp)
    elseif(_dir STREQUAL "phy")
        set(_libs openuwb_phy_hrp)
    elseif(_dir STREQUAL "sic")
        set(_libs openuwb_sic)
    elseif(_dir STREQUAL "sensing")
        set(_libs openuwb_sensing)
    elseif(_dir STREQUAL "twr")
        set(_libs openuwb_twr_core)
    elseif(_dir STREQUAL "runtime")
        set(_libs openuwb_runtime)
    elseif(_dir STREQUAL "radio")
        set(_libs openuwb_radio)
    else()
        set(_libs openuwb_phy_hrp)
    endif()
    set(${out_var} "${_libs}" PARENT_SCOPE)
endfunction()

# Which unit QA belongs to the current UWB_DEV_SCOPE.
function(uwb_qa_in_scope qa_source out_var)
    get_filename_component(_name "${qa_source}" NAME)
    get_filename_component(_dir "${qa_source}" DIRECTORY)
    get_filename_component(_dir "${_dir}" NAME)
    if(UWB_DEV_SCOPE STREQUAL "all")
        # Honor the fine-grained switches even in the all scope.
        set(_in_scope TRUE)
        if(_dir STREQUAL "sic" AND NOT UWB_ENABLE_SIC)
            set(_in_scope FALSE)
        elseif(_dir STREQUAL "sensing" AND NOT UWB_ENABLE_RADAR)
            set(_in_scope FALSE)
        elseif(_dir STREQUAL "twr" AND NOT UWB_ENABLE_TWR)
            set(_in_scope FALSE)
        elseif((_dir STREQUAL "runtime" OR _dir STREQUAL "radio")
               AND NOT UWB_ENABLE_RUNTIME)
            set(_in_scope FALSE)
        endif()
        set(${out_var} ${_in_scope} PARENT_SCOPE)
        return()
    endif()
    set(_in_scope FALSE)
    if(UWB_DEV_SCOPE STREQUAL "rx")
        if(_dir STREQUAL "support" OR _dir STREQUAL "dsp" OR _dir STREQUAL "phy"
           OR _dir STREQUAL "sic" OR _dir STREQUAL "runtime"
           OR _dir STREQUAL "radio")
            set(_in_scope TRUE)
        endif()
    elseif(UWB_DEV_SCOPE STREQUAL "tx")
        if(_name STREQUAL "qa_uwb_hrp_mod_core.cc")
            set(_in_scope TRUE)
        endif()
    elseif(UWB_DEV_SCOPE STREQUAL "radar")
        if(_dir STREQUAL "sensing" OR _dir STREQUAL "runtime"
           OR _dir STREQUAL "radio")
            set(_in_scope TRUE)
        endif()
    elseif(UWB_DEV_SCOPE STREQUAL "twr")
        if(_dir STREQUAL "twr")
            set(_in_scope TRUE)
        endif()
    endif()
    set(${out_var} ${_in_scope} PARENT_SCOPE)
endfunction()

function(uwb_register_qa qa_name qa_source)
    cmake_parse_arguments(ARG "NO_DEFAULT_TEST" "" "LINK_LIBS" ${ARGN})
    if(NOT ARG_LINK_LIBS)
        uwb_qa_link_libs_for("${qa_source}" ARG_LINK_LIBS)
    endif()
    add_executable("uwb_${qa_name}" "${qa_source}")
    target_link_libraries("uwb_${qa_name}"
        PRIVATE ${ARG_LINK_LIBS} Boost::unit_test_framework)
    set_target_properties("uwb_${qa_name}" PROPERTIES
        COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK;BOOST_TEST_MAIN")
    # Keep the PCH prototype's -m/-O flags identical to the consuming TU.
    uwb_simd_apply("uwb_${qa_name}")
    if(UWB_ENABLE_PCH AND TARGET uwb_qa_pch_plain)
        target_precompile_headers("uwb_${qa_name}" REUSE_FROM uwb_qa_pch_plain)
    endif()
    if(NOT ARG_NO_DEFAULT_TEST)
        add_test(NAME "uwb_${qa_name}" COMMAND "uwb_${qa_name}"
            WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
    endif()
endfunction()

########################################################################
# Native unit QA.
#
# The list keeps the historical registration order (it was the order of the
# old uwb_core_qa_sources list) so a `ctest -N` diff against the Phase 0
# baseline only shows the tree move, never a reshuffle.
# qa_uwb_twr_timestamp.cc is deliberately absent: it has never had an
# add_executable/add_test, and registering it would change the test count.
########################################################################
set(UWB_UNIT_QA_SOURCES
    ${UWB_UNIT_DIR}/dsp/qa_uwb_autocorr_gate.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_cfo_derotate.cc
    ${UWB_UNIT_DIR}/support/qa_uwb_defaults.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_demod_core_timing_cfo_sfd.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_demod_core_cir_pipeline.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_cfo_ml.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_demod_core_kernel_grid.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_demod_core_fec_joint_ml.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_demod_core_heavy_grid.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_demod_sic_avx512.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_discrete_cir_estimator.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_discrete_cir_reconstruct.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_fixed_tone_cancel_core.cc
    ${UWB_UNIT_DIR}/radio/qa_uwb_fixed_tone_cancel_source.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_fullpacket_cfo_diag.cc
    ${UWB_UNIT_DIR}/runtime/qa_uwb_fullpacket_fast_rx.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_fullpacket_payload_cache.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_fullpacket_sic_core.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_hrp_mod_core.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_hrp_mod_rx_loopback.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_refactor_contracts.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_preamble_sic_core_basics.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_preamble_sic_core_model.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_preamble_sic_core_perf.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_preamble_sic_core_heavy.cc
    ${UWB_UNIT_DIR}/sensing/qa_uwb_radar_cir_core.cc
    ${UWB_UNIT_DIR}/sensing/qa_uwb_radar_cir_estimator.cc
    ${UWB_UNIT_DIR}/sensing/qa_uwb_radar_sfd_core.cc
    ${UWB_UNIT_DIR}/sensing/qa_uwb_radar_timing_core.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_rational_resampler_core.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_rational_resampler_core_65_32.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_rational_resampler_m2a.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_resampler_avx512.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_resampler_burst_opt.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_resampler_transposed.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_rx_fused_cmf.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_rx_fused_phr.cc
    ${UWB_UNIT_DIR}/phy/qa_uwb_rx_phase_sfo.cc
    ${UWB_UNIT_DIR}/runtime/qa_uwb_scheduled_dump_io.cc
    ${UWB_UNIT_DIR}/sic/qa_uwb_sic_unified.cc
    ${UWB_UNIT_DIR}/runtime/qa_uwb_target_lock_core.cc
    ${UWB_UNIT_DIR}/runtime/qa_uwb_timed_rx_plan.cc
    ${UWB_UNIT_DIR}/dsp/qa_uwb_tone_canceller.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_config.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_config_parity.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_frame.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_math.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_native_phy.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_phy_impairment.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_tof_input.cc
    ${UWB_UNIT_DIR}/runtime/qa_uwb_uhd_sc16_pipeline.cc
    # Phase 1 (1C): the header-only diagnostic logger. Appended last so the
    # historical test numbering above is untouched.
    ${UWB_UNIT_DIR}/support/qa_uwb_log.cc
    # Phase 2: the reserved HRP revision vocabulary and the runtime config
    # placeholders.
    ${UWB_UNIT_DIR}/phy/qa_hrp_standard_rev.cc
    ${UWB_UNIT_DIR}/runtime/qa_openuwb_runtime_config.cc
    # SIC live processor QA. Appended last so the historical test numbering
    # above is untouched.
    ${UWB_UNIT_DIR}/runtime/qa_uwb_live_packet_processor_sic.cc
    # Timing refactor (2026-11-21): the support-layer wall-clock primitives.
    # Appended after the historical numbering for the same reason.
    ${UWB_UNIT_DIR}/support/qa_uwb_timing.cc
    # Timing Stage-1 refactor (2026-11-24): K=4 fold fast path, strict gates,
    # K=64 fallback and the warm-path zero-alloc contract. Appended last.
    ${UWB_UNIT_DIR}/phy/qa_uwb_timing_refactor.cc
)

set(UWB_SHARDED_QA_BINARIES
    qa_uwb_demod_core_heavy_grid.cc
    qa_uwb_preamble_sic_core_basics.cc
    qa_uwb_preamble_sic_core_model.cc
    qa_uwb_preamble_sic_core_perf.cc)

foreach(qa_source IN LISTS UWB_UNIT_QA_SOURCES)
    # NAME keeps the ".cc" suffix: it is part of every CTest name.
    get_filename_component(qa_name "${qa_source}" NAME)
    uwb_qa_in_scope("${qa_source}" _in_scope)
    if(_in_scope)
        if("${qa_name}" IN_LIST UWB_SHARDED_QA_BINARIES)
            uwb_register_qa("${qa_name}" "${qa_source}" NO_DEFAULT_TEST)
        else()
            uwb_register_qa("${qa_name}" "${qa_source}")
        endif()
    endif()
endforeach()

########################################################################
# M1-B pure protocol-core tests: link only the protocol and test transport,
# without the UHD runtime.  They need no working directory (no file I/O).
########################################################################
function(uwb_add_pure_twr_test test_name test_source)
    add_executable(${test_name} ${test_source})
    target_include_directories(${test_name}
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    target_link_libraries(${test_name}
        PRIVATE uwb_twr_core uwb_twr_fake_link Boost::unit_test_framework)
    set_target_properties(${test_name} PROPERTIES
        COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK;BOOST_TEST_MAIN")
    uwb_simd_apply(${test_name})
    if(UWB_ENABLE_PCH AND TARGET uwb_qa_pch_plain)
        target_precompile_headers(${test_name} REUSE_FROM uwb_qa_pch_plain)
    endif()
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

if(UWB_ENABLE_TWR)
uwb_add_pure_twr_test(uwb_qa_uwb_twr_protocol_time.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_protocol_time.cc)
uwb_add_pure_twr_test(uwb_qa_uwb_twr_fake_link.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_fake_link.cc)
# The endpoint core unit QA and the two-endpoint end-to-end matrix both
# include twr_m1b_test_support.h from their own directory.
uwb_add_pure_twr_test(uwb_qa_uwb_twr_core.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_core.cc)
uwb_add_pure_twr_test(uwb_qa_uwb_twr_m1b_e2e.cc
    ${UWB_UNIT_DIR}/twr/qa_uwb_twr_m1b_e2e.cc)
endif()

# UHD-gated runtime QA: only compiled when UHD is available; verifies the
# device-unavailable contract without any USRP (skips when one exists).
# This is the only test that disappears in the no-uhd preset.
if(UWB_HAVE_UHD AND UWB_ENABLE_RADIO)
    uwb_register_qa(qa_uwb_uhd_backend_device.cc
        ${UWB_UNIT_DIR}/radio/qa_uwb_uhd_backend_device.cc)
    # Device discovery can take several seconds on transport timeouts.
    set_tests_properties(uwb_qa_uwb_uhd_backend_device.cc PROPERTIES
        TIMEOUT 300)
endif()

########################################################################
# Fixture roots.  UWB_TESTDATA_DIR must stay a *quoted string literal*:
# qa_uwb_uhd_sc16_pipeline.cc concatenates it with adjacent literals.
########################################################################
function(uwb_qa_set_testdata_dir target)
    if(TARGET ${target})
        target_compile_definitions(${target} PRIVATE
            "UWB_TESTDATA_DIR=\"${CMAKE_SOURCE_DIR}/testdata\"")
    endif()
endfunction()

# Radar QAs and the SC16 pipeline read generated fixtures.
foreach(_t IN ITEMS
        uwb_qa_uwb_radar_sfd_core.cc
        uwb_qa_uwb_radar_cir_estimator.cc
        uwb_qa_uwb_uhd_sc16_pipeline.cc
        uwb_qa_uwb_radar_timing_core.cc
        uwb_qa_uwb_radar_cir_core.cc
        uwb_qa_uwb_rational_resampler_core_65_32.cc
        uwb_qa_uwb_rational_resampler_m2a.cc
        uwb_qa_uwb_live_packet_processor_sic.cc
        uwb_qa_uwb_resampler_transposed.cc
        uwb_qa_uwb_twr_native_phy.cc
        uwb_qa_uwb_hrp_mod_core.cc
        uwb_qa_uwb_hrp_mod_rx_loopback.cc
        uwb_qa_uwb_preamble_sic_core_basics.cc
        uwb_qa_uwb_preamble_sic_core_model.cc
        uwb_qa_uwb_preamble_sic_core_perf.cc
        uwb_qa_uwb_twr_config_parity.cc
        uwb_qa_uwb_twr_frame.cc
        uwb_qa_uwb_twr_math.cc)
    uwb_qa_set_testdata_dir(${_t})
endforeach()

# Every QA target now receives the production ISA flags inside
# uwb_register_qa().  This list only selects the PCH prototype that also
# precompiles <immintrin.h>.
set(UWB_QA_PCH_SIMD_TARGETS
    uwb_qa_uwb_fullpacket_fast_rx.cc
    uwb_qa_uwb_live_packet_processor_sic.cc
    uwb_qa_uwb_rx_fused_cmf.cc
    uwb_qa_uwb_rx_fused_phr.cc
    uwb_qa_uwb_uhd_sc16_pipeline.cc
    uwb_qa_uwb_rational_resampler_core.cc
    uwb_qa_uwb_resampler_burst_opt.cc
    uwb_qa_uwb_resampler_transposed.cc
    uwb_qa_uwb_rational_resampler_core_65_32.cc
    uwb_qa_uwb_cfo_ml.cc)

# Stage 1: SIMD-chain rebinding.  The nine QA targets above were first
# bound to the plain PCH by uwb_register_qa(); they are compiled through
# the production ISA chain, so rebind them to the SIMD PCH.
if(UWB_ENABLE_PCH AND TARGET uwb_qa_pch_simd)
    foreach(_t IN LISTS UWB_QA_PCH_SIMD_TARGETS)
        if(TARGET ${_t})
            target_precompile_headers(${_t} REUSE_FROM uwb_qa_pch_simd)
        endif()
    endforeach()
endif()

########################################################################
# Integration tests.
#
# install_consumer/ proves an installed prefix links and runs with no source
# tree, no build tree and no GNU Radio/PMT dependency.  The shell case pins the
# pipeline runner's exit status contract.  The Python tests pin the TWR config
# parity and the offline two-endpoint demo/verifier.
########################################################################
set(UWB_NATIVE_INSTALL_PREFIX "" CACHE PATH
    "Installed UWB prefix for the consumer test (empty skips the test)")

if(UNIX AND UWB_DEV_SCOPE STREQUAL "all")
    add_test(NAME uwb_qa_native_install_consumer
        COMMAND ${UWB_INTEGRATION_DIR}/install_consumer/run_native_install_consumer.sh
                --prefix "${UWB_NATIVE_INSTALL_PREFIX}"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(uwb_qa_native_install_consumer PROPERTIES
        TIMEOUT 300
        SKIP_RETURN_CODE 77
        RUN_SERIAL TRUE)
    message(STATUS "UWB CTest: registered uwb_qa_native_install_consumer "
                   "(prefix='${UWB_NATIVE_INSTALL_PREFIX}'; empty => SKIP)")
endif()

if(TARGET openuwb_pipeline_runner)
    add_test(NAME uwb_qa_sc16_pipeline_runner_exit
        COMMAND sh ${UWB_INTEGRATION_DIR}/qa_uwb_uhd_sc16_pipeline_runner_exit.sh
                --runner $<TARGET_FILE:uwb_uhd_sc16_pipeline_runner>
                --fixture ${CMAKE_SOURCE_DIR}/testdata/uhd_sc16/sim_preamble_sc16_737m28_rate.bin
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(uwb_qa_sc16_pipeline_runner_exit PROPERTIES TIMEOUT 120)
endif()

if(Python3_Interpreter_FOUND AND UWB_ENABLE_TWR)
    # The parity suite derives the repository root from its own location; if
    # that ever regresses, three cases silently become SKIPs, so run it with
    # -v through the test and check the per-case output, not just the exit code.
    add_test(NAME uwb_qa_twr_config_py
        COMMAND ${Python3_EXECUTABLE} ${UWB_INTEGRATION_DIR}/test_twr_config.py
        WORKING_DIRECTORY ${UWB_INTEGRATION_DIR})
    set_property(TEST uwb_qa_twr_config_py APPEND PROPERTY
        ENVIRONMENT
            "PYTHONPATH="
            "PYTHONDONTWRITEBYTECODE=1"
            "UWB_TWR_PYTHON_QA=ctest")
    set_tests_properties(uwb_qa_twr_config_py PROPERTIES
        TIMEOUT 300
        # A hand-edited checkout or a stale installed uwb package must not make
        # this test's result depend on run order, so it never runs alongside
        # anything else.
        RUN_SERIAL TRUE)
    message(STATUS "UWB CTest: registered uwb_qa_twr_config_py "
                   "(tests/integration/test_twr_config.py)")

    add_test(NAME uwb_qa_twr_fake_demo_py
        COMMAND ${Python3_EXECUTABLE} ${UWB_INTEGRATION_DIR}/test_twr_fake_demo.py
        WORKING_DIRECTORY ${UWB_INTEGRATION_DIR})
    set_property(TEST uwb_qa_twr_fake_demo_py APPEND PROPERTY
        ENVIRONMENT
            "PYTHONPATH="
            "PYTHONDONTWRITEBYTECODE=1"
            "UWB_TWR_FAKE_DEMO=$<TARGET_FILE:twr_fake_demo>")
    set_tests_properties(uwb_qa_twr_fake_demo_py PROPERTIES
        TIMEOUT 300 RUN_SERIAL TRUE)

    # The independent verifier re-derives the ToF from the frame fields and the
    # transport model without calling the C++ at all.  It needs a real
    # request/output pair, generated in a scratch dir on the fly.
    add_test(NAME uwb_qa_twr_m1b_verify_py
        COMMAND sh -c
            "d=$(mktemp -d ${CMAKE_CURRENT_BINARY_DIR}/m1b_verify.XXXXXX) && \
${Python3_EXECUTABLE} ${UWB_INTEGRATION_DIR}/twr_fake_demo.py --scenario ds_unity \
    --request-out $d/req.json --output $d/out.json && \
${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/twr/verify_m1_b_protocol.py \
    --request $d/req.json --output $d/out.json; rc=$?; rm -rf $d; exit $rc")
    set_property(TEST uwb_qa_twr_m1b_verify_py APPEND PROPERTY
        ENVIRONMENT
            "PYTHONPATH="
            "PYTHONDONTWRITEBYTECODE=1"
            "UWB_TWR_FAKE_DEMO=$<TARGET_FILE:twr_fake_demo>")
    set_tests_properties(uwb_qa_twr_m1b_verify_py PROPERTIES
        TIMEOUT 300 RUN_SERIAL TRUE)

    # M2-A: run the offline native-rate demo and re-verify the manifest from
    # the frame bytes without the C++ (the demo binary stays in apps/).
    add_test(NAME uwb_qa_twr_m2a_native_phy_verify_py
        COMMAND sh -c
            "d=$(mktemp -d ${CMAKE_CURRENT_BINARY_DIR}/m2a_verify.XXXXXX) && \
$<TARGET_FILE:twr_native_phy_demo> --frame poll --native-rate 737280000 --iq cf32 \
    --output $d/out.json && \
${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/twr/verify_m2_a_native_phy.py \
    --output $d/out.json --self-test; rc=$?; rm -rf $d; exit $rc")
    set_property(TEST uwb_qa_twr_m2a_native_phy_verify_py APPEND PROPERTY
        ENVIRONMENT
            "PYTHONPATH="
            "PYTHONDONTWRITEBYTECODE=1")
    set_tests_properties(uwb_qa_twr_m2a_native_phy_verify_py PROPERTIES
        TIMEOUT 300 RUN_SERIAL TRUE)
endif()

########################################################################
# Benchmarks.  They are measurement rigs, not production binaries: they are
# never installed.  The two with a CTest case keep it (a smoke gate, not a
# timing assertion).
########################################################################

# HRP TX modulator timing (P3: 64-SYNC PSDU<=127 B >= 1000 pkt/s)
if(UWB_ENABLE_BENCHMARKS OR UWB_DEV_SCOPE STREQUAL "tx")
add_executable(benchmark_hrp_mod ${UWB_BENCH_DIR}/benchmark_hrp_mod.cc)
target_link_libraries(benchmark_hrp_mod
    uwb_core
    Volk::volk
)
add_test(NAME uwb_benchmark_hrp_mod
    COMMAND benchmark_hrp_mod 4 32)
set_tests_properties(uwb_benchmark_hrp_mod PROPERTIES TIMEOUT 120)
endif()

if(UWB_ENABLE_BENCHMARKS)
add_executable(benchmark_energy_gate_formats
    ${UWB_BENCH_DIR}/benchmark_energy_gate_formats.cc)
target_link_libraries(benchmark_energy_gate_formats
    uwb_core
)

# M2-A closeout: the A13 allocation/latency observation tool (G0 A.4).
# An observation tool, not a timing gate; links only uwb_twr_phy.
if(UWB_ENABLE_TWR)
add_executable(twr_m2a_bench ${UWB_BENCH_DIR}/twr_m2a_bench.cc)
target_include_directories(twr_m2a_bench PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(twr_m2a_bench PRIVATE
    uwb_twr_phy Volk::volk Threads::Threads)
target_compile_definitions(twr_m2a_bench PRIVATE
    "UWB_TESTDATA_DIR=\"${CMAKE_SOURCE_DIR}/testdata\"")
endif()

# 2026-11-11 round 2: pure-core resampler benchmark.
add_executable(benchmark_sc16_burst_core
    ${UWB_BENCH_DIR}/benchmark_sc16_burst_core.cc)
target_link_libraries(benchmark_sc16_burst_core uwb_native Threads::Threads)
uwb_simd_apply(benchmark_sc16_burst_core)
target_compile_definitions(benchmark_sc16_burst_core PRIVATE
    "UWB_TESTDATA_DIR=\"${CMAKE_SOURCE_DIR}/testdata\"")

# 2026-11-15: fixed-frequency CW cancellation.  The benchmark measures the
# production core (two-pass vs fused, ABBA); the replay driver
# (fixed_tone_cancel_replay) stays in apps/ as an operator tool.
add_executable(benchmark_fixed_tone_cancel
    ${UWB_BENCH_DIR}/benchmark_fixed_tone_cancel.cc)
target_link_libraries(benchmark_fixed_tone_cancel uwb_native Threads::Threads)
uwb_simd_apply(benchmark_fixed_tone_cancel)

# 2026-11-11 round 2: steady-state pipeline benchmark driver.  It includes
# pipeline_bench_setup.h / pipeline_result_json.h from apps/, where they stay
# because uwb_uhd_sc16_pipeline_runner.cc uses them too.
add_executable(benchmark_sc16_pipeline_steady
    ${UWB_BENCH_DIR}/benchmark_sc16_pipeline_steady.cc)
target_link_libraries(benchmark_sc16_pipeline_steady uwb_native Threads::Threads)
target_include_directories(benchmark_sc16_pipeline_steady PRIVATE ${UWB_APPS_DIR})
uwb_simd_apply(benchmark_sc16_pipeline_steady)
target_compile_definitions(benchmark_sc16_pipeline_steady PRIVATE
    "UWB_TESTDATA_DIR=\"${CMAKE_SOURCE_DIR}/testdata\"")
if(UWB_HAVE_UHD)
    target_link_libraries(benchmark_sc16_pipeline_steady UHD::UHD)
endif()
endif() # UWB_ENABLE_BENCHMARKS

########################################################################
# Historical benchmark / diagnostic EXPERIMENT rigs (UWB_BUILD_EXPERIMENTS).
########################################################################
if(UWB_BUILD_EXPERIMENTS AND UWB_ENABLE_BENCHMARKS)
    add_executable(benchmark_fullpacket_compare
        ${UWB_BENCH_DIR}/benchmark_fullpacket_compare.cc)
    target_link_libraries(benchmark_fullpacket_compare uwb_native)
    # fullpacket_rx_pipeline.h stays in apps/ (shared with production).
    target_include_directories(benchmark_fullpacket_compare PRIVATE ${UWB_APPS_DIR})

    # Mature RX uses the existing AVX2/FMA FIR and compiler vectorization.
    # Keep the portable executable for same-source, same-library paired
    # measurement.
    add_executable(benchmark_fullpacket_compare_fast_rx
        ${UWB_BENCH_DIR}/benchmark_fullpacket_compare.cc)
    target_link_libraries(benchmark_fullpacket_compare_fast_rx uwb_native)
    target_include_directories(benchmark_fullpacket_compare_fast_rx PRIVATE
        ${UWB_APPS_DIR})
    uwb_simd_apply(benchmark_fullpacket_compare_fast_rx)

    # Isolated generic-ISA CFO replay: diagnostic only, production kernel
    # unchanged.
    add_executable(benchmark_fullpacket_cfo_diag
        ${UWB_BENCH_DIR}/benchmark_fullpacket_cfo_diag.cc)
    target_link_libraries(benchmark_fullpacket_cfo_diag uwb_native)
    target_compile_options(benchmark_fullpacket_cfo_diag PRIVATE -fPIC)

    # Offline-only copy of reconstruction, split by waveform fields for timing.
    add_executable(benchmark_fullpacket_segment_diag
        ${UWB_BENCH_DIR}/benchmark_fullpacket_segment_diag.cc
        ${UWB_BENCH_DIR}/fullpacket_segment_diagnostics.cc)
    target_link_libraries(benchmark_fullpacket_segment_diag uwb_native)
    target_compile_options(benchmark_fullpacket_segment_diag PRIVATE -fPIC)
    add_test(NAME uwb_qa_fullpacket_segment_diag
        COMMAND benchmark_fullpacket_segment_diag --self-test)
    set_tests_properties(uwb_qa_fullpacket_segment_diag PROPERTIES
        TIMEOUT 60 RUN_SERIAL TRUE)

    # Controlled tap-count sweep using the unchanged production discrete
    # backend.
    add_executable(benchmark_discrete_cir_taps
        ${UWB_BENCH_DIR}/benchmark_discrete_cir_taps.cc)
    target_link_libraries(benchmark_discrete_cir_taps uwb_native)
    target_compile_options(benchmark_discrete_cir_taps PRIVATE -fPIC)

    # Unified SIC timing benchmark.
    add_executable(bench_unified_timing
        ${UWB_BENCH_DIR}/bench_unified_timing.cc)
    target_link_libraries(bench_unified_timing uwb_native)
    uwb_simd_apply(bench_unified_timing)

    # Phase-2 operator-level benchmark: one demod/SIC kernel at a time, with a
    # bit-identity gate against the AVX2 reference for every AVX-512 candidate.
    add_executable(bench_demod_sic_kernels
        ${UWB_BENCH_DIR}/bench_demod_sic_kernels.cc)
    target_link_libraries(bench_demod_sic_kernels uwb_native)
    uwb_simd_apply(bench_demod_sic_kernels)
endif()

########################################################################
# Sharded CTest registrations (one compiled binary, parallel Boost cases).
# Registered BEFORE the label foreach() below so labels apply cleanly.
########################################################################
if(TARGET uwb_qa_uwb_demod_core_heavy_grid.cc)
    foreach(_snr IN ITEMS m10db m12db m14db noise_scalar noise_fast)
        add_test(NAME "uwb_qa_uwb_demod_core_heavy_grid.cc_${_snr}"
            COMMAND uwb_qa_uwb_demod_core_heavy_grid.cc
                    "--run_test=test_demod_core_low_snr_sfd_discrete_grid_${_snr}"
            WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
        set_tests_properties("uwb_qa_uwb_demod_core_heavy_grid.cc_${_snr}" PROPERTIES
            LABELS "heavy")
    endforeach()
endif()
if(TARGET uwb_qa_uwb_preamble_sic_core_basics.cc)
    add_test(NAME uwb_qa_uwb_preamble_sic_core_basics.cc_1
        COMMAND uwb_qa_uwb_preamble_sic_core_basics.cc
                "--run_test=test_sic_core_integer_single_path,test_sic_core_integer_four_path,test_sic_core_fractional,test_sic_core_cfo,test_sic_core_roi_narrow_vs_wide,test_sic_core_comm_outside_sensing_roi_bypasses,test_sic_core_wrong_code_no_commit,test_sic_core_bypass_bit_exact,test_sic_core_oracle_truth_isolation,test_sic_core_timings_no_double_count,test_sic_core_gain_condition_semantics,test_sic_core_fractional_selection_differential"
        WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
    add_test(NAME uwb_qa_uwb_preamble_sic_core_basics.cc_2
        COMMAND uwb_qa_uwb_preamble_sic_core_basics.cc
                "--run_test=test_sic_core_gain_fit_ill_conditioned,test_sic_core_adjacent_path_conditioning,test_sic_core_union_low_edge_negative_delay,test_sic_s1_prefix_isolation,test_sic_s1_cfo_isolation,test_sic_s1_full_support_and_tail,test_sic_s1_chunk_invariance"
        WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
endif()
if(TARGET uwb_qa_uwb_preamble_sic_core_model.cc)
    add_test(NAME uwb_qa_uwb_preamble_sic_core_model.cc_1
        COMMAND uwb_qa_uwb_preamble_sic_core_model.cc
                "--run_test=test_sic_s1_extrapolation_mismatch,test_sic_s1_bypass_untrusted_model,test_sic_s1_degenerate_nonfinite_model_bypass,test_sic_s1_prefix_guard_page"
        WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
    add_test(NAME uwb_qa_uwb_preamble_sic_core_model.cc_2
        COMMAND uwb_qa_uwb_preamble_sic_core_model.cc
                "--run_test=test_sic_core_golden_crosscheck,test_sic_core_k_sweep,test_sic_core_noise_bypass_is_clean"
        WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
endif()
if(TARGET uwb_qa_uwb_preamble_sic_core_perf.cc)
    add_test(NAME uwb_qa_uwb_preamble_sic_core_perf.cc_1
        COMMAND uwb_qa_uwb_preamble_sic_core_perf.cc
                "--run_test=test_sic_perf_e0_cache_matches_uncached,test_sic_perf_e0_k_and_flags,test_sic_perf_recon_cache_and_diag"
        WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
    add_test(NAME uwb_qa_uwb_preamble_sic_core_perf.cc_2
        COMMAND uwb_qa_uwb_preamble_sic_core_perf.cc
                "--run_test=test_sic_perf_r01_self_energy_geometry,test_sic_core_fit_reps_sweep,test_sic_core_energy_prefix_fallback"
        WORKING_DIRECTORY "${UWB_QA_WORKING_DIRECTORY}")
endif()

########################################################################
# CTest grouping labels (non-destructive: keeps all targets and order).
#   core / rx_sic : essential openuwb_rx, SIC, PHY, runtime & resampler QA
#   twr           : two-way ranging unit + Python integration QA
#   sensing       : monostatic radar CIR/SFD/timing QA
#   dsp_extra     : secondary resampler & tone-cancellation QA
#   legacy_diag   : historical diagnostics & experiment self-tests
#   infra         : contract/config/log/install/benchmark smoke QA
########################################################################
foreach(_t IN ITEMS
    uwb_qa_uwb_live_packet_processor_sic.cc
    uwb_qa_uwb_uhd_sc16_pipeline.cc
    uwb_qa_sc16_pipeline_runner_exit
    uwb_qa_uwb_sic_unified.cc
    uwb_qa_uwb_discrete_cir_estimator.cc
    uwb_qa_uwb_discrete_cir_reconstruct.cc
    uwb_qa_uwb_fullpacket_sic_core.cc
    uwb_qa_uwb_preamble_sic_core_basics.cc_1
    uwb_qa_uwb_preamble_sic_core_basics.cc_2
    uwb_qa_uwb_preamble_sic_core_model.cc_1
    uwb_qa_uwb_preamble_sic_core_model.cc_2
    uwb_qa_uwb_preamble_sic_core_perf.cc_1
    uwb_qa_uwb_preamble_sic_core_perf.cc_2
    uwb_qa_uwb_demod_core_timing_cfo_sfd.cc
    uwb_qa_uwb_demod_core_cir_pipeline.cc
    uwb_qa_uwb_cfo_ml.cc
    uwb_qa_uwb_demod_core_kernel_grid.cc
    uwb_qa_uwb_demod_core_fec_joint_ml.cc
    uwb_qa_uwb_hrp_mod_core.cc
    uwb_qa_uwb_hrp_mod_rx_loopback.cc
    uwb_qa_uwb_rx_fused_cmf.cc
    uwb_qa_uwb_rx_fused_phr.cc
    uwb_qa_uwb_rx_phase_sfo.cc
    uwb_qa_uwb_cfo_derotate.cc
    uwb_qa_uwb_target_lock_core.cc
    uwb_qa_uwb_timed_rx_plan.cc
    uwb_qa_uwb_autocorr_gate.cc
    uwb_qa_uwb_rational_resampler_core.cc
    uwb_qa_uwb_resampler_transposed.cc
    uwb_qa_uwb_timing_refactor.cc
)
    if(TEST ${_t})
        set_tests_properties(${_t} PROPERTIES LABELS "core;rx_sic")
    endif()
endforeach()

# Heavy low-SNR grid scans: compile/run in full and explicit Sweep runs, but
# stay out of the `rx_sic` fast loop.  `ctest -L heavy` opts in.
if(TEST uwb_qa_uwb_preamble_sic_core_heavy.cc)
    set_tests_properties(uwb_qa_uwb_preamble_sic_core_heavy.cc PROPERTIES
        LABELS "heavy")
endif()
foreach(_t IN ITEMS
    uwb_qa_uwb_twr_config.cc
    uwb_qa_uwb_twr_config_parity.cc
    uwb_qa_uwb_twr_frame.cc
    uwb_qa_uwb_twr_math.cc
    uwb_qa_uwb_twr_native_phy.cc
    uwb_qa_uwb_twr_phy_impairment.cc
    uwb_qa_uwb_twr_tof_input.cc
    uwb_qa_uwb_twr_protocol_time.cc
    uwb_qa_uwb_twr_fake_link.cc
    uwb_qa_uwb_twr_core.cc
    uwb_qa_uwb_twr_m1b_e2e.cc
    uwb_qa_twr_config_py
    uwb_qa_twr_fake_demo_py
    uwb_qa_twr_m1b_verify_py
    uwb_qa_twr_m2a_native_phy_verify_py
)
    if(TEST ${_t})
        set_tests_properties(${_t} PROPERTIES LABELS "twr")
    endif()
endforeach()

foreach(_t IN ITEMS
    uwb_qa_uwb_radar_cir_core.cc
    uwb_qa_uwb_radar_cir_estimator.cc
    uwb_qa_uwb_radar_sfd_core.cc
    uwb_qa_uwb_radar_timing_core.cc
)
    if(TEST ${_t})
        set_tests_properties(${_t} PROPERTIES LABELS "sensing")
    endif()
endforeach()

foreach(_t IN ITEMS
    uwb_qa_uwb_rational_resampler_core_65_32.cc
    uwb_qa_uwb_rational_resampler_m2a.cc
    uwb_qa_uwb_resampler_avx512.cc
    uwb_qa_uwb_resampler_burst_opt.cc
    uwb_qa_uwb_fixed_tone_cancel_core.cc
    uwb_qa_uwb_fixed_tone_cancel_source.cc
    uwb_qa_uwb_tone_canceller.cc
)
    if(TEST ${_t})
        set_tests_properties(${_t} PROPERTIES LABELS "dsp_extra")
    endif()
endforeach()

foreach(_t IN ITEMS
    uwb_qa_uwb_fullpacket_cfo_diag.cc
    uwb_qa_uwb_fullpacket_fast_rx.cc
    uwb_qa_uwb_fullpacket_payload_cache.cc
    uwb_qa_fullpacket_segment_diag
    uwb_qa_uwb_demod_sic_avx512.cc
)
    if(TEST ${_t})
        set_tests_properties(${_t} PROPERTIES LABELS "legacy_diag")
    endif()
endforeach()

foreach(_t IN ITEMS
    uwb_qa_uwb_defaults.cc
    uwb_qa_uwb_refactor_contracts.cc
    uwb_qa_uwb_log.cc
    uwb_qa_hrp_standard_rev.cc
    uwb_qa_openuwb_runtime_config.cc
    uwb_qa_uwb_scheduled_dump_io.cc
    uwb_benchmark_hrp_mod
    uwb_qa_native_install_consumer
)
    if(TEST ${_t})
        set_tests_properties(${_t} PROPERTIES LABELS "infra")
    endif()
endforeach()

if(TEST uwb_qa_uwb_uhd_backend_device.cc)
    set_tests_properties(uwb_qa_uwb_uhd_backend_device.cc PROPERTIES LABELS "radio_rx")
endif()
