# Copyright 2011-2020 Free Software Foundation, Inc.
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.20)
# Phase 3: the repository root is the project.  `OpenUWB` is the package name
# consumers find; targets keep the `openuwb_` prefix (see lib/CMakeLists.txt).
project(OpenUWB LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Catch an old build command that explicitly asks for the retired adapter.
if(UWB_BUILD_GNURADIO)
    message(FATAL_ERROR "The GNU Radio adapter is archived on archive/gnuradio-adapter. Configure this branch in a fresh native build directory.")
endif()

include(GNUInstallDirs)
# Runtime resources (taps, preamble template) live under share/openuwb; the
# lowercase directory is deliberate, `share/OpenUWB` would look like a casing
# accident on case-sensitive filesystems.
set(UWB_DATA_DIR ${CMAKE_INSTALL_DATADIR}/openuwb)
# Primary package: find_package(OpenUWB CONFIG REQUIRED) + openuwb::* targets.
set(OPENUWB_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/OpenUWB)
# Legacy package kept for existing consumers.
# TODO(openuwb-compat-removal): drop the uwb package once downstream migrated.
set(UWB_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/uwb)

find_package(Volk REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Python3 COMPONENTS Interpreter QUIET)

# Select the release build type by default to get optimization flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "")

option(UWB_BUILD_EXPERIMENTS
    "Build historical benchmark/diagnostic experiment apps" ON)

# Stage 1 (2026-11-21): build-base upgrades.  Generator selection moves to
# CMakePresets.json (Ninja); the options below cover the two general lifts:
# precompiled headers for the QA targets and ccache acceleration.
option(UWB_ENABLE_PCH
    "Enable precompiled headers for the QA/library/app targets" ON)
option(UWB_USE_CCACHE
    "Use ccache as the C++ compiler launcher when available" ON)

if(UWB_USE_CCACHE)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        execute_process(
            COMMAND "${CCACHE_PROGRAM}" --version
            RESULT_VARIABLE _uwb_ccache_rc
            OUTPUT_QUIET
            ERROR_QUIET)
        if(_uwb_ccache_rc EQUAL 0)
            set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
            message(STATUS "UWB: using ccache (${CCACHE_PROGRAM}) as compiler launcher")
        else()
            message(WARNING
                "UWB: found ccache (${CCACHE_PROGRAM}) but '--version' exited "
                "with ${_uwb_ccache_rc}; compiling without ccache")
        endif()
    else()
        message(STATUS "UWB: ccache not found; compiling without it")
    endif()
endif()

# ccache can cache PCH-using TUs with GCC only when the driver emits a
# preprocessed form it understands.  Without this flag every QA/lib/app TU that
# consumes the CMake PCH is reported as "uncacheable".
if(CMAKE_CXX_COMPILER_LAUNCHER AND UWB_ENABLE_PCH AND
   CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fpch-preprocess)
    message(STATUS
        "UWB: ccache + PCH: adding -fpch-preprocess to keep PCH TUs cacheable")
endif()

########################################################################
# Stage 2 (2026-11-21): functional-scope builds.
#
# UWB_DEV_SCOPE selects which domain libraries, apps and QA targets are
# configured.  `all` (the default and the delivery/CI setting) is byte-for-byte
# the previous behaviour.  The four single-domain values are for the daily
# edit/build/test loop; they force the finer UWB_ENABLE_* switches to the
# mapping documented in
# docs/performance/多Subagent并行开发方案_按功能裁剪与编译测试再加速_2026-11-21.md.
#
# Scoped builds never register install() rules and never register the
# install_consumer test: a partial tree must not be installed over a full one.
########################################################################
set(UWB_DEV_SCOPE "all" CACHE STRING
    "Functional scope: all|rx|tx|twr|radar (single-domain scopes are dev only)")
set_property(CACHE UWB_DEV_SCOPE PROPERTY STRINGS all rx tx twr radar)
string(TOLOWER "${UWB_DEV_SCOPE}" UWB_DEV_SCOPE)

option(UWB_ENABLE_TX "Configure the transmit/PHY-modulation domain" ON)
option(UWB_ENABLE_RX "Configure the receive/demodulation domain" ON)
option(UWB_ENABLE_SIC "Configure the successive-interference-cancellation domain" ON)
option(UWB_ENABLE_RADAR "Configure the monostatic-radar/sensing domain" ON)
option(UWB_ENABLE_TWR "Configure the two-way-ranging domain" ON)
option(UWB_ENABLE_BENCHMARKS "Configure the benchmark/measurement executables" ON)
# Development-only: build the domain archives as shared libraries so editing a
# library TU relinks one .so instead of every QA executable.  The delivery
# presets keep the static archives (default OFF).
option(UWB_DEV_SHARED_LIBS
    "Build openuwb_* domain libraries as shared objects (dev only)" OFF)
set(UWB_LIBRARY_TYPE STATIC)
if(UWB_DEV_SHARED_LIBS)
    set(UWB_LIBRARY_TYPE SHARED)
    # A changed .so does not require relinking every executable: symbols are
    # resolved at run time.  Without this CMake/Ninja would still relink all
    # consumers because the shared object is a link input.
    set(CMAKE_LINK_DEPENDS_NO_SHARED ON)
endif()

if(NOT UWB_DEV_SCOPE MATCHES "^(all|rx|tx|twr|radar)$")
    message(FATAL_ERROR
        "UWB_DEV_SCOPE must be all, rx, tx, twr or radar "
        "(got '${UWB_DEV_SCOPE}')")
endif()

if(NOT UWB_DEV_SCOPE STREQUAL "all")
    # Single-domain scope: lock the fine-grained switches to the table in the
    # plan so the user cannot accidentally configure a half-consistent tree.
    set(UWB_ENABLE_BENCHMARKS OFF)
    if(UWB_DEV_SCOPE STREQUAL "rx")
        set(UWB_ENABLE_TX OFF)
        set(UWB_ENABLE_RX ON)
        set(UWB_ENABLE_SIC ON)
        set(UWB_ENABLE_RADAR OFF)
        set(UWB_ENABLE_TWR OFF)
    elseif(UWB_DEV_SCOPE STREQUAL "tx")
        set(UWB_ENABLE_TX ON)
        set(UWB_ENABLE_RX OFF)
        set(UWB_ENABLE_SIC OFF)
        set(UWB_ENABLE_RADAR OFF)
        set(UWB_ENABLE_TWR OFF)
    elseif(UWB_DEV_SCOPE STREQUAL "radar")
        set(UWB_ENABLE_TX OFF)
        set(UWB_ENABLE_RX ON)
        set(UWB_ENABLE_SIC ON)
        set(UWB_ENABLE_RADAR ON)
        set(UWB_ENABLE_TWR OFF)
    elseif(UWB_DEV_SCOPE STREQUAL "twr")
        set(UWB_ENABLE_TX ON)
        set(UWB_ENABLE_RX ON)
        set(UWB_ENABLE_SIC OFF)
        set(UWB_ENABLE_RADAR OFF)
        set(UWB_ENABLE_TWR ON)
    endif()
endif()

# Internal consistency: the radar chain needs SIC for the scheduled pipeline.
# RX alone (the TWR scope builds PHY receive cores but no runtime) does not.
if(UWB_ENABLE_RADAR AND NOT (UWB_ENABLE_RX AND UWB_ENABLE_SIC))
    message(WARNING
        "UWB_ENABLE_RADAR=ON requires RX + SIC; enabling both.")
    set(UWB_ENABLE_RX ON)
    set(UWB_ENABLE_SIC ON)
endif()
# The runtime archive needs SIC; radio exists to feed it.  TWR's RX bit only
# selects the PHY receive cores, so no runtime/radio is pulled into that tree.
set(UWB_ENABLE_RUNTIME OFF)
if(UWB_ENABLE_RX AND UWB_ENABLE_SIC)
    set(UWB_ENABLE_RUNTIME ON)
endif()
set(UWB_ENABLE_RADIO ${UWB_ENABLE_RUNTIME})
# The receive applications need the runtime + radio chain.
set(UWB_ENABLE_RX_APPS OFF)
if(UWB_ENABLE_RUNTIME AND
   (UWB_DEV_SCOPE STREQUAL "all" OR UWB_DEV_SCOPE STREQUAL "rx"))
    set(UWB_ENABLE_RX_APPS ON)
endif()

# A partial tree must never produce install rules or install_consumer.
if(NOT UWB_DEV_SCOPE STREQUAL "all")
    set(CMAKE_SKIP_INSTALL_RULES TRUE)
    message(STATUS
        "UWB_DEV_SCOPE=${UWB_DEV_SCOPE}: install() rules and the "
        "install_consumer test are disabled for this development build")
endif()

message(STATUS
    "UWB scope: ${UWB_DEV_SCOPE} "
    "(TX=${UWB_ENABLE_TX} RX=${UWB_ENABLE_RX} SIC=${UWB_ENABLE_SIC} "
    "RUNTIME=${UWB_ENABLE_RUNTIME} RADAR=${UWB_ENABLE_RADAR} "
    "TWR=${UWB_ENABLE_TWR} BENCHMARKS=${UWB_ENABLE_BENCHMARKS})")

# Make sure our local CMake Modules path comes first
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Stage 2: linker selection (AUTO|lld|mold|gold|ld).  AUTO prefers the fastest
# installed linker (mold > lld > gold > ld); the probe keeps the build working
# on hosts that only ship GNU ld.  Only --fuse-ld is added, never
# -march=native or any other code-generation flag.
include(UwbLinker)

# Stage 2: a small stable STL PCH for the library/app translation units.  The
# QA PCH prototypes live in tests/CMakeLists.txt; this module covers lib/ and
# apps/ where the previous Stage 1 work had no PCH at all.
include(UwbPch)

########################################################################
# Compile-time ISA selection for the self-written DSP kernels
# (UWB_SIMD=AUTO|AVX2|AVX512|SCALAR).  See cmake/Modules/UwbSimd.cmake and
# docs/performance/AVX512内核与编译期CPU选择开发规划_2026-11-16.md.
########################################################################
include(UwbSimd)
set(_uwb_sparse_mf_avx512_default OFF)
if(UWB_SIMD_HAVE_AVX512)
    set(_uwb_sparse_mf_avx512_default ON)
endif()
option(UWB_SPARSE_MF_AVX512_16LAG
    "Use AVX-512 16-lag batching for sparse ternary matched filtering when AVX-512 is enabled"
    ${_uwb_sparse_mf_avx512_default})
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/openuwb/support/simd_config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/openuwb/support/simd_config.h
    @ONLY)

########################################################################
# openuwb::timing (support layer): the wall-clock measurement primitives.
#
# The instrumentation RAII scopes default to ON, which is the historical,
# always-on behaviour of the stage-boundary telemetry.  OFF compiles only
# ScopedNs<false> by default (a TU can still opt in with ScopedNs<true>);
# the plain elapsed_ns/elapsed_us/thread_cpu helpers are never compiled out
# because the stage-boundary contracts quote their values verbatim.
########################################################################
option(UWB_ENABLE_TIMING
    "Compile the openuwb::timing RAII scopes enabled by default" ON)
set(UWB_TIMING_ENABLED_DEFAULT 1)
if(NOT UWB_ENABLE_TIMING)
    set(UWB_TIMING_ENABLED_DEFAULT 0)
endif()
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/openuwb/support/timing_config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/openuwb/support/timing_config.h
    @ONLY)

# The generated header must be visible to every subdirectory (and to the
# installed consumer, which gets the same file from the install rule).
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

if(APPLE)
    if(NOT CMAKE_INSTALL_NAME_DIR)
        set(CMAKE_INSTALL_NAME_DIR
            ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} CACHE
            PATH "Library Install Name Destination Directory" FORCE)
    endif(NOT CMAKE_INSTALL_NAME_DIR)
    if(NOT CMAKE_INSTALL_RPATH)
        set(CMAKE_INSTALL_RPATH
            ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} CACHE
            PATH "Library Install RPath" FORCE)
    endif(NOT CMAKE_INSTALL_RPATH)
    if(NOT CMAKE_BUILD_WITH_INSTALL_RPATH)
        set(CMAKE_BUILD_WITH_INSTALL_RPATH ON CACHE
            BOOL "Do Build Using Library Install RPath" FORCE)
    endif(NOT CMAKE_BUILD_WITH_INSTALL_RPATH)
endif(APPLE)

########################################################################
# Optional UHD backend discovery (AUTO/ON/OFF
# semantics).  Sets UWB_HAVE_UHD and makes the UHD::UHD imported target.
########################################################################
set(ENABLE_UHD_BACKEND "AUTO" CACHE STRING
    "Build the UHD device backends: AUTO (build when UHD found), ON (require UHD), OFF (never)")
set_property(CACHE ENABLE_UHD_BACKEND PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${ENABLE_UHD_BACKEND}" ENABLE_UHD_BACKEND_UC)
if(NOT ENABLE_UHD_BACKEND_UC MATCHES "^(AUTO|ON|OFF)$")
    message(FATAL_ERROR
        "ENABLE_UHD_BACKEND must be AUTO, ON or OFF "
        "(got '${ENABLE_UHD_BACKEND}')")
endif()

set(UWB_HAVE_UHD 0)
if(NOT ENABLE_UHD_BACKEND_UC STREQUAL "OFF")
    find_package(UHD CONFIG QUIET)
    if(NOT UHD_FOUND)
        find_package(PkgConfig QUIET)
        if(PkgConfig_FOUND)
            pkg_check_modules(UHD_PC QUIET uhd)
        endif()
    endif()
    if(UHD_FOUND OR UHD_PC_FOUND)
        set(UWB_HAVE_UHD 1)
    elseif(ENABLE_UHD_BACKEND_UC STREQUAL "ON")
        message(FATAL_ERROR
            "ENABLE_UHD_BACKEND=ON but UHD was not found "
            "(install libuhd-dev / UHD CMake or pkg-config support). "
            "Use ENABLE_UHD_BACKEND=AUTO to skip the backend instead.")
    else()
        message(STATUS
            "UHD not found: UHD device backends disabled; the native replay "
            "chain, library and QA build and test without UHD.")
    endif()
endif()
if(UWB_HAVE_UHD AND NOT TARGET UHD::UHD)
    add_library(UHD::UHD INTERFACE IMPORTED GLOBAL)
    if(UHD_INCLUDE_DIRS)
        set_target_properties(UHD::UHD PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${UHD_INCLUDE_DIRS}")
    endif()
    set_target_properties(UHD::UHD PROPERTIES
        INTERFACE_LINK_LIBRARIES "${UHD_LIBRARIES};${UHD_PC_LIBRARIES}")
endif()
message(STATUS "UHD backends (ENABLE_UHD_BACKEND=${ENABLE_UHD_BACKEND_UC}): ${UWB_HAVE_UHD}")

if(UWB_HAVE_UHD)
    # The exported openuwb::uwb_native references UHD::UHD; the installed config
    # must define that target even when the consumer's UHD package does not
    # (some UHD packages only set UHD_LIBRARIES/UHD_INCLUDE_DIRS or pkg-config).
    # \${...} escapes keep the consumer-time ${CMAKE_CURRENT_LIST_DIR} literal.
    set(OPENUWB_CONFIG_UHD_DEP
        "include(\"\${CMAKE_CURRENT_LIST_DIR}/OpenUWBConfigUhdDep.cmake\")")
    # The legacy package lives in a sibling directory (cmake/uwb vs
    # cmake/OpenUWB), so it must reach across instead of looking beside itself.
    set(UWB_CONFIG_UHD_DEP
        "include(\"\${CMAKE_CURRENT_LIST_DIR}/../OpenUWB/OpenUWBConfigUhdDep.cmake\")")
else()
    set(OPENUWB_CONFIG_UHD_DEP "# this build has no UHD dependency")
    set(UWB_CONFIG_UHD_DEP "# this build has no UHD dependency")
endif()


configure_file(cmake/cmake_uninstall.cmake.in cmake_uninstall.cmake @ONLY)
add_custom_target(uninstall
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

add_subdirectory(include/gnuradio/uwb)
add_subdirectory(lib)
add_subdirectory(apps)
# Phase 1: the unit QA, benchmarks and integration harness live in tests/.
# Registered last so it can reference the app targets it drives
# (twr_fake_demo, twr_native_phy_demo, uwb_uhd_sc16_pipeline_runner).
add_subdirectory(tests)

include(CMakePackageConfigHelpers)
# Primary package: find_package(OpenUWB CONFIG REQUIRED) + openuwb::* targets.
configure_package_config_file(
    cmake/Modules/OpenUWBConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/OpenUWBConfig.cmake
    INSTALL_DESTINATION ${OPENUWB_CMAKE_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/OpenUWBConfig.cmake
    cmake/Modules/OpenUWBConfigUhdDep.cmake
    DESTINATION ${OPENUWB_CMAKE_DIR})

# Legacy package name.  Same exported targets, with the old uwb:: spellings
# re-created on top.
# TODO(openuwb-compat-removal): remove the uwb package and its config.
configure_package_config_file(
    cmake/Modules/uwbConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/uwbConfig.cmake
    INSTALL_DESTINATION ${UWB_CMAKE_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/uwbConfig.cmake
    DESTINATION ${UWB_CMAKE_DIR})

# Install the runtime resources the native live receiver needs (resampler
# taps + preamble template), so an INSTALLED uwb_dw3000_live does not depend
# on a source checkout, an old build tree, RPATH or PYTHONPATH.
########################################################################
install(FILES
    ${CMAKE_SOURCE_DIR}/testdata/resampler_65_48/taps_quality_minorder.txt
    ${CMAKE_SOURCE_DIR}/testdata/resampler_65_48/taps_quality_eq32.txt
    ${CMAKE_SOURCE_DIR}/testdata/resampler_65_48/taps_realtime_minorder.txt
    ${CMAKE_SOURCE_DIR}/testdata/resampler_65_48/taps_realtime_eq16.txt
    DESTINATION ${UWB_DATA_DIR}/resampler_65_48
)
install(FILES
    ${CMAKE_SOURCE_DIR}/testdata/reference_preamble_code10_998p4.cf32
    DESTINATION ${UWB_DATA_DIR}
)
