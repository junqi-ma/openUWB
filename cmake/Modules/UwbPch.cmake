# UwbPch.cmake -- lightweight precompiled headers for lib/ and apps/.
#
# Plan: docs/performance/多Subagent并行开发方案_按功能裁剪与编译测试再加速_2026-11-21.md (B2)
#
# The QA tree has its own PCH prototypes (tests/CMakeLists.txt) because it must
# match Boost.Test macros.  Library and application translation units mostly
# share the stable standard-library surface, so one small header list per
# target is enough and keeps the PCH valid across the per-target ISA flags.
#
# Every target gets its own PCH (target_precompile_headers PRIVATE), so the
# compile options/definitions always match the consuming target.  GCC rejects
# a PCH built with different -m/-D flags, hence nothing is shared globally.

set(UWB_STL_PCH_HEADERS
    <algorithm>
    <array>
    <atomic>
    <chrono>
    <cmath>
    <complex>
    <condition_variable>
    <cstddef>
    <cstdint>
    <cstdio>
    <cstdlib>
    <cstring>
    <deque>
    <fstream>
    <functional>
    <iostream>
    <limits>
    <map>
    <memory>
    <mutex>
    <numeric>
    <optional>
    <random>
    <set>
    <sstream>
    <string>
    <string_view>
    <thread>
    <tuple>
    <type_traits>
    <unordered_map>
    <utility>
    <vector>)

# Skip INTERFACE/ALIAS targets, tiny/empty stub targets (where generating a
# 3-second private .gch exceeds total compile time), and reuse one shared
# prototype (`uwb_lib_pch_simd`) across all PIC+SIMD library targets so only one
# library .gch is compiled instead of one per target.
function(uwb_apply_stl_pch target)
    if(NOT UWB_ENABLE_PCH OR NOT TARGET ${target})
        return()
    endif()
    get_target_property(_aliased ${target} ALIASED_TARGET)
    if(_aliased)
        return()
    endif()
    get_target_property(_type ${target} TYPE)
    if(_type STREQUAL "INTERFACE_LIBRARY")
        return()
    endif()

    get_target_property(_srcs ${target} SOURCES)
    set(_tu_count 0)
    if(_srcs)
        foreach(_src IN LISTS _srcs)
            if(_src MATCHES "\\.(c|cc|cpp|cxx)$" AND
               NOT _src MATCHES "openuwb_compat_empty\\.cc$")
                math(EXPR _tu_count "${_tu_count} + 1")
            endif()
        endforeach()
    endif()

    get_target_property(_is_pic ${target} POSITION_INDEPENDENT_CODE)
    if(_type MATCHES "LIBRARY" AND _is_pic AND TARGET uwb_lib_pch_simd)
        # Reusing an already-built shared .gch is free even for 1-TU libs
        # (e.g., UWB_DEV_SCOPE=tx openuwb_phy_hrp).
        if(_tu_count LESS 1)
            return()
        endif()
        target_precompile_headers(${target} REUSE_FROM uwb_lib_pch_simd)
        return()
    endif()

    # Private PCH only pays off when a target compiles at least 3 TUs (e.g.,
    # openuwb_rx with 4 TUs); single-TU apps lose ~2.5 s generating a .gch.
    if(_tu_count LESS 3)
        return()
    endif()
    target_precompile_headers(${target} PRIVATE ${UWB_STL_PCH_HEADERS})
endfunction()
