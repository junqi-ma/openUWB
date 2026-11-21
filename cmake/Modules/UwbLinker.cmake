# UwbLinker.cmake -- optional fast-linker selection.
#
# Plan: docs/performance/多Subagent并行开发方案_按功能裁剪与编译测试再加速_2026-11-21.md (B2)
#
# UWB_LINKER=AUTO|lld|mold|gold|ld selects the linker used by the C++ driver.
# AUTO probes in preference order (mold, lld, gold, ld) and falls back to the
# toolchain default, so a host without any of the optional linkers still
# builds.  This only adds -fuse-ld=<name>; it never changes code generation
# (unlike -march=native, which the project forbids).

set(UWB_LINKER "AUTO" CACHE STRING
    "Linker driver: AUTO|lld|mold|gold|ld")
set_property(CACHE UWB_LINKER PROPERTY STRINGS AUTO lld mold gold ld)
string(TOLOWER "${UWB_LINKER}" UWB_LINKER_REQUESTED)

if(NOT UWB_LINKER_REQUESTED MATCHES "^(auto|lld|mold|gold|ld)$")
    message(FATAL_ERROR
        "UWB_LINKER must be AUTO, lld, mold, gold or ld "
        "(got '${UWB_LINKER}')")
endif()

# Probe -fuse-ld=<name> by both checking the binary runs cleanly AND linking a
# minimal C++ program.  check_cxx_compiler_flag() does not run the linker on all
# generators, so a dangling wrapper script on PATH could otherwise look valid at
# configure time and fail later during ninja.
function(_uwb_linker_available name out_var)
    set(_bin "")
    if(name STREQUAL "mold")
        find_program(_uwb_ld_mold NAMES mold ld.mold)
        set(_bin "${_uwb_ld_mold}")
    elseif(name STREQUAL "lld")
        find_program(_uwb_ld_lld NAMES ld.lld lld)
        set(_bin "${_uwb_ld_lld}")
    elseif(name STREQUAL "gold")
        find_program(_uwb_ld_gold NAMES ld.gold)
        set(_bin "${_uwb_ld_gold}")
    elseif(name STREQUAL "ld")
        set(${out_var} TRUE PARENT_SCOPE)
        return()
    endif()

    if(NOT _bin OR NOT EXISTS "${_bin}")
        set(${out_var} FALSE PARENT_SCOPE)
        return()
    endif()
    execute_process(
        COMMAND "${_bin}" --version
        RESULT_VARIABLE _uwb_ld_rc
        OUTPUT_QUIET
        ERROR_QUIET)
    if(NOT _uwb_ld_rc EQUAL 0)
        set(${out_var} FALSE PARENT_SCOPE)
        return()
    endif()

    include(CheckCXXSourceCompiles)
    include(CMakePushCheckState)
    string(TOUPPER "${name}" _cache_suffix)
    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_LINK_OPTIONS "-fuse-ld=${name}")
    check_cxx_source_compiles("int main() { return 0; }"
        UWB_CC_LINKS_${_cache_suffix})
    cmake_pop_check_state()
    set(${out_var} ${UWB_CC_LINKS_${_cache_suffix}} PARENT_SCOPE)
endfunction()

set(_uwb_linker_selected "")
if(UWB_LINKER_REQUESTED STREQUAL "auto")
    foreach(_candidate IN ITEMS mold lld gold)
        _uwb_linker_available(${_candidate} _ok)
        if(_ok)
            set(_uwb_linker_selected ${_candidate})
            break()
        endif()
    endforeach()
elseif(NOT UWB_LINKER_REQUESTED STREQUAL "ld")
    _uwb_linker_available(${UWB_LINKER_REQUESTED} _ok)
    if(_ok)
        set(_uwb_linker_selected ${UWB_LINKER_REQUESTED})
    else()
        message(WARNING
            "UWB_LINKER=${UWB_LINKER_REQUESTED} is not usable by this "
            "compiler; falling back to the toolchain default linker")
    endif()
endif()

if(_uwb_linker_selected)
    add_link_options("-fuse-ld=${_uwb_linker_selected}")
    message(STATUS "UWB: using ${_uwb_linker_selected} linker (-fuse-ld)")
else()
    message(STATUS "UWB: using the toolchain default linker")
endif()
