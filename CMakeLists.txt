cmake_minimum_required(VERSION 3.16)
project(ikafssn LANGUAGES C CXX)

set(IKAFSSN_VERSION "0.1.2026.05.05")  # ikafssnindex: parallelise per-partition combine_each merge and EF dictionary encode

# --- Optional benchmark suite (Google Benchmark via FetchContent) ---
# OFF by default; ON pulls Google Benchmark v1.9.4 and adds bench/ subdir.
option(IKAFSSN_BUILD_BENCH
       "Build benchmark binaries under bench/ (off by default)"
       OFF)

# --- Optional analysis tools (Phase 5a) ---
# OFF by default; ON adds tools/ subdir (throwaway dump utilities).
option(IKAFSSN_BUILD_TOOLS
       "Build analysis tools under tools/ (Phase 5a feasibility study)"
       OFF)

# Get timezone offset in +HH:MM format for build timestamp
string(TIMESTAMP _tz "%z")
string(REGEX REPLACE "^([+-][0-9][0-9])([0-9][0-9])$" "\\1:\\2" IKAFSSN_BUILD_TZ_OFFSET "${_tz}")

configure_file(
    "${CMAKE_SOURCE_DIR}/src/core/version.hpp.in"
    "${CMAKE_BINARY_DIR}/generated/core/version.hpp"
    @ONLY
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release build: -O2 with loop unrolling (safer than default -O3,
# and ensures NTTP extraction loops are fully unrolled with immediates).
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -funroll-loops" CACHE STRING
    "C++ flags for Release builds" FORCE)

# --- Runtime SIMD dispatch ---
# Enables CPU-feature detection and ISA-specific kernel selection at runtime.
# Disabling forces the Scalar tier even on capable hardware.
option(IKAFSSN_ENABLE_SIMD
       "Enable runtime SIMD dispatch (CPU detection + ISA-specific kernels)"
       ON)
if(IKAFSSN_ENABLE_SIMD)
    add_compile_definitions(IKAFSSN_ENABLE_SIMD=1)
else()
    add_compile_definitions(IKAFSSN_ENABLE_SIMD=0)
endif()

# --- NCBI C++ Toolkit ---
set(NCBI_TOOLKIT_DIR "${CMAKE_SOURCE_DIR}/ncbi-cxx-toolkit" CACHE PATH
    "Path to NCBI C++ Toolkit install root (contains include/ and CMake-*/)")
set(NCBI_TOOLKIT_BUILD_TAG "CMake-GCC1330-Release" CACHE STRING
    "NCBI Toolkit build subdirectory name (e.g. CMake-GCC1330-Release)")

set(NCBI_TOOLKIT_INCLUDE
    "${NCBI_TOOLKIT_DIR}/include"
    "${NCBI_TOOLKIT_DIR}/${NCBI_TOOLKIT_BUILD_TAG}/inc"
)
set(NCBI_TOOLKIT_CMAKE
    "${NCBI_TOOLKIT_DIR}/${NCBI_TOOLKIT_BUILD_TAG}/cmake/ncbi-cpp-toolkit.cmake")

# Define dependency targets required by NCBI export before including it
if(NOT TARGET LMDB::lmdb_LIB)
    find_library(LMDB_LIBRARY NAMES lmdb REQUIRED)
    add_library(LMDB::lmdb_LIB UNKNOWN IMPORTED)
    set_target_properties(LMDB::lmdb_LIB PROPERTIES
        IMPORTED_LOCATION "${LMDB_LIBRARY}"
    )
endif()

if(NOT TARGET NCBI::DL_LIBS)
    add_library(NCBI::DL_LIBS INTERFACE IMPORTED)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set_target_properties(NCBI::DL_LIBS PROPERTIES
            INTERFACE_LINK_LIBRARIES "dl"
        )
    endif()
endif()

if(NOT TARGET NCBI::RT_LIBS)
    add_library(NCBI::RT_LIBS INTERFACE IMPORTED)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set_target_properties(NCBI::RT_LIBS PROPERTIES
            INTERFACE_LINK_LIBRARIES "rt"
        )
    endif()
endif()

if(NOT TARGET NCBI::MATH_LIBS)
    add_library(NCBI::MATH_LIBS INTERFACE IMPORTED)
    set_target_properties(NCBI::MATH_LIBS PROPERTIES
        INTERFACE_LINK_LIBRARIES "m"
    )
endif()

include("${NCBI_TOOLKIT_CMAKE}")

# --- Parasail (for Stage 3 alignment) ---
set(PARASAIL_DIR "${CMAKE_SOURCE_DIR}/parasail" CACHE PATH
    "Path to Parasail install root")
add_library(parasail STATIC IMPORTED)
set_target_properties(parasail PROPERTIES
    IMPORTED_LOCATION "${PARASAIL_DIR}/lib/libparasail.a"
    INTERFACE_INCLUDE_DIRECTORIES "${PARASAIL_DIR}/include"
)

# --- htslib (for SAM/BAM output) ---
set(HTSLIB_DIR "${CMAKE_SOURCE_DIR}/htslib" CACHE PATH
    "Path to htslib install root")
find_library(DEFLATE_LIB NAMES deflate REQUIRED)
add_library(hts STATIC IMPORTED)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_DIR}/lib/libhts.a"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_DIR}/include"
    INTERFACE_LINK_LIBRARIES "z;bz2;lzma;${DEFLATE_LIB};pthread;m"
)

# --- Intel TBB ---
find_package(TBB REQUIRED)

# --- Drogon (optional, for ikafssnhttpd) ---
option(BUILD_HTTPD "Build ikafssnhttpd (requires Drogon)" ON)
if(BUILD_HTTPD)
    set(DROGON_DIR "${CMAKE_SOURCE_DIR}/drogon" CACHE PATH
        "Path to Drogon install root (contains lib/cmake/Drogon/)")
    # Prepend DROGON_DIR to CMAKE_PREFIX_PATH so that DrogonConfig.cmake's
    # transitive find_dependency(Trantor) also resolves under the same prefix.
    list(PREPEND CMAKE_PREFIX_PATH "${DROGON_DIR}")
    find_package(Drogon REQUIRED)
endif()

# --- OpenSSL (for ikafssnclient checkpointing SHA256) ---
find_package(OpenSSL REQUIRED)

# --- libcurl + jsoncpp (optional, for ikafssnclient HTTP mode) ---
option(ENABLE_CLIENT_HTTP "Enable HTTP mode in ikafssnclient (requires libcurl)" ON)
if(ENABLE_CLIENT_HTTP)
    find_package(CURL REQUIRED)
    find_package(jsoncpp REQUIRED)
    # Determine the correct jsoncpp imported target.
    # - Ubuntu's libjsoncpp-dev provides jsoncpp_lib + JsonCpp::JsonCpp.
    # - Drogon's bundled Findjsoncpp.cmake only sets JSONCPP_INCLUDE_DIRS /
    #   JSONCPP_LIBRARIES variables (no imported target).  Create one manually.
    if(TARGET jsoncpp_lib)
        set(JSONCPP_TARGET jsoncpp_lib)
    elseif(TARGET JsonCpp::JsonCpp)
        set(JSONCPP_TARGET JsonCpp::JsonCpp)
    elseif(JSONCPP_LIBRARIES)
        add_library(_jsoncpp_imported UNKNOWN IMPORTED)
        set_target_properties(_jsoncpp_imported PROPERTIES
            IMPORTED_LOCATION "${JSONCPP_LIBRARIES}")
        if(JSONCPP_INCLUDE_DIRS)
            set_target_properties(_jsoncpp_imported PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${JSONCPP_INCLUDE_DIRS}")
        endif()
        set(JSONCPP_TARGET _jsoncpp_imported)
    else()
        message(FATAL_ERROR "jsoncpp found but no target or variables available")
    endif()
endif()

# --- libcurl (optional, for ikafssnretrieve remote mode) ---
option(ENABLE_REMOTE_RETRIEVE "Enable NCBI efetch in ikafssnretrieve (requires libcurl)" ON)
if(ENABLE_REMOTE_RETRIEVE)
    find_package(CURL REQUIRED)
endif()

# --- x86-simd-sort (Phase 2c, M-4) ---
# Header + per-arch dispatch library used by index_builder for the posting
# buffer sort. Header-only on its own; the per-arch object files (avx2 / skx)
# are compiled with -march=haswell / -march=skylake-avx512 respectively, and
# the dispatcher in lib/x86simdsort.cpp picks at runtime via CPU feature bits.
# Only enabled on x86_64 with IKAFSSN_ENABLE_SIMD; aarch64 / SIMD-disabled
# builds fall back to tbb::parallel_sort.
if(IKAFSSN_ENABLE_SIMD AND CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)")
    include(FetchContent)
    FetchContent_Declare(
        x86_simd_sort
        URL      https://github.com/numpy/x86-simd-sort/archive/refs/tags/v7.0.tar.gz
        URL_HASH SHA256=f1658b0562a56255b157b7ffbb81e53e78e2571db1962262e55db6ca01733907
    )
    FetchContent_MakeAvailable(x86_simd_sort)

    # Per-arch object libraries; runtime dispatch in x86simdsort.cpp picks one.
    # All four arches must be linked in because the dispatcher's __attribute__
    # ((constructor)) symbol-resolution functions reference every variant
    # unconditionally; the runtime CPU check only chooses among them.
    add_library(ikafssn_xss_avx2 OBJECT
        ${x86_simd_sort_SOURCE_DIR}/lib/x86simdsort-avx2.cpp)
    target_include_directories(ikafssn_xss_avx2 PRIVATE
        ${x86_simd_sort_SOURCE_DIR}/src
        ${x86_simd_sort_SOURCE_DIR}/lib)
    target_compile_options(ikafssn_xss_avx2 PRIVATE -march=haswell)
    set_property(TARGET ikafssn_xss_avx2 PROPERTY POSITION_INDEPENDENT_CODE ON)

    add_library(ikafssn_xss_skx OBJECT
        ${x86_simd_sort_SOURCE_DIR}/lib/x86simdsort-skx.cpp)
    target_include_directories(ikafssn_xss_skx PRIVATE
        ${x86_simd_sort_SOURCE_DIR}/src
        ${x86_simd_sort_SOURCE_DIR}/lib)
    target_compile_options(ikafssn_xss_skx PRIVATE -march=skylake-avx512)
    set_property(TARGET ikafssn_xss_skx PROPERTY POSITION_INDEPENDENT_CODE ON)

    add_library(ikafssn_xss_icl OBJECT
        ${x86_simd_sort_SOURCE_DIR}/lib/x86simdsort-icl.cpp)
    target_include_directories(ikafssn_xss_icl PRIVATE
        ${x86_simd_sort_SOURCE_DIR}/src
        ${x86_simd_sort_SOURCE_DIR}/lib)
    target_compile_options(ikafssn_xss_icl PRIVATE -march=icelake-client)
    set_property(TARGET ikafssn_xss_icl PROPERTY POSITION_INDEPENDENT_CODE ON)

    add_library(ikafssn_xss_spr OBJECT
        ${x86_simd_sort_SOURCE_DIR}/lib/x86simdsort-spr.cpp)
    target_include_directories(ikafssn_xss_spr PRIVATE
        ${x86_simd_sort_SOURCE_DIR}/src
        ${x86_simd_sort_SOURCE_DIR}/lib)
    target_compile_options(ikafssn_xss_spr PRIVATE -march=sapphirerapids)
    set_property(TARGET ikafssn_xss_spr PROPERTY POSITION_INDEPENDENT_CODE ON)

    add_library(ikafssn_xss STATIC
        ${x86_simd_sort_SOURCE_DIR}/lib/x86simdsort.cpp
        $<TARGET_OBJECTS:ikafssn_xss_avx2>
        $<TARGET_OBJECTS:ikafssn_xss_skx>
        $<TARGET_OBJECTS:ikafssn_xss_icl>
        $<TARGET_OBJECTS:ikafssn_xss_spr>
    )
    target_include_directories(ikafssn_xss PUBLIC
        ${x86_simd_sort_SOURCE_DIR}/lib)
    target_include_directories(ikafssn_xss PRIVATE
        ${x86_simd_sort_SOURCE_DIR}/src
        ${x86_simd_sort_SOURCE_DIR}/utils)
    set_property(TARGET ikafssn_xss PROPERTY POSITION_INDEPENDENT_CODE ON)

    set(IKAFSSN_HAS_X86_SIMD_SORT TRUE)
    add_compile_definitions(IKAFSSN_HAS_X86_SIMD_SORT=1)
else()
    set(IKAFSSN_HAS_X86_SIMD_SORT FALSE)
    add_compile_definitions(IKAFSSN_HAS_X86_SIMD_SORT=0)
endif()

# --- FastPFor (Phase 5d/5f: per-tier object libraries for runtime dispatch) ---
#
# FastPFor 0.4.0 ships with FASTPFOR_SIMD_MODE = "native" by default, which
# bakes the build host's -march into the static archive — incompatible with
# ikafssn's runtime SIMD dispatch policy.  Worse, FastPFor's own source
# files contain no explicit AVX-512 intrinsics (only SSE2-SSE4.2 and a
# handful of AVX2 in src/varintdecode.c); whatever AVX-512 instructions
# end up in libFastPFOR.a come purely from the compiler's auto-
# vectorization of surrounding scalar/SSE code.
#
# We side-step the upstream build entirely by populating the source tree
# (no add_subdirectory) and compiling the FastPFor TUs ourselves into 4
# OBJECT libraries, one per ISA tier (Phase 5f x86_64 4-tier ladder):
#
#   - ikafssn_pfd_sse42         -msse4.2 -mpopcnt
#   - ikafssn_pfd_avx2          -mavx2 -mbmi -mbmi2 -mfma
#   - ikafssn_pfd_avx512bw      -mavx512bw -mavx512vl + AVX2 deps
#   - ikafssn_pfd_avx512vbmi2   -mavx512vbmi2 -mavx512vbmi + above
#
# The sse42 tier is semantically equivalent to FASTPFOR_SIMD_MODE=portable
# (FastPFor's SSE4.2 baseline). The other three tiers are equivalent to
# FASTPFOR_SIMD_MODE=native with the upstream `-march=native` replaced by
# the explicit ISA flag set listed above. Because we drive the per-tier
# compile from `target_compile_options`, no `-march=native` ever reaches
# `flags.make`, so the resulting object files run on any CPU that supports
# the named tier — e.g. a GitHub Actions runner without AVX-512 VBMI2 can
# still build the avx512vbmi2 OBJECT library because we only set the ISA
# flags rather than asking the compiler to introspect the build host.
#
# AVX-512 VBMI/VBMI2 are merged into a single top tier — FastPFor's actual
# source contains no byte-level shuffle/permute that VBMI/VBMI2 would
# unlock, so the gap between the two tiers is purely auto-vectorization
# and doesn't justify a separate object library.
#
# Symbol collision (every tier compiles the same FastPForLib::* symbols)
# is avoided by remapping the namespace at compile time via
# `-DFastPForLib=FastPForLib_<tier>`.  ikafssn_pfd_dispatch.cpp has
# forward declarations that name the per-tier wrapper namespaces and
# selects the right entry points from current_simd_cap() at first use.
include(FetchContent)
FetchContent_Declare(
    fastpfor
    URL      https://github.com/lemire/FastPFor/archive/refs/tags/v0.4.0.tar.gz
    URL_HASH SHA256=94c15a3e182c4f4ac0551cdd4c1afdb8465d337293a637c12483d0da026d75af
)
FetchContent_GetProperties(fastpfor)
if(NOT fastpfor_POPULATED)
    # Populate sources but do not run FastPFor's own CMakeLists (which
    # would build a "native" libFastPFOR archive we can't use).
    FetchContent_Populate(fastpfor)
endif()

set(FASTPFOR_SRC_DIR "${fastpfor_SOURCE_DIR}")
set(FASTPFOR_INC_DIR "${FASTPFOR_SRC_DIR}/headers")

# --- SIMDe (Phase 5h: AArch64 NEON support for FastPFor) ---
#
# FastPFor's headers/common.h has a native dispatch path for aarch64 that
# requires <simde/x86/sse4.1.h>; SIMDe (https://github.com/simd-everywhere/simde)
# is the upstream-blessed SSE/AVX → NEON emulation layer.  We populate the
# SIMDe sources only on aarch64 builds and feed them into the per-tier
# OBJECT library defined below.  The exact commit is the one FastPFor
# 0.4.0's cmake_modules/simde.cmake pins (kept verbatim so any future
# FastPFor bump can simply re-sync).
#
# SVE / SVE2 capable CPUs are NOT given a separate FastPFor tier object:
# SIMDe translates SSE intrinsics into NEON only, and a SVE-native bitpacker
# would be a hand-coded effort orthogonal to the SIMDe path.  At runtime
# pfd_codec.cpp's dispatcher routes any aarch64 SimdCap >= NEON (which
# includes SVE / SVE2) to the single ikafssn_pfd_neon tier object.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    FetchContent_Declare(
        simde
        URL      https://github.com/simd-everywhere/simde/archive/c6ddddc4a5bee9913b60de6757227aa078192663.tar.gz
        URL_HASH SHA256=8e83d4812476598a77a4a3c633c37f277040a6bbf08fd1a9667f4655b4deff01
    )
    FetchContent_GetProperties(simde)
    if(NOT simde_POPULATED)
        FetchContent_Populate(simde)
    endif()
    set(SIMDE_INC_DIR "${simde_SOURCE_DIR}")
endif()

# FastPFor TUs we compile per tier.  Phase 5e trimmed this list down to
# only the bitpacker primitives — codecfactory.cpp / streamvbyte.c /
# varintdecode.c / horizontalbitpacking.cpp are no longer used now that
# pfd_codec_tier.cpp drives `simdpackwithoutmask` / `simdunpack` directly
# under a custom block layout (FOR-within-block + varint tail).
set(IKAFSSN_FASTPFOR_SOURCES
    ${FASTPFOR_SRC_DIR}/src/bitpacking.cpp
    ${FASTPFOR_SRC_DIR}/src/bitpackingaligned.cpp
    ${FASTPFOR_SRC_DIR}/src/bitpackingunaligned.cpp
    ${FASTPFOR_SRC_DIR}/src/simdunalignedbitpacking.cpp
    ${FASTPFOR_SRC_DIR}/src/simdbitpacking.cpp
)

# Helper: define ONE per-tier OBJECT library that bundles the FastPFor
# sources + the ikafssn pfd wrapper TU, all under a tier-specific
# namespace (FastPForLib_<tier>) and -march flag set.
function(ikafssn_add_pfd_tier_object TIER)
    set(_options)
    set(_one)
    set(_multi MARCH_FLAGS INCLUDE_DIRS EXTRA_DEFINES)
    cmake_parse_arguments(ARG "${_options}" "${_one}" "${_multi}" ${ARGN})

    set(_target ikafssn_pfd_${TIER})
    add_library(${_target} OBJECT
        ${IKAFSSN_FASTPFOR_SOURCES}
        ${CMAKE_SOURCE_DIR}/src/index/pfd_codec_tier.cpp
    )
    target_include_directories(${_target} PRIVATE
        ${FASTPFOR_INC_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/generated
        ${ARG_INCLUDE_DIRS}
    )
    target_compile_options(${_target} PRIVATE
        ${ARG_MARCH_FLAGS}
        # Silence a few warnings that the FastPFor sources hit at the
        # warning levels we use; we don't own the upstream code.
        -Wno-unused-parameter
        -Wno-unused-variable
        -Wno-unused-but-set-variable
        -Wno-unused-function
        -Wno-conversion
        -Wno-sign-compare
        -Wno-effc++
        -Wno-shadow
    )
    target_compile_definitions(${_target} PRIVATE
        FastPForLib=FastPForLib_${TIER}
        IKAFSSN_PFD_TIER_NAME=${TIER}
        ${ARG_EXTRA_DEFINES}
    )
    set_property(TARGET ${_target} PROPERTY POSITION_INDEPENDENT_CODE ON)
endfunction()

# The four x86_64 PForDelta tiers (Phase 5f ladder). On non-x86 builds
# (aarch64 etc.) FastPFor's bitpacker primitives — simdpackwithoutmask /
# simdunpack — are SSE2-only and the surrounding -m<x86 ISA> flags are
# rejected by GCC, so we simply do not build these targets there. The
# pfd_codec.cpp dispatcher is `#if`-guarded to match (Phase 5h).
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)")
    set(IKAFSSN_HAS_X86_PFD TRUE)
    ikafssn_add_pfd_tier_object(sse42
        MARCH_FLAGS -msse4.2 -mpopcnt)
    ikafssn_add_pfd_tier_object(avx2
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2)
    ikafssn_add_pfd_tier_object(avx512bw
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2
                    -mavx512f -mavx512bw -mavx512vl -mavx512cd -mavx512dq)
    ikafssn_add_pfd_tier_object(avx512vbmi2
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2
                    -mavx512f -mavx512bw -mavx512vl -mavx512cd -mavx512dq
                    -mavx512vbmi -mavx512vbmi2)
else()
    set(IKAFSSN_HAS_X86_PFD FALSE)
endif()

# AArch64 PForDelta tier (Phase 5h proper).  Single tier object covering
# NEON / SVE / SVE2 capable CPUs alike — SIMDe's SSE → NEON translation
# is the only emulation path, and adding hand-coded SVE bitpackers would
# be orthogonal work not handled here.  The MARCH_FLAGS slot is left empty
# so the compiler's default ISA (or whatever cmake/aarch64.toolchain.cmake
# / IKAFSSN_AARCH64_MARCH supplies via CMAKE_CXX_FLAGS_INIT) carries
# through unchanged.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    set(IKAFSSN_HAS_ARM_PFD TRUE)
    ikafssn_add_pfd_tier_object(neon
        INCLUDE_DIRS  ${SIMDE_INC_DIR}
        EXTRA_DEFINES SIMDE_ENABLE_NATIVE_ALIASES)
else()
    set(IKAFSSN_HAS_ARM_PFD FALSE)
endif()

# --- Phase 5i: per-tier SIMD dedup kernel ---
#
# src/index/seq_id_dedup_tier.cpp is compiled once per ISA tier into
# ikafssn_dedup_<tier> OBJECT libraries.  The runtime dispatcher in
# src/index/seq_id_dedup.cpp picks the right tier at first use based on
# current_simd_cap(), mirroring the .kix/.kpx codec dispatcher.  No
# FastPFor or SIMDe dependency — the kernel is a self-contained scalar
# loop that the compiler auto-vectorises under the per-tier -m<arch>
# flag set.
function(ikafssn_add_dedup_tier_object TIER)
    set(_options)
    set(_one)
    set(_multi MARCH_FLAGS)
    cmake_parse_arguments(ARG "${_options}" "${_one}" "${_multi}" ${ARGN})

    set(_target ikafssn_dedup_${TIER})
    add_library(${_target} OBJECT
        ${CMAKE_SOURCE_DIR}/src/index/seq_id_dedup_tier.cpp
    )
    target_include_directories(${_target} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/generated
    )
    target_compile_options(${_target} PRIVATE ${ARG_MARCH_FLAGS})
    target_compile_definitions(${_target} PRIVATE
        IKAFSSN_DEDUP_TIER_NAME=${TIER}
    )
    set_property(TARGET ${_target} PROPERTY POSITION_INDEPENDENT_CODE ON)
endfunction()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)")
    ikafssn_add_dedup_tier_object(sse42
        MARCH_FLAGS -msse4.2 -mpopcnt)
    ikafssn_add_dedup_tier_object(avx2
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2)
    ikafssn_add_dedup_tier_object(avx512bw
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2
                    -mavx512f -mavx512bw -mavx512vl -mavx512cd -mavx512dq)
    ikafssn_add_dedup_tier_object(avx512vbmi2
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2
                    -mavx512f -mavx512bw -mavx512vl -mavx512cd -mavx512dq
                    -mavx512vbmi -mavx512vbmi2)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    ikafssn_add_dedup_tier_object(neon)
endif()

# --- Phase 7b: per-tier Elias-Fano dictionary codec ---
#
# src/index/ef_codec_tier.cpp is compiled once per ISA tier into the
# ikafssn_ef_<tier> OBJECT libraries.  Each compilation gets its own
# IKAFSSN_EF_TIER_NAME so the per-tier symbols land in disjoint
# namespaces (ikafssn_ef_sse42 / ikafssn_ef_avx2 / ...).  The
# AVX2-and-above builds enable BMI2 (-mbmi2) so the EF select1
# inner-word loop uses PDEP+TZCNT.  No FastPFor / SIMDe dependency
# (the codec is self-contained).
function(ikafssn_add_ef_tier_object TIER)
    set(_options)
    set(_one)
    set(_multi MARCH_FLAGS)
    cmake_parse_arguments(ARG "${_options}" "${_one}" "${_multi}" ${ARGN})

    set(_target ikafssn_ef_${TIER})
    add_library(${_target} OBJECT
        ${CMAKE_SOURCE_DIR}/src/index/ef_codec_tier.cpp
    )
    target_include_directories(${_target} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_BINARY_DIR}/generated
    )
    target_compile_options(${_target} PRIVATE ${ARG_MARCH_FLAGS})
    target_compile_definitions(${_target} PRIVATE
        IKAFSSN_EF_TIER_NAME=${TIER}
    )
    set_property(TARGET ${_target} PROPERTY POSITION_INDEPENDENT_CODE ON)
endfunction()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)")
    ikafssn_add_ef_tier_object(sse42
        MARCH_FLAGS -msse4.2 -mpopcnt)
    ikafssn_add_ef_tier_object(avx2
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2)
    ikafssn_add_ef_tier_object(avx512bw
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2
                    -mavx512f -mavx512bw -mavx512vl -mavx512cd -mavx512dq)
    ikafssn_add_ef_tier_object(avx512vbmi2
        MARCH_FLAGS -mavx2 -mfma -mbmi -mbmi2
                    -mavx512f -mavx512bw -mavx512vl -mavx512cd -mavx512dq
                    -mavx512vbmi -mavx512vbmi2)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    ikafssn_add_ef_tier_object(neon)
endif()

# --- Subdirectories ---
add_subdirectory(src)

# --- Install rules ---
include(GNUInstallDirs)

install(TARGETS
    ikafssnindex ikafssnsearch ikafssnserver ikafssnclient ikafssnretrieve ikafssninfo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
if(BUILD_HTTPD)
    install(TARGETS ikafssnhttpd
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

install(FILES
    doc/systemd/ikafssnserver@.service
    doc/systemd/ikafssnhttpd@.service
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/systemd/system
)
install(FILES
    doc/systemd/README.md
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ikafssn/systemd
)

# --- CPack packaging ---
set(CPACK_PACKAGE_NAME "ikafssn")
set(CPACK_PACKAGE_VERSION "${IKAFSSN_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "K-mer-based alignment-free similarity search for nucleotide sequences")
set(CPACK_PACKAGE_CONTACT "Akifumi S. Tanabe")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)

# RPM-specific
set(CPACK_RPM_PACKAGE_LICENSE "Apache-2.0")
set(CPACK_RPM_PACKAGE_GROUP "Applications/Engineering")
set(CPACK_RPM_PACKAGE_URL "https://github.com/astanabe/ikafssn")
set(CPACK_RPM_PACKAGE_DESCRIPTION
    "ikafssn builds a complete inverted index over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.")
set(CPACK_RPM_PACKAGE_AUTOREQPROV ON)
set(CPACK_RPM_FILE_NAME RPM-DEFAULT)

include(CPack)

enable_testing()
add_subdirectory(test)

# --- Benchmark binaries (optional) ---
if(IKAFSSN_BUILD_BENCH)
    include(FetchContent)
    set(BENCHMARK_ENABLE_TESTING     OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL     OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        benchmark
        URL      https://github.com/google/benchmark/archive/refs/tags/v1.9.4.tar.gz
        URL_HASH SHA256=b334658edd35efcf06a99d9be21e4e93e092bd5f95074c1673d5c8705d95c104
    )
    FetchContent_MakeAvailable(benchmark)
    add_subdirectory(bench)
endif()

# --- Analysis tools (optional, Phase 5a) ---
if(IKAFSSN_BUILD_TOOLS)
    add_subdirectory(tools)
endif()
