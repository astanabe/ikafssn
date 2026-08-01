# 3.20 is what the conda recipe already requires, and it is the floor for
# FetchContent_MakeAvailable() honouring SOURCE_SUBDIR (see FastPFor below).
cmake_minimum_required(VERSION 3.20)

# Give the archives FetchContent downloads the timestamp of their extraction,
# so a changed URL rebuilds what depends on the contents.  Guarded because the
# policy only exists from CMake 3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

project(ikafssn LANGUAGES C CXX)

set(IKAFSSN_VERSION "0.1.2026.08.02")

# --- Optional benchmark suite (Google Benchmark via FetchContent) ---
# OFF by default; ON pulls Google Benchmark v1.9.4 and adds bench/ subdir.
option(IKAFSSN_BUILD_BENCH
       "Build benchmark binaries under bench/ (off by default)"
       OFF)

# --- Optional analysis tools ---
# OFF by default; ON adds tools/ subdir (posting / index inspection utilities).
option(IKAFSSN_BUILD_TOOLS
       "Build analysis tools under tools/"
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

# --- Compressed I/O codecs (used by ikafssn_io for FASTA / TSV / JSON / FASTA
#     output paths with .gz / .bz2 / .xz / .zst suffixes) ---
find_package(ZLIB    REQUIRED)
find_package(BZip2   REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET libzstd)

# --- Drogon (optional, for ikafssnhttpd) ---
option(BUILD_HTTPD "Build ikafssnhttpd (requires Drogon)" ON)
if(BUILD_HTTPD)
    set(DROGON_DIR "${CMAKE_SOURCE_DIR}/drogon" CACHE PATH
        "Path to Drogon install root (contains lib/cmake/Drogon/)")
    # Prepend DROGON_DIR to CMAKE_PREFIX_PATH so that DrogonConfig.cmake's
    # transitive find_dependency(Trantor) also resolves under the same prefix.
    list(PREPEND CMAKE_PREFIX_PATH "${DROGON_DIR}")
    find_package(Drogon REQUIRED)
    # SQLite3 backs the ikafssnhttpd job store (async REST polling).
    # FindSQLite3.cmake exposes the imported target as `SQLite::SQLite3`
    # on CMake >= 3.14; older systems and some distros only set the
    # SQLite3_LIBRARY / SQLite3_INCLUDE_DIRS variables, so we synthesise
    # the imported target if needed.
    find_package(SQLite3 REQUIRED)
    # Some upstream module (NCBI's bundled cmake exports) creates an empty
    # SQLite::SQLite3 imported target before find_package(SQLite3) runs, in
    # which case FindSQLite3's `if(NOT TARGET SQLite::SQLite3) ... endif()`
    # block never populates IMPORTED_LOCATION.  Resolve the library path
    # by hand via find_library and force-fill the target so generation
    # does not fail with CMP0111 ("missing location property").
    if(NOT SQLite3_LIBRARY)
        find_library(SQLite3_LIBRARY NAMES sqlite3 sqlite)
    endif()
    if(NOT SQLite3_INCLUDE_DIR)
        find_path(SQLite3_INCLUDE_DIR NAMES sqlite3.h)
    endif()
    if(NOT TARGET SQLite::SQLite3)
        add_library(SQLite::SQLite3 UNKNOWN IMPORTED)
    endif()
    set_target_properties(SQLite::SQLite3 PROPERTIES
        IMPORTED_LOCATION "${SQLite3_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SQLite3_INCLUDE_DIR}"
    )
endif()

# --- OpenSSL (for ikafssnclient checkpointing SHA256) ---
find_package(OpenSSL REQUIRED)

# --- libcurl + jsoncpp (optional, for ikafssnclient HTTP mode) ---
option(ENABLE_CLIENT_HTTP "Enable HTTP mode in ikafssnclient (requires libcurl)" ON)
if(ENABLE_CLIENT_HTTP)
    find_package(CURL REQUIRED)
    # jsoncpp ships a config package on every platform we build for.  The
    # module-mode fallback is Drogon's FindJsoncpp.cmake, which find_package
    # (Drogon) puts on the module path; asking for it under the name its file
    # actually uses keeps CMake from matching a module whose case differs.
    find_package(jsoncpp CONFIG QUIET)
    if(NOT jsoncpp_FOUND)
        find_package(Jsoncpp MODULE REQUIRED)
    endif()
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

# libzstd is now an unconditional project-wide dependency (see the codec
# block above).  ikafssnclient / ikafssnhttpd reuse the same PkgConfig::ZSTD
# imported target for their async-job persistence layer.

# --- x86-simd-sort ---
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

# --- FastPFor (per-tier object libraries for runtime dispatch) ---
#
# FastPFor ships with FASTPFOR_SIMD_MODE = "native" by default,
# which bakes the build host's -march into the static archive —
# incompatible with ikafssn's runtime SIMD dispatch policy.  We
# side-step the upstream build entirely by populating the source tree
# (no add_subdirectory) and compiling the FastPFor TUs ourselves into
# four OBJECT libraries, one per x86_64 ISA tier:
#
#   - ikafssn_pfd_sse42         -msse4.2 -mpopcnt
#   - ikafssn_pfd_avx2          -mavx2 -mbmi -mbmi2 -mfma
#   - ikafssn_pfd_avx512bw      -mavx512bw -mavx512vl + AVX2 deps
#   - ikafssn_pfd_avx512vbmi2   -mavx512vbmi2 -mavx512vbmi + above
#
# Driving the per-tier compile from `target_compile_options` ensures
# `-march=native` never reaches `flags.make`, so the resulting object
# files run on any CPU that supports the named tier — a build host
# without AVX-512 VBMI2 can still produce the avx512vbmi2 OBJECT
# library.
#
# Symbol collision (every tier compiles the same FastPForLib::* symbols)
# is avoided by remapping the namespace at compile time via
# `-DFastPForLib=FastPForLib_<tier>`.  pfd_codec.cpp's dispatcher
# selects the right per-tier entry points at first use.
#
# On aarch64, FastPFor's headers/common.h pulls in headers/fastpfor_neon.h,
# which maps the SSE intrinsics it uses onto native NEON.  That header
# targets NEON only, so SVE / SVE2 capable CPUs get no tier object of their
# own: pfd_codec.cpp's dispatcher routes any aarch64 SimdCap >= NEON to the
# single ikafssn_pfd_neon tier object.
include(FetchContent)
FetchContent_Declare(
    fastpfor
    URL      https://github.com/fast-pack/FastPFOR/archive/refs/tags/v0.5.0.tar.gz
    URL_HASH SHA256=d72f61fdef11d5293a4506350eaba3d293b94e955e2b65c7b2558c2cdad8b32c
    # SOURCE_SUBDIR names a directory the archive does not contain, so
    # FetchContent_MakeAvailable() downloads and extracts without calling
    # add_subdirectory(): FastPFor's own CMakeLists would build a "native"
    # libFastPFOR archive we can't use.
    SOURCE_SUBDIR ikafssn-no-add-subdirectory
)
FetchContent_MakeAvailable(fastpfor)

set(FASTPFOR_SRC_DIR "${fastpfor_SOURCE_DIR}")
set(FASTPFOR_INC_DIR "${FASTPFOR_SRC_DIR}/headers")

# FastPFor TUs we compile per tier.  Only the bitpacker primitives are
# included — pfd_codec_tier.cpp drives `simdpackwithoutmask` /
# `simdunpack` directly under a custom block layout (FOR-within-block
# + varint tail).
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
    set(_multi MARCH_FLAGS)
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
    )
    set_property(TARGET ${_target} PROPERTY POSITION_INDEPENDENT_CODE ON)
endfunction()

# The four x86_64 PForDelta tiers. On non-x86 builds (aarch64 etc.)
# FastPFor's bitpacker primitives — simdpackwithoutmask / simdunpack —
# are SSE2-only and the surrounding -m<x86 ISA> flags are rejected by
# GCC, so these targets are not built there.  pfd_codec.cpp's
# dispatcher is `#if`-guarded to match.
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

# AArch64 PForDelta tier.  A single tier object covers NEON / SVE /
# SVE2 capable CPUs alike, since FastPFor's NEON header targets NEON
# only.  The MARCH_FLAGS slot is left empty so the compiler's default
# ISA (or whatever cmake/aarch64.toolchain.cmake / IKAFSSN_AARCH64_MARCH
# supplies via CMAKE_CXX_FLAGS_INIT) carries through unchanged.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(aarch64|arm64)")
    set(IKAFSSN_HAS_ARM_PFD TRUE)
    ikafssn_add_pfd_tier_object(neon)
else()
    set(IKAFSSN_HAS_ARM_PFD FALSE)
endif()

# --- Per-tier SIMD dedup kernel ---
#
# src/index/seq_id_dedup_tier.cpp is compiled once per ISA tier into
# the ikafssn_dedup_<tier> OBJECT libraries.  The runtime dispatcher
# in src/index/seq_id_dedup.cpp picks the right tier at first use
# based on current_simd_cap(), mirroring the .kix/.kpx codec
# dispatcher.  No FastPFor dependency — the kernel is a self-contained
# scalar loop that the compiler auto-vectorises under the per-tier
# -m<arch> flag set.
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

# --- Per-tier Elias-Fano dictionary codec ---
#
# src/index/ef_codec_tier.cpp is compiled once per ISA tier into the
# ikafssn_ef_<tier> OBJECT libraries.  Each compilation gets its own
# IKAFSSN_EF_TIER_NAME so the per-tier symbols land in disjoint
# namespaces (ikafssn_ef_sse42 / ikafssn_ef_avx2 / ...).  The
# AVX2-and-above builds enable BMI2 (-mbmi2) so the EF select1
# inner-word loop uses PDEP+TZCNT.  No FastPFor dependency (the codec
# is self-contained).
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

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    install(FILES
        doc/systemd/ikafssnserver@.service
        doc/systemd/ikafssnhttpd@.service
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/systemd/system
    )
    install(FILES
        doc/systemd/README.md
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ikafssn/systemd
    )
endif()

# --- CPack packaging ---
set(CPACK_PACKAGE_NAME "ikafssn")
set(CPACK_PACKAGE_VERSION "${IKAFSSN_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "K-mer-based alignment-free similarity search for nucleotide sequences")
set(CPACK_PACKAGE_CONTACT "Akifumi S. Tanabe")
# The RPM generator strips debug info via its own debug-info splitting;
# the DEB generator does not, so without this flag the .deb retains
# DWARF sections from statically linked dependencies (htslib's
# ./configure defaults to -O2 -g) and balloons several times over.
set(CPACK_STRIP_FILES TRUE)
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

# --- Analysis tools (optional) ---
if(IKAFSSN_BUILD_TOOLS)
    add_subdirectory(tools)
endif()
