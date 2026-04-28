# Cross-compile toolchain file for AArch64 Linux (Phase 5h-prelim).
#
# Usage on an x86_64 host with the standard Debian/Ubuntu cross packages
# installed (`sudo apt install g++-aarch64-linux-gnu qemu-user-static
# binfmt-support`):
#
#   cmake -B build-arm64 \
#       -DCMAKE_TOOLCHAIN_FILE=cmake/aarch64.toolchain.cmake \
#       -DCMAKE_BUILD_TYPE=Release \
#       -DBUILD_HTTPD=OFF -DBUILD_CLIENT=OFF \
#       -DENABLE_REMOTE_RETRIEVE=OFF
#   cmake --build build-arm64 --target test_kmer_encoding test_extract_for_mask_batch
#   ctest --test-dir build-arm64 -L simd_variant -R "test_kmer_encoding|test_extract_for_mask_batch"
#
# CMAKE_CROSSCOMPILING_EMULATOR routes test binaries through qemu-aarch64
# transparently so CTest can run them on the x86_64 host.
#
# This file targets armv8-a + NEON + crypto (Cortex-A72 / Pi 4B baseline).
# Override IKAFSSN_AARCH64_MARCH on the cmake command line for newer
# CPU features:
#   -DIKAFSSN_AARCH64_MARCH=armv9-a+sve2  (CIX CP8180 / Cortex-A720+A520)
#   -DIKAFSSN_AARCH64_MARCH=armv8.2-a+sve (Neoverse V1, A64FX, etc.)
#
# NCBI C++ Toolkit / Drogon / Parasail / htslib are *not* cross-compiled
# by this file; binaries that depend on those libraries (ikafssnindex,
# ikafssnsearch, test_builder, etc.) cannot be cross-built without a
# pre-built aarch64 sysroot.  Use this file only for the SIMD-kernel
# subset (test_kmer_encoding, test_extract_for_mask_batch,
# test_packed_kmer_scanner_*, etc.) — anything heavier should be built
# natively on the AArch64 target instead.

set(CMAKE_SYSTEM_NAME      Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)

set(CMAKE_C_COMPILER   aarch64-linux-gnu-gcc)
set(CMAKE_CXX_COMPILER aarch64-linux-gnu-g++)

set(CMAKE_FIND_ROOT_PATH /usr/aarch64-linux-gnu)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

if(NOT DEFINED IKAFSSN_AARCH64_MARCH)
    set(IKAFSSN_AARCH64_MARCH "armv8-a+crc+crypto")
endif()
set(CMAKE_C_FLAGS_INIT   "-march=${IKAFSSN_AARCH64_MARCH}")
set(CMAKE_CXX_FLAGS_INIT "-march=${IKAFSSN_AARCH64_MARCH}")

# Run aarch64 binaries on the x86_64 host via qemu-user.  Used by CTest.
find_program(IKAFSSN_QEMU_AARCH64 qemu-aarch64-static)
if(IKAFSSN_QEMU_AARCH64)
    # `-cpu max` enables NEON + every optional ARMv8/v9 feature QEMU
    # knows about, including SVE2.  Override per-invocation if you want
    # to test a narrower CPU profile.
    set(CMAKE_CROSSCOMPILING_EMULATOR
        "${IKAFSSN_QEMU_AARCH64};-cpu;max,sve=on,sve2=on")
endif()
