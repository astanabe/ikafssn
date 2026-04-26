#!/bin/bash
#
# conda-build script for ikafssn-nohttpd.
#
# Builds three vendored static libraries inside the conda-build sandbox so the
# resulting binaries only request glibc symbols up to c_stdlib_version (2.34
# on Linux) and macOS 11.0 on osx-arm64:
#   - Parasail 2.6.2  (Stage 3 alignment, with degmatch CIGAR/score patch)
#   - htslib   1.23.1 (SAM/BAM output, libcurl/gcs/s3 disabled)
#   - NCBI C++ Toolkit 30.2.0 (BLAST DB seqdb_reader + blastdb_format)
#
# Then builds ikafssn itself with -DBUILD_HTTPD=OFF and installs into ${PREFIX}.
#
# All other dependencies (TBB, LMDB, libsqlite, libcurl, jsoncpp, OpenSSL,
# zlib, bzip2, xz, libdeflate) are taken from ${PREFIX} (conda-forge host
# packages). The conda-forge cross compiler chain is used so that the host
# system glibc never leaks into the binaries.

set -euo pipefail

NPROC="${NPROC:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}"

PARASAIL_VER="2.6.2"
HTSLIB_VER="1.23.1"
NCBI_VER="30.2.0"

PARASAIL_PREFIX="${SRC_DIR}/_parasail"
HTSLIB_PREFIX="${SRC_DIR}/_htslib"
NCBI_PREFIX="${SRC_DIR}/_ncbi-cxx-toolkit"

# ---- Parasail (static) ----

cd "${SRC_DIR}"
curl --retry 5 --retry-delay 2 -fsSL -o parasail.tar.gz \
  "https://github.com/jeffdaily/parasail/archive/refs/tags/v${PARASAIL_VER}.tar.gz"
mkdir -p parasail-src
tar xf parasail.tar.gz -C parasail-src --strip-components=1

cd parasail-src
patch -p1 < "${SRC_DIR}/patches/parasail-degmatch-cigar-score.patch"
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PARASAIL_PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
make -j"${NPROC}" parasail
mkdir -p "${PARASAIL_PREFIX}/lib" "${PARASAIL_PREFIX}/include"
cp libparasail.a "${PARASAIL_PREFIX}/lib/"
cp -r ../parasail.h ../parasail "${PARASAIL_PREFIX}/include/"

# ---- htslib (static) ----

cd "${SRC_DIR}"
curl --retry 5 --retry-delay 2 -fsSL -o htslib.tar.bz2 \
  "https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2"
mkdir -p htslib-src
tar xf htslib.tar.bz2 -C htslib-src --strip-components=1

cd htslib-src
autoreconf -i
./configure \
  --prefix="${HTSLIB_PREFIX}" \
  --disable-libcurl --disable-gcs --disable-s3 \
  CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS:-}" \
  LDFLAGS="-L${PREFIX}/lib ${LDFLAGS:-}"
make -j"${NPROC}"
make install

# ---- NCBI C++ Toolkit (static) ----

cd "${SRC_DIR}"
curl --retry 5 --retry-delay 2 -fsSL -o ncbi.tar.gz \
  "https://github.com/ncbi/ncbi-cxx-toolkit-public/archive/refs/tags/release/${NCBI_VER}.tar.gz"
mkdir -p ncbi-src
tar xf ncbi.tar.gz -C ncbi-src --strip-components=1

cd ncbi-src
patch -p1 < "${SRC_DIR}/patches/ncbi-cxx-toolkit-seqdb-madvise-random.patch"

# Compute the sysroot path used by the conda-build toolchain.
#   - Linux: ${BUILD_PREFIX}/${HOST}/sysroot ships glibc 2.34 stub libraries.
#     The host system glibc (e.g. Ubuntu 24.04's 2.39 with C23 strtol overload)
#     must never be consulted.
#   - macOS: there is no separate sysroot directory; conda-build instead
#     exports CONDA_BUILD_SYSROOT pointing at the appropriate macOS SDK.
SYSROOT=""
if [ "$(uname)" = "Darwin" ]; then
  if [ -n "${CONDA_BUILD_SYSROOT:-}" ] && [ -d "${CONDA_BUILD_SYSROOT}" ]; then
    SYSROOT="${CONDA_BUILD_SYSROOT}"
  fi
else
  if [ -n "${BUILD_PREFIX:-}" ] && [ -n "${HOST:-}" ] && \
     [ -d "${BUILD_PREFIX}/${HOST}/sysroot" ]; then
    SYSROOT="${BUILD_PREFIX}/${HOST}/sysroot"
  fi
fi

# NCBI's cmake-cfg-unix.sh derives CC_NAME from `$CC --version`:
#   - On Linux it uppercases the first token. conda-build supplies
#     CC=x86_64-conda-linux-gnu-cc which would yield
#     CC_NAME=X86_64-CONDA-LINUX-GNU-CC and cause cmkTool.sh to skip its
#     GCC-specific configuration (cmGCC.sh), producing a broken CMAKE_ARGS.
#   - On macOS it takes the second token of conda-forge's clang banner
#     ("arm64-apple-darwin*-clang version X") which yields "version" and
#     fails the same way.
# Both cases are avoided by exposing plain `gcc`/`g++` (Linux) or
# `clang`/`clang++` (macOS) symlinks: cmake-cfg-unix.sh's case statement
# matches `clang*` on basename, and `gcc --version` prints "gcc (...) X"
# whose first token is "gcc". Symlinks (rather than wrapper scripts) are
# important on Linux: argv[0]="gcc" is set automatically, and GCC walks
# argv[0]'s realpath to locate its libexec/cc1, so a wrapper that uses
# `exec -a gcc x86_64-conda-linux-gnu-cc` would break that lookup.
if [ "$(uname)" = "Darwin" ]; then
  NCBI_C_BASENAME="clang"
  NCBI_CXX_BASENAME="clang++"
else
  NCBI_C_BASENAME="gcc"
  NCBI_CXX_BASENAME="g++"
fi

NCBI_SHIM_BIN="${SRC_DIR}/_ncbi_compiler_shim"
mkdir -p "${NCBI_SHIM_BIN}"
# Resolve to absolute paths so the symlinks are valid regardless of CWD.
NCBI_SHIM_CC=$(command -v "${CC}")
NCBI_SHIM_CXX=$(command -v "${CXX}")
ln -sf "${NCBI_SHIM_CC}"  "${NCBI_SHIM_BIN}/${NCBI_C_BASENAME}"
ln -sf "${NCBI_SHIM_CC}"  "${NCBI_SHIM_BIN}/cc"
ln -sf "${NCBI_SHIM_CXX}" "${NCBI_SHIM_BIN}/${NCBI_CXX_BASENAME}"
ln -sf "${NCBI_SHIM_CXX}" "${NCBI_SHIM_BIN}/c++"

NCBI_SAVED_CMAKE_ARGS="${CMAKE_ARGS:-}"
NCBI_SAVED_PATH="${PATH}"
NCBI_SAVED_CC="${CC:-}"
NCBI_SAVED_CXX="${CXX:-}"
NCBI_SAVED_CFLAGS="${CFLAGS:-}"
NCBI_SAVED_CXXFLAGS="${CXXFLAGS:-}"

export PATH="${NCBI_SHIM_BIN}:${PATH}"
export CC="${NCBI_C_BASENAME}"
export CXX="${NCBI_CXX_BASENAME}"

# conda-build exports CMAKE_ARGS pre-populated with -DCMAKE_AR, -DCMAKE_RANLIB,
# --sysroot, -B build, etc. NCBI's cmake-cfg-unix.sh appends to $CMAKE_ARGS
# without initializing it, which would inject those flags into the cmake
# invocation (and their unquoted spaces break the eval call). Clear it.
unset CMAKE_ARGS

# Append C++20 to CXXFLAGS so the toolkit itself compiles in C++20 mode
# (matching ikafssn).
export CXXFLAGS="${CXXFLAGS:-} -std=c++20"

# Pass the sysroot through compiler/linker flags so the host system headers
# never leak in (on Linux this prevents glibc-2.39 __isoc23_strtol from being
# pulled in; on macOS this points clang at the conda-forge SDK).
ncbi_extra_cmake_args=()
if [ -n "${SYSROOT}" ]; then
  if [ "$(uname)" = "Darwin" ]; then
    ncbi_extra_cmake_args+=(
      "-DCMAKE_OSX_SYSROOT=${SYSROOT}"
    )
    if [ -n "${MACOSX_DEPLOYMENT_TARGET:-}" ]; then
      ncbi_extra_cmake_args+=(
        "-DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET}"
      )
    fi
  else
    ncbi_extra_cmake_args+=(
      "-DCMAKE_SYSROOT=${SYSROOT}"
      "-DCMAKE_C_FLAGS=--sysroot=${SYSROOT}"
      "-DCMAKE_CXX_FLAGS=--sysroot=${SYSROOT}"
      "-DCMAKE_EXE_LINKER_FLAGS=--sysroot=${SYSROOT}"
      "-DCMAKE_SHARED_LINKER_FLAGS=--sysroot=${SYSROOT}"
    )
  fi
fi

./cmake-configure \
  --without-debug \
  --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" \
  --with-install="${NCBI_PREFIX}" \
  "${ncbi_extra_cmake_args[@]}"

# cmake-configure creates a directory like CMake-GCC1240-Release/build for the
# detected compiler. Glob and enter it.
build_dir=$(echo CMake-*/build)
cd "${build_dir}"
make -j"${NPROC}"
make install

# Detect the build tag (e.g. CMake-GCC1240-Release) the install populated.
ncbi_build_tag=$(basename "${NCBI_PREFIX}"/CMake-*/)
echo "Detected NCBI build tag: ${ncbi_build_tag}"

# Restore conda-build's compiler/flag environment for the ikafssn build,
# which is a normal CMake project that expects the standard cross-compile
# variables (CC=x86_64-conda-linux-gnu-cc, populated CMAKE_ARGS, etc.).
export PATH="${NCBI_SAVED_PATH}"
export CC="${NCBI_SAVED_CC}"
export CXX="${NCBI_SAVED_CXX}"
export CFLAGS="${NCBI_SAVED_CFLAGS}"
export CXXFLAGS="${NCBI_SAVED_CXXFLAGS}"
if [ -n "${NCBI_SAVED_CMAKE_ARGS}" ]; then
  export CMAKE_ARGS="${NCBI_SAVED_CMAKE_ARGS}"
fi

# ---- ikafssn ----

cd "${SRC_DIR}"
ikafssn_extra_cmake_args=()
if [ -n "${SYSROOT}" ]; then
  if [ "$(uname)" = "Darwin" ]; then
    ikafssn_extra_cmake_args+=("-DCMAKE_OSX_SYSROOT=${SYSROOT}")
    if [ -n "${MACOSX_DEPLOYMENT_TARGET:-}" ]; then
      ikafssn_extra_cmake_args+=(
        "-DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET}"
      )
    fi
  else
    ikafssn_extra_cmake_args+=("-DCMAKE_SYSROOT=${SYSROOT}")
  fi
fi
cmake -S . -B build-ikafssn \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DNCBI_TOOLKIT_DIR="${NCBI_PREFIX}" \
  -DNCBI_TOOLKIT_BUILD_TAG="${ncbi_build_tag}" \
  -DPARASAIL_DIR="${PARASAIL_PREFIX}" \
  -DHTSLIB_DIR="${HTSLIB_PREFIX}" \
  -DBUILD_HTTPD=OFF \
  "${ikafssn_extra_cmake_args[@]}"
cmake --build build-ikafssn -j"${NPROC}"
cmake --install build-ikafssn

# Strip systemd unit files from the conda package — these are .deb/.rpm
# specific and have no meaning inside a conda environment.
rm -rf "${PREFIX}/lib/systemd" "${PREFIX}/share/ikafssn/systemd"
