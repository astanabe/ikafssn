cmake_minimum_required(VERSION 3.16)
project(ikafssn LANGUAGES C CXX)

set(IKAFSSN_VERSION "0.1.2026.02.28")
configure_file(
    "${CMAKE_SOURCE_DIR}/src/core/version.hpp.in"
    "${CMAKE_BINARY_DIR}/generated/core/version.hpp"
    @ONLY
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
    find_package(Drogon REQUIRED)
endif()

# --- OpenSSL (for ikafssnclient checkpointing SHA256) ---
find_package(OpenSSL REQUIRED)

# --- libcurl + jsoncpp (optional, for ikafssnclient HTTP mode) ---
option(ENABLE_CLIENT_HTTP "Enable HTTP mode in ikafssnclient (requires libcurl)" ON)
if(ENABLE_CLIENT_HTTP)
    find_package(CURL REQUIRED)
    find_package(jsoncpp REQUIRED)
    # Determine the correct jsoncpp imported target name.
    # Ubuntu provides jsoncpp_lib; some platforms only provide JsonCpp::JsonCpp.
    if(TARGET jsoncpp_lib)
        set(JSONCPP_TARGET jsoncpp_lib)
    elseif(TARGET JsonCpp::JsonCpp)
        set(JSONCPP_TARGET JsonCpp::JsonCpp)
    else()
        message(FATAL_ERROR "jsoncpp found but no known CMake target (jsoncpp_lib or JsonCpp::JsonCpp) available")
    endif()
endif()

# --- libcurl (optional, for ikafssnretrieve remote mode) ---
option(ENABLE_REMOTE_RETRIEVE "Enable NCBI efetch in ikafssnretrieve (requires libcurl)" ON)
if(ENABLE_REMOTE_RETRIEVE)
    find_package(CURL REQUIRED)
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

# --- CPack DEB packaging ---
set(CPACK_PACKAGE_NAME "ikafssn")
set(CPACK_PACKAGE_VERSION "${IKAFSSN_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "K-mer-based alignment-free similarity search for nucleotide sequences")
set(CPACK_PACKAGE_CONTACT "Akifumi S. Tanabe")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
include(CPack)

enable_testing()
add_subdirectory(test)
