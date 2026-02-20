cmake_minimum_required(VERSION 3.16)
project(ikafssn VERSION 0.1.0 LANGUAGES C CXX)

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
    set_target_properties(NCBI::DL_LIBS PROPERTIES
        INTERFACE_LINK_LIBRARIES "dl"
    )
endif()

if(NOT TARGET NCBI::RT_LIBS)
    add_library(NCBI::RT_LIBS INTERFACE IMPORTED)
    set_target_properties(NCBI::RT_LIBS PROPERTIES
        INTERFACE_LINK_LIBRARIES "rt"
    )
endif()

if(NOT TARGET NCBI::MATH_LIBS)
    add_library(NCBI::MATH_LIBS INTERFACE IMPORTED)
    set_target_properties(NCBI::MATH_LIBS PROPERTIES
        INTERFACE_LINK_LIBRARIES "m"
    )
endif()

include("${NCBI_TOOLKIT_CMAKE}")

# --- Intel TBB ---
find_package(TBB REQUIRED)

# --- Drogon (optional, for ikafssnhttpd) ---
option(BUILD_HTTPD "Build ikafssnhttpd (requires Drogon)" ON)
if(BUILD_HTTPD)
    find_package(Drogon REQUIRED)
endif()

# --- libcurl + jsoncpp (optional, for ikafssnclient HTTP mode) ---
option(ENABLE_CLIENT_HTTP "Enable HTTP mode in ikafssnclient (requires libcurl)" ON)
if(ENABLE_CLIENT_HTTP)
    find_package(CURL REQUIRED)
    find_package(jsoncpp REQUIRED)
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
    doc/systemd/README.md
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ikafssn/systemd
)

enable_testing()
add_subdirectory(test)
