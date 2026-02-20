cmake_minimum_required(VERSION 3.16)
project(ikafssn VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- NCBI C++ Toolkit ---
set(NCBI_TOOLKIT_DIR "${CMAKE_SOURCE_DIR}/ncbi-toolkit")
set(NCBI_TOOLKIT_INCLUDE
    "${NCBI_TOOLKIT_DIR}/include"
    "${NCBI_TOOLKIT_DIR}/CMake-GCC1330-Release/inc"
)
set(NCBI_TOOLKIT_CMAKE "${NCBI_TOOLKIT_DIR}/CMake-GCC1330-Release/cmake/ncbi-cpp-toolkit.cmake")

# Define dependency targets required by NCBI export before including it
if(NOT TARGET LMDB::lmdb_LIB)
    add_library(LMDB::lmdb_LIB UNKNOWN IMPORTED)
    set_target_properties(LMDB::lmdb_LIB PROPERTIES
        IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/liblmdb.so"
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

# --- Subdirectories ---
add_subdirectory(src)

enable_testing()
add_subdirectory(test)
