# Core library (header-only for now, but define as INTERFACE target)
add_library(ikafssn_core INTERFACE)
target_include_directories(ikafssn_core INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Index library (reader/writer without NCBI dependency)
add_library(ikafssn_index STATIC
    index/ksx_writer.cpp
    index/ksx_reader.cpp
    index/kix_writer.cpp
    index/kix_reader.cpp
    index/kpx_writer.cpp
    index/kpx_reader.cpp
)
target_link_libraries(ikafssn_index PUBLIC ikafssn_core ikafssn_io)

# Index builder library (requires NCBI via blastdb)
add_library(ikafssn_builder STATIC
    index/index_builder.cpp
)
target_include_directories(ikafssn_builder PRIVATE
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_link_libraries(ikafssn_builder PUBLIC
    ikafssn_core ikafssn_index ikafssn_blastdb
)

# IO library (without NCBI dependency)
add_library(ikafssn_io STATIC
    io/mmap_file.cpp
)
target_link_libraries(ikafssn_io PUBLIC ikafssn_core)

# IO library for BLAST DB access (requires NCBI C++ Toolkit)
add_library(ikafssn_blastdb STATIC
    io/blastdb_reader.cpp
)
target_include_directories(ikafssn_blastdb PRIVATE
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_link_libraries(ikafssn_blastdb PUBLIC ikafssn_core
    seqdb xobjutil xobjmgr xncbi xser xutil
)

# Utility library
add_library(ikafssn_util STATIC
    util/cli_parser.cpp
    util/size_parser.cpp
)
target_link_libraries(ikafssn_util PUBLIC ikafssn_core)

# --- Command executables ---

# ikafssnindex
add_executable(ikafssnindex
    ikafssnindex/main.cpp
)
target_include_directories(ikafssnindex PRIVATE
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_link_libraries(ikafssnindex PRIVATE
    ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
    TBB::tbb
)
