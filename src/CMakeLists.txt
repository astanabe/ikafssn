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
    TBB::tbb
)

# IO library (without NCBI dependency)
add_library(ikafssn_io STATIC
    io/mmap_file.cpp
    io/fasta_reader.cpp
    io/seqidlist_reader.cpp
    io/result_writer.cpp
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

# Search library (no NCBI dependency)
add_library(ikafssn_search STATIC
    search/oid_filter.cpp
    search/stage1_filter.cpp
    search/diagonal_filter.cpp
    search/stage2_chaining.cpp
    search/volume_searcher.cpp
)
target_link_libraries(ikafssn_search PUBLIC ikafssn_core ikafssn_index ikafssn_io)

# Utility library
add_library(ikafssn_util STATIC
    util/cli_parser.cpp
    util/size_parser.cpp
    util/socket_utils.cpp
)
target_link_libraries(ikafssn_util PUBLIC ikafssn_core)

# Protocol library (frame + messages + serializer)
add_library(ikafssn_protocol STATIC
    protocol/frame.cpp
    protocol/serializer.cpp
)
target_link_libraries(ikafssn_protocol PUBLIC ikafssn_core)

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

# ikafssnsearch
add_executable(ikafssnsearch
    ikafssnsearch/main.cpp
)
target_link_libraries(ikafssnsearch PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_util ikafssn_core
    TBB::tbb
)

# ikafssnserver
add_library(ikafssn_server_lib STATIC
    ikafssnserver/server.cpp
    ikafssnserver/connection_handler.cpp
    ikafssnserver/request_processor.cpp
)
target_link_libraries(ikafssn_server_lib PUBLIC
    ikafssn_search ikafssn_protocol ikafssn_index ikafssn_io ikafssn_util ikafssn_core
    TBB::tbb
)

add_executable(ikafssnserver
    ikafssnserver/main.cpp
)
target_link_libraries(ikafssnserver PRIVATE
    ikafssn_server_lib
)

# ikafssnclient
add_library(ikafssn_client_lib STATIC
    ikafssnclient/socket_client.cpp
)
target_link_libraries(ikafssn_client_lib PUBLIC
    ikafssn_protocol ikafssn_util ikafssn_core
)

add_executable(ikafssnclient
    ikafssnclient/main.cpp
)
target_link_libraries(ikafssnclient PRIVATE
    ikafssn_client_lib ikafssn_io ikafssn_util ikafssn_core
)
