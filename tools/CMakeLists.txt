# tools/ — analysis utilities for inspecting the on-disk index.
# Only built when configured with -DIKAFSSN_BUILD_TOOLS=ON.

# dump_postings: read .kix/.kpx and emit posting TSV for analyze_postings.py.
add_executable(dump_postings dump_postings.cpp)
target_link_libraries(dump_postings PRIVATE
    ikafssn_index ikafssn_io ikafssn_util ikafssn_core)
target_include_directories(dump_postings PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")

# index_compress_stats: walk .kix/.kpx and report byte breakdown.
add_executable(index_compress_stats index_compress_stats.cpp)
target_link_libraries(index_compress_stats PRIVATE
    ikafssn_index ikafssn_io ikafssn_util ikafssn_core)
target_include_directories(index_compress_stats PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")
