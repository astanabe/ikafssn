# tools/ — throwaway analysis utilities used during Phase 5a (FastPFOR
# feasibility study).  Only built when configured with -DIKAFSSN_BUILD_TOOLS=ON.

# dump_postings: read v3 .kix/.kpx and emit posting TSV for analyze_postings.py
add_executable(dump_postings dump_postings.cpp)
target_link_libraries(dump_postings PRIVATE
    ikafssn_index ikafssn_io ikafssn_util ikafssn_core)
target_include_directories(dump_postings PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")

# index_compress_stats: walk v8 .kix/.kpx binary and report byte breakdown.
add_executable(index_compress_stats index_compress_stats.cpp)
target_link_libraries(index_compress_stats PRIVATE
    ikafssn_index ikafssn_io ikafssn_util ikafssn_core)
target_include_directories(index_compress_stats PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")
