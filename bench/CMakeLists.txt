# bench/ — micro-benchmarks (Google Benchmark) and end-to-end shell drivers.
#
# Built only when configured with -DIKAFSSN_BUILD_BENCH=ON. Google Benchmark is
# fetched by the parent CMakeLists.txt via FetchContent and exposed as the
# `benchmark::benchmark` target.

# Phase 1a — FASTA toupper SIMD
add_executable(bench_fasta_toupper bench_fasta_toupper.cpp)
target_link_libraries(bench_fasta_toupper PRIVATE
    ikafssn_io ikafssn_util benchmark::benchmark)
target_include_directories(bench_fasta_toupper PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")
