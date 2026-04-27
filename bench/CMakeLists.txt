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

# Phase 1b — ncbi2na unpack SIMD
add_executable(bench_ncbi2na_unpack bench_ncbi2na_unpack.cpp)
target_link_libraries(bench_ncbi2na_unpack PRIVATE
    ikafssn_util benchmark::benchmark)
target_include_directories(bench_ncbi2na_unpack PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")

# Phase 2b — Stage 1 batch SIMD update
add_executable(bench_stage1 bench_stage1.cpp)
target_link_libraries(bench_stage1 PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_util benchmark::benchmark)
target_include_directories(bench_stage1 PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")

# Phase 2c — index_builder posting-buffer SIMD sort
add_executable(bench_parallel_sort bench_parallel_sort.cpp)
target_link_libraries(bench_parallel_sort PRIVATE
    ikafssn_index ikafssn_io ikafssn_util benchmark::benchmark)
target_include_directories(bench_parallel_sort PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")

# Phase 3b — k-mer revcomp batch SIMD
add_executable(bench_kmer_revcomp bench_kmer_revcomp.cpp)
target_link_libraries(bench_kmer_revcomp PRIVATE
    ikafssn_util benchmark::benchmark)
target_include_directories(bench_kmer_revcomp PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}")

