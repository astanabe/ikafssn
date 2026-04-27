# bench/ — micro-benchmarks (Google Benchmark) and end-to-end shell drivers.
#
# Built only when configured with -DIKAFSSN_BUILD_BENCH=ON. Google Benchmark is
# fetched by the parent CMakeLists.txt via FetchContent and exposed as the
# `benchmark::benchmark` target.
#
# Phase 1a/1b targets (added in their respective commits):
#   - bench_fasta_toupper        (Phase 1a)
#   - bench_ncbi2na_unpack       (Phase 1b)
