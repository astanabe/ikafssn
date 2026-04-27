# ikafssn micro-benchmarks (Google Benchmark)

## Build

    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DIKAFSSN_BUILD_BENCH=ON
    make -j$(nproc)

Google Benchmark v1.9.4 is fetched by the configure step (`FetchContent`) and
linked into each `bench/bench_*.cpp` target. With
`-DIKAFSSN_BUILD_BENCH=OFF` (default), Google Benchmark is **not** fetched and
no bench targets are produced.

## Run

Each binary supports the standard Google Benchmark CLI options:

    ./bench/bench_fasta_toupper                                      # all tiers x sizes
    ./bench/bench_fasta_toupper --benchmark_filter='avx2'
    ./bench/bench_fasta_toupper --benchmark_format=json > result.json

`bench_common.hpp::apply_force_tier()` calls `force_simd_cap(<tier>)` inside
each fixture; if the host CPU cannot satisfy the requested tier the iteration
is reported via `state.SkipWithError(...)` instead of falling back silently.

End-to-end timing of `ikafssnindex` (full index construction across SIMD
tiers) is delegated to a shell driver:

    bench/run_e2e.sh build/

## Available targets

| Target                    | Phase | Notes                                      |
|---------------------------|-------|--------------------------------------------|
| `bench_fasta_toupper`     | 1a    | ASCII bulk to-upper (`io/text_simd`)       |
| `bench_ncbi2na_unpack`    | 1b    | ncbi2na -> 2-bit (`core/ncbi2na_unpack`)   |
