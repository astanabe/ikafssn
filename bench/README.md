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

## Available targets

| Target                  | Notes                                                  |
|-------------------------|--------------------------------------------------------|
| `bench_fasta_toupper`   | ASCII bulk to-upper (`io/text_simd`)                   |
| `bench_ncbi2na_unpack`  | ncbi2na -> 2-bit (`core/ncbi2na_unpack`)               |
| `bench_kmer_revcomp`    | `kmer_revcomp_batch` per tier vs scalar                |
| `bench_degenerate_scan` | `has_degenerate_base` per tier                         |
| `bench_extract_for_mask`| `extract_for_mask_batch` (spaced seed) per tier        |
| `bench_ef_codec`        | Elias-Fano dictionary encode / `EFDictionary::access`  |
| `bench_parallel_sort`   | index_builder posting list buffer sort per tier        |
| `bench_stage1`          | `flush_batch_simd` scatter kernel per Stage 1 width    |

## End-to-end shell drivers

These are not Google Benchmark targets; they drive the installed binaries and
need no `-DIKAFSSN_BUILD_BENCH=ON`.

| Script                 | Measures                                                    |
|------------------------|-------------------------------------------------------------|
| `run_e2e.sh`           | `ikafssnindex` full index construction, one run per SIMD tier |
| `run_e2e_search.sh`    | `ikafssnsearch` wall time, one run per SIMD tier; builds the index itself |
| `run_warm_e2e.sh`      | `ikafssnsearch` against an existing index, repeated, with per-stage attribution |
| `run_cold_e2e.sh`      | two `ikafssnsearch` builds A/B, index evicted before every run |

    bench/run_e2e.sh build/
    bench/run_e2e_search.sh build/ --format json

`run_warm_e2e.sh` is the one to reach for when comparing two builds of the
search path.  It repeats a fixed command `REPS` times per sweep point, drops
the first run so every timed run is page-cache warm, and reports the median
and minimum of

- `/usr/bin/time -v` counters (wall, user, max RSS, minor page faults), and
- every key of the `Timing run_search (s):` / `Timing overall (s):` lines that
  `ikafssnsearch` prints at `info` level (the latter prefixed `overall_`).

It also records the result file's line count and `sort | sha256sum` on every
run, so output equivalence against a baseline is a direct comparison, and runs
`perf stat` once outside the timed set so its overhead stays out of the median.
`-stage1_max_nhit_per_volume` is swept via `M_SWEEP` (the token `default` omits
the flag).  Exit code 2 — some query skipped, e.g. shorter than
`-min_query_length` — counts as success.

    BUILD=build IX=/path/to/index/prefix QUERY=queries.fasta \
    K=11 T=21 TEMPLATE_TYPE=both MODE=2 STRAND=2 MIN_SCORE=0.1 \
    M_SWEEP="100 1000 10000" REPS=5 OUTDIR=/tmp/warm \
    bench/run_warm_e2e.sh --format json

`run_cold_e2e.sh` is the opposite regime: it is for an index too large to stay
resident, where the search is I/O bound rather than CPU bound.  It evicts the
index with `posix_fadvise(POSIX_FADV_DONTNEED)` before every run — no
privileges needed, and no other process's cache is disturbed — and alternates
the two binaries so a residual cache trend hits both alike.  Without the
eviction the page cache carries over between runs and the same binary varies by
2x or more.  Alongside wall time it reports the read volume summed over the
`/proc/diskstats` devices in `DISKS` and the major fault count: a change that
leaves those two alone has not altered which posting bytes the search touches.

    BIN_A=/usr/bin/ikafssnsearch BIN_B=build/src/ikafssnsearch \
    IX=/path/to/index/prefix QUERY=queries.fasta REPS=3 \
    EXTRA_ARGS="-k 11 -t 21 -template_type both -mode 1 -strand 2" \
    bench/run_cold_e2e.sh --format summary

Output equivalence is reported as line count and sorted sha256, not a byte
comparison: with a large result set mode 1 does not fix the row order within a
query.

See the header comment of each script for the full environment variable list.
`parse_time.py` turns a single `/usr/bin/time -v` log into the same JSON
record shape, and `compare.py` diffs two Google Benchmark JSON outputs.
