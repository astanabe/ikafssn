# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ikafssn (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) builds a complete inverted index over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.

- **Primary**: https://github.com/astanabe/ikafssn
- **Secondary**: https://gitlab.com/astanabe/ikafssn

## Build Commands

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make test          # run all unit tests
ctest --test-dir build -R test_kmer_encoding   # run a single test
```

CMake options for selective builds:
- `-DBUILD_HTTPD=OFF` — skip ikafssnhttpd (requires Drogon)
- `-DBUILD_CLIENT=OFF` — skip ikafssnclient (requires libcurl for HTTP mode)
- `-DENABLE_REMOTE_RETRIEVE=OFF` — disable NCBI efetch in ikafssnretrieve

## Dependencies

This project's external dependencies and their status on this dev machine (Ubuntu 24.04.4 LTS, x86_64):

| Library | Version | Location | Install method |
|---|---|---|---|
| g++ | 13.3.0 | system | apt (pre-installed) |
| CMake | 3.28.3 | system | apt (pre-installed) |
| NCBI C++ Toolkit 30.0.0 | 30.0.0 | `./ncbi-toolkit/` | source build |
| Intel TBB (oneTBB) | 2021.11.0 | system | `sudo apt install libtbb-dev` |
| Drogon | 1.8.7 | system | `sudo apt install libdrogon-dev` |
| libcurl | 8.5.0 | system | apt (pre-installed) |
| BLAST+ | 2.12.0 | system | apt (pre-installed, for test data generation) |

### NCBI C++ Toolkit

Installed at `./ncbi-toolkit/` (project-local, built from source):
- Headers: `ncbi-toolkit/include/`
- Libraries: `ncbi-toolkit/CMake-GCC1330-Release/lib/` (static `.a`)
- CMake exports: `ncbi-toolkit/CMake-GCC1330-Release/cmake/ncbi-cpp-toolkit.cmake`
- Exported targets: `seqdb`, `blastdb_format`, `xobjutil`, `xobjmgr`, `xncbi`, `xser`, `xutil`, etc.

## Architecture

### Command binaries (separate executables, not subcommands)

Each command links only its required dependencies to allow lightweight deployment:

| Command | Purpose | Key Dependencies |
|---|---|---|
| `ikafssnindex` | Build k-mer inverted index from BLAST DB | NCBI C++ Toolkit, TBB |
| `ikafssnsearch` | Local direct search (mmap index) | NCBI C++ Toolkit, TBB |
| `ikafssnretrieve` | Extract matched subsequences | NCBI C++ Toolkit, libcurl (remote) |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) | NCBI C++ Toolkit, TBB |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver | Drogon |
| `ikafssnclient` | Client (socket or HTTP) | libcurl (HTTP mode) |
| `ikafssninfo` | Index/DB info display | NCBI C++ Toolkit (optional) |

### Shared library layers (`src/`)

- **`core/`** — Fundamental types (`Hit`, `ChainResult`), constants, k-mer 2-bit encoding/revcomp (templates parameterized on `KmerInt` = `uint16_t` for k<=8, `uint32_t` for k=9-13), LEB128 varint
- **`index/`** — Reader/writer for four index file formats (`.kix` main index, `.kpx` position data, `.ksx` sequence metadata, `.khx` build-time exclusion bitset), index builder with partition+buffer strategy
- **`search/`** — Two-stage search pipeline: Stage 1 (ID posting scan with OID filter, coverscore or matchscore) -> Stage 2 (position-aware chaining DP with diagonal filter, chainscore). Mode 1 skips Stage 2 entirely. Configurable sort_score (stage1 score or chainscore). Defaults: stage1_topn=0 (unlimited, no sort), num_results=0 (unlimited, no sort), min_stage1_score=0.5 (fractional), min_score=1 — speed-first defaults that skip sorting. Set positive stage1_topn/num_results to enable sorting. Fractional min_stage1_score (0 < P < 1) resolves per-query threshold as `ceil(Nqkmer * P) - Nhighfreq`, using `.khx` for build-time exclusion awareness
- **`protocol/`** — Length-prefixed binary protocol for client-server communication (frame header + typed messages)
- **`io/`** — BLAST DB reader (CSeqDB wrapper), FASTA reader, mmap RAII wrapper, seqidlist reader (text/binary auto-detect), result writer/reader
- **`util/`** — CLI parser, size string parser ("8G"), socket utilities, progress display, logger

### Index file formats

Per BLAST DB volume, three files are generated with naming pattern `<db_prefix>.<volume_index>.<kk>mer.{kix,kpx,ksx}` where `<kk>` is the zero-padded 2-digit k value (e.g. `nt.00.09mer.kix`, `nt.01.11mer.kpx`):

- **`.kix`** — Header + direct-address table (offsets[4^k], counts[4^k]) + delta-compressed ID postings
- **`.kpx`** — Header + pos_offsets[4^k] + delta-compressed position postings (correlated with .kix ID postings)
- **`.ksx`** — Header + seq_lengths[] + accession string table (enables standalone result display without BLAST DB)
- **`.khx`** — Header (32B, magic "KMHX") + bitset (ceil(4^k / 8) bytes). Generated only when `-max_freq_build` is used. Records which k-mers were excluded during index build (bit i = 1 means k-mer i was excluded). Used at search time for fractional `-min_stage1_score` threshold adjustment

ID and position postings are in separate files so Stage 1 never touches `.kpx`, maximizing page cache efficiency.

### Key design conventions

- **C++17**, CMake >= 3.16, little-endian only (Linux x86_64/aarch64)
- k-mer values use direct-address tables (array index = k-mer integer value), no hash maps
- Template dispatch on `KmerInt` type; runtime dispatch reads `kmer_type` from index header, then calls the correct template instantiation — no virtual calls in hot loops
- Delta encoding uses LEB128 varint; position deltas reset at sequence boundaries (detected by checking if the corresponding ID delta is non-zero)
- N-containing k-mers are skipped via an N-counter (shared logic between indexing and query scanning)
- Reverse complement: index stores forward strand only; search generates both forward and revcomp of query
- Parallelization via Intel TBB: indexing uses `parallel_for` + `combinable` for counting and partition scan, `parallel_sort` for posting buffer, `task_group` for volume-level parallelism; search uses `parallel_for_each` for (query, volume) jobs. Thread count is controlled centrally via `tbb::global_control` in each main binary
- mmap is read-only and shared across threads; `score_per_seq` arrays are thread-local
- One server process serves one BLAST DB; multiple DBs require multiple processes

## Documentation

Usage documentation is maintained in two languages:
- `doc/ikafssn.en.md` — English
- `doc/ikafssn.ja.md` — Japanese

## Test structure

Tests are in `test/` using CTest. Real SSU_eukaryote_rRNA BLAST DB at `db/SSU_eukaryote_rRNA` is used for all BLAST-DB-dependent tests. `test/scripts/setup_ssu_testdata.sh` generates derived test data (ambig DB, multi-volume DBs, queries) in `/tmp/ikafssn_ssu_test/`. Shared fixture: `test/ssu_test_fixture.hpp`.

## Development Environment Rules

- `sudo` が必要なコマンドは Claude Code から直接実行しない。ユーザーに提示して手動実行を求めること。
