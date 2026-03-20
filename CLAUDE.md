# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **重要**: ユーザーとの画面上のやり取りには**常に日本語**を使用すること。ただし、作成・編集するファイルは日本語ファイルと指定されない限り英語で記述する。このリポジトリで日本語ファイルは `doc/ikafssn.ja.md` のみ。

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
| NCBI C++ Toolkit 30.1.0 | 30.1.0 | `./ncbi-toolkit/` | source build |
| Intel TBB (oneTBB) | 2021.11.0 | system | `sudo apt install libtbb-dev` |
| Drogon | 1.8.7 | system | `sudo apt install libdrogon-dev` |
| libcurl | 8.5.0 | system | apt (pre-installed) |
| Parasail | 2.6.2 | `./parasail/` | source build |
| htslib | 1.23 | `./htslib/` | source build |
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
| `ikafssnindex` | Build k-mer inverted index from BLAST DB (`-mode 1` skips `.kpx`; `-template_type both` builds coding and optimal indexes sequentially) | NCBI C++ Toolkit, TBB |
| `ikafssnsearch` | Local direct search (mmap index), in-silico PCR primer mode | NCBI C++ Toolkit, TBB, Parasail (mode 3), htslib (SAM/BAM) |
| `ikafssnretrieve` | Extract matched subsequences | NCBI C++ Toolkit, libcurl (remote) |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) | NCBI C++ Toolkit, TBB |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver(s) (multi-backend) | Drogon |
| `ikafssnclient` | Client (socket or HTTP), in-silico PCR primer mode | libcurl (HTTP mode) |
| `ikafssninfo` | Index/DB info display (local or remote) | NCBI C++ Toolkit (local), socket/HTTP (remote) |

### Shared library layers (`src/`)

- **`core/`** — Fundamental types (`Hit`, `ChainResult`), constants, k-mer 2-bit encoding/revcomp (templates parameterized on `KmerInt` = `uint16_t` or `uint32_t`, selected by `kmer_type_for(k, t)` based on bit width: `2*k` bits), `AmbigInfo` struct + `expand_ambig_kmer_multi()` for configurable multi-position degenerate expansion (`max_expansion` parameter controls product limit), LEB128 varint, spaced seed (discontiguous megablast) support (`spaced_seed.hpp`: `TemplateType` enum, 24 template masks for k=8/9 × t=13/15/18 (derived from discontiguous MegaBLAST template design principles) and k=11/12 × t=16/18/21 (MegaBLAST-native) × coding/optimal, `get_seed_masks()`, `reverse_complement_string()` for string-level RC used by spaced seed query scanning). "Both" search mode merges separate coding and optimal indexes at search time rather than using a single combined index
- **`index/`** — Reader/writer for four index file formats (`.kix` main index, `.kpx` position data, `.ksx` sequence metadata, `.khx` build-time exclusion bitset). `.kix` and `.kpx` use format version 3: `.kix` v3 removes the counts table and uses sentinel-based offsets (`offsets[table_size+1]`, posting byte length = `offsets[kmer+1]-offsets[kmer]`); both `.kix` v3 and `.kpx` v3 support conditional uint32_t offsets when posting data < 4 GiB (`KIX_FLAG_OFFSET32` flag / `KpxHeader.offset_type` field), halving dictionary size. `KixReader`/`KpxReader` use dual pointers (`offsets32_`/`offsets64_`) with inline branch dispatch. `KixReader` provides `posting_offset()`, `posting_byte_length()`, `count_postings()` (on-demand varint decode), and `bulk_count_postings()`. Table index types use `uint32_t` (max table size = 4^12 = 16M). Index builder with partition+buffer strategy (`skip_kpx` config flag omits `.kpx` for mode 1 indexes); builder rewrites files at finalize to select optimal offset width. For spaced seeds (`-t > 0`), builder uses `scan_spaced()`; `-template_type both` builds coding and optimal indexes sequentially (one after the other). All readers expose `willneed_size()` / `apply_madvise(bool)` for budget-based madvise allocation (used by `-memory_limit` in search/server)
- **`search/`** — Three-stage search pipeline: Stage 1 (ID posting scan with OID filter, coverscore or matchscore) -> Stage 2 (position-aware chaining DP with diagonal filter, chainscore) -> Stage 3 (Parasail pairwise alignment with BLAST DB subject retrieval, alnscore/ppositive/CIGAR). Mode 1 skips Stage 2+3, mode 2 skips Stage 3. The sort key is auto-determined by mode (1=stage1, 2=chainscore, 3=alnscore). Defaults: stage1_topn=0 (unlimited, no sort), num_results=0 (unlimited, no sort), stage1_min_score=0.5 (fractional), stage2_min_score=1 — speed-first defaults that skip sorting. Set positive stage1_topn/num_results to enable sorting. Fractional stage1_min_score (0 < P < 1) resolves per-query threshold as `ceil(Nqkmer * P) - Nhighfreq`, using `.khx` for build-time exclusion awareness. `QueryKmerData` uses SoA layout (`fwd_positions`/`fwd_kmer_values` separate vectors) to eliminate pair padding overhead. `stage1_filter` takes SoA pointers `(const uint32_t* positions, const KmerInt* kmers, size_t n)`. Stage 1 uses byte-limit decoding (`SeqIdDecoder(data, end)` + `has_more()`) instead of count-based loops, matching the counts-free `.kix` v3 format. `Stage1Buffer` uses AoS `Stage1Entry<Tier>` with 3-tier template specialization (`Stage1Tier::T8` = 2B, `T16` = 4B, `T32` = 8B per entry); tier is selected from actual preprocessed k-mer position counts via `select_tier()`. Query k-mers with IUPAC degenerate bases are expanded when the product of per-position variant counts <= `max_degen_expand` (configurable, default 16 for search/server, 4 for index); Stage 1 uses per-position dedup (`last_pos` in AoS entry) to avoid inflating scores from expanded k-mers; Stage 2 deduplicates `(q_pos, s_pos)` hits before diagonal filter and uses `span` (= t for spaced seeds, k for contiguous) for chain endpoint coordinates; `max_nhit_per_subject` (default 1) controls how many non-overlapping chains are extracted per subject via greedy best-chain removal; Stage 3 uses parasail_sg semi-global alignment with configurable score matrix (DEGMATCH default, also dnafull/nuc44 via `-stage3_score_matrix`), optional traceback for CIGAR/ppositive/qseq/sseq, context extension for subject subsequence retrieval (with overlap resolution for multi-chain hits when context > 0), and volume-parallel BLAST DB prefetch. CIGAR match (`=`) / mismatch (`X`) is determined by score matrix (score > 0 = match). `OutputHit.oid`/`OutputHit.volume` carry the BLAST DB OID and volume index from Stage 2, enabling Stage 3 to directly open the correct volume reader and fetch by OID without building a full accession-to-OID map. For spaced seeds, query preprocessing uses string-level reverse complement instead of `kmer_revcomp()`, and `SearchConfig.t` propagates the template length to all pipeline stages. "Both" template_type search uses `search_volume_both()` to merge results from separate coding and optimal indexes at search time. In-silico PCR primer mode (`-primer`) is supported by both `ikafssnsearch` and `ikafssnclient`; see `io/primer_query.hpp/cpp` for the shared primer pair parsing logic
- **`protocol/`** — Length-prefixed binary protocol (v8) for client-server communication. Frame header (12B, magic + payload_size + msg_type + msg_version=8 + reserved) + typed messages. SearchRequest/Response include `db` for multi-DB routing. SearchRequest includes `max_degen_expand`, `stage2_max_nhit_per_subject`, `t`, `template_type`, and `score_matrix` (0 = server default, 1=degmatch, 2=dnafull, 3=nuc44). SearchResponse echoes `t`. KmerGroupInfo includes `t` and `template_type`. InfoResponse contains per-database metadata (`DatabaseInfo` with name, default_k, max_mode, kmer groups; `VolumeInfo` with volume_index, num_sequences, total_postings, total_bases, db) plus server-level `max_queue_size`/`queue_depth`/`max_seqs_per_req`. `ikafssnclient` uses `max_seqs_per_req` to split queries into batches before sending, avoiding oversized requests that would be partially rejected. `info_format.hpp/cpp` provides shared validation (`validate_info()` — checks db/k/mode/t/template_type against InfoResponse, returns error string with capability listing) and formatting (`format_server_info()`, `format_all_databases()`) used by ikafssnclient (pre-flight validation), ikafssnhttpd (merged-info validation), and ikafssninfo (remote display)
- **`io/`** — BLAST DB reader (CSeqDB wrapper), FASTA reader, mmap RAII wrapper, seqidlist reader (text/binary auto-detect), result writer/reader (including output format parsing/validation and SAM/BAM dispatch), `.kvx` manifest reader, volume discovery (`discover_volumes()`, `index_file_stem()`, `khx_path_for()` — shared by search/server/info/index), primer query parser (`primer_query.hpp/cpp`: `parse_primer_pairs()` parses primer pair FASTA (even number of records → fwd/rev pairs), generates concatenated query `fwd + N×insert_length + RC(rev)`, and computes per-primer k-mer position counts for threshold resolution — shared by ikafssnsearch and ikafssnclient)
- **`util/`** — CLI parser (supports multi-value options via `get_strings()` for repeatable flags like `-ix`), size string parser ("8G"), socket utilities, progress display, logger, common CLI init helpers (`check_version()` supporting `-version`/`--version` with build timestamp, `print_version_header()` for `-h` display, `format_build_timestamp()`, `make_logger()`, `resolve_threads()`, `default_memory_limit()`, `format_size()` in `common_init.hpp`), context parameter parser (`context_parser.hpp`)

### Index file formats

Per BLAST DB volume, index files are generated with naming pattern `<vol_basename>.<kk>mer.{kix,kpx,ksx}` where `<vol_basename>` is the BLAST DB volume's own basename (from `FindVolumePaths()`) and `<kk>` is the zero-padded 2-digit k value (e.g. `nt.00.09mer.kix`, `nt.01.11mer.kpx` for standard multi-volume; `foo.09mer.kix`, `bar.09mer.kix` for aggregated DBs). When spaced seeds are enabled (`-t > 0`), the naming extends to `<vol_basename>.<kk>mer.<tt>mer.<type>.{kix,kpx,ksx}` where `<tt>` is the zero-padded template length and `<type>` is `cod` (coding) or `opt` (optimal). When built with `-mode 1`, `.kpx` files are omitted:

- **`.kvx`** — Text manifest file for volume discovery (naming: `<db_base>.<kk>mer.kvx`). Contains a `DBLIST` line listing volume basenames in order. Always generated by `ikafssnindex`. Readers use `.kvx` to discover volumes
- **`.kix`** — Header + direct-address table (offsets[4^k + 1] with sentinel, uint32 or uint64 depending on posting data size) + delta-compressed ID postings. No counts table (format v3); posting count derived from offset differences
- **`.kpx`** — Header + pos_offsets[4^k] (uint32 or uint64 depending on posting data size) + delta-compressed position postings (correlated with .kix ID postings). Not generated when `ikafssnindex -mode 1` is used; `discover_volumes()` sets `has_kpx=false` for such indexes, and search/server/info commands handle the absence gracefully
- **`.ksx`** — Header + seq_lengths[] + accession string table (enables standalone result display without BLAST DB)
- **`.khx`** — Header (32B, magic "KMHX") + bitset (ceil(4^k / 8) bytes). Shared across all volumes (one per k value, naming: `<db_base>.<kk>mer.khx`). Generated only when `-max_freq_build` is set to a value other than 1 (disabled). K-mer counts are aggregated across all volumes before applying the threshold. Records which k-mers were excluded during index build (bit i = 1 means k-mer i was excluded). Used at search time for fractional `-stage1_min_score` threshold adjustment

ID and position postings are in separate files so Stage 1 never touches `.kpx`, maximizing page cache efficiency.

### Key design conventions

- **C++17**, CMake >= 3.16, little-endian only (Linux x86_64/aarch64)
- k-mer values use direct-address tables (array index = k-mer integer value), no hash maps
- Template dispatch on `KmerInt` type; CLI-side dispatch uses `kmer_type_for(k, t)` (bit-width-based: `2*k`; >16 bits → `uint32_t`); reader-side dispatch reads `kmer_type` from index header. No virtual calls in hot loops
- Delta encoding uses LEB128 varint; position deltas reset at sequence boundaries (detected by checking if the corresponding ID delta is non-zero)
- N-containing k-mers are skipped via an N-counter (shared logic between indexing and query scanning)
- Reverse complement: index stores forward strand only; search generates both forward and revcomp of query
- Parallelization via Intel TBB: indexing uses `parallel_for` + `combinable` for counting and partition scan, `parallel_sort` for posting buffer, `task_group` for volume-level parallelism; `ikafssnsearch` adaptively selects `parallel_for` (query-level, when queries > threads*2 or single volume) or `parallel_for_each` (query×volume, for few queries with multiple volumes), and uses `tbb::parallel_sort` for final cross-volume result sorting; `ikafssnserver` always uses `parallel_for` (query-level) to preserve CPU headroom, and post-processes per-query sort/truncate with `parallel_for` + `nth_element`. Stage 1 and volume-level sort+truncate use `nth_element` to select top-K before sorting, avoiding full O(N log N) sort. Thread count is controlled centrally via `tbb::global_control` in each main binary
- mmap is read-only and shared across threads; `Stage1Buffer` (AoS entries) is thread-local. `ikafssnsearch` and `ikafssnserver` apply budget-based madvise hints via `-memory_limit` (default: half of RAM): MADV_WILLNEED is allocated in priority order (khx → kix dict → kpx dict → ksx) until the budget is exhausted, then MADV_RANDOM for the remainder
- One server process can serve multiple BLAST DBs simultaneously via repeatable `-ix`/`-db` flags; each DB is identified by its basename (from `parse_index_prefix()`), stored in `DatabaseEntry`, and looked up by `db` in search requests. Per-DB state includes `kmer_groups` (vector of `KmerGroup`, keyed by (k, t, template_type) via `find_group()`), `default_k`, `default_t`, `default_template_type`, `max_mode`, resolved `SearchConfig` (with per-DB `max_freq` resolution), `Stage3Config`, and context settings. Spaced seed search resolves masks via `get_seed_masks()`. "Both" search merges separate coding and optimal indexes at search time
- `ikafssnhttpd` supports multiple backends via repeatable `-server_socket`/`-server_tcp` flags (CLI order = priority). `BackendManager` (`src/ikafssnhttpd/backend_manager.hpp/cpp`) handles backend lifecycle: init with cross-server DB validation (per shared (db, k, t, template_type) tuple: total sequences and total bases must match; k-value sets may differ across backends), merged-info aggregation (union of k-value groups across backends, max_mode = max across backends), priority+capacity routing (considering both slot availability and `max_seqs_per_req`), exclusion on failure (`-exclusion_time`), and periodic heartbeat (`-heartbeat_interval`). `HttpController` uses `BackendManager` for merged-info validation, routed search, and aggregated info JSON (with per-mode capacity and `max_seqs_per_req` in `modes` array, plus top-level `max_seqs_per_req`)

## Documentation

Usage documentation is maintained in two languages:
- `doc/ikafssn.en.md` — English
- `doc/ikafssn.ja.md` — Japanese

## Test structure

Tests are in `test/` using CTest. Real SSU_eukaryote_rRNA BLAST DB at `db/SSU_eukaryote_rRNA` is used for all BLAST-DB-dependent tests. `test/scripts/setup_ssu_testdata.sh` generates derived test data (ambig DB, multi-volume DBs, queries) in `/tmp/ikafssn_ssu_test/`. Shared fixture: `test/ssu_test_fixture.hpp`.

## Plan Mode Rules

- プランモードには**自発的に入らない**。ユーザーからの明示的な指示があった場合のみ入ること。
- 通常のワークフロー:
  1. プランモード外でユーザーと相談しながら方針を練る。
  2. 方針が十分に固まった段階で、ユーザーが以下のいずれかで計画立案を指示する:
     - 「プランモードに入って上記方針に基づいてプランを立てて提示するように」等の明示的な指示
     - `/plan` でプランモードに入った後、「上記方針に基づいてプランを立てて提示するように」等の指示
  3. この明示的な指示を受けて初めて計画立案を開始する。
- 計画を提示する際は、プランファイル（`.md`）の**絶対パス**を表示すること。

## Version Management

- バージョン番号は `CMakeLists.txt` の `IKAFSSN_VERSION` で管理される。形式は `"0.1.YYYY.MM.DD"`（本日の日付）。
- **コミット前に必ず確認**: `IKAFSSN_VERSION` の日付部分が本日の日付と一致しているか確認し、古い場合は本日の日付に更新してからコミットすること。ユーザーからの指示がなくても自主的に行う。

## Development Environment Rules

- `sudo` が必要なコマンドは Claude Code から直接実行しない。ユーザーに提示して手動実行を求めること。
