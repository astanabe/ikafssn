# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **Important**: Always use Japanese for on-screen communication with the user. However, write files in English unless they are designated as Japanese files. The only Japanese file in this repository is `doc/ikafssn.ja.md`. This file (`CLAUDE.md`) is therefore English-only — including any new sections, examples, and example values inside tables.

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
| NCBI C++ Toolkit 30.2.0 | 30.2.0 | `./ncbi-toolkit/` | source build |
| Intel TBB (oneTBB) | 2021.11.0 | system | `sudo apt install libtbb-dev` |
| Drogon | 1.9.12 | `./drogon/` | source build (static; trantor 1.5.26 placed at `drogon-1.9.12/trantor`) |
| libcurl | 8.5.0 | system | apt (pre-installed) |
| Parasail | 2.6.2 | `./parasail/` | source build |
| htslib | 1.23.1 | `./htslib/` | source build |
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
| `ikafssnindex` | Build k-mer inverted index from BLAST DB (`-mode 1` skips `.kpx`; `-template_type both` builds coding and optimal indexes sequentially; `-freq_threshold_part` controls .kpx v9 per-(kmer, seq_id) partitioning, default 8, max 255) | NCBI C++ Toolkit, TBB |
| `ikafssnsearch` | Local direct search (mmap index), in-silico PCR primer mode | NCBI C++ Toolkit, TBB, Parasail (mode 3), htslib (SAM/BAM) |
| `ikafssnretrieve` | Extract matched subsequences | NCBI C++ Toolkit, libcurl (remote) |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) | NCBI C++ Toolkit, TBB |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver(s) (multi-backend) | Drogon |
| `ikafssnclient` | Client (socket or HTTP), in-silico PCR primer mode | libcurl (HTTP mode) |
| `ikafssninfo` | Index/DB info display (local or remote) | NCBI C++ Toolkit (local), socket/HTTP (remote) |

### Shared library layers (`src/`)

- **`core/`** — Fundamental types (`Hit`, `ChainResult`), constants, k-mer 2-bit encoding/revcomp (templates parameterized on `KmerInt` = `uint16_t` or `uint32_t`, selected by `kmer_type_for(k, t)` based on bit width: `2*k` bits), `AmbigInfo` struct + `expand_ambig_kmer_multi()` for configurable multi-position degenerate expansion (`max_expansion` parameter controls product limit), LEB128 varint, spaced seed (discontiguous megablast) support (`spaced_seed.hpp`: `TemplateType` enum, 24 template masks for k=8/9 × t=13/15/18 (derived from discontiguous MegaBLAST template design principles) and k=11/12 × t=16/18/21 (MegaBLAST-native) × coding/optimal, `get_seed_masks()`, `reverse_complement_string()` for string-level RC used by spaced seed query scanning; accumulator-based extraction utilities — `ExtractionRun`/`SpacedSeedExtractor` structs decompose each mask into contiguous-1-bit runs, 24 `constexpr` extractors pre-computed at compile time, `extract_kmer_ct<Ext, KmerInt>()` NTTP function fully unrolled by the compiler with all shifts/masks as immediates at `-O2 -funroll-loops`, `extract_for_mask<KmerInt>()` switch dispatch selects the NTTP specialization per mask at ~0 cost via branch prediction, `extractor_for_mask()` returns the constexpr extractor for slow-path `kmer_bit_offset` access). `PackedKmerScanner::scan_spaced()` and `KmerScanner::scan_spaced()`/`scan_spaced_ambig()` use a sliding `uint64_t` accumulator (shift+OR per position, shared across masks) + window bitsets (`ambig_bits`/`n_bits`/`degen_bits`) for O(1) ambiguity checking, replacing the previous O(t) per-position k-mer rebuild and O(|ambig_entries|) linear ambiguity search. "Both" search mode merges separate coding and optimal indexes at search time rather than using a single combined index
- **`index/`** — Reader/writer for four index file formats (`.kix` main index, `.kpx` position data, `.ksx` sequence metadata, `.khx` build-time exclusion bitset).  All four files share format version 9 (Phase 7 — single major version policy, all index files bump together; users only have to track one number).  `.kix` v9 (magic `"KIX9"`, 96 B header) stores the **distinct seq_id** delta stream (codec body wire-compatible with v8) — intra-sequence k-mer duplicates are removed at build time by a SIMD dedup kernel (`src/index/seq_id_dedup.{hpp,cpp}` runtime dispatcher + `seq_id_dedup_tier.cpp` per-tier implementation, mirroring the FastPFor tier dispatch).  `-max_freq_build` thresholds by **the number of containing sequences**, not total occurrences.  `.kpx` v9 (magic `"KPX9"`, 64 B header) keeps the per-(kmer, seq_id) partition grouping but classifies every distinct seq_id via a **2-bit kind map** indexed against the .kix decoded distinct seq_id array (00 = short_occ1, 01 = short_occ_ge2, 10 = partition, 11 = reserved); the seq_id is therefore not stored in the .kpx posting list at all and the short bucket is split into occ=1 (positions only) and occ>=2 (u8 occ_count[] + positions) sub-buckets.  Decoding is candidate-set-driven (`open_stream_kpx_for_candidates(posting, bytes, kix_decoded, kix_count, candidates, n_candidates, scratch, out)` performs a 2-pointer merge walk of the .kix decoded array against the (sorted) candidate list, resolving each candidate's (kind, rank) pair which then indexes into per-kind decoded position buffers held in the caller-supplied `PosDecodeScratch`).  Both `.kix` and `.kpx` posting list bodies still go through `src/index/pfd_codec.{hpp,cpp}` (runtime dispatcher) + `pfd_codec_tier.cpp` (per-tier implementation): `.kix` posting list body is FastPFor's `CompositeCodec<SIMDFastPFor<4>, VariableByte>` (PForDelta with VByte exception stream); each partition group, the short_occ1 sub-bucket and the short_occ_ge2 sub-bucket use the same FOR-within-block stream — per-128-element block stores its min and the bitpacked spread `(value - min)` with an 8 B aligned header `[u32 min][u8 b][3 B pad]` followed by 16*b body bytes (Phase 6 proposal F), and the stream tail switches to a packed bit-width form `[u8 tail_count][u32 tail_min][u8 tail_b][bitpacked body]` (Phase 6 proposal D) replacing the prior varint stream.  Bitpacking primitives come from FastPFor v0.4.0 (`simdpackwithoutmask` / `simdunpack`), recompiled into per-ISA OBJECT libraries under different `FastPForLib_<tier>` namespace remappings.  On x86_64 the four tiers `ikafssn_pfd_{sse42,avx2,avx512bw,avx512vbmi2}` are built alongside the parallel Phase 7b EF tier ladder `ikafssn_ef_{sse42,avx2,avx512bw,avx512vbmi2}` (SSE4.2 is the floor, Scalar is rejected at startup with `exit(2)`).  On aarch64 single `ikafssn_pfd_neon` and `ikafssn_ef_neon` tiers are built — `pfd_codec_tier.cpp` compiles via SIMDe (https://github.com/simd-everywhere/simde) translating its SSE intrinsics to NEON, while `ef_codec_tier.cpp` is intrinsics-free and compiles directly; SVE / SVE2 capable CPUs are routed to the same NEON tier objects (Phase 5h, single-tier design — a SVE-native bitpacker would be orthogonal hand-coded work and is not implemented).  The EF dispatcher (`ef_codec.cpp`) routes inner-word select1 through BMI2 `PDEP+TZCNT` on AVX2-and-above tiers and falls back to `__builtin_popcountll` + bit-stripping `__builtin_ctzll` on SSE4.2 / NEON.  The runtime dispatchers in `pfd_codec.cpp` and `ef_codec.cpp` select the per-arch tier at first use based on `current_simd_cap()`.  Codec selection follows `format_version`; v8-or-older indexes are rejected at open with a clear "rebuild your index" message.  Phase 7a–7e replace the raw u32/u64 dictionaries with Elias-Fano blobs: `.kix` keeps sentinel-based dictionary entries (`offsets[table_size+1]`) and `.kpx` is non-sentinel-terminated (`pos_offsets[table_size]`), both encoded by `encode_dictionary_ef(offsets, D, U_raw, out)` in `dictionary_io.{hpp,cpp}` (writers) / `EFDictionary::open()` (readers).  Phase 7c+7d further deduplicate the per-posting-list headers: `.kix` is just `[u32 distinct_count][body]` (4 B; `body_words` derived from EF dictionary's `posting_byte_length`); `.kpx` is headerless — body starts at the 2-bit kind map, with distinct/partition/short1/short2 counts derived from `kix_count` + a SIMD popcount on the kind map, and `short2_position_count` derived from the running cum sum that builds the short2 offset table.  The `KIX_FLAG_OFFSET32` flag and `KpxHeader.offset_type` byte are reserved (writers force-clear / set the EF sentinel `0xFF`; readers ignore them).  `KixReader`/`KpxReader` hold an `ef::EFDictionary` member; `posting_list_offset(i)` and `pos_offset(i)` route through `dict_.access(i)`.  `KixReader::bulk_count_postings()` runs in O(table_size), walking the EF dictionary with sequential `access_pair()` calls and reading the O(1) leading `distinct_count` field of every posting list (used at index-build / -info time only — the search hot path no longer counts posting lists per query).  Table index types use `uint32_t` (max table size = 4^12 = 16M).  Index builder uses a partition+buffer strategy (`skip_kpx` config flag omits `.kpx` for mode 1 indexes); builder rewrites files at finalize to attach the EF dictionary blob.  For spaced seeds (`-t > 0`), builder uses `scan_spaced()`; `-template_type both` builds coding and optimal indexes sequentially (one after the other).  All readers expose `willneed_size()` / `apply_madvise(bool)` for budget-based madvise allocation (used by `-memory_limit` in search/server).  v8-or-older readers/writers are fully removed; users upgrading from v8 must reindex
- **`search/`** — Three-stage search pipeline: Stage 1 (ID posting list scan with OID filter, coverscore or matchscore) -> Stage 2 (position-aware chaining DP with diagonal filter, chainscore) -> Stage 3 (Parasail pairwise alignment with BLAST DB subject retrieval, alnscore/ppositive/CIGAR). Mode 1 skips Stage 2+3, mode 2 skips Stage 3. The sort key is auto-determined by mode (1=stage1, 2=chainscore, 3=alnscore). Defaults: stage1_topn=0 (unlimited, no sort), num_results=0 (unlimited, no sort), stage1_min_score=0.5 (fractional), stage2_min_score=1 — speed-first defaults that skip sorting. Set positive stage1_topn/num_results to enable sorting. Fractional stage1_min_score (0 < P < 1) resolves per-query threshold as `ceil(Nqkmer * P) - Nhighfreq` where `Nqkmer = max(0, seq_len - span + 1)` is the pure window count (`span = t` for spaced seeds, `k` otherwise; content-independent) and `Nhighfreq = #{positions with ANY-excluded k-mer via .khx} + (Nqkmer − #emitted positions)` (the second term covers windows that fail to emit a k-mer due to degenerate over-expansion or truly-invalid characters). When the resolved threshold falls below 1 on every searched strand the query is skipped with `kSkipThresholdUnreachable` rather than silently returning zero hits. Other skip reasons (`kSkipQueryTooShort` for `seq_len < span`, `kSkipDegenRejected` when `accept_qdegen=0` and the query has IUPAC degenerates, `kSkipInvalidChar` for non-IUPAC characters) are detected inside `preprocess_query` and propagated via `QueryKmerData.skip_reason`/`skip_detail`; downstream binaries emit a sentinel record (`*SKIPPED:<reason>`) in TSV / JSON (`"status": "skipped"`) / SAM (unmapped + `XR:Z`/`XD:Z` aux tags), and `result_reader` silently drops these rows for backward compatibility. `QueryKmerData` uses SoA layout (`fwd_positions`/`fwd_kmer_values` separate vectors) to eliminate pair padding overhead. `stage1_filter` takes SoA pointers `(const uint32_t* positions, const KmerInt* kmers, size_t n)`. Stage 1 uses byte-limit decoding (`SeqIdDecoder(data, end)` + `has_more()`); the `.kix` codec (`CompositeCodec<SIMDFastPFor<4>, VariableByte>` + delta cumsum) is hidden behind that decoder API so Stage 1 callers do not need to know about PForDelta blocks.  Stage 2 re-decodes the .kix posting list per k-mer with `SeqIdDecoder` and feeds `decoded_data()` / `decoded_count()` into `PosDecoder` so the .kpx 2-pointer merge walk can resolve each candidate's (kind, rank) pair.  A thread-local `pfd::PosDecodeScratch` holds the per-kind position buffers across calls so the search hot path avoids per-k-mer allocation. `Stage1Buffer` uses AoS `Stage1Entry<Tier>` with 3-tier template specialization (`Stage1Tier::T8` = 2B, `T16` = 4B, `T32` = 8B per entry); tier is selected from actual preprocessed k-mer position counts via `select_tier()`. Query k-mers with IUPAC degenerate bases are expanded when the product of per-position variant counts <= `max_degen_expand` (configurable, default 16 for search/server, 4 for index); Stage 1 uses per-position dedup (`last_pos` in AoS entry) to avoid inflating scores from expanded k-mers; Stage 2 deduplicates `(q_pos, s_pos)` hits before diagonal filter and uses `span` (= t for spaced seeds, k for contiguous) for chain endpoint coordinates; `max_nhit_per_subject` (default 1) controls how many non-overlapping chains are extracted per subject via greedy best-chain removal; Stage 3 uses parasail_sg semi-global alignment with configurable score matrix (DEGMATCH default, also dnafull/nuc44 via `-stage3_score_matrix`), optional traceback for CIGAR/ppositive/qseq/sseq, context extension for subject subsequence retrieval (with overlap resolution for multi-chain hits when context > 0), and volume-parallel BLAST DB prefetch. CIGAR match (`=`) / mismatch (`X`) is determined by score matrix (score > 0 = match). `OutputHit.oid`/`OutputHit.volume` carry the BLAST DB OID and volume index from Stage 2, enabling Stage 3 to directly open the correct volume reader and fetch by OID without building a full accession-to-OID map. For spaced seeds, query preprocessing uses string-level reverse complement instead of `kmer_revcomp()`, and `SearchConfig.t` propagates the template length to all pipeline stages. "Both" template_type search uses `search_volume_both()` to merge results from separate coding and optimal indexes at search time; Stage 1 walks both streams in `q_pos` order through a *single shared* `Stage1Buffer` so per-(sid, q_pos) dedup carries across templates (merged Stage 1 score stays in `[0, Nqkmer]` instead of `[0, 2×Nqkmer]`), and the unified threshold is `min(thr_cod, thr_opt)` so hybrid queries (one region cod-only, another opt-only) are not lost. `stage1_filter_accumulate` / `stage1_filter_finish` are the cross-template primitives; the existing `stage1_filter` is now a wrapper that calls accumulate-then-finish on a fresh-cleared buf. In-silico PCR primer mode (`-primer`) is supported by both `ikafssnsearch` and `ikafssnclient`; see `io/primer_query.hpp/cpp` for the shared primer pair parsing logic
- **`protocol/`** — Length-prefixed binary protocol (v10) for client-server communication. Frame header (12B, magic + payload_size + msg_type + msg_version=10 + reserved) + typed messages. v10 (Phase 1) replaces the legacy `QueryResult.skipped` (uint8) with `QueryResult.skip_reason` (uint8 `SkipReason` enum) + `skip_detail` (str16) + `qlen` (u32); v9 wire format is rejected by `read_frame()`. SearchRequest/Response include `db` for multi-DB routing. SearchRequest includes `max_degen_expand`, `stage2_max_nhit_per_subject`, `t`, `template_type`, and `score_matrix` (0 = server default, 1=degmatch, 2=dnafull, 3=nuc44). SearchResponse echoes `t`. KmerGroupInfo includes `t` and `template_type`. InfoResponse contains per-database metadata (`DatabaseInfo` with name, default_k, max_mode, kmer groups; `VolumeInfo` with volume_index, num_sequences, total_postings, total_bases, db) plus server-level `max_queue_size`/`queue_depth`/`max_seqs_per_req`. `ikafssnclient` uses `max_seqs_per_req` to split queries into batches before sending, avoiding oversized requests that would be partially rejected. `info_format.hpp/cpp` provides shared validation (`validate_info()` — checks db/k/mode/t/template_type against InfoResponse, returns error string with capability listing) and formatting (`format_server_info()`, `format_all_databases()`) used by ikafssnclient (pre-flight validation), ikafssnhttpd (merged-info validation), and ikafssninfo (remote display)
- **`io/`** — BLAST DB reader (CSeqDB wrapper), FASTA reader, mmap RAII wrapper, seqidlist reader (text/binary auto-detect), result writer/reader (including output format parsing/validation and SAM/BAM dispatch), `.kvx` manifest reader, volume discovery (`discover_volumes()`, `index_file_stem()`, `khx_path_for()` — shared by search/server/info/index), primer query parser (`primer_query.hpp/cpp`: `parse_primer_pairs()` parses primer pair FASTA (even number of records → fwd/rev pairs), generates concatenated query `fwd + N×insert_length + RC(rev)`, and computes per-primer k-mer position counts for threshold resolution — shared by ikafssnsearch and ikafssnclient)
- **`util/`** — CLI parser (supports multi-value options via `get_strings()` for repeatable flags like `-ix`), size string parser ("8G"), socket utilities, progress display, logger, common CLI init helpers (`check_version()` supporting `-version`/`--version` with build timestamp, `print_version_header()` for `-h` display, `format_build_timestamp()`, `make_logger()`, `resolve_threads()`, `default_memory_limit()`, `format_size()` in `common_init.hpp`), context parameter parser (`context_parser.hpp`)

### Index file formats

Per BLAST DB volume, index files are generated with naming pattern `<vol_basename>.<kk>mer.{kix,kpx,ksx}` where `<vol_basename>` is the BLAST DB volume's own basename (from `FindVolumePaths()`) and `<kk>` is the zero-padded 2-digit k value (e.g. `nt.00.09mer.kix`, `nt.01.11mer.kpx` for standard multi-volume; `foo.09mer.kix`, `bar.09mer.kix` for aggregated DBs). When spaced seeds are enabled (`-t > 0`), the naming extends to `<vol_basename>.<kk>mer.<tt>mer.<type>.{kix,kpx,ksx}` where `<tt>` is the zero-padded template length and `<type>` is `cod` (coding) or `opt` (optimal). When built with `-mode 1`, `.kpx` files are omitted:

- **`.kvx`** — Text manifest file for volume discovery (naming: `<db_base>.<kk>mer.kvx`). Header lines include `FORMAT_VERSION 9` (informational; the parser ignores unknown lines for forward compatibility) followed by `TITLE` and `DBLIST "vol0" "vol1" ...`.  Always generated by `ikafssnindex`.  Readers use `.kvx` to discover volumes
- **`.kix`** — Magic `"KIX9"`, 96 B header, **Elias-Fano dictionary blob** (`offsets[4^k + 1]` sentinel-terminated; encoded by `encode_dictionary_ef`).  Per-posting-list format: `[u32 distinct_count][u32 body[(bytes-4)/4]]` where the body is FastPFor's `CompositeCodec<SIMDFastPFor<4>, VariableByte>` output over the **distinct seq_id** delta stream `[abs_first, d1, d2, ...]` with `d_i >= 1` (cumsum at decode reconstructs absolute distinct seq_ids).  `body_words` is derived at decode time from the EF dictionary's posting_byte_length (Phase 7c dedup B; the in-blob u32 was redundant).  Intra-sequence k-mer duplicates are removed at build time by the SIMD dedup kernel (`seq_id_dedup_tier.cpp`).  PForDelta's exception stream handles the wide-value `abs_first` while the per-block bit-width tracks the typically-small consecutive deltas.  `distinct_count` is available in O(1) from the leading u32 of each posting list.  Header carries `total_distinct_postings` (sum across k-mers).  `-max_freq_build` thresholds by `distinct_count` (i.e. the number of sequences containing the k-mer).  The .kix decoded distinct seq_id array is the resolution table for the .kpx 2-bit kind map; the .kpx decoder requires the .kix decoded array as input
- **`.kpx`** — Magic `"KPX9"`, 64 B header (Phase 6 layout: per-(kmer, seq_id) partition groups + 2-bit kind map + occ=1/occ>=2 split short buckets), **Elias-Fano dictionary blob** `pos_offsets[4^k]` (no sentinel; encoded by `encode_dictionary_ef`).  Per-posting-list format (Phase 7c+7d, no fixed header — all four redundant u32 fields removed): `[2-bit kind map: ceil(distinct_count*2/8) bytes]` (00 = short_occ1, 01 = short_occ_ge2, 10 = partition, 11 = reserved — indexed by position in the .kix decoded distinct seq_id array, so the seq_id is *not* stored in .kpx); then `partition_count` partition groups in .kix sid order each `[u32 occ_count][FOR-block stream]`; then a FOR-block stream over the `short1_count` occ=1 positions; then `[u8 occ_count[short2_count]][FOR-block stream over short2_position_count positions]` for the occ>=2 sub-bucket.  `distinct_count` comes from `kix_count` (decoder parameter); `partition_count`/`short1_count`/`short2_count` come from `popcount_kinds(kind_map)`; `short2_position_count` is the cum sum that builds the short2 offset table.  Empty `.kpx` posting lists emit zero bytes.  Each FOR-block is 8 + 16*b bytes: `[u32 min][u8 b][3 B pad][16*b bytes bitpacked (value - min)]` (proposal F: 8 B aligned body within the block).  Stream tail layout (proposal D): `[u8 tail_count]` plus when `tail_count > 0` a `[u32 tail_min][u8 tail_b][bitpacked body: ceil(tail_count*tail_b/8) B]` packed bit-width stream — replaces the v7 varint tail.  A (k-mer, seq_id) cluster whose occurrence count exceeds the build-time `-freq_threshold_part` (default 8, `>` comparison; max 255 since occ_count is u8) gets its own partition group; clusters with occ=1 land in the short_occ1 sub-bucket (no occ_count bytes); the rest land in short_occ_ge2.  Decoding is candidate-set-driven and requires the **.kix decoded distinct seq_id array** as input alongside the (sorted) candidate set (`pfd::open_stream_kpx_for_candidates(posting, bytes, kix_decoded, kix_count, candidates, n_candidates, scratch, out)` returns per-candidate position vectors).  The decoder runs a 2-pointer merge walk between kix_decoded and candidates and resolves each candidate's (kind, rank) pair via the kind map; `PosDecodeScratch` is a per-thread reusable buffer that holds the per-kind decoded position arrays so the search hot path avoids per-k-mer allocation.  Header carries `total_position_count` (sum across k-mers).  Not generated when `ikafssnindex -mode 1` is used; `discover_volumes()` sets `has_kpx=false` for such indexes, and search/server/info commands handle the absence gracefully
- **`.ksx`** — Magic `"KMSX"`, format_version 9 (data layout unchanged from v3).  Header + seq_lengths[] + accession string table (enables standalone result display without BLAST DB).  Multi-defline OIDs (BLAST DBs built with `makeblastdb -parse_seqids` carrying multiple Seq_ids per OID — typically used to register identical sequences under several accessions) store every accession in a single `'\x01'`-joined string per OID.  `BlastDbReader::get_accession()` and `KsxReader::accession()` return that joined form verbatim; output writers (TSV / SAM / FASTA defline / protocol `sseqid`) emit it as-is so consumers can split on `'\x01'` themselves.  In-process callers that need individual tokens (e.g. `OidFilter` reverse map, `local_retriever` accession→OID map) use `split_accessions()` from `io/accession_utils.hpp`
- **`.khx`** — Magic `"KMHX"`, format_version 9 (data layout unchanged from v3).  32 B header + bitset (ceil(4^k / 8) bytes).  Shared across all volumes (one per k value, naming: `<db_base>.<kk>mer.khx`).  Generated only when `-max_freq_build` is set to a value other than 1 (disabled).  K-mer **distinct seq_id counts** are aggregated across all volumes before applying the threshold (aligned with .kix's distinct-count semantic).  Records which k-mers were excluded during index build (bit i = 1 means k-mer i was excluded).  Used at search time for fractional `-stage1_min_score` threshold adjustment

ID and position posting lists are stored in separate files so Stage 1 never touches `.kpx`, maximizing page cache efficiency.

### Key design conventions

- **C++20**, CMake >= 3.16, little-endian only (Linux x86_64/aarch64). Release builds use `-O2 -DNDEBUG -funroll-loops` (set via `CMAKE_CXX_FLAGS_RELEASE CACHE FORCE` in CMakeLists.txt; `-funroll-loops` ensures C++20 NTTP extraction loops are fully unrolled with immediate constants, matching hand-coded BLAST+-style performance without manual per-template functions)
- **x86_64 baseline is SSE4.2** (Nehalem 2008+, AMD Bulldozer 2011+); pre-SSE4.2 CPUs are rejected at startup with `exit(2)` from `init_simd_dispatch()`. The runtime SIMD ladder is 4 tier (SSE4.2 / AVX2 / AVX512BW / AVX512VBMI2); AVX-512 VBMI without VBMI2 (Ice Lake client) is silently demoted to AVX512BW. `IKAFSSN_FORCE_SIMD=avx512vbmi` is accepted as an alias for `avx512bw`. The `IKAFSSN_ENABLE_SIMD=OFF` build option keeps the legacy Scalar path for debug/portability builds
- **aarch64 baseline is NEON** (Armv8.0+, ASIMD); pre-NEON aarch64 CPUs are rejected at startup with `exit(2)` from `init_simd_dispatch()`.  SVE / SVE2 capable CPUs run the same `ikafssn_pfd_neon` PForDelta tier object (FastPFor side).  All per-kernel SIMD files (`text_simd`, `ncbi2na_unpack`, `kmer_revcomp_simd`, `degenerate_scan_simd`, `spaced_seed_simd`) ship a NEON dispatch on the aarch64 baseline; `text_simd` and `ncbi2na_unpack` additionally have separate SVE / SVE2 paths, and `kmer_revcomp_simd` carries SVE / SVE2 dispatch slots that currently fall back to scalar pending hand-coded kernels.  Cross-architecture `IKAFSSN_FORCE_SIMD` requests (e.g. `=sse42` on aarch64) are detected and collapsed to Scalar (= startup reject or test SKIP)
- k-mer values use direct-address tables (array index = k-mer integer value), no hash maps
- Template dispatch on `KmerInt` type; CLI-side dispatch uses `kmer_type_for(k, t)` (bit-width-based: `2*k`; >16 bits → `uint32_t`); reader-side dispatch reads `kmer_type` from index header. No virtual calls in hot loops
- ID posting lists are **distinct seq_id** delta-encoded (first absolute, then strictly-positive differences; intra-sequence duplicates removed at build time by the SIMD dedup kernel) and stored via FastPFor's `CompositeCodec<SIMDFastPFor<4>, VariableByte>` (PForDelta with VByte exception stream); position posting lists are stored as absolute values, partitioned per-(kmer, seq_id) above the build-time `freq_threshold_part` (default 8, max 255) into independent FOR-within-block streams plus an occ=1 / occ>=2 split short bucket classified by a 2-bit kind map (Phase 6), so each partition group sees only intra-sequence positions and absolute-magnitude effects no longer leak into other clusters' bit-widths.  `-max_freq_build` thresholds by **distinct seq_id count** (number of containing sequences), matching the original design intent.  Search-time high-frequency k-mer filtering has been removed — only build-time exclusion via `.khx` is in effect
- N-containing k-mers are skipped via an N-counter (shared logic between indexing and query scanning)
- Reverse complement: index stores forward strand only; search generates both forward and revcomp of query
- Parallelization via Intel TBB: indexing uses `parallel_for` + `combinable` for counting and partition scan, `parallel_sort` for the posting list buffer, `task_group` for volume-level parallelism; `ikafssnsearch` adaptively selects `parallel_for` (query-level, when queries > threads*2 or single volume) or `parallel_for_each` (query×volume, for few queries with multiple volumes), and uses `tbb::parallel_sort` for final cross-volume result sorting; `ikafssnserver` always uses `parallel_for` (query-level) to preserve CPU headroom, and post-processes per-query sort/truncate with `parallel_for` + `nth_element`. Stage 1 and volume-level sort+truncate use `nth_element` to select top-K before sorting, avoiding full O(N log N) sort. Thread count is controlled centrally via `tbb::global_control` in each main binary
- mmap is read-only and shared across threads; `Stage1Buffer` (AoS entries) is thread-local. `ikafssnsearch` and `ikafssnserver` apply budget-based madvise hints via `-memory_limit` (default: half of RAM): MADV_WILLNEED is allocated in priority order (khx → kix dict → kpx dict → ksx) until the budget is exhausted, then MADV_RANDOM for the remainder
- One server process can serve multiple BLAST DBs simultaneously via repeatable `-ix`/`-db` flags; each DB is identified by its basename (from `parse_index_prefix()`), stored in `DatabaseEntry`, and looked up by `db` in search requests. Per-DB state includes `kmer_groups` (vector of `KmerGroup`, keyed by (k, t, template_type) via `find_group()`), `default_k`, `default_t`, `default_template_type`, `max_mode`, resolved `SearchConfig`, `Stage3Config`, and context settings. Spaced seed search resolves masks via `get_seed_masks()`. "Both" search merges separate coding and optimal indexes at search time
- `ikafssnhttpd` supports multiple backends via repeatable `-server_socket`/`-server_tcp` flags (CLI order = priority). `BackendManager` (`src/ikafssnhttpd/backend_manager.hpp/cpp`) handles backend lifecycle: init with cross-server DB validation (per shared (db, k, t, template_type) tuple: total sequences and total bases must match; k-value sets may differ across backends), merged-info aggregation (union of k-value groups across backends, max_mode = max across backends), priority+capacity routing (considering both slot availability and `max_seqs_per_req`), exclusion on failure (`-exclusion_time`), and periodic heartbeat (`-heartbeat_interval`). `HttpController` uses `BackendManager` for merged-info validation, routed search, and aggregated info JSON (with per-mode capacity and `max_seqs_per_req` in `modes` array, plus top-level `max_seqs_per_req`)

## Terminology (inverted index terms)

This is a **permanent rule** that applies to every change in this repository, not just the change in front of you right now. When describing the inverted index — in any of: documentation, code identifiers, source comments, docstrings, log messages, error strings, CLI help text, wire-format / protocol field names, or conversation with the user — use the canonical full-text-search vocabulary defined below. There are **no exceptions**: code, tests, fixtures, tools, scripts, build files, and Markdown documents are all bound by the same rule.

### Canonical vocabulary

The canonical mapping (full-text-search domain) is:

| Concept | Canonical term | Definition |
|---|---|---|
| Indexed unit | **term** (interchangeable with **k-mer** in this codebase) | A single token being indexed. In ikafssn the term is a k-mer. |
| List of term IDs | **dictionary** | The list of term IDs. In ikafssn's direct-address layout, the dictionary array is indexed by the term ID (k-mer integer value) and stores the byte offset to that term's posting list inside the posting file. |
| Sequence ID where a term occurs | **posting** | A single entry inside a posting list — the ID of one sequence containing the term (in `.kpx`, also carries position information). |
| List of postings for one term | **posting list** | All postings for a given term, concatenated. |
| Concatenation of all posting lists | **posting file** | The bulk region holding every term's posting list, in dictionary order. In ikafssn this region lives inside `.kix` and `.kpx` after the header + dictionary. |

Auxiliary terms that follow from the canonical vocabulary:
- **dictionary entry** — one offset value in the dictionary array (`offsets[i]`).
- **posting list header** — the fixed-size metadata at the front of a posting list (e.g. v9 `.kix`: `[u32 distinct_count]` = 4 B; v9 `.kpx`: empty — body starts at the kind map).
- **posting list body** — the variable-size encoded content of a posting list after its header (FastPFor codec output, FOR-block stream, kind map, etc.).

A nucleotide sequence is **not** called a "document". Do not import that word from generic IR literature into this codebase. The term "sequence" (or `seq_id` for its identifier) stands on its own.

### Forbidden expressions (always rewrite)

- "offset table" → **dictionary**
- "posting data region" / "posting list region" / "posting blob" / "per-kmer payload" → **posting file**
- "per-kmer header" → **posting list header**
- "posting" used alone to mean the full per-k-mer data → **posting list**
- "payload" in the narrow sense (variable-size content of a posting list) → **posting list body**
- "document" (when the actual referent is a nucleotide sequence) → **sequence**

### Application — no exceptions

This convention is binding on **every** artefact in the repository:

- code identifiers — function names, method names, struct / class members, free variables, namespaces, file names, CMake target names;
- public API symbols and wire-format / protocol field names (rename in lockstep with the appropriate `format_version` / `msg_version` bump if the rename touches the wire format);
- log messages, error strings, and CLI help text;
- comments, docstrings, and any Markdown documentation;
- test names and fixture file names;
- conversation with the user.

When you encounter a legacy identifier or comment that uses a forbidden term, **rename it as part of the change you are already making**. Do not file the rename as a follow-up. Do not retain the legacy term "for compatibility" or "for diff readability" — those are not exceptions. The only legitimate reason to defer a rename is that the canonical term itself is genuinely ambiguous in context; in that case follow the next subsection.

### When the canonical term is ambiguous, ask

If you are about to introduce a new concept or rename an existing identifier and you find that:

- multiple terms in the canonical vocabulary could plausibly apply, or
- the concept does not fit any of the canonical terms above, or
- the established vocabulary in the surrounding code disagrees with the canonical mapping in a way that is not obvious to resolve,

**stop and ask the user before choosing a name.** Do not invent a new convention silently. Do not pick "the closest one" without checking. State the alternatives, the trade-off, and your tentative recommendation, then wait for the user's decision. The choice you make becomes part of this codebase forever, so a few seconds of clarification is always cheaper than a wrong commit.

## Documentation

Usage documentation is maintained in two languages:
- `doc/ikafssn.en.md` — English
- `doc/ikafssn.ja.md` — Japanese

## Test structure

Tests are in `test/` using CTest. Real SSU_eukaryote_rRNA BLAST DB at `db/SSU_eukaryote_rRNA` is used for all BLAST-DB-dependent tests. `test/scripts/setup_ssu_testdata.sh` generates derived test data (ambig DB, multi-volume DBs, queries) in `/tmp/ikafssn_ssu_test/`. Shared fixture: `test/ssu_test_fixture.hpp`.

## Plan Mode Rules

- **Do not enter plan mode autonomously.** Enter plan mode only when the user explicitly instructs you to do so.
- Standard workflow:
  1. Discuss the approach with the user outside of plan mode.
  2. Once the approach is well-defined, the user issues a planning instruction such as:
     - "Enter plan mode and prepare a plan based on the above approach", or
     - After entering plan mode via `/plan`, "Prepare a plan based on the above approach".
  3. Begin planning only after this explicit instruction is received.
- When presenting a plan, display the **absolute path** of the plan file (`.md`).

## Version Management

- The version number is managed via `IKAFSSN_VERSION` in `CMakeLists.txt`. Format: `"0.1.YYYY.MM.DD"` (today's date).
- **Verify before every commit**: confirm that the date portion of `IKAFSSN_VERSION` matches today's date. If outdated, update it to today's date before committing. Do this proactively, even without an explicit user instruction.

## Development Environment Rules

- Do not execute commands that require `sudo` directly from Claude Code. Present such commands to the user for manual execution.

## Git Push Policy

- **Direct push to the `main` branch is permitted.** This repository is a personal development repository with no PR-review protection configured. Claude Code's default safety guard (refusing direct push to `main`) does not apply here. When the user instructs "commit and push", execute `git push origin main` without further confirmation.
