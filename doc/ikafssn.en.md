# ikafssn User Guide

**ikafssn** (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) is a suite of tools that builds inverted indexes over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.

- **Primary repository**: <https://github.com/astanabe/ikafssn>
- **Secondary repository**: <https://gitlab.com/astanabe/ikafssn>

## Overview

ikafssn consists of seven independent command-line programs:

| Command | Purpose |
|---|---|
| `ikafssnindex` | Build k-mer inverted index from BLAST DB |
| `ikafssnsearch` | Local direct search (mmap index) |
| `ikafssnretrieve` | Extract matched subsequences |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver |
| `ikafssnclient` | Client (socket or HTTP) |
| `ikafssninfo` | Index/DB and server information display |

Each command is a standalone executable that links only its required dependencies.

## Quick Start

```bash
# 1. Build an index from a BLAST DB
ikafssnindex -db mydb -k 11 -o ./index

# 2. Search with a query FASTA
ikafssnsearch -ix ./index/mydb -query query.fasta

# 3. Extract matched subsequences
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db mydb > matches.fasta
```

## Commands

### Common options

All commands support the following options:

```
  -h, --help              Show usage information (with version header)
  -version, --version     Print version and build information
  -v, --verbose           Verbose output
```

The `-version` output format is:

```
<command>: <version>
 Package: ikafssn <version>, build <ISO 8601 timestamp>
```

### Compressed I/O

`ikafssnsearch`, `ikafssnclient`, and `ikafssnretrieve` read FASTA query
inputs and write TSV / JSON / FASTA output through transparent codec
wrappers:

- **Inputs** are auto-detected from the first few bytes of the source, so
  `.gz`, `.bz2`, `.xz`, and `.zst` FASTA files (and the same formats on
  stdin) decode without any user-side flag.  Concatenated streams
  (`cat a.fa.gz b.fa.gz > c.fa.gz`) decode correctly across all four
  codecs.
- **Outputs** select the codec from the trailing path suffix
  (`out.tsv.gz`, `out.json.zst`, `out.fa.bz2`, …).  Empty `-o` or
  `-o -` always writes uncompressed bytes to stdout regardless of any
  other setting; pipe through your own encoder if you need stdout
  compressed.
- `-compression_level <int>` controls the compression level for the
  selected codec.  Per-codec accepted ranges and defaults are: gzip 0..9
  (default 6), bzip2 1..9 (default 9), xz 0..9 (default 6), zstd
  `ZSTD_minCLevel()..ZSTD_maxCLevel()` (default 3).  The level is
  validated against the codec at parse time, before any I/O is opened.
- SAM and BAM **reject** the four compressed suffixes — BAM is already
  BGZF and we deliberately do not double-wrap.  Use `-output_format tsv` /
  `-output_format json` if you want a compressed result file.

### ikafssnindex

Build a k-mer inverted index from a BLAST database. For each volume, index files are generated: `.kix` (ID postings), `.kpx` (position postings, unless `-mode 1`), and `.ksx` (sequence metadata). When `-max_freq_build` is used, a shared `.khx` file (build-time exclusion bitset) is also generated. The `.khx` file is shared across all volumes (one per k value, not per volume).

**Indexed unit:** Each parent BLAST OID is split into one or more **fragments** at index time. A fragment is the unit registered as one internal SeqId in `.kix` / `.kpx` / `.ksx`; adjacent fragments of the same parent share `-overlap_length` bases so a chain that crosses the boundary still has at least one fragment that fully covers it. When `-min_length_split 0` is used, every parent has exactly one fragment that spans the whole parent (the degenerate / no-split layout). Search-side dedup folds fragment-relative coordinates back to **parent-relative coordinates** so downstream tools always see one canonical row per parent OID.

```
ikafssnindex [options]

Required:
  -db <path>              BLAST DB prefix
  -k <int>                K-mer length (5-15)
  -o <dir>                Output directory

Options:
  -mode <1|2|3>           Search mode the index will support (default: 1)
                          1 = Stage 1 only (skip .kpx generation, saves disk and time; default)
                          2 = Stage 1+2
                          3 = Stage 1+2+3 (same as 2 for index build)
  -min_seq_length <int>   Minimum length of a parent OID accepted into the
                          index (default: 64).  Parents shorter than this are
                          dropped at index time and never registered as
                          fragments.  Persisted in the .kix / .kpx / .ksx
                          headers; ikafssnsearch / ikafssnclient cross-check
                          their -min_query_length against this so a too-short
                          query cannot silently produce zero hits.
  -min_length_split <int> Fragment-splitting threshold (default: 50000 under
                          -mode 1, 0 under -mode 2/3).  A parent longer than
                          this is split into multiple overlapping fragments
                          using the kafssstore calcsegment2 formula
                          (DNA2 mode: ncbi4na==0xF / N runs cut a parent
                          into valid segments first, and each valid segment
                          is then carved into fragments of approximately
                          this size).  0 disables splitting (every parent
                          becomes a single fragment).
  -overlap_length <int>   Per-fragment trailing overlap, in bases (default:
                          500 under -mode 1, 0 otherwise).  Adjacent
                          fragments of the same valid segment share this
                          many bases.  Must satisfy
                          0 < overlap_length < min_length_split / 2 when
                          splitting is enabled.  This value also caps the
                          maximum query length at search time
                          (kSkipQueryTooLong is emitted for longer queries),
                          because the parent-relative dedup keys assume
                          every chain hit fits inside at most two adjacent
                          fragments.
  -memory_limit <size>    Per-volume sort buffer budget (default: half of physical RAM)
                          Accepts K, M, G suffixes.  Bounds peak RAM during
                          the postings pass by partitioning k-mer entries
                          to fit this budget per volume.
  -max_freq_build <num>   Exclude k-mers with cross-volume count above this threshold
                          1 or 1.0: disable (no exclusion, default)
                          0 < x < 1: fraction of total NSEQ across all volumes
                          > 1: absolute count threshold
                          0: not allowed (error)
                          Counts are aggregated across all volumes before filtering
  -freq_threshold_part <int>
                          .kpx per-(kmer, seq_id) partition threshold (default: 8, max: 255)
                          A (k-mer, seq_id) cluster with occurrence count > threshold
                          gets its own partition group; lower-multiplicity clusters
                          merge into a shared short bucket (improves chromosome-class
                          subject compression by decoupling per-block bit-width from
                          absolute position magnitude)
  -nthread_highfreq_filter <int>
                          Threads for cross-volume high-frequency filtering
                          (default: min(8, nthread))
  -max_degen_expand <int> Max degenerate expansion per k-mer (default: 4, max: 16, 0/1: disable)
                          Controls how many non-degenerate k-mers are generated from
                          a k-mer containing IUPAC degenerate bases. Expansion occurs
                          when the product of per-position variant counts <= this limit.
  -t <int>                Template length for spaced seeds (default: 0)
                          0: contiguous k-mers (traditional mode)
                          13, 15, 18: spaced seed template length (requires -k 8 or 9)
                          16, 18, 21: spaced seed template length (requires -k 11 or 12)
  -template_type <str>    Template type for spaced seeds (required when -t is specified)
                          coding: coding template only
                          optimal: optimal template only
                          both: build coding and optimal indexes sequentially
  -nthread <int>          Number of threads (default: all cores)
                          Parallelizes per-volume counting / partition
                          scan / sort and the per-OID metadata pass.
                          Volumes are processed one at a time; each
                          volume gets the full -nthread pool.
  -force_rebuild <0|1>    Force a full rebuild (default: 0)
                          0 = per-volume resume: existing .kix / .kpx /
                              .ksx are strictly validated (file size,
                              magic / format_version, full posting list
                              walk).  Volumes that pass are reused;
                              volumes that fail or are missing are
                              rebuilt.  .ksx.tmp / .kix.tmp / .kpx.tmp
                              from an interrupted previous run are also
                              picked up.  Post-build validation prints
                              any failing paths and exits non-zero
                              without deleting them.
                          1 = delete tmp / final files for this build's
                              parameters and rebuild every volume.
  -v, --verbose           Verbose output
```

**Examples:**

```bash
# Small DB, plenty of memory
ikafssnindex -db mydb -k 11 -o ./index -memory_limit 128G

# Large DB, limited memory, multi-threaded
ikafssnindex -db nt -k 11 -o ./nt_index -memory_limit 32G -nthread 32

# Exclude high-frequency k-mers during build (absolute)
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 50000

# Exclude k-mers appearing in >1% of total sequences across all volumes
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 0.01

# Build mode 1 index (Stage 1 only, no .kpx files)
ikafssnindex -db mydb -k 11 -o ./index -mode 1

# Build spaced seed index (k=11, t=16, coding template)
ikafssnindex -db mydb -k 11 -t 16 -template_type coding -o ./index

# Build spaced seed index with coding template only
ikafssnindex -db mydb -k 11 -t 18 -template_type coding -o ./index

# Build both coding and optimal indexes sequentially
ikafssnindex -db mydb -k 11 -t 18 -template_type both -o ./index

# Build spaced seed index with k=12, t=21
ikafssnindex -db mydb -k 12 -t 21 -o ./index

# Build spaced seed index optimized for PCR (k=8, t=13)
ikafssnindex -db mydb -k 8 -t 13 -o ./index
```

### ikafssnsearch

Local direct search command. Directly mmaps index files for searching. Does not require a running server.

```
ikafssnsearch [options]

Required:
  -ix <prefix>            Index prefix (like blastn -db)
  -query <path>           Query FASTA file (- for stdin)

Options:
  -k <int>                K-mer size to use (required if multiple k values exist)
  -o <path>               Output file (default: stdout)
  -nthread <int>          Threads for search and output formatting
                          (default: all cores)
  -memory_limit <size>    Search memory budget (default: half of RAM)
                          Pins .khx/.ksx metadata; residual caps the
                          Stage 3 posting heap.  Accepts K, M, G suffixes
  -mode <1|2|3>           Search mode (default: 1)
                          1=Stage 1 only, 2=Stage 1+2, 3=Stage 1+2+3
  -min_query_length <int> Minimum query length, in bases (default: 64).
                          Queries shorter than this are skipped with
                          kSkipQueryTooShort.  Must be >= the index's
                          min_seq_length (recorded in the .kix / .ksx
                          header); a smaller value would let queries match
                          parents the index never saw, so the integrity
                          check refuses to start.
  -db <path>              BLAST DB path for mode 3 (default: same as -ix)
  -stage1_max_nhit_per_subject <int>  Max Stage 1 candidates per parent
                          (default: 1, 0=unlimited); applied per
                          (query, parent, volume, strand)
  -stage1_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode
                          (default: 3); 1/2=top-N (no ties), 3/4=top-N + ties
  -stage1_max_nhit_per_volume <int>  Max Stage 1 candidates per
                          (query, volume, strand) (default: 0=unlimited)
  -stage1_max_nhit_in_total <int>  Max Stage 1 candidates per query across all
                          volumes / strands (default: 0=unlimited).  Setting
                          this also sets -stage1_max_nhit_per_volume to the same
                          value unless that is given explicitly (must be >= it).
  -stage1_min_score <num> Stage 1 minimum score (default: 0.5)
                          Integer (>= 1): absolute threshold
                          Fraction (0 < P < 1): proportion of query k-mers,
                            resolved per query as ceil(Nqkmer * P) - Nhighfreq
  -stage2_min_score <int> Minimum chain score (default: 0 = adaptive)
                          0 = use resolved Stage 1 threshold as minimum
                          >= 1: absolute minimum chain score
  -stage2_max_gap <int>   Chaining diagonal gap tolerance (default: 100)
  -stage2_max_lookback <int>  Chaining DP lookback window (default: 64, 0=unlimited)
  -stage2_max_nhit_per_subject <int>  Max chains per subject (default: 1, 0=unlimited)
  -stage2_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3)
                           1/2=top-N (no ties), 3/4=top-N + score ties;
                           1/3=strands merged per parent, 2/4=strands separate
  -stage2_max_nhit_in_total <int>  Max chains per query across all volumes,
                          by chainscore (default: 0=unlimited)
  -stage2_min_nhit_diag <int>  Diagonal filter min hits (default: 1)
  -context_extend <value>        Context extension for mode 3 (default: 2.0)
                          Integer: bases to extend; Decimal: query length multiplier
  -stage3_traceback <0|1> Enable traceback in mode 3 (default: 0)
  -stage3_gapopen <int>   Gap open penalty for mode 3 (default: 10)
  -stage3_gapext <int>    Gap extension penalty for mode 3 (default: 1)
  -stage3_min_ppositive <num>  Min percent positive filter for mode 3 (default: 0)
  -stage3_min_npositive <int>  Min positive-scoring positions filter for mode 3 (default: 0)
  -stage3_max_nhit_per_subject <int>  Max hits per subject, by alnscore (default: 1, 0=unlimited)
  -stage3_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3)
                           1/2=top-N (no ties), 3/4=top-N + score ties;
                           1/3=strands merged per parent, 2/4=strands separate
  -stage3_max_nhit_in_total <int>  Max hits per query across all volumes,
                          by alnscore (default: 0=unlimited)
  -stage3_score_matrix <name>  Score matrix: degmatch, dnafull, nuc44 (default: degmatch)
  -seqidlist <path>       Include only listed accessions
  -negative_seqidlist <path>  Exclude listed accessions
  -strand <-1|1|2>       Strand to search (default: 2)
                          1=plus only, -1=minus only, 2=both
  -accept_qdegen <0|1>    Accept queries with degenerate bases (default: 1)
  -max_degen_expand <int> Max degenerate expansion per k-mer (default: 16, max: 256, 0/1: disable).
                          This is a query parameter, not a filter: it also breaks ties when
                          several otherwise-identical index variants differ only in their
                          build-time max_degen_expand (the variant whose build value equals
                          this query value is preferred, else the largest is chosen).
  -t <int>                Template length for spaced seeds (default: 0 = contiguous)
                          0: contiguous k-mers (traditional mode); this is the default, NOT a wildcard
                          13, 15, 18: spaced seed template length (requires -k 8 or 9)
                          16, 18, 21: spaced seed template length (requires -k 11 or 12)
  -template_type <str>    Template type for spaced seeds (default: both; invalid with -t 0)
                          coding: use coding index only
                          optimal: use optimal index only
                          both: merge coding and optimal indexes at search time
  -min_seq_length <int>   Select the index variant built with this min_seq_length (default: any)
  -min_length_split <int> Select the index variant built with this min_length_split (default: any)
  -overlap_length <int>   Select the index variant built with this overlap_length (default: any)
  -max_freq_build <int>   Select the index variant built with this max_freq_build (default: any)
  -output_format <tsv|json|sam|bam>  Output format (default: tsv)
  -compression_level <int> Output compression level (gzip/bzip2/xz/zstd; default per codec: gzip=6, bzip2=9, xz=6, zstd=3)
                          Codec is selected by -o suffix (.gz/.bz2/.xz/.zst); SAM/BAM reject all four
  -v, --verbose           Verbose logging

Primer mode (alternative to -query):
  -primer <path>           Primer pair FASTA (even number of sequences; mutually exclusive with -query)
  -insert_length <int>     Expected insert length (required with -primer)
  -stage1_primer_score <num>  Stage 1 threshold for primer (0<v<=1: fraction, v>=2: absolute; default: 0.5)
  -stage2_primer_score_add <int>  Stage 2 threshold addon: max(Lf,Lr) + N (default: 1)
```

(Note: -query and -primer auto-detect gzip/bzip2/xz/zstd-compressed FASTA
 inputs from their leading magic bytes; no flag is required.)

The `-ix` option specifies the index prefix path (without extension), similar to `blastn -db`. For example, if `ikafssnindex -db nt -k 11 -o /data/index` generated the following files for a multi-volume BLAST DB (`nt` with volumes `nt.00`, `nt.01`):

```
/data/index/nt.00.k11.minlen64.minsplit50000.ovllen500.kix
/data/index/nt.00.k11.minlen64.minsplit50000.ovllen500.ksx
/data/index/nt.01.k11.minlen64.minsplit50000.ovllen500.kix
/data/index/nt.01.k11.minlen64.minsplit50000.ovllen500.ksx
/data/index/nt.k11.minlen64.minsplit50000.ovllen500.kvx
```

then specify `-ix /data/index/nt`. The prefix `/data/index/nt` is split into the directory `/data/index/` and the base name `nt`. Volumes are discovered via the `.kvx` manifest file, which lists the volume basenames. For aggregated databases (e.g. `combined` aggregating `foo` and `bar`), the index files would be `foo.k11.…kix`, `bar.k11.…kix`, with `combined.k11.…kvx` as the manifest.

**Index-variant selection.** An index variant is identified by the eight parameters encoded in its file name: `k`, `t`, `template_type`, `min_seq_length`, `min_length_split`, `overlap_length`, `max_freq_build`, and `max_degen_expand`. `ikafssnsearch` discovers every variant under the prefix and then narrows the set down with the options above: each unset option is a wildcard (matches any value), each given option is an exact filter. After filtering, the remaining variants must collapse to a single variant identified by the first seven parameters (`max_degen_expand` excepted) — otherwise the search aborts with an "index filter is ambiguous" error listing the options you can use to narrow it, or a "no index variant matches" error if nothing matched. So when only one variant exists, every option can be omitted; when several exist, supply just enough options to make the choice unique.

`max_degen_expand` is the one parameter allowed to remain non-unique: if several variants survive that differ only in their build-time `max_degen_expand`, the tie is broken by `-max_degen_expand` as described above. The build-time `-max_degen_expand` (from `ikafssnindex`) and the query-time `-max_degen_expand` need not match.

The template_type dimension resolves as follows:

| Options | Result |
|---|---|
| `-t 0` (or `-t` omitted) | the contiguous variant (single) |
| `-t>0 -template_type coding` | the coding variant (single) |
| `-t>0 -template_type optimal` | the optimal variant (single) |
| `-t>0 -template_type both` | the coding+optimal pair (both required; error if either is missing) |
| `-t>0` (no `-template_type`) | the coding+optimal pair if **both** exist, otherwise whichever single side exists (coding-only → coding, optimal-only → optimal) |

In the both-pair case, coding and optimal must agree on all identifying parameters except `template_type` — **including `max_degen_expand`** (no "chanpon"/mixing). The shared `max_degen_expand` is tie-broken among the values present on *both* sides (query value preferred, else the largest common value); if the two sides share no `max_degen_expand`, the search aborts rather than mixing them.

Two asymmetries are worth remembering. First, `-t` defaults to `0` (the contiguous index), which is a concrete value, not a wildcard: omitting `-t` selects the contiguous variant, it does not match spaced-seed variants. (`ikafssninfo` and `ikafssnserver`, by contrast, treat an unset `-t` as a wildcard.) Second, `-t 0` combined with `-template_type` is an error, because the contiguous index has no template type.

When `-accept_qdegen` is 0, queries containing IUPAC degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) are skipped with a `degen_rejected` skip-marker (TSV `*SKIPPED:degen_rejected`, JSON `"status": "skipped"`, SAM unmapped record with `XR:Z:degen_rejected`) and a stderr warning, and the exit code is 2. See "Skip reasons" below for the complete reason list and per-format representation. Set `-accept_qdegen 1` to allow such queries. K-mers containing exactly one degenerate base are expanded to all possible variants (e.g., R→A,G produces 2 k-mers; N→A,C,G,T produces 4) and used for search. K-mers whose per-position expansion product exceeds `-max_degen_expand` are skipped; when this occurs, a warning is emitted to stderr once per query indicating the query name and that such k-mers are ignored. (Note: this is a per-window unemit, not a whole-query skip — the rest of the query is still searched, and the unemit position is reflected in `Nhighfreq` for fractional thresholds.) In server mode (`ikafssnserver`), this warning is propagated through the protocol to `ikafssnclient`, which displays the same message. This matches the indexer's handling of subject-side degenerate bases.

`-seqidlist` and `-negative_seqidlist` are mutually exclusive. Both text (one accession per line) and binary (generated by `blastdb_aliastool -seqid_file_in`) formats are accepted, auto-detected by magic bytes.

Options that only take effect in a later stage are rejected with an error (non-zero exit) when an earlier `-mode` is selected, rather than being silently ignored. Stage 2 options (`-stage2_min_score`, `-stage2_max_gap`, `-stage2_max_lookback`, `-stage2_min_nhit_diag`, `-stage2_max_nhit_per_subject`, `-stage2_max_nhit_per_subject_mode`, `-stage2_max_nhit_in_total`) require `-mode 2` or higher; Stage 3 options (`-stage3_traceback`, `-stage3_gapopen`, `-stage3_gapext`, `-stage3_min_ppositive`, `-stage3_min_npositive`, `-stage3_score_matrix`, `-stage3_max_nhit_per_subject`, `-stage3_max_nhit_per_subject_mode`, `-stage3_max_nhit_in_total`, `-context_extend`, `-db`) require `-mode 3`. The same rule is enforced by `ikafssnclient` (CLI) and `ikafssnhttpd` (JSON request body).

**File descriptors in mode 3.** Stage 3 opens every volume of the BLAST DB at once so hits can be fetched in parallel, and each open volume costs two file descriptors that the NCBI toolkit holds for the lifetime of the mapping. `ikafssnsearch` raises its own `RLIMIT_NOFILE` soft limit to `2 × volumes + 64` before Stage 3 runs, and exits with instructions when the hard limit is too low. Modes 1 and 2 do not open the BLAST DB at all. See [OS settings for large databases](#os-settings-for-large-databases).

`ikafssnindex` likewise rejects contradictory option combinations: `-template_type` requires `-t > 0`, `-freq_threshold_part` requires `-mode 2` or `3`, and `-nthread_highfreq_filter` requires an active `-max_freq_build`. `ikafssnretrieve` rejects the efetch-only options (`-api_key`, `-batch_size`, `-max_nretry`, `-timeout`, `-range_threshold`) when run with local `-db`. For `ikafssnclient` and `ikafssninfo`, the HTTP authentication options (`-user`, `-http_user`, `-http_password`, `-netrc_file`) require `-http`, are mutually exclusive across the three methods, and `-http_password` requires `-http_user`.

**Timing output.** Alongside the per-stage counts, every run prints two timing lines to stderr. Both are at `info` level, so they appear without `-v`:

```
[INFO] Timing run_search (s): s1_open=0.512 s1_compute=12.345 s1_fold=0.031 s1_intotal=0.002 s2_open=0.220 s2a=45.678 s2b=3.210 s2_free=0.015 dedup=0.012 parent_topn=0.034 s2_intotal=0.005 total=62.064
[INFO] Timing overall (s): open_index=1.234 preprocess=0.456 run_search=62.064 convert=2.345 stage3=0.000 stage3_select=0.001 write=3.456 total=69.556
```

`Timing run_search` breaks the search orchestrator down into index volume open / close (`s1_open`, `s2_open`), Stage 1 accumulation (`s1_compute`), the mode 1 result fold (`s1_fold`), the Stage 1 and Stage 2 in-total (L) caps (`s1_intotal`, `s2_intotal`), Stage 2A and Stage 2B (`s2a`, `s2b`), release of the per-job Stage 2 state (`s2_free`), fragment-overlap dedup (`dedup`), and parent top-N selection (`parent_topn`). With `-mode 1` only the Stage 1 keys are printed. `Timing overall` covers the stages around it: index open, query preprocessing, `run_search` itself, conversion to output records, Stage 3 alignment, the Stage 3 caps, and output writing. `total` is measured rather than summed, so the gap between it and the sum of the other keys is time not attributed to any listed stage. Since the orchestrator is shared, `ikafssnserver` logs the `Timing run_search` line once per search request too.

**Examples:**

```bash
# Basic search
ikafssnsearch -ix ./index/mydb -query query.fasta -nthread 8

# Specify k-mer size (required if index contains multiple k values)
ikafssnsearch -ix ./index/mydb -k 11 -query query.fasta

# Increase sensitivity (keep more candidates per query)
ikafssnsearch -ix ./index/mydb -query query.fasta \
    -stage2_min_score 2 -stage1_max_nhit_per_subject 0 -stage1_max_nhit_in_total 2000

# Restrict to specific accessions
ikafssnsearch -ix ./index/mydb -query query.fasta -seqidlist targets.txt

# Exclude specific accessions
ikafssnsearch -ix ./index/mydb -query query.fasta -negative_seqidlist exclude.txt

# Fractional Stage 1 threshold (50% of query k-mers)
ikafssnsearch -ix ./index/mydb -query query.fasta -stage1_min_score 0.5

# Mode 3: pairwise alignment with traceback
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -stage3_max_nhit_in_total 5

# Mode 3: SAM output
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -output_format sam -o result.sam

# Mode 3: BAM output
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -output_format bam -o result.bam

# Mode 3: filter by percent positive
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -stage3_min_ppositive 90

# Mode 3: with context extension (50 bases each side)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context_extend 50 -stage3_max_nhit_in_total 5

# Mode 3: with context extension (0.1x query length each side)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context_extend 0.1 -stage3_max_nhit_in_total 5

# Pipe to ikafssnretrieve
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta
```

#### In-Silico PCR (Primer Mode)

`ikafssnsearch` and `ikafssnclient` support an in-silico PCR mode via the `-primer` option. In this mode, primer pairs are used to construct virtual PCR amplicon queries that are searched against the index. In `ikafssnclient`, the primer-to-query conversion is performed client-side, so the server requires no changes.

**Options:**

```
  -primer <path>              Primer FASTA file (mutually exclusive with -query)
  -insert_length <int>        Expected insert length (required with -primer)
  -stage1_primer_score <num>  Stage 1 threshold for primer (default: 0.5)
                              0 < v <= 1: fraction of primer-derived k-mers
                              v >= 2: absolute threshold
  -stage2_primer_score_add <int>  Stage 2 threshold addon (default: 1)
```

`-primer` and `-query` are mutually exclusive and cannot be specified together. When `-primer` is used, `-insert_length` is required. The options `-stage1_min_score` and `-stage2_min_score` cannot be used with `-primer`; use `-stage1_primer_score` and `-stage2_primer_score_add` instead.

**Primer pair processing:**

The primer FASTA file must contain an even number of sequences. Consecutive pairs of sequences form forward/reverse primer pairs (sequences 1+2 are the first pair, 3+4 are the second pair, and so on).

For each primer pair, a query sequence is constructed as:

```
fwd + N*insert_length + RC(rev)
```

Where `fwd` is the forward primer sequence, `N*insert_length` is a string of N characters of the specified length (representing the unknown insert region), and `RC(rev)` is the reverse complement of the reverse primer sequence.

**Threshold resolution logic:**

In primer mode, only the primer-derived k-mers in the constructed query contribute to matching (k-mers in the N region are skipped because they contain N). The `-stage1_primer_score` threshold is resolved as follows:

- `0 < v <= 1` (fraction): Stage 1 threshold is resolved as `ceil(Nprimer_kmer * v)` where Nprimer_kmer is the number of primer-derived k-mer positions
- `v >= 2` (absolute): used directly as an absolute threshold

For Stage 2, the minimum chain score is set to `max(Lf, Lr) + stage2_primer_score_add`, where Lf and Lr are the k-mer position counts of the forward and reverse primers respectively. Additionally, `-stage2_max_gap` is automatically set to the insert length unless explicitly overridden on the command line.

**Examples:**

```bash
# Basic In-Silico PCR
ikafssnsearch -ix ./index/mydb -primer primers.fasta -insert_length 500

# Set Stage 1 threshold to 80% of primer k-mers
ikafssnsearch -ix ./index/mydb -primer primers.fasta -insert_length 300 \
    -stage1_primer_score 0.8

# Mode 3: alignment with traceback
ikafssnsearch -ix ./index/mydb -primer primers.fasta -insert_length 500 \
    -mode 3 -stage3_traceback 1 -stage3_max_nhit_in_total 10
```

### ikafssnretrieve

Extract matched subsequences based on search results. Supports local BLAST DB extraction and remote retrieval via NCBI E-utilities (efetch).

**FASTA defline:** Each retrieved record is emitted as
`>parent_accession:start-end sseqid=<sseqid> query=<qseqid> strand=<+|-> score=<chainscore|alnscore>`,
with `start` / `end` 1-based and inclusive in **parent-relative coordinates** (the parent OID accession on the left, never a fragment-derived synthetic name). The record ID carries only the first accession so it stays a single well-formed token; the `sseqid=` field carries the value verbatim, including the `\x01`-joined multi-defline form, and is always present. Unlike the search output's `sstart` / `send`, this range always ascends: on a `-` strand row the range still names the same span in forward numbering and the emitted sequence is its reverse complement. On the local path the sequence fetch is routed through `BlastDbReader::get_subsequence(parent_oid, start, end)` so chromosome-scale parents only decode the requested window instead of the whole OID.

**Requirements of the local (`-db`) path:**

- The results file must carry a `volume` column: it tells `ikafssnretrieve` which BLAST DB volume holds each hit, so only the volumes that carry a hit are opened, one at a time. Every TSV and JSON layout ikafssn writes includes it. A results file without the column is rejected.
- The BLAST DB must be searchable by accession — a BLAST v5 database (which carries an LMDB), or a v4 database built with `makeblastdb -parse_seqids`. It must also be the database the index was built from; a different one will not resolve the accessions.
- Records are staged in a scratch file so the output keeps the input hit order regardless of the order the volumes are walked in. It is created in the directory of `-o`, or in the current directory when writing to standard output, and needs room for the uncompressed output. `TMPDIR` is deliberately not used, since it is a memory-backed tmpfs on many systems. The file is unlinked as soon as it is created, so nothing is left behind even if the process is killed.
- File descriptor use is a small constant — around five — no matter how many volumes the database has.

```
ikafssnretrieve [options]

Sequence source (one required):
  -db <path>              Local BLAST DB prefix
  -remote                 Retrieve from NCBI efetch

Input:
  -tsv <path>             Search results file (TSV format)
  (none)                  Read from stdin

Common options:
  -o <path>               Output FASTA file (default: stdout)
                          Suffix (.gz/.bz2/.xz/.zst) selects an output codec
  -context_extend <value>        Context extension (default: 2.0)
                          Integer: bases to add before/after match region
                          Decimal: multiplier of query length (qlen)
  -compression_level <int> Output compression level (defaults: gzip=6, bzip2=9, xz=6, zstd=3)
  -v, --verbose           Verbose logging

(Note: -tsv auto-detects gzip/bzip2/xz/zstd-compressed result files from
 their leading magic bytes; no flag is required.)

Remote options (-remote):
  -api_key <key>          NCBI API key (or NCBI_API_KEY env var)
  -batch_size <int>       Accessions per batch (default: 100)
  -max_nretry <int>          Max retries (default: 3)
  -timeout <int>          Request timeout in seconds (default: 30)
  -range_threshold <int>  Seq length threshold for individual fetch (default: 100000)
```

**Examples:**

```bash
# Local BLAST DB extraction (file input)
ikafssnsearch -ix ./index/mydb -query query.fasta -o results.tsv
ikafssnretrieve -db nt -tsv results.tsv -o matches.fasta

# Local BLAST DB extraction (pipe)
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# Server search results also work
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# Remote retrieval via NCBI efetch
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -remote > matches.fasta

# Remote retrieval with API key (higher throughput)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta \
    | ikafssnretrieve -remote -api_key XXXXXXXX > matches.fasta

# Include 50bp of context around match region
ikafssnretrieve -db nt -tsv results.tsv -context_extend 50

# Context as fraction of query length (0.1x each side)
ikafssnretrieve -db nt -tsv results.tsv -context_extend 0.1
```

### ikafssnserver

Search daemon. Keeps index files memory-mapped and accepts search requests via UNIX domain socket or TCP socket.

When a `-db` is given, every volume of that BLAST DB is opened at startup and stays open for the server's lifetime: mode 3 requests touch several volumes each and run concurrently, so re-opening them per request would pay the setup cost every time without ever releasing a descriptor. Each open volume costs two file descriptors, so a server needs roughly `2 × (total BLAST DB volumes) + 256`. The server raises its own `RLIMIT_NOFILE` soft limit to cover that before it loads anything; if the hard limit is too low it reports what to change and exits at startup rather than failing mid-request. See [OS settings for large databases](#os-settings-for-large-databases).

```
ikafssnserver [options]

Required:
  -ix <prefix>            Index prefix (repeatable for multi-DB)

Listener (at least one required):
  -socket <path>          UNIX domain socket path
  -tcp <host>:<port>      TCP listen address

Index-variant load filter (each unset field is a wildcard):
  -k <int>                Load only this k-mer length (default: any)
  -t <int>                Load only this template length: 0=contiguous, 13/15/16/18/21
                          (default: any; -t 0 loads the contiguous index only)
  -template_type <str>    coding, optimal, contiguous, both (default: any; invalid with -t 0)
  -min_seq_length <int>   Load only variants with this min_seq_length (default: any)
  -min_length_split <int> Load only variants with this min_length_split (default: any)
  -overlap_length <int>   Load only variants with this overlap_length (default: any)
  -max_freq_build <int>   Load only variants with this max_freq_build (default: any)
  -max_degen_expand <int> Load only variants with this max_degen_expand (default: any).
                          Load filter only; the query-side max_degen_expand is supplied by
                          the client (server fallback when a request omits it: 16).

Options:
  -nthread <int>          Worker threads (default: all cores)
  -max_queue_size <int>   Max concurrent query sequences globally (default: 1024)
  -max_nseq_per_req <int> Max sequences accepted per request (default: thread count)
  -max_concurrent_search <int>  Cap concurrent searches sharing -memory_limit's
                          residual posting_budget, bounding in-flight posting
                          heap (Stage 3) (default: 0 = unlimited)
  -pid <path>             PID file path
  -db <path>              BLAST DB path for mode 3 (repeatable, paired with -ix;
                          default: same as corresponding -ix prefix)
  -stage1_max_nhit_per_subject <int>  Default max Stage 1 candidates per parent (default: 1, 0=unlimited)
  -stage1_max_nhit_per_subject_mode <1|2|3|4>  Default per-subject selection mode (default: 3)
  -stage1_max_nhit_per_volume <int>  Default max Stage 1 candidates per (query,volume,strand) (default: 0)
  -stage1_max_nhit_in_total <int>  Default max Stage 1 candidates per query across volumes (default: 0)
  -stage1_min_score <num> Default Stage 1 minimum score (default: 0.5)
                          Integer (>= 1) or fraction (0 < P < 1)
  -stage2_min_score <int> Default minimum chain score (default: 0 = adaptive)
  -stage2_max_gap <int>   Default chaining gap tolerance (default: 100)
  -stage2_max_lookback <int>  Default chaining DP lookback window (default: 64, 0=unlimited)
  -stage2_max_nhit_per_subject <int>  Default max chains per subject (default: 1, 0=unlimited)
  -stage2_max_nhit_per_subject_mode <1|2|3|4>  Default per-subject selection mode (default: 3)
                           1/2=top-N (no ties), 3/4=top-N + score ties;
                           1/3=strands merged per parent, 2/4=strands separate
  -stage2_max_nhit_in_total <int>  Default max Stage 2 chains per query across volumes (default: 0)
  -stage2_min_nhit_diag <int> Default diagonal filter min hits (default: 1)
  -context_extend <value>        Default context extension (default: 2.0)
                          Integer: bases to extend; Decimal: query length multiplier
  -stage3_traceback <0|1> Default traceback mode (default: 0)
  -stage3_gapopen <int>   Default gap open penalty (default: 10)
  -stage3_gapext <int>    Default gap extension penalty (default: 1)
  -stage3_min_ppositive <num>  Default min percent positive (default: 0)
  -stage3_min_npositive <int>  Default min positive-scoring positions (default: 0)
  -stage3_max_nhit_per_subject <int>  Default max Stage 3 hits per subject (default: 1, 0=unlimited)
  -stage3_max_nhit_per_subject_mode <1|2|3|4>  Default per-subject selection mode (default: 3)
  -stage3_max_nhit_in_total <int>  Default max Stage 3 hits per query across volumes (default: 0)
  -stage3_score_matrix <name>  Default score matrix: degmatch, dnafull, nuc44 (default: degmatch)
  -accept_qdegen <0|1>    Default accept queries with degenerate bases (default: 1)
  -memory_limit <size>    Search memory budget (default: half of RAM)
                          Pins .khx/.ksx metadata; residual caps the
                          Stage 3 posting heap and concurrent-search pool.
                          Accepts K, M, G suffixes
  -shutdown_timeout <int> Graceful shutdown timeout in seconds (default: 30)
  -v, --verbose           Verbose logging
```

The load filter restricts which index variants the server loads, using the same wildcard semantics as `ikafssninfo` local mode (every unset field is a wildcard; an unset `-t` is a full wildcard; `-t 0` combined with `-template_type` is an error). All eight identifying parameters — **including `-max_degen_expand`** — are load filters here. `-max_degen_expand` is *only* a load filter on the server: the query-side `max_degen_expand` (k-mer expansion + tie-break) is supplied by the client, and the server applies an internal fallback of 16 only when a request omits it. Leaving `-max_degen_expand` unset loads every build-time value (so the per-request tie-break can choose); setting it restricts the loadable set, and clients can then only select among the loaded variants.

The server exposes each loaded variant as a distinct group over the wire (the info response carries all eight identifying parameters per group). A search request carries the client-resolved variant identity; the server matches `t` and `template_type` literally (0 = contiguous), the four indexing fields exactly, and — when several `max_degen_expand` variants remain — picks the one whose build value equals the request's `max_degen_expand` (else the largest). For a `template_type=both` request the server selects a coding+optimal pair that shares one `max_degen_expand` (no mixing). `-min_query_length` and the long-query cap are derived per request from the selected variant's `min_seq_length` and `overlap_length`.

**Examples:**

```bash
# Listen on UNIX socket
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -nthread 16

# Listen on TCP (remote access)
ikafssnserver -ix ./nt_index -tcp 0.0.0.0:9100 -nthread 32

# Both simultaneously
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -tcp 0.0.0.0:9100

# Daemon operation (started from systemd)
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -pid /var/run/ikafssn.pid

# Mode 3 support: specify BLAST DB and Stage 3 defaults
ikafssnserver -ix ./nt_index -db nt -socket /var/run/ikafssn.sock \
    -stage3_traceback 1 -context_extend 50

# Multi-DB: serve two databases in one process
ikafssnserver -ix ./nt_index -db nt -ix ./rs_index -db refseq_genomic \
    -socket /var/run/ikafssn.sock -nthread 32
```

**Operational characteristics:**

- One process can serve multiple BLAST DB indexes simultaneously. Specify `-ix` (and optionally `-db`) multiple times to load several databases. Each database is identified by its basename (the last path component of the `-ix` prefix) and clients must specify `-db <name>` when the server hosts more than one database.
- If `-db` is specified, the count must match the number of `-ix` flags (paired in order). Databases without a `-db` override default to the `-ix` prefix as the BLAST DB path. A database with no `-db` path supports modes 1-2 only (max_mode=2); providing `-db` enables mode 3 (max_mode=3).
- If the index prefix matches multiple index variants (differing in any of the eight identifying parameters — k, t, template_type, min_seq_length, min_length_split, overlap_length, max_freq_build, max_degen_expand), all are loaded (subject to the load filter) and each client request selects one by carrying its resolved variant identity. The load filter can restrict which variants are loaded in the first place.
- On SIGTERM/SIGINT, performs graceful shutdown: stops accepting new connections, waits for in-flight requests to complete (up to `-shutdown_timeout` seconds), then exits.
- **Per-sequence concurrency control:** The server limits concurrency at the per-sequence level, not per-connection. When a request arrives, the server attempts to acquire permits for each valid query sequence. If the global limit (`-max_queue_size`) is reached, excess sequences are returned to the client as "rejected" for retry. The `-max_nseq_per_req` option caps how many permits a single request can acquire, preventing one large request from monopolizing all slots.
- **Inter-request posting budget pool (`-max_concurrent_search`):** The default value `0` lets every concurrent request independently use the full residual `posting_budget` (the leftover of `-memory_limit` after `.khx` / `.ksx` `WILLNEED`). When set to `N >= 1`, requests share that budget through a lease/release pool: each request holds a single lease for its whole search, and at most `N` searches run at once. This bounds in-flight posting heap — consumed by the Stage 3 alignment batcher — to `posting_budget` regardless of how many connections are open, at the cost of serialising additional searches until a lease is released. Use this when sizing `-memory_limit` for a single request would be undesirable but operators still want a hard cap on per-process posting RAM.

### ikafssnhttpd

HTTP REST API daemon. Connects to one or more `ikafssnserver` instances and exposes search as an asynchronous, polling-based REST API backed by a SQLite job store and a worker pool. Uses the Drogon framework. Multiple backends can be specified for multi-database support or load balancing of same-database replicas.

On startup, it connects to all configured backends to cache their capabilities (retrying with exponential backoff for up to 30 seconds). If the same database name appears on multiple backends, cross-server validation ensures that for each shared (db, k) pair, total sequence counts and total bases are identical; mismatches cause a startup error. Backends are allowed to have different k-value sets for the same database (e.g., server A has k=10, server B has k=10 and k=11); the merged capabilities expose the union of all k-value groups. Search requests are validated against the merged capabilities (synchronous, no backend round-trip) to reject obviously invalid requests immediately. Validated requests are persisted in the job store and then dispatched by the worker pool to the best available backend based on priority and slot availability.

**Routing and health:**

- **Priority**: Backends are prioritized by CLI argument order (first = highest priority).
- **Selection**: For each search request, the backend with the highest priority and available effective capacity is selected. Effective capacity considers both slot availability (`max_queue_size - queue_depth`) and per-request cap (`max_nseq_per_req`), taking the minimum of the two. If all backends are full, the highest-priority one is used.
- **Pre-check**: Before each search, a fresh info request is sent to the selected backend to verify connectivity.
- **Exclusion**: If a backend fails to respond (connection error on info or search), it is excluded for `-exclusion_time` seconds. Excluded backends are automatically re-checked during heartbeat and re-enabled once reachable.
- **Heartbeat**: A background thread refreshes all backends' info every `-heartbeat_interval` seconds.
- **Retry in httpd**: A failed dispatch is retried up to `-max_nretry` times by the worker pool before the job is marked `failed`. Per-query rejections returned by a backend are also retried by the worker.

```
ikafssnhttpd [options]

Backend connection (at least one required; order = priority):
  -server_socket <path>      UNIX socket path to ikafssnserver
  -server_tcp <host>:<port>  TCP address of ikafssnserver

Job store (async REST):
  -sqlite_db <path>           SQLite job store path
                              (default: /var/lib/ikafssnhttpd/jobs.db)
  -query_dir <path>           Per-job query file directory.  Each queued or
                              running job's request body is persisted here as
                              one zstd file; the file is removed when the job
                              reaches a terminal status.
                              (default: /var/lib/ikafssnhttpd/queries)
  -query_compression_level <int>  Zstd compression level applied to plaintext
                              JSON uploads (Content-Type: application/json).
                              Bodies uploaded as application/zstd are stored
                              verbatim (their client-chosen level is
                              preserved) and ignore this option.
                              (default: 3, range: 1-22)
  -result_dir <path>          Per-job result file directory (one zstd file per
                              completed job, served via sendfile(2) to clients)
                              (default: /var/lib/ikafssnhttpd/results)
  -result_compression_level <int>  Zstd compression level for result files
                              (default: 3, range: 1-22)
  -retention_time <int>       Done/failed job retention in seconds (default: 86400)
                              After this duration the housekeeper purges the row
                              and removes the matching result file.
  -max_nretry <int>           Per-job dispatch retry cap (default: 3)
  -nthread_worker <int>       Backend search worker pool size (default: 4)

Options:
  -listen <host>:<port>       HTTP listen address (default: 0.0.0.0:8080)
  -path_prefix <prefix>       API path prefix (e.g., /nt)
  -nthread <int>              Drogon I/O threads (default: all cores)
  -heartbeat_interval <int>   Heartbeat interval in seconds (default: 3600)
  -exclusion_time <int>       Backend exclusion time in seconds (default: 3600)
  -pid <path>                 PID file path
  -v, --verbose               Verbose logging
```

**REST API endpoints (async, polling-based):**

| Method | Path | Description |
|---|---|---|
| POST | `/api/v1/jobs` | Submit a search job. Body is a JSON document. `ikafssnclient` always sends it as a single zstd frame with `Content-Type: application/zstd`; ad-hoc HTTP clients (curl, scripts, etc.) may send it uncompressed with `Content-Type: application/json` to skip the pre-compression step. Any other Content-Type is rejected with HTTP 415. Returns `{job_id, status}`. |
| GET  | `/api/v1/jobs/{job_id}` | Job status (`queued` / `running` / `done` / `failed` / `timeout`) |
| GET  | `/api/v1/jobs/{job_id}/result` | Fetch the serialized SearchResponse blob (only valid once status = `done`). Response body is a single zstd frame; `Content-Type: application/zstd`. Decode with `zstd -d` or `ZSTD_decompress` before deserialising. |
| GET  | `/api/v1/info` | Aggregated index information from all backends |
| GET  | `/api/v1/health` | Health check (OK if any backend is reachable) |

The `/api/v1/info` response aggregates databases from all healthy backends. For databases served by multiple backends, capacity is reported per mode in a `modes` array within each kmer_group, showing the sum of `max_queue_size`, `queue_depth`, and `max_nseq_per_req` (computed as `sum(min(available_i, per_req_i))` across backends) across all serving backends. A top-level `max_nseq_per_req` field reports the minimum across all modes.

`POST /api/v1/jobs` validates the body and persists it to the per-job query store at `<query_dir>/<ab>/<job_id>.json.zst` (sharded by the first two characters of `job_id`) before inserting the SQLite row. Bodies uploaded as `application/zstd` are stored verbatim — the client's compression level is preserved and the daemon never re-compresses; plaintext JSON uploads (`application/json`) are compressed at `-query_compression_level` before they hit disk. The worker pool then dispatches the request asynchronously to a backend. When a job finishes successfully the worker writes its serialized `SearchResponse` to a per-job zstd file at `<result_dir>/<ab>/<job_id>.bin.zst` and stamps the SQLite row to `done`; `GET /api/v1/jobs/{job_id}/result` then streams that file directly via `sendfile(2)` with `Content-Type: application/zstd`, so no in-memory copy of the result is held by `ikafssnhttpd`. The matching query file is unlinked at the same point — only queued / running jobs occupy `-query_dir`, so a backlog never inflates `jobs.db` (which now carries metadata only). Clients poll `GET /api/v1/jobs/{job_id}` until status reaches `done` and download the result. Jobs are kept for `-retention_time` seconds; the housekeeper purges expired rows from SQLite and then removes the matching result file (and any leftover query file as insurance). On startup, leftover `*.bin.zst.tmp` and `*.json.zst.tmp` files (from a hard kill mid-write) are swept once. A periodic orphan sweep also removes result and query files whose SQLite row has already disappeared (rare: only after a failed `mark_done` or a worker crash between `mark_done` and the query unlink).

**Examples:**

```bash
# Single backend via UNIX socket
ikafssnhttpd -server_socket /var/run/ikafssn.sock -listen 0.0.0.0:8080

# Single backend via TCP
ikafssnhttpd -server_tcp 10.0.1.5:9100 -listen 0.0.0.0:8080

# Multiple backends for load balancing (same DB on two servers)
ikafssnhttpd -server_tcp server1:9100 -server_tcp server2:9100 -listen :8080

# Multiple backends with different DBs
ikafssnhttpd -server_socket /var/run/nt.sock -server_socket /var/run/rs.sock -listen :8080

# Mixed: primary + failover
ikafssnhttpd -server_socket /var/run/primary.sock -server_tcp backup:9100 -listen :8080
```

### ikafssnclient

Client command. Connects to `ikafssnserver` via socket or `ikafssnhttpd` via HTTP. Output format is identical to `ikafssnsearch`. Before sending any queries, the client performs pre-flight validation by fetching server capabilities and checking that the requested database name, k-mer size, and mode are valid. It then resolves the index variant exactly the way `ikafssnsearch` does — applying the same `-k` / `-t` / `-template_type` / `-min_seq_length` / `-min_length_split` / `-overlap_length` / `-max_freq_build` filter and the same template_type resolution (single, both-required, or wildcard-with-coding/optimal-fallback; both-pairs share one `max_degen_expand`, no mixing) to the server's per-variant group list, then sending the resolved identity in the request. As with `ikafssnsearch`, the result must collapse to one variant (or one both-pair) using the seven parameters other than `max_degen_expand`; `max_degen_expand` itself is left to the server's tie-break. An ambiguous or empty selection is reported locally before any query data is transmitted. Invalid parameters produce an error with available database listings. The client uses the server's `max_nseq_per_req` and available slot count to automatically split queries into appropriately-sized batches, avoiding oversized requests that would be partially rejected. Within each batch, if the server still rejects some query sequences due to concurrency limits, the client automatically retries the rejected queries with exponential backoff (30s, 60s, 120s, 120s, ...) until all queries are processed.

**Socket/TCP mode (synchronous) — checkpointing:** When connected via `-socket` or `-tcp`, the client saves intermediate results to a temporary directory during batch processing. If the process is interrupted (e.g., Ctrl+C, network failure), re-running the same command resumes from where it left off, skipping already-completed queries. The temporary directory is named `{output}.{input}.{ix_name}.{kk}.ikafssn.tmp/` and is automatically cleaned up after successful completion. A directory-based lock prevents concurrent runs with the same parameters. Resume validation checks the search parameters, input file SHA256, and integrity of each batch file.

**HTTP mode (asynchronous, polling-based):** When connected via `-http`, the client splits the input into one job per batch and submits each via `POST /api/v1/jobs` to `ikafssnhttpd`. The request body is zstd-compressed before being sent (`Content-Type: application/zstd`) so wire bandwidth stays small even for million-base query batches; the server transparently decompresses on receipt. Each job's metadata is persisted under `~/.ikafssnclient/<group_id>/<job_id>.json` along with a zstd-compressed defline cache (`<job_id>.deflines.zst`); the client then polls `GET /api/v1/jobs/{job_id}` until the job reaches a terminal state and downloads the result. The downloaded body is already a single zstd frame so it is written verbatim to `<job_id>.result.bin.zst` (no double-compression), then decoded and merged into the user's output file. Because all state is on disk, the client is resumable across reboots — `-resume <id>` re-enters the polling loop for an existing group or job, and `-submit_only` returns the `group_id` immediately after submission so a separate process (or a later invocation of `-resume`) can collect the result. Note: `*.result.bin.zst` files written by older clients used a doubled-zstd format; they are not compatible with the new code and any partially completed groups should be discarded after upgrade.

```
ikafssnclient [options]

Connection (one required):
  -socket <path>           UNIX domain socket path
  -tcp <host>:<port>       TCP server address
  -http <url>              ikafssnhttpd URL (e.g., http://example.com:8080)
                           HTTP mode uses async REST polling.

Async REST job management (HTTP mode only; mutually exclusive with a search):
  -submit_only             Submit each batch, print the group_id, then exit
                           without polling.  Use with -resume later to collect.
  -jobs                    List all locally-tracked job groups under
                           ~/.ikafssnclient/ (no server contact)
  -job_detail <id>          Show jobs in a group, or detail of a single job
  -resume <id>             Resume polling for an existing group or single job

Required (for a fresh search):
  -query <path>            Query FASTA file (- for stdin)
  -ix <name>               Target database name on server

Primer mode (alternative to -query):
  -primer <path>           Primer pair FASTA (mutually exclusive with -query)
  -insert_length <int>     Expected insert length (required with -primer)
  -stage1_primer_score <num>  Stage 1 threshold (0<v<=1: fraction, v>=2: absolute; default: 0.5)
  -stage2_primer_score_add <int>  Stage 2 threshold addon: max(Lf,Lr) + N (default: 1)

Options:
  -o <path>                Output file (default: stdout)
  -k <int>                 K-mer size (default: server default)
  -mode <1|2|3>            Search mode (default: server default)
  -stage1_max_nhit_per_subject <int>  Max Stage 1 candidates per parent (default: server default)
  -stage1_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: server default)
  -stage1_max_nhit_per_volume <int>  Max Stage 1 candidates per (query,volume,strand) (default: server default)
  -stage1_max_nhit_in_total <int>  Max Stage 1 candidates per query across volumes (default: server default)
  -stage1_min_score <num>  Stage 1 minimum score (default: server default)
                           Integer (>= 1) or fraction (0 < P < 1)
  -stage2_min_score <int>  Minimum chain score (default: server default)
                           0 = explicitly request adaptive mode
  -stage2_max_gap <int>    Chaining gap tolerance (default: server default)
  -stage2_max_lookback <int>  Chaining DP lookback window (default: server default)
  -stage2_max_nhit_per_subject <int>  Max chains per subject (default: server default)
  -stage2_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3)
                           1/2=top-N (no ties), 3/4=top-N + score ties;
                           1/3=strands merged per parent, 2/4=strands separate
  -stage2_max_nhit_in_total <int>  Max chains per query across all volumes (requires -mode 2 or higher)
  -stage2_min_nhit_diag <int> Diagonal filter min hits (default: server default)
  -context_extend <value>         Context extension (default: 2.0)
                           Integer: bases to extend; Decimal: query length multiplier
  -stage3_traceback <0|1>  Enable traceback (default: server default)
  -stage3_gapopen <int>    Gap open penalty (default: server default)
  -stage3_gapext <int>     Gap extension penalty (default: server default)
  -stage3_min_ppositive <num> Min percent positive filter (default: server default)
  -stage3_min_npositive <int> Min positive-scoring positions filter (default: server default)
  -stage3_max_nhit_per_subject <int>  Max hits per subject (requires -mode 3)
  -stage3_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3; requires -mode 3)
  -stage3_max_nhit_in_total <int>  Max hits per query across all volumes (requires -mode 3)
  -stage3_score_matrix <name> Score matrix: degmatch, dnafull, nuc44 (default: server default)
  -seqidlist <path>        Include only listed accessions
  -negative_seqidlist <path>  Exclude listed accessions
  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: server default)
  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 1)
  -max_degen_expand <int>  Max degenerate expansion (default: 16, max: 256, 0/1: disable).
                           A query parameter identical to ikafssnsearch (NOT delegated to the
                           server); it also breaks ties among max_degen_expand variants.
  -min_query_length <int>  Minimum query length, in bases (default: 64).
                           Pre-filters queries client-side so that the
                           server's identical check never has to.  Must
                           be >= the server-reported min_seq_length
                           (surfaced via the InfoResponse); a smaller
                           value is rejected at the pre-flight stage.
  -t <int>                 Template length for spaced seeds (default: 0 = contiguous, NOT a wildcard)
                           0: contiguous k-mers; 13, 15, 18 (k=8-9); 16, 18, 21 (k=11-12)
  -template_type <str>     Template type for spaced seeds (default: both for -t>0; invalid with -t 0)
                           coding: use coding index only
                           optimal: use optimal index only
                           both: merge coding and optimal indexes at search time
  -min_seq_length <int>    Select the index variant with this min_seq_length (default: any)
  -min_length_split <int>  Select the index variant with this min_length_split (default: any)
  -overlap_length <int>    Select the index variant with this overlap_length (default: any)
  -max_freq_build <int>    Select the index variant with this max_freq_build (default: any)
  -output_format <tsv|json|sam|bam>  Output format (default: tsv)
  -compression_level <int> Output compression level (defaults: gzip=6, bzip2=9, xz=6, zstd=3)
                           Codec is selected by -o suffix (.gz/.bz2/.xz/.zst); SAM/BAM reject all four
  -v, --verbose            Verbose logging

(Note: -query and -primer auto-detect gzip/bzip2/xz/zstd-compressed FASTA
 inputs from their leading magic bytes; no flag is required.)

HTTP Authentication (HTTP mode only):
  -user <user:password>   Credentials (curl-style)
  -http_user <USER>       Username (wget-style)
  -http_password <PASS>   Password (used with -http_user)
  -netrc_file <path>      .netrc file for credentials
```

**Examples:**

```bash
# UNIX socket (local, single-DB server)
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta

# TCP (local or remote)
ikafssnclient -tcp 10.0.1.5:9100 -ix nt -query query.fasta

# HTTP
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta

# Restrict search scope
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -seqidlist targets.txt

# Pipe to ikafssnretrieve
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# Specify k-mer size
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -k 9

# HTTP with Basic Auth (curl-style)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta -user admin:secret

# HTTP with Basic Auth (wget-style)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta -http_user=admin -http_password=secret

# HTTP with .netrc credentials
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta -netrc_file ~/.netrc

# Mode 3: pairwise alignment with traceback
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -mode 3 -stage3_traceback 1

# Mode 3: SAM output
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -mode 3 -stage3_traceback 1 -output_format sam -o result.sam

# In-Silico PCR (primer mode)
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -primer primers.fasta -insert_length 500

# In-Silico PCR with custom thresholds
ikafssnclient -tcp 10.0.1.5:9100 -ix nt -primer primers.fasta -insert_length 300 \
    -stage1_primer_score 0.8

# In-Silico PCR with mode 3 alignment
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -primer primers.fasta -insert_length 500 \
    -mode 3 -stage3_traceback 1 -stage3_max_nhit_in_total 10
```

### ikafssninfo

Display index/database information. Supports two modes: **local mode** (reads index files directly) and **remote mode** (queries a running `ikafssnserver` or `ikafssnhttpd`).

```
ikafssninfo [options]

Required (one of):
  -ix <prefix>             Index prefix [local mode]
  -socket <path>           UNIX socket to ikafssnserver [remote mode]
  -tcp <host>:<port>       TCP address of ikafssnserver [remote mode]
  -http <url>              ikafssnhttpd URL [remote mode]

Local mode options (index-variant filter; each unset field is a wildcard):
  -db <path>               BLAST DB prefix (default: auto-detect from -ix)
  -k <int>                 K-mer length (default: any)
  -t <int>                 Template length: 0=contiguous, 13/15/16/18/21
                           (default: any; -t 0 selects the contiguous index only)
  -template_type <str>     coding, optimal, contiguous, both (default: any; invalid with -t 0)
  -min_seq_length <int>    Filter by min_seq_length (default: any)
  -min_length_split <int>  Filter by min_length_split (default: any)
  -overlap_length <int>    Filter by overlap_length (default: any)
  -max_freq_build <int>    Filter by max_freq_build (default: any)
  -max_degen_expand <int>  Filter by max_degen_expand (default: any)
  -stats <0|1>             Compute k-mer frequency distribution (default: 0; slow)

Remote HTTP authentication:
  -user <user:password>   Credentials (curl-style)
  -http_user <USER>       Username (wget-style)
  -http_password <PASS>   Password (used with -http_user)
  -netrc_file <path>      .netrc file for credentials

Options:
  -v, --verbose            Verbose output
```

`-ix` and remote options (`-socket`, `-tcp`, `-http`) are mutually exclusive. Only one remote option may be specified at a time.

**Local mode** reads index files directly and displays detailed statistics. When `-db` is not specified, `ikafssninfo` attempts to auto-detect the BLAST DB by checking whether the index prefix path corresponds to a valid BLAST DB. Unlike `ikafssnsearch` / `ikafssnclient`, local mode does not resolve down to a single variant: it reports **every** index variant that passes the filter, each under its own `=== Index variant: … ===` heading, and statistics are never combined across variants. The filter options use wildcard semantics — an unset `-t` matches any template length (whereas in `ikafssnsearch` an unset `-t` means contiguous), `-max_degen_expand` is an ordinary filter here (not a tie-break), and `-t 0` combined with `-template_type` is an error.

Local mode output includes:

- K-mer length (k) and integer type (uint16/uint32)
- Number of volumes
- Per-volume statistics: parent (BLAST OID) count, fragment (internal SeqId) count, fragment-length distribution (min / median / mean / max), total postings, file sizes, excluded k-mer count (if `.khx` present)
- Overall statistics: total parents, total fragments, aggregated fragment-length distribution, total postings, total index size, compression ratio
- The variant's eight identifying parameters (`k`, `t`, `template_type`, `min_seq_length`, `min_length_split`, `overlap_length`, `max_freq_build`, `max_degen_expand`) parsed from the index file name
- With `-stats 1`: k-mer frequency distribution (min, max, mean, percentiles)
- With `-db` (or auto-detected): BLAST DB title, sequence count, total bases, volume paths

The "Parents" / "Fragments" split makes it visible whether the index uses fragment splitting (`min_length_split > 0` and `Fragments > Parents`) or the degenerate one-fragment-per-parent layout (`Fragments == Parents`).

**Remote mode** queries a running server and displays its capabilities.

Remote mode output includes:

- Active/max sequence slots
- Per-database information: name, default k, max mode, and one group per index variant. Each group lists its eight identifying parameters (`k`, `t`, `template_type`, `min_seq_length`, `min_length_split`, `overlap_length`, `max_freq_build`, `max_degen_expand`) along with volume counts, sequence counts, total bases, and posting statistics
- With `-v`: per-volume details (sequence count, total bases, postings) within each variant group

**Examples:**

```bash
# Local: basic index info
ikafssninfo -ix ./index/mydb

# Local: include BLAST DB info
ikafssninfo -ix ./index/mydb -db mydb

# Local: detailed frequency distribution
ikafssninfo -ix ./index/mydb -stats 1

# Remote: query server via UNIX socket
ikafssninfo -socket /var/run/ikafssn.sock

# Remote: query server via TCP
ikafssninfo -tcp 10.0.1.5:9100

# Remote: query server via HTTP
ikafssninfo -http http://search.example.com:8080

# Remote: verbose (show per-volume details)
ikafssninfo -socket /var/run/ikafssn.sock -v

# Remote: HTTP with authentication
ikafssninfo -http http://search.example.com:8080 -user admin:secret
```

## Search Pipeline

ikafssn uses a three-stage search pipeline:

The default parameters prioritize throughput: `stage1_min_score=0.5` (fractional) filters candidates by requiring at least 50% of query k-mers to match. Each stage carries the same family of candidate-count limits — per-subject (N), per-volume (M, Stage 1 only), and in-total (L) — so the result count can be bounded at whichever stage is the terminal stage for the chosen `-mode` (mode 1 → Stage 1 L, mode 2 → Stage 2 L, mode 3 → Stage 3 L). Stage 1 keeps at most one candidate per parent by default (`stage1_max_nhit_per_subject=1`); its per-volume (`stage1_max_nhit_per_volume`) and in-total (`stage1_max_nhit_in_total`) limits are unlimited by default. Because the Stage 1 and Stage 2 limits prune candidates before the expensive Stage 3 alignment, they are the levers for bounding Stage 3 cost; the Stage 3 limits themselves are applied after alignment (by alnscore) and refine the output without reducing the alignment work.

1. **Stage 1 (Candidate Selection):** Scans ID postings for each query k-mer and accumulates a **coverscore** per sequence (the number of query k-mer positions matching the sequence). Sequences at or above `stage1_min_score` are selected as candidates. The candidate set is then capped, in the order N → M → L (all coverscore, tie-inclusive): `stage1_max_nhit_per_subject` (N) keeps the top-N per parent within each (query, volume, strand); `stage1_max_nhit_per_volume` (M) keeps the top-M per (query, volume, strand); `stage1_max_nhit_in_total` (L) keeps the top-L per query across all volumes and strands.

2. **Stage 2 (Collinear Chaining):** For each candidate, collects position-level hits from the `.kpx` file, applies a diagonal filter, and runs a chaining DP to find the best collinear chain. The chain length is reported as **chainscore**. Chains with `chainscore >= stage2_min_score` are reported. The DP inner loop is limited by `-stage2_max_lookback` (default: 64), restricting each hit to consider only the preceding B hits as potential chain predecessors. This reduces worst-case complexity from O(n²) to O(n×B) when a single query×subject pair has a very large number of hits. Set to 0 for an unlimited lookback (O(n²)). When `-stage2_max_nhit_per_subject` is greater than 1 (or 0 for unlimited), multiple non-overlapping chains are extracted per subject using greedy best-chain removal: the best chain is found and its hits are removed, then the DP is re-run on the remaining hits, repeating until the limit is reached or no chain meets `min_score`. `-stage2_max_nhit_per_subject` (N) is applied twice — once per (query, fragment, strand) during extraction, and once more per (query, parent[, strand]) after the overlap-region dedup — both scored by chainscore. `-stage2_max_nhit_per_subject_mode` (default: 3) controls both stages: modes 1/2 keep the top N chains with no tie handling, modes 3/4 keep the top N plus every chain tying the N-th chainscore; modes 1/3 merge the two strands into one per-parent group, modes 2/4 keep the strands separate. The sentinel value 0 resolves to mode 3; explicit values must be 1..4. This option requires `-mode 2` or higher. After the per-subject selection, `-stage2_max_nhit_in_total` (L) caps the chains kept per query across all volumes and strands by chainscore (tie-inclusive; 0 = unlimited), bounding the chains that enter Stage 3.

3. **Stage 3 (Pairwise Alignment):** For each Stage 2 hit, retrieves the subject subsequence from the BLAST DB (with optional context extension via `-context_extend`), and performs semi-global pairwise alignment using the Parasail library (using the score matrix specified by `-stage3_score_matrix`, default: DEGMATCH). The alignment score (**alnscore**) is computed for all hits. When `-stage3_traceback 1` is enabled, CIGAR strings, percent positive (ppositive), positive-scoring position count (npositive), negative-scoring count (nnegative), and aligned sequences (with gaps) are also computed. Hits can be filtered by `-stage3_min_ppositive` and `-stage3_min_npositive` (traceback mode only). After alignment and dedup, `-stage3_max_nhit_per_subject` (N, default 1) keeps the top-N hits per (qseqid, sseqid[, sstrand]) group by alnscore — with `-stage3_max_nhit_per_subject_mode` controlling tie / strand semantics as in Stage 1/2 — and `-stage3_max_nhit_in_total` (L) then caps the hits kept per query across all volumes by alnscore (tie-inclusive; 0 = unlimited). These limits are applied after alignment, so they refine the result set rather than reduce the alignment work. Subject subsequences are fetched hit-parallel using the full `-nthread` pool; within each TBB task, hits are visited in (volume, OID) order to preserve sequential mmap locality.

**Context extension (`-context_extend`):** an integer is a base count, a decimal is a multiplier of that hit's query length (`2` means 2 bases, `2.0` means twice the query length). The amount is applied to **each side** of the Stage 2 chain and clamped to the parent, so a 100 bp query with `-context_extend 2.0` widens a 100 bp chain to at most `200 + 100 + 200` = 500 bases. Clamping is per side and never compensated: a chain starting at the first base of the parent still extends only the configured amount downstream.

Because the alignment is semi-global, a free end gap can skip a prefix or suffix of only one of the two sequences, so whichever of the query and the subject window is shorter is consumed in full. With a multiplier the window is normally longer than the query, and the whole query is aligned (`qstart` = 1, `qend` = `qlen`) even where its ends match nothing. Use `-context_extend 0`, or a small base count, when the query may carry sequence absent from the subject (primers, adapters) and the reported boundaries should stop at the match.

**Adaptive `-stage2_min_score` (default):** When `-stage2_min_score 0` (the default), the minimum chain score is set adaptively per query to the resolved Stage 1 threshold. With fractional `-stage1_min_score` (e.g. `0.5`), this means each query gets a per-query adaptive threshold based on its k-mer composition. With absolute `-stage1_min_score`, the configured value is used. Set `-stage2_min_score` to a positive integer to override this behavior with a fixed threshold.

**Mode 1 (Stage 1 only):** When `-mode 1` is specified, Stages 2 and 3 are skipped entirely. The `.kpx` file is not accessed, saving I/O and memory. Results contain only Stage 1 scores; position fields (qstart, qend, sstart, send) and chainscore are omitted. The sort key is forced to stage1 score.

**Mode 3 (Full pipeline):** When `-mode 3` is specified, all three stages are executed. A BLAST DB is required (specified via `-db`, defaulting to the index prefix). The sort key is automatically set to alnscore. SAM/BAM output requires `-mode 3` with `-stage3_traceback 1`.

By default, both forward and reverse complement strands of the query are searched. Use `-strand 1` to search only the plus (forward) strand, or `-strand -1` to search only the minus (reverse complement) strand.

### Spaced Seed Template Masks

When spaced seeds are enabled (`-t > 0`), k-mers are extracted using discontiguous megablast-style bitmask templates. Each mask selects k positions from a window of t bases. Two template types are available for each (k, t) combination: **coding** (optimized for coding regions) and **optimal** (optimized for non-coding regions). At index time, `coding` or `optimal` is specified to build a single-template index. At search time, the **both** option merges separate coding and optimal indexes to combine their results.

All templates are derived from the design principles of discontiguous MegaBLAST templates. **Coding** templates follow a periodic "110" structure that maximizes coverage of the second codon position, with excess gaps placed at the first codon position. **Optimal** templates are designed to minimize overlap with the corresponding coding template and use a non-periodic structure; they employ bookend patterns that vary with template length (t=13: "111" at both ends; t=15: "111" at the start + "11" at the end; t=18: "111" or "11" at the start + "11" at the end). The k=11/12 templates are the native discontiguous MegaBLAST templates, while k=8/9 templates are newly designed following the same principles for shorter seed weights suited to PCR amplicon analysis.

Valid (k, t) combinations and their mask values:

| k | t | Type | Mask (binary, left-to-right = pos 0 to t−1) | Mask (hex) |
|---|---|---|---|---|
| 8 | 13 | coding | `1101100101101` | 0x1B2D |
| 8 | 13 | optimal | `1110011010011` | 0x1CD3 |
| 8 | 15 | coding | `100101100101101` | 0x4B2D |
| 8 | 15 | optimal | `111010010010011` | 0x7493 |
| 8 | 18 | coding | `100100100100101101` | 0x2492D |
| 8 | 18 | optimal | `110010001010010011` | 0x32293 |
| 9 | 13 | coding | `1101101101101` | 0x1B6D |
| 9 | 13 | optimal | `1110111010011` | 0x1DD3 |
| 9 | 15 | coding | `100101101101101` | 0x4B6D |
| 9 | 15 | optimal | `111010011010011` | 0x74D3 |
| 9 | 18 | coding | `100100101100101101` | 0x24B2D |
| 9 | 18 | optimal | `111010001010010011` | 0x3A293 |
| 11 | 16 | coding | `1101101101101101` | 0xDB6D |
| 11 | 16 | optimal | `1110010110110111` | 0xE5B7 |
| 11 | 18 | coding | `101101100101101101` | 0x2D96D |
| 11 | 18 | optimal | `111010010110010111` | 0x3A597 |
| 11 | 21 | coding | `100101100101100101101` | 0x12CB2D |
| 11 | 21 | optimal | `111010010100010010111` | 0x1D2897 |
| 12 | 16 | coding | `1111101101101101` | 0xFB6D |
| 12 | 16 | optimal | `1110110110110111` | 0xEDB7 |
| 12 | 18 | coding | `101101101101101101` | 0x2DB6D |
| 12 | 18 | optimal | `111010110010110111` | 0x3ACB7 |
| 12 | 21 | coding | `100101101101100101101` | 0x12DB2D |
| 12 | 21 | optimal | `111010010110010010111` | 0x1D2C97 |

**KmerInt type selection:** The integer type used for k-mer values is determined by the required bit width: `2*k` bits. If the bit width exceeds 16, `uint32_t` is used; otherwise `uint16_t`. This means k <= 8 uses `uint16_t` (16 bits or fewer), while k >= 9 uses `uint32_t`.

### High-Frequency K-mer Filtering

High-frequency k-mer filtering is performed exclusively at index build time via `-max_freq_build`. There is no search-time high-frequency filter — the search hot path no longer pays the cost of per-query k-mer count aggregation across volumes.

**Build-time exclusion** (`-max_freq_build`): When indexing with `-max_freq_build`, high-frequency k-mers are excluded from the index entirely. K-mer counts are aggregated across all volumes before applying the threshold, so a k-mer that is locally below the threshold in each volume but exceeds it globally will be correctly excluded. A single shared `.khx` file (per k value, not per volume) records which k-mers were excluded. When a fractional value (0 < x < 1) is specified, the threshold is resolved using the total NSEQ across all volumes. At search time, when fractional `-stage1_min_score` is used, k-mers excluded at build time are recognized from the `.khx` file and subtracted from the threshold calculation (`Nhighfreq` term in the formula below).

### Fractional Stage 1 Threshold

When `-stage1_min_score` is specified as a fraction (0 < P < 1), the threshold is resolved per query as:

```
threshold = ceil(Nqkmer * P) - Nhighfreq
```

Where:
- **Nqkmer**: pure window count, `max(0, seq_len - span + 1)` where `span = t` for spaced seeds (`-t > 0`) and `span = k` otherwise. Content-independent.
- **Nhighfreq**: number of windows that fail to contribute a usable Stage 1 hit, defined as the union of three cases:
  1. **Build-time exclusion**: a window position where *at least one* expanded k-mer is excluded by `.khx`.
  2. **Degenerate over-expansion**: a window whose IUPAC degenerate-base expansion product exceeds `-max_degen_expand` and is therefore not emitted.
  3. **Truly invalid character in window**: a window containing a character outside `[ACGT]` ∪ IUPAC. (Note: an entire query containing such a character is rejected with `kSkipInvalidChar` before threshold resolution; see "Skip reasons" below.)

  Equivalent computation used by ikafssn: `Nhighfreq = #{positions with ANY-excluded k-mer} + (Nqkmer − #emitted positions)`.

If the resolved threshold falls below 1 on every searched strand, the query is skipped with `kSkipThresholdUnreachable` (see below) instead of producing zero hits silently.

### Skip reasons

When ikafssn cannot run a query through Stages 1–3 it emits a sentinel record so consumers can see why no hits were produced. The sentinel is emitted across every output format:

| Reason (enum / string) | Meaning |
|---|---|
| `kSkipQueryTooShort` / `query_too_short` | `seq_len < span` (k-mer / spaced-seed window cannot fit), or `seq_len < -min_query_length`. |
| `kSkipDegenRejected` / `degen_rejected` | `-accept_qdegen 0` was set and the query contains IUPAC degenerate bases. |
| `kSkipInvalidChar` / `invalid_char` | Query contains a character outside `[ACGT]` ∪ IUPAC. The detail names the offending position. |
| `kSkipThresholdUnreachable` / `threshold_unreachable` | Resolved fractional threshold is `< 1` on every searched strand. |
| `kSkipQueryTooLong` / `query_too_long` | `seq_len > overlap_length` on a fragment-split index.  The parent-relative dedup keys rely on every chain hit fitting inside at most two adjacent fragments, so queries longer than the overlap are rejected up-front rather than producing per-fragment partial chains.  Indexes with `min_length_split == 0` (no splitting) report `overlap_length == 0`, which disables this check. |

Output representations:

- **TSV**: a single row per skipped query with `sseqid = "*SKIPPED:<reason>"`, `sstrand = '*'`, all numeric columns set to 0, and `qlen` populated. `result_reader` silently drops these rows so existing pipelines keep working.
- **JSON**: a query object with `"status": "skipped"`, `"skip_reason": "<reason>"`, `"skip_detail": "<detail>"`, and an empty `"hits": []`. Non-skipped queries get `"status": "ok"`.
- **SAM/BAM**: an unmapped record with `FLAG = 4`, `RNAME = *`, `POS = 0`, `CIGAR = *`, plus aux tags `XR:Z:<reason>` and `XD:Z:<detail>`. Use `samtools view -f 4` to extract them.
- **stderr**: a `Warning: query 'X' skipped: <reason> (<detail>)` line is also emitted at preprocess time.

Searches with at least one skipped query exit with code `2`.

### `-template_type both`: cross-template merging

When `-template_type both` is selected for spaced seeds, ikafssn searches the coding and optimal indexes simultaneously and merges Stage 1 scores into a single, unified per-position score in `[0, Nqkmer]`. This avoids two failure modes of a naïve per-template scoring:

1. **Cross-template double-counting**: a sequence position that matches both the coding and the optimal mask would otherwise contribute `+2` to a query position's score, inflating the dynamic range to `[0, 2 × Nqkmer]`. ikafssn walks both streams in `q_pos` order through a shared per-(sid, q_pos) deduplication buffer, so each query position contributes at most `+1`.
2. **Hybrid-query miss**: a query whose region A only matches the coding template and region B only matches the optimal template would fail an additive `thr_cod + thr_opt` threshold. ikafssn instead applies a unified threshold of `min(thr_cod, thr_opt)` so hybrid matches are retained.

The resolved Stage 1 score reported for `-template_type both` therefore behaves the same as for `coding` or `optimal` alone: it is the count of query window positions matched on the subject (after build-time `.khx` exclusion), bounded by `Nqkmer`. `Nhighfreq` is computed independently per side and used as the threshold offset.

### Score Types

ikafssn computes three types of scores:

| Score | Description | Computed in |
|---|---|---|
| **coverscore** | Number of query k-mer positions that match the subject sequence. Each query position is counted at most once per subject, however many positions of the subject its k-mer occurs at. This is the Stage 1 score. | Stage 1 |
| **chainscore** | Length (number of k-mer hits) of the best collinear chain found by the chaining DP. Requires position data from `.kpx`. | Stage 2 |
| **alnscore** | Semi-global pairwise alignment score computed by Parasail (using `-stage3_score_matrix`, default: degmatch). Requires BLAST DB for subject sequence retrieval. | Stage 3 |

- The Stage 1 score is always coverscore.
- The sort key is auto-determined by mode: mode 1 sorts by coverscore, mode 2 by chainscore, mode 3 by alnscore.
- In `-mode 1`, only Stage 1 scores are available; chainscore and alnscore are not computed.

### Stage 3 Scoring Matrix

Stage 3 pairwise alignment uses a configurable scoring matrix specified by `-stage3_score_matrix`. Three matrices are available:

| Matrix | Description |
|---|---|
| **degmatch** (default) | Assigns positive scores to degenerate base pairs that share at least one possible nucleotide identity. Suitable for primer matching and searches involving degenerate bases. |
| **dnafull** | Traditional EMBOSS DNA full matrix (created by Todd Lowe). Uses ambiguous nucleotide codes with probabilities rounded to nearest integer. Degenerate pairs that share partial overlap receive negative or low scores. |
| **nuc44** | NCBI BLAST nucleotide matrix. Similar to dnafull but with slightly different degenerate base scoring. |

**CIGAR `=`/`X` determination:** The extended CIGAR operators `=` (sequence match) and `X` (sequence mismatch) are determined based on alignment scores, not strict base identity. A position is reported as `=` when its score in the matrix is positive (score > 0), and as `X` when its score is zero or negative. This means that with the DEGMATCH matrix, a degenerate base pair like N-A (score 1) is counted as a match (`=`), while with NUC44/DNAFULL it would be counted as a mismatch (`X`).

**Column naming:** The output columns are `ppositive` (percent positive-scoring), `npositive` (number of positive-scoring positions), and `nnegative` (number of negative-scoring positions). With the DEGMATCH matrix these counts represent positive-scoring positions rather than strictly identical positions. The matching filter options are `-stage3_min_ppositive` and `-stage3_min_npositive`.

**Aligned region only:** Stage 3 aligns the query against a context-extended subject window, so the raw alignment carries terminal gaps that span the rest of that window. Those are stripped before reporting: `qstart` / `qend` / `sstart` / `send`, `cigar`, `qseq` / `sseq` and the `ppositive` denominator all describe the aligned region between the terminal gaps, matching what BLAST reports for the same alignment. `ppositive` is therefore `npositive` over the aligned length, not over the window length.

**DEGMATCH matrix (16x16):**

The full DEGMATCH scoring matrix. Rows and columns correspond to the 16 symbols: A, T, G, C, S (G/C), W (A/T), R (A/G), Y (T/C), K (G/T), M (A/C), B (G/T/C), V (A/G/C), H (A/T/C), D (A/G/T), N (A/T/G/C), and `*` (stop/invalid).

|   | A | T | G | C | S | W | R | Y | K | M | B | V | H | D | N | * |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 5 | -4 | -4 | -4 | -4 | 3 | 3 | -4 | -4 | 3 | -4 | 2 | 2 | 2 | 1 | -4 |
| **T** | -4 | 5 | -4 | -4 | -4 | 3 | -4 | 3 | 3 | -4 | 2 | -4 | 2 | 2 | 1 | -4 |
| **G** | -4 | -4 | 5 | -4 | 3 | -4 | 3 | -4 | 3 | -4 | 2 | 2 | -4 | 2 | 1 | -4 |
| **C** | -4 | -4 | -4 | 5 | 3 | -4 | -4 | 3 | -4 | 3 | 2 | 2 | 2 | -4 | 1 | -4 |
| **S** | -4 | -4 | 3 | 3 | 4 | -4 | 1 | 1 | 1 | 1 | 2 | 2 | 1 | 1 | 1 | -4 |
| **W** | 3 | 3 | -4 | -4 | -4 | 4 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 1 | -4 |
| **R** | 3 | -4 | 3 | -4 | 1 | 1 | 4 | -4 | 1 | 1 | 1 | 2 | 1 | 2 | 1 | -4 |
| **Y** | -4 | 3 | -4 | 3 | 1 | 1 | -4 | 4 | 1 | 1 | 2 | 1 | 2 | 1 | 1 | -4 |
| **K** | -4 | 3 | 3 | -4 | 1 | 1 | 1 | 1 | 4 | -4 | 2 | 1 | 1 | 2 | 1 | -4 |
| **M** | 3 | -4 | -4 | 3 | 1 | 1 | 1 | 1 | -4 | 4 | 1 | 2 | 2 | 1 | 1 | -4 |
| **B** | -4 | 2 | 2 | 2 | 2 | 1 | 1 | 2 | 2 | 1 | 3 | 1 | 1 | 1 | 1 | -4 |
| **V** | 2 | -4 | 2 | 2 | 2 | 1 | 2 | 1 | 1 | 2 | 1 | 3 | 1 | 1 | 1 | -4 |
| **H** | 2 | 2 | -4 | 2 | 1 | 2 | 1 | 2 | 1 | 2 | 1 | 1 | 3 | 1 | 1 | -4 |
| **D** | 2 | 2 | 2 | -4 | 1 | 2 | 2 | 1 | 2 | 1 | 1 | 1 | 1 | 3 | 1 | -4 |
| **N** | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | -4 |
| **\*** | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 | -4 |

Score rules:
- Standard match (A-A, T-T, G-G, C-C) = 5
- Standard mismatch (e.g. A-T, A-G) = -4
- 2-fold degenerate self (S-S, W-W, R-R, Y-Y, K-K, M-M) = 4
- 2-fold degenerate vs contained base (e.g. R-A, R-G) = 3
- 2-fold degenerate vs partial overlap with another 2-fold (e.g. R-K) = 1
- 2-fold degenerate vs non-overlapping (e.g. S-W) = -4
- 3-fold degenerate self (B-B, V-V, H-H, D-D) = 3
- 3-fold degenerate vs contained base/2-fold (e.g. B-T, B-K) = 2
- 3-fold degenerate vs partial overlap = 1
- N vs any (except `*`) = 1, N-N = 1
- `*` vs any = -4

## Output Format

**Coordinate convention:** TSV and JSON output follows BLAST `-outfmt 6`.

- All four of `qstart` / `qend` / `sstart` / `send` are **1-based and inclusive**.
- `qstart` / `qend` are **query-relative** — positions in the query sequence as given — so `qstart` < `qend` on both strands.
- `sstart` / `send` are **parent-relative** in forward numbering but run in alignment order, so a reverse-strand hit (`sstrand` = `-`) has **`sstart` > `send`**.

`sseqid` is the **parent OID's** accession (never a fragment-derived synthetic name) and `slen` is the parent OID's full length. After Stage 2 / Stage 3 dedup folds per-fragment chains together, every row in the output describes one canonical hit per `(qseqid, sseqid, sstrand, send, alnscore)` tuple.

SAM / BAM output does not use this convention: its `POS` field is a 0-based leftmost coordinate carried by htslib and is never swapped by strand (see [SAM/BAM Format](#sambam-format)).

### Tab Format (default)

**Mode 2** (default):

Tab-separated columns (the Stage 1 score column is always `coverscore`):

```
# qseqid  sseqid  sstrand  qstart  qend  qlen  sstart  send  slen  coverscore  chainscore  volume
```

**Mode 1** (`-mode 1`):

```
# qseqid  sseqid  sstrand  qlen  slen  coverscore  volume
```

**Mode 3, traceback=0** (`-mode 3`):

```
# qseqid  sseqid  sstrand  qend  qlen  send  slen  coverscore  chainscore  alnscore  volume
```

Note: `qstart` and `sstart` are omitted because accurate alignment start positions are unavailable without traceback. On the reverse strand the alignment's end is the *lower* subject coordinate, so this layout's `send` is the smaller of the two subject positions the traceback layout would report.

**Mode 3, traceback=1** (`-mode 3 -stage3_traceback 1`):

```
# qseqid  sseqid  sstrand  qstart  qend  qlen  sstart  send  slen  coverscore  chainscore  alnscore  ppositive  npositive  nnegative  cigar  qseq  sseq  volume
```

### JSON Format

`results` holds one object per `qseqid`, in the order the queries first appear in the input, and each object lists that query's hits in the order the search produced them. Query sequences that share a `qseqid` are searched independently — each gets its own hits — but they share one object here, since the object is keyed by the name.

**Mode 2** (default):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qstart": 1,
          "qend": 150,
          "qlen": 200,
          "sstart": 1001,
          "send": 1150,
          "slen": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "volume": 0
        }
      ]
    }
  ]
}
```

**Mode 1** (`-mode 1`):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qlen": 200,
          "slen": 5000,
          "coverscore": 8,
          "volume": 0
        }
      ]
    }
  ]
}
```

**Mode 3, traceback=0** (`-mode 3`):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qend": 150,
          "qlen": 200,
          "send": 1150,
          "slen": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "alnscore": 240,
          "volume": 0
        }
      ]
    }
  ]
}
```

**Mode 3, traceback=1** (`-mode 3 -stage3_traceback 1`):

```json
{
  "results": [
    {
      "qseqid": "query1",
      "hits": [
        {
          "sseqid": "NC_001234.5",
          "sstrand": "+",
          "qstart": 1,
          "qend": 150,
          "qlen": 200,
          "sstart": 1001,
          "send": 1150,
          "slen": 5000,
          "coverscore": 8,
          "chainscore": 12,
          "alnscore": 240,
          "ppositive": 95.3,
          "npositive": 143,
          "nnegative": 7,
          "cigar": "50=2X48=1I50=",
          "qseq": "ACGT...",
          "sseq": "ACGT...",
          "volume": 0
        }
      ]
    }
  ]
}
```

### SAM/BAM Format

SAM/BAM output requires `-mode 3 -stage3_traceback 1`. Use `-output_format sam` for SAM or `-output_format bam` for BAM (BAM requires `-o <path>`).

SAM records contain:
- **QNAME**: qseqid
- **FLAG**: 0 (forward) or 16 (reverse)
- **RNAME**: sseqid
- **POS**: leftmost aligned parent-relative position, 1-based. This is the alignment's *lower* subject coordinate on both strands — unlike the TSV / JSON `sstart`, which swaps with `send` on the reverse strand.
- **MAPQ**: 255
- **CIGAR**: extended CIGAR with =/X/I/D operators, always in reference orientation. Query bases outside the aligned region are reported as leading / trailing hard clips (`H`) because SEQ carries the aligned region only.
- **SEQ**: ungapped aligned query sequence in reference orientation, i.e. the reverse complement of the given query when FLAG 0x10 is set, as the SAM specification requires
- **QUAL**: * (not available)
- **Tags**: `AS:i` (alnscore), `NM:i` (nnegative), `cs:i` (chainscore), `cv:i` (coverscore)

## Deployment Architecture

### Single Machine

```
ikafssnsearch (standalone, no server needed)
    or
ikafssnserver → ikafssnclient (via UNIX socket)
```

### Multi-Machine

```
Machine A (search server):
  ikafssnserver -ix ./index/mydb -tcp 0.0.0.0:9100

Machine B (HTTP frontend):
  ikafssnhttpd -server_tcp A:9100 -listen :8080
  nginx (TLS, auth, rate limiting) → ikafssnhttpd
```

### Multiple Databases

A single `ikafssnserver` process can serve multiple databases simultaneously:

```
# Single process, multiple DBs (recommended)
ikafssnserver -ix ./nt_index -db nt -ix ./rs_index -db refseq_genomic \
    -socket /var/run/ikafssn.sock

ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta
ikafssnclient -socket /var/run/ikafssn.sock -ix refseq_genomic -query query.fasta
```

### Multi-Backend (Load Balancing / Multi-Server)

A single `ikafssnhttpd` can connect to multiple `ikafssnserver` instances:

```
# Load balancing: same DB on two servers, one httpd
ikafssnserver -ix ./nt_index -db nt -tcp 0.0.0.0:9100   # Server A
ikafssnserver -ix ./nt_index -db nt -tcp 0.0.0.0:9100   # Server B
ikafssnhttpd -server_tcp A:9100 -server_tcp B:9100 -listen :8080

# Different DBs on separate servers, unified through one httpd
ikafssnserver -ix ./nt_index -db nt -socket /var/run/nt.sock
ikafssnserver -ix ./rs_index -db refseq -socket /var/run/rs.sock
ikafssnhttpd -server_socket /var/run/nt.sock -server_socket /var/run/rs.sock -listen :8080
```

When the same database name appears on multiple backends, `ikafssnhttpd` verifies at startup that for each shared (db, k) pair, total sequence counts and total bases are identical. Backends may have different k-value sets for the same database; the merged capabilities expose the union of all k-value groups, so requests for a k-value available on only some backends are naturally routed to those backends. Requests are routed to the highest-priority backend with available effective capacity (considering both slot availability and `max_nseq_per_req`). Note that capacity values (`max_queue_size`, `queue_depth`, `max_nseq_per_req`) are shared per server across all databases served by that server.

Alternatively, separate processes with path-based HTTP routing:

```
ikafssnserver -ix ./nt_index  -socket /var/run/ikafssn_nt.sock
ikafssnserver -ix ./rs_index  -socket /var/run/ikafssn_rs.sock

ikafssnhttpd -server_socket /var/run/ikafssn_nt.sock -listen :8080 -path_prefix /nt
ikafssnhttpd -server_socket /var/run/ikafssn_rs.sock -listen :8081 -path_prefix /rs
# nginx routes /nt → :8080, /rs → :8081
```

### systemd Integration

Sample systemd unit files are provided in `doc/systemd/`. See the files for configuration details.

## OS settings for large databases

A BLAST DB that is split into many volumes can need more open file descriptors than the default limit allows. An open volume costs **two file descriptors** (`.nin` + `.nsq`), and the NCBI toolkit holds them for as long as the volume's memory mapping lives. NCBI `nt`, for example, is several hundred volumes.

Which commands need how many:

| Command | Volumes open at once | Descriptors needed |
|---|---|---|
| `ikafssnindex` | one | a small constant (under 10) |
| `ikafssnretrieve` | one | a small constant (under 10) |
| `ikafssnsearch -mode 1` / `-mode 2` | none | a small constant |
| `ikafssnsearch -mode 3` | all | `2 × volumes + 64` |
| `ikafssnserver` (with `-db`) | all | `2 × (volumes across all databases) + 256` |
| `ikafssnhttpd` | none (opens no BLAST DB) | a small constant |

`ikafssnsearch -mode 3` and `ikafssnserver` raise their own `RLIMIT_NOFILE` soft limit to the figure above before they open anything. That succeeds without any privilege as long as the **hard** limit is high enough; when it is not, they report the shortfall and exit. Raising the hard limit is an OS-level setting:

**Interactive shell** — applies to commands started from that shell:

```bash
ulimit -n 4096          # soft limit
ulimit -Hn 4096         # hard limit (lowering it is irreversible for the shell)
```

**All logins on the machine** — `/etc/security/limits.d/99-ikafssn.conf` (read by PAM at login):

```
*  soft  nofile  4096
*  hard  nofile  65536
```

**systemd services** — `ulimit` and `limits.d` do not apply to units; set the limit in the unit's `[Service]` section (or via a drop-in under `/etc/systemd/system/<unit>.d/`):

```ini
[Service]
LimitNOFILE=65536
```

**Kernel ceilings** — `fs.nr_open` caps what any single process may be given, and `fs.file-max` caps the whole system. Raise them in `/etc/sysctl.d/99-ikafssn.conf` if the per-process limit you need exceeds them:

```
fs.nr_open = 1048576
fs.file-max = 2097152
```

## Index File Formats

All fragment-indexing parameters that change the *content* of an index
are encoded directly in the file name, so multiple indexes built with
different parameters can co-exist in the same output directory.

Per BLAST DB volume, three files are generated:

```
<vol_basename>.k<k>[.t<t>][.minlen<X>][.minsplit<X>][.ovllen<X>]
              [.maxfreq<X>][.maxexpand<X>][.cod|.opt].(kix|kpx|ksx)
```

The shared per-DB manifest (and, when `-max_freq_build` is used, the
exclusion bitset) follow the same convention without the volume index:

```
<db_base>.k<k>[.t<t>][.minlen<X>][.minsplit<X>][.ovllen<X>]
          [.maxfreq<X>][.maxexpand<X>][.cod|.opt].(kvx|khx)
```

Suffix-omission rules:

| Suffix | Omitted when |
|---|---|
| `k<X>`         | never (always emitted, no zero-padding) |
| `t<X>`         | `t == 0` |
| `minlen<X>`    | `min_seq_length == 0` |
| `minsplit<X>`  | `min_length_split == 0` |
| `ovllen<X>`    | `overlap_length == 0` (synchronized with `minsplit`) |
| `maxfreq<X>`   | `max_freq_build == 1` (the resolved absolute threshold) |
| `maxexpand<X>` | `max_degen_expand` is `0` or `1` |
| `cod` / `opt`  | `t == 0` |

Examples:

- Default `-mode 1` (`-min_seq_length 64 -min_length_split 50000 -overlap_length 500`):
  - `nt.00.k11.minlen64.minsplit50000.ovllen500.kix`
  - `nt.k11.minlen64.minsplit50000.ovllen500.kvx`
- Mode 1 with cross-volume filter (`-max_freq_build 1000 -max_degen_expand 4`):
  - `nt.00.k11.minlen64.minsplit50000.ovllen500.maxfreq1000.maxexpand4.kix`
  - `nt.k11.minlen64.minsplit50000.ovllen500.maxfreq1000.maxexpand4.khx`
- Mode 2 / 3 (`-min_length_split 0 -overlap_length 0`):
  - `nt.00.k11.minlen64.kix`
- Spaced seed (`-k 11 -t 16 -template_type coding -mode 1`):
  - `nt.00.k11.t16.minlen64.minsplit50000.ovllen500.cod.kix`
  - `nt.k11.t16.minlen64.minsplit50000.ovllen500.cod.kvx`

Fractional `-max_freq_build` (e.g. `0.001`) is resolved to an absolute
threshold at the end of the metadata-collection pass, and that
absolute value is what appears in the file name, so the same fraction
applied to a different fragment set produces a different file name.

The `.khx` file contains a 32-byte header (magic "KMHX", format version,
k) followed by a bitset of `ceil(4^k / 8)` bytes. Bit *i* = 1 indicates
that k-mer *i* was excluded during index build based on cross-volume
aggregated counts.

ID and position postings are stored in separate files so that Stage 1
filtering never touches `.kpx`, maximizing page cache efficiency.

**Multi-accession deflines:** When the source BLAST DB was built with `makeblastdb -parse_seqids` and carries multi-defline records (the NCBI convention for registering identical sequences under several accessions, separated by `\x01` / `^A` in the FASTA defline), `ikafssnindex` preserves **all** accessions for each OID. The `.ksx` accession string for such OIDs contains every accession joined by `\x01`, and search output emits the same `\x01`-joined string in the `sseqid` column / SAM RNAME / FASTA defline / protocol `sseqid` field. Downstream consumers should split on `\x01` to recover individual accessions. The `-seqidlist` filter and `ikafssnretrieve` accept either the full `\x01`-joined form or any individual constituent accession.

**Index format version:** The current index format is **v11** for every index file (`.kix`, `.kpx`, `.ksx`, `.khx`, `.kvx`). Layout summary:

- **`.kix` / `.kpx` magic** is `KIX11` / `KPX11`; the fragment-indexing parameters that used to live in the headers (`min_seq_length`, `min_length_split`, `overlap_length`, plus the resolved `max_freq_build` / `max_degen_expand`) are now encoded in the file name and parsed once at volume discovery time.
- **`.ksx` two-stage layout:** the sequence-metadata file records each parent OID's `(parent_length, blast_oid, accession)` in a parent table, followed by a fragment table that maps every internal SeqId to `(parent_idx, fragment_start, fragment_end)`. Convenience accessors (`KsxReader::seq_length` / `accession`) take a SeqId and resolve to the matching parent. Magic is `KMSX`.
- **`.kix` body:** Elias-Fano dictionary; each posting list is `[u32 distinct_count]` followed by a FastPFor `CompositeCodec<SIMDFastPFor<4>, VariableByte>` body over the distinct seq_id delta stream.
- **`.kpx` body:** Elias-Fano dictionary; each posting list starts directly at a 2-bit kind map (no per-list header) and continues with the partition / short_occ1 / short_occ_ge2 sub-buckets, all bit-packed as frame-of-reference blocks. The decoder is driven by the candidate set, so only the blocks holding a wanted position are unpacked.
- **`.kvx`:** manifest text format with a `FORMAT_VERSION` line set to 11.

The fragment splitter is a port of kafssstore's `split_long_sequence_positions` (ncbi4na==0xF cuts, calcsegment2 formula). When `-min_length_split 0`, every parent is registered as a single fragment that spans the whole parent — i.e. one fragment per OID, and the file name carries no `minsplit` / `ovllen` suffix.

Indexes whose `format_version` does not match are rejected at open with a message such as:

```
KixReader: index format version mismatch (got 10, expected 11). Please rebuild with the current ikafssnindex.
```

(and analogous messages for `.kpx` / `.ksx` / `.khx` / `.kvx`). Rebuild with the current `ikafssnindex`:

```
# Default (mode 1, fragment-aware)
ikafssnindex -db <BLAST_DB_prefix> -k <k> -o <out_dir>

# Disable splitting (one fragment per parent)
ikafssnindex -db <BLAST_DB_prefix> -k <k> -o <out_dir> \
    -min_length_split 0 -overlap_length 0
```

For NCBI nt-scale databases (~700 volumes at k=12 t=21), expect tens of hours of rebuild time on a 32-core host. With `-min_length_split` non-zero the on-disk size grows in proportion to the average overlap length per parent, and chromosome-scale parents are split into multiple fragments of bounded length — which bounds Stage 2's per-subject memory.

## Installation

### macOS 26 Tahoe (Homebrew)

On macOS 26 (Tahoe) with Apple Silicon (aarch64), install via the Homebrew Tap:

```bash
brew tap astanabe/ikafssn
brew install ikafssn
```

This installs pre-built bottles when available, or builds from source as a fallback.

### Conda (channel)

ikafssn is distributed through a dedicated Conda channel hosted at `https://conda.ikafssn.org`. Pre-built `.conda` packages are available for `linux-64`, `linux-aarch64`, and `osx-arm64`. The Linux variants target glibc 2.34 or newer (RHEL 9, Ubuntu 22.04, Debian 12 and later); the `osx-arm64` variant targets macOS 11 or newer.

Install into an existing miniforge / miniconda environment with `mamba` (or `conda`); runtime dependencies are pulled from `conda-forge`:

```bash
mamba install -c https://conda.ikafssn.org -c conda-forge ikafssn
```

Or to create a dedicated environment:

```bash
mamba create -n ikafssn -c https://conda.ikafssn.org -c conda-forge ikafssn
```

The channel itself only hosts `repodata.json` index files (via GitHub Pages on [astanabe/conda-ikafssn](https://github.com/astanabe/conda-ikafssn)); the `.conda` binaries are streamed directly from this project's [GitHub Releases](https://github.com/astanabe/ikafssn/releases) page using the [CEP-15 `base_url`](https://github.com/conda/ceps/blob/main/cep-0015.md) mechanism.

Requirements:
- `conda` 23.7+ or `micromamba` 2.0+ (older clients silently ignore `base_url` and fail with HTTP 404 on package download).

### Ubuntu (APT channel)

ikafssn is distributed through a signed APT channel hosted at `https://deb.ikafssn.org`. The channel always reflects the latest release; older versions are not retained.

| Suite | Ubuntu release | Architectures |
|---|---|---|
| `jammy` | Ubuntu 22.04 LTS | `amd64`, `arm64` |
| `noble` | Ubuntu 24.04 LTS | `amd64`, `arm64` |
| `resolute` | Ubuntu 26.04 LTS | `amd64`, `arm64` |

One-time setup (run as root or via `sudo`):

```bash
sudo install -d -m 0755 /etc/apt/keyrings
curl -fsSL https://deb.ikafssn.org/ikafssn-archive-keyring.asc \
  | sudo gpg --dearmor -o /etc/apt/keyrings/ikafssn-archive-keyring.gpg
echo "deb [signed-by=/etc/apt/keyrings/ikafssn-archive-keyring.gpg] https://deb.ikafssn.org/ $(lsb_release -cs) main" \
  | sudo tee /etc/apt/sources.list.d/ikafssn.list
sudo apt update
sudo apt install ikafssn
```

Subsequent upgrades:

```bash
sudo apt update && sudo apt upgrade ikafssn
```

The channel hosts the `.deb` binaries themselves and the per-suite `Packages` / `Release` / `InRelease` / `Release.gpg` metadata. The public key fingerprint can also be verified out-of-band against `ikafssn-archive-keyring.asc` committed at the root of [astanabe/ikafssn](https://github.com/astanabe/ikafssn).

### Enterprise Linux (DNF / YUM channel)

ikafssn is distributed through a signed DNF channel hosted at `https://rpm.ikafssn.org`. The channel always reflects the latest release; older versions are not retained.

| `releasever` | Distribution | Architectures |
|---|---|---|
| `9` | AlmaLinux / Rocky Linux / RHEL 9 | `x86_64`, `aarch64` |
| `10` | AlmaLinux / Rocky Linux / RHEL 10 | `x86_64`, `aarch64` |

One-time setup (run as root or via `sudo`):

```bash
sudo tee /etc/yum.repos.d/ikafssn.repo > /dev/null <<'EOF'
[ikafssn]
name=ikafssn
baseurl=https://rpm.ikafssn.org/el$releasever/$basearch/
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://rpm.ikafssn.org/ikafssn-archive-keyring.asc
EOF
sudo dnf install ikafssn
```

Subsequent upgrades:

```bash
sudo dnf upgrade ikafssn
```

Both the package signature (`gpgcheck=1`) and the repository metadata signature (`repo_gpgcheck=1`) are verified against the channel's public key.

### Direct download (fallback)

If installing through Homebrew / Conda / APT / DNF is not an option, the same `.deb`, `.rpm`, `.conda`, and bottle artifacts can be downloaded directly from the [GitHub Releases](https://github.com/astanabe/ikafssn/releases) page.

Package naming conventions:

```
ikafssn_<version>_ubuntu-<ubuntu_ver>_<arch>.deb         # Ubuntu 22.04 / 24.04 / 26.04, amd64 / arm64
ikafssn-<version>.el<el_ver>.<arch>.rpm                   # EL 9 / 10, x86_64 / aarch64
ikafssn_<version>_<conda_subdir>.conda                    # linux-64 / linux-aarch64 / osx-arm64
ikafssn-<version>.arm64_tahoe.bottle.tar.gz               # macOS 26 (Tahoe) arm64 Homebrew bottle
```

Install a downloaded `.deb`:

```bash
sudo apt install ./ikafssn_<version>_ubuntu-<ubuntu_ver>_<arch>.deb
```

Install a downloaded `.rpm`:

```bash
sudo dnf install ./ikafssn-<version>.el<el_ver>.<arch>.rpm
```

### Verify installation

```bash
ikafssnindex -version
```

## Building from Source

### Dependencies

- C++20 compiler (GCC >= 11, Clang >= 13)
- CMake >= 3.16
- NCBI C++ Toolkit (for BLAST DB access)
- Intel TBB (for parallelization)
- Parasail >= 2.6 (for Stage 3 pairwise alignment)
- htslib >= 1.17 (for SAM/BAM output)
- zlib, libbz2, liblzma, libzstd (for transparent compressed I/O); pkg-config required for libzstd discovery
- Drogon (for ikafssnhttpd, optional)
- libcurl (for HTTP client mode and remote retrieval, optional)

CMake downloads FastPFor (the `.kix` / `.kpx` posting list codec) and, on x86_64,
x86-simd-sort at configure time, so the build host needs network access unless the
CMake FetchContent cache is already populated.  Nothing else has to be installed for them.

### CPU requirements

- **x86_64**: SSE4.2 minimum (Intel Nehalem 2008+, AMD Bulldozer 2011+). The runtime SIMD dispatcher additionally targets AVX2, AVX-512 BW, and AVX-512 VBMI2 when present. CPUs without SSE4.2 are rejected at startup with `exit(2)`.
- **aarch64**: NEON (ASIMD) minimum (Armv8.0+). Pre-NEON aarch64 CPUs are rejected at startup with `exit(2)`. SVE / SVE2 capable CPUs use the NEON-tier FastPFor codec object (the dispatcher routes any tier ≥ NEON to a single NEON OBJECT library, built on FastPFor's native NEON header); per-kernel SIMD files (toupper, ncbi2na unpack, k-mer revcomp, degenerate scan, spaced-seed) keep their separate SVE / SVE2 dispatches.

### Installing Dependencies

Install the required packages (excluding NCBI C++ Toolkit) with the following commands.

**Ubuntu Server 22.04 / 24.04 / 26.04:**

```bash
sudo apt install build-essential cmake pkg-config libtbb-dev liblmdb-dev libsqlite3-dev \
    libcurl4-openssl-dev libjsoncpp-dev
sudo apt install zlib1g-dev libbz2-dev liblzma-dev libzstd-dev libdeflate-dev autoconf \
    libssl-dev uuid-dev
```

`pkg-config` is required by CMake's `pkg_check_modules` call for libzstd. The second line installs dependencies required for building Parasail and htslib from source, plus `libssl-dev` and `uuid-dev` which are needed for building the NCBI C++ Toolkit and Drogon from source. If ikafssnhttpd is not needed, build with `-DBUILD_HTTPD=OFF` and `uuid-dev` may be omitted (`libssl-dev` is still required by the NCBI C++ Toolkit).

**AlmaLinux 9 / Rocky Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y epel-release
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y oracle-epel-release-el9
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**AlmaLinux 10 / Rocky Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ pkgconfig tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libzstd-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

On EL9, `jsoncpp-devel` requires EPEL and `lmdb-devel` requires CRB. On EL10, both are in CRB so EPEL is not needed for these packages. `pkgconfig` is required by CMake's `pkg_check_modules` call for libzstd. The second-to-last line of each block installs dependencies required for building Parasail and htslib from source. The last line installs dependencies needed to build the NCBI C++ Toolkit and Drogon from source. If ikafssnhttpd is not needed, build with `-DBUILD_HTTPD=OFF` and `libuuid-devel` may be omitted (`openssl-devel` is still required by the NCBI C++ Toolkit).

**macOS (Homebrew):**

```bash
brew install cmake pkg-config tbb lmdb sqlite3 curl jsoncpp \
    xz zstd libdeflate autoconf automake libtool openssl@3
```

`openssl@3` is required for building the NCBI C++ Toolkit and Drogon from source (Drogon is built separately, see below). On macOS, use `make -j$(sysctl -n hw.ncpu)` instead of `make -j$(nproc)` in the build steps below.

### Parasail

ikafssn uses the Parasail library for Stage 3 pairwise alignment. By default, CMake looks for Parasail at `./parasail` relative to the source root. If Parasail is installed elsewhere, specify the path with `-DPARASAIL_DIR`.

To download, build, and install Parasail into `./parasail`, run the following from the ikafssn source root:

```bash
curl -L -o parasail-2.6.2.tar.gz \
    https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz
tar xf parasail-2.6.2.tar.gz
cd parasail-2.6.2
patch -p1 < ../patches/parasail-degmatch-cigar-score.patch
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/parasail" \
    -DBUILD_SHARED_LIBS=OFF
make -j$(nproc)
make install
cd ../..
```

On macOS with CMake >= 4.0, add `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` to the `cmake` command above.

### htslib

ikafssn uses htslib for SAM/BAM output. By default, CMake looks for htslib at `./htslib` relative to the source root. If htslib is installed elsewhere, specify the path with `-DHTSLIB_DIR`.

To download, build, and install htslib into `./htslib`, run the following from the ikafssn source root:

```bash
curl -L -o htslib-1.23.1.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2
tar xf htslib-1.23.1.tar.bz2
cd htslib-1.23.1
autoreconf -i
./configure --prefix="$(realpath ..)/htslib" --disable-libcurl --disable-gcs --disable-s3
make -j$(nproc)
make install
cd ..
```

On macOS, add Homebrew paths so that configure can find xz (lzma) and libdeflate:

```bash
./configure --prefix="$(realpath ..)/htslib" --disable-libcurl --disable-gcs --disable-s3 \
    CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib"
```

### NCBI C++ Toolkit

ikafssn requires a pre-built NCBI C++ Toolkit. By default, CMake looks for the toolkit at `./ncbi-cxx-toolkit` relative to the source root. If the toolkit is installed elsewhere, specify the path with `-DNCBI_TOOLKIT_DIR`.

The build subdirectory name within the toolkit (e.g. `CMake-GCC1330-Release`) is auto-detected by default but can be overridden with `-DNCBI_TOOLKIT_BUILD_TAG` if needed.

To download, build, and install the toolkit into `./ncbi-cxx-toolkit`, run the following from the ikafssn source root:

```bash
curl -L -o ncbi-cxx-toolkit-public-release-30.2.0.tar.gz \
    https://github.com/ncbi/ncbi-cxx-toolkit-public/archive/refs/tags/release/30.2.0.tar.gz
tar xf ncbi-cxx-toolkit-public-release-30.2.0.tar.gz
cd ncbi-cxx-toolkit-public-release-30.2.0
patch -p1 < ../patches/ncbi-cxx-toolkit-seqdb-madvise-random.patch
export CXXFLAGS="-std=c++20"
./cmake-configure \
    --without-debug \
    --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" \
    --with-install="$(realpath ..)/ncbi-cxx-toolkit"
cd CMake-GCC*/build
make -j$(nproc)
make install
cd ../../..
```

Only the libraries required by ikafssn (`seqdb`, `blastdb_format`, and their dependencies) are built. The full toolkit build is not necessary.

On macOS, the Homebrew include path must be visible to the compiler (for `lmdb.h`), and the build directory glob pattern differs:

```bash
export CFLAGS="-I$(brew --prefix)/include"
export CXXFLAGS="-I$(brew --prefix)/include -std=c++20"
patch -p1 < ../patches/ncbi-cxx-toolkit-seqdb-madvise-random.patch
./cmake-configure \
    --without-debug \
    --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" \
    --with-install="$(realpath ..)/ncbi-cxx-toolkit"
cd CMake-Clang*/build
make -j$(sysctl -n hw.ncpu)
make install
cd ../../..
```

### Drogon

ikafssnhttpd uses the Drogon HTTP framework. By default, CMake looks for Drogon at `./drogon` relative to the source root. If Drogon is installed elsewhere, specify the path with `-DDROGON_DIR`. If `-DBUILD_HTTPD=OFF` is used, Drogon is not required.

The Drogon release tarball does not include the trantor submodule contents, so the corresponding trantor release must be downloaded separately (Drogon 1.9.12 references trantor v1.5.26). To download, build, and install Drogon and trantor into `./drogon`, run the following from the ikafssn source root:

```bash
curl -L -o drogon-1.9.12.tar.gz \
    https://github.com/drogonframework/drogon/archive/refs/tags/v1.9.12.tar.gz
tar xf drogon-1.9.12.tar.gz
curl -L -o trantor-1.5.26.tar.gz \
    https://github.com/an-tao/trantor/archive/refs/tags/v1.5.26.tar.gz
tar xf trantor-1.5.26.tar.gz
rmdir drogon-1.9.12/trantor
mv trantor-1.5.26 drogon-1.9.12/trantor
cd drogon-1.9.12
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/drogon" \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_CTL=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_ORM=OFF \
    -DBUILD_BROTLI=OFF \
    -DBUILD_YAML_CONFIG=OFF \
    -DBUILD_DOC=OFF
make -j$(nproc)
make install
cd ../..
```

Disabling ORM (PostgreSQL/MySQL/SQLite/Redis), brotli, yaml-cpp, `drogon_ctl`, examples, and documentation generation reduces the link-time dependency surface to just jsoncpp, libuuid, zlib, and OpenSSL.

On macOS, Homebrew's `openssl@3` is keg-only, so `OPENSSL_ROOT_DIR` must be specified explicitly:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$(realpath ../..)/drogon" \
    -DOPENSSL_ROOT_DIR="$(brew --prefix)/opt/openssl@3" \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DBUILD_CTL=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_ORM=OFF \
    -DBUILD_BROTLI=OFF \
    -DBUILD_YAML_CONFIG=OFF \
    -DBUILD_DOC=OFF
make -j$(sysctl -n hw.ncpu)
make install
```

### Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make test
```

If the NCBI C++ Toolkit is installed at a non-default location:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DNCBI_TOOLKIT_DIR=/path/to/ncbi-cxx-toolkit
```

### Installation

```bash
sudo make install
```

By default, executables are installed to `/usr/local/bin` and systemd unit files to `/usr/local/share/ikafssn/systemd`. To change the install prefix:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/ikafssn
make -j$(nproc)
sudo make install
```

In this example, executables are installed to `/opt/ikafssn/bin`.

### CMake Options

| Option | Default | Description |
|---|---|---|
| `NCBI_TOOLKIT_DIR` | `${CMAKE_SOURCE_DIR}/ncbi-cxx-toolkit` | Path to NCBI C++ Toolkit install root |
| `NCBI_TOOLKIT_BUILD_TAG` | `CMake-GCC1330-Release` | Toolkit build subdirectory name |
| `PARASAIL_DIR` | `${CMAKE_SOURCE_DIR}/parasail` | Path to Parasail install root |
| `HTSLIB_DIR` | `${CMAKE_SOURCE_DIR}/htslib` | Path to htslib install root |
| `DROGON_DIR` | `${CMAKE_SOURCE_DIR}/drogon` | Path to Drogon install root (used when `BUILD_HTTPD=ON`) |
| `BUILD_HTTPD` | ON | Build ikafssnhttpd (requires Drogon) |
| `BUILD_CLIENT` | ON | Build ikafssnclient (requires libcurl for HTTP mode) |
| `ENABLE_REMOTE_RETRIEVE` | ON | Enable NCBI efetch in ikafssnretrieve |
