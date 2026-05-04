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
  BGZF and we deliberately do not double-wrap.  Use `-outfmt tab` /
  `-outfmt json` if you want a compressed result file.

### ikafssnindex

Build a k-mer inverted index from a BLAST database. For each volume, index files are generated: `.kix` (ID postings), `.kpx` (position postings, unless `-mode 1`), and `.ksx` (sequence metadata). When `-max_freq_build` is used, a shared `.khx` file (build-time exclusion bitset) is also generated. The `.khx` file is shared across all volumes (one per k value, not per volume).

```
ikafssnindex [options]

Required:
  -db <path>              BLAST DB prefix
  -k <int>                K-mer length (5-15)
  -o <dir>                Output directory

Options:
  -mode <1|2|3>           Search mode the index will support (default: 2)
                          1 = Stage 1 only (skip .kpx generation, saves disk and time)
                          2 = Stage 1+2 (default)
                          3 = Stage 1+2+3 (same as 2 for index build)
  -memory_limit <size>    Memory limit (default: half of physical RAM)
                          Accepts K, M, G suffixes
                          Partitions are auto-calculated to fit within this limit
  -max_freq_build <num>   Exclude k-mers with cross-volume count above this threshold
                          1 or 1.0: disable (no exclusion, default)
                          0 < x < 1: fraction of total NSEQ across all volumes
                          > 1: absolute count threshold
                          0: not allowed (error)
                          Counts are aggregated across all volumes before filtering
  -freq_threshold_part <int>
                          .kpx v7 per-(kmer, seq_id) partition threshold (default: 8, max: 255)
                          A (k-mer, seq_id) cluster with occurrence count > threshold
                          gets its own partition group; lower-multiplicity clusters
                          merge into a shared short bucket (improves chromosome-class
                          subject compression by decoupling per-block bit-width from
                          absolute position magnitude)
  -highfreq_filter_threads <int>
                          Threads for cross-volume high-frequency filtering
                          (default: min(8, threads))
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
  -threads <int>          Number of threads (default: all cores)
                          Parallelizes counting, partition scan, sort,
                          and volume processing
  -v, --verbose           Verbose output
```

**Examples:**

```bash
# Small DB, plenty of memory
ikafssnindex -db mydb -k 11 -o ./index -memory_limit 128G

# Large DB, limited memory, multi-threaded
ikafssnindex -db nt -k 11 -o ./nt_index -memory_limit 32G -threads 32

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
  -threads <int>          Parallel search threads (default: all cores)
  -memory_limit <size>    madvise WILLNEED budget (default: half of RAM)
                          Accepts K, M, G suffixes
  -mode <1|2|3>           Search mode (default: 2)
                          1=Stage 1 only, 2=Stage 1+2, 3=Stage 1+2+3
  -db <path>              BLAST DB path for mode 3 (default: same as -ix)
  -stage1_score <1|2>     Stage 1 score type (default: 1)
                          1=coverscore, 2=matchscore
  -stage1_topn <int>      Stage 1 candidate limit, 0=unlimited (default: 0)
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
  -stage2_min_diag_hits <int>  Diagonal filter min hits (default: 1)
  -context <value>        Context extension for mode 3 (default: 2.0)
                          Integer: bases to extend; Decimal: query length multiplier
  -stage3_traceback <0|1> Enable traceback in mode 3 (default: 0)
  -stage3_gapopen <int>   Gap open penalty for mode 3 (default: 10)
  -stage3_gapext <int>    Gap extension penalty for mode 3 (default: 1)
  -stage3_min_ppositive <num>  Min percent positive filter for mode 3 (default: 0)
  -stage3_min_npositive <int>  Min positive-scoring positions filter for mode 3 (default: 0)
  -stage3_score_matrix <name>  Score matrix: degmatch, dnafull, nuc44 (default: degmatch)
  -stage3_fetch_threads <int>  Threads for BLAST DB fetch in mode 3 (default: min(8, threads); error if exceeds -threads)
  -num_results <int>      Max results per query, 0=unlimited (default: 0)
  -seqidlist <path>       Include only listed accessions
  -negative_seqidlist <path>  Exclude listed accessions
  -strand <-1|1|2>       Strand to search (default: 2)
                          1=plus only, -1=minus only, 2=both
  -accept_qdegen <0|1>    Accept queries with degenerate bases (default: 1)
  -max_degen_expand <int> Max degenerate expansion per k-mer (default: 16, max: 256, 0/1: disable)
  -t <int>                Template length for spaced seeds (default: 0)
                          0: contiguous k-mers (traditional mode)
                          13, 15, 18: spaced seed template length (requires -k 8 or 9)
                          16, 18, 21: spaced seed template length (requires -k 11 or 12)
  -template_type <str>    Template type for spaced seeds (default: both)
                          coding: use coding index only
                          optimal: use optimal index only
                          both: merge coding and optimal indexes at search time
  -outfmt <tab|json|sam|bam>  Output format (default: tab)
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
/data/index/nt.00.11mer.kix
/data/index/nt.00.11mer.kpx
/data/index/nt.00.11mer.ksx
/data/index/nt.01.11mer.kix
/data/index/nt.01.11mer.kpx
/data/index/nt.01.11mer.ksx
/data/index/nt.11mer.kvx
```

then specify `-ix /data/index/nt`. The prefix `/data/index/nt` is split into the directory `/data/index/` and the base name `nt`. Volumes are discovered via the `.kvx` manifest file (`nt.11mer.kvx`), which lists the volume basenames. For aggregated databases (e.g. `combined` aggregating `foo` and `bar`), the index files would be `foo.11mer.kix`, `bar.11mer.kix`, with `combined.11mer.kvx` as the manifest.

If the index directory contains indexes for multiple k-mer sizes (e.g. both `nt.09mer.kvx` and `nt.11mer.kvx`), you must specify `-k` to select which one to use. If only a single k-mer size exists, `-k` can be omitted.

When `-accept_qdegen` is 0, queries containing IUPAC degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) are skipped with a `degen_rejected` skip-marker (TSV `*SKIPPED:degen_rejected`, JSON `"status": "skipped"`, SAM unmapped record with `XR:Z:degen_rejected`) and a stderr warning, and the exit code is 2. See "Skip reasons" below for the complete reason list and per-format representation. Set `-accept_qdegen 1` to allow such queries. K-mers containing exactly one degenerate base are expanded to all possible variants (e.g., R→A,G produces 2 k-mers; N→A,C,G,T produces 4) and used for search. K-mers whose per-position expansion product exceeds `-max_degen_expand` are skipped; when this occurs, a warning is emitted to stderr once per query indicating the query name and that such k-mers are ignored. (Note: this is a per-window unemit, not a whole-query skip — the rest of the query is still searched, and the unemit position is reflected in `Nhighfreq` for fractional thresholds.) In server mode (`ikafssnserver`), this warning is propagated through the protocol to `ikafssnclient`, which displays the same message. This matches the indexer's handling of subject-side degenerate bases.

`-seqidlist` and `-negative_seqidlist` are mutually exclusive. Both text (one accession per line) and binary (generated by `blastdb_aliastool -seqid_file_in`) formats are accepted, auto-detected by magic bytes.

**Examples:**

```bash
# Basic search
ikafssnsearch -ix ./index/mydb -query query.fasta -threads 8

# Specify k-mer size (required if index contains multiple k values)
ikafssnsearch -ix ./index/mydb -k 11 -query query.fasta

# Increase sensitivity
ikafssnsearch -ix ./index/mydb -query query.fasta \
    -stage2_min_score 2 -stage1_topn 2000

# Restrict to specific accessions
ikafssnsearch -ix ./index/mydb -query query.fasta -seqidlist targets.txt

# Exclude specific accessions
ikafssnsearch -ix ./index/mydb -query query.fasta -negative_seqidlist exclude.txt

# Fractional Stage 1 threshold (50% of query k-mers)
ikafssnsearch -ix ./index/mydb -query query.fasta -stage1_min_score 0.5

# Mode 3: pairwise alignment with traceback
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -num_results 5

# Mode 3: SAM output
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -outfmt sam -o result.sam

# Mode 3: BAM output
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -outfmt bam -o result.bam

# Mode 3: filter by percent positive
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -stage3_traceback 1 -stage3_min_ppositive 90

# Mode 3: with context extension (50 bases each side)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context 50 -num_results 5

# Mode 3: with context extension (0.1x query length each side)
ikafssnsearch -ix ./index/mydb -query query.fasta -mode 3 -context 0.1 -num_results 5

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
    -mode 3 -stage3_traceback 1 -num_results 10
```

### ikafssnretrieve

Extract matched subsequences based on search results. Supports local BLAST DB extraction and remote retrieval via NCBI E-utilities (efetch).

```
ikafssnretrieve [options]

Sequence source (one required):
  -db <path>              Local BLAST DB prefix
  -remote                 Retrieve from NCBI efetch

Input:
  -results <path>         Search results file (tab format)
  (none)                  Read from stdin

Common options:
  -o <path>               Output FASTA file (default: stdout)
                          Suffix (.gz/.bz2/.xz/.zst) selects an output codec
  -context <value>        Context extension (default: 2.0)
                          Integer: bases to add before/after match region
                          Decimal: multiplier of query length (qlen)
  -compression_level <int> Output compression level (defaults: gzip=6, bzip2=9, xz=6, zstd=3)
  -v, --verbose           Verbose logging

(Note: -results auto-detects gzip/bzip2/xz/zstd-compressed result files from
 their leading magic bytes; no flag is required.)

Remote options (-remote):
  -api_key <key>          NCBI API key (or NCBI_API_KEY env var)
  -batch_size <int>       Accessions per batch (default: 100)
  -retries <int>          Max retries (default: 3)
  -timeout <int>          Request timeout in seconds (default: 30)
  -range_threshold <int>  Seq length threshold for individual fetch (default: 100000)
```

**Examples:**

```bash
# Local BLAST DB extraction (file input)
ikafssnsearch -ix ./index/mydb -query query.fasta -o results.tsv
ikafssnretrieve -db nt -results results.tsv -o matches.fasta

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
ikafssnretrieve -db nt -results results.tsv -context 50

# Context as fraction of query length (0.1x each side)
ikafssnretrieve -db nt -results results.tsv -context 0.1
```

### ikafssnserver

Search daemon. Keeps index files memory-mapped and accepts search requests via UNIX domain socket or TCP socket.

```
ikafssnserver [options]

Required:
  -ix <prefix>            Index prefix (repeatable for multi-DB)

Listener (at least one required):
  -socket <path>          UNIX domain socket path
  -tcp <host>:<port>      TCP listen address

Options:
  -threads <int>          Worker threads (default: all cores)
  -max_queue_size <int>   Max concurrent query sequences globally (default: 1024)
  -max_seqs_per_req <int> Max sequences accepted per request (default: thread count)
  -pid <path>             PID file path
  -db <path>              BLAST DB path for mode 3 (repeatable, paired with -ix;
                          default: same as corresponding -ix prefix)
  -stage1_topn <int>      Default Stage 1 candidate limit (default: 0)
  -stage1_min_score <num> Default Stage 1 minimum score (default: 0.5)
                          Integer (>= 1) or fraction (0 < P < 1)
  -stage2_min_score <int> Default minimum chain score (default: 0 = adaptive)
  -stage2_max_gap <int>   Default chaining gap tolerance (default: 100)
  -stage2_max_lookback <int>  Default chaining DP lookback window (default: 64, 0=unlimited)
  -stage2_max_nhit_per_subject <int>  Default max chains per subject (default: 1, 0=unlimited)
  -stage2_min_diag_hits <int> Default diagonal filter min hits (default: 1)
  -context <value>        Default context extension (default: 2.0)
                          Integer: bases to extend; Decimal: query length multiplier
  -stage3_traceback <0|1> Default traceback mode (default: 0)
  -stage3_gapopen <int>   Default gap open penalty (default: 10)
  -stage3_gapext <int>    Default gap extension penalty (default: 1)
  -stage3_min_ppositive <num>  Default min percent positive (default: 0)
  -stage3_min_npositive <int>  Default min positive-scoring positions (default: 0)
  -stage3_score_matrix <name>  Default score matrix: degmatch, dnafull, nuc44 (default: degmatch)
  -stage3_fetch_threads <int>  Threads for BLAST DB fetch (default: min(8, threads))
  -num_results <int>      Default max results per query (default: 0)
  -accept_qdegen <0|1>    Default accept queries with degenerate bases (default: 1)
  -max_degen_expand <int> Max degenerate expansion per k-mer (default: 16, max: 256, 0/1: disable)
  -memory_limit <size>    madvise WILLNEED budget (default: half of RAM)
                          Accepts K, M, G suffixes
  -shutdown_timeout <int> Graceful shutdown timeout in seconds (default: 30)
  -v, --verbose           Verbose logging
```

**Examples:**

```bash
# Listen on UNIX socket
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -threads 16

# Listen on TCP (remote access)
ikafssnserver -ix ./nt_index -tcp 0.0.0.0:9100 -threads 32

# Both simultaneously
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -tcp 0.0.0.0:9100

# Daemon operation (started from systemd)
ikafssnserver -ix ./nt_index -socket /var/run/ikafssn.sock -pid /var/run/ikafssn.pid

# Mode 3 support: specify BLAST DB and Stage 3 defaults
ikafssnserver -ix ./nt_index -db nt -socket /var/run/ikafssn.sock \
    -stage3_traceback 1 -context 50

# Multi-DB: serve two databases in one process
ikafssnserver -ix ./nt_index -db nt -ix ./rs_index -db refseq_genomic \
    -socket /var/run/ikafssn.sock -threads 32
```

**Operational characteristics:**

- One process can serve multiple BLAST DB indexes simultaneously. Specify `-ix` (and optionally `-db`) multiple times to load several databases. Each database is identified by its basename (the last path component of the `-ix` prefix) and clients must specify `-db <name>` when the server hosts more than one database.
- If `-db` is specified, the count must match the number of `-ix` flags (paired in order). Databases without a `-db` override default to the `-ix` prefix as the BLAST DB path. A database with no `-db` path supports modes 1-2 only (max_mode=2); providing `-db` enables mode 3 (max_mode=3).
- If the index prefix matches indexes for multiple k-mer sizes, all are loaded and clients can specify k per request.
- On SIGTERM/SIGINT, performs graceful shutdown: stops accepting new connections, waits for in-flight requests to complete (up to `-shutdown_timeout` seconds), then exits.
- **Per-sequence concurrency control:** The server limits concurrency at the per-sequence level, not per-connection. When a request arrives, the server attempts to acquire permits for each valid query sequence. If the global limit (`-max_queue_size`) is reached, excess sequences are returned to the client as "rejected" for retry. The `-max_seqs_per_req` option caps how many permits a single request can acquire, preventing one large request from monopolizing all slots.

### ikafssnhttpd

HTTP REST API daemon. Connects to one or more `ikafssnserver` instances and exposes search as an HTTP API. Uses the Drogon framework. Multiple backends can be specified for multi-database support or load balancing of same-database replicas.

On startup, it connects to all configured backends to cache their capabilities (retrying with exponential backoff for up to 30 seconds). If the same database name appears on multiple backends, cross-server validation ensures that for each shared (db, k) pair, total sequence counts and total bases are identical; mismatches cause a startup error. Backends are allowed to have different k-value sets for the same database (e.g., server A has k=10, server B has k=10 and k=11); the merged capabilities expose the union of all k-value groups. Search requests are validated against the merged capabilities (synchronous, no backend round-trip) to reject obviously invalid requests immediately, then routed to the best available backend based on priority and slot availability.

**Routing and health:**

- **Priority**: Backends are prioritized by CLI argument order (first = highest priority).
- **Selection**: For each search request, the backend with the highest priority and available effective capacity is selected. Effective capacity considers both slot availability (`max_queue_size - queue_depth`) and per-request cap (`max_seqs_per_req`), taking the minimum of the two. If all backends are full, the highest-priority one is used.
- **Pre-check**: Before each search, a fresh info request is sent to the selected backend to verify connectivity.
- **Exclusion**: If a backend fails to respond (connection error on info or search), it is excluded for `-exclusion_time` seconds. Excluded backends are automatically re-checked during heartbeat and re-enabled once reachable.
- **Heartbeat**: A background thread refreshes all backends' info every `-heartbeat_interval` seconds.
- **No retry in httpd**: If a search request fails after backend selection, the error is returned to the client. `ikafssnclient` handles retry of rejected queries.

```
ikafssnhttpd [options]

Backend connection (at least one required; order = priority):
  -server_socket <path>      UNIX socket path to ikafssnserver
  -server_tcp <host>:<port>  TCP address of ikafssnserver

Options:
  -listen <host>:<port>       HTTP listen address (default: 0.0.0.0:8080)
  -path_prefix <prefix>       API path prefix (e.g., /nt)
  -threads <int>              Drogon I/O threads (default: all cores)
  -heartbeat_interval <int>   Heartbeat interval in seconds (default: 3600)
  -exclusion_time <int>       Backend exclusion time in seconds (default: 3600)
  -pid <path>                 PID file path
  -v, --verbose               Verbose logging
```

**REST API endpoints:**

| Method | Path | Description |
|---|---|---|
| POST | `/api/v1/search` | Search request (query sequences in JSON body) |
| GET | `/api/v1/info` | Aggregated index information from all backends |
| GET | `/api/v1/health` | Health check (OK if any backend is reachable) |

The `/api/v1/info` response aggregates databases from all healthy backends. For databases served by multiple backends, capacity is reported per mode in a `modes` array within each kmer_group, showing the sum of `max_queue_size`, `queue_depth`, and `max_seqs_per_req` (computed as `sum(min(available_i, per_req_i))` across backends) across all serving backends. A top-level `max_seqs_per_req` field reports the minimum across all modes.

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

Client command. Connects to `ikafssnserver` via socket or `ikafssnhttpd` via HTTP. Output format is identical to `ikafssnsearch`. Before sending any queries, the client performs pre-flight validation by fetching server capabilities and checking that the requested database name, k-mer size, and mode are valid. Invalid parameters produce an error with available database listings before any query data is transmitted. The client uses the server's `max_seqs_per_req` and available slot count to automatically split queries into appropriately-sized batches, avoiding oversized requests that would be partially rejected. Within each batch, if the server still rejects some query sequences due to concurrency limits, the client automatically retries the rejected queries with exponential backoff (30s, 60s, 120s, 120s, ...) until all queries are processed.

**Checkpointing:** The client automatically saves intermediate results to a temporary directory during batch processing. If the process is interrupted (e.g., Ctrl+C, network failure), re-running the same command resumes from where it left off, skipping already-completed queries. The temporary directory is named `{output}.{input}.{ix_name}.{kk}.ikafssn.tmp/` and is automatically cleaned up after successful completion. A directory-based lock prevents concurrent runs with the same parameters. Resume validation checks the search parameters, input file SHA256, and integrity of each batch file.

```
ikafssnclient [options]

Connection (one required):
  -socket <path>           UNIX domain socket path
  -tcp <host>:<port>       TCP server address
  -http <url>              ikafssnhttpd URL (e.g., http://example.com:8080)

Required:
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
  -stage1_score <1|2>      Stage 1 score type (default: server default)
  -stage1_topn <int>       Stage 1 candidate limit (default: server default)
  -stage1_min_score <num>  Stage 1 minimum score (default: server default)
                           Integer (>= 1) or fraction (0 < P < 1)
  -stage2_min_score <int>  Minimum chain score (default: server default)
                           0 = explicitly request adaptive mode
  -stage2_max_gap <int>    Chaining gap tolerance (default: server default)
  -stage2_max_lookback <int>  Chaining DP lookback window (default: server default)
  -stage2_max_nhit_per_subject <int>  Max chains per subject (default: server default)
  -stage2_min_diag_hits <int> Diagonal filter min hits (default: server default)
  -context <value>         Context extension (default: 2.0)
                           Integer: bases to extend; Decimal: query length multiplier
  -stage3_traceback <0|1>  Enable traceback (default: server default)
  -stage3_gapopen <int>    Gap open penalty (default: server default)
  -stage3_gapext <int>     Gap extension penalty (default: server default)
  -stage3_min_ppositive <num> Min percent positive filter (default: server default)
  -stage3_min_npositive <int> Min positive-scoring positions filter (default: server default)
  -stage3_score_matrix <name> Score matrix: degmatch, dnafull, nuc44 (default: server default)
  -num_results <int>       Max results per query (default: server default)
  -seqidlist <path>        Include only listed accessions
  -negative_seqidlist <path>  Exclude listed accessions
  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: server default)
  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 1)
  -max_degen_expand <int>  Max degenerate expansion (default: server default, max: 256)
  -t <int>                 Template length for spaced seeds (default: server default)
                           0: contiguous k-mers; 13, 15, 18 (k=8-9); 16, 18, 21 (k=11-12)
  -template_type <str>     Template type for spaced seeds (default: server default)
                           coding: use coding index only
                           optimal: use optimal index only
                           both: merge coding and optimal indexes at search time
  -outfmt <tab|json|sam|bam>  Output format (default: tab)
  -compression_level <int> Output compression level (defaults: gzip=6, bzip2=9, xz=6, zstd=3)
                           Codec is selected by -o suffix (.gz/.bz2/.xz/.zst); SAM/BAM reject all four
  -v, --verbose            Verbose logging

(Note: -query and -primer auto-detect gzip/bzip2/xz/zstd-compressed FASTA
 inputs from their leading magic bytes; no flag is required.)

HTTP Authentication (HTTP mode only):
  --user <user:password>   Credentials (curl-style)
  --http-user <USER>       Username (wget-style)
  --http-password <PASS>   Password (used with --http-user)
  --netrc-file <path>      .netrc file for credentials
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
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta --user admin:secret

# HTTP with Basic Auth (wget-style)
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta --http-user=admin --http-password=secret

# HTTP with .netrc credentials
ikafssnclient -http http://search.example.com:8080 -ix nt -query query.fasta --netrc-file ~/.netrc

# Mode 3: pairwise alignment with traceback
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -mode 3 -stage3_traceback 1

# Mode 3: SAM output
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -query query.fasta -mode 3 -stage3_traceback 1 -outfmt sam -o result.sam

# In-Silico PCR (primer mode)
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -primer primers.fasta -insert_length 500

# In-Silico PCR with custom thresholds
ikafssnclient -tcp 10.0.1.5:9100 -ix nt -primer primers.fasta -insert_length 300 \
    -stage1_primer_score 0.8

# In-Silico PCR with mode 3 alignment
ikafssnclient -socket /var/run/ikafssn.sock -ix nt -primer primers.fasta -insert_length 500 \
    -mode 3 -stage3_traceback 1 -num_results 10
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

Local mode options:
  -db <path>               BLAST DB prefix (default: auto-detect from -ix)

Remote HTTP authentication:
  --user <user:password>   Credentials (curl-style)
  --http-user <USER>       Username (wget-style)
  --http-password <PASS>   Password (used with --http-user)
  --netrc-file <path>      .netrc file for credentials

Options:
  -v, --verbose            Verbose output
```

`-ix` and remote options (`-socket`, `-tcp`, `-http`) are mutually exclusive. Only one remote option may be specified at a time.

**Local mode** reads index files directly and displays detailed statistics. When `-db` is not specified, `ikafssninfo` attempts to auto-detect the BLAST DB by checking whether the index prefix path corresponds to a valid BLAST DB.

Local mode output includes:

- K-mer length (k) and integer type (uint16/uint32)
- Number of volumes
- Per-volume statistics: sequence count, total postings, file sizes, excluded k-mer count (if `.khx` present)
- Overall statistics: total sequences, total postings, total index size, compression ratio
- With `-v`: k-mer frequency distribution (min, max, mean, percentiles)
- With `-db` (or auto-detected): BLAST DB title, sequence count, total bases, volume paths

**Remote mode** queries a running server and displays its capabilities.

Remote mode output includes:

- Active/max sequence slots
- Per-database information: name, default k, max mode, k-mer groups with volume counts, sequence counts, total bases, and posting statistics
- With `-v`: per-volume details (sequence count, total bases, postings) within each k-mer group

**Examples:**

```bash
# Local: basic index info
ikafssninfo -ix ./index/mydb

# Local: include BLAST DB info
ikafssninfo -ix ./index/mydb -db mydb

# Local: detailed frequency distribution
ikafssninfo -ix ./index/mydb -v

# Remote: query server via UNIX socket
ikafssninfo -socket /var/run/ikafssn.sock

# Remote: query server via TCP
ikafssninfo -tcp 10.0.1.5:9100

# Remote: query server via HTTP
ikafssninfo -http http://search.example.com:8080

# Remote: verbose (show per-volume details)
ikafssninfo -socket /var/run/ikafssn.sock -v

# Remote: HTTP with authentication
ikafssninfo -http http://search.example.com:8080 --user admin:secret
```

## Search Pipeline

ikafssn uses a three-stage search pipeline:

The default parameters prioritize throughput: `stage1_topn=0` and `num_results=0` disable sorting, and `stage1_min_score=0.5` (fractional) filters candidates by requiring at least 50% of query k-mers to match. To get ranked output, set positive values for `-stage1_topn` and/or `-num_results`, which triggers sorting but may reduce speed for large result sets.

1. **Stage 1 (Candidate Selection):** Scans ID postings for each query k-mer and accumulates scores per sequence. Two score types are available: **coverscore** (number of distinct query k-mers matching the sequence) and **matchscore** (total k-mer position matches). Sequences exceeding `stage1_min_score` are selected as candidates. When `stage1_topn > 0`, candidates are sorted by score and truncated. When `stage1_topn = 0` (default), all qualifying candidates are returned without sorting.

2. **Stage 2 (Collinear Chaining):** For each candidate, collects position-level hits from the `.kpx` file, applies a diagonal filter, and runs a chaining DP to find the best collinear chain. The chain length is reported as **chainscore**. Chains with `chainscore >= stage2_min_score` are reported. The DP inner loop is limited by `-stage2_max_lookback` (default: 64), restricting each hit to consider only the preceding B hits as potential chain predecessors. This reduces worst-case complexity from O(n²) to O(n×B) when a single query×subject pair has a very large number of hits. Set to 0 for unlimited (original O(n²) behavior). When `-stage2_max_nhit_per_subject` is greater than 1 (or 0 for unlimited), multiple non-overlapping chains are extracted per subject using greedy best-chain removal: the best chain is found and its hits are removed, then the DP is re-run on the remaining hits, repeating until the limit is reached or no chain meets `min_score`.

3. **Stage 3 (Pairwise Alignment):** For each Stage 2 hit, retrieves the subject subsequence from the BLAST DB (with optional context extension via `-context`), and performs semi-global pairwise alignment using the Parasail library (using the score matrix specified by `-stage3_score_matrix`, default: DEGMATCH). The alignment score (**alnscore**) is computed for all hits. When `-stage3_traceback 1` is enabled, CIGAR strings, percent positive (ppositive), positive-scoring position count (npositive), negative-scoring count (nnegative), and aligned sequences (with gaps) are also computed. Hits can be filtered by `-stage3_min_ppositive` and `-stage3_min_npositive` (traceback mode only). Subject sequences are pre-fetched in parallel across BLAST DB volumes controlled by `-stage3_fetch_threads`.

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
| `kSkipQueryTooShort` / `query_too_short` | `seq_len < span` (k-mer / spaced-seed window cannot fit). |
| `kSkipDegenRejected` / `degen_rejected` | `-accept_qdegen 0` was set and the query contains IUPAC degenerate bases. |
| `kSkipInvalidChar` / `invalid_char` | Query contains a character outside `[ACGT]` ∪ IUPAC. The detail names the offending position. |
| `kSkipThresholdUnreachable` / `threshold_unreachable` | Resolved fractional threshold is `< 1` on every searched strand. |

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
| **coverscore** | Number of distinct query k-mers that match the subject sequence. Each query k-mer is counted at most once per subject, regardless of how many positions it matches. | Stage 1 |
| **matchscore** | Total number of (query k-mer, subject position) matches. A single query k-mer matching multiple positions in the subject counts multiple times. | Stage 1 |
| **chainscore** | Length (number of k-mer hits) of the best collinear chain found by the chaining DP. Requires position data from `.kpx`. | Stage 2 |
| **alnscore** | Semi-global pairwise alignment score computed by Parasail (using `-stage3_score_matrix`, default: degmatch). Requires BLAST DB for subject sequence retrieval. | Stage 3 |

- `-stage1_score` selects which score type Stage 1 uses (1=coverscore, 2=matchscore). This affects candidate ranking and the stage1 score reported in output.
- The sort key is auto-determined by mode: mode 1 sorts by stage1 score, mode 2 by chainscore, mode 3 by alnscore.
- In `-mode 1`, only Stage 1 scores are available; chainscore and alnscore are not computed.

### Stage 3 Scoring Matrix

Stage 3 pairwise alignment uses a configurable scoring matrix specified by `-stage3_score_matrix`. Three matrices are available:

| Matrix | Description |
|---|---|
| **degmatch** (default) | Assigns positive scores to degenerate base pairs that share at least one possible nucleotide identity. Suitable for primer matching and searches involving degenerate bases. |
| **dnafull** | Traditional EMBOSS DNA full matrix (created by Todd Lowe). Uses ambiguous nucleotide codes with probabilities rounded to nearest integer. Degenerate pairs that share partial overlap receive negative or low scores. |
| **nuc44** | NCBI BLAST nucleotide matrix. Similar to dnafull but with slightly different degenerate base scoring. |

**CIGAR `=`/`X` determination:** The extended CIGAR operators `=` (sequence match) and `X` (sequence mismatch) are determined based on alignment scores, not strict base identity. A position is reported as `=` when its score in the matrix is positive (score > 0), and as `X` when its score is zero or negative. This means that with the DEGMATCH matrix, a degenerate base pair like N-A (score 1) is counted as a match (`=`), while with NUC44/DNAFULL it would be counted as a mismatch (`X`).

**Column name changes:** The output columns previously named `pident` (percent identity), `nident` (number of identical positions), and `mismatch` (number of mismatches) have been renamed to `ppositive` (percent positive-scoring), `npositive` (number of positive-scoring positions), and `nnegative` (number of negative-scoring positions). This reflects the fact that with the DEGMATCH matrix, these counts represent positive-scoring positions rather than strictly identical positions. The corresponding filter options have been renamed from `-stage3_min_pident`/`-stage3_min_nident` to `-stage3_min_ppositive`/`-stage3_min_npositive`.

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

### Tab Format (default)

**Mode 2** (default):

Tab-separated columns, where `coverscore` is replaced by `matchscore` when `-stage1_score 2`:

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

Note: `qstart` and `sstart` are omitted because accurate alignment start positions are unavailable without traceback.

**Mode 3, traceback=1** (`-mode 3 -stage3_traceback 1`):

```
# qseqid  sseqid  sstrand  qstart  qend  qlen  sstart  send  slen  coverscore  chainscore  alnscore  ppositive  npositive  nnegative  cigar  qseq  sseq  volume
```

### JSON Format

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
          "qstart": 0,
          "qend": 150,
          "qlen": 200,
          "sstart": 1000,
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
          "qstart": 0,
          "qend": 150,
          "qlen": 200,
          "sstart": 1000,
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

SAM/BAM output requires `-mode 3 -stage3_traceback 1`. Use `-outfmt sam` for SAM or `-outfmt bam` for BAM (BAM requires `-o <path>`).

SAM records contain:
- **QNAME**: qseqid
- **FLAG**: 0 (forward) or 16 (reverse)
- **RNAME**: sseqid
- **POS**: sstart + 1 (1-based)
- **MAPQ**: 255
- **CIGAR**: extended CIGAR with =/X/I/D operators
- **SEQ**: ungapped query sequence
- **QUAL**: * (not available)
- **Tags**: `AS:i` (alnscore), `NM:i` (nnegative), `cs:i` (chainscore), `cv:i` (coverscore), `ms:i` (matchscore)

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

When the same database name appears on multiple backends, `ikafssnhttpd` verifies at startup that for each shared (db, k) pair, total sequence counts and total bases are identical. Backends may have different k-value sets for the same database; the merged capabilities expose the union of all k-value groups, so requests for a k-value available on only some backends are naturally routed to those backends. Requests are routed to the highest-priority backend with available effective capacity (considering both slot availability and `max_seqs_per_req`). Note that capacity values (`max_queue_size`, `queue_depth`, `max_seqs_per_req`) are shared per server across all databases served by that server.

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

## Index File Formats

Per BLAST DB volume, three files are generated using the BLAST DB volume's own basename:

```
<vol_basename>.<kk>mer.kix   — ID postings (direct-address table + delta-compressed)
<vol_basename>.<kk>mer.kpx   — Position postings (delta-compressed)
<vol_basename>.<kk>mer.ksx   — Sequence metadata (lengths + accessions)
```

A `.kvx` manifest file is always generated for volume discovery:

```
<db_base>.<kk>mer.kvx        — Volume manifest (text, lists volume basenames)
```

When `-max_freq_build` is used, a shared exclusion bitset file is also generated (one per k value, shared across all volumes):

```
<db_base>.<kk>mer.khx        — Build-time exclusion bitset (shared across volumes)
```

Examples:
- Standard multi-volume (`nt` with volumes `nt.00`, `nt.01`): `nt.00.11mer.kix`, `nt.01.11mer.kpx`, `nt.11mer.kvx`, `nt.11mer.khx`
- Aggregated (`combined` with volumes `foo`, `bar`): `foo.11mer.kix`, `bar.11mer.kix`, `combined.11mer.kvx`

The `.khx` file contains a 32-byte header (magic "KMHX", format version, k) followed by a bitset of `ceil(4^k / 8)` bytes. Bit *i* = 1 indicates that k-mer *i* was excluded during index build based on cross-volume aggregated counts.

ID and position postings are stored in separate files so that Stage 1 filtering never touches `.kpx`, maximizing page cache efficiency.

### Spaced Seed Index File Naming

When spaced seeds are enabled (`-t > 0`), the file naming includes the template length and type:

```
<vol_basename>.<kk>mer.<tt>mer.<type>.kix
<vol_basename>.<kk>mer.<tt>mer.<type>.kpx
<vol_basename>.<kk>mer.<tt>mer.<type>.ksx
<db_base>.<kk>mer.<tt>mer.<type>.kvx
<db_base>.<kk>mer.<tt>mer.<type>.khx
```

Where `<tt>` is the zero-padded template length and `<type>` is `cod` (coding) or `opt` (optimal).

Examples:
- Spaced seed with k=11, t=16, coding: `nt.00.11mer.16mer.cod.kix`
- Spaced seed with k=11, t=16, optimal: `nt.00.11mer.16mer.opt.kix`
- Spaced seed with k=12, t=21, coding: `nt.00.12mer.21mer.cod.kix`
- Manifest: `nt.11mer.16mer.cod.kvx`

**Multi-accession deflines:** When the source BLAST DB was built with `makeblastdb -parse_seqids` and carries multi-defline records (the NCBI convention for registering identical sequences under several accessions, separated by `\x01` / `^A` in the FASTA defline), `ikafssnindex` preserves **all** accessions for each OID. The `.ksx` accession string for such OIDs contains every accession joined by `\x01`, and search output emits the same `\x01`-joined string in the `sseqid` column / SAM RNAME / FASTA defline / protocol `sseqid` field. Downstream consumers should split on `\x01` to recover individual accessions. The `-seqidlist` filter and `ikafssnretrieve` accept either the full `\x01`-joined form or any individual constituent accession.

**Index format version:** The current index format is version 9 for all index files (`.kix`, `.kpx`, `.ksx`, `.khx`). Key changes from earlier versions:

- **`.kix` v9 (Phase 7):** The dictionary is stored as an Elias-Fano blob in place of the raw `u32`/`u64` `offsets[]` array (4.6× smaller on average across NCBI nt_v4). Each posting list header is now just `[u32 distinct_count]` (4 B) — `body_words` was removed (Phase 7c dedup B; derived from the EF dictionary's `posting_byte_length`). Body encoding is unchanged: FastPFor's `CompositeCodec<SIMDFastPFor<4>, VariableByte>` (PForDelta with VByte exception stream) over the **distinct seq_id** delta stream `[abs_first, d1, d2, ...]` with `d_i >= 1`. The legacy `KIX_FLAG_OFFSET32` flag (0x04) is reserved (writers force-clear it; readers ignore it).
- **`.kpx` v9 (Phase 7):** The `pos_offsets` dictionary is also Elias-Fano. Per-posting-list, all four redundant header `u32` fields were dropped — the body starts directly at the 2-bit kind map. `distinct_count` is taken from the `kix_count` decoder parameter (Phase 7c dedup A); `partition_count`, `short1_count`, `short2_count` are derived from a SIMD popcount of the kind map (Phase 7d dedup C); `short2_position_count` is derived as the cumulative sum that builds the short2 offset table from the `u8 occ_count[]` array (Phase 7d dedup D). Empty `.kpx` posting lists emit zero bytes. The `KpxHeader.offset_type` byte is reserved (writers set the EF sentinel `0xFF`; readers ignore it).
- **`.ksx` / `.khx` v9 (Phase 7):** Data layouts are unchanged from v8; the `format_version` field bumps for family-wide alignment (single-major-version policy). Magic strings stay `KMSX` / `KMHX`.
- **`.kvx` v9 (Phase 7):** Manifest text format is unchanged; the `FORMAT_VERSION` line bumps to 9.

### Migration (v8 → v9)

ikafssn 0.1.2026.05.03+ requires v9 indexes. v8 indexes are rejected at open with the message:

```
KixReader: index format version mismatch (got 8, expected 9). Please rebuild with the current ikafssnindex.
```

(and analogous messages for `.kpx` / `.ksx` / `.khx`). To migrate, rebuild the index from the BLAST DB:

```
ikafssnindex -db <BLAST_DB_prefix> -k <k> -o <out_dir>
```

For NCBI nt-scale databases (~700 volumes at k=12 t=21), expect tens of hours of rebuild time on a 32-core host. v9 indexes are 25–35 % smaller on disk than v8 (the dictionary EF blob is ~4.6× smaller, posting list headers shrink by 4–24 B per non-empty k-mer); RAM/page-cache savings are substantially larger because the dictionary and posting list headers sit on the Stage 1 hot path.

Indexes built before v8 are still rejected with the same "rebuild your index" message. Rebuild after upgrading.

## Installation

### Ubuntu (.deb package)

Pre-built `.deb` packages are available for Ubuntu 22.04 and 24.04 (amd64 and arm64). Download the appropriate package from the [GitHub Releases](https://github.com/astanabe/ikafssn/releases) page.

Package naming convention:

```
ikafssn_<version>_ubuntu-<ubuntu_ver>_<arch>.deb
```

Install the package:

```bash
sudo apt install ./ikafssn_<version>_ubuntu-<ubuntu_ver>_<arch>.deb
```

### Enterprise Linux (.rpm package)

Pre-built `.rpm` packages are available for AlmaLinux / RHEL / Rocky Linux 9 and 10 (x86_64 and aarch64) from the [GitHub Releases](https://github.com/astanabe/ikafssn/releases) page.

Package naming convention:

```
ikafssn-<version>.el<el_ver>.<arch>.rpm
```

where `<el_ver>` is `9` or `10` and `<arch>` is `x86_64` or `aarch64`.

Install the package:

```bash
sudo dnf install ./ikafssn-<version>.el<el_ver>.<arch>.rpm
```

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
- Drogon (for ikafssnhttpd, optional)
- libcurl (for HTTP client mode and remote retrieval, optional)

### CPU requirements

- **x86_64**: SSE4.2 minimum (Intel Nehalem 2008+, AMD Bulldozer 2011+). The runtime SIMD dispatcher additionally targets AVX2, AVX-512 BW, and AVX-512 VBMI2 when present. CPUs without SSE4.2 are rejected at startup with `exit(2)`.
- **aarch64**: NEON (ASIMD) minimum (Armv8.0+). Pre-NEON aarch64 CPUs are rejected at startup with `exit(2)`. SVE / SVE2 capable CPUs use the NEON-tier FastPFor codec object (the dispatcher routes any tier ≥ NEON to a single SIMDe-translated NEON OBJECT library); per-kernel SIMD files (toupper, ncbi2na unpack, k-mer revcomp, degenerate scan, spaced-seed) keep their separate SVE / SVE2 dispatches.

### Installing Dependencies

Install the required packages (excluding NCBI C++ Toolkit) with the following commands.

**Ubuntu Server 22.04 / 24.04:**

```bash
sudo apt install build-essential cmake libtbb-dev liblmdb-dev libsqlite3-dev \
    libcurl4-openssl-dev libjsoncpp-dev
sudo apt install zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev autoconf \
    libssl-dev uuid-dev
```

The second line installs dependencies required for building Parasail and htslib from source, plus `libssl-dev` and `uuid-dev` which are needed for building the NCBI C++ Toolkit and Drogon from source. If ikafssnhttpd is not needed, build with `-DBUILD_HTTPD=OFF` and `uuid-dev` may be omitted (`libssl-dev` is still required by the NCBI C++ Toolkit).

**AlmaLinux 9 / Rocky Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y epel-release
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y oracle-epel-release-el9
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**AlmaLinux 10 / Rocky Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

**Oracle Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y zlib-devel bzip2-devel xz-devel libdeflate-devel autoconf
sudo dnf install -y libuuid-devel openssl-devel
```

On EL9, `jsoncpp-devel` requires EPEL and `lmdb-devel` requires CRB. On EL10, both are in CRB so EPEL is not needed for these packages. The second-to-last line of each block installs dependencies required for building Parasail and htslib from source. The last line installs dependencies needed to build the NCBI C++ Toolkit and Drogon from source. If ikafssnhttpd is not needed, build with `-DBUILD_HTTPD=OFF` and `libuuid-devel` may be omitted (`openssl-devel` is still required by the NCBI C++ Toolkit).

**macOS (Homebrew):**

```bash
brew install cmake tbb lmdb sqlite3 curl jsoncpp \
    xz libdeflate autoconf automake libtool openssl@3
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
