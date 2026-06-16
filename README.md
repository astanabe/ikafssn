# ikafssn

**ikafssn** (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) is a suite of command-line tools that builds inverted indexes over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.

A conference talk on the methods used in ikafssn is available on YouTube: <https://youtu.be/ppoTB6MHsqY> (spoken in Japanese, slides in English).

## Features

- Builds a k-mer inverted index directly from NCBI BLAST databases
- Supports both contiguous k-mers and discontiguous megablast-style spaced seeds (`-t 13/15/18` with k=8 or 9 for PCR-optimized search, `-t 16/18/21` with k=11 or 12, coding/optimal/both templates)
- Three-stage search pipeline: fast candidate filtering (Stage 1), position-aware collinear chaining (Stage 2), and Parasail pairwise alignment with CIGAR/percent identity output (Stage 3), with configurable mode selection (1/2/3)
- Client-server architecture with UNIX/TCP socket and HTTP REST API support, with multi-database serving from a single process
- Handles IUPAC ambiguous bases by expanding degenerate k-mers during indexing and search (configurable expansion limit)
- Reports queries that cannot be searched as skip-markers in TSV / JSON / SAM (with reason and detail), so no query is silently dropped
- Transparent compressed I/O: FASTA inputs auto-detect gzip / bzip2 / xz / zstd from leading magic bytes; TSV / JSON / FASTA outputs select the codec by `.gz` / `.bz2` / `.xz` / `.zst` suffix (SAM/BAM excluded)
- Parallel indexing and search via Intel TBB
- Lightweight per-command executables, each linking only its required dependencies

## Comparison with BLAST (blastn)

ikafssn and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (blastn) both search NCBI BLAST databases, but differ in indexing strategy, search algorithm, and deployment model:

| Aspect | ikafssn | blastn |
|---|---|---|
| **Primary use case** | K-mer-based similarity search across large nucleotide collections | Local alignment search for homologous sequences |
| **Search algorithm** | Seed-chain-align: k-mer inverted index → candidate filtering → collinear chaining → pairwise alignment | Seed-extend: word seeds → ungapped extension → gapped alignment (Smith-Waterman) |
| **Database format** | NCBI BLAST DB (reads via C++ Toolkit) | NCBI BLAST DB (native) |
| **Index structure** | Pre-built k-mer inverted index (direct-address table, 4^k entries) stored on disk | No pre-built database index; per-query lookup table built from query words, database scanned sequentially |
| **Pipeline modes & index content** | Default Mode 1: Stage 1 candidate filtering only, over a presence-only index (`.kix`, sequence IDs without positions); long parents are split into overlapping fragments at index-build time to preserve locality without storing positions. Mode 2 (optional) adds collinear chaining and Mode 3 (optional) adds pairwise alignment, both over a position-bearing index (`.kpx`) | No pre-built database index and no staged search modes; each query runs the full seed → ungapped extension → gapped alignment pipeline |
| **Seeding** | All k-mers indexed exhaustively (contiguous or spaced seeds); high-frequency filtering at index build time only (optional `-max_freq_build`) | Exact word matches (default word size 28 for megablast, 11 for blastn); discontiguous megablast templates |
| **Scoring model** | K-mer match count (Stage 1), chain length (Stage 2), semi-global alignment score (Stage 3) | E-value based on local alignment score (bit score) with statistical significance model |
| **Alignment** | Parasail semi-global, 1-piece affine gap (optional Stage 3) | BLAST gapped extension, affine gap penalties, with X-drop heuristic |
| **Hits per subject** | Configurable via `-stage2_max_nhit_per_subject` (default: 1; 0 = unlimited) | Multiple HSPs per subject by default |
| **Ambiguous bases** | IUPAC degenerate expansion in both index and query (configurable limit) | Not expanded; seeds at ambiguous positions are skipped, extension penalizes via scoring matrix |
| **Client-server mode** | Built-in: `ikafssnserver` + `ikafssnclient` (UNIX/TCP socket), `ikafssnhttpd` (HTTP REST) | Not built-in (requires external wrappers or NCBI cloud BLAST) |
| **Parallelization** | Intel TBB (`parallel_for`, `parallel_sort`); multi-query and multi-volume parallelism | OpenMP-based multi-threading (`-num_threads`) |
| **Output format** | TSV, JSON, SAM, BAM | Tabular, XML, ASN.1, HTML, SAM, and others (`-outfmt`) |
| **Expected query sequences** | Short to moderate-length sequences (PCR amplicons, marker genes) | Any length (short queries to full chromosomes) |

## Comparison with minimap2

ikafssn and [minimap2](https://github.com/lh3/minimap2) both follow the seed-chain-align paradigm, but differ in design goals and key algorithmic choices:

| Aspect | ikafssn | minimap2 |
|---|---|---|
| **Primary use case** | Database search: find all similar entries across a large sequence collection | Read mapping: map reads to a single reference genome |
| **Database input format** | NCBI BLAST DB | FASTA, FASTQ |
| **Query input format** | FASTA | FASTA, FASTQ |
| **Expected query sequences** | Short to moderate-length sequences such as PCR amplicons and marker genes (hundreds to a few thousand bases) | Genomic reads (Illumina short reads, PacBio/ONT long reads) and assembled sequences (contigs, chromosomes) |
| **Expected query quality** | High-accuracy sequences with few or no sequencing errors (e.g., Sanger, error-corrected consensus) | Designed to tolerate high error rates (ONT ~5–15%, PacBio CLR ~10–15%; also handles HiFi and Illumina) |
| **K-mer size** | k = 5–15 (`uint16_t` for k ≤ 8, `uint32_t` for k ≥ 9; selected by `2*k` bit width) | k = 1–28 (default 15) |
| **Seeding** | All k-mers indexed in a direct-address table (4^k entries) | Minimizers (subsampled k-mers) indexed in a hash table |
| **Pipeline modes & index content** | Default Mode 1: Stage 1 candidate filtering only, over a presence-only index (`.kix`, sequence IDs without positions); long parents are split into overlapping fragments at index-build time to preserve locality without storing positions. Mode 2 (optional) adds collinear chaining and Mode 3 (optional) adds pairwise alignment, both over a position-bearing index (`.kpx`) | Single seed-chain-align pipeline with no presence-only mode; the minimizer index always stores positions (reference ID, position, strand) |
| **Candidate filtering** | Explicit Stage 1: scan ID posting lists to score and filter candidates before chaining | No separate filtering stage; seeds go directly to chaining |
| **Chaining DP score** | Chain length (anchor count) with a diagonal-deviation constraint (`max_gap`) | Estimated matching bases minus a gap penalty with logarithmic distance term |
| **Hits per subject** | Configurable via `-stage2_max_nhit_per_subject` (default: 1 best chain; set >1 or 0 for unlimited) | Multiple chains per subject (primary, secondary, supplementary) |
| **Default result limit** | Unlimited (returns all candidates above threshold) | 1 primary + up to 5 secondary per query (`-N 5`, `-p 0.8`) |
| **Alignment** | Parasail semi-global, 1-piece affine gap (optional Stage 3) | KSW2 with SIMD, 2-piece affine gap and Z-drop heuristic |
| **Ambiguous bases** | IUPAC degenerate expansion (configurable limit) | Not supported (N-containing k-mers skipped) |
| **Output format** | TSV, JSON, SAM, BAM | PAF, SAM |
| **Index partitioning** | BLAST DB multi-volume with `.kvx` manifest | `-I` batch partitioning with `--split-prefix` |

## Comparison with MMseqs2

ikafssn and [MMseqs2](https://github.com/soedinglab/MMseqs2) (nucleotide-vs-nucleotide search, `--search-type 3`) both prune candidates with a k-mer prefilter before alignment, but differ in index design, default pipeline, and intended workload. MMseqs2 is primarily a protein search and clustering suite; nucleotide-vs-nucleotide search is a secondary mode that is most effective for highly similar sequences:

| Aspect | ikafssn | MMseqs2 (`--search-type 3`) |
|---|---|---|
| **Primary use case** | K-mer-based similarity search across large nucleotide collections | Many-against-many protein search and clustering; nucleotide-vs-nucleotide as a secondary mode |
| **Search algorithm** | Seed-chain-align: k-mer inverted index → candidate filtering → collinear chaining → pairwise alignment | Seed-prefilter-align: k-mer match → ungapped diagonal scoring → Smith-Waterman |
| **Database input format** | NCBI BLAST DB (reads via C++ Toolkit) | Own MMseqs2 DB built by `createdb` from FASTA/FASTQ (BLAST DB not read directly) |
| **Index structure** | Pre-built on-disk k-mer inverted index (direct-address table, 4^k entries), always accessed via mmap | Optional on-disk index built by `createindex` for repeated searches (loaded into memory or mmapped, `--db-load-mode`); without it, the k-mer index is built in memory at search time — but unlike blastn, a target-side k-mer index is always used |
| **Pipeline modes & index content** | Default Mode 1: Stage 1 candidate filtering only, over a presence-only index (`.kix`, sequence IDs without positions); long parents are split into overlapping fragments at index-build time to preserve locality without storing positions. Mode 2 (optional) adds collinear chaining and Mode 3 (optional) adds pairwise alignment, both over a position-bearing index (`.kpx`) | Single prefilter → alignment pipeline with no presence-only mode; the k-mer index always stores positions (`seqId` + `position_j`), required for diagonal scoring |
| **Seeding** | All k-mers indexed exhaustively (contiguous or spaced seeds); high-frequency filtering at index build time only (optional `-max_freq_build`) | k-mers (default k=15 for nucleotides), optional spaced seeds (`--spaced-kmer-mode`); per query position, similar k-mers are derived from a substitution matrix and looked up together, controlled by sensitivity `-s` |
| **Candidate filtering** | Explicit Stage 1: scan ID posting lists to score and filter candidates; an optional per-diagonal hit-count filter and collinear chaining handle positional consistency in Stage 2 | Prefilter requiring two consecutive k-mer matches on the same diagonal, scored by an ungapped diagonal alignment threshold |
| **Alignment** | Parasail semi-global, 1-piece affine gap (Mode 3, optional) | Vectorized (striped SIMD) Smith-Waterman local alignment |
| **Scoring model** | K-mer match count (Stage 1), chain length (Stage 2), semi-global alignment score (Stage 3) | Prefilter diagonal score; alignment score, sequence identity, and E-value |
| **Ambiguous bases** | IUPAC degenerate expansion in both index and query (configurable limit) | No IUPAC degenerate expansion |
| **Client-server mode** | Built-in: `ikafssnserver` + `ikafssnclient` (UNIX/TCP socket), `ikafssnhttpd` (HTTP REST) | Not built-in (a separate desktop / local web-server app exists) |
| **Parallelization** | Intel TBB (`parallel_for`, `parallel_sort`); multi-query and multi-volume parallelism | OpenMP multi-threading + MPI distribution; GPU acceleration (NVIDIA Ampere or newer) |
| **Output format** | TSV, JSON, SAM, BAM | BLAST m8 tabular (via `convertalis`), and other MMseqs2 formats |
| **Clustering** | Not provided (search only) | Core feature (`cluster` / `linclust`), mature mainly for proteins; limited for nucleotides |
| **Expected query sequences** | Short to moderate-length sequences (PCR amplicons, marker genes) | Best for highly similar nucleotide sequences; less sensitive for divergent ones |

## Commands

| Command | Purpose |
|---|---|
| `ikafssnindex` | Build k-mer inverted index from BLAST DB |
| `ikafssnsearch` | Local direct search (mmap index) |
| `ikafssnretrieve` | Extract matched subsequences (local or NCBI efetch) |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver |
| `ikafssnclient` | Client (socket or HTTP) |
| `ikafssninfo` | Index/DB and server information display |

## Quick Start

```bash
# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Index a BLAST DB
ikafssnindex -db mydb -k 11 -o ./index

# Search
ikafssnsearch -ix ./index -query query.fasta

# Retrieve matched subsequences
ikafssnsearch -ix ./index -query query.fasta | ikafssnretrieve -db mydb > matches.fasta
```

## Dependencies

- C++20 compiler (GCC >= 11, Clang >= 13)
- CMake >= 3.16
- NCBI C++ Toolkit (for BLAST DB access)
- Intel TBB (for parallelization)
- Parasail >= 2.6 (for Stage 3 pairwise alignment)
- htslib >= 1.17 (for SAM/BAM output)
- zlib, libbz2, liblzma, libzstd (for transparent compressed I/O); pkg-config required for libzstd discovery
- Drogon (optional, for ikafssnhttpd)
- libcurl (optional, for HTTP client mode and remote retrieval)

## Documentation

For detailed usage, command-line options, deployment architecture, and index file format specifications, see:

- **English**: [doc/ikafssn.en.md](doc/ikafssn.en.md)
- **Japanese**: [doc/ikafssn.ja.md](doc/ikafssn.ja.md)

## Repository

- **Primary**: <https://github.com/astanabe/ikafssn>
- **Secondary**: <https://gitlab.com/astanabe/ikafssn>

## License

This project is licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for details.
