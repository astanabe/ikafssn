# ikafssn

**ikafssn** (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) is a suite of command-line tools that builds inverted indexes over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.

A conference talk on the methods used in ikafssn is available on YouTube: <https://youtu.be/ppoTB6MHsqY> (spoken in Japanese, slides in English).

## Features

- Builds a k-mer inverted index directly from NCBI BLAST databases
- Supports both contiguous k-mers and discontiguous megablast-style spaced seeds (`-t 13/15/18` with k=8 or 9 for PCR-optimized search, `-t 16/18/21` with k=11 or 12, coding/optimal/both templates)
- Three-stage search pipeline: fast candidate filtering (Stage 1), position-aware collinear chaining (Stage 2), and Parasail pairwise alignment with CIGAR/percent identity output (Stage 3), with configurable mode selection (1/2/3)
- Client-server architecture with UNIX/TCP socket and HTTP REST API support, with multi-database serving from a single process
- Handles IUPAC ambiguous bases by expanding degenerate k-mers during indexing and search (configurable expansion limit)
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
| **Seeding** | All k-mers indexed exhaustively (contiguous or spaced seeds); high-frequency filtering at build time (optional `-max_freq_build`) and/or search time | Exact word matches (default word size 28 for megablast, 11 for blastn); discontiguous megablast templates |
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
| **K-mer size** | k = 5–16 (`uint16_t` for k ≤ 8, `uint32_t` for k ≥ 9; selected by `2*k` bit width) | k = 1–28 (default 15) |
| **Seeding** | All k-mers indexed in a direct-address table (4^k entries) | Minimizers (subsampled k-mers) indexed in a hash table |
| **Candidate filtering** | Explicit Stage 1: scan ID posting lists to score and filter candidates before chaining | No separate filtering stage; seeds go directly to chaining |
| **Chaining DP score** | Chain length (anchor count) with a diagonal-deviation constraint (`max_gap`) | Estimated matching bases minus a gap penalty with logarithmic distance term |
| **Hits per subject** | Configurable via `-stage2_max_nhit_per_subject` (default: 1 best chain; set >1 or 0 for unlimited) | Multiple chains per subject (primary, secondary, supplementary) |
| **Default result limit** | Unlimited (returns all candidates above threshold) | 1 primary + up to 5 secondary per query (`-N 5`, `-p 0.8`) |
| **Alignment** | Parasail semi-global, 1-piece affine gap (optional Stage 3) | KSW2 with SIMD, 2-piece affine gap and Z-drop heuristic |
| **Ambiguous bases** | IUPAC degenerate expansion (configurable limit) | Not supported (N-containing k-mers skipped) |
| **Output format** | TSV, JSON, SAM, BAM | PAF, SAM |
| **Index partitioning** | BLAST DB multi-volume with `.kvx` manifest | `-I` batch partitioning with `--split-prefix` |

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

- C++17 compiler (GCC >= 10, Clang >= 12)
- CMake >= 3.16
- NCBI C++ Toolkit (for BLAST DB access)
- Intel TBB (for parallelization)
- Parasail >= 2.6 (for Stage 3 pairwise alignment)
- htslib >= 1.17 (for SAM/BAM output)
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
