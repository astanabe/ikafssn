# ikafssn

**ikafssn** (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) is a suite of command-line tools that builds inverted indexes over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.

## Features

- Builds a k-mer inverted index directly from NCBI BLAST databases
- Two-stage search pipeline: fast candidate filtering (Stage 1) followed by position-aware collinear chaining (Stage 2), with configurable scoring (coverscore/matchscore/chainscore) and optional Stage 1-only mode
- Client-server architecture with UNIX/TCP socket and HTTP REST API support, with multi-database serving from a single process
- Handles IUPAC ambiguous bases by expanding single-ambiguity k-mers during indexing
- Parallel indexing and search via Intel TBB
- Lightweight per-command executables, each linking only its required dependencies

## Commands

| Command | Purpose |
|---|---|
| `ikafssnindex` | Build k-mer inverted index from BLAST DB |
| `ikafssnsearch` | Local direct search (mmap index) |
| `ikafssnretrieve` | Extract matched subsequences (local or NCBI efetch) |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver |
| `ikafssnclient` | Client (socket or HTTP) |
| `ikafssninfo` | Index/DB information display |

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

This project is licensed under the GNU General Public License v2.0. See [LICENSE](LICENSE) for details.
