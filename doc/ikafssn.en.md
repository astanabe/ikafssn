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
| `ikafssninfo` | Index/DB information display |

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

### ikafssnindex

Build a k-mer inverted index from a BLAST database. For each volume, three index files are generated: `.kix` (ID postings), `.kpx` (position postings), and `.ksx` (sequence metadata). When `-max_freq_build` is used, a `.khx` file (build-time exclusion bitset) is also generated.

```
ikafssnindex [options]

Required:
  -db <path>              BLAST DB prefix
  -k <int>                K-mer length (5-13)
  -o <dir>                Output directory

Options:
  -memory_limit <size>    Memory limit (default: half of physical RAM)
                          Accepts K, M, G suffixes
                          Partitions are auto-calculated to fit within this limit
  -max_freq_build <num>   Exclude k-mers with count above this threshold
                          >= 1: absolute count threshold
                          0 < x < 1: fraction of NSEQ per volume
                          (default: 0 = no exclusion)
  -openvol <int>          Max volumes processed simultaneously (default: 1)
                          Controls peak memory usage for multi-volume DBs
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

# Large DB, allow 2 volumes to be processed simultaneously
ikafssnindex -db nt -k 11 -o ./nt_index -openvol 2

# Exclude high-frequency k-mers during build (absolute)
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 50000

# Exclude k-mers appearing in >1% of sequences per volume
ikafssnindex -db nt -k 11 -o ./nt_index -max_freq_build 0.01
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
  -min_score <int>        Minimum score (default: 1)
  -max_gap <int>          Chaining diagonal gap tolerance (default: 100)
  -max_freq <int>         High-frequency k-mer skip threshold (default: auto)
  -min_diag_hits <int>    Diagonal filter min hits (default: 2)
  -stage1_topn <int>      Stage 1 candidate limit, 0=unlimited (default: 0)
  -min_stage1_score <num> Stage 1 minimum score (default: 0.5)
                          Integer (>= 1): absolute threshold
                          Fraction (0 < P < 1): proportion of query k-mers,
                            resolved per query as ceil(Nqkmer * P) - Nhighfreq
  -num_results <int>      Max results per query, 0=unlimited (default: 0)
  -mode <1|2>             Search mode (default: 2)
                          1=Stage 1 only, 2=Stage 1 + Stage 2
  -stage1_score <1|2>     Stage 1 score type (default: 1)
                          1=coverscore, 2=matchscore
  -sort_score <1|2>       Sort key for results (default: 2)
                          1=stage1 score, 2=chainscore
  -seqidlist <path>       Include only listed accessions
  -negative_seqidlist <path>  Exclude listed accessions
  -strand <-1|1|2>       Strand to search (default: 2)
                          1=plus only, -1=minus only, 2=both
  -accept_qdegen <0|1>    Accept queries with degenerate bases (default: 0)
  -outfmt <tab|json>      Output format (default: tab)
  -v, --verbose           Verbose logging
```

The `-ix` option specifies the index prefix path (without extension), similar to `blastn -db`. For example, if `ikafssnindex -db nt -k 11 -o /data/index` generated the following files:

```
/data/index/nt.00.11mer.kix
/data/index/nt.00.11mer.kpx
/data/index/nt.00.11mer.ksx
/data/index/nt.01.11mer.kix
/data/index/nt.01.11mer.kpx
/data/index/nt.01.11mer.ksx
```

then specify `-ix /data/index/nt`. The prefix `/data/index/nt` is split into the directory `/data/index/` and the base name `nt`, and all files matching `nt.<vol>.<k>mer.k?x` in that directory are discovered automatically.

If the index directory contains indexes for multiple k-mer sizes (e.g. both `nt.00.09mer.kix` and `nt.00.11mer.kix`), you must specify `-k` to select which one to use. If only a single k-mer size exists, `-k` can be omitted.

When `-accept_qdegen` is 0 (default), queries containing IUPAC degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) are skipped with a warning, and the exit code is 2. Set `-accept_qdegen 1` to allow such queries (degenerate positions are treated as invalid and skipped during k-mer scanning).

`-seqidlist` and `-negative_seqidlist` are mutually exclusive. Both text (one accession per line) and binary (generated by `blastdb_aliastool -seqid_file_in`) formats are accepted, auto-detected by magic bytes.

**Examples:**

```bash
# Basic search
ikafssnsearch -ix ./index/mydb -query query.fasta -threads 8

# Specify k-mer size (required if index contains multiple k values)
ikafssnsearch -ix ./index/mydb -k 11 -query query.fasta

# Increase sensitivity
ikafssnsearch -ix ./index/mydb -query query.fasta \
    -min_score 2 -stage1_topn 2000 -max_freq 50000

# Restrict to specific accessions
ikafssnsearch -ix ./index/mydb -query query.fasta -seqidlist targets.txt

# Exclude specific accessions
ikafssnsearch -ix ./index/mydb -query query.fasta -negative_seqidlist exclude.txt

# Fractional Stage 1 threshold (50% of query k-mers)
ikafssnsearch -ix ./index/mydb -query query.fasta -min_stage1_score 0.5

# Pipe to ikafssnretrieve
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -db nt > matches.fasta
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
  -context <int>          Bases to add before/after match region (default: 0)
  -v, --verbose           Verbose logging

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
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# Remote retrieval via NCBI efetch
ikafssnsearch -ix ./index/mydb -query query.fasta | ikafssnretrieve -remote > matches.fasta

# Remote retrieval with API key (higher throughput)
ikafssnclient -http http://search.example.com:8080 -query query.fasta \
    | ikafssnretrieve -remote -api_key XXXXXXXX > matches.fasta

# Include 50bp of context around match region
ikafssnretrieve -db nt -results results.tsv -context 50
```

### ikafssnserver

Search daemon. Keeps index files memory-mapped and accepts search requests via UNIX domain socket or TCP socket.

```
ikafssnserver [options]

Required:
  -ix <prefix>            Index prefix (like blastn -db)

Listener (at least one required):
  -socket <path>          UNIX domain socket path
  -tcp <host>:<port>      TCP listen address

Options:
  -threads <int>          Worker threads (default: all cores)
  -max_query <int>        Max concurrent query sequences globally (default: 1024)
  -max_seqs_per_req <int> Max sequences accepted per request (default: thread count)
  -pid <path>             PID file path
  -min_score <int>        Default minimum chain score (default: 1)
  -max_gap <int>          Default chaining gap tolerance (default: 100)
  -max_freq <int>         Default high-freq k-mer skip threshold (default: auto)
  -min_diag_hits <int>    Default diagonal filter min hits (default: 2)
  -stage1_topn <int>      Default Stage 1 candidate limit (default: 0)
  -min_stage1_score <num> Default Stage 1 minimum score (default: 0.5)
                          Integer (>= 1) or fraction (0 < P < 1)
  -num_results <int>      Default max results per query (default: 0)
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
```

**Operational characteristics:**

- One process serves one BLAST DB index. To serve multiple DBs, start separate processes.
- If the index prefix matches indexes for multiple k-mer sizes, all are loaded and clients can specify k per request.
- On SIGTERM/SIGINT, performs graceful shutdown: stops accepting new connections, waits for in-flight requests to complete (up to `-shutdown_timeout` seconds), then exits.
- **Per-sequence concurrency control:** The server limits concurrency at the per-sequence level, not per-connection. When a request arrives, the server attempts to acquire permits for each valid query sequence. If the global limit (`-max_query`) is reached, excess sequences are returned to the client as "rejected" for retry. The `-max_seqs_per_req` option caps how many permits a single request can acquire, preventing one large request from monopolizing all slots.

### ikafssnhttpd

HTTP REST API daemon. Connects to `ikafssnserver` and exposes search as an HTTP API. Uses the Drogon framework.

```
ikafssnhttpd [options]

Backend connection (one required):
  -server_socket <path>      UNIX socket path to ikafssnserver
  -server_tcp <host>:<port>  TCP address of ikafssnserver

Options:
  -listen <host>:<port>  HTTP listen address (default: 0.0.0.0:8080)
  -path_prefix <prefix>  API path prefix (e.g., /nt)
  -threads <int>         Drogon I/O threads (default: all cores)
  -pid <path>            PID file path
  -v, --verbose          Verbose logging
```

**REST API endpoints:**

| Method | Path | Description |
|---|---|---|
| POST | `/api/v1/search` | Search request (query sequences in JSON body) |
| GET | `/api/v1/info` | Index information |
| GET | `/api/v1/health` | Health check |

**Examples:**

```bash
# Connect to local ikafssnserver via UNIX socket
ikafssnhttpd -server_socket /var/run/ikafssn.sock -listen 0.0.0.0:8080

# Connect to remote ikafssnserver via TCP
ikafssnhttpd -server_tcp 10.0.1.5:9100 -listen 0.0.0.0:8080

# Multiple DBs with path prefix (pair with nginx routing)
ikafssnhttpd -server_socket /var/run/ikafssn_nt.sock -listen :8080 -path_prefix /nt
ikafssnhttpd -server_socket /var/run/ikafssn_rs.sock -listen :8081 -path_prefix /rs
```

### ikafssnclient

Client command. Connects to `ikafssnserver` via socket or `ikafssnhttpd` via HTTP. Output format is identical to `ikafssnsearch`. If the server rejects some query sequences due to concurrency limits, the client automatically retries the rejected queries with exponential backoff (30s, 60s, 120s, 120s, ...) until all queries are processed.

```
ikafssnclient [options]

Connection (one required):
  -socket <path>           UNIX domain socket path
  -tcp <host>:<port>       TCP server address
  -http <url>              ikafssnhttpd URL (e.g., http://example.com:8080)

Required:
  -query <path>            Query FASTA file (- for stdin)

Options:
  -o <path>                Output file (default: stdout)
  -k <int>                 K-mer size (default: server default)
  -min_score <int>         Minimum score (default: server default)
  -max_gap <int>           Chaining gap tolerance (default: server default)
  -max_freq <int>          High-freq k-mer skip threshold (default: server default)
  -min_diag_hits <int>     Diagonal filter min hits (default: server default)
  -stage1_topn <int>       Stage 1 candidate limit (default: server default)
  -min_stage1_score <num>  Stage 1 minimum score (default: server default)
                           Integer (>= 1) or fraction (0 < P < 1)
  -num_results <int>       Max results per query (default: server default)
  -mode <1|2>              Search mode (default: server default)
  -stage1_score <1|2>      Stage 1 score type (default: server default)
  -sort_score <1|2>        Sort key for results (default: server default)
  -seqidlist <path>        Include only listed accessions
  -negative_seqidlist <path>  Exclude listed accessions
  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: server default)
  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 0)
  -outfmt <tab|json>       Output format (default: tab)
  -v, --verbose            Verbose logging
```

**Examples:**

```bash
# UNIX socket (local)
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta

# TCP (local or remote)
ikafssnclient -tcp 10.0.1.5:9100 -query query.fasta

# HTTP
ikafssnclient -http http://search.example.com:8080 -query query.fasta

# Restrict search scope
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta -seqidlist targets.txt

# Pipe to ikafssnretrieve
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta | ikafssnretrieve -db nt > matches.fasta

# Specify k-mer size
ikafssnclient -socket /var/run/ikafssn.sock -query query.fasta -k 9
```

### ikafssninfo

Display index statistics and optionally BLAST DB information.

```
ikafssninfo [options]

Required:
  -ix <prefix>            Index prefix (like blastn -db)

Options:
  -db <path>              BLAST DB prefix (also show DB info)
  -v, --verbose           Verbose output (k-mer frequency distribution, etc.)
```

**Output includes:**

- K-mer length (k) and integer type (uint16/uint32)
- Number of volumes
- Per-volume statistics: sequence count, total postings, file sizes, excluded k-mer count (if `.khx` present)
- Overall statistics: total sequences, total postings, total index size, compression ratio
- With `-v`: k-mer frequency distribution (min, max, mean, percentiles)
- With `-db`: BLAST DB title, sequence count, total bases, volume paths

**Examples:**

```bash
# Basic index info
ikafssninfo -ix ./index/mydb

# Include BLAST DB info
ikafssninfo -ix ./index/mydb -db mydb

# Detailed frequency distribution
ikafssninfo -ix ./index/mydb -v
```

## Search Pipeline

ikafssn uses a two-stage search pipeline:

The default parameters prioritize throughput: `stage1_topn=0` and `num_results=0` disable sorting, and `min_stage1_score=0.5` (fractional) filters candidates by requiring at least 50% of query k-mers to match. To get ranked output, set positive values for `-stage1_topn` and/or `-num_results`, which triggers sorting but may reduce speed for large result sets.

1. **Stage 1 (Candidate Selection):** Scans ID postings for each query k-mer and accumulates scores per sequence. Two score types are available: **coverscore** (number of distinct query k-mers matching the sequence) and **matchscore** (total k-mer position matches). Sequences exceeding `min_stage1_score` are selected as candidates. When `stage1_topn > 0`, candidates are sorted by score and truncated. When `stage1_topn = 0` (default), all qualifying candidates are returned without sorting.

2. **Stage 2 (Collinear Chaining):** For each candidate, collects position-level hits from the `.kpx` file, applies a diagonal filter, and runs a chaining DP to find the best collinear chain. The chain length is reported as **chainscore**. Chains with `chainscore >= min_score` are reported.

**Mode 1 (Stage 1 only):** When `-mode 1` is specified, Stage 2 is skipped entirely. The `.kpx` file is not accessed, saving I/O and memory. Results contain only Stage 1 scores; position fields (q_start, q_end, s_start, s_end) and chainscore are omitted. In mode 1, `-min_score` applies to the Stage 1 score, and `-sort_score` is forced to 1 (stage1 score).

By default, both forward and reverse complement strands of the query are searched. Use `-strand 1` to search only the plus (forward) strand, or `-strand -1` to search only the minus (reverse complement) strand.

### High-Frequency K-mer Filtering

K-mers occurring more than `max_freq` times are skipped during both stages. When `max_freq` is not specified (default), it is auto-calculated as:

```
max_freq = mean_count * 10    (clamped to [1000, 100000])
where mean_count = total_postings / 4^k
```

This is computed per volume from the `.kix` header.

**Build-time exclusion** (`-max_freq_build`): When indexing with `-max_freq_build`, high-frequency k-mers are excluded from the index entirely. A `.khx` file records which k-mers were excluded. At search time, when fractional `-min_stage1_score` is used, k-mers excluded at build time are recognized from the `.khx` file and subtracted from the threshold calculation.

### Fractional Stage 1 Threshold

When `-min_stage1_score` is specified as a fraction (0 < P < 1), the threshold is resolved per query as:

```
threshold = ceil(Nqkmer * P) - Nhighfreq
```

Where:
- **Nqkmer**: number of query k-mers (distinct for coverscore, total positions for matchscore)
- **Nhighfreq**: number of query k-mers that are excluded, combining:
  - Search-time exclusion: k-mers with count > `max_freq`
  - Build-time exclusion: k-mers marked in `.khx` (if present)

If the resolved threshold is <= 0, the query strand is skipped with a warning.

### Score Types

ikafssn computes three types of scores:

| Score | Description | Computed in |
|---|---|---|
| **coverscore** | Number of distinct query k-mers that match the subject sequence. Each query k-mer is counted at most once per subject, regardless of how many positions it matches. | Stage 1 |
| **matchscore** | Total number of (query k-mer, subject position) matches. A single query k-mer matching multiple positions in the subject counts multiple times. | Stage 1 |
| **chainscore** | Length (number of k-mer hits) of the best collinear chain found by the chaining DP. Requires position data from `.kpx`. | Stage 2 |

- `-stage1_score` selects which score type Stage 1 uses (1=coverscore, 2=matchscore). This affects candidate ranking and the stage1 score reported in output.
- `-sort_score` selects which score is used for final result sorting (1=stage1 score, 2=chainscore).
- In `-mode 1`, only Stage 1 scores are available; chainscore is not computed.

## Output Format

### Tab Format (default)

**Mode 2** (default):

Tab-separated columns, where `coverscore` is replaced by `matchscore` when `-stage1_score 2`:

```
# query_id  accession  strand  q_start  q_end  s_start  s_end  coverscore  chainscore  volume
```

**Mode 1** (`-mode 1`):

```
# query_id  accession  strand  coverscore  volume
```

### JSON Format

**Mode 2** (default):

```json
{
  "results": [
    {
      "query_id": "query1",
      "hits": [
        {
          "accession": "NC_001234.5",
          "strand": "+",
          "q_start": 0,
          "q_end": 150,
          "s_start": 1000,
          "s_end": 1150,
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
      "query_id": "query1",
      "hits": [
        {
          "accession": "NC_001234.5",
          "strand": "+",
          "coverscore": 8,
          "volume": 0
        }
      ]
    }
  ]
}
```

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

Per BLAST DB volume, three files are generated (plus an optional fourth):

```
<db_prefix>.<volume_index>.<kk>mer.kix   — ID postings (direct-address table + delta-compressed)
<db_prefix>.<volume_index>.<kk>mer.kpx   — Position postings (delta-compressed)
<db_prefix>.<volume_index>.<kk>mer.ksx   — Sequence metadata (lengths + accessions)
<db_prefix>.<volume_index>.<kk>mer.khx   — Build-time exclusion bitset (only when -max_freq_build is used)
```

Example: `nt.00.11mer.kix`, `nt.01.11mer.kpx`

The `.khx` file contains a 32-byte header (magic "KMHX", format version, k) followed by a bitset of `ceil(4^k / 8)` bytes. Bit *i* = 1 indicates that k-mer *i* was excluded during index build.

ID and position postings are stored in separate files so that Stage 1 filtering never touches `.kpx`, maximizing page cache efficiency.

## Building from Source

### Dependencies

- C++17 compiler (GCC >= 10, Clang >= 12)
- CMake >= 3.16
- NCBI C++ Toolkit (for BLAST DB access)
- Intel TBB (for parallelization)
- Drogon (for ikafssnhttpd, optional)
- libcurl (for HTTP client mode and remote retrieval, optional)

### Installing Dependencies

Install the required packages (excluding NCBI C++ Toolkit) with the following commands.

**Ubuntu Server 24.04:**

```bash
sudo apt install build-essential cmake libtbb-dev liblmdb-dev libsqlite3-dev \
    libcurl4-openssl-dev libjsoncpp-dev
sudo apt install libdrogon-dev uuid-dev libmariadb-dev libyaml-cpp-dev
```

The second line installs Drogon and its additional dependencies that are not automatically pulled in by `libdrogon-dev` on Ubuntu. If ikafssnhttpd is not needed, omit the second line and build with `-DBUILD_HTTPD=OFF`.

**AlmaLinux 9 / Rocky Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y epel-release
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y libuuid-devel openssl-devel zlib-devel
```

**Oracle Linux 9:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf install -y oracle-epel-release-el9
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y libuuid-devel openssl-devel zlib-devel
```

**AlmaLinux 10 / Rocky Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y libuuid-devel openssl-devel zlib-devel
```

**Oracle Linux 10:**

```bash
sudo dnf config-manager --set-enabled crb
sudo dnf group install -y "Development Tools"
sudo dnf install -y cmake gcc-c++ tbb-devel lmdb-devel sqlite-devel \
    libcurl-devel jsoncpp-devel
sudo dnf install -y libuuid-devel openssl-devel zlib-devel
```

On EL9, `jsoncpp-devel` requires EPEL and `lmdb-devel` requires CRB. On EL10, both are in CRB so EPEL is not needed for these packages. Drogon is not packaged for EL9/EL10; the last line of each block installs dependencies needed to build Drogon from source. If ikafssnhttpd is not needed, omit the last line and build with `-DBUILD_HTTPD=OFF`.

### NCBI C++ Toolkit

ikafssn requires a pre-built NCBI C++ Toolkit. By default, CMake looks for the toolkit at `./ncbi-cxx-toolkit` relative to the source root. If the toolkit is installed elsewhere, specify the path with `-DNCBI_TOOLKIT_DIR`.

The build subdirectory name within the toolkit (e.g. `CMake-GCC1330-Release`) is auto-detected by default but can be overridden with `-DNCBI_TOOLKIT_BUILD_TAG` if needed.

To download, build, and install the toolkit into `./ncbi-cxx-toolkit`, run the following from the ikafssn source root:

```bash
curl -L -o ncbi-cxx-toolkit-30.0.0.tar.gz \
    https://github.com/ncbi/ncbi-cxx-toolkit-public/archive/refs/tags/release/30.0.0.tar.gz
tar xf ncbi-cxx-toolkit-30.0.0.tar.gz
cd ncbi-cxx-toolkit-public-release-30.0.0
./cmake-configure \
    --without-debug \
    --with-projects="objtools/blast/seqdb_reader;objtools/blast/blastdb_format" \
    --with-install="$(realpath ..)/ncbi-cxx-toolkit"
cd CMake-GCC*/build
make -j$(nproc)
make install
cd ../..
```

Only the libraries required by ikafssn (`seqdb`, `blastdb_format`, and their dependencies) are built. The full toolkit build is not necessary.

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
| `BUILD_HTTPD` | ON | Build ikafssnhttpd (requires Drogon) |
| `BUILD_CLIENT` | ON | Build ikafssnclient (requires libcurl for HTTP mode) |
| `ENABLE_REMOTE_RETRIEVE` | ON | Enable NCBI efetch in ikafssnretrieve |
