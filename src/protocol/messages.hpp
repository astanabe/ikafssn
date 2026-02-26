#pragma once

#include <climits>
#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Seqidlist filter mode in protocol messages
enum class SeqidlistMode : uint8_t {
    kNone     = 0,  // no filter
    kInclude  = 1,  // seqidlist (include only)
    kExclude  = 2,  // negative_seqidlist (exclude)
};

// Query sequence within a search request
struct QueryEntry {
    std::string query_id;
    std::string sequence;
};

// Search request message (client -> server)
struct SearchRequest {
    uint8_t  k = 0;                  // 0 = server default
    uint16_t stage2_min_score = 0;   // 0 = server default (see has_stage2_min_score)
    uint16_t stage2_max_gap = 0;     // 0 = server default
    uint32_t stage1_max_freq = 0;    // 0 = server default
    uint8_t  stage2_min_diag_hits = 0; // 0 = server default
    uint16_t stage1_topn = 0;        // 0 = server default
    uint16_t stage1_min_score = 0;   // 0 = server default
    uint16_t num_results = 0;        // 0 = server default
    uint16_t stage1_min_score_frac_x10000 = 0; // P * 10000, 0 = use integer field
    uint16_t stage1_max_freq_frac_x10000 = 0;  // P * 10000, 0 = use integer field
    SeqidlistMode seqidlist_mode = SeqidlistMode::kNone;
    uint8_t  mode = 0;              // 0 = server default, 1 = stage1 only, 2 = stage1+stage2
    uint8_t  stage1_score = 0;      // 0 = server default, 1 = coverscore, 2 = matchscore
    uint8_t  accept_qdegen = 1;     // 0 = reject degenerate queries, 1 = accept
    int8_t   strand = 0;            // 0 = server default, 1 = plus, -1 = minus, 2 = both
    uint8_t  has_stage2_min_score = 0; // 1 = stage2_min_score was explicitly set by client
    uint16_t stage2_max_lookback = 0;  // 0 = server default
    // Stage 3 parameters (mode 3)
    uint8_t  stage3_traceback = 0;                // 0 = score only, 1 = full traceback
    int16_t  stage3_gapopen = INT16_MIN;          // INT16_MIN = server default
    int16_t  stage3_gapext = INT16_MIN;           // INT16_MIN = server default
    uint16_t stage3_min_pident_x100 = 0;          // pident * 100; 0 = no filter
    uint32_t stage3_min_nident = 0;               // 0 = no filter
    uint32_t context_abs = 0;                     // absolute context bases (when frac=0)
    uint16_t context_frac_x10000 = 0;             // ratio * 10000 (when > 0, ratio mode)
    std::string db_name;                           // target database name (empty = error)
    std::vector<std::string> seqids;
    std::vector<QueryEntry> queries;
};

// A single hit in the search response
struct ResponseHit {
    std::string accession;
    uint8_t  strand;     // 0 = '+', 1 = '-'
    uint32_t q_start;
    uint32_t q_end;
    uint32_t q_length = 0;    // query full sequence length
    uint32_t s_start;
    uint32_t s_end;
    uint32_t s_length = 0;    // subject full sequence length
    uint16_t score;           // chainscore
    uint16_t stage1_score;
    uint16_t volume;
    // Stage 3 fields (populated only when mode=3)
    int32_t  alnscore = 0;
    uint32_t nident = 0;
    uint32_t nmismatch = 0;
    uint16_t pident_x100 = 0;  // percent identity * 100
    std::string cigar;          // CIGAR string (traceback only)
    std::string q_seq;          // aligned query (traceback only)
    std::string s_seq;          // aligned subject (traceback only)
};

// Per-query warning flags (bitmask)
enum QueryWarning : uint8_t {
    kWarnMultiDegen = 0x01,  // k-mers with 2+ degenerate bases were skipped
};

// Per-query result in the search response
struct QueryResult {
    std::string query_id;
    std::vector<ResponseHit> hits;
    uint8_t skipped = 0;   // 0 = normal, 1 = skipped (degenerate bases)
    uint8_t warnings = 0;  // bitmask of QueryWarning flags
};

// Search response message (server -> client)
struct SearchResponse {
    uint8_t  status = 0;  // 0 = success
    uint8_t  k = 0;
    uint8_t  mode = 2;              // 1 = stage1 only, 2 = stage1+stage2, 3 = stage1+stage2+stage3
    uint8_t  stage1_score = 1;      // 1 = coverscore, 2 = matchscore
    uint8_t  stage3_traceback = 0;  // echo back: 1 = traceback fields populated
    std::string db_name;            // echo back database name
    std::vector<QueryResult> results;
    std::vector<std::string> rejected_query_ids;  // queries rejected due to concurrency limit
};

// Error response message (server -> client)
struct ErrorResponse {
    uint32_t error_code = 0;
    std::string message;
};

// Health request: empty payload
struct HealthRequest {};

// Health response: simple status
struct HealthResponse {
    uint8_t status = 0; // 0 = OK
};

// Info request: empty payload
struct InfoRequest {};

// Per-volume info in the info response
struct VolumeInfo {
    uint16_t volume_index;
    uint32_t num_sequences;
    uint64_t total_postings;
    std::string db_name;
};

// Per-k group info in the info response
struct KmerGroupInfo {
    uint8_t  k;
    uint8_t  kmer_type; // 0 = uint16, 1 = uint32
    std::vector<VolumeInfo> volumes;
};

// Per-database info in the info response
struct DatabaseInfo {
    std::string name;
    uint8_t default_k = 0;
    uint8_t max_mode = 2;   // 1=stage1 only, 2=stage1+2, 3=stage1+2+3
    std::vector<KmerGroupInfo> groups;
};

// Info response message (server -> client)
struct InfoResponse {
    uint8_t  status = 0;           // 0 = success
    uint8_t  default_k = 0;       // global default (first DB's)
    int32_t  max_active_sequences = 0;
    int32_t  active_sequences = 0; // current in-use
    std::vector<DatabaseInfo> databases;
};

} // namespace ikafssn
