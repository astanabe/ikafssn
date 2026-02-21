#pragma once

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
    uint16_t min_score = 0;          // 0 = server default (see has_min_score)
    uint16_t max_gap = 0;            // 0 = server default
    uint32_t max_freq = 0;           // 0 = server default
    uint8_t  min_diag_hits = 0;      // 0 = server default
    uint16_t stage1_topn = 0;        // 0 = server default
    uint16_t min_stage1_score = 0;   // 0 = server default
    uint16_t num_results = 0;        // 0 = server default
    uint16_t min_stage1_score_frac_x10000 = 0; // P * 10000, 0 = use integer field
    uint16_t max_freq_frac_x10000 = 0;  // P * 10000, 0 = use integer max_freq field
    SeqidlistMode seqidlist_mode = SeqidlistMode::kNone;
    uint8_t  mode = 0;              // 0 = server default, 1 = stage1 only, 2 = stage1+stage2
    uint8_t  stage1_score_type = 0; // 0 = server default, 1 = coverscore, 2 = matchscore
    uint8_t  sort_score = 0;        // 0 = server default, 1 = stage1 score, 2 = chainscore
    uint8_t  accept_qdegen = 0;     // 0 = reject degenerate queries, 1 = accept
    int8_t   strand = 0;           // 0 = server default, 1 = plus, -1 = minus, 2 = both
    uint8_t  has_min_score = 0;    // 1 = min_score was explicitly set by client
    uint16_t chain_max_lookback = 0;  // 0 = server default
    std::vector<std::string> seqids;
    std::vector<QueryEntry> queries;
};

// A single hit in the search response
struct ResponseHit {
    std::string accession;
    uint8_t  strand;     // 0 = '+', 1 = '-'
    uint32_t q_start;
    uint32_t q_end;
    uint32_t s_start;
    uint32_t s_end;
    uint16_t score;
    uint16_t stage1_score;
    uint16_t volume;
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
    uint8_t  mode = 2;              // 1 = stage1 only, 2 = stage1+stage2
    uint8_t  stage1_score_type = 1; // 1 = coverscore, 2 = matchscore
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

// Info response message (server -> client)
struct InfoResponse {
    uint8_t  status = 0;  // 0 = success
    uint8_t  default_k = 0;
    std::vector<KmerGroupInfo> groups;
};

} // namespace ikafssn
