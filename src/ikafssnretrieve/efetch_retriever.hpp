#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

struct EfetchOptions {
    std::string api_key;           // NCBI API key (empty = none)
    uint32_t batch_size = 100;     // accessions per batch request
    uint32_t retries = 3;          // max retry count
    uint32_t timeout_sec = 30;     // request timeout
    uint32_t range_threshold = 100000;  // seq_length threshold for individual fetch
    uint32_t context = 0;          // bases to add before/after match region
};

// Information about a single fetch target, grouped from hits.
struct FetchTarget {
    std::string accession;
    uint32_t s_start;
    uint32_t s_end;
    uint32_t seq_length_estimate;  // estimated from s_end
    // Original hit info for FASTA header
    std::string qseqid;
    char strand;
    uint32_t score;
};

// Build efetch URL for batch retrieval (full sequences).
std::string build_efetch_url_batch(const std::vector<std::string>& accessions,
                                   const std::string& api_key);

// Build efetch URL for individual retrieval with range (1-based inclusive).
std::string build_efetch_url_range(const std::string& accession,
                                   uint32_t seq_start, uint32_t seq_stop,
                                   const std::string& api_key);

// Parse a multi-FASTA efetch response into (accession -> sequence) pairs.
// Extracts accession from the defline (first word after '>').
struct EfetchRecord {
    std::string accession;
    std::string sequence;
};
std::vector<EfetchRecord> parse_efetch_response(const std::string& response);

// Get rate limit sleep duration in milliseconds.
uint32_t rate_limit_sleep_ms(bool has_api_key);

// Check if an HTTP status code is retryable.
bool is_retryable_http_status(long status_code);

// Check if an HTTP status code means the accession should be skipped.
bool is_skip_http_status(long status_code);

#ifdef IKAFSSN_ENABLE_REMOTE

// Retrieve matched subsequences from NCBI efetch.
// Writes FASTA records to the output stream.
// Returns the number of successfully retrieved sequences.
uint32_t retrieve_remote(const std::vector<OutputHit>& hits,
                         const EfetchOptions& opts,
                         std::ostream& out);

#endif // IKAFSSN_ENABLE_REMOTE

} // namespace ikafssn
