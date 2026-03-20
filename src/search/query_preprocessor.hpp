#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class KixReader;
class KhxReader;
struct SearchConfig;

// Pre-processed query k-mer data with global high-freq filtering applied.
// Generated once per query before the volume loop.
template <typename KmerInt>
struct QueryKmerData {
    std::vector<uint32_t> fwd_positions;  // query positions (high-freq removed)
    std::vector<KmerInt>  fwd_kmer_values; // k-mer values (high-freq removed)
    std::vector<uint32_t> rc_positions;
    std::vector<KmerInt>  rc_kmer_values;
    uint32_t resolved_threshold_fwd = 0;  // resolved Stage 1 absolute threshold (fwd)
    uint32_t resolved_threshold_rc = 0;   // resolved Stage 1 absolute threshold (rc)
    uint32_t effective_min_score_fwd = 0;  // for Stage 2 (fwd)
    uint32_t effective_min_score_rc = 0;   // for Stage 2 (rc)
    bool has_multi_degen = false;  // true if any k-mer had 2+ degenerate bases (skipped)
};

// Pre-process a query sequence: extract k-mers, determine global high-freq
// k-mers across all volumes, filter them out, and resolve per-strand thresholds.
//
// all_kix: pointers to KixReaders for ALL volumes (for global count aggregation).
// khx: nullable pointer to shared KhxReader for build-time exclusion info.
template <typename KmerInt>
QueryKmerData<KmerInt> preprocess_query(
    const std::string& query_seq, int k,
    const std::vector<const KixReader*>& all_kix,
    const KhxReader* khx,
    const SearchConfig& config,
    uint8_t t = 0,
    const std::vector<uint32_t>& masks = {});

extern template QueryKmerData<uint16_t> preprocess_query<uint16_t>(
    const std::string&, int,
    const std::vector<const KixReader*>&,
    const KhxReader*,
    const SearchConfig&,
    uint8_t,
    const std::vector<uint32_t>&);
extern template QueryKmerData<uint32_t> preprocess_query<uint32_t>(
    const std::string&, int,
    const std::vector<const KixReader*>&,
    const KhxReader*,
    const SearchConfig&,
    uint8_t,
    const std::vector<uint32_t>&);

} // namespace ikafssn
