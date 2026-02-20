#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KixReader;
class OidFilter;

struct Stage1Config {
    uint32_t max_freq = 0;          // 0 = auto-calculate
    uint32_t stage1_topn = 500;
    uint32_t min_stage1_score = 2;
};

// Compute effective max_freq from config and index statistics.
// If config_max_freq > 0, returns it unchanged.
// Otherwise, computes: mean_count * 10, clamped to [1000, 100000].
uint32_t compute_effective_max_freq(uint32_t config_max_freq,
                                    uint64_t total_postings,
                                    uint64_t table_size);

// Run Stage 1 filtering: scan ID postings for each query k-mer,
// accumulate hit counts per seq_id, select top candidates.
//
// Template parameter KmerInt: uint16_t or uint32_t depending on k.
// query_kmers: pairs of (q_pos, kmer_value) from KmerScanner.
// Returns candidate seq_ids sorted by score descending.
template <typename KmerInt>
std::vector<SeqId> stage1_filter(
    const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config);

// Explicit instantiations declared
extern template std::vector<SeqId> stage1_filter<uint16_t>(
    const std::vector<std::pair<uint32_t, uint16_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&);
extern template std::vector<SeqId> stage1_filter<uint32_t>(
    const std::vector<std::pair<uint32_t, uint32_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&);

} // namespace ikafssn
