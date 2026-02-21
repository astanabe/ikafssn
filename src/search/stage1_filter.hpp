#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KixReader;
class OidFilter;

struct Stage1Buffer {
    std::vector<uint32_t> score_per_seq;
    std::vector<uint32_t> last_scored_pos;  // per-position dedup for degenerate expansion
    std::vector<uint32_t> dirty;

    void ensure_capacity(uint32_t num_seqs) {
        if (score_per_seq.size() < num_seqs) {
            score_per_seq.resize(num_seqs, 0);
            last_scored_pos.resize(num_seqs, UINT32_MAX);
        }
    }

    void clear_dirty() {
        for (uint32_t idx : dirty) {
            score_per_seq[idx] = 0;
            last_scored_pos[idx] = UINT32_MAX;
        }
        dirty.clear();
    }
};

struct Stage1Candidate {
    SeqId id;
    uint32_t score;
};

struct Stage1Config {
    uint32_t max_freq = 0;          // 0 = auto-calculate
    uint32_t stage1_topn = 0;       // 0 = unlimited (skip sort)
    uint32_t min_stage1_score = 1;
    uint8_t  stage1_score_type = 1; // 1 = coverscore, 2 = matchscore
};

// Compute effective max_freq from config and index statistics.
// If config_max_freq > 0, returns it unchanged.
// Otherwise, computes: mean_count * 10, clamped to [1000, 100000].
uint32_t compute_effective_max_freq(uint32_t config_max_freq,
                                    uint64_t total_postings,
                                    uint64_t table_size);

// Run Stage 1 filtering: scan ID postings for each query k-mer,
// accumulate scores per seq_id, select top candidates.
//
// stage1_score_type: 1 = coverscore (distinct query k-mers matching per seq),
//                    2 = matchscore (total k-mer position matches per seq).
// stage1_topn: 0 = unlimited (return all candidates, skip sort).
//
// Template parameter KmerInt: uint16_t or uint32_t depending on k.
// query_kmers: pairs of (q_pos, kmer_value) from KmerScanner.
// Returns candidates with scores, sorted by score desc (unless topn=0).
template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer* buf = nullptr);

// Explicit instantiations declared
extern template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const std::vector<std::pair<uint32_t, uint16_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);
extern template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const std::vector<std::pair<uint32_t, uint32_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);

} // namespace ikafssn
