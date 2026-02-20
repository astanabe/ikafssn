#include "search/stage1_filter.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace ikafssn {

template <typename KmerInt>
std::vector<SeqId> stage1_filter(
    const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config) {

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0) return {};

    // Compute effective max_freq
    uint32_t max_freq = config.max_freq;
    if (max_freq == 0) {
        // Auto: mean_count * 10, clamped to [1000, 100000]
        double mean = static_cast<double>(kix.total_postings()) /
                      static_cast<double>(kix.table_size());
        max_freq = static_cast<uint32_t>(mean * 10.0);
        if (max_freq < 1000) max_freq = 1000;
        if (max_freq > 100000) max_freq = 100000;
    }

    // score_per_seq: hit count per OID
    std::vector<uint32_t> score_per_seq(num_seqs, 0);

    const uint64_t* offsets = kix.offsets();
    const uint32_t* counts = kix.counts();
    const uint8_t* posting_data = kix.posting_data();

    for (const auto& [q_pos, kmer] : query_kmers) {
        uint64_t kmer_idx = static_cast<uint64_t>(kmer);
        uint32_t cnt = counts[kmer_idx];
        if (cnt == 0 || cnt > max_freq) continue;

        SeqIdDecoder decoder(posting_data + offsets[kmer_idx]);
        for (uint32_t i = 0; i < cnt; i++) {
            SeqId sid = decoder.next();
            if (!filter.pass(sid)) continue;
            score_per_seq[sid]++;
        }
    }

    // Collect candidates with min_stage1_score
    struct Candidate {
        SeqId id;
        uint32_t score;
    };
    std::vector<Candidate> candidates;
    for (uint32_t oid = 0; oid < num_seqs; oid++) {
        if (score_per_seq[oid] >= config.min_stage1_score) {
            candidates.push_back({oid, score_per_seq[oid]});
        }
    }

    // Sort by score descending
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) {
                  return a.score > b.score;
              });

    // Truncate to stage1_topn
    if (candidates.size() > config.stage1_topn) {
        candidates.resize(config.stage1_topn);
    }

    // Extract seq_ids
    std::vector<SeqId> result;
    result.reserve(candidates.size());
    for (const auto& c : candidates) {
        result.push_back(c.id);
    }
    return result;
}

// Explicit template instantiations
template std::vector<SeqId> stage1_filter<uint16_t>(
    const std::vector<std::pair<uint32_t, uint16_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&);
template std::vector<SeqId> stage1_filter<uint32_t>(
    const std::vector<std::pair<uint32_t, uint32_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&);

} // namespace ikafssn
