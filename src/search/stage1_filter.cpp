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

uint32_t compute_effective_max_freq(uint32_t config_max_freq,
                                    uint64_t total_postings,
                                    uint64_t table_size) {
    if (config_max_freq > 0) return config_max_freq;
    double mean = static_cast<double>(total_postings) /
                  static_cast<double>(table_size);
    uint32_t max_freq = static_cast<uint32_t>(mean * 10.0);
    if (max_freq < 1000) max_freq = 1000;
    if (max_freq > 100000) max_freq = 100000;
    return max_freq;
}

template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer* buf) {

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0) return {};

    const uint64_t* offsets = kix.offsets();
    const uint32_t* counts = kix.counts();
    const uint8_t* posting_data = kix.posting_data();

    const bool use_coverscore = (config.stage1_score_type == 1);

    if (buf) {
        // Reusable buffer path: dirty list avoids full-array alloc/clear
        buf->ensure_capacity(num_seqs);

        for (const auto& [q_pos, kmer] : query_kmers) {
            uint64_t kmer_idx = static_cast<uint64_t>(kmer);
            uint32_t cnt = counts[kmer_idx];
            if (cnt == 0) continue;

            SeqIdDecoder decoder(posting_data + offsets[kmer_idx]);
            for (uint32_t i = 0; i < cnt; i++) {
                SeqId sid = decoder.next();
                if (use_coverscore && !decoder.was_new_seq()) continue;
                if (!filter.pass(sid)) continue;
                if (buf->score_per_seq[sid] == 0) buf->dirty.push_back(sid);
                if (buf->last_scored_pos[sid] != q_pos) {
                    buf->score_per_seq[sid]++;
                    buf->last_scored_pos[sid] = q_pos;
                }
            }
        }

        // Collect candidates from dirty list only
        std::vector<Stage1Candidate> candidates;
        for (uint32_t sid : buf->dirty) {
            if (buf->score_per_seq[sid] >= config.min_stage1_score) {
                candidates.push_back({sid, buf->score_per_seq[sid]});
            }
        }

        buf->clear_dirty();

        if (config.stage1_topn == 0) {
            return candidates;
        }

        auto cmp = [](const Stage1Candidate& a, const Stage1Candidate& b) {
            return a.score > b.score;
        };

        if (candidates.size() > config.stage1_topn) {
            std::nth_element(candidates.begin(),
                             candidates.begin() + config.stage1_topn,
                             candidates.end(), cmp);
            candidates.resize(config.stage1_topn);
        }
        std::sort(candidates.begin(), candidates.end(), cmp);

        return candidates;
    }

    // Fallback: allocate local buffer (existing behavior)
    std::vector<uint32_t> score_per_seq(num_seqs, 0);
    std::vector<uint32_t> last_scored_pos(num_seqs, UINT32_MAX);

    for (const auto& [q_pos, kmer] : query_kmers) {
        uint64_t kmer_idx = static_cast<uint64_t>(kmer);
        uint32_t cnt = counts[kmer_idx];
        if (cnt == 0) continue;

        SeqIdDecoder decoder(posting_data + offsets[kmer_idx]);
        for (uint32_t i = 0; i < cnt; i++) {
            SeqId sid = decoder.next();
            if (use_coverscore && !decoder.was_new_seq()) continue;
            if (!filter.pass(sid)) continue;
            if (last_scored_pos[sid] != q_pos) {
                score_per_seq[sid]++;
                last_scored_pos[sid] = q_pos;
            }
        }
    }

    // Collect candidates with min_stage1_score
    std::vector<Stage1Candidate> candidates;
    for (uint32_t oid = 0; oid < num_seqs; oid++) {
        if (score_per_seq[oid] >= config.min_stage1_score) {
            candidates.push_back({oid, score_per_seq[oid]});
        }
    }

    if (config.stage1_topn == 0) {
        // Unlimited: return all candidates, skip sort
        return candidates;
    }

    auto cmp = [](const Stage1Candidate& a, const Stage1Candidate& b) {
        return a.score > b.score;
    };

    if (candidates.size() > config.stage1_topn) {
        std::nth_element(candidates.begin(),
                         candidates.begin() + config.stage1_topn,
                         candidates.end(), cmp);
        candidates.resize(config.stage1_topn);
    }
    std::sort(candidates.begin(), candidates.end(), cmp);

    return candidates;
}

// Explicit template instantiations
template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const std::vector<std::pair<uint32_t, uint16_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);
template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const std::vector<std::pair<uint32_t, uint32_t>>&,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);

} // namespace ikafssn
