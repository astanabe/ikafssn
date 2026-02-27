#include "search/volume_searcher.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/posting_decoder.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"
#include "search/query_preprocessor.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "core/kmer_encoding.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>

namespace ikafssn {

template <typename KmerInt>
static std::vector<std::pair<uint32_t, KmerInt>>
extract_kmers(const std::string& seq, int k) {
    std::vector<std::pair<uint32_t, KmerInt>> kmers;
    KmerScanner<KmerInt> scanner(k);
    scanner.scan_ambig(seq.data(), seq.size(),
        [&](uint32_t pos, KmerInt kmer) {
            kmers.emplace_back(pos, kmer);
        },
        [&](uint32_t pos, KmerInt base_kmer, uint8_t ncbi4na, int bit_offset) {
            expand_ambig_kmer<KmerInt>(base_kmer, ncbi4na, bit_offset,
                [&](KmerInt expanded) {
                    kmers.emplace_back(pos, expanded);
                });
        });
    return kmers;
}

// Stage 1 only: return candidates as ChainResult with stage1_score, no chaining.
static std::vector<ChainResult>
stage1_only_results(const std::vector<Stage1Candidate>& candidates,
                    bool is_reverse,
                    uint32_t min_score) {
    std::vector<ChainResult> results;
    for (const auto& c : candidates) {
        if (c.score < min_score) continue;
        ChainResult cr{};
        cr.seq_id = c.id;
        cr.chainscore = 0;
        cr.stage1_score = c.score;
        cr.q_start = 0;
        cr.q_end = 0;
        cr.s_start = 0;
        cr.s_end = 0;
        cr.is_reverse = is_reverse;
        results.push_back(cr);
    }
    return results;
}

// New search_one_strand: takes pre-resolved threshold and effective_min_score.
// High-freq k-mers have already been removed from query_kmers by the preprocessor.
template <typename KmerInt>
static std::vector<ChainResult>
search_one_strand_preprocessed(
    const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
    int k,
    bool is_reverse,
    const KixReader& kix,
    const KpxReader& kpx,
    const OidFilter& filter,
    const SearchConfig& config,
    uint32_t resolved_threshold,
    uint32_t effective_min_score,
    Stage1Buffer* buf) {

    if (resolved_threshold == 0) return {};  // threshold was non-positive

    // Stage 1: candidate selection with pre-resolved threshold
    Stage1Config stage1_config = config.stage1;
    stage1_config.min_stage1_score = resolved_threshold;

    auto candidates = stage1_filter(query_kmers, kix, filter, stage1_config, buf);
    if (candidates.empty()) return {};

    // Mode 1: Stage 1 only — return candidates directly
    if (config.mode == 1) {
        return stage1_only_results(candidates, is_reverse, effective_min_score);
    }

    // Build candidate set for fast lookup and score map
    std::unordered_set<SeqId> candidate_set;
    std::unordered_map<SeqId, uint32_t> stage1_scores;
    candidate_set.reserve(candidates.size());
    stage1_scores.reserve(candidates.size());
    for (const auto& c : candidates) {
        candidate_set.insert(c.id);
        stage1_scores[c.id] = c.score;
    }

    // Stage 2: collect hits for candidates
    // High-freq k-mers already removed from query_kmers, so skip only cnt==0
    std::unordered_map<SeqId, std::vector<Hit>> hits_per_seq;

    const uint64_t* offsets = kix.offsets();
    const uint32_t* counts = kix.counts();
    const uint8_t* id_data = kix.posting_data();
    const uint64_t* pos_offsets = kpx.pos_offsets();
    const uint8_t* pos_data = kpx.posting_data();

    for (const auto& [q_pos, kmer] : query_kmers) {
        uint64_t kmer_idx = static_cast<uint64_t>(kmer);
        uint32_t cnt = counts[kmer_idx];
        if (cnt == 0) continue;

        SeqIdDecoder id_decoder(id_data + offsets[kmer_idx]);
        PosDecoder pos_decoder(pos_data + pos_offsets[kmer_idx]);

        for (uint32_t i = 0; i < cnt; i++) {
            SeqId sid = id_decoder.next();
            uint32_t s_pos = pos_decoder.next(id_decoder.was_new_seq());

            if (candidate_set.count(sid)) {
                hits_per_seq[sid].push_back({q_pos, s_pos});
            }
        }
    }

    // Chain hits for each candidate, using effective_min_score
    Stage2Config stage2_config = config.stage2;
    stage2_config.min_score = effective_min_score;

    std::vector<ChainResult> results;
    for (const auto& c : candidates) {
        auto it = hits_per_seq.find(c.id);
        if (it == hits_per_seq.end()) continue;

        ChainResult cr = chain_hits(it->second, c.id, k, is_reverse, stage2_config);
        if (cr.chainscore >= effective_min_score) {
            cr.stage1_score = c.score;
            results.push_back(cr);
        }
    }

    return results;
}

// Sort and truncate helper.
static void sort_and_truncate(SearchResult& result, const SearchConfig& config) {
    if (config.num_results > 0) {
        if (config.sort_score == 1) {
            std::sort(result.hits.begin(), result.hits.end(),
                      [](const ChainResult& a, const ChainResult& b) {
                          return a.stage1_score > b.stage1_score;
                      });
        } else {
            std::sort(result.hits.begin(), result.hits.end(),
                      [](const ChainResult& a, const ChainResult& b) {
                          return a.chainscore > b.chainscore;
                      });
        }

        if (result.hits.size() > config.num_results) {
            result.hits.resize(config.num_results);
        }
    }
}

// Search a single volume using pre-processed QueryKmerData with globally resolved thresholds.
template <typename KmerInt>
SearchResult search_volume(
    const std::string& query_id,
    const QueryKmerData<KmerInt>& qdata,
    int k,
    const KixReader& kix,
    const KpxReader& kpx,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config,
    Stage1Buffer* buf) {

    SearchResult result;
    result.query_id = query_id;

    // Search forward strand
    if (config.strand == 2 || config.strand == 1) {
        auto fwd_results = search_one_strand_preprocessed(
            qdata.fwd_kmers, k, false, kix, kpx, filter, config,
            qdata.resolved_threshold_fwd, qdata.effective_min_score_fwd, buf);
        result.hits.insert(result.hits.end(), fwd_results.begin(), fwd_results.end());
    }

    // Search reverse complement
    if (config.strand == 2 || config.strand == -1) {
        auto rc_results = search_one_strand_preprocessed(
            qdata.rc_kmers, k, true, kix, kpx, filter, config,
            qdata.resolved_threshold_rc, qdata.effective_min_score_rc, buf);
        result.hits.insert(result.hits.end(), rc_results.begin(), rc_results.end());
    }

    sort_and_truncate(result, config);
    return result;
}

// Explicit template instantiations
template SearchResult search_volume<uint16_t>(
    const std::string&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer*);
template SearchResult search_volume<uint32_t>(
    const std::string&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer*);

} // namespace ikafssn
