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
    scanner.scan(seq.data(), seq.size(), [&](uint32_t pos, KmerInt kmer) {
        kmers.emplace_back(pos, kmer);
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
        cr.score = 0;
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
    uint32_t effective_min_score) {

    if (resolved_threshold == 0) return {};  // threshold was non-positive

    // Stage 1: candidate selection with pre-resolved threshold
    Stage1Config stage1_config = config.stage1;
    stage1_config.min_stage1_score = resolved_threshold;

    auto candidates = stage1_filter(query_kmers, kix, filter, stage1_config);
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
        if (cr.score >= effective_min_score) {
            cr.stage1_score = c.score;
            results.push_back(cr);
        }
    }

    return results;
}

// Legacy search_one_strand: per-volume threshold resolution (backward-compatible).
template <typename KmerInt>
static std::vector<ChainResult>
search_one_strand(const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
                  int k,
                  bool is_reverse,
                  const KixReader& kix,
                  const KpxReader& kpx,
                  const OidFilter& filter,
                  const SearchConfig& config,
                  const KhxReader* khx) {

    // Resolve fractional min_stage1_score if active
    Stage1Config stage1_config = config.stage1;
    if (config.min_stage1_score_frac > 0) {
        const bool use_coverscore = (stage1_config.stage1_score_type == 1);

        // Compute Nqkmer: number of query k-mers
        uint32_t Nqkmer;
        if (use_coverscore) {
            // Count distinct k-mer values in query
            std::unordered_set<uint64_t> distinct;
            for (const auto& [pos, kmer] : query_kmers) {
                distinct.insert(static_cast<uint64_t>(kmer));
            }
            Nqkmer = static_cast<uint32_t>(distinct.size());
        } else {
            // matchscore: total k-mer positions
            Nqkmer = static_cast<uint32_t>(query_kmers.size());
        }

        // Compute effective max_freq
        uint32_t max_freq = compute_effective_max_freq(
            stage1_config.max_freq, kix.total_postings(), kix.table_size());

        // Count Nhighfreq from query k-mers
        const uint32_t* counts = kix.counts();
        uint32_t Nhighfreq;
        if (use_coverscore) {
            // Count distinct high-freq k-mer values
            std::unordered_set<uint64_t> highfreq_distinct;
            for (const auto& [pos, kmer] : query_kmers) {
                uint64_t kmer_idx = static_cast<uint64_t>(kmer);
                bool is_highfreq = (counts[kmer_idx] > max_freq) ||
                    (khx != nullptr && khx->is_excluded(kmer_idx));
                if (is_highfreq) {
                    highfreq_distinct.insert(kmer_idx);
                }
            }
            Nhighfreq = static_cast<uint32_t>(highfreq_distinct.size());
        } else {
            // Count all occurrences of high-freq k-mers
            Nhighfreq = 0;
            for (const auto& [pos, kmer] : query_kmers) {
                uint64_t kmer_idx = static_cast<uint64_t>(kmer);
                bool is_highfreq = (counts[kmer_idx] > max_freq) ||
                    (khx != nullptr && khx->is_excluded(kmer_idx));
                if (is_highfreq) {
                    Nhighfreq++;
                }
            }
        }

        // Compute threshold
        int32_t threshold = static_cast<int32_t>(
            std::ceil(static_cast<double>(Nqkmer) * config.min_stage1_score_frac))
            - static_cast<int32_t>(Nhighfreq);

        if (threshold <= 0) {
            // Threshold is non-positive: skip this strand
            std::fprintf(stderr,
                "Warning: fractional min_stage1_score threshold <= 0 "
                "(strand=%s, Nqkmer=%u, Nhighfreq=%u, P=%.4f)\n",
                is_reverse ? "rc" : "fwd", Nqkmer, Nhighfreq,
                config.min_stage1_score_frac);
            return {};
        }

        stage1_config.min_stage1_score = static_cast<uint32_t>(threshold);
    }

    // Stage 1: candidate selection
    auto candidates = stage1_filter(query_kmers, kix, filter, stage1_config);
    if (candidates.empty()) return {};

    // Mode 1: Stage 1 only — return candidates directly
    if (config.mode == 1) {
        return stage1_only_results(candidates, is_reverse,
                                   config.stage2.min_score);
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

    // Compute effective max_freq (same as Stage 1)
    uint32_t max_freq = compute_effective_max_freq(
        config.stage1.max_freq, kix.total_postings(), kix.table_size());

    // Stage 2: collect hits for candidates
    std::unordered_map<SeqId, std::vector<Hit>> hits_per_seq;

    const uint64_t* offsets = kix.offsets();
    const uint32_t* counts = kix.counts();
    const uint8_t* id_data = kix.posting_data();
    const uint64_t* pos_offsets = kpx.pos_offsets();
    const uint8_t* pos_data = kpx.posting_data();

    for (const auto& [q_pos, kmer] : query_kmers) {
        uint64_t kmer_idx = static_cast<uint64_t>(kmer);
        uint32_t cnt = counts[kmer_idx];
        if (cnt == 0 || cnt > max_freq) continue;

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

    // Chain hits for each candidate
    std::vector<ChainResult> results;
    for (const auto& c : candidates) {
        auto it = hits_per_seq.find(c.id);
        if (it == hits_per_seq.end()) continue;

        ChainResult cr = chain_hits(it->second, c.id, k, is_reverse, config.stage2);
        if (cr.score >= config.stage2.min_score) {
            cr.stage1_score = c.score;
            results.push_back(cr);
        }
    }

    return results;
}

// Sort and truncate helper (shared by both search_volume overloads).
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
                          return a.score > b.score;
                      });
        }

        if (result.hits.size() > config.num_results) {
            result.hits.resize(config.num_results);
        }
    }
}

// Legacy search_volume: backward-compatible wrapper.
// Builds single-element all_kix, calls preprocess_query, delegates to new overload.
template <typename KmerInt>
SearchResult search_volume(
    const std::string& query_id,
    const std::string& query_seq,
    int k,
    const KixReader& kix,
    const KpxReader& kpx,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config,
    const KhxReader* khx) {

    SearchResult result;
    result.query_id = query_id;

    // Extract k-mers from forward strand
    auto fwd_kmers = extract_kmers<KmerInt>(query_seq, k);

    // Search forward strand
    if (config.strand == 2 || config.strand == 1) {
        auto fwd_results = search_one_strand(fwd_kmers, k, false, kix, kpx, filter, config, khx);
        result.hits.insert(result.hits.end(), fwd_results.begin(), fwd_results.end());
    }

    // Search reverse complement
    if (config.strand == 2 || config.strand == -1) {
        std::vector<std::pair<uint32_t, KmerInt>> rc_kmers;
        rc_kmers.reserve(fwd_kmers.size());
        for (const auto& [pos, kmer] : fwd_kmers) {
            KmerInt rc = kmer_revcomp(kmer, k);
            rc_kmers.emplace_back(pos, rc);
        }

        auto rc_results = search_one_strand(rc_kmers, k, true, kix, kpx, filter, config, khx);
        result.hits.insert(result.hits.end(), rc_results.begin(), rc_results.end());
    }

    sort_and_truncate(result, config);
    return result;
}

// New search_volume: uses pre-processed QueryKmerData with globally resolved thresholds.
template <typename KmerInt>
SearchResult search_volume(
    const std::string& query_id,
    const QueryKmerData<KmerInt>& qdata,
    int k,
    const KixReader& kix,
    const KpxReader& kpx,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config) {

    SearchResult result;
    result.query_id = query_id;

    // Search forward strand
    if (config.strand == 2 || config.strand == 1) {
        auto fwd_results = search_one_strand_preprocessed(
            qdata.fwd_kmers, k, false, kix, kpx, filter, config,
            qdata.resolved_threshold_fwd, qdata.effective_min_score_fwd);
        result.hits.insert(result.hits.end(), fwd_results.begin(), fwd_results.end());
    }

    // Search reverse complement
    if (config.strand == 2 || config.strand == -1) {
        auto rc_results = search_one_strand_preprocessed(
            qdata.rc_kmers, k, true, kix, kpx, filter, config,
            qdata.resolved_threshold_rc, qdata.effective_min_score_rc);
        result.hits.insert(result.hits.end(), rc_results.begin(), rc_results.end());
    }

    sort_and_truncate(result, config);
    return result;
}

// Explicit template instantiations for legacy overload
template SearchResult search_volume<uint16_t>(
    const std::string&, const std::string&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, const KhxReader*);
template SearchResult search_volume<uint32_t>(
    const std::string&, const std::string&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, const KhxReader*);

// Explicit template instantiations for new overload
template SearchResult search_volume<uint16_t>(
    const std::string&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&);
template SearchResult search_volume<uint32_t>(
    const std::string&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&);

} // namespace ikafssn
