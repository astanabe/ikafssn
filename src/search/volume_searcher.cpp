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
#include "core/spaced_seed.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>

namespace ikafssn {

template <typename KmerInt>
static std::vector<std::pair<uint32_t, KmerInt>>
extract_kmers(const std::string& seq, int k, int max_expansion = 16) {
    std::vector<std::pair<uint32_t, KmerInt>> kmers;
    KmerScanner<KmerInt> scanner(k);
    scanner.scan_ambig(seq.data(), seq.size(),
        [&](uint32_t pos, KmerInt kmer) {
            kmers.emplace_back(pos, kmer);
        },
        [&](uint32_t pos, KmerInt base_kmer, const AmbigInfo* infos, int count) {
            expand_ambig_kmer_multi<KmerInt>(base_kmer, infos, count,
                [&](KmerInt expanded) {
                    kmers.emplace_back(pos, expanded);
                });
        },
        nullptr,
        max_expansion);
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
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    int k,
    bool is_reverse,
    const KixReader& kix,
    const KpxReader& kpx,
    const OidFilter& filter,
    const SearchConfig& config,
    uint32_t resolved_threshold,
    uint32_t effective_min_score,
    Stage1Buffer* buf) {

    if (resolved_threshold == 0 || n_kmers == 0) return {};

    // Stage 1: candidate selection with pre-resolved threshold
    Stage1Config stage1_config = config.stage1;
    stage1_config.min_stage1_score = resolved_threshold;

    auto candidates = stage1_filter(positions, kmers, n_kmers, kix, filter, stage1_config, buf);
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
    std::unordered_map<SeqId, std::vector<Hit>> hits_per_seq;

    const uint8_t* id_data = kix.posting_data();
    const uint8_t* pos_data = kpx.posting_data();

    for (size_t qi = 0; qi < n_kmers; qi++) {
        uint32_t q_pos = positions[qi];
        auto kmer_idx = kmers[qi];
        auto off = kix.posting_offset(kmer_idx);
        auto end_off = kix.posting_offset(kmer_idx + 1);
        if (off == end_off) continue;

        SeqIdDecoder id_decoder(id_data + off, id_data + end_off);
        PosDecoder pos_decoder(pos_data + kpx.pos_offset(kmer_idx));

        while (id_decoder.has_more()) {
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

        auto chains = chain_hits(it->second, c.id, seed_span(config.t, k), is_reverse, stage2_config);
        for (auto& cr : chains) {
            cr.stage1_score = c.score;
            results.push_back(cr);
        }
    }

    return results;
}

// Sort and truncate helper.
static void sort_and_truncate(SearchResult& result, const SearchConfig& config) {
    if (config.num_results > 0) {
        auto cmp = (config.sort_score == 1)
            ? [](const ChainResult& a, const ChainResult& b) {
                  return a.stage1_score > b.stage1_score;
              }
            : [](const ChainResult& a, const ChainResult& b) {
                  return a.chainscore > b.chainscore;
              };

        if (result.hits.size() > config.num_results) {
            std::nth_element(result.hits.begin(),
                             result.hits.begin() + config.num_results,
                             result.hits.end(), cmp);
            result.hits.resize(config.num_results);
        }
        std::sort(result.hits.begin(), result.hits.end(), cmp);
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
            qdata.fwd_positions.data(), qdata.fwd_kmer_values.data(),
            qdata.fwd_positions.size(),
            k, false, kix, kpx, filter, config,
            qdata.resolved_threshold_fwd, qdata.effective_min_score_fwd, buf);
        result.hits.insert(result.hits.end(), fwd_results.begin(), fwd_results.end());
    }

    // Search reverse complement
    if (config.strand == 2 || config.strand == -1) {
        auto rc_results = search_one_strand_preprocessed(
            qdata.rc_positions.data(), qdata.rc_kmer_values.data(),
            qdata.rc_positions.size(),
            k, true, kix, kpx, filter, config,
            qdata.resolved_threshold_rc, qdata.effective_min_score_rc, buf);
        result.hits.insert(result.hits.end(), rc_results.begin(), rc_results.end());
    }

    sort_and_truncate(result, config);
    return result;
}

// Collect position hits for Stage 2 from one index set.
template <typename KmerInt>
static void collect_position_hits(
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    const KixReader& kix, const KpxReader& kpx,
    const std::unordered_set<SeqId>& candidate_set,
    std::unordered_map<SeqId, std::vector<Hit>>& hits_per_seq) {

    const uint8_t* id_data = kix.posting_data();
    const uint8_t* pos_data = kpx.posting_data();

    for (size_t qi = 0; qi < n_kmers; qi++) {
        uint32_t q_pos = positions[qi];
        auto kmer_idx = kmers[qi];
        auto off = kix.posting_offset(kmer_idx);
        auto end_off = kix.posting_offset(kmer_idx + 1);
        if (off == end_off) continue;

        SeqIdDecoder id_decoder(id_data + off, id_data + end_off);
        PosDecoder pos_decoder(pos_data + kpx.pos_offset(kmer_idx));

        while (id_decoder.has_more()) {
            SeqId sid = id_decoder.next();
            uint32_t s_pos = pos_decoder.next(id_decoder.was_new_seq());

            if (candidate_set.count(sid)) {
                hits_per_seq[sid].push_back({q_pos, s_pos});
            }
        }
    }
}

// Search a single volume using merged coding+optimal ("both" mode).
template <typename KmerInt>
static std::vector<ChainResult>
search_one_strand_both(
    const uint32_t* pos_cod, const KmerInt* kmers_cod, size_t n_cod,
    const uint32_t* pos_opt, const KmerInt* kmers_opt, size_t n_opt,
    int k,
    bool is_reverse,
    const KixReader& kix_cod, const KpxReader& kpx_cod,
    const KixReader& kix_opt, const KpxReader& kpx_opt,
    const OidFilter& filter,
    const SearchConfig& config,
    uint32_t resolved_threshold_cod,
    uint32_t resolved_threshold_opt,
    uint32_t effective_min_score,
    Stage1Buffer* buf_cod,
    Stage1Buffer* buf_opt) {

    if (n_cod == 0 && n_opt == 0) return {};

    // Stage 1: run independently on coding and optimal
    Stage1Config s1cfg = config.stage1;
    s1cfg.min_stage1_score = 1;     // collect all, merge later
    s1cfg.stage1_topn = 0;          // no truncation per-template

    std::vector<Stage1Candidate> cand_cod, cand_opt;
    if (n_cod > 0) {
        cand_cod = stage1_filter(pos_cod, kmers_cod, n_cod, kix_cod, filter, s1cfg, buf_cod);
    }
    if (n_opt > 0) {
        cand_opt = stage1_filter(pos_opt, kmers_opt, n_opt, kix_opt, filter, s1cfg, buf_opt);
    }

    // Merge: sum scores per SeqId
    std::unordered_map<SeqId, uint32_t> merged_scores;
    merged_scores.reserve(cand_cod.size() + cand_opt.size());
    for (const auto& c : cand_cod) merged_scores[c.id] += c.score;
    for (const auto& c : cand_opt) merged_scores[c.id] += c.score;

    // Apply combined threshold
    uint32_t combined_threshold = resolved_threshold_cod + resolved_threshold_opt;
    if (combined_threshold == 0) return {};

    // Mode 1: Stage 1 only — return merged candidates directly
    if (config.mode == 1) {
        std::vector<ChainResult> results;
        for (const auto& [sid, score] : merged_scores) {
            if (score < combined_threshold) continue;
            ChainResult cr{};
            cr.seq_id = sid;
            cr.chainscore = 0;
            cr.stage1_score = score;
            cr.is_reverse = is_reverse;
            results.push_back(cr);
        }
        return results;
    }

    // Filter candidates
    std::unordered_set<SeqId> candidate_set;
    std::unordered_map<SeqId, uint32_t> stage1_scores;
    for (const auto& [sid, score] : merged_scores) {
        if (score >= combined_threshold) {
            candidate_set.insert(sid);
            stage1_scores[sid] = score;
        }
    }
    if (candidate_set.empty()) return {};

    // Apply stage1_topn if set
    if (config.stage1.stage1_topn > 0 && candidate_set.size() > config.stage1.stage1_topn) {
        std::vector<Stage1Candidate> sorted_cands;
        sorted_cands.reserve(candidate_set.size());
        for (auto sid : candidate_set) {
            sorted_cands.push_back({sid, stage1_scores[sid]});
        }
        auto cmp = [](const Stage1Candidate& a, const Stage1Candidate& b) {
            return a.score > b.score;
        };
        std::nth_element(sorted_cands.begin(),
                         sorted_cands.begin() + config.stage1.stage1_topn,
                         sorted_cands.end(), cmp);
        sorted_cands.resize(config.stage1.stage1_topn);
        candidate_set.clear();
        for (const auto& c : sorted_cands) candidate_set.insert(c.id);
    }

    // Stage 2: collect position hits from both indexes
    std::unordered_map<SeqId, std::vector<Hit>> hits_per_seq;

    collect_position_hits(pos_cod, kmers_cod, n_cod, kix_cod, kpx_cod,
                          candidate_set, hits_per_seq);
    collect_position_hits(pos_opt, kmers_opt, n_opt, kix_opt, kpx_opt,
                          candidate_set, hits_per_seq);

    // Chain hits
    Stage2Config stage2_config = config.stage2;
    stage2_config.min_score = effective_min_score;

    std::vector<ChainResult> results;
    for (auto sid : candidate_set) {
        auto it = hits_per_seq.find(sid);
        if (it == hits_per_seq.end()) continue;

        auto chains = chain_hits(it->second, sid, seed_span(config.t, k), is_reverse, stage2_config);
        for (auto& cr : chains) {
            cr.stage1_score = stage1_scores[sid];
            results.push_back(cr);
        }
    }

    return results;
}

template <typename KmerInt>
SearchResult search_volume_both(
    const std::string& query_id,
    const QueryKmerData<KmerInt>& qdata_cod,
    const QueryKmerData<KmerInt>& qdata_opt,
    int k,
    const KixReader& kix_cod, const KpxReader& kpx_cod,
    const KixReader& kix_opt, const KpxReader& kpx_opt,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config,
    Stage1Buffer* buf_cod,
    Stage1Buffer* buf_opt) {

    SearchResult result;
    result.query_id = query_id;

    // Search forward strand
    if (config.strand == 2 || config.strand == 1) {
        auto fwd_results = search_one_strand_both(
            qdata_cod.fwd_positions.data(), qdata_cod.fwd_kmer_values.data(),
            qdata_cod.fwd_positions.size(),
            qdata_opt.fwd_positions.data(), qdata_opt.fwd_kmer_values.data(),
            qdata_opt.fwd_positions.size(),
            k, false,
            kix_cod, kpx_cod, kix_opt, kpx_opt,
            filter, config,
            qdata_cod.resolved_threshold_fwd, qdata_opt.resolved_threshold_fwd,
            std::max(qdata_cod.effective_min_score_fwd, qdata_opt.effective_min_score_fwd),
            buf_cod, buf_opt);
        result.hits.insert(result.hits.end(), fwd_results.begin(), fwd_results.end());
    }

    // Search reverse complement
    if (config.strand == 2 || config.strand == -1) {
        auto rc_results = search_one_strand_both(
            qdata_cod.rc_positions.data(), qdata_cod.rc_kmer_values.data(),
            qdata_cod.rc_positions.size(),
            qdata_opt.rc_positions.data(), qdata_opt.rc_kmer_values.data(),
            qdata_opt.rc_positions.size(),
            k, true,
            kix_cod, kpx_cod, kix_opt, kpx_opt,
            filter, config,
            qdata_cod.resolved_threshold_rc, qdata_opt.resolved_threshold_rc,
            std::max(qdata_cod.effective_min_score_rc, qdata_opt.effective_min_score_rc),
            buf_cod, buf_opt);
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

template SearchResult search_volume_both<uint16_t>(
    const std::string&,
    const QueryKmerData<uint16_t>&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer*, Stage1Buffer*);
template SearchResult search_volume_both<uint32_t>(
    const std::string&,
    const QueryKmerData<uint32_t>&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer*, Stage1Buffer*);

} // namespace ikafssn
