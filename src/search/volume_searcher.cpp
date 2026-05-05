#include "search/volume_searcher.hpp"
#include "search/parallel_search.hpp"
#include "search/oid_filter.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"
#include "search/query_preprocessor.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "core/spaced_seed.hpp"

#include <algorithm>
#include <iterator>

namespace ikafssn {

// Drain a JobState by running Stage 2B for every candidate in the order
// produced by Stage 1.  Used by the volume-level wrappers below; the
// parallel orchestrator drives Stage 2B as one global parallel_for over
// (ext_job, sid) pairs.
static std::vector<ChainResult>
drain_stage2b_single_template(const JobState& state) {
    if (state.mode1_only) return state.mode1_results;
    std::vector<ChainResult> results;
    for (const auto& c : state.candidates) {
        auto chains = stage2b_one_subject(c.id, c.score, state);
        if (!chains.empty()) {
            results.insert(results.end(),
                           std::make_move_iterator(chains.begin()),
                           std::make_move_iterator(chains.end()));
        }
    }
    return results;
}

static std::vector<ChainResult>
drain_stage2b_both_template(const JobState& state) {
    if (state.mode1_only) return state.mode1_results;
    std::vector<ChainResult> results;
    for (SeqId sid : state.sorted_candidate_sids) {
        auto sit = state.stage1_scores.find(sid);
        uint32_t score = (sit != state.stage1_scores.end()) ? sit->second : 0;
        auto chains = stage2b_one_subject(sid, score, state);
        if (!chains.empty()) {
            results.insert(results.end(),
                           std::make_move_iterator(chains.begin()),
                           std::make_move_iterator(chains.end()));
        }
    }
    return results;
}

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
    Stage1Buffer& buf) {

    SearchResult result;
    result.query_id = query_id;

    if (config.strand == 2 || config.strand == 1) {
        JobState state;
        stage1_one_strand_single(
            qdata.fwd_positions.data(), qdata.fwd_kmer_values.data(),
            qdata.fwd_positions.size(),
            k, false, kix, filter, config,
            qdata.resolved_threshold_fwd, qdata.effective_min_score_fwd,
            buf, state);
        stage2a_one_strand_single(
            qdata.fwd_positions.data(), qdata.fwd_kmer_values.data(),
            qdata.fwd_positions.size(),
            kix, kpx, state);
        auto fwd_results = drain_stage2b_single_template(state);
        result.hits.insert(result.hits.end(),
                           std::make_move_iterator(fwd_results.begin()),
                           std::make_move_iterator(fwd_results.end()));
    }

    if (config.strand == 2 || config.strand == -1) {
        JobState state;
        stage1_one_strand_single(
            qdata.rc_positions.data(), qdata.rc_kmer_values.data(),
            qdata.rc_positions.size(),
            k, true, kix, filter, config,
            qdata.resolved_threshold_rc, qdata.effective_min_score_rc,
            buf, state);
        stage2a_one_strand_single(
            qdata.rc_positions.data(), qdata.rc_kmer_values.data(),
            qdata.rc_positions.size(),
            kix, kpx, state);
        auto rc_results = drain_stage2b_single_template(state);
        result.hits.insert(result.hits.end(),
                           std::make_move_iterator(rc_results.begin()),
                           std::make_move_iterator(rc_results.end()));
    }

    sort_and_truncate(result, config);
    return result;
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
    Stage1Buffer& buf) {

    SearchResult result;
    result.query_id = query_id;

    auto unify_min = [](uint32_t a, uint32_t b) -> uint32_t {
        if (a > 0 && b > 0) return std::min(a, b);
        return std::max(a, b);
    };

    if (config.strand == 2 || config.strand == 1) {
        JobState state;
        stage1_one_strand_both(
            qdata_cod.fwd_positions.data(), qdata_cod.fwd_kmer_values.data(),
            qdata_cod.fwd_positions.size(),
            qdata_opt.fwd_positions.data(), qdata_opt.fwd_kmer_values.data(),
            qdata_opt.fwd_positions.size(),
            k, false,
            kix_cod, kix_opt,
            filter, config,
            qdata_cod.resolved_threshold_fwd, qdata_opt.resolved_threshold_fwd,
            unify_min(qdata_cod.effective_min_score_fwd,
                      qdata_opt.effective_min_score_fwd),
            buf, state);
        stage2a_one_strand_both(
            qdata_cod.fwd_positions.data(), qdata_cod.fwd_kmer_values.data(),
            qdata_cod.fwd_positions.size(),
            qdata_opt.fwd_positions.data(), qdata_opt.fwd_kmer_values.data(),
            qdata_opt.fwd_positions.size(),
            kix_cod, kpx_cod, kix_opt, kpx_opt, state);
        auto fwd_results = drain_stage2b_both_template(state);
        result.hits.insert(result.hits.end(),
                           std::make_move_iterator(fwd_results.begin()),
                           std::make_move_iterator(fwd_results.end()));
    }

    if (config.strand == 2 || config.strand == -1) {
        JobState state;
        stage1_one_strand_both(
            qdata_cod.rc_positions.data(), qdata_cod.rc_kmer_values.data(),
            qdata_cod.rc_positions.size(),
            qdata_opt.rc_positions.data(), qdata_opt.rc_kmer_values.data(),
            qdata_opt.rc_positions.size(),
            k, true,
            kix_cod, kix_opt,
            filter, config,
            qdata_cod.resolved_threshold_rc, qdata_opt.resolved_threshold_rc,
            unify_min(qdata_cod.effective_min_score_rc,
                      qdata_opt.effective_min_score_rc),
            buf, state);
        stage2a_one_strand_both(
            qdata_cod.rc_positions.data(), qdata_cod.rc_kmer_values.data(),
            qdata_cod.rc_positions.size(),
            qdata_opt.rc_positions.data(), qdata_opt.rc_kmer_values.data(),
            qdata_opt.rc_positions.size(),
            kix_cod, kpx_cod, kix_opt, kpx_opt, state);
        auto rc_results = drain_stage2b_both_template(state);
        result.hits.insert(result.hits.end(),
                           std::make_move_iterator(rc_results.begin()),
                           std::make_move_iterator(rc_results.end()));
    }

    sort_and_truncate(result, config);
    return result;
}

template SearchResult search_volume<uint16_t>(
    const std::string&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer&);
template SearchResult search_volume<uint32_t>(
    const std::string&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer&);

template SearchResult search_volume_both<uint16_t>(
    const std::string&,
    const QueryKmerData<uint16_t>&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer&);
template SearchResult search_volume_both<uint32_t>(
    const std::string&,
    const QueryKmerData<uint32_t>&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer&);

} // namespace ikafssn
