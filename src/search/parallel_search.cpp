#include "search/parallel_search.hpp"

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
#include "index/pfd_codec.hpp"
#include "core/spaced_seed.hpp"

#include <algorithm>
#include <climits>
#include <cstdint>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>

namespace ikafssn {

// Per-thread scratch for the .kpx candidate-set decoder.  The buffers grow
// monotonically across calls so the search hot path avoids per-k-mer
// allocation.  Each TBB worker holds its own instance.
static inline pfd::PosDecodeScratch& tls_pos_scratch() {
    thread_local pfd::PosDecodeScratch scratch;
    return scratch;
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

// Collect position hits for Stage 2A from one index set.
template <typename KmerInt>
static void collect_position_hits(
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    const KixReader& kix, const KpxReader& kpx,
    const std::vector<SeqId>& candidate_sids,
    std::unordered_map<SeqId, std::vector<Hit>>& hits_per_seq) {

    const uint8_t* kix_data = kix.posting_file();
    const uint8_t* pos_data = kpx.posting_file();
    pfd::PosDecodeScratch& scratch = tls_pos_scratch();

    for (size_t qi = 0; qi < n_kmers; qi++) {
        uint32_t q_pos = positions[qi];
        auto kmer_idx = kmers[qi];

        uint64_t kix_off, kix_len;
        kix.posting_list_range(kmer_idx, kix_off, kix_len);
        if (kix_len == 0) continue;

        SeqIdDecoder kix_dec(kix_data + kix_off,
                             kix_data + kix_off + kix_len);
        kix_dec.ensure_decoded();

        PosDecoder pos_decoder(pos_data + kpx.pos_offset(kmer_idx),
                               pos_data + kpx.posting_file_size(),
                               kix_dec.decoded_data(),
                               kix_dec.decoded_count(),
                               candidate_sids.data(),
                               candidate_sids.size(),
                               &scratch);
        pos_decoder.for_each_candidate(
            [&](uint32_t sid, const std::vector<uint32_t>& positions_v) {
                auto& bucket = hits_per_seq[sid];
                bucket.reserve(bucket.size() + positions_v.size());
                for (uint32_t spos : positions_v) {
                    bucket.push_back({q_pos, spos});
                }
            });
    }
}

template <typename KmerInt>
void phase_a_one_strand_single(
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    int k,
    bool is_reverse,
    const KixReader& kix,
    const KpxReader& kpx,
    const OidFilter& filter,
    const SearchConfig& config,
    uint32_t resolved_threshold,
    uint32_t effective_min_score,
    Stage1Buffer& buf,
    PhaseAState& state) {

    state.is_reverse = is_reverse;
    state.both_mode = false;
    state.effective_min_score = effective_min_score;
    state.span = seed_span(config.t, k);
    state.stage2_config = config.stage2;
    state.stage2_config.min_score = effective_min_score;

    if (resolved_threshold == 0 || n_kmers == 0) return;

    Stage1Config stage1_config = config.stage1;
    stage1_config.min_stage1_score = resolved_threshold;

    state.candidates =
        stage1_filter(positions, kmers, n_kmers, kix, filter, stage1_config, buf);
    if (state.candidates.empty()) return;

    if (config.mode == 1) {
        state.mode1_only = true;
        state.mode1_results =
            stage1_only_results(state.candidates, is_reverse, effective_min_score);
        return;
    }

    std::vector<SeqId> candidate_sids;
    candidate_sids.reserve(state.candidates.size());
    state.stage1_scores.reserve(state.candidates.size());
    for (const auto& c : state.candidates) {
        candidate_sids.push_back(c.id);
        state.stage1_scores[c.id] = c.score;
    }
    std::sort(candidate_sids.begin(), candidate_sids.end());

    collect_position_hits(positions, kmers, n_kmers, kix, kpx,
                          candidate_sids, state.hits_per_seq);
}

template <typename KmerInt>
void phase_a_one_strand_both(
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
    Stage1Buffer& buf,
    PhaseAState& state) {

    state.is_reverse = is_reverse;
    state.both_mode = true;
    state.effective_min_score = effective_min_score;
    state.span = seed_span(config.t, k);
    state.stage2_config = config.stage2;
    state.stage2_config.min_score = effective_min_score;

    if (n_cod == 0 && n_opt == 0) return;

    uint32_t unified_threshold = 0;
    if (resolved_threshold_cod > 0 && resolved_threshold_opt > 0) {
        unified_threshold = std::min(resolved_threshold_cod, resolved_threshold_opt);
    } else {
        unified_threshold = std::max(resolved_threshold_cod, resolved_threshold_opt);
    }
    if (unified_threshold == 0) return;

    Stage1Config s1cfg_acc = config.stage1;
    s1cfg_acc.min_stage1_score = 1;
    s1cfg_acc.stage1_topn = 0;

    {
        size_t ic = 0, io = 0;
        while (ic < n_cod || io < n_opt) {
            uint32_t p_c = (ic < n_cod) ? pos_cod[ic] : UINT32_MAX;
            uint32_t p_o = (io < n_opt) ? pos_opt[io] : UINT32_MAX;
            uint32_t cur = (p_c < p_o) ? p_c : p_o;

            if (p_c == cur) {
                size_t end = ic;
                while (end < n_cod && pos_cod[end] == cur) end++;
                stage1_filter_accumulate(pos_cod + ic, kmers_cod + ic,
                                         end - ic, kix_cod, filter,
                                         s1cfg_acc, buf);
                ic = end;
            }
            if (p_o == cur) {
                size_t end = io;
                while (end < n_opt && pos_opt[end] == cur) end++;
                stage1_filter_accumulate(pos_opt + io, kmers_opt + io,
                                         end - io, kix_opt, filter,
                                         s1cfg_acc, buf);
                io = end;
            }
        }
    }

    Stage1Config s1cfg_fin = config.stage1;
    s1cfg_fin.min_stage1_score = unified_threshold;
    state.candidates = stage1_filter_finish(buf, s1cfg_fin);

    if (state.candidates.empty()) return;

    if (config.mode == 1) {
        state.mode1_only = true;
        state.mode1_results =
            stage1_only_results(state.candidates, is_reverse, effective_min_score);
        return;
    }

    state.sorted_candidate_sids.reserve(state.candidates.size());
    state.stage1_scores.reserve(state.candidates.size());
    for (const auto& c : state.candidates) {
        state.sorted_candidate_sids.push_back(c.id);
        state.stage1_scores[c.id] = c.score;
    }
    std::sort(state.sorted_candidate_sids.begin(),
              state.sorted_candidate_sids.end());

    collect_position_hits(pos_cod, kmers_cod, n_cod, kix_cod, kpx_cod,
                          state.sorted_candidate_sids, state.hits_per_seq);
    collect_position_hits(pos_opt, kmers_opt, n_opt, kix_opt, kpx_opt,
                          state.sorted_candidate_sids, state.hits_per_seq);
}

std::vector<ChainResult>
phase_b_one_subject(SeqId sid, uint32_t stage1_score, const PhaseAState& state) {
    auto it = state.hits_per_seq.find(sid);
    if (it == state.hits_per_seq.end()) return {};
    auto chains = chain_hits(it->second, sid, state.span,
                             state.is_reverse, state.stage2_config);
    for (auto& cr : chains) cr.stage1_score = stage1_score;
    return chains;
}

void sort_and_truncate(SearchResult& result, const SearchConfig& config) {
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

// One element of the flat (job, strand, sid) tuple stream consumed by
// Phase B's parallel_for.  Stage1 score is captured at flatten time so
// Phase B does not need to look up `state.stage1_scores` again.
namespace {
struct FlatBSlot {
    size_t job_idx;
    uint8_t strand_idx;  // 0 = fwd, 1 = rc
    SeqId sid;
    uint32_t stage1_score;
};
}  // namespace

template <typename KmerInt>
std::vector<OrchestratorHit> run_search_jobs(
    const std::vector<QueryBundle<KmerInt>>& queries,
    const std::vector<VolumeBundle<KmerInt>>& volumes,
    const std::vector<std::pair<size_t, size_t>>& jobs,
    int k,
    const SearchConfig& config,
    bool both_mode,
    tbb::task_arena& arena,
    tbb::enumerable_thread_specific<Stage1Buffer>& tls_bufs) {

    // Per-job, per-strand Phase A state.  Indexed [job_idx][strand_idx].
    struct JobStates { PhaseAState fwd; PhaseAState rc; };
    std::vector<JobStates> job_states(jobs.size());

    if (jobs.empty()) return {};

    const bool run_fwd = (config.strand == 2 || config.strand == 1);
    const bool run_rc  = (config.strand == 2 || config.strand == -1);

    auto unify_min = [](uint32_t a, uint32_t b) -> uint32_t {
        if (a > 0 && b > 0) return std::min(a, b);
        return std::max(a, b);
    };

    // Phase A: one task per job.
    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, jobs.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                auto& buf = tls_bufs.local();
                for (size_t ji = range.begin(); ji != range.end(); ++ji) {
                    const auto [qi, vi] = jobs[ji];
                    const auto& qb = queries[qi];
                    const auto& vb = volumes[vi];
                    auto& js = job_states[ji];

                    if (both_mode) {
                        const auto& qd_cod = *qb.qdata_primary;
                        const auto& qd_opt = *qb.qdata_secondary;
                        if (run_fwd) {
                            phase_a_one_strand_both<KmerInt>(
                                qd_cod.fwd_positions.data(),
                                qd_cod.fwd_kmer_values.data(),
                                qd_cod.fwd_positions.size(),
                                qd_opt.fwd_positions.data(),
                                qd_opt.fwd_kmer_values.data(),
                                qd_opt.fwd_positions.size(),
                                k, false,
                                *vb.kix, *vb.kpx, *vb.kix_opt, *vb.kpx_opt,
                                *vb.filter, config,
                                qd_cod.resolved_threshold_fwd,
                                qd_opt.resolved_threshold_fwd,
                                unify_min(qd_cod.effective_min_score_fwd,
                                          qd_opt.effective_min_score_fwd),
                                buf, js.fwd);
                        }
                        if (run_rc) {
                            phase_a_one_strand_both<KmerInt>(
                                qd_cod.rc_positions.data(),
                                qd_cod.rc_kmer_values.data(),
                                qd_cod.rc_positions.size(),
                                qd_opt.rc_positions.data(),
                                qd_opt.rc_kmer_values.data(),
                                qd_opt.rc_positions.size(),
                                k, true,
                                *vb.kix, *vb.kpx, *vb.kix_opt, *vb.kpx_opt,
                                *vb.filter, config,
                                qd_cod.resolved_threshold_rc,
                                qd_opt.resolved_threshold_rc,
                                unify_min(qd_cod.effective_min_score_rc,
                                          qd_opt.effective_min_score_rc),
                                buf, js.rc);
                        }
                    } else {
                        const auto& qd = *qb.qdata_primary;
                        if (run_fwd) {
                            phase_a_one_strand_single<KmerInt>(
                                qd.fwd_positions.data(),
                                qd.fwd_kmer_values.data(),
                                qd.fwd_positions.size(),
                                k, false,
                                *vb.kix, *vb.kpx, *vb.filter, config,
                                qd.resolved_threshold_fwd,
                                qd.effective_min_score_fwd,
                                buf, js.fwd);
                        }
                        if (run_rc) {
                            phase_a_one_strand_single<KmerInt>(
                                qd.rc_positions.data(),
                                qd.rc_kmer_values.data(),
                                qd.rc_positions.size(),
                                k, true,
                                *vb.kix, *vb.kpx, *vb.filter, config,
                                qd.resolved_threshold_rc,
                                qd.effective_min_score_rc,
                                buf, js.rc);
                        }
                    }
                }
            });
    });

    // Mode 1 fast path: chain_hits is skipped, so collect mode1_results
    // directly without building the (job, strand, sid) flat list.
    if (config.mode == 1) {
        std::vector<OrchestratorHit> results;
        for (size_t ji = 0; ji < jobs.size(); ++ji) {
            const auto [qi, vi] = jobs[ji];
            const auto& vb = volumes[vi];
            auto& js = job_states[ji];
            for (auto* st : { &js.fwd, &js.rc }) {
                for (auto& cr : st->mode1_results) {
                    OrchestratorHit oh;
                    oh.query_idx = qi;
                    oh.volume_idx = vi;
                    oh.volume_index = vb.volume_index;
                    oh.cr = std::move(cr);
                    results.push_back(std::move(oh));
                }
                st->mode1_results.clear();
            }
        }
        return results;
    }

    // Build flat (job, strand, sid) work list for Phase B.  Iteration order
    // mirrors the wrapper helpers in volume_searcher.cpp: single-template
    // walks state.candidates in stage1_filter output order; both-mode walks
    // sorted_candidate_sids in ascending seq_id order.
    std::vector<FlatBSlot> flat;
    {
        size_t reserve = 0;
        for (auto& js : job_states) {
            reserve += js.fwd.candidates.size() + js.rc.candidates.size();
        }
        flat.reserve(reserve);
    }
    for (size_t ji = 0; ji < jobs.size(); ++ji) {
        auto& js = job_states[ji];
        PhaseAState* arms[2] = { &js.fwd, &js.rc };
        for (uint8_t s = 0; s < 2; ++s) {
            const PhaseAState* st = arms[s];
            if (st->mode1_only) continue;  // mode 1 already handled
            if (st->candidates.empty()) continue;
            if (st->both_mode) {
                for (SeqId sid : st->sorted_candidate_sids) {
                    auto sit = st->stage1_scores.find(sid);
                    uint32_t sc = (sit != st->stage1_scores.end()) ? sit->second : 0;
                    flat.push_back({ ji, s, sid, sc });
                }
            } else {
                for (const auto& c : st->candidates) {
                    flat.push_back({ ji, s, c.id, c.score });
                }
            }
        }
    }

    if (flat.empty()) return {};

    // Phase B: flat parallel_for over (job, strand, sid) tuples.  TBB
    // work-stealing handles fairness across queries with very different
    // candidate counts.
    tbb::combinable<std::vector<OrchestratorHit>> tls_results;
    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, flat.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                auto& local = tls_results.local();
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    const auto& slot = flat[i];
                    const auto [qi, vi] = jobs[slot.job_idx];
                    const auto& vb = volumes[vi];
                    const auto& js = job_states[slot.job_idx];
                    const PhaseAState& st = (slot.strand_idx == 0) ? js.fwd : js.rc;
                    auto chains = phase_b_one_subject(slot.sid, slot.stage1_score, st);
                    for (auto& cr : chains) {
                        OrchestratorHit oh;
                        oh.query_idx = qi;
                        oh.volume_idx = vi;
                        oh.volume_index = vb.volume_index;
                        oh.cr = std::move(cr);
                        local.push_back(std::move(oh));
                    }
                }
            });
    });

    std::vector<OrchestratorHit> results;
    tls_results.combine_each([&results](std::vector<OrchestratorHit>& local) {
        results.insert(results.end(),
                       std::make_move_iterator(local.begin()),
                       std::make_move_iterator(local.end()));
    });
    return results;
}

// Explicit template instantiations.
template void phase_a_one_strand_single<uint16_t>(
    const uint32_t*, const uint16_t*, size_t, int, bool,
    const KixReader&, const KpxReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);
template void phase_a_one_strand_single<uint32_t>(
    const uint32_t*, const uint32_t*, size_t, int, bool,
    const KixReader&, const KpxReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);

template void phase_a_one_strand_both<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const uint32_t*, const uint16_t*, size_t,
    int, bool,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);
template void phase_a_one_strand_both<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const uint32_t*, const uint32_t*, size_t,
    int, bool,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);

template std::vector<OrchestratorHit> run_search_jobs<uint16_t>(
    const std::vector<QueryBundle<uint16_t>>&,
    const std::vector<VolumeBundle<uint16_t>>&,
    const std::vector<std::pair<size_t, size_t>>&,
    int, const SearchConfig&, bool,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);
template std::vector<OrchestratorHit> run_search_jobs<uint32_t>(
    const std::vector<QueryBundle<uint32_t>>&,
    const std::vector<VolumeBundle<uint32_t>>&,
    const std::vector<std::pair<size_t, size_t>>&,
    int, const SearchConfig&, bool,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);

} // namespace ikafssn
