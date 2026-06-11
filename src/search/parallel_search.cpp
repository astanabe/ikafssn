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

    // Hoisted out of the qi loop so the StreamCtx::decoded buffer is
    // reused across all per-thread probes.
    thread_local SeqIdDecoder kix_dec;

    for (size_t qi = 0; qi < n_kmers; qi++) {
        uint32_t q_pos = positions[qi];
        auto kmer_idx = kmers[qi];

        uint64_t kix_off, kix_len;
        kix.posting_list_range(kmer_idx, kix_off, kix_len);
        if (kix_len == 0) continue;

        kix_dec.reset(kix_data + kix_off,
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
void stage1_one_strand_single(
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    int k,
    bool is_reverse,
    const KixReader& kix,
    const OidFilter& filter,
    const SearchConfig& config,
    uint32_t resolved_threshold,
    uint32_t effective_min_score,
    Stage1Buffer& buf,
    JobState& state) {

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

    state.stage1_scores.reserve(state.candidates.size());
    for (const auto& c : state.candidates) {
        state.stage1_scores[c.id] = c.score;
    }
}

template <typename KmerInt>
void stage2a_one_strand_single(
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    const KixReader& kix, const KpxReader& kpx,
    JobState& state) {

    if (state.mode1_only) return;
    if (state.candidates.empty()) return;

    std::vector<SeqId> candidate_sids;
    candidate_sids.reserve(state.candidates.size());
    for (const auto& c : state.candidates) {
        candidate_sids.push_back(c.id);
    }
    std::sort(candidate_sids.begin(), candidate_sids.end());

    collect_position_hits(positions, kmers, n_kmers, kix, kpx,
                          candidate_sids, state.hits_per_seq);
}

template <typename KmerInt>
void stage1_one_strand_both(
    const uint32_t* pos_cod, const KmerInt* kmers_cod, size_t n_cod,
    const uint32_t* pos_opt, const KmerInt* kmers_opt, size_t n_opt,
    int k,
    bool is_reverse,
    const KixReader& kix_cod, const KixReader& kix_opt,
    const OidFilter& filter,
    const SearchConfig& config,
    uint32_t resolved_threshold_cod,
    uint32_t resolved_threshold_opt,
    uint32_t effective_min_score,
    Stage1Buffer& buf,
    JobState& state) {

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

    // Count distinct merged query positions (Nq) so the score cutoff can be
    // engaged.  A sequence's final score is at most the number of distinct
    // positions it matches, so once score + remaining < unified_threshold it can
    // never reach the threshold and its scatter is skipped.
    uint32_t Nq = 0;
    {
        size_t ic = 0, io = 0;
        while (ic < n_cod || io < n_opt) {
            uint32_t p_c = (ic < n_cod) ? pos_cod[ic] : UINT32_MAX;
            uint32_t p_o = (io < n_opt) ? pos_opt[io] : UINT32_MAX;
            uint32_t cur = (p_c < p_o) ? p_c : p_o;
            while (ic < n_cod && pos_cod[ic] == cur) ic++;
            while (io < n_opt && pos_opt[io] == cur) io++;
            Nq++;
        }
    }
    s1cfg_acc.cutoff_threshold = unified_threshold;

    {
        size_t ic = 0, io = 0;
        uint32_t g = 0;   // position groups consumed so far
        while (ic < n_cod || io < n_opt) {
            uint32_t p_c = (ic < n_cod) ? pos_cod[ic] : UINT32_MAX;
            uint32_t p_o = (io < n_opt) ? pos_opt[io] : UINT32_MAX;
            uint32_t cur = (p_c < p_o) ? p_c : p_o;

            // Positions still to consume, including the current one.
            s1cfg_acc.cutoff_remaining = Nq - g;

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
            g++;
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
}

template <typename KmerInt>
void stage2a_one_strand_both(
    const uint32_t* pos_cod, const KmerInt* kmers_cod, size_t n_cod,
    const uint32_t* pos_opt, const KmerInt* kmers_opt, size_t n_opt,
    const KixReader& kix_cod, const KpxReader& kpx_cod,
    const KixReader& kix_opt, const KpxReader& kpx_opt,
    JobState& state) {

    if (state.mode1_only) return;
    if (state.candidates.empty()) return;

    collect_position_hits(pos_cod, kmers_cod, n_cod, kix_cod, kpx_cod,
                          state.sorted_candidate_sids, state.hits_per_seq);
    collect_position_hits(pos_opt, kmers_opt, n_opt, kix_opt, kpx_opt,
                          state.sorted_candidate_sids, state.hits_per_seq);
}

std::vector<ChainResult>
stage2b_one_subject(SeqId sid, uint32_t stage1_score, const JobState& state) {
    auto it = state.hits_per_seq.find(sid);
    if (it == state.hits_per_seq.end()) return {};
    auto chains = chain_hits(it->second, sid, state.span,
                             state.is_reverse, state.stage2_config);
    for (auto& cr : chains) cr.stage1_score = stage1_score;
    return chains;
}

namespace {

// Fan out one ext_job into the right Stage 1 entry point.  Used by
// run_stage1_jobs.
template <typename KmerInt>
void dispatch_stage1(
    const ExtJob& ej,
    const QueryBundle<KmerInt>& qb,
    const VolumeBundle<KmerInt>& vb,
    int k, const SearchConfig& config, bool both_mode,
    Stage1Buffer& buf,
    JobState& state) {

    const bool is_reverse = (ej.strand_idx != 0);

    if (both_mode) {
        const auto& qd_cod = *qb.qdata_primary;
        const auto& qd_opt = *qb.qdata_secondary;
        if (!is_reverse) {
            auto unify = [](uint32_t a, uint32_t b) {
                if (a > 0 && b > 0) return std::min(a, b);
                return std::max(a, b);
            };
            stage1_one_strand_both<KmerInt>(
                qd_cod.fwd_positions.data(),
                qd_cod.fwd_kmer_values.data(),
                qd_cod.fwd_positions.size(),
                qd_opt.fwd_positions.data(),
                qd_opt.fwd_kmer_values.data(),
                qd_opt.fwd_positions.size(),
                k, false,
                *vb.kix, *vb.kix_opt,
                *vb.filter, config,
                qd_cod.resolved_threshold_fwd,
                qd_opt.resolved_threshold_fwd,
                unify(qd_cod.effective_min_score_fwd,
                      qd_opt.effective_min_score_fwd),
                buf, state);
        } else {
            auto unify = [](uint32_t a, uint32_t b) {
                if (a > 0 && b > 0) return std::min(a, b);
                return std::max(a, b);
            };
            stage1_one_strand_both<KmerInt>(
                qd_cod.rc_positions.data(),
                qd_cod.rc_kmer_values.data(),
                qd_cod.rc_positions.size(),
                qd_opt.rc_positions.data(),
                qd_opt.rc_kmer_values.data(),
                qd_opt.rc_positions.size(),
                k, true,
                *vb.kix, *vb.kix_opt,
                *vb.filter, config,
                qd_cod.resolved_threshold_rc,
                qd_opt.resolved_threshold_rc,
                unify(qd_cod.effective_min_score_rc,
                      qd_opt.effective_min_score_rc),
                buf, state);
        }
    } else {
        const auto& qd = *qb.qdata_primary;
        if (!is_reverse) {
            stage1_one_strand_single<KmerInt>(
                qd.fwd_positions.data(),
                qd.fwd_kmer_values.data(),
                qd.fwd_positions.size(),
                k, false,
                *vb.kix, *vb.filter, config,
                qd.resolved_threshold_fwd,
                qd.effective_min_score_fwd,
                buf, state);
        } else {
            stage1_one_strand_single<KmerInt>(
                qd.rc_positions.data(),
                qd.rc_kmer_values.data(),
                qd.rc_positions.size(),
                k, true,
                *vb.kix, *vb.filter, config,
                qd.resolved_threshold_rc,
                qd.effective_min_score_rc,
                buf, state);
        }
    }
}

template <typename KmerInt>
void dispatch_stage2a(
    const ExtJob& ej,
    const QueryBundle<KmerInt>& qb,
    const VolumeBundle<KmerInt>& vb,
    bool both_mode,
    JobState& state) {

    const bool is_reverse = (ej.strand_idx != 0);

    if (both_mode) {
        const auto& qd_cod = *qb.qdata_primary;
        const auto& qd_opt = *qb.qdata_secondary;
        if (!is_reverse) {
            stage2a_one_strand_both<KmerInt>(
                qd_cod.fwd_positions.data(),
                qd_cod.fwd_kmer_values.data(),
                qd_cod.fwd_positions.size(),
                qd_opt.fwd_positions.data(),
                qd_opt.fwd_kmer_values.data(),
                qd_opt.fwd_positions.size(),
                *vb.kix, *vb.kpx,
                *vb.kix_opt, *vb.kpx_opt,
                state);
        } else {
            stage2a_one_strand_both<KmerInt>(
                qd_cod.rc_positions.data(),
                qd_cod.rc_kmer_values.data(),
                qd_cod.rc_positions.size(),
                qd_opt.rc_positions.data(),
                qd_opt.rc_kmer_values.data(),
                qd_opt.rc_positions.size(),
                *vb.kix, *vb.kpx,
                *vb.kix_opt, *vb.kpx_opt,
                state);
        }
    } else {
        const auto& qd = *qb.qdata_primary;
        if (!is_reverse) {
            stage2a_one_strand_single<KmerInt>(
                qd.fwd_positions.data(),
                qd.fwd_kmer_values.data(),
                qd.fwd_positions.size(),
                *vb.kix, *vb.kpx, state);
        } else {
            stage2a_one_strand_single<KmerInt>(
                qd.rc_positions.data(),
                qd.rc_kmer_values.data(),
                qd.rc_positions.size(),
                *vb.kix, *vb.kpx, state);
        }
    }
}

}  // namespace

template <typename KmerInt>
void run_stage1_jobs(
    const std::vector<size_t>& batch_indices,
    const std::vector<ExtJob>& ext_jobs,
    const std::vector<QueryBundle<KmerInt>>& queries,
    const std::vector<VolumeBundle<KmerInt>>& volumes,
    int k,
    const SearchConfig& config,
    bool both_mode,
    std::vector<JobState>& states,
    tbb::task_arena& arena,
    tbb::enumerable_thread_specific<Stage1Buffer>& tls_bufs) {

    if (batch_indices.empty()) return;

    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, batch_indices.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                auto& buf = tls_bufs.local();
                for (size_t idx = range.begin(); idx != range.end(); ++idx) {
                    size_t ji = batch_indices[idx];
                    const auto& ej = ext_jobs[ji];
                    const auto& qb = queries[ej.qi];
                    const auto& vb = volumes[ej.vi];
                    dispatch_stage1<KmerInt>(ej, qb, vb, k, config, both_mode,
                                              buf, states[ji]);
                }
            });
    });
}

template <typename KmerInt>
void run_stage2a_jobs(
    const std::vector<size_t>& batch_indices,
    const std::vector<ExtJob>& ext_jobs,
    const std::vector<QueryBundle<KmerInt>>& queries,
    const std::vector<VolumeBundle<KmerInt>>& volumes,
    bool both_mode,
    std::vector<JobState>& states,
    tbb::task_arena& arena) {

    if (batch_indices.empty()) return;

    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, batch_indices.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t idx = range.begin(); idx != range.end(); ++idx) {
                    size_t ji = batch_indices[idx];
                    const auto& ej = ext_jobs[ji];
                    const auto& qb = queries[ej.qi];
                    const auto& vb = volumes[ej.vi];
                    dispatch_stage2a<KmerInt>(ej, qb, vb, both_mode, states[ji]);
                }
            });
    });
}

namespace {
// One element of the flat (ext_job, sid) tuple stream consumed by Stage 2B's
// parallel_for.  Stage1 score is captured at flatten time so Stage 2B does
// not need to look up `state.stage1_scores` again.
struct FlatBSlot {
    size_t   ext_job_idx;
    SeqId    sid;
    uint32_t stage1_score;
};
}  // namespace

void run_stage2b_jobs(
    const std::vector<size_t>& batch_indices,
    const std::vector<ExtJob>& ext_jobs,
    const std::vector<JobState>& states,
    const std::vector<uint16_t>& volume_indices,
    tbb::task_arena& arena,
    std::vector<OrchestratorHit>& out_results) {

    if (batch_indices.empty()) return;

    std::vector<FlatBSlot> flat;
    {
        size_t reserve = 0;
        for (size_t ji : batch_indices) reserve += states[ji].candidates.size();
        flat.reserve(reserve);
    }
    for (size_t ji : batch_indices) {
        const JobState& st = states[ji];
        if (st.mode1_only) continue;
        if (st.candidates.empty()) continue;
        if (st.both_mode) {
            for (SeqId sid : st.sorted_candidate_sids) {
                auto sit = st.stage1_scores.find(sid);
                uint32_t sc = (sit != st.stage1_scores.end()) ? sit->second : 0;
                flat.push_back({ ji, sid, sc });
            }
        } else {
            for (const auto& c : st.candidates) {
                flat.push_back({ ji, c.id, c.score });
            }
        }
    }

    if (flat.empty()) return;

    tbb::combinable<std::vector<OrchestratorHit>> tls_results;
    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, flat.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                auto& local = tls_results.local();
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    const auto& slot = flat[i];
                    const ExtJob& ej = ext_jobs[slot.ext_job_idx];
                    const JobState& st = states[slot.ext_job_idx];
                    auto chains = stage2b_one_subject(slot.sid, slot.stage1_score, st);
                    for (auto& cr : chains) {
                        OrchestratorHit oh;
                        oh.query_idx = ej.qi;
                        oh.volume_idx = ej.vi;
                        oh.volume_index = volume_indices[ej.vi];
                        oh.cr = std::move(cr);
                        local.push_back(std::move(oh));
                    }
                }
            });
    });

    tls_results.combine_each([&out_results](std::vector<OrchestratorHit>& local) {
        out_results.insert(out_results.end(),
                           std::make_move_iterator(local.begin()),
                           std::make_move_iterator(local.end()));
    });
}

// Explicit template instantiations.
template void stage1_one_strand_single<uint16_t>(
    const uint32_t*, const uint16_t*, size_t, int, bool,
    const KixReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, JobState&);
template void stage1_one_strand_single<uint32_t>(
    const uint32_t*, const uint32_t*, size_t, int, bool,
    const KixReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, JobState&);

template void stage2a_one_strand_single<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const KpxReader&, JobState&);
template void stage2a_one_strand_single<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const KpxReader&, JobState&);

template void stage1_one_strand_both<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const uint32_t*, const uint16_t*, size_t,
    int, bool,
    const KixReader&, const KixReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, JobState&);
template void stage1_one_strand_both<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const uint32_t*, const uint32_t*, size_t,
    int, bool,
    const KixReader&, const KixReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, JobState&);

template void stage2a_one_strand_both<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    JobState&);
template void stage2a_one_strand_both<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    JobState&);

template void run_stage1_jobs<uint16_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint16_t>>&,
    const std::vector<VolumeBundle<uint16_t>>&,
    int, const SearchConfig&, bool,
    std::vector<JobState>&,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);
template void run_stage1_jobs<uint32_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint32_t>>&,
    const std::vector<VolumeBundle<uint32_t>>&,
    int, const SearchConfig&, bool,
    std::vector<JobState>&,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);

template void run_stage2a_jobs<uint16_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint16_t>>&,
    const std::vector<VolumeBundle<uint16_t>>&,
    bool,
    std::vector<JobState>&,
    tbb::task_arena&);
template void run_stage2a_jobs<uint32_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint32_t>>&,
    const std::vector<VolumeBundle<uint32_t>>&,
    bool,
    std::vector<JobState>&,
    tbb::task_arena&);

} // namespace ikafssn
