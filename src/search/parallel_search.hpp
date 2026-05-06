#pragma once

// Stage-split parallel search building blocks.
//
// Each (query, volume, strand) job runs in three stages:
//   Stage 1  — candidate selection from the .kix posting file.  Output is
//              `JobState::candidates` and `JobState::stage1_scores`; in mode 1
//              the per-strand result list `JobState::mode1_results` is also
//              produced.
//   Stage 2A — collect (q_pos, s_pos) hits from the .kpx posting file for
//              every candidate that survived Stage 1.  Populates
//              `JobState::hits_per_seq`.
//   Stage 2B — per-subject chain_hits().  Pure function over read-only
//              JobState; runs as a flat parallel_for over (ext_job, sid)
//              tuples.
//
// `search_orchestrator.{hpp,cpp}` drives the three stages with a global
// volume-batched WILLNEED window: Stage 1 has its own batch loop, and
// Stage 2 has a single batch loop where each iteration runs Stage 2A then
// Stage 2B for the batch — the batch's per-ext_job transient state
// (hits_per_seq etc.) is freed before the next batch begins, so peak
// memory tracks `posting_budget` rather than the total volume × query
// fan-out.  The .kix / .kpx readers are opened / closed per batch so the
// kernel-level page cache is bounded by the configured memory budget.

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/types.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"
#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/task_arena.h>

namespace ikafssn {

class KixReader;
class KpxReader;
class KsxReader;
class OidFilter;

// One element of the (query, volume, strand) work axis.  `strand_idx` is
// 0 for the forward strand and 1 for the reverse-complement strand.
struct ExtJob {
    size_t  qi = 0;          // index into the QueryBundle vector
    size_t  vi = 0;          // index into the VolumeBundle vector
    uint8_t strand_idx = 0;  // 0 = fwd, 1 = rc
};

// Per-(ext_job) intermediate state populated by Stage 1 / Stage 2A and
// consumed by Stage 2B.  When `mode1_only` is true the chain_hits step is
// skipped and `mode1_results` carries the final per-strand results.
struct JobState {
    std::vector<Stage1Candidate> candidates;
    std::vector<SeqId> sorted_candidate_sids;  // populated only in both-mode
    std::unordered_map<SeqId, std::vector<Hit>> hits_per_seq;
    std::unordered_map<SeqId, uint32_t> stage1_scores;
    bool is_reverse = false;
    bool mode1_only = false;
    bool both_mode = false;
    uint32_t effective_min_score = 0;
    int span = 0;  // seed span (t for spaced seeds, k for contiguous)
    Stage2Config stage2_config;
    std::vector<ChainResult> mode1_results;
};

// Stage 1 — single-template path.  Populates state.candidates,
// state.stage1_scores, and (when mode == 1) state.mode1_results.
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
    JobState& state);

extern template void stage1_one_strand_single<uint16_t>(
    const uint32_t*, const uint16_t*, size_t, int, bool,
    const KixReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, JobState&);
extern template void stage1_one_strand_single<uint32_t>(
    const uint32_t*, const uint32_t*, size_t, int, bool,
    const KixReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, JobState&);

// Stage 2A — single-template path.  Collects (q_pos, s_pos) hits for the
// candidates produced by stage1_one_strand_single().  Populates
// state.hits_per_seq.
template <typename KmerInt>
void stage2a_one_strand_single(
    const uint32_t* positions, const KmerInt* kmers, size_t n_kmers,
    const KixReader& kix, const KpxReader& kpx,
    JobState& state);

extern template void stage2a_one_strand_single<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const KpxReader&, JobState&);
extern template void stage2a_one_strand_single<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const KpxReader&, JobState&);

// Stage 1 — both-mode (cross-template) path.  Stage 1 accumulate-then-finish
// runs against a shared Stage1Buffer so per-(sid, q_pos) dedup carries
// across templates.
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
    JobState& state);

extern template void stage1_one_strand_both<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const uint32_t*, const uint16_t*, size_t,
    int, bool,
    const KixReader&, const KixReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, JobState&);
extern template void stage1_one_strand_both<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const uint32_t*, const uint32_t*, size_t,
    int, bool,
    const KixReader&, const KixReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, JobState&);

// Stage 2A — both-mode path.  Hits from coding and optimal indexes are
// merged into one hits_per_seq.
template <typename KmerInt>
void stage2a_one_strand_both(
    const uint32_t* pos_cod, const KmerInt* kmers_cod, size_t n_cod,
    const uint32_t* pos_opt, const KmerInt* kmers_opt, size_t n_opt,
    const KixReader& kix_cod, const KpxReader& kpx_cod,
    const KixReader& kix_opt, const KpxReader& kpx_opt,
    JobState& state);

extern template void stage2a_one_strand_both<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    JobState&);
extern template void stage2a_one_strand_both<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    JobState&);

// Stage 2B — run chain_hits() for one subject in `state`.  Pure function over
// read-only state; safe to call concurrently across distinct (state, sid)
// pairs.  Returns chain results for the subject (empty if no hits or no
// chain meets min_score).
std::vector<ChainResult>
stage2b_one_subject(SeqId sid, uint32_t stage1_score, const JobState& state);

// Sort and truncate a SearchResult per the search config.  Exposed so the
// volume-level wrappers and the orchestrator's per-volume / per-query post
// pass share one implementation.
void sort_and_truncate(SearchResult& result, const SearchConfig& config);

// Volume-side bundle — pointers are non-owning and caller-managed.
// `kix_opt` / `kpx_opt` are nullptr in single-template mode and point to the
// optimal-side readers in both-mode.
template <typename KmerInt>
struct VolumeBundle {
    const KixReader* kix = nullptr;
    const KpxReader* kpx = nullptr;
    const KsxReader* ksx = nullptr;
    const OidFilter* filter = nullptr;
    uint16_t volume_index = 0;
    const KixReader* kix_opt = nullptr;
    const KpxReader* kpx_opt = nullptr;
};

// Query-side bundle — `qdata_secondary` is null in single-template mode and
// points to the optimal-side QueryKmerData in both-mode (qdata_primary then
// holds the coding side).
template <typename KmerInt>
struct QueryBundle {
    const std::string* query_id = nullptr;
    const QueryKmerData<KmerInt>* qdata_primary = nullptr;
    const QueryKmerData<KmerInt>* qdata_secondary = nullptr;
};

// Boundary representation produced by Stage 2B.  Caller maps each element
// back to its output type (OutputHit for ikafssnsearch, ResponseHit for
// ikafssnserver).
struct OrchestratorHit {
    size_t query_idx = 0;     // index into the QueryBundle vector
    size_t volume_idx = 0;    // index into the VolumeBundle vector
    uint16_t volume_index = 0;
    ChainResult cr;
};

// Stage 1 orchestrator — runs `stage1_one_strand_*` for every ext_job in
// `batch_indices` (indices into `ext_jobs`).  States in `states` not
// referenced by `batch_indices` are left untouched.  `volumes` must have
// entries populated for every vi referenced by `batch_indices`.
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
    tbb::enumerable_thread_specific<Stage1Buffer>& tls_bufs);

extern template void run_stage1_jobs<uint16_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint16_t>>&,
    const std::vector<VolumeBundle<uint16_t>>&,
    int, const SearchConfig&, bool,
    std::vector<JobState>&,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);
extern template void run_stage1_jobs<uint32_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint32_t>>&,
    const std::vector<VolumeBundle<uint32_t>>&,
    int, const SearchConfig&, bool,
    std::vector<JobState>&,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);

// Stage 2A orchestrator — runs `stage2a_one_strand_*` for every ext_job in
// `batch_indices`.  Skips ext_jobs whose Stage 1 produced no candidates or
// that ran in mode 1 only.
template <typename KmerInt>
void run_stage2a_jobs(
    const std::vector<size_t>& batch_indices,
    const std::vector<ExtJob>& ext_jobs,
    const std::vector<QueryBundle<KmerInt>>& queries,
    const std::vector<VolumeBundle<KmerInt>>& volumes,
    bool both_mode,
    std::vector<JobState>& states,
    tbb::task_arena& arena);

extern template void run_stage2a_jobs<uint16_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint16_t>>&,
    const std::vector<VolumeBundle<uint16_t>>&,
    bool,
    std::vector<JobState>&,
    tbb::task_arena&);
extern template void run_stage2a_jobs<uint32_t>(
    const std::vector<size_t>&,
    const std::vector<ExtJob>&,
    const std::vector<QueryBundle<uint32_t>>&,
    const std::vector<VolumeBundle<uint32_t>>&,
    bool,
    std::vector<JobState>&,
    tbb::task_arena&);

// Stage 2B orchestrator — flat parallel_for over (ext_job, sid) tuples for
// the ext_jobs in `batch_indices`.  `volume_indices` provides the wire-level
// volume_index per ext_job so the returned OrchestratorHits carry the
// correct value.  Chains are appended to `out_results` (caller-owned
// accumulator) so per-batch invocations can fold into a single result vector
// without intermediate copies.
void run_stage2b_jobs(
    const std::vector<size_t>& batch_indices,
    const std::vector<ExtJob>& ext_jobs,
    const std::vector<JobState>& states,
    const std::vector<uint16_t>& volume_indices,
    tbb::task_arena& arena,
    std::vector<OrchestratorHit>& out_results);

} // namespace ikafssn
