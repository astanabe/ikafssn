#pragma once

// Parallel search orchestrator.
//
// Each (query, volume) job is split into two phases:
//   Phase A — Stage 1 candidate selection + Stage 2A per-(kmer) hit
//             collection.  Output is `PhaseAState` per strand.
//   Phase B — Stage 2B chain_hits().  Pure function over read-only
//             PhaseAState; runs as a flat parallel_for over all
//             (job, strand, sid) tuples across the request.
//
// The orchestrator runs Phase A as a parallel_for over jobs and then a
// second parallel_for over the union of (job, strand, sid) tuples.  TBB
// work-stealing inside the second phase naturally allocates more threads
// to queries with more candidates, replacing the prior adaptive
// query-level vs (query, volume)-level scheduling heuristic.

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

// Per-(query, volume, strand) intermediate state populated by Phase A and
// consumed by Phase B.  When `mode1_only` is true the chain_hits step is
// skipped and `mode1_results` carries the final per-strand results.
struct PhaseAState {
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

// Phase A — single-template path.  Stage 1 + Stage 2A hit collection.
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
    PhaseAState& state);

extern template void phase_a_one_strand_single<uint16_t>(
    const uint32_t*, const uint16_t*, size_t, int, bool,
    const KixReader&, const KpxReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);
extern template void phase_a_one_strand_single<uint32_t>(
    const uint32_t*, const uint32_t*, size_t, int, bool,
    const KixReader&, const KpxReader&, const OidFilter&,
    const SearchConfig&, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);

// Phase A — both-mode (cross-template) path.  Stage 1 accumulate-then-finish
// runs against a shared Stage1Buffer so per-(sid, q_pos) dedup carries
// across templates.  Stage 2A hits from coding and optimal indexes are
// merged into one hits_per_seq.
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
    PhaseAState& state);

extern template void phase_a_one_strand_both<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const uint32_t*, const uint16_t*, size_t,
    int, bool,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);
extern template void phase_a_one_strand_both<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const uint32_t*, const uint32_t*, size_t,
    int, bool,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const OidFilter&, const SearchConfig&,
    uint32_t, uint32_t, uint32_t, Stage1Buffer&, PhaseAState&);

// Phase B — run chain_hits() for one subject in `state`.  Pure function over
// read-only state; safe to call concurrently across distinct (state, sid)
// pairs.  Returns chain results for the subject (empty if no hits or no
// chain meets min_score).  The caller assigns sid's stage1 score.
std::vector<ChainResult>
phase_b_one_subject(SeqId sid, uint32_t stage1_score, const PhaseAState& state);

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

// Boundary representation produced by `run_search_jobs`. Caller maps each
// element back to its output type (OutputHit for ikafssnsearch,
// ResponseHit for ikafssnserver).
struct OrchestratorHit {
    size_t query_idx = 0;     // index into the QueryBundle vector
    size_t volume_idx = 0;    // index into the VolumeBundle vector
    uint16_t volume_index = 0;
    ChainResult cr;
};

// Top-level orchestrator: Phase A parallel_for over (q, v), then flat
// parallel_for over (job, strand, sid).  Returns a flat vector of
// OrchestratorHit values aggregated across all jobs.
//
// `tls_bufs` is shared by reference so the caller controls its lifetime;
// `arena` is the TBB task arena both phases run inside.
template <typename KmerInt>
std::vector<OrchestratorHit> run_search_jobs(
    const std::vector<QueryBundle<KmerInt>>& queries,
    const std::vector<VolumeBundle<KmerInt>>& volumes,
    const std::vector<std::pair<size_t, size_t>>& jobs,
    int k,
    const SearchConfig& config,
    bool both_mode,
    tbb::task_arena& arena,
    tbb::enumerable_thread_specific<Stage1Buffer>& tls_bufs);

extern template std::vector<OrchestratorHit> run_search_jobs<uint16_t>(
    const std::vector<QueryBundle<uint16_t>>&,
    const std::vector<VolumeBundle<uint16_t>>&,
    const std::vector<std::pair<size_t, size_t>>&,
    int, const SearchConfig&, bool,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);
extern template std::vector<OrchestratorHit> run_search_jobs<uint32_t>(
    const std::vector<QueryBundle<uint32_t>>&,
    const std::vector<VolumeBundle<uint32_t>>&,
    const std::vector<std::pair<size_t, size_t>>&,
    int, const SearchConfig&, bool,
    tbb::task_arena&,
    tbb::enumerable_thread_specific<Stage1Buffer>&);

} // namespace ikafssn
