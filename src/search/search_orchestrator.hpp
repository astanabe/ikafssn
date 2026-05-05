#pragma once

// Top-level search orchestrator shared by ikafssnsearch (one-shot) and
// ikafssnserver (per-request).  Drives the three search stages with a
// global volume-batched WILLNEED window:
//
//   Stage 1   — one batch loop over volumes (.kix only).
//   Stage 2A  — one batch loop over volumes (.kix + .kpx).  Only mode 2/3.
//   Stage 2B  — one global parallel_for over (ext_job, sid) tuples.
//
// .kix / .kpx readers are opened / closed per batch so the kernel-level
// page cache footprint is bounded by the configured posting_budget.
// Persistent madvise WILLNEED is reserved for .khx / .ksx by the caller;
// `posting_budget` is the residual budget the orchestrator may spend on
// kix / kpx posting bodies inside one batch.

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "io/volume_discovery.hpp"
#include "search/oid_filter.hpp"
#include "search/parallel_search.hpp"
#include "search/volume_searcher.hpp"

namespace ikafssn {

class KixReader;
class KpxReader;
class KsxReader;
class Logger;

struct VolumeMeta {
    DiscoveredVolume files;
    size_t kix_posting_size = 0;  // posting_file_size (mode-aware cost)
    size_t kpx_posting_size = 0;  // 0 in mode 1
    size_t kix_full_size = 0;     // willneed_size_full
    size_t kpx_full_size = 0;     // 0 in mode 1
    uint16_t volume_index = 0;
    uint32_t num_sequences = 0;
};

template <typename KmerInt>
struct RunSearchInputs {
    // Single-template: only volumes_cod is used (volumes_opt empty).
    // Both-mode: volumes_cod[i] and volumes_opt[i] are the cod/opt sides
    // of volume i (parallel arrays).
    std::vector<VolumeMeta> volumes_cod;
    std::vector<VolumeMeta> volumes_opt;
    std::vector<const KsxReader*> ksx_per_volume;  // size = volumes_cod.size()
    std::vector<OidFilter> oid_filters;            // size = volumes_cod.size()
    const std::vector<QueryBundle<KmerInt>>* queries = nullptr;
    const std::vector<uint8_t>* query_skip_reason = nullptr;
    SearchConfig config;
    bool both_mode = false;
    int k = 0;
    int nthread = 1;
    uint64_t posting_budget = 0;
    Logger* logger = nullptr;
    uint32_t max_num_seqs = 0;
    Stage1Tier tier = Stage1Tier::T32;
};

template <typename KmerInt>
std::vector<OrchestratorHit> run_search(const RunSearchInputs<KmerInt>& in);

extern template std::vector<OrchestratorHit>
run_search<uint16_t>(const RunSearchInputs<uint16_t>&);
extern template std::vector<OrchestratorHit>
run_search<uint32_t>(const RunSearchInputs<uint32_t>&);

}  // namespace ikafssn
