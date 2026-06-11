#pragma once

// Top-level search orchestrator shared by ikafssnsearch (one-shot) and
// ikafssnserver (per-request).  Drives the three search stages with a
// volume-batched group loop:
//
//   Stage 1  — one group loop over volumes (.kix only).
//   Stage 2  — one group loop over volumes (.kix + .kpx).  Only mode 2/3.
//              Each iteration runs Stage 2A then Stage 2B for the group
//              and frees the group's per-ext_job transient JobState
//              (hits_per_seq etc.) before the next group begins, so peak
//              heap usage stays bounded by one group rather than the total
//              volume × query fan-out.
//
// Volumes are bundled into a group until the group's ext_job count reaches
// the thread target (`nthread`): many-query runs get one volume per group,
// few-query runs bundle volumes to saturate the arena.  .kix / .kpx readers
// are opened / closed per group, so the kernel-level page cache footprint is
// bounded by one group's worth of mappings.  Persistent madvise WILLNEED is
// reserved for .khx / .ksx by the caller.

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "io/volume_discovery.hpp"
#include "search/oid_filter.hpp"
#include "search/parallel_search.hpp"
#include "search/search_config.hpp"

namespace ikafssn {

class KixReader;
class KpxReader;
class KsxReader;
class Logger;

struct VolumeMeta {
    DiscoveredVolume files;
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
    Logger* logger = nullptr;
    uint32_t max_num_seqs = 0;
    Stage1Width width = Stage1Width::T32;
};

template <typename KmerInt>
std::vector<OrchestratorHit> run_search(const RunSearchInputs<KmerInt>& in);

extern template std::vector<OrchestratorHit>
run_search<uint16_t>(const RunSearchInputs<uint16_t>&);
extern template std::vector<OrchestratorHit>
run_search<uint32_t>(const RunSearchInputs<uint32_t>&);

}  // namespace ikafssn
