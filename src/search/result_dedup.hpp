#pragma once

// v10 (Phase 3): result dedup helpers for fragment-based indexing.
//
// Fragment splitting introduces overlap regions between adjacent fragments
// of the same parent OID; a single chain can therefore be reported once
// per fragment that covers the same parent-relative region.  The two
// helpers below collapse those duplicates back to one row.
//
//   dedup_stage2_orchestrator_hits()
//       Applied just before run_search() returns in modes 2 and 3.  Keys
//       hits by (query_idx, volume_idx, parent_idx, is_reverse,
//       parent-relative sstart, parent-relative send, chainscore).  Stage 2
//       coordinates are mapped to parent-relative space via the matching
//       KsxReader::fragment_start().
//
//   dedup_stage3_output_hits()
//       Applied in callers after run_stage3() (mode 3 only).  Keys hits by
//       (qseqid, sseqid, sstrand, send, alnscore); the plan assumes
//       traceback=0 so qstart / qend / sstart can shift freely between
//       duplicates.  Skip-marker rows (skip_reason != 0) are not
//       deduplicated and are preserved in their original order.

#include <vector>

namespace ikafssn {

class KsxReader;
struct OrchestratorHit;
struct OutputHit;

void dedup_stage2_orchestrator_hits(
    std::vector<OrchestratorHit>& hits,
    const std::vector<const KsxReader*>& ksx_per_volume);

void dedup_stage3_output_hits(std::vector<OutputHit>& hits);

}  // namespace ikafssn
