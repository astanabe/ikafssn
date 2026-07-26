#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

struct Stage2Config {
    uint32_t max_gap = 100;         // max diagonal deviation between consecutive chain hits
    uint32_t min_nhit_diag = 1;     // diagonal filter threshold (1 = disabled)
    uint32_t min_score = 0;         // minimum chain score to report (0 = adaptive)
    uint32_t chain_max_lookback = 64; // chaining DP lookback window (0 = unlimited O(n²))
    uint32_t max_nhit_per_subject = 1; // max chains per subject (0 = unlimited)
    // Selection mode for max_nhit_per_subject, 1..4.  The mode picks the
    // tie behaviour of the per-fragment chain extraction (this file) and the
    // parent grouping of the parent-level selector (select_parent_topn):
    //   1 = take-N (no ties),        parent grouping merges strands
    //   2 = take-N (no ties),        parent grouping splits strands
    //   3 = top-N + score ties,      parent grouping merges strands
    //   4 = top-N + score ties,      parent grouping splits strands
    // CLI / server resolve the sentinel 0 (auto) to 3 before reaching here.
    uint8_t max_nhit_per_subject_mode = 3;
    // L: max chains per query across all volumes / strands (0 = unlimited).
    // Applied after the per-subject selector, by chainscore, keeping the top-L
    // plus every chain tying the L-th chainscore.  Bounds the chains that enter
    // Stage 3 alignment.
    uint32_t max_nhit_in_total = 0;
};

// Run Stage 2 chaining on hits for a single candidate sequence.
// 1. Sort hits by q_pos (then s_pos) and drop duplicate (q_pos, s_pos) pairs
// 2. Apply diagonal filter
// 3. Run chaining DP over a lookback window
// 4. Traceback best chain
// 5. If max_nhit_per_subject > 1 (or 0=unlimited), remove used hits and repeat.
//    In a tie-inclusive mode (max_nhit_per_subject_mode 3/4) extraction
//    continues past N while the next chain ties the N-th chain's chainscore.
//
// Chains are written to `out` (cleared first); it ends up empty when no chain
// passes min_score.  Working buffers live in per-thread scratch, so reusing
// one `out` across calls makes the steady state allocation-free.
void chain_hits(const std::vector<Hit>& hits,
                SeqId seq_id,
                int span,
                bool is_reverse,
                const Stage2Config& config,
                std::vector<ChainResult>& out);

} // namespace ikafssn
