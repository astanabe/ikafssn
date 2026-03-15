#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

struct Stage2Config {
    uint32_t max_gap = 100;         // max diagonal deviation between consecutive chain hits
    uint32_t min_diag_hits = 1;     // diagonal filter threshold (1 = disabled)
    uint32_t min_score = 0;         // minimum chain score to report (0 = adaptive)
    uint32_t chain_max_lookback = 64; // chaining DP lookback window (0 = unlimited O(n²))
    uint32_t max_nhit_per_subject = 1; // max chains per subject (0 = unlimited)
};

// Run Stage 2 chaining on hits for a single candidate sequence.
// 1. Apply diagonal filter
// 2. Sort hits by q_pos (then s_pos)
// 3. Run O(n^2) chaining DP
// 4. Traceback best chain
// 5. If max_nhit_per_subject > 1 (or 0=unlimited), remove used hits and repeat
//
// Returns vector of ChainResult (empty if no chain passes min_score).
std::vector<ChainResult> chain_hits(const std::vector<Hit>& hits,
                                    SeqId seq_id,
                                    int k,
                                    bool is_reverse,
                                    const Stage2Config& config);

} // namespace ikafssn
