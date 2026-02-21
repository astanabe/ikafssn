#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

struct Stage2Config {
    uint32_t max_gap = 100;         // max diagonal deviation between consecutive chain hits
    uint32_t min_diag_hits = 2;     // diagonal filter threshold
    uint32_t min_score = 0;         // minimum chain score to report (0 = adaptive)
};

// Run Stage 2 chaining on hits for a single candidate sequence.
// 1. Apply diagonal filter
// 2. Sort hits by q_pos (then s_pos)
// 3. Run O(n^2) chaining DP
// 4. Traceback best chain
//
// Returns ChainResult (or empty with score=0 if no chain passes min_score).
ChainResult chain_hits(const std::vector<Hit>& hits,
                       SeqId seq_id,
                       int k,
                       bool is_reverse,
                       const Stage2Config& config);

} // namespace ikafssn
