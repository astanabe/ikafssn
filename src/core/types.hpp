#pragma once

#include <cstdint>

namespace ikafssn {

using SeqId = uint32_t;
using SeqPos = uint32_t;

struct Hit {
    SeqPos q_pos;
    SeqPos s_pos;
};

struct ChainResult {
    SeqId  seq_id;
    uint32_t score;
    uint32_t stage1_score = 0;
    SeqPos q_start;
    SeqPos q_end;
    SeqPos s_start;
    SeqPos s_end;
    bool   is_reverse;
};

// Returns 0 for k <= 8 (uint16_t), 1 for k >= 9 (uint32_t)
inline constexpr uint8_t kmer_type_for_k(int k) {
    return k >= 9 ? 1 : 0;
}

} // namespace ikafssn
