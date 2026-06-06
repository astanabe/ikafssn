#pragma once

#include "core/types.hpp"
#include "search/stage1_filter.hpp"

#include <cstdint>
#include <vector>

namespace ikafssn {

// Batch-update kernel for the Stage 1 hot loop. SoA scores / last_pos buffers.
//
// Semantics (must hold for every width, bit-exact with scalar):
//   for each sid in sid_batch[0..count):
//     if (Cutoff && scores[sid] + remaining < cutoff_T) continue;  // dead, skip
//     if (scores[sid] == 0) dirty.push_back(sid);
//     if (last_pos[sid] != q_pos) { scores[sid]++; last_pos[sid] = q_pos; }
//
// `dirty` is updated in place; `scores_void` and `last_pos_void` are width-erased
// SoA pointers (same width as Stage1WidthTraits<Width>::ScoreT/PosT).
//
// The cutoff prunes sids that can no longer reach the final threshold `cutoff_T`:
// `remaining` is an upper bound on the additional score a sid can still gain (the
// number of query positions left to process), so a sid with
// `scores[sid] + remaining < cutoff_T` is provably below threshold and is skipped.
// Pruning a dead sid never changes the candidate set (its final score stays below
// cutoff_T whether or not it is pruned) and never affects any other sid.  A
// `cutoff_T <= 1` disables the cutoff (no sid can be dead when the threshold is 1).
//
// `flush_batch_simd<Width>` dispatches to the SIMD kernel matching the
// active CPU tier (or the scalar fallback when SIMD is disabled), selecting the
// cutoff-enabled kernel when `cutoff_T > 1`.
template <Stage1Width Width>
void flush_batch_simd(
    const SeqId* sid_batch,
    int count,
    typename Stage1WidthTraits<Width>::PosT q_pos,
    void* scores_void,
    void* last_pos_void,
    std::vector<uint32_t>& dirty,
    uint32_t remaining,
    uint32_t cutoff_T);

extern template void flush_batch_simd<Stage1Width::T8> (
    const SeqId*, int, uint8_t,  void*, void*, std::vector<uint32_t>&,
    uint32_t, uint32_t);
extern template void flush_batch_simd<Stage1Width::T16>(
    const SeqId*, int, uint16_t, void*, void*, std::vector<uint32_t>&,
    uint32_t, uint32_t);
extern template void flush_batch_simd<Stage1Width::T32>(
    const SeqId*, int, uint32_t, void*, void*, std::vector<uint32_t>&,
    uint32_t, uint32_t);

} // namespace ikafssn
