#pragma once

#include "core/types.hpp"
#include "search/stage1_filter.hpp"

#include <cstdint>
#include <vector>

namespace ikafssn {

// Batch-update kernel for the Stage 1 hot loop. SoA scores / last_pos buffers.
//
// Semantics (must hold for every tier, bit-exact with scalar):
//   for each sid in sid_batch[0..count):
//     if (scores[sid] == 0) dirty.push_back(sid);
//     if (last_pos[sid] != q_pos) { scores[sid]++; last_pos[sid] = q_pos; }
//
// `dirty` is updated in place; `scores_void` and `last_pos_void` are tier-erased
// SoA pointers (same width as Stage1TierTraits<Tier>::ScoreT/PosT).
//
// `flush_batch_simd<Tier>` dispatches to the SIMD kernel matching the
// active CPU tier (or the scalar fallback when SIMD is disabled).
//
// `cutoff_remaining` / `cutoff_threshold` carry the C9 score cutoff: a lane
// whose gathered score + cutoff_remaining < cutoff_threshold is dead (can never
// reach the final threshold) and is skipped entirely (no score++/last_pos/dirty).
// The default 0 / 0 disables the cutoff, so existing callers are unchanged.
template <Stage1Tier Tier>
void flush_batch_simd(
    const SeqId* sid_batch,
    int count,
    typename Stage1TierTraits<Tier>::PosT q_pos,
    void* scores_void,
    void* last_pos_void,
    std::vector<uint32_t>& dirty,
    uint32_t cutoff_remaining = 0,
    uint32_t cutoff_threshold = 0);

extern template void flush_batch_simd<Stage1Tier::T8> (
    const SeqId*, int, uint8_t,  void*, void*, std::vector<uint32_t>&, uint32_t, uint32_t);
extern template void flush_batch_simd<Stage1Tier::T16>(
    const SeqId*, int, uint16_t, void*, void*, std::vector<uint32_t>&, uint32_t, uint32_t);
extern template void flush_batch_simd<Stage1Tier::T32>(
    const SeqId*, int, uint32_t, void*, void*, std::vector<uint32_t>&, uint32_t, uint32_t);

} // namespace ikafssn
