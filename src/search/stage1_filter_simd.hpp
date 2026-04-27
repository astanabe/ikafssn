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
// Phase 2a: scalar fallback only — `flush_batch_simd<Tier>` always dispatches
// to the scalar implementation. Phase 2b will add SIMD-accelerated kernels.
template <Stage1Tier Tier>
void flush_batch_simd(
    const SeqId* sid_batch,
    int count,
    typename Stage1TierTraits<Tier>::PosT q_pos,
    void* scores_void,
    void* last_pos_void,
    std::vector<uint32_t>& dirty);

extern template void flush_batch_simd<Stage1Tier::T8> (
    const SeqId*, int, uint8_t,  void*, void*, std::vector<uint32_t>&);
extern template void flush_batch_simd<Stage1Tier::T16>(
    const SeqId*, int, uint16_t, void*, void*, std::vector<uint32_t>&);
extern template void flush_batch_simd<Stage1Tier::T32>(
    const SeqId*, int, uint32_t, void*, void*, std::vector<uint32_t>&);

} // namespace ikafssn
