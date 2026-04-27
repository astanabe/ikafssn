#include "search/stage1_filter_simd.hpp"

namespace ikafssn {

namespace {
template <Stage1Tier Tier>
inline void flush_scalar(const SeqId* sid_batch, int count,
                         typename Stage1TierTraits<Tier>::PosT q_pos,
                         typename Stage1TierTraits<Tier>::ScoreT* scores,
                         typename Stage1TierTraits<Tier>::PosT* last_pos,
                         std::vector<uint32_t>& dirty) {
    for (int i = 0; i < count; i++) {
        SeqId sid = sid_batch[i];
        if (scores[sid] == 0) dirty.push_back(sid);
        if (last_pos[sid] != q_pos) {
            scores[sid]++;
            last_pos[sid] = q_pos;
        }
    }
}
} // namespace

template <Stage1Tier Tier>
void flush_batch_simd(const SeqId* sid_batch, int count,
                      typename Stage1TierTraits<Tier>::PosT q_pos,
                      void* scores_void, void* last_pos_void,
                      std::vector<uint32_t>& dirty) {
    using ScoreT = typename Stage1TierTraits<Tier>::ScoreT;
    using PosT   = typename Stage1TierTraits<Tier>::PosT;
    auto* scores   = reinterpret_cast<ScoreT*>(scores_void);
    auto* last_pos = reinterpret_cast<PosT*>(last_pos_void);

    // Phase 2a: scalar only. Phase 2b will add per-tier SIMD kernels here.
    flush_scalar<Tier>(sid_batch, count, q_pos, scores, last_pos, dirty);
}

template void flush_batch_simd<Stage1Tier::T8> (
    const SeqId*, int, uint8_t,  void*, void*, std::vector<uint32_t>&);
template void flush_batch_simd<Stage1Tier::T16>(
    const SeqId*, int, uint16_t, void*, void*, std::vector<uint32_t>&);
template void flush_batch_simd<Stage1Tier::T32>(
    const SeqId*, int, uint32_t, void*, void*, std::vector<uint32_t>&);

} // namespace ikafssn
