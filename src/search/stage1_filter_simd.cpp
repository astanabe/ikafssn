#include "search/stage1_filter_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>

#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace ikafssn {

namespace {

template <Stage1Tier Tier>
void flush_scalar_impl(const SeqId* sid_batch, int count,
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

#if defined(__x86_64__)
// Full 16-lane AVX-512 batch update for the T32 tier. Uses CD's vpconflictd to
// detect duplicates within the batch — duplicate lanes (those with at least one
// lower-index lane sharing the same sid) fall back to scalar so that the
// post-scatter state is observed, which preserves the scalar iteration order
// semantics:
//   for each sid in sid_batch:
//     if (scores[sid] == 0) dirty.push_back(sid);
//     if (last_pos[sid] != q_pos) { scores[sid]++; last_pos[sid] = q_pos; }
//
// The kernel handles arbitrary count (>= 16 entry path required by dispatcher,
// but the loop and tail naturally absorb any size) so it is safe to call from
// tests with non-multiple-of-16 lengths.
__attribute__((target("avx512cd,avx512bw,avx512f")))
void flush_avx512_t32(const SeqId* sid_batch, int count,
                      uint32_t q_pos,
                      uint32_t* scores, uint32_t* last_pos,
                      std::vector<uint32_t>& dirty) {
    const __m512i qp   = _mm512_set1_epi32(static_cast<int32_t>(q_pos));
    const __m512i ones = _mm512_set1_epi32(1);
    const __m512i zero = _mm512_setzero_si512();

    int i = 0;
    while (i + 16 <= count) {
        const __m512i sids = _mm512_loadu_si512(
            reinterpret_cast<const __m512i*>(sid_batch + i));

        // For each lane, conflict[j] is a bitmap of lower lanes sharing the
        // same value. unique_mask[j] = (conflict[j] == 0).
        const __m512i conflict = _mm512_conflict_epi32(sids);
        const __mmask16 unique_mask = _mm512_cmpeq_epi32_mask(conflict, zero);

        // Gather last_pos[sid] for unique lanes.
        const __m512i lp = _mm512_mask_i32gather_epi32(
            zero, unique_mask, sids, last_pos, 4);
        // update_mask = unique & (last_pos != q_pos)
        const __mmask16 update_mask =
            _mm512_mask_cmpneq_epi32_mask(unique_mask, lp, qp);

        // Gather scores[sid] for unique lanes (needed for both ++ and == 0).
        const __m512i sc_old = _mm512_mask_i32gather_epi32(
            zero, unique_mask, sids, scores, 4);
        // dirty push (unique lanes where scores_old == 0).
        const __mmask16 zero_mask =
            _mm512_mask_cmpeq_epi32_mask(unique_mask, sc_old, zero);

        // sc_new[update lanes] = sc_old + 1.
        const __m512i sc_new = _mm512_mask_add_epi32(sc_old, update_mask,
                                                     sc_old, ones);

        // Scatter: scores and last_pos for update_mask lanes. unique_mask
        // guarantees no within-batch index collisions on the scattered indices.
        _mm512_mask_i32scatter_epi32(scores,   update_mask, sids, sc_new, 4);
        _mm512_mask_i32scatter_epi32(last_pos, update_mask, sids, qp,     4);

        // Push dirty entries in lane order (matches scalar iteration order).
        if (zero_mask != 0) {
            uint32_t bm = static_cast<uint32_t>(zero_mask);
            while (bm) {
                int bit = __builtin_ctz(bm);
                dirty.push_back(static_cast<uint32_t>(sid_batch[i + bit]));
                bm &= bm - 1;
            }
        }

        // Conflict lanes (later occurrences of an earlier sid in the same
        // chunk) — process scalar so they observe the post-scatter state.
        const __mmask16 conflict_mask =
            static_cast<__mmask16>((~static_cast<uint32_t>(unique_mask)) & 0xFFFFu);
        if (conflict_mask != 0) {
            uint32_t bm = static_cast<uint32_t>(conflict_mask);
            while (bm) {
                int bit = __builtin_ctz(bm);
                SeqId sid = sid_batch[i + bit];
                if (scores[sid] == 0) dirty.push_back(sid);
                if (last_pos[sid] != q_pos) {
                    scores[sid]++;
                    last_pos[sid] = q_pos;
                }
                bm &= bm - 1;
            }
        }

        i += 16;
    }

    // Tail (count - i < 16): scalar.
    for (; i < count; i++) {
        SeqId sid = sid_batch[i];
        if (scores[sid] == 0) dirty.push_back(sid);
        if (last_pos[sid] != q_pos) {
            scores[sid]++;
            last_pos[sid] = q_pos;
        }
    }
}
#endif // __x86_64__

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

#if defined(__x86_64__)
    if constexpr (Tier == Stage1Tier::T32) {
        // AVX-512 BW/VBMI/VBMI2 share one kernel — VBMI's vpermb and
        // VBMI2's compress/expand do not help when all values are uint32_t
        // and updates are gather+conflict+scatter shaped. AVX2 lacks
        // scatter, so it falls through to scalar.
        if (count >= 16) {
            SimdCap cap = current_simd_cap();
            if (cap >= SimdCap::AVX512BW) {
                flush_avx512_t32(sid_batch, count, q_pos, scores, last_pos, dirty);
                return;
            }
        }
    }
#endif

    flush_scalar_impl<Tier>(sid_batch, count, q_pos, scores, last_pos, dirty);
}

template void flush_batch_simd<Stage1Tier::T8> (
    const SeqId*, int, uint8_t,  void*, void*, std::vector<uint32_t>&);
template void flush_batch_simd<Stage1Tier::T16>(
    const SeqId*, int, uint16_t, void*, void*, std::vector<uint32_t>&);
template void flush_batch_simd<Stage1Tier::T32>(
    const SeqId*, int, uint32_t, void*, void*, std::vector<uint32_t>&);

} // namespace ikafssn
