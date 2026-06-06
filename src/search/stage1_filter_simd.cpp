#include "search/stage1_filter_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>

#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace ikafssn {

namespace {

template <Stage1Width Width, bool Cutoff>
void flush_scalar_impl(const SeqId* sid_batch, int count,
                       typename Stage1WidthTraits<Width>::PosT q_pos,
                       typename Stage1WidthTraits<Width>::ScoreT* scores,
                       typename Stage1WidthTraits<Width>::PosT* last_pos,
                       std::vector<uint32_t>& dirty,
                       uint32_t remaining, uint32_t cutoff_T) {
    for (int i = 0; i < count; i++) {
        SeqId sid = sid_batch[i];
        if constexpr (Cutoff) {
            if (static_cast<uint32_t>(scores[sid]) + remaining < cutoff_T) continue;
        }
        if (scores[sid] == 0) dirty.push_back(sid);
        if (last_pos[sid] != q_pos) {
            scores[sid]++;
            last_pos[sid] = q_pos;
        }
    }
    (void)remaining;
    (void)cutoff_T;
}

#if defined(__x86_64__)
// Full 16-lane AVX-512 batch update for the T32 width. Uses CD's vpconflictd to
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
//
// `Cutoff` (NTTP) adds the dead-sid pruning described in the header: a unique
// lane whose `scores[sid] + remaining < cutoff_T` can never reach the final
// threshold, so its last_pos gather, both scatters, and dirty push are all
// elided.  When `Cutoff == false` the kernel emits the identical instruction
// sequence as the cutoff-free path (the if constexpr branches collapse away).
template <bool Cutoff>
__attribute__((target("avx512cd,avx512bw,avx512f")))
void flush_avx512_t32(const SeqId* sid_batch, int count,
                      uint32_t q_pos,
                      uint32_t* scores, uint32_t* last_pos,
                      std::vector<uint32_t>& dirty,
                      uint32_t remaining, uint32_t cutoff_T) {
    const __m512i qp   = _mm512_set1_epi32(static_cast<int32_t>(q_pos));
    const __m512i ones = _mm512_set1_epi32(1);
    const __m512i zero = _mm512_setzero_si512();
    [[maybe_unused]] const __m512i rem  =
        _mm512_set1_epi32(static_cast<int32_t>(remaining));
    [[maybe_unused]] const __m512i tvec =
        _mm512_set1_epi32(static_cast<int32_t>(cutoff_T));

    int i = 0;
    while (i + 16 <= count) {
        const __m512i sids = _mm512_loadu_si512(
            reinterpret_cast<const __m512i*>(sid_batch + i));

        // For each lane, conflict[j] is a bitmap of lower lanes sharing the
        // same value. unique_mask[j] = (conflict[j] == 0).
        const __m512i conflict = _mm512_conflict_epi32(sids);
        const __mmask16 unique_mask = _mm512_cmpeq_epi32_mask(conflict, zero);

        // alive_mask drops dead lanes when Cutoff is on; otherwise it is
        // exactly unique_mask so the remaining instruction sequence matches
        // the cutoff-free path.
        __mmask16 alive_mask = unique_mask;
        __m512i sc_old = zero;
        if constexpr (Cutoff) {
            // Gather scores[sid] for unique lanes (needed up front for the
            // dead test), then dead = (sc_old + remaining < cutoff_T).
            sc_old = _mm512_mask_i32gather_epi32(zero, unique_mask, sids, scores, 4);
            const __m512i sc_plus = _mm512_add_epi32(sc_old, rem);
            const __mmask16 dead =
                _mm512_mask_cmplt_epi32_mask(unique_mask, sc_plus, tvec);
            alive_mask = static_cast<__mmask16>(
                static_cast<uint32_t>(unique_mask) & ~static_cast<uint32_t>(dead));
        }

        // Gather last_pos[sid] for alive lanes.
        const __m512i lp = _mm512_mask_i32gather_epi32(
            zero, alive_mask, sids, last_pos, 4);
        // update_mask = alive & (last_pos != q_pos)
        const __mmask16 update_mask =
            _mm512_mask_cmpneq_epi32_mask(alive_mask, lp, qp);

        // Gather scores[sid] for alive lanes (needed for both ++ and == 0).
        // The cutoff path already gathered sc_old for unique lanes above and
        // its alive-lane values are still valid, so only the cutoff-free path
        // needs the gather here.
        if constexpr (!Cutoff) {
            sc_old = _mm512_mask_i32gather_epi32(zero, alive_mask, sids, scores, 4);
        }
        // dirty push (alive lanes where scores_old == 0).
        const __mmask16 zero_mask =
            _mm512_mask_cmpeq_epi32_mask(alive_mask, sc_old, zero);

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
                if constexpr (Cutoff) {
                    if (scores[sid] + remaining < cutoff_T) { bm &= bm - 1; continue; }
                }
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
        if constexpr (Cutoff) {
            if (scores[sid] + remaining < cutoff_T) continue;
        }
        if (scores[sid] == 0) dirty.push_back(sid);
        if (last_pos[sid] != q_pos) {
            scores[sid]++;
            last_pos[sid] = q_pos;
        }
    }
    (void)remaining;
    (void)cutoff_T;
}
#endif // __x86_64__

} // namespace

template <Stage1Width Width>
void flush_batch_simd(const SeqId* sid_batch, int count,
                      typename Stage1WidthTraits<Width>::PosT q_pos,
                      void* scores_void, void* last_pos_void,
                      std::vector<uint32_t>& dirty,
                      uint32_t remaining, uint32_t cutoff_T) {
    using ScoreT = typename Stage1WidthTraits<Width>::ScoreT;
    using PosT   = typename Stage1WidthTraits<Width>::PosT;
    auto* scores   = reinterpret_cast<ScoreT*>(scores_void);
    auto* last_pos = reinterpret_cast<PosT*>(last_pos_void);

    // cutoff_T <= 1 disables the cutoff: with threshold 1 no sid can be dead.
    const bool cutoff = (cutoff_T > 1);

#if defined(__x86_64__)
    if constexpr (Width == Stage1Width::T32) {
        // AVX-512 BW/VBMI/VBMI2 share one kernel — VBMI's vpermb and
        // VBMI2's compress/expand do not help when all values are uint32_t
        // and updates are gather+conflict+scatter shaped. AVX2 lacks
        // scatter, so it falls through to scalar.
        if (count >= 16) {
            SimdCap cap = current_simd_cap();
            if (cap >= SimdCap::AVX512BW) {
                if (cutoff) {
                    flush_avx512_t32</*Cutoff=*/true>(sid_batch, count, q_pos,
                                                      scores, last_pos, dirty,
                                                      remaining, cutoff_T);
                } else {
                    flush_avx512_t32</*Cutoff=*/false>(sid_batch, count, q_pos,
                                                       scores, last_pos, dirty,
                                                       remaining, cutoff_T);
                }
                return;
            }
        }
    }
#endif

    if (cutoff) {
        flush_scalar_impl<Width, /*Cutoff=*/true>(sid_batch, count, q_pos,
                                                 scores, last_pos, dirty,
                                                 remaining, cutoff_T);
    } else {
        flush_scalar_impl<Width, /*Cutoff=*/false>(sid_batch, count, q_pos,
                                                  scores, last_pos, dirty,
                                                  remaining, cutoff_T);
    }
}

template void flush_batch_simd<Stage1Width::T8> (
    const SeqId*, int, uint8_t,  void*, void*, std::vector<uint32_t>&,
    uint32_t, uint32_t);
template void flush_batch_simd<Stage1Width::T16>(
    const SeqId*, int, uint16_t, void*, void*, std::vector<uint32_t>&,
    uint32_t, uint32_t);
template void flush_batch_simd<Stage1Width::T32>(
    const SeqId*, int, uint32_t, void*, void*, std::vector<uint32_t>&,
    uint32_t, uint32_t);

} // namespace ikafssn
