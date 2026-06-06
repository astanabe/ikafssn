#include "search/stage1_filter.hpp"
#include "search/stage1_filter_simd.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/decode_cache.hpp"
#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <type_traits>

namespace ikafssn {

namespace {
template <Stage1Width Width>
inline std::size_t width_byte_size(uint32_t num_seqs) noexcept {
    using ScoreT = typename Stage1WidthTraits<Width>::ScoreT;
    return static_cast<std::size_t>(num_seqs) * sizeof(ScoreT);
}

template <Stage1Width Width>
inline void reset_all_typed(Stage1Buffer& buf) {
    using PosT = typename Stage1WidthTraits<Width>::PosT;
    auto* s = score_ptr<Width>(buf);
    auto* p = last_pos_ptr<Width>(buf);
    std::memset(s, 0, static_cast<std::size_t>(buf.capacity) * sizeof(typename Stage1WidthTraits<Width>::ScoreT));
    constexpr PosT sentinel = std::numeric_limits<PosT>::max();
    for (uint32_t i = 0; i < buf.capacity; i++) p[i] = sentinel;
}
} // namespace

void Stage1Buffer::ensure_capacity(uint32_t num_seqs) {
    if (capacity >= num_seqs) return;
    capacity = num_seqs;
    std::size_t bytes = 0;
    switch (width) {
    case Stage1Width::T8:  bytes = width_byte_size<Stage1Width::T8> (num_seqs); break;
    case Stage1Width::T16: bytes = width_byte_size<Stage1Width::T16>(num_seqs); break;
    case Stage1Width::T32: bytes = width_byte_size<Stage1Width::T32>(num_seqs); break;
    }
    score_data.resize(bytes);
    last_pos_data.resize(bytes);
    reset_all();
}

void Stage1Buffer::reset_all() {
    switch (width) {
    case Stage1Width::T8:  reset_all_typed<Stage1Width::T8> (*this); break;
    case Stage1Width::T16: reset_all_typed<Stage1Width::T16>(*this); break;
    case Stage1Width::T32: reset_all_typed<Stage1Width::T32>(*this); break;
    }
}

template <Stage1Width Width>
void Stage1Buffer::clear_dirty_typed() {
    using PosT = typename Stage1WidthTraits<Width>::PosT;
    using ScoreT = typename Stage1WidthTraits<Width>::ScoreT;
    // When the dirty fraction is high enough, a sequential memset/reset over
    // the whole buffer beats random per-index writes (avoids touching the
    // dirty index list twice and stays linear in cache lines).
    if (dirty.size() * 8 > capacity) {
        reset_all_typed<Width>(*this);
    } else {
        auto* s = score_ptr<Width>(*this);
        auto* p = last_pos_ptr<Width>(*this);
        constexpr PosT sentinel = std::numeric_limits<PosT>::max();
        for (uint32_t idx : dirty) {
            s[idx] = ScoreT{0};
            p[idx] = sentinel;
        }
    }
    dirty.clear();
}

template void Stage1Buffer::clear_dirty_typed<Stage1Width::T8>();
template void Stage1Buffer::clear_dirty_typed<Stage1Width::T16>();
template void Stage1Buffer::clear_dirty_typed<Stage1Width::T32>();

// Internal: accumulate per-(kix, q_pos) updates into buf without clearing dirty.
// Used by both stage1_filter (single template) and stage1_filter_accumulate
// (cross-template "both" mode). HasFilter is a compile-time switch that lets
// the no-filter case drop the per-sid pass() branch entirely.
template <typename KmerInt, Stage1Width Width, bool HasFilter>
static void stage1_filter_accumulate_impl(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf,
    uint32_t remaining) {

    using PosT = typename Stage1WidthTraits<Width>::PosT;

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0 || n == 0) return;

    const uint8_t* posting_file = kix.posting_file();
    const bool use_coverscore = (config.stage1_score_type == 1);
    const DecodedKmerCache* dc = kix.decode_cache();

    // The cutoff is engaged only when the caller passes a positive `remaining`
    // (cost-ordered position groups on a cached volume).  cutoff_T is the
    // resolved final threshold; flush_batch_simd treats cutoff_T <= 1 as off.
    const uint32_t cutoff_T = (remaining > 0) ? config.min_stage1_score : 0;

    buf.ensure_capacity(num_seqs);
    auto* scores   = score_ptr<Width>(buf);
    auto* last_pos = last_pos_ptr<Width>(buf);

    constexpr int kBatch = 16;
    SeqId sid_batch[kBatch];
    int batch_count = 0;

    // Cap a single flush_batch_simd call's count to stay within `int` even
    // for a conserved k-mer present in billions of sequences.  Cached sids
    // are distinct within one posting list, so splitting at any boundary
    // gives the identical per-sid result as a single call.
    constexpr uint32_t kFlushChunk = 1u << 24;

    constexpr int kDecBatch = SeqIdDecoder::kMaxBatch;
    static_assert(kDecBatch >= 1, "kDecBatch must be positive");
    SeqId   raw_sids[kDecBatch];
    uint8_t was_new[kDecBatch];

    // Hoisted out of the qi loop so the underlying StreamCtx::decoded
    // capacity (which can grow to millions of u32 for conserved k-mers)
    // is reused across all probes on this thread instead of being
    // allocated and freed per probe.
    thread_local SeqIdDecoder decoder;

    for (size_t qi = 0; qi < n; qi++) {
        auto q_pos = static_cast<PosT>(positions[qi]);
        auto kmer_idx = kmers[qi];

        if (dc) {
            auto lk = dc->lookup(static_cast<uint32_t>(kmer_idx));
            if (lk.sids != nullptr) {
                // Cached probe: lk.sids is the same ascending-distinct array
                // open_stream_kix would decode, so the flush semantics are
                // independent of chunk boundaries and the result is identical
                // to the uncached path.  An empty posting list matches the
                // uncached `byte_len == 0 => continue`.
                if (lk.count == 0) continue;
                if constexpr (!HasFilter) {
                    // The cache stores distinct seq_ids, so the
                    // `use_coverscore && !was_new` skip is always a no-op
                    // here; there are no intra-list duplicates to drop.
                    for (uint32_t base = 0; base < lk.count; base += kFlushChunk) {
                        int chunk = static_cast<int>(
                            std::min<uint32_t>(kFlushChunk, lk.count - base));
                        flush_batch_simd<Width>(lk.sids + base, chunk, q_pos,
                                               scores, last_pos, buf.dirty,
                                               remaining, cutoff_T);
                    }
                } else {
                    for (uint32_t i = 0; i < lk.count; i++) {
                        SeqId sid = lk.sids[i];
                        if (!filter.pass(sid)) continue;
                        sid_batch[batch_count++] = sid;
                        if (batch_count == kBatch) {
                            flush_batch_simd<Width>(sid_batch, kBatch, q_pos,
                                                   scores, last_pos, buf.dirty,
                                                   remaining, cutoff_T);
                            batch_count = 0;
                        }
                    }
                    if (batch_count > 0) {
                        flush_batch_simd<Width>(sid_batch, batch_count, q_pos,
                                               scores, last_pos, buf.dirty,
                                               remaining, cutoff_T);
                        batch_count = 0;
                    }
                }
                continue;
            }
            // Lookup miss (invariant: every probed k-mer is in the filled
            // unique set, so this is only an insurance path) falls through
            // to the on-the-fly decode below.
        }

        uint64_t off, byte_len;
        kix.posting_list_range(kmer_idx, off, byte_len);
        if (byte_len == 0) continue;

        decoder.reset(posting_file + off, posting_file + off + byte_len);
        while (decoder.has_more()) {
            int n_dec = decoder.next_batch(raw_sids, was_new, kDecBatch);
            if (n_dec == 0) break;
            for (int i = 0; i < n_dec; i++) {
                SeqId sid = raw_sids[i];
                if (use_coverscore && !was_new[i]) continue;
                if constexpr (HasFilter) {
                    if (!filter.pass(sid)) continue;
                }
                sid_batch[batch_count++] = sid;
                if (batch_count == kBatch) {
                    flush_batch_simd<Width>(sid_batch, kBatch, q_pos,
                                           scores, last_pos, buf.dirty,
                                           remaining, cutoff_T);
                    batch_count = 0;
                }
            }
        }
        if (batch_count > 0) {
            flush_batch_simd<Width>(sid_batch, batch_count, q_pos,
                                   scores, last_pos, buf.dirty,
                                   remaining, cutoff_T);
            batch_count = 0;
        }
    }
}

// Dispatch by filter mode to one of the HasFilter specializations.
template <typename KmerInt, Stage1Width Width>
static inline void dispatch_accumulate_by_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf,
    uint32_t remaining) {
    if (filter.mode() == OidFilterMode::kNone) {
        stage1_filter_accumulate_impl<KmerInt, Width, /*HasFilter=*/false>(
            positions, kmers, n, kix, filter, config, buf, remaining);
    } else {
        stage1_filter_accumulate_impl<KmerInt, Width, /*HasFilter=*/true>(
            positions, kmers, n, kix, filter, config, buf, remaining);
    }
}

// Internal: harvest candidates from buf and clear dirty.
template <Stage1Width Width>
static std::vector<Stage1Candidate> stage1_filter_finish_impl(
    Stage1Buffer& buf, const Stage1Config& config) {
    using PosT = typename Stage1WidthTraits<Width>::PosT;
    using ScoreT = typename Stage1WidthTraits<Width>::ScoreT;
    auto* scores = score_ptr<Width>(buf);

    std::vector<Stage1Candidate> candidates;

    if (config.stage1_topn == 0) {
        // Fused path: in the unlimited-topn case we know the dirty list
        // will be consumed and reset right away, so walk it once and do
        // the per-index reset in the same iteration (instead of running
        // a separate clear_dirty_typed pass over the same indices).
        // The bulk-reset fallback in clear_dirty_typed is still desirable
        // when dirty is dense, so honour the same threshold here.
        if (buf.dirty.size() * 8 > buf.capacity) {
            for (uint32_t sid : buf.dirty) {
                if (scores[sid] >= config.min_stage1_score) {
                    candidates.push_back({sid, static_cast<uint32_t>(scores[sid])});
                }
            }
            reset_all_typed<Width>(buf);
        } else {
            auto* last_pos = last_pos_ptr<Width>(buf);
            constexpr PosT sentinel = std::numeric_limits<PosT>::max();
            for (uint32_t sid : buf.dirty) {
                ScoreT s = scores[sid];
                if (s >= config.min_stage1_score) {
                    candidates.push_back({sid, static_cast<uint32_t>(s)});
                }
                scores[sid] = ScoreT{0};
                last_pos[sid] = sentinel;
            }
        }
        buf.dirty.clear();
        return candidates;
    }

    for (uint32_t sid : buf.dirty) {
        if (scores[sid] >= config.min_stage1_score) {
            candidates.push_back({sid, static_cast<uint32_t>(scores[sid])});
        }
    }

    buf.clear_dirty_typed<Width>();

    // Top-N is ties-inclusive: keep every candidate whose score is at least the
    // N-th largest score, so candidates tied at the cutoff are all retained
    // (the result may exceed stage1_topn).  This makes the candidate set a pure
    // function of the scores, independent of the order in which positions were
    // accumulated (the cost-ordered reorder only permutes the dirty list).
    auto by_score = [](const Stage1Candidate& a, const Stage1Candidate& b) {
        return a.score > b.score;
    };
    if (candidates.size() > config.stage1_topn) {
        std::nth_element(candidates.begin(),
                         candidates.begin() + (config.stage1_topn - 1),
                         candidates.end(), by_score);
        uint32_t threshold = candidates[config.stage1_topn - 1].score;
        candidates.erase(
            std::remove_if(candidates.begin(), candidates.end(),
                           [&](const Stage1Candidate& c) { return c.score < threshold; }),
            candidates.end());
    }
    // Deterministic final order: score descending, then sid ascending.
    std::sort(candidates.begin(), candidates.end(),
              [](const Stage1Candidate& a, const Stage1Candidate& b) {
                  return a.score != b.score ? a.score > b.score : a.id < b.id;
              });
    return candidates;
}

// Internal implementation with KmerInt + Width template dispatch.
template <typename KmerInt, Stage1Width Width>
static std::vector<Stage1Candidate> stage1_filter_impl(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {

    if (n == 0) return {};
    dispatch_accumulate_by_filter<KmerInt, Width>(
        positions, kmers, n, kix, filter, config, buf, /*remaining=*/0);
    return stage1_filter_finish_impl<Width>(buf, config);
}

// Public dispatch: selects width from buffer.
template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {

    switch (buf.width) {
    case Stage1Width::T8:
        return stage1_filter_impl<KmerInt, Stage1Width::T8>(
            positions, kmers, n, kix, filter, config, buf);
    case Stage1Width::T16:
        return stage1_filter_impl<KmerInt, Stage1Width::T16>(
            positions, kmers, n, kix, filter, config, buf);
    case Stage1Width::T32:
    default:
        return stage1_filter_impl<KmerInt, Stage1Width::T32>(
            positions, kmers, n, kix, filter, config, buf);
    }
}

// Cross-template accumulation entry point (no dirty clear).
template <typename KmerInt>
void stage1_filter_accumulate(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf,
    uint32_t remaining) {

    switch (buf.width) {
    case Stage1Width::T8:
        dispatch_accumulate_by_filter<KmerInt, Stage1Width::T8>(
            positions, kmers, n, kix, filter, config, buf, remaining);
        break;
    case Stage1Width::T16:
        dispatch_accumulate_by_filter<KmerInt, Stage1Width::T16>(
            positions, kmers, n, kix, filter, config, buf, remaining);
        break;
    case Stage1Width::T32:
    default:
        dispatch_accumulate_by_filter<KmerInt, Stage1Width::T32>(
            positions, kmers, n, kix, filter, config, buf, remaining);
        break;
    }
}

// Harvest candidates after one or more accumulate calls.
std::vector<Stage1Candidate> stage1_filter_finish(
    Stage1Buffer& buf, const Stage1Config& config) {
    switch (buf.width) {
    case Stage1Width::T8:
        return stage1_filter_finish_impl<Stage1Width::T8>(buf, config);
    case Stage1Width::T16:
        return stage1_filter_finish_impl<Stage1Width::T16>(buf, config);
    case Stage1Width::T32:
    default:
        return stage1_filter_finish_impl<Stage1Width::T32>(buf, config);
    }
}

// Explicit template instantiations
template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);
template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);

template void stage1_filter_accumulate<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&, uint32_t);
template void stage1_filter_accumulate<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&, uint32_t);

} // namespace ikafssn
