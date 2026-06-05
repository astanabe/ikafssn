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
template <Stage1Tier Tier>
inline std::size_t tier_byte_size(uint32_t num_seqs) noexcept {
    using ScoreT = typename Stage1TierTraits<Tier>::ScoreT;
    return static_cast<std::size_t>(num_seqs) * sizeof(ScoreT);
}

template <Stage1Tier Tier>
inline void reset_all_typed(Stage1Buffer& buf) {
    using PosT = typename Stage1TierTraits<Tier>::PosT;
    auto* s = score_ptr<Tier>(buf);
    auto* p = last_pos_ptr<Tier>(buf);
    std::memset(s, 0, static_cast<std::size_t>(buf.capacity) * sizeof(typename Stage1TierTraits<Tier>::ScoreT));
    constexpr PosT sentinel = std::numeric_limits<PosT>::max();
    for (uint32_t i = 0; i < buf.capacity; i++) p[i] = sentinel;
}
} // namespace

void Stage1Buffer::ensure_capacity(uint32_t num_seqs) {
    if (capacity >= num_seqs) return;
    capacity = num_seqs;
    std::size_t bytes = 0;
    switch (tier) {
    case Stage1Tier::T8:  bytes = tier_byte_size<Stage1Tier::T8> (num_seqs); break;
    case Stage1Tier::T16: bytes = tier_byte_size<Stage1Tier::T16>(num_seqs); break;
    case Stage1Tier::T32: bytes = tier_byte_size<Stage1Tier::T32>(num_seqs); break;
    }
    score_data.resize(bytes);
    last_pos_data.resize(bytes);
    reset_all();
}

void Stage1Buffer::reset_all() {
    switch (tier) {
    case Stage1Tier::T8:  reset_all_typed<Stage1Tier::T8> (*this); break;
    case Stage1Tier::T16: reset_all_typed<Stage1Tier::T16>(*this); break;
    case Stage1Tier::T32: reset_all_typed<Stage1Tier::T32>(*this); break;
    }
}

template <Stage1Tier Tier>
void Stage1Buffer::clear_dirty_typed() {
    using PosT = typename Stage1TierTraits<Tier>::PosT;
    using ScoreT = typename Stage1TierTraits<Tier>::ScoreT;
    // When the dirty fraction is high enough, a sequential memset/reset over
    // the whole buffer beats random per-index writes (avoids touching the
    // dirty index list twice and stays linear in cache lines).
    if (dirty.size() * 8 > capacity) {
        reset_all_typed<Tier>(*this);
    } else {
        auto* s = score_ptr<Tier>(*this);
        auto* p = last_pos_ptr<Tier>(*this);
        constexpr PosT sentinel = std::numeric_limits<PosT>::max();
        for (uint32_t idx : dirty) {
            s[idx] = ScoreT{0};
            p[idx] = sentinel;
        }
    }
    dirty.clear();
}

template void Stage1Buffer::clear_dirty_typed<Stage1Tier::T8>();
template void Stage1Buffer::clear_dirty_typed<Stage1Tier::T16>();
template void Stage1Buffer::clear_dirty_typed<Stage1Tier::T32>();

// Internal: accumulate per-(kix, q_pos) updates into buf without clearing dirty.
// Used by both stage1_filter (single template) and stage1_filter_accumulate
// (cross-template "both" mode). HasFilter is a compile-time switch that lets
// the no-filter case drop the per-sid pass() branch entirely.
template <typename KmerInt, Stage1Tier Tier, bool HasFilter>
static void stage1_filter_accumulate_impl(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {

    using PosT = typename Stage1TierTraits<Tier>::PosT;

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0 || n == 0) return;

    const uint8_t* posting_file = kix.posting_file();
    const bool use_coverscore = (config.stage1_score_type == 1);
    const DecodedKmerCache* dc = kix.decode_cache();

    buf.ensure_capacity(num_seqs);
    auto* scores   = score_ptr<Tier>(buf);
    auto* last_pos = last_pos_ptr<Tier>(buf);

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
                        flush_batch_simd<Tier>(lk.sids + base, chunk, q_pos,
                                               scores, last_pos, buf.dirty);
                    }
                } else {
                    for (uint32_t i = 0; i < lk.count; i++) {
                        SeqId sid = lk.sids[i];
                        if (!filter.pass(sid)) continue;
                        sid_batch[batch_count++] = sid;
                        if (batch_count == kBatch) {
                            flush_batch_simd<Tier>(sid_batch, kBatch, q_pos,
                                                   scores, last_pos, buf.dirty);
                            batch_count = 0;
                        }
                    }
                    if (batch_count > 0) {
                        flush_batch_simd<Tier>(sid_batch, batch_count, q_pos,
                                               scores, last_pos, buf.dirty);
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
                    flush_batch_simd<Tier>(sid_batch, kBatch, q_pos,
                                           scores, last_pos, buf.dirty);
                    batch_count = 0;
                }
            }
        }
        if (batch_count > 0) {
            flush_batch_simd<Tier>(sid_batch, batch_count, q_pos,
                                   scores, last_pos, buf.dirty);
            batch_count = 0;
        }
    }
}

// Dispatch by filter mode to one of the HasFilter specializations.
template <typename KmerInt, Stage1Tier Tier>
static inline void dispatch_accumulate_by_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {
    if (filter.mode() == OidFilterMode::kNone) {
        stage1_filter_accumulate_impl<KmerInt, Tier, /*HasFilter=*/false>(
            positions, kmers, n, kix, filter, config, buf);
    } else {
        stage1_filter_accumulate_impl<KmerInt, Tier, /*HasFilter=*/true>(
            positions, kmers, n, kix, filter, config, buf);
    }
}

// Internal: harvest candidates from buf and clear dirty.
template <Stage1Tier Tier>
static std::vector<Stage1Candidate> stage1_filter_finish_impl(
    Stage1Buffer& buf, const Stage1Config& config) {
    using PosT = typename Stage1TierTraits<Tier>::PosT;
    using ScoreT = typename Stage1TierTraits<Tier>::ScoreT;
    auto* scores = score_ptr<Tier>(buf);

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
            reset_all_typed<Tier>(buf);
        } else {
            auto* last_pos = last_pos_ptr<Tier>(buf);
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

    buf.clear_dirty_typed<Tier>();

    auto cmp = [](const Stage1Candidate& a, const Stage1Candidate& b) {
        return a.score > b.score;
    };
    if (candidates.size() > config.stage1_topn) {
        std::nth_element(candidates.begin(),
                         candidates.begin() + config.stage1_topn,
                         candidates.end(), cmp);
        candidates.resize(config.stage1_topn);
    }
    std::sort(candidates.begin(), candidates.end(), cmp);
    return candidates;
}

// Internal implementation with KmerInt + Tier template dispatch.
template <typename KmerInt, Stage1Tier Tier>
static std::vector<Stage1Candidate> stage1_filter_impl(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {

    if (n == 0) return {};
    dispatch_accumulate_by_filter<KmerInt, Tier>(
        positions, kmers, n, kix, filter, config, buf);
    return stage1_filter_finish_impl<Tier>(buf, config);
}

// Public dispatch: selects tier from buffer.
template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {

    switch (buf.tier) {
    case Stage1Tier::T8:
        return stage1_filter_impl<KmerInt, Stage1Tier::T8>(
            positions, kmers, n, kix, filter, config, buf);
    case Stage1Tier::T16:
        return stage1_filter_impl<KmerInt, Stage1Tier::T16>(
            positions, kmers, n, kix, filter, config, buf);
    case Stage1Tier::T32:
    default:
        return stage1_filter_impl<KmerInt, Stage1Tier::T32>(
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
    Stage1Buffer& buf) {

    switch (buf.tier) {
    case Stage1Tier::T8:
        dispatch_accumulate_by_filter<KmerInt, Stage1Tier::T8>(
            positions, kmers, n, kix, filter, config, buf);
        break;
    case Stage1Tier::T16:
        dispatch_accumulate_by_filter<KmerInt, Stage1Tier::T16>(
            positions, kmers, n, kix, filter, config, buf);
        break;
    case Stage1Tier::T32:
    default:
        dispatch_accumulate_by_filter<KmerInt, Stage1Tier::T32>(
            positions, kmers, n, kix, filter, config, buf);
        break;
    }
}

// Harvest candidates after one or more accumulate calls.
std::vector<Stage1Candidate> stage1_filter_finish(
    Stage1Buffer& buf, const Stage1Config& config) {
    switch (buf.tier) {
    case Stage1Tier::T8:
        return stage1_filter_finish_impl<Stage1Tier::T8>(buf, config);
    case Stage1Tier::T16:
        return stage1_filter_finish_impl<Stage1Tier::T16>(buf, config);
    case Stage1Tier::T32:
    default:
        return stage1_filter_finish_impl<Stage1Tier::T32>(buf, config);
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
    Stage1Buffer&);
template void stage1_filter_accumulate<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);

} // namespace ikafssn
