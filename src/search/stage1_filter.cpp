#include "search/stage1_filter.hpp"
#include "search/stage1_filter_simd.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/hit_limits.hpp"
#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <type_traits>
#include <utility>

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
    Stage1Buffer& buf) {

    using PosT = typename Stage1WidthTraits<Width>::PosT;

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0 || n == 0) return;

    const uint8_t* posting_file = kix.posting_file();

    buf.ensure_capacity(num_seqs);
    auto* scores   = score_ptr<Width>(buf);
    auto* last_pos = last_pos_ptr<Width>(buf);

    const uint32_t cutoff_rem = config.cutoff_remaining;
    const uint32_t cutoff_thr = config.cutoff_threshold;

    // Hoisted out of the qi loop so the underlying StreamCtx::decoded
    // capacity (which can grow to millions of u32 for conserved k-mers)
    // is reused across all probes on this thread instead of being
    // allocated and freed per probe.
    thread_local SeqIdDecoder decoder;

    for (size_t qi = 0; qi < n; qi++) {
        auto q_pos = static_cast<PosT>(positions[qi]);
        auto kmer_idx = kmers[qi];
        uint64_t off, byte_len;
        kix.posting_list_range(kmer_idx, off, byte_len);
        if (byte_len == 0) continue;

        decoder.reset(posting_file + off, posting_file + off + byte_len);
        decoder.ensure_decoded();
        const SeqId* decoded = decoder.decoded_data();
        const std::size_t n_dec = decoder.decoded_count();

        if constexpr (HasFilter) {
            // Excluded seq_ids have to be squeezed out, so survivors are
            // staged in a small fixed buffer.
            constexpr int kBatch = 16;
            SeqId sid_batch[kBatch];
            int batch_count = 0;
            for (std::size_t i = 0; i < n_dec; i++) {
                SeqId sid = decoded[i];
                if (!filter.pass(sid)) continue;
                sid_batch[batch_count++] = sid;
                if (batch_count == kBatch) {
                    flush_batch_simd<Width>(sid_batch, kBatch, q_pos,
                                           scores, last_pos, buf.dirty,
                                           cutoff_rem, cutoff_thr);
                    batch_count = 0;
                }
            }
            if (batch_count > 0) {
                flush_batch_simd<Width>(sid_batch, batch_count, q_pos,
                                       scores, last_pos, buf.dirty,
                                       cutoff_rem, cutoff_thr);
            }
        } else {
            // .kix postings are distinct per k-mer, so coverscore counts a sid
            // at most once per k-mer and the decoded span goes straight to the
            // scatter kernel, which does the per-query-position dedup
            // (last_pos != q_pos).  Split because decoded_count() is a
            // std::size_t while the kernel takes an int count.
            constexpr std::size_t kMaxBatch = std::size_t{1} << 20;
            for (std::size_t i = 0; i < n_dec; i += kMaxBatch) {
                const std::size_t len = std::min(kMaxBatch, n_dec - i);
                flush_batch_simd<Width>(decoded + i, static_cast<int>(len),
                                       q_pos, scores, last_pos, buf.dirty,
                                       cutoff_rem, cutoff_thr);
            }
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
    Stage1Buffer& buf) {
    if (filter.mode() == OidFilterMode::kNone) {
        stage1_filter_accumulate_impl<KmerInt, Width, /*HasFilter=*/false>(
            positions, kmers, n, kix, filter, config, buf);
    } else {
        stage1_filter_accumulate_impl<KmerInt, Width, /*HasFilter=*/true>(
            positions, kmers, n, kix, filter, config, buf);
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

    // The candidate-count limits (M/N/L) are applied downstream, so finish
    // returns every candidate at or above min_stage1_score.  We know the dirty
    // list is consumed and reset right away, so walk it once and do the
    // per-index reset in the same iteration.  The bulk-reset fallback in
    // clear_dirty_typed is still desirable when dirty is dense, so honour the
    // same threshold here.
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
        positions, kmers, n, kix, filter, config, buf);
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
    Stage1Buffer& buf) {

    switch (buf.width) {
    case Stage1Width::T8:
        dispatch_accumulate_by_filter<KmerInt, Stage1Width::T8>(
            positions, kmers, n, kix, filter, config, buf);
        break;
    case Stage1Width::T16:
        dispatch_accumulate_by_filter<KmerInt, Stage1Width::T16>(
            positions, kmers, n, kix, filter, config, buf);
        break;
    case Stage1Width::T32:
    default:
        dispatch_accumulate_by_filter<KmerInt, Stage1Width::T32>(
            positions, kmers, n, kix, filter, config, buf);
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

// --- Stage 1 candidate-count limits ---------------------------------------

void stage1_limit_topk(std::vector<Stage1Candidate>& candidates, uint32_t k,
                       bool tie_inclusive) {
    if (k == 0 || candidates.size() <= k) return;
    // One group, so the scores go straight into the single-group selection.
    std::vector<uint32_t> scores;
    scores.reserve(candidates.size());
    for (const auto& c : candidates) scores.push_back(c.score);
    const uint32_t kth = in_total_threshold(scores, k);

    std::vector<Stage1Candidate> out;
    out.reserve(k);
    if (tie_inclusive) {
        for (const auto& c : candidates)
            if (c.score >= kth) out.push_back(c);
    } else {
        uint32_t gt = 0;
        for (const auto& c : candidates) if (c.score > kth) ++gt;
        uint32_t need_eq = (k > gt) ? (k - gt) : 0;
        for (const auto& c : candidates) {
            if (c.score > kth) {
                out.push_back(c);
            } else if (c.score == kth && need_eq > 0) {
                out.push_back(c);
                --need_eq;
            }
        }
    }
    candidates = std::move(out);
}

void stage1_limit_per_parent(std::vector<Stage1Candidate>& candidates, uint32_t n,
                             bool tie_inclusive, const uint32_t* parent_index) {
    if (n == 0 || parent_index == nullptr || candidates.empty()) return;

    const size_t cnt = candidates.size();
    // The group key is the parent, kept as a callback so no per-candidate key
    // array is materialised.
    std::vector<char> keep;
    group_topn_keep(
        cnt, n, tie_inclusive,
        [&](uint32_t a, uint32_t b) {
            return parent_index[candidates[a].id] < parent_index[candidates[b].id];
        },
        [&](uint32_t a, uint32_t b) {
            return parent_index[candidates[a].id] == parent_index[candidates[b].id];
        },
        [&](uint32_t a) { return candidates[a].score; },
        keep);

    std::vector<Stage1Candidate> out;
    out.reserve(cnt);
    for (size_t t = 0; t < cnt; ++t)
        if (keep[t]) out.push_back(candidates[t]);
    candidates = std::move(out);
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
