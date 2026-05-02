#include "search/stage1_filter.hpp"
#include "search/stage1_filter_simd.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
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
    auto* s = score_ptr<Tier>(*this);
    auto* p = last_pos_ptr<Tier>(*this);
    constexpr PosT sentinel = std::numeric_limits<PosT>::max();
    for (uint32_t idx : dirty) {
        s[idx] = ScoreT{0};
        p[idx] = sentinel;
    }
    dirty.clear();
}

template void Stage1Buffer::clear_dirty_typed<Stage1Tier::T8>();
template void Stage1Buffer::clear_dirty_typed<Stage1Tier::T16>();
template void Stage1Buffer::clear_dirty_typed<Stage1Tier::T32>();

// Internal implementation with KmerInt + Tier template dispatch.
template <typename KmerInt, Stage1Tier Tier>
static std::vector<Stage1Candidate> stage1_filter_impl(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf) {

    using PosT = typename Stage1TierTraits<Tier>::PosT;

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0 || n == 0) return {};

    const uint8_t* posting_file = kix.posting_file();
    const bool use_coverscore = (config.stage1_score_type == 1);

    buf.ensure_capacity(num_seqs);
    auto* scores   = score_ptr<Tier>(buf);
    auto* last_pos = last_pos_ptr<Tier>(buf);

    // Batch sids per (kix, q_pos) so flush_batch_simd<Tier> can apply gather/
    // scatter SIMD updates in one shot. q_pos changes per qi iteration, so the
    // partial batch must be flushed at the end of each qi to keep the scalar
    // semantics of comparing last_pos[sid] against the current q_pos.
    constexpr int kBatch = 16;
    SeqId sid_batch[kBatch];
    int batch_count = 0;

    // Phase 3a: decode varints in batches via SeqIdDecoder::next_batch() so
    // the inner loop can amortize the per-call dispatch overhead and keep the
    // input bytes hot in cache. Per-element coverscore / OID-filter skip is
    // applied scalar-style on the small fixed-size buffer.
    constexpr int kDecBatch = SeqIdDecoder::kMaxBatch;
    static_assert(kDecBatch >= 1, "kDecBatch must be positive");
    SeqId   raw_sids[kDecBatch];
    uint8_t was_new[kDecBatch];

    for (size_t qi = 0; qi < n; qi++) {
        auto q_pos = static_cast<PosT>(positions[qi]);
        auto kmer_idx = kmers[qi];
        // Phase 7d hot-path: one EF access_pair instead of two access calls.
        uint64_t off, byte_len;
        kix.posting_list_range(kmer_idx, off, byte_len);
        if (byte_len == 0) continue;

        SeqIdDecoder decoder(posting_file + off, posting_file + off + byte_len);
        while (decoder.has_more()) {
            int n_dec = decoder.next_batch(raw_sids, was_new, kDecBatch);
            if (n_dec == 0) break;
            for (int i = 0; i < n_dec; i++) {
                SeqId sid = raw_sids[i];
                if (use_coverscore && !was_new[i]) continue;
                if (!filter.pass(sid)) continue;
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

    std::vector<Stage1Candidate> candidates;
    for (uint32_t sid : buf.dirty) {
        if (scores[sid] >= config.min_stage1_score) {
            candidates.push_back({sid, static_cast<uint32_t>(scores[sid])});
        }
    }

    buf.clear_dirty_typed<Tier>();

    if (config.stage1_topn == 0) return candidates;

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

// Explicit template instantiations
template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);
template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);

} // namespace ikafssn
