#include "search/stage1_filter.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace ikafssn {

uint32_t compute_effective_max_freq(uint32_t config_max_freq,
                                    uint64_t total_postings,
                                    uint32_t table_size) {
    if (config_max_freq > 0) return config_max_freq;
    double mean = static_cast<double>(total_postings) /
                  static_cast<double>(table_size);
    uint32_t max_freq = static_cast<uint32_t>(mean * 10.0);
    if (max_freq < 1000) max_freq = 1000;
    if (max_freq > 100000) max_freq = 100000;
    return max_freq;
}

// Internal implementation with KmerInt + Tier template dispatch.
template <typename KmerInt, Stage1Tier Tier>
static std::vector<Stage1Candidate> stage1_filter_impl(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer* buf) {

    using Entry = Stage1Entry<Tier>;
    using ScoreT = decltype(Entry::score);
    using PosT = decltype(Entry::last_pos);
    constexpr PosT SENTINEL = std::numeric_limits<PosT>::max();

    uint32_t num_seqs = kix.num_sequences();
    if (num_seqs == 0 || n == 0) return {};

    const uint8_t* posting_data = kix.posting_data();
    const bool use_coverscore = (config.stage1_score_type == 1);

    if (buf) {
        buf->ensure_capacity(num_seqs);
        auto* entries = reinterpret_cast<Entry*>(buf->data.data());

        for (size_t qi = 0; qi < n; qi++) {
            auto q_pos = static_cast<PosT>(positions[qi]);
            auto kmer_idx = kmers[qi];
            auto off = kix.posting_offset(kmer_idx);
            auto end_off = kix.posting_offset(kmer_idx + 1);
            if (off == end_off) continue;

            SeqIdDecoder decoder(posting_data + off, posting_data + end_off);
            while (decoder.has_more()) {
                SeqId sid = decoder.next();
                if (use_coverscore && !decoder.was_new_seq()) continue;
                if (!filter.pass(sid)) continue;
                if (entries[sid].score == 0) buf->dirty.push_back(sid);
                if (entries[sid].last_pos != q_pos) {
                    entries[sid].score++;
                    entries[sid].last_pos = q_pos;
                }
            }
        }

        std::vector<Stage1Candidate> candidates;
        for (uint32_t sid : buf->dirty) {
            if (entries[sid].score >= config.min_stage1_score) {
                candidates.push_back({sid, static_cast<uint32_t>(entries[sid].score)});
            }
        }

        buf->clear_dirty_typed<Tier>();

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

    // Fallback: allocate local T32 buffer (always safe)
    std::vector<Stage1Entry<Stage1Tier::T32>> local_entries(num_seqs);
    for (auto& e : local_entries) {
        e.score = 0;
        e.last_pos = UINT32_MAX;
    }

    for (size_t qi = 0; qi < n; qi++) {
        uint32_t q_pos = positions[qi];
        auto kmer_idx = kmers[qi];
        auto off = kix.posting_offset(kmer_idx);
        auto end_off = kix.posting_offset(kmer_idx + 1);
        if (off == end_off) continue;

        SeqIdDecoder decoder(posting_data + off, posting_data + end_off);
        while (decoder.has_more()) {
            SeqId sid = decoder.next();
            if (use_coverscore && !decoder.was_new_seq()) continue;
            if (!filter.pass(sid)) continue;
            if (local_entries[sid].last_pos != q_pos) {
                local_entries[sid].score++;
                local_entries[sid].last_pos = q_pos;
            }
        }
    }

    std::vector<Stage1Candidate> candidates;
    for (uint32_t oid = 0; oid < num_seqs; oid++) {
        if (local_entries[oid].score >= config.min_stage1_score) {
            candidates.push_back({oid, local_entries[oid].score});
        }
    }

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

// Public dispatch: selects tier from buffer (or uses T32 fallback).
template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer* buf) {

    Stage1Tier tier = buf ? buf->tier : Stage1Tier::T32;
    switch (tier) {
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

// Explicit template instantiations (2 KmerInt types × dispatch internally)
template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);
template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);

} // namespace ikafssn
