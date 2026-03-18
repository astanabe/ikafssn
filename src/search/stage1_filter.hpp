#pragma once

#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KixReader;
class OidFilter;

// Tier selection for Stage1Buffer: controls entry size per sequence.
enum class Stage1Tier : uint8_t { T8 = 0, T16 = 1, T32 = 2 };

// AoS entry: score + last_scored_pos packed together.
template <Stage1Tier> struct Stage1Entry;

template <> struct Stage1Entry<Stage1Tier::T8> {
    uint8_t score;
    uint8_t last_pos;
};

template <> struct Stage1Entry<Stage1Tier::T16> {
    uint16_t score;
    uint16_t last_pos;
};

template <> struct Stage1Entry<Stage1Tier::T32> {
    uint32_t score;
    uint32_t last_pos;
};

static_assert(sizeof(Stage1Entry<Stage1Tier::T8>)  == 2, "T8 entry must be 2 bytes");
static_assert(sizeof(Stage1Entry<Stage1Tier::T16>) == 4, "T16 entry must be 4 bytes");
static_assert(sizeof(Stage1Entry<Stage1Tier::T32>) == 8, "T32 entry must be 8 bytes");

// Type-erased Stage1Buffer. Internally stores AoS entries at the selected tier.
struct Stage1Buffer {
    std::vector<uint8_t> data;       // raw storage for Stage1Entry<Tier>[]
    std::vector<uint32_t> dirty;     // dirty list of modified seq IDs
    Stage1Tier tier = Stage1Tier::T32;
    uint32_t capacity = 0;           // num_seqs capacity

    void ensure_capacity(uint32_t num_seqs) {
        if (capacity >= num_seqs) return;
        capacity = num_seqs;
        switch (tier) {
        case Stage1Tier::T8:
            data.resize(num_seqs * sizeof(Stage1Entry<Stage1Tier::T8>));
            break;
        case Stage1Tier::T16:
            data.resize(num_seqs * sizeof(Stage1Entry<Stage1Tier::T16>));
            break;
        case Stage1Tier::T32:
            data.resize(num_seqs * sizeof(Stage1Entry<Stage1Tier::T32>));
            break;
        }
        std::memset(data.data(), 0, data.size());
        // Set sentinel values for last_pos
        reset_all();
    }

    // Reset all entries to zero score and sentinel last_pos.
    void reset_all() {
        switch (tier) {
        case Stage1Tier::T8: {
            auto* e = reinterpret_cast<Stage1Entry<Stage1Tier::T8>*>(data.data());
            for (uint32_t i = 0; i < capacity; i++) {
                e[i].score = 0;
                e[i].last_pos = std::numeric_limits<uint8_t>::max();
            }
            break;
        }
        case Stage1Tier::T16: {
            auto* e = reinterpret_cast<Stage1Entry<Stage1Tier::T16>*>(data.data());
            for (uint32_t i = 0; i < capacity; i++) {
                e[i].score = 0;
                e[i].last_pos = std::numeric_limits<uint16_t>::max();
            }
            break;
        }
        case Stage1Tier::T32: {
            auto* e = reinterpret_cast<Stage1Entry<Stage1Tier::T32>*>(data.data());
            for (uint32_t i = 0; i < capacity; i++) {
                e[i].score = 0;
                e[i].last_pos = std::numeric_limits<uint32_t>::max();
            }
            break;
        }
        }
    }

    template <Stage1Tier Tier>
    void clear_dirty_typed() {
        using Entry = Stage1Entry<Tier>;
        using PosT = decltype(Entry::last_pos);
        auto* entries = reinterpret_cast<Entry*>(data.data());
        for (uint32_t idx : dirty) {
            entries[idx].score = 0;
            entries[idx].last_pos = std::numeric_limits<PosT>::max();
        }
        dirty.clear();
    }
};

// Determine the optimal tier based on max query k-mer position count and max position value.
inline Stage1Tier select_tier(uint32_t max_kmer_positions, uint32_t max_position_value) {
    uint32_t limit = std::max(max_kmer_positions, max_position_value);
    if (limit < 255) return Stage1Tier::T8;
    if (limit < 65535) return Stage1Tier::T16;
    return Stage1Tier::T32;
}

struct Stage1Candidate {
    SeqId id;
    uint32_t score;
};

struct Stage1Config {
    static constexpr uint32_t MAX_FREQ_DISABLED = UINT32_MAX;

    uint32_t max_freq = 0;
    uint32_t stage1_topn = 0;
    uint32_t min_stage1_score = 1;
    uint8_t  stage1_score_type = 1;
};

uint32_t compute_effective_max_freq(uint32_t config_max_freq,
                                    uint64_t total_postings,
                                    uint32_t table_size);

template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer* buf = nullptr);

extern template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);
extern template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer*);

} // namespace ikafssn
