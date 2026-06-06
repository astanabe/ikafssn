#pragma once

#include <cstdint>
#include <type_traits>
#include <vector>

#include "core/types.hpp"
#include "util/aligned_buffer.hpp"

namespace ikafssn {

class KixReader;
class OidFilter;

// Score/position width selection for Stage1Buffer: controls entry size per sequence.
enum class Stage1Width : uint8_t { T8 = 0, T16 = 1, T32 = 2 };

// Determine the position type (PosT) and score type (ScoreT) for a width.
// last_pos and score share the same bit width.
template <Stage1Width Width> struct Stage1WidthTraits;
template <> struct Stage1WidthTraits<Stage1Width::T8>  { using PosT = uint8_t;  using ScoreT = uint8_t;  };
template <> struct Stage1WidthTraits<Stage1Width::T16> { using PosT = uint16_t; using ScoreT = uint16_t; };
template <> struct Stage1WidthTraits<Stage1Width::T32> { using PosT = uint32_t; using ScoreT = uint32_t; };

// SoA Stage 1 buffer: score and last_pos kept in separate cache-line aligned
// heap buffers (64-byte alignment for AVX-512 gather/scatter friendly access).
// The underlying byte storage is reinterpret_cast to the width-specific element
// type via the score_ptr / last_pos_ptr helpers below.
struct Stage1Buffer {
    AlignedBuffer<uint8_t> score_data;     // width-erased: bytes for score[]
    AlignedBuffer<uint8_t> last_pos_data;  // width-erased: bytes for last_pos[]
    std::vector<uint32_t>  dirty;          // sentinel-marked seq IDs touched in this query batch
    Stage1Width width = Stage1Width::T32;
    uint32_t   capacity = 0;               // num_seqs capacity

    Stage1Buffer() = default;
    Stage1Buffer(Stage1Buffer&&) noexcept = default;
    Stage1Buffer& operator=(Stage1Buffer&&) noexcept = default;
    Stage1Buffer(const Stage1Buffer&) = delete;
    Stage1Buffer& operator=(const Stage1Buffer&) = delete;

    void ensure_capacity(uint32_t num_seqs);
    void reset_all();

    template <Stage1Width Width>
    void clear_dirty_typed();
};

// SoA accessors: header-inline so the hot path is fully visible to the optimizer.
template <Stage1Width Width>
inline typename Stage1WidthTraits<Width>::ScoreT* score_ptr(Stage1Buffer& buf) {
    using ScoreT = typename Stage1WidthTraits<Width>::ScoreT;
    return reinterpret_cast<ScoreT*>(buf.score_data.data());
}

template <Stage1Width Width>
inline typename Stage1WidthTraits<Width>::PosT* last_pos_ptr(Stage1Buffer& buf) {
    using PosT = typename Stage1WidthTraits<Width>::PosT;
    return reinterpret_cast<PosT*>(buf.last_pos_data.data());
}

extern template void Stage1Buffer::clear_dirty_typed<Stage1Width::T8>();
extern template void Stage1Buffer::clear_dirty_typed<Stage1Width::T16>();
extern template void Stage1Buffer::clear_dirty_typed<Stage1Width::T32>();

// Determine the optimal width based on max query k-mer position count and max position value.
inline Stage1Width select_width(uint32_t max_kmer_positions, uint32_t max_position_value) {
    uint32_t limit = std::max(max_kmer_positions, max_position_value);
    if (limit < 255) return Stage1Width::T8;
    if (limit < 65535) return Stage1Width::T16;
    return Stage1Width::T32;
}

struct Stage1Candidate {
    SeqId id;
    uint32_t score;
};

struct Stage1Config {
    uint32_t stage1_topn = 0;
    uint32_t min_stage1_score = 1;
    uint8_t  stage1_score_type = 1;
};

template <typename KmerInt>
std::vector<Stage1Candidate> stage1_filter(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf);

extern template std::vector<Stage1Candidate> stage1_filter<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);
extern template std::vector<Stage1Candidate> stage1_filter<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&);

// Cross-template accumulation ("both" mode): accumulate score / last_pos
// updates into buf without clearing dirty. Multiple calls share the same
// buf to merge scores from independent (kix, kmers) streams; per-position
// dedup (`last_pos[sid] != q_pos`) carries across calls so that two
// templates hitting the same query position contribute exactly +1 to the
// score.
//
// Use config.stage1_score_type only; stage1_topn is applied later by
// stage1_filter_finish().
//
// `remaining` engages the Stage 1 score cutoff: it is an upper bound on the
// additional score any sequence can still gain from the query positions left to
// process (the number of remaining distinct positions).  When `remaining > 0`
// the flush kernel prunes sequences whose `score + remaining < config.min_stage1_score`
// (the resolved final threshold), since they can never reach it.  Pass
// `remaining == 0` to disable the cutoff and preserve the original behaviour
// (used for uncached volumes and callers that do not order positions by cost).
template <typename KmerInt>
void stage1_filter_accumulate(
    const uint32_t* positions, const KmerInt* kmers, size_t n,
    const KixReader& kix,
    const OidFilter& filter,
    const Stage1Config& config,
    Stage1Buffer& buf,
    uint32_t remaining = 0);

extern template void stage1_filter_accumulate<uint16_t>(
    const uint32_t*, const uint16_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&, uint32_t);
extern template void stage1_filter_accumulate<uint32_t>(
    const uint32_t*, const uint32_t*, size_t,
    const KixReader&, const OidFilter&, const Stage1Config&,
    Stage1Buffer&, uint32_t);

// Collect candidates from buf populated by one or more
// stage1_filter_accumulate calls, then clear dirty so the buf is reusable.
// `min_stage1_score` and `stage1_topn` from config are applied here.
std::vector<Stage1Candidate> stage1_filter_finish(
    Stage1Buffer& buf, const Stage1Config& config);

} // namespace ikafssn
