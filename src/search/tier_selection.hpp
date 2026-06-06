#pragma once

// Tier selection helper shared by ikafssnsearch and ikafssnserver.
// Both pick `Stage1Tier` from the maximum per-strand k-mer position count
// observed across all preprocessed queries; the ladder is computed by
// `select_tier()` in stage1_filter.hpp.

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "search/stage1_filter.hpp"

namespace ikafssn {

// Walk a container of preprocessed-query wrappers (each exposing a `.qdata`
// member that is a QueryKmerData<KmerInt>) and update the running max of
// `max(fwd_positions.size(), rc_positions.size())` across all entries.
//
// Templated on the container so callers can pass either ikafssnsearch's
// `vector<PreprocessedQuery<KmerInt>>` or ikafssnserver's
// `vector<PreprocessedQuery{16,32}>` without repeating the body.
template <typename Container>
inline uint32_t accumulate_max_kmer_positions(uint32_t cur_max,
                                              const Container& c) {
    for (const auto& elem : c) {
        const auto& d = elem.qdata;
        cur_max = std::max(cur_max,
            static_cast<uint32_t>(std::max(d.fwd_positions.size(),
                                           d.rc_positions.size())));
    }
    return cur_max;
}

} // namespace ikafssn
