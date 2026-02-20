#include "search/diagonal_filter.hpp"

#include <unordered_map>

namespace ikafssn {

std::vector<Hit> diagonal_filter(const std::vector<Hit>& hits,
                                 uint32_t min_diag_hits) {
    if (min_diag_hits <= 1) return hits; // no filtering needed

    // Count hits per diagonal
    std::unordered_map<int32_t, uint32_t> diag_counts;
    diag_counts.reserve(hits.size());
    for (const auto& hit : hits) {
        int32_t diag = static_cast<int32_t>(hit.s_pos) - static_cast<int32_t>(hit.q_pos);
        diag_counts[diag]++;
    }

    // Keep hits on diagonals with enough counts
    std::vector<Hit> filtered;
    filtered.reserve(hits.size());
    for (const auto& hit : hits) {
        int32_t diag = static_cast<int32_t>(hit.s_pos) - static_cast<int32_t>(hit.q_pos);
        if (diag_counts[diag] >= min_diag_hits) {
            filtered.push_back(hit);
        }
    }

    return filtered;
}

} // namespace ikafssn
