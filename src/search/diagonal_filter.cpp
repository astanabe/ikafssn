#include "search/diagonal_filter.hpp"

#include <algorithm>

namespace ikafssn {

void diagonal_filter(std::vector<Hit>& hits,
                     uint32_t min_nhit_diag,
                     std::vector<int32_t>& diag_scratch) {
    if (min_nhit_diag <= 1) return;
    const size_t n = hits.size();
    if (n == 0) return;

    // Every hit's diagonal, sorted, so equal_range below counts occurrences.
    diag_scratch.clear();
    diag_scratch.reserve(n);
    for (const Hit& hit : hits) {
        diag_scratch.push_back(static_cast<int32_t>(hit.s_pos) -
                               static_cast<int32_t>(hit.q_pos));
    }
    std::sort(diag_scratch.begin(), diag_scratch.end());

    size_t w = 0;
    for (size_t r = 0; r < n; r++) {
        const int32_t diag = static_cast<int32_t>(hits[r].s_pos) -
                             static_cast<int32_t>(hits[r].q_pos);
        const auto range =
            std::equal_range(diag_scratch.begin(), diag_scratch.end(), diag);
        if (static_cast<uint32_t>(range.second - range.first) >= min_nhit_diag) {
            hits[w++] = hits[r];
        }
    }
    hits.resize(w);
}

} // namespace ikafssn
