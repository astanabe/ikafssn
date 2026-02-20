#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

// Apply diagonal filter: keep only hits on diagonals that have
// at least min_diag_hits occurrences.
// diagonal = s_pos - q_pos (can be negative, stored as int32_t).
std::vector<Hit> diagonal_filter(const std::vector<Hit>& hits,
                                 uint32_t min_diag_hits);

} // namespace ikafssn
