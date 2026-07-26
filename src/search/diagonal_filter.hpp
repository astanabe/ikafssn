#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

// Apply diagonal filter in place: keep only hits on diagonals that have
// at least min_nhit_diag occurrences.  Input order is preserved, so a
// (q_pos, s_pos)-ordered input stays ordered.  min_nhit_diag <= 1 is a no-op.
// diagonal = s_pos - q_pos (can be negative, stored as int32_t).
// `diag_scratch` is a caller-owned working buffer; its contents on entry are
// ignored and it grows to the largest hit count seen so far.
void diagonal_filter(std::vector<Hit>& hits,
                     uint32_t min_nhit_diag,
                     std::vector<int32_t>& diag_scratch);

} // namespace ikafssn
