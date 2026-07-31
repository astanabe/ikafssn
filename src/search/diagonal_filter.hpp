#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

// Keep only the hits whose diagonal (s_pos - q_pos, signed) carries at least
// min_nhit_diag hits; <= 1 is a no-op.  Filters in place and preserves input
// order, so a (q_pos, s_pos)-ordered input stays ordered.  `diag_scratch` is a
// caller-owned working buffer, overwritten on entry.
void diagonal_filter(std::vector<Hit>& hits,
                     uint32_t min_nhit_diag,
                     std::vector<int32_t>& diag_scratch);

} // namespace ikafssn
