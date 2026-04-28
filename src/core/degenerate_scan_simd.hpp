#pragma once

#include <cstddef>

namespace ikafssn {

// True if [data, data+n) contains any IUPAC degenerate base
// (R/Y/S/W/K/M/B/D/H/V/N, upper or lower case).
//
// Bit-exact with the scalar reference: degenerate_base_table()[c] != 0 for
// any c in the buffer. Internally tier-dispatched on current_simd_cap();
// scalar fallback covers small n and the SIMD-disabled build.
bool has_degenerate_base(const char* data, std::size_t n) noexcept;

} // namespace ikafssn
