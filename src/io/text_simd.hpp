#pragma once

#include <cstddef>
#include <cstdint>

namespace ikafssn {

// In-place ASCII bulk to-upper. Maps 'a'..'z' to 'A'..'Z' (-0x20); leaves all
// other byte values untouched. Locale-independent (does not call std::locale).
//
// Internally dispatches on the runtime SIMD capability (current_simd_cap()).
// Inputs smaller than the per-tier kSimdMinBytes_* thresholds, and the tail
// remainder of any SIMD chunk loop, fall back to a scalar pass.
//
// Bit-exact equivalent to the scalar reference implementation:
//   for (size_t i = 0; i < n; ++i)
//       if (p[i] >= 'a' && p[i] <= 'z') p[i] -= 0x20;
void toupper_inplace_ascii(std::uint8_t* p, std::size_t n) noexcept;

} // namespace ikafssn
