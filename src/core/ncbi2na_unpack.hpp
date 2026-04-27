#pragma once

#include <cstddef>
#include <cstdint>

namespace ikafssn {

// Bulk-unpack ncbi2na (MSB-first packed, 4 bases/byte) into a flat array of
// 2-bit base codes (one base per output byte, value in {0, 1, 2, 3}).
//
// Layout:
//   packed[start_base_pos / 4 .. (start_base_pos + count + 3) / 4)  (input)
//   out_bases[0 .. count)                                           (output)
//
// `start_base_pos` may have any 4-residue (the head misaligned region is
// processed by the scalar fallback inside the dispatcher). The internal SIMD
// kernels are entered only with start_base_pos % 4 == 0.
//
// Bit-exact equivalent to the scalar reference:
//   uint8_t b = packed[(start_base_pos + i) / 4];
//   out[i] = (b >> (6 - 2 * ((start_base_pos + i) % 4))) & 0x03;
void unpack_ncbi2na_chunk(const char*    packed,
                          std::uint32_t  start_base_pos,
                          std::uint32_t  count,
                          std::uint8_t*  out_bases) noexcept;

// Suggested per-call chunk size for callers prefetching SIMD-friendly windows.
inline constexpr std::uint32_t kNcbi2naUnpackChunkSize = 64;

} // namespace ikafssn
