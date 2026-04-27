#pragma once

#include <cstddef>
#include <cstdint>

namespace ikafssn {

// Decode up to max_count consecutive LEB128 uint32_t varints from [in, in_end).
// Writes decoded values into out_values[] and the consumed byte count into
// *out_consumed (so the caller can advance its pointer).
//
// Returns the number of varints actually decoded (0..max_count). Stops early
// when reaching in_end mid-varint or after consuming all available bytes.
//
// Pure decode: no delta arithmetic, no sequence-boundary detection. The caller
// (SeqIdDecoder / PosDecoder) is responsible for those.
//
// Tier dispatch happens internally at the function-entry boundary. The caller
// pays one runtime branch on current_simd_cap() per call.
//
// max_count must be <= kVarintMaxBatch (16) so per-tier kernels can keep
// fixed-size temporary buffers on the stack.
inline constexpr int kVarintMaxBatch = 16;

int varint_decode_batch(const std::uint8_t* in,
                        const std::uint8_t* in_end,
                        std::uint32_t* out_values,
                        std::size_t* out_consumed,
                        int max_count) noexcept;

} // namespace ikafssn
