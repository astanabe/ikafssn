#pragma once

#include "core/spaced_seed.hpp"

#include <cstddef>
#include <cstdint>

namespace ikafssn {

// Batch extraction of spaced-seed k-mers from a sliding accumulator stream.
//
// Both APIs are bit-exact with extract_kmer_ct() / extract_for_mask() applied
// element-wise; tier dispatch happens internally via current_simd_cap(); the
// scalar fallback covers small n and the SIMD-disabled build. (Phase 4b-1, M-3.)

// NTTP-specialized batch extraction. The extractor is known at compile time,
// allowing the run-application loop to unroll into immediate-shift SIMD ops.
template <SpacedSeedExtractor Ext, typename KmerInt>
void extract_kmer_ct_batch(const std::uint64_t* accums, std::size_t n,
                           KmerInt* out) noexcept;

// Mask-dispatched batch extraction. Mirrors extract_for_mask() shape; the
// switch on `mask` selects the matching NTTP specialization.
template <typename KmerInt>
void extract_for_mask_batch(const std::uint64_t* accums, std::size_t n,
                            std::uint32_t mask, KmerInt* out) noexcept;

extern template void extract_for_mask_batch<std::uint16_t>(
    const std::uint64_t*, std::size_t, std::uint32_t, std::uint16_t*) noexcept;
extern template void extract_for_mask_batch<std::uint32_t>(
    const std::uint64_t*, std::size_t, std::uint32_t, std::uint32_t*) noexcept;

} // namespace ikafssn
