#pragma once

#include <cstddef>
#include <cstdint>

namespace ikafssn {

// Batch reverse-complement: out[i] = kmer_revcomp(in[i], k) for i in [0, n).
//
// Bit-exact with the scalar kmer_revcomp() inline template. Tier dispatch
// happens internally; scalar fallback covers small n and the SIMD-disabled
// build. in/out may alias (in == out is valid; the kernel reads each element
// before writing it back).
//
// Specializations are provided for KmerInt = uint16_t (k <= 8) and uint32_t
// (k >= 9). Other widths are not used by the project.
template <typename KmerInt>
void kmer_revcomp_batch(const KmerInt* in, KmerInt* out, std::size_t n,
                        int k) noexcept;

extern template void kmer_revcomp_batch<std::uint16_t>(
    const std::uint16_t*, std::uint16_t*, std::size_t, int) noexcept;
extern template void kmer_revcomp_batch<std::uint32_t>(
    const std::uint32_t*, std::uint32_t*, std::size_t, int) noexcept;

} // namespace ikafssn
