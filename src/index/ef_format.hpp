#pragma once

// Phase 7 — Elias-Fano (EF) succinct dictionary header.
//
// On-disk layout written by encode_dictionary_ef() and consumed by
// EFDictionary::open():
//
//   [EFHeader  : 32 B]
//   [lower bits: ((D * l + 63) / 64) * 8 bytes]
//   [upper bits: ((upper_bits + 63) / 64) * 8 bytes]
//   [select samples: select_count * 8 bytes]
//
// The dictionary stores D non-strictly monotonic offsets in
// [0, U_raw].  The encoder applies the canonical strict-monotonicity
// transform ``ef[i] = offsets[i] + i`` (strictly increasing) before EF
// encoding, so the EF universe is U_ef = U_raw + D.  ``access(i)``
// reverses the transform by subtracting i.
//
// l = floor(log2(U_ef / D))   (capped at 63; 0 when D == 0 or U_ef <= D)
//
// Lower bits hold the low l bits of every ef-value, packed contiguously
// little-endian in u64 words.  Upper bits hold the unary-coded gaps of
// the high (b - l) bits where b = ceil(log2(U_ef)): for each i we set
// the bit at position (ef[i] >> l) + i, so the i-th set bit appears at
// position (high_i + i).  Total upper-bit length = D + max_high + 1.
//
// Select samples are bit positions of the i-th set bit for i = 0,
// kSelectStep, 2*kSelectStep, ... (one u64 per sample).  ``access(i)``
// looks up sample (i / kSelectStep), then word-scan-forward to locate
// the (i % kSelectStep)-th additional set bit.

#include <cstdint>

namespace ikafssn::ef {

// Sample one bit position per kSelectStep set bits.  Fixed at compile
// time so encoder and decoder agree without storing it in the header.
inline constexpr uint32_t kSelectStep = 64;

inline constexpr char EF_MAGIC[4] = {'E', 'F', '0', '1'};

#pragma pack(push, 1)
struct EFHeader {
    char     magic[4];        // 0x00: "EF01"
    uint8_t  l;               // 0x04: lower bits per entry (0..63)
    uint8_t  reserved0[3];    // 0x05
    uint64_t D;               // 0x08: number of entries
    uint64_t U;               // 0x10: U_ef (strict-monotonic upper bound)
    uint32_t upper_bits;      // 0x18: total upper-bit length
    uint32_t select_count;    // 0x1C: number of select samples
};
#pragma pack(pop)

static_assert(sizeof(EFHeader) == 32, "EFHeader must be 32 bytes");

} // namespace ikafssn::ef
