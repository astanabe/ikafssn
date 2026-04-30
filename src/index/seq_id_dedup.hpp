#pragma once

// Phase 5i — SIMD dedup kernel for the .kix v7 build pipeline.
//
// Given a sorted seq_id array (with duplicates from intra-sequence k-mer
// repetitions), produce:
//   - distinct_sid_out[0..d-1]: unique seq_ids in input order
//   - occ_count_out[0..d-1]:    per-seq_id occurrence count (u8; <= 255)
//   - return value:             d = number of distinct seq_ids
//
// The caller is expected to have ensured occurrence runs are <= 255 by
// capping `freq_threshold_part` at 255 (CLI validation, see
// src/ikafssnindex/main.cpp).  The kernel does *not* enforce that
// constraint — it simply accumulates whatever count it sees and would
// silently truncate above 255.  Index-builder users always pass a
// pre-bounded threshold, so the assertion is enforced there.
//
// Implementation: per-ISA OBJECT libraries selected at runtime via the
// same mechanism as src/index/pfd_codec.cpp.  Currently a scalar fall-
// back lives in each tier file; future work can add gather/compress-
// store specialisations.

#include <cstdint>
#include <cstddef>

namespace ikafssn::seq_id_dedup {

// Run the dedup kernel.
//
//   sid              sorted seq_id input, length n
//   n                input length
//   distinct_sid_out output buffer, must hold at least n elements
//   occ_count_out    output buffer, must hold at least n bytes
//
// Returns the number of distinct seq_ids written to the out buffers.
std::uint32_t dedup_seq_ids(const std::uint32_t* sid, std::uint32_t n,
                            std::uint32_t* distinct_sid_out,
                            std::uint8_t*  occ_count_out);

// Diagnostic name of the active ISA tier ("sse42" / "avx2" /
// "avx512bw" / "avx512vbmi2" on x86_64; "neon" on aarch64).
const char* active_tier_name();

} // namespace ikafssn::seq_id_dedup
