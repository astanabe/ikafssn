#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace ikafssn {

// Posting buffer entry used by index_builder before delta-compressing into the
// .kix / .kpx streams. The struct layout is shared between index_builder.cpp
// and parallel_sort_dispatch.cpp.
struct TempEntry {
    uint32_t kmer_value;
    uint32_t seq_id;
    uint32_t pos;
};
static_assert(sizeof(TempEntry) == 12, "TempEntry must be 12 bytes");

// Sort `buf` by (kmer_value, seq_id, pos). Equivalent to the previous
// `tbb::parallel_sort + 3-key comparator` but on x86_64 with AVX2+ uses
// Intel x86-simd-sort's `keyvalue_qsort<uint64_t, uint32_t>` to sort the
// (kmer_value << 32 | seq_id) keys against `pos` values, followed by a local
// pos-only sort within each (kmer, seq_id) block. Output is bit-exact with
// the TBB path because TempEntry tuples are unique per index build.
void parallel_sort_temp_entries(std::vector<TempEntry>& buf);

// Returns true when `parallel_sort_temp_entries` may dispatch to the SIMD
// path on this build + runtime configuration. Note: `parallel_sort_temp_entries`
// also gates on a per-call size threshold (large partitions go through TBB
// regardless), so a true return here does not guarantee the SIMD path runs
// for any specific call.
[[nodiscard]] bool parallel_sort_simd_active();

// Per-entry overhead (bytes) the index_builder uses to size partitions
// against `memory_limit`. The SIMD path needs an extra uint64_t key +
// uint32_t val array (12 extra bytes), but it only activates for small
// partitions (n <= 128K), so the dominant cost is always the TempEntry
// itself. Returning sizeof(TempEntry) keeps partition sizing simple and the
// SIMD-path scratch arrays comfortably below memory_limit on any realistic
// configuration.
[[nodiscard]] inline std::size_t parallel_sort_entry_overhead() {
    return sizeof(TempEntry);
}

} // namespace ikafssn
