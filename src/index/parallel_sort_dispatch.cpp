#include "index/parallel_sort_dispatch.hpp"
#include "util/simd_dispatch.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include <tbb/parallel_sort.h>

#if IKAFSSN_HAS_X86_SIMD_SORT
// Public dispatched API. The library's per-arch implementations
// (x86simdsort-skx.cpp, x86simdsort-avx2.cpp) are linked in by CMake; the
// runtime selector inside libxss uses __builtin_cpu_supports. We additionally
// gate at our SimdCap so that IKAFSSN_FORCE_SIMD=scalar/sse42 routes to TBB.
#include "x86simdsort.h"
#endif

namespace ikafssn {

namespace {
// Size threshold above which the SIMD path is bypassed in favor of
// `tbb::parallel_sort`. x86-simd-sort runs single-threaded; once n is large
// enough that TBB's split overhead is amortized across cores, the multi-
// thread parallel sort beats single-thread SIMD. Empirically (Zen5 9950X,
// 32 logical cores, see bench/bench_parallel_sort) the crossover sits below
// 16K entries — at n=1024 the SIMD path is ~12% faster, at n=16K it is
// ~18% slower. 8192 keeps the SIMD path exercised by tests for the regime
// where it wins, without regressing real ikafssnindex partitions which are
// typically much larger.
constexpr std::size_t kSimdSortMaxN = 8192;
} // namespace

bool parallel_sort_simd_active() {
#if IKAFSSN_HAS_X86_SIMD_SORT
    return current_simd_cap() >= SimdCap::AVX2;
#else
    return false;
#endif
}

namespace {

void tbb_sort_fallback(std::vector<TempEntry>& buf) {
    tbb::parallel_sort(buf.begin(), buf.end(),
        [](const TempEntry& a, const TempEntry& b) {
            if (a.kmer_value != b.kmer_value) return a.kmer_value < b.kmer_value;
            if (a.seq_id     != b.seq_id)     return a.seq_id < b.seq_id;
            return a.pos < b.pos;
        });
}

#if IKAFSSN_HAS_X86_SIMD_SORT
// Key-value SIMD sort:
//   1. key[i] = (uint64_t(kmer) << 32) | seq_id
//   2. val[i] = pos
//   3. x86simdsort::keyvalue_qsort<uint64_t, uint32_t>(key, val, n)
//   4. For each contiguous block of identical keys (same kmer + seq_id),
//      std::sort the val slice to recover (kmer, seq_id, pos) total order.
// keyvalue_qsort is unstable for duplicate keys, so step 4 is required to
// make the (pos) tail match the TBB path.
void simd_sort_keyvalue(std::vector<TempEntry>& buf) {
    const std::size_t n = buf.size();
    std::vector<uint64_t> key(n);
    std::vector<uint32_t> val(n);
    for (std::size_t i = 0; i < n; ++i) {
        key[i] = (static_cast<uint64_t>(buf[i].kmer_value) << 32)
               |  static_cast<uint64_t>(buf[i].seq_id);
        val[i] = buf[i].pos;
    }

    x86simdsort::keyvalue_qsort<uint64_t, uint32_t>(key.data(), val.data(), n);

    // Local pos sort within each (kmer, seq_id) block. Each block is small
    // (typically a few entries — multiple positions of the same k-mer in one
    // sequence) so std::sort is cheap; the parallelism win came from step 3.
    std::size_t i = 0;
    while (i < n) {
        std::size_t j = i + 1;
        while (j < n && key[j] == key[i]) ++j;
        if (j - i > 1) std::sort(val.data() + i, val.data() + j);
        i = j;
    }

    // Reassemble TempEntry in place.
    for (std::size_t i = 0; i < n; ++i) {
        buf[i].kmer_value = static_cast<uint32_t>(key[i] >> 32);
        buf[i].seq_id     = static_cast<uint32_t>(key[i] & 0xFFFFFFFFu);
        buf[i].pos        = val[i];
    }
}
#endif // IKAFSSN_HAS_X86_SIMD_SORT

} // namespace

void parallel_sort_temp_entries(std::vector<TempEntry>& buf) {
    if (buf.size() < 2) return;

#if IKAFSSN_HAS_X86_SIMD_SORT
    if (parallel_sort_simd_active() && buf.size() <= kSimdSortMaxN) {
        simd_sort_keyvalue(buf);
        return;
    }
#endif

    tbb_sort_fallback(buf);
}

} // namespace ikafssn
