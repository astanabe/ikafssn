#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <tbb/task_arena.h>

namespace ikafssn {

class KixReader;

// Per-volume Stage 1 decode cache.
//
// Stage 1 parallelises over (query x volume x strand) jobs and re-decodes
// the same .kix posting list once per query that probes a given k-mer.  For
// a batch of similar queries that share conserved k-mers this re-decodes the
// same posting list up to (number of queries) times inside one volume.
//
// DecodedKmerCache materialises every unique k-mer's posting list for one
// volume exactly once into a single arena and lets all queries read it.  The
// decoded seq_id arrays are byte-identical to open_stream_kix's output, so a
// cached probe produces numerically identical Stage 1 scores.
//
// fill() sizes the cache from the posting list headers, rejects an
// over-budget or corrupt volume, then decodes every posting list into the
// arena in parallel.  The orchestrator applies the volume's WILLNEED before
// fill() so both the header reads and the decode ride the prefetch, and
// charges the decoded byte size against the same (soft) budget as the
// WILLNEED page cache.  Because -memory_limit is a soft limit, the batch may
// overshoot it by up to one volume's worth; a volume whose decode alone
// exceeds the budget, or any corrupt posting list, leaves the cache empty so
// Stage 1 falls back to the on-the-fly decode path for that volume.
class DecodedKmerCache {
public:
    struct Lookup {
        const uint32_t* sids;  // ascending distinct seq_ids; nullptr => miss
        uint32_t        count; // number of seq_ids (0 for a present-but-empty
                               // posting list)
    };

    // Build the cache for `unique` (sorted ascending) over the open volume
    // `kix`: read the per-k-mer headers to lay out the arena, reject the
    // volume if the decoded size exceeds `budget` (no decode performed) or any
    // posting list is corrupt, otherwise decode every posting list into the
    // arena via `arena` (disjoint slots => race-free).  Returns true on a
    // usable cache (size then available via decoded_bytes()); false leaves the
    // cache empty for the caller to run the volume uncached.
    template <typename KmerInt>
    bool fill(const KixReader& kix, const std::vector<KmerInt>& unique,
              tbb::task_arena& arena, size_t budget);

    // Decoded arena byte size (valid after a successful fill()).
    size_t decoded_bytes() const { return storage_.size() * sizeof(uint32_t); }

    // Free the decoded arena and layout state, dropping resident heap.
    void release();

    // Resolve a probed k-mer to its decoded seq_id array.  Returns
    // {nullptr, 0} for a k-mer absent from the filled unique set (the caller
    // then decodes on the fly); otherwise {ptr, count} where count may be 0.
    Lookup lookup(uint32_t kmer) const;

private:
    std::vector<uint32_t> unique_;   // promoted sorted unique k-mer values
    std::vector<uint32_t> storage_;  // arena of decoded seq_ids
    std::vector<size_t>   offsets_;  // size unique_.size()+1; prefix sums
};

}  // namespace ikafssn
