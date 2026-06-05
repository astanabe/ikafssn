#include "search/decode_cache.hpp"

#include "index/kix_reader.hpp"
#include "index/pfd_codec.hpp"

#include <algorithm>
#include <atomic>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace ikafssn {

template <typename KmerInt>
bool DecodedKmerCache::fill(const KixReader& kix,
                            const std::vector<KmerInt>& unique,
                            tbb::task_arena& arena, size_t budget) {
    const size_t n = unique.size();
    unique_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        unique_[i] = static_cast<uint32_t>(unique[i]);
    }
    offsets_.assign(n + 1, 0);

    const uint8_t* posting_file = kix.posting_file();

    // Pass 1 (serial, cheap): read each k-mer's distinct_count from its
    // posting list header and prefix-sum into offsets_.  The volume's WILLNEED
    // is already applied, so these header reads ride the prefetch.
    size_t total = 0;
    for (size_t i = 0; i < n; ++i) {
        offsets_[i] = total;
        uint64_t off, len;
        kix.posting_list_range(unique_[i], off, len);
        total += pfd::posting_count(posting_file + off, len);
    }
    offsets_[n] = total;

    // Reject an over-budget volume before decoding (no wasted Pass 2); the
    // caller runs it uncached.
    if (total > budget / sizeof(uint32_t)) {
        release();
        return false;
    }

    storage_.resize(total);

    // Pass 2 (parallel): decode each k-mer's posting list into its disjoint
    // arena slot.  Slots never overlap, so the writes are race-free.
    std::atomic<bool> corrupt{false};
    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, n),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    const size_t slot = offsets_[i];
                    const size_t cnt = offsets_[i + 1] - slot;
                    if (cnt == 0) continue;
                    uint64_t off, len;
                    kix.posting_list_range(unique_[i], off, len);
                    size_t got = pfd::decode_kix_into(posting_file + off, len,
                                                      storage_.data() + slot, cnt);
                    if (got != cnt) corrupt.store(true, std::memory_order_relaxed);
                }
            });
    });

    if (corrupt.load(std::memory_order_relaxed)) {
        release();
        return false;
    }
    return true;
}

void DecodedKmerCache::release() {
    std::vector<uint32_t>().swap(storage_);
    std::vector<uint32_t>().swap(unique_);
    std::vector<size_t>().swap(offsets_);
}

DecodedKmerCache::Lookup DecodedKmerCache::lookup(uint32_t kmer) const {
    auto it = std::lower_bound(unique_.begin(), unique_.end(), kmer);
    if (it == unique_.end() || *it != kmer) {
        return {nullptr, 0};
    }
    const size_t i = static_cast<size_t>(it - unique_.begin());
    return {storage_.data() + offsets_[i],
            static_cast<uint32_t>(offsets_[i + 1] - offsets_[i])};
}

template bool DecodedKmerCache::fill<uint16_t>(
    const KixReader&, const std::vector<uint16_t>&, tbb::task_arena&, size_t);
template bool DecodedKmerCache::fill<uint32_t>(
    const KixReader&, const std::vector<uint32_t>&, tbb::task_arena&, size_t);

}  // namespace ikafssn
