#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include "io/mmap_file.hpp"
#include "index/ef_codec.hpp"
#include "index/kix_format.hpp"

namespace ikafssn {

class DecodedKmerCache;

class KixReader {
public:
    bool open(const std::string& path);
    void close();

    const KixHeader& header() const { return *header_; }
    int k() const { return header_->k; }
    uint8_t kmer_type() const { return header_->kmer_type; }
    uint32_t num_sequences() const { return header_->num_sequences; }
    uint64_t total_distinct_postings() const { return header_->total_distinct_postings; }
    uint8_t t() const { return header_->t; }
    uint8_t template_type() const { return header_->template_type; }
    uint32_t table_size() const { return table_size_; }

    // Raw pointer to the start of the ID posting file
    const uint8_t* posting_file() const { return posting_file_; }
    size_t posting_file_size() const { return posting_file_size_; }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

    // Byte size of the whole mapping (dictionary head + posting file).  Used
    // to charge a volume against the madvise WILLNEED budget.
    size_t willneed_size_full() const;

    // Hint that the posting body will be touched at random (page cache is
    // preserved but readahead is suppressed).  The dictionary head is
    // re-pinned WILLNEED + HUGEPAGE.  Used by the fallback batch path
    // when a single volume's posting body is larger than the WILLNEED
    // budget allows.
    void   apply_madvise_posting_random();

    // Dictionary WILLNEED + HUGEPAGE, posting body MADV_RANDOM (suppress
    // default readahead).  A subsequent apply_madvise_posting_ranges()
    // actively pre-faults the requested ranges regardless of the hint.
    void   apply_madvise_dict_only();

    // Byte size of the dictionary region (header + Elias-Fano blob).
    // The first dict_size() bytes of the mapping are the dictionary;
    // everything beyond is the posting body.
    size_t dict_size() const;

    // Apply MADV_WILLNEED to a set of posting-body ranges. Each pair is
    // (offset, length) in posting-file coordinates — the same coordinate
    // space returned by posting_list_range. Used by the orchestrator's
    // range-WILLNEED path so the kernel only prefetches the posting
    // lists this batch will actually probe.
    void apply_madvise_posting_ranges(
        const std::vector<std::pair<uint64_t, uint64_t>>& ranges);

    // Get posting list byte offset for a k-mer
    uint64_t posting_list_offset(uint32_t kmer) const {
        return dict_.access(kmer);
    }

    // Byte length of posting list for a k-mer.  The Elias-Fano dictionary
    // resolves both endpoints in a single pair-fetch.
    uint64_t posting_list_byte_length(uint32_t kmer) const {
        uint64_t s, e;
        dict_.access_pair(kmer, s, e);
        return e - s;
    }

    // Hot-path helper: fetch (offset, byte_length) for a k-mer in one
    // EF access_pair, halving dispatcher round-trips on the search path.
    void posting_list_range(uint32_t kmer,
                            uint64_t& offset,
                            uint64_t& byte_length) const {
        uint64_t s, e;
        dict_.access_pair(kmer, s, e);
        offset = s;
        byte_length = e - s;
    }

    // Bulk count all postings: for each k-mer, returns the number of
    // distinct sequences containing it (intra-sequence duplicates are
    // removed at build time).  Used by ikafssninfo and index_filter for
    // build-time exclusion via .khx; not used on the search hot path.
    std::vector<uint32_t> bulk_count_postings() const;

    // Optional per-volume decode cache shared read-only by every Stage 1
    // probe on this volume.  The orchestrator publishes the cache on the
    // mutable reader before entering the parallel Stage 1 region (publish-
    // once, full barrier => happens-before) and clears it afterwards, so
    // the hot path reads decode_cache() without a data race.  When set,
    // stage1_filter looks each k-mer up here before falling back to an
    // on-the-fly posting-list decode.
    void set_decode_cache(const DecodedKmerCache* c) { decode_cache_ = c; }
    const DecodedKmerCache* decode_cache() const { return decode_cache_; }

private:
    MmapFile mmap_;
    const KixHeader* header_ = nullptr;
    ef::EFDictionary dict_;
    const uint8_t* posting_file_ = nullptr;
    size_t posting_file_size_ = 0;
    uint32_t table_size_ = 0;
    const DecodedKmerCache* decode_cache_ = nullptr;
};

} // namespace ikafssn
