#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include "io/mmap_file.hpp"
#include "index/ef_codec.hpp"
#include "index/kix_format.hpp"

namespace ikafssn {

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

    // Byte size of the whole mapping (dictionary head + the entire posting
    // file).  Used by ikafssnsearch to charge a per-volume batch against the
    // madvise WILLNEED budget — including the posting list body — so the
    // batched search loop pre-faults a whole volume before issuing
    // query × volume jobs, then releases it for the next batch.
    size_t willneed_size_full() const;

    // Dictionary WILLNEED + HUGEPAGE, posting body MADV_RANDOM (suppress
    // default readahead).  Stage 1 / Stage 2A access the posting body at
    // random, so demand faults track real need without readahead over-fetch.
    void   apply_madvise_dict_only();

    // Byte size of the dictionary region (header + Elias-Fano blob).
    // The first dict_size() bytes of the mapping are the dictionary;
    // everything beyond is the posting body.
    size_t dict_size() const;

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

private:
    MmapFile mmap_;
    const KixHeader* header_ = nullptr;
    ef::EFDictionary dict_;
    const uint8_t* posting_file_ = nullptr;
    size_t posting_file_size_ = 0;
    uint32_t table_size_ = 0;
};

} // namespace ikafssn
