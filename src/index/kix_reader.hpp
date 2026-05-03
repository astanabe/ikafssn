#pragma once

#include <cstdint>
#include <string>
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
    // Phase 7a: dictionary is Elias-Fano; the legacy flag is no longer
    // consulted at read time (kept on the header for byte-stability).
    bool is_offset32() const { return false; }

    // Raw pointer to the start of the ID posting file
    const uint8_t* posting_file() const { return posting_file_; }
    size_t posting_file_size() const { return posting_file_size_; }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

    // Get posting list byte offset for a k-mer
    uint64_t posting_list_offset(uint32_t kmer) const {
        return dict_.access(kmer);
    }

    // Byte length of posting list for a k-mer.  The Elias-Fano dictionary
    // resolves both endpoints in a single pair-fetch (a single select1 +
    // adjacent low-bits read in 7b SIMD; two sequential accesses in the
    // 7a scalar PoC).
    uint64_t posting_list_byte_length(uint32_t kmer) const {
        uint64_t s, e;
        dict_.access_pair(kmer, s, e);
        return e - s;
    }

    // Phase 7d hot-path helper: fetch (offset, byte_length) for a k-mer
    // in one EF access_pair, halving the dispatcher round-trips that
    // Stage 1 / Stage 2 paid by calling posting_list_offset(kmer) +
    // posting_list_offset(kmer+1) (or _byte_length + _offset).
    void posting_list_range(uint32_t kmer,
                            uint64_t& offset,
                            uint64_t& byte_length) const {
        uint64_t s, e;
        dict_.access_pair(kmer, s, e);
        offset = s;
        byte_length = e - s;
    }

    // Bulk count all postings: for each k-mer, returns the number of
    // distinct sequences containing it (v7: intra-sequence duplicates are
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
