#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include "io/mmap_file.hpp"
#include "index/ef_codec.hpp"
#include "index/kpx_format.hpp"

namespace ikafssn {

class KpxReader {
public:
    bool open(const std::string& path);
    void close();
    bool is_open() const { return mmap_.is_open(); }

    const KpxHeader& header() const { return *header_; }
    int k() const { return header_->k; }
    uint8_t t() const { return header_->t; }
    uint8_t template_type() const { return header_->template_type; }
    uint64_t total_position_count() const { return header_->total_position_count; }
    uint32_t table_size() const { return table_size_; }

    // Raw pointer to the start of the position posting file
    const uint8_t* posting_file() const { return posting_file_; }
    size_t posting_file_size() const { return posting_file_size_; }

    // pos_offset(kmer) resolves through the EF dictionary.  .kpx has
    // no sentinel — callers needing the byte length of the last k-mer's
    // posting list use posting_file_size() as a loose upper bound (the
    // .kpx posting list bodies are self-delimiting).
    uint64_t pos_offset(uint32_t kmer) const {
        return pos_dict_.access(kmer);
    }

    // Hot-path helper: fan out (lo, hi) for a k-mer's .kpx posting list
    // slice in one EF access_pair.  hi is set to posting_file_size()
    // when kmer is the last entry.
    void pos_offset_range(uint32_t kmer,
                          uint64_t& lo, uint64_t& hi) const {
        if (kmer + 1u < table_size_) {
            pos_dict_.access_pair(kmer, lo, hi);
        } else {
            lo = pos_dict_.access(kmer);
            hi = posting_file_size_;
        }
    }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

    // Byte size of the whole mapping (pos_offsets dictionary head + .kpx
    // posting file).  Used to charge a volume against the madvise budget in
    // Stage 2 / Stage 3 (mode 1 never reads .kpx).
    size_t willneed_size_full() const;

    // Dictionary WILLNEED + HUGEPAGE, posting body MADV_RANDOM (suppress
    // default readahead).  Stage 2A reads the posting body at random, so
    // demand faults track real need without readahead over-fetch.
    void   apply_madvise_dict_only();

    // Byte size of the dictionary region (header + Elias-Fano blob).
    size_t dict_size() const;

private:
    MmapFile mmap_;
    const KpxHeader* header_ = nullptr;
    ef::EFDictionary pos_dict_;
    const uint8_t* posting_file_ = nullptr;
    size_t posting_file_size_ = 0;
    uint32_t table_size_ = 0;
};

} // namespace ikafssn
