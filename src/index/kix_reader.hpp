#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "io/mmap_file.hpp"
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
    uint64_t total_postings() const { return header_->total_postings; }
    uint8_t t() const { return header_->t; }
    uint8_t template_type() const { return header_->template_type; }
    uint32_t table_size() const { return table_size_; }
    bool is_offset32() const { return offset32_; }

    // Raw pointer to the start of ID posting section
    const uint8_t* posting_data() const { return posting_data_; }
    size_t posting_data_size() const { return posting_data_size_; }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

    // Get posting byte offset for a k-mer
    uint64_t posting_offset(uint32_t kmer) const {
        if (offset32_) return offsets32_[kmer];
        return offsets64_[kmer];
    }

    // Byte length of posting data for a k-mer
    uint64_t posting_byte_length(uint32_t kmer) const {
        return posting_offset(kmer + 1) - posting_offset(kmer);
    }

    // Count postings for a k-mer (on-demand varint decode).
    uint32_t count_postings(uint32_t kmer) const;

    // Bulk count all postings. Returns counts[table_size].
    std::vector<uint32_t> bulk_count_postings() const;

private:
    MmapFile mmap_;
    const KixHeader* header_ = nullptr;
    const uint64_t* offsets64_ = nullptr;
    const uint32_t* offsets32_ = nullptr;
    bool offset32_ = false;
    const uint8_t* posting_data_ = nullptr;
    size_t posting_data_size_ = 0;
    uint32_t table_size_ = 0;
};

} // namespace ikafssn
