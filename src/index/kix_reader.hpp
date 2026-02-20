#pragma once

#include <cstdint>
#include <string>
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
    uint64_t table_size() const { return table_size_; }

    // Direct access to tables
    const uint64_t* offsets() const { return offsets_; }
    const uint32_t* counts() const { return counts_; }

    // Raw pointer to the start of ID posting section
    const uint8_t* posting_data() const { return posting_data_; }
    size_t posting_data_size() const { return posting_data_size_; }

    // Convenience: get posting offset and count for a k-mer
    uint64_t posting_offset(uint64_t kmer) const { return offsets_[kmer]; }
    uint32_t posting_count(uint64_t kmer) const { return counts_[kmer]; }

private:
    MmapFile mmap_;
    const KixHeader* header_ = nullptr;
    const uint64_t* offsets_ = nullptr;
    const uint32_t* counts_ = nullptr;
    const uint8_t* posting_data_ = nullptr;
    size_t posting_data_size_ = 0;
    uint64_t table_size_ = 0;
};

} // namespace ikafssn
