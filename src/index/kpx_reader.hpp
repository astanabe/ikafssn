#pragma once

#include <cstdint>
#include <string>
#include "io/mmap_file.hpp"
#include "index/kpx_format.hpp"

namespace ikafssn {

class KpxReader {
public:
    bool open(const std::string& path);
    void close();

    const KpxHeader& header() const { return *header_; }
    int k() const { return header_->k; }
    uint64_t total_postings() const { return header_->total_postings; }
    uint64_t table_size() const { return table_size_; }

    const uint64_t* pos_offsets() const { return pos_offsets_; }

    // Raw pointer to position posting data
    const uint8_t* posting_data() const { return posting_data_; }
    size_t posting_data_size() const { return posting_data_size_; }

    uint64_t pos_offset(uint64_t kmer) const { return pos_offsets_[kmer]; }

private:
    MmapFile mmap_;
    const KpxHeader* header_ = nullptr;
    const uint64_t* pos_offsets_ = nullptr;
    const uint8_t* posting_data_ = nullptr;
    size_t posting_data_size_ = 0;
    uint64_t table_size_ = 0;
};

} // namespace ikafssn
