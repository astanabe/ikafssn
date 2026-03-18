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
    bool is_open() const { return mmap_.is_open(); }

    const KpxHeader& header() const { return *header_; }
    int k() const { return header_->k; }
    uint8_t t() const { return header_->t; }
    uint8_t template_type() const { return header_->template_type; }
    uint64_t total_postings() const { return header_->total_postings; }
    uint32_t table_size() const { return table_size_; }
    bool is_offset32() const { return offset32_; }

    // Raw pointer to position posting data
    const uint8_t* posting_data() const { return posting_data_; }
    size_t posting_data_size() const { return posting_data_size_; }

    uint64_t pos_offset(uint32_t kmer) const {
        if (offset32_) return pos_offsets32_[kmer];
        return pos_offsets64_[kmer];
    }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

private:
    MmapFile mmap_;
    const KpxHeader* header_ = nullptr;
    const uint64_t* pos_offsets64_ = nullptr;
    const uint32_t* pos_offsets32_ = nullptr;
    bool offset32_ = false;
    const uint8_t* posting_data_ = nullptr;
    size_t posting_data_size_ = 0;
    uint32_t table_size_ = 0;
};

} // namespace ikafssn
