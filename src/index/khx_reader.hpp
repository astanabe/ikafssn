#pragma once

#include <cstdint>
#include <string>
#include "io/mmap_file.hpp"

namespace ikafssn {

class KhxReader {
public:
    bool open(const std::string& path);
    void close();

    bool is_open() const { return mmap_.is_open(); }

    int k() const { return k_; }
    uint8_t t() const { return t_; }
    uint8_t template_type() const { return template_type_; }

    // Check if a k-mer was excluded during index build.
    bool is_excluded(uint64_t kmer_idx) const {
        return (bitset_[kmer_idx / 8] >> (kmer_idx % 8)) & 1;
    }

    // Count total number of excluded k-mers.
    uint64_t count_excluded() const;

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

private:
    MmapFile mmap_;
    int k_ = 0;
    uint8_t t_ = 0;
    uint8_t template_type_ = 0;
    const uint8_t* bitset_ = nullptr;
    uint64_t tbl_size_ = 0;
};

} // namespace ikafssn
