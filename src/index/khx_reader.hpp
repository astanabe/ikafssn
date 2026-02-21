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

    // Check if a k-mer was excluded during index build.
    bool is_excluded(uint64_t kmer_idx) const {
        return (bitset_[kmer_idx / 8] >> (kmer_idx % 8)) & 1;
    }

    // Count total number of excluded k-mers.
    uint64_t count_excluded() const;

private:
    MmapFile mmap_;
    int k_ = 0;
    const uint8_t* bitset_ = nullptr;
    uint64_t tbl_size_ = 0;
};

} // namespace ikafssn
