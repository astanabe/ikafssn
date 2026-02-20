#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include "io/mmap_file.hpp"

namespace ikafssn {

class KsxReader {
public:
    bool open(const std::string& path);
    void close();

    uint32_t num_sequences() const { return num_sequences_; }
    uint32_t seq_length(uint32_t oid) const;
    std::string_view accession(uint32_t oid) const;

private:
    MmapFile mmap_;
    uint32_t num_sequences_ = 0;
    const uint32_t* seq_lengths_ = nullptr;
    const uint32_t* acc_offsets_ = nullptr;
    const char* acc_strings_ = nullptr;
};

} // namespace ikafssn
