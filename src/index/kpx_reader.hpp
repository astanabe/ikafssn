#pragma once

#include <cstdint>
#include <string>
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
    // Phase 7e: pos_offsets dictionary is Elias-Fano; the legacy
    // offset_type byte is preserved on the header as reserved/sentinel
    // (set to 0xFF on EF writes by the writer).
    bool is_offset32() const { return false; }

    // Raw pointer to the start of the position posting file
    const uint8_t* posting_file() const { return posting_file_; }
    size_t posting_file_size() const { return posting_file_size_; }

    // Phase 7e: pos_offset(kmer) resolves through the EF dictionary.
    // .kpx has no sentinel — callers needing the byte length of the
    // last k-mer's posting list use posting_file_size() as a loose
    // upper bound (the .kpx posting list bodies are self-delimiting).
    uint64_t pos_offset(uint32_t kmer) const {
        return pos_dict_.access(kmer);
    }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

private:
    MmapFile mmap_;
    const KpxHeader* header_ = nullptr;
    ef::EFDictionary pos_dict_;
    const uint8_t* posting_file_ = nullptr;
    size_t posting_file_size_ = 0;
    uint32_t table_size_ = 0;
};

} // namespace ikafssn
