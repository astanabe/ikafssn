#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Writes a .kix file (format version 3):
// 1. Header
// 2. offsets[table_size + 1]  (sentinel at end = total posting data bytes)
// 3. Delta-compressed ID postings
class KixWriter {
public:
    KixWriter(int k, uint8_t kmer_type);

    // Set metadata (optional, call before write)
    void set_volume_info(uint16_t volume_index, uint16_t total_volumes);
    void set_db(const std::string& name);
    void set_num_sequences(uint32_t n);
    void set_flags(uint32_t flags);

    // Add a posting list for a k-mer. postings must be sorted by seq_id.
    // Caller must call this for k-mers in ascending order (0, 1, 2, ..., 4^k-1).
    // Empty posting lists should be added with count=0 / empty vector.
    void add_posting_list(uint32_t kmer_value, const std::vector<uint32_t>& seq_ids);

    // Finalize and write to file. Returns true on success.
    bool write(const std::string& path);

private:
    int k_;
    uint8_t kmer_type_;
    uint32_t num_sequences_ = 0;
    uint16_t volume_index_ = 0;
    uint16_t total_volumes_ = 1;
    uint32_t flags_ = 0;
    std::string db_;

    uint32_t table_size_;

    // Accumulated data: offsets has table_size_ + 1 entries
    std::vector<uint64_t> offsets_;
    std::vector<uint8_t> posting_data_; // all delta-compressed postings concatenated
    uint64_t total_postings_ = 0;
};

} // namespace ikafssn
