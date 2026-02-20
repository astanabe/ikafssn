#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Writes a .kix file using Method B (post-hoc offset construction):
// 1. Write header
// 2. Reserve table area (offsets + counts) with zeros
// 3. Write delta-compressed ID postings k-mer by k-mer
// 4. Seek back and overwrite tables with actual values
class KixWriter {
public:
    KixWriter(int k, uint8_t kmer_type);

    // Set metadata (optional, call before write)
    void set_volume_info(uint16_t volume_index, uint16_t total_volumes);
    void set_db_name(const std::string& name);
    void set_num_sequences(uint32_t n);
    void set_flags(uint32_t flags);

    // Add a posting list for a k-mer. postings must be sorted by seq_id.
    // Caller must call this for k-mers in ascending order (0, 1, 2, ..., 4^k-1).
    // Empty posting lists should be added with count=0 / empty vector.
    void add_posting_list(uint64_t kmer_value, const std::vector<uint32_t>& seq_ids);

    // Finalize and write to file. Returns true on success.
    bool write(const std::string& path);

private:
    int k_;
    uint8_t kmer_type_;
    uint32_t num_sequences_ = 0;
    uint16_t volume_index_ = 0;
    uint16_t total_volumes_ = 1;
    uint32_t flags_ = 0;
    std::string db_name_;

    uint64_t table_size_;

    // Accumulated data
    std::vector<uint64_t> offsets_;
    std::vector<uint32_t> counts_;
    std::vector<uint8_t> posting_data_; // all delta-compressed postings concatenated
    uint64_t total_postings_ = 0;
};

} // namespace ikafssn
