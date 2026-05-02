#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Writes a .kix file (format version 7):
// 1. Header
// 2. dictionary entries `offsets[table_size + 1]`  (sentinel at end = total
//    posting file bytes)
// 3. Per-k-mer posting list:
//      [u32 distinct_count][u32 body_words][u32 body[body_words]]
//    over the distinct seq_id delta stream encoded with FastPFor.
class KixWriter {
public:
    KixWriter(int k, uint8_t kmer_type);

    // Set metadata (optional, call before write)
    void set_volume_info(uint16_t volume_index, uint16_t total_volumes);
    void set_db(const std::string& name);
    void set_num_sequences(uint32_t n);
    void set_flags(uint32_t flags);

    // Add a posting list for a k-mer.  In v7 `seq_ids` must be **sorted
    // and distinct** (the on-disk posting list stores distinct seq_ids only;
    // intra-sequence k-mer duplicates are removed at build time).
    // Caller must call this for k-mers in ascending order.
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
    // All per-k-mer posting lists concatenated, delta-compressed.
    std::vector<uint8_t> posting_file_;
    uint64_t total_distinct_postings_ = 0;
};

} // namespace ikafssn
