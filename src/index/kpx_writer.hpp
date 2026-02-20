#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Writes a .kpx file.
// Position postings are delta-compressed with reset at sequence boundaries.
// The caller must provide both seq_ids (to detect boundaries) and positions.
class KpxWriter {
public:
    explicit KpxWriter(int k);

    struct PostingEntry {
        uint32_t seq_id;
        uint32_t pos;
    };

    // Add position posting list for a k-mer.
    // entries must be in the same order as the corresponding ID postings in .kix.
    // Called for k-mers in ascending order.
    void add_posting_list(uint64_t kmer_value, const std::vector<PostingEntry>& entries);

    // Write the .kpx file. Returns true on success.
    bool write(const std::string& path) const;

    uint64_t total_postings() const { return total_postings_; }

private:
    int k_;
    uint64_t table_size_;
    std::vector<uint64_t> pos_offsets_;
    std::vector<uint8_t> posting_data_;
    uint64_t total_postings_ = 0;
};

} // namespace ikafssn
