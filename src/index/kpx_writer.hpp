#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Writes a .kpx file (v7 layout).
// Position postings are stored as absolute values; the codec partitions
// each per-kmer posting on (seq_id) so that high-multiplicity clusters
// get their own partition group while the rest go into a self-describing
// short bucket (see src/index/pfd_codec.hpp).
class KpxWriter {
public:
    // freq_threshold_part: per-(kmer, seq_id) occurrence count strictly
    // greater than this becomes a partition group.  Default 8; max 255
    // (the short bucket stores per-seq_id occurrence counts as u8).
    explicit KpxWriter(int k, uint32_t freq_threshold_part = 8);

    struct PostingEntry {
        uint32_t seq_id;
        uint32_t pos;
    };

    // Add position posting list for a k-mer.
    // entries must be sorted by (seq_id, pos): all positions for one seq_id
    // are contiguous, ordered by absolute position.  Called for k-mers in
    // ascending order.
    void add_posting_list(uint32_t kmer_value, const std::vector<PostingEntry>& entries);

    // Write the .kpx file. Returns true on success.
    bool write(const std::string& path) const;

    uint64_t total_position_count() const { return total_position_count_; }

private:
    int k_;
    uint32_t table_size_;
    uint32_t freq_threshold_part_;
    std::vector<uint64_t> pos_offsets_;
    std::vector<uint8_t> posting_data_;
    uint64_t total_position_count_ = 0;
};

} // namespace ikafssn
