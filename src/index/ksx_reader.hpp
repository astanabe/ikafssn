#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include "io/mmap_file.hpp"

namespace ikafssn {

// Two-stage .ksx reader (v10 layout — see src/index/ksx_format.hpp).
//
// External callers continue to use the simple SeqId-based API:
//   - seq_length(seq_id):   fragment length (= frag_end - frag_start + 1)
//   - accession(seq_id):    accession string of the fragment's parent OID
// On top of that the reader exposes the new fragment-aware accessors
// (parent_index / fragment_start / fragment_end / blast_oid / parent_length /
//  num_parents / num_sequences / min_seq_length / min_length_split /
//  overlap_length).
class KsxReader {
public:
    bool open(const std::string& path);
    void close();

    bool is_open() const { return mmap_.is_open(); }
    uint32_t num_sequences() const { return num_sequences_; }
    uint32_t num_parents()   const { return num_parents_; }

    // Convenience accessors keyed by internal SeqId (= fragment number).
    uint32_t seq_length(uint32_t seq_id) const;
    std::string_view accession(uint32_t seq_id) const;

    // Fragment-aware accessors.
    uint32_t parent_index(uint32_t seq_id)   const { return fragment_parent_idx_[seq_id]; }
    uint32_t fragment_start(uint32_t seq_id) const { return fragment_start_[seq_id]; }
    uint32_t fragment_end(uint32_t seq_id)   const { return fragment_end_[seq_id]; }

    // Parent-keyed accessors.
    uint32_t parent_length(uint32_t parent_idx) const { return parent_lengths_[parent_idx]; }
    uint32_t blast_oid(uint32_t parent_idx)     const { return parent_blast_oids_[parent_idx]; }
    std::string_view parent_accession(uint32_t parent_idx) const;

    // Fragment-indexing header values.
    uint32_t min_seq_length()   const { return min_seq_length_; }
    uint32_t min_length_split() const { return min_length_split_; }
    uint32_t overlap_length()   const { return overlap_length_; }

    // madvise budget API
    size_t willneed_size() const;
    void apply_madvise(bool willneed);

private:
    MmapFile mmap_;
    uint32_t num_sequences_   = 0;
    uint32_t num_parents_     = 0;
    uint32_t min_seq_length_  = 0;
    uint32_t min_length_split_= 0;
    uint32_t overlap_length_  = 0;

    const uint32_t* parent_lengths_       = nullptr;
    const uint32_t* parent_blast_oids_    = nullptr;
    const uint32_t* parent_acc_offsets_   = nullptr;
    const char*     parent_acc_strings_   = nullptr;

    const uint32_t* fragment_parent_idx_  = nullptr;
    const uint32_t* fragment_start_       = nullptr;
    const uint32_t* fragment_end_         = nullptr;
};

} // namespace ikafssn
