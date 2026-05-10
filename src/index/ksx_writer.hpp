#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// Two-stage .ksx writer (see src/index/ksx_format.hpp).  Parents are
// registered first via add_parent(); each parent then receives one or
// more fragments via add_fragment().  When fragment splitting is
// disabled (`min_length_split == 0`) the index builder emits exactly
// one fragment per parent spanning the whole parent.
class KsxWriter {
public:
    // Register a parent OID.  Returns the parent index (0-based, in
    // registration order) for use with add_fragment().
    uint32_t add_parent(uint32_t blast_oid,
                        uint32_t parent_length,
                        const std::string& accession);

    // Register a fragment derived from `parent_idx`.  start/end are
    // 1-based, parent-relative coordinates (inclusive).  Must satisfy
    // 1 <= start <= end <= parent_length.  Returns the internal SeqId
    // (0-based) assigned to this fragment.
    uint32_t add_fragment(uint32_t parent_idx,
                          uint32_t frag_start,
                          uint32_t frag_end);

    void set_min_seq_length(uint32_t v)   { min_seq_length_ = v; }
    void set_min_length_split(uint32_t v) { min_length_split_ = v; }
    void set_overlap_length(uint32_t v)   { overlap_length_ = v; }

    // Write the .ksx file. Returns true on success.
    bool write(const std::string& path) const;

    uint32_t num_sequences() const {
        return static_cast<uint32_t>(fragment_parent_idx_.size());
    }
    uint32_t num_parents() const {
        return static_cast<uint32_t>(parent_lengths_.size());
    }

private:
    // Parent table
    std::vector<uint32_t>    parent_lengths_;
    std::vector<uint32_t>    parent_blast_oids_;
    std::vector<std::string> parent_accessions_;

    // Fragment table (size = num_sequences)
    std::vector<uint32_t> fragment_parent_idx_;
    std::vector<uint32_t> fragment_start_;
    std::vector<uint32_t> fragment_end_;

    uint32_t min_seq_length_   = 0;
    uint32_t min_length_split_ = 0;
    uint32_t overlap_length_   = 0;
};

} // namespace ikafssn
