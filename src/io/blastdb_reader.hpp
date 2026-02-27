#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <memory>

namespace ikafssn {

// Wrapper around NCBI CSeqDB for reading BLAST nucleotide databases.
// Provides a simplified interface for k-mer index construction.
class BlastDbReader {
public:
    BlastDbReader();
    ~BlastDbReader();

    // Non-copyable
    BlastDbReader(const BlastDbReader&) = delete;
    BlastDbReader& operator=(const BlastDbReader&) = delete;

    // Move
    BlastDbReader(BlastDbReader&&) noexcept;
    BlastDbReader& operator=(BlastDbReader&&) noexcept;

    // Open a BLAST DB volume by path.
    // Returns true on success, false on error (message to stderr).
    bool open(const std::string& db_path);

    // Close the database.
    void close();

    bool is_open() const;

    // Number of sequences in this volume.
    uint32_t num_sequences() const;

    // Sequence length in bases for given OID.
    uint32_t seq_length(uint32_t oid) const;

    // Get sequence as ACGTN string. N replaces all ambiguous bases.
    // Returns empty string on error.
    std::string get_sequence(uint32_t oid) const;

    // Raw sequence data from BLAST DB (ncbi2na packed + ambiguity data).
    // Pointers are into mmap region; call ret_raw_sequence() when done.
    struct RawSequence {
        const char* ncbi2na_data;  // ncbi2na packed data pointer (mmap)
        int ncbi2na_bytes;         // ncbi2na data byte length
        const char* ambig_data;    // ambiguity data pointer (ncbi2na_data + ncbi2na_bytes)
        int ambig_bytes;           // ambiguity data byte length
        uint32_t seq_length;       // number of bases
    };

    // Get raw ncbi2na packed data + ambiguity data for given OID.
    // Must call ret_raw_sequence() after use.
    RawSequence get_raw_sequence(uint32_t oid) const;

    // Release raw sequence buffer obtained from get_raw_sequence().
    void ret_raw_sequence(const RawSequence& raw) const;

    // Get primary accession for given OID.
    // Returns empty string if not available.
    std::string get_accession(uint32_t oid) const;

    // Get DB title.
    std::string get_title() const;

    // Static: find all volume paths for a DB prefix.
    static std::vector<std::string> find_volume_paths(const std::string& db);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ikafssn
