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

    // Total residue count across every sequence in this volume.
    uint64_t total_length() const;

    // Memory-map access pattern hint passed through to the kernel
    // (CSeqDB::SetMMapStrategy under the hood).  NOT scoped to this reader's
    // volume: CSeqDBAtlas is a process-wide singleton and applies the hint
    // to every file it has mapped, so one call covers all open volumes.
    //   kNormal     — OS default; readahead enabled.
    //   kRandom     — disable readahead (sparse OID lookups).
    //   kSequential — aggressive readahead (whole-volume scans).
    //   kWillNeed   — pre-fault pages into the page cache up front.
    //   kDontNeed   — release pages from the page cache.
    //   kNoReuse    — hint that data will be touched only once.
    enum class MMapStrategy {
        kNormal,
        kRandom,
        kSequential,
        kWillNeed,
        kDontNeed,
        kNoReuse,
    };
    void set_mmap_strategy(MMapStrategy s) const;

    // Sequence length in bases for given OID.
    uint32_t seq_length(uint32_t oid) const;

    // Get subsequence [start, end] (0-based inclusive) as an IUPAC string,
    // decoding only the requested range from ncbi2na.
    // end is clamped to seq_length-1 if out of range.
    // Returns empty string if start > end (after clamping) or on error.
    std::string get_subsequence(uint32_t oid, uint32_t start, uint32_t end) const;

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

    // Collect every accession registered for the given OID into `out`,
    // joined by '\x01' (BLAST's multi-defline separator).  The joined form
    // is preserved in `.ksx` and emitted as-is in the `sseqid` output
    // field; consumers split on '\x01' (see io/accession_utils.hpp's
    // `split_accessions`).
    // Returns false when the lookup failed, distinguishing an unreadable
    // volume from an unlabelled OID (true with an empty `out`).
    bool get_accession(uint32_t oid, std::string& out) const;

    // Reverse of get_accession(): collect the OIDs of this volume carrying
    // `accession`.  Needs the accession index `makeblastdb -parse_seqids`
    // builds (v5: `.nos` beside the LMDB; v4: the `.nsi` / `.nsd` string
    // ISAM).  Returns false when the lookup failed; an absent accession
    // yields true with an empty `oids`.
    bool accession_to_oids(const std::string& accession,
                           std::vector<uint32_t>& oids) const;

    // Get DB title.
    std::string get_title() const;

    // Static: find all volume paths for a DB prefix.
    static std::vector<std::string> find_volume_paths(const std::string& db);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Result of extract_context_subseq().
struct ContextSubseq {
    std::string seq;        // decoded bases (empty if the OID is empty,
                            // the range is invalid, or decoding failed)
    uint32_t ext_start = 0; // 0-based parent-relative start actually used
    uint32_t ext_end = 0;   // 0-based parent-relative end actually used (exclusive)
    uint32_t seq_len = 0;   // parent OID length in bases
};

// Extract the match region [sstart, send) (0-based parent-relative half-open,
// the convention every hit coordinate uses internally) expanded by `ctx` bases
// on each side, clamped to the parent OID, and decode it via partial ncbi2na
// unpack.  Shared by Stage 3 alignment and ikafssnretrieve.
ContextSubseq extract_context_subseq(const BlastDbReader& reader, uint32_t oid,
                                     uint32_t sstart, uint32_t send, uint32_t ctx);

} // namespace ikafssn
