#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class BlastDbReader;
class Logger;

// Configuration for index building.
struct IndexBuilderConfig {
    int k = 11;                         // k-mer length
    uint64_t memory_limit = uint64_t(8) << 30; // per-volume memory budget (default: 8 GB)
    bool keep_tmp = false;              // true: keep .tmp files (skip rename to final)
    int threads = 1;                    // threads (counting + partition scan + sort)
    bool verbose = false;
    bool skip_kpx = false;              // true: skip .kpx generation (mode 1 index)
    int max_degen_expand = 4;           // max degenerate expansion per k-mer (0/1: disable)
    uint8_t t = 0;                      // template length (0=contiguous, 16/18/21=spaced)
    uint8_t template_type = 0;          // TemplateType enum value
    // Per-(kmer, seq_id) partition threshold for .kpx. A (k-mer, seq_id)
    // cluster with occurrence count > threshold is split into its own
    // partition group; smaller clusters merge into a shared short bucket.
    // Default 8 mirrors -max_freq_build's `>` comparison semantics.
    uint32_t freq_threshold_part = 8;

    // Fragment-indexing parameters.
    uint32_t min_seq_length    = 64;  // shorter sequences are skipped at index time
    uint32_t min_length_split  = 0;   // 0 = splitting disabled
    uint32_t overlap_length    = 0;   // 0 = splitting disabled
};

// Per-fragment metadata captured during the metadata pass.  Consumed
// by build_postings() to drive the per-fragment k-mer scanner.
struct VolumeMetadata {
    uint32_t num_sequences = 0;             // total fragment count for this volume
    uint32_t num_parents = 0;
    std::vector<uint32_t> seq_id_to_blast_oid;
    std::vector<uint32_t> seq_id_to_frag_start;
    std::vector<uint32_t> seq_id_to_frag_end;
};

// Collect parent / fragment metadata, write .ksx.tmp, and populate
// `out` for the postings pass.  Returns true on success.
bool build_metadata(BlastDbReader& db,
                    const IndexBuilderConfig& config,
                    const std::string& output_prefix,
                    VolumeMetadata& out,
                    const Logger& logger);

// Scan k-mers per fragment, emit per-partition posting lists, and
// write .kix.tmp / .kpx.tmp.  `volume_index` / `total_volumes` /
// `db_name` are recorded in the .kix header for cross-volume
// bookkeeping.
template <typename KmerInt>
bool build_postings(BlastDbReader& db,
                    const IndexBuilderConfig& config,
                    const std::string& output_prefix,
                    const VolumeMetadata& meta,
                    uint16_t volume_index,
                    uint16_t total_volumes,
                    const std::string& db_name,
                    const Logger& logger);

// Convenience wrapper: run build_metadata + build_postings for one
// volume and rename the resulting .ksx.tmp / .kix.tmp / .kpx.tmp to
// their final names.  Used by tests to stand up a single-volume index.
template <typename KmerInt>
bool build_index(BlastDbReader& db,
                 const IndexBuilderConfig& config,
                 const std::string& output_prefix,
                 uint16_t volume_index,
                 uint16_t total_volumes,
                 const std::string& db_name,
                 const Logger& logger);

} // namespace ikafssn
