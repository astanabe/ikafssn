#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KsxReader;
struct FastaRecord;
struct OrchestratorHit;

struct OutputHit {
    std::string qseqid;
    // Parent OID accession.  When the index uses fragment splitting
    // (overlap_length > 0), this is the parent accession common to all
    // fragments of a single OID, NOT a per-fragment accession.
    std::string sseqid;
    char sstrand;         // '+' or '-'
    uint32_t qstart;
    uint32_t qend;
    // Parent-OID-relative subject coordinates (1-based).  The
    // orchestrator writes Stage 2 chains in fragment-relative space and
    // the (orchestrator -> output) boundary shifts them by the matching
    // KsxReader::fragment_start() so downstream tools (Stage 3, output
    // writers, ikafssnretrieve) all see one canonical coordinate system
    // per parent OID.
    uint32_t sstart;
    uint32_t send;
    uint32_t coverscore = 0;
    uint32_t matchscore = 0;
    uint32_t chainscore = 0;
    uint16_t volume;
    uint32_t oid = 0;        // internal: BLAST DB OID (not written to output)

    // Stage 3 fields (populated only when mode == 3)
    int32_t alnscore = 0;
    std::string cigar;
    uint32_t npositive = 0;
    uint32_t nnegative = 0;
    double ppositive = 0.0;
    std::string qseq;        // aligned query (with gaps, traceback only)
    std::string sseq;        // aligned subject (with gaps, traceback only)
    uint32_t qlen = 0;       // query full sequence length
    // Subject full sequence length = parent OID length.  Always derived
    // from KsxReader::parent_length(parent_idx); never the (smaller)
    // fragment length.
    uint32_t slen = 0;

    // Skip marker: when skip_reason != 0 this OutputHit represents a query
    // that was not searched (no Stage 1/2/3 ran). It carries qseqid + qlen
    // and is emitted as a sentinel row in TSV / JSON / SAM so consumers can
    // see why the query produced no hits.
    uint8_t     skip_reason = 0;   // ikafssn::SkipReason (0 = normal hit)
    std::string skip_detail;       // human-readable detail
};

enum class OutputFormat { kTsv, kJson, kSam, kBam };

// Parse an output format string ("tsv", "json", "sam", "bam").
// Returns true on success. On failure, out is unchanged and error_msg is set.
bool parse_output_format(const std::string& str, OutputFormat& out,
                         std::string& error_msg);

// Validate that the output format is compatible with mode/traceback settings.
// Returns true if valid. On failure, sets error_msg.
bool validate_output_format(OutputFormat fmt, uint8_t mode, bool traceback,
                            const std::string& output_path,
                            std::string& error_msg);

// Write results in TSV (tab-delimited) format.
void write_results_tsv(std::ostream& out,
                       const std::vector<OutputHit>& hits,
                       uint8_t mode = 2,
                       uint8_t stage1_score_type = 1,
                       bool stage3_traceback = false);

// Write results in JSON format.
void write_results_json(std::ostream& out,
                        const std::vector<OutputHit>& hits,
                        uint8_t mode = 2,
                        uint8_t stage1_score_type = 1,
                        bool stage3_traceback = false);

// Write per-query JSON objects without the outer {"results": [...]} wrapper.
// Used for checkpoint batch intermediate files.
void write_results_json_fragment(std::ostream& out,
                                  const std::vector<OutputHit>& hits,
                                  uint8_t mode = 2,
                                  uint8_t stage1_score_type = 1,
                                  bool stage3_traceback = false);

// Write results in the specified format (tsv or json).
void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt,
                   uint8_t mode = 2,
                   uint8_t stage1_score_type = 1,
                   bool stage3_traceback = false);

// -------------------------------------------------------------------------
// Mode 1 parallel TSV / JSON writers.
//
// Mode 1's output (k-mer-only "did this query hit this subject" rows) is
// the only path where the OrchestratorHit -> OutputHit boundary is a
// dominant single-thread cost: nt_v4-scale runs spend a significant
// fraction of wall on the dump.  These writers consume OrchestratorHit
// directly, format per-thread chunk strings via std::to_chars, and let
// the caller drive the (serial) compressed output sink.
//
// Mode 2 / Mode 3 go through the OutputHit-based writers above, which
// carry chain / Stage 3 / SAM fields.
// -------------------------------------------------------------------------
struct Mode1ParallelInputs {
    const std::vector<OrchestratorHit>* hits = nullptr;
    const std::vector<FastaRecord>*     queries = nullptr;
    // Indexed by OrchestratorHit::volume_idx (the bundle index, not
    // volume_index).  Each entry resolves accession / seq_length for the
    // sids produced by that volume.  Caller owns the readers.
    std::vector<const KsxReader*>       ksx_per_volume;
    const std::vector<uint8_t>*         skip_reason = nullptr;
    const std::vector<std::string>*     skip_detail = nullptr;
    uint8_t  stage1_score_type = 1;
    int      nthread = 1;
};

// TSV mode 1 parallel writer: header is emitted serially, hit rows are
// formatted into per-thread chunk buffers (parallel_for over the hits
// vector range) and written to `out` in chunk order.  Skip rows are
// appended after hits.  Caller is responsible for sorting `hits` ahead
// of time when -nresult > 0 (the mode-1 sort is by query_idx + score).
void write_results_tsv_mode1_parallel(std::ostream& out,
                                       const Mode1ParallelInputs& in);

// JSON mode 1 parallel writer: chunk boundaries align to query_idx
// boundaries (each chunk is a contiguous range of queries).  Caller must
// have sorted `hits` by query_idx so that all rows for one query land in
// one chunk.  Skip queries are emitted as objects with "status":
// "skipped" mirroring the existing OutputHit writer.
void write_results_json_mode1_parallel(std::ostream& out,
                                        const Mode1ParallelInputs& in);

} // namespace ikafssn
