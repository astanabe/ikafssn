#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KsxReader;

struct OutputHit {
    std::string query_id;
    std::string accession;
    char strand;          // '+' or '-'
    uint32_t q_start;
    uint32_t q_end;
    uint32_t s_start;
    uint32_t s_end;
    uint32_t coverscore = 0;
    uint32_t matchscore = 0;
    uint32_t chainscore = 0;
    uint16_t volume;

    // Stage 3 fields (populated only when mode == 3)
    int32_t alnscore = 0;
    std::string cigar;
    uint32_t nident = 0;
    uint32_t nmismatch = 0;
    double pident = 0.0;
    std::string q_seq;        // aligned query (with gaps, traceback only)
    std::string s_seq;        // aligned subject (with gaps, traceback only)
    uint32_t q_length = 0;    // query full sequence length
    uint32_t s_length = 0;    // subject full sequence length (for SAM @SQ)
};

enum class OutputFormat { kTab, kJson, kSam, kBam };

// Parse an output format string ("tab", "json", "sam", "bam").
// Returns true on success. On failure, out is unchanged and error_msg is set.
bool parse_output_format(const std::string& str, OutputFormat& out,
                         std::string& error_msg);

// Validate that the output format is compatible with mode/traceback settings.
// Returns true if valid. On failure, sets error_msg.
bool validate_output_format(OutputFormat fmt, uint8_t mode, bool traceback,
                            const std::string& output_path,
                            std::string& error_msg);

// Write results in tab-delimited format.
void write_results_tab(std::ostream& out,
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

// Write results in the specified format (tab or json).
void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt,
                   uint8_t mode = 2,
                   uint8_t stage1_score_type = 1,
                   bool stage3_traceback = false);

} // namespace ikafssn
