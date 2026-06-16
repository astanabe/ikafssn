#pragma once

#include "core/spaced_seed.hpp"
#include "util/cli_parser.hpp"

#include <cstdint>
#include <string>

namespace ikafssn {

// Parse "-max_degen_expand" with range check [0, 256].
// `default_val` is used when the flag is absent. Returns false on out-of-range.
bool parse_max_degen_expand(const CliParser& cli, int default_val,
                            uint16_t& out, std::string& err);

// Validate "-stage3_score_matrix" value is one of degmatch/dnafull/nuc44.
// Sets `out` to the raw string. When the flag is absent, leaves `out` untouched
// and returns true (caller decides default behaviour).
bool parse_score_matrix(const CliParser& cli, std::string& out, std::string& err);

// Map a validated score-matrix name to its protocol code:
// 1=degmatch, 2=dnafull, 3=nuc44, 0=empty/unknown.
uint8_t score_matrix_code(const std::string& name);

// Parse "-t" template length, accepting only {0, 13, 15, 16, 18, 21}.
bool parse_spaced_seed_t(const CliParser& cli, uint8_t& out, std::string& err);

// Parse "-template_type" string ("coding"/"optimal"/"both").
// When `allow_contiguous` is false (typical for search-side), kContiguous is rejected.
// `default_value` is returned when the flag is absent.
bool parse_template_type_cli(const CliParser& cli, TemplateType default_value,
                             bool allow_contiguous,
                             TemplateType& out, std::string& err);

// Reject the combination "-t 0 -template_type ...": the contiguous index has
// no template type, so pairing it with an explicit -template_type is an error
// in every tool.  Returns true (and leaves `err` untouched) otherwise.
bool validate_t_template_type_combo(const CliParser& cli, std::string& err);

// Cross-check primer-mode options shared by ikafssnsearch and ikafssnclient:
// rejects combinations with -stage1_min_score / -stage2_min_score, and
// requires -insert_length when -primer is given.
bool validate_primer_mode_options(const CliParser& cli, std::string& err);

// Validate the HTTP authentication options shared by ikafssnclient and
// ikafssninfo: -user / -http_user / -http_password / -netrc_file.
//   * When the transport is not HTTP (socket / TCP), any of them is an error.
//   * On HTTP, the three methods (-user, -http_user[+-http_password],
//     -netrc_file) are mutually exclusive, and -http_password requires
//     -http_user.
// Returns true (leaving `err` untouched) when the options are consistent.
bool validate_http_auth_options(const CliParser& cli, bool is_http_transport,
                                std::string& err);

} // namespace ikafssn
