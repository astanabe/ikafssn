#pragma once

#include <string>

#include "protocol/messages.hpp"

namespace ikafssn {

// Format database capability listing for error messages.
// "  db    k=9 (mode 1-2), k=11 (mode 1-3)\n" per database.
std::string format_all_databases(const InfoResponse& info);

// Full server info display for ikafssninfo remote mode.
std::string format_server_info(const InfoResponse& info, bool verbose);

// Validate db/k/mode against server capabilities.
// check_slots=true adds soft slot check (active >= max -> error).
// Returns "" on success, or error message string on failure.
// Error messages include full server capability listing.
std::string validate_info(const InfoResponse& info,
                          const std::string& db,
                          uint8_t k, uint8_t mode,
                          bool check_slots);

} // namespace ikafssn
