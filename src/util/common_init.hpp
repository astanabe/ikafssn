#pragma once

#include "util/cli_parser.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <thread>

// IKAFSSN_VERSION is defined in the generated core/version.hpp.
// Callers must include core/version.hpp before using check_version().

namespace ikafssn {

// Print "<cmd_name> <version>" to stderr if --version is present.
// Returns true if --version was handled (caller should return 0).
inline bool check_version(const CliParser& cli, const char* cmd_name) {
    if (cli.has("--version")) {
        std::fprintf(stderr, "%s %s\n", cmd_name, IKAFSSN_VERSION);
        return true;
    }
    return false;
}

// Create a Logger from -v / --verbose flags.
inline Logger make_logger(const CliParser& cli) {
    bool verbose = cli.has("-v") || cli.has("--verbose");
    return Logger(verbose ? Logger::kDebug : Logger::kInfo);
}

// Resolve thread count from CLI (0 or negative → hardware_concurrency).
inline int resolve_threads(const CliParser& cli,
                           const std::string& key = "-threads") {
    int n = cli.get_int(key, 0);
    if (n <= 0) {
        n = static_cast<int>(std::thread::hardware_concurrency());
        if (n <= 0) n = 1;
    }
    return n;
}

} // namespace ikafssn
