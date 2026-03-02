#pragma once

#include "util/cli_parser.hpp"
#include "util/logger.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <thread>
#include <unistd.h>

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

// Detect physical memory and return half of it (minimum 1 GB).
inline uint64_t default_memory_limit() {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page_size > 0) {
        uint64_t half = static_cast<uint64_t>(pages) * static_cast<uint64_t>(page_size) / 2;
        if (half >= (uint64_t(1) << 30)) return half;
    }
    return uint64_t(1) << 30; // fallback: 1 GB
}

// Format a byte size as a human-readable string with K/M/G suffix.
inline std::string format_size(uint64_t bytes) {
    if (bytes >= (uint64_t(1) << 30) && bytes % (uint64_t(1) << 30) == 0)
        return std::to_string(bytes >> 30) + "G";
    if (bytes >= (uint64_t(1) << 20) && bytes % (uint64_t(1) << 20) == 0)
        return std::to_string(bytes >> 20) + "M";
    if (bytes >= (uint64_t(1) << 10) && bytes % (uint64_t(1) << 10) == 0)
        return std::to_string(bytes >> 10) + "K";
    return std::to_string(bytes);
}

} // namespace ikafssn
