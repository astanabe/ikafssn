#pragma once

#include "util/cli_parser.hpp"
#include "util/logger.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <unistd.h>

// IKAFSSN_VERSION and IKAFSSN_BUILD_TZ_OFFSET are defined in the generated
// core/version.hpp.  Callers must include core/version.hpp before this header.

namespace ikafssn {

// Format __DATE__ ("Mmm dd yyyy") and __TIME__ ("HH:MM:SS") into ISO 8601
// extended format with the timezone offset captured at CMake configure time.
inline std::string format_build_timestamp(const char* date, const char* time_str) {
    static const char* months[] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };
    int month = 1;
    for (int i = 0; i < 12; i++) {
        if (std::strncmp(date, months[i], 3) == 0) {
            month = i + 1;
            break;
        }
    }
    int day = (date[4] == ' ') ? (date[5] - '0')
                                : ((date[4] - '0') * 10 + (date[5] - '0'));
    int year = std::atoi(date + 7);

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%s%s",
                  year, month, day, time_str, IKAFSSN_BUILD_TZ_OFFSET);
    return std::string(buf);
}

// Print version info if -version or --version is present.
// Uses __DATE__ and __TIME__ default arguments so the build timestamp reflects
// the compile time of the caller's translation unit.
// Returns true if version was handled (caller should return 0).
inline bool check_version(const CliParser& cli, const char* cmd_name,
                           const char* build_date = __DATE__,
                           const char* build_time = __TIME__) {
    if (cli.has("-version") || cli.has("--version")) {
        std::string ts = format_build_timestamp(build_date, build_time);
        std::fprintf(stderr, "%s: %s\n Package: ikafssn %s, build %s\n",
                     cmd_name, IKAFSSN_VERSION, IKAFSSN_VERSION, ts.c_str());
        return true;
    }
    return false;
}

// Print "cmd_name: version\n\n" header for -h / usage display.
inline void print_version_header(const char* cmd_name) {
    std::fprintf(stderr, "%s: %s\n\n", cmd_name, IKAFSSN_VERSION);
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
