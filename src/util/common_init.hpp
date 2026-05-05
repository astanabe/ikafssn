#pragma once

#include "util/cli_parser.hpp"
#include "util/logger.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <unistd.h>

// IKAFSSN_VERSION and IKAFSSN_BUILD_TZ_OFFSET are defined in the generated
// core/version.hpp.  Callers must include core/version.hpp before this header.

namespace ikafssn {

// Format __DATE__ ("Mmm dd yyyy") and __TIME__ ("HH:MM:SS") into ISO 8601
// extended format with the timezone offset captured at CMake configure time.
std::string format_build_timestamp(const char* date, const char* time_str);

// Print version info if -version or --version is present.
// Returns true if version was handled (caller should return 0).
// build_date / build_time default to the caller's __DATE__ / __TIME__ so the
// stamp reflects each translation unit's compile time.
bool check_version(const CliParser& cli, const char* cmd_name,
                   const char* build_date, const char* build_time);

inline bool check_version(const CliParser& cli, const char* cmd_name) {
    return check_version(cli, cmd_name, __DATE__, __TIME__);
}

// Print "cmd_name: version\n" header for -h / usage display.  The caller
// is responsible for any blank line that should follow before the next
// section (e.g. "Usage: ...").
void print_version_header(const char* cmd_name);

// Create a Logger from -v / --verbose flags.
Logger make_logger(const CliParser& cli);

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

// Format a byte size as a compact string with K/M/G suffix.
// Only adds a suffix when the value is an exact multiple of the unit;
// otherwise returns the raw byte count. Suitable for echoing config values.
inline std::string format_size(uint64_t bytes) {
    if (bytes >= (uint64_t(1) << 30) && bytes % (uint64_t(1) << 30) == 0)
        return std::to_string(bytes >> 30) + "G";
    if (bytes >= (uint64_t(1) << 20) && bytes % (uint64_t(1) << 20) == 0)
        return std::to_string(bytes >> 20) + "M";
    if (bytes >= (uint64_t(1) << 10) && bytes % (uint64_t(1) << 10) == 0)
        return std::to_string(bytes >> 10) + "K";
    return std::to_string(bytes);
}

// Format a byte size for human display: "1.5 GiB", "512 KiB", "42 B".
// Always emits a unit; one decimal digit for KiB and above.
std::string format_size_human(uint64_t bytes);

} // namespace ikafssn
