#include "util/common_init.hpp"
#include "core/version.hpp"

#include <cstdio>
#include <cstring>
#include <cstdlib>

namespace ikafssn {

std::string format_build_timestamp(const char* date, const char* time_str) {
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

bool check_version(const CliParser& cli, const char* cmd_name,
                   const char* build_date, const char* build_time) {
    if (cli.has("-version") || cli.has("--version")) {
        std::string ts = format_build_timestamp(build_date, build_time);
        std::fprintf(stderr, "%s: %s\n Package: ikafssn %s, build %s\n",
                     cmd_name, IKAFSSN_VERSION, IKAFSSN_VERSION, ts.c_str());
        return true;
    }
    return false;
}

void print_version_header(const char* cmd_name) {
    std::fprintf(stderr, "%s: %s\n\n", cmd_name, IKAFSSN_VERSION);
}

Logger make_logger(const CliParser& cli) {
    bool verbose = cli.has("-v") || cli.has("--verbose");
    return Logger(verbose ? Logger::kDebug : Logger::kInfo);
}

std::string format_size_human(uint64_t bytes) {
    auto split_decimal = [](uint64_t value, uint64_t unit) {
        return std::to_string(value / unit) + "."
             + std::to_string((value % unit) * 10 / unit);
    };
    if (bytes >= uint64_t(1) << 30) return split_decimal(bytes, uint64_t(1) << 30) + " GiB";
    if (bytes >= uint64_t(1) << 20) return split_decimal(bytes, uint64_t(1) << 20) + " MiB";
    if (bytes >= uint64_t(1) << 10) return split_decimal(bytes, uint64_t(1) << 10) + " KiB";
    return std::to_string(bytes) + " B";
}

} // namespace ikafssn
