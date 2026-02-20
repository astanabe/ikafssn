#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

// Parse a size string with optional suffix (K, M, G).
// Returns parsed size in bytes. Returns 0 on parse error.
// Examples: "8G" -> 8589934592, "512M" -> 536870912, "1024" -> 1024
inline uint64_t parse_size_string(const std::string& s) {
    if (s.empty()) return 0;

    char* end = nullptr;
    double val = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || val < 0) return 0;

    uint64_t multiplier = 1;
    if (end && *end != '\0') {
        switch (*end) {
            case 'K': case 'k': multiplier = uint64_t(1) << 10; break;
            case 'M': case 'm': multiplier = uint64_t(1) << 20; break;
            case 'G': case 'g': multiplier = uint64_t(1) << 30; break;
            default: return 0; // unknown suffix
        }
    }
    return static_cast<uint64_t>(val * static_cast<double>(multiplier));
}

} // namespace ikafssn
