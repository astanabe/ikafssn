#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

struct ContextParam {
    bool is_ratio = false;
    double ratio = 0.0;
    uint32_t abs = 0;
};

// Parse -context value.
// If the string contains '.', interprets as ratio; otherwise as absolute.
// Returns false and sets error_msg on negative values.
inline bool parse_context(const std::string& value, ContextParam& out,
                          std::string& error_msg) {
    out = {};
    if (value.find('.') != std::string::npos) {
        out.is_ratio = true;
        out.ratio = std::stod(value);
        if (out.ratio < 0) {
            error_msg = "Error: -context ratio must be >= 0";
            return false;
        }
    } else {
        int v = std::stoi(value);
        if (v < 0) {
            error_msg = "Error: -context must be >= 0";
            return false;
        }
        out.abs = static_cast<uint32_t>(v);
    }
    return true;
}

} // namespace ikafssn
