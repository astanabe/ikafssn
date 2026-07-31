#pragma once

// Append-style text formatting helpers for the output writers.
//
// The writers build rows in a std::string instead of going through the
// stream, so per-value formatting has to be cheap and independent of the
// stream's locale.  std::to_chars is both.

#include <charconv>
#include <cstdint>
#include <string>
#include <string_view>

namespace ikafssn {

inline void append_uint(std::string& buf, uint64_t v) {
    char tmp[24];
    auto r = std::to_chars(tmp, tmp + sizeof(tmp), v);
    buf.append(tmp, static_cast<size_t>(r.ptr - tmp));
}

inline void append_int(std::string& buf, int64_t v) {
    char tmp[24];
    auto r = std::to_chars(tmp, tmp + sizeof(tmp), v);
    buf.append(tmp, static_cast<size_t>(r.ptr - tmp));
}

// Six significant digits, matching what std::ostream's default precision
// produces for a double (its num_put always formats in the C locale).
inline void append_double_g6(std::string& buf, double v) {
    char tmp[32];
    auto r = std::to_chars(tmp, tmp + sizeof(tmp), v,
                           std::chars_format::general, 6);
    buf.append(tmp, static_cast<size_t>(r.ptr - tmp));
}

inline void append_str(std::string& buf, std::string_view s) {
    buf.append(s.data(), s.size());
}

// Append `s` as a quoted JSON string, escaping the characters that JSON
// forbids raw.  Control characters other than \n / \r / \t are passed
// through unchanged; the values written here never contain them.
inline void json_escape_into(std::string& buf, std::string_view s) {
    buf.push_back('"');
    for (char c : s) {
        switch (c) {
            case '"':  buf.append("\\\""); break;
            case '\\': buf.append("\\\\"); break;
            case '\n': buf.append("\\n");  break;
            case '\r': buf.append("\\r");  break;
            case '\t': buf.append("\\t");  break;
            default:   buf.push_back(c);   break;
        }
    }
    buf.push_back('"');
}

} // namespace ikafssn
