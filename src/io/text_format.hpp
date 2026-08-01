#pragma once

// Append-style text formatting helpers for the output writers.
//
// The writers build rows in a std::string instead of going through the
// stream, so per-value formatting has to be cheap and independent of the
// stream's locale.  std::to_chars is both.

#include <charconv>
#include <cstdint>
#include <cstdio>
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
// This one goes through snprintf rather than std::to_chars because libc++
// marks the floating-point to_chars unavailable below macOS 13.3; "%.6g" is
// the conversion to_chars(general, 6) is specified against, and nothing here
// calls setlocale, so both give the same bytes on every target.
inline void append_double_g6(std::string& buf, double v) {
    char tmp[32];
    const int n = std::snprintf(tmp, sizeof(tmp), "%.6g", v);
    if (n > 0) buf.append(tmp, static_cast<size_t>(n));
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
