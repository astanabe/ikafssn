#pragma once

// Helper for splitting multi-accession strings.
//
// `BlastDbReader::get_accession()` and `KsxReader::accession()` may
// return a string containing multiple accessions joined by '\x01'
// (BLAST's native multi-defline separator).  Output writers emit the
// raw '\x01'-joined form so downstream consumers can split on '\x01'
// themselves; this header only provides the split helper used by
// internal call sites that need to walk individual tokens (e.g. the
// OID filter's reverse map, and the leading token ikafssnretrieve looks
// an accession up by and puts in its FASTA deflines).

#include <string_view>
#include <vector>

namespace ikafssn {

inline constexpr char kAccessionSeparator = '\x01';

// Split a '\x01'-joined accession string into individual tokens
// (zero-copy std::string_view).  Empty tokens are dropped.
inline std::vector<std::string_view>
split_accessions(std::string_view s) {
    std::vector<std::string_view> out;
    out.reserve(1);
    std::size_t start = 0;
    while (start <= s.size()) {
        std::size_t end = s.find(kAccessionSeparator, start);
        if (end == std::string_view::npos) end = s.size();
        if (end > start) out.emplace_back(s.data() + start, end - start);
        if (end == s.size()) break;
        start = end + 1;
    }
    return out;
}

} // namespace ikafssn
