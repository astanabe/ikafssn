#include "io/seqidlist_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

namespace ikafssn {

// Binary seqidlist magic: BLAST binary seqidlist files start with specific bytes.
// The blastdb_aliastool binary format uses ASN.1 encoding.
// We detect binary vs text by checking if first bytes are printable text or not.

static bool is_binary_seqidlist(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return false;

    // Read first 4 bytes
    unsigned char buf[4] = {};
    file.read(reinterpret_cast<char*>(buf), 4);
    auto n = file.gcount();
    if (n == 0) return false;

    // Binary seqidlist files produced by blastdb_aliastool start with
    // ASN.1 BER encoding. The first byte is typically 0x30 (SEQUENCE)
    // or other non-printable values. Text files start with printable
    // chars (accessions) or '>' or '#'.
    for (int i = 0; i < n; i++) {
        unsigned char c = buf[i];
        if (c == '\n' || c == '\r' || c == '\t') continue;
        if (c == '>' || c == '#') return false; // text markers
        if (c < 0x20 || c > 0x7E) return true;  // non-printable = binary
    }
    return false;
}

static std::vector<std::string> read_text_seqidlist(const std::string& path) {
    std::vector<std::string> result;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::fprintf(stderr, "read_seqidlist: cannot open %s\n", path.c_str());
        return result;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Remove trailing \r
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        // Skip empty / whitespace-only lines
        if (line.empty()) continue;
        bool all_space = true;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                all_space = false;
                break;
            }
        }
        if (all_space) continue;

        // Skip comment lines
        if (line[0] == '#') continue;

        // Trim leading '>' if present
        size_t start = 0;
        if (line[0] == '>') start = 1;

        // Trim leading/trailing whitespace
        while (start < line.size() && std::isspace(static_cast<unsigned char>(line[start])))
            start++;
        size_t end = line.size();
        while (end > start && std::isspace(static_cast<unsigned char>(line[end - 1])))
            end--;

        if (start < end) {
            // Take first whitespace-delimited token as accession
            size_t tok_end = start;
            while (tok_end < end && !std::isspace(static_cast<unsigned char>(line[tok_end])))
                tok_end++;
            result.emplace_back(line.substr(start, tok_end - start));
        }
    }

    return result;
}

static std::vector<std::string> read_binary_seqidlist(const std::string& path) {
    // Binary seqidlist from blastdb_aliastool stores GI numbers or Seq-ids
    // in ASN.1 format. For our purposes, we read the raw GI/accession data.
    // This is a simplified parser for the common case.
    // If the binary format cannot be parsed, fall back to warning.
    std::fprintf(stderr, "read_seqidlist: binary seqidlist format detected for %s\n", path.c_str());
    std::fprintf(stderr, "read_seqidlist: binary format parsing not yet fully implemented, "
                         "please use text format\n");
    return {};
}

std::vector<std::string> read_seqidlist(const std::string& path) {
    if (is_binary_seqidlist(path)) {
        return read_binary_seqidlist(path);
    }
    return read_text_seqidlist(path);
}

} // namespace ikafssn
