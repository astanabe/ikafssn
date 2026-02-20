#include "io/seqidlist_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>

namespace ikafssn {

// BLAST v5 binary seqidlist format:
//   Byte 0:       0x00 (null byte, marks v5 binary format)
//   Bytes 1-8:    uint64 file_size
//   Bytes 9-16:   uint64 num_ids
//   Bytes 17-20:  uint32 title_length
//   Next bytes:   title string (title_length bytes)
//   Next 1 byte:  create_date string length (uint8)
//   Next bytes:   create_date string
//   Next 8 bytes: uint64 db_vol_length (0 if not associated with DB)
//   If db_vol_length > 0:
//     Next 1 byte:  db_create_date string length (uint8)
//     Next bytes:   db_create_date string
//     Next 4 bytes: uint32 db_vol_names string length
//     Next bytes:   db_vol_names string
//   DATA SECTION (num_ids entries):
//     1 byte:  id_len (uint8). If 0xFF, next 4 bytes are uint32 actual length.
//     id_len bytes: ID string (accession)

static bool is_binary_seqidlist(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return false;

    // BLAST v5 binary seqidlist starts with a null byte (0x00)
    unsigned char first_byte = 0xFF;
    file.read(reinterpret_cast<char*>(&first_byte), 1);
    if (file.gcount() < 1) return false;

    return first_byte == 0x00;
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

// Helper: read a little-endian uint64 from a byte buffer
static uint64_t read_u64(const unsigned char* p) {
    uint64_t v = 0;
    for (int i = 7; i >= 0; i--)
        v = (v << 8) | p[i];
    return v;
}

// Helper: read a little-endian uint32 from a byte buffer
static uint32_t read_u32(const unsigned char* p) {
    uint32_t v = 0;
    for (int i = 3; i >= 0; i--)
        v = (v << 8) | p[i];
    return v;
}

static std::vector<std::string> read_binary_seqidlist(const std::string& path) {
    std::vector<std::string> result;

    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        std::fprintf(stderr, "read_seqidlist: cannot open %s\n", path.c_str());
        return result;
    }

    auto file_size = file.tellg();
    if (file_size < 21) { // minimum header size
        std::fprintf(stderr, "read_seqidlist: binary file too small: %s\n", path.c_str());
        return result;
    }

    file.seekg(0);
    std::vector<unsigned char> buf(static_cast<size_t>(file_size));
    file.read(reinterpret_cast<char*>(buf.data()), file_size);

    const unsigned char* p = buf.data();
    const unsigned char* end = p + buf.size();

    // Byte 0: null byte (already verified by is_binary_seqidlist)
    p += 1;

    // Bytes 1-8: file_size (uint64)
    if (p + 8 > end) return result;
    // uint64_t header_file_size = read_u64(p); // for validation
    p += 8;

    // Bytes 9-16: num_ids (uint64)
    if (p + 8 > end) return result;
    uint64_t num_ids = read_u64(p);
    p += 8;

    // Bytes 17-20: title_length (uint32)
    if (p + 4 > end) return result;
    uint32_t title_len = read_u32(p);
    p += 4;

    // Title string
    if (p + title_len > end) return result;
    p += title_len;

    // Create date string length (1 byte) + string
    if (p + 1 > end) return result;
    uint8_t date_len = *p;
    p += 1;
    if (p + date_len > end) return result;
    p += date_len;

    // DB volume length (uint64)
    if (p + 8 > end) return result;
    uint64_t db_vol_length = read_u64(p);
    p += 8;

    if (db_vol_length > 0) {
        // DB create date string length (1 byte) + string
        if (p + 1 > end) return result;
        uint8_t db_date_len = *p;
        p += 1;
        if (p + db_date_len > end) return result;
        p += db_date_len;

        // DB volume names string length (uint32) + string
        if (p + 4 > end) return result;
        uint32_t vol_names_len = read_u32(p);
        p += 4;
        if (p + vol_names_len > end) return result;
        p += vol_names_len;
    }

    // DATA SECTION: read num_ids entries
    result.reserve(static_cast<size_t>(num_ids));
    for (uint64_t i = 0; i < num_ids; i++) {
        if (p + 1 > end) {
            std::fprintf(stderr, "read_seqidlist: binary file truncated at entry %llu/%llu: %s\n",
                         (unsigned long long)i, (unsigned long long)num_ids, path.c_str());
            break;
        }

        uint8_t first_byte = *p;
        p += 1;

        uint32_t id_len;
        if (first_byte == 0xFF) {
            // Long ID: next 4 bytes are actual length
            if (p + 4 > end) break;
            id_len = read_u32(p);
            p += 4;
        } else {
            id_len = first_byte;
        }

        if (p + id_len > end) break;
        result.emplace_back(reinterpret_cast<const char*>(p), id_len);
        p += id_len;
    }

    return result;
}

std::vector<std::string> read_seqidlist(const std::string& path) {
    if (is_binary_seqidlist(path)) {
        return read_binary_seqidlist(path);
    }
    return read_text_seqidlist(path);
}

} // namespace ikafssn
