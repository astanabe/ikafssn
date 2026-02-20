#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

namespace ikafssn {

struct AmbiguityEntry {
    uint32_t position;   // 0-based base offset
    uint32_t run_length; // number of consecutive ambiguous bases
    uint8_t  ncbi4na;    // ncbi4na value (0-15)
};

class AmbiguityParser {
public:
    // Parse ambiguity data and return entries sorted by position.
    // ambig_data: pointer to ambiguity data
    // ambig_bytes: byte length of ambiguity data
    static std::vector<AmbiguityEntry> parse(const char* ambig_data, int ambig_bytes) {
        std::vector<AmbiguityEntry> entries;
        if (!ambig_data || ambig_bytes < 4) return entries;

        const auto* data = reinterpret_cast<const uint8_t*>(ambig_data);

        // Read header word (big-endian)
        uint32_t header = read_be32(data);
        bool new_format = (header & 0x80000000u) != 0;

        if (new_format) {
            // New format: 8 bytes/entry
            // Header bits 30-0 = number of Int4 words (= num_entries * 2)
            uint32_t num_words = header & 0x7FFFFFFFu;
            uint32_t num_entries = num_words / 2;
            int required = 4 + static_cast<int>(num_entries) * 8;
            if (ambig_bytes < required) return entries;

            entries.reserve(num_entries);
            for (uint32_t i = 0; i < num_entries; i++) {
                const uint8_t* p = data + 4 + i * 8;
                uint32_t w0 = read_be32(p);
                uint32_t w1 = read_be32(p + 4);

                AmbiguityEntry e;
                e.ncbi4na    = static_cast<uint8_t>((w0 >> 28) & 0xF);
                e.run_length = ((w0 >> 16) & 0xFFF) + 1;
                e.position   = w1;
                entries.push_back(e);
            }
        } else {
            // Old format: 4 bytes/entry
            // Header bits 30-0 = number of entries
            uint32_t num_entries = header & 0x7FFFFFFFu;
            int required = 4 + static_cast<int>(num_entries) * 4;
            if (ambig_bytes < required) return entries;

            entries.reserve(num_entries);
            for (uint32_t i = 0; i < num_entries; i++) {
                uint32_t word = read_be32(data + 4 + i * 4);

                AmbiguityEntry e;
                e.ncbi4na    = static_cast<uint8_t>((word >> 28) & 0xF);
                e.run_length = ((word >> 24) & 0xF) + 1;
                e.position   = word & 0x00FFFFFFu;
                entries.push_back(e);
            }
        }

        // Sort by position (normally already sorted, but guarantee it)
        std::sort(entries.begin(), entries.end(),
            [](const AmbiguityEntry& a, const AmbiguityEntry& b) {
                return a.position < b.position;
            });

        return entries;
    }

private:
    static uint32_t read_be32(const uint8_t* p) {
        return (static_cast<uint32_t>(p[0]) << 24) |
               (static_cast<uint32_t>(p[1]) << 16) |
               (static_cast<uint32_t>(p[2]) <<  8) |
               (static_cast<uint32_t>(p[3]));
    }
};

} // namespace ikafssn
