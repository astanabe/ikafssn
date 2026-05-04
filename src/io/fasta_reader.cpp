#include "io/fasta_reader.hpp"
#include "io/compressed_stream.hpp"
#include "io/text_simd.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

namespace ikafssn {

static void finish_record(std::vector<FastaRecord>& records,
                          std::string& cur_id, std::string& cur_seq) {
    if (!cur_id.empty()) {
        toupper_inplace_ascii(reinterpret_cast<std::uint8_t*>(cur_seq.data()),
                              cur_seq.size());
        records.push_back({std::move(cur_id), std::move(cur_seq)});
    }
    cur_id.clear();
    cur_seq.clear();
}

std::vector<FastaRecord> read_fasta_stream(std::istream& in) {
    std::vector<FastaRecord> records;
    std::string line;
    std::string cur_id;
    std::string cur_seq;

    while (std::getline(in, line)) {
        // Remove trailing \r if present (Windows line endings)
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        if (line.empty())
            continue;

        if (line[0] == '>') {
            finish_record(records, cur_id, cur_seq);
            // Extract ID: first word after '>'
            size_t start = 1;
            while (start < line.size() && std::isspace(static_cast<unsigned char>(line[start])))
                start++;
            size_t end = start;
            while (end < line.size() && !std::isspace(static_cast<unsigned char>(line[end])))
                end++;
            cur_id = line.substr(start, end - start);
        } else if (line[0] == ';') {
            // Comment line, skip
            continue;
        } else {
            cur_seq += line;
        }
    }

    finish_record(records, cur_id, cur_seq);
    return records;
}

std::vector<FastaRecord> read_fasta(const std::string& path) {
    std::string err;
    auto in = open_input_compressed(path, err);
    if (!in) {
        std::fprintf(stderr, "read_fasta: %s\n", err.c_str());
        return {};
    }
    return read_fasta_stream(*in.stream);
}

} // namespace ikafssn
