#include "io/fasta_reader.hpp"

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
        // Convert sequence to uppercase
        for (auto& c : cur_seq)
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        records.push_back({std::move(cur_id), std::move(cur_seq)});
    }
    cur_id.clear();
    cur_seq.clear();
}

static std::vector<FastaRecord> read_fasta_stream(std::istream& in) {
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
    if (path == "-") {
        return read_fasta_stream(std::cin);
    }

    std::ifstream file(path);
    if (!file.is_open()) {
        std::fprintf(stderr, "read_fasta: cannot open %s\n", path.c_str());
        return {};
    }
    return read_fasta_stream(file);
}

} // namespace ikafssn
