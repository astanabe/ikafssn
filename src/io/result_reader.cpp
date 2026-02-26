#include "io/result_reader.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

namespace ikafssn {

static bool parse_line(const std::string& line, OutputHit& hit) {
    // Split into tab-separated fields
    std::vector<std::string> fields;
    std::string::size_type start = 0;
    while (true) {
        auto pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }

    if (fields.size() < 7) return false;

    hit.query_id = fields[0];
    hit.accession = fields[1];

    if (fields[2].size() != 1 || (fields[2][0] != '+' && fields[2][0] != '-'))
        return false;
    hit.strand = fields[2][0];

    try {
        // Supported formats (with q_len/s_len):
        //   7 fields (mode 1):  query_id accession strand q_len s_len stage1_score volume
        //  12 fields (mode 2):  query_id accession strand q_start q_end q_len s_start s_end s_len stage1_score chainscore volume
        //  13 fields (mode 3, no traceback): + alnscore
        //  20 fields (mode 3, traceback): + alnscore pident nident nmismatch cigar q_seq s_seq

        if (fields.size() >= 20) {
            // Mode 3, traceback (20 fields)
            hit.q_start = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.q_end = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.s_start = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.s_end = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.stage1_score = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.score = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.alnscore = static_cast<int32_t>(std::stol(fields[11]));
            hit.pident = std::stod(fields[12]);
            hit.nident = static_cast<uint32_t>(std::stoul(fields[13]));
            hit.nmismatch = static_cast<uint32_t>(std::stoul(fields[14]));
            hit.cigar = fields[15];
            hit.q_seq = fields[16];
            hit.s_seq = fields[17];
            hit.volume = static_cast<uint16_t>(std::stoul(fields[18]));
        } else if (fields.size() >= 13) {
            // Mode 3, no traceback (13 fields)
            hit.q_start = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.q_end = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.s_start = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.s_end = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.stage1_score = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.score = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.alnscore = static_cast<int32_t>(std::stol(fields[11]));
            hit.volume = static_cast<uint16_t>(std::stoul(fields[12]));
        } else if (fields.size() >= 12) {
            // Mode 2 (12 fields)
            hit.q_start = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.q_end = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.s_start = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.s_end = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.stage1_score = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.score = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.volume = static_cast<uint16_t>(std::stoul(fields[11]));
        } else if (fields.size() >= 7) {
            // Mode 1 (7 fields)
            hit.q_start = 0;
            hit.q_end = 0;
            hit.s_start = 0;
            hit.s_end = 0;
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.stage1_score = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.score = 0;
            hit.volume = static_cast<uint16_t>(std::stoul(fields[6]));
        }
    } catch (...) {
        return false;
    }

    return true;
}

std::vector<OutputHit> read_results_tab(std::istream& in) {
    std::vector<OutputHit> results;
    std::string line;
    int line_num = 0;

    while (std::getline(in, line)) {
        line_num++;
        // Remove trailing \r
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        // Skip empty lines and comment/header lines
        if (line.empty() || line[0] == '#')
            continue;

        OutputHit hit;
        if (parse_line(line, hit)) {
            results.push_back(std::move(hit));
        } else {
            std::fprintf(stderr, "result_reader: skipping invalid line %d\n", line_num);
        }
    }

    return results;
}

std::vector<OutputHit> read_results_tab(const std::string& path) {
    if (path == "-") {
        return read_results_tab(std::cin);
    }

    std::ifstream file(path);
    if (!file.is_open()) {
        std::fprintf(stderr, "result_reader: cannot open %s\n", path.c_str());
        return {};
    }
    return read_results_tab(file);
}

} // namespace ikafssn
