#include "io/result_reader.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace ikafssn {

// Split a string by tab delimiter
static std::vector<std::string> split_tabs(const std::string& line) {
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
    return fields;
}

// Build column-name-to-index map from a header line (leading "# " stripped).
static std::unordered_map<std::string, size_t>
build_column_map(const std::string& header_body) {
    auto cols = split_tabs(header_body);
    std::unordered_map<std::string, size_t> m;
    for (size_t i = 0; i < cols.size(); i++) {
        m[cols[i]] = i;
    }
    return m;
}

// Helper: get a field as string; returns fallback if column absent or index out of range.
static const std::string& field_str(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& cmap,
    const std::string& name,
    const std::string& fallback) {
    auto it = cmap.find(name);
    if (it == cmap.end() || it->second >= fields.size()) return fallback;
    return fields[it->second];
}

static uint32_t field_u32(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& cmap,
    const std::string& name) {
    auto it = cmap.find(name);
    if (it == cmap.end() || it->second >= fields.size()) return 0;
    return static_cast<uint32_t>(std::stoul(fields[it->second]));
}

static int32_t field_i32(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& cmap,
    const std::string& name) {
    auto it = cmap.find(name);
    if (it == cmap.end() || it->second >= fields.size()) return 0;
    return static_cast<int32_t>(std::stol(fields[it->second]));
}

static uint16_t field_u16(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& cmap,
    const std::string& name) {
    auto it = cmap.find(name);
    if (it == cmap.end() || it->second >= fields.size()) return 0;
    return static_cast<uint16_t>(std::stoul(fields[it->second]));
}

static double field_dbl(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& cmap,
    const std::string& name) {
    auto it = cmap.find(name);
    if (it == cmap.end() || it->second >= fields.size()) return 0.0;
    return std::stod(fields[it->second]);
}

static char field_char(
    const std::vector<std::string>& fields,
    const std::unordered_map<std::string, size_t>& cmap,
    const std::string& name,
    char fallback) {
    auto it = cmap.find(name);
    if (it == cmap.end() || it->second >= fields.size()) return fallback;
    const auto& s = fields[it->second];
    return s.empty() ? fallback : s[0];
}

// Parse a data line using the column map built from the header.
static bool parse_line_with_map(
    const std::string& line,
    const std::unordered_map<std::string, size_t>& cmap,
    OutputHit& hit) {
    auto fields = split_tabs(line);

    // Required columns: query_id, accession, strand
    static const std::string empty;
    const auto& qid = field_str(fields, cmap, "query_id", empty);
    const auto& acc = field_str(fields, cmap, "accession", empty);
    if (qid.empty() || acc.empty()) return false;

    char strand = field_char(fields, cmap, "strand", '\0');
    if (strand != '+' && strand != '-') return false;

    try {
        hit.query_id = qid;
        hit.accession = acc;
        hit.strand = strand;

        hit.q_start = field_u32(fields, cmap, "q_start");
        hit.q_end = field_u32(fields, cmap, "q_end");
        hit.q_length = field_u32(fields, cmap, "q_len");
        hit.s_start = field_u32(fields, cmap, "s_start");
        hit.s_end = field_u32(fields, cmap, "s_end");
        hit.s_length = field_u32(fields, cmap, "s_len");

        hit.coverscore = field_u32(fields, cmap, "coverscore");
        hit.matchscore = field_u32(fields, cmap, "matchscore");
        hit.chainscore = field_u32(fields, cmap, "chainscore");
        hit.alnscore = field_i32(fields, cmap, "alnscore");
        hit.pident = field_dbl(fields, cmap, "pident");
        hit.nident = field_u32(fields, cmap, "nident");
        hit.nmismatch = field_u32(fields, cmap, "nmismatch");

        hit.cigar = field_str(fields, cmap, "cigar", empty);
        hit.q_seq = field_str(fields, cmap, "q_seq", empty);
        hit.s_seq = field_str(fields, cmap, "s_seq", empty);

        hit.volume = field_u16(fields, cmap, "volume");
    } catch (...) {
        return false;
    }

    return true;
}

// Legacy field-count-based parser (fallback when no header line is present).
static bool parse_line_legacy(const std::string& line, OutputHit& hit) {
    auto fields = split_tabs(line);

    if (fields.size() < 7) return false;

    hit.query_id = fields[0];
    hit.accession = fields[1];

    if (fields[2].size() != 1 || (fields[2][0] != '+' && fields[2][0] != '-'))
        return false;
    hit.strand = fields[2][0];

    try {
        if (fields.size() >= 20) {
            hit.q_start = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.q_end = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.s_start = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.s_end = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.chainscore = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.alnscore = static_cast<int32_t>(std::stol(fields[11]));
            hit.pident = std::stod(fields[12]);
            hit.nident = static_cast<uint32_t>(std::stoul(fields[13]));
            hit.nmismatch = static_cast<uint32_t>(std::stoul(fields[14]));
            hit.cigar = fields[15];
            hit.q_seq = fields[16];
            hit.s_seq = fields[17];
            hit.volume = static_cast<uint16_t>(std::stoul(fields[18]));
        } else if (fields.size() >= 13) {
            hit.q_start = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.q_end = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.s_start = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.s_end = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.chainscore = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.alnscore = static_cast<int32_t>(std::stol(fields[11]));
            hit.volume = static_cast<uint16_t>(std::stoul(fields[12]));
        } else if (fields.size() >= 12) {
            hit.q_start = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.q_end = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.s_start = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.s_end = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.chainscore = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.volume = static_cast<uint16_t>(std::stoul(fields[11]));
        } else if (fields.size() >= 7) {
            hit.q_start = 0;
            hit.q_end = 0;
            hit.s_start = 0;
            hit.s_end = 0;
            hit.q_length = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.s_length = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.chainscore = 0;
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

    // Collect all lines; remember the last header line for column map.
    std::string last_header;
    std::vector<std::pair<int, std::string>> data_lines; // (line_num, line)

    while (std::getline(in, line)) {
        line_num++;
        // Remove trailing \r
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        if (line.empty()) continue;

        if (line[0] == '#') {
            last_header = line;
            continue;
        }

        data_lines.emplace_back(line_num, std::move(line));
    }

    // Determine parse strategy: header-based or legacy fallback
    bool use_header = false;
    std::unordered_map<std::string, size_t> cmap;

    if (!last_header.empty() && last_header.size() > 2 && last_header[1] == ' ') {
        // Strip "# " prefix and build column map
        auto header_body = last_header.substr(2);
        cmap = build_column_map(header_body);
        // Require at least the three mandatory columns
        if (cmap.count("query_id") && cmap.count("accession") && cmap.count("strand")) {
            use_header = true;
        }
    }

    for (const auto& [lnum, dline] : data_lines) {
        OutputHit hit;
        bool ok;
        if (use_header) {
            ok = parse_line_with_map(dline, cmap, hit);
        } else {
            ok = parse_line_legacy(dline, hit);
        }

        if (ok) {
            results.push_back(std::move(hit));
        } else {
            std::fprintf(stderr, "result_reader: skipping invalid line %d\n", lnum);
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
