#include "io/result_reader.hpp"
#include "io/compressed_stream.hpp"
#include "io/output_coords.hpp"

#include <algorithm>
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

// Turn the BLAST-convention coordinates `hit` holds as written into the
// internal form.  A row carrying only the two end coordinates (mode 3 without
// traceback) cannot be inverted — the other end of each interval was never
// emitted — so it keeps its ends and leaves both starts at 0.
static void to_internal_hit(OutputHit& hit, bool has_start_coords) {
    if (!has_start_coords) {
        hit.qstart = 0;
        hit.sstart = 0;
        return;
    }
    // Clamp so a malformed row cannot underflow past zero.
    const uint32_t qstart = std::max(hit.qstart, 1u);
    const uint32_t sstart = std::max(hit.sstart, 1u);
    const uint32_t send   = std::max(hit.send, 1u);
    to_internal_coords(qstart, hit.qend, sstart, send,
                       hit.qlen, hit.sstrand == '-',
                       hit.qstart, hit.qend, hit.sstart, hit.send);
}

// Returns true if the sseqid is a skip-marker sentinel emitted by
// result_writer for queries that were not searched.  Two prefixes
// exist: "*SKIPPED:" (server-produced) and "*FAILED:" (client-produced
// async-job failure).  Both must be silently dropped on read.
static bool is_skip_sseqid(const std::string& s) {
    if (s.size() >= 9 && s.compare(0, 9, "*SKIPPED:") == 0) return true;
    if (s.size() >= 8 && s.compare(0, 8, "*FAILED:") == 0)  return true;
    return false;
}

// Parse a data line using the column map built from the header.
// Returns true if the line was consumed; out-param `is_skip_marker` is set
// when the line is a skip sentinel (caller should drop it without reporting
// an error).
static bool parse_line_with_map(
    const std::string& line,
    const std::unordered_map<std::string, size_t>& cmap,
    OutputHit& hit,
    bool& is_skip_marker) {
    is_skip_marker = false;
    auto fields = split_tabs(line);

    // Required columns: qseqid, sseqid, sstrand
    static const std::string empty;
    const auto& qid = field_str(fields, cmap, "qseqid", empty);
    const auto& acc = field_str(fields, cmap, "sseqid", empty);
    if (qid.empty() || acc.empty()) return false;

    if (is_skip_sseqid(acc)) {
        // Drop skip-marker rows silently — readers consuming this file
        // are looking for real hits. The skip information is preserved
        // in the original file for human inspection.
        is_skip_marker = true;
        return true;
    }

    char sstrand = field_char(fields, cmap, "sstrand", '\0');
    if (sstrand != '+' && sstrand != '-') return false;

    try {
        hit.qseqid = qid;
        hit.sseqid = acc;
        hit.sstrand = sstrand;

        hit.qstart = field_u32(fields, cmap, "qstart");
        hit.qend = field_u32(fields, cmap, "qend");
        hit.qlen = field_u32(fields, cmap, "qlen");
        hit.sstart = field_u32(fields, cmap, "sstart");
        hit.send = field_u32(fields, cmap, "send");
        hit.slen = field_u32(fields, cmap, "slen");

        hit.coverscore = field_u32(fields, cmap, "coverscore");
        hit.chainscore = field_u32(fields, cmap, "chainscore");
        hit.alnscore = field_i32(fields, cmap, "alnscore");
        hit.ppositive = field_dbl(fields, cmap, "ppositive");
        hit.npositive = field_u32(fields, cmap, "npositive");
        hit.nnegative = field_u32(fields, cmap, "nnegative");

        hit.cigar = field_str(fields, cmap, "cigar", empty);
        hit.qseq = field_str(fields, cmap, "qseq", empty);
        hit.sseq = field_str(fields, cmap, "sseq", empty);

        hit.volume = field_u16(fields, cmap, "volume");

        to_internal_hit(hit, cmap.count("qstart") != 0 && cmap.count("sstart") != 0);
    } catch (...) {
        return false;
    }

    return true;
}

// Field-count-based parser, used as a fallback when no `# ` header
// line is present in the input TSV.
static bool parse_line_no_header(const std::string& line, OutputHit& hit) {
    auto fields = split_tabs(line);

    if (fields.size() < 7) return false;

    hit.qseqid = fields[0];
    hit.sseqid = fields[1];

    if (fields[2].size() != 1 || (fields[2][0] != '+' && fields[2][0] != '-'))
        return false;
    hit.sstrand = fields[2][0];

    bool has_start_coords = true;
    try {
        if (fields.size() >= 20) {
            hit.qstart = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.qend = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.qlen = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.sstart = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.send = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.slen = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.chainscore = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.alnscore = static_cast<int32_t>(std::stol(fields[11]));
            hit.ppositive = std::stod(fields[12]);
            hit.npositive = static_cast<uint32_t>(std::stoul(fields[13]));
            hit.nnegative = static_cast<uint32_t>(std::stoul(fields[14]));
            hit.cigar = fields[15];
            hit.qseq = fields[16];
            hit.sseq = fields[17];
            hit.volume = static_cast<uint16_t>(std::stoul(fields[18]));
        } else if (fields.size() >= 13) {
            hit.qstart = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.qend = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.qlen = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.sstart = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.send = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.slen = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.chainscore = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.alnscore = static_cast<int32_t>(std::stol(fields[11]));
            hit.volume = static_cast<uint16_t>(std::stoul(fields[12]));
        } else if (fields.size() >= 12) {
            hit.qstart = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.qend = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.qlen = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.sstart = static_cast<uint32_t>(std::stoul(fields[6]));
            hit.send = static_cast<uint32_t>(std::stoul(fields[7]));
            hit.slen = static_cast<uint32_t>(std::stoul(fields[8]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[9]));
            hit.chainscore = static_cast<uint32_t>(std::stoul(fields[10]));
            hit.volume = static_cast<uint16_t>(std::stoul(fields[11]));
        } else if (fields.size() >= 7) {
            hit.qstart = 0;
            hit.qend = 0;
            hit.sstart = 0;
            hit.send = 0;
            hit.qlen = static_cast<uint32_t>(std::stoul(fields[3]));
            hit.slen = static_cast<uint32_t>(std::stoul(fields[4]));
            hit.coverscore = static_cast<uint32_t>(std::stoul(fields[5]));
            hit.chainscore = 0;
            hit.volume = static_cast<uint16_t>(std::stoul(fields[6]));
            has_start_coords = false;
        }
        to_internal_hit(hit, has_start_coords);
    } catch (...) {
        return false;
    }

    return true;
}

std::vector<OutputHit> read_results_tsv(std::istream& in) {
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

    // Determine parse strategy: header-based or no-header fallback
    bool use_header = false;
    std::unordered_map<std::string, size_t> cmap;

    if (!last_header.empty() && last_header.size() > 2 && last_header[1] == ' ') {
        // Strip "# " prefix and build column map
        auto header_body = last_header.substr(2);
        cmap = build_column_map(header_body);
        // Require at least the three mandatory columns
        if (cmap.count("qseqid") && cmap.count("sseqid") && cmap.count("sstrand")) {
            use_header = true;
        }
    }

    for (const auto& [lnum, dline] : data_lines) {
        OutputHit hit;
        bool ok;
        bool is_skip_marker = false;
        if (use_header) {
            ok = parse_line_with_map(dline, cmap, hit, is_skip_marker);
        } else {
            ok = parse_line_no_header(dline, hit);
        }

        if (is_skip_marker) {
            // Silently drop skip sentinel rows.
            continue;
        }

        if (ok) {
            results.push_back(std::move(hit));
        } else {
            std::fprintf(stderr, "result_reader: skipping invalid line %d\n", lnum);
        }
    }

    return results;
}

std::vector<OutputHit> read_results_tsv(const std::string& path) {
    std::string err;
    auto in = open_input_compressed(path, err);
    if (!in) {
        std::fprintf(stderr, "result_reader: %s\n", err.c_str());
        return {};
    }
    return read_results_tsv(*in.stream);
}

} // namespace ikafssn
