#include "io/result_writer.hpp"

#include <map>

namespace ikafssn {

static const char* stage1_score_name(uint8_t stage1_score_type) {
    return (stage1_score_type == 2) ? "matchscore" : "coverscore";
}

void write_results_tab(std::ostream& out,
                       const std::vector<OutputHit>& hits,
                       uint8_t mode,
                       uint8_t stage1_score_type,
                       bool stage3_traceback) {
    const char* s1name = stage1_score_name(stage1_score_type);

    if (mode == 1) {
        out << "# qseqid\tsseqid\tsstrand\tqlen\tslen\t" << s1name << "\tvolume\n";
        for (const auto& h : hits) {
            out << h.qseqid << '\t'
                << h.sseqid << '\t'
                << h.sstrand << '\t'
                << h.qlen << '\t'
                << h.slen << '\t'
                << ((stage1_score_type == 2) ? h.matchscore : h.coverscore) << '\t'
                << h.volume << '\n';
        }
    } else if (mode == 3 && stage3_traceback) {
        out << "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\t"
            << s1name << "\tchainscore\talnscore\tpident\tnident\tmismatch\tcigar\tqseq\tsseq\tvolume\n";
        for (const auto& h : hits) {
            out << h.qseqid << '\t'
                << h.sseqid << '\t'
                << h.sstrand << '\t'
                << h.qstart << '\t'
                << h.qend << '\t'
                << h.qlen << '\t'
                << h.sstart << '\t'
                << h.send << '\t'
                << h.slen << '\t'
                << ((stage1_score_type == 2) ? h.matchscore : h.coverscore) << '\t'
                << h.chainscore << '\t'
                << h.alnscore << '\t'
                << h.pident << '\t'
                << h.nident << '\t'
                << h.mismatch << '\t'
                << h.cigar << '\t'
                << h.qseq << '\t'
                << h.sseq << '\t'
                << h.volume << '\n';
        }
    } else if (mode == 3) {
        out << "# qseqid\tsseqid\tsstrand\tqend\tqlen\tsend\tslen\t"
            << s1name << "\tchainscore\talnscore\tvolume\n";
        for (const auto& h : hits) {
            out << h.qseqid << '\t'
                << h.sseqid << '\t'
                << h.sstrand << '\t'
                << h.qend << '\t'
                << h.qlen << '\t'
                << h.send << '\t'
                << h.slen << '\t'
                << ((stage1_score_type == 2) ? h.matchscore : h.coverscore) << '\t'
                << h.chainscore << '\t'
                << h.alnscore << '\t'
                << h.volume << '\n';
        }
    } else {
        out << "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\t"
            << s1name << "\tchainscore\tvolume\n";
        for (const auto& h : hits) {
            out << h.qseqid << '\t'
                << h.sseqid << '\t'
                << h.sstrand << '\t'
                << h.qstart << '\t'
                << h.qend << '\t'
                << h.qlen << '\t'
                << h.sstart << '\t'
                << h.send << '\t'
                << h.slen << '\t'
                << ((stage1_score_type == 2) ? h.matchscore : h.coverscore) << '\t'
                << h.chainscore << '\t'
                << h.volume << '\n';
        }
    }
}

static void json_escape(std::ostream& out, const std::string& s) {
    out << '"';
    for (char c : s) {
        switch (c) {
            case '"':  out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n";  break;
            case '\r': out << "\\r";  break;
            case '\t': out << "\\t";  break;
            default:   out << c;      break;
        }
    }
    out << '"';
}

// Inner helper: write per-query JSON objects.
// If is_fragment is true, objects are emitted without outer wrapper and
// with trailing comma after each (caller handles the last comma).
static void write_results_json_inner(std::ostream& out,
                                      const std::vector<OutputHit>& hits,
                                      uint8_t mode,
                                      uint8_t stage1_score_type,
                                      bool stage3_traceback,
                                      bool is_fragment) {
    const char* s1name = stage1_score_name(stage1_score_type);

    // Group hits by qseqid (preserve order of first appearance)
    std::vector<std::string> query_order;
    std::map<std::string, std::vector<const OutputHit*>> by_query;
    for (const auto& h : hits) {
        if (by_query.find(h.qseqid) == by_query.end()) {
            query_order.push_back(h.qseqid);
        }
        by_query[h.qseqid].push_back(&h);
    }

    for (size_t qi = 0; qi < query_order.size(); qi++) {
        const auto& qid = query_order[qi];
        const auto& qhits = by_query[qid];
        out << "    {\n      \"qseqid\": ";
        json_escape(out, qid);
        out << ",\n      \"hits\": [\n";
        for (size_t hi = 0; hi < qhits.size(); hi++) {
            const auto* h = qhits[hi];
            out << "        {\n";
            out << "          \"sseqid\": "; json_escape(out, h->sseqid); out << ",\n";
            out << "          \"sstrand\": \"" << h->sstrand << "\",\n";
            if (mode == 2 || (mode == 3 && stage3_traceback)) {
                out << "          \"qstart\": " << h->qstart << ",\n";
                out << "          \"qend\": " << h->qend << ",\n";
            } else if (mode == 3) {
                out << "          \"qend\": " << h->qend << ",\n";
            }
            out << "          \"qlen\": " << h->qlen << ",\n";
            if (mode == 2 || (mode == 3 && stage3_traceback)) {
                out << "          \"sstart\": " << h->sstart << ",\n";
                out << "          \"send\": " << h->send << ",\n";
            } else if (mode == 3) {
                out << "          \"send\": " << h->send << ",\n";
            }
            out << "          \"slen\": " << h->slen << ",\n";
            out << "          \"" << s1name << "\": " << ((stage1_score_type == 2) ? h->matchscore : h->coverscore) << ",\n";
            if (mode != 1) {
                out << "          \"chainscore\": " << h->chainscore << ",\n";
            }
            if (mode == 3) {
                out << "          \"alnscore\": " << h->alnscore << ",\n";
                if (stage3_traceback) {
                    out << "          \"pident\": " << h->pident << ",\n";
                    out << "          \"nident\": " << h->nident << ",\n";
                    out << "          \"mismatch\": " << h->mismatch << ",\n";
                    out << "          \"cigar\": "; json_escape(out, h->cigar); out << ",\n";
                    out << "          \"qseq\": "; json_escape(out, h->qseq); out << ",\n";
                    out << "          \"sseq\": "; json_escape(out, h->sseq); out << ",\n";
                }
            }
            out << "          \"volume\": " << h->volume << "\n";
            out << "        }";
            if (hi + 1 < qhits.size()) out << ',';
            out << '\n';
        }
        out << "      ]\n    }";
        if (is_fragment || qi + 1 < query_order.size()) out << ',';
        out << '\n';
    }
}

void write_results_json(std::ostream& out,
                        const std::vector<OutputHit>& hits,
                        uint8_t mode,
                        uint8_t stage1_score_type,
                        bool stage3_traceback) {
    out << "{\n  \"results\": [\n";
    write_results_json_inner(out, hits, mode, stage1_score_type,
                              stage3_traceback, false);
    out << "  ]\n}\n";
}

void write_results_json_fragment(std::ostream& out,
                                  const std::vector<OutputHit>& hits,
                                  uint8_t mode,
                                  uint8_t stage1_score_type,
                                  bool stage3_traceback) {
    write_results_json_inner(out, hits, mode, stage1_score_type,
                              stage3_traceback, true);
}

void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt,
                   uint8_t mode,
                   uint8_t stage1_score_type,
                   bool stage3_traceback) {
    switch (fmt) {
        case OutputFormat::kTab:
            write_results_tab(out, hits, mode, stage1_score_type, stage3_traceback);
            break;
        case OutputFormat::kJson:
            write_results_json(out, hits, mode, stage1_score_type, stage3_traceback);
            break;
        case OutputFormat::kSam:
        case OutputFormat::kBam:
            // SAM/BAM handled separately via sam_writer
            break;
    }
}

bool validate_output_format(OutputFormat fmt, uint8_t mode, bool traceback,
                            const std::string& output_path,
                            std::string& error_msg) {
    if ((fmt == OutputFormat::kSam || fmt == OutputFormat::kBam) &&
        (mode != 3 || !traceback)) {
        error_msg = "Error: SAM/BAM output requires -mode 3 and -stage3_traceback 1";
        return false;
    }
    if (fmt == OutputFormat::kBam && output_path.empty()) {
        error_msg = "Error: BAM output requires -o <path>";
        return false;
    }
    return true;
}

bool parse_output_format(const std::string& str, OutputFormat& out,
                         std::string& error_msg) {
    if (str == "tab") {
        out = OutputFormat::kTab;
    } else if (str == "json") {
        out = OutputFormat::kJson;
    } else if (str == "sam") {
        out = OutputFormat::kSam;
    } else if (str == "bam") {
        out = OutputFormat::kBam;
    } else {
        error_msg = "Error: unknown output format '" + str + "'";
        return false;
    }
    return true;
}

} // namespace ikafssn
