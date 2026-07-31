#include "io/result_writer.hpp"
#include "io/compressed_stream.hpp"
#include "io/text_format.hpp"
#include "protocol/messages.hpp"

#include <map>
#include <string>

namespace ikafssn {

// Stage 1 score is always coverscore (the number of distinct query k-mers
// matching the sequence).
static constexpr const char* kStage1ScoreName = "coverscore";

namespace {

// Rows are formatted into a std::string and handed to the stream in blocks
// of at least this size.
constexpr std::size_t kBlockSize = 1u << 20;

inline void write_block(std::ostream& out, std::string& buf) {
    if (!buf.empty()) {
        out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
        buf.clear();
    }
}

inline void write_block_if_full(std::ostream& out, std::string& buf) {
    if (buf.size() >= kBlockSize) write_block(out, buf);
}

} // namespace

void write_results_tsv(std::ostream& out,
                       const std::vector<OutputHit>& hits,
                       uint8_t mode,
                       bool stage3_traceback) {
    const char* s1name = kStage1ScoreName;
    std::string buf;
    buf.reserve(kBlockSize + 4096);

    if (mode == 1) {
        buf.append("# qseqid\tsseqid\tsstrand\tqlen\tslen\t");
        buf.append(s1name);
        buf.append("\tvolume\n");
        for (const auto& h : hits) {
            if (h.skip_reason != 0) {
                append_str(buf, h.qseqid);                              buf.push_back('\t');
                append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                buf.append("\t*\t");
                append_uint(buf, h.qlen);
                buf.append("\t0\t0\t0\n");
                write_block_if_full(out, buf);
                continue;
            }
            append_str(buf, h.qseqid);        buf.push_back('\t');
            append_str(buf, h.sseqid);        buf.push_back('\t');
            buf.push_back(h.sstrand);         buf.push_back('\t');
            append_uint(buf, h.qlen);         buf.push_back('\t');
            append_uint(buf, h.slen);         buf.push_back('\t');
            append_uint(buf, h.coverscore);   buf.push_back('\t');
            append_uint(buf, h.volume);       buf.push_back('\n');
            write_block_if_full(out, buf);
        }
    } else if (mode == 3 && stage3_traceback) {
        buf.append("# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\t");
        buf.append(s1name);
        buf.append("\tchainscore\talnscore\tppositive\tnpositive\tnnegative\tcigar\tqseq\tsseq\tvolume\n");
        for (const auto& h : hits) {
            if (h.skip_reason != 0) {
                append_str(buf, h.qseqid);                              buf.push_back('\t');
                append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                buf.append("\t*\t0\t0\t");
                append_uint(buf, h.qlen);
                buf.append("\t0\t0\t0\t0\t0\t0\t0\t0\t0\t*\t*\t*\t0\n");
                write_block_if_full(out, buf);
                continue;
            }
            append_str(buf, h.qseqid);           buf.push_back('\t');
            append_str(buf, h.sseqid);           buf.push_back('\t');
            buf.push_back(h.sstrand);            buf.push_back('\t');
            append_uint(buf, h.qstart);          buf.push_back('\t');
            append_uint(buf, h.qend);            buf.push_back('\t');
            append_uint(buf, h.qlen);            buf.push_back('\t');
            append_uint(buf, h.sstart);          buf.push_back('\t');
            append_uint(buf, h.send);            buf.push_back('\t');
            append_uint(buf, h.slen);            buf.push_back('\t');
            append_uint(buf, h.coverscore);      buf.push_back('\t');
            append_uint(buf, h.chainscore);      buf.push_back('\t');
            append_int(buf, h.alnscore);         buf.push_back('\t');
            append_double_g6(buf, h.ppositive);  buf.push_back('\t');
            append_uint(buf, h.npositive);       buf.push_back('\t');
            append_uint(buf, h.nnegative);       buf.push_back('\t');
            append_str(buf, h.cigar);            buf.push_back('\t');
            append_str(buf, h.qseq);             buf.push_back('\t');
            append_str(buf, h.sseq);             buf.push_back('\t');
            append_uint(buf, h.volume);          buf.push_back('\n');
            write_block_if_full(out, buf);
        }
    } else if (mode == 3) {
        buf.append("# qseqid\tsseqid\tsstrand\tqend\tqlen\tsend\tslen\t");
        buf.append(s1name);
        buf.append("\tchainscore\talnscore\tvolume\n");
        for (const auto& h : hits) {
            if (h.skip_reason != 0) {
                append_str(buf, h.qseqid);                              buf.push_back('\t');
                append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                buf.append("\t*\t0\t");
                append_uint(buf, h.qlen);
                buf.append("\t0\t0\t0\t0\t0\t0\n");
                write_block_if_full(out, buf);
                continue;
            }
            append_str(buf, h.qseqid);        buf.push_back('\t');
            append_str(buf, h.sseqid);        buf.push_back('\t');
            buf.push_back(h.sstrand);         buf.push_back('\t');
            append_uint(buf, h.qend);         buf.push_back('\t');
            append_uint(buf, h.qlen);         buf.push_back('\t');
            append_uint(buf, h.send);         buf.push_back('\t');
            append_uint(buf, h.slen);         buf.push_back('\t');
            append_uint(buf, h.coverscore);   buf.push_back('\t');
            append_uint(buf, h.chainscore);   buf.push_back('\t');
            append_int(buf, h.alnscore);      buf.push_back('\t');
            append_uint(buf, h.volume);       buf.push_back('\n');
            write_block_if_full(out, buf);
        }
    } else {
        buf.append("# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\t");
        buf.append(s1name);
        buf.append("\tchainscore\tvolume\n");
        for (const auto& h : hits) {
            if (h.skip_reason != 0) {
                append_str(buf, h.qseqid);                              buf.push_back('\t');
                append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                buf.append("\t*\t0\t0\t");
                append_uint(buf, h.qlen);
                buf.append("\t0\t0\t0\t0\t0\t0\n");
                write_block_if_full(out, buf);
                continue;
            }
            append_str(buf, h.qseqid);        buf.push_back('\t');
            append_str(buf, h.sseqid);        buf.push_back('\t');
            buf.push_back(h.sstrand);         buf.push_back('\t');
            append_uint(buf, h.qstart);       buf.push_back('\t');
            append_uint(buf, h.qend);         buf.push_back('\t');
            append_uint(buf, h.qlen);         buf.push_back('\t');
            append_uint(buf, h.sstart);       buf.push_back('\t');
            append_uint(buf, h.send);         buf.push_back('\t');
            append_uint(buf, h.slen);         buf.push_back('\t');
            append_uint(buf, h.coverscore);   buf.push_back('\t');
            append_uint(buf, h.chainscore);   buf.push_back('\t');
            append_uint(buf, h.volume);       buf.push_back('\n');
            write_block_if_full(out, buf);
        }
    }

    write_block(out, buf);
}

// Inner helper: write per-query JSON objects.
// If is_fragment is true, objects are emitted without outer wrapper and
// with trailing comma after each (caller handles the last comma).
static void write_results_json_inner(std::ostream& out,
                                      const std::vector<OutputHit>& hits,
                                      uint8_t mode,
                                      bool stage3_traceback,
                                      bool is_fragment) {
    const char* s1name = kStage1ScoreName;
    std::string buf;
    buf.reserve(kBlockSize + 4096);

    // Group hits by qseqid (preserve order of first appearance)
    std::vector<std::string> query_order;
    std::map<std::string, std::vector<const OutputHit*>> by_query;
    std::map<std::string, const OutputHit*> by_query_skip; // first skip marker per query
    for (const auto& h : hits) {
        if (by_query.find(h.qseqid) == by_query.end() &&
            by_query_skip.find(h.qseqid) == by_query_skip.end()) {
            query_order.push_back(h.qseqid);
        }
        if (h.skip_reason != 0) {
            if (by_query_skip.find(h.qseqid) == by_query_skip.end())
                by_query_skip[h.qseqid] = &h;
        } else {
            by_query[h.qseqid].push_back(&h);
        }
    }

    for (size_t qi = 0; qi < query_order.size(); qi++) {
        const auto& qid = query_order[qi];
        auto skip_it = by_query_skip.find(qid);
        if (skip_it != by_query_skip.end() && by_query[qid].empty()) {
            const OutputHit* sk = skip_it->second;
            const bool failed = (sk->skip_reason == kFailHttpJob);
            buf.append("    {\n      \"qseqid\": ");
            json_escape_into(buf, qid);
            buf.append(",\n      \"qlen\": ");
            append_uint(buf, sk->qlen);
            if (failed) {
                buf.append(",\n      \"status\": \"failed\"");
                buf.append(",\n      \"reason\": ");
                json_escape_into(buf, sk->skip_detail);
            } else {
                buf.append(",\n      \"status\": \"skipped\"");
                buf.append(",\n      \"skip_reason\": ");
                json_escape_into(buf, skip_reason_str(sk->skip_reason));
                buf.append(",\n      \"skip_detail\": ");
                json_escape_into(buf, sk->skip_detail);
            }
            buf.append(",\n      \"hits\": []\n    }");
            if (is_fragment || qi + 1 < query_order.size()) buf.push_back(',');
            buf.push_back('\n');
            write_block_if_full(out, buf);
            continue;
        }
        const auto& qhits = by_query[qid];
        buf.append("    {\n      \"qseqid\": ");
        json_escape_into(buf, qid);
        buf.append(",\n      \"status\": \"ok\"");
        buf.append(",\n      \"hits\": [\n");
        for (size_t hi = 0; hi < qhits.size(); hi++) {
            const auto* h = qhits[hi];
            buf.append("        {\n");
            buf.append("          \"sseqid\": ");
            json_escape_into(buf, h->sseqid);
            buf.append(",\n");
            buf.append("          \"sstrand\": \"");
            buf.push_back(h->sstrand);
            buf.append("\",\n");
            if (mode == 2 || (mode == 3 && stage3_traceback)) {
                buf.append("          \"qstart\": ");
                append_uint(buf, h->qstart);
                buf.append(",\n          \"qend\": ");
                append_uint(buf, h->qend);
                buf.append(",\n");
            } else if (mode == 3) {
                buf.append("          \"qend\": ");
                append_uint(buf, h->qend);
                buf.append(",\n");
            }
            buf.append("          \"qlen\": ");
            append_uint(buf, h->qlen);
            buf.append(",\n");
            if (mode == 2 || (mode == 3 && stage3_traceback)) {
                buf.append("          \"sstart\": ");
                append_uint(buf, h->sstart);
                buf.append(",\n          \"send\": ");
                append_uint(buf, h->send);
                buf.append(",\n");
            } else if (mode == 3) {
                buf.append("          \"send\": ");
                append_uint(buf, h->send);
                buf.append(",\n");
            }
            buf.append("          \"slen\": ");
            append_uint(buf, h->slen);
            buf.append(",\n          \"");
            buf.append(s1name);
            buf.append("\": ");
            append_uint(buf, h->coverscore);
            buf.append(",\n");
            if (mode != 1) {
                buf.append("          \"chainscore\": ");
                append_uint(buf, h->chainscore);
                buf.append(",\n");
            }
            if (mode == 3) {
                buf.append("          \"alnscore\": ");
                append_int(buf, h->alnscore);
                buf.append(",\n");
                if (stage3_traceback) {
                    buf.append("          \"ppositive\": ");
                    append_double_g6(buf, h->ppositive);
                    buf.append(",\n          \"npositive\": ");
                    append_uint(buf, h->npositive);
                    buf.append(",\n          \"nnegative\": ");
                    append_uint(buf, h->nnegative);
                    buf.append(",\n          \"cigar\": ");
                    json_escape_into(buf, h->cigar);
                    buf.append(",\n          \"qseq\": ");
                    json_escape_into(buf, h->qseq);
                    buf.append(",\n          \"sseq\": ");
                    json_escape_into(buf, h->sseq);
                    buf.append(",\n");
                }
            }
            buf.append("          \"volume\": ");
            append_uint(buf, h->volume);
            buf.append("\n        }");
            if (hi + 1 < qhits.size()) buf.push_back(',');
            buf.push_back('\n');
            write_block_if_full(out, buf);
        }
        buf.append("      ]\n    }");
        if (is_fragment || qi + 1 < query_order.size()) buf.push_back(',');
        buf.push_back('\n');
        write_block_if_full(out, buf);
    }

    write_block(out, buf);
}

void write_results_json(std::ostream& out,
                        const std::vector<OutputHit>& hits,
                        uint8_t mode,
                        bool stage3_traceback) {
    out << "{\n  \"results\": [\n";
    write_results_json_inner(out, hits, mode, stage3_traceback, false);
    out << "  ]\n}\n";
}

void write_results_json_fragment(std::ostream& out,
                                  const std::vector<OutputHit>& hits,
                                  uint8_t mode,
                                  bool stage3_traceback) {
    write_results_json_inner(out, hits, mode, stage3_traceback, true);
}

void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt,
                   uint8_t mode,
                   bool stage3_traceback) {
    switch (fmt) {
        case OutputFormat::kTsv:
            write_results_tsv(out, hits, mode, stage3_traceback);
            break;
        case OutputFormat::kJson:
            write_results_json(out, hits, mode, stage3_traceback);
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
    if (fmt == OutputFormat::kSam || fmt == OutputFormat::kBam) {
        if (detect_format_from_extension(output_path) != CompressionFormat::kNone) {
            error_msg = "Error: SAM/BAM output does not support compression "
                        "suffix; remove .gz/.bz2/.xz/.zst or switch to "
                        "-output_format tsv/json";
            return false;
        }
    }
    return true;
}

// Mode 1 parallel TSV / JSON writers live in
// src/search/result_writer_mode1.cpp (ikafssn_search) because they
// resolve sseqid / slen via KsxReader; ikafssn_io does not depend on
// ikafssn_index.

bool parse_output_format(const std::string& str, OutputFormat& out,
                         std::string& error_msg) {
    if (str == "tsv") {
        out = OutputFormat::kTsv;
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
