#include "io/result_writer.hpp"
#include "io/compressed_stream.hpp"
#include "io/output_coords.hpp"
#include "io/text_format.hpp"
#include "protocol/messages.hpp"
#include "util/query_order.hpp"

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace ikafssn {

// Stage 1 score is always coverscore (the number of distinct query k-mers
// matching the sequence).
static constexpr const char* kStage1ScoreName = "coverscore";

namespace {

// Formatted rows are handed to the stream once they reach this size.
constexpr std::size_t kFlushSize = 1u << 20;

// Work formatted per wave, counted in rows (TSV) or hits (JSON).  The
// parallel path holds one wave of output in memory, so bounding it bounds
// the extra footprint.
constexpr std::size_t kWaveWeight = 64 * 1024;

// Below this much total work the arena setup costs more than the
// formatting it parallelises, so the serial path stays faster.
constexpr std::size_t kParallelMinWeight = 256 * 1024;

inline void flush_buf(std::ostream& out, std::string& buf) {
    if (!buf.empty()) {
        out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
        buf.clear();
    }
}

// Format items [begin, end) and write them to `out` in item order.
// `format_item(buf, i)` appends item i, `weight_of(i)` is its share of the
// wave budget and `total_weight` their sum.  A wave is formatted into
// per-chunk buffers concurrently and written in chunk order, so the bytes
// do not depend on the thread count.
template <typename FormatItem, typename WeightOf>
void write_items(std::ostream& out, std::size_t begin, std::size_t end,
                 std::size_t total_weight, int nthread,
                 FormatItem format_item, WeightOf weight_of) {
    if (begin >= end) return;

    const int threads = std::max(1, nthread);
    if (threads == 1 || total_weight < kParallelMinWeight) {
        std::string buf;
        buf.reserve(kFlushSize + 4096);
        for (std::size_t i = begin; i < end; ++i) {
            format_item(buf, i);
            if (buf.size() >= kFlushSize) flush_buf(out, buf);
        }
        flush_buf(out, buf);
        return;
    }

    const std::size_t nchunks = static_cast<std::size_t>(threads) * 4;
    std::vector<std::string> bufs(nchunks);
    tbb::task_arena arena(threads);

    for (std::size_t wlo = begin; wlo < end;) {
        // Always take at least one item, so a single oversized item still
        // makes progress.
        std::size_t whi = wlo, weight = 0;
        do {
            weight += weight_of(whi);
            ++whi;
        } while (whi < end && weight < kWaveWeight);

        const std::size_t n = whi - wlo;
        const std::size_t nc = std::min(n, nchunks);
        arena.execute([&] {
            tbb::parallel_for(
                tbb::blocked_range<std::size_t>(0, nc),
                [&](const tbb::blocked_range<std::size_t>& r) {
                    for (std::size_t ci = r.begin(); ci != r.end(); ++ci) {
                        // flush_buf() empties the buffer but keeps its
                        // capacity, so waves reuse the allocations.
                        std::string& buf = bufs[ci];
                        const std::size_t lo = wlo + n * ci / nc;
                        const std::size_t hi = wlo + n * (ci + 1) / nc;
                        for (std::size_t i = lo; i < hi; ++i)
                            format_item(buf, i);
                    }
                });
        });
        for (std::size_t ci = 0; ci < nc; ++ci) flush_buf(out, bufs[ci]);
        wlo = whi;
    }
}

// Every TSV row costs the same share of the wave budget.
inline std::size_t one_row(std::size_t) { return 1; }

inline OutputCoords blast_coords(const OutputHit& h) {
    return to_output_coords(h.qstart, h.qend, h.sstart, h.send,
                            h.qlen, h.sstrand == '-');
}

} // namespace

void write_results_tsv(std::ostream& out,
                       const std::vector<OutputHit>& hits,
                       uint8_t mode,
                       bool stage3_traceback,
                       int nthread) {
    const char* s1name = kStage1ScoreName;

    std::string header;
    if (mode == 1) {
        header.append("# qseqid\tsseqid\tsstrand\tqlen\tslen\t");
        header.append(s1name);
        header.append("\tvolume\n");
    } else if (mode == 3 && stage3_traceback) {
        header.append("# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\t");
        header.append(s1name);
        header.append("\tchainscore\talnscore\tppositive\tnpositive\tnnegative\tcigar\tqseq\tsseq\tvolume\n");
    } else if (mode == 3) {
        header.append("# qseqid\tsseqid\tsstrand\tqend\tqlen\tsend\tslen\t");
        header.append(s1name);
        header.append("\tchainscore\talnscore\tvolume\n");
    } else {
        header.append("# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\t");
        header.append(s1name);
        header.append("\tchainscore\tvolume\n");
    }
    out.write(header.data(), static_cast<std::streamsize>(header.size()));

    if (mode == 1) {
        write_items(out, 0, hits.size(), hits.size(), nthread,
            [&](std::string& buf, std::size_t i) {
                const auto& h = hits[i];
                if (h.skip_reason != 0) {
                    append_str(buf, h.qseqid);                              buf.push_back('\t');
                    append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                    buf.append("\t*\t");
                    append_uint(buf, h.qlen);
                    buf.append("\t0\t0\t0\n");
                    return;
                }
                append_str(buf, h.qseqid);        buf.push_back('\t');
                append_str(buf, h.sseqid);        buf.push_back('\t');
                buf.push_back(h.sstrand);         buf.push_back('\t');
                append_uint(buf, h.qlen);         buf.push_back('\t');
                append_uint(buf, h.slen);         buf.push_back('\t');
                append_uint(buf, h.coverscore);   buf.push_back('\t');
                append_uint(buf, h.volume);       buf.push_back('\n');
            }, one_row);
    } else if (mode == 3 && stage3_traceback) {
        write_items(out, 0, hits.size(), hits.size(), nthread,
            [&](std::string& buf, std::size_t i) {
                const auto& h = hits[i];
                if (h.skip_reason != 0) {
                    append_str(buf, h.qseqid);                              buf.push_back('\t');
                    append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                    buf.append("\t*\t0\t0\t");
                    append_uint(buf, h.qlen);
                    buf.append("\t0\t0\t0\t0\t0\t0\t0\t0\t0\t*\t*\t*\t0\n");
                    return;
                }
                const OutputCoords c = blast_coords(h);
                append_str(buf, h.qseqid);           buf.push_back('\t');
                append_str(buf, h.sseqid);           buf.push_back('\t');
                buf.push_back(h.sstrand);            buf.push_back('\t');
                append_uint(buf, c.qstart);          buf.push_back('\t');
                append_uint(buf, c.qend);            buf.push_back('\t');
                append_uint(buf, h.qlen);            buf.push_back('\t');
                append_uint(buf, c.sstart);          buf.push_back('\t');
                append_uint(buf, c.send);            buf.push_back('\t');
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
            }, one_row);
    } else if (mode == 3) {
        write_items(out, 0, hits.size(), hits.size(), nthread,
            [&](std::string& buf, std::size_t i) {
                const auto& h = hits[i];
                if (h.skip_reason != 0) {
                    append_str(buf, h.qseqid);                              buf.push_back('\t');
                    append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                    buf.append("\t*\t0\t");
                    append_uint(buf, h.qlen);
                    buf.append("\t0\t0\t0\t0\t0\t0\n");
                    return;
                }
                const OutputCoords c = blast_coords(h);
                append_str(buf, h.qseqid);        buf.push_back('\t');
                append_str(buf, h.sseqid);        buf.push_back('\t');
                buf.push_back(h.sstrand);         buf.push_back('\t');
                append_uint(buf, c.qend);         buf.push_back('\t');
                append_uint(buf, h.qlen);         buf.push_back('\t');
                append_uint(buf, c.send);         buf.push_back('\t');
                append_uint(buf, h.slen);         buf.push_back('\t');
                append_uint(buf, h.coverscore);   buf.push_back('\t');
                append_uint(buf, h.chainscore);   buf.push_back('\t');
                append_int(buf, h.alnscore);      buf.push_back('\t');
                append_uint(buf, h.volume);       buf.push_back('\n');
            }, one_row);
    } else {
        write_items(out, 0, hits.size(), hits.size(), nthread,
            [&](std::string& buf, std::size_t i) {
                const auto& h = hits[i];
                if (h.skip_reason != 0) {
                    append_str(buf, h.qseqid);                              buf.push_back('\t');
                    append_skipped_sseqid(buf, h.skip_reason, h.skip_detail);
                    buf.append("\t*\t0\t0\t");
                    append_uint(buf, h.qlen);
                    buf.append("\t0\t0\t0\t0\t0\t0\n");
                    return;
                }
                const OutputCoords c = blast_coords(h);
                append_str(buf, h.qseqid);        buf.push_back('\t');
                append_str(buf, h.sseqid);        buf.push_back('\t');
                buf.push_back(h.sstrand);         buf.push_back('\t');
                append_uint(buf, c.qstart);       buf.push_back('\t');
                append_uint(buf, c.qend);         buf.push_back('\t');
                append_uint(buf, h.qlen);         buf.push_back('\t');
                append_uint(buf, c.sstart);       buf.push_back('\t');
                append_uint(buf, c.send);         buf.push_back('\t');
                append_uint(buf, h.slen);         buf.push_back('\t');
                append_uint(buf, h.coverscore);   buf.push_back('\t');
                append_uint(buf, h.chainscore);   buf.push_back('\t');
                append_uint(buf, h.volume);       buf.push_back('\n');
            }, one_row);
    }
}

// Inner helper: write per-query JSON objects.
// If is_fragment is true, objects are emitted without outer wrapper and
// with trailing comma after each (caller handles the last comma).
static void write_results_json_inner(std::ostream& out,
                                      const std::vector<OutputHit>& hits,
                                      uint8_t mode,
                                      bool stage3_traceback,
                                      bool is_fragment,
                                      int nthread) {
    const char* s1name = kStage1ScoreName;

    // Group hits by qseqid, keeping the order of first appearance.  Each
    // query owns the slice [query_begin[qi], query_begin[qi] + hit_count[qi])
    // of `query_hits`; a query's hits need not be contiguous in `hits`.
    QueryOrder qorder;
    std::vector<uint32_t> group_of(hits.size());
    std::vector<uint32_t> hit_count;                // non-skip hits per query
    std::vector<const OutputHit*> skip_of;          // first skip marker per query
    for (std::size_t i = 0; i < hits.size(); ++i) {
        const uint32_t g = qorder.id_of(hits[i].qseqid);
        group_of[i] = g;
        if (g == hit_count.size()) {
            hit_count.push_back(0);
            skip_of.push_back(nullptr);
        }
        if (hits[i].skip_reason != 0) {
            if (skip_of[g] == nullptr) skip_of[g] = &hits[i];
        } else {
            ++hit_count[g];
        }
    }

    const std::size_t nqueries = hit_count.size();
    std::vector<uint32_t> query_begin(nqueries + 1, 0);
    for (std::size_t g = 0; g < nqueries; ++g)
        query_begin[g + 1] = query_begin[g] + hit_count[g];
    std::vector<const OutputHit*> query_hits(query_begin[nqueries]);
    {
        std::vector<uint32_t> fill(query_begin.begin(), query_begin.end() - 1);
        for (std::size_t i = 0; i < hits.size(); ++i) {
            if (hits[i].skip_reason != 0) continue;
            query_hits[fill[group_of[i]]++] = &hits[i];
        }
    }

    write_items(out, 0, nqueries, hits.size(), nthread,
        [&](std::string& buf, std::size_t qi) {
        const std::string_view qid = qorder.qseqids()[qi];
        if (skip_of[qi] != nullptr && hit_count[qi] == 0) {
            const OutputHit* sk = skip_of[qi];
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
            if (is_fragment || qi + 1 < nqueries) buf.push_back(',');
            buf.push_back('\n');
            return;
        }
        const std::size_t nhits = hit_count[qi];
        const OutputHit* const* qhits = query_hits.data() + query_begin[qi];
        buf.append("    {\n      \"qseqid\": ");
        json_escape_into(buf, qid);
        buf.append(",\n      \"status\": \"ok\"");
        buf.append(",\n      \"hits\": [\n");
        for (std::size_t hi = 0; hi < nhits; hi++) {
            const auto* h = qhits[hi];
            const OutputCoords c = blast_coords(*h);
            buf.append("        {\n");
            buf.append("          \"sseqid\": ");
            json_escape_into(buf, h->sseqid);
            buf.append(",\n");
            buf.append("          \"sstrand\": \"");
            buf.push_back(h->sstrand);
            buf.append("\",\n");
            if (mode == 2 || (mode == 3 && stage3_traceback)) {
                buf.append("          \"qstart\": ");
                append_uint(buf, c.qstart);
                buf.append(",\n          \"qend\": ");
                append_uint(buf, c.qend);
                buf.append(",\n");
            } else if (mode == 3) {
                buf.append("          \"qend\": ");
                append_uint(buf, c.qend);
                buf.append(",\n");
            }
            buf.append("          \"qlen\": ");
            append_uint(buf, h->qlen);
            buf.append(",\n");
            if (mode == 2 || (mode == 3 && stage3_traceback)) {
                buf.append("          \"sstart\": ");
                append_uint(buf, c.sstart);
                buf.append(",\n          \"send\": ");
                append_uint(buf, c.send);
                buf.append(",\n");
            } else if (mode == 3) {
                buf.append("          \"send\": ");
                append_uint(buf, c.send);
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
            if (hi + 1 < nhits) buf.push_back(',');
            buf.push_back('\n');
        }
        buf.append("      ]\n    }");
        if (is_fragment || qi + 1 < nqueries) buf.push_back(',');
        buf.push_back('\n');
        },
        [&](std::size_t qi) { return hit_count[qi] + std::size_t{1}; });
}

void write_results_json(std::ostream& out,
                        const std::vector<OutputHit>& hits,
                        uint8_t mode,
                        bool stage3_traceback,
                        int nthread) {
    out << "{\n  \"results\": [\n";
    write_results_json_inner(out, hits, mode, stage3_traceback, false, nthread);
    out << "  ]\n}\n";
}

void write_results_json_fragment(std::ostream& out,
                                  const std::vector<OutputHit>& hits,
                                  uint8_t mode,
                                  bool stage3_traceback,
                                  int nthread) {
    write_results_json_inner(out, hits, mode, stage3_traceback, true, nthread);
}

void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt,
                   uint8_t mode,
                   bool stage3_traceback,
                   int nthread) {
    switch (fmt) {
        case OutputFormat::kTsv:
            write_results_tsv(out, hits, mode, stage3_traceback, nthread);
            break;
        case OutputFormat::kJson:
            write_results_json(out, hits, mode, stage3_traceback, nthread);
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
