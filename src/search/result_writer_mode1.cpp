// Mode 1 parallel TSV / JSON writers.
//
// These live in ikafssn_search rather than ikafssn_io because they
// resolve sseqid / slen via KsxReader and consume OrchestratorHit;
// ikafssn_io does not depend on ikafssn_index or ikafssn_search.  The
// declarations live in io/result_writer.hpp so callers in main.cpp see
// them alongside the OutputHit-based writers.

#include "io/result_writer.hpp"

#include "io/fasta_reader.hpp"
#include "io/text_format.hpp"
#include "index/ksx_reader.hpp"
#include "protocol/messages.hpp"
#include "search/parallel_search.hpp"
#include "search/parent_hit.hpp"

#include <algorithm>
#include <string>
#include <string_view>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace ikafssn {

namespace {

// Stage 1 score is always coverscore.
constexpr const char* kMode1Stage1ScoreName = "coverscore";

void compute_chunk_ranges(size_t n, size_t nchunks,
                          std::vector<std::pair<size_t,size_t>>& out) {
    out.clear();
    if (nchunks == 0) nchunks = 1;
    out.reserve(nchunks);
    size_t per = n / nchunks;
    size_t rem = n % nchunks;
    size_t cur = 0;
    for (size_t i = 0; i < nchunks; ++i) {
        size_t take = per + (i < rem ? 1 : 0);
        out.emplace_back(cur, cur + take);
        cur += take;
    }
}

} // namespace

void write_results_tsv_mode1_parallel(std::ostream& out,
                                       const Mode1ParallelInputs& in) {
    const char* s1name = kMode1Stage1ScoreName;

    out << "# qseqid\tsseqid\tsstrand\tqlen\tslen\t" << s1name << "\tvolume\n";

    const auto& hits     = *in.hits;
    const auto& queries  = *in.queries;
    const size_t n_hits  = hits.size();

    const int nthread = std::max(1, in.nthread);
    const size_t nchunks = static_cast<size_t>(nthread) * 4;

    std::vector<std::pair<size_t,size_t>> ranges;
    compute_chunk_ranges(n_hits, nchunks, ranges);
    std::vector<std::string> chunk_bufs(ranges.size());

    if (n_hits > 0) {
        tbb::task_arena arena(nthread);
        arena.execute([&] {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, ranges.size()),
                [&](const tbb::blocked_range<size_t>& rng) {
                    for (size_t ci = rng.begin(); ci != rng.end(); ++ci) {
                        size_t lo = ranges[ci].first;
                        size_t hi = ranges[ci].second;
                        if (lo == hi) continue;
                        std::string& buf = chunk_bufs[ci];
                        buf.reserve((hi - lo) * 96);
                        for (size_t i = lo; i < hi; ++i) {
                            const auto& oh = hits[i];
                            const auto& cr = oh.cr;
                            const auto* ksx = in.ksx_per_volume[oh.volume_idx];
                            const std::string& qid = queries[oh.query_idx].id;
                            // Report parent accession + parent length so
                            // mode-1 output is fragment-agnostic.
                            const ParentName pn =
                                resolve_parent_name(*ksx, cr.seq_id);
                            std::string_view sseqid = pn.sseqid;
                            uint32_t qlen = static_cast<uint32_t>(
                                queries[oh.query_idx].sequence.size());
                            uint32_t slen = pn.slen;

                            append_str(buf, qid);          buf.push_back('\t');
                            append_str(buf, sseqid);       buf.push_back('\t');
                            buf.push_back(cr.is_reverse ? '-' : '+');
                            buf.push_back('\t');
                            append_uint(buf, qlen);        buf.push_back('\t');
                            append_uint(buf, slen);        buf.push_back('\t');
                            append_uint(buf, cr.stage1_score);
                            buf.push_back('\t');
                            append_uint(buf, oh.volume_index);
                            buf.push_back('\n');
                        }
                    }
                });
        });
    }

    for (auto& b : chunk_bufs) {
        if (!b.empty())
            out.write(b.data(), static_cast<std::streamsize>(b.size()));
    }

    if (in.skip_reason && in.skip_detail) {
        std::string sbuf;
        for (size_t qi = 0; qi < queries.size(); ++qi) {
            uint8_t reason = (*in.skip_reason)[qi];
            if (reason == 0) continue;
            uint32_t qlen = static_cast<uint32_t>(queries[qi].sequence.size());
            append_str(sbuf, queries[qi].id);   sbuf.push_back('\t');
            append_skipped_sseqid(sbuf, reason, (*in.skip_detail)[qi]);
            sbuf.push_back('\t');
            sbuf.push_back('*');                sbuf.push_back('\t');
            append_uint(sbuf, qlen);            sbuf.push_back('\t');
            sbuf.push_back('0');                sbuf.push_back('\t');
            sbuf.push_back('0');                sbuf.push_back('\t');
            sbuf.push_back('0');                sbuf.push_back('\n');
        }
        if (!sbuf.empty())
            out.write(sbuf.data(), static_cast<std::streamsize>(sbuf.size()));
    }
}

void write_results_json_mode1_parallel(std::ostream& out,
                                        const Mode1ParallelInputs& in) {
    const char* s1name = kMode1Stage1ScoreName;

    const auto& hits    = *in.hits;
    const auto& queries = *in.queries;
    const size_t n_q    = queries.size();
    const int nthread   = std::max(1, in.nthread);

    // hit_begin[q] = first hit index whose query_idx >= q.  Caller has
    // sorted `hits` by query_idx so each query's run is contiguous.
    std::vector<size_t> hit_begin(n_q + 1, hits.size());
    {
        size_t i = 0;
        for (size_t q = 0; q < n_q; ++q) {
            while (i < hits.size() && hits[i].query_idx < q) ++i;
            hit_begin[q] = i;
        }
        hit_begin[n_q] = hits.size();
    }

    const size_t nchunks_target = static_cast<size_t>(nthread) * 4;
    const size_t nchunks = std::min(std::max<size_t>(1, n_q),
                                     std::max<size_t>(1, nchunks_target));
    std::vector<std::pair<size_t,size_t>> qranges;
    compute_chunk_ranges(n_q, nchunks, qranges);
    std::vector<std::string> chunk_bufs(qranges.size());

    auto format_one_query = [&](std::string& buf, size_t qi,
                                 bool is_last_overall) {
        const std::string& qid = queries[qi].id;
        uint32_t qlen = static_cast<uint32_t>(queries[qi].sequence.size());
        uint8_t reason = in.skip_reason ? (*in.skip_reason)[qi] : 0;
        size_t lo = hit_begin[qi];
        size_t hi = hit_begin[qi + 1];

        if (reason != 0 && lo == hi) {
            const bool failed = (reason == kFailHttpJob);
            buf.append("    {\n      \"qseqid\": ");
            json_escape_into(buf, qid);
            buf.append(",\n      \"qlen\": ");
            append_uint(buf, qlen);
            if (failed) {
                buf.append(",\n      \"status\": \"failed\"");
                buf.append(",\n      \"reason\": ");
                json_escape_into(buf, in.skip_detail ? (*in.skip_detail)[qi] : "");
            } else {
                buf.append(",\n      \"status\": \"skipped\"");
                buf.append(",\n      \"skip_reason\": ");
                json_escape_into(buf, std::string_view(skip_reason_str(reason)));
                buf.append(",\n      \"skip_detail\": ");
                json_escape_into(buf, in.skip_detail ? (*in.skip_detail)[qi] : "");
            }
            buf.append(",\n      \"hits\": []\n    }");
            if (!is_last_overall) buf.push_back(',');
            buf.push_back('\n');
            return;
        }

        buf.append("    {\n      \"qseqid\": ");
        json_escape_into(buf, qid);
        buf.append(",\n      \"status\": \"ok\"");
        buf.append(",\n      \"hits\": [\n");
        for (size_t i = lo; i < hi; ++i) {
            const auto& oh = hits[i];
            const auto& cr = oh.cr;
            const auto* ksx = in.ksx_per_volume[oh.volume_idx];
            // Report parent accession + parent length so mode-1 JSON
            // output mirrors the TSV (fragment-agnostic).
            const ParentName pn = resolve_parent_name(*ksx, cr.seq_id);
            std::string_view sseqid = pn.sseqid;
            uint32_t slen = pn.slen;

            buf.append("        {\n          \"sseqid\": ");
            json_escape_into(buf, sseqid);
            buf.append(",\n          \"sstrand\": \"");
            buf.push_back(cr.is_reverse ? '-' : '+');
            buf.append("\",\n          \"qlen\": ");
            append_uint(buf, qlen);
            buf.append(",\n          \"slen\": ");
            append_uint(buf, slen);
            buf.append(",\n          \"");
            buf.append(s1name);
            buf.append("\": ");
            append_uint(buf, cr.stage1_score);
            buf.append(",\n          \"volume\": ");
            append_uint(buf, oh.volume_index);
            buf.append("\n        }");
            if (i + 1 < hi) buf.push_back(',');
            buf.push_back('\n');
        }
        buf.append("      ]\n    }");
        if (!is_last_overall) buf.push_back(',');
        buf.push_back('\n');
    };

    const size_t last_q = (n_q > 0) ? (n_q - 1) : 0;

    if (n_q > 0) {
        tbb::task_arena arena(nthread);
        arena.execute([&] {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, qranges.size()),
                [&](const tbb::blocked_range<size_t>& rng) {
                    for (size_t ci = rng.begin(); ci != rng.end(); ++ci) {
                        size_t qlo = qranges[ci].first;
                        size_t qhi = qranges[ci].second;
                        if (qlo == qhi) continue;
                        std::string& buf = chunk_bufs[ci];
                        buf.reserve((qhi - qlo) * 256);
                        for (size_t qi = qlo; qi < qhi; ++qi) {
                            format_one_query(buf, qi, qi == last_q);
                        }
                    }
                });
        });
    }

    out << "{\n  \"results\": [\n";
    for (auto& b : chunk_bufs) {
        if (!b.empty())
            out.write(b.data(), static_cast<std::streamsize>(b.size()));
    }
    out << "  ]\n}\n";
}

} // namespace ikafssn
