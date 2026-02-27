#include "io/sam_writer.hpp"
#include "core/version.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <htslib/sam.h>
#include <htslib/hts.h>

namespace ikafssn {

// Parse a parasail-style extended CIGAR string (e.g., "3=2I1X4D14=")
// into htslib uint32_t array.
static std::vector<uint32_t> parse_cigar_to_htslib(const std::string& cigar_str) {
    std::vector<uint32_t> cigar;
    size_t i = 0;
    while (i < cigar_str.size()) {
        uint32_t len = 0;
        while (i < cigar_str.size() && cigar_str[i] >= '0' && cigar_str[i] <= '9') {
            len = len * 10 + (cigar_str[i] - '0');
            i++;
        }
        if (i >= cigar_str.size()) break;
        char op_char = cigar_str[i++];
        int op;
        switch (op_char) {
            case '=': op = BAM_CEQUAL; break;
            case 'X': op = BAM_CDIFF;  break;
            case 'I': op = BAM_CINS;   break;
            case 'D': op = BAM_CDEL;   break;
            case 'M': op = BAM_CMATCH; break;
            default:  op = BAM_CMATCH; break;
        }
        cigar.push_back(bam_cigar_gen(len, op));
    }
    return cigar;
}

// Extract the raw query sequence (without gaps) from aligned q_seq,
// or from the original query_id. For SAM, we need the ungapped sequence.
static std::string ungap_sequence(const std::string& aligned_seq) {
    std::string seq;
    seq.reserve(aligned_seq.size());
    for (char c : aligned_seq) {
        if (c != '-') seq += c;
    }
    return seq;
}

static void write_sam_bam_impl(const std::string& output_path,
                                const std::vector<OutputHit>& hits,
                                uint8_t stage1_score_type,
                                bool is_bam) {
    const char* mode = is_bam ? "wb" : "w";
    const char* path = output_path.c_str();
    if (!is_bam && (output_path.empty() || output_path == "-")) {
        path = "-";
    }

    samFile* fp = sam_open(path, mode);
    if (!fp) return;

    sam_hdr_t* hdr = sam_hdr_init();

    // @HD
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "unsorted", NULL);

    // @SQ: collect unique (sseqid, slen) pairs, ordered by first appearance
    std::vector<std::string> sq_order;
    std::map<std::string, uint32_t> sq_lengths;
    for (const auto& h : hits) {
        if (sq_lengths.find(h.sseqid) == sq_lengths.end()) {
            sq_order.push_back(h.sseqid);
            sq_lengths[h.sseqid] = h.slen;
        }
    }
    for (const auto& acc : sq_order) {
        std::string len_str = std::to_string(sq_lengths[acc]);
        sam_hdr_add_line(hdr, "SQ", "SN", acc.c_str(), "LN", len_str.c_str(), NULL);
    }

    // @PG
    sam_hdr_add_line(hdr, "PG", "ID", "ikafssnsearch", "VN", IKAFSSN_VERSION, NULL);

    sam_hdr_write(fp, hdr);

    // Write records
    bam1_t* b = bam_init1();

    for (const auto& h : hits) {
        uint16_t flag = 0;
        if (h.sstrand == '-') flag |= BAM_FREVERSE;

        int32_t tid = sam_hdr_name2tid(hdr, h.sseqid.c_str());
        hts_pos_t pos = static_cast<hts_pos_t>(h.sstart); // 0-based in htslib

        // Parse CIGAR
        auto cigar = parse_cigar_to_htslib(h.cigar);

        // Get ungapped query sequence
        std::string seq = ungap_sequence(h.qseq);

        // Set the record
        bam_set1(b,
                 h.qseqid.size(), h.qseqid.c_str(),
                 flag, tid, pos, 255,
                 cigar.size(), cigar.data(),
                 -1, -1, 0,
                 seq.size(), seq.c_str(), nullptr,
                 0);

        // Add aux tags
        // AS:i:alnscore
        {
            int32_t as = h.alnscore;
            bam_aux_append(b, "AS", 'i', sizeof(int32_t),
                           reinterpret_cast<const uint8_t*>(&as));
        }
        // NM:i:mismatch
        {
            int32_t nm = static_cast<int32_t>(h.mismatch);
            bam_aux_append(b, "NM", 'i', sizeof(int32_t),
                           reinterpret_cast<const uint8_t*>(&nm));
        }
        // cs:i:chainscore
        {
            int32_t cs = static_cast<int32_t>(h.chainscore);
            bam_aux_append(b, "cs", 'i', sizeof(int32_t),
                           reinterpret_cast<const uint8_t*>(&cs));
        }
        // cv:i:coverscore
        {
            int32_t cv = static_cast<int32_t>(h.coverscore);
            bam_aux_append(b, "cv", 'i', sizeof(int32_t),
                           reinterpret_cast<const uint8_t*>(&cv));
        }
        // ms:i:matchscore
        {
            int32_t ms = static_cast<int32_t>(h.matchscore);
            bam_aux_append(b, "ms", 'i', sizeof(int32_t),
                           reinterpret_cast<const uint8_t*>(&ms));
        }

        sam_write1(fp, hdr, b);
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);
}

void write_results_sam(const std::string& output_path,
                       const std::vector<OutputHit>& hits,
                       uint8_t stage1_score_type) {
    write_sam_bam_impl(output_path, hits, stage1_score_type, false);
}

void write_results_bam(const std::string& output_path,
                       const std::vector<OutputHit>& hits,
                       uint8_t stage1_score_type) {
    write_sam_bam_impl(output_path, hits, stage1_score_type, true);
}

bool merge_sam_files(const std::vector<std::string>& batch_paths,
                     const std::string& output_path, bool as_bam) {
    if (batch_paths.empty()) return true;

    // Phase 1: Build merged header by reading all input files' headers.
    // Collect @SQ lines (union by SN, ordered by first appearance).
    std::vector<std::string> sq_order;
    std::map<std::string, uint32_t> sq_lengths;
    sam_hdr_t* first_hdr = nullptr;

    for (const auto& path : batch_paths) {
        samFile* in = sam_open(path.c_str(), "r");
        if (!in) continue;

        sam_hdr_t* hdr = sam_hdr_read(in);
        if (!hdr) { sam_close(in); continue; }

        if (!first_hdr) {
            first_hdr = sam_hdr_dup(hdr);
        }

        int nsq = sam_hdr_nref(hdr);
        for (int i = 0; i < nsq; i++) {
            const char* sn = sam_hdr_tid2name(hdr, i);
            uint64_t ln = sam_hdr_tid2len(hdr, i);
            std::string name(sn);
            if (sq_lengths.find(name) == sq_lengths.end()) {
                sq_order.push_back(name);
                sq_lengths[name] = static_cast<uint32_t>(ln);
            }
        }

        sam_hdr_destroy(hdr);
        sam_close(in);
    }

    if (!first_hdr) return false;

    // Build merged header
    sam_hdr_t* merged = sam_hdr_init();

    // Copy @HD from first_hdr
    if (sam_hdr_find_line_id(first_hdr, "HD", NULL, NULL, NULL) == 0) {
        sam_hdr_add_line(merged, "HD", "VN", "1.6", "SO", "unsorted", NULL);
    }

    // Add @SQ lines
    for (const auto& sn : sq_order) {
        std::string ln_str = std::to_string(sq_lengths[sn]);
        sam_hdr_add_line(merged, "SQ", "SN", sn.c_str(), "LN", ln_str.c_str(), NULL);
    }

    // Copy @PG from first_hdr
    {
        int npg = sam_hdr_count_lines(first_hdr, "PG");
        for (int i = 0; i < npg; i++) {
            // Extract PG ID and VN from first header
            kstring_t ks = KS_INITIALIZE;
            if (sam_hdr_find_line_pos(first_hdr, "PG", i, &ks) == 0) {
                sam_hdr_add_lines(merged, ks.s, ks.l);
            }
            ks_free(&ks);
        }
    }

    sam_hdr_destroy(first_hdr);

    // Phase 2: Open output and write merged header
    const char* mode = as_bam ? "wb" : "w";
    const char* out_path = output_path.c_str();
    if (!as_bam && (output_path.empty() || output_path == "-")) {
        out_path = "-";
    }

    samFile* out = sam_open(out_path, mode);
    if (!out) { sam_hdr_destroy(merged); return false; }

    sam_hdr_write(out, merged);

    // Phase 3: Read each batch file and remap tids to merged header
    bam1_t* b = bam_init1();
    for (const auto& path : batch_paths) {
        samFile* in = sam_open(path.c_str(), "r");
        if (!in) continue;

        sam_hdr_t* batch_hdr = sam_hdr_read(in);
        if (!batch_hdr) { sam_close(in); continue; }

        // Build tid remapping: batch_tid -> merged_tid
        int batch_nsq = sam_hdr_nref(batch_hdr);
        std::vector<int32_t> tid_map(batch_nsq, -1);
        for (int i = 0; i < batch_nsq; i++) {
            const char* sn = sam_hdr_tid2name(batch_hdr, i);
            int32_t merged_tid = sam_hdr_name2tid(merged, sn);
            tid_map[i] = merged_tid;
        }

        while (sam_read1(in, batch_hdr, b) >= 0) {
            // Remap tid
            if (b->core.tid >= 0 && b->core.tid < batch_nsq) {
                b->core.tid = tid_map[b->core.tid];
            }
            // Remap mate tid (should be -1 for our data, but handle anyway)
            if (b->core.mtid >= 0 && b->core.mtid < batch_nsq) {
                b->core.mtid = tid_map[b->core.mtid];
            }
            sam_write1(out, merged, b);
        }

        sam_hdr_destroy(batch_hdr);
        sam_close(in);
    }

    bam_destroy1(b);
    sam_hdr_destroy(merged);
    sam_close(out);
    return true;
}

bool write_all_results(const std::string& output_path,
                       const std::vector<OutputHit>& hits,
                       OutputFormat fmt,
                       uint8_t mode,
                       uint8_t stage1_score_type,
                       bool stage3_traceback) {
    if (fmt == OutputFormat::kSam) {
        write_results_sam(output_path.empty() ? "-" : output_path,
                          hits, stage1_score_type);
    } else if (fmt == OutputFormat::kBam) {
        write_results_bam(output_path, hits, stage1_score_type);
    } else if (output_path.empty()) {
        write_results(std::cout, hits, fmt, mode, stage1_score_type,
                      stage3_traceback);
    } else {
        std::ofstream out(output_path);
        if (!out.is_open()) {
            std::fprintf(stderr, "Error: cannot open output file %s\n",
                         output_path.c_str());
            return false;
        }
        write_results(out, hits, fmt, mode, stage1_score_type,
                      stage3_traceback);
    }
    return true;
}

} // namespace ikafssn
