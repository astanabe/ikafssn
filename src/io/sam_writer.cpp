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

    // @SQ: collect unique (accession, s_length) pairs, ordered by first appearance
    std::vector<std::string> sq_order;
    std::map<std::string, uint32_t> sq_lengths;
    for (const auto& h : hits) {
        if (sq_lengths.find(h.accession) == sq_lengths.end()) {
            sq_order.push_back(h.accession);
            sq_lengths[h.accession] = h.s_length;
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
        if (h.strand == '-') flag |= BAM_FREVERSE;

        int32_t tid = sam_hdr_name2tid(hdr, h.accession.c_str());
        hts_pos_t pos = static_cast<hts_pos_t>(h.s_start); // 0-based in htslib

        // Parse CIGAR
        auto cigar = parse_cigar_to_htslib(h.cigar);

        // Get ungapped query sequence
        std::string seq = ungap_sequence(h.q_seq);

        // Set the record
        bam_set1(b,
                 h.query_id.size(), h.query_id.c_str(),
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
        // NM:i:nmismatch
        {
            int32_t nm = static_cast<int32_t>(h.nmismatch);
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
