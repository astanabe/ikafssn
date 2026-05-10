#include "test_util.hpp"
#include "io/result_writer.hpp"
#include "io/result_reader.hpp"
#include "io/sam_writer.hpp"
#include "protocol/messages.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <htslib/sam.h>

using namespace ikafssn;

// Build an OutputHit representing a skipped query.
static OutputHit make_skip_hit(const std::string& qid, uint8_t reason,
                                const std::string& detail, uint32_t qlen) {
    OutputHit h;
    h.qseqid = qid;
    h.qlen = qlen;
    h.skip_reason = reason;
    h.skip_detail = detail;
    h.sstrand = '*';
    return h;
}

// Build a normal hit.
static OutputHit make_normal_hit(const std::string& qid,
                                  const std::string& sseqid,
                                  uint32_t coverscore) {
    OutputHit h;
    h.qseqid = qid;
    h.sseqid = sseqid;
    h.sstrand = '+';
    h.qstart = 1; h.qend = 50; h.qlen = 100;
    h.sstart = 1; h.send = 50; h.slen = 200;
    h.coverscore = coverscore;
    h.chainscore = coverscore;
    h.volume = 0;
    return h;
}

// TSV emits a sentinel "*SKIPPED:reason" row for skipped queries.
static void test_tsv_skip_row() {
    std::fprintf(stderr, "-- test_tsv_skip_row\n");

    std::vector<OutputHit> hits;
    hits.push_back(make_normal_hit("q1", "ACC1", 50));
    hits.push_back(make_skip_hit("q2", kSkipDegenRejected,
                                  "query contains IUPAC degenerate bases", 60));
    hits.push_back(make_normal_hit("q3", "ACC2", 30));

    std::ostringstream oss;
    write_results_tsv(oss, hits, /*mode=*/2, /*stage1_score_type=*/1,
                       /*traceback=*/false);
    std::string out = oss.str();

    CHECK(out.find("*SKIPPED:degen_rejected") != std::string::npos);
    CHECK(out.find("q1\tACC1") != std::string::npos);
    CHECK(out.find("q3\tACC2") != std::string::npos);

    // Skip row's strand is '*'
    auto skip_pos = out.find("*SKIPPED:degen_rejected");
    auto line_end = out.find('\n', skip_pos);
    auto line_begin = out.rfind('\n', skip_pos) + 1;
    std::string skip_line = out.substr(line_begin, line_end - line_begin);
    CHECK(skip_line.find("\t*\t") != std::string::npos);
}

// result_reader silently drops *SKIPPED:* rows so existing
// pipelines stay compatible.
static void test_tsv_reader_drops_skip_rows() {
    std::fprintf(stderr, "-- test_tsv_reader_drops_skip_rows\n");

    std::vector<OutputHit> hits;
    hits.push_back(make_normal_hit("q1", "ACC1", 50));
    hits.push_back(make_skip_hit("q2", kSkipDegenRejected, "...", 60));
    hits.push_back(make_normal_hit("q3", "ACC2", 30));

    std::ostringstream oss;
    write_results_tsv(oss, hits, 2, 1, false);

    std::istringstream iss(oss.str());
    auto parsed = read_results_tsv(iss);

    // Only normal hits round-trip; skip row is dropped.
    CHECK_EQ(parsed.size(), 2u);
    CHECK(parsed[0].qseqid == "q1" || parsed[1].qseqid == "q1");
    CHECK(parsed[0].qseqid == "q3" || parsed[1].qseqid == "q3");
    for (const auto& h : parsed) {
        CHECK(h.skip_reason == 0);
    }
}

// JSON output uses "status": "skipped" with reason and detail
// for skipped queries; normal queries have "status": "ok" and a hits array.
static void test_json_skip_status() {
    std::fprintf(stderr, "-- test_json_skip_status\n");

    std::vector<OutputHit> hits;
    hits.push_back(make_normal_hit("q1", "ACC1", 50));
    hits.push_back(make_skip_hit("q2", kSkipInvalidChar,
                                  "first invalid char '*' at pos 5", 60));

    std::ostringstream oss;
    write_results_json(oss, hits, 2, 1, false);
    std::string out = oss.str();

    CHECK(out.find("\"status\": \"skipped\"") != std::string::npos);
    CHECK(out.find("\"skip_reason\": \"invalid_char\"") != std::string::npos);
    CHECK(out.find("\"qlen\": 60") != std::string::npos);
    CHECK(out.find("\"status\": \"ok\"") != std::string::npos);
}

// SAM output emits an unmapped record (FLAG=4) with XR:Z and XD:Z
// aux tags carrying the skip reason and detail. Verify by reading the file
// back through htslib.
static void test_sam_skip_unmapped() {
    std::fprintf(stderr, "-- test_sam_skip_unmapped\n");
    std::string tmpdir = test_tmpdir("/tmp/test_skip_output");
    std::filesystem::create_directories(tmpdir);
    std::string sam_path = tmpdir + "/skip.sam";

    std::vector<OutputHit> hits;
    auto skip = make_skip_hit("qskip", kSkipThresholdUnreachable,
                               "Nqkmer=10 threshold=-3", 12);
    hits.push_back(skip);

    write_results_sam(sam_path, hits, 1);

    // Read back via htslib
    samFile* fp = sam_open(sam_path.c_str(), "r");
    CHECK(fp != nullptr);
    sam_hdr_t* hdr = sam_hdr_read(fp);
    CHECK(hdr != nullptr);

    bam1_t* b = bam_init1();
    int n = 0;
    bool found_unmapped = false;
    bool found_xr = false;
    bool found_xd = false;
    while (sam_read1(fp, hdr, b) >= 0) {
        n++;
        if (b->core.flag & BAM_FUNMAP) {
            found_unmapped = true;
            uint8_t* xr = bam_aux_get(b, "XR");
            uint8_t* xd = bam_aux_get(b, "XD");
            if (xr && bam_aux2Z(xr) &&
                std::string(bam_aux2Z(xr)) == "threshold_unreachable") {
                found_xr = true;
            }
            if (xd && bam_aux2Z(xd) &&
                std::string(bam_aux2Z(xd)).find("Nqkmer") != std::string::npos) {
                found_xd = true;
            }
        }
    }
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);

    CHECK_EQ(n, 1);
    CHECK(found_unmapped);
    CHECK(found_xr);
    CHECK(found_xd);

    std::filesystem::remove_all(tmpdir);
}

int main() {
    test_tsv_skip_row();
    test_tsv_reader_drops_skip_rows();
    test_json_skip_status();
    test_sam_skip_unmapped();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
