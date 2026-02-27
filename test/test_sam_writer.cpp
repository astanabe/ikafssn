#include "test_util.hpp"

#include "io/sam_writer.hpp"
#include "io/result_writer.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <htslib/sam.h>
#include <htslib/hts.h>

using namespace ikafssn;

static std::string g_test_dir;

static std::vector<OutputHit> make_test_hits() {
    std::vector<OutputHit> hits;

    OutputHit h1;
    h1.qseqid = "query1";
    h1.sseqid = "FJ876973";
    h1.sstrand = '+';
    h1.qstart = 0;
    h1.qend = 19;
    h1.sstart = 100;
    h1.send = 119;
    h1.chainscore = 50;
    h1.coverscore = 30;
    h1.volume = 0;
    h1.alnscore = 100;
    h1.cigar = "20=";
    h1.nident = 20;
    h1.mismatch = 0;
    h1.pident = 100.0;
    h1.qseq = "ACGTACGTACGTACGTACGT";
    h1.sseq = "ACGTACGTACGTACGTACGT";
    h1.slen = 1800;
    hits.push_back(h1);

    OutputHit h2;
    h2.qseqid = "query1";
    h2.sseqid = "GQ912721";
    h2.sstrand = '-';
    h2.qstart = 0;
    h2.qend = 14;
    h2.sstart = 200;
    h2.send = 216;
    h2.chainscore = 30;
    h2.coverscore = 20;
    h2.volume = 0;
    h2.alnscore = 60;
    h2.cigar = "10=2D5=";
    h2.nident = 15;
    h2.mismatch = 0;
    h2.pident = 88.2;
    h2.qseq = "ACGTACGTAC--ACGTC";
    h2.sseq = "ACGTACGTACTTACGTC";
    h2.slen = 1700;
    hits.push_back(h2);

    return hits;
}

static void test_sam_basic() {
    std::fprintf(stderr, "-- test_sam_basic\n");

    auto hits = make_test_hits();
    std::string sam_path = g_test_dir + "/test_output.sam";

    write_results_sam(sam_path, hits, 1);

    // Read back and verify
    std::ifstream in(sam_path);
    CHECK(in.is_open());

    std::string line;
    bool has_hd = false, has_sq = false, has_pg = false;
    int record_count = 0;

    while (std::getline(in, line)) {
        if (line.substr(0, 3) == "@HD") has_hd = true;
        if (line.substr(0, 3) == "@SQ") has_sq = true;
        if (line.substr(0, 3) == "@PG") has_pg = true;
        if (!line.empty() && line[0] != '@') record_count++;
    }

    CHECK(has_hd);
    CHECK(has_sq);
    CHECK(has_pg);
    CHECK(record_count == 2);
}

static void test_sam_reverse_flag() {
    std::fprintf(stderr, "-- test_sam_reverse_flag\n");

    auto hits = make_test_hits();
    std::string sam_path = g_test_dir + "/test_reverse.sam";

    write_results_sam(sam_path, hits, 1);

    // Parse with htslib to check flags
    samFile* fp = sam_open(sam_path.c_str(), "r");
    CHECK(fp != nullptr);

    sam_hdr_t* hdr = sam_hdr_read(fp);
    CHECK(hdr != nullptr);

    bam1_t* b = bam_init1();
    int rec = 0;
    while (sam_read1(fp, hdr, b) >= 0) {
        if (rec == 0) {
            // First hit is forward
            CHECK((b->core.flag & BAM_FREVERSE) == 0);
        } else if (rec == 1) {
            // Second hit is reverse
            CHECK((b->core.flag & BAM_FREVERSE) != 0);
        }
        rec++;
    }
    CHECK(rec == 2);

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);
}

static void test_bam_roundtrip() {
    std::fprintf(stderr, "-- test_bam_roundtrip\n");

    auto hits = make_test_hits();
    std::string bam_path = g_test_dir + "/test_roundtrip.bam";

    write_results_bam(bam_path, hits, 1);

    // Read back with htslib
    samFile* fp = sam_open(bam_path.c_str(), "r");
    CHECK(fp != nullptr);

    sam_hdr_t* hdr = sam_hdr_read(fp);
    CHECK(hdr != nullptr);

    bam1_t* b = bam_init1();
    int rec_count = 0;
    while (sam_read1(fp, hdr, b) >= 0) {
        rec_count++;
        // Check AS tag exists
        uint8_t* as = bam_aux_get(b, "AS");
        CHECK(as != nullptr);
        if (as) {
            int32_t as_val = bam_aux2i(as);
            CHECK(as_val > 0);
        }
    }
    CHECK(rec_count == 2);

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);
}

static void test_cigar_encoding() {
    std::fprintf(stderr, "-- test_cigar_encoding\n");

    auto hits = make_test_hits();
    std::string sam_path = g_test_dir + "/test_cigar.sam";

    write_results_sam(sam_path, hits, 1);

    // Parse with htslib to check CIGAR
    samFile* fp = sam_open(sam_path.c_str(), "r");
    CHECK(fp != nullptr);

    sam_hdr_t* hdr = sam_hdr_read(fp);
    CHECK(hdr != nullptr);

    bam1_t* b = bam_init1();
    int rec = 0;
    while (sam_read1(fp, hdr, b) >= 0) {
        if (rec == 0) {
            // First hit: "20=" -> 1 CIGAR op
            CHECK(b->core.n_cigar == 1);
            if (b->core.n_cigar >= 1) {
                uint32_t* cigar = bam_get_cigar(b);
                CHECK(bam_cigar_op(cigar[0]) == BAM_CEQUAL);
                CHECK(bam_cigar_oplen(cigar[0]) == 20);
            }
        } else if (rec == 1) {
            // Second hit: "10=2I5=" -> 3 CIGAR ops
            CHECK(b->core.n_cigar == 3);
        }
        rec++;
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    sam_close(fp);
}

int main() {
    // Create temp directory
    g_test_dir = "/tmp/ikafssn_test_sam_writer";
    std::filesystem::create_directories(g_test_dir);

    test_sam_basic();
    test_sam_reverse_flag();
    test_bam_roundtrip();
    test_cigar_encoding();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
