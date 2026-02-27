#include "test_util.hpp"
#include "io/result_reader.hpp"
#include "io/result_writer.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace ikafssn;

static std::string g_test_dir;

static void write_file(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    f << content;
}

static void test_basic_parse() {
    std::fprintf(stderr, "-- test_basic_parse\n");

    std::string path = g_test_dir + "/basic.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "query1\tACC001\t+\t0\t49\t500\t100\t149\t2000\t5\t15\t0\n"
        "query1\tACC002\t-\t10\t39\t500\t200\t229\t3000\t3\t10\t1\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 2u);

    CHECK(results[0].qseqid == "query1");
    CHECK(results[0].sseqid == "ACC001");
    CHECK_EQ(results[0].sstrand, '+');
    CHECK_EQ(results[0].qstart, 0u);
    CHECK_EQ(results[0].qend, 49u);
    CHECK_EQ(results[0].qlen, 500u);
    CHECK_EQ(results[0].sstart, 100u);
    CHECK_EQ(results[0].send, 149u);
    CHECK_EQ(results[0].slen, 2000u);
    CHECK_EQ(results[0].coverscore, 5u);
    CHECK_EQ(results[0].chainscore, 15u);
    CHECK_EQ(results[0].volume, 0u);

    CHECK(results[1].qseqid == "query1");
    CHECK(results[1].sseqid == "ACC002");
    CHECK_EQ(results[1].sstrand, '-');
    CHECK_EQ(results[1].qstart, 10u);
    CHECK_EQ(results[1].qend, 39u);
    CHECK_EQ(results[1].qlen, 500u);
    CHECK_EQ(results[1].sstart, 200u);
    CHECK_EQ(results[1].send, 229u);
    CHECK_EQ(results[1].slen, 3000u);
    CHECK_EQ(results[1].coverscore, 3u);
    CHECK_EQ(results[1].chainscore, 10u);
    CHECK_EQ(results[1].volume, 1u);
}

static void test_skip_header_and_blank() {
    std::fprintf(stderr, "-- test_skip_header_and_blank\n");

    std::string path = g_test_dir + "/header.tsv";
    write_file(path,
        "# comment line\n"
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "\n"
        "query1\tACC001\t+\t0\t49\t500\t100\t149\t2000\t5\t15\t0\n"
        "\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 1u);
    CHECK(results[0].sseqid == "ACC001");
}

static void test_invalid_lines() {
    std::fprintf(stderr, "-- test_invalid_lines\n");

    std::string path = g_test_dir + "/invalid.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "query1\tACC001\t+\t0\t49\t500\t100\t149\t2000\t5\t15\t0\n"
        "too_few_fields\tACC002\n"
        "query2\tACC003\tX\t0\t49\t500\t100\t149\t2000\t5\t15\t0\n"   // bad strand
        "query3\tACC004\t+\tabc\t49\t500\t100\t149\t2000\t5\t15\t0\n"  // bad number
        "query4\tACC005\t-\t5\t55\t500\t300\t350\t2000\t8\t20\t2\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].sseqid == "ACC001");
    CHECK(results[1].sseqid == "ACC005");
}

static void test_empty_input() {
    std::fprintf(stderr, "-- test_empty_input\n");

    std::string path = g_test_dir + "/empty.tsv";
    write_file(path, "");

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 0u);
}

static void test_header_only() {
    std::fprintf(stderr, "-- test_header_only\n");

    std::string path = g_test_dir + "/header_only.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tsstart\tsend\tscore\tvolume\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 0u);
}

static void test_stream_interface() {
    std::fprintf(stderr, "-- test_stream_interface\n");

    std::istringstream iss(
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "q1\tA1\t+\t0\t10\t100\t20\t30\t200\t3\t5\t0\n"
        "q2\tA2\t-\t5\t15\t100\t25\t35\t200\t4\t8\t1\n"
    );

    auto results = read_results_tab(iss);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].qseqid == "q1");
    CHECK(results[1].qseqid == "q2");
}

static void test_roundtrip() {
    std::fprintf(stderr, "-- test_roundtrip\n");

    // Create hits, write them, then read back
    std::vector<OutputHit> hits;
    hits.push_back({"queryA", "ACC100", '+', 0, 99, 500, 599, 25, 8, 0});
    hits.push_back({"queryA", "ACC200", '-', 10, 89, 1000, 1079, 18, 5, 1});
    hits.push_back({"queryB", "ACC300", '+', 0, 49, 0, 49, 12, 3, 0});

    std::ostringstream oss;
    write_results_tab(oss, hits);

    std::istringstream iss(oss.str());
    auto read_back = read_results_tab(iss);

    CHECK_EQ(read_back.size(), 3u);
    for (size_t i = 0; i < 3; i++) {
        CHECK(read_back[i].qseqid == hits[i].qseqid);
        CHECK(read_back[i].sseqid == hits[i].sseqid);
        CHECK_EQ(read_back[i].sstrand, hits[i].sstrand);
        CHECK_EQ(read_back[i].qstart, hits[i].qstart);
        CHECK_EQ(read_back[i].qend, hits[i].qend);
        CHECK_EQ(read_back[i].sstart, hits[i].sstart);
        CHECK_EQ(read_back[i].send, hits[i].send);
        CHECK_EQ(read_back[i].coverscore, hits[i].coverscore);
        CHECK_EQ(read_back[i].chainscore, hits[i].chainscore);
        CHECK_EQ(read_back[i].volume, hits[i].volume);
    }
}

static void test_roundtrip_mode3_no_traceback() {
    std::fprintf(stderr, "-- test_roundtrip_mode3_no_traceback\n");

    // Mode 3, no traceback: qstart/sstart should NOT be in output
    std::vector<OutputHit> hits;
    OutputHit h;
    h.qseqid = "qryM3";
    h.sseqid = "ACC_M3";
    h.sstrand = '+';
    h.qstart = 5;    // will NOT be written
    h.qend = 95;
    h.qlen = 100;
    h.sstart = 200;  // will NOT be written
    h.send = 800;
    h.slen = 5000;
    h.coverscore = 10;
    h.chainscore = 50;
    h.alnscore = 120;
    h.volume = 2;
    hits.push_back(h);

    std::ostringstream oss;
    write_results_tab(oss, hits, /*mode=*/3, /*stage1_score_type=*/1, /*stage3_traceback=*/false);

    // Verify qstart/sstart are not in header
    std::string output = oss.str();
    CHECK(output.find("qstart") == std::string::npos);
    CHECK(output.find("sstart") == std::string::npos);
    CHECK(output.find("qend") != std::string::npos);
    CHECK(output.find("send") != std::string::npos);

    // Read back
    std::istringstream iss(output);
    auto read_back = read_results_tab(iss);
    CHECK_EQ(read_back.size(), 1u);
    CHECK(read_back[0].qseqid == "qryM3");
    CHECK(read_back[0].sseqid == "ACC_M3");
    CHECK_EQ(read_back[0].sstrand, '+');
    CHECK_EQ(read_back[0].qstart, 0u);   // not present -> default 0
    CHECK_EQ(read_back[0].qend, 95u);
    CHECK_EQ(read_back[0].qlen, 100u);
    CHECK_EQ(read_back[0].sstart, 0u);   // not present -> default 0
    CHECK_EQ(read_back[0].send, 800u);
    CHECK_EQ(read_back[0].slen, 5000u);
    CHECK_EQ(read_back[0].coverscore, 10u);
    CHECK_EQ(read_back[0].chainscore, 50u);
    CHECK_EQ(read_back[0].alnscore, 120);
    CHECK_EQ(read_back[0].volume, 2u);
}

static void test_roundtrip_mode3_traceback() {
    std::fprintf(stderr, "-- test_roundtrip_mode3_traceback\n");

    // Mode 3, with traceback: qstart/sstart SHOULD be present
    std::vector<OutputHit> hits;
    OutputHit h;
    h.qseqid = "qryTB";
    h.sseqid = "ACC_TB";
    h.sstrand = '-';
    h.qstart = 3;
    h.qend = 97;
    h.qlen = 100;
    h.sstart = 150;
    h.send = 750;
    h.slen = 3000;
    h.coverscore = 7;
    h.chainscore = 40;
    h.alnscore = 200;
    h.pident = 95.5;
    h.nident = 90;
    h.mismatch = 4;
    h.cigar = "50M2I48M";
    h.qseq = "ACGT";
    h.sseq = "ACGT";
    h.volume = 1;
    hits.push_back(h);

    std::ostringstream oss;
    write_results_tab(oss, hits, /*mode=*/3, /*stage1_score_type=*/1, /*stage3_traceback=*/true);

    std::string output = oss.str();
    CHECK(output.find("qstart") != std::string::npos);
    CHECK(output.find("sstart") != std::string::npos);

    std::istringstream iss(output);
    auto read_back = read_results_tab(iss);
    CHECK_EQ(read_back.size(), 1u);
    CHECK_EQ(read_back[0].qstart, 3u);
    CHECK_EQ(read_back[0].qend, 97u);
    CHECK_EQ(read_back[0].sstart, 150u);
    CHECK_EQ(read_back[0].send, 750u);
    CHECK_EQ(read_back[0].alnscore, 200);
    CHECK(read_back[0].cigar == "50M2I48M");
    CHECK(read_back[0].qseq == "ACGT");
    CHECK(read_back[0].sseq == "ACGT");
}

static void test_header_reordered_columns() {
    std::fprintf(stderr, "-- test_header_reordered_columns\n");

    // Header with columns in a different order
    std::string path = g_test_dir + "/reordered.tsv";
    write_file(path,
        "# sseqid\tqseqid\tvolume\tsstrand\tcoverscore\tqlen\tslen\n"
        "ACC_R1\tqR1\t3\t+\t12\t400\t5000\n"
        "ACC_R2\tqR2\t0\t-\t8\t300\t4000\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 2u);

    CHECK(results[0].qseqid == "qR1");
    CHECK(results[0].sseqid == "ACC_R1");
    CHECK_EQ(results[0].sstrand, '+');
    CHECK_EQ(results[0].coverscore, 12u);
    CHECK_EQ(results[0].qlen, 400u);
    CHECK_EQ(results[0].slen, 5000u);
    CHECK_EQ(results[0].volume, 3u);
    // Missing columns default to 0
    CHECK_EQ(results[0].qstart, 0u);
    CHECK_EQ(results[0].sstart, 0u);
    CHECK_EQ(results[0].chainscore, 0u);

    CHECK(results[1].qseqid == "qR2");
    CHECK(results[1].sseqid == "ACC_R2");
    CHECK_EQ(results[1].sstrand, '-');
    CHECK_EQ(results[1].coverscore, 8u);
    CHECK_EQ(results[1].volume, 0u);
}

static void test_header_matchscore() {
    std::fprintf(stderr, "-- test_header_matchscore\n");

    // Test with matchscore instead of coverscore
    std::string path = g_test_dir + "/matchscore.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqlen\tslen\tmatchscore\tvolume\n"
        "qM1\tACC_MS\t+\t100\t2000\t42\t0\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 1u);
    CHECK_EQ(results[0].matchscore, 42u);
}

static void test_legacy_no_header() {
    std::fprintf(stderr, "-- test_legacy_no_header\n");

    // No valid header line — should fall back to legacy field-count parser
    std::istringstream iss(
        "q1\tA1\t+\t0\t10\t100\t20\t30\t200\t3\t5\t0\n"
        "q2\tA2\t-\t5\t15\t100\t25\t35\t200\t4\t8\t1\n"
    );

    auto results = read_results_tab(iss);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].qseqid == "q1");
    CHECK_EQ(results[0].qstart, 0u);
    CHECK_EQ(results[0].qend, 10u);
    CHECK(results[1].qseqid == "q2");
    CHECK_EQ(results[1].qstart, 5u);
    CHECK_EQ(results[1].qend, 15u);
}

static void test_windows_line_endings() {
    std::fprintf(stderr, "-- test_windows_line_endings\n");

    std::string path = g_test_dir + "/crlf.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\r\n"
        "q1\tA1\t+\t0\t10\t100\t20\t30\t200\t3\t5\t0\r\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 1u);
    CHECK(results[0].qseqid == "q1");
    CHECK(results[0].sseqid == "A1");
}

int main() {
    g_test_dir = "/tmp/ikafssn_result_reader_test";
    std::filesystem::create_directories(g_test_dir);

    test_basic_parse();
    test_skip_header_and_blank();
    test_invalid_lines();
    test_empty_input();
    test_header_only();
    test_stream_interface();
    test_roundtrip();
    test_roundtrip_mode3_no_traceback();
    test_roundtrip_mode3_traceback();
    test_header_reordered_columns();
    test_header_matchscore();
    test_legacy_no_header();
    test_windows_line_endings();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
