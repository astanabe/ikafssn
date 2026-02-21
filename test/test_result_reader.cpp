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
        "# query_id\taccession\tstrand\tq_start\tq_end\ts_start\ts_end\tcoverscore\tchainscore\tvolume\n"
        "query1\tACC001\t+\t0\t49\t100\t149\t5\t15\t0\n"
        "query1\tACC002\t-\t10\t39\t200\t229\t3\t10\t1\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 2u);

    CHECK(results[0].query_id == "query1");
    CHECK(results[0].accession == "ACC001");
    CHECK_EQ(results[0].strand, '+');
    CHECK_EQ(results[0].q_start, 0u);
    CHECK_EQ(results[0].q_end, 49u);
    CHECK_EQ(results[0].s_start, 100u);
    CHECK_EQ(results[0].s_end, 149u);
    CHECK_EQ(results[0].stage1_score, 5u);
    CHECK_EQ(results[0].score, 15u);
    CHECK_EQ(results[0].volume, 0u);

    CHECK(results[1].query_id == "query1");
    CHECK(results[1].accession == "ACC002");
    CHECK_EQ(results[1].strand, '-');
    CHECK_EQ(results[1].q_start, 10u);
    CHECK_EQ(results[1].q_end, 39u);
    CHECK_EQ(results[1].s_start, 200u);
    CHECK_EQ(results[1].s_end, 229u);
    CHECK_EQ(results[1].stage1_score, 3u);
    CHECK_EQ(results[1].score, 10u);
    CHECK_EQ(results[1].volume, 1u);
}

static void test_skip_header_and_blank() {
    std::fprintf(stderr, "-- test_skip_header_and_blank\n");

    std::string path = g_test_dir + "/header.tsv";
    write_file(path,
        "# comment line\n"
        "# another comment\n"
        "\n"
        "query1\tACC001\t+\t0\t49\t100\t149\t5\t15\t0\n"
        "\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 1u);
    CHECK(results[0].accession == "ACC001");
}

static void test_invalid_lines() {
    std::fprintf(stderr, "-- test_invalid_lines\n");

    std::string path = g_test_dir + "/invalid.tsv";
    write_file(path,
        "# header\n"
        "query1\tACC001\t+\t0\t49\t100\t149\t5\t15\t0\n"
        "too_few_fields\tACC002\n"
        "query2\tACC003\tX\t0\t49\t100\t149\t5\t15\t0\n"   // bad strand
        "query3\tACC004\t+\tabc\t49\t100\t149\t5\t15\t0\n"  // bad number
        "query4\tACC005\t-\t5\t55\t300\t350\t8\t20\t2\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].accession == "ACC001");
    CHECK(results[1].accession == "ACC005");
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
        "# query_id\taccession\tstrand\tq_start\tq_end\ts_start\ts_end\tscore\tvolume\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 0u);
}

static void test_stream_interface() {
    std::fprintf(stderr, "-- test_stream_interface\n");

    std::istringstream iss(
        "# header\n"
        "q1\tA1\t+\t0\t10\t20\t30\t3\t5\t0\n"
        "q2\tA2\t-\t5\t15\t25\t35\t4\t8\t1\n"
    );

    auto results = read_results_tab(iss);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].query_id == "q1");
    CHECK(results[1].query_id == "q2");
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
        CHECK(read_back[i].query_id == hits[i].query_id);
        CHECK(read_back[i].accession == hits[i].accession);
        CHECK_EQ(read_back[i].strand, hits[i].strand);
        CHECK_EQ(read_back[i].q_start, hits[i].q_start);
        CHECK_EQ(read_back[i].q_end, hits[i].q_end);
        CHECK_EQ(read_back[i].s_start, hits[i].s_start);
        CHECK_EQ(read_back[i].s_end, hits[i].s_end);
        CHECK_EQ(read_back[i].stage1_score, hits[i].stage1_score);
        CHECK_EQ(read_back[i].score, hits[i].score);
        CHECK_EQ(read_back[i].volume, hits[i].volume);
    }
}

static void test_windows_line_endings() {
    std::fprintf(stderr, "-- test_windows_line_endings\n");

    std::string path = g_test_dir + "/crlf.tsv";
    write_file(path,
        "# header\r\n"
        "q1\tA1\t+\t0\t10\t20\t30\t3\t5\t0\r\n"
    );

    auto results = read_results_tab(path);
    CHECK_EQ(results.size(), 1u);
    CHECK(results[0].query_id == "q1");
    CHECK(results[0].accession == "A1");
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
    test_windows_line_endings();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
