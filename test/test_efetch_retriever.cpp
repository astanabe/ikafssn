#include "test_util.hpp"
#include "ikafssnretrieve/efetch_retriever.hpp"

#include <string>
#include <vector>

using namespace ikafssn;

static void test_build_url_batch() {
    std::fprintf(stderr, "-- test_build_url_batch\n");

    std::vector<std::string> accs = {"ACC001", "ACC002", "ACC003"};
    std::string url = build_efetch_url_batch(accs, "");

    CHECK(url.find("db=nuccore") != std::string::npos);
    CHECK(url.find("rettype=fasta") != std::string::npos);
    CHECK(url.find("id=ACC001,ACC002,ACC003") != std::string::npos);
    CHECK(url.find("api_key") == std::string::npos);
}

static void test_build_url_batch_with_key() {
    std::fprintf(stderr, "-- test_build_url_batch_with_key\n");

    std::vector<std::string> accs = {"ACC001"};
    std::string url = build_efetch_url_batch(accs, "mykey123");

    CHECK(url.find("id=ACC001") != std::string::npos);
    CHECK(url.find("api_key=mykey123") != std::string::npos);
}

static void test_build_url_range() {
    std::fprintf(stderr, "-- test_build_url_range\n");

    std::string url = build_efetch_url_range("ACC001", 100, 500, "");

    CHECK(url.find("id=ACC001") != std::string::npos);
    CHECK(url.find("seq_start=100") != std::string::npos);
    CHECK(url.find("seq_stop=500") != std::string::npos);
    CHECK(url.find("api_key") == std::string::npos);
}

static void test_build_url_range_with_key() {
    std::fprintf(stderr, "-- test_build_url_range_with_key\n");

    std::string url = build_efetch_url_range("ACC001", 1, 1000, "key456");

    CHECK(url.find("seq_start=1") != std::string::npos);
    CHECK(url.find("seq_stop=1000") != std::string::npos);
    CHECK(url.find("api_key=key456") != std::string::npos);
}

static void test_parse_response_single() {
    std::fprintf(stderr, "-- test_parse_response_single\n");

    std::string response =
        ">ACC001.1 Some description\n"
        "ATCGATCGATCG\n"
        "TTTTAAAA\n";

    auto records = parse_efetch_response(response);
    CHECK_EQ(records.size(), 1u);
    CHECK(records[0].accession == "ACC001");
    CHECK(records[0].sequence == "ATCGATCGATCGTTTTAAAA");
}

static void test_parse_response_multiple() {
    std::fprintf(stderr, "-- test_parse_response_multiple\n");

    std::string response =
        ">ACC001.1 First sequence\n"
        "AAAA\n"
        ">ACC002.2 Second sequence\n"
        "CCCC\n"
        ">ACC003.1 Third\n"
        "GGGG\n"
        "TTTT\n";

    auto records = parse_efetch_response(response);
    CHECK_EQ(records.size(), 3u);
    CHECK(records[0].accession == "ACC001");
    CHECK(records[0].sequence == "AAAA");
    CHECK(records[1].accession == "ACC002");
    CHECK(records[1].sequence == "CCCC");
    CHECK(records[2].accession == "ACC003");
    CHECK(records[2].sequence == "GGGGTTTT");
}

static void test_parse_response_gi_format() {
    std::fprintf(stderr, "-- test_parse_response_gi_format\n");

    std::string response =
        ">gi|12345|gb|ACC001.1| Some description\n"
        "ATCG\n";

    auto records = parse_efetch_response(response);
    CHECK_EQ(records.size(), 1u);
    CHECK(records[0].accession == "ACC001");
    CHECK(records[0].sequence == "ATCG");
}

static void test_parse_response_empty() {
    std::fprintf(stderr, "-- test_parse_response_empty\n");

    auto records = parse_efetch_response("");
    CHECK_EQ(records.size(), 0u);
}

static void test_parse_response_lowercase() {
    std::fprintf(stderr, "-- test_parse_response_lowercase\n");

    std::string response =
        ">ACC001.1 test\n"
        "atcg\n"
        "NNNN\n";

    auto records = parse_efetch_response(response);
    CHECK_EQ(records.size(), 1u);
    CHECK(records[0].sequence == "ATCGNNNN");
}

static void test_rate_limit() {
    std::fprintf(stderr, "-- test_rate_limit\n");

    CHECK_EQ(rate_limit_sleep_ms(false), 334u);
    CHECK_EQ(rate_limit_sleep_ms(true), 100u);
}

static void test_retryable_status() {
    std::fprintf(stderr, "-- test_retryable_status\n");

    CHECK(is_retryable_http_status(429));
    CHECK(is_retryable_http_status(503));
    CHECK(!is_retryable_http_status(200));
    CHECK(!is_retryable_http_status(400));
    CHECK(!is_retryable_http_status(404));
    CHECK(!is_retryable_http_status(500));
}

static void test_skip_status() {
    std::fprintf(stderr, "-- test_skip_status\n");

    CHECK(is_skip_http_status(400));
    CHECK(is_skip_http_status(404));
    CHECK(!is_skip_http_status(200));
    CHECK(!is_skip_http_status(429));
    CHECK(!is_skip_http_status(503));
}

int main() {
    test_build_url_batch();
    test_build_url_batch_with_key();
    test_build_url_range();
    test_build_url_range_with_key();
    test_parse_response_single();
    test_parse_response_multiple();
    test_parse_response_gi_format();
    test_parse_response_empty();
    test_parse_response_lowercase();
    test_rate_limit();
    test_retryable_status();
    test_skip_status();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
