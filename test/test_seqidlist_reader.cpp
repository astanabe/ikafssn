#include "test_util.hpp"
#include "io/seqidlist_reader.hpp"

#include <fstream>
#include <string>
#include <filesystem>

using namespace ikafssn;

static std::string g_test_dir;

static void write_file(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    f << content;
}

static void test_basic_text() {
    std::fprintf(stderr, "-- test_basic_text\n");

    std::string path = g_test_dir + "/basic.txt";
    write_file(path, "ACC001\nACC002\nACC003\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 3u);
    CHECK(result[0] == "ACC001");
    CHECK(result[1] == "ACC002");
    CHECK(result[2] == "ACC003");
}

static void test_trim_angle_bracket() {
    std::fprintf(stderr, "-- test_trim_angle_bracket\n");

    std::string path = g_test_dir + "/bracket.txt";
    write_file(path, ">ACC001\n>ACC002\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
    CHECK(result[0] == "ACC001");
    CHECK(result[1] == "ACC002");
}

static void test_skip_blank_lines() {
    std::fprintf(stderr, "-- test_skip_blank_lines\n");

    std::string path = g_test_dir + "/blank.txt";
    write_file(path, "ACC001\n\n   \n\nACC002\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
}

static void test_skip_comments() {
    std::fprintf(stderr, "-- test_skip_comments\n");

    std::string path = g_test_dir + "/comments.txt";
    write_file(path, "# this is a comment\nACC001\n# another comment\nACC002\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
}

static void test_first_token_only() {
    std::fprintf(stderr, "-- test_first_token_only\n");

    std::string path = g_test_dir + "/tokens.txt";
    write_file(path, "ACC001 some extra text\nACC002\tmore\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
    CHECK(result[0] == "ACC001");
    CHECK(result[1] == "ACC002");
}

static void test_empty_file() {
    std::fprintf(stderr, "-- test_empty_file\n");

    std::string path = g_test_dir + "/empty.txt";
    write_file(path, "");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 0u);
}

static void test_nonexistent_file() {
    std::fprintf(stderr, "-- test_nonexistent_file\n");

    auto result = read_seqidlist("/nonexistent/path/file.txt");
    CHECK_EQ(result.size(), 0u);
}

int main() {
    g_test_dir = "/tmp/ikafssn_seqidlist_test";
    std::filesystem::create_directories(g_test_dir);

    test_basic_text();
    test_trim_angle_bracket();
    test_skip_blank_lines();
    test_skip_comments();
    test_first_token_only();
    test_empty_file();
    test_nonexistent_file();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
