#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <regex>
#include <string>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_index_dir;
static std::string g_build_dir;

static void setup() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_index_dir = "/tmp/ikafssn_info_test_index";
    g_build_dir = std::string(SOURCE_DIR) + "/build/src";

    // Create output directory
    std::filesystem::create_directories(g_index_dir);

    // Build a small index for testing using the build_index function
    BlastDbReader db;
    if (!db.open(g_testdb_path)) {
        std::fprintf(stderr, "SKIP: cannot open test BLAST DB at %s\n",
                     g_testdb_path.c_str());
        std::exit(0);
    }

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 7;
    config.partitions = 1;
    config.buffer_size = uint64_t(1) << 30;
    config.max_freq_build = 0;

    std::string prefix = g_index_dir + "/testdb.00.07mer";
    bool ok = build_index<uint16_t>(db, config, prefix, 0, 1, "testdb", logger);
    if (!ok) {
        std::fprintf(stderr, "FAIL: index build failed\n");
        std::exit(1);
    }
    db.close();
}

static void cleanup() {
    std::error_code ec;
    std::filesystem::remove_all(g_index_dir, ec);
}

// Run ikafssninfo and capture stdout
static std::string run_ikafssninfo(const std::string& args) {
    std::string cmd = g_build_dir + "/ikafssninfo " + args + " 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return {};
    std::string result;
    char buf[4096];
    while (fgets(buf, sizeof(buf), pipe)) {
        result += buf;
    }
    pclose(pipe);
    return result;
}

static int run_ikafssninfo_exit(const std::string& args, bool capture_stderr = false) {
    std::string cmd = g_build_dir + "/ikafssninfo " + args;
    if (!capture_stderr) cmd += " 2>/dev/null";
    else cmd += " 2>&1";
    int ret = std::system(cmd.c_str());
    if (WIFEXITED(ret)) return WEXITSTATUS(ret);
    return -1;
}

static void test_basic_info() {
    std::fprintf(stderr, "-- test_basic_info\n");

    std::string output = run_ikafssninfo("-ix " + g_index_dir);
    CHECK(!output.empty());

    // Check that key information is present
    CHECK(output.find("K-mer length (k):  7") != std::string::npos);
    CHECK(output.find("K-mer integer type: uint16") != std::string::npos);
    CHECK(output.find("Number of volumes: 1") != std::string::npos);
    CHECK(output.find("Volume 0:") != std::string::npos);
    CHECK(output.find("Sequences:") != std::string::npos);
    CHECK(output.find("Total postings:") != std::string::npos);
    CHECK(output.find(".kix:") != std::string::npos);
    CHECK(output.find(".kpx:") != std::string::npos);
    CHECK(output.find(".ksx:") != std::string::npos);
    CHECK(output.find("Overall Statistics") != std::string::npos);
    CHECK(output.find("Total sequences:") != std::string::npos);
    CHECK(output.find("Compression:") != std::string::npos);
}

static void test_verbose_info() {
    std::fprintf(stderr, "-- test_verbose_info\n");

    std::string output = run_ikafssninfo("-ix " + g_index_dir + " -v");
    CHECK(!output.empty());

    // Verbose mode should show frequency distribution
    CHECK(output.find("K-mer frequency distribution:") != std::string::npos);
    CHECK(output.find("Non-empty k-mers:") != std::string::npos);
    CHECK(output.find("Min count:") != std::string::npos);
    CHECK(output.find("Max count:") != std::string::npos);
    CHECK(output.find("Mean count:") != std::string::npos);
    CHECK(output.find("Percentiles:") != std::string::npos);
    CHECK(output.find("Aggregated K-mer Frequency Distribution") != std::string::npos);
}

static void test_db_info() {
    std::fprintf(stderr, "-- test_db_info\n");

    std::string output = run_ikafssninfo("-ix " + g_index_dir + " -db " + g_testdb_path);
    CHECK(!output.empty());

    // DB info section should be present
    CHECK(output.find("BLAST DB Information") != std::string::npos);
    CHECK(output.find("DB prefix:") != std::string::npos);
    CHECK(output.find("DB sequences:") != std::string::npos);
    CHECK(output.find("DB total bases:") != std::string::npos);
}

static void test_consistency_with_reader() {
    std::fprintf(stderr, "-- test_consistency_with_reader\n");

    // Verify that the info displayed by ikafssninfo matches
    // what we read directly from the index files.
    std::regex kix_pattern(R"((.+)\.(\d+)\.(\d+)mer\.kix)");

    uint32_t found_seqs = 0;
    uint64_t found_postings = 0;

    for (const auto& entry : std::filesystem::directory_iterator(g_index_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        std::smatch m;
        if (std::regex_match(fname, m, kix_pattern)) {
            KixReader kix;
            CHECK(kix.open(entry.path().string()));
            CHECK_EQ(kix.k(), 7);
            found_seqs += kix.num_sequences();
            found_postings += kix.total_postings();
            kix.close();
        }
    }

    // Now run ikafssninfo and verify the totals match
    std::string output = run_ikafssninfo("-ix " + g_index_dir);
    std::string seq_str = "Total sequences:   " + std::to_string(found_seqs);
    CHECK(output.find(seq_str) != std::string::npos);
}

static void test_help() {
    std::fprintf(stderr, "-- test_help\n");

    std::string cmd = g_build_dir + "/ikafssninfo --help 2>&1";
    FILE* pipe = popen(cmd.c_str(), "r");
    CHECK(pipe != nullptr);
    std::string result;
    char buf[4096];
    while (fgets(buf, sizeof(buf), pipe)) {
        result += buf;
    }
    int ret = pclose(pipe);
    CHECK(WIFEXITED(ret) && WEXITSTATUS(ret) == 0);
    CHECK(result.find("Usage:") != std::string::npos);
    CHECK(result.find("-ix") != std::string::npos);
    CHECK(result.find("-db") != std::string::npos);
}

static void test_missing_ix() {
    std::fprintf(stderr, "-- test_missing_ix\n");

    // Should return non-zero exit code when -ix is missing
    CHECK(run_ikafssninfo_exit("") != 0);
}

static void test_nonexistent_dir() {
    std::fprintf(stderr, "-- test_nonexistent_dir\n");

    CHECK(run_ikafssninfo_exit("-ix /nonexistent/path") != 0);
}

int main() {
    setup();

    test_basic_info();
    test_verbose_info();
    test_db_info();
    test_consistency_with_reader();
    test_help();
    test_missing_ix();
    test_nonexistent_dir();

    cleanup();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
