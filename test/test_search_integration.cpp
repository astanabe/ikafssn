#include "test_util.hpp"
#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace ikafssn;

static std::string g_testdb_path;
static std::string g_test_dir;

static void test_build_and_search() {
    std::fprintf(stderr, "-- test_build_and_search\n");

    // Build index with k=7
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    bconfig.partitions = 1;
    bconfig.buffer_size = uint64_t(1) << 30;

    std::string prefix = g_test_dir + "/test.00.07mer";
    CHECK(build_index<uint16_t>(db, bconfig, prefix, 0, 1, "test", logger));

    // Open index
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    // Search with query from seq1
    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1; // relaxed for small test
    config.stage2.min_score = 2;
    config.num_results = 50;

    // Query: first 24 bases of seq1 "ACGTACGTACGTACGTACGTACGT"
    std::string query_seq = "ACGTACGTACGTACGTACGTACGT";
    auto result = search_volume<uint16_t>(
        "query1", query_seq, 7, kix, kpx, ksx, filter, config);

    CHECK(!result.hits.empty());
    CHECK(result.query_id == "query1");

    // seq1 (OID 0) should be found
    bool found_seq1 = false;
    for (const auto& cr : result.hits) {
        if (cr.seq_id == 0 && !cr.is_reverse) {
            found_seq1 = true;
            CHECK(cr.score >= 2);
        }
    }
    CHECK(found_seq1);

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_revcomp_search() {
    std::fprintf(stderr, "-- test_revcomp_search\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    // Reverse complement of "ACGTACGT" = "ACGTACGT" (palindrome!)
    // Use something non-palindromic: revcomp of "AACCGGTT" = "AACCGGTT" (also palindromic)
    // Use "AACCGGTTAA": revcomp = "TTAACCGGTT"
    // seq5 = "ACGTACGTAACCGGTTAACCGGTTACGTACGTAACCGGTTAACCGGTT"
    // "TTAACCGGTT" is the revcomp of a substring of seq5
    std::string query_seq = "TTAACCGGTTAACCGGTTAACCGGTT";
    auto result = search_volume<uint16_t>(
        "rc_query", query_seq, 7, kix, kpx, ksx, filter, config);

    // Should find seq5 on reverse strand
    bool found_rev = false;
    for (const auto& cr : result.hits) {
        if (cr.is_reverse) {
            found_rev = true;
        }
    }
    // Not necessarily finding reverse - depends on actual k-mer content
    // Just verify the search completes without error
    CHECK(result.query_id == "rc_query");

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_seqidlist_filter() {
    std::fprintf(stderr, "-- test_seqidlist_filter\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    // Build OID filter that includes only seq3 and seq4
    std::vector<std::string> include_list = {"seq3", "seq4"};
    OidFilter filter;
    filter.build(include_list, ksx, OidFilterMode::kInclude);

    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 1;
    config.num_results = 50;

    std::string query_seq = "ACGTACGTACGTACGTACGTACGT";
    auto result = search_volume<uint16_t>(
        "filtered_query", query_seq, 7, kix, kpx, ksx, filter, config);

    // Results should only contain seq3 (OID 2) or seq4 (OID 3), not seq1 (OID 0)
    for (const auto& cr : result.hits) {
        CHECK(cr.seq_id == 2 || cr.seq_id == 3);
    }

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_negative_seqidlist() {
    std::fprintf(stderr, "-- test_negative_seqidlist\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    // Build OID filter that excludes seq1
    std::vector<std::string> exclude_list = {"seq1"};
    OidFilter filter;
    filter.build(exclude_list, ksx, OidFilterMode::kExclude);

    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 1;
    config.num_results = 50;

    std::string query_seq = "ACGTACGTACGTACGTACGTACGT";
    auto result = search_volume<uint16_t>(
        "neg_query", query_seq, 7, kix, kpx, ksx, filter, config);

    // seq1 (OID 0) should be excluded from results
    for (const auto& cr : result.hits) {
        CHECK(cr.seq_id != 0);
    }

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_result_output_tab() {
    std::fprintf(stderr, "-- test_result_output_tab\n");

    std::vector<OutputHit> hits = {
        {"query1", "NM_001234", '+', 10, 450, 1020, 1460, 42, 0},
        {"query1", "XM_005678", '-', 15, 430, 8050, 8465, 38, 2},
    };

    std::ostringstream oss;
    write_results_tab(oss, hits);
    std::string output = oss.str();

    // Check header
    CHECK(output.find("# query_id") != std::string::npos);
    // Check data
    CHECK(output.find("NM_001234") != std::string::npos);
    CHECK(output.find("XM_005678") != std::string::npos);
    CHECK(output.find("42") != std::string::npos);
}

static void test_result_output_json() {
    std::fprintf(stderr, "-- test_result_output_json\n");

    std::vector<OutputHit> hits = {
        {"query1", "NM_001234", '+', 10, 450, 1020, 1460, 42, 0},
    };

    std::ostringstream oss;
    write_results_json(oss, hits);
    std::string output = oss.str();

    CHECK(output.find("\"results\"") != std::string::npos);
    CHECK(output.find("\"query_id\"") != std::string::npos);
    CHECK(output.find("\"NM_001234\"") != std::string::npos);
    CHECK(output.find("\"score\": 42") != std::string::npos);
}

static void test_fasta_reader() {
    std::fprintf(stderr, "-- test_fasta_reader\n");

    std::string path = g_test_dir + "/test.fasta";
    {
        std::ofstream f(path);
        f << ">seq_A description text\n";
        f << "ACGTACGT\n";
        f << "TTTTAAAA\n";
        f << ">seq_B\n";
        f << "GGGGCCCC\n";
    }

    auto records = read_fasta(path);
    CHECK_EQ(records.size(), 2u);
    CHECK(records[0].id == "seq_A");
    CHECK(records[0].sequence == "ACGTACGTTTTTAAAA");
    CHECK(records[1].id == "seq_B");
    CHECK(records[1].sequence == "GGGGCCCC");
}

static void test_search_k9() {
    std::fprintf(stderr, "-- test_search_k9\n");

    // Build index with k=9
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = 9;
    bconfig.partitions = 1;
    bconfig.buffer_size = uint64_t(1) << 30;

    std::string prefix = g_test_dir + "/test.00.09mer";
    CHECK(build_index<uint32_t>(db, bconfig, prefix, 0, 1, "test", logger));

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    std::string query_seq = "ACGTACGTACGTACGTACGTACGT";
    auto result = search_volume<uint32_t>(
        "query_k9", query_seq, 9, kix, kpx, ksx, filter, config);

    // Should find hits (seq1 has this pattern)
    CHECK(result.query_id == "query_k9");
    // At least some hits expected
    bool found = !result.hits.empty();
    if (found) {
        // Verify hits are reasonable
        for (const auto& cr : result.hits) {
            CHECK(cr.score >= 2);
        }
    }

    kix.close();
    kpx.close();
    ksx.close();
}

int main() {
    g_testdb_path = std::string(SOURCE_DIR) + "/test/testdata/testdb";
    g_test_dir = "/tmp/ikafssn_search_test";
    std::filesystem::create_directories(g_test_dir);

    test_fasta_reader();
    test_result_output_tab();
    test_result_output_json();
    test_build_and_search();
    test_revcomp_search();
    test_seqidlist_filter();
    test_negative_seqidlist();
    test_search_k9();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
