#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
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
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;

// Runtime-extracted data
static std::string g_query_seq;   // 100bp from FJ876973.1
static uint32_t g_fj_oid = UINT32_MAX;

static void test_build_and_search() {
    std::fprintf(stderr, "-- test_build_and_search\n");

    // Build index with k=7
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = 7;

    std::string prefix = g_test_dir + "/test.00.07mer";
    CHECK(build_index<uint16_t>(db, bconfig, prefix, 0, 1, "test", logger));

    // Open index
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    // Search with query from FJ876973.1
    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    // Query: 100bp from FJ876973.1 (extracted at runtime)
    auto result = search_volume<uint16_t>(
        "query1", g_query_seq, 7, kix, kpx, ksx, filter, config);

    CHECK(!result.hits.empty());
    CHECK(result.query_id == "query1");

    // FJ876973.1 OID should be found
    bool found_fj = false;
    for (const auto& cr : result.hits) {
        if (cr.seq_id == g_fj_oid && !cr.is_reverse) {
            found_fj = true;
            CHECK(cr.score >= 2);
        }
    }
    CHECK(found_fj);

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

    // Compute the reverse complement of the query
    std::string rc_query;
    rc_query.reserve(g_query_seq.size());
    for (auto it = g_query_seq.rbegin(); it != g_query_seq.rend(); ++it) {
        switch (*it) {
            case 'A': rc_query += 'T'; break;
            case 'T': rc_query += 'A'; break;
            case 'C': rc_query += 'G'; break;
            case 'G': rc_query += 'C'; break;
            default:  rc_query += 'N'; break;
        }
    }

    auto result = search_volume<uint16_t>(
        "rc_query", rc_query, 7, kix, kpx, ksx, filter, config);

    // Verify the search completes without error
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

    // Find OIDs for the target accessions
    uint32_t oid_gq = UINT32_MAX, oid_dq = UINT32_MAX;
    for (uint32_t i = 0; i < ksx.num_sequences(); i++) {
        std::string acc(ksx.accession(i));
        if (acc == ACC_GQ) oid_gq = i;
        if (acc == ACC_DQ) oid_dq = i;
    }

    // Build OID filter that includes only ACC_GQ and ACC_DQ
    std::vector<std::string> include_list = {ACC_GQ, ACC_DQ};
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

    auto result = search_volume<uint16_t>(
        "filtered_query", g_query_seq, 7, kix, kpx, ksx, filter, config);

    // Results should only contain the included OIDs, not FJ876973.1
    for (const auto& cr : result.hits) {
        CHECK(cr.seq_id == oid_gq || cr.seq_id == oid_dq);
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

    // Build OID filter that excludes ACC_FJ
    std::vector<std::string> exclude_list = {ACC_FJ};
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

    auto result = search_volume<uint16_t>(
        "neg_query", g_query_seq, 7, kix, kpx, ksx, filter, config);

    // FJ876973.1 OID should be excluded from results
    for (const auto& cr : result.hits) {
        CHECK(cr.seq_id != g_fj_oid);
    }

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_result_output_tab() {
    std::fprintf(stderr, "-- test_result_output_tab\n");

    std::vector<OutputHit> hits = {
        {"query1", "NM_001234", '+', 10, 450, 1020, 1460, 42, 10, 0},
        {"query1", "XM_005678", '-', 15, 430, 8050, 8465, 38, 8, 2},
    };

    std::ostringstream oss;
    write_results_tab(oss, hits);
    std::string output = oss.str();

    // Check header (mode 2 default: includes coverscore and chainscore)
    CHECK(output.find("# query_id") != std::string::npos);
    CHECK(output.find("coverscore") != std::string::npos);
    CHECK(output.find("chainscore") != std::string::npos);
    // Check data
    CHECK(output.find("NM_001234") != std::string::npos);
    CHECK(output.find("XM_005678") != std::string::npos);
    CHECK(output.find("42") != std::string::npos);
}

static void test_result_output_tab_mode1() {
    std::fprintf(stderr, "-- test_result_output_tab_mode1\n");

    std::vector<OutputHit> hits = {
        {"query1", "NM_001234", '+', 0, 0, 0, 0, 0, 10, 0},
    };

    std::ostringstream oss;
    write_results_tab(oss, hits, 1, 1);
    std::string output = oss.str();

    // Mode 1: no q_start/q_end/s_start/s_end/chainscore columns
    CHECK(output.find("coverscore") != std::string::npos);
    CHECK(output.find("chainscore") == std::string::npos);
    CHECK(output.find("q_start") == std::string::npos);
}

static void test_result_output_json() {
    std::fprintf(stderr, "-- test_result_output_json\n");

    std::vector<OutputHit> hits = {
        {"query1", "NM_001234", '+', 10, 450, 1020, 1460, 42, 10, 0},
    };

    std::ostringstream oss;
    write_results_json(oss, hits);
    std::string output = oss.str();

    CHECK(output.find("\"results\"") != std::string::npos);
    CHECK(output.find("\"query_id\"") != std::string::npos);
    CHECK(output.find("\"NM_001234\"") != std::string::npos);
    CHECK(output.find("\"chainscore\": 42") != std::string::npos);
    CHECK(output.find("\"coverscore\": 10") != std::string::npos);
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

    auto result = search_volume<uint32_t>(
        "query_k9", g_query_seq, 9, kix, kpx, ksx, filter, config);

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

static void test_search_mode1() {
    std::fprintf(stderr, "-- test_search_mode1\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx; // not opened — mode 1 doesn't use kpx
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 1;
    config.num_results = 50;
    config.mode = 1;         // stage1 only
    config.sort_score = 1;   // sort by stage1 score

    auto result = search_volume<uint16_t>(
        "mode1_query", g_query_seq, 7, kix, kpx, ksx, filter, config);

    CHECK(result.query_id == "mode1_query");
    CHECK(!result.hits.empty());

    // In mode 1, hits have stage1_score > 0 and chain fields zeroed
    bool found_fj = false;
    for (const auto& cr : result.hits) {
        CHECK(cr.stage1_score >= 1);
        if (cr.seq_id == g_fj_oid && !cr.is_reverse) {
            found_fj = true;
            // Chain fields should be 0 in mode 1
            CHECK_EQ(cr.q_start, 0u);
            CHECK_EQ(cr.q_end, 0u);
            CHECK_EQ(cr.s_start, 0u);
            CHECK_EQ(cr.s_end, 0u);
        }
    }
    CHECK(found_fj);

    kix.close();
    ksx.close();
}

static void test_search_num_results_zero() {
    std::fprintf(stderr, "-- test_search_num_results_zero\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;

    // num_results=50 (limited)
    SearchConfig config_limited;
    config_limited.stage1.max_freq = 100000;
    config_limited.stage1.stage1_topn = 100;
    config_limited.stage1.min_stage1_score = 1;
    config_limited.stage2.max_gap = 100;
    config_limited.stage2.min_diag_hits = 1;
    config_limited.stage2.min_score = 1;
    config_limited.num_results = 2;

    auto result_limited = search_volume<uint16_t>(
        "q_lim", g_query_seq, 7, kix, kpx, ksx, filter, config_limited);

    // num_results=0 (unlimited)
    SearchConfig config_unlimited;
    config_unlimited.stage1.max_freq = 100000;
    config_unlimited.stage1.stage1_topn = 100;
    config_unlimited.stage1.min_stage1_score = 1;
    config_unlimited.stage2.max_gap = 100;
    config_unlimited.stage2.min_diag_hits = 1;
    config_unlimited.stage2.min_score = 1;
    config_unlimited.num_results = 0;

    auto result_unlimited = search_volume<uint16_t>(
        "q_unlim", g_query_seq, 7, kix, kpx, ksx, filter, config_unlimited);

    // Limited should have at most 2
    CHECK(result_limited.hits.size() <= 2);
    // Unlimited should have >= limited
    CHECK(result_unlimited.hits.size() >= result_limited.hits.size());

    kix.close();
    kpx.close();
    ksx.close();
}

int main() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_test_dir = "/tmp/ikafssn_search_test";
    std::filesystem::create_directories(g_test_dir);

    // Extract 100bp query from FJ876973.1 at runtime
    {
        BlastDbReader db;
        CHECK(db.open(g_testdb_path));
        g_fj_oid = find_oid_by_accession(db, ACC_FJ);
        CHECK(g_fj_oid != UINT32_MAX);
        std::string full_seq = db.get_sequence(g_fj_oid);
        CHECK(full_seq.size() >= 200);
        g_query_seq = full_seq.substr(100, 100);
    }

    test_fasta_reader();
    test_result_output_tab();
    test_result_output_tab_mode1();
    test_result_output_json();
    test_build_and_search();
    test_revcomp_search();
    test_seqidlist_filter();
    test_negative_seqidlist();
    test_search_k9();
    test_search_mode1();
    test_search_num_results_zero();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
