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
#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"
#include "protocol/messages.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "core/spaced_seed.hpp"
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
    Stage1Buffer buf;
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

    // Query: 100bp from FJ876973.1 (extracted at runtime)
    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config, buf);

    CHECK(!result.hits.empty());
    CHECK(result.query_id == "query1");

    // FJ876973.1 OID should be found
    bool found_fj = false;
    for (const auto& cr : result.hits) {
        if (cr.seq_id == g_fj_oid && !cr.is_reverse) {
            found_fj = true;
            CHECK(cr.chainscore >= 2);
        }
    }
    CHECK(found_fj);

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_revcomp_search() {
    Stage1Buffer buf;
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

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

    auto qdata = preprocess_query<uint16_t>(rc_query, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "rc_query", qdata, 7, kix, kpx, ksx, filter, config, buf);

    // Verify the search completes without error
    CHECK(result.query_id == "rc_query");

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_seqidlist_filter() {
    Stage1Buffer buf;
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 1;
    config.nresult = 50;

    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "filtered_query", qdata, 7, kix, kpx, ksx, filter, config, buf);

    // Results should only contain the included OIDs, not FJ876973.1
    for (const auto& cr : result.hits) {
        CHECK(cr.seq_id == oid_gq || cr.seq_id == oid_dq);
    }

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_negative_seqidlist() {
    Stage1Buffer buf;
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 1;
    config.nresult = 50;

    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "neg_query", qdata, 7, kix, kpx, ksx, filter, config, buf);

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
        {"query1", "NM_001234", '+', 10, 450, 1020, 1460, 10, 0, 42, 0},
        {"query1", "XM_005678", '-', 15, 430, 8050, 8465, 8, 0, 38, 2},
    };

    std::ostringstream oss;
    write_results_tsv(oss, hits);
    std::string output = oss.str();

    // Check header (mode 2 default: includes coverscore and chainscore)
    CHECK(output.find("# qseqid") != std::string::npos);
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
        {"query1", "NM_001234", '+', 0, 0, 0, 0, 10, 0, 0, 0},
    };

    std::ostringstream oss;
    write_results_tsv(oss, hits, 1, 1);
    std::string output = oss.str();

    // Mode 1: no qstart/qend/sstart/send/chainscore columns
    CHECK(output.find("coverscore") != std::string::npos);
    CHECK(output.find("chainscore") == std::string::npos);
    CHECK(output.find("qstart") == std::string::npos);
}

static void test_result_output_json() {
    std::fprintf(stderr, "-- test_result_output_json\n");

    std::vector<OutputHit> hits = {
        {"query1", "NM_001234", '+', 10, 450, 1020, 1460, 10, 0, 42, 0},
    };

    std::ostringstream oss;
    write_results_json(oss, hits);
    std::string output = oss.str();

    CHECK(output.find("\"results\"") != std::string::npos);
    CHECK(output.find("\"qseqid\"") != std::string::npos);
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
    Stage1Buffer buf;
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

    auto qdata = preprocess_query<uint32_t>(g_query_seq, 9, nullptr, config);
    auto result = search_volume<uint32_t>(
        "query_k9", qdata, 9, kix, kpx, ksx, filter, config, buf);

    // Should find hits (seq1 has this pattern)
    CHECK(result.query_id == "query_k9");
    // At least some hits expected
    bool found = !result.hits.empty();
    if (found) {
        // Verify hits are reasonable
        for (const auto& cr : result.hits) {
            CHECK(cr.chainscore >= 2);
        }
    }

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_search_mode1() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_search_mode1\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx; // not opened — mode 1 doesn't use kpx
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 1;
    config.nresult = 50;
    config.mode = 1;         // stage1 only
    config.sort_score = 1;   // sort by stage1 score

    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "mode1_query", qdata, 7, kix, kpx, ksx, filter, config, buf);

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

static void test_search_nresult_zero() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_search_nresult_zero\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;

    // nresult=50 (limited)
    SearchConfig config_limited;
    config_limited.stage1.stage1_topn = 100;
    config_limited.stage1.min_stage1_score = 1;
    config_limited.stage2.max_gap = 100;
    config_limited.stage2.min_nhit_diag = 1;
    config_limited.stage2.min_score = 1;
    config_limited.nresult = 2;

    auto qdata_lim = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config_limited);
    auto result_limited = search_volume<uint16_t>(
        "q_lim", qdata_lim, 7, kix, kpx, ksx, filter, config_limited, buf);

    // nresult=0 (unlimited)
    SearchConfig config_unlimited;
    config_unlimited.stage1.stage1_topn = 100;
    config_unlimited.stage1.min_stage1_score = 1;
    config_unlimited.stage2.max_gap = 100;
    config_unlimited.stage2.min_nhit_diag = 1;
    config_unlimited.stage2.min_score = 1;
    config_unlimited.nresult = 0;

    auto qdata_unlim = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config_unlimited);
    auto result_unlimited = search_volume<uint16_t>(
        "q_unlim", qdata_unlim, 7, kix, kpx, ksx, filter, config_unlimited, buf);

    // Limited should have at most 2
    CHECK(result_limited.hits.size() <= 2);
    // Unlimited should have >= limited
    CHECK(result_unlimited.hits.size() >= result_limited.hits.size());

    kix.close();
    kpx.close();
    ksx.close();
}

// regression: build coding + optimal indexes for k=9, t=15, run
// search_volume_both, and verify (1) FJ876973.1 is found via the unified
// Stage 1 path, and (2) the merged Stage 1 score stays in [0, Nqkmer]
// (i.e. no double-counting from the two templates).
static void test_search_both_template() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_search_both_template\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);

    // Build coding-side spaced seed index (k=9, t=15).
    IndexBuilderConfig bcfg_cod;
    bcfg_cod.k = 9;
    bcfg_cod.t = 15;
    bcfg_cod.template_type = static_cast<uint8_t>(TemplateType::kCoding);
    std::string prefix_cod = g_test_dir + "/test.00.09mer.15mer.cod";
    CHECK(build_index<uint32_t>(db, bcfg_cod, prefix_cod, 0, 1, "test", logger));

    // Build optimal-side spaced seed index (k=9, t=15).
    IndexBuilderConfig bcfg_opt;
    bcfg_opt.k = 9;
    bcfg_opt.t = 15;
    bcfg_opt.template_type = static_cast<uint8_t>(TemplateType::kOptimal);
    std::string prefix_opt = g_test_dir + "/test.00.09mer.15mer.opt";
    CHECK(build_index<uint32_t>(db, bcfg_opt, prefix_opt, 0, 1, "test", logger));

    KixReader kix_cod, kix_opt;
    KpxReader kpx_cod, kpx_opt;
    KsxReader ksx;
    CHECK(kix_cod.open(prefix_cod + ".kix"));
    CHECK(kpx_cod.open(prefix_cod + ".kpx"));
    CHECK(kix_opt.open(prefix_opt + ".kix"));
    CHECK(kpx_opt.open(prefix_opt + ".kpx"));
    CHECK(ksx.open(prefix_cod + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 1;
    config.nresult = 50;
    config.t = 15;

    const auto masks_cod = get_seed_masks(9, 15, TemplateType::kCoding);
    const auto masks_opt = get_seed_masks(9, 15, TemplateType::kOptimal);

    auto qdata_cod = preprocess_query<uint32_t>(g_query_seq, 9, nullptr, config,
                                                15, masks_cod);
    auto qdata_opt = preprocess_query<uint32_t>(g_query_seq, 9, nullptr, config,
                                                15, masks_opt);

    auto result = search_volume_both<uint32_t>(
        "both_query", qdata_cod, qdata_opt, 9,
        kix_cod, kpx_cod, kix_opt, kpx_opt, ksx, filter, config, buf);

    CHECK(result.query_id == "both_query");
    CHECK(!result.hits.empty());

    // invariant: per-position dedup keeps Stage 1 score in [0, Nqkmer].
    // Nqkmer = seq_len - t + 1 (window count). A naive cod-then-opt
    // accumulation would inflate scores up to 2*Nqkmer; the q_pos-merge walk
    // in search_one_strand_both keeps the merged score capped at Nqkmer.
    uint32_t Nqkmer = static_cast<uint32_t>(g_query_seq.size()) - 15 + 1;
    bool found_fj = false;
    for (const auto& cr : result.hits) {
        CHECK(cr.stage1_score <= Nqkmer);
        if (cr.seq_id == g_fj_oid) found_fj = true;
    }
    CHECK(found_fj);

    kix_cod.close(); kpx_cod.close();
    kix_opt.close(); kpx_opt.close();
    ksx.close();
}

// v10: -min_query_length integration.  Verify that
//   1. queries shorter than config.min_query_length are skipped with
//      kSkipQueryTooShort and produce a clear skip_detail message;
//   2. queries at or above the threshold preprocess normally.
//
// The integrity check between -min_query_length and the index's
// min_seq_length is performed by ikafssnsearch's open_volumes lambda
// before any preprocessing runs, so we exercise it here by inspecting the
// .kix's min_seq_length field directly (the lambda's failure path
// terminates the binary, which is not unit-test-friendly).
static void test_min_query_length_skip() {
    std::fprintf(stderr, "-- test_min_query_length_skip\n");

    std::string prefix = g_test_dir + "/test.00.07mer";

    // Index built with default min_seq_length=64 in test_build_and_search.
    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
    CHECK_EQ(kix.min_seq_length(), 64u);
    kix.close();

    // Query length 50 < min_query_length 64 -> kSkipQueryTooShort.
    SearchConfig config;
    config.min_query_length = 64;
    std::string short_query(50, 'A');
    auto qdata_short = preprocess_query<uint16_t>(short_query, 7, nullptr, config);
    CHECK_EQ(qdata_short.skip_reason, kSkipQueryTooShort);
    CHECK(!qdata_short.skip_detail.empty());

    // Query length 100 >= min_query_length -> preprocesses normally.
    auto qdata_ok = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    CHECK_EQ(qdata_ok.skip_reason, 0u);
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
    test_search_nresult_zero();
    test_search_both_template();
    test_min_query_length_skip();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
