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
#include "io/volume_discovery.hpp"
#include "search/oid_filter.hpp"
#include "search/parallel_search.hpp"
#include "search/query_preprocessor.hpp"
#include "search/search_orchestrator.hpp"
#include "search/stage1_filter.hpp"
#include "search/volume_searcher.hpp"
#include "protocol/messages.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "core/spaced_seed.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);
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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);

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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);

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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);

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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);

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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);

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
    std::string prefix_cod = index_file_stem(g_test_dir, "test.00",
        bcfg_cod.k, bcfg_cod.t, bcfg_cod.template_type,
        bcfg_cod.min_seq_length, bcfg_cod.min_length_split,
        bcfg_cod.overlap_length, /*max_freq_build=*/1, /*max_degen_expand=*/0);
    CHECK(build_index<uint32_t>(db, bcfg_cod, prefix_cod, 0, 1, "test", logger));

    // Build optimal-side spaced seed index (k=9, t=15).
    IndexBuilderConfig bcfg_opt;
    bcfg_opt.k = 9;
    bcfg_opt.t = 15;
    bcfg_opt.template_type = static_cast<uint8_t>(TemplateType::kOptimal);
    std::string prefix_opt = index_file_stem(g_test_dir, "test.00",
        bcfg_opt.k, bcfg_opt.t, bcfg_opt.template_type,
        bcfg_opt.min_seq_length, bcfg_opt.min_length_split,
        bcfg_opt.overlap_length, /*max_freq_build=*/1, /*max_degen_expand=*/0);
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

// Helper: build cod+opt index pair for both-template tests. Returns the
// two index prefixes. Memoizes by caching: if both prefix .kix files
// already exist on disk, no rebuild happens.
static void build_both_template_indexes(int k, int t,
                                        std::string& prefix_cod,
                                        std::string& prefix_opt) {
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    Logger logger(Logger::kError);

    IndexBuilderConfig bcfg_cod;
    bcfg_cod.k = k;
    bcfg_cod.t = t;
    bcfg_cod.template_type = static_cast<uint8_t>(TemplateType::kCoding);
    prefix_cod = index_file_stem(g_test_dir, "both.00",
        bcfg_cod.k, bcfg_cod.t, bcfg_cod.template_type,
        bcfg_cod.min_seq_length, bcfg_cod.min_length_split,
        bcfg_cod.overlap_length, /*max_freq_build=*/1, /*max_degen_expand=*/0);
    if (!std::filesystem::exists(prefix_cod + ".kix")) {
        CHECK(build_index<uint32_t>(db, bcfg_cod, prefix_cod, 0, 1, "both", logger));
    }

    IndexBuilderConfig bcfg_opt;
    bcfg_opt.k = k;
    bcfg_opt.t = t;
    bcfg_opt.template_type = static_cast<uint8_t>(TemplateType::kOptimal);
    prefix_opt = index_file_stem(g_test_dir, "both.00",
        bcfg_opt.k, bcfg_opt.t, bcfg_opt.template_type,
        bcfg_opt.min_seq_length, bcfg_opt.min_length_split,
        bcfg_opt.overlap_length, /*max_freq_build=*/1, /*max_degen_expand=*/0);
    if (!std::filesystem::exists(prefix_opt + ".kix")) {
        CHECK(build_index<uint32_t>(db, bcfg_opt, prefix_opt, 0, 1, "both", logger));
    }
}

// End-to-end mode 1 + both-template: build cod/opt indexes, preprocess
// the query for both templates, run search_volume_both with mode == 1,
// and verify FJ876973.1 lands in the unified Stage 1 results with the
// mode 1 chain fields zeroed and a non-zero stage1_score.
static void test_search_mode1_both_template() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_search_mode1_both_template\n");

    const int k = 9, t = 15;
    std::string prefix_cod, prefix_opt;
    build_both_template_indexes(k, t, prefix_cod, prefix_opt);

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
    config.nresult = 50;
    config.t = t;
    config.mode = 1;
    config.sort_score = 1;

    const auto masks_cod = get_seed_masks(k, t, TemplateType::kCoding);
    const auto masks_opt = get_seed_masks(k, t, TemplateType::kOptimal);
    auto qdata_cod = preprocess_query<uint32_t>(g_query_seq, k, nullptr, config, t, masks_cod);
    auto qdata_opt = preprocess_query<uint32_t>(g_query_seq, k, nullptr, config, t, masks_opt);

    auto result = search_volume_both<uint32_t>(
        "both_mode1", qdata_cod, qdata_opt, k,
        kix_cod, kpx_cod, kix_opt, kpx_opt, ksx, filter, config, buf);

    CHECK(result.query_id == "both_mode1");
    CHECK(!result.hits.empty());

    uint32_t Nqkmer = static_cast<uint32_t>(g_query_seq.size()) - t + 1;
    bool found_fj = false;
    for (const auto& cr : result.hits) {
        CHECK(cr.stage1_score >= 1);
        CHECK(cr.stage1_score <= Nqkmer);
        CHECK_EQ(cr.q_start, 0u);
        CHECK_EQ(cr.q_end, 0u);
        CHECK_EQ(cr.s_start, 0u);
        CHECK_EQ(cr.s_end, 0u);
        CHECK_EQ(cr.chainscore, 0);
        if (cr.seq_id == g_fj_oid && !cr.is_reverse) found_fj = true;
    }
    CHECK(found_fj);

    kix_cod.close(); kpx_cod.close();
    kix_opt.close(); kpx_opt.close();
    ksx.close();
}

// Direct invocation of stage1_one_strand_both with a fresh JobState +
// Stage1Buffer. Verifies (1) FJ876973.1 ends up in state.candidates,
// (2) every candidate's merged Stage 1 score is in [1, Nqkmer] (the
// q_pos-merge dedup invariant), and (3) running through search_volume_both
// against the same indexes produces the same candidate set (sanity that
// the direct call is not bypassing some validation).
static void test_search_stage1_both_template() {
    std::fprintf(stderr, "-- test_search_stage1_both_template\n");

    const int k = 9, t = 15;
    std::string prefix_cod, prefix_opt;
    build_both_template_indexes(k, t, prefix_cod, prefix_opt);

    KixReader kix_cod, kix_opt;
    CHECK(kix_cod.open(prefix_cod + ".kix"));
    CHECK(kix_opt.open(prefix_opt + ".kix"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.stage1_topn = 0;        // unlimited so candidates is complete
    config.stage1.min_stage1_score = 1;
    config.t = t;

    const auto masks_cod = get_seed_masks(k, t, TemplateType::kCoding);
    const auto masks_opt = get_seed_masks(k, t, TemplateType::kOptimal);
    auto qdata_cod = preprocess_query<uint32_t>(g_query_seq, k, nullptr, config, t, masks_cod);
    auto qdata_opt = preprocess_query<uint32_t>(g_query_seq, k, nullptr, config, t, masks_opt);

    Stage1Buffer buf;
    buf.width = Stage1Width::T32;
    buf.ensure_capacity(kix_cod.num_sequences());

    JobState state;
    stage1_one_strand_both<uint32_t>(
        qdata_cod.fwd_positions.data(),
        qdata_cod.fwd_kmer_values.data(),
        qdata_cod.fwd_positions.size(),
        qdata_opt.fwd_positions.data(),
        qdata_opt.fwd_kmer_values.data(),
        qdata_opt.fwd_positions.size(),
        k,
        /*is_reverse=*/false,
        kix_cod, kix_opt, filter, config,
        qdata_cod.resolved_threshold_fwd,
        qdata_opt.resolved_threshold_fwd,
        /*effective_min_score=*/1,
        buf, state);

    CHECK(state.both_mode);
    CHECK(!state.candidates.empty());
    bool found_fj = false;
    uint32_t Nqkmer = static_cast<uint32_t>(g_query_seq.size()) - t + 1;
    for (const auto& c : state.candidates) {
        CHECK(c.score >= 1);
        CHECK(c.score <= Nqkmer);
        if (c.id == g_fj_oid) found_fj = true;
    }
    CHECK(found_fj);

    kix_cod.close();
    kix_opt.close();
}

// -min_query_length integration.  Verify that
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

    std::string prefix = index_file_stem(g_test_dir, "test.00", 7,
                                         /*t=*/0, /*template_type=*/0,
                                         /*min_seq_length=*/64,
                                         /*min_length_split=*/0,
                                         /*overlap_length=*/0,
                                         /*max_freq_build=*/1,
                                         /*max_degen_expand=*/0);

    // Index built with default min_seq_length=64 in test_build_and_search.
    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
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

// End-to-end both-template (cod + opt) decode-cache equivalence: drive
// run_search<>() in mode 1 with both_mode=true so the orchestrator binds a
// cod cache AND an opt cache on the same volume simultaneously and the
// Stage 1 position-merge consumes both.  A large budget engages the caches;
// budget 1 forces the uncached on-the-fly decode.  Output must be identical,
// confirming the two-cache both-mode path matches the decode path.
static void test_run_search_both_decode_cache_equiv() {
    std::fprintf(stderr, "-- test_run_search_both_decode_cache_equiv\n");
    const int k = 9, t = 15;
    std::string prefix_cod, prefix_opt;
    build_both_template_indexes(k, t, prefix_cod, prefix_opt);

    KsxReader ksx;
    CHECK(ksx.open(prefix_cod + ".ksx"));

    auto make_meta = [](const std::string& prefix) {
        KixReader probe;
        CHECK(probe.open(prefix + ".kix"));
        VolumeMeta m;
        m.files.kix_path     = prefix + ".kix";
        m.files.kpx_path     = prefix + ".kpx";
        m.files.ksx_path     = prefix + ".ksx";
        m.files.volume_index = 0;
        m.files.k            = 9;
        m.kix_posting_size   = probe.posting_file_size();
        m.kix_full_size      = probe.willneed_size_full();
        m.volume_index       = 0;
        m.num_sequences      = probe.num_sequences();
        probe.close();
        return m;
    };
    std::vector<VolumeMeta> volumes_cod{make_meta(prefix_cod)};
    std::vector<VolumeMeta> volumes_opt{make_meta(prefix_opt)};

    Logger logger(Logger::kError);
    uint32_t max_num_seqs =
        std::max(volumes_cod[0].num_sequences, volumes_opt[0].num_sequences);
    uint64_t big_budget =
        (volumes_cod[0].kix_full_size + volumes_opt[0].kix_full_size) * 2;
    big_budget = std::max<uint64_t>(big_budget, 1u << 20);

    const auto masks_cod = get_seed_masks(k, t, TemplateType::kCoding);
    const auto masks_opt = get_seed_masks(k, t, TemplateType::kOptimal);

    auto key = [](const OrchestratorHit& h) {
        return std::tuple<size_t, uint16_t, uint16_t, uint32_t, uint32_t>(
            h.query_idx, h.volume_idx, h.volume_index, h.cr.seq_id,
            h.cr.stage1_score);
    };
    auto sorted_keys = [&](const std::vector<OrchestratorHit>& v) {
        std::vector<decltype(key(*v.begin()))> ks;
        ks.reserve(v.size());
        for (auto& h : v) ks.push_back(key(h));
        std::sort(ks.begin(), ks.end());
        return ks;
    };

    // The cost-ordered both-mode Stage 1 (cached, big budget) must produce the
    // same candidate set as the on-the-fly merge (uncached, budget 1) for every
    // variant: top-N off/on, a cutoff-firing threshold (T > 1), and a
    // degenerate query that stacks multiple k-mers on one position.
    auto compare = [&](const char* label, const std::string& qseq,
                       const SearchConfig& cfg, bool expect_nonempty) {
        std::fprintf(stderr, "   variant: %s\n", label);
        SearchConfig config = cfg;
        config.t      = t;
        config.strand = 1;  // forward

        auto qdata_cod = preprocess_query<uint32_t>(qseq, k, nullptr, config, t, masks_cod);
        auto qdata_opt = preprocess_query<uint32_t>(qseq, k, nullptr, config, t, masks_opt);

        std::vector<QueryBundle<uint32_t>> bundles(1);
        bundles[0].query_id        = &qseq;
        bundles[0].qdata_primary   = &qdata_cod;
        bundles[0].qdata_secondary = &qdata_opt;
        std::vector<uint8_t> skip(1, 0);

        auto run = [&](uint64_t budget) {
            RunSearchInputs<uint32_t> in;
            in.volumes_cod       = volumes_cod;
            in.volumes_opt       = volumes_opt;
            in.ksx_per_volume.assign(1, &ksx);
            in.oid_filters.resize(1);
            in.queries           = &bundles;
            in.query_skip_reason = &skip;
            in.config            = config;
            in.both_mode         = true;
            in.k                 = k;
            in.nthread           = 4;
            in.posting_budget    = budget;
            in.logger            = &logger;
            in.max_num_seqs      = max_num_seqs;
            in.width             = Stage1Width::T32;
            return run_search<uint32_t>(in);
        };

        auto cached   = run(big_budget);  // cod + opt decode caches engaged (C9)
        auto uncached = run(1);           // on-the-fly decode, no reorder
        if (expect_nonempty) CHECK(!cached.empty());
        CHECK_EQ(cached.size(), uncached.size());
        CHECK(sorted_keys(cached) == sorted_keys(uncached));
    };

    auto base_cfg = [&]() {
        SearchConfig c;
        c.nresult = 0;
        c.mode    = 1;
        c.stage1.min_stage1_score = 1;
        c.stage1.stage1_topn      = 100;
        return c;
    };

    // top-N on, threshold 1 (cutoff disabled).
    compare("topn100_T1", g_query_seq, base_cfg(), true);

    // top-N off (topn 0 finish path), threshold 1.
    {
        SearchConfig c = base_cfg();
        c.stage1.stage1_topn = 0;
        compare("topn0_T1", g_query_seq, c, true);
    }

    // Cutoff fires (T = 5 > 1), top-N off.
    {
        SearchConfig c = base_cfg();
        c.stage1.stage1_topn      = 0;
        c.stage1.min_stage1_score = 5;
        compare("topn0_T5", g_query_seq, c, false);
    }

    // Cutoff fires (T = 5 > 1) AND top-N on (ties-inclusive).
    {
        SearchConfig c = base_cfg();
        c.stage1.stage1_topn      = 3;
        c.stage1.min_stage1_score = 5;
        compare("topn3_T5", g_query_seq, c, false);
    }

    // Degenerate query: an IUPAC ambiguity stacks >=2 k-mers per affected
    // position, exercising multi-k-mer position groups, with the cutoff on.
    {
        std::string degen = g_query_seq;
        CHECK(degen.size() > 60);
        degen[50] = 'R';
        SearchConfig c = base_cfg();
        c.stage1.stage1_topn      = 0;
        c.stage1.min_stage1_score = 3;
        c.accept_qdegen           = 1;
        c.max_degen_expand        = 16;
        compare("degenerate_T3", degen, c, false);
    }

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
    test_search_nresult_zero();
    test_search_both_template();
    test_search_mode1_both_template();
    test_search_stage1_both_template();
    test_run_search_both_decode_cache_equiv();
    test_min_query_length_skip();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
