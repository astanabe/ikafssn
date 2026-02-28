#include "test_util.hpp"
#include "ssu_test_fixture.hpp"

#include "search/stage3_alignment.hpp"
#include "search/volume_searcher.hpp"
#include "search/query_preprocessor.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;
static std::string g_query_seq;

static void test_stage3_pipeline() {
    std::fprintf(stderr, "-- test_stage3_pipeline (build -> search -> align)\n");

    // Build index with k=7
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    // Extract query: 100bp from FJ876973.1
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    CHECK(fj_seq.size() >= 200);
    g_query_seq = fj_seq.substr(50, 100);
    db.close();

    Logger logger(Logger::kDebug);
    IndexBuilderConfig bconfig;
    bconfig.k = 7;

    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/s3test.00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    // Open index
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    // Stage 1+2 search
    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 50;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 10;
    config.mode = 2;

    std::vector<const KixReader*> all_kix = {&kix};
    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, all_kix, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config);

    CHECK(!result.hits.empty());

    // Convert to OutputHit
    std::vector<OutputHit> all_hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        all_hits.push_back(oh);
    }

    CHECK(!all_hits.empty());

    // Prepare query records for Stage 3
    std::vector<FastaRecord> queries;
    queries.push_back({"query1", g_query_seq});

    // Run Stage 3 with traceback
    Stage3Config s3config;
    s3config.traceback = true;
    s3config.gapopen = 10;
    s3config.gapext = 1;
    s3config.fetch_threads = 1;

    auto filtered = run_stage3(all_hits, queries, g_testdb_path, s3config,
                               false, 0.0, 0, logger);

    CHECK(!filtered.empty());

    // Verify alignment results
    for (const auto& h : filtered) {
        // alnscore should be positive for real matches
        CHECK(h.alnscore > 0);
        // CIGAR should be non-empty
        CHECK(!h.cigar.empty());
        // pident should be reasonable (at least some identity)
        CHECK(h.pident > 0.0);
        // nident should be positive
        CHECK(h.nident > 0);
        // Coordinates should be within sequence bounds
        CHECK(h.slen > 0);
        CHECK(h.sstart < h.slen);
        CHECK(h.send < h.slen);
    }

    // The exact match hit (FJ876973.1) should have high pident
    bool found_high_pident = false;
    for (const auto& h : filtered) {
        if (h.sseqid == ACC_FJ && h.sstrand == '+') {
            CHECK(h.pident > 90.0);
            found_high_pident = true;
        }
    }
    CHECK(found_high_pident);
}

static void test_stage3_score_only() {
    std::fprintf(stderr, "-- test_stage3_score_only (traceback=0)\n");

    Logger logger(Logger::kError);

    // Build a quick index
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    std::string query = fj_seq.substr(50, 100);
    db.close();

    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/s3test2.00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 50;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 5;

    std::vector<const KixReader*> all_kix = {&kix};
    auto qdata = preprocess_query<uint16_t>(query, 7, all_kix, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config);
    CHECK(!result.hits.empty());

    std::vector<OutputHit> all_hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        all_hits.push_back(oh);
    }

    std::vector<FastaRecord> queries;
    queries.push_back({"query1", query});

    // Score-only mode (no traceback)
    Stage3Config s3config;
    s3config.traceback = false;
    s3config.fetch_threads = 1;

    auto filtered = run_stage3(all_hits, queries, g_testdb_path, s3config,
                               false, 0.0, 0, logger);
    CHECK(!filtered.empty());

    for (const auto& h : filtered) {
        CHECK(h.alnscore > 0);
        // CIGAR should be empty (no traceback)
        CHECK(h.cigar.empty());
        // pident should be 0 (not computed without traceback)
        CHECK(h.pident == 0.0);
    }
}

static void test_stage3_context() {
    std::fprintf(stderr, "-- test_stage3_context\n");

    Logger logger(Logger::kError);

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    std::string query = fj_seq.substr(100, 50);
    db.close();

    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/s3test3.00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 50;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 5;

    std::vector<const KixReader*> all_kix = {&kix};
    auto qdata = preprocess_query<uint16_t>(query, 7, all_kix, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config);

    if (result.hits.empty()) {
        std::fprintf(stderr, "  (no hits, skipping context test)\n");
        return;
    }

    std::vector<OutputHit> all_hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        all_hits.push_back(oh);
    }

    std::vector<FastaRecord> queries;
    queries.push_back({"query1", query});

    // With integer context
    Stage3Config s3config;
    s3config.traceback = true;
    s3config.fetch_threads = 1;

    auto filtered = run_stage3(all_hits, queries, g_testdb_path, s3config,
                               false, 0.0, 50, logger);  // 50bp context
    CHECK(!filtered.empty());

    for (const auto& h : filtered) {
        CHECK(h.alnscore > 0);
        CHECK(!h.cigar.empty());
    }
}

int main() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_test_dir = "/tmp/ikafssn_test_stage3";
    std::filesystem::create_directories(g_test_dir);

    test_stage3_pipeline();
    test_stage3_score_only();
    test_stage3_context();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
