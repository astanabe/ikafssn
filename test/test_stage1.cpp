#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "search/stage1_filter.hpp"
#include "search/oid_filter.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/ksx_reader.hpp"
#include "io/blastdb_reader.hpp"
#include "core/kmer_encoding.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <filesystem>
#include <string>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_index_dir;

// Runtime-extracted query and OID for FJ876973.1
static std::string g_query_seq;
static uint32_t g_fj_oid = UINT32_MAX;

static bool build_test_index() {
    BlastDbReader db;
    if (!db.open(g_testdb_path)) return false;

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 7;

    std::string prefix = g_index_dir + "/test.00.07mer";
    return build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger);
}

static void test_stage1_basic() {
    std::fprintf(stderr, "-- test_stage1_basic\n");

    KixReader kix;
    std::string kix_path = g_index_dir + "/test.00.07mer.kix";
    CHECK(kix.open(kix_path));

    // Query: 100bp from FJ876973.1 (extracted at runtime)
    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter; // no filter
    Stage1Config config;
    config.max_freq = 100000;
    config.stage1_topn = 100;
    config.min_stage1_score = 1;

    auto candidates = stage1_filter(query_kmers, kix, filter, config);

    // Should find FJ876973.1 OID in candidates
    CHECK(!candidates.empty());
    bool found_fj = false;
    for (const auto& c : candidates) {
        if (c.id == g_fj_oid) found_fj = true;
    }
    CHECK(found_fj);

    kix.close();
}

static void test_stage1_max_freq_skip() {
    std::fprintf(stderr, "-- test_stage1_max_freq_skip\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    // Without frequency filter
    OidFilter filter;
    Stage1Config config_unrestricted;
    config_unrestricted.max_freq = 100000;
    config_unrestricted.stage1_topn = 100;
    config_unrestricted.min_stage1_score = 1;
    auto candidates_all = stage1_filter(query_kmers, kix, filter, config_unrestricted);

    // With max_freq=1 (very restrictive) -> fewer candidates
    Stage1Config config_restricted;
    config_restricted.max_freq = 1;
    config_restricted.stage1_topn = 100;
    config_restricted.min_stage1_score = 1;
    auto candidates_filtered = stage1_filter(query_kmers, kix, filter, config_restricted);

    CHECK(candidates_filtered.size() <= candidates_all.size());

    kix.close();
}

static void test_stage1_min_score() {
    std::fprintf(stderr, "-- test_stage1_min_score\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter;
    Stage1Config config;
    config.max_freq = 100000;
    config.stage1_topn = 100;
    config.min_stage1_score = 999999; // Very high threshold

    auto candidates = stage1_filter(query_kmers, kix, filter, config);
    CHECK(candidates.empty());

    kix.close();
}

static void test_stage1_with_oid_filter() {
    std::fprintf(stderr, "-- test_stage1_with_oid_filter\n");

    KixReader kix;
    KsxReader ksx;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));
    CHECK(ksx.open(g_index_dir + "/test.00.07mer.ksx"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    // Filter: include only ACC_GQ and ACC_DQ (exclude FJ876973.1)
    OidFilter filter;
    filter.build({ACC_GQ, ACC_DQ}, ksx, OidFilterMode::kInclude);

    Stage1Config config;
    config.max_freq = 100000;
    config.stage1_topn = 100;
    config.min_stage1_score = 1;

    auto candidates = stage1_filter(query_kmers, kix, filter, config);

    // FJ876973.1 OID should NOT be in results
    for (const auto& c : candidates) {
        CHECK(c.id != g_fj_oid);
    }

    kix.close();
    ksx.close();
}

static void test_stage1_coverscore_vs_matchscore() {
    std::fprintf(stderr, "-- test_stage1_coverscore_vs_matchscore\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter;

    // Coverscore mode (default, type=1)
    Stage1Config config_cover;
    config_cover.max_freq = 100000;
    config_cover.stage1_topn = 100;
    config_cover.min_stage1_score = 1;
    config_cover.stage1_score_type = 1;
    auto cover_candidates = stage1_filter(query_kmers, kix, filter, config_cover);

    // Matchscore mode (type=2)
    Stage1Config config_match;
    config_match.max_freq = 100000;
    config_match.stage1_topn = 100;
    config_match.min_stage1_score = 1;
    config_match.stage1_score_type = 2;
    auto match_candidates = stage1_filter(query_kmers, kix, filter, config_match);

    // Both should find FJ876973.1
    bool cover_found = false, match_found = false;
    uint32_t cover_score = 0, match_score = 0;
    for (const auto& c : cover_candidates) {
        if (c.id == g_fj_oid) { cover_found = true; cover_score = c.score; }
    }
    for (const auto& c : match_candidates) {
        if (c.id == g_fj_oid) { match_found = true; match_score = c.score; }
    }
    CHECK(cover_found);
    CHECK(match_found);

    // Matchscore >= coverscore always (matchscore counts all postings,
    // coverscore counts only distinct query k-mers)
    CHECK(match_score >= cover_score);

    kix.close();
}

static void test_stage1_topn_zero() {
    std::fprintf(stderr, "-- test_stage1_topn_zero\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter;

    // topn=0 (unlimited)
    Stage1Config config_unlimited;
    config_unlimited.max_freq = 100000;
    config_unlimited.stage1_topn = 0;
    config_unlimited.min_stage1_score = 1;
    auto unlimited = stage1_filter(query_kmers, kix, filter, config_unlimited);

    // topn=2 (limited)
    Stage1Config config_limited;
    config_limited.max_freq = 100000;
    config_limited.stage1_topn = 2;
    config_limited.min_stage1_score = 1;
    auto limited = stage1_filter(query_kmers, kix, filter, config_limited);

    // Unlimited should return >= limited count
    CHECK(unlimited.size() >= limited.size());

    // Limited should return at most 2
    CHECK(limited.size() <= 2);

    kix.close();
}

int main() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_index_dir = "/tmp/ikafssn_stage1_test";
    std::filesystem::create_directories(g_index_dir);

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

    CHECK(build_test_index());

    test_stage1_basic();
    test_stage1_max_freq_skip();
    test_stage1_min_score();
    test_stage1_with_oid_filter();
    test_stage1_coverscore_vs_matchscore();
    test_stage1_topn_zero();

    std::filesystem::remove_all(g_index_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
