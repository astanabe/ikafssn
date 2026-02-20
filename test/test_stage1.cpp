#include "test_util.hpp"
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

static std::string g_testdb_path;
static std::string g_index_dir;

static bool build_test_index() {
    BlastDbReader db;
    if (!db.open(g_testdb_path)) return false;

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 7;
    config.partitions = 1;
    config.buffer_size = uint64_t(1) << 30;

    std::string prefix = g_index_dir + "/test.00.07mer";
    return build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger);
}

static void test_stage1_basic() {
    std::fprintf(stderr, "-- test_stage1_basic\n");

    KixReader kix;
    std::string kix_path = g_index_dir + "/test.00.07mer.kix";
    CHECK(kix.open(kix_path));

    // Query: "ACGTACGTACGTACGT" (from seq1)
    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    std::string query = "ACGTACGTACGTACGT";
    scanner.scan(query.data(), query.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter; // no filter
    Stage1Config config;
    config.max_freq = 100000;
    config.stage1_topn = 10;
    config.min_stage1_score = 1;

    auto candidates = stage1_filter(query_kmers, kix, filter, config);

    // Should find at least seq1 (OID 0)
    CHECK(!candidates.empty());
    bool found_seq1 = false;
    for (SeqId id : candidates) {
        if (id == 0) found_seq1 = true;
    }
    CHECK(found_seq1);

    kix.close();
}

static void test_stage1_max_freq_skip() {
    std::fprintf(stderr, "-- test_stage1_max_freq_skip\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    std::string query = "ACGTACGTACGTACGT";
    scanner.scan(query.data(), query.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter;
    Stage1Config config;
    config.max_freq = 1; // Very restrictive
    config.stage1_topn = 10;
    config.min_stage1_score = 1;

    auto candidates = stage1_filter(query_kmers, kix, filter, config);
    // With max_freq=1, many k-mers skipped
    CHECK(candidates.size() <= 5);

    kix.close();
}

static void test_stage1_min_score() {
    std::fprintf(stderr, "-- test_stage1_min_score\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<std::pair<uint32_t, uint16_t>> query_kmers;
    KmerScanner<uint16_t> scanner(7);
    std::string query = "ACGTACGTACGTACGT";
    scanner.scan(query.data(), query.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    OidFilter filter;
    Stage1Config config;
    config.max_freq = 100000;
    config.stage1_topn = 10;
    config.min_stage1_score = 999; // Very high threshold

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
    std::string query = "ACGTACGTACGTACGT";
    scanner.scan(query.data(), query.size(), [&](uint32_t pos, uint16_t kmer) {
        query_kmers.emplace_back(pos, kmer);
    });

    // Filter: only allow seq3 and seq4
    OidFilter filter;
    filter.build({"seq3", "seq4"}, ksx, OidFilterMode::kInclude);

    Stage1Config config;
    config.max_freq = 100000;
    config.stage1_topn = 10;
    config.min_stage1_score = 1;

    auto candidates = stage1_filter(query_kmers, kix, filter, config);

    // seq1 (OID 0) should NOT be in results
    for (SeqId id : candidates) {
        CHECK(id != 0); // seq1 excluded
    }

    kix.close();
    ksx.close();
}

int main() {
    g_testdb_path = std::string(SOURCE_DIR) + "/test/testdata/testdb";
    g_index_dir = "/tmp/ikafssn_stage1_test";
    std::filesystem::create_directories(g_index_dir);

    CHECK(build_test_index());

    test_stage1_basic();
    test_stage1_max_freq_skip();
    test_stage1_min_score();
    test_stage1_with_oid_filter();

    std::filesystem::remove_all(g_index_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
