#include "test_util.hpp"
#include "index/khx_writer.hpp"
#include "index/khx_reader.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.khx";

static void test_basic_roundtrip() {
    std::fprintf(stderr, "-- test_khx_basic_roundtrip\n");

    const int k = 5;
    const uint64_t tbl = table_size(k); // 4^5 = 1024
    Logger logger(Logger::kError);

    // Build counts: some above threshold, most below
    std::vector<uint32_t> counts(tbl, 10);
    counts[0] = 100;    // excluded
    counts[1] = 100;    // excluded
    counts[100] = 200;  // excluded
    counts[1023] = 50;  // excluded

    uint64_t freq_threshold = 20;

    CHECK(write_khx(TEST_FILE, k, counts, freq_threshold, logger));

    // Read back
    KhxReader reader;
    CHECK(reader.open(TEST_FILE));
    CHECK_EQ(reader.k(), k);

    // Verify excluded bits
    CHECK(reader.is_excluded(0));
    CHECK(reader.is_excluded(1));
    CHECK(reader.is_excluded(100));
    CHECK(reader.is_excluded(1023));

    // Verify non-excluded bits
    CHECK(!reader.is_excluded(2));
    CHECK(!reader.is_excluded(50));
    CHECK(!reader.is_excluded(500));

    // Count excluded
    CHECK_EQ(reader.count_excluded(), 4u);

    reader.close();
    std::remove(TEST_FILE);
}

static void test_no_exclusions() {
    std::fprintf(stderr, "-- test_khx_no_exclusions\n");

    const int k = 5;
    const uint64_t tbl = table_size(k);
    Logger logger(Logger::kError);

    std::vector<uint32_t> counts(tbl, 10);
    uint64_t freq_threshold = 100; // no counts exceed this

    CHECK(write_khx(TEST_FILE, k, counts, freq_threshold, logger));

    KhxReader reader;
    CHECK(reader.open(TEST_FILE));
    CHECK_EQ(reader.count_excluded(), 0u);

    for (uint64_t i = 0; i < tbl; i++) {
        CHECK(!reader.is_excluded(i));
    }

    reader.close();
    std::remove(TEST_FILE);
}

static void test_all_excluded() {
    std::fprintf(stderr, "-- test_khx_all_excluded\n");

    const int k = 5;
    const uint64_t tbl = table_size(k);
    Logger logger(Logger::kError);

    std::vector<uint32_t> counts(tbl, 100);
    uint64_t freq_threshold = 0; // all counts exceed this

    CHECK(write_khx(TEST_FILE, k, counts, freq_threshold, logger));

    KhxReader reader;
    CHECK(reader.open(TEST_FILE));
    CHECK_EQ(reader.count_excluded(), tbl);

    for (uint64_t i = 0; i < tbl; i++) {
        CHECK(reader.is_excluded(i));
    }

    reader.close();
    std::remove(TEST_FILE);
}

static void test_larger_k() {
    std::fprintf(stderr, "-- test_khx_larger_k\n");

    const int k = 7;
    const uint64_t tbl = table_size(k); // 4^7 = 16384
    Logger logger(Logger::kError);

    std::vector<uint32_t> counts(tbl, 5);
    // Exclude a specific pattern
    for (uint64_t i = 0; i < tbl; i += 3) {
        counts[i] = 50;
    }
    uint64_t freq_threshold = 10;

    CHECK(write_khx(TEST_FILE, k, counts, freq_threshold, logger));

    KhxReader reader;
    CHECK(reader.open(TEST_FILE));
    CHECK_EQ(reader.k(), k);

    // Verify pattern
    for (uint64_t i = 0; i < tbl; i++) {
        if (i % 3 == 0) {
            CHECK(reader.is_excluded(i));
        } else {
            CHECK(!reader.is_excluded(i));
        }
    }

    reader.close();
    std::remove(TEST_FILE);
}

int main() {
    test_basic_roundtrip();
    test_no_exclusions();
    test_all_excluded();
    test_larger_k();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
