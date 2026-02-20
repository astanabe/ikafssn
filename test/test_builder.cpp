#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/kmer_encoding.hpp"
#include "core/varint.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

using namespace ikafssn;
using namespace ssu_fixture;

// Test data paths
static std::string g_testdb_path;
static std::string g_output_dir;

// Decode ID posting list from raw data
static std::vector<uint32_t> decode_id_postings(
    const uint8_t* data, uint64_t offset, uint32_t count) {
    std::vector<uint32_t> result;
    if (count == 0) return result;
    result.reserve(count);

    const uint8_t* p = data + offset;
    uint32_t prev_id = 0;
    for (uint32_t i = 0; i < count; i++) {
        uint32_t delta;
        p += varint_decode(p, delta);
        uint32_t id = (i == 0) ? delta : prev_id + delta;
        result.push_back(id);
        prev_id = id;
    }
    return result;
}

// Decode pos posting list from raw data (needs seq_ids to detect boundaries)
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, uint64_t offset, uint32_t count,
    const std::vector<uint32_t>& seq_ids) {
    std::vector<uint32_t> result;
    if (count == 0) return result;
    result.reserve(count);

    const uint8_t* p = data + offset;
    uint32_t prev_pos = 0;
    uint32_t prev_id = UINT32_MAX;
    for (uint32_t i = 0; i < count; i++) {
        uint32_t val;
        p += varint_decode(p, val);
        bool new_seq = (seq_ids[i] != prev_id);
        uint32_t pos = new_seq ? val : prev_pos + val;
        result.push_back(pos);
        prev_pos = pos;
        prev_id = seq_ids[i];
    }
    return result;
}

static void test_build_and_verify_ksx() {
    std::fprintf(stderr, "-- test_build_and_verify_ksx\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t num_seqs = db.num_sequences();
    CHECK(num_seqs > 0);

    // Build index with k=7
    Logger logger(Logger::kError); // quiet
    IndexBuilderConfig config;
    config.k = 7;
    config.verbose = false;

    std::string prefix = g_output_dir + "/test.00.07mer";
    CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger));

    // Verify .ksx
    KsxReader ksx;
    CHECK(ksx.open(prefix + ".ksx"));
    CHECK_EQ(ksx.num_sequences(), num_seqs);

    // Verify sequence lengths match BLAST DB
    for (uint32_t i = 0; i < num_seqs; i++) {
        CHECK_EQ(ksx.seq_length(i), db.seq_length(i));
    }

    // Verify accessions are non-empty
    for (uint32_t i = 0; i < num_seqs; i++) {
        auto acc = ksx.accession(i);
        CHECK(!acc.empty());
    }

    ksx.close();
}

static void test_build_and_verify_kix_kpx() {
    std::fprintf(stderr, "-- test_build_and_verify_kix_kpx\n");

    // Index should already be built from previous test
    std::string prefix = g_output_dir + "/test.00.07mer";

    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
    KpxReader kpx;
    CHECK(kpx.open(prefix + ".kpx"));

    // Verify headers
    CHECK_EQ(kix.k(), 7);
    CHECK_EQ(kix.kmer_type(), 0u); // uint16_t for k=7
    CHECK(kix.num_sequences() > 0);
    CHECK_EQ(kpx.k(), 7);
    CHECK_EQ(kix.total_postings(), kpx.total_postings());

    uint64_t tbl_sz = table_size(7); // 4^7 = 16384
    CHECK_EQ(kix.table_size(), tbl_sz);
    CHECK_EQ(kpx.table_size(), tbl_sz);

    // Verify total postings via counts
    uint64_t sum_counts = 0;
    for (uint64_t i = 0; i < tbl_sz; i++) {
        sum_counts += kix.counts()[i];
    }
    CHECK_EQ(sum_counts, kix.total_postings());

    // Verify the posting lists are correctly sorted
    // For each kmer with postings, decode and check seq_ids are non-decreasing
    uint32_t kmers_checked = 0;
    for (uint64_t kmer = 0; kmer < tbl_sz; kmer++) {
        uint32_t cnt = kix.counts()[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids = decode_id_postings(
            kix.posting_data(), kix.offsets()[kmer], cnt);
        CHECK_EQ(ids.size(), static_cast<size_t>(cnt));

        // IDs must be non-decreasing
        for (size_t j = 1; j < ids.size(); j++) {
            CHECK(ids[j] >= ids[j - 1]);
        }

        // All IDs must be valid OIDs
        for (uint32_t id : ids) {
            CHECK(id < kix.num_sequences());
        }

        // Decode positions and verify they are valid
        std::vector<uint32_t> positions = decode_pos_postings(
            kpx.posting_data(), kpx.pos_offsets()[kmer], cnt, ids);
        CHECK_EQ(positions.size(), static_cast<size_t>(cnt));

        kmers_checked++;
    }
    CHECK(kmers_checked > 0);

    kix.close();
    kpx.close();
}

static void test_known_kmer_in_index() {
    std::fprintf(stderr, "-- test_known_kmer_in_index\n");

    // Find FJ876973.1 OID and extract first 7bp as target k-mer
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t target_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(target_oid != UINT32_MAX);

    std::string full_seq = db.get_sequence(target_oid);
    CHECK(full_seq.size() >= 7);
    std::string first7 = full_seq.substr(0, 7);

    KmerScanner<uint16_t> scanner(7);
    uint16_t target_kmer = 0;
    bool found = false;
    scanner.scan(first7.data(), first7.size(), [&](uint32_t /*pos*/, uint16_t kmer) {
        target_kmer = kmer;
        found = true;
    });
    CHECK(found);

    // Now check the index
    std::string prefix = g_output_dir + "/test.00.07mer";
    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));

    uint32_t cnt = kix.counts()[target_kmer];
    CHECK(cnt > 0);

    std::vector<uint32_t> ids = decode_id_postings(
        kix.posting_data(), kix.offsets()[target_kmer], cnt);

    // FJ876973.1 OID should be in the posting list
    bool has_target = false;
    for (uint32_t id : ids) {
        if (id == target_oid) has_target = true;
    }
    CHECK(has_target);

    kix.close();
}

static void test_build_k9_uint32() {
    std::fprintf(stderr, "-- test_build_k9_uint32\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 9;

    std::string prefix = g_output_dir + "/test.00.09mer";
    CHECK(build_index<uint32_t>(db, config, prefix, 0, 1, "test", logger));

    // Verify header
    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
    CHECK_EQ(kix.k(), 9);
    CHECK_EQ(kix.kmer_type(), 1u); // uint32_t for k>=9
    CHECK(kix.num_sequences() > 0);
    CHECK_EQ(kix.table_size(), table_size(9)); // 4^9 = 262144

    // Verify counts sum
    uint64_t sum = 0;
    for (uint64_t i = 0; i < kix.table_size(); i++) {
        sum += kix.counts()[i];
    }
    CHECK_EQ(sum, kix.total_postings());

    kix.close();
}

static void test_build_with_max_freq_build() {
    std::fprintf(stderr, "-- test_build_with_max_freq_build\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);

    // Build without exclusion
    IndexBuilderConfig config1;
    config1.k = 7;
    config1.max_freq_build = 0;

    std::string prefix1 = g_output_dir + "/freq_test1.00.07mer";
    CHECK(build_index<uint16_t>(db, config1, prefix1, 0, 1, "test", logger));

    // Build with max_freq_build=3 (exclude kmers with count > 3)
    IndexBuilderConfig config2;
    config2.k = 7;
    config2.max_freq_build = 3;

    std::string prefix2 = g_output_dir + "/freq_test2.00.07mer";
    CHECK(build_index<uint16_t>(db, config2, prefix2, 0, 1, "test", logger));

    // The filtered index should have fewer or equal total postings
    KixReader kix1, kix2;
    CHECK(kix1.open(prefix1 + ".kix"));
    CHECK(kix2.open(prefix2 + ".kix"));

    CHECK(kix2.total_postings() <= kix1.total_postings());

    // Verify that high-frequency kmers in kix2 have been zeroed
    for (uint64_t kmer = 0; kmer < kix2.table_size(); kmer++) {
        if (kix1.counts()[kmer] > 3) {
            CHECK_EQ(kix2.counts()[kmer], 0u);
        }
    }

    kix1.close();
    kix2.close();
}

static void test_build_with_memory_limits() {
    std::fprintf(stderr, "-- test_build_with_memory_limits\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);

    // Build with large memory_limit (-> 1 partition)
    IndexBuilderConfig config1;
    config1.k = 7;

    std::string prefix1 = g_output_dir + "/part1.00.07mer";
    CHECK(build_index<uint16_t>(db, config1, prefix1, 0, 1, "test", logger));

    // Build with tiny memory_limit to force multiple partitions
    IndexBuilderConfig config4;
    config4.k = 7;
    config4.memory_limit = uint64_t(256) << 10; // 256 KB -> forces multiple partitions

    std::string prefix4 = g_output_dir + "/part4.00.07mer";
    CHECK(build_index<uint16_t>(db, config4, prefix4, 0, 1, "test", logger));

    // Both should produce identical total_postings and counts
    KixReader kix1, kix4;
    CHECK(kix1.open(prefix1 + ".kix"));
    CHECK(kix4.open(prefix4 + ".kix"));

    CHECK_EQ(kix1.total_postings(), kix4.total_postings());

    for (uint64_t kmer = 0; kmer < kix1.table_size(); kmer++) {
        CHECK_EQ(kix1.counts()[kmer], kix4.counts()[kmer]);
    }

    // Verify posting data matches for each kmer
    for (uint64_t kmer = 0; kmer < kix1.table_size(); kmer++) {
        uint32_t cnt = kix1.counts()[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids1 = decode_id_postings(
            kix1.posting_data(), kix1.offsets()[kmer], cnt);
        std::vector<uint32_t> ids4 = decode_id_postings(
            kix4.posting_data(), kix4.offsets()[kmer], cnt);

        CHECK_EQ(ids1.size(), ids4.size());
        for (size_t j = 0; j < ids1.size(); j++) {
            CHECK_EQ(ids1[j], ids4[j]);
        }
    }

    kix1.close();
    kix4.close();
}

static void test_build_parallel_scan() {
    std::fprintf(stderr, "-- test_build_parallel_scan\n");

    // Build with threads=1 and threads=2, verify identical results.
    // This tests that the parallel partition scan produces correct output.
    Logger logger(Logger::kError);

    // Build with 1 thread
    {
        BlastDbReader db;
        CHECK(db.open(g_testdb_path));

        IndexBuilderConfig config;
        config.k = 7;
        config.memory_limit = uint64_t(256) << 10; // 256 KB -> forces multiple partitions
        config.threads = 1;

        std::string prefix = g_output_dir + "/parscan_st.00.07mer";
        CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger));
    }

    // Build with 2 threads
    {
        BlastDbReader db;
        CHECK(db.open(g_testdb_path));

        IndexBuilderConfig config;
        config.k = 7;
        config.memory_limit = uint64_t(256) << 10; // 256 KB -> forces multiple partitions
        config.threads = 2;

        std::string prefix = g_output_dir + "/parscan_mt.00.07mer";
        CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger));
    }

    // Compare kix files
    KixReader kix_st, kix_mt;
    CHECK(kix_st.open(g_output_dir + "/parscan_st.00.07mer.kix"));
    CHECK(kix_mt.open(g_output_dir + "/parscan_mt.00.07mer.kix"));

    CHECK_EQ(kix_st.table_size(), kix_mt.table_size());
    CHECK_EQ(kix_st.total_postings(), kix_mt.total_postings());

    // Compare counts
    bool counts_match = true;
    for (uint64_t i = 0; i < kix_st.table_size(); i++) {
        if (kix_st.counts()[i] != kix_mt.counts()[i]) {
            counts_match = false;
            std::fprintf(stderr, "  counts mismatch at kmer %lu: st=%u mt=%u\n",
                         (unsigned long)i, kix_st.counts()[i], kix_mt.counts()[i]);
            break;
        }
    }
    CHECK(counts_match);

    // Compare posting data for each kmer
    for (uint64_t kmer = 0; kmer < kix_st.table_size(); kmer++) {
        uint32_t cnt = kix_st.counts()[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids_st = decode_id_postings(
            kix_st.posting_data(), kix_st.offsets()[kmer], cnt);
        std::vector<uint32_t> ids_mt = decode_id_postings(
            kix_mt.posting_data(), kix_mt.offsets()[kmer], cnt);

        CHECK_EQ(ids_st.size(), ids_mt.size());
        for (size_t j = 0; j < ids_st.size(); j++) {
            if (ids_st[j] != ids_mt[j]) {
                std::fprintf(stderr, "  posting mismatch at kmer %lu, idx %zu\n",
                             (unsigned long)kmer, j);
                CHECK_EQ(ids_st[j], ids_mt[j]);
                break;
            }
        }
    }

    // Compare kpx files
    KpxReader kpx_st, kpx_mt;
    CHECK(kpx_st.open(g_output_dir + "/parscan_st.00.07mer.kpx"));
    CHECK(kpx_mt.open(g_output_dir + "/parscan_mt.00.07mer.kpx"));

    CHECK_EQ(kpx_st.total_postings(), kpx_mt.total_postings());

    for (uint64_t kmer = 0; kmer < kix_st.table_size(); kmer++) {
        uint32_t cnt = kix_st.counts()[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids = decode_id_postings(
            kix_st.posting_data(), kix_st.offsets()[kmer], cnt);
        std::vector<uint32_t> pos_st = decode_pos_postings(
            kpx_st.posting_data(), kpx_st.pos_offsets()[kmer], cnt, ids);
        std::vector<uint32_t> pos_mt = decode_pos_postings(
            kpx_mt.posting_data(), kpx_mt.pos_offsets()[kmer], cnt, ids);

        CHECK_EQ(pos_st.size(), pos_mt.size());
        for (size_t j = 0; j < pos_st.size(); j++) {
            if (pos_st[j] != pos_mt[j]) {
                std::fprintf(stderr, "  pos mismatch at kmer %lu, idx %zu\n",
                             (unsigned long)kmer, j);
                CHECK_EQ(pos_st[j], pos_mt[j]);
                break;
            }
        }
    }

    kix_st.close(); kix_mt.close();
    kpx_st.close(); kpx_mt.close();
}

int main(int argc, char* argv[]) {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_output_dir = "/tmp/ikafssn_builder_test";

    // Create output directory
    std::filesystem::create_directories(g_output_dir);

    test_build_and_verify_ksx();
    test_build_and_verify_kix_kpx();
    test_known_kmer_in_index();
    test_build_k9_uint32();
    test_build_with_max_freq_build();
    test_build_with_memory_limits();
    test_build_parallel_scan();

    // Clean up
    std::filesystem::remove_all(g_output_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
