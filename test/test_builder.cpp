#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "index/index_builder.hpp"
#include "index/index_filter.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/kmer_encoding.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/posting_decoder.hpp"
#include "util/logger.hpp"
#include "util/simd_dispatch.hpp"

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

// Decode ID posting list via the v4 SeqIdDecoder.
static std::vector<uint32_t> decode_id_postings(
    const uint8_t* data, uint64_t offset, uint64_t byte_len) {
    std::vector<uint32_t> result;
    if (byte_len == 0) return result;
    SeqIdDecoder dec(data + offset, data + offset + byte_len);
    while (dec.has_more()) result.push_back(dec.next());
    return result;
}

// Decode pos posting list via PosDecoder.  Candidate-set-driven: pass
// a sorted distinct seq_id list (which equals the .kix decoded distinct
// seq_id array for this k-mer) as both the kix_decoded view and the
// candidate set, then concatenate the per-candidate position vectors
// in candidate order.
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, uint64_t offset, uint64_t byte_len,
    const std::vector<uint32_t>& distinct_sids) {
    std::vector<uint32_t> result;
    if (byte_len == 0 || distinct_sids.empty()) return result;
    pfd::PosDecodeScratch scratch;
    PosDecoder dec(data + offset, data + offset + byte_len,
                   distinct_sids.data(), distinct_sids.size(),
                   distinct_sids.data(), distinct_sids.size(),
                   &scratch);
    for (size_t i = 0; i < distinct_sids.size(); i++) {
        const auto& v = dec.positions_for(i);
        result.insert(result.end(), v.begin(), v.end());
    }
    return result;
}

// Compute the .kpx per-kmer position count by summing partition
// occ_counts + short1_count + sum(short2 u8 occ_count[]).  The body
// starts directly at the 2-bit kind map and partition / short1 /
// short2 counts are derived via popcount of the kind map.
static uint32_t kpx_position_count(const KixReader& kix,
                                   const KpxReader& kpx,
                                   uint32_t kmer) {
    const uint8_t* p = kpx.posting_file() + kpx.pos_offset(kmer);
    const uint64_t kix_byte_len = kix.posting_list_byte_length(kmer);
    const uint32_t distinct_count = pfd::posting_count(
        kix.posting_file() + kix.posting_list_offset(kmer), kix_byte_len);
    if (distinct_count == 0) return 0;

    // Derive the per-kind counts from the kind map (mirrors the
    // pfd_codec_tier popcount_kinds helper but kept inline for the test).
    const std::size_t kind_map_bytes = (std::size_t(distinct_count) * 2 + 7) / 8;
    uint32_t partition_count = 0, short1_count = 0, short2_count = 0;
    for (uint32_t k = 0; k < distinct_count; k++) {
        uint8_t kind = static_cast<uint8_t>(
            (p[k >> 2] >> ((k & 3) * 2)) & 0x03);
        switch (kind) {
            case 0: short1_count++; break;
            case 1: short2_count++; break;
            case 2: partition_count++; break;
        }
    }

    const uint8_t* gp = p + kind_map_bytes;
    uint32_t partition_pos = 0;
    for (uint32_t g = 0; g < partition_count; g++) {
        uint32_t gcnt;
        std::memcpy(&gcnt, gp, sizeof(uint32_t));
        gp += sizeof(uint32_t);
        partition_pos += gcnt;
        // Skip the FOR-block stream: full blocks (8 + 16*b each) + tail.
        const uint32_t num_blocks = gcnt / 128u;
        const uint32_t tail_count = gcnt % 128u;
        for (uint32_t b = 0; b < num_blocks; b++) {
            uint8_t bw = gp[4];
            gp += 8 + std::size_t(16) * bw;
        }
        const uint8_t got_tail = *gp++;
        (void)got_tail;
        if (tail_count > 0) {
            gp += sizeof(uint32_t);          // tail_min
            uint8_t tail_b = *gp++;
            const std::size_t body_bits = std::size_t(tail_count) * tail_b;
            gp += (body_bits + 7) / 8;
        }
    }

    // Skip short1 FOR stream — its position count is just short1_count
    // (one position per cluster).  No need to walk it for this aggregate.
    // Then read short2's u8 occ_count[] and horizontal-sum.
    if (short1_count > 0) {
        const uint32_t num_blocks = short1_count / 128u;
        const uint32_t tail_count = short1_count % 128u;
        for (uint32_t b = 0; b < num_blocks; b++) {
            uint8_t bw = gp[4];
            gp += 8 + std::size_t(16) * bw;
        }
        const uint8_t got_tail = *gp++;
        (void)got_tail;
        if (tail_count > 0) {
            gp += sizeof(uint32_t);
            uint8_t tail_b = *gp++;
            const std::size_t body_bits = std::size_t(tail_count) * tail_b;
            gp += (body_bits + 7) / 8;
        }
    }

    uint32_t short2_pos = 0;
    for (uint32_t i = 0; i < short2_count; i++) short2_pos += gp[i];

    return partition_pos + short1_count + short2_pos;
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
    // v7: distinct count <= position count (intra-sequence k-mer
    // duplicates are removed from the .kix stream but kept in .kpx).
    CHECK(kix.total_distinct_postings() <= kpx.total_position_count());

    uint64_t tbl_sz = table_size(7); // 4^7 = 16384
    CHECK_EQ(kix.table_size(), tbl_sz);
    CHECK_EQ(kpx.table_size(), tbl_sz);

    // Verify total distinct postings via counts
    auto counts = kix.bulk_count_postings();
    uint64_t sum_counts = 0;
    for (uint32_t i = 0; i < tbl_sz; i++) {
        sum_counts += counts[i];
    }
    CHECK_EQ(sum_counts, kix.total_distinct_postings());

    // Verify the posting lists are correctly sorted
    // For each kmer with postings, decode and check seq_ids are non-decreasing
    uint32_t kmers_checked = 0;
    for (uint32_t kmer = 0; kmer < tbl_sz; kmer++) {
        uint32_t cnt = counts[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids = decode_id_postings(
            kix.posting_file(), kix.posting_list_offset(kmer),
            kix.posting_list_byte_length(kmer));
        CHECK_EQ(ids.size(), static_cast<size_t>(cnt));

        // IDs must be non-decreasing
        for (size_t j = 1; j < ids.size(); j++) {
            CHECK(ids[j] >= ids[j - 1]);
        }

        // All IDs must be valid OIDs
        for (uint32_t id : ids) {
            CHECK(id < kix.num_sequences());
        }

        // Decode positions via v7 candidate-set API: ids is already a
        // sorted distinct seq_id list (decoded from .kix v7).  Total
        // positions returned equal the .kpx per-kmer position_count.
        std::vector<uint32_t> positions = decode_pos_postings(
            kpx.posting_file(), kpx.pos_offset(kmer),
            kpx.posting_file_size() - kpx.pos_offset(kmer),
            ids);
        CHECK_EQ(positions.size(),
                 static_cast<size_t>(kpx_position_count(kix, kpx, kmer)));

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

    auto cnts_target = kix.bulk_count_postings();
    uint32_t cnt = cnts_target[target_kmer];
    CHECK(cnt > 0);

    std::vector<uint32_t> ids = decode_id_postings(
        kix.posting_file(), kix.posting_list_offset(target_kmer),
        kix.posting_list_byte_length(target_kmer));

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
    auto counts_k9 = kix.bulk_count_postings();
    uint64_t sum = 0;
    for (uint32_t i = 0; i < kix.table_size(); i++) {
        sum += counts_k9[i];
    }
    CHECK_EQ(sum, kix.total_distinct_postings());

    kix.close();
}

static void test_build_with_max_freq_build() {
    std::fprintf(stderr, "-- test_build_with_max_freq_build\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);

    // Build without exclusion (reference)
    IndexBuilderConfig config1;
    config1.k = 7;

    std::string prefix1 = g_output_dir + "/freq_test1.00.07mer";
    CHECK(build_index<uint16_t>(db, config1, prefix1, 0, 1, "test", logger));

    // Build with keep_tmp=true, then apply cross-volume filter with threshold=3
    IndexBuilderConfig config2;
    config2.k = 7;
    config2.keep_tmp = true;

    std::string prefix2 = g_output_dir + "/freq_test2.00.07mer";
    CHECK(build_index<uint16_t>(db, config2, prefix2, 0, 1, "test", logger));

    // Apply cross-volume filtering (single volume, threshold=3)
    std::string khx_path2 = g_output_dir + "/freq_test2.07mer.khx";
    CHECK(filter_volumes_cross_volume({prefix2}, {prefix2}, khx_path2, 7, 3, 1, logger));

    // The filtered index should have fewer or equal total postings
    KixReader kix1, kix2;
    CHECK(kix1.open(prefix1 + ".kix"));
    CHECK(kix2.open(prefix2 + ".kix"));

    CHECK(kix2.total_distinct_postings() <= kix1.total_distinct_postings());

    // Verify that high-frequency kmers in kix2 have been zeroed
    auto counts1 = kix1.bulk_count_postings();
    auto counts2 = kix2.bulk_count_postings();
    for (uint32_t kmer = 0; kmer < kix2.table_size(); kmer++) {
        if (counts1[kmer] > 3) {
            CHECK_EQ(counts2[kmer], 0u);
        }
    }

    kix1.close();
    kix2.close();

    // Build with keep_tmp=true, then apply cross-volume filter with fractional threshold
    uint32_t nseqs = db.num_sequences();
    uint64_t frac_threshold = static_cast<uint64_t>(0.5 * nseqs);
    if (frac_threshold == 0) frac_threshold = 1;

    IndexBuilderConfig config3;
    config3.k = 7;
    config3.keep_tmp = true;

    std::string prefix3 = g_output_dir + "/freq_test3.00.07mer";
    CHECK(build_index<uint16_t>(db, config3, prefix3, 0, 1, "test", logger));

    std::string khx_path3 = g_output_dir + "/freq_test3.07mer.khx";
    CHECK(filter_volumes_cross_volume({prefix3}, {prefix3}, khx_path3, 7, frac_threshold, 1, logger));

    // Re-open kix1 for comparison with fractional-threshold build
    CHECK(kix1.open(prefix1 + ".kix"));

    KixReader kix3;
    CHECK(kix3.open(prefix3 + ".kix"));
    CHECK(kix3.total_distinct_postings() <= kix1.total_distinct_postings());

    auto counts1b = kix1.bulk_count_postings();
    auto counts3 = kix3.bulk_count_postings();
    for (uint32_t kmer = 0; kmer < kix3.table_size(); kmer++) {
        if (counts1b[kmer] > frac_threshold) {
            CHECK_EQ(counts3[kmer], 0u);
        }
    }
    kix1.close();
    kix3.close();
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

    CHECK_EQ(kix1.total_distinct_postings(), kix4.total_distinct_postings());

    auto counts_m1 = kix1.bulk_count_postings();
    auto counts_m4 = kix4.bulk_count_postings();
    for (uint32_t kmer = 0; kmer < kix1.table_size(); kmer++) {
        CHECK_EQ(counts_m1[kmer], counts_m4[kmer]);
    }

    // Verify posting list bytes match for each kmer
    for (uint32_t kmer = 0; kmer < kix1.table_size(); kmer++) {
        uint32_t cnt = counts_m1[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids1 = decode_id_postings(
            kix1.posting_file(), kix1.posting_list_offset(kmer),
            kix1.posting_list_byte_length(kmer));
        std::vector<uint32_t> ids4 = decode_id_postings(
            kix4.posting_file(), kix4.posting_list_offset(kmer),
            kix4.posting_list_byte_length(kmer));

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
    CHECK_EQ(kix_st.total_distinct_postings(), kix_mt.total_distinct_postings());

    // Compare counts
    auto counts_st = kix_st.bulk_count_postings();
    auto counts_mt = kix_mt.bulk_count_postings();
    bool counts_match = true;
    for (uint32_t i = 0; i < kix_st.table_size(); i++) {
        if (counts_st[i] != counts_mt[i]) {
            counts_match = false;
            std::fprintf(stderr, "  counts mismatch at kmer %lu: st=%u mt=%u\n",
                         (unsigned long)i, counts_st[i], counts_mt[i]);
            break;
        }
    }
    CHECK(counts_match);

    // Compare posting list bytes for each kmer
    for (uint32_t kmer = 0; kmer < kix_st.table_size(); kmer++) {
        uint32_t cnt = counts_st[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids_st = decode_id_postings(
            kix_st.posting_file(), kix_st.posting_list_offset(kmer),
            kix_st.posting_list_byte_length(kmer));
        std::vector<uint32_t> ids_mt = decode_id_postings(
            kix_mt.posting_file(), kix_mt.posting_list_offset(kmer),
            kix_mt.posting_list_byte_length(kmer));

        CHECK_EQ(ids_st.size(), ids_mt.size());
        for (size_t j = 0; j < ids_st.size(); j++) {
            if (ids_st[j] != ids_mt[j]) {
                std::fprintf(stderr, "  posting list mismatch at kmer %lu, idx %zu\n",
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

    CHECK_EQ(kpx_st.total_position_count(), kpx_mt.total_position_count());

    for (uint32_t kmer = 0; kmer < kix_st.table_size(); kmer++) {
        uint32_t cnt = counts_st[kmer];
        if (cnt == 0) continue;

        std::vector<uint32_t> ids = decode_id_postings(
            kix_st.posting_file(), kix_st.posting_list_offset(kmer),
            kix_st.posting_list_byte_length(kmer));
        std::vector<uint32_t> pos_st = decode_pos_postings(
            kpx_st.posting_file(), kpx_st.pos_offset(kmer),
            kpx_st.posting_file_size() - kpx_st.pos_offset(kmer),
            ids);
        std::vector<uint32_t> pos_mt = decode_pos_postings(
            kpx_mt.posting_file(), kpx_mt.pos_offset(kmer),
            kpx_mt.posting_file_size() - kpx_mt.pos_offset(kmer),
            ids);

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

// v10: -min_seq_length filter.  Build the same DB twice with different
// thresholds and verify that the higher threshold produces a strictly
// smaller / equal index (more sequences excluded).
static void test_min_seq_length_filter() {
    std::fprintf(stderr, "-- test_min_seq_length_filter\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    const uint32_t db_nseq = db.num_sequences();

    // Build with default min_seq_length=64 first.
    Logger logger(Logger::kError);
    IndexBuilderConfig cfg64;
    cfg64.k = 7;
    cfg64.min_seq_length = 64;
    std::string prefix64 = g_output_dir + "/test.minlen64.07mer";
    CHECK(build_index<uint16_t>(db, cfg64, prefix64, 0, 1, "test", logger));

    KsxReader ksx64;
    CHECK(ksx64.open(prefix64 + ".ksx"));
    const uint32_t kept64 = ksx64.num_sequences();
    // 1 fragment per parent, so num_parents == num_sequences.
    CHECK_EQ(ksx64.num_parents(), kept64);
    // Every fragment spans the whole parent.
    for (uint32_t i = 0; i < kept64; i++) {
        CHECK_EQ(ksx64.fragment_start(i), 1u);
        CHECK_EQ(ksx64.fragment_end(i), ksx64.parent_length(i));
        CHECK_EQ(ksx64.parent_index(i), i);
    }
    ksx64.close();

    // .kix mirrors the fragment count from .ksx.
    KixReader kix64;
    CHECK(kix64.open(prefix64 + ".kix"));
    CHECK_EQ(kix64.num_sequences(), kept64);
    kix64.close();
    KpxReader kpx64;
    CHECK(kpx64.open(prefix64 + ".kpx"));
    kpx64.close();

    // Build again with a very large min_seq_length so most / all sequences
    // are filtered out.  10 Mbp exceeds every SSU rRNA gene length.
    IndexBuilderConfig cfg_big;
    cfg_big.k = 7;
    cfg_big.min_seq_length = 10'000'000u;
    std::string prefix_big = g_output_dir + "/test.minlen10M.07mer";
    CHECK(build_index<uint16_t>(db, cfg_big, prefix_big, 0, 1, "test", logger));

    KsxReader ksx_big;
    CHECK(ksx_big.open(prefix_big + ".ksx"));
    const uint32_t kept_big = ksx_big.num_sequences();
    CHECK(kept_big <= kept64);
    CHECK(kept_big < db_nseq);
    ksx_big.close();
}

// v10: degenerate fragment layout — 1 parent = 1 fragment, parent_length
// equals fragment length, blast_oid is monotonically increasing across
// surviving parents (i.e. the original BLAST DB OID after the
// min_seq_length filter).
static void test_degenerate_fragment_layout() {
    std::fprintf(stderr, "-- test_degenerate_fragment_layout\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig cfg;
    cfg.k = 7;
    cfg.min_seq_length = 64;
    std::string prefix = g_output_dir + "/test.degenfrag.07mer";
    CHECK(build_index<uint16_t>(db, cfg, prefix, 0, 1, "test", logger));

    KsxReader ksx;
    CHECK(ksx.open(prefix + ".ksx"));

    CHECK_EQ(ksx.num_parents(), ksx.num_sequences());
    uint32_t prev_blast_oid = 0;
    bool first = true;
    for (uint32_t pidx = 0; pidx < ksx.num_parents(); pidx++) {
        uint32_t blast_oid = ksx.blast_oid(pidx);
        if (!first) CHECK(blast_oid > prev_blast_oid);
        prev_blast_oid = blast_oid;
        first = false;
        // Parent length matches BLAST DB.
        CHECK_EQ(ksx.parent_length(pidx), db.seq_length(blast_oid));
        // Accession matches.
        CHECK(ksx.parent_accession(pidx) == db.get_accession(blast_oid));
    }
    ksx.close();
}

int main(int argc, char* argv[]) {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_output_dir = test_tmpdir("/tmp/ikafssn_builder_test");

    // Create output directory
    std::filesystem::create_directories(g_output_dir);

    test_build_and_verify_ksx();
    test_build_and_verify_kix_kpx();
    test_known_kmer_in_index();
    test_build_k9_uint32();
    test_build_with_max_freq_build();
    test_build_with_memory_limits();
    test_build_parallel_scan();
    test_min_seq_length_filter();
    test_degenerate_fragment_layout();

    // Clean up
    std::filesystem::remove_all(g_output_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
