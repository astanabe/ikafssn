#include "test_util.hpp"
#include "kpx_writer.hpp"
#include "index/kpx_reader.hpp"
#include "core/config.hpp"
#include "search/posting_decoder.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.kpx";

// Decode position posting lists from a v8 .kpx posting list via the
// candidate-set API.  Builds the kix_decoded array as the sorted union
// of `candidates` (treating it as the full distinct seq_id list for the
// k-mer — for these tests the writer adds postings only for those sids,
// so the .kix decoded array equals the candidate set).  Returns a flat
// vector of positions ordered by the candidate seq_id list (positions
// for candidates[0] first, then [1], etc.).
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, const uint8_t* end,
    const std::vector<uint32_t>& candidates) {
    std::vector<uint32_t> distinct_sorted = candidates;
    std::sort(distinct_sorted.begin(), distinct_sorted.end());
    distinct_sorted.erase(std::unique(distinct_sorted.begin(),
                                       distinct_sorted.end()),
                          distinct_sorted.end());

    pfd::PosDecodeScratch scratch;
    PosDecoder decoder(data, end,
                       distinct_sorted.data(), distinct_sorted.size(),
                       distinct_sorted.data(), distinct_sorted.size(),
                       &scratch);
    std::vector<uint32_t> result;
    for (size_t i = 0; i < distinct_sorted.size(); i++) {
        const auto& v = decoder.positions_for(i);
        result.insert(result.end(), v.begin(), v.end());
    }
    return result;
}

static void test_single_seq() {
    int k = 7;
    uint32_t ts = table_size(k);

    // Single sequence, multiple positions for one k-mer (occ>=2 path).
    std::vector<KpxWriter::PostingEntry> entries = {
        {5, 10}, {5, 20}, {5, 30}, {5, 100}
    };
    std::vector<uint32_t> candidates = {5};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 42) {
                writer.add_posting_list(i, entries);
            } else {
                writer.add_posting_list(i, {});
            }
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.k(), k);
        CHECK_EQ(reader.total_position_count(), 4u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(42),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), 4u);
        CHECK_EQ(positions[0], 10u);
        CHECK_EQ(positions[1], 20u);
        CHECK_EQ(positions[2], 30u);
        CHECK_EQ(positions[3], 100u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

// All occ=1 clusters → entire posting list goes through the short_occ1
// sub-bucket (no u8 occ_count[] bytes, no partition groups).
static void test_all_occ1() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> candidates;
    for (uint32_t sid = 0; sid < 50; sid++) {
        entries.push_back({sid, 100u + sid});
        candidates.push_back(sid);
    }

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 17) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 50u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(17),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), 50u);
        for (uint32_t i = 0; i < 50; i++) {
            CHECK_EQ(positions[i], 100u + i);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// Mixed occ=1 / occ>=2 (no partition since default threshold = 8 and all
// occ counts <= 8).  Exercises both short sub-buckets.
static void test_mixed_short_buckets() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries = {
        {1, 5},                 // occ=1
        {3, 10}, {3, 20},       // occ=2 → short_occ_ge2
        {7, 99},                // occ=1
        {9, 1}, {9, 2}, {9, 3}, // occ=3
        {11, 50},               // occ=1
    };
    std::vector<uint32_t> candidates = {1, 3, 7, 9, 11};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 5) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 8u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(5),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), 8u);
        CHECK_EQ(positions[0], 5u);    // sid=1 occ=1
        CHECK_EQ(positions[1], 10u);   // sid=3
        CHECK_EQ(positions[2], 20u);
        CHECK_EQ(positions[3], 99u);   // sid=7
        CHECK_EQ(positions[4], 1u);    // sid=9
        CHECK_EQ(positions[5], 2u);
        CHECK_EQ(positions[6], 3u);
        CHECK_EQ(positions[7], 50u);   // sid=11
        reader.close();
    }
    std::remove(TEST_FILE);
}

// Partition-only case: every cluster exceeds the threshold.  Hits the
// partition path with 0 short_occ1 / short_occ_ge2 sub-buckets.
static void test_partition_only() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> candidates;
    for (uint32_t sid : {3u, 4u, 5u}) {
        for (uint32_t pos = 0; pos < 12; pos++) {
            entries.push_back({sid, 1000u + sid * 100 + pos});
        }
        candidates.push_back(sid);
    }

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);  // 12 > 8
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 21) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 36u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(21),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), 36u);
        for (size_t j = 0; j < entries.size(); j++) {
            CHECK_EQ(positions[j], entries[j].pos);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// FOR-block stream tail size sweep: 0, 1, 127 elements in the tail
// (alongside one full block of 128).  Single partition group with
// occurrence counts up to 255 (the u8 cap inherited from the dedup
// kernel).  For larger sweeps the short_occ1 sub-bucket is used in
// test_short1_full_blocks below.
static void test_partition_tail_sizes() {
    int k = 5;
    uint32_t ts = table_size(k);

    auto build_and_check = [&](uint32_t total_count) {
        std::vector<KpxWriter::PostingEntry> entries;
        entries.reserve(total_count);
        for (uint32_t pos = 0; pos < total_count; pos++) {
            entries.push_back({42u, 5000u + pos * 7u});
        }
        std::vector<uint32_t> candidates = {42};

        KpxWriter writer(k, /*freq_threshold_part=*/8);  // total_count > 8
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 9) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));

        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), total_count);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(9),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), static_cast<size_t>(total_count));
        for (uint32_t i = 0; i < total_count; i++) {
            CHECK_EQ(positions[i], 5000u + i * 7u);
        }
        reader.close();
        std::remove(TEST_FILE);
    };

    build_and_check(128);   // 1 full block, tail_count = 0
    build_and_check(129);   // 1 full block, tail_count = 1
    build_and_check(255);   // 1 full block, tail_count = 127
}

// Regression: a single (k-mer, seq_id) cluster with > 255 occurrences.
// occ_count must stay wide through the dedup → encoder pipeline: a u8
// would saturate the run at 255, drop the rest of the positions, and
// corrupt the encoder's pos_cursor walk.
static void test_partition_above_255() {
    int k = 5;
    uint32_t ts = table_size(k);

    const uint32_t occ = 1000;  // well above the short bucket's u8 cap
    std::vector<KpxWriter::PostingEntry> entries;
    entries.reserve(occ);
    for (uint32_t pos = 0; pos < occ; pos++) {
        entries.push_back({77u, 100u + pos * 11u});
    }
    std::vector<uint32_t> candidates = {77};

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);  // occ=1000 → partition
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 19) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), occ);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(19),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), static_cast<size_t>(occ));
        for (uint32_t i = 0; i < occ; i++) {
            CHECK_EQ(positions[i], 100u + i * 11u);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// 2+ full FOR blocks via the short_occ1 sub-bucket: each distinct sid
// contributes exactly 1 position, so the per-(kmer, seq_id) u8 occ cap
// does not bound the FOR stream length.
static void test_short1_full_blocks() {
    int k = 5;
    uint32_t ts = table_size(k);

    const uint32_t n = 256;  // 2 full FOR blocks, tail_count = 0
    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> candidates;
    entries.reserve(n);
    candidates.reserve(n);
    for (uint32_t sid = 0; sid < n; sid++) {
        entries.push_back({sid, 100u + sid * 3u});
        candidates.push_back(sid);
    }

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 11) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), n);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(11),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), static_cast<size_t>(n));
        for (uint32_t i = 0; i < n; i++) {
            CHECK_EQ(positions[i], 100u + i * 3u);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// distinct_count == 0 empty posting list (no k-mer matched, but writer still
// emits the all-zero blob so offsets stay valid).
static void test_empty_posting() {
    int k = 5;
    uint32_t ts = table_size(k);

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 0u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(0),
            reader.posting_file() + reader.posting_file_size(),
            {0});
        CHECK_EQ(positions.size(), 0u);
        reader.close();
    }
    std::remove(TEST_FILE);
}

// All candidates miss: posting list has data, but caller-supplied candidate
// set does not include any of its sids → every per-candidate bucket
// must come back empty.
static void test_all_candidates_miss() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries = {
        {1, 10}, {2, 20}, {3, 30}, {3, 31},
    };

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 99) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        // The actual distinct sids in the posting list are {1, 2, 3}; the
        // .kix decoded array must reflect that.  Candidates {7, 8, 9}
        // intentionally have no overlap.
        std::vector<uint32_t> kix_decoded = {1, 2, 3};
        std::vector<uint32_t> candidates  = {7, 8, 9};
        pfd::PosDecodeScratch scratch;
        PosDecoder decoder(reader.posting_file() + reader.pos_offset(99),
                           reader.posting_file() + reader.posting_file_size(),
                           kix_decoded.data(), kix_decoded.size(),
                           candidates.data(), candidates.size(),
                           &scratch);
        for (size_t i = 0; i < candidates.size(); i++) {
            CHECK_EQ(decoder.positions_for(i).size(), 0u);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

static void test_multiple_kmers() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries_a = {
        {0, 100}, {0, 200}, {1, 50},
    };
    std::vector<uint32_t> cand_a = {0, 1};

    std::vector<KpxWriter::PostingEntry> entries_b = {
        {2, 0}, {2, 1000},
    };
    std::vector<uint32_t> cand_b = {2};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 0) {
                writer.add_posting_list(i, entries_a);
            } else if (i == 10) {
                writer.add_posting_list(i, entries_b);
            } else {
                writer.add_posting_list(i, {});
            }
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 5u);

        auto pos_a = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(0),
            reader.posting_file() + reader.posting_file_size(),
            cand_a);
        CHECK_EQ(pos_a.size(), 3u);
        CHECK_EQ(pos_a[0], 100u);
        CHECK_EQ(pos_a[1], 200u);
        CHECK_EQ(pos_a[2], 50u);

        auto pos_b = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(10),
            reader.posting_file() + reader.posting_file_size(),
            cand_b);
        CHECK_EQ(pos_b.size(), 2u);
        CHECK_EQ(pos_b[0], 0u);
        CHECK_EQ(pos_b[1], 1000u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

// Partition / short bucket boundary across thresholds.  A 20-occurrence
// run on seq_id=10 plus 1 occurrence each on {11, 12, 13} (occ=1 →
// short_occ1 sub-bucket).
//   threshold = 8   -> seq_id=10 partitioned, others short_occ1
//   threshold = 16  -> seq_id=10 partitioned (20 > 16), others short_occ1
//   threshold = 100 -> all clusters short_occ1, no partition groups
static void test_partition_group_threshold() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    for (uint32_t pos = 0; pos < 20; pos++) {
        entries.push_back({10, 1000 + pos});
    }
    for (uint32_t sid : {11u, 12u, 13u}) {
        entries.push_back({sid, 7});
    }
    const uint32_t expected_count = static_cast<uint32_t>(entries.size());
    std::vector<uint32_t> candidates = {10, 11, 12, 13};

    auto run_threshold = [&](uint32_t threshold) {
        KpxWriter writer(k, threshold);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 33) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));

        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), expected_count);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(33),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), static_cast<size_t>(expected_count));
        for (size_t j = 0; j < entries.size(); j++) {
            CHECK_EQ(positions[j], entries[j].pos);
        }

        reader.close();
        std::remove(TEST_FILE);
    };

    run_threshold(8);    // seq_id=10 cluster (20 > 8) → partition
    run_threshold(16);   // seq_id=10 cluster (20 > 16) → partition
    run_threshold(100);  // nothing exceeds 100 → all short
}

int main() {
    test_single_seq();
    test_all_occ1();
    test_mixed_short_buckets();
    test_partition_only();
    test_partition_tail_sizes();
    test_partition_above_255();
    test_short1_full_blocks();
    test_empty_posting();
    test_all_candidates_miss();
    test_multiple_kmers();
    test_partition_group_threshold();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
