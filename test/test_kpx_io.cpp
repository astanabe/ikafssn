#include "test_util.hpp"
#include "index/kpx_writer.hpp"
#include "index/kpx_reader.hpp"
#include "core/config.hpp"
#include "search/posting_decoder.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.kpx";

// Decode position postings from a v7 .kpx posting blob via the
// candidate-set API.  Returns a flat vector of positions ordered by the
// candidate seq_id list (positions for candidates[0] first, then [1], etc.).
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, const uint8_t* end,
    const std::vector<uint32_t>& candidates) {
    std::vector<uint32_t> distinct_sorted = candidates;
    std::sort(distinct_sorted.begin(), distinct_sorted.end());
    distinct_sorted.erase(std::unique(distinct_sorted.begin(),
                                       distinct_sorted.end()),
                          distinct_sorted.end());

    PosDecoder decoder(data, end,
                       distinct_sorted.data(), distinct_sorted.size());
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

    // Single sequence, multiple positions for one k-mer
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
            reader.posting_data() + reader.pos_offset(42),
            reader.posting_data() + reader.posting_data_size(),
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

static void test_multiple_seqs() {
    int k = 5;
    uint32_t ts = table_size(k); // 1024

    // Multiple sequences with arbitrary positions; default threshold = 8
    // so all clusters land in the short bucket.
    std::vector<KpxWriter::PostingEntry> entries = {
        {0, 10}, {0, 20},
        {1, 5},  {1, 15}, {1, 25},
        {3, 100},
    };
    std::vector<uint32_t> candidates = {0, 1, 3};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 7) {
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
        CHECK_EQ(reader.total_position_count(), 6u);

        auto positions = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(7),
            reader.posting_data() + reader.posting_data_size(),
            candidates);

        CHECK_EQ(positions.size(), 6u);
        CHECK_EQ(positions[0], 10u);
        CHECK_EQ(positions[1], 20u);
        CHECK_EQ(positions[2], 5u);
        CHECK_EQ(positions[3], 15u);
        CHECK_EQ(positions[4], 25u);
        CHECK_EQ(positions[5], 100u);

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
            reader.posting_data() + reader.pos_offset(0),
            reader.posting_data() + reader.posting_data_size(),
            cand_a);
        CHECK_EQ(pos_a.size(), 3u);
        CHECK_EQ(pos_a[0], 100u);
        CHECK_EQ(pos_a[1], 200u);
        CHECK_EQ(pos_a[2], 50u);

        auto pos_b = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(10),
            reader.posting_data() + reader.posting_data_size(),
            cand_b);
        CHECK_EQ(pos_b.size(), 2u);
        CHECK_EQ(pos_b[0], 0u);
        CHECK_EQ(pos_b[1], 1000u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

// Phase 5g-2 / 5i: exercise the partition / short-bucket boundary across
// a few thresholds.  We build a single k-mer with a 20-occurrence run on
// seq_id=10 plus 1 occurrence each on seq_id={11, 12, 13}.
//   threshold = 8   -> seq_id=10 partitioned, others short
//   threshold = 16  -> seq_id=10 partitioned (20 > 16), others short
//   threshold = 100 -> all clusters short, no partition groups
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
            reader.posting_data() + reader.pos_offset(33),
            reader.posting_data() + reader.posting_data_size(),
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
    test_multiple_seqs();
    test_multiple_kmers();
    test_partition_group_threshold();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
