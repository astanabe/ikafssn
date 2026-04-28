#include "test_util.hpp"
#include "index/kpx_writer.hpp"
#include "index/kpx_reader.hpp"
#include "core/config.hpp"
#include "search/posting_decoder.hpp"

#include <cstdio>
#include <vector>
#include <string>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.kpx";

// Decode position postings from a v6 .kpx posting blob.  v6 needs the
// caller-side seq_id stream so that the partition+short merge can walk
// lock-step with it; tests pass the per-(kmer) seq_id list explicitly.
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, const uint8_t* end,
    const std::vector<uint32_t>& seq_ids,
    uint32_t expected_count) {

    std::vector<uint32_t> result;
    PosDecoder decoder(data, end, seq_ids.data(), seq_ids.size());
    while (decoder.has_more() && result.size() < expected_count) {
        result.push_back(decoder.next());
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
    std::vector<uint32_t> seq_ids = {5, 5, 5, 5};

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
        CHECK_EQ(reader.total_postings(), 4u);

        auto positions = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(42),
            reader.posting_data() + reader.posting_data_size(),
            seq_ids,
            4);

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
    std::vector<uint32_t> seq_ids = {0, 0, 1, 1, 1, 3};

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
        CHECK_EQ(reader.total_postings(), 6u);

        auto positions = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(7),
            reader.posting_data() + reader.posting_data_size(),
            seq_ids,
            6);

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

    // Two k-mers with different posting patterns
    std::vector<KpxWriter::PostingEntry> entries_a = {
        {0, 100}, {0, 200}, {1, 50},
    };
    std::vector<uint32_t> sids_a = {0, 0, 1};

    std::vector<KpxWriter::PostingEntry> entries_b = {
        {2, 0}, {2, 1000},
    };
    std::vector<uint32_t> sids_b = {2, 2};

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
        CHECK_EQ(reader.total_postings(), 5u);

        // k-mer 0
        auto pos_a = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(0),
            reader.posting_data() + reader.posting_data_size(),
            sids_a,
            3);
        CHECK_EQ(pos_a.size(), 3u);
        CHECK_EQ(pos_a[0], 100u);
        CHECK_EQ(pos_a[1], 200u);
        CHECK_EQ(pos_a[2], 50u);

        // k-mer 10
        auto pos_b = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(10),
            reader.posting_data() + reader.posting_data_size(),
            sids_b,
            2);
        CHECK_EQ(pos_b.size(), 2u);
        CHECK_EQ(pos_b[0], 0u);
        CHECK_EQ(pos_b[1], 1000u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

// Phase 5g-2: exercise the partition / short-bucket boundary across a
// few thresholds.  We build a single k-mer with a 20-occurrence run on
// seq_id=10 plus 1 occurrence each on seq_id={11, 12, 13}.
//   threshold = 8   -> seq_id=10 partitioned, others short
//   threshold = 16  -> seq_id=10 partitioned (20 > 16), others short
//   threshold = 100 -> all clusters short, no partition groups
static void test_partition_group_threshold() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> seq_ids;
    for (uint32_t pos = 0; pos < 20; pos++) {
        entries.push_back({10, 1000 + pos});
        seq_ids.push_back(10);
    }
    for (uint32_t sid : {11u, 12u, 13u}) {
        entries.push_back({sid, 7});
        seq_ids.push_back(sid);
    }
    const uint32_t expected_count = static_cast<uint32_t>(entries.size());

    auto run_threshold = [&](uint32_t threshold) {
        KpxWriter writer(k, threshold);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 33) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));

        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_postings(), expected_count);

        auto positions = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(33),
            reader.posting_data() + reader.posting_data_size(),
            seq_ids,
            expected_count);

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
