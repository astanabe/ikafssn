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

// Decode position postings from a v4 .kpx posting blob.  v4 stores absolute
// positions (no within-seq delta reset), so the legacy id_deltas argument
// from v3 is no longer required.
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, const uint8_t* end, uint32_t expected_count) {

    std::vector<uint32_t> result;
    PosDecoder decoder(data, end);
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

    // Multiple sequences with arbitrary positions.  v4 encodes absolute
    // positions, so no delta-reset bookkeeping is needed.
    std::vector<KpxWriter::PostingEntry> entries = {
        {0, 10}, {0, 20},
        {1, 5},  {1, 15}, {1, 25},
        {3, 100},
    };

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
    std::vector<KpxWriter::PostingEntry> entries_b = {
        {2, 0}, {2, 1000},
    };

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
            3);
        CHECK_EQ(pos_a.size(), 3u);
        CHECK_EQ(pos_a[0], 100u);
        CHECK_EQ(pos_a[1], 200u);
        CHECK_EQ(pos_a[2], 50u);

        // k-mer 10
        auto pos_b = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(10),
            reader.posting_data() + reader.posting_data_size(),
            2);
        CHECK_EQ(pos_b.size(), 2u);
        CHECK_EQ(pos_b[0], 0u);
        CHECK_EQ(pos_b[1], 1000u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

int main() {
    test_single_seq();
    test_multiple_seqs();
    test_multiple_kmers();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
