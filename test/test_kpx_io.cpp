#include "test_util.hpp"
#include "index/kpx_writer.hpp"
#include "index/kpx_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <cstdio>
#include <vector>
#include <string>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.kpx";

// Decode position postings given the corresponding ID deltas.
// id_deltas[i]: delta of seq_id at position i (0 = same seq, >0 = new seq).
// Returns decoded positions.
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, uint32_t count,
    const std::vector<uint32_t>& id_deltas) {

    std::vector<uint32_t> result;
    result.reserve(count);
    size_t offset = 0;

    for (uint32_t i = 0; i < count; i++) {
        uint32_t val;
        offset += varint_decode(data + offset, val);
        if (i == 0 || id_deltas[i] != 0) {
            // First entry or new sequence: raw pos
            result.push_back(val);
        } else {
            // Same sequence: delta from previous pos
            result.push_back(result.back() + val);
        }
    }
    return result;
}

static void test_single_seq() {
    int k = 7;
    uint64_t ts = table_size(k);

    // Single sequence, multiple positions for one k-mer
    // seq_id=5, positions: 10, 20, 30, 100
    std::vector<KpxWriter::PostingEntry> entries = {
        {5, 10}, {5, 20}, {5, 30}, {5, 100}
    };

    {
        KpxWriter writer(k);
        for (uint64_t i = 0; i < ts; i++) {
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

        // Decode positions for k-mer 42
        // All same seq_id -> id_deltas = {raw, 0, 0, 0}
        // First id delta is conceptually "raw" (nonzero since it's the first)
        std::vector<uint32_t> id_deltas = {5, 0, 0, 0}; // first is raw (nonzero)

        auto positions = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(42),
            4, id_deltas);

        CHECK_EQ(positions.size(), 4u);
        CHECK_EQ(positions[0], 10u);
        CHECK_EQ(positions[1], 20u);
        CHECK_EQ(positions[2], 30u);
        CHECK_EQ(positions[3], 100u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_seq_boundary_reset() {
    int k = 5;
    uint64_t ts = table_size(k); // 1024

    // Multiple sequences with position delta reset at boundaries
    // seq 0: pos 10, 20
    // seq 1: pos 5, 15       <-- delta resets here
    // seq 1: pos 25           <-- delta from 15
    // seq 3: pos 100          <-- delta resets here
    std::vector<KpxWriter::PostingEntry> entries = {
        {0, 10}, {0, 20},
        {1, 5},  {1, 15}, {1, 25},
        {3, 100},
    };

    // Corresponding ID deltas (as would be in .kix)
    // raw(0), delta(0), delta(1), delta(0), delta(0), delta(2)
    std::vector<uint32_t> id_deltas = {0, 0, 1, 0, 0, 2};
    // But the first entry always uses raw, so id_deltas[0] is "raw" regardless.
    // In our decoder, i==0 -> raw. For i>0, id_deltas[i] != 0 -> raw.

    {
        KpxWriter writer(k);
        for (uint64_t i = 0; i < ts; i++) {
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
            6, id_deltas);

        CHECK_EQ(positions.size(), 6u);
        CHECK_EQ(positions[0], 10u);   // raw
        CHECK_EQ(positions[1], 20u);   // delta 10 from 10
        CHECK_EQ(positions[2], 5u);    // raw (new seq)
        CHECK_EQ(positions[3], 15u);   // delta 10 from 5
        CHECK_EQ(positions[4], 25u);   // delta 10 from 15
        CHECK_EQ(positions[5], 100u);  // raw (new seq)

        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_multiple_kmers() {
    int k = 5;
    uint64_t ts = table_size(k);

    // Two k-mers with different posting patterns
    std::vector<KpxWriter::PostingEntry> entries_a = {
        {0, 100}, {0, 200}, {1, 50},
    };
    std::vector<KpxWriter::PostingEntry> entries_b = {
        {2, 0}, {2, 1000},
    };

    {
        KpxWriter writer(k);
        for (uint64_t i = 0; i < ts; i++) {
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
        std::vector<uint32_t> id_deltas_a = {0, 0, 1};
        auto pos_a = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(0),
            3, id_deltas_a);
        CHECK_EQ(pos_a[0], 100u);
        CHECK_EQ(pos_a[1], 200u);
        CHECK_EQ(pos_a[2], 50u); // raw (new seq)

        // k-mer 10
        std::vector<uint32_t> id_deltas_b = {2, 0};
        auto pos_b = decode_pos_postings(
            reader.posting_data() + reader.pos_offset(10),
            2, id_deltas_b);
        CHECK_EQ(pos_b[0], 0u);
        CHECK_EQ(pos_b[1], 1000u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

int main() {
    test_single_seq();
    test_seq_boundary_reset();
    test_multiple_kmers();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
