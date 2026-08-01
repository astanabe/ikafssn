#include "test_util.hpp"
#include "kix_writer.hpp"
#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"
#include "search/seq_id_decoder.hpp"

#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

using namespace ikafssn;

// Per-tier path: the SIMD variants of this test can run concurrently.
static const std::string TEST_FILE = test_tmpdir("/tmp/test_ikafssn") + ".kix";

// Collect a whole .kix posting list into a vector.
static std::vector<uint32_t> decode_id_postings(const uint8_t* data, uint64_t byte_len) {
    SeqIdDecoder decoder(data, data + byte_len);
    decoder.ensure_decoded();
    return std::vector<uint32_t>(decoder.decoded_data(),
                                 decoder.decoded_data() + decoder.decoded_count());
}

static void test_k7_uint16() {
    int k = 7;
    uint8_t kmer_type = 0; // uint16_t
    uint32_t ts = table_size(k); // 4^7 = 16384

    // Create synthetic posting list bytes for a few k-mers.  v7 .kix stores
    // distinct seq_ids only — duplicates are not allowed at the writer
    // boundary (intra-sequence k-mer duplicates are removed by the
    // builder's SIMD dedup kernel).
    std::vector<std::vector<uint32_t>> postings(ts);
    postings[0] = {0, 1, 2, 5, 10};           // k-mer 0
    postings[1] = {3, 7};                       // k-mer 1
    postings[100] = {0, 1, 2, 3, 4, 5};          // 6 distinct seq_ids
    postings[ts - 1] = {999};                    // last k-mer

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(1000);
        writer.set_db("testdb");
        writer.set_volume_info(0, 1);
        writer.set_flags(KIX_FLAG_HAS_KSX);

        for (uint32_t i = 0; i < ts; i++) {
            writer.add_posting_list(i, postings[i]);
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KixReader reader;
        CHECK(reader.open(TEST_FILE));

        // Check header
        CHECK_EQ(reader.k(), k);
        CHECK_EQ(reader.kmer_type(), kmer_type);
        CHECK_EQ(reader.num_sequences(), 1000u);
        CHECK_EQ(reader.table_size(), ts);
        CHECK(std::memcmp(reader.header().magic, KIX_MAGIC, sizeof(KIX_MAGIC)) == 0);
        CHECK(std::string(reader.header().db, reader.header().db_len) == "testdb");

        // Check byte-lengths (non-empty k-mers have > 0 byte length)
        CHECK(reader.posting_list_byte_length(0) > 0);
        CHECK(reader.posting_list_byte_length(1) > 0);
        CHECK_EQ(reader.posting_list_byte_length(2), 0u);
        CHECK(reader.posting_list_byte_length(100) > 0);
        CHECK(reader.posting_list_byte_length(ts - 1) > 0);

        // Check counts via bulk_count_postings
        auto counts = reader.bulk_count_postings();
        CHECK_EQ(counts[0], 5u);
        CHECK_EQ(counts[1], 2u);
        CHECK_EQ(counts[2], 0u);
        CHECK_EQ(counts[100], 6u);
        CHECK_EQ(counts[ts - 1], 1u);

        // Decode and verify postings
        auto decoded0 = decode_id_postings(
            reader.posting_file() + reader.posting_list_offset(0),
            reader.posting_list_byte_length(0));
        CHECK_EQ(decoded0.size(), 5u);
        CHECK_EQ(decoded0[0], 0u);
        CHECK_EQ(decoded0[1], 1u);
        CHECK_EQ(decoded0[2], 2u);
        CHECK_EQ(decoded0[3], 5u);
        CHECK_EQ(decoded0[4], 10u);

        auto decoded1 = decode_id_postings(
            reader.posting_file() + reader.posting_list_offset(1),
            reader.posting_list_byte_length(1));
        CHECK_EQ(decoded1.size(), 2u);
        CHECK_EQ(decoded1[0], 3u);
        CHECK_EQ(decoded1[1], 7u);

        auto decoded100 = decode_id_postings(
            reader.posting_file() + reader.posting_list_offset(100),
            reader.posting_list_byte_length(100));
        CHECK_EQ(decoded100.size(), 6u);
        CHECK_EQ(decoded100[0], 0u);
        CHECK_EQ(decoded100[1], 1u);
        CHECK_EQ(decoded100[2], 2u);
        CHECK_EQ(decoded100[3], 3u);
        CHECK_EQ(decoded100[4], 4u);
        CHECK_EQ(decoded100[5], 5u);

        auto decoded_last = decode_id_postings(
            reader.posting_file() + reader.posting_list_offset(ts - 1),
            reader.posting_list_byte_length(ts - 1));
        CHECK_EQ(decoded_last.size(), 1u);
        CHECK_EQ(decoded_last[0], 999u);

        // Check total distinct postings
        CHECK_EQ(reader.total_distinct_postings(), 5u + 2u + 6u + 1u);

        reader.close();
    }

    std::remove(TEST_FILE.c_str());
}

static void test_k9_uint32() {
    int k = 9;
    uint8_t kmer_type = 1; // uint32_t
    uint32_t ts = table_size(k); // 4^9 = 262144

    // Sparse: only a few k-mers have postings
    std::vector<std::pair<uint32_t, std::vector<uint32_t>>> sparse_postings = {
        {0, {0}},
        {1000, {5, 10, 15, 20}},
        {ts - 1, {100, 200}},
    };

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(300);

        for (uint32_t i = 0; i < ts; i++) {
            std::vector<uint32_t> ids;
            for (auto& [kmer, posts] : sparse_postings) {
                if (kmer == i) {
                    ids = posts;
                    break;
                }
            }
            writer.add_posting_list(i, ids);
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KixReader reader;
        CHECK(reader.open(TEST_FILE));

        CHECK_EQ(reader.k(), k);
        CHECK_EQ(reader.kmer_type(), 1u);
        CHECK_EQ(reader.table_size(), ts);

        // Verify posting counts via bulk decode
        auto counts = reader.bulk_count_postings();
        CHECK_EQ(counts[0], 1u);
        CHECK_EQ(counts[1], 0u);
        CHECK_EQ(counts[1000], 4u);
        CHECK_EQ(counts[ts - 1], 2u);

        // Decode k-mer 1000
        auto decoded = decode_id_postings(
            reader.posting_file() + reader.posting_list_offset(1000),
            reader.posting_list_byte_length(1000));
        CHECK_EQ(decoded.size(), 4u);
        CHECK_EQ(decoded[0], 5u);
        CHECK_EQ(decoded[1], 10u);
        CHECK_EQ(decoded[2], 15u);
        CHECK_EQ(decoded[3], 20u);

        reader.close();
    }

    std::remove(TEST_FILE.c_str());
}

static void test_empty_postings() {
    int k = 5;
    uint8_t kmer_type = 0;
    uint32_t ts = table_size(k); // 4^5 = 1024

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(0);
        // All posting lists empty
        for (uint32_t i = 0; i < ts; i++) {
            writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KixReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_distinct_postings(), 0u);
        for (uint32_t i = 0; i < ts; i++) {
            CHECK_EQ(reader.posting_list_byte_length(i), 0u);
        }
        reader.close();
    }

    std::remove(TEST_FILE.c_str());
}

// The decoder reconstructs deltas a SIMD register at a time and finishes
// the remainder element by element, so lengths around a register width or
// a 128-element codec block catch a wrong lane count or a mis-seeded
// running total.  Uneven gaps make every seq_id depend on all before it.
static void test_delta_reconstruction_lengths() {
    int k = 7;
    uint8_t kmer_type = 0;
    uint32_t ts = table_size(k);

    const std::vector<uint32_t> lengths = {1, 2, 3, 4, 5, 15, 16, 17,
                                           127, 128, 129, 1000};

    // Distinct, strictly increasing seq_ids with an uneven gap pattern.
    std::vector<std::vector<uint32_t>> postings(ts);
    for (size_t li = 0; li < lengths.size(); li++) {
        std::vector<uint32_t>& ids = postings[li];
        uint32_t sid = static_cast<uint32_t>(li);
        for (uint32_t i = 0; i < lengths[li]; i++) {
            ids.push_back(sid);
            sid += 1 + ((i * 37 + li * 11) % 97);
        }
    }

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(200000);
        for (uint32_t i = 0; i < ts; i++) {
            writer.add_posting_list(i, postings[i]);
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KixReader reader;
        CHECK(reader.open(TEST_FILE));
        for (size_t li = 0; li < lengths.size(); li++) {
            auto decoded = decode_id_postings(
                reader.posting_file() + reader.posting_list_offset(li),
                reader.posting_list_byte_length(li));
            CHECK_EQ(decoded.size(), size_t(lengths[li]));
            CHECK(decoded == postings[li]);
        }
        reader.close();
    }

    std::remove(TEST_FILE.c_str());
}

int main() {
    test_k7_uint16();
    test_k9_uint32();
    test_empty_postings();
    test_delta_reconstruction_lengths();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
