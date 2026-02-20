#include "test_util.hpp"
#include "index/kix_writer.hpp"
#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.kix";

// Decode a delta-compressed posting list from raw bytes
static std::vector<uint32_t> decode_id_postings(const uint8_t* data, uint32_t count) {
    std::vector<uint32_t> result;
    result.reserve(count);
    size_t offset = 0;
    uint32_t prev = 0;
    for (uint32_t i = 0; i < count; i++) {
        uint32_t delta;
        offset += varint_decode(data + offset, delta);
        if (i == 0) {
            prev = delta; // first is raw
        } else {
            prev += delta;
        }
        result.push_back(prev);
    }
    return result;
}

static void test_k7_uint16() {
    int k = 7;
    uint8_t kmer_type = 0; // uint16_t
    uint64_t ts = table_size(k); // 4^7 = 16384

    // Create synthetic posting data for a few k-mers
    std::vector<std::vector<uint32_t>> postings(ts);
    postings[0] = {0, 1, 2, 5, 10};           // k-mer 0
    postings[1] = {3, 7};                       // k-mer 1
    postings[100] = {0, 0, 0, 1, 1, 2};         // repeated IDs (delta=0)
    postings[ts - 1] = {999};                    // last k-mer

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(1000);
        writer.set_db_name("testdb");
        writer.set_volume_info(0, 1);
        writer.set_flags(KIX_FLAG_HAS_KSX);

        for (uint64_t i = 0; i < ts; i++) {
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
        CHECK(std::memcmp(reader.header().magic, KIX_MAGIC, 4) == 0);
        CHECK(std::string(reader.header().db_name, reader.header().db_name_len) == "testdb");

        // Check counts
        CHECK_EQ(reader.posting_count(0), 5u);
        CHECK_EQ(reader.posting_count(1), 2u);
        CHECK_EQ(reader.posting_count(2), 0u);
        CHECK_EQ(reader.posting_count(100), 6u);
        CHECK_EQ(reader.posting_count(ts - 1), 1u);

        // Decode and verify postings
        auto decoded0 = decode_id_postings(
            reader.posting_data() + reader.posting_offset(0),
            reader.posting_count(0));
        CHECK_EQ(decoded0.size(), 5u);
        CHECK_EQ(decoded0[0], 0u);
        CHECK_EQ(decoded0[1], 1u);
        CHECK_EQ(decoded0[2], 2u);
        CHECK_EQ(decoded0[3], 5u);
        CHECK_EQ(decoded0[4], 10u);

        auto decoded1 = decode_id_postings(
            reader.posting_data() + reader.posting_offset(1),
            reader.posting_count(1));
        CHECK_EQ(decoded1.size(), 2u);
        CHECK_EQ(decoded1[0], 3u);
        CHECK_EQ(decoded1[1], 7u);

        auto decoded100 = decode_id_postings(
            reader.posting_data() + reader.posting_offset(100),
            reader.posting_count(100));
        CHECK_EQ(decoded100.size(), 6u);
        CHECK_EQ(decoded100[0], 0u);
        CHECK_EQ(decoded100[1], 0u);
        CHECK_EQ(decoded100[2], 0u);
        CHECK_EQ(decoded100[3], 1u);
        CHECK_EQ(decoded100[4], 1u);
        CHECK_EQ(decoded100[5], 2u);

        auto decoded_last = decode_id_postings(
            reader.posting_data() + reader.posting_offset(ts - 1),
            reader.posting_count(ts - 1));
        CHECK_EQ(decoded_last.size(), 1u);
        CHECK_EQ(decoded_last[0], 999u);

        // Check total postings
        CHECK_EQ(reader.total_postings(), 5u + 2u + 6u + 1u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_k9_uint32() {
    int k = 9;
    uint8_t kmer_type = 1; // uint32_t
    uint64_t ts = table_size(k); // 4^9 = 262144

    // Sparse: only a few k-mers have postings
    std::vector<std::pair<uint64_t, std::vector<uint32_t>>> sparse_postings = {
        {0, {0}},
        {1000, {5, 10, 15, 20}},
        {ts - 1, {100, 200}},
    };

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(300);

        for (uint64_t i = 0; i < ts; i++) {
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

        // Verify posting counts
        CHECK_EQ(reader.posting_count(0), 1u);
        CHECK_EQ(reader.posting_count(1), 0u);
        CHECK_EQ(reader.posting_count(1000), 4u);
        CHECK_EQ(reader.posting_count(ts - 1), 2u);

        // Decode k-mer 1000
        auto decoded = decode_id_postings(
            reader.posting_data() + reader.posting_offset(1000),
            reader.posting_count(1000));
        CHECK_EQ(decoded.size(), 4u);
        CHECK_EQ(decoded[0], 5u);
        CHECK_EQ(decoded[1], 10u);
        CHECK_EQ(decoded[2], 15u);
        CHECK_EQ(decoded[3], 20u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_empty_postings() {
    int k = 5;
    uint8_t kmer_type = 0;
    uint64_t ts = table_size(k); // 4^5 = 1024

    {
        KixWriter writer(k, kmer_type);
        writer.set_num_sequences(0);
        // All posting lists empty
        for (uint64_t i = 0; i < ts; i++) {
            writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KixReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_postings(), 0u);
        for (uint64_t i = 0; i < ts; i++) {
            CHECK_EQ(reader.posting_count(i), 0u);
        }
        reader.close();
    }

    std::remove(TEST_FILE);
}

int main() {
    test_k7_uint16();
    test_k9_uint32();
    test_empty_postings();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
