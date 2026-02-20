#include "test_util.hpp"
#include "index/ksx_writer.hpp"
#include "index/ksx_reader.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.ksx";

static void test_basic_roundtrip() {
    // Write 5 sequences
    struct SeqInfo {
        uint32_t length;
        std::string accession;
    };
    std::vector<SeqInfo> seqs = {
        {1000, "NC_000001.11"},
        {500,  "NC_000002.12"},
        {2500, "NM_001301717.2"},
        {100,  "XR_001"},
        {9999, "AB123456789"},
    };

    {
        KsxWriter writer;
        for (auto& s : seqs) {
            writer.add_sequence(s.length, s.accession);
        }
        CHECK_EQ(writer.num_sequences(), 5u);
        CHECK(writer.write(TEST_FILE));
    }

    // Read back
    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_sequences(), 5u);

        for (uint32_t i = 0; i < seqs.size(); i++) {
            CHECK_EQ(reader.seq_length(i), seqs[i].length);
            CHECK(reader.accession(i) == seqs[i].accession);
        }
        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_empty_accession() {
    {
        KsxWriter writer;
        writer.add_sequence(100, "");
        writer.add_sequence(200, "ACC2");
        writer.add_sequence(300, "");
        CHECK(writer.write(TEST_FILE));
    }

    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_sequences(), 3u);

        CHECK_EQ(reader.seq_length(0), 100u);
        CHECK(reader.accession(0) == "");

        CHECK_EQ(reader.seq_length(1), 200u);
        CHECK(reader.accession(1) == "ACC2");

        CHECK_EQ(reader.seq_length(2), 300u);
        CHECK(reader.accession(2) == "");

        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_long_accession() {
    std::string long_acc(200, 'X');

    {
        KsxWriter writer;
        writer.add_sequence(42, long_acc);
        CHECK(writer.write(TEST_FILE));
    }

    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_sequences(), 1u);
        CHECK_EQ(reader.seq_length(0), 42u);
        CHECK(reader.accession(0) == long_acc);
        reader.close();
    }

    std::remove(TEST_FILE);
}

int main() {
    test_basic_roundtrip();
    test_empty_accession();
    test_long_accession();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
