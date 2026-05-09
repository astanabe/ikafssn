#include "test_util.hpp"
#include "index/ksx_writer.hpp"
#include "index/ksx_reader.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.ksx";

// Phase 1 degenerate layout: one fragment per parent spanning the whole
// parent (frag_start=1, frag_end=parent_length).  Verifies that the
// SeqId-based read API still returns parent-equivalent answers and that
// the v10 fragment-aware accessors expose the underlying parent info.
static void test_basic_roundtrip() {
    struct SeqInfo {
        uint32_t blast_oid;
        uint32_t length;
        std::string accession;
    };
    std::vector<SeqInfo> seqs = {
        {0, 1000, "NC_000001.11"},
        {1, 500,  "NC_000002.12"},
        {2, 2500, "NM_001301717.2"},
        {3, 100,  "XR_001"},
        {4, 9999, "AB123456789"},
    };

    {
        KsxWriter writer;
        writer.set_min_seq_length(64);
        writer.set_min_length_split(0);
        writer.set_overlap_length(0);
        for (auto& s : seqs) {
            uint32_t parent_idx = writer.add_parent(s.blast_oid, s.length, s.accession);
            writer.add_fragment(parent_idx, 1, s.length);
        }
        CHECK_EQ(writer.num_sequences(), 5u);
        CHECK_EQ(writer.num_parents(), 5u);
        CHECK(writer.write(TEST_FILE));
    }

    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_sequences(), 5u);
        CHECK_EQ(reader.num_parents(), 5u);
        CHECK_EQ(reader.min_seq_length(), 64u);
        CHECK_EQ(reader.min_length_split(), 0u);
        CHECK_EQ(reader.overlap_length(), 0u);

        for (uint32_t i = 0; i < seqs.size(); i++) {
            CHECK_EQ(reader.seq_length(i), seqs[i].length);
            CHECK(reader.accession(i) == seqs[i].accession);
            CHECK_EQ(reader.parent_index(i), i);
            CHECK_EQ(reader.fragment_start(i), 1u);
            CHECK_EQ(reader.fragment_end(i), seqs[i].length);
            CHECK_EQ(reader.parent_length(i), seqs[i].length);
            CHECK_EQ(reader.blast_oid(i), seqs[i].blast_oid);
            CHECK(reader.parent_accession(i) == seqs[i].accession);
        }
        reader.close();
    }

    std::remove(TEST_FILE);
}

static void test_empty_accession() {
    {
        KsxWriter writer;
        writer.set_min_seq_length(64);
        uint32_t p0 = writer.add_parent(0, 100, "");      writer.add_fragment(p0, 1, 100);
        uint32_t p1 = writer.add_parent(1, 200, "ACC2");  writer.add_fragment(p1, 1, 200);
        uint32_t p2 = writer.add_parent(2, 300, "");      writer.add_fragment(p2, 1, 300);
        CHECK(writer.write(TEST_FILE));
    }

    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_sequences(), 3u);
        CHECK_EQ(reader.num_parents(), 3u);

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
        writer.set_min_seq_length(64);
        uint32_t p = writer.add_parent(0, 42, long_acc);
        writer.add_fragment(p, 1, 42);
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

// Multi-fragment per parent: not used by the Phase 1 builder but the
// v10 wire format is already laid out for it.  Verifies that the reader
// resolves SeqId -> parent / start / end correctly when several
// fragments share a parent.
static void test_multi_fragment_per_parent() {
    {
        KsxWriter writer;
        writer.set_min_seq_length(64);
        writer.set_min_length_split(50000);
        writer.set_overlap_length(500);
        // Parent 0: full length 100000, split into 2 fragments [1..50500] and [50001..100000].
        uint32_t p0 = writer.add_parent(0, 100000, "BIG_CHR");
        writer.add_fragment(p0, 1, 50500);          // SeqId 0
        writer.add_fragment(p0, 50001, 100000);     // SeqId 1
        // Parent 1: short, single fragment.
        uint32_t p1 = writer.add_parent(1, 200, "SHORT_ACC");
        writer.add_fragment(p1, 1, 200);            // SeqId 2

        CHECK(writer.write(TEST_FILE));
    }

    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_sequences(), 3u);
        CHECK_EQ(reader.num_parents(), 2u);
        CHECK_EQ(reader.min_seq_length(), 64u);
        CHECK_EQ(reader.min_length_split(), 50000u);
        CHECK_EQ(reader.overlap_length(), 500u);

        // SeqId 0 -> first fragment of parent 0
        CHECK_EQ(reader.parent_index(0), 0u);
        CHECK_EQ(reader.fragment_start(0), 1u);
        CHECK_EQ(reader.fragment_end(0), 50500u);
        CHECK_EQ(reader.seq_length(0), 50500u);
        CHECK(reader.accession(0) == "BIG_CHR");
        CHECK_EQ(reader.parent_length(0), 100000u);
        CHECK_EQ(reader.blast_oid(0), 0u);

        // SeqId 1 -> second fragment of parent 0
        CHECK_EQ(reader.parent_index(1), 0u);
        CHECK_EQ(reader.fragment_start(1), 50001u);
        CHECK_EQ(reader.fragment_end(1), 100000u);
        CHECK_EQ(reader.seq_length(1), 50000u);
        CHECK(reader.accession(1) == "BIG_CHR");

        // SeqId 2 -> only fragment of parent 1
        CHECK_EQ(reader.parent_index(2), 1u);
        CHECK_EQ(reader.fragment_start(2), 1u);
        CHECK_EQ(reader.fragment_end(2), 200u);
        CHECK(reader.accession(2) == "SHORT_ACC");
        CHECK_EQ(reader.parent_length(1), 200u);
        CHECK_EQ(reader.blast_oid(1), 1u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

int main() {
    test_basic_roundtrip();
    test_empty_accession();
    test_long_accession();
    test_multi_fragment_per_parent();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
