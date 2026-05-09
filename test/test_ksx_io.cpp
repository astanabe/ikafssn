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

// Phase 5 round-trip: multiple parents that all carry several fragments,
// exercising the cross-parent boundary as well as 3+ fragments per parent.
// The Phase 1 builder still emits the degenerate 1:1 layout, so this test
// stresses the wire format / reader resolution end-to-end against the kind
// of layout the Phase 2 builder writes once -min_length_split is non-zero.
static void test_multi_parent_multi_fragment() {
    {
        KsxWriter writer;
        writer.set_min_seq_length(64);
        writer.set_min_length_split(1500);
        writer.set_overlap_length(200);

        // Parent 0 (BLAST OID 7): 3 fragments mirroring the LONGCHR_1
        // expected split (5171 bp, min=1500, ovl=200).
        uint32_t p0 = writer.add_parent(7, 5171, "LONGCHR_1");
        writer.add_fragment(p0, 1,    1791);  // SeqId 0
        writer.add_fragment(p0, 1592, 3381);  // SeqId 1
        writer.add_fragment(p0, 3182, 4771);  // SeqId 2

        // Parent 1 (BLAST OID 11): 2 fragments, asymmetric lengths.
        uint32_t p1 = writer.add_parent(11, 4000, "MID_CHR");
        writer.add_fragment(p1, 1,    2200);  // SeqId 3
        writer.add_fragment(p1, 2001, 4000);  // SeqId 4

        // Parent 2 (BLAST OID 12): single fragment.
        uint32_t p2 = writer.add_parent(12, 800, "SHORT");
        writer.add_fragment(p2, 1, 800);     // SeqId 5

        // Parent 3 (BLAST OID 21): 4 fragments — exercises the case where
        // a parent has more than 3 children of unequal length.
        uint32_t p3 = writer.add_parent(21, 8000, "LONGCHR_2");
        writer.add_fragment(p3, 1,    2200);  // SeqId 6
        writer.add_fragment(p3, 2001, 4200);  // SeqId 7
        writer.add_fragment(p3, 4001, 6200);  // SeqId 8
        writer.add_fragment(p3, 6001, 8000);  // SeqId 9

        CHECK_EQ(writer.num_parents(), 4u);
        CHECK_EQ(writer.num_sequences(), 10u);
        CHECK(writer.write(TEST_FILE));
    }

    {
        KsxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.num_parents(), 4u);
        CHECK_EQ(reader.num_sequences(), 10u);
        CHECK_EQ(reader.min_length_split(), 1500u);
        CHECK_EQ(reader.overlap_length(), 200u);

        // Per-fragment expectations: parent_idx, frag_start, frag_end,
        // frag_length, parent_length, parent_accession, blast_oid.
        struct FragExp {
            uint32_t parent_idx;
            uint32_t fstart;
            uint32_t fend;
            uint32_t fragment_length;
            uint32_t parent_length;
            const char* parent_accession;
            uint32_t blast_oid;
        };
        FragExp exp[10] = {
            {0, 1,    1791, 1791, 5171, "LONGCHR_1", 7},
            {0, 1592, 3381, 1790, 5171, "LONGCHR_1", 7},
            {0, 3182, 4771, 1590, 5171, "LONGCHR_1", 7},
            {1, 1,    2200, 2200, 4000, "MID_CHR",   11},
            {1, 2001, 4000, 2000, 4000, "MID_CHR",   11},
            {2, 1,     800,  800,  800, "SHORT",     12},
            {3, 1,    2200, 2200, 8000, "LONGCHR_2", 21},
            {3, 2001, 4200, 2200, 8000, "LONGCHR_2", 21},
            {3, 4001, 6200, 2200, 8000, "LONGCHR_2", 21},
            {3, 6001, 8000, 2000, 8000, "LONGCHR_2", 21},
        };
        for (uint32_t sid = 0; sid < 10; ++sid) {
            CHECK_EQ(reader.parent_index(sid),   exp[sid].parent_idx);
            CHECK_EQ(reader.fragment_start(sid), exp[sid].fstart);
            CHECK_EQ(reader.fragment_end(sid),   exp[sid].fend);
            CHECK_EQ(reader.seq_length(sid),     exp[sid].fragment_length);
            CHECK(reader.accession(sid) == exp[sid].parent_accession);
        }
        // Parent-keyed accessors round-trip independently.
        for (uint32_t pidx = 0; pidx < reader.num_parents(); ++pidx) {
            // Find the first sid that maps to this parent, and reuse its
            // expected parent_length / parent_accession / blast_oid.
            for (uint32_t sid = 0; sid < 10; ++sid) {
                if (exp[sid].parent_idx != pidx) continue;
                CHECK_EQ(reader.parent_length(pidx), exp[sid].parent_length);
                CHECK_EQ(reader.blast_oid(pidx),     exp[sid].blast_oid);
                CHECK(reader.parent_accession(pidx) == exp[sid].parent_accession);
                break;
            }
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

int main() {
    test_basic_roundtrip();
    test_empty_accession();
    test_long_accession();
    test_multi_fragment_per_parent();
    test_multi_parent_multi_fragment();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
