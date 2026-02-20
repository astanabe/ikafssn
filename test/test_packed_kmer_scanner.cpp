#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "core/ambiguity_parser.hpp"
#include "core/packed_kmer_scanner.hpp"
#include "core/kmer_encoding.hpp"
#include "io/blastdb_reader.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "core/varint.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <filesystem>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_ambigdb_path;
static std::string g_output_dir;

// ==================== AmbiguityParser tests ====================

static void test_ambiguity_parser_empty() {
    std::fprintf(stderr, "-- test_ambiguity_parser_empty\n");
    auto entries = AmbiguityParser::parse(nullptr, 0);
    CHECK_EQ(entries.size(), size_t(0));
}

static void test_ambiguity_parser_old_format() {
    std::fprintf(stderr, "-- test_ambiguity_parser_old_format\n");

    // Construct old format ambiguity data:
    // Header: bit31=0, bits30-0 = num_entries = 2
    // Entry 0: ncbi4na=15 (N), run_len-1=0 (run=1), position=7
    // Entry 1: ncbi4na=5 (R=A|G), run_len-1=2 (run=3), position=20

    uint8_t data[12]; // 4 header + 2 * 4 bytes
    // Header: 0x00000002 (2 entries, old format)
    data[0] = 0x00; data[1] = 0x00; data[2] = 0x00; data[3] = 0x02;

    // Entry 0: ncbi4na=15, run_len-1=0, pos=7
    // word = (15 << 28) | (0 << 24) | 7 = 0xF0000007
    data[4] = 0xF0; data[5] = 0x00; data[6] = 0x00; data[7] = 0x07;

    // Entry 1: ncbi4na=5, run_len-1=2, pos=20
    // word = (5 << 28) | (2 << 24) | 20 = 0x52000014
    data[8] = 0x52; data[9] = 0x00; data[10] = 0x00; data[11] = 0x14;

    auto entries = AmbiguityParser::parse(reinterpret_cast<const char*>(data), 12);
    CHECK_EQ(entries.size(), size_t(2));

    CHECK_EQ(entries[0].ncbi4na, uint8_t(15));
    CHECK_EQ(entries[0].run_length, uint32_t(1));
    CHECK_EQ(entries[0].position, uint32_t(7));

    CHECK_EQ(entries[1].ncbi4na, uint8_t(5));
    CHECK_EQ(entries[1].run_length, uint32_t(3));
    CHECK_EQ(entries[1].position, uint32_t(20));
}

static void test_ambiguity_parser_new_format() {
    std::fprintf(stderr, "-- test_ambiguity_parser_new_format\n");

    // Construct new format ambiguity data:
    // Header: bit31=1, bits30-0 = num_words = 2 (1 entry * 2 words)
    // Entry 0: ncbi4na=15 (N), run_len-1=99 (run=100), position=1000000

    uint8_t data[12]; // 4 header + 1 * 8 bytes
    // Header: 0x80000002
    data[0] = 0x80; data[1] = 0x00; data[2] = 0x00; data[3] = 0x02;

    // Word0: ncbi4na=15, run_len-1=99, padding=0
    // w0 = (15 << 28) | (99 << 16) | 0 = 0xF0630000
    data[4] = 0xF0; data[5] = 0x63; data[6] = 0x00; data[7] = 0x00;

    // Word1: position=1000000 = 0x000F4240
    data[8] = 0x00; data[9] = 0x0F; data[10] = 0x42; data[11] = 0x40;

    auto entries = AmbiguityParser::parse(reinterpret_cast<const char*>(data), 12);
    CHECK_EQ(entries.size(), size_t(1));
    CHECK_EQ(entries[0].ncbi4na, uint8_t(15));
    CHECK_EQ(entries[0].run_length, uint32_t(100));
    CHECK_EQ(entries[0].position, uint32_t(1000000));
}

static void test_ambiguity_parser_sort() {
    std::fprintf(stderr, "-- test_ambiguity_parser_sort\n");

    // Entries deliberately out of order: pos=20 before pos=7
    uint8_t data[12];
    data[0] = 0x00; data[1] = 0x00; data[2] = 0x00; data[3] = 0x02;

    // Entry 0: ncbi4na=5, run_len-1=0, pos=20
    data[4] = 0x50; data[5] = 0x00; data[6] = 0x00; data[7] = 0x14;

    // Entry 1: ncbi4na=15, run_len-1=0, pos=7
    data[8] = 0xF0; data[9] = 0x00; data[10] = 0x00; data[11] = 0x07;

    auto entries = AmbiguityParser::parse(reinterpret_cast<const char*>(data), 12);
    CHECK_EQ(entries.size(), size_t(2));
    // After sort, position 7 should come first
    CHECK_EQ(entries[0].position, uint32_t(7));
    CHECK_EQ(entries[1].position, uint32_t(20));
}

// ==================== PackedKmerScanner tests ====================

// Helper: encode a DNA string to ncbi2na packed format
static std::vector<char> encode_ncbi2na(const char* seq, uint32_t len) {
    uint32_t num_bytes = (len + 3) / 4;
    std::vector<char> packed(num_bytes, 0);
    for (uint32_t i = 0; i < len; i++) {
        uint8_t code;
        switch (seq[i]) {
            case 'A': case 'a': code = 0; break;
            case 'C': case 'c': code = 1; break;
            case 'G': case 'g': code = 2; break;
            case 'T': case 't': code = 3; break;
            default: code = 0; break; // placeholder for ambig
        }
        int byte_idx = i / 4;
        int bit_shift = 6 - 2 * (i % 4);
        packed[byte_idx] |= static_cast<char>(code << bit_shift);
    }
    return packed;
}

static void test_packed_scanner_matches_kmer_scanner() {
    std::fprintf(stderr, "-- test_packed_scanner_matches_kmer_scanner\n");

    // Test with a clean sequence (no ambiguity)
    const char* seq = "ACGTACGTACGTACGTACGT";
    uint32_t len = 20;
    int k = 7;

    // Reference: KmerScanner results
    std::vector<std::pair<uint32_t, uint16_t>> ref_kmers;
    KmerScanner<uint16_t> ref_scanner(k);
    ref_scanner.scan(seq, len, [&](uint32_t pos, uint16_t kmer) {
        ref_kmers.push_back({pos, kmer});
    });

    // PackedKmerScanner results
    auto packed = encode_ncbi2na(seq, len);
    std::vector<AmbiguityEntry> empty_ambig;
    std::vector<std::pair<uint32_t, uint16_t>> packed_kmers;
    int ambig_calls = 0;

    PackedKmerScanner<uint16_t> packed_scanner(k);
    packed_scanner.scan(packed.data(), len, empty_ambig,
        [&](uint32_t pos, uint16_t kmer) {
            packed_kmers.push_back({pos, kmer});
        },
        [&](uint32_t, uint16_t, uint8_t, int) {
            ambig_calls++;
        });

    CHECK_EQ(ambig_calls, 0);
    CHECK_EQ(packed_kmers.size(), ref_kmers.size());
    for (size_t i = 0; i < ref_kmers.size() && i < packed_kmers.size(); i++) {
        CHECK_EQ(packed_kmers[i].first, ref_kmers[i].first);
        CHECK_EQ(packed_kmers[i].second, ref_kmers[i].second);
    }
}

static void test_packed_scanner_single_ambig() {
    std::fprintf(stderr, "-- test_packed_scanner_single_ambig\n");

    // Sequence: ACGTACG?ACGT (12 bases), ? at position 7 is R (A|G, ncbi4na=5)
    const char* raw_seq = "ACGTACGAACGT"; // placeholder A at pos 7
    uint32_t len = 12;
    int k = 5;

    auto packed = encode_ncbi2na(raw_seq, len);

    // Ambiguity entry: R at position 7
    std::vector<AmbiguityEntry> ambig = {{7, 1, 5}}; // pos=7, run=1, ncbi4na=5 (R=A|G)

    std::vector<std::pair<uint32_t, uint16_t>> normal_kmers;
    struct AmbigKmer {
        uint32_t pos;
        uint16_t base_kmer;
        uint8_t ncbi4na;
        int bit_offset;
    };
    std::vector<AmbigKmer> ambig_kmers;

    PackedKmerScanner<uint16_t> scanner(k);
    scanner.scan(packed.data(), len, ambig,
        [&](uint32_t pos, uint16_t kmer) {
            normal_kmers.push_back({pos, kmer});
        },
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            ambig_kmers.push_back({pos, base_kmer, ncbi4na, bit_offset});
        });

    // K-mers containing position 7 are at positions 3,4,5,6,7
    // All other k-mers should be normal
    // Positions 0,1,2 are normal (window [0,4], [1,5], [2,6] - don't include pos 7)
    // Position 3: window [3,7] includes pos 7 -> ambig
    // ...
    // Position 7: window [7,11] includes pos 7 -> ambig
    // Total: 3 normal + 5 ambig = 8 k-mers (len - k + 1 = 8)
    CHECK_EQ(normal_kmers.size() + ambig_kmers.size(), size_t(8));
    CHECK_EQ(normal_kmers.size(), size_t(3));
    CHECK_EQ(ambig_kmers.size(), size_t(5));

    // Verify ambig k-mers have ncbi4na = 5 (R)
    for (const auto& ak : ambig_kmers) {
        CHECK_EQ(ak.ncbi4na, uint8_t(5));
    }
}

static void test_packed_scanner_double_ambig_skip() {
    std::fprintf(stderr, "-- test_packed_scanner_double_ambig_skip\n");

    // Sequence: ACNNACGT (8 bases), N at positions 2 and 3
    const char* raw_seq = "ACAAACGT"; // placeholders at pos 2,3
    uint32_t len = 8;
    int k = 5;

    auto packed = encode_ncbi2na(raw_seq, len);

    // Two Ns at positions 2 and 3 (run of 2)
    std::vector<AmbiguityEntry> ambig = {{2, 2, 15}}; // pos=2, run=2, ncbi4na=15 (N)

    std::vector<std::pair<uint32_t, uint16_t>> normal_kmers;
    int ambig_calls = 0;

    PackedKmerScanner<uint16_t> scanner(k);
    scanner.scan(packed.data(), len, ambig,
        [&](uint32_t pos, uint16_t kmer) {
            normal_kmers.push_back({pos, kmer});
        },
        [&](uint32_t, uint16_t, uint8_t, int) {
            ambig_calls++;
        });

    // k=5, len=8, total possible k-mers = 4
    // pos 0: [0,4] includes pos 2,3 -> 2 ambig -> skip
    // pos 1: [1,5] includes pos 2,3 -> 2 ambig -> skip
    // pos 2: [2,6] includes pos 2,3 -> 2 ambig -> skip
    // pos 3: [3,7] includes pos 3 -> 1 ambig -> ambig callback
    CHECK_EQ(normal_kmers.size(), size_t(0));
    CHECK_EQ(ambig_calls, 1);
}

static void test_packed_scanner_ambig_expansion() {
    std::fprintf(stderr, "-- test_packed_scanner_ambig_expansion\n");

    // Verify that ambig_callback provides correct bit_offset for expansion.
    // Sequence: ACRTG (5 bases), R at position 2 (A|G, ncbi4na=5)
    const char* raw_seq = "ACATG"; // placeholder A at pos 2
    uint32_t len = 5;
    int k = 5;

    auto packed = encode_ncbi2na(raw_seq, len);
    std::vector<AmbiguityEntry> ambig = {{2, 1, 5}}; // R at pos 2

    uint32_t got_pos = UINT32_MAX;
    uint16_t got_base_kmer = 0;
    uint8_t got_ncbi4na = 0;
    int got_bit_offset = -1;

    PackedKmerScanner<uint16_t> scanner(k);
    scanner.scan(packed.data(), len, ambig,
        [](uint32_t, uint16_t) {},
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            got_pos = pos;
            got_base_kmer = base_kmer;
            got_ncbi4na = ncbi4na;
            got_bit_offset = bit_offset;
        });

    CHECK_EQ(got_pos, uint32_t(0));
    CHECK_EQ(got_ncbi4na, uint8_t(5)); // R
    // Position 2 in a 5-mer: rightmost base is pos 4, so
    // bases_from_right = 4 - 2 = 2, bit_offset = 4
    CHECK_EQ(got_bit_offset, 4);

    // Expand and verify we get both A and G variants
    std::set<uint16_t> expanded;
    uint16_t clear_mask = ~(uint16_t(0x03) << got_bit_offset);
    uint16_t cleared = got_base_kmer & clear_mask;
    for (uint8_t b = 0; b < 4; b++) {
        if (got_ncbi4na & (1u << b)) {
            expanded.insert(cleared | (uint16_t(b) << got_bit_offset));
        }
    }
    CHECK_EQ(expanded.size(), size_t(2)); // A and G

    // Verify expanded k-mers match "ACATG" and "ACGTG"
    uint16_t kmer_acatg = 0, kmer_acgtg = 0;
    KmerScanner<uint16_t> ref(k);
    ref.scan("ACATG", 5, [&](uint32_t, uint16_t km) { kmer_acatg = km; });
    ref.scan("ACGTG", 5, [&](uint32_t, uint16_t km) { kmer_acgtg = km; });
    CHECK(expanded.count(kmer_acatg) == 1);
    CHECK(expanded.count(kmer_acgtg) == 1);
}

// Test all IUPAC ambiguity codes: 2-base (M,R,W,S,Y,K), 3-base (V,H,D,B), 4-base (N)
static void test_packed_scanner_all_iupac_expansion() {
    std::fprintf(stderr, "-- test_packed_scanner_all_iupac_expansion\n");

    // For each ambiguity code, place it in the middle of a 5-mer "AC?GT"
    // and verify the expansion count matches the number of bases represented.

    // ncbi4na -> expected expansion count
    // 3=M(AC)=2, 5=R(AG)=2, 9=W(AT)=2, 6=S(CG)=2, 10=Y(CT)=2, 12=K(GT)=2
    // 7=V(ACG)=3, 11=H(ACT)=3, 13=D(AGT)=3, 14=B(CGT)=3
    // 15=N(ACGT)=4
    struct TestCase {
        uint8_t ncbi4na;
        int expected_count;
    };
    TestCase cases[] = {
        { 3, 2}, { 5, 2}, { 9, 2}, { 6, 2}, {10, 2}, {12, 2},  // 2-base
        { 7, 3}, {11, 3}, {13, 3}, {14, 3},                      // 3-base
        {15, 4},                                                   // 4-base
    };

    const char* raw_seq = "ACAGT"; // placeholder A at pos 2
    uint32_t len = 5;
    int k = 5;
    auto packed = encode_ncbi2na(raw_seq, len);

    for (const auto& tc : cases) {
        std::vector<AmbiguityEntry> ambig = {{2, 1, tc.ncbi4na}};

        int ambig_calls = 0;
        std::set<uint16_t> expanded_kmers;

        PackedKmerScanner<uint16_t> scanner(k);
        scanner.scan(packed.data(), len, ambig,
            [](uint32_t, uint16_t) {},
            [&](uint32_t, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
                ambig_calls++;
                uint16_t clear_mask = ~(uint16_t(0x03) << bit_offset);
                uint16_t cleared = base_kmer & clear_mask;
                for (uint8_t b = 0; b < 4; b++) {
                    if (ncbi4na & (1u << b)) {
                        expanded_kmers.insert(cleared | (uint16_t(b) << bit_offset));
                    }
                }
            });

        CHECK_EQ(ambig_calls, 1);
        CHECK_EQ(static_cast<int>(expanded_kmers.size()), tc.expected_count);

        // Verify each expanded k-mer matches the expected base substitution
        for (uint16_t ekm : expanded_kmers) {
            // Extract the base at bit_offset=4 (pos 2 in a 5-mer)
            uint8_t base_at_pos2 = (ekm >> 4) & 0x03;
            // This base must be one of the bases in ncbi4na
            CHECK((tc.ncbi4na & (1u << base_at_pos2)) != 0);
        }
    }
}

// Test run length > 1 ambiguous regions
static void test_packed_scanner_run_length() {
    std::fprintf(stderr, "-- test_packed_scanner_run_length\n");

    // Sequence: ACGTNNNNNACGT (13 bases), N run at positions 4-8 (run=5)
    const char* raw_seq = "ACGTAAAAACGT"; // placeholders
    // Actually 13 bases but our placeholder is 12. Let me fix:
    // "ACGT" + "AAAAA" + "ACGT" = 13 chars
    const char* raw_seq2 = "ACGTAAAAACGT";
    // Wait, that's 12 chars. Let me count: A-C-G-T-A-A-A-A-A-C-G-T = 12
    // I need 13: "ACGTAAAAACGTA" but that doesn't match...
    // Let me just use the correct length:
    // positions: 0=A, 1=C, 2=G, 3=T, 4-8=NNNNN, 9=A, 10=C, 11=G, 12=T
    const char* raw_data = "ACGTAAAAACGT"; // 12 chars for 13 bases? No...
    // "ACGT" (4) + "AAAAA" (5) + "ACGT" (4) = 13 chars
    const char* raw13 = "ACGTAAAAACGT";
    // A(0)C(1)G(2)T(3)A(4)A(5)A(6)A(7)A(8)C(9)G(10)T(11) = 12 chars. Hmm.
    // Let me use 12 bases instead: ACGTNNNNNACG, run at 4-8
    // Or more simply: use "ACGTAAAACGT" = 11 bases, run at 4-7 (run=4)
    // Actually let me simplify: 10 bases, run of 3 Ns

    // ACNNNACGTA (10 bases), N run at positions 2-4 (run=3)
    const char* seq10 = "ACAAAACGTA";
    uint32_t len = 10;
    int k = 5;

    auto packed = encode_ncbi2na(seq10, len);
    std::vector<AmbiguityEntry> ambig = {{2, 3, 15}}; // pos=2, run=3, ncbi4na=15 (N)

    std::vector<std::pair<uint32_t, uint16_t>> normal_kmers;
    int ambig_calls = 0;
    int skip_count = 0;

    PackedKmerScanner<uint16_t> scanner(k);
    scanner.scan(packed.data(), len, ambig,
        [&](uint32_t pos, uint16_t kmer) {
            normal_kmers.push_back({pos, kmer});
        },
        [&](uint32_t, uint16_t, uint8_t, int) {
            ambig_calls++;
        });

    // k=5, len=10, possible k-mers = 6 (positions 0-5)
    // pos 0: [0,4] includes 2,3,4 -> 3 ambig -> skip
    // pos 1: [1,5] includes 2,3,4 -> 3 ambig -> skip
    // pos 2: [2,6] includes 2,3,4 -> 3 ambig -> skip
    // pos 3: [3,7] includes 3,4 -> 2 ambig -> skip
    // pos 4: [4,8] includes 4 -> 1 ambig -> ambig callback
    // pos 5: [5,9] no ambig -> normal
    int total = static_cast<int>(normal_kmers.size()) + ambig_calls;
    CHECK_EQ(total, 2);  // 1 ambig + 1 normal
    CHECK_EQ(static_cast<int>(normal_kmers.size()), 1);
    CHECK_EQ(ambig_calls, 1);
}

// Test sequence length not a multiple of 4 (partial last byte)
static void test_packed_scanner_partial_byte() {
    std::fprintf(stderr, "-- test_packed_scanner_partial_byte\n");

    // Sequence lengths 5,6,7 (not multiples of 4) with no ambiguity
    for (uint32_t len = 5; len <= 7; len++) {
        const char* seqs[] = {"ACGTA", "ACGTAC", "ACGTACG"};
        const char* seq = seqs[len - 5];
        int k = 5;

        auto packed = encode_ncbi2na(seq, len);
        std::vector<AmbiguityEntry> empty_ambig;

        // Reference
        std::vector<std::pair<uint32_t, uint16_t>> ref_kmers;
        KmerScanner<uint16_t> rs(k);
        rs.scan(seq, len, [&](uint32_t pos, uint16_t kmer) {
            ref_kmers.push_back({pos, kmer});
        });

        // PackedKmerScanner
        std::vector<std::pair<uint32_t, uint16_t>> packed_kmers;
        PackedKmerScanner<uint16_t> ps(k);
        ps.scan(packed.data(), len, empty_ambig,
            [&](uint32_t pos, uint16_t kmer) {
                packed_kmers.push_back({pos, kmer});
            },
            [](uint32_t, uint16_t, uint8_t, int) {});

        CHECK_EQ(packed_kmers.size(), ref_kmers.size());
        for (size_t i = 0; i < ref_kmers.size() && i < packed_kmers.size(); i++) {
            CHECK_EQ(packed_kmers[i].first, ref_kmers[i].first);
            CHECK_EQ(packed_kmers[i].second, ref_kmers[i].second);
        }
    }
}

// Test long ambiguous run (16 consecutive Ns)
static void test_packed_scanner_long_ambig_run() {
    std::fprintf(stderr, "-- test_packed_scanner_long_ambig_run\n");

    // "ACGT" + 16*N + "ACGT" = 24 bases
    // In ncbi2na, N positions have placeholder bases.
    const char* raw_seq = "ACGTAAAAAAAAAAAAAAAACGT";
    // A(0)C(1)G(2)T(3) A*16(4-19) A(20)C(21)G(22)T(23) = 24
    // Wait, that's "ACGT" + 16*A + "ACGT" = 4+16+4 = 24. Correct.
    uint32_t len = 24;
    int k = 5;

    auto packed = encode_ncbi2na(raw_seq, len);
    // N run at positions 4-19 (run=16)
    std::vector<AmbiguityEntry> ambig = {{4, 16, 15}}; // ncbi4na=15 (N)

    std::vector<std::pair<uint32_t, uint16_t>> normal_kmers;
    int ambig_calls = 0;

    PackedKmerScanner<uint16_t> scanner(k);
    scanner.scan(packed.data(), len, ambig,
        [&](uint32_t pos, uint16_t kmer) {
            normal_kmers.push_back({pos, kmer});
        },
        [&](uint32_t, uint16_t, uint8_t, int) {
            ambig_calls++;
        });

    // Window positions 0..19 (20 total k-mers)
    // pos 0: [0,4] includes pos 4 -> 1 ambig
    // pos 1-15: [1,5]..[15,19] includes 2+ ambig positions -> skip
    // pos 16: [16,20] includes pos 16-19 (4 ambig) -> skip
    //   Actually pos 16: window [16,20], ambig at 16,17,18,19 -> 4 -> skip
    // pos 17: [17,21] ambig at 17,18,19 -> 3 -> skip
    // pos 18: [18,22] ambig at 18,19 -> 2 -> skip
    // pos 19: [19,23] ambig at 19 -> 1 ambig
    // So: normal at none from the Ns area, ambig_calls at pos 0 and 19
    // Clean k-mers: none that have 0 ambig in window
    // Wait, we also need to check if there are k-mers entirely before pos 4:
    // pos 0: [0,4] -> pos 4 is ambig -> 1 ambig (ambig_call)
    // No k-mer starts before 0 for k=5. So no fully clean k-mers before the N run.
    // After the N run: pos 19: [19,23] includes pos 19 (ambig) -> 1 ambig
    // No clean k-mers after either.

    // Let's verify: only 2 ambig callbacks, 0 normal
    CHECK_EQ(static_cast<int>(normal_kmers.size()), 0);
    CHECK_EQ(ambig_calls, 2);
}

// ==================== Integration tests with BLAST DB ====================

static void test_blastdb_raw_sequence() {
    std::fprintf(stderr, "-- test_blastdb_raw_sequence\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    // Find FJ876973.1 and verify raw sequence access
    uint32_t oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(oid != UINT32_MAX);

    auto raw = db.get_raw_sequence(oid);
    CHECK(raw.ncbi2na_data != nullptr);
    CHECK_EQ(raw.seq_length, db.seq_length(oid));
    CHECK(raw.ncbi2na_bytes > 0);

    // Verify we can decode bases correctly
    for (uint32_t i = 0; i < raw.seq_length; i++) {
        uint8_t code = ncbi2na_base_at(raw.ncbi2na_data, i);
        CHECK(code < 4);
    }

    // Verify get_sequence returns valid bases
    std::string seq = db.get_sequence(oid);
    CHECK_EQ(seq.size(), static_cast<size_t>(raw.seq_length));
    for (char c : seq) {
        CHECK(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
    }

    db.ret_raw_sequence(raw);
}

static void test_blastdb_raw_sequence_with_ambig() {
    std::fprintf(stderr, "-- test_blastdb_raw_sequence_with_ambig\n");

    // Open the ambig DB where FJ876973.1 has R injected at position 100
    BlastDbReader db;
    CHECK(db.open(g_ambigdb_path));

    uint32_t oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(oid != UINT32_MAX);

    auto raw = db.get_raw_sequence(oid);
    CHECK(raw.ncbi2na_data != nullptr);
    CHECK_EQ(raw.seq_length, db.seq_length(oid));

    // Check ambiguity data exists (R was injected at position 100)
    auto ambig = AmbiguityParser::parse(raw.ambig_data, raw.ambig_bytes);
    CHECK(ambig.size() > 0);

    // Verify get_sequence shows R at position 100
    std::string seq = db.get_sequence(oid);
    CHECK(seq.size() > 100);
    CHECK(seq[100] == 'R');

    db.ret_raw_sequence(raw);
}

static void test_packed_scanner_matches_kmer_scanner_on_blastdb() {
    std::fprintf(stderr, "-- test_packed_scanner_matches_kmer_scanner_on_blastdb\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    int k = 7;

    // Compare on the first 20 non-ambiguous sequences to keep runtime bounded
    uint32_t limit = std::min(db.num_sequences(), uint32_t(20));
    for (uint32_t oid = 0; oid < limit; oid++) {
        std::string seq = db.get_sequence(oid);
        bool has_ambig = false;
        for (char c : seq) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                has_ambig = true;
                break;
            }
        }
        if (has_ambig) continue;

        auto raw = db.get_raw_sequence(oid);
        auto ambig = AmbiguityParser::parse(raw.ambig_data, raw.ambig_bytes);

        std::vector<std::pair<uint32_t, uint16_t>> packed_kmers;
        PackedKmerScanner<uint16_t> ps(k);
        ps.scan(raw.ncbi2na_data, raw.seq_length, ambig,
            [&](uint32_t pos, uint16_t kmer) {
                packed_kmers.push_back({pos, kmer});
            },
            [](uint32_t, uint16_t, uint8_t, int) {});

        std::vector<std::pair<uint32_t, uint16_t>> ref_kmers;
        KmerScanner<uint16_t> rs(k);
        rs.scan(seq.data(), seq.size(),
            [&](uint32_t pos, uint16_t kmer) {
                ref_kmers.push_back({pos, kmer});
            });

        CHECK_EQ(packed_kmers.size(), ref_kmers.size());
        for (size_t i = 0; i < ref_kmers.size() && i < packed_kmers.size(); i++) {
            CHECK_EQ(packed_kmers[i].first, ref_kmers[i].first);
            CHECK_EQ(packed_kmers[i].second, ref_kmers[i].second);
        }

        db.ret_raw_sequence(raw);
    }
}

static void test_ambig_db_index_build() {
    std::fprintf(stderr, "-- test_ambig_db_index_build\n");

    BlastDbReader db;
    CHECK(db.open(g_ambigdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 5;

    std::string prefix = g_output_dir + "/ambig.00.05mer";
    CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "ambig", logger));

    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
    CHECK_EQ(kix.k(), 5);
    CHECK(kix.total_postings() > 0);

    // Verify counts sum matches total_postings
    uint64_t sum = 0;
    for (uint64_t i = 0; i < kix.table_size(); i++) {
        sum += kix.counts()[i];
    }
    CHECK_EQ(sum, kix.total_postings());

    kix.close();
}

static void test_ambig_expansion_in_index() {
    std::fprintf(stderr, "-- test_ambig_expansion_in_index\n");

    // In ssu_ambigdb, FJ876973.1 has R (A|G) injected at position 100.
    // The 5-mer spanning position 100 should be expanded into A and G variants.

    BlastDbReader db;
    CHECK(db.open(g_ambigdb_path));

    uint32_t oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(oid != UINT32_MAX);

    // Read the sequence to get the bases around position 100
    std::string seq = db.get_sequence(oid);
    CHECK(seq.size() > 104);
    CHECK(seq[100] == 'R');

    // Extract the 5-mer at position 96..100: 4 clean bases + R
    // Expand: replace R(pos 100) with A and G
    std::string kmer_base = seq.substr(96, 5);  // has 'R' at index 4
    std::string kmer_a = kmer_base;
    std::string kmer_g = kmer_base;
    kmer_a[4] = 'A';
    kmer_g[4] = 'G';

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 5;

    std::string prefix = g_output_dir + "/ambig_exp.00.05mer";
    CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "ambig", logger));

    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));

    // Encode both expansion variants
    uint16_t kval_a = 0, kval_g = 0;
    KmerScanner<uint16_t> ref(5);
    ref.scan(kmer_a.data(), 5, [&](uint32_t, uint16_t km) { kval_a = km; });
    ref.scan(kmer_g.data(), 5, [&](uint32_t, uint16_t km) { kval_g = km; });

    // Both should have non-zero counts in the index
    CHECK(kix.counts()[kval_a] > 0);
    CHECK(kix.counts()[kval_g] > 0);

    kix.close();
}

// Decode ID posting list
static std::vector<uint32_t> decode_id_postings(
    const uint8_t* data, uint64_t offset, uint32_t count) {
    std::vector<uint32_t> result;
    if (count == 0) return result;
    result.reserve(count);
    const uint8_t* p = data + offset;
    uint32_t prev_id = 0;
    for (uint32_t i = 0; i < count; i++) {
        uint32_t delta;
        p += varint_decode(p, delta);
        uint32_t id = (i == 0) ? delta : prev_id + delta;
        result.push_back(id);
        prev_id = id;
    }
    return result;
}

static void test_ssu_db_kmer_check() {
    std::fprintf(stderr, "-- test_ssu_db_kmer_check\n");

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    // Find FJ876973.1 and extract first 7bp
    uint32_t target_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(target_oid != UINT32_MAX);

    std::string full_seq = db.get_sequence(target_oid);
    CHECK(full_seq.size() >= 7);
    std::string first7 = full_seq.substr(0, 7);

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 7;

    std::string prefix = g_output_dir + "/compat.00.07mer";
    CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger));

    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
    CHECK_EQ(kix.k(), 7);
    CHECK(kix.num_sequences() > 0);

    uint16_t target_kmer = 0;
    KmerScanner<uint16_t> ref(7);
    ref.scan(first7.data(), 7, [&](uint32_t, uint16_t km) { target_kmer = km; });

    CHECK(kix.counts()[target_kmer] > 0);

    auto ids = decode_id_postings(
        kix.posting_data(), kix.offsets()[target_kmer],
        kix.counts()[target_kmer]);
    bool has_target = false;
    for (uint32_t id : ids) {
        if (id == target_oid) has_target = true;
    }
    CHECK(has_target);

    kix.close();
}

// Test index build with odd-length sequence (DQ235612.1 = 1809bp, odd)
static void test_ambig_db_odd_length() {
    std::fprintf(stderr, "-- test_ambig_db_odd_length\n");

    BlastDbReader db;
    CHECK(db.open(g_ambigdb_path));

    // DQ235612.1 is 1809bp (odd length)
    uint32_t oid = find_oid_by_accession(db, ACC_DQ);
    CHECK(oid != UINT32_MAX);

    uint32_t len = db.seq_length(oid);
    CHECK(len % 2 == 1); // Verify odd length

    // Build index and verify it succeeds
    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 5;

    std::string prefix = g_output_dir + "/odd.00.05mer";
    CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "ambig", logger));

    KixReader kix;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kix.total_postings() > 0);

    uint64_t sum = 0;
    for (uint64_t i = 0; i < kix.table_size(); i++) {
        sum += kix.counts()[i];
    }
    CHECK_EQ(sum, kix.total_postings());

    kix.close();
}

int main(int argc, char* argv[]) {
    check_ssu_available();
    check_derived_data_ready();

    g_testdb_path = ssu_db_prefix();
    g_ambigdb_path = ambig_db_prefix();
    g_output_dir = "/tmp/ikafssn_packed_scanner_test";
    std::filesystem::create_directories(g_output_dir);

    // Unit tests: AmbiguityParser
    test_ambiguity_parser_empty();
    test_ambiguity_parser_old_format();
    test_ambiguity_parser_new_format();
    test_ambiguity_parser_sort();

    // Unit tests: PackedKmerScanner
    test_packed_scanner_matches_kmer_scanner();
    test_packed_scanner_single_ambig();
    test_packed_scanner_double_ambig_skip();
    test_packed_scanner_ambig_expansion();
    test_packed_scanner_all_iupac_expansion();
    test_packed_scanner_run_length();
    test_packed_scanner_partial_byte();
    test_packed_scanner_long_ambig_run();

    // Integration tests
    test_blastdb_raw_sequence();
    test_blastdb_raw_sequence_with_ambig();
    test_packed_scanner_matches_kmer_scanner_on_blastdb();
    test_ambig_db_index_build();
    test_ambig_expansion_in_index();
    test_ssu_db_kmer_check();
    test_ambig_db_odd_length();

    // Clean up
    std::filesystem::remove_all(g_output_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
