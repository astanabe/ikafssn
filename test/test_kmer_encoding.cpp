#include "test_util.hpp"
#include "core/kmer_encoding.hpp"
#include <vector>
#include <string>

using namespace ikafssn;

static void test_base_encoding() {
    CHECK_EQ(encode_base('A'), 0);
    CHECK_EQ(encode_base('C'), 1);
    CHECK_EQ(encode_base('G'), 2);
    CHECK_EQ(encode_base('T'), 3);
    CHECK_EQ(encode_base('a'), 0);
    CHECK_EQ(encode_base('c'), 1);
    CHECK_EQ(encode_base('g'), 2);
    CHECK_EQ(encode_base('t'), 3);
    CHECK_EQ(encode_base('N'), BASE_ENCODE_INVALID);
    CHECK_EQ(encode_base('n'), BASE_ENCODE_INVALID);
    CHECK_EQ(encode_base('X'), BASE_ENCODE_INVALID);
}

static void test_known_kmer() {
    // "ACGT" (k=4): A=00, C=01, G=10, T=11 -> 0b00011011 = 0x1B = 27
    uint16_t kmer = 0;
    const char* seq = "ACGT";
    for (int i = 0; i < 4; i++) {
        kmer = ((kmer << 2) | encode_base(seq[i])) & kmer_mask<uint16_t>(4);
    }
    CHECK_EQ(kmer, 0x1B);
}

static void test_revcomp_involution_u16() {
    // revcomp(revcomp(x)) == x for all valid k-mers
    for (int k = MIN_K; k <= 8; k++) {
        uint16_t mask = kmer_mask<uint16_t>(k);
        // Test a sample of values
        for (uint32_t v = 0; v < table_size(k) && v < 1024; v++) {
            uint16_t kmer = static_cast<uint16_t>(v);
            uint16_t rc1 = kmer_revcomp<uint16_t>(kmer, k);
            uint16_t rc2 = kmer_revcomp<uint16_t>(rc1, k);
            CHECK_EQ(rc2, kmer);
            CHECK((rc1 & mask) == rc1); // no bits outside mask
        }
    }
}

static void test_revcomp_involution_u32() {
    for (int k = 9; k <= MAX_K; k++) {
        uint32_t mask = kmer_mask<uint32_t>(k);
        // Test a sample of values
        for (uint32_t v = 0; v < 1024; v++) {
            uint32_t rc1 = kmer_revcomp<uint32_t>(v, k);
            uint32_t rc2 = kmer_revcomp<uint32_t>(rc1, k);
            CHECK_EQ(rc2, v);
            CHECK((rc1 & mask) == rc1);
        }
        // Test high values near max
        uint64_t ts = table_size(k);
        for (uint64_t v = ts - 1024; v < ts; v++) {
            uint32_t kmer = static_cast<uint32_t>(v);
            uint32_t rc1 = kmer_revcomp<uint32_t>(kmer, k);
            uint32_t rc2 = kmer_revcomp<uint32_t>(rc1, k);
            CHECK_EQ(rc2, kmer);
        }
    }
}

static void test_revcomp_known() {
    // "ACGT" -> revcomp = "ACGT" (palindrome for k=4)
    // A=00,C=01,G=10,T=11 -> 0b00011011
    // revcomp: complement = ~0b00011011 = 0b11100100, reverse 2-bit pairs:
    // T=11,G=10,C=01,A=00 -> complement -> A=00,C=01,G=10,T=11 = ACGT
    uint16_t acgt = 0x1B;
    CHECK_EQ(kmer_revcomp<uint16_t>(acgt, 4), acgt);

    // "AAAA" (k=4) -> 0b00000000 = 0
    // revcomp = "TTTT" = 0b11111111 = 255
    CHECK_EQ(kmer_revcomp<uint16_t>(uint16_t(0), 4), uint16_t(0xFF));

    // "TTTT" -> revcomp = "AAAA"
    CHECK_EQ(kmer_revcomp<uint16_t>(uint16_t(0xFF), 4), uint16_t(0));
}

static void test_scanner_basic() {
    // Scan "ACGTACGT" with k=5 -> should produce 4 k-mers
    std::vector<std::pair<uint32_t, uint16_t>> results;
    KmerScanner<uint16_t> scanner(5);
    std::string seq = "ACGTACGT";
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint16_t kmer) {
        results.push_back({pos, kmer});
    });
    CHECK_EQ(results.size(), 4u); // 8 - 5 + 1 = 4
    // Positions should be 0, 1, 2, 3
    for (size_t i = 0; i < results.size(); i++) {
        CHECK_EQ(results[i].first, static_cast<uint32_t>(i));
    }
}

static void test_scanner_with_n() {
    // "ACNGTACGT" with k=5: N at pos 2 should skip k-1=4 bases after it
    std::vector<std::pair<uint32_t, uint16_t>> results;
    KmerScanner<uint16_t> scanner(5);
    std::string seq = "ACNGTACGT";
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint16_t kmer) {
        results.push_back({pos, kmer});
    });
    // After N at pos 2, need 4 more valid bases. Next valid bases: G(3),T(4),A(5),C(6),G(7),T(8)
    // First valid k-mer after N starts at pos 3+4-4 = pos 3? Let's think carefully:
    // pos 0: A (n_count was 4, now 3)
    // pos 1: C (n_count 2)
    // pos 2: N -> n_count = 4
    // pos 3: G (n_count 3)
    // pos 4: T (n_count 2)
    // pos 5: A (n_count 1)
    // pos 6: C (n_count 0) -> still decrementing
    // pos 7: G (n_count was 0, but we decremented to -1? No, n_count > 0 check.
    // Wait: at pos 6, n_count was 1, enc valid, kmer updated, n_count > 0 -> n_count-- (now 0), continue
    // At pos 7: G, enc valid, kmer updated, n_count == 0, so callback at pos 7-5+1 = 3
    // Wait, that's not right. Let me re-read the scanner:
    // n_count starts at k-1 = 4
    // pos 0: A, enc=0, kmer updated, n_count=4>0, n_count-- -> 3
    // pos 1: C, enc=1, kmer updated, n_count=3>0, n_count-- -> 2
    // pos 2: N, enc=0xFF, n_count=4, kmer=0
    // pos 3: G, enc=2, kmer updated, n_count=4>0, n_count-- -> 3
    // pos 4: T, enc=3, kmer updated, n_count=3>0, n_count-- -> 2
    // pos 5: A, enc=0, kmer updated, n_count=2>0, n_count-- -> 1
    // pos 6: C, enc=1, kmer updated, n_count=1>0, n_count-- -> 0
    // pos 7: G, enc=2, kmer updated, n_count=0, callback(7-5+1=3, kmer)
    // pos 8: T, enc=3, kmer updated, n_count=0, callback(8-5+1=4, kmer)
    // So 2 k-mers at positions 3 and 4
    CHECK_EQ(results.size(), 2u);
    CHECK_EQ(results[0].first, 3u);
    CHECK_EQ(results[1].first, 4u);
}

static void test_scanner_k8_boundary() {
    // k=8 with uint16_t (max for uint16_t)
    std::string seq = "ACGTACGTACGT"; // len=12, produces 12-8+1=5 k-mers
    std::vector<std::pair<uint32_t, uint16_t>> results;
    KmerScanner<uint16_t> scanner(8);
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint16_t kmer) {
        results.push_back({pos, kmer});
    });
    CHECK_EQ(results.size(), 5u);
}

static void test_scanner_k9_u32() {
    // k=9 with uint32_t
    std::string seq = "ACGTACGTACGT"; // len=12, produces 12-9+1=4 k-mers
    std::vector<std::pair<uint32_t, uint32_t>> results;
    KmerScanner<uint32_t> scanner(9);
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint32_t kmer) {
        results.push_back({pos, kmer});
    });
    CHECK_EQ(results.size(), 4u);
}

static void test_scanner_k13_u32() {
    // k=13 with uint32_t
    std::string seq = "ACGTACGTACGTACGT"; // len=16, produces 16-13+1=4 k-mers
    std::vector<std::pair<uint32_t, uint32_t>> results;
    KmerScanner<uint32_t> scanner(13);
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint32_t kmer) {
        results.push_back({pos, kmer});
    });
    CHECK_EQ(results.size(), 4u);
    // Check all k-mers are within valid range
    uint32_t mask = kmer_mask<uint32_t>(13);
    for (auto& [pos, kmer] : results) {
        CHECK((kmer & mask) == kmer);
    }
}

static void test_scanner_k16_u32() {
    // k=16 with uint32_t (max supported, uses all 32 bits)
    std::string seq = "ACGTACGTACGTACGTACGT"; // len=20, produces 20-16+1=5 k-mers
    std::vector<std::pair<uint32_t, uint32_t>> results;
    KmerScanner<uint32_t> scanner(16);
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint32_t kmer) {
        results.push_back({pos, kmer});
    });
    CHECK_EQ(results.size(), 5u);
    // k=16 uses full uint32_t: mask = 0xFFFFFFFF
    uint32_t mask = kmer_mask<uint32_t>(16);
    CHECK_EQ(mask, uint32_t(0xFFFFFFFF));
    for (auto& [pos, kmer] : results) {
        CHECK((kmer & mask) == kmer);
    }
    // Verify revcomp involution at k=16
    for (auto& [pos, kmer] : results) {
        uint32_t rc1 = kmer_revcomp<uint32_t>(kmer, 16);
        uint32_t rc2 = kmer_revcomp<uint32_t>(rc1, 16);
        CHECK_EQ(rc2, kmer);
    }
}

static void test_contains_degenerate_base() {
    // ACGT only should return false
    CHECK(!contains_degenerate_base("ACGT"));
    CHECK(!contains_degenerate_base("acgt"));
    CHECK(!contains_degenerate_base("ACGTACGTACGT"));
    CHECK(!contains_degenerate_base(""));

    // Each IUPAC ambiguity code should return true
    CHECK(contains_degenerate_base("ACGTR"));  // R
    CHECK(contains_degenerate_base("ACGTY"));  // Y
    CHECK(contains_degenerate_base("ACGTS"));  // S
    CHECK(contains_degenerate_base("ACGTW"));  // W
    CHECK(contains_degenerate_base("ACGTK"));  // K
    CHECK(contains_degenerate_base("ACGTM"));  // M
    CHECK(contains_degenerate_base("ACGTB"));  // B
    CHECK(contains_degenerate_base("ACGTD"));  // D
    CHECK(contains_degenerate_base("ACGTH"));  // H
    CHECK(contains_degenerate_base("ACGTV"));  // V
    CHECK(contains_degenerate_base("ACGTN"));  // N

    // Lowercase
    CHECK(contains_degenerate_base("acgtr"));
    CHECK(contains_degenerate_base("acgty"));
    CHECK(contains_degenerate_base("acgts"));
    CHECK(contains_degenerate_base("acgtw"));
    CHECK(contains_degenerate_base("acgtk"));
    CHECK(contains_degenerate_base("acgtm"));
    CHECK(contains_degenerate_base("acgtb"));
    CHECK(contains_degenerate_base("acgtd"));
    CHECK(contains_degenerate_base("acgth"));
    CHECK(contains_degenerate_base("acgtv"));
    CHECK(contains_degenerate_base("acgtn"));

    // Single degenerate character
    CHECK(contains_degenerate_base("N"));
    CHECK(contains_degenerate_base("R"));
}

static void test_degenerate_ncbi4na() {
    // All 11 IUPAC degenerate codes
    CHECK_EQ(degenerate_ncbi4na('R'), 0x05u); // A|G
    CHECK_EQ(degenerate_ncbi4na('Y'), 0x0Au); // C|T
    CHECK_EQ(degenerate_ncbi4na('S'), 0x06u); // G|C
    CHECK_EQ(degenerate_ncbi4na('W'), 0x09u); // A|T
    CHECK_EQ(degenerate_ncbi4na('K'), 0x0Cu); // G|T
    CHECK_EQ(degenerate_ncbi4na('M'), 0x03u); // A|C
    CHECK_EQ(degenerate_ncbi4na('B'), 0x0Eu); // C|G|T
    CHECK_EQ(degenerate_ncbi4na('D'), 0x0Du); // A|G|T
    CHECK_EQ(degenerate_ncbi4na('H'), 0x0Bu); // A|C|T
    CHECK_EQ(degenerate_ncbi4na('V'), 0x07u); // A|C|G
    CHECK_EQ(degenerate_ncbi4na('N'), 0x0Fu); // A|C|G|T

    // Lowercase
    CHECK_EQ(degenerate_ncbi4na('r'), 0x05u);
    CHECK_EQ(degenerate_ncbi4na('y'), 0x0Au);
    CHECK_EQ(degenerate_ncbi4na('n'), 0x0Fu);

    // Normal bases return 0 (not degenerate)
    CHECK_EQ(degenerate_ncbi4na('A'), 0u);
    CHECK_EQ(degenerate_ncbi4na('C'), 0u);
    CHECK_EQ(degenerate_ncbi4na('G'), 0u);
    CHECK_EQ(degenerate_ncbi4na('T'), 0u);
    CHECK_EQ(degenerate_ncbi4na('a'), 0u);
}

static void test_scan_ambig_no_degen() {
    // Pure ACGT sequence should produce same results as scan()
    std::string seq = "ACGTACGT";
    int k = 5;

    std::vector<std::pair<uint32_t, uint16_t>> scan_results;
    KmerScanner<uint16_t> scanner(k);
    scanner.scan(seq.c_str(), seq.size(), [&](uint32_t pos, uint16_t kmer) {
        scan_results.push_back({pos, kmer});
    });

    std::vector<std::pair<uint32_t, uint16_t>> ambig_results;
    int ambig_count = 0;
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) {
            ambig_results.push_back({pos, kmer});
        },
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            ambig_count++;
        });

    CHECK_EQ(ambig_count, 0);
    CHECK_EQ(ambig_results.size(), scan_results.size());
    for (size_t i = 0; i < scan_results.size(); i++) {
        CHECK_EQ(ambig_results[i].first, scan_results[i].first);
        CHECK_EQ(ambig_results[i].second, scan_results[i].second);
    }
}

static void test_scan_ambig_single_R() {
    // "ACGTR" k=5 -> single degenerate base R at position 4
    std::string seq = "ACGTR";
    int k = 5;

    int normal_count = 0;
    int ambig_count = 0;
    uint8_t got_ncbi4na = 0;

    KmerScanner<uint16_t> scanner(k);
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) { normal_count++; },
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            ambig_count++;
            got_ncbi4na = ncbi4na;
            CHECK_EQ(pos, 0u); // k-mer starts at position 0
            CHECK_EQ(bit_offset, 0); // R is at the last (rightmost) position
        });

    CHECK_EQ(normal_count, 0);
    CHECK_EQ(ambig_count, 1);
    CHECK_EQ(got_ncbi4na, 0x05u); // R = A|G
}

static void test_scan_ambig_two_degen_skip() {
    // "ACRSW" k=5 -> 2 degenerate bases in the window -> skip
    std::string seq = "ACRSW";
    int k = 5;

    int normal_count = 0;
    int ambig_count = 0;

    KmerScanner<uint16_t> scanner(k);
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) { normal_count++; },
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            ambig_count++;
        });

    CHECK_EQ(normal_count, 0);
    CHECK_EQ(ambig_count, 0);
}

static void test_scan_ambig_expand_correct_kmers() {
    // "ACGTR" k=5: R -> A,G
    // ACGTA: A=00,C=01,G=10,T=11,A=00 -> 0b0001101100 = 0x6C
    // ACGTG: A=00,C=01,G=10,T=11,G=10 -> 0b0001101110 = 0x6E
    std::string seq = "ACGTR";
    int k = 5;

    std::vector<uint16_t> expanded;
    KmerScanner<uint16_t> scanner(k);
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) {},
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            expand_ambig_kmer<uint16_t>(base_kmer, ncbi4na, bit_offset,
                [&](uint16_t exp) { expanded.push_back(exp); });
        });

    CHECK_EQ(expanded.size(), 2u);
    // R = A|G -> bit0=A, bit2=G
    // A=00 at bit_offset=0: 0x6C (ACGTA)
    // G=10 at bit_offset=0: 0x6E (ACGTG)
    CHECK_EQ(expanded[0], uint16_t(0x6C));
    CHECK_EQ(expanded[1], uint16_t(0x6E));
}

static void test_scan_ambig_N_expansion() {
    // "ACGTN" k=5 -> N expands to 4 variants (A,C,G,T)
    std::string seq = "ACGTN";
    int k = 5;

    std::vector<uint16_t> expanded;
    KmerScanner<uint16_t> scanner(k);
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) {},
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            CHECK_EQ(ncbi4na, 0x0Fu); // N = A|C|G|T
            expand_ambig_kmer<uint16_t>(base_kmer, ncbi4na, bit_offset,
                [&](uint16_t exp) { expanded.push_back(exp); });
        });

    CHECK_EQ(expanded.size(), 4u);
}

static void test_scan_ambig_degen_at_start() {
    // "RACGT" k=5 -> R at start (position 0 of the window)
    // bit_offset should be (k-1)*2 = 8 for the first base in the window
    std::string seq = "RACGT";
    int k = 5;

    int ambig_count = 0;
    KmerScanner<uint16_t> scanner(k);
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) {},
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            ambig_count++;
            CHECK_EQ(pos, 0u);
            CHECK_EQ(bit_offset, 8); // (5-1-0)*2 = 8... but depends on slot mapping
        });

    CHECK_EQ(ambig_count, 1);
}

static void test_expand_ambig_kmer_shared() {
    // Verify the shared expand_ambig_kmer produces correct results
    // base_kmer: 5-mer with placeholder A (00) at bit_offset=4 (bits 5-4)
    // Encoding: bits 9-8=A(00), 7-6=C(01), 5-4=A(00, placeholder), 3-2=G(10), 1-0=T(11)
    // = 0b00_01_00_10_11 = 0x4B = 75
    uint16_t base_kmer = 0x4B;
    uint8_t ncbi4na = 0x05; // R = A|G
    int bit_offset = 4;

    std::vector<uint16_t> results;
    expand_ambig_kmer<uint16_t>(base_kmer, ncbi4na, bit_offset,
        [&](uint16_t exp) { results.push_back(exp); });

    CHECK_EQ(results.size(), 2u);
    // A(00) at bit_offset=4: cleared | (0<<4) = 0x4B (placeholder was already A)
    CHECK_EQ(results[0], uint16_t(0x4B));
    // G(10) at bit_offset=4: cleared | (2<<4) = 0x4B & ~0x30 | 0x20 = 0x6B
    CHECK_EQ(results[1], uint16_t(0x6B));
}

static void test_scan_ambig_sliding_window() {
    // Test that degenerate base correctly exits the window as it slides
    // "RACGTACGT" k=5:
    // pos 0: RACGT -> 1 degen (R) -> ambig_callback
    // pos 1: ACGTA -> 0 degen -> normal callback
    // pos 2: CGTAC -> 0 degen -> normal callback
    // pos 3: GTACG -> 0 degen -> normal callback
    // pos 4: TACGT -> 0 degen -> normal callback
    std::string seq = "RACGTACGT";
    int k = 5;

    int normal_count = 0;
    int ambig_count = 0;

    KmerScanner<uint16_t> scanner(k);
    scanner.scan_ambig(seq.c_str(), seq.size(),
        [&](uint32_t pos, uint16_t kmer) { normal_count++; },
        [&](uint32_t pos, uint16_t base_kmer, uint8_t ncbi4na, int bit_offset) {
            ambig_count++;
        });

    CHECK_EQ(ambig_count, 1); // only pos 0
    CHECK_EQ(normal_count, 4); // pos 1,2,3,4
}

int main() {
    test_base_encoding();
    test_known_kmer();
    test_revcomp_involution_u16();
    test_revcomp_involution_u32();
    test_revcomp_known();
    test_scanner_basic();
    test_scanner_with_n();
    test_scanner_k8_boundary();
    test_scanner_k9_u32();
    test_scanner_k13_u32();
    test_scanner_k16_u32();
    test_contains_degenerate_base();
    test_degenerate_ncbi4na();
    test_scan_ambig_no_degen();
    test_scan_ambig_single_R();
    test_scan_ambig_two_degen_skip();
    test_scan_ambig_expand_correct_kmers();
    test_scan_ambig_N_expansion();
    test_scan_ambig_degen_at_start();
    test_expand_ambig_kmer_shared();
    test_scan_ambig_sliding_window();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
