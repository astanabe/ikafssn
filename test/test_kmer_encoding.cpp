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
    // k=13 with uint32_t (max supported)
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
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
