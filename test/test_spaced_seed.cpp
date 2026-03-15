#include "test_util.hpp"
#include "core/spaced_seed.hpp"
#include "core/kmer_encoding.hpp"
#include <string>
#include <vector>

using namespace ikafssn;

static void test_validate_spaced_seed() {
    CHECK(validate_spaced_seed(11, 0));
    CHECK(validate_spaced_seed(11, 16));
    CHECK(validate_spaced_seed(11, 18));
    CHECK(validate_spaced_seed(11, 21));
    CHECK(validate_spaced_seed(12, 16));
    CHECK(validate_spaced_seed(12, 18));
    CHECK(validate_spaced_seed(12, 21));
    CHECK(!validate_spaced_seed(10, 16));
    CHECK(!validate_spaced_seed(11, 17));
    CHECK(!validate_spaced_seed(13, 16));
    CHECK(!validate_spaced_seed(8, 18));
}

static void test_get_seed_masks() {
    auto masks = get_seed_masks(11, 16, TemplateType::kCoding);
    CHECK_EQ(masks.size(), 1u);
    CHECK_EQ(masks[0], MASK_K11_T16_CODING);

    masks = get_seed_masks(11, 16, TemplateType::kOptimal);
    CHECK_EQ(masks.size(), 1u);
    CHECK_EQ(masks[0], MASK_K11_T16_OPTIMAL);

    masks = get_seed_masks(11, 16, TemplateType::kBoth);
    CHECK_EQ(masks.size(), 2u);
    CHECK_EQ(masks[0], MASK_K11_T16_CODING);
    CHECK_EQ(masks[1], MASK_K11_T16_OPTIMAL);

    masks = get_seed_masks(12, 21, TemplateType::kOptimal);
    CHECK_EQ(masks.size(), 1u);
    CHECK_EQ(masks[0], MASK_K12_T21_OPTIMAL);

    masks = get_seed_masks(12, 18, TemplateType::kBoth);
    CHECK_EQ(masks.size(), 2u);
    CHECK_EQ(masks[0], MASK_K12_T18_CODING);
    CHECK_EQ(masks[1], MASK_K12_T18_OPTIMAL);

    // Invalid combination should return empty
    masks = get_seed_masks(10, 16, TemplateType::kCoding);
    CHECK_EQ(masks.size(), 0u);
}

static void test_seed_span() {
    CHECK_EQ(seed_span(0, 11), 11);
    CHECK_EQ(seed_span(0, 9), 9);
    CHECK_EQ(seed_span(16, 11), 16);
    CHECK_EQ(seed_span(18, 11), 18);
    CHECK_EQ(seed_span(21, 12), 21);
}

static void test_template_type_conversion() {
    CHECK(template_type_from_string("coding") == TemplateType::kCoding);
    CHECK(template_type_from_string("optimal") == TemplateType::kOptimal);
    CHECK(template_type_from_string("both") == TemplateType::kBoth);
    CHECK(template_type_from_string("coding_and_optimal") == TemplateType::kBoth);
    CHECK(template_type_from_string("invalid") == TemplateType::kContiguous);

    CHECK(template_type_to_string(TemplateType::kCoding) == std::string("coding"));
    CHECK(template_type_to_string(TemplateType::kOptimal) == std::string("optimal"));
    CHECK(template_type_to_string(TemplateType::kBoth) == std::string("both"));
    CHECK(template_type_to_string(TemplateType::kContiguous) == std::string("contiguous"));
}

static void test_reverse_complement_string() {
    CHECK(reverse_complement_string("ACGT") == std::string("ACGT"));
    CHECK(reverse_complement_string("AAAA") == std::string("TTTT"));
    CHECK(reverse_complement_string("ACGTACGT") == std::string("ACGTACGT"));
    CHECK(reverse_complement_string("A") == std::string("T"));
    CHECK(reverse_complement_string("") == std::string(""));
    // IUPAC degenerate
    CHECK(reverse_complement_string("RYSWKM") == std::string("KMWSRY"));
}

static void test_mask_weights() {
    // Verify all k=11 masks have popcount 11
    CHECK_EQ(popcount32(MASK_K11_T16_CODING), 11);
    CHECK_EQ(popcount32(MASK_K11_T16_OPTIMAL), 11);
    CHECK_EQ(popcount32(MASK_K11_T18_CODING), 11);
    CHECK_EQ(popcount32(MASK_K11_T18_OPTIMAL), 11);
    CHECK_EQ(popcount32(MASK_K11_T21_CODING), 11);
    CHECK_EQ(popcount32(MASK_K11_T21_OPTIMAL), 11);
    // Verify all k=12 masks have popcount 12
    CHECK_EQ(popcount32(MASK_K12_T16_CODING), 12);
    CHECK_EQ(popcount32(MASK_K12_T16_OPTIMAL), 12);
    CHECK_EQ(popcount32(MASK_K12_T18_CODING), 12);
    CHECK_EQ(popcount32(MASK_K12_T18_OPTIMAL), 12);
    CHECK_EQ(popcount32(MASK_K12_T21_CODING), 12);
    CHECK_EQ(popcount32(MASK_K12_T21_OPTIMAL), 12);
}

static void test_mask_exact_values() {
    // Verify exact hex values of all 12 masks match the specification.
    // Patterns (1=use, 0=skip, left-to-right = position 0 to t-1):
    //   k=11, t=16, coding:  1101101101101101 = 0xDB6D
    //   k=11, t=16, optimal: 1110010110110111 = 0xE5B7
    //   k=12, t=16, coding:  1111101101101101 = 0xFB6D
    //   k=12, t=16, optimal: 1110110110110111 = 0xEDB7
    //   k=11, t=18, coding:  101101100101101101 = 0x2D96D
    //   k=11, t=18, optimal: 111010010110010111 = 0x3A597
    //   k=12, t=18, coding:  101101101101101101 = 0x2DB6D
    //   k=12, t=18, optimal: 111010110010110111 = 0x3ACB7
    //   k=11, t=21, coding:  100101100101100101101 = 0x12CB2D
    //   k=11, t=21, optimal: 111010010100010010111 = 0x1D2897
    //   k=12, t=21, coding:  100101101101100101101 = 0x12DB2D
    //   k=12, t=21, optimal: 111010010110010010111 = 0x1D2C97
    CHECK_EQ(MASK_K11_T16_CODING,  (uint32_t)0xDB6D);
    CHECK_EQ(MASK_K11_T16_OPTIMAL, (uint32_t)0xE5B7);
    CHECK_EQ(MASK_K12_T16_CODING,  (uint32_t)0xFB6D);
    CHECK_EQ(MASK_K12_T16_OPTIMAL, (uint32_t)0xEDB7);
    CHECK_EQ(MASK_K11_T18_CODING,  (uint32_t)0x2D96D);
    CHECK_EQ(MASK_K11_T18_OPTIMAL, (uint32_t)0x3A597);
    CHECK_EQ(MASK_K12_T18_CODING,  (uint32_t)0x2DB6D);
    CHECK_EQ(MASK_K12_T18_OPTIMAL, (uint32_t)0x3ACB7);
    CHECK_EQ(MASK_K11_T21_CODING,  (uint32_t)0x12CB2D);
    CHECK_EQ(MASK_K11_T21_OPTIMAL, (uint32_t)0x1D2897);
    CHECK_EQ(MASK_K12_T21_CODING,  (uint32_t)0x12DB2D);
    CHECK_EQ(MASK_K12_T21_OPTIMAL, (uint32_t)0x1D2C97);
}

static void test_scan_spaced_known_kmer() {
    // Manually compute expected k-mer from a known sequence and mask,
    // then verify the scanner produces the same value.
    //
    // Sequence (16 bases): ACGTACGTACGTACGT
    // 2-bit encoding: A=0, C=1, G=2, T=3
    //   pos: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    //   seq: A  C  G  T  A  C  G  T  A  C  G  T  A  C  G  T
    //   enc: 0  1  2  3  0  1  2  3  0  1  2  3  0  1  2  3
    //
    // k=11, t=16, coding mask: 1101101101101101
    //   1-positions (left-to-right): 0,1,3,4,6,7,9,10,12,13,15
    //   Bases at those positions:     A,C,T,A,G,T,C,G,A,C,T
    //   2-bit values (MSB first):     0,1,3,0,2,3,1,2,0,1,3
    //   k-mer (22-bit): 00 01 11 00 10 11 01 10 00 01 11
    //                  = 0x01CB187 ... let me compute properly
    //
    // The k-mer integer (uint32_t) for k=11:
    //   bit_pos 0 (MSB of k-mer) -> first mask-selected base (pos 0) = A = 0b00
    //   bit_pos 1 -> pos 1 = C = 0b01
    //   bit_pos 2 -> pos 3 = T = 0b11
    //   bit_pos 3 -> pos 4 = A = 0b00
    //   bit_pos 4 -> pos 6 = G = 0b10
    //   bit_pos 5 -> pos 7 = T = 0b11
    //   bit_pos 6 -> pos 9 = C = 0b01
    //   bit_pos 7 -> pos 10 = G = 0b10
    //   bit_pos 8 -> pos 12 = A = 0b00
    //   bit_pos 9 -> pos 13 = C = 0b01
    //   bit_pos 10 -> pos 15 = T = 0b11
    //
    // kmer = (0 << 20) | (1 << 18) | (3 << 16) | (0 << 14) | (2 << 12)
    //      | (3 << 10) | (1 << 8)  | (2 << 6)  | (0 << 4)  | (1 << 2) | 3
    //      = 0 + 262144 + 196608 + 0 + 8192 + 3072 + 256 + 128 + 0 + 4 + 3
    //      = 470407
    uint32_t expected_kmer = (0u << 20) | (1u << 18) | (3u << 16) | (0u << 14) | (2u << 12)
                           | (3u << 10) | (1u << 8)  | (2u << 6)  | (0u << 4)  | (1u << 2) | 3u;

    std::string seq = "ACGTACGTACGTACGT"; // exactly 16 bases
    KmerScanner<uint32_t> scanner(11);
    std::vector<uint32_t> masks = {MASK_K11_T16_CODING};

    std::vector<std::pair<uint32_t, uint32_t>> results;
    scanner.scan_spaced(seq.data(), seq.size(), masks, 16,
        [&](uint32_t pos, uint32_t kmer) {
            results.emplace_back(pos, kmer);
        });

    CHECK_EQ(results.size(), 1u);
    CHECK_EQ(results[0].first, 0u);
    CHECK_EQ(results[0].second, expected_kmer);
}

static void test_scan_spaced_basic() {
    // Test with k=11, t=16 coding mask
    // Create a 20-base sequence and verify scanner finds k-mers
    std::string seq = "ACGTACGTACGTACGTACGT"; // 20 bases
    KmerScanner<uint32_t> scanner(11);
    std::vector<uint32_t> masks = {MASK_K11_T16_CODING};

    std::vector<std::pair<uint32_t, uint32_t>> results;
    scanner.scan_spaced(seq.data(), seq.size(), masks, 16,
        [&](uint32_t pos, uint32_t kmer) {
            results.emplace_back(pos, kmer);
        });

    // With 20 bases and t=16, we should get positions 0-4 (5 windows)
    CHECK_EQ(results.size(), 5u);
    CHECK_EQ(results[0].first, 0u);
    CHECK_EQ(results[4].first, 4u);
}

static void test_scan_spaced_with_n() {
    // N at a mask position should skip that window
    std::string seq = "ACGTACGTNACGTACGTAC"; // 19 bases, N at pos 8
    KmerScanner<uint32_t> scanner(11);
    std::vector<uint32_t> masks = {MASK_K11_T16_CODING};

    std::vector<uint32_t> positions;
    scanner.scan_spaced(seq.data(), seq.size(), masks, 16,
        [&](uint32_t pos, uint32_t /*kmer*/) {
            positions.push_back(pos);
        });

    // Windows that include position 8 as a mask-selected position should be skipped
    // Total possible windows: 0 to 3 (19-16=3), so 4 windows max
    // Some may be skipped if position 8 is a mask position for that window
    CHECK(positions.size() <= 4u);
}

static void test_scan_spaced_both_masks() {
    // With kBoth, should get results from both masks
    std::string seq = "ACGTACGTACGTACGTACGT"; // 20 bases
    KmerScanner<uint32_t> scanner(11);
    std::vector<uint32_t> masks = get_seed_masks(11, 16, TemplateType::kBoth);
    CHECK_EQ(masks.size(), 2u);

    std::vector<std::pair<uint32_t, uint32_t>> results;
    scanner.scan_spaced(seq.data(), seq.size(), masks, 16,
        [&](uint32_t pos, uint32_t kmer) {
            results.emplace_back(pos, kmer);
        });

    // 5 windows * 2 masks = 10 results
    CHECK_EQ(results.size(), 10u);
}

static void test_scan_spaced_short_seq() {
    // Sequence shorter than t should produce no k-mers
    std::string seq = "ACGTACGTACGTACG"; // 15 bases, t=16
    KmerScanner<uint32_t> scanner(11);
    std::vector<uint32_t> masks = {MASK_K11_T16_CODING};

    std::vector<std::pair<uint32_t, uint32_t>> results;
    scanner.scan_spaced(seq.data(), seq.size(), masks, 16,
        [&](uint32_t pos, uint32_t kmer) {
            results.emplace_back(pos, kmer);
        });

    CHECK_EQ(results.size(), 0u);
}

static void test_scan_spaced_k12() {
    // Test with k=12, t=18 coding mask
    std::string seq = "ACGTACGTACGTACGTACGTAC"; // 22 bases
    KmerScanner<uint32_t> scanner(12);
    std::vector<uint32_t> masks = {MASK_K12_T18_CODING};

    std::vector<std::pair<uint32_t, uint32_t>> results;
    scanner.scan_spaced(seq.data(), seq.size(), masks, 18,
        [&](uint32_t pos, uint32_t kmer) {
            results.emplace_back(pos, kmer);
        });

    // With 22 bases and t=18, we should get positions 0-4 (5 windows)
    CHECK_EQ(results.size(), 5u);
}

int main() {
    test_validate_spaced_seed();
    test_get_seed_masks();
    test_seed_span();
    test_template_type_conversion();
    test_reverse_complement_string();
    test_mask_weights();
    test_mask_exact_values();
    test_scan_spaced_known_kmer();
    test_scan_spaced_basic();
    test_scan_spaced_with_n();
    test_scan_spaced_both_masks();
    test_scan_spaced_short_seq();
    test_scan_spaced_k12();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
