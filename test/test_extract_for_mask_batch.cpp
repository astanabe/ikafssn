// Phase 4b-1 — bit-exact tests for the SIMD batch spaced-seed extractor.
//
// All 24 canonical masks are verified against the scalar reference
// (extract_for_mask) at the exact n boundaries that drive each tier's chunk
// width (2 / 4 / 8) plus a "many chunks" size. The test_kmer_encoding-style
// test_util.hpp framework is reused; init_simd_dispatch() is called in main
// so the IKAFSSN_FORCE_SIMD env var can pin a specific tier.

#include "test_util.hpp"

#include "core/spaced_seed.hpp"
#include "core/spaced_seed_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <random>
#include <type_traits>
#include <vector>

using namespace ikafssn;

template <typename KmerInt>
static void check_mask(uint32_t mask, std::mt19937& rng) {
    // Cover lane edges (1, 2-1, 2, 2+1, 4-1, 4, 4+1, 8-1, 8, 8+1) and a
    // larger size that runs the chunk loop many times.
    for (std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{2},
                          std::size_t{3}, std::size_t{4}, std::size_t{5},
                          std::size_t{7}, std::size_t{8}, std::size_t{9},
                          std::size_t{15}, std::size_t{16}, std::size_t{17},
                          std::size_t{63}, std::size_t{64}, std::size_t{65},
                          std::size_t{1024}}) {
        std::vector<std::uint64_t> a(n);
        std::vector<KmerInt> got(n), gold(n);
        for (auto& v : a) {
            // 64-bit random, but only the low 2*t bits matter so any pattern
            // exercises every run.
            v = (std::uint64_t)rng() ^ ((std::uint64_t)rng() << 32);
        }
        for (std::size_t i = 0; i < n; ++i) {
            gold[i] = extract_for_mask<KmerInt>(a[i], mask);
        }
        extract_for_mask_batch<KmerInt>(a.data(), n, mask, got.data());
        for (std::size_t i = 0; i < n; ++i) {
            CHECK_EQ(got[i], gold[i]);
        }
    }
}

static void test_all_masks_bit_exact() {
    std::mt19937 rng(0xBEEFu);
    // K=8 (uint16_t)
    check_mask<std::uint16_t>(MASK_K8_T13_CODING, rng);
    check_mask<std::uint16_t>(MASK_K8_T13_OPTIMAL, rng);
    check_mask<std::uint16_t>(MASK_K8_T15_CODING, rng);
    check_mask<std::uint16_t>(MASK_K8_T15_OPTIMAL, rng);
    check_mask<std::uint16_t>(MASK_K8_T18_CODING, rng);
    check_mask<std::uint16_t>(MASK_K8_T18_OPTIMAL, rng);
    // K=9, 11, 12 (uint32_t)
    check_mask<std::uint32_t>(MASK_K9_T13_CODING, rng);
    check_mask<std::uint32_t>(MASK_K9_T13_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K9_T15_CODING, rng);
    check_mask<std::uint32_t>(MASK_K9_T15_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K9_T18_CODING, rng);
    check_mask<std::uint32_t>(MASK_K9_T18_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K11_T16_CODING, rng);
    check_mask<std::uint32_t>(MASK_K11_T16_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K11_T18_CODING, rng);
    check_mask<std::uint32_t>(MASK_K11_T18_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K11_T21_CODING, rng);
    check_mask<std::uint32_t>(MASK_K11_T21_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K12_T16_CODING, rng);
    check_mask<std::uint32_t>(MASK_K12_T16_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K12_T18_CODING, rng);
    check_mask<std::uint32_t>(MASK_K12_T18_OPTIMAL, rng);
    check_mask<std::uint32_t>(MASK_K12_T21_CODING, rng);
    check_mask<std::uint32_t>(MASK_K12_T21_OPTIMAL, rng);
}

static void test_n_zero_guard() {
    std::uint16_t out16[1] = {0xDEAD};
    std::uint32_t out32[1] = {0xC0FFEE};
    extract_for_mask_batch<std::uint16_t>(nullptr, 0, MASK_K8_T13_CODING,
                                           out16);
    extract_for_mask_batch<std::uint32_t>(nullptr, 0, MASK_K11_T16_CODING,
                                           out32);
    CHECK_EQ(out16[0], 0xDEAD);
    CHECK_EQ(out32[0], 0xC0FFEEu);
}

template <typename KmerInt>
static void check_unaligned(uint32_t mask, std::mt19937& rng) {
    constexpr std::size_t n = 17;  // enough to cross every tier's chunk size
    // Allocate larger backing buffers and shift the pointer by 1 element.
    std::vector<std::uint64_t> a_buf(n + 1);
    std::vector<KmerInt> got_buf(n + 1, KmerInt{0xAA});
    std::vector<KmerInt> gold(n);
    for (auto& v : a_buf) v = (std::uint64_t)rng() ^ ((std::uint64_t)rng() << 32);
    const std::uint64_t* a = a_buf.data() + 1;
    KmerInt* got = got_buf.data() + 1;
    for (std::size_t i = 0; i < n; ++i) {
        gold[i] = extract_for_mask<KmerInt>(a[i], mask);
    }
    extract_for_mask_batch<KmerInt>(a, n, mask, got);
    for (std::size_t i = 0; i < n; ++i) CHECK_EQ(got[i], gold[i]);
}

static void test_unaligned_inputs() {
    std::mt19937 rng(0xACCE55u);
    check_unaligned<std::uint16_t>(MASK_K8_T18_CODING, rng);
    check_unaligned<std::uint32_t>(MASK_K11_T16_CODING, rng);
    check_unaligned<std::uint32_t>(MASK_K12_T21_OPTIMAL, rng);
}

static void test_kBoth_two_masks_separate() {
    // "Both" template_type: caller invokes the batch entry point twice with
    // the coding and optimal masks separately. The two streams must be
    // independently bit-exact.
    std::mt19937 rng(0xDA7Au);
    constexpr std::size_t n = 64;
    std::vector<std::uint64_t> a(n);
    std::vector<std::uint32_t> got_c(n), got_o(n), gold_c(n), gold_o(n);
    for (auto& v : a) v = (std::uint64_t)rng() ^ ((std::uint64_t)rng() << 32);
    for (std::size_t i = 0; i < n; ++i) {
        gold_c[i] = extract_for_mask<std::uint32_t>(a[i], MASK_K11_T16_CODING);
        gold_o[i] = extract_for_mask<std::uint32_t>(a[i], MASK_K11_T16_OPTIMAL);
    }
    extract_for_mask_batch<std::uint32_t>(a.data(), n, MASK_K11_T16_CODING,
                                           got_c.data());
    extract_for_mask_batch<std::uint32_t>(a.data(), n, MASK_K11_T16_OPTIMAL,
                                           got_o.data());
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_EQ(got_c[i], gold_c[i]);
        CHECK_EQ(got_o[i], gold_o[i]);
    }
}

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();
    test_all_masks_bit_exact();
    test_n_zero_guard();
    test_unaligned_inputs();
    test_kBoth_two_masks_separate();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
