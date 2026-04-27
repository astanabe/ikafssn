#include "test_util.hpp"
#include "core/ncbi2na_unpack.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using namespace ikafssn;

// Reference implementation (single source of truth for expected output).
static void scalar_unpack(const char* packed, std::uint32_t start,
                          std::uint32_t count, std::uint8_t* out) {
    for (std::uint32_t i = 0; i < count; ++i) {
        std::uint32_t pos = start + i;
        std::uint8_t  b   = static_cast<std::uint8_t>(packed[pos >> 2]);
        out[i] = static_cast<std::uint8_t>(
            (b >> (6 - 2 * (pos & 3))) & 0x03);
    }
}

static void check_one(const std::vector<char>& packed,
                      std::uint32_t start, std::uint32_t count,
                      const char* label) {
    std::vector<std::uint8_t> got(count, 0xFFu);
    std::vector<std::uint8_t> exp(count, 0xFFu);
    unpack_ncbi2na_chunk(packed.data(), start, count, got.data());
    scalar_unpack(packed.data(), start, count, exp.data());
    if (got != exp) {
        std::fprintf(stderr,
            "FAIL: %s start=%u count=%u — first mismatch at",
            label, start, count);
        for (std::uint32_t i = 0; i < count; ++i) {
            if (got[i] != exp[i]) {
                std::fprintf(stderr,
                    " i=%u got=%u expect=%u (packed[%u]=0x%02X)\n",
                    i, got[i], exp[i], (start + i) >> 2,
                    static_cast<std::uint8_t>(packed[(start + i) >> 2]));
                break;
            }
        }
        g_fail_count++;
    } else {
        g_pass_count++;
    }
}

static void test_special_byte_patterns() {
    // Build a packed buffer with controlled bytes; exercise count and start
    // covering all 4 residues.
    static constexpr std::uint8_t kPatterns[] = {0x00, 0xFF, 0xAA, 0x55,
                                                 0x1B, 0xE4, 0x39, 0xC6};
    std::vector<char> buf(4096);
    for (auto p : kPatterns) {
        for (auto& b : buf) b = static_cast<char>(p);
        for (std::uint32_t start_off : {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u}) {
            for (std::uint32_t count :
                 {0u, 1u, 3u, 4u, 7u, 8u, 15u, 16u, 31u, 32u,
                  63u, 64u, 65u, 127u, 128u, 1024u, 4000u}) {
                if (start_off + count > buf.size() * 4) continue;
                check_one(buf, start_off, count, "uniform pattern");
            }
        }
    }
}

static void test_random_inputs() {
    std::mt19937 rng(0xCAFEBABE);
    std::vector<char> buf(8192);
    for (auto& b : buf) b = static_cast<char>(rng() & 0xFF);
    for (std::uint32_t start_off : {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 11u, 17u}) {
        for (std::uint32_t count :
             {0u, 1u, 4u, 16u, 32u, 64u, 100u, 256u, 1023u, 1024u, 1025u,
              4096u, 16384u}) {
            if (start_off + count > buf.size() * 4) continue;
            check_one(buf, start_off, count, "random");
        }
    }
}

static void test_periodic_byte_table() {
    // 256-byte periodic input covers every possible packed byte.
    std::vector<char> buf(2048);
    for (std::size_t i = 0; i < buf.size(); ++i)
        buf[i] = static_cast<char>(i & 0xFF);
    for (std::uint32_t start_off : {0u, 1u, 2u, 3u}) {
        for (std::uint32_t count : {64u, 128u, 256u, 1024u, 4096u, 8000u}) {
            if (start_off + count > buf.size() * 4) continue;
            check_one(buf, start_off, count, "periodic 256");
        }
    }
}

static void test_zero_count() {
    std::vector<char> buf(64, 0x55);
    std::uint8_t sentinel[1] = {0xAA};
    unpack_ncbi2na_chunk(buf.data(), 0, 0, sentinel);
    CHECK_EQ(sentinel[0], static_cast<std::uint8_t>(0xAA));
}

static void test_chunk_boundary() {
    // Cross the 64-base SIMD chunk boundary at every offset to make sure the
    // tail/head split logic is correct.
    std::vector<char> buf(4096);
    std::mt19937 rng(0xDEADBEEF);
    for (auto& b : buf) b = static_cast<char>(rng() & 0xFF);
    for (std::uint32_t start = 0; start < 64; ++start) {
        for (std::uint32_t count : {60u, 64u, 65u, 66u, 67u, 68u,
                                    127u, 128u, 129u, 256u}) {
            check_one(buf, start, count, "chunk boundary");
        }
    }
}

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    test_zero_count();
    test_special_byte_patterns();
    test_random_inputs();
    test_periodic_byte_table();
    test_chunk_boundary();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
