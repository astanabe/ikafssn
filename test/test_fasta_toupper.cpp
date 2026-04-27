#include "test_util.hpp"
#include "io/text_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

using namespace ikafssn;

// Scalar reference: the contract toupper_inplace_ascii() must match bit-exactly.
static void scalar_toupper(std::uint8_t* p, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        std::uint8_t c = p[i];
        if (c >= 'a' && c <= 'z') p[i] = static_cast<std::uint8_t>(c - 0x20);
    }
}

// Run toupper on `input`, compare against scalar golden, fail if mismatch.
// `head_skew` shifts the start pointer to exercise misaligned inputs.
static void check_pattern(const std::vector<std::uint8_t>& input,
                          std::size_t head_skew,
                          const char* label) {
    if (head_skew >= input.size() && !input.empty()) head_skew = 0;
    std::size_t n = input.size() - head_skew;

    std::vector<std::uint8_t> got(input.begin() + head_skew, input.end());
    std::vector<std::uint8_t> expect = got;
    scalar_toupper(expect.data(), expect.size());
    toupper_inplace_ascii(got.data(), n);

    if (got != expect) {
        std::fprintf(stderr,
            "FAIL: %s (n=%zu skew=%zu) — first mismatch at",
            label, n, head_skew);
        for (std::size_t i = 0; i < n; ++i) {
            if (got[i] != expect[i]) {
                std::fprintf(stderr,
                    " i=%zu got=0x%02X expect=0x%02X\n",
                    i, got[i], expect[i]);
                break;
            }
        }
        g_fail_count++;
    } else {
        g_pass_count++;
    }
}

static void test_sizes_lowercase_only() {
    static constexpr std::size_t sizes[] = {
        0, 1, 8, 15, 16, 31, 32, 63, 64, 127, 128, 255, 256,
        511, 512, 1023, 1024, 65536, 1u << 20
    };
    for (std::size_t n : sizes) {
        std::vector<std::uint8_t> buf(n, 'a');
        for (std::size_t i = 0; i < n; ++i)
            buf[i] = static_cast<std::uint8_t>('a' + (i % 26));
        check_pattern(buf, 0, "lowercase only");
    }
}

static void test_sizes_uppercase_only() {
    static constexpr std::size_t sizes[] = {0, 1, 16, 64, 255, 256, 1024, 65536};
    for (std::size_t n : sizes) {
        std::vector<std::uint8_t> buf(n);
        for (std::size_t i = 0; i < n; ++i)
            buf[i] = static_cast<std::uint8_t>('A' + (i % 26));
        check_pattern(buf, 0, "uppercase only");
    }
}

static void test_mixed_with_iupac_and_specials() {
    // Mixed ASCII printable + IUPAC degenerate + digits + symbols + NUL +
    // 0x80..0xFF (must remain unchanged).
    std::mt19937 rng(0x1234);
    static constexpr std::size_t sizes[] = {
        16, 32, 63, 64, 128, 256, 1024, 65535, 1u << 20
    };
    const char* extras = "0123456789-_.,;:[](){}!@#$%^&*+=<>? \t\r\n\0RYWSKMBDHVNryswskmbdhvn";
    std::size_t n_extras = 64; // includes NUL

    for (std::size_t n : sizes) {
        std::vector<std::uint8_t> buf(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::uint32_t r = rng();
            int kind = r & 0x7;
            std::uint8_t c;
            if (kind < 4) {
                // ASCII letter, mixed case
                c = static_cast<std::uint8_t>(((r >> 8) & 1)
                        ? 'A' + ((r >> 16) % 26)
                        : 'a' + ((r >> 16) % 26));
            } else if (kind == 4) {
                c = static_cast<std::uint8_t>(extras[(r >> 16) % n_extras]);
            } else if (kind == 5) {
                // High byte: must NOT be modified.
                c = static_cast<std::uint8_t>(0x80 | ((r >> 16) & 0x7F));
            } else {
                // Around the [a-z] boundary: '`' = 0x60, '{' = 0x7B
                c = static_cast<std::uint8_t>(0x60 + ((r >> 16) % 30));
            }
            buf[i] = c;
        }
        check_pattern(buf, 0, "mixed");
    }
}

static void test_misaligned_starts() {
    std::mt19937 rng(0xC0FFEE);
    for (std::size_t skew : {1u, 3u, 5u, 7u, 9u, 11u, 13u, 17u}) {
        for (std::size_t total : {64u, 128u, 256u, 1024u, 65536u}) {
            std::vector<std::uint8_t> buf(total + skew);
            for (auto& b : buf)
                b = static_cast<std::uint8_t>(rng() & 0xFF);
            check_pattern(buf, skew, "misaligned");
        }
    }
}

static void test_high_byte_invariant() {
    // All bytes in 0x80..0xFF: nothing should change.
    std::vector<std::uint8_t> buf(4096);
    for (std::size_t i = 0; i < buf.size(); ++i)
        buf[i] = static_cast<std::uint8_t>(0x80 | (i & 0x7F));
    auto golden = buf;
    toupper_inplace_ascii(buf.data(), buf.size());
    CHECK(buf == golden);
}

static void test_boundary_chars() {
    // Exhaustive check on 256 bytes covering 0x00..0xFF — only 'a'..'z' must
    // change, all other bytes must be preserved.
    std::vector<std::uint8_t> buf(256);
    for (std::size_t i = 0; i < 256; ++i)
        buf[i] = static_cast<std::uint8_t>(i);
    std::vector<std::uint8_t> golden(256);
    for (std::size_t i = 0; i < 256; ++i)
        golden[i] = static_cast<std::uint8_t>(
            (i >= 'a' && i <= 'z') ? (i - 0x20) : i);
    toupper_inplace_ascii(buf.data(), buf.size());
    CHECK(buf == golden);
}

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    test_sizes_lowercase_only();
    test_sizes_uppercase_only();
    test_mixed_with_iupac_and_specials();
    test_misaligned_starts();
    test_high_byte_invariant();
    test_boundary_chars();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
