#include "test_util.hpp"
#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

// Round-trip helper: encode `offsets` (size D) with U_raw, open the
// resulting blob, and verify access(i) returns offsets[i] for all i.
void check_roundtrip(const std::vector<std::uint64_t>& offsets,
                     std::uint64_t U_raw,
                     const char* label) {
    std::vector<std::uint8_t> blob;
    std::size_t n = ef::encode_dictionary_ef(offsets.data(),
                                             offsets.size(),
                                             U_raw, blob);
    CHECK_EQ(n, blob.size());
    CHECK(blob.size() >= sizeof(ef::EFHeader));

    ef::EFDictionary dict;
    bool ok = dict.open(blob.data(), blob.size());
    if (!ok) {
        std::fprintf(stderr, "FAIL [%s]: EFDictionary::open returned false\n", label);
        return;
    }
    CHECK_EQ(dict.size(), offsets.size());

    for (std::size_t i = 0; i < offsets.size(); ++i) {
        std::uint64_t got = dict.access(static_cast<std::uint32_t>(i));
        if (got != offsets[i]) {
            std::fprintf(stderr,
                "FAIL [%s]: access(%zu) = %llu, want %llu\n",
                label, i,
                static_cast<unsigned long long>(got),
                static_cast<unsigned long long>(offsets[i]));
            // Bump fail counter once per mismatch.
            CHECK_EQ(got, offsets[i]);
            // Stop early so a wrong codec doesn't spam thousands of lines.
            return;
        }
    }
}

void test_empty() {
    std::vector<std::uint64_t> offsets;
    std::vector<std::uint8_t> blob;
    std::size_t n = ef::encode_dictionary_ef(offsets.data(), 0, 0, blob);
    CHECK_EQ(n, sizeof(ef::EFHeader));

    ef::EFDictionary dict;
    bool ok = dict.open(blob.data(), blob.size());
    CHECK(ok);
    CHECK_EQ(dict.size(), 0u);
}

void test_single_zero() {
    std::vector<std::uint64_t> offsets = {0};
    check_roundtrip(offsets, 1, "single_zero");
}

void test_single_large() {
    std::vector<std::uint64_t> offsets = {1234567ULL};
    check_roundtrip(offsets, 2000000ULL, "single_large");
}

void test_all_zero_deltas() {
    // Non-strictly monotonic input: all entries equal.  The strict-
    // monotonicity transform turns this into a perfectly increasing
    // sequence that EF must round-trip exactly.
    std::vector<std::uint64_t> offsets(64, 0ULL);
    check_roundtrip(offsets, 1, "all_zero_deltas");
}

void test_dense_monotonic() {
    // Cover the case l == 0 (U/D == 1).
    std::vector<std::uint64_t> offsets;
    for (std::uint64_t i = 0; i < 256; ++i) offsets.push_back(i);
    check_roundtrip(offsets, 256, "dense_monotonic");
}

void test_sparse_monotonic() {
    // Wide gaps stress the upper-bit unary stream.
    std::vector<std::uint64_t> offsets;
    for (std::uint64_t i = 0; i < 1000; ++i) offsets.push_back(i * 12345ULL);
    check_roundtrip(offsets, 12345ULL * 1000ULL + 1, "sparse_monotonic");
}

void test_runs_and_jumps() {
    // Mix runs of repeats with sparse jumps; representative of real
    // .kix dictionaries where many adjacent k-mers share an offset
    // (empty posting lists) interleaved with large jumps for hot
    // k-mers.
    std::vector<std::uint64_t> offsets;
    std::uint64_t off = 0;
    for (int i = 0; i < 1024; ++i) {
        offsets.push_back(off);
        if (i % 3 == 0) off += 100;     // small step
        else if (i % 17 == 0) off += 50000; // big jump
        // else: same offset (run)
    }
    check_roundtrip(offsets, off + 1, "runs_and_jumps");
}

void test_select_step_boundaries() {
    // Boundary lengths around kSelectStep (64): D = 1, 63, 64, 65,
    // 127, 128, 129, etc.  Catches off-by-one in select-sample
    // indexing.
    for (std::size_t D : {std::size_t(1), std::size_t(63), std::size_t(64),
                          std::size_t(65), std::size_t(127), std::size_t(128),
                          std::size_t(129), std::size_t(255), std::size_t(256),
                          std::size_t(257)}) {
        std::vector<std::uint64_t> offsets;
        offsets.reserve(D);
        for (std::size_t i = 0; i < D; ++i) {
            offsets.push_back(i * 7ULL);
        }
        char buf[64];
        std::snprintf(buf, sizeof(buf), "select_boundary_D=%zu", D);
        check_roundtrip(offsets, 7ULL * D + 1, buf);
    }
}

void test_random_strict_monotonic() {
    std::mt19937_64 rng(0xC0FFEEULL);
    std::vector<std::uint64_t> offsets;
    std::uint64_t off = 0;
    for (int i = 0; i < 5000; ++i) {
        off += (rng() % 100);
        offsets.push_back(off);
    }
    check_roundtrip(offsets, off + 1, "random_5000");
}

void test_random_with_repeats() {
    std::mt19937_64 rng(0xBEEFULL);
    std::vector<std::uint64_t> offsets;
    std::uint64_t off = 0;
    for (int i = 0; i < 5000; ++i) {
        if ((rng() & 3) != 0) {
            // 75 % chance: keep offset (empty posting list)
        } else {
            off += (rng() % 1000) + 1;
        }
        offsets.push_back(off);
    }
    check_roundtrip(offsets, off + 1, "random_with_repeats_5000");
}

void test_simulated_kix_k9_layout() {
    // Simulate a .kix dictionary for k=9 (4^9 = 262144 entries +
    // sentinel).  Posting list sizes drawn from a heavy-tailed
    // distribution so the offset stream has both runs and large
    // jumps.
    std::mt19937_64 rng(0xABCDULL);
    constexpr std::size_t M = 4096;  // small subset to keep test runtime sane
    std::vector<std::uint64_t> offsets;
    offsets.reserve(M + 1);
    std::uint64_t off = 0;
    for (std::size_t i = 0; i < M; ++i) {
        offsets.push_back(off);
        // Heavy-tailed: most are tiny, a few are huge.
        std::uint32_t r = rng() & 0xFF;
        std::uint64_t step;
        if (r < 200) step = (rng() & 0x3F);                 // < 64 B
        else if (r < 250) step = 100 + (rng() & 0x3FF);     // 100..1124
        else step = 10000 + (rng() & 0xFFFF);               // big
        off += step;
    }
    offsets.push_back(off);  // sentinel
    check_roundtrip(offsets, off + 1, "sim_kix_k9_subset");
}

void test_access_pair() {
    std::vector<std::uint64_t> offsets;
    for (std::uint64_t i = 0; i < 100; ++i) offsets.push_back(i * 13ULL);
    std::vector<std::uint8_t> blob;
    ef::encode_dictionary_ef(offsets.data(), offsets.size(),
                             100 * 13ULL, blob);
    ef::EFDictionary dict;
    CHECK(dict.open(blob.data(), blob.size()));
    for (std::uint32_t i = 0; i + 1 < offsets.size(); ++i) {
        std::uint64_t s, e;
        dict.access_pair(i, s, e);
        CHECK_EQ(s, offsets[i]);
        CHECK_EQ(e, offsets[i + 1]);
    }
}

void test_open_rejects_bad_magic() {
    std::vector<std::uint8_t> blob(sizeof(ef::EFHeader), 0);
    blob[0] = 'X'; blob[1] = 'Y'; blob[2] = 'Z'; blob[3] = '0';
    ef::EFDictionary dict;
    CHECK(!dict.open(blob.data(), blob.size()));
}

void test_open_rejects_truncated() {
    std::vector<std::uint64_t> offsets;
    for (std::uint64_t i = 0; i < 200; ++i) offsets.push_back(i);
    std::vector<std::uint8_t> blob;
    ef::encode_dictionary_ef(offsets.data(), offsets.size(), 200, blob);
    ef::EFDictionary dict;
    CHECK(!dict.open(blob.data(), blob.size() - 1));
}

void test_active_tier_name_present() {
    const char* name = ef::active_tier_name();
    CHECK(name != nullptr);
    CHECK(*name != '\0');
}

} // anonymous namespace

int main() {
    test_empty();
    test_single_zero();
    test_single_large();
    test_all_zero_deltas();
    test_dense_monotonic();
    test_sparse_monotonic();
    test_runs_and_jumps();
    test_select_step_boundaries();
    test_random_strict_monotonic();
    test_random_with_repeats();
    test_simulated_kix_k9_layout();
    test_access_pair();
    test_open_rejects_bad_magic();
    test_open_rejects_truncated();
    test_active_tier_name_present();

    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
