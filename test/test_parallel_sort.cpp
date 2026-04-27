#include "test_util.hpp"
#include "index/parallel_sort_dispatch.hpp"
#include "util/simd_dispatch.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

std::vector<TempEntry> golden_sort(std::vector<TempEntry> v) {
    std::sort(v.begin(), v.end(), [](const TempEntry& a, const TempEntry& b) {
        if (a.kmer_value != b.kmer_value) return a.kmer_value < b.kmer_value;
        if (a.seq_id     != b.seq_id)     return a.seq_id < b.seq_id;
        return a.pos < b.pos;
    });
    return v;
}

bool entries_equal(const std::vector<TempEntry>& a, const std::vector<TempEntry>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].kmer_value != b[i].kmer_value) return false;
        if (a[i].seq_id     != b[i].seq_id)     return false;
        if (a[i].pos        != b[i].pos)        return false;
    }
    return true;
}

void run(std::vector<TempEntry> input, const char* label) {
    auto golden = golden_sort(input);
    auto dut    = input;
    parallel_sort_temp_entries(dut);
    bool ok = entries_equal(golden, dut);
    if (!ok) {
        std::fprintf(stderr, "FAIL[%s]: size=%zu\n", label, input.size());
        std::size_t mismatches = 0;
        for (std::size_t i = 0; i < std::min(golden.size(), dut.size()); ++i) {
            if (golden[i].kmer_value != dut[i].kmer_value
             || golden[i].seq_id     != dut[i].seq_id
             || golden[i].pos        != dut[i].pos) {
                if (mismatches < 5) {
                    std::fprintf(stderr,
                        "  [%zu] golden=(%u,%u,%u) dut=(%u,%u,%u)\n", i,
                        golden[i].kmer_value, golden[i].seq_id, golden[i].pos,
                        dut[i].kmer_value, dut[i].seq_id, dut[i].pos);
                }
                ++mismatches;
            }
        }
        std::fprintf(stderr, "  total mismatches: %zu\n", mismatches);
    }
    CHECK(ok);
}

void test_basic() {
    run({}, "empty");
    run({{1, 2, 3}}, "single");
    run({{2, 0, 0}, {1, 0, 0}}, "size2_kmer_swap");
    run({{1, 2, 0}, {1, 1, 0}}, "size2_seqid_swap");
    run({{1, 0, 2}, {1, 0, 1}}, "size2_pos_swap");
}

void test_random_sizes() {
    std::mt19937 rng(7);
    for (std::size_t n : {2u, 3u, 7u, 16u, 17u, 64u, 100u, 1000u, 10000u}) {
        std::vector<TempEntry> v;
        v.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            v.push_back({rng() & 0xFFFFFFu, rng() & 0xFFFFu, rng() & 0xFFFFu});
        }
        char label[64];
        std::snprintf(label, sizeof(label), "rand_n=%zu", n);
        run(std::move(v), label);
    }
}

void test_duplicates_within_kmer() {
    // Same (kmer, seq_id), multiple pos — exercise the local pos-sort step.
    std::vector<TempEntry> v;
    for (int p = 19; p >= 0; --p) v.push_back({100, 5, static_cast<uint32_t>(p)});
    for (int p = 0; p < 20; ++p)  v.push_back({100, 7, static_cast<uint32_t>(p ^ 0x5)});
    for (int p = 0; p < 5; ++p)   v.push_back({99,  3, static_cast<uint32_t>(p)});
    run(std::move(v), "dup_kmer_seqid");
}

void test_already_sorted() {
    std::vector<TempEntry> v;
    for (uint32_t k = 0; k < 100; ++k) {
        for (uint32_t s = 0; s < 5; ++s) {
            for (uint32_t p = 0; p < 5; ++p) {
                v.push_back({k, s, p});
            }
        }
    }
    run(std::move(v), "sorted");
}

void test_reverse_sorted() {
    std::vector<TempEntry> v;
    for (uint32_t k = 100; k > 0; --k) {
        for (uint32_t s = 5; s > 0; --s) {
            for (uint32_t p = 5; p > 0; --p) {
                v.push_back({k - 1, s - 1, p - 1});
            }
        }
    }
    run(std::move(v), "reverse");
}

void test_large_random() {
    std::mt19937 rng(42);
    std::vector<TempEntry> v;
    v.reserve(1 << 17);  // 128K entries
    for (std::size_t i = 0; i < (1u << 17); ++i) {
        v.push_back({rng() & 0xFFFFFFu, rng() & 0xFFFFu, rng() & 0xFFFFu});
    }
    run(std::move(v), "large_rand_128K");
}

} // namespace

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    std::fprintf(stderr, "active sort path: %s\n",
                 parallel_sort_simd_active() ? "x86-simd-sort" : "tbb::parallel_sort");

    test_basic();
    test_random_sizes();
    test_duplicates_within_kmer();
    test_already_sorted();
    test_reverse_sorted();
    test_large_random();

    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
