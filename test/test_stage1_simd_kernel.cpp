#include "test_util.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage1_filter_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

template <Stage1Width Width>
void scalar_golden(const SeqId* sid_batch, int count,
                   typename Stage1WidthTraits<Width>::PosT q_pos,
                   typename Stage1WidthTraits<Width>::ScoreT* scores,
                   typename Stage1WidthTraits<Width>::PosT* last_pos,
                   std::vector<uint32_t>& dirty) {
    for (int i = 0; i < count; i++) {
        SeqId sid = sid_batch[i];
        if (scores[sid] == 0) dirty.push_back(sid);
        if (last_pos[sid] != q_pos) {
            scores[sid]++;
            last_pos[sid] = q_pos;
        }
    }
}

template <Stage1Width Width>
void verify_case(const std::vector<SeqId>& sids,
                 typename Stage1WidthTraits<Width>::PosT q_pos,
                 const std::vector<typename Stage1WidthTraits<Width>::ScoreT>& scores_init,
                 const std::vector<typename Stage1WidthTraits<Width>::PosT>& last_pos_init,
                 const char* label) {
    using ScoreT = typename Stage1WidthTraits<Width>::ScoreT;
    using PosT   = typename Stage1WidthTraits<Width>::PosT;

    auto scores_g   = scores_init;
    auto last_pos_g = last_pos_init;
    std::vector<uint32_t> dirty_g;
    scalar_golden<Width>(sids.data(), static_cast<int>(sids.size()), q_pos,
                        scores_g.data(), last_pos_g.data(), dirty_g);

    auto scores_s   = scores_init;
    auto last_pos_s = last_pos_init;
    std::vector<uint32_t> dirty_s;
    flush_batch_simd<Width>(sids.data(), static_cast<int>(sids.size()), q_pos,
                           static_cast<void*>(scores_s.data()),
                           static_cast<void*>(last_pos_s.data()),
                           dirty_s);

    bool ok = (scores_g == scores_s) && (last_pos_g == last_pos_s) && (dirty_g == dirty_s);
    if (!ok) {
        std::fprintf(stderr, "FAIL[%s]: scores match=%d, last_pos match=%d, dirty match=%d\n",
                     label,
                     static_cast<int>(scores_g == scores_s),
                     static_cast<int>(last_pos_g == last_pos_s),
                     static_cast<int>(dirty_g == dirty_s));
        if (dirty_g != dirty_s) {
            std::fprintf(stderr, "  scalar dirty (%zu): ", dirty_g.size());
            for (auto v : dirty_g) std::fprintf(stderr, "%u ", v);
            std::fprintf(stderr, "\n  simd   dirty (%zu): ", dirty_s.size());
            for (auto v : dirty_s) std::fprintf(stderr, "%u ", v);
            std::fprintf(stderr, "\n");
        }
    }
    CHECK(ok);
    (void)label;
    (void)q_pos;  // suppress unused warning under some compilers
    using NoWarn = ScoreT;
    using NoWarn2 = PosT;
    (void)sizeof(NoWarn);
    (void)sizeof(NoWarn2);
}

void test_t32_basic() {
    constexpr uint32_t NUM_SEQS = 1000;
    constexpr uint32_t SENT = std::numeric_limits<uint32_t>::max();

    auto fresh = [&](const std::vector<SeqId>& sids, uint32_t q_pos, const char* label) {
        std::vector<uint32_t> scores(NUM_SEQS, 0);
        std::vector<uint32_t> lp(NUM_SEQS, SENT);
        verify_case<Stage1Width::T32>(sids, q_pos, scores, lp, label);
    };

    fresh({}, 42, "t32:empty");
    fresh({0}, 42, "t32:single0");
    fresh({500}, 42, "t32:single500");

    std::mt19937 rng(42);
    for (int n : {1, 7, 8, 9, 15, 16, 17, 31, 32, 64, 100}) {
        std::vector<SeqId> sids;
        sids.reserve(n);
        for (int i = 0; i < n; i++) sids.push_back(rng() % NUM_SEQS);
        char label[64];
        std::snprintf(label, sizeof(label), "t32:rand_n=%d", n);
        fresh(sids, 100, label);
    }

    // Duplicates within and across chunks.
    fresh({5, 10, 5, 20, 10, 30, 5}, 50, "t32:dup_small");
    fresh(std::vector<SeqId>(16, 7), 100, "t32:all_same_16");
    fresh(std::vector<SeqId>(33, 11), 100, "t32:all_same_33");

    // Aligned 16 unique
    {
        std::vector<SeqId> sids;
        for (int i = 0; i < 16; i++) sids.push_back(i * 10);
        fresh(sids, 100, "t32:16_unique");
    }

    // Cross-chunk duplicates
    {
        std::vector<SeqId> sids;
        for (int i = 0; i < 32; i++) sids.push_back((i * 7 + 3) % 50);
        fresh(sids, 99, "t32:32_dup_cross");
    }

    // q_pos boundary values
    fresh({1, 2, 3}, 0, "t32:qpos_0");
    fresh({1, 2, 3}, SENT, "t32:qpos_sent");

    // Pre-populated state: scores[5] != 0, last_pos[5] != q_pos
    {
        std::vector<uint32_t> scores(NUM_SEQS, 0);
        std::vector<uint32_t> lp(NUM_SEQS, SENT);
        scores[5] = 7; lp[5] = 99;
        std::vector<SeqId> sids = {5, 10, 5};
        verify_case<Stage1Width::T32>(sids, 100, scores, lp, "t32:prepop_scoreNZ");
    }

    // last_pos == q_pos case (no update)
    {
        std::vector<uint32_t> scores(NUM_SEQS, 0);
        std::vector<uint32_t> lp(NUM_SEQS, SENT);
        scores[5] = 1; lp[5] = 100;
        std::vector<SeqId> sids = {5};
        verify_case<Stage1Width::T32>(sids, 100, scores, lp, "t32:lp_eq_qpos");
    }

    // Mixed: lane 0 unique with score==0, lane 1 conflict with lane 0
    {
        std::vector<uint32_t> scores(NUM_SEQS, 0);
        std::vector<uint32_t> lp(NUM_SEQS, SENT);
        std::vector<SeqId> sids;
        for (int i = 0; i < 16; i++) sids.push_back(i % 4);  // 4 unique sids, each x4
        verify_case<Stage1Width::T32>(sids, 100, scores, lp, "t32:4unique_x4");
    }

    // Larger random with overlaps
    {
        std::mt19937 rng2(7);
        std::vector<SeqId> sids;
        sids.reserve(200);
        for (int i = 0; i < 200; i++) sids.push_back(rng2() % 64);  // many duplicates
        std::vector<uint32_t> scores(NUM_SEQS, 0);
        std::vector<uint32_t> lp(NUM_SEQS, SENT);
        verify_case<Stage1Width::T32>(sids, 1234, scores, lp, "t32:200_dense_dup");
    }
}

void test_t8_t16() {
    // T8 / T16 always go through scalar; verify the dispatcher still produces
    // correct output (same as the scalar reference).
    constexpr uint32_t NUM_SEQS = 200;

    {
        std::vector<uint8_t> scores(NUM_SEQS, 0);
        std::vector<uint8_t> lp(NUM_SEQS, std::numeric_limits<uint8_t>::max());
        std::vector<SeqId> sids = {3, 7, 3, 100, 50, 7, 7};
        verify_case<Stage1Width::T8>(sids, 50, scores, lp, "t8:basic");
    }
    {
        std::vector<uint16_t> scores(NUM_SEQS, 0);
        std::vector<uint16_t> lp(NUM_SEQS, std::numeric_limits<uint16_t>::max());
        std::vector<SeqId> sids = {3, 7, 3, 100, 50, 7, 7};
        verify_case<Stage1Width::T16>(sids, 50, scores, lp, "t16:basic");
    }
}

} // namespace

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    test_t32_basic();
    test_t8_t16();

    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
