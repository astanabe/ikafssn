#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "search/stage1_filter.hpp"
#include "search/hit_limits.hpp"
#include "search/parallel_search.hpp"
#include "volume_search_helper.hpp"
#include "search/query_preprocessor.hpp"
#include "search/oid_filter.hpp"
#include "index/index_builder.hpp"
#include "index/index_filter.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "io/blastdb_reader.hpp"
#include "core/kmer_encoding.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"
#include "util/simd_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <string>
#include <unordered_set>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_index_dir;

// Runtime-extracted query and OID for FJ876973.1
static std::string g_query_seq;
static uint32_t g_fj_oid = UINT32_MAX;

static bool build_test_index() {
    BlastDbReader db;
    if (!db.open(g_testdb_path)) return false;

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 7;

    std::string prefix = g_index_dir + "/test.00.07mer";
    return build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger);
}

static void test_stage1_basic() {
    std::fprintf(stderr, "-- test_stage1_basic\n");

    KixReader kix;
    std::string kix_path = g_index_dir + "/test.00.07mer.kix";
    CHECK(kix.open(kix_path));

    // Query: 100bp from FJ876973.1 (extracted at runtime)
    std::vector<uint32_t> positions;
    std::vector<uint16_t> kmer_values;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        positions.push_back(pos);
        kmer_values.push_back(kmer);
    });

    OidFilter filter; // no filter
    Stage1Config config;
    config.min_stage1_score = 1;

    Stage1Buffer buf;
    auto candidates = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, filter, config, buf);

    // Should find FJ876973.1 OID in candidates
    CHECK(!candidates.empty());
    bool found_fj = false;
    for (const auto& c : candidates) {
        if (c.id == g_fj_oid) found_fj = true;
    }
    CHECK(found_fj);

    kix.close();
}

static void test_stage1_min_score() {
    std::fprintf(stderr, "-- test_stage1_min_score\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<uint32_t> positions;
    std::vector<uint16_t> kmer_values;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        positions.push_back(pos);
        kmer_values.push_back(kmer);
    });

    OidFilter filter;
    Stage1Config config;
    config.min_stage1_score = 999999; // Very high threshold

    Stage1Buffer buf;
    auto candidates = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, filter, config, buf);
    CHECK(candidates.empty());

    kix.close();
}

static void test_stage1_with_oid_filter() {
    std::fprintf(stderr, "-- test_stage1_with_oid_filter\n");

    KixReader kix;
    KsxReader ksx;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));
    CHECK(ksx.open(g_index_dir + "/test.00.07mer.ksx"));

    std::vector<uint32_t> positions;
    std::vector<uint16_t> kmer_values;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        positions.push_back(pos);
        kmer_values.push_back(kmer);
    });

    // Filter: include only ACC_GQ and ACC_DQ (exclude FJ876973.1)
    OidFilter filter;
    filter.build({ACC_GQ, ACC_DQ}, ksx, OidFilterMode::kInclude);

    Stage1Config config;
    config.min_stage1_score = 1;

    Stage1Buffer buf;
    auto candidates = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, filter, config, buf);

    // FJ876973.1 OID should NOT be in results
    for (const auto& c : candidates) {
        CHECK(c.id != g_fj_oid);
    }

    kix.close();
    ksx.close();
}

// Stage 1 candidate-count selectors (synthetic; no index needed).
static void test_stage1_limit_topk() {
    std::fprintf(stderr, "-- test_stage1_limit_topk\n");
    // scores (desc): 5,5,4,3,3,2
    auto make = []() {
        return std::vector<Stage1Candidate>{
            {10, 5}, {11, 5}, {12, 4}, {13, 3}, {14, 3}, {15, 2}};
    };
    // tie-inclusive k=3: 3rd score is 4 -> keep score>=4 => {5,5,4}
    { auto c = make(); stage1_limit_topk(c, 3, true); CHECK(c.size() == 3); }
    // tie-inclusive k=4: 4th score is 3 -> keep score>=3 => {5,5,4,3,3}
    { auto c = make(); stage1_limit_topk(c, 4, true); CHECK(c.size() == 5); }
    // exact k=4: keep exactly 4
    { auto c = make(); stage1_limit_topk(c, 4, false); CHECK(c.size() == 4); }
    // k=0 and k>=size are no-ops
    { auto c = make(); stage1_limit_topk(c, 0, true);  CHECK(c.size() == 6); }
    { auto c = make(); stage1_limit_topk(c, 10, true); CHECK(c.size() == 6); }
}

static void test_stage1_limit_per_parent() {
    std::fprintf(stderr, "-- test_stage1_limit_per_parent\n");
    // sid -> parent: {0,0,0,1,1,1}
    std::vector<uint32_t> parent = {0, 0, 0, 1, 1, 1};
    auto make = []() {
        return std::vector<Stage1Candidate>{
            {0, 5}, {1, 4}, {2, 4}, {3, 9}, {4, 1}, {5, 1}};
    };
    // n=1 tie-inclusive: parent0 top1=5, parent1 top1=9 -> 2 kept
    { auto c = make(); stage1_limit_per_parent(c, 1, true, parent.data());
      CHECK(c.size() == 2); }
    // n=2 tie-inclusive: parent0 {5,4,4} (2nd=4, ties) -> 3; parent1 {9,1,1} -> 3
    { auto c = make(); stage1_limit_per_parent(c, 2, true, parent.data());
      CHECK(c.size() == 6); }
    // n=2 exact: parent0 keeps {5,4}, parent1 keeps {9,1} -> 4
    { auto c = make(); stage1_limit_per_parent(c, 2, false, parent.data());
      CHECK(c.size() == 4); }
    // n=0 and null parent_index are no-ops
    { auto c = make(); stage1_limit_per_parent(c, 0, true, parent.data());
      CHECK(c.size() == 6); }
    { auto c = make(); stage1_limit_per_parent(c, 1, true, nullptr);
      CHECK(c.size() == 6); }
}

// Without tie inclusion the selector still has to name exactly which of the
// tying candidates survive: the ones that arrived first.
static void test_stage1_limit_per_parent_tie_order() {
    std::fprintf(stderr, "-- test_stage1_limit_per_parent_tie_order\n");
    std::vector<uint32_t> parent = {0, 0, 0};
    std::vector<Stage1Candidate> c{{0, 7}, {1, 7}, {2, 7}};
    stage1_limit_per_parent(c, 2, false, parent.data());
    CHECK_EQ(c.size(), 2u);
    if (c.size() != 2) return;
    CHECK_EQ(c[0].id, SeqId{0});
    CHECK_EQ(c[1].id, SeqId{1});
}

// Stage 1 in-total (L) limit.  The threshold is per query, so it aggregates
// over every ext_job of that query (both strands, every volume), and the prune
// must leave each job's candidate list in ascending SeqId order — the .kpx
// candidate-set decoder in Stage 2A depends on it.
static void test_limit_candidates_in_total() {
    std::fprintf(stderr, "-- test_limit_candidates_in_total\n");

    // Two ext_jobs for query 0 (one per strand) and one for query 1.
    const std::vector<ExtJob> ext_jobs = {
        {/*qi=*/0, /*vi=*/0, /*strand_idx=*/0},
        {/*qi=*/0, /*vi=*/0, /*strand_idx=*/1},
        {/*qi=*/1, /*vi=*/0, /*strand_idx=*/0},
    };
    auto make_states = []() {
        std::vector<JobState> st(3);
        st[0].candidates = {{10, 5}, {20, 3}, {30, 7}};
        st[1].candidates = {{5, 9}, {40, 1}};
        st[2].candidates = {{1, 4}, {2, 4}};
        return st;
    };

    // L == 0 is a no-op.
    {
        auto st = make_states();
        limit_candidates_in_total(st, ext_jobs, 0);
        CHECK_EQ(st[0].candidates.size(), 3u);
        CHECK_EQ(st[1].candidates.size(), 2u);
        CHECK_EQ(st[2].candidates.size(), 2u);
    }
    // L above a query's candidate count keeps that query untouched.
    {
        auto st = make_states();
        limit_candidates_in_total(st, ext_jobs, 100);
        CHECK_EQ(st[0].candidates.size(), 3u);
        CHECK_EQ(st[1].candidates.size(), 2u);
        CHECK_EQ(st[2].candidates.size(), 2u);
    }
    // L == 3 over query 0's five candidates (9,7,5,3,1): the threshold is 5,
    // so the first job keeps 7 and 5 and the second keeps 9.  Query 1 has
    // fewer than three candidates and is left alone.
    {
        auto st = make_states();
        limit_candidates_in_total(st, ext_jobs, 3);
        CHECK_EQ(st[0].candidates.size(), 2u);
        if (st[0].candidates.size() == 2) {
            // Ascending SeqId, as Stage 2A requires.
            CHECK_EQ(st[0].candidates[0].id, SeqId{10});
            CHECK_EQ(st[0].candidates[1].id, SeqId{30});
        }
        CHECK_EQ(st[1].candidates.size(), 1u);
        if (st[1].candidates.size() == 1)
            CHECK_EQ(st[1].candidates[0].id, SeqId{5});
        CHECK_EQ(st[2].candidates.size(), 2u);
    }
    // L == 1 is tie-inclusive: query 1's two tying candidates both survive,
    // while query 0 keeps only its single best.
    {
        auto st = make_states();
        limit_candidates_in_total(st, ext_jobs, 1);
        CHECK_EQ(st[0].candidates.size(), 0u);
        CHECK_EQ(st[1].candidates.size(), 1u);
        if (st[1].candidates.size() == 1)
            CHECK_EQ(st[1].candidates[0].id, SeqId{5});
        CHECK_EQ(st[2].candidates.size(), 2u);
    }
}

static std::string g_maxfreq_index_dir;

static bool build_maxfreq_index() {
    BlastDbReader db;
    if (!db.open(g_testdb_path)) return false;

    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 7;
    config.keep_tmp = true;

    std::string prefix = g_maxfreq_index_dir + "/test.00.07mer";
    if (!build_index<uint16_t>(db, config, prefix, 0, 1, "test", logger))
        return false;

    // Resolve fractional threshold: 0.5 * nseq
    uint32_t nseq = db.num_sequences();
    uint64_t freq_threshold = static_cast<uint64_t>(0.5 * nseq);
    if (freq_threshold == 0) freq_threshold = 1;

    // Apply cross-volume filter (single volume)
    std::string khx_path = g_maxfreq_index_dir + "/test.07mer.khx";
    return filter_volumes_cross_volume({prefix}, {prefix}, khx_path, 7, freq_threshold, 1, logger);
}

static void test_stage1_fractional_threshold() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_stage1_fractional_threshold\n");

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));
    CHECK(kpx.open(g_index_dir + "/test.00.07mer.kpx"));
    CHECK(ksx.open(g_index_dir + "/test.00.07mer.ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.min_stage1_score_frac = 0.5; // 50% of query k-mers

    // Search with fractional threshold (should produce results)
    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "test_query", qdata, 7,
        kix, kpx, ksx, filter, config, buf);

    // The fractional threshold resolves to ceil(Nqkmer * 0.5)
    // With a 100bp query, k=7: 94 k-mer positions, many distinct values
    // At 50%, threshold should be ~47 (for coverscore)
    // With exact match, FJ876973.1 should have high score but threshold is high
    // At minimum, the search should complete without error
    // (result may or may not contain hits depending on threshold)

    // Verify that with a lower fraction, we get more results
    SearchConfig config_low;
    config_low.stage1.min_stage1_score = 1;
    config_low.min_stage1_score_frac = 0.05; // 5% of query k-mers

    auto qdata_low = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config_low);
    auto result_low = search_volume<uint16_t>(
        "test_query", qdata_low, 7,
        kix, kpx, ksx, filter, config_low, buf);

    // Lower fraction should produce at least as many results as higher fraction
    CHECK(result_low.hits.size() >= result.hits.size());

    kix.close();
    kpx.close();
    ksx.close();
}

static void test_stage1_fractional_with_highfreq() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_stage1_fractional_with_highfreq\n");

    // Build index with max_freq_build to create .khx file
    CHECK(build_maxfreq_index());

    // Verify .khx was created (shared path, no volume index)
    std::string khx_path = g_maxfreq_index_dir + "/test.07mer.khx";
    KhxReader khx;
    CHECK(khx.open(khx_path));
    uint64_t excluded = khx.count_excluded();
    // With max_freq_build=0.5, some k-mers should be excluded
    // (exact count depends on data)
    std::fprintf(stderr, "   .khx excluded k-mers: %lu\n",
                 static_cast<unsigned long>(excluded));

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(g_maxfreq_index_dir + "/test.00.07mer.kix"));
    CHECK(kpx.open(g_maxfreq_index_dir + "/test.00.07mer.kpx"));
    CHECK(ksx.open(g_maxfreq_index_dir + "/test.00.07mer.ksx"));

    OidFilter filter;

    // Search with fractional threshold and khx awareness
    // Use P=0.3 so threshold stays positive after Nhighfreq subtraction:
    // ceil(Nqkmer * 0.3) - Nhighfreq > 0
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.min_stage1_score_frac = 0.3; // 30% of query k-mers

    auto qdata_with = preprocess_query<uint16_t>(g_query_seq, 7, &khx, config);
    auto result_with_khx = search_volume<uint16_t>(
        "test_query", qdata_with, 7,
        kix, kpx, ksx, filter, config, buf);

    auto qdata_without = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result_without_khx = search_volume<uint16_t>(
        "test_query", qdata_without, 7,
        kix, kpx, ksx, filter, config, buf);

    // With khx, the Nhighfreq subtraction makes the effective threshold lower
    // (because more k-mers are recognized as excluded), so we should get
    // at least as many results as without khx
    CHECK(result_with_khx.hits.size() >= result_without_khx.hits.size());

    khx.close();
    kix.close();
    kpx.close();
    ksx.close();
}

// Cover the HasFilter NTTP dispatch: a default-constructed OidFilter (kNone)
// must produce identical results to an explicit kInclude filter that lists
// every OID. The two paths exercise different compiled functions in
// stage1_filter_accumulate_impl.
static void test_stage1_filter_none_vs_pass_all() {
    std::fprintf(stderr, "-- test_stage1_filter_none_vs_pass_all\n");

    KixReader kix;
    KsxReader ksx;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));
    CHECK(ksx.open(g_index_dir + "/test.00.07mer.ksx"));

    std::vector<uint32_t> positions;
    std::vector<uint16_t> kmer_values;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        positions.push_back(pos);
        kmer_values.push_back(kmer);
    });

    OidFilter none_filter;  // kNone — HasFilter=false specialization

    // Build a kInclude filter listing every accession in this index volume.
    OidFilter pass_all;
    std::vector<std::string> all_accs;
    for (uint32_t oid = 0; oid < ksx.num_sequences(); oid++) {
        std::string_view acc = ksx.accession(oid);
        if (!acc.empty()) all_accs.emplace_back(acc);
    }
    pass_all.build(all_accs, ksx, OidFilterMode::kInclude);

    Stage1Config config;
    config.min_stage1_score = 1;

    Stage1Buffer buf_none, buf_all;
    auto a = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, none_filter, config, buf_none);
    auto b = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, pass_all, config, buf_all);

    CHECK(a.size() == b.size());
    auto sort_by_id = [](std::vector<Stage1Candidate>& v) {
        std::sort(v.begin(), v.end(),
                  [](const Stage1Candidate& x, const Stage1Candidate& y) { return x.id < y.id; });
    };
    sort_by_id(a);
    sort_by_id(b);
    for (size_t i = 0; i < a.size(); i++) {
        CHECK(a[i].id == b[i].id);
        CHECK(a[i].score == b[i].score);
    }

    kix.close();
    ksx.close();
}

// Reach the high-dirty branch of clear_dirty_typed by marking more than
// capacity/8 sids dirty. Verifies that after clearing, every score / last_pos
// slot has returned to the post-reset_all values (zero score, sentinel pos).
static void test_clear_dirty_bulk_reset() {
    std::fprintf(stderr, "-- test_clear_dirty_bulk_reset\n");

    Stage1Buffer buf;
    buf.width = Stage1Width::T32;
    buf.ensure_capacity(64);  // small capacity so 9 dirty entries trip the threshold

    using PosT   = Stage1WidthTraits<Stage1Width::T32>::PosT;
    using ScoreT = Stage1WidthTraits<Stage1Width::T32>::ScoreT;
    auto* scores   = score_ptr<Stage1Width::T32>(buf);
    auto* last_pos = last_pos_ptr<Stage1Width::T32>(buf);

    constexpr PosT sentinel = std::numeric_limits<PosT>::max();

    // 9 of 64 -> dirty.size()*8 = 72 > capacity (64): bulk path taken.
    for (uint32_t i = 0; i < 9; i++) {
        uint32_t sid = i * 5;
        scores[sid]   = static_cast<ScoreT>(i + 1);
        last_pos[sid] = static_cast<PosT>(i);
        buf.dirty.push_back(sid);
    }

    buf.clear_dirty_typed<Stage1Width::T32>();
    CHECK(buf.dirty.empty());
    for (uint32_t i = 0; i < buf.capacity; i++) {
        CHECK(scores[i] == ScoreT{0});
        CHECK(last_pos[i] == sentinel);
    }
}

// Drive the thread_local SeqIdDecoder through three consecutive runs with
// different query inputs on the same thread. If reset() forgot any field
// (decoded_, ctx_.count), the result of run 2 or 3 would differ from the
// result computed in isolation by a fresh call.
static void test_seq_id_decoder_reset_across_runs() {
    std::fprintf(stderr, "-- test_seq_id_decoder_reset_across_runs\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    auto scan = [&](const std::string& seq,
                    std::vector<uint32_t>& positions,
                    std::vector<uint16_t>& kmers) {
        KmerScanner<uint16_t> scanner(7);
        scanner.scan(seq.data(), seq.size(), [&](uint32_t pos, uint16_t kmer) {
            positions.push_back(pos);
            kmers.push_back(kmer);
        });
    };

    std::vector<uint32_t> p1, p2, p3;
    std::vector<uint16_t> k1, k2, k3;
    scan(g_query_seq, p1, k1);
    scan(g_query_seq.substr(0, 50), p2, k2);
    scan(g_query_seq.substr(20), p3, k3);

    OidFilter filter;
    Stage1Config config;
    config.min_stage1_score = 1;

    auto canonical = [&](const std::vector<uint32_t>& p,
                         const std::vector<uint16_t>& k) {
        Stage1Buffer buf;
        return stage1_filter(p.data(), k.data(), p.size(), kix, filter, config, buf);
    };

    // Reference values from independent buffers.
    auto r1 = canonical(p1, k1);
    auto r2 = canonical(p2, k2);
    auto r3 = canonical(p3, k3);

    // Now run the same three calls back-to-back on a shared buffer so the
    // thread_local SeqIdDecoder is exercised across runs.
    Stage1Buffer shared;
    auto s1 = stage1_filter(p1.data(), k1.data(), p1.size(), kix, filter, config, shared);
    auto s2 = stage1_filter(p2.data(), k2.data(), p2.size(), kix, filter, config, shared);
    auto s3 = stage1_filter(p3.data(), k3.data(), p3.size(), kix, filter, config, shared);

    auto compare = [](std::vector<Stage1Candidate> a, std::vector<Stage1Candidate> b) {
        CHECK(a.size() == b.size());
        auto sort_by_id = [](std::vector<Stage1Candidate>& v) {
            std::sort(v.begin(), v.end(),
                      [](const Stage1Candidate& x, const Stage1Candidate& y) { return x.id < y.id; });
        };
        sort_by_id(a);
        sort_by_id(b);
        for (size_t i = 0; i < a.size(); i++) {
            CHECK(a[i].id == b[i].id);
            CHECK(a[i].score == b[i].score);
        }
    };
    compare(r1, s1);
    compare(r2, s2);
    compare(r3, s3);

    kix.close();
}

// Reproduce the finish path by calling stage1_filter twice.  The second call
// must produce the same candidates as the first, proving that finish leaves no
// leftover score / last_pos / dirty state.
static void test_stage1_finish_no_side_effects() {
    std::fprintf(stderr, "-- test_stage1_finish_no_side_effects\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    std::vector<uint32_t> positions;
    std::vector<uint16_t> kmer_values;
    KmerScanner<uint16_t> scanner(7);
    scanner.scan(g_query_seq.data(), g_query_seq.size(), [&](uint32_t pos, uint16_t kmer) {
        positions.push_back(pos);
        kmer_values.push_back(kmer);
    });

    OidFilter filter;
    Stage1Config config;
    config.min_stage1_score = 1;

    Stage1Buffer buf;
    auto a = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, filter, config, buf);
    auto b = stage1_filter(positions.data(), kmer_values.data(), positions.size(), kix, filter, config, buf);

    CHECK(a.size() == b.size());
    auto sort_by_id = [](std::vector<Stage1Candidate>& v) {
        std::sort(v.begin(), v.end(),
                  [](const Stage1Candidate& x, const Stage1Candidate& y) { return x.id < y.id; });
    };
    sort_by_id(a);
    sort_by_id(b);
    for (size_t i = 0; i < a.size(); i++) {
        CHECK(a[i].id == b[i].id);
        CHECK(a[i].score == b[i].score);
    }
    kix.close();
}

static void test_adaptive_min_score() {
    std::fprintf(stderr, "-- test_adaptive_min_score\n");

    KixReader kix;
    CHECK(kix.open(g_index_dir + "/test.00.07mer.kix"));

    // Test adaptive min_score: min_score=0 with fractional threshold
    // effective_min_score should equal resolved Stage 1 threshold
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.min_score = 0; // adaptive
    config.min_stage1_score_frac = 0.3; // 30% of query k-mers

    auto qdata = preprocess_query<uint16_t>(
        g_query_seq, 7, nullptr, config);

    // Adaptive: effective_min_score should equal resolved threshold
    CHECK(qdata.effective_min_score_fwd == qdata.resolved_threshold_fwd);
    CHECK(qdata.effective_min_score_rc == qdata.resolved_threshold_rc);
    CHECK(qdata.resolved_threshold_fwd > 0);
    CHECK(qdata.resolved_threshold_rc > 0);

    // Test explicit min_score overrides adaptive
    SearchConfig config_explicit;
    config_explicit.stage1.min_stage1_score = 1;
    config_explicit.stage2.min_score = 5; // explicit
    config_explicit.min_stage1_score_frac = 0.3;

    auto qdata_explicit = preprocess_query<uint16_t>(
        g_query_seq, 7, nullptr, config_explicit);

    CHECK(qdata_explicit.effective_min_score_fwd == 5);
    CHECK(qdata_explicit.effective_min_score_rc == 5);

    kix.close();
}

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_index_dir = test_tmpdir("/tmp/ikafssn_stage1_test");
    g_maxfreq_index_dir = test_tmpdir("/tmp/ikafssn_stage1_maxfreq_test");
    std::filesystem::create_directories(g_index_dir);
    std::filesystem::create_directories(g_maxfreq_index_dir);

    // Extract 100bp query from FJ876973.1 at runtime
    {
        BlastDbReader db;
        CHECK(db.open(g_testdb_path));
        g_fj_oid = find_oid_by_accession(db, ACC_FJ);
        CHECK(g_fj_oid != UINT32_MAX);
        std::string full_seq = db.get_subsequence(g_fj_oid, 0, db.seq_length(g_fj_oid) - 1);
        CHECK(full_seq.size() >= 200);
        g_query_seq = full_seq.substr(100, 100);
    }

    CHECK(build_test_index());

    test_stage1_basic();
    test_stage1_min_score();
    test_stage1_with_oid_filter();
    test_stage1_limit_topk();
    test_stage1_limit_per_parent();
    test_stage1_limit_per_parent_tie_order();
    test_limit_candidates_in_total();
    test_stage1_fractional_threshold();
    test_stage1_fractional_with_highfreq();
    test_clear_dirty_bulk_reset();
    test_stage1_finish_no_side_effects();
    test_seq_id_decoder_reset_across_runs();
    test_stage1_filter_none_vs_pass_all();
    test_adaptive_min_score();

    std::filesystem::remove_all(g_index_dir);
    std::filesystem::remove_all(g_maxfreq_index_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
