#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <mutex>
#include <string>
#include <vector>

#include <tbb/parallel_for_each.h>
#include <tbb/task_arena.h>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_multivol_a_path;
static std::string g_multivol_b_path;
static std::string g_test_dir;

// Runtime-extracted queries from the multi-volume DBs
static std::string g_query_fj;  // 100bp from FJ876973.1 (in multivol_a)
static std::string g_query_gq;  // 100bp from GQ912721.1 (in multivol_b)

// Build indexes for two separate BLAST DBs as volume 0 and volume 1
// in the same output directory. This simulates a multi-volume DB.
static void build_multivolume_index(int k, const std::string& db_base) {
    Logger logger(Logger::kError);

    // Volume 0: multivol_a
    {
        BlastDbReader db;
        CHECK(db.open(g_multivol_a_path));

        IndexBuilderConfig config;
        config.k = k;

        char kk_str[8];
        std::snprintf(kk_str, sizeof(kk_str), "%02d", k);
        std::string prefix = g_test_dir + "/" + db_base + ".00." +
                             std::string(kk_str) + "mer";

        if (k < K_TYPE_THRESHOLD) {
            CHECK(build_index<uint16_t>(db, config, prefix, 0, 2, db_base, logger));
        } else {
            CHECK(build_index<uint32_t>(db, config, prefix, 0, 2, db_base, logger));
        }
    }

    // Volume 1: multivol_b
    {
        BlastDbReader db;
        CHECK(db.open(g_multivol_b_path));

        IndexBuilderConfig config;
        config.k = k;

        char kk_str[8];
        std::snprintf(kk_str, sizeof(kk_str), "%02d", k);
        std::string prefix = g_test_dir + "/" + db_base + ".01." +
                             std::string(kk_str) + "mer";

        if (k < K_TYPE_THRESHOLD) {
            CHECK(build_index<uint16_t>(db, config, prefix, 1, 2, db_base, logger));
        } else {
            CHECK(build_index<uint32_t>(db, config, prefix, 1, 2, db_base, logger));
        }
    }
}

// Helper: search all volumes sequentially (reference implementation)
static std::vector<OutputHit> search_sequential(
        const std::string& db_base, int k,
        const std::vector<FastaRecord>& queries,
        const SearchConfig& config) {
    char kk_str[8];
    std::snprintf(kk_str, sizeof(kk_str), "%02d", k);

    std::vector<OutputHit> all_hits;

    for (int vi = 0; vi < 2; vi++) {
        char vol_str[8];
        std::snprintf(vol_str, sizeof(vol_str), "%02d", vi);
        std::string prefix = g_test_dir + "/" + db_base + "." +
                             std::string(vol_str) + "." +
                             std::string(kk_str) + "mer";

        KixReader kix;
        KpxReader kpx;
        KsxReader ksx;
        CHECK(kix.open(prefix + ".kix"));
        CHECK(kpx.open(prefix + ".kpx"));
        CHECK(ksx.open(prefix + ".ksx"));

        OidFilter filter;

        for (const auto& query : queries) {
            SearchResult sr;
            if (k < K_TYPE_THRESHOLD) {
                sr = search_volume<uint16_t>(
                    query.id, query.sequence, k, kix, kpx, ksx, filter, config);
            } else {
                sr = search_volume<uint32_t>(
                    query.id, query.sequence, k, kix, kpx, ksx, filter, config);
            }

            for (const auto& cr : sr.hits) {
                OutputHit oh;
                oh.query_id = sr.query_id;
                oh.accession = std::string(ksx.accession(cr.seq_id));
                oh.strand = cr.is_reverse ? '-' : '+';
                oh.q_start = cr.q_start;
                oh.q_end = cr.q_end;
                oh.s_start = cr.s_start;
                oh.s_end = cr.s_end;
                oh.score = cr.score;
                oh.volume = static_cast<uint16_t>(vi);
                all_hits.push_back(oh);
            }
        }

        kix.close();
        kpx.close();
        ksx.close();
    }

    // Sort by (query_id, score desc, accession, volume)
    std::sort(all_hits.begin(), all_hits.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  if (a.query_id != b.query_id) return a.query_id < b.query_id;
                  if (a.score != b.score) return a.score > b.score;
                  if (a.accession != b.accession) return a.accession < b.accession;
                  return a.volume < b.volume;
              });

    return all_hits;
}

// Helper: search all volumes in parallel using TBB
static std::vector<OutputHit> search_parallel(
        const std::string& db_base, int k,
        const std::vector<FastaRecord>& queries,
        const SearchConfig& config,
        int num_threads) {
    char kk_str[8];
    std::snprintf(kk_str, sizeof(kk_str), "%02d", k);

    // Pre-open volumes
    struct VolumeData {
        KixReader kix;
        KpxReader kpx;
        KsxReader ksx;
        uint16_t volume_index;
    };

    std::vector<VolumeData> vols(2);
    for (int vi = 0; vi < 2; vi++) {
        char vol_str[8];
        std::snprintf(vol_str, sizeof(vol_str), "%02d", vi);
        std::string prefix = g_test_dir + "/" + db_base + "." +
                             std::string(vol_str) + "." +
                             std::string(kk_str) + "mer";
        CHECK(vols[vi].kix.open(prefix + ".kix"));
        CHECK(vols[vi].kpx.open(prefix + ".kpx"));
        CHECK(vols[vi].ksx.open(prefix + ".ksx"));
        vols[vi].volume_index = static_cast<uint16_t>(vi);
    }

    // Build job list
    struct Job { size_t query_idx; size_t volume_idx; };
    std::vector<Job> jobs;
    for (size_t qi = 0; qi < queries.size(); qi++) {
        for (size_t vi = 0; vi < 2; vi++) {
            jobs.push_back({qi, vi});
        }
    }

    std::vector<OutputHit> all_hits;
    std::mutex mutex;

    tbb::task_arena arena(num_threads);
    arena.execute([&] {
        tbb::parallel_for_each(jobs.begin(), jobs.end(),
            [&](const Job& job) {
                const auto& query = queries[job.query_idx];
                const auto& vd = vols[job.volume_idx];
                OidFilter filter;

                SearchResult sr;
                if (k < K_TYPE_THRESHOLD) {
                    sr = search_volume<uint16_t>(
                        query.id, query.sequence, k,
                        vd.kix, vd.kpx, vd.ksx, filter, config);
                } else {
                    sr = search_volume<uint32_t>(
                        query.id, query.sequence, k,
                        vd.kix, vd.kpx, vd.ksx, filter, config);
                }

                if (!sr.hits.empty()) {
                    std::vector<OutputHit> local;
                    for (const auto& cr : sr.hits) {
                        OutputHit oh;
                        oh.query_id = sr.query_id;
                        oh.accession = std::string(vd.ksx.accession(cr.seq_id));
                        oh.strand = cr.is_reverse ? '-' : '+';
                        oh.q_start = cr.q_start;
                        oh.q_end = cr.q_end;
                        oh.s_start = cr.s_start;
                        oh.s_end = cr.s_end;
                        oh.score = cr.score;
                        oh.volume = vd.volume_index;
                        local.push_back(oh);
                    }
                    std::lock_guard<std::mutex> lock(mutex);
                    all_hits.insert(all_hits.end(), local.begin(), local.end());
                }
            });
    });

    // Sort by (query_id, score desc, accession, volume)
    std::sort(all_hits.begin(), all_hits.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  if (a.query_id != b.query_id) return a.query_id < b.query_id;
                  if (a.score != b.score) return a.score > b.score;
                  if (a.accession != b.accession) return a.accession < b.accession;
                  return a.volume < b.volume;
              });

    for (auto& vd : vols) {
        vd.kix.close();
        vd.kpx.close();
        vd.ksx.close();
    }

    return all_hits;
}

static void test_multivolume_search() {
    std::fprintf(stderr, "-- test_multivolume_search\n");

    const int k = 7;
    const std::string db_base = "mvtest";

    build_multivolume_index(k, db_base);

    // Queries extracted from sequences in both volumes
    std::vector<FastaRecord> queries = {
        {"query_fj", g_query_fj}
    };

    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    auto results = search_sequential(db_base, k, queries, config);

    // Should find hits from volume 0 (mvseq1 has ACGT pattern)
    bool found_vol0 = false;
    // Should find hits from volume 1 (mvseq5 has ACGT pattern)
    bool found_vol1 = false;

    for (const auto& h : results) {
        if (h.volume == 0) found_vol0 = true;
        if (h.volume == 1) found_vol1 = true;
    }

    CHECK(found_vol0);
    CHECK(found_vol1);
    CHECK(!results.empty());
}

static void test_parallel_equals_sequential() {
    std::fprintf(stderr, "-- test_parallel_equals_sequential\n");

    const int k = 7;
    const std::string db_base = "mvtest";

    // Index already built by test_multivolume_search
    std::vector<FastaRecord> queries = {
        {"query_fj", g_query_fj},
        {"query_gq", g_query_gq}
    };

    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    auto seq_results = search_sequential(db_base, k, queries, config);
    auto par_results = search_parallel(db_base, k, queries, config, 4);

    // Both should produce the same number of results
    CHECK_EQ(seq_results.size(), par_results.size());

    // Verify each result matches
    for (size_t i = 0; i < seq_results.size() && i < par_results.size(); i++) {
        CHECK(seq_results[i].query_id == par_results[i].query_id);
        CHECK(seq_results[i].accession == par_results[i].accession);
        CHECK(seq_results[i].strand == par_results[i].strand);
        CHECK_EQ(seq_results[i].score, par_results[i].score);
        CHECK_EQ(seq_results[i].volume, par_results[i].volume);
    }
}

static void test_result_merge_ordering() {
    std::fprintf(stderr, "-- test_result_merge_ordering\n");

    const int k = 7;
    const std::string db_base = "mvtest";

    std::vector<FastaRecord> queries = {
        {"query_fj", g_query_fj}
    };

    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    auto results = search_parallel(db_base, k, queries, config, 2);

    // Verify results are sorted by score descending within each query
    for (size_t i = 1; i < results.size(); i++) {
        if (results[i].query_id == results[i - 1].query_id) {
            CHECK(results[i].score <= results[i - 1].score);
        }
    }
}

static void test_parallel_counting_pass() {
    std::fprintf(stderr, "-- test_parallel_counting_pass\n");

    // Build index with k=7 using parallel counting (threads=2)
    // and verify it produces the same result as threads=1
    Logger logger(Logger::kError);
    BlastDbReader db;
    CHECK(db.open(g_multivol_a_path));

    // Build with 1 thread
    {
        IndexBuilderConfig config;
        config.k = 7;
        config.threads = 1;

        std::string prefix = g_test_dir + "/pcnt_st.00.07mer";
        CHECK(build_index<uint16_t>(db, config, prefix, 0, 1, "pcnt", logger));
    }

    // Build with 2 threads (need to re-open DB as it may have internal state)
    BlastDbReader db2;
    CHECK(db2.open(g_multivol_a_path));
    {
        IndexBuilderConfig config;
        config.k = 7;
        config.threads = 2;

        std::string prefix = g_test_dir + "/pcnt_mt.00.07mer";
        CHECK(build_index<uint16_t>(db2, config, prefix, 0, 1, "pcnt", logger));
    }

    // Compare kix files: offsets and counts should be identical
    KixReader kix_st, kix_mt;
    CHECK(kix_st.open(g_test_dir + "/pcnt_st.00.07mer.kix"));
    CHECK(kix_mt.open(g_test_dir + "/pcnt_mt.00.07mer.kix"));

    CHECK_EQ(kix_st.table_size(), kix_mt.table_size());
    CHECK_EQ(kix_st.total_postings(), kix_mt.total_postings());

    // Compare counts arrays
    bool counts_match = true;
    for (uint64_t i = 0; i < kix_st.table_size(); i++) {
        if (kix_st.counts()[i] != kix_mt.counts()[i]) {
            counts_match = false;
            std::fprintf(stderr, "  counts mismatch at kmer %lu: st=%u mt=%u\n",
                         (unsigned long)i, kix_st.counts()[i], kix_mt.counts()[i]);
            break;
        }
    }
    CHECK(counts_match);

    // Search both and verify same results
    KpxReader kpx_st, kpx_mt;
    KsxReader ksx_st, ksx_mt;
    CHECK(kpx_st.open(g_test_dir + "/pcnt_st.00.07mer.kpx"));
    CHECK(kpx_mt.open(g_test_dir + "/pcnt_mt.00.07mer.kpx"));
    CHECK(ksx_st.open(g_test_dir + "/pcnt_st.00.07mer.ksx"));
    CHECK(ksx_mt.open(g_test_dir + "/pcnt_mt.00.07mer.ksx"));

    OidFilter filter;
    SearchConfig sconfig;
    sconfig.stage1.max_freq = 100000;
    sconfig.stage1.stage1_topn = 100;
    sconfig.stage1.min_stage1_score = 1;
    sconfig.stage2.max_gap = 100;
    sconfig.stage2.min_diag_hits = 1;
    sconfig.stage2.min_score = 2;
    sconfig.num_results = 50;

    auto sr_st = search_volume<uint16_t>(
        "q", g_query_fj, 7, kix_st, kpx_st, ksx_st, filter, sconfig);
    auto sr_mt = search_volume<uint16_t>(
        "q", g_query_fj, 7, kix_mt, kpx_mt, ksx_mt, filter, sconfig);

    CHECK_EQ(sr_st.hits.size(), sr_mt.hits.size());
    for (size_t i = 0; i < sr_st.hits.size() && i < sr_mt.hits.size(); i++) {
        CHECK_EQ(sr_st.hits[i].seq_id, sr_mt.hits[i].seq_id);
        CHECK_EQ(sr_st.hits[i].score, sr_mt.hits[i].score);
    }

    kix_st.close(); kix_mt.close();
    kpx_st.close(); kpx_mt.close();
    ksx_st.close(); ksx_mt.close();
}

static void test_multivolume_k9() {
    std::fprintf(stderr, "-- test_multivolume_k9\n");

    const int k = 9;
    const std::string db_base = "mvtest9";

    build_multivolume_index(k, db_base);

    std::vector<FastaRecord> queries = {
        {"query_fj9", g_query_fj}
    };

    SearchConfig config;
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    auto seq_results = search_sequential(db_base, k, queries, config);
    auto par_results = search_parallel(db_base, k, queries, config, 2);

    CHECK_EQ(seq_results.size(), par_results.size());

    for (size_t i = 0; i < seq_results.size() && i < par_results.size(); i++) {
        CHECK(seq_results[i].query_id == par_results[i].query_id);
        CHECK_EQ(seq_results[i].score, par_results[i].score);
    }
}

int main() {
    check_ssu_available();
    check_derived_data_ready();

    g_multivol_a_path = multivol_a_prefix();
    g_multivol_b_path = multivol_b_prefix();
    g_test_dir = "/tmp/ikafssn_multivolume_test";
    std::filesystem::create_directories(g_test_dir);

    // Extract queries from the multi-volume DBs at runtime
    {
        BlastDbReader db_a;
        CHECK(db_a.open(g_multivol_a_path));
        uint32_t oid_fj = find_oid_by_accession(db_a, ACC_FJ);
        CHECK(oid_fj != UINT32_MAX);
        std::string seq = db_a.get_sequence(oid_fj);
        CHECK(seq.size() >= 200);
        g_query_fj = seq.substr(100, 100);
    }
    {
        BlastDbReader db_b;
        CHECK(db_b.open(g_multivol_b_path));
        uint32_t oid_gq = find_oid_by_accession(db_b, ACC_GQ);
        CHECK(oid_gq != UINT32_MAX);
        std::string seq = db_b.get_sequence(oid_gq);
        CHECK(seq.size() >= 200);
        g_query_gq = seq.substr(50, 100);
    }

    test_multivolume_search();
    test_parallel_equals_sequential();
    test_result_merge_ordering();
    test_parallel_counting_pass();
    test_multivolume_k9();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
