#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "io/volume_discovery.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"
#include "search/search_orchestrator.hpp"
#include "search/stage1_filter.hpp"
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

// Volume basenames for the two test BLAST DBs
static const char* VOL_BASENAME_A = "ssu_multivol_a";
static const char* VOL_BASENAME_B = "ssu_multivol_b";

// Build indexes for two separate BLAST DBs as volume 0 and volume 1
// in the same output directory, using volume basenames (new naming convention).
// Also writes a .kvx manifest.
static void build_multivolume_index(int k, const std::string& db_base) {
    Logger logger(Logger::kError);

    char kk_str[8];
    std::snprintf(kk_str, sizeof(kk_str), "%02d", k);

    auto vol_stem = [&](const std::string& base) {
        IndexBuilderConfig config;
        config.k = k;
        return index_file_stem(g_test_dir, base, config.k,
                               config.t, config.template_type,
                               config.min_seq_length,
                               config.min_length_split,
                               config.overlap_length,
                               /*max_freq_build=*/1,
                               /*max_degen_expand=*/0);
    };

    // Volume 0: multivol_a (using volume basename)
    {
        BlastDbReader db;
        CHECK(db.open(g_multivol_a_path));

        IndexBuilderConfig config;
        config.k = k;

        std::string prefix = vol_stem(VOL_BASENAME_A);

        if (k < K_TYPE_THRESHOLD) {
            CHECK(build_index<uint16_t>(db, config, prefix, 0, 2, db_base, logger));
        } else {
            CHECK(build_index<uint32_t>(db, config, prefix, 0, 2, db_base, logger));
        }
    }

    // Volume 1: multivol_b (using volume basename)
    {
        BlastDbReader db;
        CHECK(db.open(g_multivol_b_path));

        IndexBuilderConfig config;
        config.k = k;

        std::string prefix = vol_stem(VOL_BASENAME_B);

        if (k < K_TYPE_THRESHOLD) {
            CHECK(build_index<uint16_t>(db, config, prefix, 1, 2, db_base, logger));
        } else {
            CHECK(build_index<uint32_t>(db, config, prefix, 1, 2, db_base, logger));
        }
    }

    // Write .kvx manifest
    {
        std::string kvx_path = vol_stem(db_base) + ".kvx";
        FILE* fp = std::fopen(kvx_path.c_str(), "w");
        CHECK(fp != nullptr);
        std::fprintf(fp, "#\n# ikafssn index volume manifest\n#\n");
        std::fprintf(fp, "TITLE %s\n", db_base.c_str());
        std::fprintf(fp, "DBLIST \"%s\" \"%s\"\n", VOL_BASENAME_A, VOL_BASENAME_B);
        std::fclose(fp);
    }
}

// Helper: search all volumes sequentially (reference implementation)
static std::vector<OutputHit> search_sequential(
        const std::string& db_base, int k,
        const std::vector<FastaRecord>& queries,
        const SearchConfig& config) {
    auto vols = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(vols.size(), 2u);
    std::vector<OutputHit> all_hits;

    for (size_t vi = 0; vi < 2; vi++) {
        const std::string& kix_path = vols[vi].kix_path;
        const std::string& kpx_path = vols[vi].kpx_path;
        const std::string& ksx_path = vols[vi].ksx_path;

        KixReader kix;
        KpxReader kpx;
        KsxReader ksx;
        CHECK(kix.open(kix_path));
        CHECK(kpx.open(kpx_path));
        CHECK(ksx.open(ksx_path));

        OidFilter filter;

        for (const auto& query : queries) {
            SearchResult sr;
            if (k < K_TYPE_THRESHOLD) {
    Stage1Buffer buf;
                auto qdata = preprocess_query<uint16_t>(query.sequence, k, nullptr, config);
                sr = search_volume<uint16_t>(
                    query.id, qdata, k, kix, kpx, ksx, filter, config, buf);
            } else {
    Stage1Buffer buf;
                auto qdata = preprocess_query<uint32_t>(query.sequence, k, nullptr, config);
                sr = search_volume<uint32_t>(
                    query.id, qdata, k, kix, kpx, ksx, filter, config, buf);
            }

            for (const auto& cr : sr.hits) {
                OutputHit oh;
                oh.qseqid = sr.query_id;
                oh.sseqid = std::string(ksx.accession(cr.seq_id));
                oh.sstrand = cr.is_reverse ? '-' : '+';
                oh.qstart = cr.q_start;
                oh.qend = cr.q_end;
                oh.sstart = cr.s_start;
                oh.send = cr.s_end;
                oh.chainscore = cr.chainscore;
                oh.volume = static_cast<uint16_t>(vi);
                all_hits.push_back(oh);
            }
        }

        kix.close();
        kpx.close();
        ksx.close();
    }

    // Sort by (qseqid, chainscore desc, sseqid, volume)
    std::sort(all_hits.begin(), all_hits.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  if (a.qseqid != b.qseqid) return a.qseqid < b.qseqid;
                  if (a.chainscore != b.chainscore) return a.chainscore > b.chainscore;
                  if (a.sseqid != b.sseqid) return a.sseqid < b.sseqid;
                  return a.volume < b.volume;
              });

    return all_hits;
}

// Helper: search all volumes in parallel using TBB
static std::vector<OutputHit> search_parallel(
        const std::string& db_base, int k,
        const std::vector<FastaRecord>& queries,
        const SearchConfig& config,
        int nthread) {
    auto discovered = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(discovered.size(), 2u);

    // Pre-open volumes
    struct VolumeData {
        KixReader kix;
        KpxReader kpx;
        KsxReader ksx;
        uint16_t volume_index;
    };

    std::vector<VolumeData> vols(2);
    for (size_t vi = 0; vi < 2; vi++) {
        CHECK(vols[vi].kix.open(discovered[vi].kix_path));
        CHECK(vols[vi].kpx.open(discovered[vi].kpx_path));
        CHECK(vols[vi].ksx.open(discovered[vi].ksx_path));
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

    tbb::task_arena arena(nthread);
    arena.execute([&] {
        tbb::parallel_for_each(jobs.begin(), jobs.end(),
            [&](const Job& job) {
                const auto& query = queries[job.query_idx];
                const auto& vd = vols[job.volume_idx];
                OidFilter filter;

                SearchResult sr;
                if (k < K_TYPE_THRESHOLD) {
    Stage1Buffer buf;
                    auto qdata = preprocess_query<uint16_t>(query.sequence, k, nullptr, config);
                    sr = search_volume<uint16_t>(
                        query.id, qdata, k,
                        vd.kix, vd.kpx, vd.ksx, filter, config, buf);
                } else {
    Stage1Buffer buf;
                    auto qdata = preprocess_query<uint32_t>(query.sequence, k, nullptr, config);
                    sr = search_volume<uint32_t>(
                        query.id, qdata, k,
                        vd.kix, vd.kpx, vd.ksx, filter, config, buf);
                }

                if (!sr.hits.empty()) {
                    std::vector<OutputHit> local;
                    for (const auto& cr : sr.hits) {
                        OutputHit oh;
                        oh.qseqid = sr.query_id;
                        oh.sseqid = std::string(vd.ksx.accession(cr.seq_id));
                        oh.sstrand = cr.is_reverse ? '-' : '+';
                        oh.qstart = cr.q_start;
                        oh.qend = cr.q_end;
                        oh.sstart = cr.s_start;
                        oh.send = cr.s_end;
                        oh.chainscore = cr.chainscore;
                        oh.volume = vd.volume_index;
                        local.push_back(oh);
                    }
                    std::lock_guard<std::mutex> lock(mutex);
                    all_hits.insert(all_hits.end(), local.begin(), local.end());
                }
            });
    });

    // Sort by (qseqid, chainscore desc, sseqid, volume)
    std::sort(all_hits.begin(), all_hits.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  if (a.qseqid != b.qseqid) return a.qseqid < b.qseqid;
                  if (a.chainscore != b.chainscore) return a.chainscore > b.chainscore;
                  if (a.sseqid != b.sseqid) return a.sseqid < b.sseqid;
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

    auto seq_results = search_sequential(db_base, k, queries, config);
    auto par_results = search_parallel(db_base, k, queries, config, 4);

    // Both should produce the same number of results
    CHECK_EQ(seq_results.size(), par_results.size());

    // Verify each result matches
    for (size_t i = 0; i < seq_results.size() && i < par_results.size(); i++) {
        CHECK(seq_results[i].qseqid == par_results[i].qseqid);
        CHECK(seq_results[i].sseqid == par_results[i].sseqid);
        CHECK(seq_results[i].sstrand == par_results[i].sstrand);
        CHECK_EQ(seq_results[i].chainscore, par_results[i].chainscore);
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

    auto results = search_parallel(db_base, k, queries, config, 2);

    // Verify results are sorted by score descending within each query
    for (size_t i = 1; i < results.size(); i++) {
        if (results[i].qseqid == results[i - 1].qseqid) {
            CHECK(results[i].chainscore <= results[i - 1].chainscore);
        }
    }
}

static void test_parallel_counting_pass() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_parallel_counting_pass\n");

    // Build index with k=7 using parallel counting (threads=2)
    // and verify it produces the same result as threads=1
    Logger logger(Logger::kError);
    BlastDbReader db;
    CHECK(db.open(g_multivol_a_path));

    auto pcnt_stem = [&](const std::string& base) {
        IndexBuilderConfig cfg;
        cfg.k = 7;
        return index_file_stem(g_test_dir, base, cfg.k,
                               cfg.t, cfg.template_type,
                               cfg.min_seq_length,
                               cfg.min_length_split,
                               cfg.overlap_length,
                               /*max_freq_build=*/1,
                               /*max_degen_expand=*/0);
    };
    std::string st_prefix = pcnt_stem("pcnt_st.00");
    std::string mt_prefix = pcnt_stem("pcnt_mt.00");

    // Build with 1 thread
    {
        IndexBuilderConfig config;
        config.k = 7;
        config.threads = 1;

        CHECK(build_index<uint16_t>(db, config, st_prefix, 0, 1, "pcnt", logger));
    }

    // Build with 2 threads (need to re-open DB as it may have internal state)
    BlastDbReader db2;
    CHECK(db2.open(g_multivol_a_path));
    {
        IndexBuilderConfig config;
        config.k = 7;
        config.threads = 2;

        CHECK(build_index<uint16_t>(db2, config, mt_prefix, 0, 1, "pcnt", logger));
    }

    // Compare kix files: offsets and counts should be identical
    KixReader kix_st, kix_mt;
    CHECK(kix_st.open(st_prefix + ".kix"));
    CHECK(kix_mt.open(mt_prefix + ".kix"));

    CHECK_EQ(kix_st.table_size(), kix_mt.table_size());
    CHECK_EQ(kix_st.total_distinct_postings(), kix_mt.total_distinct_postings());

    // Compare posting list byte lengths (replaces counts comparison)
    bool counts_match = true;
    for (uint32_t i = 0; i < kix_st.table_size(); i++) {
        if (kix_st.posting_list_byte_length(i) != kix_mt.posting_list_byte_length(i)) {
            counts_match = false;
            std::fprintf(stderr, "  byte_length mismatch at kmer %u: st=%lu mt=%lu\n",
                         i, (unsigned long)kix_st.posting_list_byte_length(i),
                         (unsigned long)kix_mt.posting_list_byte_length(i));
            break;
        }
    }
    CHECK(counts_match);

    // Search both and verify same results
    KpxReader kpx_st, kpx_mt;
    KsxReader ksx_st, ksx_mt;
    CHECK(kpx_st.open(st_prefix + ".kpx"));
    CHECK(kpx_mt.open(mt_prefix + ".kpx"));
    CHECK(ksx_st.open(st_prefix + ".ksx"));
    CHECK(ksx_mt.open(mt_prefix + ".ksx"));

    OidFilter filter;
    SearchConfig sconfig;
    sconfig.stage1.stage1_topn = 100;
    sconfig.stage1.min_stage1_score = 1;
    sconfig.stage2.max_gap = 100;
    sconfig.stage2.min_nhit_diag = 1;
    sconfig.stage2.min_score = 2;
    sconfig.nresult = 50;

    auto qdata_st = preprocess_query<uint16_t>(g_query_fj, 7, nullptr, sconfig);
    auto sr_st = search_volume<uint16_t>(
        "q", qdata_st, 7, kix_st, kpx_st, ksx_st, filter, sconfig, buf);
    auto qdata_mt = preprocess_query<uint16_t>(g_query_fj, 7, nullptr, sconfig);
    auto sr_mt = search_volume<uint16_t>(
        "q", qdata_mt, 7, kix_mt, kpx_mt, ksx_mt, filter, sconfig, buf);

    CHECK_EQ(sr_st.hits.size(), sr_mt.hits.size());
    for (size_t i = 0; i < sr_st.hits.size() && i < sr_mt.hits.size(); i++) {
        CHECK_EQ(sr_st.hits[i].seq_id, sr_mt.hits[i].seq_id);
        CHECK_EQ(sr_st.hits[i].chainscore, sr_mt.hits[i].chainscore);
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
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 50;

    auto seq_results = search_sequential(db_base, k, queries, config);
    auto par_results = search_parallel(db_base, k, queries, config, 2);

    CHECK_EQ(seq_results.size(), par_results.size());

    for (size_t i = 0; i < seq_results.size() && i < par_results.size(); i++) {
        CHECK(seq_results[i].qseqid == par_results[i].qseqid);
        CHECK_EQ(seq_results[i].chainscore, par_results[i].chainscore);
    }
}

// Drive run_search<>() in mode 2 with two different posting_budget
// values: one large enough that Stage 1 + Stage 2 share a single batch,
// one tiny enough that each volume goes through a Stage 2 Tier1
// fallback (range planner emits per-volume single-vol batches). With
// the range-WILLNEED planner the per-volume kpx ranges must match the
// single-batch result.
static void test_mode2_batched_equals_single() {
    std::fprintf(stderr, "-- test_mode2_batched_equals_single (Phase 5)\n");
    const int k = 7;
    const std::string db_base = "mvtest";

    auto discovered = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(discovered.size(), 2u);

    std::vector<KsxReader> ksxs(2);
    std::vector<VolumeMeta> volumes_cod(2);
    for (size_t vi = 0; vi < 2; vi++) {
        const auto& vf = discovered[vi];
        CHECK(ksxs[vi].open(vf.ksx_path));
        KixReader kix_probe;
        KpxReader kpx_probe;
        CHECK(kix_probe.open(vf.kix_path));
        CHECK(kpx_probe.open(vf.kpx_path));
        volumes_cod[vi].files            = vf;
        volumes_cod[vi].kix_posting_size = kix_probe.posting_file_size();
        volumes_cod[vi].kpx_posting_size = kpx_probe.posting_file_size();
        volumes_cod[vi].kix_full_size    = kix_probe.willneed_size_full();
        volumes_cod[vi].kpx_full_size    = kpx_probe.willneed_size_full();
        volumes_cod[vi].volume_index     = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences    = kix_probe.num_sequences();
        kix_probe.close();
        kpx_probe.close();
    }

    SearchConfig config;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.nresult = 0;
    config.mode = 2;

    auto qdata = preprocess_query<uint16_t>(g_query_fj, k, nullptr, config);
    std::vector<QueryBundle<uint16_t>> bundles(1);
    bundles[0].query_id = &g_query_fj;
    bundles[0].qdata_primary = &qdata;
    std::vector<uint8_t> skip_reason(1, 0);

    Logger logger(Logger::kError);
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    auto run_with_budget = [&](uint64_t budget) {
        RunSearchInputs<uint16_t> in;
        in.volumes_cod      = volumes_cod;
        in.ksx_per_volume.resize(2);
        in.oid_filters.resize(2);
        for (size_t vi = 0; vi < 2; vi++) in.ksx_per_volume[vi] = &ksxs[vi];
        in.queries           = &bundles;
        in.query_skip_reason = &skip_reason;
        in.config            = config;
        in.both_mode         = false;
        in.k                 = k;
        in.nthread           = 2;
        in.posting_budget    = budget;
        in.logger            = &logger;
        in.max_num_seqs      = max_num_seqs;
        in.tier              = Stage1Tier::T32;
        return run_search<uint16_t>(in);
    };

    uint64_t big_budget = 0;
    for (auto& v : volumes_cod) big_budget += v.kix_full_size + v.kpx_full_size;
    big_budget = std::max<uint64_t>(big_budget * 2, 1u << 20);
    auto big = run_with_budget(big_budget);

    auto small = run_with_budget(1);  // forces single-volume Tier1 batches

    auto key = [](const OrchestratorHit& h) {
        return std::tuple<size_t, uint16_t, uint16_t, uint32_t, int32_t,
                          uint32_t, uint32_t, uint32_t, uint32_t>(
            h.query_idx, h.volume_idx, h.volume_index,
            h.cr.seq_id, h.cr.chainscore,
            h.cr.q_start, h.cr.q_end, h.cr.s_start, h.cr.s_end);
    };
    auto sorted_keys = [&](const std::vector<OrchestratorHit>& v) {
        std::vector<decltype(key(v.front()))> ks;
        ks.reserve(v.size());
        for (auto& h : v) ks.push_back(key(h));
        std::sort(ks.begin(), ks.end());
        return ks;
    };

    CHECK_EQ(big.size(), small.size());
    CHECK(sorted_keys(big) == sorted_keys(small));

    for (auto& ksx : ksxs) ksx.close();
}

// Drive run_search<>() directly with a multi-volume index, in mode 1, with
// two different posting_budget values: one large enough to fit both
// volumes in a single Stage 1 batch, one small enough to force two
// batches. The orchestrator's per-batch mode1_results fold (Phase 4a / C5-a)
// must produce identical output across both batch shapes.
static void test_mode1_batched_equals_single() {
    std::fprintf(stderr, "-- test_mode1_batched_equals_single (Phase 4a)\n");
    const int k = 7;
    const std::string db_base = "mvtest";  // index built by earlier test

    auto discovered = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(discovered.size(), 2u);

    // Per-volume readers (opened once for sizing + ksx, orchestrator
    // reopens kix internally per batch).
    std::vector<KsxReader> ksxs(2);
    std::vector<VolumeMeta> volumes_cod(2);
    for (size_t vi = 0; vi < 2; vi++) {
        const auto& vf = discovered[vi];
        CHECK(ksxs[vi].open(vf.ksx_path));
        KixReader probe;
        CHECK(probe.open(vf.kix_path));
        volumes_cod[vi].files            = vf;
        volumes_cod[vi].kix_posting_size = probe.posting_file_size();
        volumes_cod[vi].kix_full_size    = probe.willneed_size_full();
        volumes_cod[vi].volume_index     = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences    = probe.num_sequences();
        probe.close();
    }

    SearchConfig config;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.nresult = 0;
    config.mode = 1;

    auto qdata = preprocess_query<uint16_t>(g_query_fj, k, nullptr, config);
    std::vector<QueryBundle<uint16_t>> bundles(1);
    bundles[0].query_id = &g_query_fj;  // string pointer (test-only lifetime ok)
    bundles[0].qdata_primary = &qdata;
    std::vector<uint8_t> skip_reason(1, 0);

    Logger logger(Logger::kError);
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    auto run_with_budget = [&](uint64_t budget) {
        RunSearchInputs<uint16_t> in;
        in.volumes_cod      = volumes_cod;  // copy: orchestrator reads only
        in.ksx_per_volume.resize(2);
        in.oid_filters.resize(2);
        for (size_t vi = 0; vi < 2; vi++) in.ksx_per_volume[vi] = &ksxs[vi];
        in.queries           = &bundles;
        in.query_skip_reason = &skip_reason;
        in.config            = config;
        in.both_mode         = false;
        in.k                 = k;
        in.nthread           = 2;
        in.posting_budget    = budget;
        in.logger            = &logger;
        in.max_num_seqs      = max_num_seqs;
        in.tier              = Stage1Tier::T32;
        return run_search<uint16_t>(in);
    };

    // Large budget: single batch over both volumes.
    uint64_t single_budget = 0;
    for (auto& v : volumes_cod) single_budget += v.kix_full_size;
    single_budget = std::max<uint64_t>(single_budget * 2, 1u << 20);
    auto big = run_with_budget(single_budget);

    // Small budget: each volume's kix_full_size on its own exceeds the
    // budget, so plan_stage1_batches emits one (tier1) batch per volume.
    uint64_t small_budget = 1;  // any positive < per-volume size
    auto small = run_with_budget(small_budget);

    auto key = [](const OrchestratorHit& h) {
        return std::tuple<size_t, uint16_t, uint16_t, uint32_t, uint32_t>(
            h.query_idx, h.volume_idx, h.volume_index,
            h.cr.seq_id, h.cr.stage1_score);
    };
    auto sorted_keys = [&](const std::vector<OrchestratorHit>& v) {
        std::vector<decltype(key(v.front()))> ks;
        ks.reserve(v.size());
        for (auto& h : v) ks.push_back(key(h));
        std::sort(ks.begin(), ks.end());
        return ks;
    };

    CHECK_EQ(big.size(), small.size());
    CHECK(sorted_keys(big) == sorted_keys(small));

    // Stress the parallel_scan fold path with many queries on the same
    // volume so each batch produces enough candidates for the prefix-sum
    // / parallel write to exercise multiple TBB workers. Run the same
    // configuration twice and require byte-identical output to confirm
    // the fold is deterministic.
    {
        std::vector<std::string> queries_storage(16, g_query_fj);
        std::vector<QueryKmerData<uint16_t>> qdatas;
        qdatas.reserve(queries_storage.size());
        for (auto& q : queries_storage) {
            qdatas.push_back(preprocess_query<uint16_t>(q, k, nullptr, config));
        }
        std::vector<QueryBundle<uint16_t>> many_bundles(queries_storage.size());
        for (size_t i = 0; i < queries_storage.size(); i++) {
            many_bundles[i].query_id      = &queries_storage[i];
            many_bundles[i].qdata_primary = &qdatas[i];
        }
        std::vector<uint8_t> skip(queries_storage.size(), 0);

        auto run = [&]() {
            RunSearchInputs<uint16_t> in;
            in.volumes_cod      = volumes_cod;
            in.ksx_per_volume.resize(2);
            in.oid_filters.resize(2);
            for (size_t vi = 0; vi < 2; vi++) in.ksx_per_volume[vi] = &ksxs[vi];
            in.queries           = &many_bundles;
            in.query_skip_reason = &skip;
            in.config            = config;
            in.both_mode         = false;
            in.k                 = k;
            in.nthread           = 4;
            in.posting_budget    = single_budget;
            in.logger            = &logger;
            in.max_num_seqs      = max_num_seqs;
            in.tier              = Stage1Tier::T32;
            return run_search<uint16_t>(in);
        };

        auto a = run();
        auto b = run();
        CHECK_EQ(a.size(), b.size());
        CHECK(sorted_keys(a) == sorted_keys(b));
    }

    for (auto& ksx : ksxs) ksx.close();
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
    test_mode1_batched_equals_single();
    test_mode2_batched_equals_single();
    test_result_merge_ordering();
    test_parallel_counting_pass();
    test_multivolume_k9();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
