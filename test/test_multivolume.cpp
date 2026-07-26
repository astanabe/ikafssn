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
#include "volume_search_helper.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

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
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

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
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

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
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

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
    sconfig.stage1.min_stage1_score = 1;
    sconfig.stage2.max_gap = 100;
    sconfig.stage2.min_nhit_diag = 1;
    sconfig.stage2.min_score = 2;

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
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

    auto seq_results = search_sequential(db_base, k, queries, config);
    auto par_results = search_parallel(db_base, k, queries, config, 2);

    CHECK_EQ(seq_results.size(), par_results.size());

    for (size_t i = 0; i < seq_results.size() && i < par_results.size(); i++) {
        CHECK(seq_results[i].qseqid == par_results[i].qseqid);
        CHECK_EQ(seq_results[i].chainscore, par_results[i].chainscore);
    }
}

// Drive run_search<>() in mode 2 with two different nthread values over a
// single-query, two-volume index.  With nthread large the query cannot
// fill the thread pool, so both volumes bundle into one group; with
// nthread == 1 the thread target is reached after one volume, so each
// volume runs in its own group.  The two grouping shapes must produce an
// identical sorted key set.
static void test_mode2_batched_equals_single() {
    std::fprintf(stderr, "-- test_mode2_batched_equals_single\n");
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
        volumes_cod[vi].volume_index     = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences    = kix_probe.num_sequences();
        kix_probe.close();
        kpx_probe.close();
    }

    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.mode = 2;

    auto qdata = preprocess_query<uint16_t>(g_query_fj, k, nullptr, config);
    std::vector<QueryBundle<uint16_t>> bundles(1);
    bundles[0].query_id = &g_query_fj;
    bundles[0].qdata_primary = &qdata;
    std::vector<uint8_t> skip_reason(1, 0);

    Logger logger(Logger::kError);
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    auto run_with_nthread = [&](int nthread) {
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
        in.nthread           = nthread;
        in.logger            = &logger;
        in.max_num_seqs      = max_num_seqs;
        in.width              = Stage1Width::T32;
        return run_search<uint16_t>(in);
    };

    // nthread large: the lone query can't fill the pool, so both volumes
    // bundle into one group.
    auto big = run_with_nthread(8);

    // nthread == 1: the thread target is met after one volume, so each
    // volume runs in its own group.
    auto small = run_with_nthread(1);

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

// bench/run_warm_e2e.sh attributes a run by parsing the key=value pairs of the
// "Timing run_search (s):" line, so pin the key set: renaming one key would
// silently drop a metric from every measurement instead of failing loudly.
static void test_timing_line_keys() {
    std::fprintf(stderr, "-- test_timing_line_keys\n");
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
        CHECK(kix_probe.open(vf.kix_path));
        volumes_cod[vi].files         = vf;
        volumes_cod[vi].volume_index  = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences = kix_probe.num_sequences();
        kix_probe.close();
    }

    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

    auto qdata = preprocess_query<uint16_t>(g_query_fj, k, nullptr, config);
    std::vector<QueryBundle<uint16_t>> bundles(1);
    bundles[0].query_id = &g_query_fj;
    bundles[0].qdata_primary = &qdata;
    std::vector<uint8_t> skip_reason(1, 0);

    Logger info_logger(Logger::kInfo);
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    // Run run_search with stderr redirected to a file and return the last
    // Timing line it wrote.
    auto captured_timing_line = [&](uint8_t mode) {
        RunSearchInputs<uint16_t> in;
        in.volumes_cod = volumes_cod;
        in.ksx_per_volume.resize(2);
        in.oid_filters.resize(2);
        for (size_t vi = 0; vi < 2; vi++) in.ksx_per_volume[vi] = &ksxs[vi];
        in.queries           = &bundles;
        in.query_skip_reason = &skip_reason;
        in.config            = config;
        in.config.mode       = mode;
        in.both_mode         = false;
        in.k                 = k;
        in.nthread           = 2;
        in.logger            = &info_logger;
        in.max_num_seqs      = max_num_seqs;
        in.width             = Stage1Width::T32;

        const std::string log_path =
            g_test_dir + "/timing_mode" + std::to_string(mode) + ".log";
        std::fflush(stderr);
        int saved_fd = dup(STDERR_FILENO);
        int log_fd = open(log_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        CHECK(saved_fd >= 0);
        CHECK(log_fd >= 0);
        CHECK(dup2(log_fd, STDERR_FILENO) >= 0);
        close(log_fd);

        run_search<uint16_t>(in);

        std::fflush(stderr);
        CHECK(dup2(saved_fd, STDERR_FILENO) >= 0);
        close(saved_fd);

        std::ifstream f(log_path);
        std::string line, timing;
        while (std::getline(f, line)) {
            if (line.find("Timing run_search (s):") != std::string::npos)
                timing = line;
        }
        return timing;
    };

    const std::string mode2 = captured_timing_line(2);
    CHECK(!mode2.empty());
    for (const char* key : {"s1_open=", "s1_compute=", "s1_fold=", "s1_intotal=",
                            "s2_open=", "s2a=", "s2b=", "s2_free=", "dedup=",
                            "parent_topn=", "s2_intotal=", "total="}) {
        CHECK(mode2.find(key) != std::string::npos);
    }

    // The mode 1 early return reports the Stage 1 keys only.
    const std::string mode1 = captured_timing_line(1);
    CHECK(!mode1.empty());
    for (const char* key : {"s1_open=", "s1_compute=", "s1_fold=", "s1_intotal=",
                            "total="}) {
        CHECK(mode1.find(key) != std::string::npos);
    }
    CHECK(mode1.find("s2a=") == std::string::npos);

    for (auto& ksx : ksxs) ksx.close();
}

// Drive run_search<>() in mode 2 with the Stage 2 in-total (L) cap and check
// the result against the tie-inclusive top-L reduction (by chainscore)
// recomputed from the unlimited baseline.
static void test_stage2_in_total_limit() {
    std::fprintf(stderr, "-- test_stage2_in_total_limit\n");
    const int k = 7;
    const std::string db_base = "mvtest";  // index built by an earlier test

    auto discovered = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(discovered.size(), 2u);

    std::vector<KsxReader> ksxs(2);
    std::vector<VolumeMeta> volumes_cod(2);
    for (size_t vi = 0; vi < 2; vi++) {
        const auto& vf = discovered[vi];
        CHECK(ksxs[vi].open(vf.ksx_path));
        KixReader probe;
        CHECK(probe.open(vf.kix_path));
        volumes_cod[vi].files         = vf;
        volumes_cod[vi].volume_index  = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences = probe.num_sequences();
        probe.close();
    }
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    Logger logger(Logger::kError);
    auto run = [&](const SearchConfig& cfg) {
        auto qdata = preprocess_query<uint16_t>(g_query_fj, k, nullptr, cfg);
        std::vector<QueryBundle<uint16_t>> bundles(1);
        bundles[0].query_id      = &g_query_fj;
        bundles[0].qdata_primary = &qdata;
        std::vector<uint8_t> skip(1, 0);
        RunSearchInputs<uint16_t> in;
        in.volumes_cod = volumes_cod;
        in.ksx_per_volume.resize(2);
        in.oid_filters.resize(2);
        for (size_t vi = 0; vi < 2; vi++) in.ksx_per_volume[vi] = &ksxs[vi];
        in.queries           = &bundles;
        in.query_skip_reason = &skip;
        in.config            = cfg;
        in.both_mode         = false;
        in.k                 = k;
        in.nthread           = 4;
        in.logger            = &logger;
        in.max_num_seqs      = max_num_seqs;
        in.width             = Stage1Width::T32;
        return run_search<uint16_t>(in);
    };

    // Unlimited baseline: keep every chain (no per-subject or in-total cap).
    SearchConfig base;
    base.stage1.min_stage1_score = 1;
    base.stage2.max_gap = 100;
    base.stage2.min_nhit_diag = 1;
    base.stage2.min_score = 2;
    base.stage2.max_nhit_per_subject = 0;
    base.mode = 2;
    auto baseline = run(base);
    CHECK(!baseline.empty());

    using Key = std::tuple<size_t, uint16_t, bool, uint32_t, uint32_t, uint32_t>;
    auto hkey = [](const OrchestratorHit& h) {
        return Key(h.query_idx, h.volume_idx, h.cr.is_reverse, h.cr.seq_id,
                   h.cr.s_start, h.cr.s_end);
    };
    auto keyset = [&](const std::vector<OrchestratorHit>& v) {
        std::set<Key> s;
        for (auto& h : v) s.insert(hkey(h));
        return s;
    };

    // Tie-inclusive top-L per query, by chainscore, recomputed from baseline.
    auto expected = [&](uint32_t L) {
        std::map<size_t, std::vector<std::pair<uint32_t, Key>>> by_q;
        for (auto& h : baseline)
            by_q[h.query_idx].push_back({h.cr.chainscore, hkey(h)});
        std::set<Key> keep;
        for (auto& kv : by_q) {
            auto& v = kv.second;
            if (v.size() <= L) {
                for (auto& p : v) keep.insert(p.second);
                continue;
            }
            std::vector<uint32_t> sc;
            sc.reserve(v.size());
            for (auto& p : v) sc.push_back(p.first);
            std::nth_element(sc.begin(), sc.begin() + (L - 1), sc.end(),
                             std::greater<uint32_t>());
            uint32_t thr = sc[L - 1];
            for (auto& p : v)
                if (p.first >= thr) keep.insert(p.second);
        }
        return keep;
    };

    for (uint32_t L : {1u, 2u, 3u}) {
        SearchConfig cfg = base;
        cfg.stage2.max_nhit_in_total = L;
        auto limited = run(cfg);
        CHECK(limited.size() <= baseline.size());
        CHECK(keyset(limited) == expected(L));
    }

    for (auto& ksx : ksxs) ksx.close();
}

// Drive run_search<>() directly with a multi-volume index, in mode 1, with
// two different nthread values over a single-query, two-volume index: a
// large nthread bundles both volumes into one Stage 1 group, while
// nthread == 1 reaches the thread target after one volume and runs each
// volume in its own group.  The orchestrator's per-group mode1_results
// fold must produce identical output across both grouping shapes.
static void test_mode1_batched_equals_single() {
    std::fprintf(stderr, "-- test_mode1_batched_equals_single\n");
    const int k = 7;
    const std::string db_base = "mvtest";  // index built by earlier test

    auto discovered = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(discovered.size(), 2u);

    // Per-volume readers (opened once for ksx + num_sequences; the
    // orchestrator reopens kix internally per group).
    std::vector<KsxReader> ksxs(2);
    std::vector<VolumeMeta> volumes_cod(2);
    for (size_t vi = 0; vi < 2; vi++) {
        const auto& vf = discovered[vi];
        CHECK(ksxs[vi].open(vf.ksx_path));
        KixReader probe;
        CHECK(probe.open(vf.kix_path));
        volumes_cod[vi].files            = vf;
        volumes_cod[vi].volume_index     = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences    = probe.num_sequences();
        probe.close();
    }

    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.mode = 1;

    auto qdata = preprocess_query<uint16_t>(g_query_fj, k, nullptr, config);
    std::vector<QueryBundle<uint16_t>> bundles(1);
    bundles[0].query_id = &g_query_fj;  // string pointer (test-only lifetime ok)
    bundles[0].qdata_primary = &qdata;
    std::vector<uint8_t> skip_reason(1, 0);

    Logger logger(Logger::kError);
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    auto run_with_nthread = [&](int nthread) {
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
        in.nthread           = nthread;
        in.logger            = &logger;
        in.max_num_seqs      = max_num_seqs;
        in.width              = Stage1Width::T32;
        return run_search<uint16_t>(in);
    };

    // nthread large: the lone query can't fill the pool, so both volumes
    // bundle into one Stage 1 group.
    auto big = run_with_nthread(8);

    // nthread == 1: the thread target is met after one volume, so each
    // volume runs in its own group.
    auto small = run_with_nthread(1);

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
            in.logger            = &logger;
            in.max_num_seqs      = max_num_seqs;
            in.width              = Stage1Width::T32;
            return run_search<uint16_t>(in);
        };

        auto a = run();
        auto b = run();
        CHECK_EQ(a.size(), b.size());
        CHECK(sorted_keys(a) == sorted_keys(b));
    }

    for (auto& ksx : ksxs) ksx.close();
}

// Drive run_search<>() in mode 1 with the Stage 1 candidate-count limits
// (per-subject N, per-volume M, in-total L) and check each result against the
// top-K+ties reduction recomputed from the unlimited baseline.
static void test_stage1_nhit_limits() {
    std::fprintf(stderr, "-- test_stage1_nhit_limits\n");
    const int k = 7;
    const std::string db_base = "mvtest";  // index built by an earlier test

    auto discovered = discover_volumes(g_test_dir + "/" + db_base, k);
    CHECK_EQ(discovered.size(), 2u);

    std::vector<KsxReader> ksxs(2);
    std::vector<VolumeMeta> volumes_cod(2);
    for (size_t vi = 0; vi < 2; vi++) {
        const auto& vf = discovered[vi];
        CHECK(ksxs[vi].open(vf.ksx_path));
        KixReader probe;
        CHECK(probe.open(vf.kix_path));
        volumes_cod[vi].files         = vf;
        volumes_cod[vi].volume_index  = static_cast<uint16_t>(vi);
        volumes_cod[vi].num_sequences = probe.num_sequences();
        probe.close();
    }
    uint32_t max_num_seqs = 0;
    for (auto& v : volumes_cod) max_num_seqs = std::max(max_num_seqs, v.num_sequences);

    // Several copies of the query so the per-group candidate lists are large
    // enough for the limits to bite.
    std::vector<std::string> qstore(8, g_query_fj);
    Logger logger(Logger::kError);

    auto run = [&](const SearchConfig& cfg) {
        std::vector<QueryKmerData<uint16_t>> qdatas;
        qdatas.reserve(qstore.size());
        for (auto& q : qstore)
            qdatas.push_back(preprocess_query<uint16_t>(q, k, nullptr, cfg));
        std::vector<QueryBundle<uint16_t>> bundles(qstore.size());
        for (size_t i = 0; i < qstore.size(); i++) {
            bundles[i].query_id      = &qstore[i];
            bundles[i].qdata_primary = &qdatas[i];
        }
        std::vector<uint8_t> skip(qstore.size(), 0);
        RunSearchInputs<uint16_t> in;
        in.volumes_cod = volumes_cod;
        in.ksx_per_volume.resize(2);
        in.oid_filters.resize(2);
        for (size_t vi = 0; vi < 2; vi++) in.ksx_per_volume[vi] = &ksxs[vi];
        in.queries           = &bundles;
        in.query_skip_reason = &skip;
        in.config            = cfg;
        in.both_mode         = false;
        in.k                 = k;
        in.nthread           = 4;
        in.logger            = &logger;
        in.max_num_seqs      = max_num_seqs;
        in.width             = Stage1Width::T32;
        return run_search<uint16_t>(in);
    };

    SearchConfig base;
    base.stage1.min_stage1_score = 1;
    base.mode = 1;
    auto baseline = run(base);
    CHECK(!baseline.empty());

    using Key = std::tuple<size_t, uint16_t, bool, uint32_t>;  // q, vol, strand, seqid
    auto hkey = [](const OrchestratorHit& h) {
        return Key(h.query_idx, h.volume_idx, h.cr.is_reverse, h.cr.seq_id);
    };
    auto keyset = [&](const std::vector<OrchestratorHit>& v) {
        std::set<Key> s;
        for (auto& h : v) s.insert(hkey(h));
        return s;
    };

    // Recompute the tie-inclusive top-K survivors from the baseline, grouping
    // hits by `groupid`.
    auto expected = [&](uint32_t K, auto groupid) {
        std::map<std::string, std::vector<std::pair<uint32_t, Key>>> groups;
        for (auto& h : baseline)
            groups[groupid(h)].push_back({h.cr.stage1_score, hkey(h)});
        std::set<Key> keep;
        for (auto& kv : groups) {
            auto& v = kv.second;
            if (v.size() <= K) {
                for (auto& p : v) keep.insert(p.second);
                continue;
            }
            std::vector<uint32_t> sc;
            sc.reserve(v.size());
            for (auto& p : v) sc.push_back(p.first);
            std::sort(sc.begin(), sc.end(), std::greater<uint32_t>());
            uint32_t kth = sc[K - 1];
            for (auto& p : v)
                if (p.first >= kth) keep.insert(p.second);
        }
        return keep;
    };

    auto gid_total = [](const OrchestratorHit& h) {
        return std::to_string(h.query_idx);
    };
    auto gid_vol = [](const OrchestratorHit& h) {
        return std::to_string(h.query_idx) + ":" +
               std::to_string(h.volume_idx) + ":" +
               std::to_string(h.cr.is_reverse);
    };
    auto gid_parent = [&](const OrchestratorHit& h) {
        uint32_t p = ksxs[h.volume_idx].parent_index(h.cr.seq_id);
        return gid_vol(h) + ":" + std::to_string(p);
    };

    // L (in-total): per query, across volumes / strands.
    for (uint32_t L : {1u, 2u}) {
        SearchConfig c = base;
        c.stage1.max_nhit_in_total = L;
        auto got = keyset(run(c));
        CHECK(got == expected(L, gid_total));
        CHECK(got.size() <= baseline.size());
    }
    // M (per-volume): per (query, volume, strand).
    for (uint32_t M : {1u, 2u}) {
        SearchConfig c = base;
        c.stage1.max_nhit_per_volume = M;
        CHECK(keyset(run(c)) == expected(M, gid_vol));
    }
    // N (per-subject), mode 3 = tie-inclusive: per (query, parent, volume, strand).
    {
        SearchConfig c = base;
        c.stage1.max_nhit_per_subject = 1;
        c.stage1.max_nhit_per_subject_mode = 3;
        CHECK(keyset(run(c)) == expected(1, gid_parent));
    }
    // Limits at/above the candidate count are no-ops.
    {
        SearchConfig c = base;
        c.stage1.max_nhit_in_total = 100000;
        c.stage1.max_nhit_per_volume = 100000;
        c.stage1.max_nhit_per_subject = 100000;
        CHECK(keyset(run(c)) == keyset(baseline));
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
    test_timing_line_keys();
    test_stage1_nhit_limits();
    test_stage2_in_total_limit();
    test_result_merge_ordering();
    test_parallel_counting_pass();
    test_multivolume_k9();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
