#include "test_util.hpp"
#include "ssu_test_fixture.hpp"

#include "search/stage3_alignment.hpp"
#include "volume_search_helper.hpp"
#include "search/query_preprocessor.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/output_coords.hpp"
#include "io/result_writer.hpp"
#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "core/spaced_seed.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;
static std::string g_query_seq;

static void test_stage3_pipeline() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_stage3_pipeline (build -> search -> align)\n");

    // Build index with k=7
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    // Extract query: 100bp from FJ876973.1
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    CHECK(fj_seq.size() >= 200);
    g_query_seq = fj_seq.substr(50, 100);
    db.close();

    Logger logger(Logger::kDebug);
    IndexBuilderConfig bconfig;
    bconfig.k = 7;

    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/s3test.00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    // Open index
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    // Stage 1+2 search
    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.mode = 2;

    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config, buf);

    CHECK(!result.hits.empty());

    // Convert to OutputHit
    std::vector<OutputHit> all_hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        all_hits.push_back(oh);
    }

    CHECK(!all_hits.empty());

    // Prepare query records for Stage 3
    std::vector<FastaRecord> queries;
    queries.push_back({"query1", g_query_seq});

    // Run Stage 3 with traceback
    Stage3Config s3config;
    s3config.traceback = true;
    s3config.gapopen = 10;
    s3config.gapext = 1;

    auto filtered = run_stage3(all_hits, queries, g_testdb_path, s3config,
                               false, 0.0, 0, logger);

    CHECK(!filtered.empty());

    // Verify alignment results
    for (const auto& h : filtered) {
        // alnscore should be positive for real matches
        CHECK(h.alnscore > 0);
        // CIGAR should be non-empty
        CHECK(!h.cigar.empty());
        // ppositive should be reasonable (at least some identity)
        CHECK(h.ppositive > 0.0);
        // npositive should be positive
        CHECK(h.npositive > 0);
        // Coordinates should be within sequence bounds
        CHECK(h.slen > 0);
        CHECK(h.sstart < h.slen);
        CHECK(h.send < h.slen);
    }

    // The exact match hit (FJ876973.1) should have high ppositive
    bool found_high_ppositive = false;
    for (const auto& h : filtered) {
        if (h.sseqid == ACC_FJ && h.sstrand == '+') {
            CHECK(h.ppositive > 90.0);
            found_high_ppositive = true;
        }
    }
    CHECK(found_high_ppositive);
}

static void test_stage3_score_only() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_stage3_score_only (traceback=0)\n");

    Logger logger(Logger::kError);

    // Build a quick index
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    std::string query = fj_seq.substr(50, 100);
    db.close();

    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/s3test2.00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

    auto qdata = preprocess_query<uint16_t>(query, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config, buf);
    CHECK(!result.hits.empty());

    std::vector<OutputHit> all_hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        all_hits.push_back(oh);
    }

    std::vector<FastaRecord> queries;
    queries.push_back({"query1", query});

    // Score-only mode (no traceback)
    Stage3Config s3config;
    s3config.traceback = false;

    auto filtered = run_stage3(all_hits, queries, g_testdb_path, s3config,
                               false, 0.0, 0, logger);
    CHECK(!filtered.empty());

    for (const auto& h : filtered) {
        CHECK(h.alnscore > 0);
        // CIGAR should be empty (no traceback)
        CHECK(h.cigar.empty());
        // ppositive should be 0 (not computed without traceback)
        CHECK(h.ppositive == 0.0);
    }
}

static void test_stage3_context() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_stage3_context\n");

    Logger logger(Logger::kError);

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    std::string query = fj_seq.substr(100, 50);
    db.close();

    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/s3test3.00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

    auto qdata = preprocess_query<uint16_t>(query, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config, buf);

    if (result.hits.empty()) {
        std::fprintf(stderr, "  (no hits, skipping context test)\n");
        return;
    }

    std::vector<OutputHit> all_hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        all_hits.push_back(oh);
    }

    std::vector<FastaRecord> queries;
    queries.push_back({"query1", query});

    // With integer context
    Stage3Config s3config;
    s3config.traceback = true;

    auto filtered = run_stage3(all_hits, queries, g_testdb_path, s3config,
                               false, 0.0, 50, logger);  // 50bp context
    CHECK(!filtered.empty());

    for (const auto& h : filtered) {
        CHECK(h.alnscore > 0);
        CHECK(!h.cigar.empty());
    }
}

// ---------------------------------------------------------------------------
// Stage 3 batching regression: a small posting_budget should partition hits
// into multiple batches without changing the filtered output.
// ---------------------------------------------------------------------------
static void build_hits_for_batching(std::vector<OutputHit>& all_hits,
                                    std::vector<FastaRecord>& queries,
                                    std::string& query_str,
                                    Logger& logger,
                                    const std::string& prefix_suffix)
{
    Stage1Buffer buf;

    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t fj_oid = find_oid_by_accession(db, ACC_FJ);
    CHECK(fj_oid != UINT32_MAX);
    std::string fj_seq = db.get_sequence(fj_oid);
    CHECK(fj_seq.size() >= 200);
    query_str = fj_seq.substr(50, 100);
    db.close();

    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    BlastDbReader db2;
    CHECK(db2.open(g_testdb_path));
    std::string prefix = g_test_dir + "/" + prefix_suffix + ".00.07mer";
    CHECK(build_index<uint16_t>(db2, bconfig, prefix, 0, 1, "test", logger));
    db2.close();

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.mode = 2;

    auto qdata = preprocess_query<uint16_t>(query_str, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config, buf);

    CHECK(!result.hits.empty());

    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = 0;
        oh.oid = cr.seq_id;
        oh.qlen = static_cast<uint32_t>(query_str.size());
        all_hits.push_back(oh);
    }

    queries.push_back({"query1", query_str});
}

namespace {
struct BatchingHitKey {
    std::string qseqid;
    std::string sseqid;
    char        sstrand;
    uint32_t    sstart;
    uint32_t    send;
    int32_t     alnscore;

    bool operator<(const BatchingHitKey& o) const {
        if (qseqid != o.qseqid) return qseqid < o.qseqid;
        if (sseqid != o.sseqid) return sseqid < o.sseqid;
        if (sstrand != o.sstrand) return sstrand < o.sstrand;
        if (sstart != o.sstart) return sstart < o.sstart;
        if (send != o.send) return send < o.send;
        return alnscore < o.alnscore;
    }
    bool operator==(const BatchingHitKey& o) const {
        return qseqid == o.qseqid && sseqid == o.sseqid && sstrand == o.sstrand
            && sstart == o.sstart && send == o.send && alnscore == o.alnscore;
    }
};
} // namespace

static std::vector<BatchingHitKey> to_keys(const std::vector<OutputHit>& hits) {
    std::vector<BatchingHitKey> keys;
    keys.reserve(hits.size());
    for (const auto& h : hits) {
        keys.push_back({h.qseqid, h.sseqid, h.sstrand, h.sstart, h.send, h.alnscore});
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

static void test_stage3_batching_equivalence() {
    std::fprintf(stderr, "-- test_stage3_batching_equivalence (small budget == unbatched)\n");
    Logger logger(Logger::kError);

    std::vector<OutputHit> hits_a;
    std::vector<FastaRecord> queries_a;
    std::string query_a;
    build_hits_for_batching(hits_a, queries_a, query_a, logger, "s3batch_a");

    std::vector<OutputHit> hits_b = hits_a;  // independent copy for second run
    std::vector<FastaRecord> queries_b = queries_a;

    Stage3Config cfg;
    cfg.traceback = true;
    cfg.gapopen = 10;
    cfg.gapext = 1;

    // Unbatched run.
    cfg.posting_budget = 0;
    auto filtered_a = run_stage3(hits_a, queries_a, g_testdb_path, cfg,
                                 false, 0.0, 0, logger);

    // Force >=3 batches: pick a budget that holds roughly one third of the
    // groups.  Since all hits share (qseqid, sseqid) only when sseqid
    // matches, with SSU SSU the result mostly has many distinct sseqids =>
    // many groups => budget = (sum / 3) splits them.
    uint64_t total_cost = 0;
    for (const auto& h : hits_b) {
        uint64_t subseq = static_cast<uint64_t>(h.send - h.sstart + 1);
        uint64_t tb = 2ull * std::max<uint64_t>(h.qlen, subseq);
        total_cost += subseq + tb + 256;
    }
    cfg.posting_budget = total_cost / 4;  // small enough to force splitting
    if (cfg.posting_budget == 0) cfg.posting_budget = 1024;

    auto filtered_b = run_stage3(hits_b, queries_b, g_testdb_path, cfg,
                                 false, 0.0, 0, logger);

    CHECK(filtered_a.size() == filtered_b.size());
    auto ka = to_keys(filtered_a);
    auto kb = to_keys(filtered_b);
    CHECK(ka == kb);
    if (ka != kb) {
        std::fprintf(stderr, "  filtered_a: %zu hits, filtered_b: %zu hits\n",
                     filtered_a.size(), filtered_b.size());
    }
}

static void test_stage3_oversize_group_solo() {
    std::fprintf(stderr, "-- test_stage3_oversize_group_solo (single group exceeds budget)\n");
    Logger logger(Logger::kError);

    std::vector<OutputHit> hits;
    std::vector<FastaRecord> queries;
    std::string query;
    build_hits_for_batching(hits, queries, query, logger, "s3oversize");

    // Find a (qseqid, sseqid, sstrand) with at least one hit and copy that
    // group repeatedly so the group cost dwarfs the budget.  The
    // overlap-resolution loop expects sstart-sorted hits within a group,
    // which the planner normalizes anyway.
    if (hits.empty()) {
        std::fprintf(stderr, "  (no hits, skipping oversize test)\n");
        return;
    }
    // All synthesized clones share qseqid/sseqid/sstrand of hits[0].  We
    // make the group dominant by adding many copies with identical
    // coordinates — they will collide in overlap resolution but that is
    // exactly what should also happen in the unbatched path, and what we
    // care about here is that the planner emits a solo batch.
    OutputHit base = hits[0];
    std::vector<OutputHit> baseline = hits;  // for equivalence comparison
    for (int i = 0; i < 8; i++) hits.push_back(base);
    baseline = hits;  // baseline now matches the inflated input

    Stage3Config cfg;
    cfg.traceback = false;
    cfg.gapopen = 10;
    cfg.gapext = 1;

    cfg.posting_budget = 0;
    auto unbatched = run_stage3(hits, queries, g_testdb_path, cfg,
                                false, 0.0, 0, logger);

    cfg.posting_budget = 1;  // any positive group cost > 1 forces a solo batch
    auto batched = run_stage3(baseline, queries, g_testdb_path, cfg,
                              false, 0.0, 0, logger);

    CHECK(unbatched.size() == batched.size());
    CHECK(to_keys(unbatched) == to_keys(batched));
}

// Reverse-strand chains reach Stage 3 only once the query k-mer positions are
// query-strand-relative, so this covers the alignment of a reverse-complement
// query against the forward subject end to end.
static void test_stage3_reverse_strand() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_stage3_reverse_strand\n");

    std::string prefix = g_test_dir + "/s3test.00.07mer";
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;
    config.mode = 2;

    Logger logger(Logger::kError);
    Stage3Config s3config;
    s3config.traceback = true;
    s3config.gapopen = 10;
    s3config.gapext = 1;

    // Align one query and return its ACC_FJ hit.  `found` distinguishes a
    // missing hit from a default-constructed one.
    auto align_one = [&](const std::string& query, char want_strand,
                         bool& found) -> OutputHit {
        found = false;
        auto qdata = preprocess_query<uint16_t>(query, 7, nullptr, config);
        auto result = search_volume<uint16_t>(
            "q", qdata, 7, kix, kpx, ksx, filter, config, buf);

        std::vector<OutputHit> hits;
        for (const auto& cr : result.hits) {
            if (std::string(ksx.accession(cr.seq_id)) != ACC_FJ) continue;
            if ((cr.is_reverse ? '-' : '+') != want_strand) continue;
            OutputHit oh;
            oh.qseqid = "q";
            oh.sseqid = ACC_FJ;
            oh.sstrand = want_strand;
            oh.qstart = cr.q_start;
            oh.qend = cr.q_end;
            oh.sstart = cr.s_start;
            oh.send = cr.s_end;
            oh.chainscore = cr.chainscore;
            oh.coverscore = cr.stage1_score;
            oh.volume = 0;
            oh.oid = cr.seq_id;
            oh.qlen = static_cast<uint32_t>(query.size());
            hits.push_back(std::move(oh));
        }
        if (hits.empty()) return {};

        std::vector<FastaRecord> queries{{"q", query}};
        auto filtered = run_stage3(hits, queries, g_testdb_path, s3config,
                                   false, 0.0, 0, logger);
        if (filtered.empty()) return {};
        found = true;
        return filtered.front();
    };

    bool found_fwd = false, found_rev = false;
    OutputHit fwd = align_one(g_query_seq, '+', found_fwd);
    OutputHit rev = align_one(reverse_complement_string(g_query_seq), '-',
                              found_rev);

    CHECK(found_fwd);
    CHECK(found_rev);
    if (!found_fwd || !found_rev) return;

    // The reverse-complement query is aligned against the same forward
    // subject, so the alignment itself is identical.
    CHECK(rev.ppositive > 90.0);
    CHECK(rev.npositive > 0);
    CHECK(!rev.cigar.empty());
    CHECK(rev.cigar == fwd.cigar);
    CHECK(rev.qseq == fwd.qseq);
    CHECK(rev.sseq == fwd.sseq);
    CHECK_EQ(rev.sstart, fwd.sstart);
    CHECK_EQ(rev.send, fwd.send);
    CHECK(rev.sstart < rev.slen);
    CHECK(rev.send <= rev.slen);

    // On output the reverse hit keeps ascending query coordinates and its
    // subject coordinates descend.
    OutputCoords c = to_output_coords(rev.qstart, rev.qend, rev.sstart,
                                      rev.send, rev.qlen, true);
    CHECK(c.qstart <= c.qend);
    CHECK(c.sstart > c.send);
}

int main() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_test_dir = "/tmp/ikafssn_test_stage3";
    std::filesystem::create_directories(g_test_dir);

    test_stage3_pipeline();
    test_stage3_score_only();
    test_stage3_context();
    test_stage3_batching_equivalence();
    test_stage3_oversize_group_solo();
    test_stage3_reverse_strand();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
