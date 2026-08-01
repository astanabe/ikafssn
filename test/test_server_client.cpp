#include "test_util.hpp"
#include "ssu_test_fixture.hpp"

#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "io/volume_discovery.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"
#include "volume_search_helper.hpp"
#include "ikafssnserver/server.hpp"
#include "ikafssnserver/connection_handler.hpp"
#include "ikafssnserver/request_processor.hpp"
#include "ikafssnclient/socket_client.hpp"
#include "protocol/frame.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"
#include "core/config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>

#include <tbb/task_arena.h>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;

// Send a health-check request to the server and read its response (test-only;
// production health checks go through ikafssnhttpd's BackendClient).
static bool socket_health_check(int fd, HealthResponse& resp) {
    HealthRequest hreq;
    auto payload = serialize(hreq);
    if (!write_frame(fd, MsgType::kHealthRequest, payload)) {
        return false;
    }
    FrameHeader hdr;
    std::vector<uint8_t> resp_payload;
    if (!read_frame(fd, hdr, resp_payload)) {
        return false;
    }
    if (static_cast<MsgType>(hdr.msg_type) != MsgType::kHealthResponse) {
        return false;
    }
    return deserialize(resp_payload, resp);
}

// Build test index and return the index prefix path
static std::string build_test_index(int k) {
    std::string ix_dir = g_test_dir + "/sc_index";
    std::filesystem::create_directories(ix_dir);

    BlastDbReader db;
    if (!db.open(g_testdb_path)) {
        std::fprintf(stderr, "Cannot open test BLAST DB at %s\n", g_testdb_path.c_str());
        return {};
    }

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = k;

    auto stem_for = [&](const std::string& base) {
        return index_file_stem(ix_dir, base, bconfig.k,
                               bconfig.t, bconfig.template_type,
                               bconfig.min_seq_length,
                               bconfig.min_length_split,
                               bconfig.overlap_length,
                               /*max_freq_build=*/1,
                               /*max_degen_expand=*/0);
    };
    std::string prefix = stem_for("test.00");

    bool ok;
    if (k < K_TYPE_THRESHOLD) {
        ok = build_index<uint16_t>(db, bconfig, prefix, 0, 1, "test", logger);
    } else {
        ok = build_index<uint32_t>(db, bconfig, prefix, 0, 1, "test", logger);
    }

    if (!ok) {
        std::fprintf(stderr, "Failed to build index for k=%d\n", k);
        return {};
    }

    // Write .kvx manifest (normally done by ikafssnindex main)
    std::string kvx_path = stem_for("test") + ".kvx";
    FILE* fp = std::fopen(kvx_path.c_str(), "w");
    if (fp) {
        std::fprintf(fp, "#\n# ikafssn index volume manifest\n#\n");
        std::fprintf(fp, "TITLE test\n");
        std::fprintf(fp, "DBLIST \"test.00\"\n");
        std::fclose(fp);
    }

    return ix_dir + "/test";
}

// Derive DB name from index prefix (same logic as server)
static std::string db_from_prefix(const std::string& ix_prefix) {
    return parse_index_prefix(ix_prefix).db;
}

// Test: direct local search produces results, and server produces same results
static void test_server_client_search() {
    std::fprintf(stderr, "-- test_server_client_search\n");

    int k = 7;
    std::string ix_prefix = build_test_index(k);
    CHECK(!ix_prefix.empty());

    std::string db = db_from_prefix(ix_prefix);

    // Read query from derived test data
    std::string query_fasta = queries_path();
    auto queries = read_fasta(query_fasta);
    CHECK(!queries.empty());

    // --- Direct local search for reference results ---
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    auto vols = discover_volumes(ix_prefix, k);
    CHECK(!vols.empty());
    CHECK(kix.open(vols[0].kix_path));
    CHECK(kpx.open(vols[0].kpx_path));
    CHECK(ksx.open(vols[0].ksx_path));

    // Build search config matching server defaults.
    // The local reference uses search_volume(), which chains per
    // (fragment, strand) and does not run the orchestrator's parent-level
    // top-N selection.  Mode 2 (strands separate, no score ties) makes the
    // parent selection a no-op on this non-fragmented index, so the two
    // paths stay directly comparable.
    SearchConfig config;
    config.stage2.min_score = 1;
    config.stage2.max_nhit_per_subject_mode = 2;
    OidFilter no_filter;

    std::vector<SearchResult> local_results;
    for (const auto& q : queries) {
    Stage1Buffer buf;
        auto qdata = preprocess_query<uint16_t>(q.sequence, k, nullptr, config);
        auto sr = search_volume<uint16_t>(q.id, qdata, k,
                                          kix, kpx, ksx, no_filter, config, buf);
        local_results.push_back(sr);
    }

    // --- Set up server in a thread ---
    std::string sock_path = g_test_dir + "/test_server.sock";
    // Remove stale socket
    ::unlink(sock_path.c_str());

    // Load server index
    Server server;
    Logger logger(Logger::kError);
    ServerConfig server_config;
    server_config.search_config = config;
    CHECK(server.load_database(ix_prefix, ix_prefix, server_config, logger));
    CHECK(server.default_k() == k);

    // Start a listening socket
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    // Accept one connection in a background thread and process it
    tbb::task_arena arena(1);
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server, server_config, arena, logger);
        }
    });

    // Give server thread time to reach accept()
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    // --- Client connects and searches ---
    int fd = unix_connect(sock_path);
    CHECK(fd >= 0);

    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    req.db = db;
    req.min_seq_length = 64;  // resolved variant identity (default builder min_seq_length)
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    CHECK(socket_search(fd, req, resp));
    close_fd(fd);

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());

    // --- Verify results match ---
    CHECK(resp.status == 0);
    CHECK_EQ(resp.k, static_cast<uint8_t>(k));
    CHECK_EQ(resp.results.size(), local_results.size());

    for (size_t qi = 0; qi < resp.results.size() && qi < local_results.size(); qi++) {
        const auto& server_qr = resp.results[qi];
        const auto& local_sr = local_results[qi];
        CHECK(server_qr.qseqid == local_sr.query_id);
        CHECK_EQ(server_qr.hits.size(), local_sr.hits.size());

        // Build comparable lists sorted by (sseqid, sstart, qstart) to avoid ordering issues
        struct HitKey {
            std::string sseqid;
            uint32_t qstart, qend, sstart, send;
            uint16_t chainscore;
            bool is_reverse;
            bool operator<(const HitKey& o) const {
                if (sseqid != o.sseqid) return sseqid < o.sseqid;
                if (sstart != o.sstart) return sstart < o.sstart;
                return qstart < o.qstart;
            }
        };

        std::vector<HitKey> server_sorted, local_sorted;
        for (const auto& sh : server_qr.hits) {
            server_sorted.push_back({sh.sseqid, sh.qstart, sh.qend,
                                     sh.sstart, sh.send, sh.chainscore, sh.sstrand == 1});
        }
        for (const auto& lh : local_sr.hits) {
            local_sorted.push_back({std::string(ksx.accession(lh.seq_id)),
                                    lh.q_start, lh.q_end, lh.s_start, lh.s_end,
                                    static_cast<uint16_t>(lh.chainscore), lh.is_reverse});
        }
        std::sort(server_sorted.begin(), server_sorted.end());
        std::sort(local_sorted.begin(), local_sorted.end());

        for (size_t hi = 0; hi < server_sorted.size() && hi < local_sorted.size(); hi++) {
            CHECK(server_sorted[hi].sseqid == local_sorted[hi].sseqid);
            CHECK_EQ(server_sorted[hi].qstart, local_sorted[hi].qstart);
            CHECK_EQ(server_sorted[hi].qend, local_sorted[hi].qend);
            CHECK_EQ(server_sorted[hi].sstart, local_sorted[hi].sstart);
            CHECK_EQ(server_sorted[hi].send, local_sorted[hi].send);
            CHECK_EQ(server_sorted[hi].chainscore, local_sorted[hi].chainscore);
            CHECK(server_sorted[hi].is_reverse == local_sorted[hi].is_reverse);
        }
    }
}

// ---------------------------------------------------------------------------
// Request processing without a socket: process_search_request() is where the
// orchestrator hits are turned into response hits and routed back to their
// QueryResult, so these drive it directly.
// ---------------------------------------------------------------------------

// Two query sequences of different lengths, so a hit's qlen identifies the
// query it was aligned against.
struct TwoQueries {
    std::string a;   // 200 bp of ACC_FJ
    std::string b;   // 150 bp of ACC_GQ
};

static TwoQueries make_two_queries() {
    TwoQueries q;
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t fj = find_oid_by_accession(db, ACC_FJ);
    uint32_t gq = find_oid_by_accession(db, ACC_GQ);
    CHECK(fj != UINT32_MAX);
    CHECK(gq != UINT32_MAX);
    if (fj == UINT32_MAX || gq == UINT32_MAX) return q;
    std::string fj_seq = db.get_sequence(fj);
    std::string gq_seq = db.get_sequence(gq);
    CHECK(fj_seq.size() >= 400);
    CHECK(gq_seq.size() >= 400);
    db.close();
    q.a = fj_seq.substr(100, 200);
    q.b = gq_seq.substr(100, 150);
    return q;
}

// Load the test index with a BLAST DB attached, so mode 3 is available.
static bool load_server_with_db(Server& server, ServerConfig& cfg,
                                std::string& dbname) {
    std::string ix_prefix = build_test_index(7);
    CHECK(!ix_prefix.empty());
    if (ix_prefix.empty()) return false;
    cfg.search_config.stage2.min_score = 1;
    cfg.stage3_config.traceback = true;
    cfg.context_is_ratio = false;
    cfg.context_abs = 0;
    Logger logger(Logger::kError);
    if (!server.load_database(ix_prefix, g_testdb_path, cfg, logger)) return false;
    dbname = db_from_prefix(ix_prefix);
    return true;
}

static SearchRequest make_request(const std::string& dbname, uint8_t mode) {
    SearchRequest req;
    req.k = 7;
    req.db = dbname;
    req.min_seq_length = 64;
    req.mode = mode;
    req.stage3_traceback = 1;
    return req;
}

// Every hit must reach the QueryResult of the query it was searched with; the
// per-hit qlen is the query's own length, so it names that query.
static void check_hits_routed(const SearchResponse& resp) {
    for (const auto& qr : resp.results) {
        for (const auto& h : qr.hits) CHECK_EQ(h.qlen, qr.qlen);
    }
}

static void test_server_mode3_roundtrip() {
    std::fprintf(stderr, "-- test_server_mode3_roundtrip\n");

    Server server;
    ServerConfig cfg;
    std::string dbname;
    if (!load_server_with_db(server, cfg, dbname)) { CHECK(false); return; }
    const DatabaseEntry* entry = server.find_database(dbname);
    CHECK(entry != nullptr);
    if (!entry) return;

    TwoQueries q = make_two_queries();
    SearchRequest req = make_request(dbname, 3);
    req.queries.push_back({"qa", q.a});
    req.queries.push_back({"qb", q.b});

    tbb::task_arena arena(1);
    SearchResponse resp = process_search_request(req, *entry, server, arena);

    CHECK_EQ(resp.status, uint8_t{0});
    CHECK_EQ(resp.results.size(), 2u);
    if (resp.results.size() != 2) return;
    CHECK(resp.results[0].qseqid == "qa");
    CHECK(resp.results[1].qseqid == "qb");
    check_hits_routed(resp);

    for (const auto& qr : resp.results) {
        CHECK(!qr.hits.empty());
        for (const auto& h : qr.hits) {
            // Stage 3 ran: every hit carries an alignment.
            CHECK(h.alnscore != 0);
            CHECK(!h.cigar.empty());
        }
    }
}

// Two queries may share a qseqid; each still gets its own results, so a
// name-keyed route back would drop one of them entirely.
static void test_server_mode3_duplicate_qseqid() {
    std::fprintf(stderr, "-- test_server_mode3_duplicate_qseqid\n");

    Server server;
    ServerConfig cfg;
    std::string dbname;
    if (!load_server_with_db(server, cfg, dbname)) { CHECK(false); return; }
    const DatabaseEntry* entry = server.find_database(dbname);
    CHECK(entry != nullptr);
    if (!entry) return;

    TwoQueries q = make_two_queries();
    SearchRequest req = make_request(dbname, 3);
    req.queries.push_back({"dup", q.a});
    req.queries.push_back({"dup", q.b});

    tbb::task_arena arena(1);
    SearchResponse resp = process_search_request(req, *entry, server, arena);

    CHECK_EQ(resp.status, uint8_t{0});
    CHECK_EQ(resp.results.size(), 2u);
    if (resp.results.size() != 2) return;
    CHECK(resp.results[0].qseqid == "dup");
    CHECK(resp.results[1].qseqid == "dup");
    CHECK(!resp.results[0].hits.empty());
    CHECK(!resp.results[1].hits.empty());
    check_hits_routed(resp);
}

// The mode 2 hits of one query arrive in the order the Stage 2 dedup left
// them: one contiguous run per subject, ascending (strand, sstart) inside it.
static void test_server_mode2_hit_order() {
    std::fprintf(stderr, "-- test_server_mode2_hit_order\n");

    Server server;
    ServerConfig cfg;
    std::string dbname;
    if (!load_server_with_db(server, cfg, dbname)) { CHECK(false); return; }
    const DatabaseEntry* entry = server.find_database(dbname);
    CHECK(entry != nullptr);
    if (!entry) return;

    TwoQueries q = make_two_queries();
    SearchRequest req = make_request(dbname, 2);
    req.queries.push_back({"qa", q.a});
    req.queries.push_back({"qb", q.b});

    tbb::task_arena arena(1);
    SearchResponse resp = process_search_request(req, *entry, server, arena);

    CHECK_EQ(resp.status, uint8_t{0});
    CHECK_EQ(resp.results.size(), 2u);
    if (resp.results.size() != 2) return;
    check_hits_routed(resp);

    for (const auto& qr : resp.results) {
        CHECK(!qr.hits.empty());
        std::vector<std::string> seen;
        for (size_t i = 0; i < qr.hits.size(); ++i) {
            const auto& h = qr.hits[i];
            if (i == 0 || h.sseqid != qr.hits[i - 1].sseqid) {
                // A subject must not come back after another one intervened.
                CHECK(std::find(seen.begin(), seen.end(), h.sseqid) == seen.end());
                seen.push_back(h.sseqid);
                continue;
            }
            const auto& prev = qr.hits[i - 1];
            CHECK(prev.sstrand <= h.sstrand);
            if (prev.sstrand == h.sstrand) CHECK(prev.sstart <= h.sstart);
        }
    }
}

// Test: health check
static void test_health_check() {
    std::fprintf(stderr, "-- test_health_check\n");

    std::string sock_path = g_test_dir + "/test_health.sock";
    ::unlink(sock_path.c_str());

    int k = 7;
    std::string ix_prefix = g_test_dir + "/sc_index/test";

    Server server;
    Logger logger(Logger::kError);
    ServerConfig server_config;
    CHECK(server.load_database(ix_prefix, ix_prefix, server_config, logger));

    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    tbb::task_arena arena(1);
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server, server_config, arena, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    int fd = unix_connect(sock_path);
    CHECK(fd >= 0);

    HealthResponse hresp;
    CHECK(socket_health_check(fd, hresp));
    CHECK_EQ(hresp.status, 0);

    close_fd(fd);
    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());
}

// Test: seqidlist filtering through server
static void test_seqidlist_filter_via_server() {
    std::fprintf(stderr, "-- test_seqidlist_filter_via_server\n");

    int k = 7;
    std::string ix_prefix = g_test_dir + "/sc_index/test";
    std::string db = db_from_prefix(ix_prefix);
    std::string sock_path = g_test_dir + "/test_seqidlist.sock";
    ::unlink(sock_path.c_str());

    // Read queries
    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    // Get an accession from the ksx to use as seqidlist filter
    KsxReader ksx;
    auto vols = discover_volumes(ix_prefix, k);
    CHECK(!vols.empty());
    CHECK(ksx.open(vols[0].ksx_path));

    std::string target_acc;
    if (ksx.num_sequences() > 0) {
        target_acc = std::string(ksx.accession(0));
    }
    CHECK(!target_acc.empty());

    Server server;
    Logger logger(Logger::kError);
    ServerConfig server_config;
    CHECK(server.load_database(ix_prefix, ix_prefix, server_config, logger));

    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    tbb::task_arena arena(1);
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server, server_config, arena, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    int fd = unix_connect(sock_path);
    CHECK(fd >= 0);

    // Build request with seqidlist (include mode)
    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    req.db = db;
    req.min_seq_length = 64;  // resolved variant identity (default builder min_seq_length)
    req.seqidlist_mode = SeqidlistMode::kInclude;
    req.seqids = {target_acc};
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    CHECK(socket_search(fd, req, resp));
    close_fd(fd);

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());

    CHECK(resp.status == 0);

    // All hits should have the target accession
    for (const auto& qr : resp.results) {
        for (const auto& hit : qr.hits) {
            CHECK(hit.sseqid == target_acc);
        }
    }
}

int main() {
    check_ssu_available();
    check_derived_data_ready();

    g_testdb_path = ssu_db_prefix();

    // Create temp directory
    g_test_dir = "/tmp/ikafssn_test_server_client_" + std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_server_client_search();
    test_server_mode3_roundtrip();
    test_server_mode3_duplicate_qseqid();
    test_server_mode2_hit_order();
    test_health_check();
    test_seqidlist_filter_via_server();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
