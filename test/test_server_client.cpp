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
#include "search/volume_searcher.hpp"
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
    SearchConfig config;
    config.stage2.min_score = 1;
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

// Helper: run two concurrent search clients against a single server fixture
// configured with a given max_concurrent_search.  Returns the peak number of
// concurrent leases observed inside the server's BudgetPool, plus both
// responses for invariant checks.
struct PoolConcurrencyOutcome {
    int peak_leases = 0;
    SearchResponse resp_a;
    SearchResponse resp_b;
    bool a_ok = false;
    bool b_ok = false;
};

static PoolConcurrencyOutcome run_two_concurrent_searches(
    const std::string& sock_path,
    const std::string& ix_prefix,
    const std::string& db,
    int k,
    const std::vector<FastaRecord>& queries,
    int max_concurrent_search)
{
    PoolConcurrencyOutcome out;

    ::unlink(sock_path.c_str());

    Server server;
    Logger logger(Logger::kError);
    ServerConfig server_config;
    server_config.max_concurrent_search = max_concurrent_search;
    // memory_limit=64MiB leaves a comfortable residual posting_budget after
    // khx/ksx WILLNEED on the small SSU test index; the concurrency
    // invariant is what we're asserting on, not any specific numeric budget.
    server_config.memory_limit = 64ull << 20;

    CHECK(server.load_database(ix_prefix, ix_prefix, server_config, logger));

    // Manually drive the budget + pool configure steps that Server::run()
    // would normally invoke (we don't want a full accept loop here — just
    // two concurrent request_processor invocations sharing the pool).
    server.apply_madvise_budget(server_config.memory_limit, logger);
    {
        uint64_t pb = server.posting_budget();
        uint64_t floor = (max_concurrent_search > 0)
            ? std::max<uint64_t>(1ull << 20,
                                  pb / static_cast<uint64_t>(max_concurrent_search))
            : 0;
        server.pool().configure(pb, floor);
    }

    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    // Two parallel server-side handlers.  The accept loop runs on its own
    // thread to feed each connection into a dedicated worker thread so both
    // request_processor invocations can overlap inside the BudgetPool.
    tbb::task_arena arena(2);
    std::atomic<int> accepted{0};
    std::vector<std::thread> workers;
    std::thread accept_thread([&] {
        for (int i = 0; i < 2; ++i) {
            int client_fd = accept_connection(listen_fd);
            if (client_fd < 0) break;
            accepted.fetch_add(1);
            workers.emplace_back([&, client_fd] {
                handle_connection(client_fd, server, server_config, arena, logger);
            });
        }
    });

    // Give the accept thread a moment to reach accept().
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    auto run_client = [&](SearchResponse& resp, bool& ok) {
        int fd = unix_connect(sock_path);
        if (fd < 0) { ok = false; return; }
        SearchRequest req;
        req.k = static_cast<uint8_t>(k);
        req.db = db;
        for (const auto& q : queries) {
            req.queries.push_back({q.id, q.sequence});
        }
        ok = socket_search(fd, req, resp);
        close_fd(fd);
    };

    std::thread ta([&] { run_client(out.resp_a, out.a_ok); });
    std::thread tb([&] { run_client(out.resp_b, out.b_ok); });
    ta.join();
    tb.join();

    accept_thread.join();
    for (auto& w : workers) w.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());

    out.peak_leases = server.pool().peak_leases();
    return out;
}

// Test: -max_concurrent_search 1 serialises two concurrent searches at the
// budget-bound stages; peak_leases must be 1 and both clients still succeed.
static void test_budget_pool_serialized() {
    std::fprintf(stderr, "-- test_budget_pool_serialized\n");

    int k = 7;
    std::string ix_prefix = g_test_dir + "/sc_index/test";
    std::string db = db_from_prefix(ix_prefix);

    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    std::string sock_path = g_test_dir + "/test_pool_serial.sock";
    auto out = run_two_concurrent_searches(
        sock_path, ix_prefix, db, k, queries, /*max_concurrent_search=*/1);

    CHECK(out.a_ok);
    CHECK(out.b_ok);
    CHECK_EQ(static_cast<uint64_t>(out.resp_a.status), uint64_t{0});
    CHECK_EQ(static_cast<uint64_t>(out.resp_b.status), uint64_t{0});
    CHECK_EQ(static_cast<uint64_t>(out.peak_leases), uint64_t{1});
    CHECK_EQ(out.resp_a.results.size(), out.resp_b.results.size());
}

// Test: -max_concurrent_search 2 lets two requests overlap; peak_leases
// reaches 2.
static void test_budget_pool_overlapping() {
    std::fprintf(stderr, "-- test_budget_pool_overlapping\n");

    int k = 7;
    std::string ix_prefix = g_test_dir + "/sc_index/test";
    std::string db = db_from_prefix(ix_prefix);

    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    std::string sock_path = g_test_dir + "/test_pool_overlap.sock";
    auto out = run_two_concurrent_searches(
        sock_path, ix_prefix, db, k, queries, /*max_concurrent_search=*/2);

    CHECK(out.a_ok);
    CHECK(out.b_ok);
    CHECK_EQ(static_cast<uint64_t>(out.resp_a.status), uint64_t{0});
    CHECK_EQ(static_cast<uint64_t>(out.resp_b.status), uint64_t{0});
    // peak may be 1 or 2 depending on scheduling, but must never exceed 2.
    CHECK(out.peak_leases >= 1 && out.peak_leases <= 2);
}

int main() {
    check_ssu_available();
    check_derived_data_ready();

    g_testdb_path = ssu_db_prefix();

    // Create temp directory
    g_test_dir = "/tmp/ikafssn_test_server_client_" + std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_server_client_search();
    test_health_check();
    test_seqidlist_filter_via_server();
    test_budget_pool_serialized();
    test_budget_pool_overlapping();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
