#include "ikafssnhttpd/http_controller.hpp"
#include "ikafssnhttpd/backend_manager.hpp"
#include "ikafssnhttpd/backend_client.hpp"
#include "ikafssnhttpd/job_store.hpp"
#include "ikafssnhttpd/job_worker.hpp"
#include "ikafssnhttpd/job_housekeeper.hpp"
#include "core/version.hpp"
#include "util/cli_parser.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "util/logger.hpp"
#include "util/socket_utils.hpp"

#include <drogon/HttpAppFramework.h>
#include <trantor/utils/Logger.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <thread>

#include <unistd.h>

using namespace ikafssn;

static void print_usage(const char* prog) {
    print_version_header("ikafssnhttpd");
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Backend connection (at least one required; order = priority):\n"
        "  -server_socket <path>      UNIX socket path to ikafssnserver\n"
        "  -server_tcp <host>:<port>  TCP address of ikafssnserver\n"
        "\n"
        "Job store (async REST):\n"
        "  -db <path>                  SQLite job DB path\n"
        "                              (default: /var/lib/ikafssnhttpd/jobs.db)\n"
        "  -retention_time <int>       Done/failed retention in seconds (default: 86400)\n"
        "  -max_nretry <int>           Per-job retry cap (default: 3)\n"
        "  -nthread_worker <int>       Backend search worker pool size (default: 4)\n"
        "\n"
        "Options:\n"
        "  -listen <host>:<port>       HTTP listen address (default: 0.0.0.0:8080)\n"
        "  -path_prefix <prefix>       API path prefix (e.g., /nt)\n"
        "  -nthread <int>              Drogon I/O threads (default: all cores)\n"
        "  -heartbeat_interval <int>   Heartbeat interval in seconds (default: 3600)\n"
        "  -exclusion_time <int>       Backend exclusion time in seconds (default: 3600)\n"
        "  -pid <path>                 PID file path\n"
        "  -v, --verbose               Verbose logging\n",
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (check_version(cli, "ikafssnhttpd")) return 0;

    if (cli.has("-h") || cli.has("-help")) {
        print_usage(argv[0]);
        return 0;
    }

    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);
    bool verbose = logger.verbose();

    auto manager = std::make_shared<BackendManager>();
    int backend_count = 0;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-server_socket") == 0 && i + 1 < argc) {
            std::string addr = argv[++i];
            manager->add_backend(BackendMode::kUnix, addr);
            logger.info("Backend %d: UNIX socket %s", backend_count, addr.c_str());
            backend_count++;
        } else if (std::strcmp(argv[i], "-server_tcp") == 0 && i + 1 < argc) {
            std::string addr = argv[++i];
            manager->add_backend(BackendMode::kTcp, addr);
            logger.info("Backend %d: TCP %s", backend_count, addr.c_str());
            backend_count++;
        }
    }

    if (backend_count == 0) {
        std::fprintf(stderr,
            "Error: at least one -server_socket or -server_tcp is required\n");
        print_usage(argv[0]);
        return 1;
    }

    int heartbeat_interval = cli.get_int("-heartbeat_interval", 3600);
    int exclusion_time = cli.get_int("-exclusion_time", 3600);
    manager->set_exclusion_time(exclusion_time);

    logger.info("Connecting to %d backend(s)...", backend_count);
    if (!manager->init(30, logger)) {
        std::fprintf(stderr,
            "Error: Failed to initialize backends after 30 seconds. "
            "Ensure ikafssnserver(s) are running.\n");
        return 1;
    }
    logger.info("All reachable backends initialized successfully.");

    manager->start_heartbeat(heartbeat_interval, logger);

    std::string db_path = cli.get_string("-db", "/var/lib/ikafssnhttpd/jobs.db");
    int retention_time = cli.get_int("-retention_time", 86400);
    int max_nretry     = cli.get_int("-max_nretry", 3);
    int nthread_worker = cli.get_int("-nthread_worker", 4);

    JobStore store;
    {
        std::string err;
        if (!store.open(db_path, err)) {
            std::fprintf(stderr, "Error: %s\n", err.c_str());
            manager->stop_heartbeat();
            return 1;
        }
    }
    logger.info("JobStore opened at %s (retention=%ds, max_nretry=%d, workers=%d)",
                db_path.c_str(), retention_time, max_nretry, nthread_worker);

    {
        std::string err;
        int64_t requeued = store.requeue_orphans(err);
        if (!err.empty()) {
            logger.error("requeue_orphans failed: %s", err.c_str());
        } else if (requeued > 0) {
            logger.info("Re-queued %lld orphan job(s) on startup",
                        static_cast<long long>(requeued));
        }
    }

    JobWorker worker(store, manager, logger, max_nretry);
    worker.start(nthread_worker);

    JobHousekeeper housekeeper(store, logger);
    housekeeper.start(retention_time);

    HttpController controller(manager, store, worker);
    std::string path_prefix = cli.get_string("-path_prefix");
    controller.register_routes(path_prefix);

    std::string listen_addr = cli.get_string("-listen", "0.0.0.0:8080");
    std::string host;
    uint16_t port;
    if (!parse_host_port(listen_addr, host, port)) {
        std::fprintf(stderr,
            "Error: invalid listen address '%s' (expected host:port)\n",
            listen_addr.c_str());
        housekeeper.stop();
        worker.stop();
        manager->stop_heartbeat();
        return 1;
    }

    int threads = resolve_threads(cli);

    drogon::app()
        .addListener(host, port)
        .setThreadNum(static_cast<size_t>(threads))
        .setLogLevel(verbose ? trantor::Logger::kDebug
                             : trantor::Logger::kWarn);

    std::string pid_file = cli.get_string("-pid");
    if (!pid_file.empty()) {
        FILE* f = std::fopen(pid_file.c_str(), "w");
        if (f) {
            std::fprintf(f, "%d\n", ::getpid());
            std::fclose(f);
        }
    }

    logger.info("Starting HTTP server on %s:%u (threads: %d, backends: %d)",
                host.c_str(), port, threads, backend_count);
    if (!path_prefix.empty()) {
        logger.info("API path prefix: %s", path_prefix.c_str());
    }
    logger.info("Heartbeat interval: %d seconds, exclusion time: %d seconds",
                heartbeat_interval, exclusion_time);

    drogon::app().run();

    housekeeper.stop();
    worker.stop();
    manager->stop_heartbeat();

    if (!pid_file.empty()) {
        std::remove(pid_file.c_str());
    }

    return 0;
}
