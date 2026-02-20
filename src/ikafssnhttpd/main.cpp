#include "ikafssnhttpd/http_controller.hpp"
#include "ikafssnhttpd/backend_client.hpp"
#include "util/cli_parser.hpp"
#include "util/logger.hpp"
#include "util/socket_utils.hpp"

#include <drogon/HttpAppFramework.h>
#include <trantor/utils/Logger.h>

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

#include <unistd.h>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Backend connection (one required):\n"
        "  -server_socket <path>      UNIX socket path to ikafssnserver\n"
        "  -server_tcp <host>:<port>  TCP address of ikafssnserver\n"
        "\n"
        "Options:\n"
        "  -listen <host>:<port>  HTTP listen address (default: 0.0.0.0:8080)\n"
        "  -path_prefix <prefix>  API path prefix (e.g., /nt)\n"
        "  -threads <int>         Drogon I/O threads (default: 4)\n"
        "  -pid <path>            PID file path\n"
        "  -v, --verbose          Verbose logging\n",
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (cli.has("-h") || cli.has("--help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (!cli.has("-server_socket") && !cli.has("-server_tcp")) {
        std::fprintf(stderr,
            "Error: one of -server_socket or -server_tcp is required\n");
        print_usage(argv[0]);
        return 1;
    }

    bool verbose = cli.has("-v") || cli.has("--verbose");
    Logger logger(verbose ? Logger::kDebug : Logger::kInfo);

    // Create backend client
    BackendMode mode;
    std::string backend_addr;
    if (cli.has("-server_socket")) {
        mode = BackendMode::kUnix;
        backend_addr = cli.get_string("-server_socket");
        logger.info("Backend: UNIX socket %s", backend_addr.c_str());
    } else {
        mode = BackendMode::kTcp;
        backend_addr = cli.get_string("-server_tcp");
        logger.info("Backend: TCP %s", backend_addr.c_str());
    }

    auto backend = std::make_shared<BackendClient>(mode, backend_addr);

    // Create HTTP controller and register routes
    HttpController controller(backend);
    std::string path_prefix = cli.get_string("-path_prefix");
    controller.register_routes(path_prefix);

    // Parse listen address
    std::string listen_addr = cli.get_string("-listen", "0.0.0.0:8080");
    std::string host;
    uint16_t port;
    if (!parse_host_port(listen_addr, host, port)) {
        std::fprintf(stderr,
            "Error: invalid listen address '%s' (expected host:port)\n",
            listen_addr.c_str());
        return 1;
    }

    int threads = cli.get_int("-threads", 4);

    // Configure Drogon
    drogon::app()
        .addListener(host, port)
        .setThreadNum(static_cast<size_t>(threads))
        .setLogLevel(verbose ? trantor::Logger::kDebug
                             : trantor::Logger::kWarn);

    // PID file
    std::string pid_file = cli.get_string("-pid");
    if (!pid_file.empty()) {
        FILE* f = std::fopen(pid_file.c_str(), "w");
        if (f) {
            std::fprintf(f, "%d\n", ::getpid());
            std::fclose(f);
        }
    }

    logger.info("Starting HTTP server on %s:%u (threads: %d)",
                host.c_str(), port, threads);
    if (!path_prefix.empty()) {
        logger.info("API path prefix: %s", path_prefix.c_str());
    }

    // Run Drogon (blocks until shutdown via SIGTERM/SIGINT)
    drogon::app().run();

    // Cleanup PID file
    if (!pid_file.empty()) {
        std::remove(pid_file.c_str());
    }

    return 0;
}
