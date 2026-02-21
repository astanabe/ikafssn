#include "ikafssnserver/server.hpp"
#include "util/cli_parser.hpp"

#include <csignal>
#include <cstdio>

using namespace ikafssn;

static Server* g_server = nullptr;

static void signal_handler(int /*sig*/) {
    if (g_server) {
        g_server->request_shutdown();
    }
}

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Required:\n"
        "  -ix <dir>                Index directory\n"
        "\n"
        "Listener (at least one required):\n"
        "  -socket <path>           UNIX domain socket path\n"
        "  -tcp <host>:<port>       TCP listen address\n"
        "\n"
        "Options:\n"
        "  -threads <int>           Worker threads (default: all cores)\n"
        "  -pid <path>              PID file path\n"
        "  -min_score <int>         Default minimum chain score (default: 1)\n"
        "  -max_gap <int>           Default chaining gap tolerance (default: 100)\n"
        "  -max_freq <int>          Default high-freq k-mer skip threshold (default: auto)\n"
        "  -min_diag_hits <int>     Default diagonal filter min hits (default: 2)\n"
        "  -stage1_topn <int>       Default Stage 1 candidate limit (default: 0)\n"
        "  -min_stage1_score <num>  Default Stage 1 minimum score; integer or 0<P<1 fraction (default: 0.5)\n"
        "  -num_results <int>       Default max results per query (default: 0)\n"
        "  -shutdown_timeout <int>  Graceful shutdown timeout in seconds (default: 30)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (cli.has("-h") || cli.has("--help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (!cli.has("-ix")) {
        std::fprintf(stderr, "Error: -ix is required\n");
        print_usage(argv[0]);
        return 1;
    }

    if (!cli.has("-socket") && !cli.has("-tcp")) {
        std::fprintf(stderr, "Error: at least one of -socket or -tcp is required\n");
        print_usage(argv[0]);
        return 1;
    }

    ServerConfig config;
    config.ix_dir = cli.get_string("-ix");
    config.unix_socket_path = cli.get_string("-socket");
    config.tcp_addr = cli.get_string("-tcp");
    config.pid_file = cli.get_string("-pid");
    config.num_threads = cli.get_int("-threads", 0);
    config.shutdown_timeout = cli.get_int("-shutdown_timeout", 30);

    bool verbose = cli.has("-v") || cli.has("--verbose");
    config.log_level = verbose ? Logger::kDebug : Logger::kInfo;

    config.search_config.stage1.max_freq =
        static_cast<uint32_t>(cli.get_int("-max_freq", 0));
    config.search_config.stage1.stage1_topn =
        static_cast<uint32_t>(cli.get_int("-stage1_topn", 0));
    {
        double min_s1 = cli.get_double("-min_stage1_score", 0.5);
        if (min_s1 > 0 && min_s1 < 1.0) {
            config.search_config.min_stage1_score_frac = min_s1;
        } else {
            config.search_config.stage1.min_stage1_score = static_cast<uint32_t>(min_s1);
        }
    }
    config.search_config.stage2.max_gap =
        static_cast<uint32_t>(cli.get_int("-max_gap", 100));
    config.search_config.stage2.min_diag_hits =
        static_cast<uint32_t>(cli.get_int("-min_diag_hits", 2));
    config.search_config.stage2.min_score =
        static_cast<uint32_t>(cli.get_int("-min_score", 1));
    config.search_config.num_results =
        static_cast<uint32_t>(cli.get_int("-num_results", 0));

    Server server;
    g_server = &server;

    // Install signal handlers
    struct sigaction sa;
    sa.sa_handler = signal_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGTERM, &sa, nullptr);
    sigaction(SIGINT, &sa, nullptr);

    int ret = server.run(config);
    g_server = nullptr;
    return ret;
}
