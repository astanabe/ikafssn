#include "ikafssnserver/server.hpp"
#include "core/version.hpp"
#include "util/cli_parser.hpp"
#include "util/common_init.hpp"
#include "util/context_parser.hpp"

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
        "  -ix <prefix>             Index prefix (like blastn -db)\n"
        "\n"
        "Listener (at least one required):\n"
        "  -socket <path>           UNIX domain socket path\n"
        "  -tcp <host>:<port>       TCP listen address\n"
        "\n"
        "Options:\n"
        "  -threads <int>           Worker threads (default: all cores)\n"
        "  -max_query <int>         Max concurrent query sequences globally (default: 1024)\n"
        "  -max_seqs_per_req <int>  Max sequences accepted per request (default: thread count)\n"
        "  -pid <path>              PID file path\n"
        "  -db <path>               BLAST DB path for mode 3 (default: same as -ix)\n"
        "  -stage2_min_score <int>  Default minimum chain score (default: 0 = adaptive)\n"
        "  -stage2_max_gap <int>    Default chaining gap tolerance (default: 100)\n"
        "  -stage2_max_lookback <int>  Default chaining DP lookback window (default: 64, 0=unlimited)\n"
        "  -stage1_max_freq <num>   Default high-freq k-mer skip threshold (default: 0.5)\n"
        "                           0 < x < 1: fraction of total NSEQ across all volumes\n"
        "                           >= 1: absolute count threshold; 0 = auto\n"
        "  -stage2_min_diag_hits <int>  Default diagonal filter min hits (default: 1)\n"
        "  -stage1_topn <int>       Default Stage 1 candidate limit (default: 0)\n"
        "  -stage1_min_score <num>  Default Stage 1 minimum score; integer or 0<P<1 fraction (default: 0.5)\n"
        "  -num_results <int>       Default max results per query (default: 0)\n"
        "  -accept_qdegen <0|1>     Default accept queries with degenerate bases (default: 1)\n"
        "  -context <value>         Default context extension (int=bases, decimal=ratio, default: 0)\n"
        "  -stage3_traceback <0|1>  Default traceback mode (default: 0)\n"
        "  -stage3_gapopen <int>    Default gap open penalty (default: 10)\n"
        "  -stage3_gapext <int>     Default gap extension penalty (default: 1)\n"
        "  -stage3_min_pident <num> Default min percent identity (default: 0)\n"
        "  -stage3_min_nident <int> Default min identical bases (default: 0)\n"
        "  -stage3_fetch_threads <int>  Threads for BLAST DB fetch (default: min(8, threads))\n"
        "  -shutdown_timeout <int>  Graceful shutdown timeout in seconds (default: 30)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (check_version(cli, "ikafssnserver")) return 0;

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
    config.ix_prefix = cli.get_string("-ix");
    config.unix_socket_path = cli.get_string("-socket");
    config.tcp_addr = cli.get_string("-tcp");
    config.pid_file = cli.get_string("-pid");
    config.num_threads = cli.get_int("-threads", 0);
    config.max_query = cli.get_int("-max_query", 0);
    if (config.max_query < 0) {
        std::fprintf(stderr, "Error: -max_query must be >= 0\n");
        return 1;
    }
    config.max_seqs_per_req = cli.get_int("-max_seqs_per_req", 0);
    if (config.max_seqs_per_req < 0) {
        std::fprintf(stderr, "Error: -max_seqs_per_req must be >= 0\n");
        return 1;
    }
    config.shutdown_timeout = cli.get_int("-shutdown_timeout", 30);

    config.log_level = make_logger(cli).level();

    // Search config (renamed options)
    config.max_freq_raw = cli.get_double("-stage1_max_freq", 0.5);
    config.search_config.stage1.stage1_topn =
        static_cast<uint32_t>(cli.get_int("-stage1_topn", 0));
    {
        double min_s1 = cli.get_double("-stage1_min_score", 0.5);
        if (min_s1 > 0 && min_s1 < 1.0) {
            config.search_config.min_stage1_score_frac = min_s1;
        } else {
            config.search_config.stage1.min_stage1_score = static_cast<uint32_t>(min_s1);
        }
    }
    config.search_config.stage2.max_gap =
        static_cast<uint32_t>(cli.get_int("-stage2_max_gap", 100));
    config.search_config.stage2.chain_max_lookback =
        static_cast<uint32_t>(cli.get_int("-stage2_max_lookback", 64));
    config.search_config.stage2.min_diag_hits =
        static_cast<uint32_t>(cli.get_int("-stage2_min_diag_hits", 1));
    config.search_config.stage2.min_score =
        static_cast<uint32_t>(cli.get_int("-stage2_min_score", 0));
    config.search_config.num_results =
        static_cast<uint32_t>(cli.get_int("-num_results", 0));
    config.search_config.accept_qdegen =
        static_cast<uint8_t>(cli.get_int("-accept_qdegen", 1));

    // BLAST DB path for mode 3 (default: same as index prefix)
    config.db_path = cli.get_string("-db", config.ix_prefix);

    // Stage 3 config
    config.stage3_config.gapopen = cli.get_int("-stage3_gapopen", 10);
    config.stage3_config.gapext = cli.get_int("-stage3_gapext", 1);
    config.stage3_config.traceback = (cli.get_int("-stage3_traceback", 0) != 0);
    config.stage3_config.min_pident = cli.get_double("-stage3_min_pident", 0.0);
    config.stage3_config.min_nident = static_cast<uint32_t>(cli.get_int("-stage3_min_nident", 0));

    // Resolve fetch_threads after num_threads is known
    int num_threads_resolved = resolve_threads(cli);
    if (cli.has("-stage3_fetch_threads")) {
        config.stage3_config.fetch_threads = cli.get_int("-stage3_fetch_threads", 8);
        if (config.stage3_config.fetch_threads > num_threads_resolved) {
            std::fprintf(stderr,
                "Error: -stage3_fetch_threads (%d) exceeds -threads (%d)\n",
                config.stage3_config.fetch_threads, num_threads_resolved);
            return 1;
        }
    } else {
        config.stage3_config.fetch_threads = std::min(8, num_threads_resolved);
    }

    // Parse -context: integer (bases) or decimal (query length multiplier)
    {
        ContextParam ctx_param;
        std::string err;
        if (!parse_context(cli.get_string("-context", "0"), ctx_param, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        config.context_is_ratio = ctx_param.is_ratio;
        config.context_ratio = ctx_param.ratio;
        config.context_abs = ctx_param.abs;
    }

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
