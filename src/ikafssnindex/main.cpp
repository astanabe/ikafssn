#include "io/blastdb_reader.hpp"
#include "index/index_builder.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "util/cli_parser.hpp"
#include "util/size_parser.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <filesystem>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n\n"
        "Required:\n"
        "  -db <path>             BLAST DB prefix\n"
        "  -k <int>               k-mer length (%d-%d)\n"
        "  -o <dir>               Output directory\n\n"
        "Options:\n"
        "  -buffer_size <size>    Buffer size (default: 8G)\n"
        "                         Accepts K, M, G suffixes\n"
        "  -partitions <int>      Number of partitions (default: 4)\n"
        "                         Powers of 2 recommended\n"
        "  -max_freq_build <int>  Exclude k-mers with count > threshold\n"
        "                         (default: 0 = no exclusion)\n"
        "  -threads <int>         Threads for counting pass (default: 1)\n"
        "  -v, --verbose          Verbose output\n",
        prog, MIN_K, MAX_K);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    // Check for help
    if (cli.has("-h") || cli.has("--help") || argc < 2) {
        print_usage(argv[0]);
        return (argc < 2) ? 1 : 0;
    }

    // Required arguments
    std::string db_path = cli.get_string("-db");
    int k = cli.get_int("-k", 0);
    std::string out_dir = cli.get_string("-o");

    if (db_path.empty()) {
        std::fprintf(stderr, "Error: -db is required\n");
        print_usage(argv[0]);
        return 1;
    }
    if (k == 0) {
        std::fprintf(stderr, "Error: -k is required\n");
        print_usage(argv[0]);
        return 1;
    }
    if (out_dir.empty()) {
        std::fprintf(stderr, "Error: -o is required\n");
        print_usage(argv[0]);
        return 1;
    }

    // Validate k
    if (k < MIN_K || k > MAX_K) {
        std::fprintf(stderr, "Error: k must be between %d and %d\n", MIN_K, MAX_K);
        return 1;
    }

    // Optional arguments
    std::string buf_size_str = cli.get_string("-buffer_size", "8G");
    uint64_t buffer_size = parse_size_string(buf_size_str);
    if (buffer_size == 0) {
        std::fprintf(stderr, "Error: invalid -buffer_size '%s'\n", buf_size_str.c_str());
        return 1;
    }

    int partitions = cli.get_int("-partitions", 4);
    if (partitions < 1 || (partitions & (partitions - 1)) != 0) {
        std::fprintf(stderr, "Error: -partitions must be a power of 2 (got %d)\n", partitions);
        return 1;
    }

    uint64_t max_freq_build = 0;
    if (cli.has("-max_freq_build")) {
        max_freq_build = static_cast<uint64_t>(cli.get_int("-max_freq_build", 0));
    }

    bool verbose = cli.has("-v") || cli.has("--verbose");

    Logger logger(verbose ? Logger::kDebug : Logger::kInfo);

    // Create output directory
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::fprintf(stderr, "Error: cannot create output directory '%s': %s\n",
                     out_dir.c_str(), ec.message().c_str());
        return 1;
    }

    // Find volume paths
    std::vector<std::string> vol_paths = BlastDbReader::find_volume_paths(db_path);
    if (vol_paths.empty()) {
        // Single-volume DB: try opening directly
        vol_paths.push_back(db_path);
    }

    logger.info("Database: %s (%zu volume(s))", db_path.c_str(), vol_paths.size());
    logger.info("Parameters: k=%d, buffer=%s, partitions=%d",
                k, buf_size_str.c_str(), partitions);

    // Extract DB base name from path
    std::string db_base = std::filesystem::path(db_path).filename().string();

    int threads = cli.get_int("-threads", 1);
    if (threads < 1) threads = 1;

    // Build config
    IndexBuilderConfig config;
    config.k = k;
    config.buffer_size = buffer_size;
    config.partitions = partitions;
    config.max_freq_build = max_freq_build;
    config.threads = threads;
    config.verbose = verbose;

    uint16_t total_volumes = static_cast<uint16_t>(vol_paths.size());

    // Process each volume sequentially
    for (uint16_t vi = 0; vi < total_volumes; vi++) {
        logger.info("=== Volume %d/%d: %s ===", vi + 1, total_volumes,
                    vol_paths[vi].c_str());

        BlastDbReader db;
        if (!db.open(vol_paths[vi])) {
            std::fprintf(stderr, "Error: cannot open volume '%s'\n",
                         vol_paths[vi].c_str());
            return 1;
        }

        // Output prefix: <out_dir>/<db_base>.<vol_idx>.<kk>mer
        char kk_str[8];
        std::snprintf(kk_str, sizeof(kk_str), "%02d", k);
        char vol_str[8];
        std::snprintf(vol_str, sizeof(vol_str), "%02d", vi);
        std::string prefix = out_dir + "/" + db_base + "." +
                             vol_str + "." + kk_str + "mer";

        bool ok;
        if (k < K_TYPE_THRESHOLD) {
            ok = build_index<uint16_t>(db, config, prefix,
                                        vi, total_volumes, db_base, logger);
        } else {
            ok = build_index<uint32_t>(db, config, prefix,
                                        vi, total_volumes, db_base, logger);
        }

        if (!ok) {
            std::fprintf(stderr, "Error: index build failed for volume %d\n", vi);
            return 1;
        }
    }

    logger.info("All volumes completed successfully.");
    return 0;
}
