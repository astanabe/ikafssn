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
#include <unistd.h>
#include <vector>
#include <filesystem>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <tbb/global_control.h>
#include <tbb/task_group.h>

using namespace ikafssn;

// Detect physical memory and return half of it (minimum 1 GB).
static uint64_t default_memory_limit() {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page_size > 0) {
        uint64_t half = static_cast<uint64_t>(pages) * static_cast<uint64_t>(page_size) / 2;
        if (half >= (uint64_t(1) << 30)) return half;
    }
    return uint64_t(1) << 30; // fallback: 1 GB
}

static void print_usage(const char* prog, const std::string& default_mem) {
    std::fprintf(stderr,
        "Usage: %s [options]\n\n"
        "Required:\n"
        "  -db <path>             BLAST DB prefix\n"
        "  -k <int>               k-mer length (%d-%d)\n"
        "  -o <dir>               Output directory\n\n"
        "Options:\n"
        "  -memory_limit <size>   Memory limit (default: %s = half of RAM)\n"
        "                         Accepts K, M, G suffixes\n"
        "  -max_freq_build <num>  Exclude k-mers with count > threshold\n"
        "                         >= 1: absolute count threshold\n"
        "                         0 < x < 1: fraction of NSEQ per volume\n"
        "                         (default: 0 = no exclusion)\n"
        "  -openvol <int>         Max volumes processed simultaneously\n"
        "                         (default: 1)\n"
        "  -threads <int>         Number of threads (default: all cores)\n"
        "  -v, --verbose          Verbose output\n",
        prog, MIN_K, MAX_K, default_mem.c_str());
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    // Compute default memory limit string for help display
    uint64_t default_mem = default_memory_limit();
    std::string default_mem_str;
    if (default_mem >= (uint64_t(1) << 30) && default_mem % (uint64_t(1) << 30) == 0)
        default_mem_str = std::to_string(default_mem >> 30) + "G";
    else
        default_mem_str = std::to_string(default_mem >> 20) + "M";

    // Check for help
    if (cli.has("-h") || cli.has("--help") || argc < 2) {
        print_usage(argv[0], default_mem_str);
        return (argc < 2) ? 1 : 0;
    }

    // Required arguments
    std::string db_path = cli.get_string("-db");
    int k = cli.get_int("-k", 0);
    std::string out_dir = cli.get_string("-o");

    if (db_path.empty()) {
        std::fprintf(stderr, "Error: -db is required\n");
        print_usage(argv[0], default_mem_str);
        return 1;
    }
    if (k == 0) {
        std::fprintf(stderr, "Error: -k is required\n");
        print_usage(argv[0], default_mem_str);
        return 1;
    }
    if (out_dir.empty()) {
        std::fprintf(stderr, "Error: -o is required\n");
        print_usage(argv[0], default_mem_str);
        return 1;
    }

    // Validate k
    if (k < MIN_K || k > MAX_K) {
        std::fprintf(stderr, "Error: k must be between %d and %d\n", MIN_K, MAX_K);
        return 1;
    }

    // Optional arguments
    uint64_t memory_limit;
    std::string mem_limit_str;
    if (cli.has("-memory_limit")) {
        mem_limit_str = cli.get_string("-memory_limit");
        memory_limit = parse_size_string(mem_limit_str);
        if (memory_limit == 0) {
            std::fprintf(stderr, "Error: invalid -memory_limit '%s'\n", mem_limit_str.c_str());
            return 1;
        }
    } else {
        memory_limit = default_mem;
        mem_limit_str = default_mem_str;
    }

    double max_freq_build = 0;
    if (cli.has("-max_freq_build")) {
        max_freq_build = cli.get_double("-max_freq_build", 0);
        if (max_freq_build < 0) {
            std::fprintf(stderr, "Error: -max_freq_build must be >= 0\n");
            return 1;
        }
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

    int threads = cli.get_int("-threads", 0);
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
        if (threads <= 0) threads = 1;
    }

    int openvol = cli.get_int("-openvol", 1);
    if (openvol < 1) openvol = 1;

    logger.info("Database: %s (%zu volume(s))", db_path.c_str(), vol_paths.size());
    logger.info("Parameters: k=%d, memory_limit=%s, openvol=%d, threads=%d",
                k, mem_limit_str.c_str(), openvol, threads);

    // Extract DB base name from path
    std::string db_base = std::filesystem::path(db_path).filename().string();

    // Centralized TBB thread control
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, threads);

    // Build config (per-volume memory budget = total limit / openvol)
    IndexBuilderConfig config;
    config.k = k;
    config.memory_limit = memory_limit / static_cast<uint64_t>(openvol);
    config.max_freq_build = max_freq_build;
    config.threads = threads;
    config.verbose = verbose;

    uint16_t total_volumes = static_cast<uint16_t>(vol_paths.size());

    // Process volumes via TBB task_group with concurrency limited by -openvol.
    // The main thread gates submission: it waits until a slot is available
    // before submitting the next volume task, so at most openvol volumes
    // are active simultaneously.
    std::atomic<bool> any_error{false};
    std::vector<std::string> error_messages(total_volumes);
    std::mutex log_mutex;

    int max_active = std::min(openvol, static_cast<int>(total_volumes));
    std::mutex vol_mutex;
    std::condition_variable vol_cv;
    int active_volumes = 0;

    tbb::task_group tg;
    for (uint16_t vi = 0; vi < total_volumes; vi++) {
        // Wait until a slot is available
        {
            std::unique_lock<std::mutex> lock(vol_mutex);
            vol_cv.wait(lock, [&] { return active_volumes < max_active; });
            active_volumes++;
        }

        if (any_error.load(std::memory_order_relaxed)) break;

        tg.run([&, vi]() {
            if (any_error.load(std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(vol_mutex);
                active_volumes--;
                vol_cv.notify_one();
                return;
            }

            {
                std::lock_guard<std::mutex> lock(log_mutex);
                logger.info("=== Volume %d/%d: %s ===", vi + 1, total_volumes,
                            vol_paths[vi].c_str());
            }

            BlastDbReader db;
            if (!db.open(vol_paths[vi])) {
                error_messages[vi] = "cannot open volume '" + vol_paths[vi] + "'";
                any_error.store(true, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(vol_mutex);
                active_volumes--;
                vol_cv.notify_one();
                return;
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
                error_messages[vi] = "index build failed for volume " + std::to_string(vi);
                any_error.store(true, std::memory_order_relaxed);
            }

            {
                std::lock_guard<std::mutex> lock(vol_mutex);
                active_volumes--;
            }
            vol_cv.notify_one();
        });
    }
    tg.wait();

    if (any_error.load()) {
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            if (!error_messages[vi].empty()) {
                std::fprintf(stderr, "Error: %s\n", error_messages[vi].c_str());
            }
        }
        return 1;
    }

    logger.info("All volumes completed successfully.");
    return 0;
}
