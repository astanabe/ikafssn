#include "io/blastdb_reader.hpp"
#include "io/volume_discovery.hpp"
#include "index/index_builder.hpp"
#include "index/index_filter.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "util/cli_parser.hpp"
#include "util/common_init.hpp"
#include "util/size_parser.hpp"
#include "util/logger.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
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
        "  -mode <1|2|3>          Search mode the index will support (default: 2)\n"
        "                         1 = Stage 1 only (skip .kpx generation)\n"
        "                         2 = Stage 1+2 (default)\n"
        "                         3 = Stage 1+2+3 (same as 2 for index)\n"
        "  -memory_limit <size>   Memory limit (default: %s = half of RAM)\n"
        "                         Accepts K, M, G suffixes\n"
        "  -max_freq_build <num>  Exclude k-mers with cross-volume count > threshold\n"
        "                         >= 1: absolute count threshold\n"
        "                         0 < x < 1: fraction of total NSEQ across all volumes\n"
        "                         Counts are aggregated across all volumes before filtering\n"
        "                         (default: 0 = no exclusion)\n"
        "  -highfreq_filter_threads <int>\n"
        "                         Threads for cross-volume filtering (default: min(8, threads))\n"
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

    if (check_version(cli, "ikafssnindex")) return 0;

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

    // Parse -mode (1, 2, or 3; default 2)
    int index_mode = cli.get_int("-mode", 2);
    if (index_mode < 1 || index_mode > 3) {
        std::fprintf(stderr, "Error: -mode must be 1, 2, or 3\n");
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

    Logger logger = make_logger(cli);
    bool verbose = logger.verbose();

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

    int threads = resolve_threads(cli);

    int highfreq_filter_threads;
    if (cli.has("-highfreq_filter_threads")) {
        highfreq_filter_threads = cli.get_int("-highfreq_filter_threads", 8);
        if (highfreq_filter_threads < 1) {
            std::fprintf(stderr, "Error: -highfreq_filter_threads must be >= 1\n");
            return 1;
        }
        if (highfreq_filter_threads > threads) {
            std::fprintf(stderr,
                "Error: -highfreq_filter_threads (%d) exceeds -threads (%d)\n",
                highfreq_filter_threads, threads);
            return 1;
        }
    } else {
        highfreq_filter_threads = std::min(8, threads);
    }

    int openvol = cli.get_int("-openvol", 1);
    if (openvol < 1) openvol = 1;

    logger.info("Database: %s (%zu volume(s))", db_path.c_str(), vol_paths.size());
    logger.info("Parameters: k=%d, mode=%d, memory_limit=%s, openvol=%d, threads=%d",
                k, index_mode, mem_limit_str.c_str(), openvol, threads);

    // Extract DB base name from path
    std::string db_base = std::filesystem::path(db_path).filename().string();

    // Centralized TBB thread control
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, threads);

    // Build config (per-volume memory budget = total limit / openvol)
    IndexBuilderConfig config;
    config.k = k;
    config.memory_limit = memory_limit / static_cast<uint64_t>(openvol);
    config.threads = threads;
    config.verbose = verbose;
    config.skip_kpx = (index_mode == 1);
    // When max_freq_build is active, keep .tmp files for cross-volume filtering
    config.keep_tmp = (max_freq_build > 0);

    // Resolve fractional -max_freq_build to absolute threshold
    uint64_t freq_threshold = 0;
    if (max_freq_build > 0) {
        if (max_freq_build < 1.0) {
            uint64_t total_nseq = 0;
            for (const auto& vp : vol_paths) {
                BlastDbReader tmp_db;
                if (!tmp_db.open(vp)) {
                    std::fprintf(stderr, "Error: cannot open volume '%s' for NSEQ count\n",
                                 vp.c_str());
                    return 1;
                }
                total_nseq += tmp_db.num_sequences();
            }
            double resolved = std::ceil(max_freq_build * total_nseq);
            if (resolved < 1.0) resolved = 1.0;
            freq_threshold = static_cast<uint64_t>(resolved);
            logger.info("-max_freq_build=%.6g (fraction of total NSEQ=%lu) -> threshold=%lu",
                        max_freq_build, static_cast<unsigned long>(total_nseq),
                        static_cast<unsigned long>(freq_threshold));
        } else {
            freq_threshold = static_cast<uint64_t>(max_freq_build);
        }
    }

    uint16_t total_volumes = static_cast<uint16_t>(vol_paths.size());

    // Extract per-volume basenames from BLAST DB volume paths
    std::vector<std::string> vol_basenames(total_volumes);
    for (uint16_t vi = 0; vi < total_volumes; vi++) {
        vol_basenames[vi] = std::filesystem::path(vol_paths[vi]).filename().string();
    }

    // Check for duplicate basenames
    {
        std::set<std::string> seen;
        for (const auto& bn : vol_basenames) {
            if (!seen.insert(bn).second) {
                std::fprintf(stderr, "Error: duplicate volume basename '%s'\n", bn.c_str());
                return 1;
            }
        }
    }

    // Pre-compute per-volume output prefixes
    std::vector<std::string> vol_prefixes(total_volumes);
    for (uint16_t vi = 0; vi < total_volumes; vi++) {
        vol_prefixes[vi] = index_file_stem(out_dir, vol_basenames[vi], k);
    }

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

            const std::string& prefix = vol_prefixes[vi];

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

    // Write .kvx manifest
    {
        std::string kvx_path = index_file_stem(out_dir, db_base, k) + ".kvx";
        FILE* fp = std::fopen(kvx_path.c_str(), "w");
        if (!fp) {
            std::fprintf(stderr, "Error: cannot write %s\n", kvx_path.c_str());
            return 1;
        }
        std::fprintf(fp, "#\n# ikafssn index volume manifest\n#\n");
        std::fprintf(fp, "TITLE %s\n", db_base.c_str());
        std::fprintf(fp, "DBLIST");
        for (const auto& bn : vol_basenames) {
            std::fprintf(fp, " \"%s\"", bn.c_str());
        }
        std::fprintf(fp, "\n");
        std::fclose(fp);
        logger.info("Wrote volume manifest: %s", kvx_path.c_str());
    }

    // Post-build cross-volume frequency filtering
    if (max_freq_build > 0) {
        std::string khx_path = khx_path_for(out_dir, db_base, k);

        if (!filter_volumes_cross_volume(vol_prefixes, khx_path, k,
                                         freq_threshold, highfreq_filter_threads,
                                         logger)) {
            std::fprintf(stderr, "Error: cross-volume filtering failed\n");
            return 1;
        }
    }

    logger.info("All volumes completed successfully.");
    return 0;
}
