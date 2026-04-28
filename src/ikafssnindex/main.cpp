#include "io/blastdb_reader.hpp"
#include "io/volume_discovery.hpp"
#include "index/index_builder.hpp"
#include "index/index_filter.hpp"
#include "index/kix_format.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"
#include "core/spaced_seed.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "util/cli_parser.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "util/size_parser.hpp"
#include "util/logger.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

// Result of validating existing index files against a BLAST DB volume.
enum class IndexStatus {
    kNotFound,   // No complete index files exist (needs building)
    kValid,      // Index exists and matches BLAST DB
    kMismatch    // Index exists but metadata doesn't match BLAST DB
};

struct IndexValidation {
    IndexStatus status = IndexStatus::kNotFound;
    std::string detail;  // mismatch description (for error messages)
};

// Validate existing index files for a volume against the BLAST DB.
// Checks: file existence, header validity, num_sequences, total_bases.
// suffix: "" for final files (.kix), ".tmp" for temp files (.kix.tmp).
static IndexValidation validate_existing_index(
    const std::string& output_prefix,
    const std::string& vol_path,
    int k, uint8_t t, uint8_t template_type,
    bool skip_kpx,
    const char* suffix = "") {

    IndexValidation result;
    std::string kix_path = output_prefix + ".kix" + suffix;
    std::string ksx_path = output_prefix + ".ksx" + suffix;

    // Check .kix exists and read header
    FILE* fp = std::fopen(kix_path.c_str(), "rb");
    if (!fp) return result; // kNotFound
    KixHeader kix_hdr{};
    bool ok = (std::fread(&kix_hdr, sizeof(kix_hdr), 1, fp) == 1);
    std::fclose(fp);
    if (!ok) return result; // kNotFound (truncated file)

    if (std::memcmp(kix_hdr.magic, KIX_MAGIC, 4) != 0) return result;
    if (kix_hdr.format_version != KIX_FORMAT_VERSION) return result;

    // Check .ksx exists, read header + seq_lengths to compute total_bases
    fp = std::fopen(ksx_path.c_str(), "rb");
    if (!fp) return result; // kNotFound
    KsxHeader ksx_hdr{};
    ok = (std::fread(&ksx_hdr, sizeof(ksx_hdr), 1, fp) == 1);
    if (!ok) { std::fclose(fp); return result; }

    if (std::memcmp(ksx_hdr.magic, KSX_MAGIC, 4) != 0) { std::fclose(fp); return result; }
    if (ksx_hdr.format_version != KSX_FORMAT_VERSION) { std::fclose(fp); return result; }

    uint32_t ksx_nseq = ksx_hdr.num_sequences;
    std::vector<uint32_t> seq_lengths(ksx_nseq);
    if (ksx_nseq > 0) {
        size_t read_count = std::fread(seq_lengths.data(), sizeof(uint32_t), ksx_nseq, fp);
        if (read_count != ksx_nseq) { std::fclose(fp); return result; }
    }
    std::fclose(fp);

    uint64_t ksx_total_bases = 0;
    for (uint32_t i = 0; i < ksx_nseq; i++) {
        ksx_total_bases += seq_lengths[i];
    }

    // Check .kpx exists (if required)
    if (!skip_kpx) {
        std::string kpx_path = output_prefix + ".kpx" + suffix;
        if (!std::filesystem::exists(kpx_path)) return result; // kNotFound
    }

    // At this point, all index files exist with valid headers.
    // Now compare against BLAST DB.
    BlastDbReader db;
    if (!db.open(vol_path)) {
        result.status = IndexStatus::kMismatch;
        result.detail = "cannot open BLAST DB '" + vol_path + "' for validation";
        return result;
    }
    uint32_t db_nseq = db.num_sequences();
    uint64_t db_total_bases = db.total_length();

    // Check k, t, template_type
    if (kix_hdr.k != static_cast<uint8_t>(k) ||
        kix_hdr.t != t ||
        kix_hdr.template_type != template_type) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: index parameters (k=%d, t=%d, template_type=%d) do not match "
            "requested (k=%d, t=%d, template_type=%d)",
            kix_path.c_str(),
            kix_hdr.k, kix_hdr.t, kix_hdr.template_type,
            k, t, template_type);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    // Check num_sequences
    if (kix_hdr.num_sequences != db_nseq || ksx_nseq != db_nseq) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: num_sequences=%u but BLAST DB '%s' has %u",
            kix_path.c_str(), kix_hdr.num_sequences,
            vol_path.c_str(), db_nseq);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    // Check total_bases
    if (ksx_total_bases != db_total_bases) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: total_bases=%lu but BLAST DB '%s' has %lu",
            ksx_path.c_str(), static_cast<unsigned long>(ksx_total_bases),
            vol_path.c_str(), static_cast<unsigned long>(db_total_bases));
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    result.status = IndexStatus::kValid;
    return result;
}

// Remove .tmp files for a volume prefix.
static void cleanup_tmp_files(const std::string& output_prefix, bool skip_kpx,
                              const Logger& logger) {
    auto remove_tmp = [&logger](const std::string& path) {
        if (std::filesystem::exists(path)) {
            std::remove(path.c_str());
            logger.info("Removed incomplete temp file: %s", path.c_str());
        }
    };
    remove_tmp(output_prefix + ".kix.tmp");
    remove_tmp(output_prefix + ".ksx.tmp");
    if (!skip_kpx) {
        remove_tmp(output_prefix + ".kpx.tmp");
    }
}

static void print_usage(const char* prog, const std::string& default_mem) {
    print_version_header("ikafssnindex");
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
        "  -rebuild <0|1>         Rebuild mode (default: 0)\n"
        "                         0 = skip volumes with valid existing indexes\n"
        "                         1 = always rebuild all volumes from scratch\n"
        "  -memory_limit <size>   Memory limit (default: %s = half of RAM)\n"
        "                         Accepts K, M, G suffixes\n"
        "  -max_freq_build <num>  Exclude k-mers with cross-volume count > threshold\n"
        "                         1 or 1.0: disable (no exclusion, default)\n"
        "                         0 < x < 1: fraction of total NSEQ across all volumes\n"
        "                         > 1: absolute count threshold\n"
        "                         0: not allowed (error)\n"
        "                         Counts are aggregated across all volumes before filtering\n"
        "  -highfreq_filter_threads <int>\n"
        "                         Threads for cross-volume filtering (default: min(8, threads))\n"
        "  -max_degen_expand <int>  Max degenerate expansion per k-mer (default: 4, max: 16, 0/1: disable)\n"
        "  -t <int>               Template length for spaced seeds\n"
        "                         0: contiguous k-mers (default)\n"
        "                         13, 15, 18: requires -k 8 or 9\n"
        "                         16, 18, 21: requires -k 11 or 12\n"
        "  -template_type <str>   Template type: coding, optimal, or both (required with -t)\n"
        "                         both: builds coding and optimal indexes sequentially\n"
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
    std::string default_mem_str = format_size(default_mem);

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

    double max_freq_build = 1.0; // default: disabled (no exclusion)
    if (cli.has("-max_freq_build")) {
        max_freq_build = cli.get_double("-max_freq_build", 1.0);
        if (max_freq_build == 0) {
            std::fprintf(stderr, "Error: -max_freq_build=0 is not allowed. "
                "Use 1 to disable or specify a threshold (0 < x < 1 for fraction, > 1 for absolute)\n");
            return 1;
        }
        if (max_freq_build < 0) {
            std::fprintf(stderr, "Error: -max_freq_build must be > 0\n");
            return 1;
        }
    }

    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);
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

    int rebuild_mode = cli.get_int("-rebuild", 0);
    if (rebuild_mode != 0 && rebuild_mode != 1) {
        std::fprintf(stderr, "Error: -rebuild must be 0 or 1\n");
        return 1;
    }

    int max_degen_expand = cli.get_int("-max_degen_expand", 4);
    if (max_degen_expand < 0 || max_degen_expand > 16) {
        std::fprintf(stderr, "Error: -max_degen_expand must be between 0 and 16\n");
        return 1;
    }

    int cli_t = cli.get_int("-t", 0);
    if (cli_t != 0 && cli_t != 13 && cli_t != 15 && cli_t != 16 && cli_t != 18 && cli_t != 21) {
        std::fprintf(stderr, "Error: -t must be 0, 13, 15, 16, 18, or 21\n");
        return 1;
    }
    uint8_t spaced_t = static_cast<uint8_t>(cli_t);

    TemplateType spaced_type = TemplateType::kContiguous;
    if (cli.has("-template_type")) {
        spaced_type = template_type_from_string(cli.get_string("-template_type"));
        if (spaced_type == TemplateType::kContiguous) {
            std::fprintf(stderr, "Error: -template_type must be coding, optimal, or both\n");
            return 1;
        }
    }

    if (spaced_t > 0) {
        if (!cli.has("-template_type")) {
            std::fprintf(stderr,
                "Error: -template_type (coding, optimal, or both) is required when -t is specified\n");
            return 1;
        }
        if (!validate_spaced_seed(k, spaced_t)) {
            std::fprintf(stderr, "Error: -t %d is not valid for -k %d\n", spaced_t, k);
            return 1;
        }
    }

    // Determine which template types to build.
    // "both" builds coding and optimal indexes sequentially.
    std::vector<TemplateType> build_types;
    if (spaced_type == TemplateType::kBoth) {
        build_types = {TemplateType::kCoding, TemplateType::kOptimal};
    } else {
        build_types = {spaced_type};
    }

    logger.info("Database: %s (%zu volume(s))", db_path.c_str(), vol_paths.size());
    if (spaced_t > 0) {
        std::string type_display = (spaced_type == TemplateType::kBoth)
            ? "both (coding+optimal)" : template_type_to_string(spaced_type);
        logger.info("Parameters: k=%d, t=%d, template_type=%s, mode=%d, memory_limit=%s, openvol=%d, threads=%d",
                    k, spaced_t, type_display.c_str(),
                    index_mode, mem_limit_str.c_str(), openvol, threads);
    } else {
        logger.info("Parameters: k=%d, mode=%d, memory_limit=%s, openvol=%d, threads=%d",
                    k, index_mode, mem_limit_str.c_str(), openvol, threads);
    }

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
    config.max_degen_expand = max_degen_expand;
    config.t = spaced_t;
    config.template_type = static_cast<uint8_t>(spaced_type);
    // When max_freq_build is active (not 1.0 = disabled), keep .tmp files for cross-volume filtering
    bool freq_filter_active = (max_freq_build != 1.0);
    config.keep_tmp = freq_filter_active;

    // Resolve fractional -max_freq_build to absolute threshold
    uint64_t freq_threshold = 0;
    if (freq_filter_active) {
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

    // Build each template type (for "both", builds coding then optimal sequentially)
    for (TemplateType cur_type : build_types) {
        uint8_t cur_tt = static_cast<uint8_t>(cur_type);

        if (build_types.size() > 1) {
            logger.info("========== Building %s template ==========",
                        template_type_to_string(cur_type).c_str());
        }

        // Update config for this template type
        config.template_type = cur_tt;

        // Pre-compute per-volume output prefixes for this template type
        std::vector<std::string> vol_prefixes(total_volumes);
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            vol_prefixes[vi] = index_file_stem(out_dir, vol_basenames[vi], k, spaced_t, cur_tt);
        }

        // --- Resume check: determine which volumes to skip ---
        std::vector<bool> skip_volume(total_volumes, false);
        bool build_skipped = false;  // set when all volumes are skipped

        if (rebuild_mode == 0) {
            // Validate each volume's existing index
            std::vector<IndexValidation> validations(total_volumes);
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                validations[vi] = validate_existing_index(vol_prefixes[vi], vol_paths[vi],
                    k, spaced_t, cur_tt, config.skip_kpx);
            }

            // Check for mismatches — error out immediately
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (validations[vi].status == IndexStatus::kMismatch) {
                    std::fprintf(stderr,
                        "Error: existing index does not match BLAST DB:\n  %s\n"
                        "Delete the mismatched index file(s) or re-run with -rebuild 1\n",
                        validations[vi].detail.c_str());
                    return 1;
                }
            }

            if (freq_filter_active) {
                // With cross-volume filtering: skip only if ALL volumes have
                // valid final files AND .khx exists (entire pipeline completed).
                bool all_final_valid = true;
                for (uint16_t vi = 0; vi < total_volumes; vi++) {
                    if (validations[vi].status != IndexStatus::kValid) {
                        all_final_valid = false;
                        break;
                    }
                }
                std::string khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt);
                if (all_final_valid && std::filesystem::exists(khx_path)) {
                    logger.info("All %d volumes have valid indexes and .khx exists; "
                                "skipping build and filter (-rebuild 0)",
                                total_volumes);
                    build_skipped = true;
                } else {
                    // Not all final+.khx: check .tmp files for per-volume resume.
                    // Filter needs .tmp from ALL volumes, but we can reuse valid .tmp
                    // and only rebuild the missing/invalid ones.
                    for (uint16_t vi = 0; vi < total_volumes; vi++) {
                        auto tmp_val = validate_existing_index(vol_prefixes[vi], vol_paths[vi],
                            k, spaced_t, cur_tt, config.skip_kpx, ".tmp");
                        if (tmp_val.status == IndexStatus::kValid) {
                            skip_volume[vi] = true;
                            logger.info("Volume %d/%d (%s): valid .tmp index exists, "
                                        "skipping build (-rebuild 0)",
                                        vi + 1, total_volumes, vol_basenames[vi].c_str());
                        } else {
                            // Invalid or missing .tmp — clean up and rebuild
                            cleanup_tmp_files(vol_prefixes[vi], config.skip_kpx, logger);
                        }
                    }
                }
            } else {
                // Without cross-volume filtering: per-volume resume.
                for (uint16_t vi = 0; vi < total_volumes; vi++) {
                    if (validations[vi].status == IndexStatus::kValid) {
                        skip_volume[vi] = true;
                        logger.info("Volume %d/%d (%s): valid index exists, skipping (-rebuild 0)",
                                    vi + 1, total_volumes, vol_basenames[vi].c_str());
                    } else {
                        cleanup_tmp_files(vol_prefixes[vi], config.skip_kpx, logger);
                    }
                }

                // Check if all volumes are already done
                bool all_skipped = true;
                for (uint16_t vi = 0; vi < total_volumes; vi++) {
                    if (!skip_volume[vi]) { all_skipped = false; break; }
                }
                if (all_skipped) {
                    logger.info("All %d volumes have valid indexes; skipping build (-rebuild 0)",
                                total_volumes);
                    build_skipped = true;
                }
            }
        } else {
            // rebuild_mode == 1: clean up any .tmp files
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                cleanup_tmp_files(vol_prefixes[vi], config.skip_kpx, logger);
            }
        }

        if (!build_skipped) {
        // Process volumes via TBB task_group with concurrency limited by -openvol.
        {
        std::atomic<bool> any_error{false};
        std::vector<std::string> error_messages(total_volumes);
        std::mutex log_mutex;

        int max_active = std::min(openvol, static_cast<int>(total_volumes));
        std::mutex vol_mutex;
        std::condition_variable vol_cv;
        int active_volumes = 0;

        tbb::task_group tg;
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            if (skip_volume[vi]) continue;

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
                if (kmer_type_for(k, spaced_t) == 0) {
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
        } // end volume processing block
        } // end if (!build_skipped)

        // Write .kvx manifest for this template type
        {
            std::string kvx_path = index_file_stem(out_dir, db_base, k, spaced_t, cur_tt) + ".kvx";
            FILE* fp = std::fopen(kvx_path.c_str(), "w");
            if (!fp) {
                std::fprintf(stderr, "Error: cannot write %s\n", kvx_path.c_str());
                return 1;
            }
            std::fprintf(fp, "#\n# ikafssn index volume manifest\n#\n");
            std::fprintf(fp, "FORMAT_VERSION 5\n");
            std::fprintf(fp, "TITLE %s\n", db_base.c_str());
            std::fprintf(fp, "DBLIST");
            for (const auto& bn : vol_basenames) {
                std::fprintf(fp, " \"%s\"", bn.c_str());
            }
            std::fprintf(fp, "\n");
            std::fclose(fp);
            logger.info("Wrote volume manifest: %s", kvx_path.c_str());
        }

        // Post-build cross-volume frequency filtering for this template type
        if (freq_filter_active && !build_skipped) {
            std::string khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt);

            if (!filter_volumes_cross_volume(vol_prefixes, khx_path, k,
                                             freq_threshold, highfreq_filter_threads,
                                             logger)) {
                std::fprintf(stderr, "Error: cross-volume filtering failed\n");
                return 1;
            }
        }

        if (build_types.size() > 1) {
            logger.info("========== %s template completed ==========",
                        template_type_to_string(cur_type).c_str());
        }
    } // end for each build_type

    logger.info("All volumes completed successfully.");
    return 0;
}
