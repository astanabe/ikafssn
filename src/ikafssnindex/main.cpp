#include "io/blastdb_reader.hpp"
#include "io/volume_discovery.hpp"
#include "index/index_builder.hpp"
#include "index/index_filter.hpp"
#include "index/kix_format.hpp"
#include "index/ksx_reader.hpp"
#include "index/volume_validator.hpp"
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
#include <set>
#include <string>
#include <vector>
#include <filesystem>

#include <tbb/global_control.h>
using namespace ikafssn;

// Per-volume resume state.  Strict validation failures drop a volume
// one level (after deleting the bad files): kComplete → kPostingsTmp
// → kMetadataTmp → kNone.
enum class VolumeState {
    kComplete,     // final .ksx / .kix (and .kpx) pass strict validation
    kPostingsTmp,  // .ksx.tmp / .kix.tmp (and .kpx.tmp) pass strict validation
    kMetadataTmp,  // .ksx.tmp passes strict validation; no usable .kix.tmp
    kNone          // nothing reusable (or -force_rebuild 1)
};

enum class KhxState { kValid, kNeedsRebuild };

static void remove_if_exists(const std::string& path, const Logger& logger) {
    if (std::filesystem::exists(path)) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        if (ec) {
            logger.warn("Could not remove %s: %s",
                        path.c_str(), ec.message().c_str());
        } else {
            logger.info("Removed: %s", path.c_str());
        }
    }
}

static void remove_volume_tmp(const std::string& prefix, bool skip_kpx,
                              const Logger& logger) {
    remove_if_exists(prefix + ".ksx.tmp", logger);
    remove_if_exists(prefix + ".kix.tmp", logger);
    if (!skip_kpx) remove_if_exists(prefix + ".kpx.tmp", logger);
}

static void remove_volume_final(const std::string& prefix, bool skip_kpx,
                                const Logger& logger) {
    remove_if_exists(prefix + ".ksx", logger);
    remove_if_exists(prefix + ".kix", logger);
    if (!skip_kpx) remove_if_exists(prefix + ".kpx", logger);
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
        "  -mode <1|2|3>          Search mode the index will support (default: 1)\n"
        "                         1 = Stage 1 only (skip .kpx generation, default)\n"
        "                         2 = Stage 1+2\n"
        "                         3 = Stage 1+2+3 (same as 2 for index)\n"
        "  -min_seq_length <int>  Minimum sequence length; shorter sequences are skipped\n"
        "                         (default: 64)\n"
        "  -min_length_split <int>\n"
        "                         Minimum split length; long sequences are split into\n"
        "                         overlapping fragments of this size (default: 50000\n"
        "                         for -mode 1, 0 for -mode 2/3).  0 disables splitting.\n"
        "                         When non-zero, must be in [1000, 1000000].\n"
        "  -overlap_length <int>  Overlap between adjacent fragments, in bases\n"
        "                         (default: 500 for -mode 1, 0 for -mode 2/3).\n"
        "                         Must be < min_length_split / 2.\n"
        "                         -min_length_split and -overlap_length must both be 0\n"
        "                         (splitting disabled) or both non-zero.\n"
        "  -force_rebuild <0|1>   Force a full rebuild (default: 0)\n"
        "                         0 = per-volume resume: validate existing files\n"
        "                             and rebuild only the missing or corrupted\n"
        "                             volumes\n"
        "                         1 = delete tmp / final files for this build's\n"
        "                             parameters and rebuild every volume\n"
        "  -memory_limit <size>   Per-volume sort buffer budget (default: %s = half of RAM)\n"
        "                         Bounds peak RAM during the postings pass by\n"
        "                         partitioning k-mer entries to fit this budget.\n"
        "                         Accepts K, M, G suffixes\n"
        "  -max_freq_build <num>  Exclude k-mers with cross-volume count > threshold\n"
        "                         1 or 1.0: disable (no exclusion, default)\n"
        "                         0 < x < 1: fraction of total NSEQ across all volumes\n"
        "                         > 1: absolute count threshold\n"
        "                         0: not allowed (error)\n"
        "                         Counts are aggregated across all volumes before filtering\n"
        "  -freq_threshold_part <int>\n"
        "                         .kpx per-(kmer, seq_id) partition threshold (default: 8, max: 255)\n"
        "                         A (k-mer, seq_id) cluster with occurrence count > threshold\n"
        "                         is split into its own partition group; lower-multiplicity\n"
        "                         clusters merge into a shared short bucket\n"
        "  -nthread_highfreq_filter <int>\n"
        "                         Threads for cross-volume filtering (default: min(8, nthread))\n"
        "  -max_degen_expand <int>  Max degenerate expansion per k-mer (default: 4, max: 16, 0/1: disable)\n"
        "  -t <int>               Template length for spaced seeds\n"
        "                         0: contiguous k-mers (default)\n"
        "                         13, 15, 18: requires -k 8 or 9\n"
        "                         16, 18, 21: requires -k 11 or 12\n"
        "  -template_type <str>   Template type: coding, optimal, or both (required with -t)\n"
        "                         both: builds coding and optimal indexes sequentially\n"
        "  -nthread <int>         Number of threads (default: all cores)\n"
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
    if (cli.has("-h") || cli.has("-help") || argc < 2) {
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

    // -mode 1 (Stage 1 only) skips .kpx generation and is the
    // recommended default; -mode 2/3 build the full index.
    int index_mode = cli.get_int("-mode", 1);
    if (index_mode < 1 || index_mode > 3) {
        std::fprintf(stderr, "Error: -mode must be 1, 2, or 3\n");
        return 1;
    }

    // Parents shorter than -min_seq_length are skipped at index time.
    // The value is encoded in the file name and consulted by
    // ikafssnsearch / ikafssnclient to validate -min_query_length.
    uint32_t min_seq_length = 64;
    if (cli.has("-min_seq_length")) {
        int v = cli.get_int("-min_seq_length", 0);
        if (v <= 0) {
            std::fprintf(stderr,
                "Error: -min_seq_length must be a positive integer (got %d)\n", v);
            return 1;
        }
        min_seq_length = static_cast<uint32_t>(v);
    }

    // -min_length_split / -overlap_length are coupled: either both 0
    // (splitting disabled) or both non-zero.  When non-zero,
    // min_length_split is in [1000, 1000000] and overlap_length must
    // be < min_length_split / 2.  Defaults: 50000 / 500 for -mode 1,
    // 0 / 0 for -mode 2/3.
    uint32_t min_length_split;
    uint32_t overlap_length;
    if (index_mode == 1) {
        min_length_split = 50000;
        overlap_length = 500;
    } else {
        min_length_split = 0;
        overlap_length = 0;
    }
    if (cli.has("-min_length_split")) {
        long long v = cli.get_int("-min_length_split", -1);
        if (v < 0) {
            std::fprintf(stderr,
                "Error: -min_length_split must be 0 (disabled) or in [1000, 1000000] "
                "(got %lld)\n", v);
            return 1;
        }
        if (v != 0 && (v < 1000 || v > 1000000)) {
            std::fprintf(stderr,
                "Error: -min_length_split must be 0 (disabled) or in [1000, 1000000] "
                "(got %lld)\n", v);
            return 1;
        }
        min_length_split = static_cast<uint32_t>(v);
    }
    if (cli.has("-overlap_length")) {
        long long v = cli.get_int("-overlap_length", -1);
        if (v < 0) {
            std::fprintf(stderr,
                "Error: -overlap_length must be a non-negative integer (got %lld)\n", v);
            return 1;
        }
        overlap_length = static_cast<uint32_t>(v);
    }
    if ((min_length_split == 0) != (overlap_length == 0)) {
        std::fprintf(stderr,
            "Error: -min_length_split (%u) and -overlap_length (%u) must both be 0 "
            "(splitting disabled) or both non-zero (splitting enabled)\n",
            min_length_split, overlap_length);
        return 1;
    }
    if (min_length_split != 0 && overlap_length >= min_length_split / 2) {
        std::fprintf(stderr,
            "Error: -overlap_length (%u) must be < -min_length_split / 2 "
            "(%u / 2 = %u)\n",
            overlap_length, min_length_split, min_length_split / 2);
        return 1;
    }
    // Worst-case tail fragment after splitting must still contain at
    // least min_seq_length fresh bases (with two overlap regions), so
    // the parent-level min_seq_length filter stands in for a
    // fragment-level filter.
    if (min_length_split != 0) {
        const uint64_t required = static_cast<uint64_t>(min_seq_length)
                                + uint64_t{2} * static_cast<uint64_t>(overlap_length);
        if (static_cast<uint64_t>(min_length_split) < required) {
            std::fprintf(stderr,
                "Error: -min_length_split (%u) must be >= -min_seq_length + "
                "2 * -overlap_length (%u + 2 * %u = %llu)\n",
                min_length_split, min_seq_length, overlap_length,
                static_cast<unsigned long long>(required));
            return 1;
        }
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

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::fprintf(stderr, "Error: cannot create output directory '%s': %s\n",
                     out_dir.c_str(), ec.message().c_str());
        return 1;
    }

    std::vector<std::string> vol_paths = BlastDbReader::find_volume_paths(db_path);
    if (vol_paths.empty()) {
        // Single-volume DB: open directly.
        vol_paths.push_back(db_path);
    }

    int threads = resolve_threads(cli);

    int nthread_highfreq_filter;
    if (cli.has("-nthread_highfreq_filter")) {
        nthread_highfreq_filter = cli.get_int("-nthread_highfreq_filter", 8);
        if (nthread_highfreq_filter < 1) {
            std::fprintf(stderr, "Error: -nthread_highfreq_filter must be >= 1\n");
            return 1;
        }
        if (nthread_highfreq_filter > threads) {
            std::fprintf(stderr,
                "Error: -nthread_highfreq_filter (%d) exceeds -nthread (%d)\n",
                nthread_highfreq_filter, threads);
            return 1;
        }
    } else {
        nthread_highfreq_filter = std::min(8, threads);
    }

    int force_rebuild = cli.get_int("-force_rebuild", 0);
    if (force_rebuild != 0 && force_rebuild != 1) {
        std::fprintf(stderr, "Error: -force_rebuild must be 0 or 1\n");
        return 1;
    }

    int max_degen_expand = cli.get_int("-max_degen_expand", 4);
    if (max_degen_expand < 0 || max_degen_expand > 16) {
        std::fprintf(stderr, "Error: -max_degen_expand must be between 0 and 16\n");
        return 1;
    }

    int freq_threshold_part = cli.get_int("-freq_threshold_part", 8);
    if (freq_threshold_part < 0) {
        std::fprintf(stderr, "Error: -freq_threshold_part must be >= 0\n");
        return 1;
    }
    if (freq_threshold_part > 255) {
        std::fprintf(stderr,
            "Error: -freq_threshold_part must be <= 255 "
            "(short-bucket occurrence counts are u8 in the .kpx layout)\n");
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
        logger.info("Parameters: k=%d, t=%d, template_type=%s, mode=%d, "
                    "min_seq_length=%u, min_length_split=%u, overlap_length=%u, "
                    "memory_limit=%s, threads=%d",
                    k, spaced_t, type_display.c_str(),
                    index_mode, min_seq_length, min_length_split, overlap_length,
                    mem_limit_str.c_str(), threads);
    } else {
        logger.info("Parameters: k=%d, mode=%d, min_seq_length=%u, "
                    "min_length_split=%u, overlap_length=%u, "
                    "memory_limit=%s, threads=%d",
                    k, index_mode, min_seq_length,
                    min_length_split, overlap_length,
                    mem_limit_str.c_str(), threads);
    }

    std::string db_base = std::filesystem::path(db_path).filename().string();

    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, threads);

    IndexBuilderConfig config;
    config.k = k;
    config.memory_limit = memory_limit;
    config.threads = threads;
    config.verbose = verbose;
    config.skip_kpx = (index_mode == 1);
    config.max_degen_expand = max_degen_expand;
    config.t = spaced_t;
    config.template_type = static_cast<uint8_t>(spaced_type);
    config.freq_threshold_part = static_cast<uint32_t>(freq_threshold_part);
    config.min_seq_length = min_seq_length;
    config.min_length_split = min_length_split;
    config.overlap_length = overlap_length;
    // When max_freq_build is active (not 1.0 = disabled), keep .tmp files
    // for cross-volume filtering so the filter pass can read them.
    bool freq_filter_active = (max_freq_build != 1.0);
    config.keep_tmp = freq_filter_active;

    uint16_t total_volumes = static_cast<uint16_t>(vol_paths.size());

    std::vector<std::string> vol_basenames(total_volumes);
    for (uint16_t vi = 0; vi < total_volumes; vi++) {
        vol_basenames[vi] = std::filesystem::path(vol_paths[vi]).filename().string();
    }
    {
        std::set<std::string> seen;
        for (const auto& bn : vol_basenames) {
            if (!seen.insert(bn).second) {
                std::fprintf(stderr, "Error: duplicate volume basename '%s'\n", bn.c_str());
                return 1;
            }
        }
    }

    // -template_type=both builds coding then optimal sequentially.
    for (TemplateType cur_type : build_types) {
        uint8_t cur_tt = static_cast<uint8_t>(cur_type);

        if (build_types.size() > 1) {
            logger.info("========== Building %s template ==========",
                        template_type_to_string(cur_type).c_str());
        }

        config.template_type = cur_tt;

        // Tmp prefix omits max_freq_build / max_degen_expand because
        // the absolute threshold is only known after the metadata pass.
        auto tmp_prefix_for = [&](const std::string& vb) {
            return index_file_stem(out_dir, vb, k, spaced_t, cur_tt,
                                   min_seq_length, min_length_split, overlap_length,
                                   /*max_freq_build=*/1,
                                   /*max_degen_expand=*/0);
        };
        std::vector<std::string> vol_prefixes_tmp(total_volumes);
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            vol_prefixes_tmp[vi] = tmp_prefix_for(vol_basenames[vi]);
        }

        std::vector<std::string> vol_prefixes_final(total_volumes);
        std::string khx_path;
        std::string kvx_path;
        auto refresh_final_paths = [&](uint64_t freq_threshold) {
            auto final_prefix_for = [&](const std::string& vb) {
                return index_file_stem(out_dir, vb, k, spaced_t, cur_tt,
                                       min_seq_length, min_length_split, overlap_length,
                                       freq_threshold,
                                       static_cast<uint32_t>(max_degen_expand));
            };
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                vol_prefixes_final[vi] = final_prefix_for(vol_basenames[vi]);
            }
            khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt,
                                    min_seq_length, min_length_split, overlap_length,
                                    freq_threshold,
                                    static_cast<uint32_t>(max_degen_expand));
            kvx_path = index_file_stem(out_dir, db_base, k, spaced_t, cur_tt,
                                       min_seq_length, min_length_split, overlap_length,
                                       freq_threshold,
                                       static_cast<uint32_t>(max_degen_expand))
                      + ".kvx";
        };
        // Provisional value equals the resolved threshold when
        // -max_freq_build is disabled (=1) or an absolute integer (>1).
        // For a fractional value the resolved threshold needs
        // total_fragments and is substituted by the peek pass below.
        uint64_t freq_threshold_provisional =
            freq_filter_active && max_freq_build >= 1.0
                ? static_cast<uint64_t>(max_freq_build)
                : 1;
        refresh_final_paths(freq_threshold_provisional);

        if (freq_filter_active && max_freq_build < 1.0) {
            uint64_t total_fragments_peek = 0;
            bool peek_ok = true;
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                std::string ksx_path;
                std::string ksx_tmp = vol_prefixes_tmp[vi] + ".ksx.tmp";
                if (std::filesystem::exists(ksx_tmp)) {
                    ksx_path = ksx_tmp;
                } else {
                    std::error_code ec;
                    for (auto& entry : std::filesystem::directory_iterator(out_dir, ec)) {
                        if (!entry.is_regular_file()) continue;
                        std::string p = entry.path().string();
                        if (p.size() < 4 || p.substr(p.size() - 4) != ".ksx") continue;
                        IndexFilenameParts parts;
                        if (!parse_index_filename(p, parts)) continue;
                        // parts.has_vol is false for single-volume DBs (the
                        // basename has no numeric suffix), so we match on
                        // vol_basename alone to cover both single and
                        // multi-volume layouts.
                        if (parts.vol_basename != vol_basenames[vi]) continue;
                        if (parts.k != k) continue;
                        if (parts.t != spaced_t) continue;
                        if (parts.template_type != cur_tt) continue;
                        if (parts.min_seq_length != min_seq_length) continue;
                        if (parts.min_length_split != min_length_split) continue;
                        if (parts.overlap_length != overlap_length) continue;
                        if (parts.max_degen_expand !=
                            static_cast<uint32_t>(max_degen_expand)) continue;
                        ksx_path = p;
                        break;
                    }
                }
                if (ksx_path.empty()) {
                    peek_ok = false;
                    break;
                }
                KsxReader reader;
                if (!reader.open(ksx_path)) {
                    peek_ok = false;
                    break;
                }
                total_fragments_peek += reader.num_sequences();
                reader.close();
            }
            if (peek_ok && total_fragments_peek > 0) {
                double resolved = std::ceil(max_freq_build *
                                            static_cast<double>(total_fragments_peek));
                if (resolved < 1.0) resolved = 1.0;
                uint64_t peek_threshold = static_cast<uint64_t>(resolved);
                if (peek_threshold != freq_threshold_provisional) {
                    freq_threshold_provisional = peek_threshold;
                    refresh_final_paths(freq_threshold_provisional);
                    logger.info(
                        "Resolved -max_freq_build=%.6g against existing "
                        "metadata (total_fragments=%lu) -> threshold=%lu",
                        max_freq_build,
                        static_cast<unsigned long>(total_fragments_peek),
                        static_cast<unsigned long>(peek_threshold));
                }
            }
        }

        // -force_rebuild 1: wipe tmp / final files before classification.
        if (force_rebuild == 1) {
            logger.info("-force_rebuild 1: deleting existing tmp / final files");
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                remove_volume_tmp(vol_prefixes_tmp[vi], config.skip_kpx, logger);
                remove_volume_final(vol_prefixes_final[vi], config.skip_kpx, logger);
                remove_if_exists(vol_prefixes_final[vi] + ".kvx", logger);
            }
            if (freq_filter_active) {
                remove_if_exists(khx_path, logger);
            }
        }

        // Classify each volume.  Validation failures emit a warning,
        // delete the offending files, and downgrade to the next state.
        std::vector<VolumeState> state(total_volumes, VolumeState::kNone);
        if (force_rebuild == 0) {
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                const std::string& fp = vol_prefixes_final[vi];
                const std::string& tp = vol_prefixes_tmp[vi];
                const bool has_final =
                    std::filesystem::exists(fp + ".ksx") &&
                    std::filesystem::exists(fp + ".kix") &&
                    (config.skip_kpx || std::filesystem::exists(fp + ".kpx"));
                if (has_final) {
                    if (validate_volume_final_strict(fp, config.skip_kpx, logger)) {
                        state[vi] = VolumeState::kComplete;
                        continue;
                    }
                    logger.warn(
                        "Existing final files for volume %u failed strict "
                        "validation; deleting and rebuilding this volume.",
                        static_cast<unsigned>(vi));
                    remove_volume_final(fp, config.skip_kpx, logger);
                    if (freq_filter_active) {
                        // Stale .khx no longer matches the surviving
                        // volumes; drop it so the filter pass regenerates.
                        remove_if_exists(khx_path, logger);
                    }
                }

                const bool has_tmp_full =
                    std::filesystem::exists(tp + ".ksx.tmp") &&
                    std::filesystem::exists(tp + ".kix.tmp") &&
                    (config.skip_kpx || std::filesystem::exists(tp + ".kpx.tmp"));
                if (has_tmp_full) {
                    if (validate_volume_tmp_strict(tp, config.skip_kpx, logger)) {
                        state[vi] = VolumeState::kPostingsTmp;
                        continue;
                    }
                    logger.warn(
                        "Existing tmp postings files for volume %u failed "
                        "strict validation; deleting and rebuilding this volume.",
                        static_cast<unsigned>(vi));
                    remove_if_exists(tp + ".kix.tmp", logger);
                    if (!config.skip_kpx) remove_if_exists(tp + ".kpx.tmp", logger);
                    // .ksx.tmp may still be salvageable; re-evaluate below.
                }

                const bool has_ksx_tmp = std::filesystem::exists(tp + ".ksx.tmp");
                if (has_ksx_tmp) {
                    if (validate_ksx_file_strict(tp + ".ksx.tmp", logger)) {
                        state[vi] = VolumeState::kMetadataTmp;
                        continue;
                    }
                    logger.warn(
                        "Existing .ksx.tmp for volume %u failed strict "
                        "validation; deleting and rebuilding this volume.",
                        static_cast<unsigned>(vi));
                    remove_if_exists(tp + ".ksx.tmp", logger);
                }

                state[vi] = VolumeState::kNone;
            }
        }

        // Cross-volume frequency counts (and .khx) are only consistent
        // when every volume is rebuilt together: if any volume is not
        // kComplete or .khx fails strict validation, demote every
        // volume to kNone so the filter pass sees identical inputs
        // across volumes.
        KhxState khx_state = KhxState::kNeedsRebuild;
        if (freq_filter_active) {
            bool all_complete = true;
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (state[vi] != VolumeState::kComplete) {
                    all_complete = false;
                    break;
                }
            }
            bool khx_valid = all_complete &&
                validate_khx_file_strict(khx_path, k, spaced_t, cur_tt, logger);
            if (all_complete && khx_valid) {
                khx_state = KhxState::kValid;
            } else {
                if (all_complete && !khx_valid) {
                    logger.warn(
                        ".khx (%s) failed strict validation; demoting all "
                        "volumes to a full rebuild so the cross-volume "
                        "frequency filter sees identical inputs.",
                        khx_path.c_str());
                }
                for (uint16_t vi = 0; vi < total_volumes; vi++) {
                    if (state[vi] == VolumeState::kComplete) {
                        logger.info(
                            "Demoting volume %u to full rebuild "
                            "(-max_freq_build active requires every "
                            "volume to participate in the filter pass).",
                            static_cast<unsigned>(vi));
                        remove_volume_final(vol_prefixes_final[vi],
                                            config.skip_kpx, logger);
                    }
                    state[vi] = VolumeState::kNone;
                }
                remove_if_exists(khx_path, logger);
            }
        }

        // Open BLAST DB volumes once; shared by metadata and postings.
        std::vector<BlastDbReader> dbs(total_volumes);
        bool need_dbs = false;
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            if (state[vi] == VolumeState::kNone ||
                state[vi] == VolumeState::kMetadataTmp) {
                need_dbs = true;
                break;
            }
        }
        if (need_dbs) {
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (!dbs[vi].open(vol_paths[vi])) {
                    std::fprintf(stderr, "Error: cannot open volume '%s'\n",
                                 vol_paths[vi].c_str());
                    return 1;
                }
            }
        }

        // Metadata pass.  kNone: run build_metadata (writes .ksx.tmp).
        // kMetadataTmp: reuse existing .ksx.tmp.  kPostingsTmp /
        // kComplete: fetch num_sequences only, and only when
        // freq_filter_active needs total_fragments.
        std::vector<VolumeMetadata> volume_meta(total_volumes);
        int n_complete_reuse = 0, n_postings_tmp_reuse = 0;
        {
            int n_to_build = 0, n_to_reuse_ksx = 0;
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                switch (state[vi]) {
                    case VolumeState::kNone:        n_to_build++; break;
                    case VolumeState::kMetadataTmp: n_to_reuse_ksx++; break;
                    case VolumeState::kPostingsTmp: n_postings_tmp_reuse++; break;
                    case VolumeState::kComplete:    n_complete_reuse++; break;
                }
            }
            if (n_to_build > 0 || n_to_reuse_ksx > 0) {
                logger.info("=== Collecting metadata for %d of %u volume(s) "
                            "(build=%d, reuse-ksx=%d, threads=%d) ===",
                            n_to_build + n_to_reuse_ksx, total_volumes,
                            n_to_build, n_to_reuse_ksx, threads);
            } else if (n_complete_reuse + n_postings_tmp_reuse > 0) {
                logger.info("Reusing %d complete + %d postings-tmp volume(s); "
                            "no metadata pass needed.",
                            n_complete_reuse, n_postings_tmp_reuse);
            }
        }
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            if (state[vi] == VolumeState::kNone) {
                logger.info("Metadata volume %u/%u: %s",
                            vi + 1, total_volumes, vol_paths[vi].c_str());
                dbs[vi].set_mmap_strategy(BlastDbReader::MMapStrategy::kNormal);
                bool ok = build_metadata(dbs[vi], config, vol_prefixes_tmp[vi],
                                         volume_meta[vi], logger);
                dbs[vi].set_mmap_strategy(BlastDbReader::MMapStrategy::kDontNeed);
                if (!ok) {
                    std::fprintf(stderr,
                        "Error: metadata build failed for volume %u\n",
                        static_cast<unsigned>(vi));
                    return 1;
                }
            } else if (state[vi] == VolumeState::kMetadataTmp) {
                // Reconstruct VolumeMetadata from the existing .ksx.tmp
                // so the postings pass sees the same tables it would
                // have got from build_metadata.
                KsxReader reader;
                if (!reader.open(vol_prefixes_tmp[vi] + ".ksx.tmp")) {
                    std::fprintf(stderr,
                        "Error: cannot re-open %s.ksx.tmp\n",
                        vol_prefixes_tmp[vi].c_str());
                    return 1;
                }
                const uint32_t N = reader.num_sequences();
                volume_meta[vi].num_sequences = N;
                volume_meta[vi].num_parents   = reader.num_parents();
                volume_meta[vi].seq_id_to_blast_oid.resize(N);
                volume_meta[vi].seq_id_to_frag_start.resize(N);
                volume_meta[vi].seq_id_to_frag_end.resize(N);
                for (uint32_t i = 0; i < N; ++i) {
                    volume_meta[vi].seq_id_to_blast_oid[i] =
                        reader.blast_oid(reader.parent_index(i));
                    volume_meta[vi].seq_id_to_frag_start[i] = reader.fragment_start(i);
                    volume_meta[vi].seq_id_to_frag_end[i]   = reader.fragment_end(i);
                }
                reader.close();
            } else if (state[vi] == VolumeState::kPostingsTmp ||
                       state[vi] == VolumeState::kComplete) {
                // The freq-filter recompute below needs num_sequences
                // for total_fragments; the full arrays are unused here.
                if (freq_filter_active) {
                    const std::string ksx = (state[vi] == VolumeState::kPostingsTmp)
                        ? (vol_prefixes_tmp[vi] + ".ksx.tmp")
                        : (vol_prefixes_final[vi] + ".ksx");
                    KsxReader reader;
                    if (reader.open(ksx)) {
                        volume_meta[vi].num_sequences = reader.num_sequences();
                        reader.close();
                    }
                }
            }
        }

        // Resolve max_freq_build against the total fragment count.
        uint64_t total_fragments = 0;
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            total_fragments += volume_meta[vi].num_sequences;
        }
        uint64_t freq_threshold = freq_threshold_provisional;
        if (freq_filter_active) {
            if (max_freq_build < 1.0) {
                double resolved = std::ceil(max_freq_build *
                                            static_cast<double>(total_fragments));
                if (resolved < 1.0) resolved = 1.0;
                freq_threshold = static_cast<uint64_t>(resolved);
                logger.info(
                    "-max_freq_build=%.6g (fraction of total fragments=%lu) "
                    "-> threshold=%lu",
                    max_freq_build,
                    static_cast<unsigned long>(total_fragments),
                    static_cast<unsigned long>(freq_threshold));
            } else {
                freq_threshold = static_cast<uint64_t>(max_freq_build);
            }
            // Refresh the final prefixes if the resolved threshold
            // differs from the provisional one.
            if (freq_threshold != freq_threshold_provisional) {
                refresh_final_paths(freq_threshold);
            }
        }

        // Postings pass for kNone / kMetadataTmp volumes.
        {
            int n_to_build = 0;
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (state[vi] == VolumeState::kNone ||
                    state[vi] == VolumeState::kMetadataTmp) n_to_build++;
            }
            if (n_to_build > 0) {
                logger.info("=== Writing postings for %d volume(s) "
                            "(threads=%d) ===", n_to_build, threads);
            }
        }
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            if (state[vi] != VolumeState::kNone &&
                state[vi] != VolumeState::kMetadataTmp) continue;
            logger.info("Postings volume %u/%u: %s",
                        vi + 1, total_volumes, vol_paths[vi].c_str());
            dbs[vi].set_mmap_strategy(BlastDbReader::MMapStrategy::kNormal);
            bool ok;
            if (kmer_type_for(k, spaced_t) == 0) {
                ok = build_postings<uint16_t>(dbs[vi], config,
                    vol_prefixes_tmp[vi], volume_meta[vi],
                    vi, total_volumes, db_base, logger);
            } else {
                ok = build_postings<uint32_t>(dbs[vi], config,
                    vol_prefixes_tmp[vi], volume_meta[vi],
                    vi, total_volumes, db_base, logger);
            }
            dbs[vi].set_mmap_strategy(BlastDbReader::MMapStrategy::kDontNeed);
            if (!ok) {
                std::fprintf(stderr,
                    "Error: postings build failed for volume %u\n",
                    static_cast<unsigned>(vi));
                return 1;
            }
        }

        // .kvx manifest; re-written every run.
        {
            FILE* fp = std::fopen(kvx_path.c_str(), "w");
            if (!fp) {
                std::fprintf(stderr, "Error: cannot write %s\n", kvx_path.c_str());
                return 1;
            }
            std::fprintf(fp, "#\n# ikafssn index volume manifest\n#\n");
            std::fprintf(fp, "FORMAT_VERSION %u\n", KIX_FORMAT_VERSION);
            std::fprintf(fp, "TITLE %s\n", db_base.c_str());
            std::fprintf(fp, "DBLIST");
            for (const auto& bn : vol_basenames) {
                std::fprintf(fp, " \"%s\"", bn.c_str());
            }
            std::fprintf(fp, "\n");
            std::fclose(fp);
            logger.info("Wrote volume manifest: %s", kvx_path.c_str());
        }

        // Filter (freq_filter_active) or .tmp → final rename.
        if (freq_filter_active) {
            bool all_complete = true;
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (state[vi] != VolumeState::kComplete) {
                    all_complete = false;
                    break;
                }
            }
            if (!all_complete || khx_state != KhxState::kValid) {
                if (!filter_volumes_cross_volume(vol_prefixes_tmp,
                                                 vol_prefixes_final,
                                                 khx_path, k,
                                                 freq_threshold,
                                                 nthread_highfreq_filter,
                                                 logger)) {
                    std::fprintf(stderr, "Error: cross-volume filtering failed\n");
                    return 1;
                }
            }
        } else {
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (state[vi] == VolumeState::kComplete) continue;
                const std::string& tp = vol_prefixes_tmp[vi];
                const std::string& fp = vol_prefixes_final[vi];
                auto rename_one = [&](const char* ext) {
                    std::string from = tp + ext + ".tmp";
                    std::string to   = fp + ext;
                    if (std::rename(from.c_str(), to.c_str()) != 0) {
                        logger.error("Failed to rename %s -> %s",
                                     from.c_str(), to.c_str());
                        return false;
                    }
                    return true;
                };
                if (!rename_one(".ksx")) return 1;
                if (!rename_one(".kix")) return 1;
                if (!config.skip_kpx) {
                    if (!rename_one(".kpx")) return 1;
                }
            }
        }

        // Post-build validation.  On failure, list the offending paths
        // and exit non-zero; never delete files automatically.
        std::vector<std::string> failed_files;
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            const std::string& prefix = vol_prefixes_final[vi];
            const std::string ksx = prefix + ".ksx";
            const std::string kix = prefix + ".kix";
            const std::string kpx = config.skip_kpx ? std::string{} : (prefix + ".kpx");
            if (!validate_ksx_file_strict(ksx, logger)) {
                failed_files.push_back(ksx);
            }
            if (!validate_kix_kpx_strict(kix, kpx, logger)) {
                failed_files.push_back(kix);
                if (!kpx.empty()) failed_files.push_back(kpx);
            } else {
                logger.info("Validated volume: %s", prefix.c_str());
            }
        }
        if (freq_filter_active) {
            if (!validate_khx_file_strict(khx_path, k, spaced_t, cur_tt, logger)) {
                failed_files.push_back(khx_path);
            } else {
                logger.info("Validated .khx: %s", khx_path.c_str());
            }
        }
        if (!failed_files.empty()) {
            std::fprintf(stderr,
                "Error: post-build validation failed for the following file(s):\n");
            for (const auto& p : failed_files) {
                std::fprintf(stderr, "  %s\n", p.c_str());
            }
            std::fprintf(stderr,
                "Delete these files and re-run ikafssnindex to rebuild them.\n");
            return 1;
        }

        if (build_types.size() > 1) {
            logger.info("========== %s template completed ==========",
                        template_type_to_string(cur_type).c_str());
        }
    }

    logger.info("All volumes completed successfully.");
    return 0;
}
