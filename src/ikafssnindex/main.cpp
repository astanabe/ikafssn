#include "io/blastdb_reader.hpp"
#include "io/volume_discovery.hpp"
#include "index/index_builder.hpp"
#include "index/index_filter.hpp"
#include "index/kix_format.hpp"
#include "index/kix_reader.hpp"
#include "index/ksx_format.hpp"
#include "index/khx_format.hpp"
#include "index/khx_writer.hpp"
#include "core/config.hpp"
#include "core/spaced_seed.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "index/parallel_sort_dispatch.hpp"
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
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <filesystem>

#include <atomic>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

using namespace ikafssn;

// Per-volume memory-cost estimate fed to plan_volume_batches().
struct VolumeCost {
    uint16_t volume_index;
    uint64_t cost;
};

// Concurrency batch: every listed volume runs in parallel; batches are
// processed sequentially.  `oversized` flags single-volume batches
// whose lone volume exceeds memory_budget — those receive the full
// budget and rely on per-volume partitioning to fit.
struct VolumeBatch {
    std::vector<uint16_t> vols;
    uint64_t total_cost = 0;
    bool oversized = false;
};

// Greedily pack volumes into batches whose summed cost stays within
// memory_budget.
static std::vector<VolumeBatch>
plan_volume_batches(const std::vector<VolumeCost>& costs,
                    uint64_t memory_budget) {
    std::vector<VolumeBatch> batches;
    VolumeBatch cur;

    auto flush = [&]() {
        if (!cur.vols.empty()) {
            batches.push_back(std::move(cur));
            cur = VolumeBatch{};
        }
    };

    for (const auto& vc : costs) {
        const uint64_t c = std::max<uint64_t>(vc.cost, 1);
        if (c > memory_budget) {
            flush();
            VolumeBatch single;
            single.vols.push_back(vc.volume_index);
            single.total_cost = c;
            single.oversized = true;
            batches.push_back(std::move(single));
            continue;
        }
        if (cur.total_cost + c > memory_budget && !cur.vols.empty()) {
            flush();
        }
        cur.vols.push_back(vc.volume_index);
        cur.total_cost += c;
    }
    flush();
    return batches;
}

// Result of validating an existing .khx file.
enum class IndexStatus {
    kNotFound,   // .khx is missing or fails the size / header checks
    kValid,      // .khx is present and structurally consistent
};

struct IndexValidation {
    IndexStatus status = IndexStatus::kNotFound;
};

// Return on-disk size of a file, or -1 if it cannot be stat'd.
static int64_t file_size_or_neg(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return -1;
    return static_cast<int64_t>(st.st_size);
}

// Verify the on-disk .ksx size matches the two-stage layout:
//   header + parent_lengths[P] + parent_blast_oids[P] + parent_acc_offsets[P+1]
//          + accession-string-table[acc_bytes]
//          + fragment_parent_idx[N] + fragment_start[N] + fragment_end[N]
// where P = num_parents, N = num_sequences, and acc_bytes is recorded in
// parent_acc_offsets[P].
static bool verify_ksx_file_size(const std::string& ksx_path,
                                 const Logger& logger) {
    FILE* fp = std::fopen(ksx_path.c_str(), "rb");
    if (!fp) return false;
    KsxHeader hdr{};
    if (std::fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        std::fclose(fp);
        return false;
    }
    if (std::memcmp(hdr.magic, KSX_MAGIC, sizeof(KSX_MAGIC)) != 0
        || hdr.format_version != KSX_FORMAT_VERSION) {
        std::fclose(fp);
        return false;
    }
    const uint32_t num_sequences = hdr.num_sequences;
    const uint32_t num_parents   = hdr.num_parents;

    // Read parent_acc_offsets[num_parents] (the sentinel = total acc bytes).
    const long acc_offsets_pos = static_cast<long>(sizeof(KsxHeader))
                               + static_cast<long>(sizeof(uint32_t)) * num_parents
                               + static_cast<long>(sizeof(uint32_t)) * num_parents
                               + static_cast<long>(sizeof(uint32_t)) * num_parents;
    if (std::fseek(fp, acc_offsets_pos, SEEK_SET) != 0) {
        std::fclose(fp);
        return false;
    }
    uint32_t string_table_bytes = 0;
    if (std::fread(&string_table_bytes, sizeof(uint32_t), 1, fp) != 1) {
        std::fclose(fp);
        return false;
    }
    std::fclose(fp);

    int64_t expected = static_cast<int64_t>(sizeof(KsxHeader))
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_parents          // parent_lengths
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_parents          // parent_blast_oids
                     + static_cast<int64_t>(sizeof(uint32_t)) * (num_parents + 1)    // parent_acc_offsets
                     + static_cast<int64_t>(string_table_bytes)                       // accession strings
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_sequences        // fragment_parent_idx
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_sequences        // fragment_start
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_sequences;       // fragment_end
    int64_t actual = file_size_or_neg(ksx_path);
    if (actual != expected) {
        logger.warn("validate: %s file size %ld != expected %ld",
                    ksx_path.c_str(),
                    static_cast<long>(actual),
                    static_cast<long>(expected));
        return false;
    }
    return true;
}

// Open the .kix via KixReader and verify that its on-disk posting file
// region exactly matches the EF dictionary's sentinel offset (the
// expected end of the posting file).  Catches truncation past the
// dictionary as well as oversized files.  Returns true on match.
static bool verify_kix_file_size(const std::string& kix_path,
                                 const Logger& logger) {
    KixReader kix;
    if (!kix.open(kix_path)) return false;
    uint64_t expected_post = kix.posting_list_offset(kix.table_size());
    uint64_t actual_post = kix.posting_file_size();
    if (expected_post != actual_post) {
        logger.warn("validate: %s posting file size %lu != EF sentinel %lu",
                    kix_path.c_str(),
                    static_cast<unsigned long>(actual_post),
                    static_cast<unsigned long>(expected_post));
        return false;
    }
    return true;
}

// Standalone validation for a .khx file.  Checks magic, format_version,
// the (k, t, template_type) tuple matches what we are about to build,
// and the on-disk file size equals header + ceil(4^k / 8).
// Returns kValid if usable, kNotFound otherwise (callers treat the
// .khx as missing and trigger a regeneration / rebuild path).
static IndexValidation validate_khx_standalone(
    const std::string& khx_path,
    int k, uint8_t t, uint8_t template_type,
    const Logger& logger) {

    IndexValidation result;
    FILE* fp = std::fopen(khx_path.c_str(), "rb");
    if (!fp) return result; // kNotFound

    KhxHeader hdr{};
    bool ok = (std::fread(&hdr, sizeof(hdr), 1, fp) == 1);
    std::fclose(fp);
    if (!ok) return result;

    if (std::memcmp(hdr.magic, KHX_MAGIC, 4) != 0) {
        logger.warn("validate: %s has invalid magic", khx_path.c_str());
        return result;
    }
    if (hdr.format_version != KHX_FORMAT_VERSION) {
        logger.warn("validate: %s format_version=%u != expected %u",
                    khx_path.c_str(), hdr.format_version, KHX_FORMAT_VERSION);
        return result;
    }
    if (hdr.k != static_cast<uint8_t>(k) ||
        hdr.t != t ||
        hdr.template_type != template_type) {
        logger.warn("validate: %s (k=%u, t=%u, template_type=%u) does not match "
                    "requested (k=%d, t=%u, template_type=%u)",
                    khx_path.c_str(), hdr.k, hdr.t, hdr.template_type,
                    k, t, template_type);
        return result;
    }

    int64_t expected = static_cast<int64_t>(sizeof(KhxHeader))
                     + static_cast<int64_t>((table_size(k) + 7) / 8);
    int64_t actual = file_size_or_neg(khx_path);
    if (actual != expected) {
        logger.warn("validate: %s file size %ld != expected %ld",
                    khx_path.c_str(),
                    static_cast<long>(actual),
                    static_cast<long>(expected));
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
        "  -no_validate           Skip the post-build structural validation pass\n"
        "                         (default-on; walks each k-mer's .kpx posting\n"
        "                         list and checks its byte length against the\n"
        "                         EF dictionary)\n"
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

    // Post-build structural validation is on by default; -no_validate
    // opts out (e.g. when building a known-good index for benchmarking).
    const bool run_validate = !cli.has("-no_validate");

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

        // Per-volume "tmp" prefixes omit max_freq_build / max_degen_expand
        // because the absolute threshold is only known after the metadata
        // pass.  The final prefix encodes the resolved values and is used
        // for the rename / filter output.
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

        // Metadata pass: collect parent / fragment tables and write the
        // .ksx.tmp for every volume.  BlastDbReader is mmap-based and
        // CSeqDB is thread-safe, so each reader is opened once here and
        // re-used by the postings pass.
        std::vector<VolumeMetadata> volume_meta(total_volumes);
        std::vector<BlastDbReader> dbs(total_volumes);
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            if (!dbs[vi].open(vol_paths[vi])) {
                std::fprintf(stderr, "Error: cannot open volume '%s'\n",
                             vol_paths[vi].c_str());
                return 1;
            }
        }

        if (rebuild_mode == 1) {
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                cleanup_tmp_files(vol_prefixes_tmp[vi], config.skip_kpx, logger);
            }
        }

        // ~256 B per OID covers the per-OID accession + frags vector
        // with headroom for nt-class deflines.
        constexpr uint64_t META_BYTES_PER_OID = 256;
        std::vector<VolumeCost> meta_costs(total_volumes);
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            meta_costs[vi] = {vi,
                static_cast<uint64_t>(dbs[vi].num_sequences()) * META_BYTES_PER_OID};
        }
        auto meta_batches = plan_volume_batches(meta_costs, memory_limit);
        logger.info("=== Collecting metadata for %u volume(s) "
                    "(threads=%d, %zu batch(es)) ===",
                    total_volumes, threads, meta_batches.size());

        std::atomic<bool> meta_failed{false};
        std::atomic<uint16_t> meta_failed_vi{0};
        for (size_t bi = 0; bi < meta_batches.size(); bi++) {
            if (meta_failed.load()) break;
            const auto& batch = meta_batches[bi];
            const int batch_concurrency = std::min(
                static_cast<int>(batch.vols.size()), std::max(1, threads));
            tbb::task_arena meta_arena(batch_concurrency);
            meta_arena.execute([&] {
                tbb::parallel_for(size_t(0), batch.vols.size(),
                    [&](size_t bvi) {
                        if (meta_failed.load(std::memory_order_relaxed)) return;
                        uint16_t vi = batch.vols[bvi];
                        logger.info("Metadata batch %zu/%zu, volume %u: %s",
                                    bi + 1, meta_batches.size(),
                                    vi + 1, vol_paths[vi].c_str());
                        if (!build_metadata(dbs[vi], config, vol_prefixes_tmp[vi],
                                            volume_meta[vi], logger)) {
                            uint16_t expected = 0;
                            meta_failed_vi.compare_exchange_strong(expected, vi);
                            meta_failed.store(true, std::memory_order_relaxed);
                        }
                    });
            });
        }
        if (meta_failed.load()) {
            std::fprintf(stderr,
                "Error: metadata build failed for volume %u\n",
                static_cast<unsigned>(meta_failed_vi.load()));
            return 1;
        }

        uint64_t total_fragments = 0;
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            total_fragments += volume_meta[vi].num_sequences;
        }

        // Resolve max_freq_build against the total fragment count.
        uint64_t freq_threshold = 1;
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
        }

        auto final_prefix_for = [&](const std::string& vb) {
            return index_file_stem(out_dir, vb, k, spaced_t, cur_tt,
                                   min_seq_length, min_length_split, overlap_length,
                                   freq_threshold,
                                   static_cast<uint32_t>(max_degen_expand));
        };
        std::vector<std::string> vol_prefixes_final(total_volumes);
        for (uint16_t vi = 0; vi < total_volumes; vi++) {
            vol_prefixes_final[vi] = final_prefix_for(vol_basenames[vi]);
        }
        std::string khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt,
                                            min_seq_length, min_length_split, overlap_length,
                                            freq_threshold,
                                            static_cast<uint32_t>(max_degen_expand));
        std::string kvx_path = index_file_stem(out_dir, db_base, k, spaced_t, cur_tt,
                                                min_seq_length, min_length_split, overlap_length,
                                                freq_threshold,
                                                static_cast<uint32_t>(max_degen_expand))
                              + ".kvx";

        // Skip the postings / filter / rename steps when -rebuild=0 and
        // all final files already exist.  The .ksx.tmp from the metadata
        // pass is harmless (it is left for the next run's resume path).
        bool build_skipped = false;
        if (rebuild_mode == 0) {
            bool all_final_present = true;
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                std::string kix = vol_prefixes_final[vi] + ".kix";
                std::string ksx = vol_prefixes_final[vi] + ".ksx";
                std::string kpx = vol_prefixes_final[vi] + ".kpx";
                if (!std::filesystem::exists(kix) ||
                    !std::filesystem::exists(ksx) ||
                    (!config.skip_kpx && !std::filesystem::exists(kpx))) {
                    all_final_present = false;
                    break;
                }
            }
            bool khx_present = !freq_filter_active ||
                std::filesystem::exists(khx_path);
            if (all_final_present && khx_present) {
                logger.info("All %d volumes have valid final indexes; "
                            "skipping rebuild (-rebuild 0)",
                            total_volumes);
                build_skipped = true;

                // Deep validation walk so corruption is detected before
                // we report success.
                if (run_validate) {
                    for (size_t vi = 0; vi < total_volumes; vi++) {
                        const std::string& prefix = vol_prefixes_final[vi];
                        std::string kix_path = prefix + ".kix";
                        std::string ksx_path = prefix + ".ksx";
                        std::string kpx_path = config.skip_kpx
                            ? std::string{} : (prefix + ".kpx");
                        if (!verify_kix_file_size(kix_path, logger) ||
                            !verify_ksx_file_size(ksx_path, logger) ||
                            !validate_volume(kix_path, kpx_path, nullptr, logger)) {
                            std::fprintf(stderr,
                                "Error: existing index validation failed for %s; "
                                "delete and re-run with -rebuild 1\n",
                                prefix.c_str());
                            return 1;
                        }
                    }
                    if (freq_filter_active) {
                        auto khx_val = validate_khx_standalone(
                            khx_path, k, spaced_t, cur_tt, logger);
                        if (khx_val.status != IndexStatus::kValid) {
                            std::fprintf(stderr,
                                "Error: existing .khx validation failed for %s; "
                                "delete and re-run with -rebuild 1\n",
                                khx_path.c_str());
                            return 1;
                        }
                    }
                }
            }
        }

        // Postings pass: per-volume cost is the worst-case TempEntry
        // buffer (one entry per base × overhead).  Each batch divides
        // memory_limit equally among its concurrent volumes;
        // build_postings partitions to fit its own share.
        if (!build_skipped) {
            const uint64_t entry_overhead = parallel_sort_entry_overhead();
            std::vector<VolumeCost> post_costs(total_volumes);
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                post_costs[vi] = {vi, dbs[vi].total_length() * entry_overhead};
            }
            auto post_batches = plan_volume_batches(post_costs, memory_limit);
            logger.info("=== Writing postings for %u volume(s) "
                        "(threads=%d, %zu batch(es)) ===",
                        total_volumes, threads, post_batches.size());

            std::atomic<bool> post_failed{false};
            std::atomic<uint16_t> post_failed_vi{0};
            for (size_t bi = 0; bi < post_batches.size(); bi++) {
                if (post_failed.load()) break;
                const auto& batch = post_batches[bi];
                const int batch_concurrency = std::min(
                    static_cast<int>(batch.vols.size()), std::max(1, threads));
                IndexBuilderConfig batch_config = config;
                if (batch.oversized || batch_concurrency <= 1) {
                    batch_config.memory_limit = memory_limit;
                } else {
                    batch_config.memory_limit = memory_limit /
                        static_cast<uint64_t>(batch_concurrency);
                }
                tbb::task_arena post_arena(batch_concurrency);
                post_arena.execute([&] {
                    tbb::parallel_for(size_t(0), batch.vols.size(),
                        [&](size_t bvi) {
                            if (post_failed.load(std::memory_order_relaxed)) return;
                            uint16_t vi = batch.vols[bvi];
                            logger.info("Postings batch %zu/%zu, volume %u: %s",
                                        bi + 1, post_batches.size(),
                                        vi + 1, vol_paths[vi].c_str());
                            bool ok;
                            if (kmer_type_for(k, spaced_t) == 0) {
                                ok = build_postings<uint16_t>(dbs[vi], batch_config,
                                    vol_prefixes_tmp[vi], volume_meta[vi],
                                    vi, total_volumes, db_base, logger);
                            } else {
                                ok = build_postings<uint32_t>(dbs[vi], batch_config,
                                    vol_prefixes_tmp[vi], volume_meta[vi],
                                    vi, total_volumes, db_base, logger);
                            }
                            if (!ok) {
                                uint16_t expected = 0;
                                post_failed_vi.compare_exchange_strong(expected, vi);
                                post_failed.store(true, std::memory_order_relaxed);
                            }
                        });
                });
            }
            if (post_failed.load()) {
                std::fprintf(stderr,
                    "Error: postings build failed for volume %u\n",
                    static_cast<unsigned>(post_failed_vi.load()));
                return 1;
            }
        }

        // Write .kvx manifest for this template type.
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

        // Cross-volume frequency filtering: consumes .tmp, emits final
        // .kix / .kpx / .ksx and the shared .khx.
        if (!build_skipped && freq_filter_active) {
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

        // Rename .tmp -> final when no cross-volume filter pass ran.
        if (!build_skipped && !freq_filter_active) {
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
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

        // Structural validation on the finalised .kix / .kpx pair, plus
        // .kix / .ksx / .khx file size checks.
        if (run_validate && !build_skipped) {
            for (size_t vi = 0; vi < vol_prefixes_final.size(); vi++) {
                const std::string& prefix = vol_prefixes_final[vi];
                std::string kix_path = prefix + ".kix";
                std::string ksx_path = prefix + ".ksx";
                std::string kpx_path = config.skip_kpx ? std::string{} : (prefix + ".kpx");
                if (!verify_kix_file_size(kix_path, logger)) {
                    std::fprintf(stderr,
                        "Error: post-build .kix file size check failed for %s\n",
                        kix_path.c_str());
                    return 1;
                }
                if (!verify_ksx_file_size(ksx_path, logger)) {
                    std::fprintf(stderr,
                        "Error: post-build .ksx file size check failed for %s\n",
                        ksx_path.c_str());
                    return 1;
                }
                if (!validate_volume(kix_path, kpx_path, nullptr, logger)) {
                    std::fprintf(stderr,
                        "Error: post-build validation failed for volume %s\n",
                        prefix.c_str());
                    return 1;
                }
                logger.info("Validated volume: %s", prefix.c_str());
            }
            if (freq_filter_active) {
                auto khx_val = validate_khx_standalone(khx_path, k, spaced_t,
                                                       cur_tt, logger);
                if (khx_val.status != IndexStatus::kValid) {
                    std::fprintf(stderr,
                        "Error: post-build .khx validation failed for %s\n",
                        khx_path.c_str());
                    return 1;
                }
                logger.info("Validated .khx: %s", khx_path.c_str());
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
