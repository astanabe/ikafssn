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

#include <tbb/global_control.h>

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

// Return on-disk size of a file, or -1 if it cannot be stat'd.
static int64_t file_size_or_neg(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return -1;
    return static_cast<int64_t>(st.st_size);
}

// Verify the on-disk .ksx size matches the v10 two-stage layout:
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

// Validate existing index files for a volume against the BLAST DB.
// Checks (in order): file existence, header validity, BLAST DB
// metadata match (num_sequences / total_bases / k / t / template_type /
// min_seq_length), .kix and .ksx on-disk file size, and (if deep == true)
// a structural walk of the (.kix, .kpx) pair via validate_volume() so
// truncated / FOR-stream-corrupt files left behind by a previous crash
// are caught before we try to reuse them.
//
// suffix: "" for final files (.kix), ".tmp" for temp files (.kix.tmp).
static IndexValidation validate_existing_index(
    const std::string& output_prefix,
    const std::string& vol_path,
    int k, uint8_t t, uint8_t template_type,
    uint32_t min_seq_length,
    uint32_t min_length_split,
    uint32_t overlap_length,
    bool skip_kpx,
    const Logger& logger,
    bool deep,
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

    if (std::memcmp(kix_hdr.magic, KIX_MAGIC, sizeof(KIX_MAGIC)) != 0) return result;
    if (kix_hdr.format_version != KIX_FORMAT_VERSION) return result;

    // Check .ksx exists, read header + parent_lengths to compute total_bases
    // (in v10 the parent table — not the fragment table — carries the
    // full per-OID lengths that should sum to the BLAST DB total).
    fp = std::fopen(ksx_path.c_str(), "rb");
    if (!fp) return result; // kNotFound
    KsxHeader ksx_hdr{};
    ok = (std::fread(&ksx_hdr, sizeof(ksx_hdr), 1, fp) == 1);
    if (!ok) { std::fclose(fp); return result; }

    if (std::memcmp(ksx_hdr.magic, KSX_MAGIC, sizeof(KSX_MAGIC)) != 0) { std::fclose(fp); return result; }
    if (ksx_hdr.format_version != KSX_FORMAT_VERSION) { std::fclose(fp); return result; }

    uint32_t ksx_nseq     = ksx_hdr.num_sequences;
    uint32_t ksx_nparents = ksx_hdr.num_parents;
    std::vector<uint32_t> parent_lengths(ksx_nparents);
    if (ksx_nparents > 0) {
        size_t read_count = std::fread(parent_lengths.data(), sizeof(uint32_t), ksx_nparents, fp);
        if (read_count != ksx_nparents) { std::fclose(fp); return result; }
    }
    std::fclose(fp);

    uint64_t ksx_total_bases = 0;
    for (uint32_t i = 0; i < ksx_nparents; i++) {
        ksx_total_bases += parent_lengths[i];
    }

    // Check .kpx exists (if required)
    std::string kpx_path = skip_kpx ? std::string{} : (output_prefix + ".kpx" + suffix);
    if (!skip_kpx) {
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

    // Check min_seq_length matches what we are about to build with.  The
    // .ksx parent table was filtered by this value, so reusing the index
    // requires the same filter.
    if (kix_hdr.min_seq_length != min_seq_length ||
        ksx_hdr.min_seq_length != min_seq_length) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: min_seq_length (kix=%u, ksx=%u) does not match requested %u",
            kix_path.c_str(),
            kix_hdr.min_seq_length, ksx_hdr.min_seq_length,
            min_seq_length);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    // Check min_length_split / overlap_length match what we are about to
    // build with.  The .ksx fragment table is keyed off these values.
    if (kix_hdr.min_length_split != min_length_split ||
        ksx_hdr.min_length_split != min_length_split ||
        kix_hdr.overlap_length   != overlap_length   ||
        ksx_hdr.overlap_length   != overlap_length) {
        char buf[320];
        std::snprintf(buf, sizeof(buf),
            "%s: split parameters (kix min_length_split=%u, overlap_length=%u; "
            "ksx min_length_split=%u, overlap_length=%u) do not match requested "
            "(min_length_split=%u, overlap_length=%u)",
            kix_path.c_str(),
            kix_hdr.min_length_split, kix_hdr.overlap_length,
            ksx_hdr.min_length_split, ksx_hdr.overlap_length,
            min_length_split, overlap_length);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    // Recompute BLAST DB stats with the same min_seq_length filter so we
    // can compare them against the .ksx-recorded counts.  Sequences whose
    // length is below min_seq_length are excluded from both num_seqs and
    // total_bases.
    uint32_t db_filtered_nseq = 0;
    uint64_t db_filtered_total_bases = 0;
    for (uint32_t oid = 0; oid < db_nseq; oid++) {
        uint32_t slen = db.seq_length(oid);
        if (slen < min_seq_length) continue;
        db_filtered_nseq++;
        db_filtered_total_bases += slen;
    }

    // Check num_parents (= number of OIDs that survived the
    // min_seq_length filter) against the BLAST DB.  Fragment count
    // (num_sequences) is not directly comparable when splitting is on,
    // since it depends on N-runs in each parent's sequence; we only
    // require num_sequences == kix_hdr.num_sequences (internal
    // consistency between the .kix and .ksx files).
    if (ksx_nparents != db_filtered_nseq) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: num_parents=%u but BLAST DB '%s' has %u (after min_seq_length=%u filter)",
            ksx_path.c_str(), ksx_nparents,
            vol_path.c_str(), db_filtered_nseq, min_seq_length);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }
    if (kix_hdr.num_sequences != ksx_nseq) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: num_sequences=%u but %s reports %u",
            kix_path.c_str(), kix_hdr.num_sequences,
            ksx_path.c_str(), ksx_nseq);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    // Check total_bases (sum of parent lengths).
    if (ksx_total_bases != db_filtered_total_bases) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "%s: total_bases=%lu but BLAST DB '%s' has %lu (after min_seq_length=%u filter)",
            ksx_path.c_str(), static_cast<unsigned long>(ksx_total_bases),
            vol_path.c_str(), static_cast<unsigned long>(db_filtered_total_bases),
            min_seq_length);
        result.status = IndexStatus::kMismatch;
        result.detail = buf;
        return result;
    }

    // File size sanity (catches header-only-correct truncation that the
    // BLAST-DB-metadata check above can't see).  Treated as kNotFound so
    // the caller falls back to per-volume rebuild instead of erroring out.
    if (!verify_kix_file_size(kix_path, logger)) return result;
    if (!verify_ksx_file_size(ksx_path, logger)) return result;

    // Deep structural walk of the (.kix, .kpx) pair.  validate_volume
    // walks every k-mer's posting list and verifies the byte length
    // recorded in the .kpx EF dictionary matches the bytes actually
    // consumed by the kind map + partition groups + short FOR streams,
    // catching truncation / FOR-stream corruption that the size check
    // above can't see for .kpx.
    if (deep) {
        if (!validate_volume(kix_path, kpx_path, nullptr, logger)) {
            return result; // kNotFound
        }
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
        "                         .kpx v8 per-(kmer, seq_id) partition threshold (default: 8, max: 255)\n"
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
        "                         (Phase 7d default-on validation walks each\n"
        "                         k-mer's .kpx posting list and checks the\n"
        "                         byte length against the EF dictionary)\n"
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

    // Parse -mode (1, 2, or 3; default 1).
    // Phase 2 (fragment indexing): -mode 1 is now the recommended default
    // because fragment splitting brings .kpx-less search into the practical
    // range for nt-class databases.
    int index_mode = cli.get_int("-mode", 1);
    if (index_mode < 1 || index_mode > 3) {
        std::fprintf(stderr, "Error: -mode must be 1, 2, or 3\n");
        return 1;
    }

    // Parse -min_seq_length (default 64).  Sequences shorter than this
    // are skipped at index time; the value is persisted in the .kix /
    // .kpx / .ksx headers and consulted by ikafssnsearch / ikafssnclient
    // to validate -min_query_length.
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

    // Parse -min_length_split / -overlap_length (Phase 2).
    //
    // Defaults depend on -mode:
    //   - -mode 1: min_length_split = 50000, overlap_length = 500
    //   - -mode 2/3: min_length_split = 0, overlap_length = 0 (no splitting)
    //
    // The two flags are coupled: either both 0 (splitting disabled) or both
    // non-zero.  When non-zero, min_length_split must be in [1000, 1000000]
    // and overlap_length must be < min_length_split / 2.
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
    // Coupling: both 0 or both non-zero.
    if ((min_length_split == 0) != (overlap_length == 0)) {
        std::fprintf(stderr,
            "Error: -min_length_split (%u) and -overlap_length (%u) must both be 0 "
            "(splitting disabled) or both non-zero (splitting enabled)\n",
            min_length_split, overlap_length);
        return 1;
    }
    // overlap_length must leave at least min_length_split / 2 of fresh bases
    // per fragment, otherwise the kafsss split formula degenerates.
    if (min_length_split != 0 && overlap_length >= min_length_split / 2) {
        std::fprintf(stderr,
            "Error: -overlap_length (%u) must be < -min_length_split / 2 "
            "(%u / 2 = %u)\n",
            overlap_length, min_length_split, min_length_split / 2);
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

    // Phase 7d: post-build structural validation is on by default; -no_validate
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
            "(short-bucket occurrence counts are u8 in the v8 .kpx layout)\n");
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

    // Extract DB base name from path
    std::string db_base = std::filesystem::path(db_path).filename().string();

    // Centralized TBB thread control
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, threads);

    // Build config
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
            // Validate each volume's existing index (deep walk + size checks
            // included so truncated / FOR-stream-corrupt files left behind by
            // a previous crash are caught before we try to reuse them).
            std::vector<IndexValidation> validations(total_volumes);
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                validations[vi] = validate_existing_index(vol_prefixes[vi], vol_paths[vi],
                    k, spaced_t, cur_tt, min_seq_length,
                    min_length_split, overlap_length,
                    config.skip_kpx, logger, /*deep=*/true);
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
                // With cross-volume filtering: skip the build pipeline only
                // when ALL volumes have valid final files AND the shared
                // .khx is valid (entire pipeline completed).  Otherwise
                // fall back to per-volume .tmp resume.
                //
                // Note: .khx is now written before filter_one_volume() runs
                // (see filter_volumes_cross_volume), so a crash mid-pipeline
                // either leaves .tmp files behind (recoverable) or, after
                // all volumes are filtered, has .khx already on disk.  The
                // "all final + .khx missing" combination therefore only
                // arises from manual deletion of .khx; in that case the
                // resume path falls through here, finds no .tmp, and triggers
                // a full rebuild — which is the only option since the
                // pre-filter distinct-seq_id counts cannot be recovered
                // from the post-filter .kix files.
                bool all_final_valid = true;
                for (uint16_t vi = 0; vi < total_volumes; vi++) {
                    if (validations[vi].status != IndexStatus::kValid) {
                        all_final_valid = false;
                        break;
                    }
                }
                std::string khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt);
                bool khx_valid = false;
                if (all_final_valid) {
                    auto khx_val = validate_khx_standalone(khx_path, k, spaced_t, cur_tt, logger);
                    khx_valid = (khx_val.status == IndexStatus::kValid);
                }
                if (all_final_valid && khx_valid) {
                    logger.info("All %d volumes have valid indexes and .khx is valid; "
                                "skipping build and filter (-rebuild 0)",
                                total_volumes);
                    build_skipped = true;
                } else {
                    if (all_final_valid && !khx_valid) {
                        logger.warn("All %d volumes valid but .khx missing or invalid; "
                                    "cannot regenerate .khx from final files; "
                                    "rebuilding from BLAST DB",
                                    total_volumes);
                    }
                    // Per-volume .tmp resume.  Filter needs .tmp from ALL
                    // volumes, but we can reuse any valid .tmp and only
                    // rebuild the missing/invalid ones.
                    for (uint16_t vi = 0; vi < total_volumes; vi++) {
                        auto tmp_val = validate_existing_index(vol_prefixes[vi], vol_paths[vi],
                            k, spaced_t, cur_tt, min_seq_length,
                            min_length_split, overlap_length,
                            config.skip_kpx, logger, /*deep=*/true, ".tmp");
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
            // Process volumes sequentially.  Volume-level concurrency was
            // dropped along with the legacy -openvol option: per-volume
            // partition processing is already TBB-parallel end-to-end
            // (counting reduce, sort, per-k-mer encode), so an additional
            // outer task_group provided no extra parallelism while
            // doubling peak memory and forcing -memory_limit / openvol
            // accounting on every caller.
            for (uint16_t vi = 0; vi < total_volumes; vi++) {
                if (skip_volume[vi]) continue;

                logger.info("=== Volume %d/%d: %s ===", vi + 1, total_volumes,
                            vol_paths[vi].c_str());

                BlastDbReader db;
                if (!db.open(vol_paths[vi])) {
                    std::fprintf(stderr,
                        "Error: cannot open volume '%s'\n",
                        vol_paths[vi].c_str());
                    return 1;
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
                    std::fprintf(stderr,
                        "Error: index build failed for volume %u\n",
                        static_cast<unsigned>(vi));
                    return 1;
                }
            }
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

        // Post-build cross-volume frequency filtering for this template type
        if (freq_filter_active && !build_skipped) {
            std::string khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt);

            if (!filter_volumes_cross_volume(vol_prefixes, khx_path, k,
                                             freq_threshold, nthread_highfreq_filter,
                                             logger)) {
                std::fprintf(stderr, "Error: cross-volume filtering failed\n");
                return 1;
            }
        }

        // Phase 7d: structural validation on the just-finalised .kix / .kpx
        // pair (catches silent kind-map / FOR-stream corruption that the
        // v9 dedup'd headers can no longer detect by redundancy), plus
        // .kix / .ksx / .khx file size checks (catch truncation that the
        // structural walk and metadata checks can't see).
        if (run_validate && !build_skipped) {
            for (size_t vi = 0; vi < vol_prefixes.size(); vi++) {
                const std::string& prefix = vol_prefixes[vi];
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
                std::string khx_path = khx_path_for(out_dir, db_base, k, spaced_t, cur_tt);
                auto khx_val = validate_khx_standalone(khx_path, k, spaced_t, cur_tt, logger);
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
