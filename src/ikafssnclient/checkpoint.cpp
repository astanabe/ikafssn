#include "ikafssnclient/checkpoint.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>

#include <openssl/evp.h>

#include "io/sam_writer.hpp"

namespace ikafssn {

// ---------------------------------------------------------------------------
// SHA256 utilities
// ---------------------------------------------------------------------------

static std::string bytes_to_hex(const unsigned char* data, unsigned int len) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < len; i++)
        oss << std::setw(2) << static_cast<unsigned>(data[i]);
    return oss.str();
}

std::string sha256_string(const std::string& data) {
    unsigned char hash[EVP_MAX_MD_SIZE];
    unsigned int hash_len = 0;
    EVP_Digest(data.data(), data.size(), hash, &hash_len,
               EVP_sha256(), nullptr);
    return bytes_to_hex(hash, hash_len);
}

std::string sha256_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) return "";

    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr);

    char buf[8192];
    while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
        EVP_DigestUpdate(ctx, buf, static_cast<size_t>(in.gcount()));
    }

    unsigned char hash[EVP_MAX_MD_SIZE];
    unsigned int hash_len = 0;
    EVP_DigestFinal_ex(ctx, hash, &hash_len);
    EVP_MD_CTX_free(ctx);

    return bytes_to_hex(hash, hash_len);
}

// ---------------------------------------------------------------------------
// DbStats
// ---------------------------------------------------------------------------

DbStats resolve_db_stats(const InfoResponse& info, const std::string& db_name,
                          uint8_t k) {
    DbStats stats;
    stats.db_name = db_name;
    for (const auto& db : info.databases) {
        if (db.name == db_name) {
            for (const auto& grp : db.groups) {
                if (grp.k == k) {
                    for (const auto& vol : grp.volumes) {
                        stats.total_sequences += vol.num_sequences;
                        stats.total_bases += vol.total_bases;
                    }
                    break;
                }
            }
            break;
        }
    }
    return stats;
}

// ---------------------------------------------------------------------------
// build_options_text
// ---------------------------------------------------------------------------

std::string build_options_text(const SearchRequest& req, const DbStats& stats,
                                uint8_t resolved_k, OutputFormat outfmt,
                                const std::string& seqidlist_sha256,
                                const std::string& neg_seqidlist_sha256) {
    std::ostringstream oss;
    oss << "k=" << static_cast<int>(req.k) << "\n";
    oss << "resolved_k=" << static_cast<int>(resolved_k) << "\n";
    oss << "mode=" << static_cast<int>(req.mode) << "\n";
    oss << "stage1_score=" << static_cast<int>(req.stage1_score) << "\n";
    oss << "stage1_topn=" << req.stage1_topn << "\n";
    oss << "stage1_min_score=" << req.stage1_min_score << "\n";
    oss << "stage1_min_score_frac_x10000=" << req.stage1_min_score_frac_x10000 << "\n";
    oss << "stage1_max_freq=" << req.stage1_max_freq << "\n";
    oss << "stage1_max_freq_frac_x10000=" << req.stage1_max_freq_frac_x10000 << "\n";
    oss << "stage2_min_score=" << req.stage2_min_score << "\n";
    oss << "has_stage2_min_score=" << static_cast<int>(req.has_stage2_min_score) << "\n";
    oss << "stage2_max_gap=" << req.stage2_max_gap << "\n";
    oss << "stage2_max_lookback=" << req.stage2_max_lookback << "\n";
    oss << "stage2_min_diag_hits=" << static_cast<int>(req.stage2_min_diag_hits) << "\n";
    oss << "num_results=" << req.num_results << "\n";
    oss << "accept_qdegen=" << static_cast<int>(req.accept_qdegen) << "\n";
    oss << "strand=" << static_cast<int>(req.strand) << "\n";
    oss << "stage3_traceback=" << static_cast<int>(req.stage3_traceback) << "\n";
    oss << "stage3_gapopen=" << req.stage3_gapopen << "\n";
    oss << "stage3_gapext=" << req.stage3_gapext << "\n";
    oss << "stage3_min_pident_x100=" << req.stage3_min_pident_x100 << "\n";
    oss << "stage3_min_nident=" << req.stage3_min_nident << "\n";
    oss << "context_abs=" << req.context_abs << "\n";
    oss << "context_frac_x10000=" << req.context_frac_x10000 << "\n";
    oss << "seqidlist_mode=" << static_cast<int>(req.seqidlist_mode) << "\n";
    oss << "db=" << req.db << "\n";
    oss << "db_total_sequences=" << stats.total_sequences << "\n";
    oss << "db_total_bases=" << stats.total_bases << "\n";
    const char* fmt_str = "tab";
    switch (outfmt) {
        case OutputFormat::kTab:  fmt_str = "tab";  break;
        case OutputFormat::kJson: fmt_str = "json"; break;
        case OutputFormat::kSam:  fmt_str = "sam";  break;
        case OutputFormat::kBam:  fmt_str = "bam";  break;
    }
    oss << "outfmt=" << fmt_str << "\n";
    oss << "seqidlist_sha256=" << seqidlist_sha256 << "\n";
    oss << "neg_seqidlist_sha256=" << neg_seqidlist_sha256 << "\n";
    return oss.str();
}

// ---------------------------------------------------------------------------
// LockGuard
// ---------------------------------------------------------------------------

LockGuard::LockGuard(const std::string& lock_dir)
    : lock_dir_(lock_dir), locked_(false) {
    if (::mkdir(lock_dir_.c_str(), 0700) == 0) {
        locked_ = true;
    }
}

LockGuard::~LockGuard() {
    release();
}

LockGuard::LockGuard(LockGuard&& other) noexcept
    : lock_dir_(std::move(other.lock_dir_)), locked_(other.locked_) {
    other.locked_ = false;
}

LockGuard& LockGuard::operator=(LockGuard&& other) noexcept {
    if (this != &other) {
        release();
        lock_dir_ = std::move(other.lock_dir_);
        locked_ = other.locked_;
        other.locked_ = false;
    }
    return *this;
}

void LockGuard::release() {
    if (locked_) {
        ::rmdir(lock_dir_.c_str());
        locked_ = false;
    }
}

// ---------------------------------------------------------------------------
// Checkpoint helpers
// ---------------------------------------------------------------------------

static std::string basename_of(const std::string& path) {
    auto pos = path.find_last_of('/');
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

static bool file_exists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

static bool dir_exists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

static std::string read_file_string(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) return "";
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

static bool write_file_string(const std::string& path,
                               const std::string& content) {
    std::ofstream out(path, std::ios::binary);
    if (!out.is_open()) return false;
    out << content;
    return out.good();
}

static void remove_recursive(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return;

    if (S_ISDIR(st.st_mode)) {
        DIR* d = ::opendir(path.c_str());
        if (d) {
            struct dirent* ent;
            while ((ent = ::readdir(d)) != nullptr) {
                std::string name(ent->d_name);
                if (name == "." || name == "..") continue;
                remove_recursive(path + "/" + name);
            }
            ::closedir(d);
        }
        ::rmdir(path.c_str());
    } else {
        ::unlink(path.c_str());
    }
}

// ---------------------------------------------------------------------------
// Checkpoint
// ---------------------------------------------------------------------------

Checkpoint::Checkpoint(const Config& cfg, const Logger& logger)
    : cfg_(cfg), logger_(logger) {
    // Build temp dir name:
    // {output_part}.{input_part}.{ix_name}.{kk}.ikafssn.tmp/
    std::string output_part = cfg_.output_path.empty() ? "stdout" : cfg_.output_path;
    std::string input_part = (cfg_.input_path == "-") ? "stdin" : basename_of(cfg_.input_path);

    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", cfg_.resolved_k);

    temp_dir_ = output_part + "." + input_part + "." + cfg_.ix_name + "." + kk + ".ikafssn.tmp";
}

bool Checkpoint::exists() const {
    return dir_exists(temp_dir_);
}

std::string Checkpoint::options_path() const { return temp_dir_ + "/options.txt"; }
std::string Checkpoint::options_sha_path() const { return temp_dir_ + "/options.txt.sha256"; }
std::string Checkpoint::input_sha_path() const {
    std::string input_base = (cfg_.input_path == "-") ? "stdin" : basename_of(cfg_.input_path);
    return temp_dir_ + "/" + input_base + ".sha256";
}
std::string Checkpoint::stdin_fasta_path() const { return temp_dir_ + "/stdin.fasta"; }
std::string Checkpoint::meta_path() const { return temp_dir_ + "/meta.txt"; }
std::string Checkpoint::meta_sha_path() const { return temp_dir_ + "/meta.txt.sha256"; }
std::string Checkpoint::lock_dir_path() const { return temp_dir_ + "/lock"; }

std::string Checkpoint::result_extension() const {
    switch (cfg_.outfmt) {
        case OutputFormat::kJson: return ".json";
        case OutputFormat::kSam:
        case OutputFormat::kBam:  return ".sam";
        default:                  return ".txt";
    }
}

std::string Checkpoint::batch_seqid_path(int n) const {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "/batch_%03d.seqid", n);
    return temp_dir_ + buf;
}

std::string Checkpoint::batch_seqid_sha_path(int n) const {
    return batch_seqid_path(n) + ".sha256";
}

std::string Checkpoint::batch_result_path(int n) const {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "/batch_%03d", n);
    return temp_dir_ + buf + result_extension();
}

std::string Checkpoint::batch_result_sha_path(int n) const {
    return batch_result_path(n) + ".sha256";
}

bool Checkpoint::write_with_sha(const std::string& path,
                                 const std::string& content) {
    if (!write_file_string(path, content)) return false;
    std::string sha = sha256_string(content);
    return write_file_string(path + ".sha256", sha + "\n");
}

bool Checkpoint::validate_sha(const std::string& path) {
    if (!file_exists(path) || !file_exists(path + ".sha256")) return false;
    std::string expected = read_file_string(path + ".sha256");
    // Trim whitespace
    while (!expected.empty() && (expected.back() == '\n' || expected.back() == '\r'))
        expected.pop_back();
    std::string actual = sha256_file(path);
    return !actual.empty() && actual == expected;
}

bool Checkpoint::acquire_lock(LockGuard& guard) {
    guard = LockGuard(lock_dir_path());
    if (!guard.locked()) {
        logger_.error("Cannot acquire checkpoint lock: %s (another process may be running)",
                      lock_dir_path().c_str());
        return false;
    }
    return true;
}

bool Checkpoint::initialize(const std::string& options_text,
                             const std::string& input_sha256,
                             const std::string& stdin_content) {
    if (::mkdir(temp_dir_.c_str(), 0755) != 0 && errno != EEXIST) {
        logger_.error("Cannot create checkpoint directory: %s: %s",
                      temp_dir_.c_str(), std::strerror(errno));
        return false;
    }
    logger_.info("Created checkpoint directory: %s", temp_dir_.c_str());

    // Write options.txt
    if (!write_with_sha(options_path(), options_text)) {
        logger_.error("Failed to write options.txt");
        return false;
    }

    // Write input SHA256
    if (!write_file_string(input_sha_path(), input_sha256 + "\n")) {
        logger_.error("Failed to write input SHA256");
        return false;
    }

    // Save stdin content if applicable
    if (cfg_.input_path == "-" && !stdin_content.empty()) {
        if (!write_file_string(stdin_fasta_path(), stdin_content)) {
            logger_.error("Failed to write stdin.fasta");
            return false;
        }
    }

    return true;
}

bool Checkpoint::resume(const std::string& options_text,
                         const std::string& input_sha256,
                         std::unordered_set<std::string>& completed_seqids,
                         int& next_batch_num) {
    completed_seqids.clear();
    next_batch_num = 0;

    // Validate options.txt
    if (!validate_sha(options_path())) {
        logger_.warn("options.txt SHA256 validation failed");
        return false;
    }
    std::string saved_options = read_file_string(options_path());
    if (saved_options != options_text) {
        logger_.warn("Search options have changed since last run");
        return false;
    }

    // Validate input SHA256
    std::string saved_input_sha = read_file_string(input_sha_path());
    while (!saved_input_sha.empty() &&
           (saved_input_sha.back() == '\n' || saved_input_sha.back() == '\r'))
        saved_input_sha.pop_back();
    if (saved_input_sha != input_sha256) {
        logger_.warn("Input file has changed since last run");
        return false;
    }

    // Scan batch files
    for (int n = 0; ; n++) {
        std::string seqid_path = batch_seqid_path(n);
        if (!file_exists(seqid_path)) break;

        // Validate seqid file
        if (!validate_sha(seqid_path)) {
            logger_.warn("Batch %d seqid file validation failed, truncating", n);
            // Delete this batch and all subsequent
            for (int d = n; ; d++) {
                std::string sp = batch_seqid_path(d);
                if (!file_exists(sp)) break;
                ::unlink(sp.c_str());
                ::unlink(batch_seqid_sha_path(d).c_str());
                std::string rp = batch_result_path(d);
                if (file_exists(rp)) {
                    ::unlink(rp.c_str());
                    ::unlink(batch_result_sha_path(d).c_str());
                }
            }
            break;
        }

        // Check if result file exists and is valid
        std::string result_path = batch_result_path(n);
        if (!file_exists(result_path) || !validate_sha(result_path)) {
            // Delete this batch's seqid and result files
            ::unlink(seqid_path.c_str());
            ::unlink(batch_seqid_sha_path(n).c_str());
            if (file_exists(result_path)) {
                ::unlink(result_path.c_str());
                ::unlink(batch_result_sha_path(n).c_str());
            }
            // Also delete all subsequent batches
            for (int d = n + 1; ; d++) {
                std::string sp = batch_seqid_path(d);
                if (!file_exists(sp)) break;
                ::unlink(sp.c_str());
                ::unlink(batch_seqid_sha_path(d).c_str());
                std::string rp = batch_result_path(d);
                if (file_exists(rp)) {
                    ::unlink(rp.c_str());
                    ::unlink(batch_result_sha_path(d).c_str());
                }
            }
            break;
        }

        // Both valid - add seqids to completed set
        std::string seqid_content = read_file_string(seqid_path);
        std::istringstream iss(seqid_content);
        std::string line;
        while (std::getline(iss, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) completed_seqids.insert(line);
        }

        next_batch_num = n + 1;
    }

    logger_.info("Resuming: %zu completed queries, next batch %d",
                 completed_seqids.size(), next_batch_num);
    return true;
}

bool Checkpoint::write_batch_seqids(int batch_num,
                                     const std::vector<std::string>& seqids) {
    std::ostringstream oss;
    for (const auto& id : seqids) {
        oss << id << "\n";
    }
    return write_with_sha(batch_seqid_path(batch_num), oss.str());
}

bool Checkpoint::write_batch_results(int batch_num,
                                      const std::vector<OutputHit>& hits,
                                      uint8_t mode, uint8_t stage1_score,
                                      bool stage3_traceback) {
    std::string path = batch_result_path(batch_num);

    if (cfg_.outfmt == OutputFormat::kSam || cfg_.outfmt == OutputFormat::kBam) {
        // Write SAM text (always text for intermediate)
        write_results_sam(path, hits, stage1_score);
        // Compute and write SHA256
        std::string sha = sha256_file(path);
        write_file_string(path + ".sha256", sha + "\n");
        return true;
    }

    std::ostringstream oss;
    if (cfg_.outfmt == OutputFormat::kJson) {
        write_results_json_fragment(oss, hits, mode, stage1_score,
                                     stage3_traceback);
    } else {
        // Tab format: write with header for each batch
        write_results_tab(oss, hits, mode, stage1_score, stage3_traceback);
    }
    return write_with_sha(path, oss.str());
}

bool Checkpoint::write_response_meta(uint8_t mode, uint8_t stage1_score,
                                      bool stage3_traceback) {
    std::ostringstream oss;
    oss << "mode=" << static_cast<int>(mode) << "\n";
    oss << "stage1_score=" << static_cast<int>(stage1_score) << "\n";
    oss << "stage3_traceback=" << (stage3_traceback ? 1 : 0) << "\n";
    return write_with_sha(meta_path(), oss.str());
}

bool Checkpoint::read_response_meta(uint8_t& mode, uint8_t& stage1_score,
                                     bool& stage3_traceback) {
    if (!validate_sha(meta_path())) return false;

    std::string content = read_file_string(meta_path());
    std::istringstream iss(content);
    std::string line;
    while (std::getline(iss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        if (key == "mode") mode = static_cast<uint8_t>(std::stoi(val));
        else if (key == "stage1_score") stage1_score = static_cast<uint8_t>(std::stoi(val));
        else if (key == "stage3_traceback") stage3_traceback = (std::stoi(val) != 0);
    }
    return true;
}

bool Checkpoint::merge_results(const std::string& output_path,
                                uint8_t mode, uint8_t stage1_score,
                                bool stage3_traceback) {
    // Collect batch result paths
    std::vector<std::string> batch_paths;
    for (int n = 0; ; n++) {
        std::string path = batch_result_path(n);
        if (!file_exists(path)) break;
        batch_paths.push_back(path);
    }

    if (batch_paths.empty()) {
        logger_.info("No batch results to merge");
        // Still produce empty output for tab/json
        if (cfg_.outfmt == OutputFormat::kTab || cfg_.outfmt == OutputFormat::kJson) {
            std::vector<OutputHit> empty;
            if (output_path.empty()) {
                write_results(std::cout, empty, cfg_.outfmt, mode, stage1_score,
                              stage3_traceback);
            } else {
                std::ofstream out(output_path);
                write_results(out, empty, cfg_.outfmt, mode, stage1_score,
                              stage3_traceback);
            }
        }
        return true;
    }

    if (cfg_.outfmt == OutputFormat::kSam || cfg_.outfmt == OutputFormat::kBam) {
        bool as_bam = (cfg_.outfmt == OutputFormat::kBam);
        std::string out_path = output_path.empty() ? "-" : output_path;
        return merge_sam_files(batch_paths, out_path, as_bam);
    }

    // Tab or JSON merge
    auto open_output = [&]() -> std::ostream* {
        if (output_path.empty()) return &std::cout;
        static std::ofstream ofs;
        ofs.open(output_path);
        return ofs.is_open() ? &ofs : nullptr;
    };

    std::ostream* out = open_output();
    if (!out) {
        logger_.error("Cannot open output file: %s", output_path.c_str());
        return false;
    }

    if (cfg_.outfmt == OutputFormat::kJson) {
        *out << "{\n  \"results\": [\n";
        for (size_t i = 0; i < batch_paths.size(); i++) {
            std::string content = read_file_string(batch_paths[i]);
            // Content has trailing commas after each query object.
            // For the last batch, remove the trailing comma before the final newline.
            if (i + 1 == batch_paths.size() && !content.empty()) {
                // Find last comma and remove it
                auto last_comma = content.rfind(',');
                if (last_comma != std::string::npos) {
                    content.erase(last_comma, 1);
                }
            }
            *out << content;
        }
        *out << "  ]\n}\n";
    } else {
        // Tab format
        for (size_t i = 0; i < batch_paths.size(); i++) {
            std::string content = read_file_string(batch_paths[i]);
            if (i == 0) {
                *out << content;
            } else {
                // Skip header line(s) starting with '#'
                std::istringstream iss(content);
                std::string line;
                while (std::getline(iss, line)) {
                    if (!line.empty() && line[0] == '#') continue;
                    *out << line << "\n";
                }
            }
        }
    }

    return true;
}

void Checkpoint::cleanup() {
    if (!temp_dir_.empty() && dir_exists(temp_dir_)) {
        remove_recursive(temp_dir_);
        logger_.info("Removed checkpoint directory: %s", temp_dir_.c_str());
    }
}

} // namespace ikafssn
