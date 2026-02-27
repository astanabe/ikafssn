#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "io/result_writer.hpp"
#include "protocol/messages.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// SHA256 utilities (using OpenSSL EVP)
std::string sha256_file(const std::string& path);
std::string sha256_string(const std::string& data);

// DB statistics resolved from InfoResponse
struct DbStats {
    std::string db_name;
    uint64_t total_sequences = 0;
    uint64_t total_bases = 0;
};

DbStats resolve_db_stats(const InfoResponse& info, const std::string& db_name,
                          uint8_t k);

// Build a canonical options text for checkpoint validation.
std::string build_options_text(const SearchRequest& req, const DbStats& stats,
                                uint8_t resolved_k, OutputFormat outfmt,
                                const std::string& seqidlist_sha256,
                                const std::string& neg_seqidlist_sha256);

// RAII lock guard for directory-based locking
class LockGuard {
public:
    LockGuard() = default;
    explicit LockGuard(const std::string& lock_dir);
    ~LockGuard();
    LockGuard(const LockGuard&) = delete;
    LockGuard& operator=(const LockGuard&) = delete;
    LockGuard(LockGuard&& other) noexcept;
    LockGuard& operator=(LockGuard&& other) noexcept;
    bool locked() const { return locked_; }
    void release();

private:
    std::string lock_dir_;
    bool locked_ = false;
};

class Checkpoint {
public:
    struct Config {
        std::string output_path;   // "" = stdout
        std::string input_path;    // "-" = stdin
        std::string ix_name;       // -ix value (DB name)
        uint8_t resolved_k;        // actual k value (never 0)
        OutputFormat outfmt;
    };

    Checkpoint(const Config& cfg, const Logger& logger);

    const std::string& temp_dir() const { return temp_dir_; }
    bool exists() const;

    // First run: create temp dir, save options.txt + input SHA256
    bool initialize(const std::string& options_text,
                    const std::string& input_sha256,
                    const std::string& stdin_content);

    // Resume: validate, return completed seqids and next batch number.
    // Returns false on validation failure (caller should cleanup + initialize).
    bool resume(const std::string& options_text,
                const std::string& input_sha256,
                std::unordered_set<std::string>& completed_seqids,
                int& next_batch_num);

    bool acquire_lock(LockGuard& guard);

    // Write batch seqid list (before sending)
    bool write_batch_seqids(int batch_num,
                            const std::vector<std::string>& seqids);

    // Write batch results (after receiving)
    bool write_batch_results(int batch_num,
                             const std::vector<OutputHit>& hits,
                             uint8_t mode, uint8_t stage1_score,
                             bool stage3_traceback);

    // Save/load response metadata (mode, stage1_score, stage3_traceback)
    bool write_response_meta(uint8_t mode, uint8_t stage1_score,
                             bool stage3_traceback);
    bool read_response_meta(uint8_t& mode, uint8_t& stage1_score,
                            bool& stage3_traceback);

    // Merge all batch results to final output
    bool merge_results(const std::string& output_path,
                       uint8_t mode, uint8_t stage1_score,
                       bool stage3_traceback);

    void cleanup();  // remove temp directory

private:
    Config cfg_;
    const Logger& logger_;
    std::string temp_dir_;

    // File path helpers
    std::string options_path() const;
    std::string options_sha_path() const;
    std::string input_sha_path() const;
    std::string stdin_fasta_path() const;
    std::string meta_path() const;
    std::string meta_sha_path() const;
    std::string lock_dir_path() const;
    std::string batch_seqid_path(int n) const;
    std::string batch_seqid_sha_path(int n) const;
    std::string batch_result_path(int n) const;
    std::string batch_result_sha_path(int n) const;
    std::string result_extension() const;

    // Write file + its SHA256 sidecar
    bool write_with_sha(const std::string& path, const std::string& content);
    // Validate file against its SHA256 sidecar
    bool validate_sha(const std::string& path);
};

} // namespace ikafssn
