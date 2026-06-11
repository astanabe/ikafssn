#pragma once

#include <atomic>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ikafssnserver/budget_pool.hpp"
#include "ikafssnserver/request_processor.hpp"
#include "search/stage3_alignment.hpp"
#include "search/search_config.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// Per-database entry holding loaded indexes and resolved config
struct DatabaseEntry {
    std::string name;                       // DB name (basename of ix_prefix)
    std::string ix_prefix;                  // original (for logging)
    std::string db_path;                    // BLAST DB path (empty = max_mode 2)
    std::vector<KmerGroup> kmer_groups;
    int default_k = 0;                      // largest k for this DB
    uint8_t default_t = 0;
    uint8_t default_template_type = 0;
    uint8_t max_mode = 2;                   // 2 or 3
    SearchConfig resolved_search_config;    // resolved per-DB
    Stage3Config stage3_config;
    bool context_is_ratio = true;
    double context_ratio = 2.0;
    uint32_t context_abs = 0;

    // Prefer the group whose max_degen_expand matches the query value; on a
    // tie, prefer the larger max_degen_expand.
    static const KmerGroup* tie_break_expand(const KmerGroup* a,
                                             const KmerGroup* b,
                                             uint32_t query_expand) {
        if (!a) return b;
        if (!b) return a;
        bool aq = a->max_degen_expand == query_expand;
        bool bq = b->max_degen_expand == query_expand;
        if (aq != bq) return aq ? a : b;
        return a->max_degen_expand >= b->max_degen_expand ? a : b;
    }

    // Select the single group that matches the request's 7-field identity
    // (k / t / template_type + the four indexing fields), tie-breaking the
    // 8th (max_degen_expand) toward query_expand.  Returns nullptr if no
    // group matches.
    const KmerGroup* find_group(int k, uint8_t t, uint8_t template_type,
                                uint32_t min_seq_length,
                                uint32_t min_length_split,
                                uint32_t overlap_length,
                                uint64_t max_freq_build,
                                uint32_t query_expand) const {
        const KmerGroup* best = nullptr;
        for (const auto& g : kmer_groups) {
            if (g.k == k && g.t == t && g.template_type == template_type &&
                g.min_seq_length == min_seq_length &&
                g.min_length_split == min_length_split &&
                g.overlap_length == overlap_length &&
                g.max_freq_build == max_freq_build) {
                best = tie_break_expand(best, &g, query_expand);
            }
        }
        return best;
    }

    // Select a coding+optimal pair sharing the 6-field key and ONE common
    // max_degen_expand (tie-broken toward query_expand, else the largest value
    // present on both sides).  No max_degen_expand mixing between the sides.
    // Returns {nullptr, nullptr} if no pair shares a max_degen_expand.
    std::pair<const KmerGroup*, const KmerGroup*>
    find_both_groups(int k, uint8_t t,
                     uint32_t min_seq_length, uint32_t min_length_split,
                     uint32_t overlap_length, uint64_t max_freq_build,
                     uint32_t query_expand) const {
        auto matches6 = [&](const KmerGroup& g, uint8_t tt) {
            return g.k == k && g.t == t && g.template_type == tt &&
                   g.min_seq_length == min_seq_length &&
                   g.min_length_split == min_length_split &&
                   g.overlap_length == overlap_length &&
                   g.max_freq_build == max_freq_build;
        };
        std::set<uint32_t> cod_e, opt_e;
        for (const auto& g : kmer_groups) {
            if (matches6(g, 1)) cod_e.insert(g.max_degen_expand);
            if (matches6(g, 2)) opt_e.insert(g.max_degen_expand);
        }
        std::vector<uint32_t> common;
        for (uint32_t e : cod_e)
            if (opt_e.count(e)) common.push_back(e);
        if (common.empty()) return {nullptr, nullptr};
        uint32_t chosen = choose_max_degen_expand(common, query_expand);
        const KmerGroup* cod = nullptr;
        const KmerGroup* opt = nullptr;
        for (const auto& g : kmer_groups) {
            if (matches6(g, 1) && g.max_degen_expand == chosen) cod = &g;
            if (matches6(g, 2) && g.max_degen_expand == chosen) opt = &g;
        }
        return {cod, opt};
    }
};

struct ServerConfig {
    struct DbEntry {
        std::string ix_prefix;
        std::string db_path;  // empty = defaults to ix_prefix
    };
    std::vector<DbEntry> db_entries;
    std::string unix_socket_path;
    std::string tcp_addr;           // "host:port"
    std::string pid_file;
    int nthread = 0;            // 0 = auto-detect
    int max_queue_size = 0;         // 0 = default (1024). Max total in-flight sequences globally
    int max_nseq_per_req = 0;       // 0 = default (same as resolved thread count). Per-request cap
    int max_concurrent_search = 0;  // 0 = unlimited (pass-through). When >= 1, requests share
                                    // posting_budget so at most N run the budget-bound stages.
    int shutdown_timeout = 30;      // seconds
    uint64_t memory_limit = 0;     // 0 = auto (half of RAM)
    SearchConfig search_config;
    // Index-variant load filter (ikafssninfo-local semantics: an unset -t is a
    // wildcard).  max_degen_expand is deliberately not filtered here so every
    // build-time max_degen_expand variant stays loadable and the per-request
    // tie-break can choose among them.
    VariantFilter variant_filter;
    Logger::Level log_level = Logger::kInfo;
    // Stage 3 / BLAST DB config
    Stage3Config stage3_config;     // default stage3 config
    bool context_is_ratio = true;
    double context_ratio = 2.0;
    uint32_t context_abs = 0;
};

class Server {
public:
    Server() = default;
    ~Server();

    // Load a single database's indexes.
    bool load_database(const std::string& ix_prefix, const std::string& db_path,
                       const ServerConfig& config, const Logger& logger);

    // Find a database by name. Returns nullptr if not found.
    const DatabaseEntry* find_database(const std::string& name) const;

    // Get all loaded databases.
    const std::vector<DatabaseEntry>& databases() const { return databases_; }

    // Run the server (blocking). Returns exit code.
    int run(const ServerConfig& config);

    // Request graceful shutdown (called from signal handler).
    void request_shutdown();

    // Get the default k value (first DB's default_k).
    int default_k() const;

    // Get max queue size limit.
    int max_queue_size() const { return max_queue_size_; }

    // Get per-request sequence cap.
    int max_nseq_per_req() const { return max_nseq_per_req_; }

    // Get current queue depth.
    int queue_depth() const { return queue_depth_; }

    // Non-blocking: try to acquire up to n permits (capped by max_nseq_per_req_).
    // Returns count actually acquired.
    int try_acquire_sequences(int n);

    // Release n permits.
    void release_sequences(int n);

    // Per-request posting budget (residual of -memory_limit after khx/ksx).
    uint64_t posting_budget() const { return posting_budget_; }

    // Acquire a lease against the inter-request posting-budget pool.  The
    // pool is configured by run() based on -max_concurrent_search.  In
    // pass-through mode (default) this returns lease(posting_budget())
    // immediately and never blocks.
    BudgetLease acquire_posting_budget(uint64_t min, uint64_t max) {
        return pool_.acquire(min, max);
    }

private:
    // Apply persistent madvise WILLNEED for .khx + .ksx within `budget` and
    // expose the residual as posting_budget_.  Invoked by run().
    void apply_madvise_budget(uint64_t budget, const Logger& logger);

    std::vector<DatabaseEntry> databases_;
    std::unordered_map<std::string, size_t> db_index_;
    std::atomic<bool> shutdown_requested_{false};
    std::vector<int> listen_fds_;

    uint64_t posting_budget_ = 0;
    BudgetPool pool_;

    std::mutex seq_mutex_;
    int queue_depth_ = 0;
    int max_queue_size_ = 1024;  // from -max_queue_size, default 1024; overridden in run()
    int max_nseq_per_req_ = 1024;      // from -max_nseq_per_req, default = threads; overridden in run()

    void accept_loop(int listen_fd, const ServerConfig& config, const Logger& logger);
    void write_pid_file(const std::string& path);
};

} // namespace ikafssn
