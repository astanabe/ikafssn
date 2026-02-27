#pragma once

#include <atomic>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "ikafssnserver/request_processor.hpp"
#include "search/stage3_alignment.hpp"
#include "search/volume_searcher.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// Per-database entry holding loaded indexes and resolved config
struct DatabaseEntry {
    std::string name;                       // DB name (basename of ix_prefix)
    std::string ix_prefix;                  // original (for logging)
    std::string db_path;                    // BLAST DB path (empty = max_mode 2)
    std::map<int, KmerGroup> kmer_groups;
    int default_k = 0;                      // largest k for this DB
    uint8_t max_mode = 2;                   // 2 or 3
    SearchConfig resolved_search_config;    // max_freq resolved per-DB
    Stage3Config stage3_config;
    bool context_is_ratio = false;
    double context_ratio = 0.0;
    uint32_t context_abs = 0;
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
    int num_threads = 0;            // 0 = auto-detect
    int max_queue_size = 0;         // 0 = default (1024). Max total in-flight sequences globally
    int max_seqs_per_req = 0;       // 0 = default (same as resolved thread count). Per-request cap
    int shutdown_timeout = 30;      // seconds
    SearchConfig search_config;
    double max_freq_raw = 0.5;          // raw -max_freq value (fraction or absolute)
    Logger::Level log_level = Logger::kInfo;
    // Stage 3 / BLAST DB config
    Stage3Config stage3_config;     // default stage3 config
    bool context_is_ratio = false;
    double context_ratio = 0.0;
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
    int max_seqs_per_req() const { return max_seqs_per_req_; }

    // Get current queue depth.
    int queue_depth() const { return queue_depth_; }

    // Non-blocking: try to acquire up to n permits (capped by max_seqs_per_req_).
    // Returns count actually acquired.
    int try_acquire_sequences(int n);

    // Release n permits.
    void release_sequences(int n);

private:
    std::vector<DatabaseEntry> databases_;
    std::unordered_map<std::string, size_t> db_index_;
    std::atomic<bool> shutdown_requested_{false};
    std::vector<int> listen_fds_;

    std::mutex seq_mutex_;
    int queue_depth_ = 0;
    int max_queue_size_ = 1024;  // from -max_queue_size, default 1024; overridden in run()
    int max_seqs_per_req_ = 1024;      // from -max_seqs_per_req, default = threads; overridden in run()

    void accept_loop(int listen_fd, const ServerConfig& config, const Logger& logger);
    void write_pid_file(const std::string& path);
};

} // namespace ikafssn
