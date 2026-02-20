#pragma once

#include <atomic>
#include <map>
#include <string>
#include <vector>

#include "ikafssnserver/request_processor.hpp"
#include "search/volume_searcher.hpp"
#include "util/logger.hpp"

namespace ikafssn {

struct ServerConfig {
    std::string ix_dir;
    std::string unix_socket_path;
    std::string tcp_addr;           // "host:port"
    std::string pid_file;
    int num_threads = 0;            // 0 = auto-detect
    int shutdown_timeout = 30;      // seconds
    SearchConfig search_config;
    Logger::Level log_level = Logger::kInfo;
};

class Server {
public:
    Server() = default;
    ~Server();

    // Load indexes from the directory.
    bool load_indexes(const std::string& ix_dir, const Logger& logger);

    // Run the server (blocking). Returns exit code.
    int run(const ServerConfig& config);

    // Request graceful shutdown (called from signal handler).
    void request_shutdown();

    // Get the default k value.
    int default_k() const { return default_k_; }

    // Get the kmer groups (for testing).
    const std::map<int, KmerGroup>& kmer_groups() const { return kmer_groups_; }

private:
    std::map<int, KmerGroup> kmer_groups_;
    int default_k_ = 0;
    std::atomic<bool> shutdown_requested_{false};
    std::vector<int> listen_fds_;

    void accept_loop(int listen_fd, const ServerConfig& config, const Logger& logger);
    void write_pid_file(const std::string& path);
};

} // namespace ikafssn
