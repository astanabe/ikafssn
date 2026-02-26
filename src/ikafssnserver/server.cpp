#include "ikafssnserver/server.hpp"
#include "ikafssnserver/connection_handler.hpp"
#include "core/config.hpp"
#include "io/volume_discovery.hpp"
#include "util/socket_utils.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <thread>

#include <poll.h>
#include <signal.h>
#include <unistd.h>

#include <tbb/task_arena.h>
#include <tbb/task_group.h>

namespace ikafssn {

Server::~Server() {
    for (int fd : listen_fds_) {
        close_fd(fd);
    }
}

bool Server::load_indexes(const std::string& ix_prefix, const Logger& logger) {
    auto discovered = discover_volumes(ix_prefix);
    if (discovered.empty()) {
        logger.error("No index files found for prefix %s", ix_prefix.c_str());
        return false;
    }

    auto prefix_parts = parse_index_prefix(ix_prefix);

    // Group by k and open index files
    for (const auto& dv : discovered) {
        auto& group = kmer_groups_[dv.k];
        group.k = dv.k;
        group.kmer_type = kmer_type_for_k(dv.k);

        ServerVolumeData svd;
        svd.volume_index = dv.volume_index;

        if (!svd.kix.open(dv.kix_path)) {
            logger.error("Cannot open %s", dv.kix_path.c_str());
            return false;
        }
        if (!svd.kpx.open(dv.kpx_path)) {
            logger.error("Cannot open %s", dv.kpx_path.c_str());
            return false;
        }
        if (!svd.ksx.open(dv.ksx_path)) {
            logger.error("Cannot open %s", dv.ksx_path.c_str());
            return false;
        }

        group.volumes.push_back(std::move(svd));
    }

    // Sort volumes within each group, then open shared .khx per group
    for (auto& [k, group] : kmer_groups_) {
        std::sort(group.volumes.begin(), group.volumes.end(),
                  [](const ServerVolumeData& a, const ServerVolumeData& b) {
                      return a.volume_index < b.volume_index;
                  });

        // Open shared .khx for this k-mer group (non-fatal if missing)
        group.khx.open(khx_path_for(prefix_parts.parent_dir, prefix_parts.db_name, k)); // non-fatal
    }

    // Default k = largest available
    default_k_ = kmer_groups_.rbegin()->first;

    logger.info("Loaded %zu k-mer group(s):", kmer_groups_.size());
    for (const auto& [k, group] : kmer_groups_) {
        logger.info("  k=%d: %zu volume(s)", k, group.volumes.size());
    }

    return true;
}

void Server::request_shutdown() {
    shutdown_requested_.store(true, std::memory_order_release);
}

int Server::try_acquire_sequences(int n) {
    std::lock_guard<std::mutex> lock(seq_mutex_);
    int capped = std::min(n, max_seqs_per_req_);
    int available = max_active_sequences_ - active_sequences_;
    int acquired = std::min(capped, std::max(0, available));
    active_sequences_ += acquired;
    return acquired;
}

void Server::release_sequences(int n) {
    std::lock_guard<std::mutex> lock(seq_mutex_);
    active_sequences_ -= n;
}

void Server::write_pid_file(const std::string& path) {
    std::ofstream f(path);
    if (f.is_open()) {
        f << ::getpid() << "\n";
    }
}

void Server::accept_loop(int listen_fd, const ServerConfig& config, const Logger& logger) {
    int num_threads = config.num_threads;
    if (num_threads <= 0) {
        num_threads = static_cast<int>(std::thread::hardware_concurrency());
        if (num_threads <= 0) num_threads = 1;
    }

    tbb::task_arena arena(num_threads);
    tbb::task_group tg;

    while (!shutdown_requested_.load(std::memory_order_acquire)) {
        // Use poll with timeout to allow shutdown check
        struct pollfd pfd;
        pfd.fd = listen_fd;
        pfd.events = POLLIN;
        int ret = ::poll(&pfd, 1, 500); // 500ms timeout

        if (ret < 0) {
            if (errno == EINTR) continue;
            logger.error("poll() failed: %s", strerror(errno));
            break;
        }

        if (ret == 0) continue; // timeout

        int client_fd = accept_connection(listen_fd);
        if (client_fd < 0) {
            if (errno == EINTR) continue;
            if (shutdown_requested_.load(std::memory_order_acquire)) break;
            logger.error("accept() failed: %s", strerror(errno));
            continue;
        }

        logger.debug("Accepted connection (fd=%d)", client_fd);

        // Dispatch to worker thread unconditionally (no connection-level gating).
        // Per-sequence gating happens inside process_search_request().
        arena.execute([&, client_fd] {
            tg.run([&, client_fd] {
                handle_connection(client_fd, kmer_groups_, default_k_,
                                  config.search_config,
                                  config.stage3_config, config.db_path,
                                  config.context_is_ratio, config.context_ratio,
                                  config.context_abs,
                                  *this, arena, logger);
            });
        });
    }

    // Wait for in-flight requests to complete
    logger.info("Waiting for in-flight requests...");
    arena.execute([&] { tg.wait(); });
}

int Server::run(const ServerConfig& config_in) {
    // Mutable copy so we can resolve -max_freq fraction
    ServerConfig config = config_in;
    Logger logger(config.log_level);

    // Load indexes
    if (!load_indexes(config.ix_prefix, logger)) {
        return 1;
    }

    // Resolve -max_freq fraction using total NSEQ from loaded indexes
    if (config.max_freq_raw > 0 && config.max_freq_raw < 1.0) {
        if (!kmer_groups_.empty()) {
            uint64_t total_nseq = 0;
            for (const auto& vol : kmer_groups_.begin()->second.volumes)
                total_nseq += vol.ksx.num_sequences();
            auto resolved = static_cast<uint32_t>(
                std::ceil(config.max_freq_raw * total_nseq));
            if (resolved == 0) resolved = 1;
            config.search_config.stage1.max_freq = resolved;
            logger.info("-stage1_max_freq=%.6g (fraction) -> threshold=%u (total_nseq=%lu)",
                        config.max_freq_raw, resolved,
                        static_cast<unsigned long>(total_nseq));
        }
    } else {
        config.search_config.stage1.max_freq =
            static_cast<uint32_t>(config.max_freq_raw);
    }

    // Write PID file if requested
    if (!config.pid_file.empty()) {
        write_pid_file(config.pid_file);
    }

    // Resolve per-sequence concurrency limits
    {
        int num_threads = config.num_threads;
        if (num_threads <= 0) {
            num_threads = static_cast<int>(std::thread::hardware_concurrency());
            if (num_threads <= 0) num_threads = 1;
        }
        max_active_sequences_ = config.max_query > 0 ? config.max_query : 1024;
        max_seqs_per_req_ = config.max_seqs_per_req > 0 ? config.max_seqs_per_req : num_threads;
        logger.info("Max concurrent sequences: %d, max per request: %d",
                     max_active_sequences_, max_seqs_per_req_);
    }

    // Set up listening sockets
    if (config.unix_socket_path.empty() && config.tcp_addr.empty()) {
        logger.error("At least one of -socket or -tcp must be specified");
        return 1;
    }

    int listen_fd = -1;

    if (!config.unix_socket_path.empty()) {
        listen_fd = unix_listen(config.unix_socket_path);
        if (listen_fd < 0) {
            logger.error("Cannot listen on UNIX socket %s", config.unix_socket_path.c_str());
            return 1;
        }
        listen_fds_.push_back(listen_fd);
        logger.info("Listening on UNIX socket: %s", config.unix_socket_path.c_str());
    }

    if (!config.tcp_addr.empty()) {
        int tcp_fd = tcp_listen(config.tcp_addr);
        if (tcp_fd < 0) {
            logger.error("Cannot listen on TCP %s", config.tcp_addr.c_str());
            return 1;
        }
        listen_fds_.push_back(tcp_fd);
        logger.info("Listening on TCP: %s", config.tcp_addr.c_str());

        // If both UNIX and TCP, use TCP as primary accept loop.
        // For simplicity in this initial implementation, we start accept loops
        // for each socket in separate threads.
        if (listen_fd >= 0) {
            // Start UNIX socket accept loop in background thread
            std::thread unix_thread([this, listen_fd, &config, &logger] {
                accept_loop(listen_fd, config, logger);
            });

            // Run TCP accept loop in main thread
            accept_loop(tcp_fd, config, logger);

            unix_thread.join();
            return 0;
        }

        listen_fd = tcp_fd;
    }

    logger.info("Server ready, default k=%d", default_k_);

    // Run accept loop (blocks until shutdown)
    accept_loop(listen_fd, config, logger);

    // Cleanup
    if (!config.unix_socket_path.empty()) {
        ::unlink(config.unix_socket_path.c_str());
    }
    if (!config.pid_file.empty()) {
        ::unlink(config.pid_file.c_str());
    }

    logger.info("Server shut down");
    return 0;
}

} // namespace ikafssn
