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

bool Server::load_database(const std::string& ix_prefix, const std::string& db_path,
                           const ServerConfig& config, const Logger& logger) {
    auto prefix_parts = parse_index_prefix(ix_prefix);
    const std::string& db_name = prefix_parts.db;

    // Check for duplicate DB name
    if (db_index_.count(db_name)) {
        logger.error("Duplicate database name '%s' (from %s)", db_name.c_str(), ix_prefix.c_str());
        return false;
    }

    auto discovered = discover_volumes(ix_prefix);
    if (discovered.empty()) {
        logger.error("No index files found for prefix %s", ix_prefix.c_str());
        return false;
    }

    DatabaseEntry entry;
    entry.name = db_name;
    entry.ix_prefix = ix_prefix;
    entry.db_path = db_path;
    entry.max_mode = db_path.empty() ? 2 : 3;

    bool all_have_kpx = true;

    // Group by k and open index files
    for (const auto& dv : discovered) {
        auto& group = entry.kmer_groups[dv.k];
        group.k = dv.k;
        group.kmer_type = kmer_type_for_k(dv.k);

        ServerVolumeData svd;
        svd.volume_index = dv.volume_index;

        if (!svd.kix.open(dv.kix_path)) {
            logger.error("Cannot open %s", dv.kix_path.c_str());
            return false;
        }
        if (dv.has_kpx) {
            if (!svd.kpx.open(dv.kpx_path)) {
                logger.error("Cannot open %s", dv.kpx_path.c_str());
                return false;
            }
        } else {
            all_have_kpx = false;
        }
        if (!svd.ksx.open(dv.ksx_path)) {
            logger.error("Cannot open %s", dv.ksx_path.c_str());
            return false;
        }

        svd.total_bases = 0;
        for (uint32_t oid = 0; oid < svd.ksx.num_sequences(); oid++) {
            svd.total_bases += svd.ksx.seq_length(oid);
        }

        group.volumes.push_back(std::move(svd));
    }

    // Restrict max_mode if .kpx files are missing (mode 1 index)
    if (!all_have_kpx) {
        entry.max_mode = 1;
        logger.info("DB '%s': .kpx files missing, max_mode restricted to 1", db_name.c_str());
    }

    // Sort volumes within each group, then open shared .khx per group
    for (auto& [k, group] : entry.kmer_groups) {
        std::sort(group.volumes.begin(), group.volumes.end(),
                  [](const ServerVolumeData& a, const ServerVolumeData& b) {
                      return a.volume_index < b.volume_index;
                  });

        // Open shared .khx for this k-mer group (non-fatal if missing)
        group.khx.open(khx_path_for(prefix_parts.parent_dir, prefix_parts.db, k));
    }

    // Default k = largest available
    entry.default_k = entry.kmer_groups.rbegin()->first;

    // Resolve search config from server config template
    entry.resolved_search_config = config.search_config;

    // Resolve -max_freq: 1/1.0 = disable, fraction -> absolute, else integer
    if (config.max_freq_raw == 1.0) {
        entry.resolved_search_config.stage1.max_freq = Stage1Config::MAX_FREQ_DISABLED;
        logger.info("DB '%s': -stage1_max_freq=1 -> high-frequency k-mer filtering disabled",
                    db_name.c_str());
    } else if (config.max_freq_raw > 0 && config.max_freq_raw < 1.0) {
        if (!entry.kmer_groups.empty()) {
            uint64_t total_nseq = 0;
            for (const auto& vol : entry.kmer_groups.begin()->second.volumes)
                total_nseq += vol.ksx.num_sequences();
            auto resolved = static_cast<uint32_t>(
                std::ceil(config.max_freq_raw * total_nseq));
            if (resolved == 0) resolved = 1;
            entry.resolved_search_config.stage1.max_freq = resolved;
            logger.info("DB '%s': -stage1_max_freq=%.6g (fraction) -> threshold=%u (total_nseq=%lu)",
                        db_name.c_str(), config.max_freq_raw, resolved,
                        static_cast<unsigned long>(total_nseq));
        }
    } else {
        entry.resolved_search_config.stage1.max_freq =
            static_cast<uint32_t>(config.max_freq_raw);
    }

    // Copy stage3/context params from ServerConfig
    entry.stage3_config = config.stage3_config;
    entry.context_is_ratio = config.context_is_ratio;
    entry.context_ratio = config.context_ratio;
    entry.context_abs = config.context_abs;

    logger.info("Loaded DB '%s' (%zu k-mer group(s)):", db_name.c_str(), entry.kmer_groups.size());
    for (const auto& [k, group] : entry.kmer_groups) {
        logger.info("  k=%d: %zu volume(s)", k, group.volumes.size());
    }

    size_t idx = databases_.size();
    databases_.push_back(std::move(entry));
    db_index_[db_name] = idx;

    return true;
}

const DatabaseEntry* Server::find_database(const std::string& name) const {
    auto it = db_index_.find(name);
    if (it == db_index_.end()) return nullptr;
    return &databases_[it->second];
}

int Server::default_k() const {
    if (databases_.empty()) return 0;
    return databases_[0].default_k;
}

void Server::request_shutdown() {
    shutdown_requested_.store(true, std::memory_order_release);
}

int Server::try_acquire_sequences(int n) {
    std::lock_guard<std::mutex> lock(seq_mutex_);
    int capped = std::min(n, max_seqs_per_req_);
    int available = max_queue_size_ - queue_depth_;
    int acquired = std::min(capped, std::max(0, available));
    queue_depth_ += acquired;
    return acquired;
}

void Server::release_sequences(int n) {
    std::lock_guard<std::mutex> lock(seq_mutex_);
    queue_depth_ -= n;
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

        // Dispatch to worker thread
        arena.execute([&, client_fd] {
            tg.run([&, client_fd] {
                handle_connection(client_fd, *this, config, arena, logger);
            });
        });
    }

    // Wait for in-flight requests to complete
    logger.info("Waiting for in-flight requests...");
    arena.execute([&] { tg.wait(); });
}

int Server::run(const ServerConfig& config_in) {
    ServerConfig config = config_in;
    Logger logger(config.log_level);

    // Load all databases
    for (const auto& db_entry : config.db_entries) {
        if (!load_database(db_entry.ix_prefix, db_entry.db_path, config, logger)) {
            return 1;
        }
    }

    if (databases_.empty()) {
        logger.error("No databases loaded");
        return 1;
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
        max_queue_size_ = config.max_queue_size > 0 ? config.max_queue_size : 1024;
        max_seqs_per_req_ = config.max_seqs_per_req > 0 ? config.max_seqs_per_req : num_threads;
        logger.info("Max concurrent sequences: %d, max per request: %d",
                     max_queue_size_, max_seqs_per_req_);
    }

    // Log total mmap count
    size_t total_mmaps = 0;
    for (const auto& db : databases_) {
        for (const auto& [k, group] : db.kmer_groups) {
            for (const auto& vol : group.volumes) {
                total_mmaps += 2; // kix + ksx always
                if (vol.kpx.is_open()) total_mmaps++;
            }
            if (group.khx.is_open()) total_mmaps++;
        }
    }
    logger.info("Total mmap'd files across %zu DB(s): %zu", databases_.size(), total_mmaps);

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
        if (listen_fd >= 0) {
            std::thread unix_thread([this, listen_fd, &config, &logger] {
                accept_loop(listen_fd, config, logger);
            });

            accept_loop(tcp_fd, config, logger);

            unix_thread.join();
            return 0;
        }

        listen_fd = tcp_fd;
    }

    logger.info("Server ready, %zu database(s), default k=%d", databases_.size(), default_k());

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
