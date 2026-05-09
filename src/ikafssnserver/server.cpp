#include "ikafssnserver/server.hpp"
#include "ikafssnserver/connection_handler.hpp"
#include "core/config.hpp"
#include "core/spaced_seed.hpp"
#include "core/version.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "io/volume_discovery.hpp"
#include "util/common_init.hpp"
#include "util/socket_utils.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
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

    // Helper to find or create a KmerGroup by (k, t, template_type)
    auto find_or_create_group = [&](int k, uint8_t t, uint8_t tt) -> KmerGroup& {
        for (auto& g : entry.kmer_groups) {
            if (g.k == k && g.t == t && g.template_type == tt) return g;
        }
        entry.kmer_groups.push_back({});
        auto& g = entry.kmer_groups.back();
        g.k = k;
        g.kmer_type = kmer_type_for(k, t);
        g.t = t;
        g.template_type = tt;
        return g;
    };

    // Group by (k, t, template_type), validate kix/kpx headers, capture sizes.
    // .kix and .kpx are opened only briefly for header validation here;
    // the search path re-opens them per-request so concurrent requests do
    // not contend on shared madvise calls against the same mappings.
    for (const auto& dv : discovered) {
        auto& group = find_or_create_group(dv.k, dv.t, dv.template_type);

        ServerVolumeData svd;
        svd.files = dv;
        svd.volume_index = dv.volume_index;

        {
            KixReader kix_probe;
            if (!kix_probe.open(dv.kix_path)) {
                logger.error("Cannot open %s", dv.kix_path.c_str());
                return false;
            }
            svd.kix_posting_size       = kix_probe.posting_file_size();
            svd.kix_full_size          = kix_probe.willneed_size_full();
            svd.num_sequences          = kix_probe.num_sequences();
            svd.total_distinct_postings = kix_probe.total_distinct_postings();
            const auto& kix_hdr = kix_probe.header();
            svd.db_name = std::string(kix_hdr.db,
                                      strnlen(kix_hdr.db, sizeof(kix_hdr.db)));
            // v10: capture the fragment-indexing triplet from the first
            // index seen, then verify every subsequent index in this DB
            // agrees.  Mismatches mean the user has mixed indexes built
            // with different filter / split parameters under the same
            // prefix — refuse to load.
            const uint32_t kix_min     = kix_probe.min_seq_length();
            const uint32_t kix_split   = kix_probe.min_length_split();
            const uint32_t kix_overlap = kix_probe.overlap_length();
            if (entry.kmer_groups.size() == 1 && group.volumes.empty()) {
                // First index seen for this DB.
                entry.min_seq_length   = kix_min;
                entry.min_length_split = kix_split;
                entry.overlap_length   = kix_overlap;
            } else if (entry.min_seq_length   != kix_min   ||
                       entry.min_length_split != kix_split ||
                       entry.overlap_length   != kix_overlap) {
                logger.error("Index %s has min_seq_length=%u min_length_split=%u "
                             "overlap_length=%u but the DB '%s' was already loaded "
                             "with %u/%u/%u",
                             dv.kix_path.c_str(),
                             kix_min, kix_split, kix_overlap,
                             db_name.c_str(),
                             entry.min_seq_length, entry.min_length_split,
                             entry.overlap_length);
                kix_probe.close();
                return false;
            }
            kix_probe.close();
        }

        if (dv.has_kpx) {
            KpxReader kpx_probe;
            if (!kpx_probe.open(dv.kpx_path)) {
                logger.error("Cannot open %s", dv.kpx_path.c_str());
                return false;
            }
            svd.kpx_posting_size = kpx_probe.posting_file_size();
            svd.kpx_full_size    = kpx_probe.willneed_size_full();
            kpx_probe.close();
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
    for (auto& group : entry.kmer_groups) {
        std::sort(group.volumes.begin(), group.volumes.end(),
                  [](const ServerVolumeData& a, const ServerVolumeData& b) {
                      return a.volume_index < b.volume_index;
                  });

        // Open shared .khx for this k-mer group (non-fatal if missing)
        group.khx.open(khx_path_for(prefix_parts.parent_dir, prefix_parts.db,
                                     group.k, group.t, group.template_type));
    }

    // Default k = largest available (groups sorted by k ascending)
    if (!entry.kmer_groups.empty()) {
        auto max_k_it = std::max_element(entry.kmer_groups.begin(), entry.kmer_groups.end(),
            [](const KmerGroup& a, const KmerGroup& b) { return a.k < b.k; });
        entry.default_k = max_k_it->k;
        entry.default_t = max_k_it->t;
        // If both coding and optimal exist for this (k, t), default to "both" search
        if (entry.find_group(max_k_it->k, max_k_it->t, 1) &&
            entry.find_group(max_k_it->k, max_k_it->t, 2)) {
            entry.default_template_type = 3; // kBoth
        } else {
            entry.default_template_type = max_k_it->template_type;
        }
    }

    // Resolve search config from server config template
    entry.resolved_search_config = config.search_config;
    // v10: the server has no -min_query_length CLI flag of its own.  It
    // adopts the loaded index's min_seq_length as the effective floor so
    // queries shorter than the index's filter are skipped server-side
    // (defence-in-depth — well-behaved clients reject these locally).
    entry.resolved_search_config.min_query_length = entry.min_seq_length;
    // v10 (Phase 3): adopt the index's overlap_length as the upper bound so
    // queries longer than overlap_length are skipped with kSkipQueryTooLong.
    // overlap_length == 0 (no fragment splitting) leaves the check disabled.
    entry.resolved_search_config.max_query_length = entry.overlap_length;

    // Copy stage3/context params from ServerConfig
    entry.stage3_config = config.stage3_config;
    entry.context_is_ratio = config.context_is_ratio;
    entry.context_ratio = config.context_ratio;
    entry.context_abs = config.context_abs;

    logger.info("Loaded DB '%s' (%zu k-mer group(s)):", db_name.c_str(), entry.kmer_groups.size());
    for (const auto& group : entry.kmer_groups) {
        if (group.t > 0)
            logger.info("  k=%d t=%d %s: %zu volume(s)", group.k, group.t,
                        template_type_to_string(static_cast<TemplateType>(group.template_type)).c_str(),
                        group.volumes.size());
        else
            logger.info("  k=%d: %zu volume(s)", group.k, group.volumes.size());
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
    pool_.shutdown();
}

int Server::try_acquire_sequences(int n) {
    std::lock_guard<std::mutex> lock(seq_mutex_);
    int capped = std::min(n, max_nseq_per_req_);
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

void Server::apply_madvise_budget(uint64_t budget, const Logger& logger) {
    uint64_t total = budget;
    auto try_willneed = [&budget](auto& reader) {
        uint64_t sz = reader.willneed_size();
        bool fits = (sz > 0 && budget >= sz);
        reader.apply_madvise(fits);
        if (fits) budget -= sz;
    };

    // Persistent madvise covers only khx + ksx (the small per-volume metadata
    // every search re-walks).  kix / kpx are opened/closed per-request so
    // concurrent requests do not contend on shared madvise calls; whatever
    // budget remains becomes the per-request posting_budget the search
    // orchestrator may spend on kix / kpx posting bodies inside one batch.

    // Priority 1: khx (one per k-mer group per DB)
    for (auto& db : databases_)
        for (auto& group : db.kmer_groups)
            try_willneed(group.khx);

    // Priority 2: ksx metadata
    for (auto& db : databases_)
        for (auto& group : db.kmer_groups)
            for (auto& vol : group.volumes)
                try_willneed(vol.ksx);

    posting_budget_ = budget;
    logger.info("madvise budget: %s used / %s total (per-request posting_budget=%s)",
                format_size(total - budget).c_str(),
                format_size(total).c_str(),
                format_size(posting_budget_).c_str());
}

void Server::accept_loop(int listen_fd, const ServerConfig& config, const Logger& logger) {
    int nthread = config.nthread;
    if (nthread <= 0) {
        nthread = static_cast<int>(std::thread::hardware_concurrency());
        if (nthread <= 0) nthread = 1;
    }

    tbb::task_arena arena(nthread);
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

    // Apply madvise budget
    {
        uint64_t budget = config.memory_limit > 0 ? config.memory_limit : default_memory_limit();
        apply_madvise_budget(budget, logger);
    }

    // Configure the inter-request posting-budget pool.  floor == 0 disables
    // blocking entirely (pass-through mode = historical behaviour: every
    // request sees the full posting_budget).  When -max_concurrent_search
    // is >= 1, residual posting_budget is divided into N slices floored at
    // 1 MiB so at most N concurrent requests run the budget-bound stages.
    {
        uint64_t floor = 0;
        if (config.max_concurrent_search > 0) {
            floor = std::max<uint64_t>(
                1ull << 20,
                posting_budget_ / static_cast<uint64_t>(config.max_concurrent_search));
            logger.info("BudgetPool: max_concurrent_search=%d, floor=%s",
                        config.max_concurrent_search, format_size(floor).c_str());
        }
        pool_.configure(posting_budget_, floor);
    }

    // Write PID file if requested
    if (!config.pid_file.empty()) {
        write_pid_file(config.pid_file);
    }

    // Resolve per-sequence concurrency limits
    {
        int nthread = config.nthread;
        if (nthread <= 0) {
            nthread = static_cast<int>(std::thread::hardware_concurrency());
            if (nthread <= 0) nthread = 1;
        }
        max_queue_size_ = config.max_queue_size > 0 ? config.max_queue_size : 1024;
        max_nseq_per_req_ = config.max_nseq_per_req > 0 ? config.max_nseq_per_req : nthread;
        logger.info("Max concurrent sequences: %d, max per request: %d",
                     max_queue_size_, max_nseq_per_req_);
    }

    // Log persistent mmap count (.khx + .ksx only; .kix / .kpx are
    // opened/closed per-request).
    size_t total_mmaps = 0;
    for (const auto& db : databases_) {
        for (const auto& group : db.kmer_groups) {
            total_mmaps += group.volumes.size();  // ksx
            if (group.khx.is_open()) total_mmaps++;
        }
    }
    logger.info("Persistent mmap'd files across %zu DB(s): %zu",
                databases_.size(), total_mmaps);

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
