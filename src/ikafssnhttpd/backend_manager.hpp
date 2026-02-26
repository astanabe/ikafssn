#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <json/json.h>

#include "ikafssnhttpd/backend_client.hpp"
#include "protocol/messages.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// A single backend server entry.
struct BackendEntry {
    int priority;                              // CLI order (0 = highest)
    BackendMode mode;
    std::string address;
    std::shared_ptr<BackendClient> client;
    InfoResponse cached_info;
    bool info_valid = false;
    enum class Status { kHealthy, kExcluded };
    Status status = Status::kHealthy;
    std::chrono::steady_clock::time_point exclusion_expiry;
    mutable std::mutex mutex;
};

// Manages multiple ikafssnserver backends for ikafssnhttpd.
// Provides health-checking, routing, and aggregated info.
class BackendManager {
public:
    // Add a backend in CLI order. priority = order of add_backend() calls.
    void add_backend(BackendMode mode, const std::string& address);

    // Initialize: connect to all backends, fetch info, validate cross-server DBs.
    // Returns false if all backends fail or cross-server validation fails.
    bool init(int timeout_seconds, const Logger& logger);

    // Start background heartbeat thread. Periodically refreshes info from all backends.
    void start_heartbeat(int interval_seconds, const Logger& logger);

    // Stop background heartbeat thread.
    void stop_heartbeat();

    // Route a search request to the best available backend.
    // Returns false on routing/search failure.
    bool route_search(const SearchRequest& req, SearchResponse& resp,
                      std::string& error_msg);

    // Return a merged InfoResponse (DB/k/mode structure only) for validate_info().
    // Capacity fields (max_active_sequences/active_sequences) are set to 0.
    InfoResponse merged_info() const;

    // Build aggregated info JSON for the /info endpoint.
    Json::Value build_info_json() const;

    // Returns true if at least one backend is healthy.
    bool any_healthy() const;

    // Check health of a specific backend (synchronous). Used by health endpoint.
    bool check_any_health(std::string& error_msg) const;

    // Set exclusion time in seconds.
    void set_exclusion_time(int seconds) { exclusion_seconds_ = seconds; }

private:
    std::vector<std::unique_ptr<BackendEntry>> backends_;

    // Capability map: (db_name, k, mode) -> list of backend indices
    struct CapKey {
        std::string db_name;
        uint8_t k;
        uint8_t mode;
        bool operator==(const CapKey& o) const {
            return db_name == o.db_name && k == o.k && mode == o.mode;
        }
    };
    struct CapKeyHash {
        size_t operator()(const CapKey& k) const {
            size_t h = std::hash<std::string>()(k.db_name);
            h ^= std::hash<uint8_t>()(k.k) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<uint8_t>()(k.mode) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_map<CapKey, std::vector<size_t>, CapKeyHash> capability_map_;

    int exclusion_seconds_ = 3600;

    // Heartbeat thread
    std::thread heartbeat_thread_;
    std::atomic<bool> heartbeat_stop_{false};
    std::mutex heartbeat_mutex_;
    std::condition_variable heartbeat_cv_;

    // Refresh a single backend's info. Returns true on success.
    bool refresh_info(size_t idx, const Logger& logger);

    // Validate that same-named DBs across backends are consistent.
    bool validate_cross_server_dbs(const Logger& logger) const;

    // Rebuild capability map from all backends' cached info.
    void rebuild_capability_map();

    // Select the best backend for a given (db_name, k, mode) query.
    // Returns backend index, or -1 if none available.
    int select_backend(const std::string& db_name, uint8_t k, uint8_t mode) const;

    // Exclude a backend for exclusion_seconds_.
    void exclude_backend(size_t idx, const Logger& logger);
};

} // namespace ikafssn
