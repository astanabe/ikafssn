#include "ikafssnhttpd/backend_manager.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <set>
#include <thread>

namespace ikafssn {

void BackendManager::add_backend(BackendMode mode, const std::string& address) {
    auto entry = std::make_unique<BackendEntry>();
    entry->priority = static_cast<int>(backends_.size());
    entry->mode = mode;
    entry->address = address;
    entry->client = std::make_shared<BackendClient>(mode, address);
    backends_.push_back(std::move(entry));
}

bool BackendManager::init(int timeout_seconds, const Logger& logger) {
    if (backends_.empty()) {
        logger.error("No backends configured");
        return false;
    }

    // Try to connect to each backend with exponential backoff
    auto start = std::chrono::steady_clock::now();
    bool any_success = false;

    for (size_t i = 0; i < backends_.size(); i++) {
        auto& be = *backends_[i];
        int delay = 1;
        bool connected = false;

        while (true) {
            InfoResponse info;
            std::string err;
            if (be.client->info(info, err)) {
                std::lock_guard<std::mutex> lock(be.mutex);
                be.cached_info = std::move(info);
                be.info_valid = true;
                be.status = BackendEntry::Status::kHealthy;
                connected = true;
                logger.info("Backend %zu (%s) connected successfully",
                            i, be.address.c_str());
                break;
            }

            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start).count();
            if (elapsed >= timeout_seconds) {
                logger.error("Backend %zu (%s) failed to connect: %s",
                             i, be.address.c_str(), err.c_str());
                break;
            }

            std::this_thread::sleep_for(std::chrono::seconds(delay));
            delay = std::min(delay * 2, 8);
        }

        if (connected) any_success = true;
    }

    if (!any_success) {
        logger.error("All backends failed to connect");
        return false;
    }

    // Validate cross-server DB consistency
    if (!validate_cross_server_dbs(logger)) {
        return false;
    }

    // Build capability map
    rebuild_capability_map();

    return true;
}

bool BackendManager::validate_cross_server_dbs(const Logger& logger) const {
    // Collect per-DB aggregate info across backends
    struct DbStats {
        std::set<uint8_t> k_values;
        uint64_t total_sequences = 0;
        uint64_t total_bases = 0;
    };

    // Map: db_name -> backend_index -> DbStats
    std::unordered_map<std::string, std::vector<std::pair<size_t, DbStats>>> db_backends;

    for (size_t i = 0; i < backends_.size(); i++) {
        auto& be = *backends_[i];
        std::lock_guard<std::mutex> lock(be.mutex);
        if (!be.info_valid) continue;

        for (const auto& db : be.cached_info.databases) {
            DbStats stats;
            for (const auto& g : db.groups) {
                stats.k_values.insert(g.k);
                for (const auto& v : g.volumes) {
                    stats.total_sequences += v.num_sequences;
                    stats.total_bases += v.total_bases;
                }
            }
            db_backends[db.name].push_back({i, stats});
        }
    }

    // Check consistency for DBs that appear on multiple backends
    for (const auto& [db_name, entries] : db_backends) {
        if (entries.size() <= 1) continue;

        const auto& ref = entries[0].second;
        for (size_t j = 1; j < entries.size(); j++) {
            const auto& other = entries[j].second;

            if (ref.k_values != other.k_values) {
                logger.error("Cross-server validation failed for DB '%s': "
                             "k-value sets differ between backend %zu and %zu",
                             db_name.c_str(), entries[0].first, entries[j].first);
                return false;
            }

            if (ref.total_sequences != other.total_sequences) {
                logger.error("Cross-server validation failed for DB '%s': "
                             "total sequences differ between backend %zu (%lu) and %zu (%lu)",
                             db_name.c_str(), entries[0].first,
                             static_cast<unsigned long>(ref.total_sequences),
                             entries[j].first,
                             static_cast<unsigned long>(other.total_sequences));
                return false;
            }

            if (ref.total_bases != other.total_bases) {
                logger.error("Cross-server validation failed for DB '%s': "
                             "total bases differ between backend %zu (%lu) and %zu (%lu)",
                             db_name.c_str(), entries[0].first,
                             static_cast<unsigned long>(ref.total_bases),
                             entries[j].first,
                             static_cast<unsigned long>(other.total_bases));
                return false;
            }
        }
    }

    return true;
}

void BackendManager::rebuild_capability_map() {
    capability_map_.clear();

    for (size_t i = 0; i < backends_.size(); i++) {
        auto& be = *backends_[i];
        std::lock_guard<std::mutex> lock(be.mutex);
        if (!be.info_valid) continue;

        for (const auto& db : be.cached_info.databases) {
            for (const auto& g : db.groups) {
                for (uint8_t m = 1; m <= db.max_mode; m++) {
                    CapKey key{db.name, g.k, m};
                    capability_map_[key].push_back(i);
                }
                // Also add k=0 (server default) entries
                CapKey key0{db.name, 0, 0};
                auto& vec = capability_map_[key0];
                if (std::find(vec.begin(), vec.end(), i) == vec.end()) {
                    vec.push_back(i);
                }
            }
        }
    }
}

int BackendManager::select_backend(const std::string& db_name,
                                    uint8_t k, uint8_t mode) const {
    // Try exact match first
    CapKey key{db_name, k, mode};
    auto it = capability_map_.find(key);

    // Fall back to wildcard (k=0, mode=0) if exact not found
    if (it == capability_map_.end()) {
        CapKey key0{db_name, 0, 0};
        it = capability_map_.find(key0);
    }

    if (it == capability_map_.end()) {
        return -1;
    }

    const auto& candidates = it->second;
    auto now = std::chrono::steady_clock::now();

    // Filter out excluded backends, sort by priority
    int best_idx = -1;
    int best_priority = INT_MAX;
    bool best_has_capacity = false;

    for (size_t ci : candidates) {
        auto& be = *backends_[ci];
        std::lock_guard<std::mutex> lock(be.mutex);

        // Skip excluded backends (unless expired)
        if (be.status == BackendEntry::Status::kExcluded) {
            if (now < be.exclusion_expiry) {
                continue;
            }
            // Exclusion expired; treat as healthy for selection
        }

        if (!be.info_valid) continue;

        // Check capacity: consider both slot availability and per-request cap
        int32_t available = be.cached_info.max_active_sequences - be.cached_info.active_sequences;
        int32_t per_req = be.cached_info.max_seqs_per_req;
        int32_t effective = std::min(std::max(static_cast<int32_t>(0), available),
                                     per_req > 0 ? per_req : INT32_MAX);
        bool has_capacity = true;
        if (be.cached_info.max_active_sequences > 0 && effective <= 0) {
            has_capacity = false;
        }

        // Prefer: has_capacity > no_capacity, then lower priority
        if (best_idx < 0 ||
            (has_capacity && !best_has_capacity) ||
            (has_capacity == best_has_capacity && be.priority < best_priority)) {
            best_idx = static_cast<int>(ci);
            best_priority = be.priority;
            best_has_capacity = has_capacity;
        }
    }

    return best_idx;
}

bool BackendManager::refresh_info(size_t idx, const Logger& logger) {
    auto& be = *backends_[idx];
    InfoResponse info;
    std::string err;
    if (!be.client->info(info, err)) {
        logger.debug("Failed to refresh info for backend %zu (%s): %s",
                     idx, be.address.c_str(), err.c_str());
        return false;
    }
    std::lock_guard<std::mutex> lock(be.mutex);
    be.cached_info = std::move(info);
    be.info_valid = true;
    return true;
}

void BackendManager::exclude_backend(size_t idx, const Logger& logger) {
    auto& be = *backends_[idx];
    std::lock_guard<std::mutex> lock(be.mutex);
    be.status = BackendEntry::Status::kExcluded;
    be.exclusion_expiry = std::chrono::steady_clock::now() +
                          std::chrono::seconds(exclusion_seconds_);
    logger.info("Backend %zu (%s) excluded for %d seconds",
                idx, be.address.c_str(), exclusion_seconds_);
}

bool BackendManager::route_search(const SearchRequest& req, SearchResponse& resp,
                                   std::string& error_msg) {
    // Try up to 3 times with re-selection on failure
    for (int attempt = 0; attempt < 3; attempt++) {
        int idx = select_backend(req.db_name, req.k, req.mode);
        if (idx < 0) {
            error_msg = "No available backend for db=" + req.db_name;
            return false;
        }

        auto& be = *backends_[static_cast<size_t>(idx)];

        // Pre-check: refresh info
        if (!refresh_info(static_cast<size_t>(idx), Logger(Logger::kError))) {
            exclude_backend(static_cast<size_t>(idx), Logger(Logger::kError));
            continue;
        }

        // Perform search
        std::string search_err;
        if (be.client->search(req, resp, search_err)) {
            return true;
        }

        // Search failed - exclude and retry
        error_msg = search_err;
        exclude_backend(static_cast<size_t>(idx), Logger(Logger::kError));
    }

    return false;
}

InfoResponse BackendManager::merged_info() const {
    InfoResponse merged;
    merged.status = 0;
    merged.max_active_sequences = 0;
    merged.active_sequences = 0;

    // Merge databases from all healthy backends
    std::unordered_map<std::string, size_t> db_index;

    for (const auto& be : backends_) {
        std::lock_guard<std::mutex> lock(be->mutex);
        if (!be->info_valid) continue;
        if (be->status == BackendEntry::Status::kExcluded) {
            auto now = std::chrono::steady_clock::now();
            if (now < be->exclusion_expiry) continue;
        }

        for (const auto& db : be->cached_info.databases) {
            auto it = db_index.find(db.name);
            if (it == db_index.end()) {
                // New DB - add it
                db_index[db.name] = merged.databases.size();
                merged.databases.push_back(db);
                if (merged.default_k == 0) {
                    merged.default_k = db.default_k;
                }
            }
            // Same DB already exists - no need to merge (validated identical)
        }
    }

    return merged;
}

Json::Value BackendManager::build_info_json() const {
    Json::Value result;
    result["status"] = "success";

    // Collect and merge databases
    struct ModeCapInfo {
        int64_t sum_max_active = 0;
        int64_t sum_active = 0;
        int64_t sum_effective_per_req = 0; // sum(min(available_i, per_req_i))
    };

    struct MergedGroup {
        uint8_t k;
        uint8_t kmer_type;
        std::vector<VolumeInfo> volumes;
        std::map<uint8_t, ModeCapInfo> mode_capacity;
    };

    struct MergedDb {
        std::string name;
        uint8_t default_k;
        uint8_t max_mode;
        std::map<uint8_t, MergedGroup> groups; // key=k
    };

    std::map<std::string, MergedDb> merged_dbs;

    for (const auto& be : backends_) {
        std::lock_guard<std::mutex> lock(be->mutex);
        if (!be->info_valid) continue;
        if (be->status == BackendEntry::Status::kExcluded) {
            auto now = std::chrono::steady_clock::now();
            if (now < be->exclusion_expiry) continue;
        }

        for (const auto& db : be->cached_info.databases) {
            auto& mdb = merged_dbs[db.name];
            if (mdb.name.empty()) {
                mdb.name = db.name;
                mdb.default_k = db.default_k;
                mdb.max_mode = db.max_mode;
            }

            for (const auto& g : db.groups) {
                auto& mg = mdb.groups[g.k];
                if (mg.k == 0) {
                    mg.k = g.k;
                    mg.kmer_type = g.kmer_type;
                    mg.volumes = g.volumes;
                }

                // Accumulate capacity per mode
                for (uint8_t m = 1; m <= db.max_mode; m++) {
                    auto& cap = mg.mode_capacity[m];
                    cap.sum_max_active += be->cached_info.max_active_sequences;
                    cap.sum_active += be->cached_info.active_sequences;
                    int32_t avail = be->cached_info.max_active_sequences
                                  - be->cached_info.active_sequences;
                    int32_t pr = be->cached_info.max_seqs_per_req;
                    cap.sum_effective_per_req += std::min(
                        std::max(static_cast<int32_t>(0), avail),
                        pr > 0 ? pr : std::max(static_cast<int32_t>(0), avail));
                }
            }
        }
    }

    Json::Value databases_arr(Json::arrayValue);
    for (const auto& [name, mdb] : merged_dbs) {
        Json::Value dbobj;
        dbobj["name"] = mdb.name;
        dbobj["default_k"] = mdb.default_k;
        dbobj["max_mode"] = mdb.max_mode;

        Json::Value groups_arr(Json::arrayValue);
        for (const auto& [k, mg] : mdb.groups) {
            Json::Value gobj;
            gobj["k"] = mg.k;
            gobj["kmer_type"] = (mg.kmer_type == 0) ? "uint16" : "uint32";

            uint64_t group_total_sequences = 0;
            uint64_t group_total_bases = 0;
            uint64_t group_total_postings = 0;

            Json::Value vols_arr(Json::arrayValue);
            for (const auto& v : mg.volumes) {
                Json::Value vobj;
                vobj["volume_index"] = v.volume_index;
                vobj["num_sequences"] = v.num_sequences;
                vobj["total_postings"] = static_cast<Json::UInt64>(v.total_postings);
                vobj["total_bases"] = static_cast<Json::UInt64>(v.total_bases);
                vobj["db_name"] = v.db_name;
                vols_arr.append(std::move(vobj));

                group_total_sequences += v.num_sequences;
                group_total_bases += v.total_bases;
                group_total_postings += v.total_postings;
            }
            gobj["volumes"] = std::move(vols_arr);
            gobj["num_volumes"] = static_cast<Json::UInt>(mg.volumes.size());
            gobj["total_sequences"] = static_cast<Json::UInt64>(group_total_sequences);
            gobj["total_bases"] = static_cast<Json::UInt64>(group_total_bases);
            gobj["total_postings"] = static_cast<Json::UInt64>(group_total_postings);

            // Modes array
            Json::Value modes_arr(Json::arrayValue);
            for (const auto& [m, cap] : mg.mode_capacity) {
                Json::Value mobj;
                mobj["mode"] = m;
                mobj["max_active_sequences"] = static_cast<Json::Int64>(cap.sum_max_active);
                mobj["active_sequences"] = static_cast<Json::Int64>(cap.sum_active);
                mobj["max_seqs_per_req"] = static_cast<Json::Int64>(cap.sum_effective_per_req);
                modes_arr.append(std::move(mobj));
            }
            gobj["modes"] = std::move(modes_arr);

            groups_arr.append(std::move(gobj));
        }
        dbobj["kmer_groups"] = std::move(groups_arr);
        databases_arr.append(std::move(dbobj));
    }
    // Top-level max_seqs_per_req: minimum across all modes of all databases
    int64_t global_max_seqs_per_req = 0;
    bool has_any_mode = false;
    for (const auto& [name, mdb] : merged_dbs) {
        for (const auto& [k, mg] : mdb.groups) {
            for (const auto& [m, cap] : mg.mode_capacity) {
                if (!has_any_mode) {
                    global_max_seqs_per_req = cap.sum_effective_per_req;
                    has_any_mode = true;
                } else {
                    global_max_seqs_per_req = std::min(global_max_seqs_per_req,
                                                       cap.sum_effective_per_req);
                }
            }
        }
    }
    result["max_seqs_per_req"] = static_cast<Json::Int64>(global_max_seqs_per_req);
    result["databases"] = std::move(databases_arr);

    return result;
}

bool BackendManager::any_healthy() const {
    auto now = std::chrono::steady_clock::now();
    for (const auto& be : backends_) {
        std::lock_guard<std::mutex> lock(be->mutex);
        if (be->status == BackendEntry::Status::kHealthy) return true;
        if (be->status == BackendEntry::Status::kExcluded &&
            now >= be->exclusion_expiry) return true;
    }
    return false;
}

bool BackendManager::check_any_health(std::string& error_msg) const {
    for (const auto& be : backends_) {
        HealthResponse hresp;
        std::string err;
        if (be->client->health_check(hresp, err)) {
            return true;
        }
    }
    error_msg = "All backends are unreachable";
    return false;
}

void BackendManager::start_heartbeat(int interval_seconds, const Logger& logger) {
    heartbeat_stop_.store(false);

    heartbeat_thread_ = std::thread([this, interval_seconds, &logger] {
        while (!heartbeat_stop_.load()) {
            {
                std::unique_lock<std::mutex> lock(heartbeat_mutex_);
                heartbeat_cv_.wait_for(lock, std::chrono::seconds(interval_seconds),
                                       [this] { return heartbeat_stop_.load(); });
            }

            if (heartbeat_stop_.load()) break;

            auto now = std::chrono::steady_clock::now();
            for (size_t i = 0; i < backends_.size(); i++) {
                auto& be = *backends_[i];

                // Check if exclusion has expired
                {
                    std::lock_guard<std::mutex> lock(be.mutex);
                    if (be.status == BackendEntry::Status::kExcluded &&
                        now >= be.exclusion_expiry) {
                        logger.info("Backend %zu (%s) exclusion expired, re-checking",
                                    i, be.address.c_str());
                    }
                }

                if (refresh_info(i, logger)) {
                    std::lock_guard<std::mutex> lock(be.mutex);
                    if (be.status == BackendEntry::Status::kExcluded &&
                        now >= be.exclusion_expiry) {
                        be.status = BackendEntry::Status::kHealthy;
                        logger.info("Backend %zu (%s) re-enabled after exclusion",
                                    i, be.address.c_str());
                    }
                }
            }

            // Rebuild capability map after refresh
            rebuild_capability_map();
        }
    });
}

void BackendManager::stop_heartbeat() {
    heartbeat_stop_.store(true);
    heartbeat_cv_.notify_all();
    if (heartbeat_thread_.joinable()) {
        heartbeat_thread_.join();
    }
}

} // namespace ikafssn
