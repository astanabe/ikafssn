#include "io/volume_discovery.hpp"
#include "io/kvx_reader.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <regex>
#include <set>

namespace ikafssn {

IndexPrefixParts parse_index_prefix(const std::string& ix_prefix) {
    std::filesystem::path p(ix_prefix);
    IndexPrefixParts parts;
    parts.parent_dir = p.parent_path().string();
    parts.db_name = p.filename().string();
    if (parts.parent_dir.empty()) parts.parent_dir = ".";
    return parts;
}

std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k) {
    char kk[8];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    return parent_dir + "/" + vol_basename + "." + kk + "mer";
}

std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db_name, int k) {
    return index_file_stem(parent_dir, db_name, k) + ".khx";
}

// Discover volumes for a single k value from a .kvx manifest.
static bool discover_from_kvx(const std::string& parent_dir,
                               const std::string& db_name,
                               int k,
                               std::vector<DiscoveredVolume>& volumes) {
    std::string stem = index_file_stem(parent_dir, db_name, k);
    std::string kvx_path = stem + ".kvx";

    auto kvx = read_kvx(kvx_path);
    if (!kvx) return false;

    bool found_any = false;
    for (uint16_t vi = 0; vi < kvx->volume_basenames.size(); vi++) {
        std::string base = index_file_stem(parent_dir, kvx->volume_basenames[vi], k);
        if (!std::filesystem::exists(base + ".kix")) continue;

        DiscoveredVolume dv;
        dv.volume_index = vi;
        dv.k = k;
        dv.kix_path = base + ".kix";
        dv.kpx_path = base + ".kpx";
        dv.ksx_path = base + ".ksx";
        volumes.push_back(std::move(dv));
        found_any = true;
    }
    return found_any;
}

// Scan directory for available k values matching the DB name.
static std::set<int> scan_k_values(const std::string& parent_dir,
                                    const std::string& db_name) {
    std::set<int> k_values;
    std::regex kvx_pattern("(\\d+)mer\\.kvx");
    std::string prefix_dot = db_name + ".";

    for (const auto& entry : std::filesystem::directory_iterator(parent_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        if (fname.size() <= prefix_dot.size() ||
            fname.compare(0, prefix_dot.size(), prefix_dot) != 0)
            continue;
        std::string suffix = fname.substr(prefix_dot.size());
        std::smatch m;
        if (std::regex_match(suffix, m, kvx_pattern)) {
            k_values.insert(std::stoi(m[1].str()));
        }
    }
    return k_values;
}

std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k) {
    auto parts = parse_index_prefix(ix_prefix);
    std::vector<DiscoveredVolume> volumes;

    if (filter_k > 0) {
        discover_from_kvx(parts.parent_dir, parts.db_name, filter_k, volumes);
    } else {
        for (int k : scan_k_values(parts.parent_dir, parts.db_name)) {
            discover_from_kvx(parts.parent_dir, parts.db_name, k, volumes);
        }
    }

    std::sort(volumes.begin(), volumes.end(),
              [](const DiscoveredVolume& a, const DiscoveredVolume& b) {
                  if (a.k != b.k) return a.k < b.k;
                  return a.volume_index < b.volume_index;
              });
    return volumes;
}

std::vector<int> discover_k_values(const std::string& ix_prefix) {
    auto parts = parse_index_prefix(ix_prefix);
    auto k_set = scan_k_values(parts.parent_dir, parts.db_name);
    return std::vector<int>(k_set.begin(), k_set.end());
}

} // namespace ikafssn
