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
    parts.db = p.filename().string();
    if (parts.parent_dir.empty()) parts.parent_dir = ".";
    return parts;
}

std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k) {
    char kk[8];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    return parent_dir + "/" + vol_basename + "." + kk + "mer";
}

std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k,
                            uint8_t t, uint8_t template_type) {
    std::string stem = index_file_stem(parent_dir, vol_basename, k);
    if (t > 0) {
        char tt[8];
        std::snprintf(tt, sizeof(tt), "%02d", static_cast<int>(t));
        std::string type_str;
        switch (template_type) {
            case 1:  type_str = "cod"; break;
            case 2:  type_str = "opt"; break;
            default: type_str = "con"; break;
        }
        stem += "." + std::string(tt) + "mer." + type_str;
    }
    return stem;
}

std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k) {
    return index_file_stem(parent_dir, db, k) + ".khx";
}

std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k,
                          uint8_t t, uint8_t template_type) {
    return index_file_stem(parent_dir, db, k, t, template_type) + ".khx";
}

// Discover volumes for a single (k, t, template_type) from a .kvx manifest.
static bool discover_from_kvx(const std::string& parent_dir,
                               const std::string& db,
                               int k,
                               uint8_t t, uint8_t template_type,
                               std::vector<DiscoveredVolume>& volumes) {
    std::string stem = index_file_stem(parent_dir, db, k, t, template_type);
    std::string kvx_path = stem + ".kvx";

    auto kvx = read_kvx(kvx_path);
    if (!kvx) return false;

    bool found_any = false;
    for (uint16_t vi = 0; vi < kvx->volume_basenames.size(); vi++) {
        std::string base = index_file_stem(parent_dir, kvx->volume_basenames[vi], k,
                                           t, template_type);
        if (!std::filesystem::exists(base + ".kix")) continue;

        DiscoveredVolume dv;
        dv.volume_index = vi;
        dv.k = k;
        dv.t = t;
        dv.template_type = template_type;
        dv.kix_path = base + ".kix";
        dv.kpx_path = base + ".kpx";
        dv.ksx_path = base + ".ksx";
        dv.has_kpx = std::filesystem::exists(dv.kpx_path);
        volumes.push_back(std::move(dv));
        found_any = true;
    }
    return found_any;
}

// Extended kvx scan result for both contiguous and spaced seed indexes.
struct KvxScanResult {
    int k;
    uint8_t t;
    uint8_t template_type;
    bool operator<(const KvxScanResult& o) const {
        if (k != o.k) return k < o.k;
        if (t != o.t) return t < o.t;
        return template_type < o.template_type;
    }
};

// Scan directory for available (k, t, template_type) combinations.
// Detects both contiguous (XXmer.kvx) and spaced (XXmer.YYmer.ZZZ.kvx) patterns.
static std::set<KvxScanResult> scan_k_values_ext(const std::string& parent_dir,
                                                   const std::string& db) {
    std::set<KvxScanResult> results;
    std::regex contiguous_pattern("(\\d+)mer\\.kvx");
    std::regex spaced_pattern("(\\d+)mer\\.(\\d+)mer\\.(cod|opt)\\.kvx");
    std::string prefix_dot = db + ".";

    for (const auto& entry : std::filesystem::directory_iterator(parent_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        if (fname.size() <= prefix_dot.size() ||
            fname.compare(0, prefix_dot.size(), prefix_dot) != 0)
            continue;
        std::string suffix = fname.substr(prefix_dot.size());
        std::smatch m;
        if (std::regex_match(suffix, m, spaced_pattern)) {
            KvxScanResult r;
            r.k = std::stoi(m[1].str());
            r.t = static_cast<uint8_t>(std::stoi(m[2].str()));
            std::string type_str = m[3].str();
            if (type_str == "cod") r.template_type = 1;
            else r.template_type = 2; // opt
            results.insert(r);
        } else if (std::regex_match(suffix, m, contiguous_pattern)) {
            results.insert({std::stoi(m[1].str()), 0, 0});
        }
    }
    return results;
}

// Backward-compatible scan returning only k values (contiguous indexes only).
static std::set<int> scan_k_values(const std::string& parent_dir,
                                    const std::string& db) {
    std::set<int> k_values;
    for (const auto& r : scan_k_values_ext(parent_dir, db)) {
        if (r.t == 0) {
            k_values.insert(r.k);
        }
    }
    return k_values;
}

std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k) {
    return discover_volumes(ix_prefix, filter_k, 0, 0);
}

std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k,
    uint8_t filter_t, uint8_t filter_template_type) {
    auto parts = parse_index_prefix(ix_prefix);
    std::vector<DiscoveredVolume> volumes;

    if (filter_k > 0) {
        discover_from_kvx(parts.parent_dir, parts.db, filter_k,
                          filter_t, filter_template_type, volumes);
    } else {
        for (const auto& r : scan_k_values_ext(parts.parent_dir, parts.db)) {
            if (filter_t > 0 && r.t != filter_t) continue;
            if (filter_template_type > 0 && r.template_type != filter_template_type) continue;
            discover_from_kvx(parts.parent_dir, parts.db, r.k,
                              r.t, r.template_type, volumes);
        }
    }

    std::sort(volumes.begin(), volumes.end(),
              [](const DiscoveredVolume& a, const DiscoveredVolume& b) {
                  if (a.k != b.k) return a.k < b.k;
                  if (a.t != b.t) return a.t < b.t;
                  if (a.template_type != b.template_type)
                      return a.template_type < b.template_type;
                  return a.volume_index < b.volume_index;
              });
    return volumes;
}

std::vector<int> discover_k_values(const std::string& ix_prefix) {
    auto parts = parse_index_prefix(ix_prefix);
    auto k_set = scan_k_values(parts.parent_dir, parts.db);
    return std::vector<int>(k_set.begin(), k_set.end());
}

} // namespace ikafssn
