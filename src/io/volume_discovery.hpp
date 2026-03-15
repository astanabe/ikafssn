#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

struct DiscoveredVolume {
    std::string kix_path;
    std::string kpx_path;
    std::string ksx_path;
    uint16_t volume_index;
    int k;
    bool has_kpx = true;
    uint8_t t = 0;              // template length (0=contiguous, 16/18/21=spaced)
    uint8_t template_type = 0;  // 0 = contiguous
};

// Split ix_prefix into parent directory and DB name.
// e.g. "path/to/nt" -> {"path/to", "nt"}
struct IndexPrefixParts {
    std::string parent_dir;
    std::string db;
};
IndexPrefixParts parse_index_prefix(const std::string& ix_prefix);

// Build the file stem for index files.
// e.g. index_file_stem("dir", "nt.00", 9) -> "dir/nt.00.09mer"
std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k);

// Build the file stem for spaced seed index files.
// t=0: "dir/vol.11mer" (delegates to contiguous overload)
// t>0: "dir/vol.11mer.16mer.both" (appends template info)
std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k,
                            uint8_t t, uint8_t template_type);

// Build the .khx file path.
// e.g. khx_path_for("dir", "nt", 9) -> "dir/nt.09mer.khx"
std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k);

// Build the .khx file path for spaced seed indexes.
std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k,
                          uint8_t t, uint8_t template_type);

// Discover index volumes from .kvx manifests.
// If filter_k > 0, only that k value. If filter_k == 0, all available k values.
// Results are sorted by (k, volume_index) ascending.
std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k = 0);

// Discover index volumes with spaced seed filter.
// If filter_t > 0, only that template length. If filter_template_type > 0, only that type.
std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k,
    uint8_t filter_t, uint8_t filter_template_type);

// Return available k values for the given index prefix.
std::vector<int> discover_k_values(const std::string& ix_prefix);

} // namespace ikafssn
