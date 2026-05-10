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

    // Fragment-indexing parameters parsed from the file name.
    uint32_t min_seq_length    = 0;
    uint32_t min_length_split  = 0;
    uint32_t overlap_length    = 0;
    uint64_t max_freq_build    = 1;  // absolute threshold (1 = disabled)
    uint32_t max_degen_expand  = 0;  // 0 / 1 = disabled
};

// Parsed components of an index file name.  vol_basename is "<DB>.<vol>"
// for per-volume files (.kix / .kpx / .ksx) or "<DB>" for per-DB files
// (.kvx / .khx); has_vol distinguishes the two.
struct IndexFilenameParts {
    std::string parent_dir;
    std::string db_basename;
    std::string vol_basename;
    bool        has_vol = false;
    int         k = 0;
    uint8_t     t = 0;
    uint8_t     template_type = 0;
    uint32_t    min_seq_length = 0;
    uint32_t    min_length_split = 0;
    uint32_t    overlap_length = 0;
    uint64_t    max_freq_build = 1;   // absolute threshold (1 = disabled)
    uint32_t    max_degen_expand = 0; // 0 / 1 = disabled
};

// Split ix_prefix into parent directory and DB name.
// e.g. "path/to/nt" -> {"path/to", "nt"}
struct IndexPrefixParts {
    std::string parent_dir;
    std::string db;
};
IndexPrefixParts parse_index_prefix(const std::string& ix_prefix);

// Build the file stem for a per-volume index file (kix/kpx/ksx).
//
// Output layout (omitted suffixes per the rules below):
//   <parent_dir>/<vol_basename>.k<k>[.t<t>].minlen<X>.minsplit<X>.ovllen<X>
//                .maxfreq<X>.maxexpand<X>[.cod|.opt]
//
// Suffix-omission rules:
//   t<X>        : t == 0
//   minlen<X>   : min_seq_length == 0
//   minsplit<X> : min_length_split == 0
//   ovllen<X>   : overlap_length == 0
//   maxfreq<X>  : max_freq_build == 1 (absolute threshold)
//   maxexpand<X>: max_degen_expand == 0 or 1
//   cod|opt     : t == 0
std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k,
                            uint8_t t,
                            uint8_t template_type,
                            uint32_t min_seq_length,
                            uint32_t min_length_split,
                            uint32_t overlap_length,
                            uint64_t max_freq_build,
                            uint32_t max_degen_expand);

// Build the .khx file path (one .khx per database; the second argument
// is the DB basename, no per-volume index inserted).
std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k,
                          uint8_t t,
                          uint8_t template_type,
                          uint32_t min_seq_length,
                          uint32_t min_length_split,
                          uint32_t overlap_length,
                          uint64_t max_freq_build,
                          uint32_t max_degen_expand);

// Parse the components encoded in an index file name (kix / kpx / ksx /
// khx / kvx).  `path` may be a full filename or just the stem.  Returns
// false if the name does not match the expected pattern.
bool parse_index_filename(const std::string& path, IndexFilenameParts& out);

// Discover index volumes from .kvx manifests.
// If filter_k > 0, only that k value. If filter_k == 0, all available k values.
// If filter_t > 0, only that template length. If filter_template_type > 0, only that type.
// Results are sorted by (k, t, template_type, volume_index) ascending.
std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k = 0,
    uint8_t filter_t = 0, uint8_t filter_template_type = 0);

} // namespace ikafssn
