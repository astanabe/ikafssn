#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

class Logger;

// Magic / format_version / file-size check against the two-stage
// parent / fragment table layout.  Accepts either a final ".ksx" or
// a ".ksx.tmp" path.
bool validate_ksx_file_strict(const std::string& ksx_path,
                              const Logger& logger);

// kix posting-file size check + full validate_volume posting list
// walk.  Empty kpx_path verifies the kix side alone (mode-1 indexes).
bool validate_kix_kpx_strict(const std::string& kix_path,
                             const std::string& kpx_path,
                             const Logger& logger);

// Magic / format_version / (k, t, template_type) tuple / file-size
// check against header + ceil(4^k / 8).
bool validate_khx_file_strict(const std::string& khx_path,
                              int k, uint8_t t, uint8_t template_type,
                              const Logger& logger);

// Strict validation of one volume's final files (prefix + .ksx /
// .kix [/ .kpx]).  skip_kpx omits the .kpx side (mode-1).
bool validate_volume_final_strict(const std::string& prefix,
                                  bool skip_kpx,
                                  const Logger& logger);

// Same as validate_volume_final_strict but for the .ksx.tmp /
// .kix.tmp / .kpx.tmp files of one volume.
bool validate_volume_tmp_strict(const std::string& prefix,
                                bool skip_kpx,
                                const Logger& logger);

} // namespace ikafssn
