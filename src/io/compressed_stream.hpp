#pragma once

#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <streambuf>
#include <string>

namespace ikafssn {

enum class CompressionFormat {
    kNone,
    kGzip,
    kBzip2,
    kXz,
    kZstd,
};

// CLI sentinel meaning "use the codec's default level".
constexpr int kCompressionLevelDefault = -1;

// Detect a compression format from a path's trailing suffix.
// Recognises ".gz" (gzip), ".bz2" (bzip2), ".xz" (xz), and ".zst" (zstd).
// Suffix matching is case-insensitive.  Returns kNone for any other suffix
// or for "-" / empty paths.
CompressionFormat detect_format_from_extension(const std::string& path);

// Detect a compression format from the leading bytes of an input stream.
// `n` is the number of valid bytes in `prefix`; values smaller than 6 are
// tolerated.  Returns kNone when no codec magic matches.
CompressionFormat detect_format_from_magic(const unsigned char* prefix,
                                           std::size_t n);

// Validate `level` against the codec's accepted range.  `level` may be
// kCompressionLevelDefault, in which case the function returns true and
// leaves the codec to translate the sentinel to its own default at init
// time.  On failure the function returns false and sets `error_msg` to a
// human-readable diagnostic such as
//   "Error: -compression_level for gzip must be in 0..9 (got 12)".
bool validate_compression_level(CompressionFormat fmt, int level,
                                std::string& error_msg);

// Owned input stream wrapper.  When `sb` is null the stream delegates to
// std::cin's rdbuf (used for "-" / empty paths); otherwise `sb` owns the
// codec / file source and is destroyed alongside `stream`.
struct OwnedIstream {
    std::unique_ptr<std::streambuf> sb;
    std::unique_ptr<std::istream>   stream;

    explicit operator bool() const { return stream != nullptr; }
};

struct OwnedOstream {
    std::unique_ptr<std::streambuf> sb;
    std::unique_ptr<std::ostream>   stream;

    explicit operator bool() const { return stream != nullptr; }
};

// Open an input stream that auto-detects the codec from the leading bytes
// of the source.  `path` may be "-" (stdin) or a regular file path.  On
// failure returns an empty OwnedIstream and sets `error_msg`.
OwnedIstream open_input_compressed(const std::string& path,
                                    std::string& error_msg);

// Open an output stream whose codec is selected by `path`'s trailing
// suffix.  Empty path or "-" routes to stdout uncompressed.
// `level` is the requested compression level or kCompressionLevelDefault
// to fall back to the codec's default.  The caller is expected to have
// already passed `level` through validate_compression_level().
OwnedOstream open_output_compressed(const std::string& path, int level,
                                     std::string& error_msg);

} // namespace ikafssn
