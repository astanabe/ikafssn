#pragma once

// Elias-Fano dictionary codec, public API.
//
// The encoder writes the on-disk EF blob described in ef_format.hpp.
// EFDictionary opens an mmap'd blob and provides O(1)-amortised
// random access via select1 + low-bits read.
//
// The body lives in ef_codec_tier.cpp, compiled per ISA tier
// (sse42 / avx2 / avx512bw / avx512vbmi2 on x86_64; neon on aarch64).

#include "index/ef_format.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace ikafssn::ef {

// Encode the dictionary into Elias-Fano form.
//
//   offsets   D non-strictly monotonic byte offsets in [0, U_raw]
//   D         number of entries (e.g. 4^k or 4^k + 1 for sentinel-
//             terminated .kix dictionaries)
//   U_raw     exclusive upper bound on raw offset values; typically
//             posting_file_.size() or offsets[D-1] + 1.  The encoder
//             does NOT scan offsets to derive this value (writer side
//             already holds posting_file_.size()).
//   out       appended-to output buffer.  On return ``out`` contains
//             [EFHeader][lower bits][upper bits][select samples].
//
// Returns the total number of bytes written into ``out`` (including
// any pre-existing ``out`` content is not counted; this is the byte
// length of just the EF blob).
std::size_t encode_dictionary_ef(const std::uint64_t* offsets,
                                 std::size_t D,
                                 std::uint64_t U_raw,
                                 std::vector<std::uint8_t>& out);

// EF dictionary reader.  Open against an mmap'd blob; access(i) returns
// the i-th raw offset.  All accessors are const + thread-safe (no
// internal mutable state).
class EFDictionary {
public:
    EFDictionary() = default;

    // Open an EF blob.  ``data`` must remain valid for the lifetime of
    // this object (typically backed by mmap).  Returns false on magic
    // mismatch or truncated blob.
    bool open(const std::uint8_t* data, std::size_t bytes);

    // Number of entries (D from the header).
    std::size_t size() const noexcept { return D_; }

    // Total byte length of the EF blob (header + lower + upper +
    // select samples).  Used by ``apply_madvise(true)`` budgeting.
    std::size_t blob_bytes() const noexcept { return blob_bytes_; }

    // Equivalent to blob_bytes(); kept for API parity with KixReader /
    // KpxReader's existing ``willneed_size()``.
    std::size_t willneed_size() const noexcept { return blob_bytes_; }

    // posix_madvise(MADV_WILLNEED) when willneed=true, MADV_RANDOM
    // otherwise.  No-op on platforms without the advice flag.
    void apply_madvise(bool willneed) const noexcept;

    // Hot-path random access.  Returns the i-th raw offset (after the
    // strict-monotonic +i transform is undone).  i must be < size();
    // out-of-range access returns 0 in release builds.
    std::uint64_t access(std::uint32_t i) const noexcept;

    // Convenience for ``access(i)`` and ``access(i+1)`` together.  The
    // common pattern in KixReader::posting_byte_length is one such pair.
    // SIMD-aware tiers may fuse the two select1 calls into one.
    void access_pair(std::uint32_t i,
                     std::uint64_t& start,
                     std::uint64_t& end) const noexcept;

private:
    EFHeader header_{};
    const std::uint8_t* data_ = nullptr;
    std::size_t blob_bytes_ = 0;

    // Decoded view into the blob (data_-relative pointers).
    const std::uint64_t* lower_ = nullptr;
    const std::uint64_t* upper_ = nullptr;
    const std::uint64_t* select_ = nullptr;

    // Cached header fields for the hot path.
    std::size_t D_ = 0;
    std::uint8_t l_ = 0;
    std::uint64_t mask_l_ = 0;
    std::uint64_t upper_bits_total_ = 0;
    std::uint32_t select_count_ = 0;
};

// Name of the active EF codec tier ("sse42" / "avx2" / "avx512bw" /
// "avx512vbmi2" on x86_64; "neon" on aarch64).  First call resolves
// the tier; subsequent calls are cached.  Useful for diagnostic logging.
const char* active_tier_name();

} // namespace ikafssn::ef
