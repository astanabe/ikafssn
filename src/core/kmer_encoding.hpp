#pragma once

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <functional>
#include <string>
#include <vector>
#include "core/config.hpp"
#include "core/degenerate_scan_simd.hpp"
#include "core/spaced_seed.hpp"
#include "core/spaced_seed_simd.hpp"

namespace ikafssn {

// 256-element LUT: char -> 2-bit encoding. 0xFF = invalid (N, etc.)
inline constexpr uint8_t BASE_ENCODE_INVALID = 0xFF;

inline const uint8_t* base_encode_table() {
    static const uint8_t table[256] = {
        // 0x00-0x3F
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        // 0x40-0x5F: @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_
        0xFF,0x00,0xFF,0x01,0xFF,0xFF,0xFF,0x02,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0x03,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        // 0x60-0x7F: `abcdefghijklmnopqrstuvwxyz{|}~DEL
        0xFF,0x00,0xFF,0x01,0xFF,0xFF,0xFF,0x02,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0x03,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        // 0x80-0xFF
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    };
    return table;
}

inline uint8_t encode_base(char c) {
    return base_encode_table()[static_cast<uint8_t>(c)];
}

// Reverse complement of a k-mer.
// Steps: ~kmer to complement all bits, then reverse 2-bit pairs, shift out unused bits.
template <typename KmerInt>
inline KmerInt kmer_revcomp(KmerInt kmer, int k) {
    constexpr int W = sizeof(KmerInt) * 8;
    KmerInt rc = ~kmer;
    // Reverse 2-bit pairs within the integer
    if constexpr (sizeof(KmerInt) == 2) {
        // uint16_t: swap bytes, then swap 2-bit pairs within bytes
        rc = (rc >> 8) | (rc << 8);
        rc = ((rc >> 4) & KmerInt(0x0F0F)) | ((rc & KmerInt(0x0F0F)) << 4);
        rc = ((rc >> 2) & KmerInt(0x3333)) | ((rc & KmerInt(0x3333)) << 2);
    } else {
        // uint32_t: swap bytes, then swap nibbles, then swap 2-bit pairs
        rc = ((rc >> 16) & 0x0000FFFF) | ((rc & 0x0000FFFF) << 16);
        rc = ((rc >>  8) & 0x00FF00FF) | ((rc & 0x00FF00FF) <<  8);
        rc = ((rc >>  4) & 0x0F0F0F0F) | ((rc & 0x0F0F0F0F) <<  4);
        rc = ((rc >>  2) & 0x33333333) | ((rc & 0x33333333) <<  2);
    }
    // Shift out unused high bits
    rc >>= (W - 2 * k);
    return rc;
}

// 256-element LUT: true for IUPAC ambiguity codes (R,Y,S,W,K,M,B,D,H,V,N)
inline const bool* degenerate_base_table() {
    static const bool table[256] = {
        // 0x00-0x3F: all false
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        // 0x40-0x5F: uppercase letters
        // @=0 A=0 B=1 C=0 D=1 E=0 F=0 G=0 H=1 I=0 J=0 K=1 L=0 M=1 N=1 O=0
        // P=0 Q=0 R=1 S=1 T=0 U=0 V=1 W=1 X=0 Y=1 Z=0 [=0 \=0 ]=0 ^=0 _=0
        0,0,1,0,1,0,0,0,1,0,0,1,0,1,1,0,
        0,0,1,1,0,0,1,1,0,1,0,0,0,0,0,0,
        // 0x60-0x7F: lowercase letters
        // `=0 a=0 b=1 c=0 d=1 e=0 f=0 g=0 h=1 i=0 j=0 k=1 l=0 m=1 n=1 o=0
        // p=0 q=0 r=1 s=1 t=0 u=0 v=1 w=1 x=0 y=1 z=0 {=0 |=0 }=0 ~=0 DEL=0
        0,0,1,0,1,0,0,0,1,0,0,1,0,1,1,0,
        0,0,1,1,0,0,1,1,0,1,0,0,0,0,0,0,
        // 0x80-0xFF: all false
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    };
    return table;
}

// Check if a sequence contains any IUPAC degenerate bases.
// SIMD-accelerated via has_degenerate_base().
inline bool contains_degenerate_base(const std::string& seq) {
    return has_degenerate_base(seq.data(), seq.size());
}

// 256-element LUT: IUPAC degenerate char -> ncbi4na bitmask (0 = not degenerate).
// Encoding: bit0=A(0x1), bit1=C(0x2), bit2=G(0x4), bit3=T(0x8)
// R=A|G=0x5, Y=C|T=0xA, S=G|C=0x6, W=A|T=0x9, K=G|T=0xC, M=A|C=0x3
// B=C|G|T=0xE, D=A|G|T=0xD, H=A|C|T=0xB, V=A|C|G=0x7, N=A|C|G|T=0xF
// Normal bases A/C/G/T return 0 (not degenerate).
inline const uint8_t* degenerate_ncbi4na_table() {
    static const uint8_t table[256] = {
        // 0x00-0x3F: all 0
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        // 0x40-0x5F: @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_
        // @=0 A=0 B=0xE C=0 D=0xD E=0 F=0 G=0 H=0xB I=0 J=0 K=0xC L=0 M=0x3 N=0xF O=0
        // P=0 Q=0 R=0x5 S=0x6 T=0 U=0 V=0x7 W=0x9 X=0 Y=0xA Z=0 [=0 \=0 ]=0 ^=0 _=0
        0,0,0x0E,0,0x0D,0,0,0,0x0B,0,0,0x0C,0,0x03,0x0F,0,
        0,0,0x05,0x06,0,0,0x07,0x09,0,0x0A,0,0,0,0,0,0,
        // 0x60-0x7F: `abcdefghijklmnopqrstuvwxyz{|}~DEL
        0,0,0x0E,0,0x0D,0,0,0,0x0B,0,0,0x0C,0,0x03,0x0F,0,
        0,0,0x05,0x06,0,0,0x07,0x09,0,0x0A,0,0,0,0,0,0,
        // 0x80-0xFF: all 0
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    };
    return table;
}

inline uint8_t degenerate_ncbi4na(char c) {
    return degenerate_ncbi4na_table()[static_cast<uint8_t>(c)];
}

// Info about one degenerate position in a k-mer window.
struct AmbigInfo {
    uint8_t ncbi4na;    // IUPAC bitmask for this position
    int bit_offset;     // 2-bit position in the k-mer integer
};

// Return the number of bases encoded in an ncbi4na bitmask (popcount of lower 4 bits).
inline int ncbi4na_expansion_count(uint8_t ncbi4na) {
    int count = 0;
    for (uint8_t b = 0; b < 4; b++) {
        if (ncbi4na & (1u << b)) count++;
    }
    return count;
}

// Expand a k-mer with multiple degenerate positions.
// Recursively iterates through each degenerate position, substituting each
// possible base, and calls action(expanded_kmer) for each combination.
template <typename KmerInt, typename Action>
inline void expand_ambig_kmer_multi(KmerInt base_kmer, const AmbigInfo* infos,
                                    int count, Action&& action) {
    if (count == 0) {
        action(base_kmer);
        return;
    }
    KmerInt clear_mask = ~(KmerInt(0x03) << infos[0].bit_offset);
    KmerInt cleared = base_kmer & clear_mask;
    for (uint8_t b = 0; b < 4; b++) {
        if (infos[0].ncbi4na & (1u << b)) {
            KmerInt variant = cleared | (KmerInt(b) << infos[0].bit_offset);
            expand_ambig_kmer_multi(variant, infos + 1, count - 1,
                                    std::forward<Action>(action));
        }
    }
}

// Expand a single ambiguous base in a k-mer and invoke action for each expansion.
// base_kmer: k-mer with placeholder ncbi2na value (0=A) at the ambiguous position
// ncbi4na: the ambiguity code bitmask (which bases it represents)
// bit_offset: 2-bit position in the k-mer integer of the ambiguous base
// action: called with each expanded KmerInt value
template <typename KmerInt, typename Action>
inline void expand_ambig_kmer(KmerInt base_kmer, uint8_t ncbi4na,
                              int bit_offset, Action&& action) {
    KmerInt clear_mask = ~(KmerInt(0x03) << bit_offset);
    KmerInt cleared = base_kmer & clear_mask;
    for (uint8_t b = 0; b < 4; b++) {
        if (ncbi4na & (1u << b)) {
            action(cleared | (KmerInt(b) << bit_offset));
        }
    }
}

// Sliding window k-mer scanner with N counter.
// Calls callback(pos, kmer) for each valid k-mer.
template <typename KmerInt>
class KmerScanner {
public:
    KmerScanner(int k) : k_(k), mask_(kmer_mask<KmerInt>(k)) {}

    // Scan a sequence of given length.
    // callback(uint32_t pos, KmerInt kmer) called for each valid k-mer.
    template <typename Callback>
    void scan(const char* seq, size_t len, Callback&& callback) const {
        if (static_cast<int>(len) < k_) return;

        KmerInt kmer = 0;
        int n_count = k_ - 1; // need k valid bases before first k-mer

        for (size_t i = 0; i < len; i++) {
            uint8_t enc = encode_base(seq[i]);
            if (enc == BASE_ENCODE_INVALID) {
                n_count = k_ - 1;
                kmer = 0;
                continue;
            }
            kmer = ((kmer << 2) | static_cast<KmerInt>(enc)) & mask_;
            if (n_count > 0) {
                n_count--;
                continue;
            }
            // Position of k-mer start = i - k + 1
            callback(static_cast<uint32_t>(i - k_ + 1), kmer);
        }
    }

    // Scan with degenerate base expansion.
    // callback(uint32_t pos, KmerInt kmer): normal k-mer (no degenerate bases)
    // ambig_callback(uint32_t pos, KmerInt base_kmer, const AmbigInfo* infos, int count):
    //   k-mer with degenerate bases whose expansion product <= max_expansion.
    //   Caller should use expand_ambig_kmer_multi().
    // K-mers whose expansion product exceeds max_expansion are skipped.
    // max_expansion: maximum expansion product (default 4). 0 or 1 disables all expansion.
    // has_multi_degen: if non-null, set to true when any k-mer is skipped due to exceeding max_expansion.
    template <typename Callback, typename AmbigCallback>
    void scan_ambig(const char* seq, size_t len,
                    Callback&& callback, AmbigCallback&& ambig_callback,
                    bool* has_multi_degen = nullptr,
                    int max_expansion = 4) const {
        if (static_cast<int>(len) < k_) return;

        const uint8_t* ncbi4na_tbl = degenerate_ncbi4na_table();

        KmerInt kmer = 0;
        int n_count = k_ - 1; // need k valid bases before first k-mer
        int degen_count = 0;   // number of degenerate bases in current window

        // Ring buffer tracking degenerate status per window slot.
        // window_degen[i % k] stores ncbi4na of the base at position i (0 = normal).
        uint8_t window_degen[MAX_K];
        std::memset(window_degen, 0, sizeof(window_degen));

        for (size_t i = 0; i < len; i++) {
            char ch = seq[i];
            uint8_t enc = encode_base(ch);
            uint8_t ncbi4na = ncbi4na_tbl[static_cast<uint8_t>(ch)];
            int slot = static_cast<int>(i % k_);

            if (enc == BASE_ENCODE_INVALID && ncbi4na == 0) {
                // Truly invalid character (not a degenerate base)
                n_count = k_ - 1;
                kmer = 0;
                degen_count = 0;
                std::memset(window_degen, 0, k_);
                continue;
            }

            // Evict the leaving slot (only when window is full, i.e. n_count == 0)
            if (n_count == 0) {
                if (window_degen[slot] != 0) {
                    degen_count--;
                }
            }

            if (ncbi4na != 0) {
                // Degenerate base: use placeholder encoding 0 (A)
                enc = 0;
                window_degen[slot] = ncbi4na;
                degen_count++;
            } else {
                // Normal base
                window_degen[slot] = 0;
            }

            kmer = ((kmer << 2) | static_cast<KmerInt>(enc)) & mask_;

            if (n_count > 0) {
                n_count--;
                continue;
            }

            uint32_t pos = static_cast<uint32_t>(i - k_ + 1);

            if (degen_count == 0) {
                callback(pos, kmer);
            } else if (max_expansion <= 1) {
                // Expansion disabled: skip all degenerate k-mers
                if (has_multi_degen) *has_multi_degen = true;
            } else {
                // Compute expansion product and collect AmbigInfo
                int product = 1;
                int info_count = 0;
                AmbigInfo infos[MAX_K];
                bool exceeded = false;
                for (int j = 0; j < k_; j++) {
                    int s = (static_cast<int>(i) + 1 + j) % k_;
                    if (window_degen[s] != 0) {
                        int ec = ncbi4na_expansion_count(window_degen[s]);
                        product *= ec;
                        if (product > max_expansion) {
                            exceeded = true;
                            break;
                        }
                        infos[info_count].ncbi4na = window_degen[s];
                        infos[info_count].bit_offset = (k_ - 1 - j) * 2;
                        info_count++;
                    }
                }
                if (!exceeded) {
                    ambig_callback(pos, kmer, infos, info_count);
                } else if (has_multi_degen) {
                    *has_multi_degen = true;
                }
            }
        }
    }

    // Scan with spaced seed templates.
    // Uses sliding accumulator + bitmask extraction for O(1) per-position
    // advance. Mask bit j (from LSB) corresponds to sequence position
    // p + (t-1-j). Per-position k-mer extraction is buffered in chunks of
    // kSpacedBatchChunk and processed via extract_for_mask_batch() at
    // flush time; callback emit order is per-position, then per-mask.
    template <typename Callback>
    void scan_spaced(const char* seq, size_t len, const std::vector<uint32_t>& masks,
                     int t, Callback&& callback) const {
        if (static_cast<int>(len) < t) return;
        const uint8_t* enc_tbl = base_encode_table();
        const int num_masks = static_cast<int>(masks.size());
        const uint64_t accum_mask = (1ULL << (2 * t)) - 1;
        const uint32_t n_window_mask = (1u << t) - 1;

        constexpr size_t kBatchChunk = 64;  // matches kNcbi2naUnpackChunkSize
        alignas(64) uint64_t accum_buf[kBatchChunk];
        alignas(64) uint32_t nbits_buf[kBatchChunk];
        alignas(64) uint32_t pos_buf[kBatchChunk];
        alignas(64) KmerInt kmer_buf[2][kBatchChunk];  // num_masks <= 2 (kBoth)
        size_t fill = 0;

        auto flush = [&]() {
            if (fill == 0) return;
            for (int mi = 0; mi < num_masks; mi++) {
                extract_for_mask_batch<KmerInt>(accum_buf, fill, masks[mi],
                                                kmer_buf[mi]);
            }
            for (size_t j = 0; j < fill; j++) {
                for (int mi = 0; mi < num_masks; mi++) {
                    if ((nbits_buf[j] & masks[mi]) == 0) {
                        callback(pos_buf[j], kmer_buf[mi][j]);
                    }
                }
            }
            fill = 0;
        };

        uint64_t accum = 0;
        uint32_t n_bits = 0;  // bit j set = window position j has invalid base
        int warmup = t - 1;

        for (size_t i = 0; i < len; i++) {
            uint8_t enc = enc_tbl[static_cast<uint8_t>(seq[i])];

            n_bits = (n_bits << 1) & n_window_mask;
            if (enc == BASE_ENCODE_INVALID) {
                n_bits |= 1u;
                enc = 0;  // placeholder
            }

            accum = ((accum << 2) | static_cast<uint64_t>(enc)) & accum_mask;

            if (warmup > 0) {
                warmup--;
                continue;
            }

            accum_buf[fill] = accum;
            nbits_buf[fill] = n_bits;
            pos_buf[fill] =
                static_cast<uint32_t>(i) - static_cast<uint32_t>(t) + 1;
            fill++;
            if (fill == kBatchChunk) flush();
        }
        flush();
    }

    // Scan with spaced seed templates, handling degenerate bases.
    // Uses sliding accumulator + dual bitsets (n_bits / degen_bits) for
    // efficient tracking. Pure-fast-path positions (no degenerate base on
    // any mask's set bits) are buffered and emitted via
    // extract_for_mask_batch() at flush time; positions that hit slow-path
    // on at least one mask drain the buffer first (so window_ncbi4na is
    // consistent at the slow position) and then take the scalar path.
    template <typename Callback, typename AmbigCallback>
    void scan_spaced_ambig(const char* seq, size_t len,
                            const std::vector<uint32_t>& masks, int t,
                            Callback&& callback, AmbigCallback&& ambig_callback,
                            bool* has_multi_degen = nullptr,
                            int max_expansion = 16) const {
        if (static_cast<int>(len) < t) return;
        const uint8_t* enc_tbl = base_encode_table();
        const uint8_t* ncbi4na_tbl = degenerate_ncbi4na_table();
        const int num_masks = static_cast<int>(masks.size());
        const uint64_t accum_mask = (1ULL << (2 * t)) - 1;
        const uint32_t window_mask = (1u << t) - 1;

        constexpr size_t kBatchChunk = 64;
        alignas(64) uint64_t accum_buf[kBatchChunk];
        alignas(64) uint32_t nbits_buf[kBatchChunk];
        alignas(64) uint32_t pos_buf[kBatchChunk];
        alignas(64) KmerInt kmer_buf[2][kBatchChunk];
        size_t fill = 0;

        auto flush = [&]() {
            if (fill == 0) return;
            for (int mi = 0; mi < num_masks; mi++) {
                extract_for_mask_batch<KmerInt>(accum_buf, fill, masks[mi],
                                                kmer_buf[mi]);
            }
            for (size_t j = 0; j < fill; j++) {
                for (int mi = 0; mi < num_masks; mi++) {
                    if ((nbits_buf[j] & masks[mi]) == 0) {
                        callback(pos_buf[j], kmer_buf[mi][j]);
                    }
                }
            }
            fill = 0;
        };

        uint64_t accum = 0;
        uint32_t n_bits = 0;      // truly invalid (N, not degenerate)
        uint32_t degen_bits = 0;   // degenerate IUPAC bases
        uint8_t window_ncbi4na[MAX_SPACED_T] = {};
        int warmup = t - 1;

        for (size_t i = 0; i < len; i++) {
            char ch = seq[i];
            uint8_t enc = enc_tbl[static_cast<uint8_t>(ch)];
            uint8_t ncbi4na = ncbi4na_tbl[static_cast<uint8_t>(ch)];

            n_bits = (n_bits << 1) & window_mask;
            degen_bits = (degen_bits << 1) & window_mask;

            if (enc == BASE_ENCODE_INVALID && ncbi4na == 0) {
                // Truly invalid (N)
                n_bits |= 1u;
                enc = 0;
            } else if (ncbi4na != 0) {
                // Degenerate IUPAC base
                degen_bits |= 1u;
                enc = 0;  // placeholder
            }

            window_ncbi4na[i % t] = ncbi4na;
            accum = ((accum << 2) | static_cast<uint64_t>(enc)) & accum_mask;

            if (warmup > 0) {
                warmup--;
                continue;
            }

            uint32_t p = static_cast<uint32_t>(i) - static_cast<uint32_t>(t) + 1;

            // If any mask hits the slow path at this position, drain the
            // buffer (so older fast-path callbacks fire in order) and then
            // process this position scalar — window_ncbi4na is read at the
            // current i, which would be stale for buffered positions.
            bool needs_slow = false;
            for (int mi = 0; mi < num_masks; mi++) {
                uint32_t mask = masks[mi];
                if ((n_bits & mask) != 0) continue;
                if ((degen_bits & mask) != 0) { needs_slow = true; break; }
            }

            if (needs_slow) {
                flush();
                for (int mi = 0; mi < num_masks; mi++) {
                    uint32_t mask = masks[mi];
                    if ((n_bits & mask) != 0) continue;

                    if ((degen_bits & mask) == 0) {
                        KmerInt kmer = extract_for_mask<KmerInt>(accum, mask);
                        callback(p, kmer);
                    } else if (max_expansion <= 1) {
                        if (has_multi_degen) *has_multi_degen = true;
                    } else {
                        KmerInt kmer = extract_for_mask<KmerInt>(accum, mask);
                        const auto& ext = extractor_for_mask(mask);
                        AmbigInfo infos[MAX_K];
                        int degen_count = 0;
                        int product = 1;
                        bool exceeded = false;

                        uint32_t affected = degen_bits & mask;
                        while (affected != 0) {
                            int j = std::countr_zero(affected);
                            affected &= affected - 1;

                            size_t seq_pos = static_cast<size_t>(p) +
                                static_cast<size_t>(t - 1 - j);
                            uint8_t a4na = window_ncbi4na[seq_pos % t];

                            int kbit_off = ext.kmer_bit_offset[j];
                            kmer &= ~(static_cast<KmerInt>(0x03) << kbit_off);

                            int ec = ncbi4na_expansion_count(a4na);
                            product *= ec;
                            if (product > max_expansion) {
                                exceeded = true;
                                break;
                            }

                            infos[degen_count].ncbi4na = a4na;
                            infos[degen_count].bit_offset = kbit_off;
                            degen_count++;
                        }

                        if (!exceeded && degen_count > 0) {
                            ambig_callback(p, kmer, infos, degen_count);
                        } else if (exceeded && has_multi_degen) {
                            *has_multi_degen = true;
                        }
                    }
                }
            } else {
                accum_buf[fill] = accum;
                nbits_buf[fill] = n_bits;
                pos_buf[fill]   = p;
                fill++;
                if (fill == kBatchChunk) flush();
            }
        }
        flush();
    }

private:
    int k_;
    KmerInt mask_;
};

} // namespace ikafssn
