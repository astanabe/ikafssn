#pragma once

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <functional>
#include <string>
#include "core/config.hpp"

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

// Check if a sequence contains any IUPAC degenerate bases
inline bool contains_degenerate_base(const std::string& seq) {
    const bool* tbl = degenerate_base_table();
    for (char c : seq) {
        if (tbl[static_cast<uint8_t>(c)]) return true;
    }
    return false;
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

    // Scan with degenerate base expansion (1-degenerate only).
    // callback(uint32_t pos, KmerInt kmer): normal k-mer (no degenerate bases)
    // ambig_callback(uint32_t pos, KmerInt base_kmer, uint8_t ncbi4na, int bit_offset):
    //   k-mer with exactly one degenerate base. Caller should use expand_ambig_kmer().
    // K-mers with 2+ degenerate bases are skipped.
    template <typename Callback, typename AmbigCallback>
    void scan_ambig(const char* seq, size_t len,
                    Callback&& callback, AmbigCallback&& ambig_callback) const {
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
            } else if (degen_count == 1) {
                // Find the single degenerate slot and compute bit_offset.
                // Window covers positions [i-k+1 .. i].
                // Iterate j=0..k-1: slot = (i+1+j) % k, bit_offset = (k-1-j)*2
                for (int j = 0; j < k_; j++) {
                    int s = (static_cast<int>(i) + 1 + j) % k_;
                    if (window_degen[s] != 0) {
                        int bit_offset = (k_ - 1 - j) * 2;
                        ambig_callback(pos, kmer, window_degen[s], bit_offset);
                        break;
                    }
                }
            }
            // degen_count >= 2: skip
        }
    }

private:
    int k_;
    KmerInt mask_;
};

} // namespace ikafssn
