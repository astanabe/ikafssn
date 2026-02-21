#pragma once

#include <cstdint>
#include <cstddef>
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

private:
    int k_;
    KmerInt mask_;
};

} // namespace ikafssn
