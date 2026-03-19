#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "core/config.hpp"

namespace ikafssn {

// Template type for spaced seed / discontiguous megablast support.
// kContiguous (0) = traditional contiguous k-mer (t=0).
// kCoding (1) = coding template (optimized for coding regions).
// kOptimal (2) = optimal template (optimized for non-coding regions).
// kBoth (3) = merged index from both coding and optimal templates.
enum class TemplateType : uint8_t {
    kContiguous = 0,
    kCoding     = 1,
    kOptimal    = 2,
    kBoth       = 3,
};

// Configuration for spaced seed indexing/searching.
struct SpacedSeedConfig {
    uint8_t t = 0;  // template length (0 = contiguous, 13/15/18 for k=8-9, 16/18/21 for k=11-12)
    TemplateType type = TemplateType::kBoth;
};

// ============================================================================
// Discontiguous megablast template masks.
//
// Each mask is a uint32_t bitmask where bit j (counting from bit 0 = rightmost)
// corresponds to position (t-1-j) in the window. A set bit means the base at
// that position contributes to the k-mer value. The number of set bits = k.
// ============================================================================

// Count set bits in a mask (constexpr popcount).
inline constexpr int popcount32(uint32_t v) {
    int c = 0;
    while (v) { c += (v & 1); v >>= 1; }
    return c;
}

// Weight 8 templates (PCR-optimized, k=8):
inline constexpr uint32_t MASK_K8_T13_CODING   = 0x165B;   // 1011001011011
inline constexpr uint32_t MASK_K8_T13_OPTIMAL  = 0x196D;   // 1100101101101
inline constexpr uint32_t MASK_K8_T15_CODING   = 0x592D;   // 101100100101101
inline constexpr uint32_t MASK_K8_T15_OPTIMAL  = 0x6965;   // 110100101100101
inline constexpr uint32_t MASK_K8_T18_CODING   = 0x25925;  // 100101100100100101
inline constexpr uint32_t MASK_K8_T18_OPTIMAL  = 0x34A25;  // 110100101000100101

// Weight 9 templates (PCR-optimized, k=9):
inline constexpr uint32_t MASK_K9_T13_CODING   = 0x1B6D;   // 1101101101101
inline constexpr uint32_t MASK_K9_T13_OPTIMAL  = 0x16DB;   // 1011011011011
inline constexpr uint32_t MASK_K9_T15_CODING   = 0x5B2D;   // 101101100101101
inline constexpr uint32_t MASK_K9_T15_OPTIMAL  = 0x6B2D;   // 110101100101101
inline constexpr uint32_t MASK_K9_T18_CODING   = 0x25965;  // 100101100101100101
inline constexpr uint32_t MASK_K9_T18_OPTIMAL  = 0x34B25;  // 110100101100100101

// Weight 11 templates:
inline constexpr uint32_t MASK_K11_T16_CODING  = 0xDB6D;   // 1101101101101101
inline constexpr uint32_t MASK_K11_T16_OPTIMAL = 0xE5B7;   // 1110010110110111
inline constexpr uint32_t MASK_K11_T18_CODING  = 0x2D96D;  // 101101100101101101
inline constexpr uint32_t MASK_K11_T18_OPTIMAL = 0x3A597;  // 111010010110010111
inline constexpr uint32_t MASK_K11_T21_CODING  = 0x12CB2D; // 100101100101100101101
inline constexpr uint32_t MASK_K11_T21_OPTIMAL = 0x1D2897; // 111010010100010010111

// Weight 12 templates:
inline constexpr uint32_t MASK_K12_T16_CODING  = 0xFB6D;   // 1111101101101101
inline constexpr uint32_t MASK_K12_T16_OPTIMAL = 0xEDB7;   // 1110110110110111
inline constexpr uint32_t MASK_K12_T18_CODING  = 0x2DB6D;  // 101101101101101101
inline constexpr uint32_t MASK_K12_T18_OPTIMAL = 0x3ACB7;  // 111010110010110111
inline constexpr uint32_t MASK_K12_T21_CODING  = 0x12DB2D; // 100101101101100101101
inline constexpr uint32_t MASK_K12_T21_OPTIMAL = 0x1D2C97; // 111010010110010010111

// Compile-time validation of mask weights.
static_assert(popcount32(MASK_K8_T13_CODING)   == 8,  "K8 T13 coding mask weight");
static_assert(popcount32(MASK_K8_T13_OPTIMAL)  == 8,  "K8 T13 optimal mask weight");
static_assert(popcount32(MASK_K8_T15_CODING)   == 8,  "K8 T15 coding mask weight");
static_assert(popcount32(MASK_K8_T15_OPTIMAL)  == 8,  "K8 T15 optimal mask weight");
static_assert(popcount32(MASK_K8_T18_CODING)   == 8,  "K8 T18 coding mask weight");
static_assert(popcount32(MASK_K8_T18_OPTIMAL)  == 8,  "K8 T18 optimal mask weight");
static_assert(popcount32(MASK_K9_T13_CODING)   == 9,  "K9 T13 coding mask weight");
static_assert(popcount32(MASK_K9_T13_OPTIMAL)  == 9,  "K9 T13 optimal mask weight");
static_assert(popcount32(MASK_K9_T15_CODING)   == 9,  "K9 T15 coding mask weight");
static_assert(popcount32(MASK_K9_T15_OPTIMAL)  == 9,  "K9 T15 optimal mask weight");
static_assert(popcount32(MASK_K9_T18_CODING)   == 9,  "K9 T18 coding mask weight");
static_assert(popcount32(MASK_K9_T18_OPTIMAL)  == 9,  "K9 T18 optimal mask weight");
static_assert(popcount32(MASK_K11_T16_CODING)  == 11, "K11 T16 coding mask weight");
static_assert(popcount32(MASK_K11_T16_OPTIMAL) == 11, "K11 T16 optimal mask weight");
static_assert(popcount32(MASK_K11_T18_CODING)  == 11, "K11 T18 coding mask weight");
static_assert(popcount32(MASK_K11_T18_OPTIMAL) == 11, "K11 T18 optimal mask weight");
static_assert(popcount32(MASK_K11_T21_CODING)  == 11, "K11 T21 coding mask weight");
static_assert(popcount32(MASK_K11_T21_OPTIMAL) == 11, "K11 T21 optimal mask weight");
static_assert(popcount32(MASK_K12_T16_CODING)  == 12, "K12 T16 coding mask weight");
static_assert(popcount32(MASK_K12_T16_OPTIMAL) == 12, "K12 T16 optimal mask weight");
static_assert(popcount32(MASK_K12_T18_CODING)  == 12, "K12 T18 coding mask weight");
static_assert(popcount32(MASK_K12_T18_OPTIMAL) == 12, "K12 T18 optimal mask weight");
static_assert(popcount32(MASK_K12_T21_CODING)  == 12, "K12 T21 coding mask weight");
static_assert(popcount32(MASK_K12_T21_OPTIMAL) == 12, "K12 T21 optimal mask weight");

// Get seed masks for a given (k, t, type) combination.
// Returns 1 mask for kCoding or kOptimal, 2 masks for kBoth.
// Caller should pass type != kContiguous and valid (k, t) combination.
inline std::vector<uint32_t> get_seed_masks(int k, uint8_t t, TemplateType type) {
    std::vector<uint32_t> masks;
    auto push = [&](uint32_t coding, uint32_t optimal) {
        if (type == TemplateType::kCoding || type == TemplateType::kBoth)
            masks.push_back(coding);
        if (type == TemplateType::kOptimal || type == TemplateType::kBoth)
            masks.push_back(optimal);
    };

    if (k == 8) {
        switch (t) {
            case 13: push(MASK_K8_T13_CODING, MASK_K8_T13_OPTIMAL); break;
            case 15: push(MASK_K8_T15_CODING, MASK_K8_T15_OPTIMAL); break;
            case 18: push(MASK_K8_T18_CODING, MASK_K8_T18_OPTIMAL); break;
            default: break;
        }
    } else if (k == 9) {
        switch (t) {
            case 13: push(MASK_K9_T13_CODING, MASK_K9_T13_OPTIMAL); break;
            case 15: push(MASK_K9_T15_CODING, MASK_K9_T15_OPTIMAL); break;
            case 18: push(MASK_K9_T18_CODING, MASK_K9_T18_OPTIMAL); break;
            default: break;
        }
    } else if (k == 11) {
        switch (t) {
            case 16: push(MASK_K11_T16_CODING, MASK_K11_T16_OPTIMAL); break;
            case 18: push(MASK_K11_T18_CODING, MASK_K11_T18_OPTIMAL); break;
            case 21: push(MASK_K11_T21_CODING, MASK_K11_T21_OPTIMAL); break;
            default: break;
        }
    } else if (k == 12) {
        switch (t) {
            case 16: push(MASK_K12_T16_CODING, MASK_K12_T16_OPTIMAL); break;
            case 18: push(MASK_K12_T18_CODING, MASK_K12_T18_OPTIMAL); break;
            case 21: push(MASK_K12_T21_CODING, MASK_K12_T21_OPTIMAL); break;
            default: break;
        }
    }
    return masks;
}

// Compute the effective table size for spaced seed indexes.
// When multiple masks are used (both mode), the table doubles to accommodate
// a 1-bit mask tag in the k-mer value that prevents cross-template matching.
inline uint32_t spaced_table_size(int k, int num_masks) {
    uint32_t base = table_size(k);
    return (num_masks > 1) ? static_cast<uint32_t>(num_masks) * base : base;
}

// Build masks and their corresponding tag values for searching/building.
// index_tt: template_type of the INDEX being used (0-3)
// search_tt: template_type requested for the SEARCH (0-3)
//
// For a "both" (3) index, each mask is tagged so coding and optimal k-mers
// occupy disjoint ranges. This allows a bot index to serve coding-only
// or optimal-only searches: the query k-mers only match the correct portion.
//
// Returns: (masks, tags) pair where tags[i] is the pre-shifted value to OR
// into k-mers from masks[i].  tags is empty when no tagging is needed.
template <typename KmerInt>
inline std::pair<std::vector<uint32_t>, std::vector<KmerInt>>
get_tagged_masks(int k, uint8_t t, uint8_t index_tt, uint8_t search_tt) {
    if (t == 0) return {{}, {}};

    auto coding_masks = get_seed_masks(k, t, TemplateType::kCoding);
    auto optimal_masks = get_seed_masks(k, t, TemplateType::kOptimal);
    uint32_t coding_mask = coding_masks.empty() ? 0 : coding_masks[0];
    uint32_t optimal_mask = optimal_masks.empty() ? 0 : optimal_masks[0];

    const KmerInt tag_optimal = static_cast<KmerInt>(1) << (2 * k);

    std::vector<uint32_t> masks;
    std::vector<KmerInt> tags;

    if (index_tt == 3) {
        // Index is "both" — tagged k-mers, table = 2*4^k
        TemplateType st = (search_tt != 0)
            ? static_cast<TemplateType>(search_tt)
            : TemplateType::kBoth;
        if (st == TemplateType::kCoding || st == TemplateType::kBoth) {
            masks.push_back(coding_mask);
            tags.push_back(KmerInt(0));  // coding tag = 0
        }
        if (st == TemplateType::kOptimal || st == TemplateType::kBoth) {
            masks.push_back(optimal_mask);
            tags.push_back(tag_optimal); // optimal tag = 1 << 2k
        }
    } else {
        // Index is single-template (coding or optimal) — no tagging
        masks = get_seed_masks(k, t, static_cast<TemplateType>(index_tt));
        // tags stays empty → no tagging
    }
    return {masks, tags};
}

// Validate spaced seed parameters.
// Returns true if valid (t=0 is always valid; t>0 requires valid (k, t) combination).
// k=8,9: t in {13,15,18}; k=11,12: t in {16,18,21}.
inline bool validate_spaced_seed(int k, uint8_t t) {
    if (t == 0) return true;
    if (k == 8 || k == 9) return (t == 13 || t == 15 || t == 18);
    if (k == 11 || k == 12) return (t == 16 || t == 18 || t == 21);
    return false;
}

// Return the effective span of the k-mer window.
// For spaced seeds (t > 0), the span is t. For contiguous k-mers, the span is k.
inline int seed_span(uint8_t t, int k) {
    return (t > 0) ? static_cast<int>(t) : k;
}

// Parse template type from CLI string.
// Returns kContiguous on unrecognized input (caller should validate).
inline TemplateType template_type_from_string(const std::string& s) {
    if (s == "coding") return TemplateType::kCoding;
    if (s == "optimal") return TemplateType::kOptimal;
    if (s == "coding_and_optimal" || s == "both") return TemplateType::kBoth;
    return TemplateType::kContiguous;
}

// Convert template type to string for file naming and display.
inline std::string template_type_to_string(TemplateType type) {
    switch (type) {
        case TemplateType::kCoding:  return "coding";
        case TemplateType::kOptimal: return "optimal";
        case TemplateType::kBoth:    return "both";
        default:                     return "contiguous";
    }
}

// Reverse complement an entire sequence string (char-level).
// Handles standard bases and IUPAC degenerate codes.
inline const char* rc_complement_table() {
    static char table[256];
    static bool init = false;
    if (!init) {
        for (int i = 0; i < 256; i++) table[i] = 'N';
        table['A'] = 'T'; table['T'] = 'A'; table['C'] = 'G'; table['G'] = 'C';
        table['a'] = 't'; table['t'] = 'a'; table['c'] = 'g'; table['g'] = 'c';
        table['R'] = 'Y'; table['Y'] = 'R'; table['S'] = 'S'; table['W'] = 'W';
        table['K'] = 'M'; table['M'] = 'K'; table['B'] = 'V'; table['V'] = 'B';
        table['D'] = 'H'; table['H'] = 'D'; table['N'] = 'N';
        table['r'] = 'y'; table['y'] = 'r'; table['s'] = 's'; table['w'] = 'w';
        table['k'] = 'm'; table['m'] = 'k'; table['b'] = 'v'; table['v'] = 'b';
        table['d'] = 'h'; table['h'] = 'd'; table['n'] = 'n';
        init = true;
    }
    return table;
}

inline std::string reverse_complement_string(const std::string& seq) {
    const char* tbl = rc_complement_table();
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); i++) {
        rc[i] = tbl[static_cast<uint8_t>(seq[seq.size() - 1 - i])];
    }
    return rc;
}

} // namespace ikafssn
