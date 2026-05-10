#pragma once

#include <cstdint>
#include <bit>
#include <string>
#include <vector>
#include "core/config.hpp"

namespace ikafssn {

// Accumulator-based spaced seed extraction utilities.
//
// These enable O(1) per-position k-mer extraction from a sliding
// accumulator. Uses C++20 NTTP (non-type template parameter) for
// compile-time specialization: at -O3 the extraction loop is fully
// unrolled with immediate constants, producing assembly equivalent to
// hand-coded per-template functions (BLAST+ style).

inline constexpr int MAX_SPACED_T = 21;         // maximum template length
inline constexpr int MAX_EXTRACTION_RUNS = 8;    // max contiguous-1-bit runs across all 24 masks

// One contiguous run of set bits in a template mask.
// Describes how to extract a contiguous group of bases from the accumulator
// and place them into the k-mer integer.
struct ExtractionRun {
    uint8_t accum_shift;   // right-shift to apply to accumulator (= 2 * j_low)
    uint8_t kmer_shift;    // left-shift for k-mer position
    uint8_t width_bits;    // number of bits to extract (= 2 * run_length)
};

// Pre-computed extraction data for a single template mask.
// Structural type: usable as C++20 NTTP.
struct SpacedSeedExtractor {
    ExtractionRun runs[MAX_EXTRACTION_RUNS];
    int8_t kmer_bit_offset[MAX_SPACED_T];  // mask bit j -> k-mer bit offset (-1 if not set)
    uint8_t num_runs;
};

// Compute extraction runs and bit-offset table for a given mask.
// Scans set-bit runs from MSB (bit t-1) to LSB (bit 0).
// constexpr: can be evaluated at compile time for known masks.
constexpr SpacedSeedExtractor compute_spaced_seed_extractor(uint32_t mask, int t, int k) {
    SpacedSeedExtractor ext{};
    int kmer_bit = 0;
    int j = t - 1;
    while (j >= 0) {
        if (mask & (1u << j)) {
            // Start of a contiguous run of set bits
            int run_len = 0;
            while (j >= 0 && (mask & (1u << j))) {
                ext.kmer_bit_offset[j] = static_cast<int8_t>((k - 1 - kmer_bit) * 2);
                run_len++;
                kmer_bit++;
                j--;
            }
            int j_low = j + 1;
            ext.runs[ext.num_runs].accum_shift = static_cast<uint8_t>(2 * j_low);
            ext.runs[ext.num_runs].kmer_shift = static_cast<uint8_t>((k - kmer_bit) * 2);
            ext.runs[ext.num_runs].width_bits = static_cast<uint8_t>(2 * run_len);
            ext.num_runs++;
        } else {
            ext.kmer_bit_offset[j] = -1;
            j--;
        }
    }
    return ext;
}

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

// Discontiguous megablast template masks.
//
// Each mask is a uint32_t bitmask where bit j (bit 0 = rightmost)
// corresponds to position (t-1-j) in the window.  A set bit means the
// base at that position contributes to the k-mer value.  The number
// of set bits = k.

// Count set bits in a mask (constexpr popcount).
inline constexpr int popcount32(uint32_t v) {
    int c = 0;
    while (v) { c += (v & 1); v >>= 1; }
    return c;
}

// Weight 8 templates (PCR-optimized, k=8, derived from discontiguous MegaBLAST design principles):
inline constexpr uint32_t MASK_K8_T13_CODING   = 0x1B2D;   // 1101100101101
inline constexpr uint32_t MASK_K8_T13_OPTIMAL  = 0x1CD3;   // 1110011010011
inline constexpr uint32_t MASK_K8_T15_CODING   = 0x4B2D;   // 100101100101101
inline constexpr uint32_t MASK_K8_T15_OPTIMAL  = 0x7493;   // 111010010010011
inline constexpr uint32_t MASK_K8_T18_CODING   = 0x2492D;  // 100100100100101101
inline constexpr uint32_t MASK_K8_T18_OPTIMAL  = 0x32293;  // 110010001010010011

// Weight 9 templates (PCR-optimized, k=9, derived from discontiguous MegaBLAST design principles):
inline constexpr uint32_t MASK_K9_T13_CODING   = 0x1B6D;   // 1101101101101
inline constexpr uint32_t MASK_K9_T13_OPTIMAL  = 0x1DD3;   // 1110111010011
inline constexpr uint32_t MASK_K9_T15_CODING   = 0x4B6D;   // 100101101101101
inline constexpr uint32_t MASK_K9_T15_OPTIMAL  = 0x74D3;   // 111010011010011
inline constexpr uint32_t MASK_K9_T18_CODING   = 0x24B2D;  // 100100101100101101
inline constexpr uint32_t MASK_K9_T18_OPTIMAL  = 0x3A293;  // 111010001010010011

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

// Pre-computed constexpr extractors for all 24 template masks.

// K=8 extractors
inline constexpr SpacedSeedExtractor EXT_K8_T13_COD  = compute_spaced_seed_extractor(MASK_K8_T13_CODING, 13, 8);
inline constexpr SpacedSeedExtractor EXT_K8_T13_OPT  = compute_spaced_seed_extractor(MASK_K8_T13_OPTIMAL, 13, 8);
inline constexpr SpacedSeedExtractor EXT_K8_T15_COD  = compute_spaced_seed_extractor(MASK_K8_T15_CODING, 15, 8);
inline constexpr SpacedSeedExtractor EXT_K8_T15_OPT  = compute_spaced_seed_extractor(MASK_K8_T15_OPTIMAL, 15, 8);
inline constexpr SpacedSeedExtractor EXT_K8_T18_COD  = compute_spaced_seed_extractor(MASK_K8_T18_CODING, 18, 8);
inline constexpr SpacedSeedExtractor EXT_K8_T18_OPT  = compute_spaced_seed_extractor(MASK_K8_T18_OPTIMAL, 18, 8);

// K=9 extractors
inline constexpr SpacedSeedExtractor EXT_K9_T13_COD  = compute_spaced_seed_extractor(MASK_K9_T13_CODING, 13, 9);
inline constexpr SpacedSeedExtractor EXT_K9_T13_OPT  = compute_spaced_seed_extractor(MASK_K9_T13_OPTIMAL, 13, 9);
inline constexpr SpacedSeedExtractor EXT_K9_T15_COD  = compute_spaced_seed_extractor(MASK_K9_T15_CODING, 15, 9);
inline constexpr SpacedSeedExtractor EXT_K9_T15_OPT  = compute_spaced_seed_extractor(MASK_K9_T15_OPTIMAL, 15, 9);
inline constexpr SpacedSeedExtractor EXT_K9_T18_COD  = compute_spaced_seed_extractor(MASK_K9_T18_CODING, 18, 9);
inline constexpr SpacedSeedExtractor EXT_K9_T18_OPT  = compute_spaced_seed_extractor(MASK_K9_T18_OPTIMAL, 18, 9);

// K=11 extractors
inline constexpr SpacedSeedExtractor EXT_K11_T16_COD = compute_spaced_seed_extractor(MASK_K11_T16_CODING, 16, 11);
inline constexpr SpacedSeedExtractor EXT_K11_T16_OPT = compute_spaced_seed_extractor(MASK_K11_T16_OPTIMAL, 16, 11);
inline constexpr SpacedSeedExtractor EXT_K11_T18_COD = compute_spaced_seed_extractor(MASK_K11_T18_CODING, 18, 11);
inline constexpr SpacedSeedExtractor EXT_K11_T18_OPT = compute_spaced_seed_extractor(MASK_K11_T18_OPTIMAL, 18, 11);
inline constexpr SpacedSeedExtractor EXT_K11_T21_COD = compute_spaced_seed_extractor(MASK_K11_T21_CODING, 21, 11);
inline constexpr SpacedSeedExtractor EXT_K11_T21_OPT = compute_spaced_seed_extractor(MASK_K11_T21_OPTIMAL, 21, 11);

// K=12 extractors
inline constexpr SpacedSeedExtractor EXT_K12_T16_COD = compute_spaced_seed_extractor(MASK_K12_T16_CODING, 16, 12);
inline constexpr SpacedSeedExtractor EXT_K12_T16_OPT = compute_spaced_seed_extractor(MASK_K12_T16_OPTIMAL, 16, 12);
inline constexpr SpacedSeedExtractor EXT_K12_T18_COD = compute_spaced_seed_extractor(MASK_K12_T18_CODING, 18, 12);
inline constexpr SpacedSeedExtractor EXT_K12_T18_OPT = compute_spaced_seed_extractor(MASK_K12_T18_OPTIMAL, 18, 12);
inline constexpr SpacedSeedExtractor EXT_K12_T21_COD = compute_spaced_seed_extractor(MASK_K12_T21_CODING, 21, 12);
inline constexpr SpacedSeedExtractor EXT_K12_T21_OPT = compute_spaced_seed_extractor(MASK_K12_T21_OPTIMAL, 21, 12);

// NTTP extraction function: compile-time specialized per extractor.
// At -O3 the loop is fully unrolled with all shifts/masks as immediates.
template <SpacedSeedExtractor Ext, typename KmerInt>
inline KmerInt extract_kmer_ct(uint64_t accum) {
    KmerInt kmer = 0;
    for (int r = 0; r < Ext.num_runs; r++) {
        kmer |= static_cast<KmerInt>(
            (accum >> Ext.runs[r].accum_shift)
            & ((1ULL << Ext.runs[r].width_bits) - 1)
        ) << Ext.runs[r].kmer_shift;
    }
    return kmer;
}

// Runtime extraction function (for custom masks or testing).
template <typename KmerInt>
inline KmerInt extract_spaced_kmer(uint64_t accum, const SpacedSeedExtractor& ext) {
    KmerInt kmer = 0;
    for (int r = 0; r < ext.num_runs; r++) {
        kmer |= static_cast<KmerInt>(
            (accum >> ext.runs[r].accum_shift)
            & ((1ULL << ext.runs[r].width_bits) - 1)
        ) << ext.runs[r].kmer_shift;
    }
    return kmer;
}

// Switch dispatch: selects the NTTP-specialized extraction for a known
// mask.  The branch predictor handles this at ~0 cost since the mask
// is constant per scan.
template <typename KmerInt>
inline KmerInt extract_for_mask(uint64_t accum, uint32_t mask) {
    switch (mask) {
        case MASK_K8_T13_CODING:   return extract_kmer_ct<EXT_K8_T13_COD, KmerInt>(accum);
        case MASK_K8_T13_OPTIMAL:  return extract_kmer_ct<EXT_K8_T13_OPT, KmerInt>(accum);
        case MASK_K8_T15_CODING:   return extract_kmer_ct<EXT_K8_T15_COD, KmerInt>(accum);
        case MASK_K8_T15_OPTIMAL:  return extract_kmer_ct<EXT_K8_T15_OPT, KmerInt>(accum);
        case MASK_K8_T18_CODING:   return extract_kmer_ct<EXT_K8_T18_COD, KmerInt>(accum);
        case MASK_K8_T18_OPTIMAL:  return extract_kmer_ct<EXT_K8_T18_OPT, KmerInt>(accum);
        case MASK_K9_T13_CODING:   return extract_kmer_ct<EXT_K9_T13_COD, KmerInt>(accum);
        case MASK_K9_T13_OPTIMAL:  return extract_kmer_ct<EXT_K9_T13_OPT, KmerInt>(accum);
        case MASK_K9_T15_CODING:   return extract_kmer_ct<EXT_K9_T15_COD, KmerInt>(accum);
        case MASK_K9_T15_OPTIMAL:  return extract_kmer_ct<EXT_K9_T15_OPT, KmerInt>(accum);
        case MASK_K9_T18_CODING:   return extract_kmer_ct<EXT_K9_T18_COD, KmerInt>(accum);
        case MASK_K9_T18_OPTIMAL:  return extract_kmer_ct<EXT_K9_T18_OPT, KmerInt>(accum);
        case MASK_K11_T16_CODING:  return extract_kmer_ct<EXT_K11_T16_COD, KmerInt>(accum);
        case MASK_K11_T16_OPTIMAL: return extract_kmer_ct<EXT_K11_T16_OPT, KmerInt>(accum);
        case MASK_K11_T18_CODING:  return extract_kmer_ct<EXT_K11_T18_COD, KmerInt>(accum);
        case MASK_K11_T18_OPTIMAL: return extract_kmer_ct<EXT_K11_T18_OPT, KmerInt>(accum);
        case MASK_K11_T21_CODING:  return extract_kmer_ct<EXT_K11_T21_COD, KmerInt>(accum);
        case MASK_K11_T21_OPTIMAL: return extract_kmer_ct<EXT_K11_T21_OPT, KmerInt>(accum);
        case MASK_K12_T16_CODING:  return extract_kmer_ct<EXT_K12_T16_COD, KmerInt>(accum);
        case MASK_K12_T16_OPTIMAL: return extract_kmer_ct<EXT_K12_T16_OPT, KmerInt>(accum);
        case MASK_K12_T18_CODING:  return extract_kmer_ct<EXT_K12_T18_COD, KmerInt>(accum);
        case MASK_K12_T18_OPTIMAL: return extract_kmer_ct<EXT_K12_T18_OPT, KmerInt>(accum);
        case MASK_K12_T21_CODING:  return extract_kmer_ct<EXT_K12_T21_COD, KmerInt>(accum);
        case MASK_K12_T21_OPTIMAL: return extract_kmer_ct<EXT_K12_T21_OPT, KmerInt>(accum);
        default: {
            // Fallback for unknown masks (e.g., testing with custom masks)
            auto ext = compute_spaced_seed_extractor(mask, std::bit_width(mask), popcount32(mask));
            return extract_spaced_kmer<KmerInt>(accum, ext);
        }
    }
}

// Return the constexpr extractor for a known mask (for slow-path kmer_bit_offset access).
inline const SpacedSeedExtractor& extractor_for_mask(uint32_t mask) {
    switch (mask) {
        case MASK_K8_T13_CODING:   return EXT_K8_T13_COD;
        case MASK_K8_T13_OPTIMAL:  return EXT_K8_T13_OPT;
        case MASK_K8_T15_CODING:   return EXT_K8_T15_COD;
        case MASK_K8_T15_OPTIMAL:  return EXT_K8_T15_OPT;
        case MASK_K8_T18_CODING:   return EXT_K8_T18_COD;
        case MASK_K8_T18_OPTIMAL:  return EXT_K8_T18_OPT;
        case MASK_K9_T13_CODING:   return EXT_K9_T13_COD;
        case MASK_K9_T13_OPTIMAL:  return EXT_K9_T13_OPT;
        case MASK_K9_T15_CODING:   return EXT_K9_T15_COD;
        case MASK_K9_T15_OPTIMAL:  return EXT_K9_T15_OPT;
        case MASK_K9_T18_CODING:   return EXT_K9_T18_COD;
        case MASK_K9_T18_OPTIMAL:  return EXT_K9_T18_OPT;
        case MASK_K11_T16_CODING:  return EXT_K11_T16_COD;
        case MASK_K11_T16_OPTIMAL: return EXT_K11_T16_OPT;
        case MASK_K11_T18_CODING:  return EXT_K11_T18_COD;
        case MASK_K11_T18_OPTIMAL: return EXT_K11_T18_OPT;
        case MASK_K11_T21_CODING:  return EXT_K11_T21_COD;
        case MASK_K11_T21_OPTIMAL: return EXT_K11_T21_OPT;
        case MASK_K12_T16_CODING:  return EXT_K12_T16_COD;
        case MASK_K12_T16_OPTIMAL: return EXT_K12_T16_OPT;
        case MASK_K12_T18_CODING:  return EXT_K12_T18_COD;
        case MASK_K12_T18_OPTIMAL: return EXT_K12_T18_OPT;
        case MASK_K12_T21_CODING:  return EXT_K12_T21_COD;
        case MASK_K12_T21_OPTIMAL: return EXT_K12_T21_OPT;
        default: {
            // Should not reach here for known masks; return first as fallback.
            static const SpacedSeedExtractor fallback{};
            return fallback;
        }
    }
}

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
