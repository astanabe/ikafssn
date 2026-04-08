#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include "core/config.hpp"
#include "core/ambiguity_parser.hpp"
#include "core/kmer_encoding.hpp"
#include "core/spaced_seed.hpp"

namespace ikafssn {

// Extracts 2-bit base code from ncbi2na packed data at given base position.
// ncbi2na packing: MSB-first, 4 bases per byte.
// Byte layout: bits 7-6 = base0, bits 5-4 = base1, bits 3-2 = base2, bits 1-0 = base3
inline uint8_t ncbi2na_base_at(const char* data, uint32_t pos) {
    uint8_t byte = static_cast<uint8_t>(data[pos >> 2]);
    return (byte >> (6 - 2 * (pos & 3))) & 0x03;
}

// Sliding window k-mer scanner that reads directly from ncbi2na packed data.
// Tracks ambiguous bases using the ambiguity entry list.
//
// callback(uint32_t pos, KmerInt kmer) - called for normal k-mers (no ambiguity)
// ambig_callback(uint32_t pos, KmerInt base_kmer, const AmbigInfo* infos, int count)
//   - called for k-mers with degenerate bases whose expansion product <= max_expansion
//   - infos: array of AmbigInfo describing each degenerate position
//   - count: number of degenerate positions
//
// K-mers whose expansion product exceeds max_expansion are skipped.
template <typename KmerInt>
class PackedKmerScanner {
public:
    PackedKmerScanner(int k) : k_(k), mask_(kmer_mask<KmerInt>(k)) {}

    template <typename Callback, typename AmbigCallback>
    void scan(const char* ncbi2na_data, uint32_t seq_length,
              const std::vector<AmbiguityEntry>& ambig_entries,
              Callback&& callback,
              AmbigCallback&& ambig_callback,
              int max_expansion = 4) const {

        if (static_cast<int>(seq_length) < k_) return;

        // Cursor into the ambiguity entry list for tracking which
        // individual ambiguous base positions enter/leave the k-mer window.
        // Each cursor tracks (entry_index, offset_within_run).
        struct Cursor {
            size_t entry_idx = 0;
            uint32_t run_offset = 0;

            uint32_t pos(const std::vector<AmbiguityEntry>& e) const {
                if (entry_idx >= e.size()) return UINT32_MAX;
                return e[entry_idx].position + run_offset;
            }

            uint8_t ncbi4na(const std::vector<AmbiguityEntry>& e) const {
                return e[entry_idx].ncbi4na;
            }

            void advance(const std::vector<AmbiguityEntry>& e) {
                if (entry_idx >= e.size()) return;
                run_offset++;
                if (run_offset >= e[entry_idx].run_length) {
                    entry_idx++;
                    run_offset = 0;
                }
            }
        };

        Cursor enter_cur;  // tracks right edge (entering positions)
        Cursor leave_cur;  // tracks left edge (leaving positions)
        int ambig_count = 0;

        // Track the single ambiguous base info when ambig_count == 1
        uint32_t single_pos = UINT32_MAX;
        uint8_t single_ncbi4na = 0;

        KmerInt kmer = 0;

        // Fill initial window: positions [0, k-1)
        for (int i = 0; i < k_ - 1; i++) {
            uint8_t code = ncbi2na_base_at(ncbi2na_data, static_cast<uint32_t>(i));
            kmer = ((kmer << 2) | static_cast<KmerInt>(code)) & mask_;

            if (enter_cur.pos(ambig_entries) == static_cast<uint32_t>(i)) {
                ambig_count++;
                single_pos = static_cast<uint32_t>(i);
                single_ncbi4na = enter_cur.ncbi4na(ambig_entries);
                enter_cur.advance(ambig_entries);
            }
        }

        // Main loop: i is the rightmost position of the k-mer window
        for (uint32_t i = static_cast<uint32_t>(k_ - 1); i < seq_length; i++) {
            uint8_t code = ncbi2na_base_at(ncbi2na_data, i);
            kmer = ((kmer << 2) | static_cast<KmerInt>(code)) & mask_;

            // Check if position i (entering right edge) is ambiguous
            if (enter_cur.pos(ambig_entries) == i) {
                ambig_count++;
                if (ambig_count == 1) {
                    single_pos = i;
                    single_ncbi4na = enter_cur.ncbi4na(ambig_entries);
                }
                enter_cur.advance(ambig_entries);
            }

            uint32_t kmer_start = i - static_cast<uint32_t>(k_) + 1;

            if (ambig_count == 0) {
                callback(kmer_start, kmer);
            } else if (max_expansion <= 1) {
                // Expansion disabled: skip all degenerate k-mers
            } else if (ambig_count == 1) {
                // Fast path: single degenerate base
                int ec = ncbi4na_expansion_count(single_ncbi4na);
                if (ec <= max_expansion) {
                    int bases_from_right = static_cast<int>(i - single_pos);
                    AmbigInfo info;
                    info.ncbi4na = single_ncbi4na;
                    info.bit_offset = bases_from_right * 2;
                    ambig_callback(kmer_start, kmer, &info, 1);
                }
            } else {
                // Multi-degen path: collect positions and check expansion product
                Cursor tmp = leave_cur;
                uint32_t win_start = kmer_start;
                // Advance tmp to first ambig position in window
                while (tmp.pos(ambig_entries) < win_start) {
                    tmp.advance(ambig_entries);
                }
                int product = 1;
                int info_count = 0;
                AmbigInfo infos[MAX_K];
                bool exceeded = false;
                uint32_t win_end = i;
                while (tmp.pos(ambig_entries) <= win_end) {
                    uint32_t apos = tmp.pos(ambig_entries);
                    uint8_t a4na = tmp.ncbi4na(ambig_entries);
                    int ec = ncbi4na_expansion_count(a4na);
                    product *= ec;
                    if (product > max_expansion) {
                        exceeded = true;
                        break;
                    }
                    int bases_from_right = static_cast<int>(i - apos);
                    infos[info_count].ncbi4na = a4na;
                    infos[info_count].bit_offset = bases_from_right * 2;
                    info_count++;
                    tmp.advance(ambig_entries);
                }
                if (!exceeded) {
                    ambig_callback(kmer_start, kmer, infos, info_count);
                }
            }

            // Check if the base leaving the window (kmer_start) was ambiguous
            if (leave_cur.pos(ambig_entries) == kmer_start) {
                ambig_count--;
                leave_cur.advance(ambig_entries);

                // If count dropped to 1, find which ambig base remains
                if (ambig_count == 1) {
                    // Scan from leave_cur to find the first ambig base
                    // in window [kmer_start+1, i+1). Bounded by k iterations.
                    Cursor tmp = leave_cur;
                    uint32_t win_start = kmer_start + 1;
                    while (tmp.pos(ambig_entries) < win_start) {
                        tmp.advance(ambig_entries);
                    }
                    single_pos = tmp.pos(ambig_entries);
                    single_ncbi4na = tmp.ncbi4na(ambig_entries);
                }
            }
        }
    }

    // Scan with spaced seed templates (packed ncbi2na data).
    // Uses sliding accumulator + bitmask extraction for O(1) per-position k-mer assembly.
    template <typename Callback, typename AmbigCallback>
    void scan_spaced(const char* ncbi2na_data, uint32_t seq_length,
                      const std::vector<AmbiguityEntry>& ambig_entries,
                      const std::vector<uint32_t>& masks, int t,
                      Callback&& callback,
                      AmbigCallback&& ambig_callback,
                      int max_expansion = 4) const {
        if (static_cast<int>(seq_length) < t) return;

        const int num_masks = static_cast<int>(masks.size());
        const uint64_t accum_mask = (1ULL << (2 * t)) - 1;
        const uint32_t ambig_window_mask = (1u << t) - 1;

        // Ambiguity enter cursor (same Cursor as contiguous scan())
        struct Cursor {
            size_t entry_idx = 0;
            uint32_t run_offset = 0;
            uint32_t pos(const std::vector<AmbiguityEntry>& e) const {
                if (entry_idx >= e.size()) return UINT32_MAX;
                return e[entry_idx].position + run_offset;
            }
            uint8_t ncbi4na(const std::vector<AmbiguityEntry>& e) const {
                return e[entry_idx].ncbi4na;
            }
            void advance(const std::vector<AmbiguityEntry>& e) {
                if (entry_idx >= e.size()) return;
                run_offset++;
                if (run_offset >= e[entry_idx].run_length) {
                    entry_idx++;
                    run_offset = 0;
                }
            }
        };

        Cursor enter_cur;
        uint64_t accum = 0;
        uint32_t ambig_bits = 0;
        uint8_t window_ncbi4na[MAX_SPACED_T] = {};
        int warmup = t - 1;

        for (uint32_t i = 0; i < seq_length; i++) {
            // 1. Advance accumulator
            uint8_t code = ncbi2na_base_at(ncbi2na_data, i);
            accum = ((accum << 2) | static_cast<uint64_t>(code)) & accum_mask;

            // 2. Shift ambig bitset (oldest falls off)
            ambig_bits = (ambig_bits << 1) & ambig_window_mask;

            // 3. Check entering position for ambiguity
            if (enter_cur.pos(ambig_entries) == i) {
                ambig_bits |= 1u;
                window_ncbi4na[i % t] = enter_cur.ncbi4na(ambig_entries);
                enter_cur.advance(ambig_entries);
            } else {
                window_ncbi4na[i % t] = 0;
            }

            // 4. Warmup
            if (warmup > 0) {
                warmup--;
                continue;
            }

            // 5. Emit k-mers
            uint32_t p = i - static_cast<uint32_t>(t) + 1;

            for (int mi = 0; mi < num_masks; mi++) {
                uint32_t mask = masks[mi];
                uint32_t affected = ambig_bits & mask;

                if (affected == 0) {
                    // Fast path: no ambiguity on set-bit positions
                    KmerInt kmer = extract_for_mask<KmerInt>(accum, mask);
                    callback(p, kmer);
                } else if (max_expansion <= 1) {
                    // skip degenerate k-mers
                } else {
                    // Slow path: some set-bit positions are ambiguous
                    KmerInt kmer = extract_for_mask<KmerInt>(accum, mask);
                    const auto& ext = extractor_for_mask(mask);
                    AmbigInfo infos[MAX_K];
                    int degen_count = 0;
                    int product = 1;
                    bool exceeded = false;

                    uint32_t bits = affected;
                    while (bits != 0) {
                        int j = std::countr_zero(bits);
                        bits &= bits - 1;

                        uint32_t seq_pos = p + static_cast<uint32_t>(t - 1 - j);
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
                    }
                }
            }
        }
    }

private:
    int k_;
    KmerInt mask_;
};

} // namespace ikafssn
