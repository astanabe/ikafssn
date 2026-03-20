#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include "core/config.hpp"
#include "core/ambiguity_parser.hpp"
#include "core/kmer_encoding.hpp"

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
    template <typename Callback, typename AmbigCallback>
    void scan_spaced(const char* ncbi2na_data, uint32_t seq_length,
                      const std::vector<AmbiguityEntry>& ambig_entries,
                      const std::vector<uint32_t>& masks, int t,
                      Callback&& callback,
                      AmbigCallback&& ambig_callback,
                      int max_expansion = 4) const {
        if (static_cast<int>(seq_length) < t) return;
        const uint32_t last_start = seq_length - static_cast<uint32_t>(t);

        for (uint32_t p = 0; p <= last_start; p++) {
            for (size_t mi = 0; mi < masks.size(); mi++) {
                uint32_t mask = masks[mi];
                KmerInt kmer = 0;
                int bit_pos = 0;
                int degen_count = 0;
                AmbigInfo infos[MAX_K];

                for (int j = t - 1; j >= 0; j--) {
                    if (!(mask & (1u << j))) continue;

                    uint32_t seq_pos = p + (static_cast<uint32_t>(t) - 1 - static_cast<uint32_t>(j));
                    uint8_t code = ncbi2na_base_at(ncbi2na_data, seq_pos);
                    int kmer_bit_offset = (k_ - 1 - bit_pos) * 2;

                    uint8_t amb_ncbi4na = 0;
                    for (const auto& ae : ambig_entries) {
                        if (seq_pos >= ae.position && seq_pos < ae.position + ae.run_length) {
                            amb_ncbi4na = ae.ncbi4na;
                            break;
                        }
                        if (ae.position > seq_pos) break;
                    }

                    if (amb_ncbi4na != 0) {
                        infos[degen_count].ncbi4na = amb_ncbi4na;
                        infos[degen_count].bit_offset = kmer_bit_offset;
                        degen_count++;
                    } else {
                        kmer |= static_cast<KmerInt>(code) << kmer_bit_offset;
                    }
                    bit_pos++;
                }

                if (degen_count == 0) {
                    callback(p, kmer);
                } else if (max_expansion <= 1) {
                    // skip
                } else {
                    int product = 1;
                    bool exceeded = false;
                    for (int d = 0; d < degen_count; d++) {
                        product *= ncbi4na_expansion_count(infos[d].ncbi4na);
                        if (product > max_expansion) {
                            exceeded = true;
                            break;
                        }
                    }
                    if (!exceeded) {
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
