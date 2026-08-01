#pragma once

// Split codec for .kix and .kpx posting lists.
//
//   .kix: FastPFor CompositeCodec<SIMDFastPFor<4>, VariableByte>
//         (PForDelta + VByte tail) over the **distinct seq_id**
//         delta stream.  Intra-sequence k-mer duplicates are removed
//         by a SIMD dedup kernel (src/index/seq_id_dedup.*) at build
//         time, so the input stream to the codec is
//         [abs_first, d1, d2, ...] with d_i >= 1 strictly.
//         Posting list layout on disk:
//           [u32 distinct_count]         — distinct seq_ids in this k-mer
//           [u32 body[(bytes-4)/4]]      — codec output, byte-unaligned
//                                          on the wire (use memcpy)
//
//   .kpx: per-(kmer, seq_id) partitioned position posting list whose
//         **decoder is driven by the .kix distinct seq_id array**.
//         Each distinct seq_id is classified into one of three kinds
//         via a 2-bit kind map; the seq_id itself is not stored in
//         the .kpx posting list (the .kix decoded array supplies the
//         resolution between rank and seq_id).
//
//             00 = short_occ1     — exactly 1 position
//             01 = short_occ_ge2  — between 2 and freq_threshold_part
//             10 = partition      — strictly more than freq_threshold_part
//             11 = reserved
//
//         Posting list layout on disk.  No fixed header: the body starts
//         at the kind map and the counts are derived at decode time
//         (partition_count / short1_count / short2_count via
//         popcount_kinds, short2_position_count by summing occ_count[]).
//           [2-bit kind map: ceil(distinct_count*2/8) bytes]
//           repeated partition_count times in .kix sid order:
//             [u32 occ_count]                     — > freq_threshold_part
//             [FOR-block stream over occ_count positions]
//           [FOR-block stream over short1_count positions]
//                                                 — absent if short1_count == 0
//           [u8 occ_count[short2_count]]          — 2..freq_threshold_part
//           [FOR-block stream over short2_position_count positions]
//                                                 — both absent if short2_count == 0
//
//            Each FOR-block is 8 + 16*b bytes:
//                [u32 min][u8 b][3 B pad][16*b bytes bitpacked (value-min)]
//            Stream tail:
//                [u8 tail_count]
//                if tail_count > 0:
//                  [u32 tail_min][u8 tail_b][bitpacked body: ceil(tail_count*tail_b/8) B]
//
//            Decoding is candidate-set-driven and needs the decoded .kix
//            distinct seq_id array: a selection pass intersects it with
//            the sorted candidate list and resolves each match to a
//            (kind, rank) pair; only the partition groups and FOR blocks
//            holding a selected rank are decoded, the rest are walked
//            past through their self-describing headers.  The result is a
//            sparse CSR over the candidates the posting list contains
//            (see open_stream_kpx_for_candidates and PosDecodeScratch).

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <vector>

namespace ikafssn::pfd {

// FastPFor SIMD block size (the codec we wrap is simdfastpfor128 =
// CompositeCodec<SIMDFastPFor<4>, VariableByte>; per-block element count is 128).
inline constexpr int kPfdBlockSize = 128;

// .kix posting list header: the leading `[u32 distinct_count]`.
// .kpx posting lists have no header, so this describes .kix only.
inline constexpr size_t kPostingListHeaderBytes = 4;

// === posting-list-level encode wrappers ===

// Encode the distinct seq_id-delta stream for a .kix posting list.
// Writes distinct_count + body into `out` (appended).  Returns the
// number of bytes written.
size_t encode_posting_kix(const uint32_t* delta_array, uint32_t count,
                          std::vector<uint8_t>& out);

// Encode the absolute-position stream for a .kpx posting list.
//   distinct_sid          length = distinct_count, sorted ascending, no dups
//   occ_count             length = distinct_count, occurrence per distinct_sid
//                         (u32 — large genomic contigs may have > 255
//                         occurrences of one k-mer in a single sequence;
//                         clusters with occ>freq_threshold_part are routed
//                         into partition groups, clusters with occ==1 enter
//                         the short_occ1 sub-bucket; the rest enter the
//                         short_occ_ge2 sub-bucket where occ <= 255 by
//                         construction since freq_threshold_part <= 255)
//   distinct_count        number of distinct seq_ids
//   abs_pos_array         length = position_count; positions are grouped by
//                         seq_id (matching distinct_sid order) and sorted
//                         ascending within each group
//   position_count        sum of occ_count
//   freq_threshold_part   max 255; per-(kmer, seq_id) occurrence count
//                         strictly greater becomes a partition group
size_t encode_posting_kpx(const uint32_t* distinct_sid,
                          const uint32_t* occ_count,
                          uint32_t distinct_count,
                          const uint32_t* abs_pos_array,
                          uint32_t position_count,
                          uint32_t freq_threshold_part,
                          std::vector<uint8_t>& out);

// === decode context (for .kix only) ===
//
// open_stream_kix materialises the entire decoded posting list into `decoded`.
// .kpx uses the candidate-set-driven API below instead.
struct StreamCtx {
    std::vector<uint32_t> decoded;
    uint32_t count = 0;     // total elements in `decoded`
};

// Initialise a StreamCtx for a .kix posting list starting at `posting_list`.
//   posting_list:  pointer to first byte of the per-k-mer posting list
//   bytes:         total byte length of the posting list
//
// Returns false on header / body size mismatch (corrupt index).
bool open_stream_kix(const uint8_t* posting_list, size_t bytes, StreamCtx& ctx);

// === .kpx candidate-set-driven decode ===
//
// Reusable scratch buffers — one per worker thread.  Every vector is reused
// across calls and its capacity grows monotonically, so the search hot path
// does no per-call allocation.
struct PosDecodeScratch {
    // One entry per candidate the posting list contains, in ascending
    // candidate-index order.  Filled by the selection pass before any
    // stream is touched; everything downstream is driven by it.
    struct Selected {
        uint32_t candidate_idx;  // index into the caller's candidate array
        uint32_t rank;           // rank of the entry within its own kind
        uint8_t  kind;           // 0 short_occ1, 1 short_occ_ge2, 2 partition
    };
    std::vector<Selected> selected;

    // Element runs the selection asked for in the FOR stream being decoded
    // (ascending, disjoint).  Rebuilt per stream.
    struct Run {
        uint32_t first;          // first element index in the stream
        uint32_t length;         // consecutive elements from `first`
    };
    std::vector<Run> selected_runs;

    // Positions of the selected entries only.  The two short buffers grow
    // to a high-water mark, so size() overstates the current call; the
    // assembly walks them in selection order instead.
    std::vector<uint32_t> selected_partition_positions;
    std::vector<uint32_t> selected_partition_offsets;  // n selected part. + 1
    std::vector<uint32_t> selected_short1_positions;   // one per selected entry
    std::vector<uint32_t> selected_short2_positions;
    std::vector<uint8_t>  short2_occ;                  // u8 occ_count[short2_count]

    // Decode result: a sparse CSR over the candidates the posting list
    // contains, in ascending candidate-index order.  The positions of
    // candidate out_candidate_idx[m] are
    // out_positions[out_offsets[m] .. out_offsets[m + 1]).
    std::vector<uint32_t> out_candidate_idx;
    std::vector<uint32_t> out_offsets;           // out_candidate_idx.size() + 1
    std::vector<uint32_t> out_positions;         // concatenated positions

    // How much of the posting list the selection let the decoder walk past,
    // for the most recent call.  A call that never reaches a stream leaves
    // that stream's counter at 0 — a posting list no candidate hits is
    // abandoned right after the selection pass, so both stream counters
    // stay 0 while skipped_kind_entries reaches distinct_count.
    uint64_t skipped_kind_entries = 0;      // not resolved one entry at a time
    uint64_t skipped_partition_groups = 0;  // group stream not decoded
    uint64_t skipped_for_blocks = 0;        // FOR block not decoded (tail counts as one)
};

// Decode a .kpx posting list against a sorted candidate seq_id array, using
// the .kix decoded distinct seq_id array (kix_decoded[0..kix_count), strictly
// increasing) to resolve position rank to seq_id.  The result lands in the
// scratch CSR described above; candidates absent from the posting list are
// not listed.  Returns false on corrupt input, in which case the CSR is
// partially filled and must not be read.
//
// Only what is actually read is validated: regions no candidate asks for
// are walked past or never visited, so a corrupt posting list that no
// candidate hits still returns true.  Nothing from it reaches the output.
bool open_stream_kpx_for_candidates(
    const uint8_t* posting_list, size_t bytes,
    const uint32_t* kix_decoded, size_t kix_count,
    const uint32_t* candidates, size_t n_candidates,
    PosDecodeScratch& scratch);

// === posting list inspection (no decode) ===

// Read the distinct_count u32 at the start of a .kix posting list.
// Returns 0 if the posting list is shorter than the 4-byte header.
inline uint32_t posting_count(const uint8_t* posting_list, size_t bytes) {
    if (bytes < kPostingListHeaderBytes) return 0;
    uint32_t cnt;
    std::memcpy(&cnt, posting_list, sizeof(uint32_t));
    return cnt;
}

// Count partition / short1 / short2 entries packed in a 2-bit kind map.
// Encoding (per pair) is
//
//   00 -> short_occ1   (incremented into *p_short1)
//   01 -> short_occ_ge2 (incremented into *p_short2)
//   10 -> partition    (incremented into *p_partition)
//   11 -> reserved (silently ignored — caller is expected to validate)
//
// km must have at least ceil(distinct_count * 2 / 8) bytes; trailing
// padding bits in the last byte must be zero (the encoder zero-fills
// the kind map before writing entries).
void popcount_kinds(const uint8_t* km,
                    uint32_t distinct_count,
                    uint32_t* p_partition,
                    uint32_t* p_short1,
                    uint32_t* p_short2) noexcept;

} // namespace ikafssn::pfd
