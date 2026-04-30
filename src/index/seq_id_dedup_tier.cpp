// Phase 5i — per-tier SIMD dedup kernel.
//
// This translation unit is compiled once per ISA tier via the
// ikafssn_dedup_<tier> OBJECT libraries declared in the top-level
// CMakeLists.txt.  Each compilation is given:
//
//   -DIKAFSSN_DEDUP_TIER_NAME=<tier>   names the per-tier inner
//                                       namespace
//   -m<arch> ...                        controls the instructions the
//                                       compiler may emit for the
//                                       scalar-with-auto-vectorisation
//                                       loops below
//
// The current implementation is a single scalar pass that the compiler
// auto-vectorises at the chosen ISA tier — sufficient for the build-
// time dedup workload which is bottlenecked by the surrounding I/O and
// the per-kmer FastPFor encoder.  Hand-coded gather/compress-store
// specialisations can be added per tier later if a dedicated benchmark
// shows a need; the per-tier OBJECT library structure is in place for
// that future work.

#include "index/seq_id_dedup.hpp"

#include <cstdint>
#include <cstddef>

#ifndef IKAFSSN_DEDUP_TIER_NAME
#error "IKAFSSN_DEDUP_TIER_NAME must be set per tier"
#endif

#define IKAFSSN_DEDUP_TIER_NS_(x) ikafssn_dedup_##x
#define IKAFSSN_DEDUP_TIER_NS(x)  IKAFSSN_DEDUP_TIER_NS_(x)

namespace ikafssn::seq_id_dedup::IKAFSSN_DEDUP_TIER_NS(IKAFSSN_DEDUP_TIER_NAME) {

std::uint32_t dedup_seq_ids(const std::uint32_t* sid, std::uint32_t n,
                            std::uint32_t* distinct_sid_out,
                            std::uint8_t*  occ_count_out) {
    if (n == 0) return 0;

    std::uint32_t d = 0;
    std::uint32_t cur_sid = sid[0];
    std::uint32_t run_len = 1;

    // Scalar pass with the surrounding ISA flags letting the compiler
    // auto-vectorise the run-length loop where it can.  Run-length
    // accumulator is bounded by freq_threshold_part (<= 255) at the
    // call site, so the u8 cast does not lose information in practice.
    for (std::uint32_t i = 1; i < n; i++) {
        const std::uint32_t v = sid[i];
        if (v == cur_sid) {
            run_len++;
        } else {
            distinct_sid_out[d] = cur_sid;
            occ_count_out[d]    = static_cast<std::uint8_t>(
                run_len > 255 ? 255 : run_len);
            d++;
            cur_sid = v;
            run_len = 1;
        }
    }
    distinct_sid_out[d] = cur_sid;
    occ_count_out[d]    = static_cast<std::uint8_t>(
        run_len > 255 ? 255 : run_len);
    d++;

    return d;
}

} // namespace ikafssn::seq_id_dedup::ikafssn_dedup_<tier>
