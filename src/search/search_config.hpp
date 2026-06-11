#pragma once

#include <cstdint>

#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"

namespace ikafssn {

struct SearchConfig {
    Stage1Config stage1;
    Stage2Config stage2;
    uint32_t nresult = 0;       // max results per query (0 = unlimited)
    uint8_t  mode = 2;          // 1 = stage1 only, 2 = stage1+stage2, 3 = stage1+stage2+stage3
                                // (CLI binaries override this default; library callers
                                // and tests get stage1+stage2 unless they opt out).
    uint8_t  sort_score = 2;    // 1 = stage1 score, 2 = chainscore, 3 = alnscore
    int8_t   strand = 2;       // 1 = plus only, -1 = minus only, 2 = both
    uint8_t  accept_qdegen = 1; // 0 = skip queries with degenerate bases, 1 = accept (expand)
    double min_stage1_score_frac = 0; // 0 = disabled, 0 < P < 1 = fractional mode
    uint16_t max_degen_expand = 16;  // max degenerate expansion per k-mer (0/1: disable)
    uint8_t  t = 0;  // template length (0 = contiguous)
    // Minimum query length.  Queries shorter than this are skipped
    // with kSkipQueryTooShort (the span-based check still runs as a
    // fallback for the t > 0 case).  0 disables the filter.
    uint32_t min_query_length = 0;
    // Maximum query length.  Set from the index's overlap_length (from
    // the .ksx header) when fragment splitting is active.  Queries
    // longer than this are skipped with kSkipQueryTooLong because dedup
    // correctness relies on every chain hit fitting inside at most two
    // adjacent fragments (i.e. the overlap window).  0 disables the
    // filter (e.g. when fragment splitting is off).
    uint32_t max_query_length = 0;
};

} // namespace ikafssn
