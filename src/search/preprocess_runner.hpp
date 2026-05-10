#pragma once

// Per-query preprocess helper shared by ikafssnsearch and ikafssnserver.
// Both binaries call `preprocess_query` for 1 or 2 template contexts per
// query and aggregate the skip / multi-degen state across templates.  This
// header exposes a single function that does that aggregation; the
// surrounding parallel_for body, warning emission, and permit bookkeeping
// stay in each binary because they differ between the two.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"

namespace ikafssn {

class KhxReader;

template <typename KmerInt>
struct PreprocessOutcome {
    bool multi_degen = false;
    uint8_t skip_reason = 0;
    std::string skip_detail;
};

// Per-template binding for one query.  `slot` is the destination
// QueryKmerData for that template; `khx` and `masks` parameterise
// preprocess_query's optional spaced-seed inputs.  An empty `masks` vector
// triggers the contiguous-k-mer path inside preprocess_query.
template <typename KmerInt>
struct PreprocessTemplateBinding {
    QueryKmerData<KmerInt>* slot = nullptr;
    const KhxReader* khx = nullptr;
    const std::vector<uint32_t>* masks = nullptr;
};

// Run preprocess_query once per template binding for a single query,
// write each result into the bound slot, and return aggregate flags
// (multi-degen across any template, first non-zero skip_reason).  The
// aggregate skip rule is "the first template that flags a skip wins".
template <typename KmerInt>
PreprocessOutcome<KmerInt> run_preprocess_one_query(
    const std::string& sequence,
    int k, uint8_t t,
    const SearchConfig& config,
    const PreprocessTemplateBinding<KmerInt>* bindings,
    size_t n_bindings) {

    PreprocessOutcome<KmerInt> out;
    for (size_t i = 0; i < n_bindings; ++i) {
        const auto& b = bindings[i];
        const std::vector<uint32_t> empty_masks;
        const std::vector<uint32_t>& masks = b.masks ? *b.masks : empty_masks;
        *b.slot = preprocess_query<KmerInt>(sequence, k, b.khx, config, t, masks);
        if (b.slot->has_multi_degen) out.multi_degen = true;
        if (b.slot->skip_reason != 0 && out.skip_reason == 0) {
            out.skip_reason = b.slot->skip_reason;
            out.skip_detail = b.slot->skip_detail;
        }
    }
    return out;
}

} // namespace ikafssn
