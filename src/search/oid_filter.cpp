#include "search/oid_filter.hpp"
#include "index/ksx_reader.hpp"

#include <cstdio>
#include <unordered_map>

namespace ikafssn {

void OidFilter::build(const std::vector<std::string>& accessions,
                      const KsxReader& ksx,
                      OidFilterMode mode) {
    mode_ = mode;
    if (mode_ == OidFilterMode::kNone || accessions.empty()) {
        mode_ = OidFilterMode::kNone;
        bitset_.clear();
        return;
    }

    uint32_t num_seqs = ksx.num_sequences();

    // Build accession → OID reverse map
    std::unordered_map<std::string, uint32_t> acc_to_oid;
    acc_to_oid.reserve(num_seqs);
    for (uint32_t oid = 0; oid < num_seqs; oid++) {
        auto acc = ksx.accession(oid);
        acc_to_oid.emplace(std::string(acc), oid);
    }

    // Resolve accessions to OIDs and build bitset
    bitset_.assign(num_seqs, false);
    for (const auto& acc : accessions) {
        auto it = acc_to_oid.find(acc);
        if (it != acc_to_oid.end()) {
            bitset_[it->second] = true;
        } else {
            std::fprintf(stderr, "OidFilter: accession '%s' not found in index\n",
                         acc.c_str());
        }
    }
}

} // namespace ikafssn
