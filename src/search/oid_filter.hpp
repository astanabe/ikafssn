#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KsxReader;

// OID filter modes
enum class OidFilterMode : uint8_t {
    kNone = 0,          // no filter, all OIDs pass
    kInclude = 1,       // only listed OIDs pass (seqidlist)
    kExclude = 2,       // listed OIDs are excluded (negative_seqidlist)
};

// OID filter: filters sequences by OID bitset.
// Built from seqidlist accessions + KsxReader accession→OID reverse map.
class OidFilter {
public:
    OidFilter() = default;

    // Build filter from accession list and KSX metadata.
    // Warns on stderr for unresolved accessions.
    void build(const std::vector<std::string>& accessions,
               const KsxReader& ksx,
               OidFilterMode mode);

    // Check if an OID passes the filter.
    bool pass(SeqId oid) const {
        if (mode_ == OidFilterMode::kNone) return true;
        if (oid >= bitset_.size()) return mode_ == OidFilterMode::kExclude;
        return mode_ == OidFilterMode::kInclude ? bitset_[oid] : !bitset_[oid];
    }

    OidFilterMode mode() const { return mode_; }

private:
    OidFilterMode mode_ = OidFilterMode::kNone;
    std::vector<bool> bitset_;
};

} // namespace ikafssn
