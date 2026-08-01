#pragma once

#include <cstdint>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ikafssn {

// Dense ids for qseqid strings in first-appearance order.  Caching the last
// key makes a run of equal qseqids O(1) with no hashing, which is the common
// case: hits reach the output writers grouped by query.
//
// The ids reference the caller's strings, so every qseqid handed to id_of()
// must outlive the QueryOrder.
class QueryOrder {
public:
    // Id of `qseqid`, assigning the next free one the first time it is seen.
    uint32_t id_of(std::string_view qseqid) {
        if (!qseqids_.empty() && qseqids_[last_id_] == qseqid) return last_id_;
        auto it = ids_.find(qseqid);
        if (it == ids_.end()) {
            const uint32_t id = static_cast<uint32_t>(qseqids_.size());
            qseqids_.push_back(qseqid);
            ids_.emplace(qseqid, id);
            last_id_ = id;
            return id;
        }
        last_id_ = it->second;
        return last_id_;
    }

    // The qseqids in first-appearance order, indexed by id.
    const std::vector<std::string_view>& qseqids() const { return qseqids_; }

private:
    std::unordered_map<std::string_view, uint32_t> ids_;
    std::vector<std::string_view> qseqids_;
    uint32_t last_id_ = 0;
};

}  // namespace ikafssn
