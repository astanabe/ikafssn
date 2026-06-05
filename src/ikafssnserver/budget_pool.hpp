#pragma once

#include <condition_variable>
#include <cstdint>
#include <mutex>

namespace ikafssn {

class BudgetPool;

// RAII handle returned by BudgetPool::acquire().  Releases the leased amount
// on destruction.  An empty / default-constructed lease is invalid.
class BudgetLease {
public:
    BudgetLease() = default;
    BudgetLease(BudgetLease&& other) noexcept;
    BudgetLease& operator=(BudgetLease&& other) noexcept;
    BudgetLease(const BudgetLease&)            = delete;
    BudgetLease& operator=(const BudgetLease&) = delete;
    ~BudgetLease();

    uint64_t value() const { return amount_; }
    bool     valid() const { return pool_ != nullptr; }
    void     release();

private:
    friend class BudgetPool;
    BudgetLease(BudgetPool* p, uint64_t amt) : pool_(p), amount_(amt) {}
    BudgetPool* pool_   = nullptr;
    uint64_t    amount_ = 0;
};

// Inter-request posting-budget pool used by ikafssnserver.  The persistent
// posting_budget computed by Server::apply_madvise_budget is the total
// capacity; concurrent requests share that capacity through lease/release.
//
// floor_min == 0 enables pass-through mode (no contention, no decrement);
// each acquire() returns lease(total) immediately and never blocks, so every
// request sees the full posting_budget.  This is the default.
//
// floor_min > 0 enables blocking mode.  acquire(min, max) blocks on a
// condition variable until at least max(min, floor_min_) bytes are free,
// then takes min(available, max).  shutdown() wakes all waiters with an
// invalid lease.
class BudgetPool {
public:
    void configure(uint64_t total, uint64_t floor_min);

    BudgetLease acquire(uint64_t min, uint64_t max);

    void shutdown();

    uint64_t total()     const { return total_; }
    uint64_t floor_min() const { return floor_min_; }

    // Test-only: peak number of concurrent leases observed since configure().
    int peak_leases() const;

private:
    friend class BudgetLease;
    void release_(uint64_t amount);

    mutable std::mutex      mu_;
    std::condition_variable cv_;
    uint64_t total_         = 0;
    uint64_t available_     = 0;
    uint64_t floor_min_     = 0;
    bool     shutdown_      = false;
    int      active_leases_ = 0;
    int      peak_leases_   = 0;
};

} // namespace ikafssn
