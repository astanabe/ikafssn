#include "ikafssnserver/budget_pool.hpp"

#include <algorithm>

namespace ikafssn {

BudgetLease::BudgetLease(BudgetLease&& other) noexcept
    : pool_(other.pool_), amount_(other.amount_) {
    other.pool_   = nullptr;
    other.amount_ = 0;
}

BudgetLease& BudgetLease::operator=(BudgetLease&& other) noexcept {
    if (this != &other) {
        release();
        pool_   = other.pool_;
        amount_ = other.amount_;
        other.pool_   = nullptr;
        other.amount_ = 0;
    }
    return *this;
}

BudgetLease::~BudgetLease() {
    release();
}

void BudgetLease::release() {
    if (pool_) {
        pool_->release_(amount_);
        pool_   = nullptr;
        amount_ = 0;
    }
}

void BudgetPool::configure(uint64_t total, uint64_t floor_min) {
    std::lock_guard<std::mutex> lock(mu_);
    total_         = total;
    available_     = total;
    floor_min_     = floor_min;
    shutdown_      = false;
    active_leases_ = 0;
    peak_leases_   = 0;
}

BudgetLease BudgetPool::acquire(uint64_t min, uint64_t max) {
    // Pass-through: no contention, no accounting against available_, and
    // every lease sees the full total.  active_leases_/peak_leases_ are
    // still tracked so tests can observe in-flight concurrency.
    if (floor_min_ == 0) {
        std::lock_guard<std::mutex> lock(mu_);
        if (shutdown_) {
            return BudgetLease();
        }
        ++active_leases_;
        peak_leases_ = std::max(peak_leases_, active_leases_);
        return BudgetLease(this, total_);
    }

    uint64_t want_min = std::max(min, floor_min_);
    if (want_min > total_) {
        // floor_min_ is bounded above by total_ at configure() time, but
        // an unusually large caller-supplied min could still exceed total_.
        // Cap it so the request is still serviceable on a fresh pool.
        want_min = total_;
    }

    std::unique_lock<std::mutex> lock(mu_);
    cv_.wait(lock, [&] {
        return shutdown_ || available_ >= want_min;
    });
    if (shutdown_) {
        return BudgetLease();
    }
    uint64_t take = std::min(available_, max);
    if (take < want_min) take = want_min;  // satisfied by predicate
    available_ -= take;
    ++active_leases_;
    peak_leases_ = std::max(peak_leases_, active_leases_);
    return BudgetLease(this, take);
}

void BudgetPool::shutdown() {
    {
        std::lock_guard<std::mutex> lock(mu_);
        shutdown_ = true;
    }
    cv_.notify_all();
}

int BudgetPool::peak_leases() const {
    std::lock_guard<std::mutex> lock(mu_);
    return peak_leases_;
}

void BudgetPool::release_(uint64_t amount) {
    {
        std::lock_guard<std::mutex> lock(mu_);
        if (active_leases_ > 0) --active_leases_;
        if (floor_min_ != 0) {
            available_ += amount;
            if (available_ > total_) available_ = total_;
        }
    }
    cv_.notify_all();
}

} // namespace ikafssn
