#include "test_util.hpp"

#include "ikafssnserver/budget_pool.hpp"

#include <atomic>
#include <chrono>
#include <random>
#include <thread>
#include <vector>

using namespace ikafssn;

// Case 1: pass-through (floor_min == 0).  Two overlapping acquires both see
// the full total; peak_leases tracks concurrent leases even though available_
// is not decremented.
static void test_passthrough_overlap() {
    std::fprintf(stderr, "-- test_passthrough_overlap\n");
    BudgetPool pool;
    pool.configure(/*total=*/1024, /*floor_min=*/0);

    BudgetLease a = pool.acquire(0, 1024);
    BudgetLease b = pool.acquire(0, 1024);
    CHECK(a.valid());
    CHECK(b.valid());
    CHECK_EQ(a.value(), uint64_t{1024});
    CHECK_EQ(b.value(), uint64_t{1024});
    CHECK_EQ(static_cast<uint64_t>(pool.peak_leases()), uint64_t{2});
}

// Case 2: serialised (floor_min == total).  A second acquire blocks until
// the first lease is released.  We assert ordering through a flag, not
// through sleep_for.
static void test_serialised_blocks() {
    std::fprintf(stderr, "-- test_serialised_blocks\n");
    BudgetPool pool;
    pool.configure(/*total=*/1024, /*floor_min=*/1024);

    BudgetLease a = pool.acquire(0, 1024);
    CHECK(a.valid());
    CHECK_EQ(a.value(), uint64_t{1024});

    std::atomic<bool> b_started{false};
    std::atomic<bool> b_acquired{false};
    std::thread b_thread([&] {
        b_started.store(true);
        BudgetLease b = pool.acquire(0, 1024);
        b_acquired.store(true);
        CHECK(b.valid());
        CHECK_EQ(b.value(), uint64_t{1024});
    });

    // Spin until thread B has actually entered acquire().
    while (!b_started.load()) std::this_thread::yield();
    // Yield a few times to let the kernel schedule B; if it had a slot
    // it would already have set b_acquired.
    for (int i = 0; i < 1000; ++i) std::this_thread::yield();
    CHECK(!b_acquired.load());

    a.release();
    b_thread.join();
    CHECK(b_acquired.load());
    CHECK_EQ(static_cast<uint64_t>(pool.peak_leases()), uint64_t{1});
}

// Case 3: acquire(min, max) returns min(available, max), but never below
// floor_min.  Use floor_min smaller than total so we can observe the cap.
static void test_max_cap_and_floor() {
    std::fprintf(stderr, "-- test_max_cap_and_floor\n");
    BudgetPool pool;
    pool.configure(/*total=*/1024, /*floor_min=*/256);

    // Caller asks for max=128 (< floor_min).  Result is floor_min.
    BudgetLease a = pool.acquire(0, 128);
    CHECK(a.valid());
    CHECK_EQ(a.value(), uint64_t{256});

    // Caller asks for max=512.  Available after a is 1024-256=768, so we
    // get 512.
    BudgetLease b = pool.acquire(0, 512);
    CHECK(b.valid());
    CHECK_EQ(b.value(), uint64_t{512});

    a.release();
    b.release();
}

// Case 4: shutdown() wakes blocked waiters with an invalid lease.
static void test_shutdown_wakes_waiters() {
    std::fprintf(stderr, "-- test_shutdown_wakes_waiters\n");
    BudgetPool pool;
    pool.configure(/*total=*/1024, /*floor_min=*/1024);

    BudgetLease a = pool.acquire(0, 1024);
    CHECK(a.valid());

    std::atomic<bool> b_returned{false};
    std::thread b_thread([&] {
        BudgetLease b = pool.acquire(0, 1024);
        // Should be invalid because shutdown fires before a releases.
        CHECK(!b.valid());
        b_returned.store(true);
    });

    // Let B park inside acquire().
    for (int i = 0; i < 1000; ++i) std::this_thread::yield();
    CHECK(!b_returned.load());

    pool.shutdown();
    b_thread.join();
    CHECK(b_returned.load());
    a.release();  // legitimate release on a still-shutdown pool
}

// Case 5: RAII destructor releases; explicit release() is idempotent.
static void test_raii_and_idempotent_release() {
    std::fprintf(stderr, "-- test_raii_and_idempotent_release\n");
    BudgetPool pool;
    pool.configure(/*total=*/1024, /*floor_min=*/1024);

    {
        BudgetLease a = pool.acquire(0, 1024);
        CHECK(a.valid());
    }
    // After destructor, another full-budget acquire must succeed without
    // blocking — proving the prior lease released.
    BudgetLease b = pool.acquire(0, 1024);
    CHECK(b.valid());
    b.release();
    b.release();  // second call must be a no-op (no double-decrement)

    BudgetLease c = pool.acquire(0, 1024);
    CHECK(c.valid());
}

// Case 6: stress test.  Random acquire/release pairs across N threads;
// in-flight sum of amount_ never exceeds total_.
static void test_invariant_under_random_ops() {
    std::fprintf(stderr, "-- test_invariant_under_random_ops\n");
    BudgetPool pool;
    const uint64_t total = 1u << 20;
    pool.configure(total, /*floor_min=*/64 * 1024);

    std::atomic<uint64_t> in_flight{0};
    std::atomic<uint64_t> peak_observed{0};
    std::atomic<bool> violated{false};

    auto worker = [&](unsigned seed) {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<uint64_t> dist(64 * 1024, total);
        for (int i = 0; i < 200; ++i) {
            uint64_t want = dist(rng);
            BudgetLease lease = pool.acquire(0, want);
            if (!lease.valid()) return;
            uint64_t v = in_flight.fetch_add(lease.value()) + lease.value();
            uint64_t prev_peak = peak_observed.load();
            while (v > prev_peak &&
                   !peak_observed.compare_exchange_weak(prev_peak, v)) {}
            if (v > total) violated.store(true);
            std::this_thread::yield();
            in_flight.fetch_sub(lease.value());
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < 4; ++i) threads.emplace_back(worker, 0xC0FFEE + i);
    for (auto& t : threads) t.join();

    CHECK(!violated.load());
    CHECK(peak_observed.load() <= total);
}

int main() {
    test_passthrough_overlap();
    test_serialised_blocks();
    test_max_cap_and_floor();
    test_shutdown_wakes_waiters();
    test_raii_and_idempotent_release();
    test_invariant_under_random_ops();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
