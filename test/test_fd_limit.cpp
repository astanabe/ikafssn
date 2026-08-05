// ensure_fd_limit() unit test.
//
// The cases run in a fixed order because lowering RLIMIT_NOFILE's hard limit
// is irreversible for an unprivileged process: the "hard limit too low" case
// must come last.

#include "test_util.hpp"
#include "util/fd_limit.hpp"

#include <cstdio>
#include <string>

#include <sys/resource.h>

using namespace ikafssn;

static rlim_t soft_limit() {
    struct rlimit rl;
    if (::getrlimit(RLIMIT_NOFILE, &rl) != 0) return 0;
    return rl.rlim_cur;
}

static rlim_t hard_limit() {
    struct rlimit rl;
    if (::getrlimit(RLIMIT_NOFILE, &rl) != 0) return 0;
    return rl.rlim_max;
}

static bool set_limits(rlim_t soft, rlim_t hard) {
    struct rlimit rl;
    rl.rlim_cur = soft;
    rl.rlim_max = hard;
    return ::setrlimit(RLIMIT_NOFILE, &rl) == 0;
}

// A request the soft limit already covers leaves the limit untouched.
static void test_already_sufficient() {
    std::fprintf(stderr, "-- test_already_sufficient\n");

    const rlim_t before = soft_limit();
    CHECK(before > 0);

    std::string err = "not cleared";
    CHECK(ensure_fd_limit(1, err));
    CHECK(err.empty());
    CHECK_EQ(soft_limit(), before);
}

// A request above the soft limit but within the hard limit raises the soft
// limit to exactly what was asked for.
static void test_raises_soft_limit() {
    std::fprintf(stderr, "-- test_raises_soft_limit\n");

    const rlim_t hard = hard_limit();
    if (hard != RLIM_INFINITY && hard < 512) {
        std::fprintf(stderr, "  hard limit %llu too low to exercise the raise\n",
                     static_cast<unsigned long long>(hard));
        return;
    }

    CHECK(set_limits(128, hard));
    CHECK_EQ(soft_limit(), static_cast<rlim_t>(128));

    std::string err = "not cleared";
    CHECK(ensure_fd_limit(512, err));
    CHECK(err.empty());
    CHECK_EQ(soft_limit(), static_cast<rlim_t>(512));
    // The hard limit is left alone.
    CHECK_EQ(hard_limit(), hard);
}

// A request above the hard limit fails and explains how to lift it.
static void test_hard_limit_too_low() {
    std::fprintf(stderr, "-- test_hard_limit_too_low\n");

    CHECK(set_limits(256, 256));

    std::string err;
    CHECK(!ensure_fd_limit(4096, err));
    CHECK(!err.empty());
    CHECK(err.find("4096") != std::string::npos);
    CHECK(err.find("ulimit -n") != std::string::npos);
    CHECK(err.find("LimitNOFILE=") != std::string::npos);
    CHECK(err.find("limits.d") != std::string::npos);
    // A failed request must not move the limit.
    CHECK_EQ(soft_limit(), static_cast<rlim_t>(256));
}

int main() {
    test_already_sufficient();
    test_raises_soft_limit();
    test_hard_limit_too_low();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
