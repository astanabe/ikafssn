#pragma once

#include <chrono>

namespace ikafssn {

// Monotonic elapsed-time probe.  Holds only a start point; callers accumulate
// the returned seconds themselves.
class Stopwatch {
public:
    Stopwatch() : start_(Clock::now()) {}

    void reset() { start_ = Clock::now(); }

    // Seconds since the last reset()/lap(), then restart.
    double lap() {
        const Clock::time_point now = Clock::now();
        const double s = std::chrono::duration<double>(now - start_).count();
        start_ = now;
        return s;
    }

    // Seconds since the last reset()/lap(), leaving the start point alone.
    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - start_).count();
    }

private:
    using Clock = std::chrono::steady_clock;
    Clock::time_point start_;
};

} // namespace ikafssn
