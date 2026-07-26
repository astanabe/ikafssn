#pragma once

#include <chrono>

namespace ikafssn {

// Monotonic elapsed-time probe used to attribute wall time to search stages.
// Callers accumulate the returned seconds in their own variables, so the class
// itself stays stateless beyond the single start point.
class Stopwatch {
public:
    using Clock = std::chrono::steady_clock;

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
    Clock::time_point start_;
};

} // namespace ikafssn
