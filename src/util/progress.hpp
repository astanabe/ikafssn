#pragma once

#include <cstdio>
#include <cstdint>
#include <string>
#include <chrono>

namespace ikafssn {

// Simple progress display on stderr.
class Progress {
public:
    Progress(const std::string& label, uint64_t total, bool enabled = true)
        : label_(label), total_(total), enabled_(enabled),
          start_(std::chrono::steady_clock::now()) {}

    void update(uint64_t current) {
        if (!enabled_ || total_ == 0) return;

        // Only print at most once per 500ms or at 100%
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            now - last_print_).count();
        if (elapsed < 500 && current < total_) return;

        last_print_ = now;
        double pct = 100.0 * static_cast<double>(current) / static_cast<double>(total_);
        auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            now - start_).count();

        std::fprintf(stderr, "\r%s: %.1f%% (%lu/%lu) [%lds]",
                     label_.c_str(), pct,
                     static_cast<unsigned long>(current),
                     static_cast<unsigned long>(total_),
                     static_cast<long>(total_elapsed));
        std::fflush(stderr);
    }

    void finish() {
        if (!enabled_) return;
        auto now = std::chrono::steady_clock::now();
        auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            now - start_).count();
        std::fprintf(stderr, "\r%s: done (%lu items, %lds)\n",
                     label_.c_str(),
                     static_cast<unsigned long>(total_),
                     static_cast<long>(total_elapsed));
        std::fflush(stderr);
    }

private:
    std::string label_;
    uint64_t total_;
    bool enabled_;
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point last_print_;
};

} // namespace ikafssn
