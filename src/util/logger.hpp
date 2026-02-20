#pragma once

#include <cstdio>
#include <cstdarg>

namespace ikafssn {

// Simple logger writing to stderr. Thread-safe if fprintf is thread-safe.
class Logger {
public:
    enum Level { kError = 0, kWarn = 1, kInfo = 2, kDebug = 3 };

    explicit Logger(Level level = kInfo) : level_(level) {}

    void set_level(Level level) { level_ = level; }
    Level level() const { return level_; }
    bool verbose() const { return level_ >= kDebug; }

    void error(const char* fmt, ...) const {
        if (level_ < kError) return;
        va_list ap;
        va_start(ap, fmt);
        log_impl("ERROR", fmt, ap);
        va_end(ap);
    }

    void warn(const char* fmt, ...) const {
        if (level_ < kWarn) return;
        va_list ap;
        va_start(ap, fmt);
        log_impl("WARN", fmt, ap);
        va_end(ap);
    }

    void info(const char* fmt, ...) const {
        if (level_ < kInfo) return;
        va_list ap;
        va_start(ap, fmt);
        log_impl("INFO", fmt, ap);
        va_end(ap);
    }

    void debug(const char* fmt, ...) const {
        if (level_ < kDebug) return;
        va_list ap;
        va_start(ap, fmt);
        log_impl("DEBUG", fmt, ap);
        va_end(ap);
    }

private:
    Level level_;

    static void log_impl(const char* tag, const char* fmt, va_list ap) {
        std::fprintf(stderr, "[%s] ", tag);
        std::vfprintf(stderr, fmt, ap);
        std::fprintf(stderr, "\n");
    }
};

} // namespace ikafssn
