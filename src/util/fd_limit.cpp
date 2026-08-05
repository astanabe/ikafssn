#include "util/fd_limit.hpp"

#include <cerrno>
#include <cstring>

#include <sys/resource.h>

namespace ikafssn {

namespace {

std::string limit_str(rlim_t v) {
    if (v == RLIM_INFINITY) return "unlimited";
    return std::to_string(static_cast<unsigned long long>(v));
}

std::string how_to_raise(uint64_t required, rlim_t hard) {
    std::string s;
    s += "Need at least ";
    s += std::to_string(static_cast<unsigned long long>(required));
    s += " file descriptor(s), but the hard limit is ";
    s += limit_str(hard);
    s += ".\n";
    s += "Raise it before running this command:\n";
    s += "  - shell session:  ulimit -n ";
    s += std::to_string(static_cast<unsigned long long>(required));
    s += "\n";
    s += "  - all logins:     add 'ulimit' entries to /etc/security/limits.d/\n"
         "                    (e.g. '<user> hard nofile ";
    s += std::to_string(static_cast<unsigned long long>(required));
    s += "')\n";
    s += "  - systemd unit:   LimitNOFILE=";
    s += std::to_string(static_cast<unsigned long long>(required));
    s += " in the [Service] section\n";
    s += "  - kernel ceiling: fs.nr_open (per process) and fs.file-max "
         "(system wide)";
    return s;
}

} // namespace

bool ensure_fd_limit(uint64_t required, std::string& error_msg) {
    error_msg.clear();

    struct rlimit rl;
    if (::getrlimit(RLIMIT_NOFILE, &rl) != 0) {
        error_msg = "getrlimit(RLIMIT_NOFILE) failed: ";
        error_msg += std::strerror(errno);
        return false;
    }

    if (rl.rlim_cur == RLIM_INFINITY ||
        static_cast<uint64_t>(rl.rlim_cur) >= required) {
        return true;
    }

    if (rl.rlim_max != RLIM_INFINITY &&
        static_cast<uint64_t>(rl.rlim_max) < required) {
        error_msg = how_to_raise(required, rl.rlim_max);
        return false;
    }

    // Raise the soft limit to exactly what is needed.  Handing the hard
    // limit through instead fails on macOS, where it is commonly
    // RLIM_INFINITY while the kernel refuses anything above
    // kern.maxfilesperproc.
    struct rlimit want = rl;
    want.rlim_cur = static_cast<rlim_t>(required);
    if (::setrlimit(RLIMIT_NOFILE, &want) != 0) {
        error_msg = "setrlimit(RLIMIT_NOFILE, ";
        error_msg += std::to_string(static_cast<unsigned long long>(required));
        error_msg += ") failed: ";
        error_msg += std::strerror(errno);
        error_msg += "\n";
        error_msg += how_to_raise(required, rl.rlim_max);
        return false;
    }
    return true;
}

} // namespace ikafssn
