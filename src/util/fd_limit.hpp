#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

// Make sure this process may hold at least `required` open file descriptors,
// raising the RLIMIT_NOFILE soft limit when it sits below that.  An open
// BLAST DB volume costs two descriptors (.nin + .nsq) for the lifetime of
// its mapping, so a command holding every volume of a large database open
// has to reserve them up front.
//
// Returns false when the hard limit is too low; `error_msg` then states the
// shortfall and the knobs that lift it, and is cleared on success.
bool ensure_fd_limit(uint64_t required, std::string& error_msg);

} // namespace ikafssn
