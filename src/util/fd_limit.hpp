#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

// Make sure this process may hold at least `required` open file descriptors,
// raising the RLIMIT_NOFILE soft limit when it sits below that.  A BLAST DB
// volume costs two descriptors while it is open (.nin + .nsq), and CSeqDB
// holds them for the lifetime of the mapping, so a command that keeps every
// volume of a large database open must reserve the descriptors up front.
//
// Returns true when the limit is already sufficient or was successfully
// raised.  Returns false when the hard limit is too low; `error_msg` then
// states the shortfall and the knobs that lift it.  `error_msg` is cleared
// on success.
bool ensure_fd_limit(uint64_t required, std::string& error_msg);

} // namespace ikafssn
