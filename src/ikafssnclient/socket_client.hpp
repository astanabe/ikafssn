#pragma once

#include <string>

#include "protocol/messages.hpp"

namespace ikafssn {

// Send a search request over a socket connection and receive the response.
// fd: connected socket file descriptor.
// Returns true on success.
bool socket_search(int fd, const SearchRequest& req, SearchResponse& resp);

// Send a health check over a socket and receive the response.
bool socket_health_check(int fd, HealthResponse& resp);

} // namespace ikafssn
