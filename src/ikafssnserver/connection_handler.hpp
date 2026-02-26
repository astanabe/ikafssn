#pragma once

#include "protocol/frame.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "util/logger.hpp"

#include <tbb/task_arena.h>

namespace ikafssn {

class Server;  // forward declaration
struct ServerConfig;  // forward declaration

// Handle a single client connection: read request, process, send response.
// Closes the connection fd when done.
void handle_connection(
    int client_fd,
    Server& server,
    const ServerConfig& config,
    tbb::task_arena& arena,
    const Logger& logger);

} // namespace ikafssn
