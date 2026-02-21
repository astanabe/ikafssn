#pragma once

#include <map>

#include "ikafssnserver/request_processor.hpp"
#include "protocol/frame.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "search/volume_searcher.hpp"
#include "util/logger.hpp"

#include <tbb/task_arena.h>

namespace ikafssn {

class Server;  // forward declaration

// Handle a single client connection: read request, process, send response.
// Closes the connection fd when done.
void handle_connection(
    int client_fd,
    const std::map<int, KmerGroup>& kmer_groups,
    int default_k,
    const SearchConfig& default_config,
    Server& server,
    tbb::task_arena& arena,
    const Logger& logger);

} // namespace ikafssn
