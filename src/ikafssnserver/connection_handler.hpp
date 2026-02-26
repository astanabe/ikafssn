#pragma once

#include <map>

#include "ikafssnserver/request_processor.hpp"
#include "protocol/frame.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "search/stage3_alignment.hpp"
#include "search/volume_searcher.hpp"
#include "util/logger.hpp"

#include <tbb/task_arena.h>

namespace ikafssn {

class Server;  // forward declaration
struct ServerConfig;  // forward declaration

// Handle a single client connection: read request, process, send response.
// Closes the connection fd when done.
void handle_connection(
    int client_fd,
    const std::map<int, KmerGroup>& kmer_groups,
    int default_k,
    const SearchConfig& default_config,
    const Stage3Config& default_stage3_config,
    const std::string& db_path,
    bool default_context_is_ratio,
    double default_context_ratio,
    uint32_t default_context_abs,
    Server& server,
    tbb::task_arena& arena,
    const Logger& logger);

} // namespace ikafssn
