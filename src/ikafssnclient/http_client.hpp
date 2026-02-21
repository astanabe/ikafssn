#pragma once

#include <string>

#include "protocol/messages.hpp"

namespace ikafssn {

// HTTP Basic Auth configuration for ikafssnclient.
struct HttpAuthConfig {
    std::string userpwd;     // "user:password" for CURLOPT_USERPWD
    std::string netrc_file;  // path for CURLOPT_NETRC_FILE
};

// Send a search request via HTTP POST to ikafssnhttpd's /api/v1/search.
// base_url: e.g., "http://example.com:8080" or "http://example.com:8080/nt"
// Returns true on success. On failure, error_msg is set.
bool http_search(const std::string& base_url, const SearchRequest& req,
                 SearchResponse& resp, std::string& error_msg,
                 const HttpAuthConfig& auth = {});

} // namespace ikafssn
