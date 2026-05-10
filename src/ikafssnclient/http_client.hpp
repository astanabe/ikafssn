#pragma once

#include <string>

#include "protocol/messages.hpp"

namespace ikafssn {

// HTTP Basic Auth configuration for ikafssnclient.
struct HttpAuthConfig {
    std::string userpwd;     // "user:password" for CURLOPT_USERPWD
    std::string netrc_file;  // path for CURLOPT_NETRC_FILE
};

// Build the JSON body for a /api/v1/jobs POST.  Exposed so
// async_http_client can reuse the same field encoding without
// duplicating it.
std::string build_request_json(const SearchRequest& req);

// Parse JSON server-info reply from /api/v1/info into InfoResponse.
bool parse_info_json(const std::string& body, InfoResponse& resp,
                     std::string& error_msg);

// Fetch server info via HTTP GET to ikafssnhttpd's /api/v1/info.
// Returns true on success. On failure, error_msg is set.
bool http_info(const std::string& base_url, InfoResponse& resp,
               std::string& error_msg, const HttpAuthConfig& auth = {});

} // namespace ikafssn
