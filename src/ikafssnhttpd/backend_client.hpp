#pragma once

#include <string>

#include "protocol/messages.hpp"

namespace ikafssn {

// Connection mode for backend server
enum class BackendMode {
    kUnix,  // UNIX domain socket
    kTcp,   // TCP socket
};

// Client for connecting to ikafssnserver backend.
// Creates a new connection for each request (stateless proxy pattern).
class BackendClient {
public:
    BackendClient(BackendMode mode, const std::string& address);

    // Send a search request and receive response.
    // On failure, sets error_msg and returns false.
    bool search(const SearchRequest& req, SearchResponse& resp,
                std::string& error_msg);

    // Send a health check and receive response.
    bool health_check(HealthResponse& resp, std::string& error_msg);

    // Send an info request and receive response.
    bool info(InfoResponse& resp, std::string& error_msg);

    const std::string& address() const { return address_; }

private:
    BackendMode mode_;
    std::string address_;

    // Create a new connection to the backend. Returns fd or -1 on error.
    int connect();
};

} // namespace ikafssn
