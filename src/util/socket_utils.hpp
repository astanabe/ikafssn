#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

// Create a UNIX domain socket, bind, and listen.
// Returns listening fd on success, -1 on error.
// Removes existing socket file if present.
int unix_listen(const std::string& path, int backlog = 128);

// Create a TCP socket, bind, and listen.
// addr format: "host:port" or ":port" (bind to all interfaces).
// Returns listening fd on success, -1 on error.
int tcp_listen(const std::string& addr, int backlog = 128);

// Accept a connection from a listening socket.
// Returns connected fd on success, -1 on error.
int accept_connection(int listen_fd);

// Connect to a UNIX domain socket.
// Returns connected fd on success, -1 on error.
int unix_connect(const std::string& path);

// Connect to a TCP address.
// addr format: "host:port".
// Returns connected fd on success, -1 on error.
int tcp_connect(const std::string& addr);

// Set socket to non-blocking mode.
bool set_nonblocking(int fd);

// Close a file descriptor safely (ignoring EINTR).
void close_fd(int fd);

// Parse "host:port" string. Returns false on invalid format.
bool parse_host_port(const std::string& addr, std::string& host, uint16_t& port);

} // namespace ikafssn
