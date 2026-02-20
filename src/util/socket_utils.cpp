#include "util/socket_utils.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <arpa/inet.h>
#include <fcntl.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

namespace ikafssn {

bool parse_host_port(const std::string& addr, std::string& host, uint16_t& port) {
    auto colon = addr.rfind(':');
    if (colon == std::string::npos) return false;

    host = addr.substr(0, colon);
    std::string port_str = addr.substr(colon + 1);
    if (port_str.empty()) return false;

    char* end = nullptr;
    long val = std::strtol(port_str.c_str(), &end, 10);
    if (*end != '\0' || val <= 0 || val > 65535) return false;
    port = static_cast<uint16_t>(val);
    return true;
}

int unix_listen(const std::string& path, int backlog) {
    int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
    if (fd < 0) return -1;

    // Remove existing socket file
    ::unlink(path.c_str());

    struct sockaddr_un addr;
    std::memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    if (path.size() >= sizeof(addr.sun_path)) {
        ::close(fd);
        return -1;
    }
    std::strncpy(addr.sun_path, path.c_str(), sizeof(addr.sun_path) - 1);

    if (::bind(fd, reinterpret_cast<struct sockaddr*>(&addr), sizeof(addr)) < 0) {
        ::close(fd);
        return -1;
    }

    if (::listen(fd, backlog) < 0) {
        ::close(fd);
        return -1;
    }

    return fd;
}

int tcp_listen(const std::string& addr, int backlog) {
    std::string host;
    uint16_t port;
    if (!parse_host_port(addr, host, port)) return -1;

    int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) return -1;

    int opt = 1;
    ::setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

    struct sockaddr_in sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sin_family = AF_INET;
    sa.sin_port = htons(port);

    if (host.empty() || host == "0.0.0.0") {
        sa.sin_addr.s_addr = INADDR_ANY;
    } else {
        if (::inet_pton(AF_INET, host.c_str(), &sa.sin_addr) != 1) {
            ::close(fd);
            return -1;
        }
    }

    if (::bind(fd, reinterpret_cast<struct sockaddr*>(&sa), sizeof(sa)) < 0) {
        ::close(fd);
        return -1;
    }

    if (::listen(fd, backlog) < 0) {
        ::close(fd);
        return -1;
    }

    return fd;
}

int accept_connection(int listen_fd) {
    int fd;
    do {
        fd = ::accept(listen_fd, nullptr, nullptr);
    } while (fd < 0 && errno == EINTR);
    return fd;
}

int unix_connect(const std::string& path) {
    int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
    if (fd < 0) return -1;

    struct sockaddr_un addr;
    std::memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    if (path.size() >= sizeof(addr.sun_path)) {
        ::close(fd);
        return -1;
    }
    std::strncpy(addr.sun_path, path.c_str(), sizeof(addr.sun_path) - 1);

    if (::connect(fd, reinterpret_cast<struct sockaddr*>(&addr), sizeof(addr)) < 0) {
        ::close(fd);
        return -1;
    }

    return fd;
}

int tcp_connect(const std::string& addr) {
    std::string host;
    uint16_t port;
    if (!parse_host_port(addr, host, port)) return -1;

    int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) return -1;

    struct sockaddr_in sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sin_family = AF_INET;
    sa.sin_port = htons(port);

    if (host.empty()) host = "127.0.0.1";
    if (::inet_pton(AF_INET, host.c_str(), &sa.sin_addr) != 1) {
        // Try hostname resolution
        struct addrinfo hints, *res;
        std::memset(&hints, 0, sizeof(hints));
        hints.ai_family = AF_INET;
        hints.ai_socktype = SOCK_STREAM;
        if (::getaddrinfo(host.c_str(), nullptr, &hints, &res) != 0) {
            ::close(fd);
            return -1;
        }
        sa.sin_addr = reinterpret_cast<struct sockaddr_in*>(res->ai_addr)->sin_addr;
        ::freeaddrinfo(res);
    }

    if (::connect(fd, reinterpret_cast<struct sockaddr*>(&sa), sizeof(sa)) < 0) {
        ::close(fd);
        return -1;
    }

    return fd;
}

bool set_nonblocking(int fd) {
    int flags = ::fcntl(fd, F_GETFL, 0);
    if (flags < 0) return false;
    return ::fcntl(fd, F_SETFL, flags | O_NONBLOCK) == 0;
}

void close_fd(int fd) {
    if (fd >= 0) {
        int ret;
        do {
            ret = ::close(fd);
        } while (ret < 0 && errno == EINTR);
    }
}

} // namespace ikafssn
