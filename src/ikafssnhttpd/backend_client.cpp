#include "ikafssnhttpd/backend_client.hpp"

#include "protocol/frame.hpp"
#include "protocol/serializer.hpp"
#include "util/socket_utils.hpp"

namespace ikafssn {

BackendClient::BackendClient(BackendMode mode, const std::string& address)
    : mode_(mode), address_(address) {}

int BackendClient::connect() {
    switch (mode_) {
    case BackendMode::kUnix:
        return unix_connect(address_);
    case BackendMode::kTcp:
        return tcp_connect(address_);
    }
    return -1;
}

bool BackendClient::search(const SearchRequest& req, SearchResponse& resp,
                           std::string& error_msg) {
    int fd = connect();
    if (fd < 0) {
        error_msg = "Failed to connect to backend server at " + address_;
        return false;
    }

    auto payload = serialize(req);
    if (!write_frame(fd, MsgType::kSearchRequest, payload)) {
        close_fd(fd);
        error_msg = "Failed to send search request to backend";
        return false;
    }

    FrameHeader hdr;
    std::vector<uint8_t> resp_payload;
    if (!read_frame(fd, hdr, resp_payload)) {
        close_fd(fd);
        error_msg = "Failed to read response from backend";
        return false;
    }
    close_fd(fd);

    MsgType type = static_cast<MsgType>(hdr.msg_type);

    if (type == MsgType::kErrorResponse) {
        ErrorResponse err;
        if (deserialize(resp_payload, err)) {
            error_msg = "Backend error " + std::to_string(err.error_code) +
                        ": " + err.message;
        } else {
            error_msg = "Backend returned unparseable error response";
        }
        return false;
    }

    if (type != MsgType::kSearchResponse) {
        error_msg = "Unexpected response type from backend";
        return false;
    }

    if (!deserialize(resp_payload, resp)) {
        error_msg = "Failed to deserialize search response";
        return false;
    }

    return true;
}

bool BackendClient::health_check(HealthResponse& resp, std::string& error_msg) {
    int fd = connect();
    if (fd < 0) {
        error_msg = "Failed to connect to backend server at " + address_;
        return false;
    }

    HealthRequest hreq;
    auto payload = serialize(hreq);
    if (!write_frame(fd, MsgType::kHealthRequest, payload)) {
        close_fd(fd);
        error_msg = "Failed to send health request to backend";
        return false;
    }

    FrameHeader hdr;
    std::vector<uint8_t> resp_payload;
    if (!read_frame(fd, hdr, resp_payload)) {
        close_fd(fd);
        error_msg = "Failed to read health response from backend";
        return false;
    }
    close_fd(fd);

    if (static_cast<MsgType>(hdr.msg_type) != MsgType::kHealthResponse) {
        error_msg = "Unexpected response type for health check";
        return false;
    }

    if (!deserialize(resp_payload, resp)) {
        error_msg = "Failed to deserialize health response";
        return false;
    }

    return true;
}

} // namespace ikafssn
