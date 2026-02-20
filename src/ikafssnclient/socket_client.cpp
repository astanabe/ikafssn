#include "ikafssnclient/socket_client.hpp"

#include "protocol/frame.hpp"
#include "protocol/serializer.hpp"

namespace ikafssn {

bool socket_search(int fd, const SearchRequest& req, SearchResponse& resp) {
    // Serialize and send request
    auto payload = serialize(req);
    if (!write_frame(fd, MsgType::kSearchRequest, payload)) {
        return false;
    }

    // Read response
    FrameHeader hdr;
    std::vector<uint8_t> resp_payload;
    if (!read_frame(fd, hdr, resp_payload)) {
        return false;
    }

    MsgType type = static_cast<MsgType>(hdr.msg_type);

    if (type == MsgType::kErrorResponse) {
        ErrorResponse err;
        if (deserialize(resp_payload, err)) {
            std::fprintf(stderr, "Server error %u: %s\n", err.error_code, err.message.c_str());
        }
        return false;
    }

    if (type != MsgType::kSearchResponse) {
        return false;
    }

    return deserialize(resp_payload, resp);
}

bool socket_health_check(int fd, HealthResponse& resp) {
    HealthRequest hreq;
    auto payload = serialize(hreq);
    if (!write_frame(fd, MsgType::kHealthRequest, payload)) {
        return false;
    }

    FrameHeader hdr;
    std::vector<uint8_t> resp_payload;
    if (!read_frame(fd, hdr, resp_payload)) {
        return false;
    }

    if (static_cast<MsgType>(hdr.msg_type) != MsgType::kHealthResponse) {
        return false;
    }

    return deserialize(resp_payload, resp);
}

} // namespace ikafssn
