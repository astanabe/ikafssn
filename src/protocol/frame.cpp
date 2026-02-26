#include "protocol/frame.hpp"

#include <cerrno>
#include <cstring>
#include <unistd.h>

namespace ikafssn {

bool write_all(int fd, const void* data, size_t n) {
    const uint8_t* p = static_cast<const uint8_t*>(data);
    size_t remaining = n;
    while (remaining > 0) {
        ssize_t w = ::write(fd, p, remaining);
        if (w < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (w == 0) return false;
        p += w;
        remaining -= static_cast<size_t>(w);
    }
    return true;
}

bool read_all(int fd, void* data, size_t n) {
    uint8_t* p = static_cast<uint8_t*>(data);
    size_t remaining = n;
    while (remaining > 0) {
        ssize_t r = ::read(fd, p, remaining);
        if (r < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (r == 0) return false; // EOF
        p += r;
        remaining -= static_cast<size_t>(r);
    }
    return true;
}

bool write_frame(int fd, MsgType type, const std::vector<uint8_t>& payload) {
    FrameHeader hdr;
    hdr.magic = FRAME_MAGIC;
    hdr.payload_length = static_cast<uint32_t>(payload.size());
    hdr.msg_type = static_cast<uint8_t>(type);
    hdr.msg_version = 3;
    hdr.reserved = 0;

    if (!write_all(fd, &hdr, sizeof(hdr))) return false;
    if (!payload.empty()) {
        if (!write_all(fd, payload.data(), payload.size())) return false;
    }
    return true;
}

bool read_frame(int fd, FrameHeader& header, std::vector<uint8_t>& payload) {
    if (!read_all(fd, &header, sizeof(header))) return false;

    if (header.magic != FRAME_MAGIC) return false;
    if (header.msg_version != 3) return false;
    if (header.payload_length > MAX_PAYLOAD_SIZE) return false;

    payload.resize(header.payload_length);
    if (header.payload_length > 0) {
        if (!read_all(fd, payload.data(), header.payload_length)) return false;
    }
    return true;
}

} // namespace ikafssn
