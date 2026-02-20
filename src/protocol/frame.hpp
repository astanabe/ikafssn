#pragma once

#include <cstdint>
#include <vector>

namespace ikafssn {

// Frame magic: "IKSV" as little-endian uint32
inline constexpr uint32_t FRAME_MAGIC = 0x56534B49; // 'I','K','S','V'

// Maximum payload size: 64 MB (sanity limit)
inline constexpr uint32_t MAX_PAYLOAD_SIZE = 64 * 1024 * 1024;

// Frame header size: 12 bytes
inline constexpr size_t FRAME_HEADER_SIZE = 12;

// Message types
enum class MsgType : uint8_t {
    // Client -> Server (0x01 - 0x7F)
    kSearchRequest  = 0x01,
    kInfoRequest    = 0x02,
    kHealthRequest  = 0x03,

    // Server -> Client (0x80 - 0xFF)
    kSearchResponse = 0x81,
    kInfoResponse   = 0x82,
    kHealthResponse = 0x83,
    kErrorResponse  = 0xFF,
};

// 12-byte frame header (little-endian, packed)
struct FrameHeader {
    uint32_t magic;           // FRAME_MAGIC
    uint32_t payload_length;  // payload size in bytes
    uint8_t  msg_type;        // MsgType
    uint8_t  msg_version;     // message version (currently 1)
    uint16_t reserved;        // 0
};
static_assert(sizeof(FrameHeader) == 12, "FrameHeader must be 12 bytes");

// Write a frame header + payload to a file descriptor.
// Returns true on success.
bool write_frame(int fd, MsgType type, const std::vector<uint8_t>& payload);

// Read a frame header + payload from a file descriptor.
// Returns true on success. Validates magic and payload size.
// On success, header and payload are filled.
bool read_frame(int fd, FrameHeader& header, std::vector<uint8_t>& payload);

// Low-level: write exactly n bytes to fd.
bool write_all(int fd, const void* data, size_t n);

// Low-level: read exactly n bytes from fd.
bool read_all(int fd, void* data, size_t n);

} // namespace ikafssn
