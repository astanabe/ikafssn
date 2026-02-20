#pragma once

#include <cstdint>
#include <cstddef>

namespace ikafssn {

// Encode a uint32_t as LEB128 into buf. Returns number of bytes written.
inline size_t varint_encode(uint32_t value, uint8_t* buf) {
    size_t n = 0;
    do {
        uint8_t byte = value & 0x7F;
        value >>= 7;
        if (value != 0)
            byte |= 0x80;
        buf[n++] = byte;
    } while (value != 0);
    return n;
}

// Decode a LEB128 uint32_t from buf. Returns number of bytes consumed.
inline size_t varint_decode(const uint8_t* buf, uint32_t& value) {
    value = 0;
    size_t n = 0;
    unsigned shift = 0;
    uint8_t byte;
    do {
        byte = buf[n++];
        value |= static_cast<uint32_t>(byte & 0x7F) << shift;
        shift += 7;
    } while (byte & 0x80);
    return n;
}

// Compute encoded size of a uint32_t in LEB128 without writing.
inline size_t varint_size(uint32_t value) {
    size_t n = 1;
    while (value >= 0x80) {
        value >>= 7;
        n++;
    }
    return n;
}

} // namespace ikafssn
