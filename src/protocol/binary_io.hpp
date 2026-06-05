#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace ikafssn {

// Append-only little-endian byte writer used by protocol serializers.
// Backed by an external buffer so callers control allocation; methods are
// thin and inline-friendly so the resulting machine code matches hand-written
// put_*() helpers.
class BinaryWriter {
public:
    explicit BinaryWriter(std::vector<uint8_t>& buf) : buf_(buf) {}

    void u8(uint8_t v)   { buf_.push_back(v); }
    void i8(int8_t v)    { buf_.push_back(static_cast<uint8_t>(v)); }

    void u16(uint16_t v) {
        buf_.push_back(static_cast<uint8_t>(v));
        buf_.push_back(static_cast<uint8_t>(v >> 8));
    }
    void i16(int16_t v) { u16(static_cast<uint16_t>(v)); }

    void u32(uint32_t v) {
        buf_.push_back(static_cast<uint8_t>(v));
        buf_.push_back(static_cast<uint8_t>(v >> 8));
        buf_.push_back(static_cast<uint8_t>(v >> 16));
        buf_.push_back(static_cast<uint8_t>(v >> 24));
    }
    void i32(int32_t v) { u32(static_cast<uint32_t>(v)); }

    void u64(uint64_t v) {
        for (int i = 0; i < 8; i++)
            buf_.push_back(static_cast<uint8_t>(v >> (i * 8)));
    }

    // Length-prefixed string with 16-bit length.
    void str16(const std::string& s) {
        u16(static_cast<uint16_t>(s.size()));
        buf_.insert(buf_.end(), s.begin(), s.end());
    }

private:
    std::vector<uint8_t>& buf_;
};

// Bounds-checked little-endian byte reader. Each get_*() returns false on
// short reads so callers can chain failures without branching on length.
class BinaryReader {
public:
    BinaryReader(const uint8_t* data, size_t size)
        : data_(data), size_(size), pos_(0) {}

    bool has(size_t n) const { return pos_ + n <= size_; }
    size_t remaining() const { return size_ - pos_; }

    bool get_u8(uint8_t& v) {
        if (!has(1)) return false;
        v = data_[pos_++];
        return true;
    }
    bool get_i8(int8_t& v) {
        uint8_t raw;
        if (!get_u8(raw)) return false;
        v = static_cast<int8_t>(raw);
        return true;
    }

    bool get_u16(uint16_t& v) {
        if (!has(2)) return false;
        v = static_cast<uint16_t>(data_[pos_]) |
            (static_cast<uint16_t>(data_[pos_ + 1]) << 8);
        pos_ += 2;
        return true;
    }
    bool get_i16(int16_t& v) {
        uint16_t raw;
        if (!get_u16(raw)) return false;
        v = static_cast<int16_t>(raw);
        return true;
    }

    bool get_u32(uint32_t& v) {
        if (!has(4)) return false;
        v = static_cast<uint32_t>(data_[pos_]) |
            (static_cast<uint32_t>(data_[pos_ + 1]) << 8) |
            (static_cast<uint32_t>(data_[pos_ + 2]) << 16) |
            (static_cast<uint32_t>(data_[pos_ + 3]) << 24);
        pos_ += 4;
        return true;
    }
    bool get_i32(int32_t& v) {
        uint32_t raw;
        if (!get_u32(raw)) return false;
        v = static_cast<int32_t>(raw);
        return true;
    }

    bool get_u64(uint64_t& v) {
        if (!has(8)) return false;
        v = 0;
        for (int i = 0; i < 8; i++)
            v |= static_cast<uint64_t>(data_[pos_ + i]) << (i * 8);
        pos_ += 8;
        return true;
    }

    bool get_str16(std::string& s) {
        uint16_t len;
        if (!get_u16(len)) return false;
        if (!has(len)) return false;
        s.assign(reinterpret_cast<const char*>(data_ + pos_), len);
        pos_ += len;
        return true;
    }

    bool skip(size_t n) {
        if (!has(n)) return false;
        pos_ += n;
        return true;
    }

private:
    const uint8_t* data_;
    size_t size_;
    size_t pos_;
};

} // namespace ikafssn
