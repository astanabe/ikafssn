#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace ikafssn {

class MmapFile {
public:
    MmapFile() = default;
    ~MmapFile();

    MmapFile(const MmapFile&) = delete;
    MmapFile& operator=(const MmapFile&) = delete;

    MmapFile(MmapFile&& other) noexcept;
    MmapFile& operator=(MmapFile&& other) noexcept;

    bool open(const std::string& path);
    void close();

    bool is_open() const { return data_ != nullptr; }
    const uint8_t* data() const { return data_; }
    size_t size() const { return size_; }

private:
    uint8_t* data_ = nullptr;
    size_t size_ = 0;
};

} // namespace ikafssn
