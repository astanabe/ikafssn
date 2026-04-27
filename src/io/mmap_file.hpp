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

    bool open(const std::string& path, bool quiet = false);
    void close();

    bool is_open() const { return data_ != nullptr; }
    const uint8_t* data() const { return data_; }
    size_t size() const { return size_; }

    // madvise hints (see madvise(2) for advice values, e.g. MADV_RANDOM)
    bool advise(int advice);
    bool advise(size_t offset, size_t length, int advice);

    // Standard pattern for indexes split into a hot dictionary head and a
    // randomly-accessed posting tail (e.g. .kix, .kpx).  When `willneed` is
    // true the dict region is marked WILLNEED + HUGEPAGE and the posting
    // region RANDOM; otherwise the entire mapping is marked RANDOM.
    void advise_dict_posting(size_t dict_size, bool willneed);

    // Standard pattern for small files accessed as a single hot region
    // (e.g. .ksx, .khx): WILLNEED for the entire mapping, or RANDOM if
    // `willneed` is false.
    void advise_all(bool willneed);

private:
    uint8_t* data_ = nullptr;
    size_t size_ = 0;
};

} // namespace ikafssn
