#include "io/mmap_file.hpp"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstdio>
#include <cstring>

namespace ikafssn {

MmapFile::~MmapFile() {
    close();
}

MmapFile::MmapFile(MmapFile&& other) noexcept
    : data_(other.data_), size_(other.size_) {
    other.data_ = nullptr;
    other.size_ = 0;
}

MmapFile& MmapFile::operator=(MmapFile&& other) noexcept {
    if (this != &other) {
        close();
        data_ = other.data_;
        size_ = other.size_;
        other.data_ = nullptr;
        other.size_ = 0;
    }
    return *this;
}

bool MmapFile::open(const std::string& path) {
    close();

    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        std::fprintf(stderr, "MmapFile: cannot open '%s': %s\n",
                     path.c_str(), strerror(errno));
        return false;
    }

    struct stat st;
    if (fstat(fd, &st) != 0) {
        std::fprintf(stderr, "MmapFile: fstat failed for '%s': %s\n",
                     path.c_str(), strerror(errno));
        ::close(fd);
        return false;
    }

    if (st.st_size == 0) {
        std::fprintf(stderr, "MmapFile: file '%s' is empty\n", path.c_str());
        ::close(fd);
        return false;
    }

    void* mapped = ::mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);

    if (mapped == MAP_FAILED) {
        std::fprintf(stderr, "MmapFile: mmap failed for '%s': %s\n",
                     path.c_str(), strerror(errno));
        return false;
    }

    data_ = static_cast<uint8_t*>(mapped);
    size_ = static_cast<size_t>(st.st_size);
    return true;
}

bool MmapFile::advise(int advice) {
    if (!data_) return false;
    return ::madvise(data_, size_, advice) == 0;
}

bool MmapFile::advise(size_t offset, size_t length, int advice) {
    if (!data_) return false;
    if (offset >= size_) return false;
    if (offset + length > size_) length = size_ - offset;
    // madvise requires page-aligned address
    size_t page_size = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    size_t aligned_offset = (offset / page_size) * page_size;
    size_t adjust = offset - aligned_offset;
    return ::madvise(data_ + aligned_offset, length + adjust, advice) == 0;
}

void MmapFile::close() {
    if (data_) {
        ::munmap(data_, size_);
        data_ = nullptr;
        size_ = 0;
    }
}

} // namespace ikafssn
