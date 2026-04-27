#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <new>
#include <type_traits>
#include <utility>

namespace ikafssn {

// POD-only aligned heap buffer. 64-byte alignment for AVX-512 gather/scatter.
// Provides the minimal subset of std::vector needed by Stage1Buffer:
//   resize(n) — allocate or reallocate to hold n elements; on grow, contents are
//               unspecified (caller must memset / fill).
//   data()    — pointer to first element
//   size()    — number of elements
//   clear()   — set size to 0 without freeing the storage
//   release() — free the underlying storage
template <typename T>
class AlignedBuffer {
    static_assert(std::is_trivially_copyable_v<T>,
                  "AlignedBuffer is for POD only");
public:
    static constexpr std::size_t kAlignment = 64;

    AlignedBuffer() noexcept = default;

    AlignedBuffer(const AlignedBuffer&) = delete;
    AlignedBuffer& operator=(const AlignedBuffer&) = delete;

    AlignedBuffer(AlignedBuffer&& other) noexcept
        : ptr_(other.ptr_), size_(other.size_), cap_(other.cap_) {
        other.ptr_ = nullptr;
        other.size_ = 0;
        other.cap_ = 0;
    }

    AlignedBuffer& operator=(AlignedBuffer&& other) noexcept {
        if (this != &other) {
            release();
            ptr_ = other.ptr_;
            size_ = other.size_;
            cap_ = other.cap_;
            other.ptr_ = nullptr;
            other.size_ = 0;
            other.cap_ = 0;
        }
        return *this;
    }

    ~AlignedBuffer() { release(); }

    void resize(std::size_t n) {
        if (n <= cap_) {
            size_ = n;
            return;
        }
        std::size_t new_cap_elems = round_up_to_alignment(n);
        std::size_t alloc_bytes = new_cap_elems * sizeof(T);
        T* p = static_cast<T*>(std::aligned_alloc(kAlignment, alloc_bytes));
        if (!p) throw std::bad_alloc();
        if (ptr_) {
            std::memcpy(p, ptr_, size_ * sizeof(T));
            std::free(ptr_);
        }
        ptr_ = p;
        size_ = n;
        cap_ = new_cap_elems;
    }

    void clear() noexcept { size_ = 0; }

    void release() noexcept {
        if (ptr_) std::free(ptr_);
        ptr_ = nullptr;
        size_ = 0;
        cap_ = 0;
    }

    [[nodiscard]] T*       data()       noexcept { return ptr_; }
    [[nodiscard]] const T* data() const noexcept { return ptr_; }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }

private:
    static constexpr std::size_t round_up_to_alignment(std::size_t n_elems) noexcept {
        // aligned_alloc requires the requested size to be a multiple of the
        // alignment. Round the byte count up; convert back to element count.
        std::size_t bytes = n_elems * sizeof(T);
        std::size_t rounded = (bytes + kAlignment - 1) & ~(kAlignment - 1);
        std::size_t elems = rounded / sizeof(T);
        // Guard against truncation when sizeof(T) does not divide kAlignment.
        if (elems * sizeof(T) < rounded) elems++;
        return elems;
    }

    T*          ptr_  = nullptr;
    std::size_t size_ = 0;
    std::size_t cap_  = 0;
};

} // namespace ikafssn
