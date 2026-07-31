#include "io/compressed_stream.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include <bzlib.h>
#include <lzma.h>
#include <zlib.h>
#include <zstd.h>

namespace ikafssn {

namespace {

// ---------------------------------------------------------------------------
// Suffix and magic detection
// ---------------------------------------------------------------------------

static std::string lower_suffix_after_last_dot(const std::string& path) {
    auto pos = path.find_last_of('.');
    if (pos == std::string::npos) return {};
    std::string s = path.substr(pos + 1);
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

// ---------------------------------------------------------------------------
// File-descriptor source with optional saved prefix.
//
// Used both for stdin and for regular files so the magic-byte detection
// path is uniform.  The streambuf reads up to 6 bytes at construction
// time, hands them to the caller for inspection, and then on underflow()
// first replays the saved prefix into the codec, followed by raw read(2)
// from the file descriptor.
// ---------------------------------------------------------------------------

class PrefixedFdSource : public std::streambuf {
public:
    PrefixedFdSource(int fd, bool owns_fd,
                     const unsigned char* prefix, std::size_t prefix_len)
        : fd_(fd), owns_fd_(owns_fd) {
        prefix_remaining_.assign(prefix, prefix + prefix_len);
        prefix_pos_ = 0;
    }

    ~PrefixedFdSource() override {
        if (owns_fd_ && fd_ >= 0) ::close(fd_);
    }

    PrefixedFdSource(const PrefixedFdSource&) = delete;
    PrefixedFdSource& operator=(const PrefixedFdSource&) = delete;

protected:
    int_type underflow() override {
        if (gptr() < egptr()) {
            return traits_type::to_int_type(*gptr());
        }
        std::size_t n = 0;
        if (prefix_pos_ < prefix_remaining_.size()) {
            std::size_t take = std::min(buf_.size(),
                                         prefix_remaining_.size() - prefix_pos_);
            std::memcpy(buf_.data(),
                        prefix_remaining_.data() + prefix_pos_, take);
            prefix_pos_ += take;
            n = take;
        } else {
            ssize_t r;
            do {
                r = ::read(fd_, buf_.data(), buf_.size());
            } while (r < 0 && errno == EINTR);
            if (r <= 0) return traits_type::eof();
            n = static_cast<std::size_t>(r);
        }
        setg(buf_.data(), buf_.data(), buf_.data() + n);
        return traits_type::to_int_type(*gptr());
    }

private:
    int fd_;
    bool owns_fd_;
    std::array<char, 64 * 1024> buf_{};
    std::vector<unsigned char> prefix_remaining_;
    std::size_t prefix_pos_ = 0;
};

// Read up to `want` bytes from `fd` into `out`.  Tolerates short reads
// (e.g. piped stdin where a single read(2) returns less than `want`).
static std::size_t read_prefix(int fd, unsigned char* out, std::size_t want) {
    std::size_t got = 0;
    while (got < want) {
        ssize_t r;
        do {
            r = ::read(fd, out + got, want - got);
        } while (r < 0 && errno == EINTR);
        if (r <= 0) break;
        got += static_cast<std::size_t>(r);
    }
    return got;
}

// ---------------------------------------------------------------------------
// Source-reading helper: pull up to `want` bytes from a streambuf into the
// codec's input window.  Returns 0 on EOF.
//
// Reads only what underflow() exposes in a single call rather than
// blocking until `want` is filled — this is important for sources like
// stdin where a short read carries useful work for the codec to process,
// and it lets the codec's own buffering interleave with the underlying
// I/O.
// ---------------------------------------------------------------------------

static std::size_t pull_from_source(std::streambuf& src,
                                    unsigned char* out, std::size_t want) {
    std::streambuf::int_type c = src.sgetc();
    if (c == std::streambuf::traits_type::eof()) return 0;
    std::streamsize avail = src.in_avail();
    if (avail <= 0) {
        out[0] = static_cast<unsigned char>(c);
        src.sbumpc();
        return 1;
    }
    std::streamsize take = std::min<std::streamsize>(avail,
                              static_cast<std::streamsize>(want));
    std::streamsize r = src.sgetn(reinterpret_cast<char*>(out), take);
    return r > 0 ? static_cast<std::size_t>(r) : 0;
}

// ---------------------------------------------------------------------------
// Per-codec input streambufs
// ---------------------------------------------------------------------------

class GzipInputStreambuf : public std::streambuf {
public:
    explicit GzipInputStreambuf(std::unique_ptr<std::streambuf> src)
        : src_(std::move(src)) {
        std::memset(&zs_, 0, sizeof(zs_));
        // 31 = gzip + zlib auto window (15 + 16).  Concatenated streams
        // are handled via inflateReset() after Z_STREAM_END.
        if (inflateInit2(&zs_, 31) != Z_OK) {
            ok_ = false;
        }
    }

    ~GzipInputStreambuf() override {
        if (ok_) inflateEnd(&zs_);
    }

    GzipInputStreambuf(const GzipInputStreambuf&) = delete;
    GzipInputStreambuf& operator=(const GzipInputStreambuf&) = delete;

protected:
    int_type underflow() override {
        if (!ok_) return traits_type::eof();
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());

        zs_.next_out = reinterpret_cast<Bytef*>(out_buf_.data());
        zs_.avail_out = static_cast<uInt>(out_buf_.size());

        while (zs_.avail_out == out_buf_.size()) {
            if (zs_.avail_in == 0) {
                std::size_t n = pull_from_source(*src_, in_buf_.data(),
                                                  in_buf_.size());
                if (n == 0) {
                    // Stream ended.  If we already produced output earlier
                    // in this call, return it; otherwise EOF.
                    if (zs_.avail_out < out_buf_.size()) break;
                    return traits_type::eof();
                }
                zs_.next_in = in_buf_.data();
                zs_.avail_in = static_cast<uInt>(n);
            }
            int r = inflate(&zs_, Z_NO_FLUSH);
            if (r == Z_STREAM_END) {
                // Concatenated stream: reset and keep going on remaining
                // input bytes (do NOT call inflateEnd + inflateInit2 — that
                // would lose zs_.next_in / zs_.avail_in).
                inflateReset(&zs_);
                if (zs_.avail_in == 0) {
                    // Look ahead — there might still be more concatenated
                    // members in the underlying source.
                    std::size_t n = pull_from_source(*src_, in_buf_.data(),
                                                      in_buf_.size());
                    if (n == 0) break;
                    zs_.next_in = in_buf_.data();
                    zs_.avail_in = static_cast<uInt>(n);
                }
                continue;
            }
            if (r != Z_OK) {
                ok_ = false;
                std::fprintf(stderr,
                    "compressed_stream: gzip inflate failed (code %d)\n", r);
                return traits_type::eof();
            }
        }

        std::size_t produced = out_buf_.size() - zs_.avail_out;
        if (produced == 0) return traits_type::eof();
        setg(out_buf_.data(), out_buf_.data(), out_buf_.data() + produced);
        return traits_type::to_int_type(*gptr());
    }

private:
    std::unique_ptr<std::streambuf> src_;
    z_stream zs_{};
    bool ok_ = true;
    std::array<unsigned char, 64 * 1024> in_buf_{};
    std::array<char, 64 * 1024> out_buf_{};
};

class Bzip2InputStreambuf : public std::streambuf {
public:
    explicit Bzip2InputStreambuf(std::unique_ptr<std::streambuf> src)
        : src_(std::move(src)) {
        if (BZ2_bzDecompressInit(&bz_, 0, 0) != BZ_OK) {
            ok_ = false;
        } else {
            inited_ = true;
        }
    }

    ~Bzip2InputStreambuf() override {
        if (inited_) BZ2_bzDecompressEnd(&bz_);
    }

    Bzip2InputStreambuf(const Bzip2InputStreambuf&) = delete;
    Bzip2InputStreambuf& operator=(const Bzip2InputStreambuf&) = delete;

protected:
    int_type underflow() override {
        if (!ok_) return traits_type::eof();
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());

        bz_.next_out = reinterpret_cast<char*>(out_buf_.data());
        bz_.avail_out = static_cast<unsigned>(out_buf_.size());

        while (bz_.avail_out == out_buf_.size()) {
            if (in_avail_ == 0) {
                std::size_t n = pull_from_source(*src_, in_buf_.data(),
                                                  in_buf_.size());
                if (n == 0) {
                    if (bz_.avail_out < out_buf_.size()) break;
                    return traits_type::eof();
                }
                in_avail_ = n;
                bz_.next_in = reinterpret_cast<char*>(in_buf_.data());
                bz_.avail_in = static_cast<unsigned>(n);
            }
            int r = BZ2_bzDecompress(&bz_);
            // bz_stream's own avail_in was decremented by BZ2_bzDecompress,
            // mirror that into our shadow counter.
            in_avail_ = static_cast<std::size_t>(bz_.avail_in);
            if (r == BZ_STREAM_END) {
                // Concatenated stream: tear down, re-init, then re-attach
                // any remaining input bytes to the freshly-initialised
                // bz_stream (BZ2_bzDecompressInit zeroes the entire
                // bz_stream, including next_in/avail_in/next_out/avail_out).
                std::size_t produced_so_far = out_buf_.size() - bz_.avail_out;
                std::size_t leftover = in_avail_;
                unsigned char* leftover_ptr = nullptr;
                if (leftover > 0) {
                    leftover_ptr = reinterpret_cast<unsigned char*>(bz_.next_in);
                }
                BZ2_bzDecompressEnd(&bz_);
                inited_ = false;
                if (BZ2_bzDecompressInit(&bz_, 0, 0) != BZ_OK) {
                    ok_ = false;
                    return traits_type::eof();
                }
                inited_ = true;
                if (leftover > 0) {
                    if (leftover_ptr != in_buf_.data()) {
                        std::memmove(in_buf_.data(), leftover_ptr, leftover);
                    }
                    bz_.next_in = reinterpret_cast<char*>(in_buf_.data());
                    bz_.avail_in = static_cast<unsigned>(leftover);
                    in_avail_ = leftover;
                } else {
                    in_avail_ = 0;
                    // Try to read another concatenated member.
                    std::size_t n = pull_from_source(*src_, in_buf_.data(),
                                                      in_buf_.size());
                    if (n == 0) {
                        // Truly exhausted — exit the outer loop.  If the
                        // outer loop's condition (avail_out unchanged)
                        // would still hold (no output produced yet), we
                        // need to break manually since we just zeroed
                        // avail_out via re-init.
                        if (produced_so_far == 0) {
                            return traits_type::eof();
                        }
                        // Restore produced output's framing and exit.
                        char* base = reinterpret_cast<char*>(out_buf_.data());
                        setg(base, base, base + produced_so_far);
                        return traits_type::to_int_type(*gptr());
                    }
                    bz_.next_in = reinterpret_cast<char*>(in_buf_.data());
                    bz_.avail_in = static_cast<unsigned>(n);
                    in_avail_ = n;
                }
                bz_.next_out = reinterpret_cast<char*>(out_buf_.data()) +
                                produced_so_far;
                bz_.avail_out = static_cast<unsigned>(out_buf_.size() -
                                                       produced_so_far);
                if (produced_so_far > 0) {
                    // We already have output; bail and return it before
                    // running another frame on the same underflow call.
                    char* base = reinterpret_cast<char*>(out_buf_.data());
                    setg(base, base, base + produced_so_far);
                    return traits_type::to_int_type(*gptr());
                }
                continue;
            }
            if (r != BZ_OK) {
                ok_ = false;
                std::fprintf(stderr,
                    "compressed_stream: bzip2 decompress failed (code %d)\n", r);
                return traits_type::eof();
            }
        }

        std::size_t produced = out_buf_.size() - bz_.avail_out;
        if (produced == 0) return traits_type::eof();
        char* base = reinterpret_cast<char*>(out_buf_.data());
        setg(base, base, base + produced);
        return traits_type::to_int_type(*gptr());
    }

private:
    std::unique_ptr<std::streambuf> src_;
    bz_stream bz_{};
    bool ok_ = true;
    bool inited_ = false;
    std::array<unsigned char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
    std::size_t in_avail_ = 0;
};

class XzInputStreambuf : public std::streambuf {
public:
    explicit XzInputStreambuf(std::unique_ptr<std::streambuf> src)
        : src_(std::move(src)) {
        std::memset(&strm_, 0, sizeof(strm_));
        if (lzma_stream_decoder(&strm_, UINT64_MAX, LZMA_CONCATENATED)
            != LZMA_OK) {
            ok_ = false;
        } else {
            inited_ = true;
        }
    }

    ~XzInputStreambuf() override {
        if (inited_) lzma_end(&strm_);
    }

    XzInputStreambuf(const XzInputStreambuf&) = delete;
    XzInputStreambuf& operator=(const XzInputStreambuf&) = delete;

protected:
    int_type underflow() override {
        if (!ok_) return traits_type::eof();
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());

        strm_.next_out = out_buf_.data();
        strm_.avail_out = out_buf_.size();

        while (strm_.avail_out == out_buf_.size()) {
            lzma_action action = LZMA_RUN;
            if (strm_.avail_in == 0) {
                std::size_t n = pull_from_source(*src_, in_buf_.data(),
                                                  in_buf_.size());
                if (n == 0) {
                    action = LZMA_FINISH;
                    strm_.next_in = nullptr;
                    strm_.avail_in = 0;
                } else {
                    strm_.next_in = in_buf_.data();
                    strm_.avail_in = n;
                }
            }
            lzma_ret r = lzma_code(&strm_, action);
            if (r == LZMA_STREAM_END) {
                // LZMA_CONCATENATED already handled multi-stream input
                // internally — we should be at EOF too.
                break;
            }
            if (r != LZMA_OK) {
                ok_ = false;
                std::fprintf(stderr,
                    "compressed_stream: xz decode failed (code %d)\n",
                    static_cast<int>(r));
                return traits_type::eof();
            }
            if (action == LZMA_FINISH && strm_.avail_in == 0 &&
                strm_.avail_out == out_buf_.size()) {
                // Ran out of input but lzma did not complete; treat as EOF.
                break;
            }
        }

        std::size_t produced = out_buf_.size() - strm_.avail_out;
        if (produced == 0) return traits_type::eof();
        char* base = reinterpret_cast<char*>(out_buf_.data());
        setg(base, base, base + produced);
        return traits_type::to_int_type(*gptr());
    }

private:
    std::unique_ptr<std::streambuf> src_;
    lzma_stream strm_ = LZMA_STREAM_INIT;
    bool ok_ = true;
    bool inited_ = false;
    std::array<unsigned char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
};

class ZstdInputStreambuf : public std::streambuf {
public:
    explicit ZstdInputStreambuf(std::unique_ptr<std::streambuf> src)
        : src_(std::move(src)) {
        ds_ = ZSTD_createDStream();
        if (!ds_) { ok_ = false; return; }
        // Cap window memory against malicious frames (1 << 27 = 128 MiB).
        ZSTD_DCtx_setParameter(ds_, ZSTD_d_windowLogMax, 27);
    }

    ~ZstdInputStreambuf() override {
        if (ds_) ZSTD_freeDStream(ds_);
    }

    ZstdInputStreambuf(const ZstdInputStreambuf&) = delete;
    ZstdInputStreambuf& operator=(const ZstdInputStreambuf&) = delete;

protected:
    int_type underflow() override {
        if (!ok_) return traits_type::eof();
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());

        ZSTD_outBuffer out{out_buf_.data(), out_buf_.size(), 0};
        while (out.pos == 0) {
            if (in_.pos == in_.size) {
                std::size_t n = pull_from_source(*src_, in_buf_.data(),
                                                  in_buf_.size());
                if (n == 0) {
                    if (out.pos > 0) break;
                    return traits_type::eof();
                }
                in_.src = in_buf_.data();
                in_.size = n;
                in_.pos = 0;
            }
            std::size_t r = ZSTD_decompressStream(ds_, &out, &in_);
            if (ZSTD_isError(r)) {
                ok_ = false;
                std::fprintf(stderr,
                    "compressed_stream: zstd decompress failed: %s\n",
                    ZSTD_getErrorName(r));
                return traits_type::eof();
            }
            // r == 0 marks frame end; ZSTD_decompressStream is happy to
            // continue into a concatenated frame on the next call.
        }

        std::size_t produced = out.pos;
        if (produced == 0) return traits_type::eof();
        char* base = reinterpret_cast<char*>(out_buf_.data());
        setg(base, base, base + produced);
        return traits_type::to_int_type(*gptr());
    }

private:
    std::unique_ptr<std::streambuf> src_;
    ZSTD_DStream* ds_ = nullptr;
    bool ok_ = true;
    std::array<unsigned char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
    ZSTD_inBuffer in_{nullptr, 0, 0};
};

// ---------------------------------------------------------------------------
// Per-codec output streambufs
// ---------------------------------------------------------------------------

class GzipOutputStreambuf : public std::streambuf {
public:
    GzipOutputStreambuf(std::unique_ptr<std::streambuf> sink, int level)
        : sink_(std::move(sink)) {
        std::memset(&zs_, 0, sizeof(zs_));
        // 31 = gzip wrapper.
        if (deflateInit2(&zs_, level, Z_DEFLATED, 31, 8,
                         Z_DEFAULT_STRATEGY) != Z_OK) {
            ok_ = false;
        } else {
            inited_ = true;
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
    }

    ~GzipOutputStreambuf() override {
        finalize();
    }

    GzipOutputStreambuf(const GzipOutputStreambuf&) = delete;
    GzipOutputStreambuf& operator=(const GzipOutputStreambuf&) = delete;

protected:
    int_type overflow(int_type ch) override {
        if (!ok_) return traits_type::eof();
        if (sync_buffer(false) < 0) return traits_type::eof();
        if (ch != traits_type::eof()) {
            *pptr() = static_cast<char>(ch);
            pbump(1);
        }
        return traits_type::not_eof(ch);
    }

    int sync() override {
        return sync_buffer(false);
    }

private:
    int sync_buffer(bool finish) {
        if (!ok_) return -1;
        std::size_t n = static_cast<std::size_t>(pptr() - pbase());
        zs_.next_in = reinterpret_cast<Bytef*>(pbase());
        zs_.avail_in = static_cast<uInt>(n);

        for (;;) {
            zs_.next_out = out_buf_.data();
            zs_.avail_out = static_cast<uInt>(out_buf_.size());
            int r = deflate(&zs_, finish ? Z_FINISH : Z_NO_FLUSH);
            std::size_t produced = out_buf_.size() - zs_.avail_out;
            if (produced > 0) {
                std::streamsize w = sink_->sputn(
                    reinterpret_cast<char*>(out_buf_.data()),
                    static_cast<std::streamsize>(produced));
                if (w != static_cast<std::streamsize>(produced)) {
                    ok_ = false;
                    return -1;
                }
            }
            if (finish) {
                if (r == Z_STREAM_END) break;
                if (r != Z_OK && r != Z_BUF_ERROR) { ok_ = false; return -1; }
            } else {
                if (zs_.avail_in == 0 && zs_.avail_out > 0) break;
                if (r != Z_OK) { ok_ = false; return -1; }
            }
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
        return 0;
    }

    void finalize() {
        if (inited_) {
            sync_buffer(true);
            deflateEnd(&zs_);
            inited_ = false;
        }
        // Flush the underlying sink's own buffer (e.g. std::filebuf).
        if (sink_) sink_->pubsync();
    }

    std::unique_ptr<std::streambuf> sink_;
    z_stream zs_{};
    bool ok_ = true;
    bool inited_ = false;
    std::array<char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
};

class Bzip2OutputStreambuf : public std::streambuf {
public:
    Bzip2OutputStreambuf(std::unique_ptr<std::streambuf> sink, int level)
        : sink_(std::move(sink)) {
        // BZ2_bzCompressInit rejects the (-1) sentinel; the caller is
        // expected to translate kCompressionLevelDefault before this call.
        if (BZ2_bzCompressInit(&bz_, level, 0, 0) != BZ_OK) {
            ok_ = false;
        } else {
            inited_ = true;
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
    }

    ~Bzip2OutputStreambuf() override {
        finalize();
    }

    Bzip2OutputStreambuf(const Bzip2OutputStreambuf&) = delete;
    Bzip2OutputStreambuf& operator=(const Bzip2OutputStreambuf&) = delete;

protected:
    int_type overflow(int_type ch) override {
        if (!ok_) return traits_type::eof();
        if (sync_buffer(BZ_RUN) < 0) return traits_type::eof();
        if (ch != traits_type::eof()) {
            *pptr() = static_cast<char>(ch);
            pbump(1);
        }
        return traits_type::not_eof(ch);
    }

    int sync() override {
        // bzip2 has no incremental flush below BZ_FINISH; we can only push
        // pending input through the codec.  This still drains pptr().
        return sync_buffer(BZ_RUN);
    }

private:
    int sync_buffer(int action) {
        if (!ok_) return -1;
        std::size_t n = static_cast<std::size_t>(pptr() - pbase());
        bz_.next_in = pbase();
        bz_.avail_in = static_cast<unsigned>(n);

        for (;;) {
            bz_.next_out = reinterpret_cast<char*>(out_buf_.data());
            bz_.avail_out = static_cast<unsigned>(out_buf_.size());
            int r = BZ2_bzCompress(&bz_, action);
            std::size_t produced = out_buf_.size() - bz_.avail_out;
            if (produced > 0) {
                std::streamsize w = sink_->sputn(
                    reinterpret_cast<char*>(out_buf_.data()),
                    static_cast<std::streamsize>(produced));
                if (w != static_cast<std::streamsize>(produced)) {
                    ok_ = false;
                    return -1;
                }
            }
            if (action == BZ_FINISH) {
                if (r == BZ_STREAM_END) break;
                if (r != BZ_FINISH_OK) { ok_ = false; return -1; }
            } else {
                if (bz_.avail_in == 0) break;
                if (r != BZ_RUN_OK) { ok_ = false; return -1; }
            }
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
        return 0;
    }

    void finalize() {
        if (inited_) {
            sync_buffer(BZ_FINISH);
            BZ2_bzCompressEnd(&bz_);
            inited_ = false;
        }
        if (sink_) sink_->pubsync();
    }

    std::unique_ptr<std::streambuf> sink_;
    bz_stream bz_{};
    bool ok_ = true;
    bool inited_ = false;
    std::array<char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
};

class XzOutputStreambuf : public std::streambuf {
public:
    XzOutputStreambuf(std::unique_ptr<std::streambuf> sink, int level)
        : sink_(std::move(sink)) {
        std::memset(&strm_, 0, sizeof(strm_));
        if (lzma_easy_encoder(&strm_, static_cast<uint32_t>(level),
                              LZMA_CHECK_CRC64) != LZMA_OK) {
            ok_ = false;
        } else {
            inited_ = true;
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
    }

    ~XzOutputStreambuf() override {
        finalize();
    }

    XzOutputStreambuf(const XzOutputStreambuf&) = delete;
    XzOutputStreambuf& operator=(const XzOutputStreambuf&) = delete;

protected:
    int_type overflow(int_type ch) override {
        if (!ok_) return traits_type::eof();
        if (sync_buffer(LZMA_RUN) < 0) return traits_type::eof();
        if (ch != traits_type::eof()) {
            *pptr() = static_cast<char>(ch);
            pbump(1);
        }
        return traits_type::not_eof(ch);
    }

    int sync() override {
        return sync_buffer(LZMA_RUN);
    }

private:
    int sync_buffer(lzma_action action) {
        if (!ok_) return -1;
        std::size_t n = static_cast<std::size_t>(pptr() - pbase());
        strm_.next_in = reinterpret_cast<uint8_t*>(pbase());
        strm_.avail_in = n;

        for (;;) {
            strm_.next_out = out_buf_.data();
            strm_.avail_out = out_buf_.size();
            lzma_ret r = lzma_code(&strm_, action);
            std::size_t produced = out_buf_.size() - strm_.avail_out;
            if (produced > 0) {
                std::streamsize w = sink_->sputn(
                    reinterpret_cast<char*>(out_buf_.data()),
                    static_cast<std::streamsize>(produced));
                if (w != static_cast<std::streamsize>(produced)) {
                    ok_ = false;
                    return -1;
                }
            }
            if (action == LZMA_FINISH) {
                if (r == LZMA_STREAM_END) break;
                if (r != LZMA_OK) { ok_ = false; return -1; }
            } else {
                if (strm_.avail_in == 0 && strm_.avail_out > 0) break;
                if (r != LZMA_OK) { ok_ = false; return -1; }
            }
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
        return 0;
    }

    void finalize() {
        if (inited_) {
            sync_buffer(LZMA_FINISH);
            lzma_end(&strm_);
            inited_ = false;
        }
        if (sink_) sink_->pubsync();
    }

    std::unique_ptr<std::streambuf> sink_;
    lzma_stream strm_ = LZMA_STREAM_INIT;
    bool ok_ = true;
    bool inited_ = false;
    std::array<char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
};

class ZstdOutputStreambuf : public std::streambuf {
public:
    ZstdOutputStreambuf(std::unique_ptr<std::streambuf> sink, int level)
        : sink_(std::move(sink)) {
        cs_ = ZSTD_createCStream();
        if (!cs_) { ok_ = false; return; }
        std::size_t r = ZSTD_CCtx_setParameter(cs_,
            ZSTD_c_compressionLevel, level);
        if (ZSTD_isError(r)) { ok_ = false; return; }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
    }

    ~ZstdOutputStreambuf() override {
        finalize();
    }

    ZstdOutputStreambuf(const ZstdOutputStreambuf&) = delete;
    ZstdOutputStreambuf& operator=(const ZstdOutputStreambuf&) = delete;

protected:
    int_type overflow(int_type ch) override {
        if (!ok_) return traits_type::eof();
        if (drain(ZSTD_e_continue) < 0) return traits_type::eof();
        if (ch != traits_type::eof()) {
            *pptr() = static_cast<char>(ch);
            pbump(1);
        }
        return traits_type::not_eof(ch);
    }

    int sync() override {
        return drain(ZSTD_e_continue);
    }

private:
    int drain(ZSTD_EndDirective mode) {
        if (!ok_) return -1;
        std::size_t n = static_cast<std::size_t>(pptr() - pbase());
        ZSTD_inBuffer in{pbase(), n, 0};

        for (;;) {
            ZSTD_outBuffer out{out_buf_.data(), out_buf_.size(), 0};
            std::size_t r = ZSTD_compressStream2(cs_, &out, &in, mode);
            if (ZSTD_isError(r)) { ok_ = false; return -1; }
            if (out.pos > 0) {
                std::streamsize w = sink_->sputn(
                    reinterpret_cast<char*>(out_buf_.data()),
                    static_cast<std::streamsize>(out.pos));
                if (w != static_cast<std::streamsize>(out.pos)) {
                    ok_ = false;
                    return -1;
                }
            }
            if (mode == ZSTD_e_end) {
                if (r == 0) break;
            } else {
                if (in.pos == in.size && r == 0) break;
            }
        }
        setp(in_buf_.data(), in_buf_.data() + in_buf_.size());
        return 0;
    }

    void finalize() {
        if (cs_) {
            drain(ZSTD_e_end);
            ZSTD_freeCStream(cs_);
            cs_ = nullptr;
        }
        if (sink_) sink_->pubsync();
    }

    std::unique_ptr<std::streambuf> sink_;
    ZSTD_CStream* cs_ = nullptr;
    bool ok_ = true;
    std::array<char, 64 * 1024> in_buf_{};
    std::array<unsigned char, 64 * 1024> out_buf_{};
};

// ---------------------------------------------------------------------------
// File-descriptor sink for output side
// ---------------------------------------------------------------------------

// The put area matters for uncompressed output: with kNone the caller's
// std::ostream writes straight into this streambuf, so without buffering
// every insertion would become one write(2).  The codec streambufs above
// hold their own 64 KiB input buffer and reach this class only in whole
// compressed blocks.
class FdSinkStreambuf final : public std::streambuf {
public:
    FdSinkStreambuf(int fd, bool owns_fd)
        : fd_(fd), owns_fd_(owns_fd), ok_(fd >= 0) {
        setp(buf_.data(), buf_.data() + buf_.size());
    }

    ~FdSinkStreambuf() override {
        // The kNone path has no other flush point: ~ostream does not sync
        // and OwnedOstream destroys the stream before this streambuf.
        drain();
        if (owns_fd_ && fd_ >= 0) ::close(fd_);
    }

    FdSinkStreambuf(const FdSinkStreambuf&) = delete;
    FdSinkStreambuf& operator=(const FdSinkStreambuf&) = delete;

protected:
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        if (!ok_ || n <= 0) return 0;
        if (n >= static_cast<std::streamsize>(buf_.size())) {
            // Buffered bytes precede this chunk in the stream, so they have
            // to reach the fd first.
            if (!drain()) return 0;
            return write_all(s, n);
        }
        std::streamsize room = epptr() - pptr();
        if (n > room && !drain()) return 0;
        std::memcpy(pptr(), s, static_cast<size_t>(n));
        pbump(static_cast<int>(n));
        return n;
    }

    int_type overflow(int_type ch) override {
        if (!drain()) return traits_type::eof();
        if (ch != traits_type::eof()) {
            *pptr() = static_cast<char>(ch);
            pbump(1);
        }
        return traits_type::not_eof(ch);
    }

    int sync() override {
        // Without this an explicit flush() — and the sink_->pubsync() the
        // codecs issue from finalize() — would leave the put area unwritten.
        return drain() ? 0 : -1;
    }

private:
    // Empty the put area onto the fd.  Resets the put area even on failure so
    // a later drain() cannot write the same bytes twice.
    bool drain() {
        if (!ok_) return false;
        std::streamsize n = pptr() - pbase();
        setp(buf_.data(), buf_.data() + buf_.size());
        if (n == 0) return true;
        return write_all(pbase(), n) == n;
    }

    std::streamsize write_all(const char* s, std::streamsize n) {
        std::streamsize wrote = 0;
        while (wrote < n) {
            ssize_t w;
            do {
                w = ::write(fd_, s + wrote, static_cast<size_t>(n - wrote));
            } while (w < 0 && errno == EINTR);
            if (w <= 0) {
                ok_ = false;
                std::fprintf(stderr, "compressed_stream: write failed: %s\n",
                             std::strerror(errno));
                break;
            }
            wrote += w;
        }
        return wrote;
    }

    int fd_;
    bool owns_fd_;
    bool ok_;
    std::array<char, 64 * 1024> buf_{};
};

} // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

CompressionFormat detect_format_from_extension(const std::string& path) {
    if (path.empty() || path == "-") return CompressionFormat::kNone;
    std::string ext = lower_suffix_after_last_dot(path);
    if (ext == "gz")  return CompressionFormat::kGzip;
    if (ext == "bz2") return CompressionFormat::kBzip2;
    if (ext == "xz")  return CompressionFormat::kXz;
    if (ext == "zst") return CompressionFormat::kZstd;
    return CompressionFormat::kNone;
}

CompressionFormat detect_format_from_magic(const unsigned char* prefix,
                                           std::size_t n) {
    if (n >= 2 && prefix[0] == 0x1F && prefix[1] == 0x8B) {
        return CompressionFormat::kGzip;
    }
    if (n >= 3 && prefix[0] == 'B' && prefix[1] == 'Z' && prefix[2] == 'h') {
        return CompressionFormat::kBzip2;
    }
    if (n >= 6 && prefix[0] == 0xFD && prefix[1] == '7' && prefix[2] == 'z' &&
        prefix[3] == 'X' && prefix[4] == 'Z' && prefix[5] == 0x00) {
        return CompressionFormat::kXz;
    }
    if (n >= 4 && prefix[0] == 0x28 && prefix[1] == 0xB5 &&
        prefix[2] == 0x2F && prefix[3] == 0xFD) {
        return CompressionFormat::kZstd;
    }
    return CompressionFormat::kNone;
}

bool validate_compression_level(CompressionFormat fmt, int level,
                                std::string& error_msg) {
    if (level == kCompressionLevelDefault) return true;
    auto fail = [&](const char* name, int lo, int hi) {
        error_msg = std::string("Error: -compression_level for ") + name +
                    " must be in " + std::to_string(lo) + ".." +
                    std::to_string(hi) + " (got " + std::to_string(level) + ")";
        return false;
    };
    switch (fmt) {
        case CompressionFormat::kNone:
            // No codec selected — silently accept (the level will be
            // ignored when writing uncompressed output).
            return true;
        case CompressionFormat::kGzip:
            if (level < 0 || level > 9) return fail("gzip", 0, 9);
            return true;
        case CompressionFormat::kBzip2:
            if (level < 1 || level > 9) return fail("bzip2", 1, 9);
            return true;
        case CompressionFormat::kXz:
            if (level < 0 || level > 9) return fail("xz", 0, 9);
            return true;
        case CompressionFormat::kZstd: {
            int lo = ZSTD_minCLevel();
            int hi = ZSTD_maxCLevel();
            if (level < lo || level > hi) return fail("zstd", lo, hi);
            return true;
        }
    }
    return true;
}

OwnedIstream open_input_compressed(const std::string& path,
                                    std::string& error_msg) {
    OwnedIstream out;

    int fd = -1;
    bool owns_fd = false;
    if (path == "-" || path.empty()) {
        fd = STDIN_FILENO;
        owns_fd = false;
    } else {
        fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) {
            error_msg = "Error: cannot open input file '" + path + "': " +
                        std::strerror(errno);
            return out;
        }
        owns_fd = true;
    }

    unsigned char prefix[6];
    std::size_t got = read_prefix(fd, prefix, sizeof(prefix));
    CompressionFormat fmt = detect_format_from_magic(prefix, got);

    auto fd_src = std::make_unique<PrefixedFdSource>(fd, owns_fd,
                                                      prefix, got);

    std::unique_ptr<std::streambuf> codec_sb;
    switch (fmt) {
        case CompressionFormat::kGzip:
            codec_sb = std::make_unique<GzipInputStreambuf>(std::move(fd_src));
            break;
        case CompressionFormat::kBzip2:
            codec_sb = std::make_unique<Bzip2InputStreambuf>(std::move(fd_src));
            break;
        case CompressionFormat::kXz:
            codec_sb = std::make_unique<XzInputStreambuf>(std::move(fd_src));
            break;
        case CompressionFormat::kZstd:
            codec_sb = std::make_unique<ZstdInputStreambuf>(std::move(fd_src));
            break;
        case CompressionFormat::kNone:
            codec_sb = std::move(fd_src);
            break;
    }

    out.stream = std::make_unique<std::istream>(codec_sb.get());
    out.sb = std::move(codec_sb);
    return out;
}

OwnedOstream open_output_compressed(const std::string& path, int level,
                                     std::string& error_msg) {
    OwnedOstream out;

    if (path.empty() || path == "-") {
        // Pure stdout passthrough — no codec, leave std::cout's rdbuf alive.
        out.stream = std::make_unique<std::ostream>(std::cout.rdbuf());
        return out;
    }

    CompressionFormat fmt = detect_format_from_extension(path);

    int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        error_msg = "Error: cannot open output file '" + path + "': " +
                    std::strerror(errno);
        return out;
    }
    auto fd_sink = std::make_unique<FdSinkStreambuf>(fd, /*owns_fd=*/true);

    // Translate the kCompressionLevelDefault sentinel for codecs that
    // refuse it.
    int eff_level = level;
    switch (fmt) {
        case CompressionFormat::kGzip:
            if (eff_level == kCompressionLevelDefault) eff_level = 6;
            break;
        case CompressionFormat::kBzip2:
            if (eff_level == kCompressionLevelDefault) eff_level = 9;
            break;
        case CompressionFormat::kXz:
            if (eff_level == kCompressionLevelDefault) eff_level = 6;
            break;
        case CompressionFormat::kZstd:
            if (eff_level == kCompressionLevelDefault) eff_level = ZSTD_CLEVEL_DEFAULT;
            break;
        case CompressionFormat::kNone:
            break;
    }

    std::unique_ptr<std::streambuf> codec_sb;
    switch (fmt) {
        case CompressionFormat::kGzip:
            codec_sb = std::make_unique<GzipOutputStreambuf>(
                std::move(fd_sink), eff_level);
            break;
        case CompressionFormat::kBzip2:
            codec_sb = std::make_unique<Bzip2OutputStreambuf>(
                std::move(fd_sink), eff_level);
            break;
        case CompressionFormat::kXz:
            codec_sb = std::make_unique<XzOutputStreambuf>(
                std::move(fd_sink), eff_level);
            break;
        case CompressionFormat::kZstd:
            codec_sb = std::make_unique<ZstdOutputStreambuf>(
                std::move(fd_sink), eff_level);
            break;
        case CompressionFormat::kNone:
            codec_sb = std::move(fd_sink);
            break;
    }

    out.stream = std::make_unique<std::ostream>(codec_sb.get());
    out.sb = std::move(codec_sb);
    return out;
}

} // namespace ikafssn
