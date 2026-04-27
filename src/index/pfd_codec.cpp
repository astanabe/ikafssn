#include "index/pfd_codec.hpp"

#include "codecfactory.h"

#include <cstring>
#include <memory>
#include <stdexcept>

namespace ikafssn::pfd {

namespace {

// FastPFor codec instances hold mutable internal scratch state inside their
// encode/decode paths (per-block exception buffers, etc.), so they are NOT
// safe to share across threads.  Use one instance per thread.
//
// `simdfastpfor128` (= CompositeCodec<SIMDFastPFor<4>, VariableByte>)
// matches the plan's block=128 spec.
std::shared_ptr<FastPForLib::IntegerCODEC>& kix_codec() {
    thread_local FastPForLib::CODECFactory factory;
    thread_local std::shared_ptr<FastPForLib::IntegerCODEC> codec =
        factory.getFromName("simdfastpfor128");
    return codec;
}

std::shared_ptr<FastPForLib::IntegerCODEC>& kpx_codec() {
    thread_local FastPForLib::CODECFactory factory;
    thread_local std::shared_ptr<FastPForLib::IntegerCODEC> codec =
        factory.getFromName("simdfastpfor128");
    return codec;
}

// Encode arr of `count` uint32 with the given codec, appending the encoded
// payload (in u32 words) plus a leading [count][payload_words] header to
// `out`.  Returns total bytes appended.
size_t encode_with(std::shared_ptr<FastPForLib::IntegerCODEC>& codec,
                   const uint32_t* arr, uint32_t count,
                   std::vector<uint8_t>& out) {
    // Pre-allocate worst-case payload buffer. FastPFor's CompositeCodec
    // never inflates by more than 8 B per int (5 B varint for the tail +
    // a small per-block overhead); 2x + 1024 is safe.
    std::vector<uint32_t> payload(static_cast<size_t>(count) * 2 + 1024);
    size_t nvalue = payload.size();
    if (count > 0) {
        codec->encodeArray(arr, count, payload.data(), nvalue);
    } else {
        nvalue = 0;
    }

    const uint32_t payload_words = static_cast<uint32_t>(nvalue);
    const size_t before = out.size();
    out.resize(before + kPostingHeaderBytes + payload_words * sizeof(uint32_t));
    uint8_t* dst = out.data() + before;
    std::memcpy(dst,                         &count,         sizeof(uint32_t));
    std::memcpy(dst + sizeof(uint32_t),      &payload_words, sizeof(uint32_t));
    if (payload_words > 0) {
        std::memcpy(dst + kPostingHeaderBytes, payload.data(),
                    payload_words * sizeof(uint32_t));
    }
    return out.size() - before;
}

// Decode a posting blob of `bytes` length using `codec`. Populates ctx.
bool decode_with(std::shared_ptr<FastPForLib::IntegerCODEC>& codec,
                 const uint8_t* posting, size_t bytes, StreamCtx& ctx) {
    ctx.decoded.clear();
    ctx.count = 0;
    ctx.pos = 0;
    if (bytes < kPostingHeaderBytes) return bytes == 0;  // empty posting OK

    uint32_t count, payload_words;
    std::memcpy(&count,         posting,                    sizeof(uint32_t));
    std::memcpy(&payload_words, posting + sizeof(uint32_t), sizeof(uint32_t));

    if (bytes < kPostingHeaderBytes + size_t(payload_words) * sizeof(uint32_t)) {
        return false;
    }
    if (count == 0) {
        ctx.count = 0;
        return true;
    }

    // The on-disk payload is uint32_t-aligned in content but the byte stream
    // is unaligned in general. memcpy into an aligned buffer for the decoder.
    std::vector<uint32_t> payload(payload_words);
    std::memcpy(payload.data(), posting + kPostingHeaderBytes,
                payload_words * sizeof(uint32_t));

    ctx.decoded.resize(count);
    size_t nvalue = ctx.decoded.size();
    codec->decodeArray(payload.data(), payload_words, ctx.decoded.data(), nvalue);
    if (nvalue != count) {
        // Size mismatch — corrupt index.
        ctx.decoded.clear();
        return false;
    }
    ctx.count = count;
    ctx.pos = 0;
    return true;
}

} // namespace

size_t encode_posting_kix(const uint32_t* delta_array, uint32_t count,
                          std::vector<uint8_t>& out) {
    return encode_with(kix_codec(), delta_array, count, out);
}

size_t encode_posting_kpx(const uint32_t* abs_pos_array, uint32_t count,
                          std::vector<uint8_t>& out) {
    return encode_with(kpx_codec(), abs_pos_array, count, out);
}

bool open_stream_kix(const uint8_t* posting, size_t bytes, StreamCtx& ctx) {
    return decode_with(kix_codec(), posting, bytes, ctx);
}

bool open_stream_kpx(const uint8_t* posting, size_t bytes, StreamCtx& ctx) {
    return decode_with(kpx_codec(), posting, bytes, ctx);
}

} // namespace ikafssn::pfd
