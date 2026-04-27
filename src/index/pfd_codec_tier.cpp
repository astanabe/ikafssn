// Phase 5d — per-tier FastPFor wrapper.
//
// This translation unit is compiled THREE TIMES (once per ISA tier) via
// the ikafssn_pfd_{avx2,avx512bw,avx512vbmi2} OBJECT libraries declared
// in the top-level CMakeLists.txt.  Each compilation is given:
//
//   -DFastPForLib=FastPForLib_<tier>     (renames the entire FastPFor
//                                          namespace at preprocessor time
//                                          so the three sets of symbols
//                                          do not collide at link time)
//   -DIKAFSSN_PFD_TIER_NAME=<tier>       (used here to define the
//                                          tier-specific ikafssn::pfd
//                                          inner namespace name)
//   -m<arch> ...                          (controls which instructions
//                                          the FastPFor TUs and this
//                                          wrapper TU are allowed to use;
//                                          AVX-512 instructions only end
//                                          up in the avx512bw and
//                                          avx512vbmi2 objects via
//                                          compiler auto-vectorization)
//
// The dispatcher in pfd_codec.cpp picks one of the three tier namespaces
// at first use based on the runtime CPU detection done by simd_dispatch.

#include "index/pfd_codec.hpp"

// FastPForLib is rewritten by the build system; we then wrap the rewritten
// namespace under a stable alias `pfor_ns` so the rest of this file does
// not need any per-tier preprocessor sprinkling.
#include "codecfactory.h"

#include <cstring>
#include <memory>
#include <vector>

#ifndef IKAFSSN_PFD_TIER_NAME
#error "IKAFSSN_PFD_TIER_NAME must be set per tier (avx2 / avx512bw / avx512vbmi2)"
#endif

#define IKAFSSN_PFD_TIER_NS_(x) ikafssn_pfd_##x
#define IKAFSSN_PFD_TIER_NS(x)  IKAFSSN_PFD_TIER_NS_(x)

namespace ikafssn::pfd::IKAFSSN_PFD_TIER_NS(IKAFSSN_PFD_TIER_NAME) {

namespace {

namespace pfor_ns = FastPForLib;

// Per-thread codec instances.  FastPFor's IntegerCODEC implementations
// hold mutable scratch state on encode/decode, so they must not be
// shared across threads.  thread_local lazy initialisation is sufficient.
std::shared_ptr<pfor_ns::IntegerCODEC>& kix_codec() {
    thread_local pfor_ns::CODECFactory factory;
    thread_local std::shared_ptr<pfor_ns::IntegerCODEC> codec =
        factory.getFromName("simdfastpfor128");
    return codec;
}

std::shared_ptr<pfor_ns::IntegerCODEC>& kpx_codec() {
    thread_local pfor_ns::CODECFactory factory;
    thread_local std::shared_ptr<pfor_ns::IntegerCODEC> codec =
        factory.getFromName("simdfastpfor128");
    return codec;
}

// Append a posting blob ([u32 count][u32 payload_words][payload]) to
// `out` using the supplied codec.  Returns bytes appended.
std::size_t encode_with(std::shared_ptr<pfor_ns::IntegerCODEC>& codec,
                        const std::uint32_t* arr, std::uint32_t count,
                        std::vector<std::uint8_t>& out) {
    // Worst-case payload growth is bounded; 2x + 1024 covers the
    // CompositeCodec's per-block overhead and the VByte tail.
    std::vector<std::uint32_t> payload(static_cast<std::size_t>(count) * 2 + 1024);
    std::size_t nvalue = payload.size();
    if (count > 0) {
        codec->encodeArray(arr, count, payload.data(), nvalue);
    } else {
        nvalue = 0;
    }

    const std::uint32_t payload_words = static_cast<std::uint32_t>(nvalue);
    const std::size_t before = out.size();
    out.resize(before + ikafssn::pfd::kPostingHeaderBytes
                      + payload_words * sizeof(std::uint32_t));
    std::uint8_t* dst = out.data() + before;
    std::memcpy(dst,                              &count,         sizeof(std::uint32_t));
    std::memcpy(dst + sizeof(std::uint32_t),      &payload_words, sizeof(std::uint32_t));
    if (payload_words > 0) {
        std::memcpy(dst + ikafssn::pfd::kPostingHeaderBytes,
                    payload.data(),
                    payload_words * sizeof(std::uint32_t));
    }
    return out.size() - before;
}

bool decode_with(std::shared_ptr<pfor_ns::IntegerCODEC>& codec,
                 const std::uint8_t* posting, std::size_t bytes,
                 ikafssn::pfd::StreamCtx& ctx) {
    ctx.decoded.clear();
    ctx.count = 0;
    ctx.pos = 0;
    if (bytes < ikafssn::pfd::kPostingHeaderBytes) return bytes == 0;

    std::uint32_t count, payload_words;
    std::memcpy(&count,         posting,                         sizeof(std::uint32_t));
    std::memcpy(&payload_words, posting + sizeof(std::uint32_t), sizeof(std::uint32_t));

    if (bytes < ikafssn::pfd::kPostingHeaderBytes
              + std::size_t(payload_words) * sizeof(std::uint32_t)) {
        return false;
    }
    if (count == 0) {
        ctx.count = 0;
        return true;
    }

    // Copy payload into an aligned uint32_t buffer; the on-disk byte
    // stream is not 4-byte aligned in general.
    std::vector<std::uint32_t> payload(payload_words);
    std::memcpy(payload.data(),
                posting + ikafssn::pfd::kPostingHeaderBytes,
                payload_words * sizeof(std::uint32_t));

    ctx.decoded.resize(count);
    std::size_t nvalue = ctx.decoded.size();
    codec->decodeArray(payload.data(), payload_words,
                       ctx.decoded.data(), nvalue);
    if (nvalue != count) {
        ctx.decoded.clear();
        return false;
    }
    ctx.count = count;
    ctx.pos = 0;
    return true;
}

} // anonymous namespace

// ===== Per-tier API (pulled in by the dispatcher in pfd_codec.cpp) =====

std::size_t encode_posting_kix(const std::uint32_t* delta_array, std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    return encode_with(kix_codec(), delta_array, count, out);
}

std::size_t encode_posting_kpx(const std::uint32_t* abs_pos_array, std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    return encode_with(kpx_codec(), abs_pos_array, count, out);
}

bool open_stream_kix(const std::uint8_t* posting, std::size_t bytes,
                     ikafssn::pfd::StreamCtx& ctx) {
    return decode_with(kix_codec(), posting, bytes, ctx);
}

bool open_stream_kpx(const std::uint8_t* posting, std::size_t bytes,
                     ikafssn::pfd::StreamCtx& ctx) {
    return decode_with(kpx_codec(), posting, bytes, ctx);
}

} // namespace ikafssn::pfd::ikafssn_pfd_<tier>
