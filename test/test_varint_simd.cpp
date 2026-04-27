#include "test_util.hpp"
#include "core/varint.hpp"
#include "index/varint_simd.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/posting_decoder.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

using namespace ikafssn;

// Encode `values` as a contiguous LEB128 stream and return the buffer.
static std::vector<std::uint8_t> encode_stream(const std::vector<std::uint32_t>& values) {
    std::vector<std::uint8_t> buf;
    buf.reserve(values.size() * 5);
    std::uint8_t tmp[5];
    for (auto v : values) {
        std::size_t n = varint_encode(v, tmp);
        for (std::size_t j = 0; j < n; j++) buf.push_back(tmp[j]);
    }
    return buf;
}

// Reference: scalar varint_decode loop, used as the bit-exact golden.
static int decode_stream_scalar(const std::uint8_t* in, const std::uint8_t* end,
                                std::uint32_t* out, std::size_t* consumed,
                                int max_count) {
    const std::uint8_t* p = in;
    int i = 0;
    for (; i < max_count && p < end; i++) {
        std::uint32_t v;
        std::size_t n = varint_decode(p, v);
        if (p + n > end) break;
        p += n;
        out[i] = v;
    }
    *consumed = static_cast<std::size_t>(p - in);
    return i;
}

static void test_basic_one_byte_values() {
    // All single-byte varints (< 128).
    std::vector<std::uint32_t> values;
    for (std::uint32_t i = 0; i < 128; i++) values.push_back(i);
    auto buf = encode_stream(values);

    std::uint32_t out[16];
    std::size_t consumed = 0;
    int n = varint_decode_batch(buf.data(), buf.data() + buf.size(),
                                out, &consumed, 16);
    CHECK_EQ(n, 16);
    for (int i = 0; i < n; i++) CHECK_EQ(out[i], values[i]);
    CHECK_EQ(consumed, 16u);
}

static void test_mixed_widths() {
    // Mix of 1, 2, 3, 4, 5-byte varints.
    std::vector<std::uint32_t> values = {
        0, 1, 127, 128, 129, 255, 16383, 16384, 1u << 20,
        1u << 21, 1u << 27, 1u << 28, UINT32_MAX, 42, 999999, 1000
    };
    auto buf = encode_stream(values);

    std::uint32_t out[16];
    std::size_t consumed = 0;
    int n = varint_decode_batch(buf.data(), buf.data() + buf.size(),
                                out, &consumed, 16);
    CHECK_EQ(n, 16);
    for (int i = 0; i < n; i++) CHECK_EQ(out[i], values[i]);
    CHECK_EQ(consumed, buf.size());
}

static void test_partial_input() {
    // Truncate input mid-varint: should stop before consuming partial byte.
    std::vector<std::uint32_t> values = {1u << 20, 999, 7};
    auto buf = encode_stream(values);

    // Truncate by 1 byte.
    for (int trunc = 1; trunc <= 3 && trunc < (int)buf.size(); trunc++) {
        std::uint32_t out[8];
        std::size_t consumed = 0;
        int n = varint_decode_batch(buf.data(), buf.data() + buf.size() - trunc,
                                    out, &consumed, 8);
        // We may decode the first 1 or 2 fully and stop on the truncated last.
        for (int i = 0; i < n; i++) CHECK_EQ(out[i], values[i]);
        CHECK(consumed <= buf.size() - trunc);
    }
}

static void test_empty_input() {
    std::uint32_t out[8];
    std::size_t consumed = 0;
    std::uint8_t buf[1] = {0x00};
    int n = varint_decode_batch(buf, buf, out, &consumed, 8);
    CHECK_EQ(n, 0);
    CHECK_EQ(consumed, 0u);
}

static void test_large_random_stream() {
    // Generate ~4096 random values mixing typical posting-delta widths.
    std::mt19937 rng(42);
    std::vector<std::uint32_t> values;
    values.reserve(4096);
    for (int i = 0; i < 4096; i++) {
        std::uint32_t r = rng();
        switch (r & 0x3) {
            case 0: values.push_back(r & 0x7Fu); break;       // 1-byte
            case 1: values.push_back(r & 0x3FFFu); break;     // 2-byte
            case 2: values.push_back(r & 0x1FFFFFu); break;   // 3-byte
            case 3: values.push_back(r); break;                // up to 5-byte
        }
    }
    auto buf = encode_stream(values);

    // Decode in chunks of 16 and verify each chunk against scalar golden.
    const std::uint8_t* p = buf.data();
    const std::uint8_t* end = buf.data() + buf.size();
    std::size_t produced = 0;
    while (p < end && produced < values.size()) {
        std::uint32_t simd_out[16];
        std::uint32_t gold_out[16];
        std::size_t simd_consumed = 0;
        std::size_t gold_consumed = 0;
        int n_simd = varint_decode_batch(p, end, simd_out, &simd_consumed, 16);
        int n_gold = decode_stream_scalar(p, end, gold_out, &gold_consumed, 16);
        CHECK_EQ(n_simd, n_gold);
        CHECK_EQ(simd_consumed, gold_consumed);
        for (int i = 0; i < n_simd; i++) {
            CHECK_EQ(simd_out[i], gold_out[i]);
            CHECK_EQ(simd_out[i], values[produced + i]);
        }
        if (n_simd == 0) break;
        produced += n_simd;
        p += simd_consumed;
    }
    CHECK_EQ(produced, values.size());
}

static void test_size_boundaries() {
    // Test inputs of exactly 16, 32, 64, 65 bytes (per-tier boundaries).
    for (int total : {16, 32, 33, 63, 64, 65, 100}) {
        std::vector<std::uint32_t> values;
        // Pack with 1-byte varints first, then 2-byte to fill remainder.
        int budget = total;
        while (budget > 0) {
            if (budget >= 1 && (values.size() % 2 == 0)) {
                values.push_back(7);
                budget -= 1;
            } else if (budget >= 2) {
                values.push_back(200);
                budget -= 2;
            } else {
                break;
            }
        }
        auto buf = encode_stream(values);
        // Decode in chunks of 8 and accumulate.
        const std::uint8_t* p = buf.data();
        const std::uint8_t* end = buf.data() + buf.size();
        std::size_t produced = 0;
        while (p < end) {
            std::uint32_t out[8];
            std::size_t consumed = 0;
            int n = varint_decode_batch(p, end, out, &consumed, 8);
            if (n == 0) break;
            for (int i = 0; i < n; i++) CHECK_EQ(out[i], values[produced + i]);
            produced += n;
            p += consumed;
        }
        CHECK_EQ(produced, values.size());
    }
}

static void test_max_count_clamping() {
    // max_count > kVarintMaxBatch should be clamped internally.
    std::vector<std::uint32_t> values;
    for (int i = 0; i < 32; i++) values.push_back(static_cast<std::uint32_t>(i));
    auto buf = encode_stream(values);

    std::uint32_t out[32];
    std::size_t consumed = 0;
    int n = varint_decode_batch(buf.data(), buf.data() + buf.size(),
                                out, &consumed, 100);
    CHECK_EQ(n, kVarintMaxBatch);
    for (int i = 0; i < n; i++) CHECK_EQ(out[i], values[i]);
}

static void test_seq_id_decoder_next_batch_matches_next() {
    // Build a delta-encoded id posting: ids = [10, 10, 12, 12, 13, 100, ...]
    // Deltas:                           [10, 0,  2,  0,  1,  87,  ...]
    // was_new_seq:                     [T,  F,  T,  F,  T,  T,    ...]
    std::vector<std::uint32_t> deltas = {
        10, 0, 2, 0, 1, 87, 0, 0, 5, 200, 1u << 20, 0, 0, 7, 1, 3, 100, 0, 1
    };
    auto buf = encode_stream(deltas);

    // Reference: decode via next() one at a time.
    std::vector<std::uint32_t> ref_sids;
    std::vector<std::uint8_t>  ref_was_new;
    {
        SeqIdDecoder dec(buf.data(), buf.data() + buf.size());
        while (dec.has_more()) {
            std::uint32_t sid = dec.next();
            ref_sids.push_back(sid);
            ref_was_new.push_back(dec.was_new_seq() ? 1 : 0);
        }
    }
    CHECK_EQ(ref_sids.size(), deltas.size());

    // Test: decode via next_batch() in chunks of various sizes.
    for (int batch : {1, 2, 8, 15, 16}) {
        SeqIdDecoder dec(buf.data(), buf.data() + buf.size());
        std::vector<std::uint32_t> got_sids;
        std::vector<std::uint8_t>  got_was_new;
        std::uint32_t out_sids[SeqIdDecoder::kMaxBatch];
        std::uint8_t  out_new[SeqIdDecoder::kMaxBatch];
        while (dec.has_more()) {
            int n = dec.next_batch(out_sids, out_new, batch);
            if (n == 0) break;
            for (int i = 0; i < n; i++) {
                got_sids.push_back(out_sids[i]);
                got_was_new.push_back(out_new[i]);
            }
        }
        CHECK_EQ(got_sids.size(), ref_sids.size());
        for (size_t i = 0; i < got_sids.size(); i++) {
            CHECK_EQ(got_sids[i], ref_sids[i]);
            CHECK_EQ(got_was_new[i], ref_was_new[i]);
        }
    }
}

static void test_pos_decoder_next_batch_matches_next() {
    // Build a position posting with delta resets driven by was_new_seq.
    // Reset semantics: when was_new_seq[i] is true, prev_pos = raw[i] (absolute);
    // otherwise prev_pos += raw[i] (delta).
    std::vector<std::uint32_t> raws = {100, 50, 30, 20, 200, 5, 10};
    std::vector<std::uint8_t>  was_new = {1, 0, 0, 0, 1, 0, 0};
    auto buf = encode_stream(raws);

    // Reference via next() one at a time.
    std::vector<std::uint32_t> ref_pos;
    {
        PosDecoder dec(buf.data());
        for (size_t i = 0; i < raws.size(); i++) {
            ref_pos.push_back(dec.next(was_new[i] != 0));
        }
    }

    // Batch via next_batch(), using both end-bounded and end-null forms.
    for (bool with_end : {false, true}) {
        PosDecoder dec = with_end
            ? PosDecoder(buf.data(), buf.data() + buf.size())
            : PosDecoder(buf.data());
        std::vector<std::uint32_t> got;
        std::uint32_t out[PosDecoder::kMaxBatch];
        size_t i = 0;
        while (i < raws.size()) {
            int want = static_cast<int>(raws.size() - i);
            if (want > PosDecoder::kMaxBatch) want = PosDecoder::kMaxBatch;
            int n = dec.next_batch(out, was_new.data() + i, want);
            if (n == 0) break;
            for (int j = 0; j < n; j++) got.push_back(out[j]);
            i += n;
        }
        CHECK_EQ(got.size(), ref_pos.size());
        for (size_t j = 0; j < got.size(); j++) CHECK_EQ(got[j], ref_pos[j]);
    }
}

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();
    test_empty_input();
    test_basic_one_byte_values();
    test_mixed_widths();
    test_partial_input();
    test_large_random_stream();
    test_size_boundaries();
    test_max_count_clamping();
    // SeqIdDecoder / PosDecoder now consume PFor-encoded posting blobs (v4),
    // not varint streams.  Coverage for the new decoder path lives in
    // test_pfd_codec.cpp; the legacy varint-driven cases below exercise an
    // input encoding that is no longer supported.
    (void)test_seq_id_decoder_next_batch_matches_next;
    (void)test_pos_decoder_next_batch_matches_next;
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
