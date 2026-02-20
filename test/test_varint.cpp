#include "test_util.hpp"
#include "core/varint.hpp"
#include <vector>

using namespace ikafssn;

static void test_roundtrip(uint32_t value) {
    uint8_t buf[5];
    size_t written = varint_encode(value, buf);
    uint32_t decoded;
    size_t consumed = varint_decode(buf, decoded);
    CHECK_EQ(decoded, value);
    CHECK_EQ(consumed, written);
    CHECK_EQ(varint_size(value), written);
}

static void test_basic_values() {
    test_roundtrip(0);
    test_roundtrip(1);
    test_roundtrip(127);
    test_roundtrip(128);
    test_roundtrip(16383);
    test_roundtrip(16384);
    test_roundtrip(UINT32_MAX);
}

static void test_encoded_sizes() {
    CHECK_EQ(varint_size(0), 1u);
    CHECK_EQ(varint_size(1), 1u);
    CHECK_EQ(varint_size(127), 1u);
    CHECK_EQ(varint_size(128), 2u);
    CHECK_EQ(varint_size(16383), 2u);
    CHECK_EQ(varint_size(16384), 3u);
    CHECK_EQ(varint_size(UINT32_MAX), 5u);
}

static void test_stream_decode() {
    // Encode multiple values into a stream, then decode sequentially
    std::vector<uint32_t> values = {0, 1, 127, 128, 255, 16384, 1000000, UINT32_MAX};
    std::vector<uint8_t> buf(values.size() * 5);

    size_t total_written = 0;
    for (auto v : values) {
        total_written += varint_encode(v, buf.data() + total_written);
    }

    size_t offset = 0;
    for (size_t i = 0; i < values.size(); i++) {
        uint32_t decoded;
        offset += varint_decode(buf.data() + offset, decoded);
        CHECK_EQ(decoded, values[i]);
    }
    CHECK_EQ(offset, total_written);
}

int main() {
    test_basic_values();
    test_encoded_sizes();
    test_stream_decode();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
