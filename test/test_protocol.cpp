#include "protocol/frame.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <unistd.h>

using namespace ikafssn;

// Helper: create a pipe pair for testing frame read/write
static bool make_pipe(int& read_fd, int& write_fd) {
    int fds[2];
    if (::pipe(fds) != 0) return false;
    read_fd = fds[0];
    write_fd = fds[1];
    return true;
}

static void test_frame_header_size() {
    std::printf("  test_frame_header_size...");
    assert(sizeof(FrameHeader) == 12);
    assert(FRAME_HEADER_SIZE == 12);
    std::printf(" OK\n");
}

static void test_frame_round_trip() {
    std::printf("  test_frame_round_trip...");

    int rfd, wfd;
    assert(make_pipe(rfd, wfd));

    std::vector<uint8_t> payload = {0x01, 0x02, 0x03, 0x04, 0x05};
    assert(write_frame(wfd, MsgType::kSearchRequest, payload));

    FrameHeader hdr;
    std::vector<uint8_t> recv_payload;
    assert(read_frame(rfd, hdr, recv_payload));

    assert(hdr.magic == FRAME_MAGIC);
    assert(hdr.payload_length == 5);
    assert(hdr.msg_type == static_cast<uint8_t>(MsgType::kSearchRequest));
    assert(hdr.msg_version == 1);
    assert(hdr.reserved == 0);
    assert(recv_payload == payload);

    ::close(rfd);
    ::close(wfd);
    std::printf(" OK\n");
}

static void test_frame_empty_payload() {
    std::printf("  test_frame_empty_payload...");

    int rfd, wfd;
    assert(make_pipe(rfd, wfd));

    std::vector<uint8_t> empty;
    assert(write_frame(wfd, MsgType::kHealthRequest, empty));

    FrameHeader hdr;
    std::vector<uint8_t> recv_payload;
    assert(read_frame(rfd, hdr, recv_payload));

    assert(hdr.payload_length == 0);
    assert(hdr.msg_type == static_cast<uint8_t>(MsgType::kHealthRequest));
    assert(recv_payload.empty());

    ::close(rfd);
    ::close(wfd);
    std::printf(" OK\n");
}

static void test_frame_invalid_magic() {
    std::printf("  test_frame_invalid_magic...");

    int rfd, wfd;
    assert(make_pipe(rfd, wfd));

    // Write a frame with bad magic
    FrameHeader bad_hdr;
    bad_hdr.magic = 0xDEADBEEF;
    bad_hdr.payload_length = 0;
    bad_hdr.msg_type = 0x01;
    bad_hdr.msg_version = 1;
    bad_hdr.reserved = 0;
    assert(write_all(wfd, &bad_hdr, sizeof(bad_hdr)));

    FrameHeader hdr;
    std::vector<uint8_t> payload;
    assert(!read_frame(rfd, hdr, payload));

    ::close(rfd);
    ::close(wfd);
    std::printf(" OK\n");
}

static void test_search_request_serialize() {
    std::printf("  test_search_request_serialize...");

    SearchRequest req;
    req.k = 9;
    req.min_score = 5;
    req.max_gap = 100;
    req.max_freq = 50000;
    req.min_diag_hits = 3;
    req.stage1_topn = 500;
    req.min_stage1_score = 2;
    req.num_results = 50;
    req.seqidlist_mode = SeqidlistMode::kInclude;
    req.seqids = {"NM_001234", "XM_005678"};
    req.queries.push_back({"query1", "ACGTACGTACGT"});
    req.queries.push_back({"query2", "TTTTAAAACCCC"});

    auto data = serialize(req);

    SearchRequest req2;
    assert(deserialize(data, req2));

    assert(req2.k == 9);
    assert(req2.min_score == 5);
    assert(req2.max_gap == 100);
    assert(req2.max_freq == 50000);
    assert(req2.min_diag_hits == 3);
    assert(req2.stage1_topn == 500);
    assert(req2.min_stage1_score == 2);
    assert(req2.num_results == 50);
    assert(req2.seqidlist_mode == SeqidlistMode::kInclude);
    assert(req2.seqids.size() == 2);
    assert(req2.seqids[0] == "NM_001234");
    assert(req2.seqids[1] == "XM_005678");
    assert(req2.queries.size() == 2);
    assert(req2.queries[0].query_id == "query1");
    assert(req2.queries[0].sequence == "ACGTACGTACGT");
    assert(req2.queries[1].query_id == "query2");
    assert(req2.queries[1].sequence == "TTTTAAAACCCC");

    std::printf(" OK\n");
}

static void test_search_request_defaults() {
    std::printf("  test_search_request_defaults...");

    SearchRequest req;
    // All defaults (zeros)
    req.queries.push_back({"q1", "ACGT"});

    auto data = serialize(req);
    SearchRequest req2;
    assert(deserialize(data, req2));

    assert(req2.k == 0);
    assert(req2.min_score == 0);
    assert(req2.max_gap == 0);
    assert(req2.max_freq == 0);
    assert(req2.seqidlist_mode == SeqidlistMode::kNone);
    assert(req2.seqids.empty());
    assert(req2.queries.size() == 1);

    std::printf(" OK\n");
}

static void test_search_response_serialize() {
    std::printf("  test_search_response_serialize...");

    SearchResponse resp;
    resp.status = 0;
    resp.k = 11;

    QueryResult qr;
    qr.query_id = "query1";

    ResponseHit hit;
    hit.accession = "NM_001234";
    hit.strand = 0; // '+'
    hit.q_start = 10;
    hit.q_end = 450;
    hit.s_start = 1020;
    hit.s_end = 1460;
    hit.score = 42;
    hit.volume = 0;
    qr.hits.push_back(hit);

    hit.accession = "XM_005678";
    hit.strand = 1; // '-'
    hit.q_start = 15;
    hit.q_end = 430;
    hit.s_start = 8050;
    hit.s_end = 8465;
    hit.score = 38;
    hit.volume = 2;
    qr.hits.push_back(hit);

    resp.results.push_back(qr);

    auto data = serialize(resp);
    SearchResponse resp2;
    assert(deserialize(data, resp2));

    assert(resp2.status == 0);
    assert(resp2.k == 11);
    assert(resp2.results.size() == 1);
    assert(resp2.results[0].query_id == "query1");
    assert(resp2.results[0].hits.size() == 2);

    const auto& h0 = resp2.results[0].hits[0];
    assert(h0.accession == "NM_001234");
    assert(h0.strand == 0);
    assert(h0.q_start == 10);
    assert(h0.q_end == 450);
    assert(h0.s_start == 1020);
    assert(h0.s_end == 1460);
    assert(h0.score == 42);
    assert(h0.volume == 0);

    const auto& h1 = resp2.results[0].hits[1];
    assert(h1.accession == "XM_005678");
    assert(h1.strand == 1);
    assert(h1.score == 38);
    assert(h1.volume == 2);

    std::printf(" OK\n");
}

static void test_error_response_serialize() {
    std::printf("  test_error_response_serialize...");

    ErrorResponse err;
    err.error_code = 404;
    err.message = "Index not found for k=15";

    auto data = serialize(err);
    ErrorResponse err2;
    assert(deserialize(data, err2));

    assert(err2.error_code == 404);
    assert(err2.message == "Index not found for k=15");

    std::printf(" OK\n");
}

static void test_health_roundtrip() {
    std::printf("  test_health_roundtrip...");

    HealthRequest hreq;
    auto req_data = serialize(hreq);
    assert(req_data.empty());

    HealthRequest hreq2;
    assert(deserialize(req_data, hreq2));

    HealthResponse hresp;
    hresp.status = 0;
    auto resp_data = serialize(hresp);
    assert(resp_data.size() == 1);

    HealthResponse hresp2;
    assert(deserialize(resp_data, hresp2));
    assert(hresp2.status == 0);

    std::printf(" OK\n");
}

static void test_search_request_with_seqidlist_exclude() {
    std::printf("  test_search_request_with_seqidlist_exclude...");

    SearchRequest req;
    req.k = 7;
    req.seqidlist_mode = SeqidlistMode::kExclude;
    req.seqids = {"ACC001", "ACC002", "ACC003"};
    req.queries.push_back({"q", "ACGTACG"});

    auto data = serialize(req);
    SearchRequest req2;
    assert(deserialize(data, req2));

    assert(req2.seqidlist_mode == SeqidlistMode::kExclude);
    assert(req2.seqids.size() == 3);
    assert(req2.seqids[2] == "ACC003");

    std::printf(" OK\n");
}

static void test_info_request_roundtrip() {
    std::printf("  test_info_request_roundtrip...");

    InfoRequest ireq;
    auto req_data = serialize(ireq);
    assert(req_data.empty());

    InfoRequest ireq2;
    assert(deserialize(req_data, ireq2));

    std::printf(" OK\n");
}

static void test_info_response_serialize() {
    std::printf("  test_info_response_serialize...");

    InfoResponse resp;
    resp.status = 0;
    resp.default_k = 11;

    KmerGroupInfo g1;
    g1.k = 7;
    g1.kmer_type = 0; // uint16
    VolumeInfo v1;
    v1.volume_index = 0;
    v1.num_sequences = 1000;
    v1.total_postings = 500000;
    v1.db_name = "testdb";
    g1.volumes.push_back(v1);
    VolumeInfo v2;
    v2.volume_index = 1;
    v2.num_sequences = 2000;
    v2.total_postings = 900000;
    v2.db_name = "testdb";
    g1.volumes.push_back(v2);
    resp.groups.push_back(g1);

    KmerGroupInfo g2;
    g2.k = 11;
    g2.kmer_type = 1; // uint32
    VolumeInfo v3;
    v3.volume_index = 0;
    v3.num_sequences = 1000;
    v3.total_postings = 450000;
    v3.db_name = "testdb";
    g2.volumes.push_back(v3);
    resp.groups.push_back(g2);

    auto data = serialize(resp);
    InfoResponse resp2;
    assert(deserialize(data, resp2));

    assert(resp2.status == 0);
    assert(resp2.default_k == 11);
    assert(resp2.groups.size() == 2);

    assert(resp2.groups[0].k == 7);
    assert(resp2.groups[0].kmer_type == 0);
    assert(resp2.groups[0].volumes.size() == 2);
    assert(resp2.groups[0].volumes[0].volume_index == 0);
    assert(resp2.groups[0].volumes[0].num_sequences == 1000);
    assert(resp2.groups[0].volumes[0].total_postings == 500000);
    assert(resp2.groups[0].volumes[0].db_name == "testdb");
    assert(resp2.groups[0].volumes[1].volume_index == 1);
    assert(resp2.groups[0].volumes[1].num_sequences == 2000);
    assert(resp2.groups[0].volumes[1].total_postings == 900000);

    assert(resp2.groups[1].k == 11);
    assert(resp2.groups[1].kmer_type == 1);
    assert(resp2.groups[1].volumes.size() == 1);
    assert(resp2.groups[1].volumes[0].total_postings == 450000);

    std::printf(" OK\n");
}

static void test_info_response_empty() {
    std::printf("  test_info_response_empty...");

    InfoResponse resp;
    resp.status = 0;
    resp.default_k = 9;
    // No groups

    auto data = serialize(resp);
    InfoResponse resp2;
    assert(deserialize(data, resp2));

    assert(resp2.status == 0);
    assert(resp2.default_k == 9);
    assert(resp2.groups.empty());

    std::printf(" OK\n");
}

static void test_frame_full_round_trip() {
    std::printf("  test_frame_full_round_trip...");

    int rfd, wfd;
    assert(make_pipe(rfd, wfd));

    // Serialize a SearchRequest and send as frame
    SearchRequest req;
    req.k = 9;
    req.queries.push_back({"q1", "ACGTACGT"});
    auto payload = serialize(req);

    assert(write_frame(wfd, MsgType::kSearchRequest, payload));

    // Read frame
    FrameHeader hdr;
    std::vector<uint8_t> recv_payload;
    assert(read_frame(rfd, hdr, recv_payload));
    assert(hdr.msg_type == static_cast<uint8_t>(MsgType::kSearchRequest));

    // Deserialize
    SearchRequest req2;
    assert(deserialize(recv_payload, req2));
    assert(req2.k == 9);
    assert(req2.queries[0].sequence == "ACGTACGT");

    ::close(rfd);
    ::close(wfd);
    std::printf(" OK\n");
}

static void test_search_request_accept_qdegen() {
    std::printf("  test_search_request_accept_qdegen...");

    SearchRequest req;
    req.k = 9;
    req.accept_qdegen = 1;
    req.queries.push_back({"q1", "ACGTRYSWKM"});

    auto data = serialize(req);
    SearchRequest req2;
    assert(deserialize(data, req2));

    assert(req2.accept_qdegen == 1);

    // Also test with accept_qdegen = 0
    req.accept_qdegen = 0;
    data = serialize(req);
    SearchRequest req3;
    assert(deserialize(data, req3));
    assert(req3.accept_qdegen == 0);

    std::printf(" OK\n");
}

static void test_search_response_skipped() {
    std::printf("  test_search_response_skipped...");

    SearchResponse resp;
    resp.status = 0;
    resp.k = 9;

    QueryResult qr1;
    qr1.query_id = "q1";
    qr1.skipped = 1;
    resp.results.push_back(qr1);

    QueryResult qr2;
    qr2.query_id = "q2";
    qr2.skipped = 0;
    ResponseHit hit;
    hit.accession = "ACC001";
    hit.strand = 0;
    hit.q_start = 0;
    hit.q_end = 100;
    hit.s_start = 0;
    hit.s_end = 100;
    hit.score = 10;
    hit.stage1_score = 5;
    hit.volume = 0;
    qr2.hits.push_back(hit);
    resp.results.push_back(qr2);

    auto data = serialize(resp);
    SearchResponse resp2;
    assert(deserialize(data, resp2));

    assert(resp2.results.size() == 2);
    assert(resp2.results[0].query_id == "q1");
    assert(resp2.results[0].skipped == 1);
    assert(resp2.results[0].hits.empty());
    assert(resp2.results[1].query_id == "q2");
    assert(resp2.results[1].skipped == 0);
    assert(resp2.results[1].hits.size() == 1);

    std::printf(" OK\n");
}

static void test_search_request_backward_compat_no_accept_qdegen() {
    std::printf("  test_search_request_backward_compat_no_accept_qdegen...");

    // Simulate old format: serialize without the accept_qdegen trailer
    // by serializing and then trimming the last byte
    SearchRequest req;
    req.k = 9;
    req.accept_qdegen = 0;
    req.queries.push_back({"q1", "ACGT"});

    auto data = serialize(req);
    // Remove the last byte (accept_qdegen trailer)
    data.pop_back();

    SearchRequest req2;
    assert(deserialize(data, req2));
    // Should default to 0
    assert(req2.accept_qdegen == 0);

    std::printf(" OK\n");
}

static void test_search_response_backward_compat_no_skipped() {
    std::printf("  test_search_response_backward_compat_no_skipped...");

    // Simulate old format: serialize a response, remove the skipped trailer
    SearchResponse resp;
    resp.status = 0;
    resp.k = 9;

    QueryResult qr;
    qr.query_id = "q1";
    qr.skipped = 1; // Will be serialized
    resp.results.push_back(qr);

    auto data = serialize(resp);
    // Remove the last byte (skipped trailer for 1 query)
    data.pop_back();

    SearchResponse resp2;
    assert(deserialize(data, resp2));
    // Should default to 0
    assert(resp2.results[0].skipped == 0);

    std::printf(" OK\n");
}

int main() {
    std::printf("test_protocol:\n");

    test_frame_header_size();
    test_frame_round_trip();
    test_frame_empty_payload();
    test_frame_invalid_magic();
    test_search_request_serialize();
    test_search_request_defaults();
    test_search_response_serialize();
    test_error_response_serialize();
    test_health_roundtrip();
    test_search_request_with_seqidlist_exclude();
    test_info_request_roundtrip();
    test_info_response_serialize();
    test_info_response_empty();
    test_frame_full_round_trip();
    test_search_request_accept_qdegen();
    test_search_response_skipped();
    test_search_request_backward_compat_no_accept_qdegen();
    test_search_response_backward_compat_no_skipped();

    std::printf("All protocol tests passed.\n");
    return 0;
}
