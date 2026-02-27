#include "protocol/serializer.hpp"

#include <cstring>
#include <limits>

namespace ikafssn {

// Helper: append little-endian integer to buffer
static void put_u8(std::vector<uint8_t>& buf, uint8_t v) {
    buf.push_back(v);
}

static void put_i8(std::vector<uint8_t>& buf, int8_t v) {
    buf.push_back(static_cast<uint8_t>(v));
}

static void put_u16(std::vector<uint8_t>& buf, uint16_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
}

static void put_i16(std::vector<uint8_t>& buf, int16_t v) {
    put_u16(buf, static_cast<uint16_t>(v));
}

static void put_u32(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
    buf.push_back(static_cast<uint8_t>(v >> 16));
    buf.push_back(static_cast<uint8_t>(v >> 24));
}

static void put_i32(std::vector<uint8_t>& buf, int32_t v) {
    put_u32(buf, static_cast<uint32_t>(v));
}

static void put_u64(std::vector<uint8_t>& buf, uint64_t v) {
    for (int i = 0; i < 8; i++) {
        buf.push_back(static_cast<uint8_t>(v >> (i * 8)));
    }
}

static void put_str16(std::vector<uint8_t>& buf, const std::string& s) {
    put_u16(buf, static_cast<uint16_t>(s.size()));
    buf.insert(buf.end(), s.begin(), s.end());
}

// Helper: read little-endian integer from buffer
class Reader {
public:
    Reader(const uint8_t* data, size_t size) : data_(data), size_(size), pos_(0) {}

    bool has(size_t n) const { return pos_ + n <= size_; }
    size_t remaining() const { return size_ - pos_; }

    bool get_u8(uint8_t& v) {
        if (!has(1)) return false;
        v = data_[pos_++];
        return true;
    }

    bool get_i8(int8_t& v) {
        uint8_t raw;
        if (!get_u8(raw)) return false;
        v = static_cast<int8_t>(raw);
        return true;
    }

    bool get_u16(uint16_t& v) {
        if (!has(2)) return false;
        v = static_cast<uint16_t>(data_[pos_]) |
            (static_cast<uint16_t>(data_[pos_ + 1]) << 8);
        pos_ += 2;
        return true;
    }

    bool get_i16(int16_t& v) {
        uint16_t raw;
        if (!get_u16(raw)) return false;
        v = static_cast<int16_t>(raw);
        return true;
    }

    bool get_u32(uint32_t& v) {
        if (!has(4)) return false;
        v = static_cast<uint32_t>(data_[pos_]) |
            (static_cast<uint32_t>(data_[pos_ + 1]) << 8) |
            (static_cast<uint32_t>(data_[pos_ + 2]) << 16) |
            (static_cast<uint32_t>(data_[pos_ + 3]) << 24);
        pos_ += 4;
        return true;
    }

    bool get_i32(int32_t& v) {
        uint32_t raw;
        if (!get_u32(raw)) return false;
        v = static_cast<int32_t>(raw);
        return true;
    }

    bool get_u64(uint64_t& v) {
        if (!has(8)) return false;
        v = 0;
        for (int i = 0; i < 8; i++) {
            v |= static_cast<uint64_t>(data_[pos_ + i]) << (i * 8);
        }
        pos_ += 8;
        return true;
    }

    bool get_str16(std::string& s) {
        uint16_t len;
        if (!get_u16(len)) return false;
        if (!has(len)) return false;
        s.assign(reinterpret_cast<const char*>(data_ + pos_), len);
        pos_ += len;
        return true;
    }

    bool skip(size_t n) {
        if (!has(n)) return false;
        pos_ += n;
        return true;
    }

private:
    const uint8_t* data_;
    size_t size_;
    size_t pos_;
};

// --- SearchRequest ---
// Wire format (all fields in natural order, no backward-compat trailer):
//   u8   k
//   u16  stage2_min_score
//   u16  stage2_max_gap
//   u32  stage1_max_freq
//   u8   stage2_min_diag_hits
//   u16  stage1_topn
//   u16  stage1_min_score
//   u16  num_results
//   u16  stage1_min_score_frac_x10000
//   u16  stage1_max_freq_frac_x10000
//   u8   seqidlist_mode
//   u8   mode
//   u8   stage1_score
//   u8   accept_qdegen
//   i8   strand
//   u8   has_stage2_min_score
//   u16  stage2_max_lookback
//   u8   stage3_traceback
//   i16  stage3_gapopen
//   i16  stage3_gapext
//   u16  stage3_min_pident_x100
//   u32  stage3_min_nident
//   u32  context_abs
//   u16  context_frac_x10000
//   str16 db
//   u32  num_seqids
//     [str16 seqid] × num_seqids
//   u16  num_queries
//     [str16 query_id, u32 seq_len, bytes seq] × num_queries

std::vector<uint8_t> serialize(const SearchRequest& req) {
    std::vector<uint8_t> buf;
    buf.reserve(256);

    put_u8(buf, req.k);
    put_u16(buf, req.stage2_min_score);
    put_u16(buf, req.stage2_max_gap);
    put_u32(buf, req.stage1_max_freq);
    put_u8(buf, req.stage2_min_diag_hits);
    put_u16(buf, req.stage1_topn);
    put_u16(buf, req.stage1_min_score);
    put_u16(buf, req.num_results);
    put_u16(buf, req.stage1_min_score_frac_x10000);
    put_u16(buf, req.stage1_max_freq_frac_x10000);
    put_u8(buf, static_cast<uint8_t>(req.seqidlist_mode));
    put_u8(buf, req.mode);
    put_u8(buf, req.stage1_score);
    put_u8(buf, req.accept_qdegen);
    put_i8(buf, req.strand);
    put_u8(buf, req.has_stage2_min_score);
    put_u16(buf, req.stage2_max_lookback);
    put_u8(buf, req.stage3_traceback);
    put_i16(buf, req.stage3_gapopen);
    put_i16(buf, req.stage3_gapext);
    put_u16(buf, req.stage3_min_pident_x100);
    put_u32(buf, req.stage3_min_nident);
    put_u32(buf, req.context_abs);
    put_u16(buf, req.context_frac_x10000);
    put_str16(buf, req.db);

    put_u32(buf, static_cast<uint32_t>(req.seqids.size()));
    for (const auto& acc : req.seqids) {
        put_str16(buf, acc);
    }

    put_u16(buf, static_cast<uint16_t>(req.queries.size()));
    for (const auto& q : req.queries) {
        put_str16(buf, q.query_id);
        put_u32(buf, static_cast<uint32_t>(q.sequence.size()));
        buf.insert(buf.end(), q.sequence.begin(), q.sequence.end());
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, SearchRequest& req) {
    Reader r(data.data(), data.size());

    if (!r.get_u8(req.k)) return false;
    if (!r.get_u16(req.stage2_min_score)) return false;
    if (!r.get_u16(req.stage2_max_gap)) return false;
    if (!r.get_u32(req.stage1_max_freq)) return false;
    if (!r.get_u8(req.stage2_min_diag_hits)) return false;
    if (!r.get_u16(req.stage1_topn)) return false;
    if (!r.get_u16(req.stage1_min_score)) return false;
    if (!r.get_u16(req.num_results)) return false;
    if (!r.get_u16(req.stage1_min_score_frac_x10000)) return false;
    if (!r.get_u16(req.stage1_max_freq_frac_x10000)) return false;

    uint8_t seqidlist_mode;
    if (!r.get_u8(seqidlist_mode)) return false;
    if (seqidlist_mode > 2) return false;
    req.seqidlist_mode = static_cast<SeqidlistMode>(seqidlist_mode);

    if (!r.get_u8(req.mode)) return false;
    if (!r.get_u8(req.stage1_score)) return false;
    if (!r.get_u8(req.accept_qdegen)) return false;
    if (!r.get_i8(req.strand)) return false;
    if (!r.get_u8(req.has_stage2_min_score)) return false;
    if (!r.get_u16(req.stage2_max_lookback)) return false;
    if (!r.get_u8(req.stage3_traceback)) return false;
    if (!r.get_i16(req.stage3_gapopen)) return false;
    if (!r.get_i16(req.stage3_gapext)) return false;
    if (!r.get_u16(req.stage3_min_pident_x100)) return false;
    if (!r.get_u32(req.stage3_min_nident)) return false;
    if (!r.get_u32(req.context_abs)) return false;
    if (!r.get_u16(req.context_frac_x10000)) return false;
    if (!r.get_str16(req.db)) return false;

    uint32_t num_seqids;
    if (!r.get_u32(num_seqids)) return false;
    req.seqids.resize(num_seqids);
    for (uint32_t i = 0; i < num_seqids; i++) {
        if (!r.get_str16(req.seqids[i])) return false;
    }

    uint16_t num_queries;
    if (!r.get_u16(num_queries)) return false;
    req.queries.resize(num_queries);
    for (uint16_t i = 0; i < num_queries; i++) {
        if (!r.get_str16(req.queries[i].query_id)) return false;
        uint32_t seq_len;
        if (!r.get_u32(seq_len)) return false;
        if (!r.has(seq_len)) return false;
        req.queries[i].sequence.assign(
            reinterpret_cast<const char*>(data.data() + (data.size() - r.remaining())),
            seq_len);
        if (!r.skip(seq_len)) return false;
    }

    return true;
}

// --- SearchResponse ---
// Wire format:
//   u8   status
//   u8   k
//   u8   mode
//   u8   stage1_score
//   u8   stage3_traceback
//   str16 db
//   u16  num_queries
//   for each query:
//     str16 query_id
//     u8    skipped
//     u8    warnings
//     u16   num_hits
//     for each hit:
//       str16  accession
//       u8     strand
//       u32    q_start
//       u32    q_end
//       u32    q_length
//       u32    s_start
//       u32    s_end
//       u32    s_length
//       u16    coverscore
//       u16    matchscore
//       u16    chainscore
//       u16    volume
//       i32    alnscore
//       u32    nident
//       u32    nmismatch
//       u16    pident_x100
//       str16  cigar
//       str16  q_seq
//       str16  s_seq
//   u16  num_rejected
//     [str16 query_id] × num_rejected

std::vector<uint8_t> serialize(const SearchResponse& resp) {
    std::vector<uint8_t> buf;
    buf.reserve(1024);

    put_u8(buf, resp.status);
    put_u8(buf, resp.k);
    put_u8(buf, resp.mode);
    put_u8(buf, resp.stage1_score);
    put_u8(buf, resp.stage3_traceback);
    put_str16(buf, resp.db);
    put_u16(buf, static_cast<uint16_t>(resp.results.size()));

    for (const auto& qr : resp.results) {
        put_str16(buf, qr.query_id);
        put_u8(buf, qr.skipped);
        put_u8(buf, qr.warnings);
        put_u16(buf, static_cast<uint16_t>(qr.hits.size()));
        for (const auto& hit : qr.hits) {
            put_str16(buf, hit.accession);
            put_u8(buf, hit.strand);
            put_u32(buf, hit.q_start);
            put_u32(buf, hit.q_end);
            put_u32(buf, hit.q_length);
            put_u32(buf, hit.s_start);
            put_u32(buf, hit.s_end);
            put_u32(buf, hit.s_length);
            put_u16(buf, hit.coverscore);
            put_u16(buf, hit.matchscore);
            put_u16(buf, hit.chainscore);
            put_u16(buf, hit.volume);
            put_i32(buf, hit.alnscore);
            put_u32(buf, hit.nident);
            put_u32(buf, hit.nmismatch);
            put_u16(buf, hit.pident_x100);
            put_str16(buf, hit.cigar);
            put_str16(buf, hit.q_seq);
            put_str16(buf, hit.s_seq);
        }
    }

    // Rejected query IDs
    put_u16(buf, static_cast<uint16_t>(resp.rejected_query_ids.size()));
    for (const auto& qid : resp.rejected_query_ids) {
        put_str16(buf, qid);
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, SearchResponse& resp) {
    Reader r(data.data(), data.size());

    if (!r.get_u8(resp.status)) return false;
    if (!r.get_u8(resp.k)) return false;
    if (!r.get_u8(resp.mode)) return false;
    if (!r.get_u8(resp.stage1_score)) return false;
    if (!r.get_u8(resp.stage3_traceback)) return false;
    if (!r.get_str16(resp.db)) return false;

    uint16_t num_queries;
    if (!r.get_u16(num_queries)) return false;
    resp.results.resize(num_queries);

    for (uint16_t qi = 0; qi < num_queries; qi++) {
        auto& qr = resp.results[qi];
        if (!r.get_str16(qr.query_id)) return false;
        if (!r.get_u8(qr.skipped)) return false;
        if (!r.get_u8(qr.warnings)) return false;

        uint16_t num_hits;
        if (!r.get_u16(num_hits)) return false;
        qr.hits.resize(num_hits);

        for (uint16_t hi = 0; hi < num_hits; hi++) {
            auto& hit = qr.hits[hi];
            if (!r.get_str16(hit.accession)) return false;
            if (!r.get_u8(hit.strand)) return false;
            if (!r.get_u32(hit.q_start)) return false;
            if (!r.get_u32(hit.q_end)) return false;
            if (!r.get_u32(hit.q_length)) return false;
            if (!r.get_u32(hit.s_start)) return false;
            if (!r.get_u32(hit.s_end)) return false;
            if (!r.get_u32(hit.s_length)) return false;
            if (!r.get_u16(hit.coverscore)) return false;
            if (!r.get_u16(hit.matchscore)) return false;
            if (!r.get_u16(hit.chainscore)) return false;
            if (!r.get_u16(hit.volume)) return false;
            if (!r.get_i32(hit.alnscore)) return false;
            if (!r.get_u32(hit.nident)) return false;
            if (!r.get_u32(hit.nmismatch)) return false;
            if (!r.get_u16(hit.pident_x100)) return false;
            if (!r.get_str16(hit.cigar)) return false;
            if (!r.get_str16(hit.q_seq)) return false;
            if (!r.get_str16(hit.s_seq)) return false;
        }
    }

    // Rejected query IDs
    uint16_t num_rejected;
    if (!r.get_u16(num_rejected)) return false;
    resp.rejected_query_ids.resize(num_rejected);
    for (uint16_t i = 0; i < num_rejected; i++) {
        if (!r.get_str16(resp.rejected_query_ids[i])) return false;
    }

    return true;
}

// --- ErrorResponse ---

std::vector<uint8_t> serialize(const ErrorResponse& err) {
    std::vector<uint8_t> buf;
    put_u32(buf, err.error_code);
    put_str16(buf, err.message);
    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, ErrorResponse& err) {
    Reader r(data.data(), data.size());
    if (!r.get_u32(err.error_code)) return false;
    if (!r.get_str16(err.message)) return false;
    return true;
}

// --- HealthRequest ---

std::vector<uint8_t> serialize(const HealthRequest& /*req*/) {
    return {};
}

bool deserialize(const std::vector<uint8_t>& /*data*/, HealthRequest& /*req*/) {
    return true;
}

// --- HealthResponse ---

std::vector<uint8_t> serialize(const HealthResponse& resp) {
    std::vector<uint8_t> buf;
    put_u8(buf, resp.status);
    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, HealthResponse& resp) {
    Reader r(data.data(), data.size());
    if (!r.get_u8(resp.status)) return false;
    return true;
}

// --- InfoRequest ---

std::vector<uint8_t> serialize(const InfoRequest& /*req*/) {
    return {};
}

bool deserialize(const std::vector<uint8_t>& /*data*/, InfoRequest& /*req*/) {
    return true;
}

// --- InfoResponse ---
// Wire format (v4):
//   u8   status
//   u8   default_k
//   i32  max_queue_size
//   i32  queue_depth
//   i32  max_seqs_per_req
//   u16  num_databases
//   for each database:
//     str16 name
//     u8    default_k
//     u8    max_mode
//     u16   num_groups
//     for each group:
//       u8  k
//       u8  kmer_type
//       u16 num_volumes
//       for each volume:
//         u16 volume_index
//         u32 num_sequences
//         u64 total_postings
//         u64 total_bases
//         str16 db

std::vector<uint8_t> serialize(const InfoResponse& resp) {
    std::vector<uint8_t> buf;
    buf.reserve(256);

    put_u8(buf, resp.status);
    put_u8(buf, resp.default_k);
    put_i32(buf, resp.max_queue_size);
    put_i32(buf, resp.queue_depth);
    put_i32(buf, resp.max_seqs_per_req);
    put_u16(buf, static_cast<uint16_t>(resp.databases.size()));

    for (const auto& db : resp.databases) {
        put_str16(buf, db.name);
        put_u8(buf, db.default_k);
        put_u8(buf, db.max_mode);
        put_u16(buf, static_cast<uint16_t>(db.groups.size()));

        for (const auto& g : db.groups) {
            put_u8(buf, g.k);
            put_u8(buf, g.kmer_type);
            put_u16(buf, static_cast<uint16_t>(g.volumes.size()));

            for (const auto& v : g.volumes) {
                put_u16(buf, v.volume_index);
                put_u32(buf, v.num_sequences);
                put_u64(buf, v.total_postings);
                put_u64(buf, v.total_bases);
                put_str16(buf, v.db);
            }
        }
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, InfoResponse& resp) {
    Reader r(data.data(), data.size());

    if (!r.get_u8(resp.status)) return false;
    if (!r.get_u8(resp.default_k)) return false;
    if (!r.get_i32(resp.max_queue_size)) return false;
    if (!r.get_i32(resp.queue_depth)) return false;
    if (!r.get_i32(resp.max_seqs_per_req)) return false;

    uint16_t num_databases;
    if (!r.get_u16(num_databases)) return false;
    resp.databases.resize(num_databases);

    for (uint16_t di = 0; di < num_databases; di++) {
        auto& db = resp.databases[di];
        if (!r.get_str16(db.name)) return false;
        if (!r.get_u8(db.default_k)) return false;
        if (!r.get_u8(db.max_mode)) return false;

        uint16_t num_groups;
        if (!r.get_u16(num_groups)) return false;
        db.groups.resize(num_groups);

        for (uint16_t gi = 0; gi < num_groups; gi++) {
            auto& g = db.groups[gi];
            if (!r.get_u8(g.k)) return false;
            if (!r.get_u8(g.kmer_type)) return false;

            uint16_t num_vols;
            if (!r.get_u16(num_vols)) return false;
            g.volumes.resize(num_vols);

            for (uint16_t vi = 0; vi < num_vols; vi++) {
                auto& v = g.volumes[vi];
                if (!r.get_u16(v.volume_index)) return false;
                if (!r.get_u32(v.num_sequences)) return false;
                if (!r.get_u64(v.total_postings)) return false;
                if (!r.get_u64(v.total_bases)) return false;
                if (!r.get_str16(v.db)) return false;
            }
        }
    }

    return true;
}

} // namespace ikafssn
