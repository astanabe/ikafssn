#include "protocol/serializer.hpp"

#include <cstring>
#include <limits>

namespace ikafssn {

// Helper: append little-endian integer to buffer
static void put_u8(std::vector<uint8_t>& buf, uint8_t v) {
    buf.push_back(v);
}

static void put_u16(std::vector<uint8_t>& buf, uint16_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
}

static void put_u32(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
    buf.push_back(static_cast<uint8_t>(v >> 16));
    buf.push_back(static_cast<uint8_t>(v >> 24));
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

    bool get_u16(uint16_t& v) {
        if (!has(2)) return false;
        v = static_cast<uint16_t>(data_[pos_]) |
            (static_cast<uint16_t>(data_[pos_ + 1]) << 8);
        pos_ += 2;
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

std::vector<uint8_t> serialize(const SearchRequest& req) {
    std::vector<uint8_t> buf;
    buf.reserve(256);

    put_u8(buf, req.k);
    put_u16(buf, req.min_score);
    put_u16(buf, req.max_gap);
    put_u32(buf, req.max_freq);
    put_u8(buf, req.min_diag_hits);
    put_u16(buf, req.stage1_topn);
    put_u16(buf, req.min_stage1_score);
    put_u16(buf, req.num_results);
    put_u8(buf, static_cast<uint8_t>(req.seqidlist_mode));
    put_u8(buf, req.mode);
    put_u8(buf, static_cast<uint8_t>(req.stage1_score_type | (req.sort_score << 4)));

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

    // Backward-compatible trailer: fractional min_stage1_score
    put_u16(buf, req.min_stage1_score_frac_x10000);

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, SearchRequest& req) {
    Reader r(data.data(), data.size());

    if (!r.get_u8(req.k)) return false;
    if (!r.get_u16(req.min_score)) return false;
    if (!r.get_u16(req.max_gap)) return false;
    if (!r.get_u32(req.max_freq)) return false;
    if (!r.get_u8(req.min_diag_hits)) return false;
    if (!r.get_u16(req.stage1_topn)) return false;
    if (!r.get_u16(req.min_stage1_score)) return false;
    if (!r.get_u16(req.num_results)) return false;

    uint8_t seqidlist_mode;
    if (!r.get_u8(seqidlist_mode)) return false;
    if (seqidlist_mode > 2) return false;
    req.seqidlist_mode = static_cast<SeqidlistMode>(seqidlist_mode);

    if (!r.get_u8(req.mode)) return false;
    uint8_t packed;
    if (!r.get_u8(packed)) return false;
    req.stage1_score_type = packed & 0x0F;
    req.sort_score = (packed >> 4) & 0x0F;

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

    // Backward-compatible trailer: fractional min_stage1_score
    if (r.remaining() >= 2) {
        r.get_u16(req.min_stage1_score_frac_x10000);
    }

    return true;
}

// --- SearchResponse ---

std::vector<uint8_t> serialize(const SearchResponse& resp) {
    std::vector<uint8_t> buf;
    buf.reserve(1024);

    put_u8(buf, resp.status);
    put_u8(buf, resp.k);
    put_u8(buf, resp.mode);
    put_u8(buf, resp.stage1_score_type);
    put_u16(buf, static_cast<uint16_t>(resp.results.size()));

    for (const auto& qr : resp.results) {
        put_str16(buf, qr.query_id);
        put_u16(buf, static_cast<uint16_t>(qr.hits.size()));
        for (const auto& hit : qr.hits) {
            put_str16(buf, hit.accession);
            put_u8(buf, hit.strand);
            put_u32(buf, hit.q_start);
            put_u32(buf, hit.q_end);
            put_u32(buf, hit.s_start);
            put_u32(buf, hit.s_end);
            put_u16(buf, hit.score);
            put_u16(buf, hit.stage1_score);
            put_u16(buf, hit.volume);
        }
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, SearchResponse& resp) {
    Reader r(data.data(), data.size());

    if (!r.get_u8(resp.status)) return false;
    if (!r.get_u8(resp.k)) return false;
    if (!r.get_u8(resp.mode)) return false;
    if (!r.get_u8(resp.stage1_score_type)) return false;

    uint16_t num_queries;
    if (!r.get_u16(num_queries)) return false;
    resp.results.resize(num_queries);

    for (uint16_t qi = 0; qi < num_queries; qi++) {
        auto& qr = resp.results[qi];
        if (!r.get_str16(qr.query_id)) return false;

        uint16_t num_hits;
        if (!r.get_u16(num_hits)) return false;
        qr.hits.resize(num_hits);

        for (uint16_t hi = 0; hi < num_hits; hi++) {
            auto& hit = qr.hits[hi];
            if (!r.get_str16(hit.accession)) return false;
            if (!r.get_u8(hit.strand)) return false;
            if (!r.get_u32(hit.q_start)) return false;
            if (!r.get_u32(hit.q_end)) return false;
            if (!r.get_u32(hit.s_start)) return false;
            if (!r.get_u32(hit.s_end)) return false;
            if (!r.get_u16(hit.score)) return false;
            if (!r.get_u16(hit.stage1_score)) return false;
            if (!r.get_u16(hit.volume)) return false;
        }
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
// Wire format:
//   u8  status
//   u8  default_k
//   u16 num_groups
//   for each group:
//     u8  k
//     u8  kmer_type
//     u16 num_volumes
//     for each volume:
//       u16 volume_index
//       u32 num_sequences
//       u64 total_postings
//       str16 db_name

std::vector<uint8_t> serialize(const InfoResponse& resp) {
    std::vector<uint8_t> buf;
    buf.reserve(256);

    put_u8(buf, resp.status);
    put_u8(buf, resp.default_k);
    put_u16(buf, static_cast<uint16_t>(resp.groups.size()));

    for (const auto& g : resp.groups) {
        put_u8(buf, g.k);
        put_u8(buf, g.kmer_type);
        put_u16(buf, static_cast<uint16_t>(g.volumes.size()));

        for (const auto& v : g.volumes) {
            put_u16(buf, v.volume_index);
            put_u32(buf, v.num_sequences);
            put_u64(buf, v.total_postings);
            put_str16(buf, v.db_name);
        }
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, InfoResponse& resp) {
    Reader r(data.data(), data.size());

    if (!r.get_u8(resp.status)) return false;
    if (!r.get_u8(resp.default_k)) return false;

    uint16_t num_groups;
    if (!r.get_u16(num_groups)) return false;
    resp.groups.resize(num_groups);

    for (uint16_t gi = 0; gi < num_groups; gi++) {
        auto& g = resp.groups[gi];
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
            if (!r.get_str16(v.db_name)) return false;
        }
    }

    return true;
}

} // namespace ikafssn
