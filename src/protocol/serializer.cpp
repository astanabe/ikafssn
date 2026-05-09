#include "protocol/serializer.hpp"
#include "protocol/binary_io.hpp"

#include <cstring>
#include <limits>

namespace ikafssn {

// --- SearchRequest ---
// Wire format (all fields in natural order, no backward-compat trailer):
//   u8   k
//   u16  stage2_min_score
//   u16  stage2_max_gap
//   u8   stage2_min_nhit_diag
//   u16  stage1_topn
//   u16  stage1_min_score
//   u16  nresult
//   u16  stage1_min_score_frac_x10000
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
//   u16  stage3_min_ppositive_x100
//   u32  stage3_min_npositive
//   u32  context_abs
//   u16  context_frac_x10000
//   u16  max_degen_expand
//   u16  stage2_max_nhit_per_subject
//   u8   t
//   u8   template_type
//   u8   score_matrix
//   str16 db
//   u32  num_seqids
//     [str16 seqid] × num_seqids
//   u16  num_queries
//     [str16 qseqid, u32 seq_len, bytes seq] × num_queries

std::vector<uint8_t> serialize(const SearchRequest& req) {
    std::vector<uint8_t> buf;
    BinaryWriter w(buf);
    buf.reserve(256);

    w.u8(req.k);
    w.u16(req.stage2_min_score);
    w.u16(req.stage2_max_gap);
    w.u8(req.stage2_min_nhit_diag);
    w.u16(req.stage1_topn);
    w.u16(req.stage1_min_score);
    w.u16(req.nresult);
    w.u16(req.stage1_min_score_frac_x10000);
    w.u8(static_cast<uint8_t>(req.seqidlist_mode));
    w.u8(req.mode);
    w.u8(req.stage1_score);
    w.u8(req.accept_qdegen);
    w.i8(req.strand);
    w.u8(req.has_stage2_min_score);
    w.u16(req.stage2_max_lookback);
    w.u8(req.stage3_traceback);
    w.i16(req.stage3_gapopen);
    w.i16(req.stage3_gapext);
    w.u16(req.stage3_min_ppositive_x100);
    w.u32(req.stage3_min_npositive);
    w.u32(req.context_abs);
    w.u16(req.context_frac_x10000);
    w.u16(req.max_degen_expand);
    w.u16(req.stage2_max_nhit_per_subject);
    w.u8(req.t);
    w.u8(req.template_type);
    w.u8(req.score_matrix);
    w.str16(req.db);

    w.u32(static_cast<uint32_t>(req.seqids.size()));
    for (const auto& acc : req.seqids) {
        w.str16(acc);
    }

    w.u16(static_cast<uint16_t>(req.queries.size()));
    for (const auto& q : req.queries) {
        w.str16(q.qseqid);
        w.u32(static_cast<uint32_t>(q.sequence.size()));
        buf.insert(buf.end(), q.sequence.begin(), q.sequence.end());
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, SearchRequest& req) {
    BinaryReader r(data.data(), data.size());

    if (!r.get_u8(req.k)) return false;
    if (!r.get_u16(req.stage2_min_score)) return false;
    if (!r.get_u16(req.stage2_max_gap)) return false;
    if (!r.get_u8(req.stage2_min_nhit_diag)) return false;
    if (!r.get_u16(req.stage1_topn)) return false;
    if (!r.get_u16(req.stage1_min_score)) return false;
    if (!r.get_u16(req.nresult)) return false;
    if (!r.get_u16(req.stage1_min_score_frac_x10000)) return false;

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
    if (!r.get_u16(req.stage3_min_ppositive_x100)) return false;
    if (!r.get_u32(req.stage3_min_npositive)) return false;
    if (!r.get_u32(req.context_abs)) return false;
    if (!r.get_u16(req.context_frac_x10000)) return false;
    if (!r.get_u16(req.max_degen_expand)) return false;
    if (!r.get_u16(req.stage2_max_nhit_per_subject)) return false;
    if (!r.get_u8(req.t)) return false;
    if (!r.get_u8(req.template_type)) return false;
    if (!r.get_u8(req.score_matrix)) return false;
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
        if (!r.get_str16(req.queries[i].qseqid)) return false;
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
//   u8   t
//   str16 db
//   u16  num_queries
//   for each query:
//     str16 qseqid
//     u8    skip_reason (0 = ok; see SkipReason enum)
//     str16 skip_detail
//     u32   qlen
//     u8    warnings
//     u16   num_hits
//     for each hit:
//       str16  sseqid
//       u8     sstrand
//       u32    qstart
//       u32    qend
//       u32    qlen
//       u32    sstart
//       u32    send
//       u32    slen
//       u16    coverscore
//       u16    matchscore
//       u16    chainscore
//       u16    volume
//       i32    alnscore
//       u32    npositive
//       u32    nnegative
//       u16    ppositive_x100
//       str16  cigar
//       str16  qseq
//       str16  sseq
//   u16  num_rejected
//     [str16 qseqid] × num_rejected

std::vector<uint8_t> serialize(const SearchResponse& resp) {
    std::vector<uint8_t> buf;
    BinaryWriter w(buf);
    buf.reserve(1024);

    w.u8(resp.status);
    w.u8(resp.k);
    w.u8(resp.mode);
    w.u8(resp.stage1_score);
    w.u8(resp.stage3_traceback);
    w.u8(resp.t);
    w.str16(resp.db);
    w.u16(static_cast<uint16_t>(resp.results.size()));

    for (const auto& qr : resp.results) {
        w.str16(qr.qseqid);
        w.u8(qr.skip_reason);
        w.str16(qr.skip_detail);
        w.u32(qr.qlen);
        w.u8(qr.warnings);
        w.u16(static_cast<uint16_t>(qr.hits.size()));
        for (const auto& hit : qr.hits) {
            w.str16(hit.sseqid);
            w.u8(hit.sstrand);
            w.u32(hit.qstart);
            w.u32(hit.qend);
            w.u32(hit.qlen);
            w.u32(hit.sstart);
            w.u32(hit.send);
            w.u32(hit.slen);
            w.u16(hit.coverscore);
            w.u16(hit.matchscore);
            w.u16(hit.chainscore);
            w.u16(hit.volume);
            w.i32(hit.alnscore);
            w.u32(hit.npositive);
            w.u32(hit.nnegative);
            w.u16(hit.ppositive_x100);
            w.str16(hit.cigar);
            w.str16(hit.qseq);
            w.str16(hit.sseq);
        }
    }

    // Rejected query IDs
    w.u16(static_cast<uint16_t>(resp.rejected_qseqids.size()));
    for (const auto& qid : resp.rejected_qseqids) {
        w.str16(qid);
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, SearchResponse& resp) {
    BinaryReader r(data.data(), data.size());

    if (!r.get_u8(resp.status)) return false;
    if (!r.get_u8(resp.k)) return false;
    if (!r.get_u8(resp.mode)) return false;
    if (!r.get_u8(resp.stage1_score)) return false;
    if (!r.get_u8(resp.stage3_traceback)) return false;
    if (!r.get_u8(resp.t)) return false;
    if (!r.get_str16(resp.db)) return false;

    uint16_t num_queries;
    if (!r.get_u16(num_queries)) return false;
    resp.results.resize(num_queries);

    for (uint16_t qi = 0; qi < num_queries; qi++) {
        auto& qr = resp.results[qi];
        if (!r.get_str16(qr.qseqid)) return false;
        if (!r.get_u8(qr.skip_reason)) return false;
        if (!r.get_str16(qr.skip_detail)) return false;
        if (!r.get_u32(qr.qlen)) return false;
        if (!r.get_u8(qr.warnings)) return false;

        uint16_t num_hits;
        if (!r.get_u16(num_hits)) return false;
        qr.hits.resize(num_hits);

        for (uint16_t hi = 0; hi < num_hits; hi++) {
            auto& hit = qr.hits[hi];
            if (!r.get_str16(hit.sseqid)) return false;
            if (!r.get_u8(hit.sstrand)) return false;
            if (!r.get_u32(hit.qstart)) return false;
            if (!r.get_u32(hit.qend)) return false;
            if (!r.get_u32(hit.qlen)) return false;
            if (!r.get_u32(hit.sstart)) return false;
            if (!r.get_u32(hit.send)) return false;
            if (!r.get_u32(hit.slen)) return false;
            if (!r.get_u16(hit.coverscore)) return false;
            if (!r.get_u16(hit.matchscore)) return false;
            if (!r.get_u16(hit.chainscore)) return false;
            if (!r.get_u16(hit.volume)) return false;
            if (!r.get_i32(hit.alnscore)) return false;
            if (!r.get_u32(hit.npositive)) return false;
            if (!r.get_u32(hit.nnegative)) return false;
            if (!r.get_u16(hit.ppositive_x100)) return false;
            if (!r.get_str16(hit.cigar)) return false;
            if (!r.get_str16(hit.qseq)) return false;
            if (!r.get_str16(hit.sseq)) return false;
        }
    }

    // Rejected query IDs
    uint16_t num_rejected;
    if (!r.get_u16(num_rejected)) return false;
    resp.rejected_qseqids.resize(num_rejected);
    for (uint16_t i = 0; i < num_rejected; i++) {
        if (!r.get_str16(resp.rejected_qseqids[i])) return false;
    }

    return true;
}

// --- ErrorResponse ---

std::vector<uint8_t> serialize(const ErrorResponse& err) {
    std::vector<uint8_t> buf;
    BinaryWriter w(buf);
    w.u32(err.error_code);
    w.str16(err.message);
    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, ErrorResponse& err) {
    BinaryReader r(data.data(), data.size());
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
    BinaryWriter w(buf);
    w.u8(resp.status);
    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, HealthResponse& resp) {
    BinaryReader r(data.data(), data.size());
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
// Wire format (v10 — adds the per-database fragment-indexing triplet
// {min_seq_length, min_length_split, overlap_length} after the existing
// per-database fields).  Older clients deserialising a v10 frame will
// stop at the wrong byte boundary; callers must rebuild client/server
// together when the format_version triplet at the top of `core/config.hpp`
// changes.
//   u8   status
//   u8   default_k
//   i32  max_queue_size
//   i32  queue_depth
//   i32  max_nseq_per_req
//   u16  num_databases
//   for each database:
//     str16 name
//     u8    default_k
//     u8    max_mode
//     u32   min_seq_length     (v10)
//     u32   min_length_split   (v10)
//     u32   overlap_length     (v10)
//     u16   num_groups
//     for each group:
//       u8  k
//       u8  kmer_type
//       u8  t
//       u8  template_type
//       u16 num_volumes
//       for each volume:
//         u16 volume_index
//         u32 num_sequences
//         u64 total_postings
//         u64 total_bases
//         str16 db

std::vector<uint8_t> serialize(const InfoResponse& resp) {
    std::vector<uint8_t> buf;
    BinaryWriter w(buf);
    buf.reserve(256);

    w.u8(resp.status);
    w.u8(resp.default_k);
    w.i32(resp.max_queue_size);
    w.i32(resp.queue_depth);
    w.i32(resp.max_nseq_per_req);
    w.u16(static_cast<uint16_t>(resp.databases.size()));

    for (const auto& db : resp.databases) {
        w.str16(db.name);
        w.u8(db.default_k);
        w.u8(db.max_mode);
        w.u32(db.min_seq_length);
        w.u32(db.min_length_split);
        w.u32(db.overlap_length);
        w.u16(static_cast<uint16_t>(db.groups.size()));

        for (const auto& g : db.groups) {
            w.u8(g.k);
            w.u8(g.kmer_type);
            w.u8(g.t);
            w.u8(g.template_type);
            w.u16(static_cast<uint16_t>(g.volumes.size()));

            for (const auto& v : g.volumes) {
                w.u16(v.volume_index);
                w.u32(v.num_sequences);
                w.u64(v.total_postings);
                w.u64(v.total_bases);
                w.str16(v.db);
            }
        }
    }

    return buf;
}

bool deserialize(const std::vector<uint8_t>& data, InfoResponse& resp) {
    BinaryReader r(data.data(), data.size());

    if (!r.get_u8(resp.status)) return false;
    if (!r.get_u8(resp.default_k)) return false;
    if (!r.get_i32(resp.max_queue_size)) return false;
    if (!r.get_i32(resp.queue_depth)) return false;
    if (!r.get_i32(resp.max_nseq_per_req)) return false;

    uint16_t num_databases;
    if (!r.get_u16(num_databases)) return false;
    resp.databases.resize(num_databases);

    for (uint16_t di = 0; di < num_databases; di++) {
        auto& db = resp.databases[di];
        if (!r.get_str16(db.name)) return false;
        if (!r.get_u8(db.default_k)) return false;
        if (!r.get_u8(db.max_mode)) return false;
        if (!r.get_u32(db.min_seq_length)) return false;
        if (!r.get_u32(db.min_length_split)) return false;
        if (!r.get_u32(db.overlap_length)) return false;

        uint16_t num_groups;
        if (!r.get_u16(num_groups)) return false;
        db.groups.resize(num_groups);

        for (uint16_t gi = 0; gi < num_groups; gi++) {
            auto& g = db.groups[gi];
            if (!r.get_u8(g.k)) return false;
            if (!r.get_u8(g.kmer_type)) return false;
            if (!r.get_u8(g.t)) return false;
            if (!r.get_u8(g.template_type)) return false;

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
