#include "io/blastdb_reader.hpp"
#include "core/ambiguity_parser.hpp"

#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>
#include <objects/seqloc/Seq_id.hpp>
#include <objects/general/Dbtag.hpp>
#include <objects/general/Object_id.hpp>

#include <cstdio>
#include <algorithm>

namespace ikafssn {

// ncbi2na 2-bit code -> ASCII character
static constexpr char ncbi2na_to_char[4] = {'A', 'C', 'G', 'T'};

// ncbi4na value -> IUPAC character
static constexpr char ncbi4na_to_iupac[16] = {
    '-',  // 0: gap
    'A',  // 1
    'C',  // 2
    'M',  // 3: A,C
    'G',  // 4
    'R',  // 5: A,G
    'S',  // 6: C,G
    'V',  // 7: A,C,G
    'T',  // 8
    'W',  // 9: A,T
    'Y',  // 10: C,T
    'H',  // 11: A,C,T
    'K',  // 12: G,T
    'D',  // 13: A,G,T
    'B',  // 14: C,G,T
    'N',  // 15: A,C,G,T
};

struct BlastDbReader::Impl {
    std::unique_ptr<ncbi::CSeqDBExpert> db;
};

BlastDbReader::BlastDbReader() : impl_(std::make_unique<Impl>()) {}

BlastDbReader::~BlastDbReader() {
    close();
}

BlastDbReader::BlastDbReader(BlastDbReader&&) noexcept = default;
BlastDbReader& BlastDbReader::operator=(BlastDbReader&&) noexcept = default;

bool BlastDbReader::open(const std::string& db_path) {
    try {
        impl_->db = std::make_unique<ncbi::CSeqDBExpert>(
            db_path, ncbi::CSeqDB::eNucleotide);
        return true;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "BlastDbReader: failed to open '%s': %s\n",
                     db_path.c_str(), e.what());
        return false;
    }
}

void BlastDbReader::close() {
    impl_->db.reset();
}

bool BlastDbReader::is_open() const {
    return impl_->db != nullptr;
}

uint32_t BlastDbReader::num_sequences() const {
    if (!impl_->db) return 0;
    return static_cast<uint32_t>(impl_->db->GetNumSeqs());
}

uint32_t BlastDbReader::seq_length(uint32_t oid) const {
    if (!impl_->db) return 0;
    return static_cast<uint32_t>(impl_->db->GetSeqLength(static_cast<int>(oid)));
}

BlastDbReader::RawSequence BlastDbReader::get_raw_sequence(uint32_t oid) const {
    RawSequence raw{};
    if (!impl_->db) return raw;

    const char* buffer = nullptr;
    int seq_len = 0;
    int ambig_len = 0;
    impl_->db->GetRawSeqAndAmbig(static_cast<int>(oid),
                                  &buffer, &seq_len, &ambig_len);
    raw.ncbi2na_data = buffer;
    raw.ncbi2na_bytes = seq_len;
    raw.ambig_data = buffer ? buffer + seq_len : nullptr;
    raw.ambig_bytes = ambig_len;
    raw.seq_length = static_cast<uint32_t>(
        impl_->db->GetSeqLength(static_cast<int>(oid)));
    return raw;
}

void BlastDbReader::ret_raw_sequence(const RawSequence& raw) const {
    if (impl_->db && raw.ncbi2na_data) {
        const char* ptr = raw.ncbi2na_data;
        impl_->db->RetSequence(&ptr);
    }
}

std::string BlastDbReader::get_sequence(uint32_t oid) const {
    if (!impl_->db) return {};

    RawSequence raw = get_raw_sequence(oid);
    if (!raw.ncbi2na_data || raw.seq_length == 0) {
        if (raw.ncbi2na_data) ret_raw_sequence(raw);
        return {};
    }

    // Decode ncbi2na packed data to ASCII
    std::string result(raw.seq_length, '\0');
    for (uint32_t i = 0; i < raw.seq_length; i++) {
        uint8_t byte = static_cast<uint8_t>(raw.ncbi2na_data[i >> 2]);
        uint8_t code = (byte >> (6 - 2 * (i & 3))) & 0x03;
        result[i] = ncbi2na_to_char[code];
    }

    // Apply ambiguity data to overwrite positions with IUPAC chars
    auto ambig_entries = AmbiguityParser::parse(raw.ambig_data, raw.ambig_bytes);
    for (const auto& entry : ambig_entries) {
        char iupac = ncbi4na_to_iupac[entry.ncbi4na];
        for (uint32_t j = 0; j < entry.run_length; j++) {
            uint32_t pos = entry.position + j;
            if (pos < raw.seq_length) {
                result[pos] = iupac;
            }
        }
    }

    ret_raw_sequence(raw);
    return result;
}

std::string BlastDbReader::get_accession(uint32_t oid) const {
    if (!impl_->db) return {};

    try {
        std::list<ncbi::CRef<ncbi::objects::CSeq_id>> ids =
            impl_->db->GetSeqIDs(static_cast<int>(oid));

        // Find the best accession: prefer textseq_id with accession
        for (const auto& id : ids) {
            if (id->IsGenbank() || id->IsEmbl() || id->IsDdbj() ||
                id->IsOther() || id->IsTpg() || id->IsTpe() || id->IsTpd()) {
                const auto* tsid = id->GetTextseq_Id();
                if (tsid && tsid->IsSetAccession()) {
                    return tsid->GetAccession();
                }
            }
        }

        // Fallback: use GetSeqIdString for the first ID
        if (!ids.empty()) {
            std::string id_str;
            ids.front()->GetLabel(&id_str, ncbi::objects::CSeq_id::eContent);
            return id_str;
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "BlastDbReader: get_accession(%u) failed: %s\n",
                     oid, e.what());
    }

    return {};
}

std::string BlastDbReader::get_title() const {
    if (!impl_->db) return {};
    return impl_->db->GetTitle();
}

std::vector<std::string> BlastDbReader::find_volume_paths(const std::string& db) {
    std::vector<std::string> paths;
    try {
        ncbi::CSeqDB::FindVolumePaths(
            db, ncbi::CSeqDB::eNucleotide, paths);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "BlastDbReader: FindVolumePaths('%s') failed: %s\n",
                     db.c_str(), e.what());
    }
    return paths;
}

} // namespace ikafssn
