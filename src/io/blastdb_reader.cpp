#include "io/blastdb_reader.hpp"

#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <objtools/blast/seqdb_reader/seqdbcommon.hpp>
#include <objects/seqloc/Seq_id.hpp>
#include <objects/general/Dbtag.hpp>
#include <objects/general/Object_id.hpp>

#include <cstdio>
#include <algorithm>

namespace ikafssn {

// NCBI NA8 encoding: single-bit flags per base.
// A=1, C=2, G=4, T=8, N=15 (all bits set), gap=0
// We map to ACGT chars; anything else becomes 'N'.
static char ncbi_na8_to_char(uint8_t val) {
    switch (val) {
        case 1:  return 'A';
        case 2:  return 'C';
        case 4:  return 'G';
        case 8:  return 'T';
        default: return 'N';
    }
}

struct BlastDbReader::Impl {
    std::unique_ptr<ncbi::CSeqDB> db;
};

BlastDbReader::BlastDbReader() : impl_(std::make_unique<Impl>()) {}

BlastDbReader::~BlastDbReader() {
    close();
}

BlastDbReader::BlastDbReader(BlastDbReader&&) noexcept = default;
BlastDbReader& BlastDbReader::operator=(BlastDbReader&&) noexcept = default;

bool BlastDbReader::open(const std::string& db_path) {
    try {
        impl_->db = std::make_unique<ncbi::CSeqDB>(
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

std::string BlastDbReader::get_sequence(uint32_t oid) const {
    if (!impl_->db) return {};

    const char* buf = nullptr;
    int len = impl_->db->GetAmbigSeq(
        static_cast<int>(oid), &buf, ncbi::kSeqDBNuclNcbiNA8);

    if (len <= 0 || !buf) {
        if (buf) impl_->db->RetAmbigSeq(&buf);
        return {};
    }

    std::string result(static_cast<size_t>(len), '\0');
    for (int i = 0; i < len; i++) {
        result[i] = ncbi_na8_to_char(static_cast<uint8_t>(buf[i]));
    }

    impl_->db->RetAmbigSeq(&buf);
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

std::vector<std::string> BlastDbReader::find_volume_paths(const std::string& db_name) {
    std::vector<std::string> paths;
    try {
        ncbi::CSeqDB::FindVolumePaths(
            db_name, ncbi::CSeqDB::eNucleotide, paths);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "BlastDbReader: FindVolumePaths('%s') failed: %s\n",
                     db_name.c_str(), e.what());
    }
    return paths;
}

} // namespace ikafssn
