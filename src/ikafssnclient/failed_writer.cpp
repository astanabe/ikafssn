#include "ikafssnclient/failed_writer.hpp"
#include "protocol/messages.hpp"

namespace ikafssn {

void synth_failed_hits(const std::vector<std::string>& deflines,
                       const std::string& reason,
                       std::vector<OutputHit>& out) {
    out.reserve(out.size() + deflines.size());
    for (const auto& def : deflines) {
        OutputHit h;
        // The cached defline is just the qseqid — first whitespace-
        // delimited token of the original FASTA header, no leading '>'.
        h.qseqid = def;
        h.skip_reason = kFailHttpJob;
        h.skip_detail = reason;
        h.sstrand = '*';
        out.push_back(std::move(h));
    }
}

} // namespace ikafssn
