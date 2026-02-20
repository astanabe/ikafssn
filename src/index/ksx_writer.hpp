#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class KsxWriter {
public:
    // Add a sequence's metadata. Must be called in OID order.
    void add_sequence(uint32_t seq_length, const std::string& accession);

    // Write the .ksx file. Returns true on success.
    bool write(const std::string& path) const;

    uint32_t num_sequences() const { return static_cast<uint32_t>(seq_lengths_.size()); }

private:
    std::vector<uint32_t> seq_lengths_;
    std::vector<std::string> accessions_;
};

} // namespace ikafssn
