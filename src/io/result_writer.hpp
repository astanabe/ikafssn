#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace ikafssn {

class KsxReader;

struct OutputHit {
    std::string query_id;
    std::string accession;
    char strand;          // '+' or '-'
    uint32_t q_start;
    uint32_t q_end;
    uint32_t s_start;
    uint32_t s_end;
    uint32_t score;
    uint16_t volume;
};

enum class OutputFormat { kTab, kJson };

// Write results in tab-delimited format.
void write_results_tab(std::ostream& out,
                       const std::vector<OutputHit>& hits);

// Write results in JSON format.
void write_results_json(std::ostream& out,
                        const std::vector<OutputHit>& hits);

// Write results in the specified format.
void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt);

} // namespace ikafssn
