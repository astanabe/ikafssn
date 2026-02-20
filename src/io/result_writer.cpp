#include "io/result_writer.hpp"

#include <map>

namespace ikafssn {

void write_results_tab(std::ostream& out,
                       const std::vector<OutputHit>& hits) {
    out << "# query_id\taccession\tstrand\tq_start\tq_end\ts_start\ts_end\tscore\tvolume\n";
    for (const auto& h : hits) {
        out << h.query_id << '\t'
            << h.accession << '\t'
            << h.strand << '\t'
            << h.q_start << '\t'
            << h.q_end << '\t'
            << h.s_start << '\t'
            << h.s_end << '\t'
            << h.score << '\t'
            << h.volume << '\n';
    }
}

static void json_escape(std::ostream& out, const std::string& s) {
    out << '"';
    for (char c : s) {
        switch (c) {
            case '"':  out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n";  break;
            case '\r': out << "\\r";  break;
            case '\t': out << "\\t";  break;
            default:   out << c;      break;
        }
    }
    out << '"';
}

void write_results_json(std::ostream& out,
                        const std::vector<OutputHit>& hits) {
    // Group hits by query_id (preserve order of first appearance)
    std::vector<std::string> query_order;
    std::map<std::string, std::vector<const OutputHit*>> by_query;
    for (const auto& h : hits) {
        if (by_query.find(h.query_id) == by_query.end()) {
            query_order.push_back(h.query_id);
        }
        by_query[h.query_id].push_back(&h);
    }

    out << "{\n  \"results\": [\n";
    for (size_t qi = 0; qi < query_order.size(); qi++) {
        const auto& qid = query_order[qi];
        const auto& qhits = by_query[qid];
        out << "    {\n      \"query_id\": ";
        json_escape(out, qid);
        out << ",\n      \"hits\": [\n";
        for (size_t hi = 0; hi < qhits.size(); hi++) {
            const auto* h = qhits[hi];
            out << "        {\n";
            out << "          \"accession\": "; json_escape(out, h->accession); out << ",\n";
            out << "          \"strand\": \"" << h->strand << "\",\n";
            out << "          \"q_start\": " << h->q_start << ",\n";
            out << "          \"q_end\": " << h->q_end << ",\n";
            out << "          \"s_start\": " << h->s_start << ",\n";
            out << "          \"s_end\": " << h->s_end << ",\n";
            out << "          \"score\": " << h->score << ",\n";
            out << "          \"volume\": " << h->volume << "\n";
            out << "        }";
            if (hi + 1 < qhits.size()) out << ',';
            out << '\n';
        }
        out << "      ]\n    }";
        if (qi + 1 < query_order.size()) out << ',';
        out << '\n';
    }
    out << "  ]\n}\n";
}

void write_results(std::ostream& out,
                   const std::vector<OutputHit>& hits,
                   OutputFormat fmt) {
    switch (fmt) {
        case OutputFormat::kTab:
            write_results_tab(out, hits);
            break;
        case OutputFormat::kJson:
            write_results_json(out, hits);
            break;
    }
}

} // namespace ikafssn
