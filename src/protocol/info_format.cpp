#include "protocol/info_format.hpp"

#include <cstdio>
#include <string>

namespace ikafssn {

std::string format_all_databases(const InfoResponse& info) {
    std::string result;
    for (const auto& db : info.databases) {
        result += "  " + db.name + "    ";
        for (size_t i = 0; i < db.groups.size(); i++) {
            if (i > 0) result += ", ";
            result += "k=" + std::to_string(db.groups[i].k)
                    + " (mode 1-" + std::to_string(db.max_mode) + ")";
        }
        result += "\n";
    }
    return result;
}

std::string validate_info(const InfoResponse& info,
                          const std::string& db_name,
                          uint8_t k, uint8_t mode,
                          bool check_slots) {
    // 1. Slot capacity check
    if (check_slots && info.max_active_sequences > 0 &&
        info.active_sequences >= info.max_active_sequences) {
        return "Error: server is at capacity ("
             + std::to_string(info.active_sequences) + "/"
             + std::to_string(info.max_active_sequences)
             + " active sequences). Try again later.";
    }

    // 2. Find target database
    const DatabaseInfo* target_db = nullptr;
    for (const auto& db : info.databases) {
        if (db.name == db_name) {
            target_db = &db;
            break;
        }
    }

    if (!target_db) {
        std::string msg = "Error: database '" + db_name
                        + "' not found on server.\nAvailable databases:\n"
                        + format_all_databases(info);
        return msg;
    }

    // 3. Check k value (k==0 means use server default, always valid)
    if (k != 0) {
        bool k_found = false;
        for (const auto& g : target_db->groups) {
            if (g.k == k) {
                k_found = true;
                break;
            }
        }
        if (!k_found) {
            std::string msg = "Error: k=" + std::to_string(k)
                            + " is not available for database '"
                            + target_db->name + "'.\nAvailable databases:\n"
                            + format_all_databases(info);
            return msg;
        }
    }

    // 4. Check mode (mode==0 means use server default, always valid)
    if (mode > 0 && mode > target_db->max_mode) {
        std::string msg = "Error: mode " + std::to_string(mode)
                        + " exceeds max mode " + std::to_string(target_db->max_mode)
                        + " for database '" + target_db->name
                        + "'.\nAvailable databases:\n"
                        + format_all_databases(info);
        return msg;
    }

    return "";
}

std::string format_server_info(const InfoResponse& info, bool verbose) {
    std::string out;
    out += "=== ikafssn Server Information ===\n\n";
    out += "Active sequences:  " + std::to_string(info.active_sequences)
         + "/" + std::to_string(info.max_active_sequences) + "\n\n";
    out += "--- Databases ---\n";

    for (const auto& db : info.databases) {
        out += "\nDatabase: " + db.name + "\n";
        out += "  Default k:       " + std::to_string(db.default_k) + "\n";
        out += "  Max mode:        " + std::to_string(db.max_mode) + "\n";
        out += "  K-mer groups:\n";

        for (const auto& g : db.groups) {
            uint64_t group_seqs = 0;
            uint64_t group_postings = 0;
            for (const auto& v : g.volumes) {
                group_seqs += v.num_sequences;
                group_postings += v.total_postings;
            }

            const char* type_str = (g.kmer_type == 0) ? "uint16" : "uint32";
            char line[256];
            std::snprintf(line, sizeof(line),
                "    k=%-3u (%s)  %zu volume(s)   %lu sequences   %lu postings\n",
                g.k, type_str, g.volumes.size(),
                static_cast<unsigned long>(group_seqs),
                static_cast<unsigned long>(group_postings));
            out += line;

            if (verbose) {
                for (const auto& v : g.volumes) {
                    char vline[256];
                    std::snprintf(vline, sizeof(vline),
                        "      Volume %u:  %u sequences  %lu postings  (%s)\n",
                        v.volume_index, v.num_sequences,
                        static_cast<unsigned long>(v.total_postings),
                        v.db_name.c_str());
                    out += vline;
                }
            }
        }
    }

    out += "\n";
    return out;
}

} // namespace ikafssn
