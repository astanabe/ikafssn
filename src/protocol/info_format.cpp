#include "protocol/info_format.hpp"

#include <cstdio>
#include <string>

#include "core/spaced_seed.hpp"

namespace ikafssn {

std::string format_all_databases(const InfoResponse& info) {
    std::string result;
    for (const auto& db : info.databases) {
        result += "  " + db.name + "    ";
        for (size_t i = 0; i < db.groups.size(); i++) {
            if (i > 0) result += ", ";
            const auto& g = db.groups[i];
            result += "k=" + std::to_string(g.k);
            if (g.t > 0) {
                result += " t=" + std::to_string(g.t)
                        + " " + template_type_to_string(static_cast<TemplateType>(g.template_type));
            }
            result += " (mode 1-" + std::to_string(db.max_mode) + ")";
        }
        result += "\n";
    }
    return result;
}

std::string validate_info(const InfoResponse& info,
                          const std::string& db,
                          uint8_t k, uint8_t mode,
                          bool check_slots,
                          uint8_t t, uint8_t template_type) {
    // 1. Slot capacity check
    if (check_slots && info.max_queue_size > 0 &&
        info.queue_depth >= info.max_queue_size) {
        return "Error: server is at capacity ("
             + std::to_string(info.queue_depth) + "/"
             + std::to_string(info.max_queue_size)
             + " active sequences). Try again later.";
    }

    // 2. Find target database
    const DatabaseInfo* target_db = nullptr;
    for (const auto& dbi : info.databases) {
        if (dbi.name == db) {
            target_db = &dbi;
            break;
        }
    }

    if (!target_db) {
        std::string msg = "Error: database '" + db
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

    // 3b. Check t/template_type (t==0 means contiguous, always valid)
    if (t != 0) {
        bool t_found = false;
        if (template_type == 3) {
            // "both" = virtual capability: coding AND optimal must both exist
            bool has_cod = false, has_opt = false;
            for (const auto& g : target_db->groups) {
                if (g.k == k && g.t == t && g.template_type == 1) has_cod = true;
                if (g.k == k && g.t == t && g.template_type == 2) has_opt = true;
            }
            t_found = has_cod && has_opt;
        } else {
            for (const auto& g : target_db->groups) {
                if (g.k == k && g.t == t &&
                    (template_type == 0 || g.template_type == template_type)) {
                    t_found = true;
                    break;
                }
            }
        }
        if (!t_found) {
            std::string msg = "Error: t=" + std::to_string(t)
                            + " is not available for database '"
                            + target_db->name + "' with k=" + std::to_string(k)
                            + ".\nAvailable databases:\n"
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
    out += "Active sequences:  " + std::to_string(info.queue_depth)
         + "/" + std::to_string(info.max_queue_size) + "\n";
    if (info.max_seqs_per_req > 0) {
        out += "Max sequences/request: " + std::to_string(info.max_seqs_per_req) + "\n";
    }
    out += "\n";
    out += "--- Databases ---\n";

    for (const auto& db : info.databases) {
        out += "\nDatabase: " + db.name + "\n";
        out += "  Default k:       " + std::to_string(db.default_k) + "\n";
        out += "  Max mode:        " + std::to_string(db.max_mode) + "\n";
        out += "  K-mer groups:\n";

        for (const auto& g : db.groups) {
            uint64_t group_seqs = 0;
            uint64_t group_bases = 0;
            uint64_t group_postings = 0;
            for (const auto& v : g.volumes) {
                group_seqs += v.num_sequences;
                group_bases += v.total_bases;
                group_postings += v.total_postings;
            }

            const char* type_str = (g.kmer_type == 0) ? "uint16" : "uint32";
            char line[512];
            if (g.t > 0) {
                std::snprintf(line, sizeof(line),
                    "    k=%-3u t=%-2u %-8s (%s)  %zu volume(s)   %lu sequences   %lu bases   %lu postings\n",
                    g.k, g.t,
                    template_type_to_string(static_cast<TemplateType>(g.template_type)).c_str(),
                    type_str, g.volumes.size(),
                    static_cast<unsigned long>(group_seqs),
                    static_cast<unsigned long>(group_bases),
                    static_cast<unsigned long>(group_postings));
            } else {
                std::snprintf(line, sizeof(line),
                    "    k=%-3u (%s)  %zu volume(s)   %lu sequences   %lu bases   %lu postings\n",
                    g.k, type_str, g.volumes.size(),
                    static_cast<unsigned long>(group_seqs),
                    static_cast<unsigned long>(group_bases),
                    static_cast<unsigned long>(group_postings));
            }
            out += line;

            if (verbose) {
                for (const auto& v : g.volumes) {
                    char vline[512];
                    std::snprintf(vline, sizeof(vline),
                        "      Volume %u:  %u sequences  %lu bases  %lu postings  (%s)\n",
                        v.volume_index, v.num_sequences,
                        static_cast<unsigned long>(v.total_bases),
                        static_cast<unsigned long>(v.total_postings),
                        v.db.c_str());
                    out += vline;
                }
            }
        }
    }

    out += "\n";
    return out;
}

} // namespace ikafssn
