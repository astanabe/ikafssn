#include "ikafssnserver/connection_handler.hpp"
#include "ikafssnserver/server.hpp"
#include "util/socket_utils.hpp"

#include <cstring>

namespace ikafssn {

static void send_error(int fd, uint32_t code, const std::string& msg) {
    ErrorResponse err;
    err.error_code = code;
    err.message = msg;
    auto payload = serialize(err);
    write_frame(fd, MsgType::kErrorResponse, payload);
}

void handle_connection(
    int client_fd,
    Server& server,
    const ServerConfig& config,
    tbb::task_arena& arena,
    const Logger& logger) {

    FrameHeader hdr;
    std::vector<uint8_t> payload;

    if (!read_frame(client_fd, hdr, payload)) {
        logger.debug("Failed to read frame from client");
        close_fd(client_fd);
        return;
    }

    MsgType type = static_cast<MsgType>(hdr.msg_type);

    switch (type) {
    case MsgType::kSearchRequest: {
        SearchRequest req;
        if (!deserialize(payload, req)) {
            send_error(client_fd, 400, "Malformed search request");
            break;
        }

        // Validate db_name
        if (req.db_name.empty()) {
            send_error(client_fd, 400, "db_name is required");
            break;
        }

        const DatabaseEntry* db = server.find_database(req.db_name);
        if (!db) {
            send_error(client_fd, 404, "Database not found: " + req.db_name);
            break;
        }

        logger.debug("Search request: db=%s, k=%d, %zu queries, %zu seqids, mode=%d",
                      req.db_name.c_str(), req.k, req.queries.size(), req.seqids.size(), req.mode);

        SearchResponse resp = process_search_request(req, *db, server, arena);

        auto resp_payload = serialize(resp);
        if (!write_frame(client_fd, MsgType::kSearchResponse, resp_payload)) {
            logger.debug("Failed to send search response");
        }
        break;
    }

    case MsgType::kHealthRequest: {
        HealthResponse hresp;
        hresp.status = 0;
        auto resp_payload = serialize(hresp);
        write_frame(client_fd, MsgType::kHealthResponse, resp_payload);
        break;
    }

    case MsgType::kInfoRequest: {
        InfoResponse iresp;
        iresp.status = 0;
        iresp.default_k = static_cast<uint8_t>(server.default_k());
        iresp.max_active_sequences = server.max_active_sequences();
        iresp.active_sequences = server.active_sequences();

        for (const auto& db : server.databases()) {
            DatabaseInfo dbi;
            dbi.name = db.name;
            dbi.default_k = static_cast<uint8_t>(db.default_k);
            dbi.max_mode = db.max_mode;

            for (const auto& [k, group] : db.kmer_groups) {
                KmerGroupInfo gi;
                gi.k = static_cast<uint8_t>(k);
                gi.kmer_type = group.kmer_type;

                for (const auto& vol : group.volumes) {
                    VolumeInfo vi;
                    vi.volume_index = vol.volume_index;
                    vi.num_sequences = vol.kix.num_sequences();
                    vi.total_postings = vol.kix.total_postings();
                    vi.total_bases = vol.total_bases;
                    const auto& kix_hdr = vol.kix.header();
                    vi.db_name = std::string(kix_hdr.db_name,
                        strnlen(kix_hdr.db_name, sizeof(kix_hdr.db_name)));
                    gi.volumes.push_back(std::move(vi));
                }

                dbi.groups.push_back(std::move(gi));
            }

            iresp.databases.push_back(std::move(dbi));
        }

        auto resp_payload = serialize(iresp);
        write_frame(client_fd, MsgType::kInfoResponse, resp_payload);
        break;
    }

    default:
        send_error(client_fd, 400, "Unknown message type");
        break;
    }

    close_fd(client_fd);
}

} // namespace ikafssn
