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
    const std::map<int, KmerGroup>& kmer_groups,
    int default_k,
    const SearchConfig& default_config,
    Server& server,
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

        logger.debug("Search request: k=%d, %zu queries, %zu seqids",
                      req.k, req.queries.size(), req.seqids.size());

        SearchResponse resp = process_search_request(
            req, kmer_groups, default_k, default_config, server, arena);

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
        iresp.default_k = static_cast<uint8_t>(default_k);

        for (const auto& [k, group] : kmer_groups) {
            KmerGroupInfo gi;
            gi.k = static_cast<uint8_t>(k);
            gi.kmer_type = group.kmer_type;

            for (const auto& vol : group.volumes) {
                VolumeInfo vi;
                vi.volume_index = vol.volume_index;
                vi.num_sequences = vol.kix.num_sequences();
                vi.total_postings = vol.kix.total_postings();
                // Extract db_name from kix header
                const auto& hdr = vol.kix.header();
                vi.db_name = std::string(hdr.db_name,
                    strnlen(hdr.db_name, sizeof(hdr.db_name)));
                gi.volumes.push_back(std::move(vi));
            }

            iresp.groups.push_back(std::move(gi));
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
