#pragma once

#include <cstdint>
#include <vector>

#include "protocol/messages.hpp"

namespace ikafssn {

// Serialize messages to byte vectors (for frame payload).

std::vector<uint8_t> serialize(const SearchRequest& req);
std::vector<uint8_t> serialize(const SearchResponse& resp);
std::vector<uint8_t> serialize(const ErrorResponse& err);
std::vector<uint8_t> serialize(const HealthRequest& req);
std::vector<uint8_t> serialize(const HealthResponse& resp);

// Deserialize messages from byte vectors.
// Returns false on malformed data.

bool deserialize(const std::vector<uint8_t>& data, SearchRequest& req);
bool deserialize(const std::vector<uint8_t>& data, SearchResponse& resp);
bool deserialize(const std::vector<uint8_t>& data, ErrorResponse& err);
bool deserialize(const std::vector<uint8_t>& data, HealthRequest& req);
bool deserialize(const std::vector<uint8_t>& data, HealthResponse& resp);

} // namespace ikafssn
