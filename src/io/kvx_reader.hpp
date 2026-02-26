#ifndef IKAFSSN_IO_KVX_READER_HPP
#define IKAFSSN_IO_KVX_READER_HPP

#include <optional>
#include <string>
#include <vector>

namespace ikafssn {

struct KvxData {
    std::string title;
    std::vector<std::string> volume_basenames;
};

// Read a .kvx manifest file.
// Returns empty optional if file doesn't exist or parse fails.
std::optional<KvxData> read_kvx(const std::string& path);

} // namespace ikafssn

#endif // IKAFSSN_IO_KVX_READER_HPP
