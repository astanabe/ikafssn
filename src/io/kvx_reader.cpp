#include "io/kvx_reader.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>

namespace ikafssn {

std::optional<KvxData> read_kvx(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) return std::nullopt;

    KvxData data;
    std::string line;

    while (std::getline(ifs, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        if (line.compare(0, 6, "TITLE ") == 0) {
            data.title = line.substr(6);
        } else if (line.compare(0, 6, "DBLIST") == 0) {
            // Parse quoted volume basenames from DBLIST line
            // Format: DBLIST "vol0" "vol1" ...
            size_t pos = 6; // skip "DBLIST"
            while (pos < line.size()) {
                // Find next opening quote
                size_t q1 = line.find('"', pos);
                if (q1 == std::string::npos) break;
                size_t q2 = line.find('"', q1 + 1);
                if (q2 == std::string::npos) break;
                data.volume_basenames.push_back(line.substr(q1 + 1, q2 - q1 - 1));
                pos = q2 + 1;
            }
        }
    }

    if (data.volume_basenames.empty()) return std::nullopt;

    return data;
}

} // namespace ikafssn
