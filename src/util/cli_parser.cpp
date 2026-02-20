#include "util/cli_parser.hpp"

#include <cstdlib>

namespace ikafssn {

CliParser::CliParser(int argc, char* argv[]) {
    if (argc > 0) {
        program_ = argv[0];
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg.size() >= 2 && arg[0] == '-') {
            // Handle --verbose or -v style flags
            std::string key = arg;

            // Check if this is a flag (no value) or key-value pair
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                opts_[key] = argv[i + 1];
                i++;
            } else {
                opts_[key] = "1";
            }
        } else {
            positional_.push_back(arg);
        }
    }
}

bool CliParser::has(const std::string& key) const {
    return opts_.count(key) > 0;
}

std::string CliParser::get_string(const std::string& key,
                                   const std::string& default_val) const {
    auto it = opts_.find(key);
    if (it != opts_.end()) return it->second;
    return default_val;
}

int CliParser::get_int(const std::string& key, int default_val) const {
    auto it = opts_.find(key);
    if (it == opts_.end()) return default_val;
    try {
        return std::stoi(it->second);
    } catch (...) {
        return default_val;
    }
}

} // namespace ikafssn
