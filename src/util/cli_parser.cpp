#include "util/cli_parser.hpp"

#include <cstdlib>

namespace ikafssn {

namespace {
// Normalize an argument key by collapsing a leading "--" to a single "-".
// This makes "-foo" and "--foo" interchangeable for callers.
inline std::string normalize_key(const std::string& key) {
    if (key.size() >= 2 && key[0] == '-' && key[1] == '-') {
        return key.substr(1);
    }
    return key;
}
} // namespace

CliParser::CliParser(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg.size() >= 2 && arg[0] == '-') {
            // Handle --key=value (also -key=value) syntax
            auto eq = arg.find('=');
            if (eq != std::string::npos) {
                std::string key = normalize_key(arg.substr(0, eq));
                opts_[key].push_back(arg.substr(eq + 1));
                continue;
            }

            std::string key = normalize_key(arg);

            // Flag (no value) or key-value pair
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                opts_[key].push_back(argv[i + 1]);
                i++;
            } else {
                opts_[key].push_back("1");
            }
        }
    }
}

bool CliParser::has(const std::string& key) const {
    return opts_.count(normalize_key(key)) > 0;
}

std::string CliParser::get_string(const std::string& key,
                                   const std::string& default_val) const {
    auto it = opts_.find(normalize_key(key));
    if (it != opts_.end() && !it->second.empty()) return it->second.back();
    return default_val;
}

std::vector<std::string> CliParser::get_strings(const std::string& key) const {
    auto it = opts_.find(normalize_key(key));
    if (it != opts_.end()) return it->second;
    return {};
}

int CliParser::get_int(const std::string& key, int default_val) const {
    auto it = opts_.find(normalize_key(key));
    if (it == opts_.end() || it->second.empty()) return default_val;
    try {
        return std::stoi(it->second.back());
    } catch (...) {
        return default_val;
    }
}

double CliParser::get_double(const std::string& key, double default_val) const {
    auto it = opts_.find(normalize_key(key));
    if (it == opts_.end() || it->second.empty()) return default_val;
    try {
        return std::stod(it->second.back());
    } catch (...) {
        return default_val;
    }
}

} // namespace ikafssn
