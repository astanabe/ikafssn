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
            // Handle --key=value syntax for double-dash args
            if (arg.size() >= 3 && arg[0] == '-' && arg[1] == '-') {
                auto eq = arg.find('=');
                if (eq != std::string::npos) {
                    opts_[arg.substr(0, eq)].push_back(arg.substr(eq + 1));
                    continue;
                }
            }

            // Handle --verbose or -v style flags
            std::string key = arg;

            // Check if this is a flag (no value) or key-value pair
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                opts_[key].push_back(argv[i + 1]);
                i++;
            } else {
                opts_[key].push_back("1");
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
    if (it != opts_.end() && !it->second.empty()) return it->second.back();
    return default_val;
}

std::vector<std::string> CliParser::get_strings(const std::string& key) const {
    auto it = opts_.find(key);
    if (it != opts_.end()) return it->second;
    return {};
}

int CliParser::get_int(const std::string& key, int default_val) const {
    auto it = opts_.find(key);
    if (it == opts_.end() || it->second.empty()) return default_val;
    try {
        return std::stoi(it->second.back());
    } catch (...) {
        return default_val;
    }
}

double CliParser::get_double(const std::string& key, double default_val) const {
    auto it = opts_.find(key);
    if (it == opts_.end() || it->second.empty()) return default_val;
    try {
        return std::stod(it->second.back());
    } catch (...) {
        return default_val;
    }
}

} // namespace ikafssn
