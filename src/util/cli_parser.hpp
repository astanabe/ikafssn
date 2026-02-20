#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace ikafssn {

// Simple command-line argument parser for -key value style arguments.
class CliParser {
public:
    CliParser(int argc, char* argv[]);

    // Check if a flag/option is present.
    bool has(const std::string& key) const;

    // Get string value for a key. Returns default_val if not found.
    std::string get_string(const std::string& key,
                           const std::string& default_val = {}) const;

    // Get integer value for a key. Returns default_val if not found or invalid.
    int get_int(const std::string& key, int default_val = 0) const;

    // Get double value for a key. Returns default_val if not found or invalid.
    double get_double(const std::string& key, double default_val = 0.0) const;

    // Get the program name (argv[0]).
    const std::string& program() const { return program_; }

    // Get positional arguments (those not preceded by a -key).
    const std::vector<std::string>& positional() const { return positional_; }

private:
    std::string program_;
    std::unordered_map<std::string, std::string> opts_;
    std::vector<std::string> positional_;
};

} // namespace ikafssn
