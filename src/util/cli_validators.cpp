#include "util/cli_validators.hpp"

namespace ikafssn {

bool parse_max_degen_expand(const CliParser& cli, int default_val,
                            uint16_t& out, std::string& err) {
    int v = cli.get_int("-max_degen_expand", default_val);
    if (v < 0 || v > 256) {
        err = "Error: -max_degen_expand must be between 0 and 256";
        return false;
    }
    out = static_cast<uint16_t>(v);
    return true;
}

bool parse_score_matrix(const CliParser& cli, std::string& out, std::string& err) {
    if (!cli.has("-stage3_score_matrix")) return true;
    std::string sm = cli.get_string("-stage3_score_matrix");
    if (sm != "degmatch" && sm != "dnafull" && sm != "nuc44") {
        err = "Error: -stage3_score_matrix must be degmatch, dnafull, or nuc44";
        return false;
    }
    out = std::move(sm);
    return true;
}

uint8_t score_matrix_code(const std::string& name) {
    if (name == "degmatch") return 1;
    if (name == "dnafull")  return 2;
    if (name == "nuc44")    return 3;
    return 0;
}

bool parse_spaced_seed_t(const CliParser& cli, uint8_t& out, std::string& err) {
    int v = cli.get_int("-t", 0);
    if (v != 0 && v != 13 && v != 15 && v != 16 && v != 18 && v != 21) {
        err = "Error: -t must be 0, 13, 15, 16, 18, or 21";
        return false;
    }
    out = static_cast<uint8_t>(v);
    return true;
}

bool parse_template_type_cli(const CliParser& cli, TemplateType default_value,
                             bool allow_contiguous,
                             TemplateType& out, std::string& err) {
    if (!cli.has("-template_type")) {
        out = default_value;
        return true;
    }
    TemplateType t = template_type_from_string(cli.get_string("-template_type"));
    if (!allow_contiguous && t == TemplateType::kContiguous) {
        err = "Error: -template_type must be coding, optimal, or both";
        return false;
    }
    out = t;
    return true;
}

bool validate_t_template_type_combo(const CliParser& cli, std::string& err) {
    if (cli.has("-template_type") && cli.has("-t") &&
        cli.get_int("-t", 0) == 0) {
        err = "Error: -template_type cannot be combined with -t 0 "
              "(the contiguous index has no template type)";
        return false;
    }
    return true;
}

bool validate_primer_mode_options(const CliParser& cli, std::string& err) {
    if (!cli.has("-primer")) return true;
    if (cli.has("-stage1_min_score")) {
        err = "Error: -stage1_min_score cannot be used with -primer; use -stage1_primer_score instead";
        return false;
    }
    if (cli.has("-stage2_min_score")) {
        err = "Error: -stage2_min_score cannot be used with -primer; use -stage2_primer_score_add instead";
        return false;
    }
    if (!cli.has("-insert_length")) {
        err = "Error: -insert_length is required with -primer";
        return false;
    }
    return true;
}

} // namespace ikafssn
