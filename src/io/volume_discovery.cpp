#include "io/volume_discovery.hpp"
#include "io/kvx_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <set>
#include <string_view>

namespace ikafssn {

IndexPrefixParts parse_index_prefix(const std::string& ix_prefix) {
    std::filesystem::path p(ix_prefix);
    IndexPrefixParts parts;
    parts.parent_dir = p.parent_path().string();
    parts.db = p.filename().string();
    if (parts.parent_dir.empty()) parts.parent_dir = ".";
    return parts;
}

namespace {

// Append a `.<key><value>` suffix to `out` if `include` is true.
void append_uint_suffix(std::string& out, const char* key, uint64_t value,
                        bool include) {
    if (!include) return;
    char buf[32];
    std::snprintf(buf, sizeof(buf), ".%s%llu", key,
                  static_cast<unsigned long long>(value));
    out += buf;
}

std::string build_stem_common(const std::string& parent_dir,
                              const std::string& base, int k,
                              uint8_t t,
                              uint8_t template_type,
                              uint32_t min_seq_length,
                              uint32_t min_length_split,
                              uint32_t overlap_length,
                              uint64_t max_freq_build,
                              uint32_t max_degen_expand) {
    std::string stem = parent_dir + "/" + base;
    char kbuf[16];
    std::snprintf(kbuf, sizeof(kbuf), ".k%d", k);
    stem += kbuf;

    append_uint_suffix(stem, "t", static_cast<uint64_t>(t), t > 0);
    append_uint_suffix(stem, "minlen", min_seq_length, min_seq_length != 0);
    append_uint_suffix(stem, "minsplit", min_length_split,
                       min_length_split != 0);
    append_uint_suffix(stem, "ovllen", overlap_length, overlap_length != 0);
    append_uint_suffix(stem, "maxfreq", max_freq_build, max_freq_build != 1);
    append_uint_suffix(stem, "maxexpand", static_cast<uint64_t>(max_degen_expand),
                       max_degen_expand > 1);

    if (t > 0) {
        const char* type_str =
            template_type == 1 ? "cod" :
            template_type == 2 ? "opt" : "con";
        stem += std::string(".") + type_str;
    }
    return stem;
}

// Match a `.<key><digits>` token against `tail`.  On match, write the
// parsed integer to `out` and return true; otherwise return false and
// leave `out` unchanged.
bool match_uint_suffix(std::string_view tail, const char* key,
                       uint64_t& out) {
    const std::size_t klen = std::strlen(key);
    if (tail.size() < klen + 2) return false;
    if (tail[0] != '.') return false;
    if (tail.compare(1, klen, key) != 0) return false;
    if (!std::isdigit(static_cast<unsigned char>(tail[1 + klen]))) return false;
    uint64_t v = 0;
    for (std::size_t i = 1 + klen; i < tail.size(); i++) {
        char c = tail[i];
        if (!std::isdigit(static_cast<unsigned char>(c))) return false;
        v = v * 10 + static_cast<uint64_t>(c - '0');
    }
    out = v;
    return true;
}

// Match the optional trailing template-type tag (.cod / .opt) at the end
// of `s` and remove it from `s`. Returns the parsed value (0 if no tag).
uint8_t consume_template_tag(std::string_view& s) {
    if (s.size() >= 4 && s.substr(s.size() - 4) == ".cod") {
        s.remove_suffix(4);
        return 1;
    }
    if (s.size() >= 4 && s.substr(s.size() - 4) == ".opt") {
        s.remove_suffix(4);
        return 2;
    }
    if (s.size() >= 4 && s.substr(s.size() - 4) == ".con") {
        s.remove_suffix(4);
        return 0;
    }
    return 0;
}

// Strip the trailing `.kix|.kpx|.ksx|.khx|.kvx` extension if present.
std::string_view strip_known_extension(std::string_view s) {
    static constexpr const char* exts[] = {
        ".kix", ".kpx", ".ksx", ".khx", ".kvx"
    };
    for (const char* e : exts) {
        std::size_t n = std::strlen(e);
        if (s.size() >= n && s.substr(s.size() - n) == e) {
            return s.substr(0, s.size() - n);
        }
    }
    return s;
}

} // namespace

std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k,
                            uint8_t t,
                            uint8_t template_type,
                            uint32_t min_seq_length,
                            uint32_t min_length_split,
                            uint32_t overlap_length,
                            uint64_t max_freq_build,
                            uint32_t max_degen_expand) {
    return build_stem_common(parent_dir, vol_basename, k, t, template_type,
                             min_seq_length, min_length_split, overlap_length,
                             max_freq_build, max_degen_expand);
}

std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k,
                          uint8_t t,
                          uint8_t template_type,
                          uint32_t min_seq_length,
                          uint32_t min_length_split,
                          uint32_t overlap_length,
                          uint64_t max_freq_build,
                          uint32_t max_degen_expand) {
    return build_stem_common(parent_dir, db, k, t, template_type,
                             min_seq_length, min_length_split, overlap_length,
                             max_freq_build, max_degen_expand) + ".khx";
}

bool parse_index_filename(const std::string& path, IndexFilenameParts& out) {
    std::filesystem::path p(path);
    out = IndexFilenameParts{};
    out.parent_dir = p.parent_path().string();
    if (out.parent_dir.empty()) out.parent_dir = ".";

    std::string fname = p.filename().string();
    std::string_view sv = strip_known_extension(fname);

    out.template_type = consume_template_tag(sv);

    // Walk encoded suffixes right-to-left.  Suffix order is accepted
    // freely (the stem builder emits a canonical order, but the parser
    // does not enforce it).  k<X> is the leftmost suffix and terminates
    // the walk.
    while (true) {
        std::size_t pos = sv.find_last_of('.');
        if (pos == std::string_view::npos) break;
        std::string_view tail = sv.substr(pos);

        uint64_t v = 0;
        if (match_uint_suffix(tail, "k", v)) {
            out.k = static_cast<int>(v);
            sv.remove_suffix(tail.size());
            break;
        } else if (match_uint_suffix(tail, "t", v)) {
            out.t = static_cast<uint8_t>(v);
            sv.remove_suffix(tail.size());
        } else if (match_uint_suffix(tail, "minlen", v)) {
            out.min_seq_length = static_cast<uint32_t>(v);
            sv.remove_suffix(tail.size());
        } else if (match_uint_suffix(tail, "minsplit", v)) {
            out.min_length_split = static_cast<uint32_t>(v);
            sv.remove_suffix(tail.size());
        } else if (match_uint_suffix(tail, "ovllen", v)) {
            out.overlap_length = static_cast<uint32_t>(v);
            sv.remove_suffix(tail.size());
        } else if (match_uint_suffix(tail, "maxfreq", v)) {
            out.max_freq_build = v;
            sv.remove_suffix(tail.size());
        } else if (match_uint_suffix(tail, "maxexpand", v)) {
            out.max_degen_expand = static_cast<uint32_t>(v);
            sv.remove_suffix(tail.size());
        } else {
            break;
        }
    }

    if (out.k == 0) return false;

    // What remains in `sv` is the base name: "<DB>.<vol>" for per-volume
    // files (numeric trailing component), "<DB>" otherwise.
    std::string base(sv);
    auto dot = base.rfind('.');
    bool has_vol = false;
    if (dot != std::string::npos) {
        std::string tail = base.substr(dot + 1);
        has_vol = !tail.empty() &&
            std::all_of(tail.begin(), tail.end(), [](unsigned char c) {
                return std::isdigit(c);
            });
    }
    out.has_vol = has_vol;
    out.vol_basename = base;
    out.db_basename = has_vol ? base.substr(0, dot) : base;
    return true;
}

namespace {

// Discover all per-volume index files registered in one .kvx manifest.
bool discover_from_kvx(const std::string& parent_dir,
                        const std::string& db,
                        const IndexFilenameParts& parts,
                        std::vector<DiscoveredVolume>& volumes) {
    std::string kvx_stem = build_stem_common(parent_dir, db,
        parts.k, parts.t, parts.template_type,
        parts.min_seq_length, parts.min_length_split, parts.overlap_length,
        parts.max_freq_build, parts.max_degen_expand);
    std::string kvx_path = kvx_stem + ".kvx";

    auto kvx = read_kvx(kvx_path);
    if (!kvx) return false;

    bool found_any = false;
    for (uint16_t vi = 0; vi < kvx->volume_basenames.size(); vi++) {
        std::string base = build_stem_common(parent_dir,
            kvx->volume_basenames[vi],
            parts.k, parts.t, parts.template_type,
            parts.min_seq_length, parts.min_length_split, parts.overlap_length,
            parts.max_freq_build, parts.max_degen_expand);
        if (!std::filesystem::exists(base + ".kix")) continue;

        DiscoveredVolume dv;
        dv.volume_index = vi;
        dv.k = parts.k;
        dv.t = parts.t;
        dv.template_type = parts.template_type;
        dv.min_seq_length = parts.min_seq_length;
        dv.min_length_split = parts.min_length_split;
        dv.overlap_length = parts.overlap_length;
        dv.max_freq_build = parts.max_freq_build;
        dv.max_degen_expand = parts.max_degen_expand;
        dv.kix_path = base + ".kix";
        dv.kpx_path = base + ".kpx";
        dv.ksx_path = base + ".ksx";
        dv.has_kpx = std::filesystem::exists(dv.kpx_path);
        volumes.push_back(std::move(dv));
        found_any = true;
    }
    return found_any;
}

// Scan the directory for .kvx manifests matching this DB basename.
std::vector<IndexFilenameParts>
scan_kvx_manifests(const std::string& parent_dir, const std::string& db) {
    std::vector<IndexFilenameParts> out;
    std::string prefix_dot = db + ".";
    std::error_code ec;
    if (!std::filesystem::exists(parent_dir, ec)) return out;
    for (const auto& entry :
         std::filesystem::directory_iterator(parent_dir, ec)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        if (fname.size() <= prefix_dot.size() ||
            fname.compare(0, prefix_dot.size(), prefix_dot) != 0)
            continue;
        if (fname.size() < 4 || fname.substr(fname.size() - 4) != ".kvx")
            continue;
        IndexFilenameParts parts;
        if (!parse_index_filename(entry.path().string(), parts)) continue;
        if (parts.has_vol) continue;       // .kvx is per-DB, never per-volume
        if (parts.db_basename != db) continue;
        out.push_back(parts);
    }
    return out;
}

} // namespace

bool variant_matches(const VariantFilter& f, const VariantFields& v) {
    if (f.has_k && v.k != f.k) return false;
    if (!f.t_is_wildcard && v.t != f.t) return false;
    if (f.has_template_type) {
        if (f.template_type == 3) {
            // "both" matches whichever spaced types exist (coding / optimal).
            if (v.template_type != 1 && v.template_type != 2) return false;
        } else if (v.template_type != f.template_type) {
            return false;
        }
    } else if (!f.t_is_wildcard && f.t > 0) {
        // A concrete non-zero -t with no -template_type selects the spaced
        // types only (coding + optimal), never the contiguous index.
        if (v.template_type == 0) return false;
    }
    if (f.has_min_seq_length   && v.min_seq_length   != f.min_seq_length)   return false;
    if (f.has_min_length_split && v.min_length_split != f.min_length_split) return false;
    if (f.has_overlap_length   && v.overlap_length   != f.overlap_length)   return false;
    if (f.has_max_freq_build   && v.max_freq_build   != f.max_freq_build)   return false;
    if (f.has_max_degen_expand && v.max_degen_expand != f.max_degen_expand) return false;
    return true;
}

bool same_variant7(const VariantFields& a, const VariantFields& b) {
    return a.k == b.k && a.t == b.t && a.template_type == b.template_type &&
           a.min_seq_length == b.min_seq_length &&
           a.min_length_split == b.min_length_split &&
           a.overlap_length == b.overlap_length &&
           a.max_freq_build == b.max_freq_build;
}

bool same_variant6(const VariantFields& a, const VariantFields& b) {
    return a.k == b.k && a.t == b.t &&
           a.min_seq_length == b.min_seq_length &&
           a.min_length_split == b.min_length_split &&
           a.overlap_length == b.overlap_length &&
           a.max_freq_build == b.max_freq_build;
}

uint32_t choose_max_degen_expand(const std::vector<uint32_t>& present,
                                 uint32_t query_value) {
    uint32_t best = 0;
    bool has_query = false;
    for (uint32_t v : present) {
        if (v == query_value) has_query = true;
        if (v > best) best = v;
    }
    return has_query ? query_value : best;
}

namespace {

constexpr const char* kAmbiguousVariantErr =
    "index filter is ambiguous: multiple index variants match; narrow the "
    "selection with -k / -t / -template_type / -min_seq_length / "
    "-min_length_split / -overlap_length / -max_freq_build";

// Select the indices of a single variant from a one-7-tuple candidate set.
std::vector<size_t> select_one_side(const std::vector<VariantFields>& cands,
                                    const std::vector<size_t>& idx,
                                    uint32_t query_max_degen_expand,
                                    std::string& err) {
    for (size_t j = 1; j < idx.size(); ++j) {
        if (!same_variant7(cands[idx[j]], cands[idx[0]])) {
            err = kAmbiguousVariantErr;
            return {};
        }
    }
    std::vector<uint32_t> present;
    present.reserve(idx.size());
    for (size_t i : idx) present.push_back(cands[i].max_degen_expand);
    uint32_t chosen = choose_max_degen_expand(present, query_max_degen_expand);
    std::vector<size_t> out;
    for (size_t i : idx)
        if (cands[i].max_degen_expand == chosen) out.push_back(i);
    return out;
}

// Resolve a coding+optimal pair with NO max_degen_expand mixing: both sides
// must agree on the 6-field pairing key and share one max_degen_expand value.
std::vector<size_t> resolve_both_pair(const std::vector<VariantFields>& cands,
                                      const std::vector<size_t>& cod,
                                      const std::vector<size_t>& opt,
                                      uint32_t query_max_degen_expand,
                                      std::string& err) {
    // Every coding and optimal candidate must share the 6-field pairing key
    // (everything but template_type and max_degen_expand).
    const VariantFields& ref = cands[cod[0]];
    for (size_t i : cod)
        if (!same_variant6(cands[i], ref)) { err = kAmbiguousVariantErr; return {}; }
    for (size_t i : opt)
        if (!same_variant6(cands[i], ref)) { err = kAmbiguousVariantErr; return {}; }

    // The shared max_degen_expand must exist on BOTH sides (no chanpon).
    std::set<uint32_t> cod_e, opt_e;
    for (size_t i : cod) cod_e.insert(cands[i].max_degen_expand);
    for (size_t i : opt) opt_e.insert(cands[i].max_degen_expand);
    std::vector<uint32_t> common;
    for (uint32_t e : cod_e)
        if (opt_e.count(e)) common.push_back(e);
    if (common.empty()) {
        err = "coding and optimal variants share no common max_degen_expand; "
              "cannot form a both-mode pair without mixing max_degen_expand";
        return {};
    }
    uint32_t chosen = choose_max_degen_expand(common, query_max_degen_expand);

    std::vector<size_t> out;
    for (size_t i : cod)
        if (cands[i].max_degen_expand == chosen) out.push_back(i);
    for (size_t i : opt)
        if (cands[i].max_degen_expand == chosen) out.push_back(i);
    return out;
}

} // namespace

std::vector<size_t> resolve_variant_indices(
    const std::vector<VariantFields>& cands, TemplateResolveMode mode,
    uint32_t query_max_degen_expand, std::string& err) {
    if (cands.empty()) {
        err = "no index variant matches the given filter";
        return {};
    }

    if (mode == TemplateResolveMode::kSingle) {
        // The filter has already restricted the candidates to one template_type.
        std::vector<size_t> all(cands.size());
        for (size_t i = 0; i < cands.size(); ++i) all[i] = i;
        return select_one_side(cands, all, query_max_degen_expand, err);
    }

    std::vector<size_t> cod, opt;
    for (size_t i = 0; i < cands.size(); ++i) {
        if (cands[i].template_type == 1) cod.push_back(i);
        else if (cands[i].template_type == 2) opt.push_back(i);
    }
    const bool have_cod = !cod.empty();
    const bool have_opt = !opt.empty();

    if (have_cod && have_opt) {
        return resolve_both_pair(cands, cod, opt, query_max_degen_expand, err);
    }

    if (mode == TemplateResolveMode::kBothRequired) {
        err = "-template_type both requires both coding and optimal index variants";
        return {};
    }

    // kBothOrSingle with only one side present: resolve that single side.
    if (have_cod) return select_one_side(cands, cod, query_max_degen_expand, err);
    if (have_opt) return select_one_side(cands, opt, query_max_degen_expand, err);
    err = "no index variant matches the given filter";
    return {};
}

std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k) {
    auto parts = parse_index_prefix(ix_prefix);
    std::vector<DiscoveredVolume> volumes;

    auto manifests = scan_kvx_manifests(parts.parent_dir, parts.db);
    for (const auto& m : manifests) {
        if (filter_k > 0 && m.k != filter_k) continue;
        discover_from_kvx(parts.parent_dir, parts.db, m, volumes);
    }

    std::sort(volumes.begin(), volumes.end(),
              [](const DiscoveredVolume& a, const DiscoveredVolume& b) {
                  if (a.k != b.k) return a.k < b.k;
                  if (a.t != b.t) return a.t < b.t;
                  if (a.template_type != b.template_type)
                      return a.template_type < b.template_type;
                  return a.volume_index < b.volume_index;
              });
    return volumes;
}

} // namespace ikafssn
