#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

struct DiscoveredVolume {
    std::string kix_path;
    std::string kpx_path;
    std::string ksx_path;
    uint16_t volume_index;
    int k;
    bool has_kpx = true;
    uint8_t t = 0;              // template length (0=contiguous, 16/18/21=spaced)
    uint8_t template_type = 0;  // 0 = contiguous

    // Fragment-indexing parameters parsed from the file name.
    uint32_t min_seq_length    = 0;
    uint32_t min_length_split  = 0;
    uint32_t overlap_length    = 0;
    uint64_t max_freq_build    = 1;  // absolute threshold (1 = disabled)
    uint32_t max_degen_expand  = 0;  // 0 / 1 = disabled
};

// Parsed components of an index file name.  vol_basename is "<DB>.<vol>"
// for per-volume files (.kix / .kpx / .ksx) or "<DB>" for per-DB files
// (.kvx / .khx); has_vol distinguishes the two.
struct IndexFilenameParts {
    std::string parent_dir;
    std::string db_basename;
    std::string vol_basename;
    bool        has_vol = false;
    int         k = 0;
    uint8_t     t = 0;
    uint8_t     template_type = 0;
    uint32_t    min_seq_length = 0;
    uint32_t    min_length_split = 0;
    uint32_t    overlap_length = 0;
    uint64_t    max_freq_build = 1;   // absolute threshold (1 = disabled)
    uint32_t    max_degen_expand = 0; // 0 / 1 = disabled
};

// Split ix_prefix into parent directory and DB name.
// e.g. "path/to/nt" -> {"path/to", "nt"}
struct IndexPrefixParts {
    std::string parent_dir;
    std::string db;
};
IndexPrefixParts parse_index_prefix(const std::string& ix_prefix);

// Build the file stem for a per-volume index file (kix/kpx/ksx).
//
// Output layout (omitted suffixes per the rules below):
//   <parent_dir>/<vol_basename>.k<k>[.t<t>].minlen<X>.minsplit<X>.ovllen<X>
//                .maxfreq<X>.maxexpand<X>[.cod|.opt]
//
// Suffix-omission rules:
//   t<X>        : t == 0
//   minlen<X>   : min_seq_length == 0
//   minsplit<X> : min_length_split == 0
//   ovllen<X>   : overlap_length == 0
//   maxfreq<X>  : max_freq_build == 1 (absolute threshold)
//   maxexpand<X>: max_degen_expand == 0 or 1
//   cod|opt     : t == 0
std::string index_file_stem(const std::string& parent_dir,
                            const std::string& vol_basename, int k,
                            uint8_t t,
                            uint8_t template_type,
                            uint32_t min_seq_length,
                            uint32_t min_length_split,
                            uint32_t overlap_length,
                            uint64_t max_freq_build,
                            uint32_t max_degen_expand);

// Build the .khx file path (one .khx per database; the second argument
// is the DB basename, no per-volume index inserted).
std::string khx_path_for(const std::string& parent_dir,
                          const std::string& db, int k,
                          uint8_t t,
                          uint8_t template_type,
                          uint32_t min_seq_length,
                          uint32_t min_length_split,
                          uint32_t overlap_length,
                          uint64_t max_freq_build,
                          uint32_t max_degen_expand);

// Parse the components encoded in an index file name (kix / kpx / ksx /
// khx / kvx).  `path` may be a full filename or just the stem.  Returns
// false if the name does not match the expected pattern.
bool parse_index_filename(const std::string& path, IndexFilenameParts& out);

// --- Shared index-variant selection ---------------------------------------
//
// An index variant is identified by 8 parameters encoded in the file name:
// k / t / template_type / min_seq_length / min_length_split / overlap_length /
// max_freq_build / max_degen_expand.  These helpers express the agreed
// selection semantics (see CLAUDE.md / the plan) so that ikafssnsearch,
// ikafssninfo, ikafssnserver and ikafssnclient all resolve variants the same
// way.  They operate on a small scalar view (VariantFields) so both the local
// DiscoveredVolume and the over-the-wire KmerGroupInfo can be fed through them.

struct VariantFields {
    int      k = 0;
    uint8_t  t = 0;
    uint8_t  template_type = 0;
    uint32_t min_seq_length = 0;
    uint32_t min_length_split = 0;
    uint32_t overlap_length = 0;
    uint64_t max_freq_build = 1;
    uint32_t max_degen_expand = 0;
};

inline VariantFields variant_fields_of(const DiscoveredVolume& v) {
    return VariantFields{v.k, v.t, v.template_type, v.min_seq_length,
                         v.min_length_split, v.overlap_length,
                         v.max_freq_build, v.max_degen_expand};
}

// A filter over the 8 identifying parameters.  Each `has_*` flag turns the
// corresponding field into an exact-match constraint; an unset flag is a
// wildcard.  `t` is special: `t_is_wildcard` selects the info/server semantics
// (unset -t matches any t), while a concrete `t` selects the search/client
// semantics (unset -t defaults to 0 = contiguous).  `max_degen_expand` is a
// filter only for ikafssninfo / ikafssnserver; ikafssnsearch / ikafssnclient
// leave `has_max_degen_expand` false and use it as a query-side tie-break.
struct VariantFilter {
    bool has_k = false;                int      k = 0;
    bool t_is_wildcard = false;        uint8_t  t = 0;
    bool has_template_type = false;    uint8_t  template_type = 0; // 0=con,1=cod,2=opt,3=both
    bool has_min_seq_length = false;   uint32_t min_seq_length = 0;
    bool has_min_length_split = false; uint32_t min_length_split = 0;
    bool has_overlap_length = false;   uint32_t overlap_length = 0;
    bool has_max_freq_build = false;   uint64_t max_freq_build = 1;
    bool has_max_degen_expand = false; uint32_t max_degen_expand = 0;
};

// True when `v` satisfies every constraint in `f`.
bool variant_matches(const VariantFilter& f, const VariantFields& v);

// Two variants are the "same" for the purpose of uniqueness when they agree on
// the 7 fields that exclude max_degen_expand.
bool same_variant7(const VariantFields& a, const VariantFields& b);

// Agreement on the 6 fields that exclude template_type and max_degen_expand
// (the both-mode coding/optimal pairing key).
bool same_variant6(const VariantFields& a, const VariantFields& b);

// Pick the max_degen_expand to use among the values present in a variant set:
// prefer the query value if it occurs, otherwise the largest present value.
uint32_t choose_max_degen_expand(const std::vector<uint32_t>& present,
                                 uint32_t query_value);

// How the template_type dimension is resolved (see resolve_variant_indices).
enum class TemplateResolveMode {
    // Candidates already share one template_type (contiguous / coding /
    // optimal).  Resolve to exactly one variant.
    kSingle,
    // Explicit `-template_type both`: a coding+optimal pair is required;
    // missing either side is an error.
    kBothRequired,
    // Template wildcard (`-t > 0` with no `-template_type`): use the
    // coding+optimal pair when both exist, otherwise the single side that
    // does (coding-only -> coding, optimal-only -> optimal).
    kBothOrSingle,
};

// Resolve filter-matched candidates to a single search target.
//
// `cands` are the VariantFields of every volume that passed variant_matches().
//
// kSingle: the candidates must collapse to exactly one variant identified by
// the 7 fields excluding max_degen_expand; the returned indices are that
// variant's volumes, with max_degen_expand tie-broken via query_max_degen_expand
// (query value preferred, else the largest present).
//
// kBothRequired / kBothOrSingle (pair case): coding and optimal must agree on
// the 6-field pairing key (k / t / min_seq_length / min_length_split /
// overlap_length / max_freq_build) AND on a common max_degen_expand value — no
// max_degen_expand mixing ("chanpon") between the two sides.  The shared
// max_degen_expand is tie-broken among the values present on BOTH sides
// (query value preferred, else the largest common one).  The returned indices
// cover the coding and optimal volumes of that single shared value.
//
// kBothOrSingle with only one side present resolves like kSingle on that side.
//
// On ambiguity, emptiness, a missing required side, or no common
// max_degen_expand, returns an empty vector and fills `err`.  The caller infers
// single-vs-pair from the template_types present in the returned indices.
std::vector<size_t> resolve_variant_indices(
    const std::vector<VariantFields>& cands, TemplateResolveMode mode,
    uint32_t query_max_degen_expand, std::string& err);

// Discover every index volume under `ix_prefix` from its .kvx manifests.
// filter_k: if > 0, only that k value; if 0, all available k values (a real
//   index always has k > 0, so 0 is unambiguously "any").  Callers that need
//   to narrow by the other variant parameters do so with variant_matches() on
//   the result (the t / template_type / indexing-field filtering lives there,
//   not here).
// Results are sorted by (k, t, template_type, volume_index) ascending.
std::vector<DiscoveredVolume> discover_volumes(
    const std::string& ix_prefix, int filter_k = 0);

} // namespace ikafssn
