#include "test_util.hpp"
#include "io/volume_discovery.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace ikafssn;

// Build a VariantFields with the 8 identifying parameters.
static VariantFields vf(int k, uint8_t t, uint8_t tt, uint32_t minlen,
                        uint32_t minsplit, uint32_t ovllen, uint64_t maxfreq,
                        uint32_t maxexpand) {
    return VariantFields{k, t, tt, minlen, minsplit, ovllen, maxfreq, maxexpand};
}

static void test_variant_matches_k() {
    std::fprintf(stderr, "-- test_variant_matches_k\n");
    VariantFilter f;                 // wildcard everything (t_is_wildcard=false, t=0)
    f.t_is_wildcard = true;          // info/server semantics: any t
    CHECK(variant_matches(f, vf(11, 0, 0, 64, 0, 0, 1, 4)));

    f.has_k = true; f.k = 11;
    CHECK(variant_matches(f, vf(11, 0, 0, 64, 0, 0, 1, 4)));
    CHECK(!variant_matches(f, vf(9, 0, 0, 64, 0, 0, 1, 4)));
}

static void test_variant_matches_t_and_template_type() {
    std::fprintf(stderr, "-- test_variant_matches_t_and_template_type\n");

    // -t 0 (search/client semantics: concrete t=0) selects contiguous only.
    {
        VariantFilter f;             // t_is_wildcard=false, t=0
        CHECK(variant_matches(f, vf(11, 0, 0, 64, 0, 0, 1, 4)));   // con
        CHECK(!variant_matches(f, vf(11, 16, 1, 64, 0, 0, 1, 4))); // cod
        CHECK(!variant_matches(f, vf(11, 16, 2, 64, 0, 0, 1, 4))); // opt
    }

    // -t 16 with no -template_type selects spaced types only (cod + opt),
    // never the contiguous index.
    {
        VariantFilter f; f.t = 16;
        CHECK(!variant_matches(f, vf(11, 0, 0, 64, 0, 0, 1, 4)));  // con excluded
        CHECK(variant_matches(f, vf(11, 16, 1, 64, 0, 0, 1, 4)));  // cod
        CHECK(variant_matches(f, vf(11, 16, 2, 64, 0, 0, 1, 4)));  // opt
        CHECK(!variant_matches(f, vf(11, 18, 1, 64, 0, 0, 1, 4))); // wrong t
    }

    // -t 16 -template_type coding (exact).
    {
        VariantFilter f; f.t = 16; f.has_template_type = true; f.template_type = 1;
        CHECK(variant_matches(f, vf(11, 16, 1, 64, 0, 0, 1, 4)));  // cod
        CHECK(!variant_matches(f, vf(11, 16, 2, 64, 0, 0, 1, 4))); // opt
    }

    // -template_type both matches whichever spaced types exist.
    {
        VariantFilter f; f.t = 16; f.has_template_type = true; f.template_type = 3;
        CHECK(variant_matches(f, vf(11, 16, 1, 64, 0, 0, 1, 4)));  // cod
        CHECK(variant_matches(f, vf(11, 16, 2, 64, 0, 0, 1, 4)));  // opt
        CHECK(!variant_matches(f, vf(11, 16, 0, 64, 0, 0, 1, 4))); // con never matches "both"
    }

    // -t wildcard (info/server, -t unset) matches con + cod + opt.
    {
        VariantFilter f; f.t_is_wildcard = true;
        CHECK(variant_matches(f, vf(11, 0, 0, 64, 0, 0, 1, 4)));
        CHECK(variant_matches(f, vf(11, 16, 1, 64, 0, 0, 1, 4)));
        CHECK(variant_matches(f, vf(11, 16, 2, 64, 0, 0, 1, 4)));
    }
}

static void test_variant_matches_indexing_fields() {
    std::fprintf(stderr, "-- test_variant_matches_indexing_fields\n");
    VariantFilter f; f.t_is_wildcard = true;
    f.has_min_seq_length = true;   f.min_seq_length = 64;
    f.has_min_length_split = true; f.min_length_split = 5000;
    f.has_overlap_length = true;   f.overlap_length = 250;
    f.has_max_freq_build = true;   f.max_freq_build = 1000;
    f.has_max_degen_expand = true; f.max_degen_expand = 8;

    CHECK(variant_matches(f, vf(11, 0, 0, 64, 5000, 250, 1000, 8)));
    CHECK(!variant_matches(f, vf(11, 0, 0, 100, 5000, 250, 1000, 8))); // minlen
    CHECK(!variant_matches(f, vf(11, 0, 0, 64, 6000, 250, 1000, 8)));  // minsplit
    CHECK(!variant_matches(f, vf(11, 0, 0, 64, 5000, 300, 1000, 8)));  // ovllen
    CHECK(!variant_matches(f, vf(11, 0, 0, 64, 5000, 250, 2000, 8)));  // maxfreq
    CHECK(!variant_matches(f, vf(11, 0, 0, 64, 5000, 250, 1000, 4)));  // maxexpand
}

static void test_choose_max_degen_expand() {
    std::fprintf(stderr, "-- test_choose_max_degen_expand\n");
    std::vector<uint32_t> present = {4, 8, 16};
    CHECK_EQ(choose_max_degen_expand(present, 8), 8u);   // query value present
    CHECK_EQ(choose_max_degen_expand(present, 32), 16u); // absent -> largest
    CHECK_EQ(choose_max_degen_expand({4, 8}, 0), 8u);    // 0 absent -> largest
}

static void test_resolve_single_unique() {
    std::fprintf(stderr, "-- test_resolve_single_unique\n");
    // Two con variants differing only in max_degen_expand (4 and 8).
    // Query maxexpand 8 must pick the maxexpand=8 variant.
    std::vector<VariantFields> cands = {
        vf(11, 0, 0, 64, 0, 0, 1, 4),
        vf(11, 0, 0, 64, 0, 0, 1, 8),
    };
    std::string err;
    auto sel = resolve_variant_indices(cands, TemplateResolveMode::kSingle, 8, err);
    CHECK(err.empty());
    CHECK_EQ(sel.size(), 1u);
    CHECK_EQ(cands[sel[0]].max_degen_expand, 8u);

    // No query match -> largest (8).
    err.clear();
    sel = resolve_variant_indices(cands, TemplateResolveMode::kSingle, 16, err);
    CHECK(err.empty());
    CHECK_EQ(sel.size(), 1u);
    CHECK_EQ(cands[sel[0]].max_degen_expand, 8u);
}

static void test_resolve_single_errors() {
    std::fprintf(stderr, "-- test_resolve_single_errors\n");
    // Empty.
    {
        std::string err;
        auto sel = resolve_variant_indices({}, TemplateResolveMode::kSingle, 16, err);
        CHECK(sel.empty());
        CHECK(!err.empty());
    }
    // Ambiguous: two distinct 7-tuples (minlen 64 vs 100).
    {
        std::vector<VariantFields> cands = {
            vf(11, 0, 0, 64, 0, 0, 1, 4),
            vf(11, 0, 0, 100, 0, 0, 1, 4),
        };
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kSingle, 16, err);
        CHECK(sel.empty());
        CHECK(err.find("ambiguous") != std::string::npos);
    }
}

static void test_resolve_both_required() {
    std::fprintf(stderr, "-- test_resolve_both_required\n");
    // cod + opt sharing the 6-tuple and a common maxexpand -> one pair.
    std::vector<VariantFields> cands = {
        vf(11, 16, 1, 64, 0, 0, 1, 4),
        vf(11, 16, 2, 64, 0, 0, 1, 4),
    };
    std::string err;
    auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothRequired, 4, err);
    CHECK(err.empty());
    CHECK_EQ(sel.size(), 2u);
    bool has_cod = false, has_opt = false;
    for (size_t i : sel) {
        if (cands[i].template_type == 1) has_cod = true;
        if (cands[i].template_type == 2) has_opt = true;
    }
    CHECK(has_cod && has_opt);

    // kBothRequired with only the coding side -> error.
    {
        std::vector<VariantFields> only_cod = {vf(11, 16, 1, 64, 0, 0, 1, 4)};
        std::string e;
        auto s = resolve_variant_indices(only_cod, TemplateResolveMode::kBothRequired, 4, e);
        CHECK(s.empty());
        CHECK(!e.empty());
    }

    // cod / opt disagree on the 6-tuple (minlen 64 vs 100) -> error.
    {
        std::vector<VariantFields> mismatch = {
            vf(11, 16, 1, 64, 0, 0, 1, 4),
            vf(11, 16, 2, 100, 0, 0, 1, 4),
        };
        std::string e;
        auto s = resolve_variant_indices(mismatch, TemplateResolveMode::kBothRequired, 4, e);
        CHECK(s.empty());
        CHECK(!e.empty());
    }
}

static void test_resolve_both_or_single() {
    std::fprintf(stderr, "-- test_resolve_both_or_single\n");
    // Both sides present -> pair.
    {
        std::vector<VariantFields> cands = {
            vf(11, 16, 1, 64, 0, 0, 1, 4),
            vf(11, 16, 2, 64, 0, 0, 1, 4),
        };
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothOrSingle, 4, err);
        CHECK(err.empty());
        CHECK_EQ(sel.size(), 2u);
    }
    // Only coding present -> single coding (NOT an error).
    {
        std::vector<VariantFields> cands = {vf(11, 16, 1, 64, 0, 0, 1, 4)};
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothOrSingle, 4, err);
        CHECK(err.empty());
        CHECK_EQ(sel.size(), 1u);
        CHECK_EQ(cands[sel[0]].template_type, 1u);
    }
    // Only optimal present -> single optimal.
    {
        std::vector<VariantFields> cands = {vf(11, 16, 2, 64, 0, 0, 1, 4)};
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothOrSingle, 4, err);
        CHECK(err.empty());
        CHECK_EQ(sel.size(), 1u);
        CHECK_EQ(cands[sel[0]].template_type, 2u);
    }
}

static void test_resolve_both_no_chanpon() {
    std::fprintf(stderr, "-- test_resolve_both_no_chanpon\n");
    // cod has {4,8}, opt has {8}.  The pair MUST use the common value 8 on
    // both sides (query 16 absent -> largest common = 8).
    {
        std::vector<VariantFields> cands = {
            vf(11, 16, 1, 64, 0, 0, 1, 4),
            vf(11, 16, 1, 64, 0, 0, 1, 8),
            vf(11, 16, 2, 64, 0, 0, 1, 8),
        };
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothRequired, 16, err);
        CHECK(err.empty());
        CHECK_EQ(sel.size(), 2u);
        for (size_t i : sel) CHECK_EQ(cands[i].max_degen_expand, 8u);  // no mixing
    }
    // Query value present on both sides -> chosen.
    {
        std::vector<VariantFields> cands = {
            vf(11, 16, 1, 64, 0, 0, 1, 4),
            vf(11, 16, 1, 64, 0, 0, 1, 8),
            vf(11, 16, 2, 64, 0, 0, 1, 4),
            vf(11, 16, 2, 64, 0, 0, 1, 8),
        };
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothRequired, 4, err);
        CHECK(err.empty());
        CHECK_EQ(sel.size(), 2u);
        for (size_t i : sel) CHECK_EQ(cands[i].max_degen_expand, 4u);
    }
    // No common max_degen_expand (cod={4}, opt={8}) -> error (would be mixing).
    {
        std::vector<VariantFields> cands = {
            vf(11, 16, 1, 64, 0, 0, 1, 4),
            vf(11, 16, 2, 64, 0, 0, 1, 8),
        };
        std::string err;
        auto sel = resolve_variant_indices(cands, TemplateResolveMode::kBothRequired, 16, err);
        CHECK(sel.empty());
        CHECK(!err.empty());
    }
}

int main() {
    test_variant_matches_k();
    test_variant_matches_t_and_template_type();
    test_variant_matches_indexing_fields();
    test_choose_max_degen_expand();
    test_resolve_single_unique();
    test_resolve_single_errors();
    test_resolve_both_required();
    test_resolve_both_or_single();
    test_resolve_both_no_chanpon();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
