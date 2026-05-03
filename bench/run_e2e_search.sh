#!/usr/bin/env bash
# bench/run_e2e_search.sh — measure ikafssnsearch wall time per SIMD tier.
#
# Usage:
#   bench/run_e2e_search.sh [build_dir] [--format json|text]
#
# Env:
#   DB     — path to BLAST DB (default: db/SSU_eukaryote_rRNA)
#   QUERY  — query FASTA (default: $TMPDIR/ikafssn_ssu_test/queries.fasta)
#   K      — k-mer length (default: 11)
#   T      — spaced seed template length (default: 0; 0 disables)
#   TEMPLATE_TYPE — coding/optimal/both (default: cod; only used when T>0)
#   IDX    — index directory (default: $TMPDIR/idx_e2e_search)
#   MODE   — search mode (default: 1)
#   COLD   — if "1", drop page cache between tier runs (requires sudo)
#   TIERS  — space-separated list of tiers to run (default: scalar sse42 avx2 avx512bw avx512vbmi avx512vbmi2)
#
# IKAFSSN_FORCE_SIMD is evaluated at process start, so each tier is timed in
# its own process invocation.
#
# JSON output schema (one object per run):
#   {
#     "schema_version": 1,
#     "db": "...", "k": 11, "t": 0, "template_type": "cod", "mode": 1,
#     "query": "...", "index": "...",
#     "host": "<hostname>", "kernel": "<uname -r>",
#     "binary": "<path>", "git_sha": "<sha or empty>",
#     "runs": [
#       { "tier": "avx2", "cold": false,
#         "wall_s": 1.23, "user_s": 0.98, "max_rss_kb": 12345,
#         "exit_code": 0 },
#       ...
#     ]
#   }

set -euo pipefail

BUILD=""
FORMAT="text"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --format)
            FORMAT=${2:-text}
            shift 2 ;;
        --format=*)
            FORMAT=${1#--format=}
            shift ;;
        *)
            if [[ -z "$BUILD" ]]; then BUILD=$1; fi
            shift ;;
    esac
done
BUILD=${BUILD:-build}

DB=${DB:-db/SSU_eukaryote_rRNA}
QUERY=${QUERY:-/tmp/ikafssn_ssu_test/queries.fasta}
K=${K:-11}
T=${T:-0}
TEMPLATE_TYPE=${TEMPLATE_TYPE:-cod}
IDX=${IDX:-/tmp/idx_e2e_search}
MODE=${MODE:-1}
COLD=${COLD:-0}
TIERS_STR=${TIERS:-"scalar sse42 avx2 avx512bw avx512vbmi avx512vbmi2"}
read -r -a TIERS <<< "$TIERS_STR"

if [[ ! -x "$BUILD/src/ikafssnsearch" ]]; then
    echo "ikafssnsearch not found at $BUILD/src/ikafssnsearch" >&2
    exit 1
fi
if [[ ! -x "$BUILD/src/ikafssnindex" ]]; then
    echo "ikafssnindex not found at $BUILD/src/ikafssnindex" >&2
    exit 1
fi
if [[ ! -f "$QUERY" ]]; then
    echo "QUERY file not found: $QUERY" >&2
    echo "Hint: run test/scripts/setup_ssu_testdata.sh to generate it." >&2
    exit 1
fi

# Build the index once (any tier — index files are bit-exact across tiers).
INDEX_BUILD_ARGS=(-db "$DB" -k "$K" -o "$IDX/" -mode "$MODE")
if [[ "$T" -gt 0 ]]; then
    INDEX_BUILD_ARGS+=(-t "$T" -template_type "$TEMPLATE_TYPE")
fi
if [[ ! -d "$IDX" ]] || [[ -z "$(ls -A "$IDX" 2>/dev/null)" ]]; then
    rm -rf "$IDX" && mkdir -p "$IDX"
    if [[ "$FORMAT" != "json" ]]; then
        echo "=== building index k=$K t=$T template_type=$TEMPLATE_TYPE mode=$MODE ==="
    fi >&2
    "$BUILD/src/ikafssnindex" "${INDEX_BUILD_ARGS[@]}" >&2
fi

# Drop the page cache once at the start when COLD=1 so the first tier
# observes a fully cold cache.  Cache state across tier runs depends on
# whether COLD=1 is also passed to subsequent runs.
drop_cache_if_cold() {
    if [[ "$COLD" == "1" ]]; then
        if command -v sudo >/dev/null 2>&1; then
            sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || true
        else
            echo "COLD=1 set but sudo not available; cache not dropped" >&2
        fi
    fi
}

SEARCH_ARGS=(-ix "$IDX/$(basename "$DB")" -k "$K" -mode "$MODE" -query "$QUERY" -o /dev/null)
if [[ "$T" -gt 0 ]]; then
    SEARCH_ARGS+=(-t "$T" -template_type "$TEMPLATE_TYPE")
fi

# Measure one tier. Echoes "wall_s user_s max_rss_kb exit_code" on stdout.
measure_tier() {
    local tier=$1
    local time_log
    time_log=$(mktemp)
    set +e
    IKAFSSN_FORCE_SIMD="$tier" /usr/bin/time -f '__BENCH__ %e %U %M' \
        "$BUILD/src/ikafssnsearch" "${SEARCH_ARGS[@]}" \
        >/dev/null 2>"$time_log"
    local exit_code=$?
    set -e
    local line
    line=$(grep -E '^__BENCH__ ' "$time_log" | tail -1 || true)
    rm -f "$time_log"
    if [[ -z "$line" ]]; then
        echo "0 0 0 $exit_code"
        return
    fi
    # __BENCH__ <wall> <user> <max_rss>
    set -- $line
    echo "$2 $3 $4 $exit_code"
}

if [[ "$FORMAT" == "json" ]]; then
    HOST=$(hostname 2>/dev/null || echo "")
    KERNEL=$(uname -r 2>/dev/null || echo "")
    GIT_SHA=$(git -C "$(dirname "$0")/.." rev-parse --short=12 HEAD 2>/dev/null || echo "")
    BIN_REAL=$(readlink -f "$BUILD/src/ikafssnsearch")
    QUERY_REAL=$(readlink -f "$QUERY")
    IDX_REAL=$(readlink -f "$IDX")
    DB_REAL=$(readlink -f "$DB" 2>/dev/null || echo "$DB")

    printf '{\n'
    printf '  "schema_version": 1,\n'
    printf '  "db": "%s",\n' "$DB_REAL"
    printf '  "k": %s,\n' "$K"
    printf '  "t": %s,\n' "$T"
    printf '  "template_type": "%s",\n' "$TEMPLATE_TYPE"
    printf '  "mode": %s,\n' "$MODE"
    printf '  "cold": %s,\n' "$([[ "$COLD" == "1" ]] && echo true || echo false)"
    printf '  "query": "%s",\n' "$QUERY_REAL"
    printf '  "index": "%s",\n' "$IDX_REAL"
    printf '  "host": "%s",\n' "$HOST"
    printf '  "kernel": "%s",\n' "$KERNEL"
    printf '  "binary": "%s",\n' "$BIN_REAL"
    printf '  "git_sha": "%s",\n' "$GIT_SHA"
    printf '  "runs": [\n'
    first=1
    for t in "${TIERS[@]}"; do
        drop_cache_if_cold
        read -r wall user rss exit_code < <(measure_tier "$t")
        if [[ $first -eq 1 ]]; then
            first=0
        else
            printf ',\n'
        fi
        printf '    {"tier": "%s", "cold": %s, "wall_s": %s, "user_s": %s, "max_rss_kb": %s, "exit_code": %s}' \
            "$t" \
            "$([[ "$COLD" == "1" ]] && echo true || echo false)" \
            "$wall" "$user" "$rss" "$exit_code"
    done
    printf '\n  ]\n'
    printf '}\n'
else
    echo "=== ikafssnsearch end-to-end on $DB (k=$K, t=$T, template_type=$TEMPLATE_TYPE, mode=$MODE) ==="
    for t in "${TIERS[@]}"; do
        drop_cache_if_cold
        printf -- "--- tier=%-15s ---\n" "$t"
        read -r wall user rss exit_code < <(measure_tier "$t")
        printf '  wall=%s s  user=%s s  max_rss=%s KB  exit=%s\n' \
            "$wall" "$user" "$rss" "$exit_code"
    done
fi
