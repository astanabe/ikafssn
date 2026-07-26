#!/usr/bin/env bash
# bench/run_warm_e2e.sh — repeated warm end-to-end ikafssnsearch measurement.
#
# Runs one fixed search command REPS times against an already-built index,
# discards the first run (page-cache warm-up) and reports the median and
# minimum of the rest.  Per run it collects /usr/bin/time -v counters and the
# per-stage breakdown that ikafssnsearch prints ("Timing run_search (s):" /
# "Timing overall (s):"), so a change can be attributed to a stage instead of
# being buried in the end-to-end wall.  The result file is hashed on every run
# so output equivalence against a baseline can be checked directly.
#
# Companion to bench/run_e2e_search.sh, which compares SIMD tiers on a
# self-built index instead.
#
# Usage:
#   bench/run_warm_e2e.sh [--format json|text]
#
# Env:
#   BUILD         — build directory (default: build)
#   IX            — index prefix; the index is NOT built here (required)
#   QUERY         — query FASTA (required)
#   DB            — BLAST DB, passed as -db (mode 3 only; default: empty)
#   K             — k-mer length (default: 11)
#   T             — spaced seed template length (default: 21; 0 disables)
#   TEMPLATE_TYPE — coding/optimal/both (default: both; used when T>0)
#   MODE          — search mode (default: 2)
#   STRAND        — -strand value (default: 2)
#   MIN_SCORE     — -stage1_min_score value (default: 0.1)
#   M_SWEEP       — -stage1_max_nhit_per_volume values to sweep
#                   (default: "100 1000 10000"; the token "default" or an
#                   empty value means one point with the flag omitted)
#   REPS          — runs per sweep point, first discarded (default: 5)
#   PERF          — 1 to add one extra perf stat run per point (default: 1)
#   OUTDIR        — where raw logs and result files land (default: mktemp -d)
#   EXTRA_ARGS    — extra ikafssnsearch arguments, word-split
#
# Exit codes 0 and 2 both count as success: 2 only means some query was
# skipped (e.g. shorter than -min_query_length), which is expected for
# tsant500.fasta.
#
# JSON output schema:
#   {
#     "schema_version": 1,
#     "index": "...", "query": "...", "db": "...",
#     "k": 11, "t": 21, "template_type": "both", "mode": 2,
#     "strand": 2, "min_score": "0.1", "reps": 5, "timed_reps": 4,
#     "host": "...", "kernel": "...", "binary": "...", "git_sha": "...",
#     "outdir": "...",
#     "points": [
#       { "m": "1000", "exit_code": 0, "lines": 12345,
#         "sha256": "...", "sha256_stable": true,
#         "perf": {"cycles": 1.0e10, "instructions": 2.0e10, ...},
#         "metrics": { "wall_s":       {"median": 1.23, "min": 1.20},
#                      "s1_compute":   {"median": 0.10, "min": 0.09},
#                      "overall_total":{"median": 1.30, "min": 1.28} } },
#       ...
#     ]
#   }
#
# Metric keys are wall_s / user_s / max_rss_kb / minor_faults, every key of
# the "Timing run_search (s):" line verbatim (s1_open, s1_compute, s1_fold,
# s1_intotal, s2_open, s2a, s2b, s2_free, dedup, parent_topn, s2_intotal,
# total), and every key of the "Timing overall (s):" line prefixed with
# "overall_".

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

FORMAT="text"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --format)    FORMAT=${2:-text}; shift 2 ;;
        --format=*)  FORMAT=${1#--format=}; shift ;;
        *)           echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

BUILD=${BUILD:-build}
IX=${IX:-}
QUERY=${QUERY:-}
DB=${DB:-}
K=${K:-11}
T=${T:-21}
TEMPLATE_TYPE=${TEMPLATE_TYPE:-both}
MODE=${MODE:-2}
STRAND=${STRAND:-2}
MIN_SCORE=${MIN_SCORE:-0.1}
M_SWEEP=${M_SWEEP-"100 1000 10000"}
REPS=${REPS:-5}
PERF=${PERF:-1}
OUTDIR=${OUTDIR:-}
EXTRA_ARGS=${EXTRA_ARGS:-}

PERF_EVENTS=${PERF_EVENTS:-cycles,instructions,dTLB-load-misses,ls_l1_d_tlb_miss.all,ls_l1_d_tlb_miss.all_l2_miss,minor-faults}

SEARCH_BIN="$BUILD/src/ikafssnsearch"
if [[ ! -x "$SEARCH_BIN" ]]; then
    echo "ikafssnsearch not found at $SEARCH_BIN" >&2
    exit 1
fi
if [[ -z "$IX" ]]; then
    echo "IX is required (index prefix; this script does not build indexes)" >&2
    exit 1
fi
if [[ -z "$QUERY" || ! -f "$QUERY" ]]; then
    echo "QUERY file not found: ${QUERY:-<unset>}" >&2
    exit 1
fi
if [[ "$REPS" -lt 2 ]]; then
    echo "REPS must be >= 2 (the first run is discarded)" >&2
    exit 1
fi

if [[ -z "$OUTDIR" ]]; then
    OUTDIR=$(mktemp -d -t ikafssn_warm_e2e.XXXXXX)
fi
mkdir -p "$OUTDIR"

if [[ -z "$M_SWEEP" ]]; then M_SWEEP="default"; fi
read -r -a M_VALUES <<< "$M_SWEEP"

BASE_ARGS=(-ix "$IX" -k "$K" -mode "$MODE" -query "$QUERY"
           -strand "$STRAND" -stage1_min_score "$MIN_SCORE")
if [[ "$T" -gt 0 ]]; then
    BASE_ARGS+=(-t "$T" -template_type "$TEMPLATE_TYPE")
fi
if [[ -n "$DB" ]]; then
    BASE_ARGS+=(-db "$DB")
fi
if [[ -n "$EXTRA_ARGS" ]]; then
    read -r -a EXTRA_ARR <<< "$EXTRA_ARGS"
    BASE_ARGS+=("${EXTRA_ARR[@]}")
fi

note() { if [[ "$FORMAT" != "json" ]]; then echo "$@" >&2; fi }

# Median and minimum per key.  Reads "key value" lines, writes
# "key median min" lines sorted by key.
agg_stats() {
    awk '
    { key = $1; n[key]++; a[key, n[key]] = $2 + 0 }
    END {
        for (k in n) {
            m = n[k];
            for (i = 1; i <= m; i++) t[i] = a[k, i];
            for (i = 2; i <= m; i++) {
                x = t[i]; j = i - 1;
                while (j >= 1 && t[j] > x) { t[j + 1] = t[j]; j-- }
                t[j + 1] = x;
            }
            med = (m % 2 == 1) ? t[(m + 1) / 2] : (t[m / 2] + t[m / 2 + 1]) / 2;
            printf "%s %.6f %.6f\n", k, med, t[1];
        }
    }' | sort -k1,1
}

# Turn one run log into "key value" lines: /usr/bin/time -v counters plus
# every key=value pair of the two Timing lines.
extract_metrics() {
    local log=$1
    python3 "$SCRIPT_DIR/parse_time.py" "$log" | python3 -c '
import json, sys
r = json.load(sys.stdin)
for k in ("wall_s", "user_s", "max_rss_kb", "minor_faults"):
    print(k, r[k])
'
    sed -n 's/.*Timing run_search (s)://p' "$log" | tail -1 \
        | tr ' ' '\n' | sed -n 's/=/ /p'
    sed -n 's/.*Timing overall (s)://p' "$log" | tail -1 \
        | tr ' ' '\n' | sed -n 's/^\([a-z0-9_]*\)=\(.*\)$/overall_\1 \2/p'
}

run_once() {
    local log=$1 out=$2; shift 2
    local rc=0
    set +e
    /usr/bin/time -v "$SEARCH_BIN" "$@" -o "$out" >/dev/null 2>"$log"
    rc=$?
    set -e
    if [[ $rc -ne 0 && $rc -ne 2 ]]; then
        echo "ikafssnsearch failed (exit $rc); see $log" >&2
        tail -20 "$log" >&2
        exit 1
    fi
    echo "$rc"
}

emit_point_json() {
    local m=$1 dir=$2 rc=$3 lines=$4 sha=$5 stable=$6
    printf '    {"m": "%s", "exit_code": %s, "lines": %s, "sha256": "%s", "sha256_stable": %s' \
        "$m" "$rc" "$lines" "$sha" "$stable"
    if [[ -s "$dir/perf.csv" ]]; then
        printf ', "perf": {'
        awk -F',' '!/^#/ && NF >= 3 && $3 != "" {
            gsub(/ /, "", $1);
            if ($1 == "<notcounted>" || $1 == "<notsupported>") next;
            if (seen++) printf ", ";
            printf "\"%s\": %s", $3, $1;
        }' "$dir/perf.csv"
        printf '}'
    fi
    printf ', "metrics": {'
    awk '{ if (NR > 1) printf ", ";
           printf "\"%s\": {\"median\": %s, \"min\": %s}", $1, $2, $3 }' \
        "$dir/stats.txt"
    printf '}}'
}

emit_point_text() {
    local m=$1 dir=$2 rc=$3 lines=$4 sha=$5 stable=$6
    printf -- "--- M=%-8s exit=%s lines=%s sha256=%s%s ---\n" \
        "$m" "$rc" "$lines" "${sha:0:16}" \
        "$([[ "$stable" == "true" ]] && echo "" || echo " (UNSTABLE)")"
    awk '{ printf "  %-16s median=%-12s min=%s\n", $1, $2, $3 }' "$dir/stats.txt"
    if [[ -s "$dir/perf.csv" ]]; then
        awk -F',' '!/^#/ && NF >= 3 && $3 != "" {
            gsub(/ /, "", $1);
            if ($1 == "<notcounted>" || $1 == "<notsupported>") next;
            printf "  perf %-16s %s\n", $3, $1;
        }' "$dir/perf.csv"
    fi
}

measure_point() {
    local m=$1
    local dir="$OUTDIR/m_$m"
    rm -rf "$dir" && mkdir -p "$dir"

    local args=("${BASE_ARGS[@]}")
    if [[ "$m" != "default" ]]; then
        args+=(-stage1_max_nhit_per_volume "$m")
    fi

    local out="$dir/out.tsv"
    local rc=0 lines=0 sha="" stable="true"
    : > "$dir/metrics.txt"

    for ((r = 1; r <= REPS; r++)); do
        note "  M=$m rep $r/$REPS"
        rc=$(run_once "$dir/rep$r.log" "$out" "${args[@]}")
        local this_lines this_sha
        this_lines=$(wc -l < "$out")
        this_sha=$(sort "$out" | sha256sum | cut -d' ' -f1)
        # The first run is the page-cache warm-up and is not timed, but its
        # output still has to match the timed runs.
        if [[ -z "$sha" ]]; then
            lines=$this_lines
            sha=$this_sha
        elif [[ "$this_sha" != "$sha" || "$this_lines" != "$lines" ]]; then
            stable="false"
        fi
        if [[ $r -gt 1 ]]; then
            extract_metrics "$dir/rep$r.log" >> "$dir/metrics.txt"
        fi
    done

    agg_stats < "$dir/metrics.txt" > "$dir/stats.txt"

    # perf runs outside the timed set so its overhead never enters the median.
    if [[ "$PERF" == "1" ]] && command -v perf >/dev/null 2>&1; then
        note "  M=$m perf run"
        set +e
        perf stat -x, -e "$PERF_EVENTS" -o "$dir/perf.csv" \
            "$SEARCH_BIN" "${args[@]}" -o "$out" \
            >/dev/null 2>"$dir/perf.log"
        local prc=$?
        set -e
        if [[ $prc -ne 0 && $prc -ne 2 ]]; then
            note "  perf stat run failed (exit $prc); see $dir/perf.log"
            : > "$dir/perf.csv"
        fi
    fi

    if [[ "$stable" != "true" ]]; then
        note "  WARNING: output differed across repetitions at M=$m"
    fi

    if [[ "$FORMAT" == "json" ]]; then
        emit_point_json "$m" "$dir" "$rc" "$lines" "$sha" "$stable"
    else
        emit_point_text "$m" "$dir" "$rc" "$lines" "$sha" "$stable"
    fi
}

if [[ "$FORMAT" == "json" ]]; then
    HOST=$(hostname 2>/dev/null || echo "")
    KERNEL=$(uname -r 2>/dev/null || echo "")
    GIT_SHA=$(git -C "$SCRIPT_DIR/.." rev-parse --short=12 HEAD 2>/dev/null || echo "")
    BIN_REAL=$(readlink -f "$SEARCH_BIN")
    QUERY_REAL=$(readlink -f "$QUERY")

    printf '{\n'
    printf '  "schema_version": 1,\n'
    printf '  "index": "%s",\n' "$IX"
    printf '  "query": "%s",\n' "$QUERY_REAL"
    printf '  "db": "%s",\n' "$DB"
    printf '  "k": %s,\n' "$K"
    printf '  "t": %s,\n' "$T"
    printf '  "template_type": "%s",\n' "$TEMPLATE_TYPE"
    printf '  "mode": %s,\n' "$MODE"
    printf '  "strand": %s,\n' "$STRAND"
    printf '  "min_score": "%s",\n' "$MIN_SCORE"
    printf '  "reps": %s,\n' "$REPS"
    printf '  "timed_reps": %s,\n' "$((REPS - 1))"
    printf '  "host": "%s",\n' "$HOST"
    printf '  "kernel": "%s",\n' "$KERNEL"
    printf '  "binary": "%s",\n' "$BIN_REAL"
    printf '  "git_sha": "%s",\n' "$GIT_SHA"
    printf '  "outdir": "%s",\n' "$OUTDIR"
    printf '  "points": [\n'
    first=1
    for m in "${M_VALUES[@]}"; do
        if [[ $first -eq 1 ]]; then first=0; else printf ',\n'; fi
        measure_point "$m"
    done
    printf '\n  ]\n'
    printf '}\n'
else
    echo "=== warm e2e: $IX (k=$K, t=$T, template_type=$TEMPLATE_TYPE, mode=$MODE,"
    echo "    strand=$STRAND, stage1_min_score=$MIN_SCORE, reps=$REPS) ==="
    echo "    query=$QUERY"
    echo "    outdir=$OUTDIR"
    for m in "${M_VALUES[@]}"; do
        measure_point "$m"
    done
fi
