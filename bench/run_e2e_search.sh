#!/usr/bin/env bash
# bench/run_e2e_search.sh — measure ikafssnsearch wall time per SIMD tier.
#
# Usage:
#   bench/run_e2e_search.sh [build_dir]
#
# Env:
#   DB     — path to BLAST DB (default: db/SSU_eukaryote_rRNA)
#   QUERY  — query FASTA (default: $TMPDIR/ikafssn_ssu_test/queries.fasta)
#   K      — k-mer length (default: 11)
#   IDX    — index directory (default: $TMPDIR/idx_e2e_search)
#   MODE   — search mode (default: 1)
#
# IKAFSSN_FORCE_SIMD is evaluated at process start, so each tier is timed in
# its own process invocation.

set -euo pipefail

BUILD=${1:-build}
DB=${DB:-db/SSU_eukaryote_rRNA}
QUERY=${QUERY:-/tmp/ikafssn_ssu_test/queries.fasta}
K=${K:-11}
IDX=${IDX:-/tmp/idx_e2e_search}
MODE=${MODE:-1}

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
if [[ ! -d "$IDX" ]] || [[ -z "$(ls -A "$IDX" 2>/dev/null)" ]]; then
    rm -rf "$IDX" && mkdir -p "$IDX"
    echo "=== building index k=$K mode=$MODE ==="
    "$BUILD/src/ikafssnindex" -ix "$DB" -k "$K" -o "$IDX/" -mode "$MODE"
fi

TIERS=(scalar sse42 avx2 avx512bw avx512vbmi avx512vbmi2)

echo "=== ikafssnsearch end-to-end on $DB (k=$K, mode=$MODE) ==="
for t in "${TIERS[@]}"; do
    printf -- "--- tier=%-15s ---\n" "$t"
    IKAFSSN_FORCE_SIMD="$t" /usr/bin/time -f '%e s wall, %U s user, %M KB max-rss' \
        "$BUILD/src/ikafssnsearch" \
            -ix "$IDX/$(basename "$DB")" -k "$K" -mode "$MODE" \
            -q "$QUERY" -o /dev/null \
        2>&1 | tail -3 || true
done
