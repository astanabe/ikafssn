#!/usr/bin/env bash
# bench/run_e2e.sh — end-to-end ikafssnindex wall-time per SIMD tier.
#
# Usage:  bench/run_e2e.sh [<build_dir> [<db_path> [<k>]]]
# Defaults: build/  db/SSU_eukaryote_rRNA  11
#
# IKAFSSN_FORCE_SIMD is consulted at process startup, so each tier requires
# a separate ikafssnindex invocation.

set -e
BUILD="${1:-build}"
DB="${2:-db/SSU_eukaryote_rRNA}"
K="${3:-11}"

INDEX_BIN="${BUILD}/src/ikafssnindex"
if [[ ! -x "${INDEX_BIN}" ]]; then
    echo "ERROR: ${INDEX_BIN} not found or not executable" >&2
    echo "       Build with -DCMAKE_BUILD_TYPE=Release first" >&2
    exit 1
fi

if [[ ! -f "${DB}.nin" && ! -f "${DB}.00.nin" ]]; then
    echo "ERROR: BLAST DB at '${DB}' not found" >&2
    exit 1
fi

# Tier list. Auto-skip tiers the host cannot satisfy (force_simd_cap downgrades
# silently to scalar; we'd produce duplicate scalar timings otherwise).
ARCH="$(uname -m)"
case "${ARCH}" in
    x86_64|i*86)
        TIERS=(scalar sse42 avx2 avx512bw avx512vbmi avx512vbmi2)
        ;;
    aarch64|arm64)
        TIERS=(scalar neon sve sve2)
        ;;
    *)
        TIERS=(scalar)
        ;;
esac

printf "=== ikafssnindex end-to-end on %s (k=%s) ===\n" "${DB}" "${K}"
printf "%-15s  %-10s  %-10s  %-12s\n" "tier" "wall(s)" "user(s)" "rss(KB)"

for tier in "${TIERS[@]}"; do
    OUT="/tmp/bench_idx_${tier}"
    rm -rf "${OUT}" && mkdir -p "${OUT}"

    # Drop page cache between runs would require root; skip and rely on
    # cumulative wall-time delta as advisory rather than absolute.
    OUTPUT=$(IKAFSSN_FORCE_SIMD="${tier}" /usr/bin/time -f '%e %U %M' \
                 "${INDEX_BIN}" -db "${DB}" -k "${K}" -o "${OUT}/" -v 0 \
                 2>&1 1>/dev/null)
    METRICS=$(printf "%s\n" "${OUTPUT}" | tail -n 1)
    WALL=$(printf "%s" "${METRICS}" | awk '{print $1}')
    USER=$(printf "%s" "${METRICS}" | awk '{print $2}')
    RSS=$(printf "%s" "${METRICS}" | awk '{print $3}')
    printf "%-15s  %-10s  %-10s  %-12s\n" "${tier}" "${WALL}" "${USER}" "${RSS}"
done
