#!/usr/bin/env bash
# bench/run_cold_e2e.sh — cold-start A/B comparison of two ikafssnsearch builds.
#
# Companion to run_warm_e2e.sh.  That script keeps the index in the page cache
# and attributes CPU time per stage; this one does the opposite — it evicts the
# index before every run so the search is I/O bound, which is the regime a
# large index (one that cannot stay resident) actually runs in.
#
# Eviction is posix_fadvise(POSIX_FADV_DONTNEED) over the index files, not
# /proc/sys/vm/drop_caches: it needs no privileges and leaves every other
# process's cache alone.  Without it the page cache carries over between runs
# and the same binary varies by 2x or more, which buries any real difference.
#
# Runs alternate A B A B ... so a residual cache trend hits both builds alike.
#
# Linux only: it needs /proc/diskstats and Python's os.posix_fadvise, and
# /usr/bin/time -v as run_warm_e2e.sh does.
#
# Usage:
#   bench/run_cold_e2e.sh [--format tsv|summary]
#
# Env:
#   BIN_A         — baseline binary (required, e.g. /usr/bin/ikafssnsearch)
#   BIN_B         — binary under test (required, e.g. build/src/ikafssnsearch)
#   IX            — index prefix; the index is NOT built here (required)
#   QUERY         — query FASTA (required)
#   EVICT_GLOB    — glob of index files to evict before each run
#                   (default: "${IX}"*)
#   DISKS         — space-separated /proc/diskstats device names to sum read
#                   sectors over.  Use the RAID members, not the md device
#                   (default: every nvme?n? and sd? line present)
#   REPS          — A/B pairs to run (default: 3).  Every run is cold, so none
#                   is discarded.
#   OUTDIR        — where per-run logs and result files land (default: mktemp -d)
#   EXTRA_ARGS    — extra ikafssnsearch arguments, word-split
#
# Exit codes 0 and 2 both count as success: 2 only means some query was
# skipped (e.g. shorter than -min_query_length).
#
# Judge a change on wall time plus the I/O counters: a change that leaves the
# read volume and the major fault count alone has not altered which posting
# bytes the search touches.  Output equivalence is reported as line count and
# sorted sha256 — at this scale mode 1 does not fix the row order within a
# query, so a byte comparison is not meaningful.

set -u

FORMAT=tsv
[ $# -gt 0 ] && case "$1" in
    --format) FORMAT=${2:-tsv} ;;
    *) echo "usage: $0 [--format tsv|summary]" >&2; exit 1 ;;
esac

: "${BIN_A:?set BIN_A to the baseline binary}"
: "${BIN_B:?set BIN_B to the binary under test}"
: "${IX:?set IX to the index prefix}"
: "${QUERY:?set QUERY to the query FASTA}"
EVICT_GLOB=${EVICT_GLOB:-"${IX}*"}
REPS=${REPS:-3}
OUTDIR=${OUTDIR:-$(mktemp -d)}
EXTRA_ARGS=${EXTRA_ARGS:-}
mkdir -p "$OUTDIR"

if [ -z "${DISKS:-}" ]; then
    DISKS=$(awk '$3 ~ /^(nvme[0-9]+n[0-9]+|sd[a-z])$/ {print $3}' /proc/diskstats)
fi

read_sectors() {
    awk -v want=" $(echo $DISKS) " \
        '{ if (index(want, " " $3 " ")) t += $6 } END { print t+0 }' /proc/diskstats
}

evict() {
    python3 - "$EVICT_GLOB" <<'PY'
import glob, os, sys
for path in sorted(glob.glob(sys.argv[1])):
    if not os.path.isfile(path):
        continue
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)
PY
}

one_run() {  # one_run <tag> <binary>
    local tag=$1 bin=$2 r0 r1
    evict
    r0=$(read_sectors)
    /usr/bin/time -v -o "$OUTDIR/$tag.time" \
        "$bin" -ix "$IX" -query "$QUERY" -o "$OUTDIR/$tag.tsv" $EXTRA_ARGS \
        > "$OUTDIR/$tag.log" 2>&1
    local rc=$?
    r1=$(read_sectors)

    local t=$OUTDIR/$tag.time
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$tag" \
        "$(awk -F'): ' '/Elapsed \(wall clock\)/{print $2}' "$t" |
          awk -F: '{ print (NF==3) ? $1*3600+$2*60+$3 : $1*60+$2 }')" \
        "$(awk -F': ' '/User time/{print $2}' "$t")" \
        "$(awk -F': ' '/System time/{print $2}' "$t")" \
        "$(awk -F': ' '/Maximum resident set size/{print $2}' "$t")" \
        "$(awk -F': ' '/Major \(requiring I\/O\) page faults/{print $2}' "$t")" \
        "$(awk -v a="$r0" -v b="$r1" 'BEGIN{printf "%.2f", (b-a)*512/1073741824}')" \
        "$(wc -l < "$OUTDIR/$tag.tsv")" \
        "$(sort "$OUTDIR/$tag.tsv" | sha256sum | cut -c1-16)" \
        "$rc"
}

RESULTS=$OUTDIR/runs.tsv
printf 'tag\twall_s\tuser_s\tsys_s\tmax_rss_kb\tmajor_faults\tread_gb\tlines\tsha256_16\texit\n' > "$RESULTS"
for i in $(seq 1 "$REPS"); do
    one_run "a_$i" "$BIN_A" >> "$RESULTS"
    one_run "b_$i" "$BIN_B" >> "$RESULTS"
done

if [ "$FORMAT" = summary ]; then
    python3 - "$RESULTS" <<'PY'
import statistics as st, sys
hdr = "tag wall_s user_s sys_s max_rss_kb major_faults read_gb lines sha256_16 exit".split()
rows = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in open(sys.argv[1])][1:]
a = [r for r in rows if r["tag"].startswith("a")]
b = [r for r in rows if r["tag"].startswith("b")]
out = {(r["lines"], r["sha256_16"], r["exit"]) for r in rows}
print("output: %s -> %s" % (sorted(out), "IDENTICAL" if len(out) == 1 else "*** DIFFERS ***"))
for key, nd in (("wall_s", 2), ("user_s", 1), ("sys_s", 1),
                ("max_rss_kb", 0), ("major_faults", 0), ("read_gb", 2)):
    x = st.median(float(r[key]) for r in a)
    y = st.median(float(r[key]) for r in b)
    print("%-13s A=%-14.*f B=%-14.*f %+6.2f%%"
          % (key, nd, x, nd, y, (y - x) / x * 100 if x else 0.0))
cores = len(__import__("os").sched_getaffinity(0))
for name, rs in (("A", a), ("B", b)):
    w = st.median(float(r["wall_s"]) for r in rs)
    c = st.median(float(r["user_s"]) + float(r["sys_s"]) for r in rs)
    print("%s CPU utilisation over %d cores: %.1f%% (under ~50%% means I/O bound)"
          % (name, cores, c / w / cores * 100))
PY
else
    cat "$RESULTS"
fi
echo "raw logs: $OUTDIR" >&2
