#!/usr/bin/env python3
"""Parse /usr/bin/time -v output and emit a single-run JSON record.

Companion utility to bench/run_e2e_search.sh.  Most of the bench harness
emits JSON directly via shell printf, but standalone profiling sessions
that invoke /usr/bin/time -v on ikafssnsearch can pipe the output through
this script to land in the same schema.

Usage:
    /usr/bin/time -v ikafssnsearch ... 2> time.log
    python3 bench/parse_time.py --tier avx2 time.log > run.json

Schema (matches bench/run_e2e_search.sh single-run object):
    {"tier": "avx2", "wall_s": 1.23, "user_s": 0.98, "max_rss_kb": 12345,
     "exit_code": 0}
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


WALL_RE = re.compile(
    r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)"
)
USER_RE = re.compile(r"User time \(seconds\):\s*([0-9.]+)")
RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*([0-9]+)")
EXIT_RE = re.compile(r"Exit status:\s*([0-9]+)")


def parse_wall(s: str) -> float:
    """Parse ``Elapsed (wall clock) time`` value (h:mm:ss or m:ss.ss)."""
    parts = s.split(":")
    if len(parts) == 3:
        h, m, sec = parts
        return int(h) * 3600 + int(m) * 60 + float(sec)
    if len(parts) == 2:
        m, sec = parts
        return int(m) * 60 + float(sec)
    return float(parts[0])


def parse_time_log(text: str) -> dict[str, float | int]:
    out: dict[str, float | int] = {
        "wall_s": 0.0,
        "user_s": 0.0,
        "max_rss_kb": 0,
        "exit_code": 0,
    }
    m = WALL_RE.search(text)
    if m:
        out["wall_s"] = parse_wall(m.group(1))
    m = USER_RE.search(text)
    if m:
        out["user_s"] = float(m.group(1))
    m = RSS_RE.search(text)
    if m:
        out["max_rss_kb"] = int(m.group(1))
    m = EXIT_RE.search(text)
    if m:
        out["exit_code"] = int(m.group(1))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("logfile", type=Path, nargs="?",
                    help="path to /usr/bin/time -v output (default: stdin)")
    ap.add_argument("--tier", default="unknown",
                    help="tier name to attach (e.g. avx2)")
    ap.add_argument("--cold", action="store_true",
                    help="mark the run as cold-cache")
    args = ap.parse_args()

    if args.logfile is None or str(args.logfile) == "-":
        text = sys.stdin.read()
    else:
        text = args.logfile.read_text()

    rec = parse_time_log(text)
    rec["tier"] = args.tier
    rec["cold"] = bool(args.cold)
    json.dump(rec, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
