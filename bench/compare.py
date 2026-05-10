#!/usr/bin/env python3
"""Compare two bench/run_e2e_search.sh JSON outputs and print a markdown table.

Usage:
    python3 bench/compare.py <baseline.json> <candidate.json>

Output columns:
    tier | baseline wall (s) | candidate wall (s) | speedup
        | baseline RSS (MB) | candidate RSS (MB)

speedup is reported as "baseline / candidate" (>1 means candidate is
faster).  Rows are kept in the tier order of the baseline input.  Tiers
present only in the candidate are appended; tiers present only in the
baseline are reported as "-".
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load_runs(path: Path) -> tuple[dict, dict[str, dict]]:
    """Return (top-level metadata, {tier: run-dict})."""
    data = json.loads(path.read_text())
    runs = {r["tier"]: r for r in data.get("runs", [])}
    return data, runs


def fmt_seconds(x: float | None) -> str:
    if x is None or x == 0:
        return "-"
    return f"{x:.2f}"


def fmt_mb(kb: int | None) -> str:
    if kb is None or kb == 0:
        return "-"
    return f"{kb / 1024:.1f}"


def fmt_speedup(a: float | None, b: float | None) -> str:
    if not a or not b:
        return "-"
    return f"{a / b:.2f}×"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("baseline", type=Path)
    ap.add_argument("candidate", type=Path)
    args = ap.parse_args()

    base_meta, base_runs = load_runs(args.baseline)
    cand_meta, cand_runs = load_runs(args.candidate)

    print(f"# Bench compare: {args.baseline.name} vs {args.candidate.name}")
    print()
    print(f"- Baseline: db={base_meta.get('db')}, k={base_meta.get('k')}, "
          f"t={base_meta.get('t')}, template_type={base_meta.get('template_type')}, "
          f"mode={base_meta.get('mode')}, cold={base_meta.get('cold')}, "
          f"git_sha={base_meta.get('git_sha')}")
    print(f"- Candidate: db={cand_meta.get('db')}, k={cand_meta.get('k')}, "
          f"t={cand_meta.get('t')}, template_type={cand_meta.get('template_type')}, "
          f"mode={cand_meta.get('mode')}, cold={cand_meta.get('cold')}, "
          f"git_sha={cand_meta.get('git_sha')}")
    print()

    if (base_meta.get("db") != cand_meta.get("db")
            or base_meta.get("k") != cand_meta.get("k")
            or base_meta.get("t") != cand_meta.get("t")
            or base_meta.get("template_type") != cand_meta.get("template_type")
            or base_meta.get("mode") != cand_meta.get("mode")):
        print("> WARNING: configurations differ between baseline and candidate.")
        print()

    print("| tier | baseline wall (s) | candidate wall (s) | speedup | baseline RSS (MB) | candidate RSS (MB) |")
    print("|---|---:|---:|---:|---:|---:|")
    seen: set[str] = set()
    for tier, base in base_runs.items():
        cand = cand_runs.get(tier)
        seen.add(tier)
        bw = base.get("wall_s")
        cw = cand.get("wall_s") if cand else None
        br = base.get("max_rss_kb")
        cr = cand.get("max_rss_kb") if cand else None
        print(f"| {tier} | {fmt_seconds(bw)} | {fmt_seconds(cw)} | "
              f"{fmt_speedup(bw, cw)} | {fmt_mb(br)} | {fmt_mb(cr)} |")
    for tier, cand in cand_runs.items():
        if tier in seen:
            continue
        cw = cand.get("wall_s")
        cr = cand.get("max_rss_kb")
        print(f"| {tier} | - | {fmt_seconds(cw)} | - | - | {fmt_mb(cr)} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
