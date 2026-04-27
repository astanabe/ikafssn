#!/usr/bin/env python3
"""Phase 5a — analyze v3 posting TSV(s) and estimate v4 (SIMD-FastPFOR* +
Simple-8b) compressed size.

Reads one or more TSV files produced by tools/dump_postings.cpp (one row per
non-empty k-mer, columns: kmer / count / seq_id_csv / pos_csv) and emits:

  * posting count distribution (1, 2-7, 8-127, 128-1023, 1024-)
  * delta distribution for .kix (seq_id deltas, including first abs)
  * delta distribution for .kpx within-seq (and abs position distribution)
  * estimated v4 byte size with FastPFOR*-style block 128 + Simple-8b tail
  * v3 vs v4 size comparison (uses --kix-bytes / --kpx-bytes for actuals)

The FastPFOR cost model is intentionally simple but matches the published
behaviour (Lemire & Boytsov 2015) closely enough for go/no-go judgment:

  - Per 128-element block:
      header:   1 byte b (bit width 0..32) + 1 byte except_count + 2 bytes pad
      body:     128 * b bits
      except :  except_count * (1 byte index + 4 bytes high-32 patch)
    where b is chosen as the 90th percentile of the block's value-bit-widths
    (a common heuristic).  Values exceeding 1<<b are coded as exceptions.

  - Tail (< 128 elements) goes through Simple-8b: each 8-byte word stores
    a selector + payload encoding 1..240 small ints.  We approximate the
    cost as ceil(N * mean_bits / 56) bytes (Simple-8b best case = 60 1-bit
    ints per word, worst case = 1 60-bit int per word; 56 = effective
    payload bits accounting for selector).

For .kix: deltas are seq_id deltas (first element is abs).
For .kpx: per the plan we use FOR-within-block on absolute positions; we
estimate b as ceil(log2(max(pos)-min(pos))) per block.
"""

import argparse
import gzip
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

PFOR_BLOCK = 128
PFOR_HEADER_BYTES = 4   # b (1) + except_count (1) + reserved (2)
KPX_FOR_HEADER_BYTES = 8  # min_position (4) + b (1) + except_count (1) + reserved (2)
EXCEPT_BYTES = 5  # 1 (index) + 4 (high-32 patch)


def open_text(path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "r")


def bit_width(v):
    if v == 0:
        return 0
    return v.bit_length()


def estimate_pfor_block(values, for_within=False):
    """Estimate compressed bytes for a 128-element block.

    values: list of 128 uint32.
    for_within: if True, subtract block min before computing widths
                (FOR-within-block, used for .kpx absolute positions).

    Heuristic: pick b = 90th percentile bit width.  Values >= (1<<b) become
    exceptions (high bits stored separately, low b bits in body).
    """
    assert len(values) == PFOR_BLOCK
    if for_within:
        mn = min(values)
        vals = [v - mn for v in values]
    else:
        vals = values
        mn = 0
    widths = sorted(bit_width(v) for v in vals)
    # 90th percentile (FastPFOR's default tradeoff)
    p90 = widths[int(0.90 * (PFOR_BLOCK - 1))]
    b = max(0, min(32, p90))
    threshold = 1 << b if b < 32 else (1 << 32)
    except_count = sum(1 for v in vals if v >= threshold)
    if except_count > 127:
        # Promote b until except_count fits
        while except_count > 127 and b < 32:
            b += 1
            threshold = 1 << b if b < 32 else (1 << 32)
            except_count = sum(1 for v in vals if v >= threshold)
    body_bits = PFOR_BLOCK * b
    body_bytes = (body_bits + 7) // 8
    hdr = KPX_FOR_HEADER_BYTES if for_within else PFOR_HEADER_BYTES
    return hdr + body_bytes + except_count * EXCEPT_BYTES, b, except_count


def estimate_simple8b(values):
    """Estimate Simple-8b bytes for a tail of < 128 ints.

    Real Simple-8b has 16 selectors; we use a coarse approximation:
    if all values are < 1<<bits, ceil(N / floor(60/bits)) words.
    """
    if not values:
        return 0
    n = len(values)
    max_v = max(values)
    bits = max(1, bit_width(max_v))
    # Simple-8b selectors: 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 60
    # Pick smallest selector that fits.
    for bw in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 60):
        if bw >= bits:
            per_word = 60 // bw
            words = math.ceil(n / per_word)
            return words * 8
    # 60-bit fallback
    return math.ceil(n / 1) * 8


def estimate_v4_kix(seq_ids):
    """Compute v4 bytes for a single k-mer's seq_id posting (delta encoded).
    Layout: 4 (count) + N x PFor block + Simple-8b tail.
    """
    n = len(seq_ids)
    if n == 0:
        return 0
    # Delta-encode (first is absolute, then differences).
    deltas = [seq_ids[0]]
    for i in range(1, n):
        deltas.append(seq_ids[i] - seq_ids[i - 1])
    body = 4  # count u32
    n_blocks = n // PFOR_BLOCK
    tail = n - n_blocks * PFOR_BLOCK
    for b in range(n_blocks):
        chunk = deltas[b * PFOR_BLOCK:(b + 1) * PFOR_BLOCK]
        sz, _, _ = estimate_pfor_block(chunk, for_within=False)
        body += sz
    body += estimate_simple8b(deltas[n_blocks * PFOR_BLOCK:])
    return body


def estimate_v4_kpx(positions):
    """Compute v4 bytes for a single k-mer's position posting (FOR-within).
    Layout: 4 (count) + N x FOR-PFor block + Simple-8b tail.
    """
    n = len(positions)
    if n == 0:
        return 0
    body = 4  # count u32
    n_blocks = n // PFOR_BLOCK
    for b in range(n_blocks):
        chunk = positions[b * PFOR_BLOCK:(b + 1) * PFOR_BLOCK]
        sz, _, _ = estimate_pfor_block(chunk, for_within=True)
        body += sz
    # Tail: Simple-8b on absolute positions (or relative to first?).
    # Use relative-to-first to keep widths small.
    tail = positions[n_blocks * PFOR_BLOCK:]
    if tail:
        mn = min(tail)
        rel = [v - mn for v in tail]
        body += 4 + estimate_simple8b(rel)  # +4 for tail base
    return body


def varint_size(v):
    if v == 0:
        return 1
    n = 0
    while v > 0:
        n += 1
        v >>= 7
    return n


def estimate_v3_kix(seq_ids):
    """v3 bytes (LEB128 of seq_id deltas, no count header)."""
    if not seq_ids:
        return 0
    total = varint_size(seq_ids[0])
    for i in range(1, len(seq_ids)):
        total += varint_size(seq_ids[i] - seq_ids[i - 1])
    return total


def estimate_v3_kpx(seq_ids, positions):
    """v3 bytes (LEB128 of position delta, reset on new seq)."""
    if not positions:
        return 0
    total = 0
    prev_sid = seq_ids[0]
    prev_pos = positions[0]
    total += varint_size(positions[0])
    for i in range(1, len(positions)):
        if seq_ids[i] != prev_sid:
            total += varint_size(positions[i])
        else:
            total += varint_size(positions[i] - prev_pos)
        prev_sid = seq_ids[i]
        prev_pos = positions[i]
    return total


def count_bucket(c):
    if c == 1:
        return "1"
    if c <= 7:
        return "2-7"
    if c <= 127:
        return "8-127"
    if c <= 1023:
        return "128-1023"
    if c <= 16383:
        return "1024-16383"
    return "16384-"


def percentile(sorted_list, p):
    if not sorted_list:
        return 0
    idx = int(p * (len(sorted_list) - 1))
    return sorted_list[idx]


def analyze_one(tsv_path):
    """Return a dict of per-file statistics."""
    count_hist = Counter()
    kix_v3 = 0
    kix_v4 = 0
    kpx_v3 = 0
    kpx_v4 = 0
    n_kmers = 0
    n_postings = 0
    sid_deltas_sample = []
    pos_within_sample = []
    pos_abs_sample = []
    block_b_dist = Counter()  # PFOR b distribution (.kix)
    block_b_dist_kpx = Counter()  # FOR b distribution (.kpx)

    SAMPLE_LIMIT = 1_000_000  # cap to keep memory bounded

    with open_text(tsv_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            kmer = int(parts[0])
            cnt = int(parts[1])
            sid_csv = parts[2]
            pos_csv = parts[3] if len(parts) >= 4 and parts[3] else ""
            seq_ids = [int(x) for x in sid_csv.split(",")] if sid_csv else []
            positions = [int(x) for x in pos_csv.split(",")] if pos_csv else []
            n_kmers += 1
            n_postings += cnt
            count_hist[count_bucket(cnt)] += 1

            kix_v3 += estimate_v3_kix(seq_ids)
            kix_v4 += estimate_v4_kix(seq_ids)
            if positions:
                kpx_v3 += estimate_v3_kpx(seq_ids, positions)
                kpx_v4 += estimate_v4_kpx(positions)

            # Sample widths for distribution analysis.
            if len(sid_deltas_sample) < SAMPLE_LIMIT:
                if seq_ids:
                    sid_deltas_sample.append(seq_ids[0])
                    for i in range(1, len(seq_ids)):
                        if len(sid_deltas_sample) >= SAMPLE_LIMIT:
                            break
                        sid_deltas_sample.append(seq_ids[i] - seq_ids[i - 1])

            if positions and len(pos_within_sample) < SAMPLE_LIMIT:
                prev_sid = seq_ids[0]
                prev_pos = positions[0]
                pos_abs_sample.append(positions[0])
                for i in range(1, len(positions)):
                    if len(pos_within_sample) >= SAMPLE_LIMIT:
                        break
                    if seq_ids[i] == prev_sid:
                        pos_within_sample.append(positions[i] - prev_pos)
                    pos_abs_sample.append(positions[i])
                    prev_sid = seq_ids[i]
                    prev_pos = positions[i]

            # Block-level b distribution (sample first 8 blocks per posting).
            n_blocks = min(8, len(seq_ids) // PFOR_BLOCK)
            for b in range(n_blocks):
                chunk = seq_ids[b * PFOR_BLOCK:(b + 1) * PFOR_BLOCK]
                # convert to deltas first
                d = [chunk[0] if b == 0 else chunk[0] - seq_ids[b * PFOR_BLOCK - 1]]
                for i in range(1, PFOR_BLOCK):
                    d.append(chunk[i] - chunk[i - 1])
                _, bw, _ = estimate_pfor_block(d, for_within=False)
                block_b_dist[bw] += 1

            n_blocks_kpx = min(8, len(positions) // PFOR_BLOCK)
            for b in range(n_blocks_kpx):
                chunk = positions[b * PFOR_BLOCK:(b + 1) * PFOR_BLOCK]
                _, bw, _ = estimate_pfor_block(chunk, for_within=True)
                block_b_dist_kpx[bw] += 1

    sid_deltas_sample.sort()
    pos_within_sample.sort()

    def stats(lst):
        if not lst:
            return {"n": 0, "mean": 0, "p50": 0, "p99": 0, "p999": 0, "max": 0}
        return {
            "n": len(lst),
            "mean": sum(lst) / len(lst),
            "p50": percentile(lst, 0.5),
            "p99": percentile(lst, 0.99),
            "p999": percentile(lst, 0.999),
            "max": lst[-1],
        }

    return {
        "tsv": str(tsv_path),
        "n_kmers": n_kmers,
        "n_postings": n_postings,
        "count_hist": dict(count_hist),
        "kix_v3_bytes_estimated": kix_v3,
        "kix_v4_bytes_estimated": kix_v4,
        "kpx_v3_bytes_estimated": kpx_v3,
        "kpx_v4_bytes_estimated": kpx_v4,
        "sid_delta_stats": stats(sid_deltas_sample),
        "pos_within_seq_delta_stats": stats(pos_within_sample),
        "block_b_kix_hist": dict(block_b_dist),
        "block_b_kpx_for_hist": dict(block_b_dist_kpx),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", nargs="+", help="TSV file(s) from dump_postings")
    ap.add_argument("--output", default=None, help="markdown output path")
    ap.add_argument("--json-output", default=None, help="raw json output path")
    ap.add_argument("--kix-bytes", type=int, default=0,
                    help="actual on-disk total .kix posting bytes (sum of "
                         "data sections) for v3 baseline reference")
    ap.add_argument("--kpx-bytes", type=int, default=0,
                    help="actual on-disk total .kpx posting bytes for v3 baseline")
    args = ap.parse_args()

    results = [analyze_one(p) for p in args.tsv]

    total_kix_v3_est = sum(r["kix_v3_bytes_estimated"] for r in results)
    total_kix_v4_est = sum(r["kix_v4_bytes_estimated"] for r in results)
    total_kpx_v3_est = sum(r["kpx_v3_bytes_estimated"] for r in results)
    total_kpx_v4_est = sum(r["kpx_v4_bytes_estimated"] for r in results)
    total_v3_est = total_kix_v3_est + total_kpx_v3_est
    total_v4_est = total_kix_v4_est + total_kpx_v4_est
    size_ratio_est = total_v4_est / total_v3_est if total_v3_est > 0 else 0.0

    actual_v3 = args.kix_bytes + args.kpx_bytes
    size_ratio_actual = (
        total_v4_est / actual_v3 if actual_v3 > 0 else 0.0
    )

    md = []
    md.append("# Phase 5a — `db/tsa_nt` posting analysis\n")
    md.append("Per-volume estimates (FastPFOR* block=128 + Simple-8b tail).\n")
    for r in results:
        md.append(f"## {r['tsv']}\n")
        md.append(f"- k-mers (non-empty): {r['n_kmers']:,}")
        md.append(f"- postings total: {r['n_postings']:,}")
        md.append("- count distribution: " + ", ".join(
            f"{k}={v:,}" for k, v in sorted(r['count_hist'].items())))
        md.append(f"- estimated v3 .kix posting bytes: {r['kix_v3_bytes_estimated']:,}")
        md.append(f"- estimated v4 .kix posting bytes: {r['kix_v4_bytes_estimated']:,}")
        md.append(f"- estimated v3 .kpx posting bytes: {r['kpx_v3_bytes_estimated']:,}")
        md.append(f"- estimated v4 .kpx posting bytes: {r['kpx_v4_bytes_estimated']:,}")
        s = r['sid_delta_stats']
        md.append(f"- seq_id delta stats: mean={s['mean']:.1f} p50={s['p50']} "
                  f"p99={s['p99']} p999={s['p999']} max={s['max']}")
        s = r['pos_within_seq_delta_stats']
        md.append(f"- pos within-seq delta stats: mean={s['mean']:.1f} p50={s['p50']} "
                  f"p99={s['p99']} p999={s['p999']} max={s['max']}")
        md.append(f"- .kix block-b histogram: {r['block_b_kix_hist']}")
        md.append(f"- .kpx FOR block-b histogram: {r['block_b_kpx_for_hist']}\n")

    md.append("## Aggregate (all volumes)\n")
    md.append(f"- Total estimated v3 posting bytes (.kix + .kpx): {total_v3_est:,}")
    md.append(f"- Total estimated v4 posting bytes (.kix + .kpx): {total_v4_est:,}")
    md.append(f"- estimated size ratio (v4 / v3): {size_ratio_est:.4f}")
    if actual_v3 > 0:
        md.append(f"- actual on-disk v3 posting bytes (.kix + .kpx): {actual_v3:,}")
        md.append(f"- estimated v4 vs actual v3 ratio: {size_ratio_actual:.4f}")

    md.append("\n### Phase 5a go/no-go gate\n")
    md.append("- Gate: (size_ratio <= 0.90) **OR** (speed_ratio <= 0.90)")
    if size_ratio_est <= 0.90:
        md.append(f"- **PASS (size)**: estimated size_ratio = {size_ratio_est:.4f} <= 0.90")
    else:
        md.append(f"- size gate not met: estimated size_ratio = {size_ratio_est:.4f}")
        md.append("  (need to also check bench_pfd_proto speed_ratio)")

    md_text = "\n".join(md) + "\n"

    if args.output:
        Path(args.output).write_text(md_text)
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(md_text)

    if args.json_output:
        Path(args.json_output).write_text(json.dumps({
            "per_file": results,
            "total_kix_v3_est": total_kix_v3_est,
            "total_kix_v4_est": total_kix_v4_est,
            "total_kpx_v3_est": total_kpx_v3_est,
            "total_kpx_v4_est": total_kpx_v4_est,
            "size_ratio_estimated": size_ratio_est,
            "actual_v3_bytes": actual_v3,
            "size_ratio_vs_actual_v3": size_ratio_actual,
        }, indent=2))
        print(f"wrote {args.json_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
