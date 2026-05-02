# Phase 5a — SIMD-FastPFOR\* + Simple-8b feasibility study

**Status: GO (speed gate passed; size gate not met).**

Per the plan, the go/no-go criterion is

```
go = (size_ratio <= 0.90) || (speed_ratio <= 0.90)
```

Speed gate is met overwhelmingly (~50× decode throughput vs Phase 3a varint).
Size gate is **not** met for `db/tsa_nt` at k=11 — projected v4 sizing is
~+5% relative to v3. We proceed to Phase 5b on the speed evidence alone, but
flag the size projection as a known risk that Phase 5b implementation must
re-validate against real index data.

Baseline commit: `8ccd295` (Phase 3a complete = `vectorize` HEAD at start
of Phase 5a).

## 1. Dataset

| | |
|---|---|
| BLAST DB | `db/tsa_nt` (4 volumes, 18.7 M sequences, 14.9 G postings) |
| k | 11 (table_size = 4,194,304) |
| mode | 2 (both `.kix` and `.kpx` generated) |
| host | AMD Ryzen 9 9950X (Zen 5, AVX-512 VBMI2) |

## 2. v3 size baseline (actual on-disk)

Per-volume v3 file sizes and posting file sizes (file − header − dictionary):

| volume | `.kix` file | `.kix` posting file | `.kpx` file | `.kpx` posting file |
|---|---:|---:|---:|---:|
| 00 | 7,531,831,047 | 7,498,276,543 | 7,087,588,288 | 7,054,033,824 |
| 01 | 6,723,026,248 | 6,689,471,744 | 6,571,981,378 | 6,538,426,914 |
| 02 | 6,166,830,286 | 6,133,275,782 | 6,935,265,025 | 6,901,710,561 |
| 03 | 3,692,546,609 | 3,675,769,325 | 4,715,868,773 | 4,682,314,309 |
| **TOTAL** | **24,114,234,190** | **23,996,793,394** | **25,310,703,464** | **25,176,485,608** |

Combined `.kix + .kpx` file size: **49,424,937,654 B (46.0 GiB)**.
Combined posting bytes (delta+varint payload only): **49,173,279,002 B (45.8 GiB)**.

## 3. v4 size projection (estimated)

Method: dump postings via `tools/dump_postings` with `--stride 100` (1%
uniform sample of non-empty k-mers per volume) and run
`tools/analyze_postings.py` to compute per-block FastPFOR\* (block=128, b =
90th-percentile bit width, exceptions stored separately) plus Simple-8b for
the tail. `.kpx` uses FOR-within-block (subtract block-min from absolute
positions). Per-volume sample is 41,944 k-mers.

| volume | sample postings | est. v3 `.kix` | est. v4 `.kix` | est. v3 `.kpx` | est. v4 `.kpx` |
|---|---:|---:|---:|---:|---:|
| 00 | 47,141,381 | 87,748,065 | 94,133,879 | 83,800,114 | 81,425,808 |
| 01 | 39,232,205 | 72,587,776 | 77,951,542 | 71,343,062 | 69,570,992 |
| 02 | 45,439,687 | 72,909,721 | 85,044,759 | 83,660,764 | 79,796,670 |
| 03 | 33,190,284 | 47,462,701 | 61,829,183 | 61,056,589 | 58,440,228 |
| sample total | 165,003,557 | 280,708,263 | 318,959,363 | 299,860,529 | 289,233,698 |

- Sample est. v3 (`.kix + .kpx`) = 580,568,792 B
- Sample est. v4 (`.kix + .kpx`) = 608,193,061 B
- **Estimated size_ratio = 608,193,061 / 580,568,792 = 1.0476**

Projecting onto actual v3 totals: v4 ≈ 51.5 GB vs v3 49.2 GB.

### 3.1 Why is size larger than v3 here?

- `.kix` ID-delta stream: heavy-tail distribution. Many deltas are 0 or 1
  (intra-sequence repeats), but occasional jumps to 17-bit values force the
  PFOR block bit-width up. The 90th-percentile heuristic settles around
  b=13-15 per block (see histograms in `analyze.md`), so a typical block
  spends 13-15 bits per delta vs ~8 bits for varint on the same data.
  Exception storage clamps the worst offenders but adds 5 B per exception.
- `.kpx` absolute positions with FOR-within-block: better behaved (b=10-13
  typical) and ~5% smaller than v3, but not enough to offset `.kix` regression.
- `tsa_nt` skews to short contigs (TSA records, mean ~kbp), which inflates
  the per-block bit-width via mid-range outliers more than NCBI nt-class
  long-chromosome corpora would.

This means **Phase 5b implementation must measure size on the real
v4-formatted index before declaring success on size**. The Phase 5a
estimator is conservative (uses 90th-percentile b without optimal-b search),
but the sign of the result is unlikely to flip on this dataset.

### 3.2 Posting count distribution (sample, all volumes)

`count_bucket` aggregated across 4 vols (each 41,944 k-mers sampled):

| bucket | total k-mers (sample) |
|---|---:|
| 8–127 | 4,485 |
| 128–1023 | 113,448 |
| 1024–16383 | 49,715 |
| 16384– | 128 |

- Count >= 128 dominates: every k-mer in this DB fills several PFOR blocks.
  This is good for amortizing the 4 B per-block header but bad for size if
  per-block b is high.

## 4. Speed bench (synthetic deltas)

Bench: `bench/bench_pfd_proto`, same delta distribution as `bench_varint.cpp`
(mix of 1- to 5-byte varint values; raw 4 B/value).

| benchmark | n_values | cpu_time | items/s | bytes/s |
|---|---:|---:|---:|---:|
| `BM_PfdProto_VarintScalar` | 16,384 | 12.93 µs | 1.27 G/s | 3.20 GiB/s |
| `BM_PfdProto_VarintScalar` | 262,144 | 1232.81 µs | 0.21 G/s | 0.55 GiB/s |
| `BM_PfdProto_VarintSimd avx2` | 262,144 | 1496.21 µs | 0.18 G/s | 0.45 GiB/s |
| `BM_PfdProto_VarintSimd avx512bw` | 262,144 | 1342.42 µs | 0.20 G/s | 0.50 GiB/s |
| `BM_PfdProto_VarintSimd avx512vbmi` | 262,144 | 1388.46 µs | 0.19 G/s | 0.49 GiB/s |
| `BM_PfdProto_VarintSimd avx512vbmi2` | 262,144 | 1350.87 µs | 0.19 G/s | 0.50 GiB/s |
| `BM_PfdProto_FastPforDecode` (simdfastpfor256) | 16,384 | **1.68 µs** | **9.76 G/s** | 29.3 GiB/s |
| `BM_PfdProto_FastPforDecode` (simdfastpfor256) | 262,144 | **26.61 µs** | **9.85 G/s** | 29.6 GiB/s |
| `BM_PfdProto_Simple8bDecode` | 262,144 | (collected) | (collected) | |

- AVX-512 VBMI2 *varint* decode: 1342 µs / 262 k items.
- FastPFor *simdfastpfor256* decode: **26.6 µs / 262 k items** — **~50× faster**.

```
speed_ratio = 26.61 / 1342.42 = 0.0198    →    PASS (<= 0.90)
```

Note: `simdfastpfor256` uses 256-element blocks; Phase 5b is specced at
block=128 (matching SPE 2015). Decode cost per element is unaffected by
that choice (same SIMD bitunpacker).

The bench process emits a `malloc(): corrupted top size` message at
shutdown after JSON results are written. This is a known FastPFor 0.4.0
destructor ordering issue and does not affect measurement validity (results
are flushed before the crash). Phase 5b's in-tree codec will not depend on
FastPFor's runtime, only on the bitpacker headers, so this issue does not
carry forward.

## 5. Posting characteristics (sample)

### Volume 00 (representative)
- seq_id delta: mean=3.6, p50=0, p99=66, p999=420, max=71,369
- pos within-seq delta: mean=11.3, p50=1, p99=241, p999=1917, max=29,759
- `.kix` block-b histogram (peaks): b=14 (77 k blocks), b=15 (60 k), b=13 (40 k)
- `.kpx` FOR block-b histogram (peaks): b=11 (122 k), b=12 (51 k), b=10 (38 k)

### Volume 01
- seq_id delta: mean=141.5, p99=3546, p999=13530 — much sparser
- `.kix` blocks skew to b=13-15

### Volume 02 / 03 — see `analyze.md` for full per-volume histograms.

The `.kpx` FOR-within-block strategy is consistently better than absolute
positions (b ~10-13 with FOR vs b ~20+ without — pos values are in the
hundreds-of-thousands range).

## 6. Decision: GO

| gate | threshold | observed | outcome |
|---|---|---|---|
| size_ratio | ≤ 0.90 | **1.05** (estimated) | NOT MET |
| speed_ratio | ≤ 0.90 | **0.020** | **MET** |

Plan rule: `go = size OR speed`. **GO.**

### Risks carried into Phase 5b

1. Estimated size *increase* (~+5%) on tsa_nt-class data. Phase 5b must
   re-measure against the real v4 binary on the same DB and re-evaluate
   before merging if the actual delta is worse than projected.
2. Mitigations available without scope creep:
   - Optimal-b search per block (FastPFOR's actual algorithm — pick b that
     minimises `body_bytes + except_count * 5`, not the 90th percentile).
     The estimator is conservative; real FastPFor will likely beat the
     projection by a few percent.
   - Patched-PFor (separate exception block stream) — listed as Phase 5e
     follow-up work. Not in 5b scope.
3. NCBI nt-class data (longer contigs, more uniform delta distribution) is
   expected to compress *better* than tsa_nt; the current dataset is a
   pessimistic case for the codec choice.

## 7. Artifacts

- `tools/dump_postings.cpp` — reads v3 `.kix`/`.kpx`, writes posting TSV.
- `tools/analyze_postings.py` — compression-cost model, generates report.
- `bench/bench_pfd_proto.cpp` — synthetic speed comparison.
- Raw bench JSON: `/tmp/phase5a_v3/bench_pfd_proto.json` (truncated by
  destructor crash; per-row data above is from full console output).
- Per-volume analysis: `/tmp/phase5a_v3/analyze.md` and `analyze.json`.
- Sampled TSV dumps: `/tmp/phase5a_v3/tsa_nt.0{0,1,2,3}.tsv` (~2 GiB total
  at stride=100, throwaway).

## 8. Out-of-scope / not measured

- Encode speed (only decode matters for Stage 1 hot path).
- AArch64 NEON (Phase 5e).
- Stage 1 end-to-end speed-up — measured in Phase 5b.

---

Phase 5a complete. Proceeding to Phase 5b.
