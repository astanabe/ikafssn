# Patches

This directory contains patches applied to third-party dependencies before building them.

## ncbi-cxx-toolkit-seqdb-mmap-strategy.patch

Adds `CSeqDB::SetMMapStrategy(EMemMapAdvise)` (and a matching
`CSeqDBAtlas::SetMMapStrategy()`) so callers can advise the kernel on
the expected access pattern for memory-mapped BLAST DB files.  The
patch also defaults `CAtlasMappedFile`'s mmap strategy to
`eMMA_Normal` (OS-decided readahead), and exposes `eMADV_NoReuse`
(`MADV_NOREUSE`) on the underlying `EMemoryAdvise` enum.

ikafssn drives this from `BlastDbReader::set_mmap_strategy()`:
- `kSequential` / `kWillNeed` for `ikafssnindex` metadata + postings
  passes and Stage 3's batched subject-sequence fetch (kafsss
  readahead, pre-fetch into page cache);
- `kDontNeed` between batches so the page cache stays bounded by
  `-memory_limit`;
- `kRandom` for callers performing sparse OID lookups, suppressing
  readahead and reducing page cache pollution.

Submitted upstream as a pull request; once merged, this patch can
be retired.

### How to apply

From the NCBI C++ Toolkit source directory (e.g. `ncbi-cxx-toolkit-public-release-30.2.0/`):

```bash
patch -p1 < /path/to/ikafssn/patches/ncbi-cxx-toolkit-seqdb-mmap-strategy.patch
```

Apply this patch after extracting the toolkit source and before running `./cmake-configure`.
