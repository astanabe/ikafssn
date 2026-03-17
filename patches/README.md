# Patches

This directory contains patches applied to third-party dependencies before building them.

## ncbi-cxx-toolkit-seqdb-madvise-random.patch

Adds `MemMapAdvise(eMMA_Random)` to `CAtlasMappedFile`'s constructor in the NCBI C++ Toolkit's SeqDB reader. This disables the kernel's default readahead (typically 128 KB per page fault) on BLAST DB memory-mapped files, reducing page cache pollution by ~32x during Stage 3 alignment. Without this patch, readahead on BLAST DB accesses evicts ikafssn's index data from the page cache, causing slow Stage 1 startup on subsequent queries.

### How to apply

From the NCBI C++ Toolkit source directory (e.g. `ncbi-cxx-toolkit-public-release-30.1.0/`):

```bash
patch -p1 < /path/to/ikafssn/patches/ncbi-cxx-toolkit-seqdb-madvise-random.patch
```

Apply this patch after extracting the toolkit source and before running `./cmake-configure`.
