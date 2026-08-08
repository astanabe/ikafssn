# Patches

This directory contains patches applied to third-party dependencies before building them.

## parasail-degmatch-cigar-score.patch

Adds the `degmatch` nucleotide score matrix (`parasail/matrices/degmatch.h`,
registered in `parasail/matrix_lookup.h`), which scores a pair of degenerate
bases positively when the two IUPAC codes share at least one possible
nucleotide.  ikafssn uses it as the default Stage 3 matrix
(`-stage3_score_matrix degmatch`) so primer-style queries and subjects
carrying degenerate bases align sensibly.

The patch also fixes `src/cigar_template.c` so the score reported alongside a
traceback matches the score of the emitted CIGAR.

### How to apply

From the Parasail source directory (e.g. `parasail-2.6.2/`):

```bash
patch -p1 < /path/to/ikafssn/patches/parasail-degmatch-cigar-score.patch
```

Apply this patch after extracting the Parasail source and before running `cmake`.
