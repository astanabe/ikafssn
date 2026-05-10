# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **Important**: Always use Japanese for on-screen communication with the user. However, write files in English unless they are designated as Japanese files. The only Japanese file in this repository is `doc/ikafssn.ja.md`. This file (`CLAUDE.md`) is therefore English-only — including any new sections, examples, and example values inside tables.

## Project Overview

ikafssn (Independent programs of K-mer-based Alignment-Free Similarity Search for Nucleotide sequences) builds a complete inverted index over NCBI BLAST DB nucleotide sequences and performs alignment-free similarity search using k-mer matching and collinear chaining.

- **Primary**: https://github.com/astanabe/ikafssn
- **Secondary**: https://gitlab.com/astanabe/ikafssn

## Build Commands

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make test          # run all unit tests
ctest --test-dir build -R test_kmer_encoding   # run a single test
```

CMake options for selective builds:
- `-DBUILD_HTTPD=OFF` — skip ikafssnhttpd (requires Drogon)
- `-DBUILD_CLIENT=OFF` — skip ikafssnclient (requires libcurl for HTTP mode)
- `-DENABLE_REMOTE_RETRIEVE=OFF` — disable NCBI efetch in ikafssnretrieve

## Command binaries

Each command is a separate executable (not a subcommand) and links only its required dependencies for lightweight deployment:

| Command | Purpose |
|---|---|
| `ikafssnindex` | Build k-mer inverted index from BLAST DB |
| `ikafssnsearch` | Local direct search (mmap index) |
| `ikafssnretrieve` | Extract matched subsequences |
| `ikafssnserver` | Search daemon (UNIX/TCP socket) |
| `ikafssnhttpd` | HTTP REST proxy to ikafssnserver(s) |
| `ikafssnclient` | Client (socket or HTTP) |
| `ikafssninfo` | Index/DB info display (local or remote) |

## Terminology (inverted index terms)

This is a **permanent rule** that applies to every change in this repository, not just the change in front of you right now. When describing the inverted index — in any of: documentation, code identifiers, source comments, docstrings, log messages, error strings, CLI help text, wire-format / protocol field names, or conversation with the user — use the canonical full-text-search vocabulary defined below. There are **no exceptions**: code, tests, fixtures, tools, scripts, build files, and Markdown documents are all bound by the same rule.

### Canonical vocabulary

The canonical mapping (full-text-search domain) is:

| Concept | Canonical term | Definition |
|---|---|---|
| Indexed unit | **term** (interchangeable with **k-mer** in this codebase) | A single token being indexed. In ikafssn the term is a k-mer. |
| List of term IDs | **dictionary** | The list of term IDs. In ikafssn's direct-address layout, the dictionary array is indexed by the term ID (k-mer integer value) and stores the byte offset to that term's posting list inside the posting file. |
| Sequence ID where a term occurs | **posting** | A single entry inside a posting list — the ID of one sequence containing the term (in `.kpx`, also carries position information). |
| List of postings for one term | **posting list** | All postings for a given term, concatenated. |
| Concatenation of all posting lists | **posting file** | The bulk region holding every term's posting list, in dictionary order. In ikafssn this region lives inside `.kix` and `.kpx` after the header + dictionary. |
| Source nucleotide sequence in the BLAST DB | **parent** (or **parent OID**) | The full-length nucleotide sequence as registered in the source BLAST DB volume — i.e. one BLAST OID and the accession that goes with it. Parents are the unit retained for accession lookup, parent length, and user-facing sstart/send coordinates. |
| Indexed slice of a parent | **fragment** | The unit registered as one internal SeqId inside `.kix` / `.kpx` / `.ksx`. A fragment is the slice of a parent produced by the fragment splitter under `-min_length_split` / `-overlap_length`; in the degenerate case (`min_length_split == 0`) every parent has exactly one fragment that spans the whole parent. Adjacent fragments of the same parent share `overlap_length` bases. The k-mer scanner indexes fragments, but the search-side dedup and downstream tools all collapse fragment-relative coordinates back to the parent. |

Auxiliary terms that follow from the canonical vocabulary:
- **dictionary entry** — one offset value in the dictionary array (`offsets[i]`).
- **posting list header** — the fixed-size metadata at the front of a posting list.
- **posting list body** — the variable-size encoded content of a posting list after its header.
- **parent-relative coordinates** — sstart / send expressed as 1-based positions within the **parent**. After Stage 2 dedup, every chain hit is reported in this coordinate space (never in fragment-relative form).
- **fragment-relative coordinates** — sstart / send expressed as 1-based positions within the **fragment**. Internal to the orchestrator: Stage 2 chains are produced in this space and shifted by `(fragment_start - 1)` at the orchestrator → output / response boundary.

A nucleotide sequence is **not** called a "document". Do not import that word from generic IR literature into this codebase. The term "sequence" stands on its own; **parent** and **fragment** are reserved for the index-side parent-OID / fragment distinction introduced in v10.

### Forbidden expressions (always rewrite)

- "offset table" → **dictionary**
- "posting data region" / "posting list region" / "posting blob" / "per-kmer payload" → **posting file**
- "per-kmer header" → **posting list header**
- "posting" used alone to mean the full per-k-mer data → **posting list**
- "payload" in the narrow sense (variable-size content of a posting list) → **posting list body**
- "document" (when the actual referent is a nucleotide sequence) → **sequence**
- "BLAST DB sequence" / "OID-level sequence" / "full-length sequence" (when distinguishing from a fragment) → **parent**
- "subsequence" / "chunk" / "window" / "tile" (when referring to a sliced indexed unit) → **fragment**
- "OID-relative coordinates" / "subject coordinates" (when emphasising parent-vs-fragment) → **parent-relative coordinates**

### Application — no exceptions

This convention is binding on **every** artefact in the repository:

- code identifiers — function names, method names, struct / class members, free variables, namespaces, file names, CMake target names;
- public API symbols and wire-format / protocol field names (rename in lockstep with the appropriate `format_version` / `msg_version` bump if the rename touches the wire format);
- log messages, error strings, and CLI help text;
- comments, docstrings, and any Markdown documentation;
- test names and fixture file names;
- conversation with the user.

When you encounter a legacy identifier or comment that uses a forbidden term, **rename it as part of the change you are already making**. Do not file the rename as a follow-up. Do not retain the legacy term "for compatibility" or "for diff readability" — those are not exceptions. The only legitimate reason to defer a rename is that the canonical term itself is genuinely ambiguous in context; in that case follow the next subsection.

### When the canonical term is ambiguous, ask

If you are about to introduce a new concept or rename an existing identifier and you find that:

- multiple terms in the canonical vocabulary could plausibly apply, or
- the concept does not fit any of the canonical terms above, or
- the established vocabulary in the surrounding code disagrees with the canonical mapping in a way that is not obvious to resolve,

**stop and ask the user before choosing a name.** Do not invent a new convention silently. Do not pick "the closest one" without checking. State the alternatives, the trade-off, and your tentative recommendation, then wait for the user's decision. The choice you make becomes part of this codebase forever, so a few seconds of clarification is always cheaper than a wrong commit.

## Documentation

Usage documentation is maintained in two languages:
- `doc/ikafssn.en.md` — English
- `doc/ikafssn.ja.md` — Japanese

## Test structure

Tests are in `test/` using CTest. Real SSU_eukaryote_rRNA BLAST DB at `db/SSU_eukaryote_rRNA` is used for all BLAST-DB-dependent tests. `test/scripts/setup_ssu_testdata.sh` generates derived test data (ambig DB, multi-volume DBs, queries) in `/tmp/ikafssn_ssu_test/`. Shared fixture: `test/ssu_test_fixture.hpp`.

## Plan Mode Rules

- **Do not enter plan mode autonomously.** Enter plan mode only when the user explicitly instructs you to do so.
- Standard workflow:
  1. Discuss the approach with the user outside of plan mode.
  2. Once the approach is well-defined, the user issues a planning instruction such as:
     - "Enter plan mode and prepare a plan based on the above approach", or
     - After entering plan mode via `/plan`, "Prepare a plan based on the above approach".
  3. Begin planning only after this explicit instruction is received.
- When presenting a plan, display the **absolute path** of the plan file (`.md`).

## Version Management

- The version number is managed via `IKAFSSN_VERSION` in `CMakeLists.txt`. Format: `"0.1.YYYY.MM.DD"` (the date the program source code was last modified).
- Update `IKAFSSN_VERSION` to today's date **only when modifying program source code** (e.g. files under `src/`, `include/`, `test/`). Do **not** update it for documentation-only changes (`doc/`, `README*`, `CLAUDE.md`, etc.) or for changes confined to the build / binary-generation system (CMake files, CI workflows, packaging).
- When a commit does include program source code changes, verify before committing that the date portion of `IKAFSSN_VERSION` matches today's date and update it if outdated.

## Development Environment Rules

- Do not execute commands that require `sudo` directly from Claude Code. Present such commands to the user for manual execution.

## Code Comments and Documentation

- **No change-history or change-rationale narration.** Comments, docstrings, log messages, and Markdown documentation must describe only the *current* implementation. Do not write things like "previously X, now Y", "v10 →v11", "this used to ...", "this was added because of issue #123", or "moved from X.cpp". Explanations of *why the current code is the way it is* are allowed even if they happen to read like change history (e.g. "we use a serial merge here because KsxWriter is thread-unsafe").
- **No "Phase X" labels in code or documentation.** Discussing a plan or its implementation in conversation may use "Phase X" freely, but committed source code, comments, log messages, CLI help text, and Markdown documents must describe each step by what it does (e.g. "metadata pass", "postings pass", "rename .tmp → final"), never by a phase number or letter.

## Git Push Policy

- **Commits may be made at Claude Code's own discretion.** Pushing to a remote, creating a release, and dispatching the GitHub Actions binary-generation workflow are **never** to be performed without an explicit user instruction.
- **Direct push to the `main` branch is permitted** when explicitly instructed. This repository is a personal development repository with no PR-review protection configured. Claude Code's default safety guard (refusing direct push to `main`) does not apply here. When the user instructs "commit and push", execute `git push origin main` without further confirmation.
- **Pre-push checklist.** When the user instructs a push, before pushing perform the following audit on every commit not yet pushed plus every uncommitted change, and commit the result before pushing:
  1. **OS-dependent headers** — confirm no Linux-only / x86-only `#include` was added that would break the macOS build.
  2. **Test updates** — make sure tests covering the changed behaviour exist and pass.
  3. **Documentation updates** — `doc/ikafssn.en.md`, `doc/ikafssn.ja.md`, and any other relevant docs reflect the change.
  4. **Binary-generation system updates** — `.github/workflows/release.yml`, `recipe/` (conda), and the formula in the `homebrew-ikafssn` repository are consistent with the change.
  5. **`IKAFSSN_VERSION`** — bump to today's date if the to-be-pushed changes touch program source code (per the Version Management rule).
  6. **Comment / identifier hygiene** — remove unused functions and identifiers, drop unnecessary or low-value comments, drop any "Phase X" wording, fix any term used with multiple inconsistent meanings, shorten long comments, and remove anything that doesn't belong (per the Code Comments and Documentation rule).
- **Plan-implementation completeness check.** When the user instructs a push that follows the implementation of a plan, also verify — without being asked — that every step of the plan has actually been implemented; if any step is missing, complete and commit it before pushing.

## Release Creation

When the user asks for "リリース作成" / "release creation" (or equivalents like "リリースして" / "リリースを作って"), follow this fixed sequence without further confirmation:

1. **Sync local main with remote.** Commit any pending source changes (verify `IKAFSSN_VERSION` matches today's date for source-touching commits per the Version Management rule) and `git push origin main`.
2. **Replace today's release if it already exists.** Look up the tag `v${IKAFSSN_VERSION}` (today's date version). If a release with that tag is already published — even if it was just created today — delete both the GitHub release and its underlying git tag (`gh release delete <tag> --cleanup-tag --yes`) before recreating it. This guarantees the new release points at the latest `main` commit.
3. **Create the release on the latest `main`.** Run `gh release create v${IKAFSSN_VERSION} --target main --title "ikafssn v${IKAFSSN_VERSION}" --notes "<short summary>"` (write a short release-notes summary derived from the recent commit log). The CI workflow expects this release tag to exist before it begins uploading assets.
4. **Trigger the binary-generation workflow.** Dispatch the `Release` workflow with the same version: `gh workflow run release.yml -f version=${IKAFSSN_VERSION}`. Do **not** run `cpack` locally — release binaries are built and uploaded only by this workflow (already noted in user memory).
5. **Report the release URL and the run URL** back to the user (`gh release view v${IKAFSSN_VERSION} --json url -q .url` and `gh run list --workflow=release.yml --limit 1 --json url -q '.[0].url'`) so they can monitor progress.
