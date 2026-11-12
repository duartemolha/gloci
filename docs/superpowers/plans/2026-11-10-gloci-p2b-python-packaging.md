# gloci P2b Python Library, Packaging and Paper Evidence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Order [owner, 2026-11-10]:** Task 1 runs right after P2a. Every other task runs after P3a and P3b,
so the Python library, docs and paper evidence cover the complete command set (spec §10, "Order").

**Goal:** Build the distribution half of spec Phase 2:
- the Python library (§8.5);
- crates.io, PyPI and bioconda packaging, with versioning and a changelog;
- the mdBook site;
- the evidence suite for a Bioinformatics Application Note: real-data benchmarks, real-data concordance, one command that regenerates every table, and Zenodo artefacts. The paper's text is not part of this plan.

**Architecture:**
- **Python reuses the CLI code path in-process.** A new `gloci::embed` module maps kwargs to an argv, parses it with the same clap `Cli`, and runs the same `commands::*` code. Inputs and outputs are in-memory endpoints registered in gloci-io, and chained steps are threads joined by bounded pipes. Results therefore equal the CLI's by construction (§8.5).
- **`gloci-py`** (PyO3 + maturin) is glue only (§8.1).
- **The book** includes the existing docs, never copies them.
- **The paper suite** extends xtask's bench and diff runners.
- No ops semantics change.

**Tech Stack:** Everything P2a uses, plus:
- PyO3 (`abi3-py310`, `extension-module`) and maturin 1.x;
- arrow-rs with `ffi`;
- pytest and mypy;
- mdBook;
- conda-build (pinned `condaforge/miniforge3` image) and `cargo-bundle-licenses`;
- GNU `time`, `gawk` and `ucsc-bigwigtobedgraph` in the bench image;
- the `tar` crate.

Use the newest releases at implementation time. PyO3 0.26 renamed `allow_threads` to `detach`.

**Spec:** `docs/superpowers/specs/2026-11-06-gloci-design.md`, covering:
- §2 (goal 6, non-goal Windows);
- §6.4;
- §8.1 and §8.5;
- §9.1–§9.3;
- §10 (P2 row, Distribution, Docs);
- §11;
- §13.

**Builds on:** P2a, fully merged, including Task 20 (benchmarks) and the final command set.
- D1–D78 stay; D80 amends D59.
- D83 closes the P1 follow-up "no CLI-only behaviour before the Python library".

**Out of scope:**
- P3;
- the paper text;
- free-threaded CPython (abi3 can't target it);
- native Windows;
- macOS x86_64 wheels (the spec names arm64; the sdist covers the rest);
- conda-forge;
- stage-2 worker pools (D77).

**P2b is done when:**
- wheels import and pass the smoke test on glibc 2.17 (x86_64, aarch64) and macOS arm64;
- every golden case that maps to Python gives CLI-identical bytes;
- a tag builds and smoke-tests every artefact and attaches it to a GitHub Release; publishing to PyPI, crates.io, GHCR and bioconda is wired but off until the owner enables each one (D80);
- both bioconda recipes pass under strict channel priority;
- the book deploys on tags;
- one full real-data run is committed, with every concordance mismatch anchored, and `cargo xtask paper` reproduces the committed tables.

## Global Constraints

P2a's Global Constraints all apply, including:
- u64 half-open coordinates;
- `rust-version = "1.91"`;
- no htslib, and one static musl binary;
- exit codes 0/1/2;
- errors that say what, where and how to fix;
- byte-identical output at every `--threads`;
- atomic output;
- record streams (§8.1);
- every difference anchored and pinned by a golden;
- CI on Linux x86_64/aarch64 and macOS aarch64.

New in P2b, from the spec:
- **§8.5:**
  - maturin abi3 wheels for manylinux x86_64/aarch64 and macOS arm64, with no runtime dependencies beyond Python;
  - in-process, releasing the GIL during work;
  - lazy pipelines whose kwargs are the CLI flag names with `-`→`_`;
  - Arrow hand-off with zero copy to pyarrow and polars, and pandas via pyarrow; columns `chrom start end` first;
  - exceptions `GlociError`, `UsageError`(`ValueError`), `DataError(.path, .line)` and `UnsortedError(DataError)`;
  - semver;
  - Python tests that reuse the golden files.
- **§10:** musl and macOS binaries, `cargo install`, a container, bioconda and PyPI; an mdBook with semantics, migration and differences.
- **§6.4:** completions and man pages come from the CLI definition.
- **§9.2:**
  - bedtools, bedops and vcftools are the binding baselines; other tools are yardsticks;
  - every tool is pinned (bedtools 2.31.1, bedops 2.4.41, vcftools 0.1.16, htslib 1.21, bcftools 1.21);
  - results are kept per release.
- **CI:**
  - every action is pinned to a 40-hex SHA;
  - CI never downloads paper datasets.

## Decisions this plan adds

| # | Decision |
|---|---|
| D79 | **The licence is MIT** (§13 left it open). bedtools is MIT, so the paper's readers know it, and every dependency is compatible. **[owner confirmed MIT in principle, 2026-11-10.]** Task 1 writes `Copyright (c) 2026 Duarte Molha`; if the owner's institution must hold the copyright, the line changes before the first public release (nothing is published before then, D80). Publish jobs fail closed without `LICENSE`. |
| D80 | **Publishing (amends D59). [revised 2026-11-10: owner: "i do not want this published in pypi just yet"]** <ul><li>A `v*` tag builds and smoke-tests the archives (with man pages and completions), wheels and sdist, and attaches them all, with `SHA256SUMS`, to a GitHub Release. That is the only thing a tag publishes by default.</li><li>The publish jobs for PyPI, crates.io and GHCR exist but each runs only when its repository variable is `true`: `GLOCI_PUBLISH_PYPI`, `GLOCI_PUBLISH_CRATES` and `GLOCI_PUBLISH_GHCR`. All are unset by default, so nothing leaves the repository until the owner turns one on.</li><li>Each publish step is idempotent, so re-running a failed release finishes it.</li><li>`workflow_dispatch` with `dry_run` publishes nothing, not even the GitHub Release.</li><li>Bioconda is not submitted in P2b. The recipes are built and tested locally and in the nightly job, and the first bioconda PR waits for the owner.</li></ul> |
| D81 | **One version** across the crates, the Python package, the conda packages, the tag and `--version`. <ul><li>Semver; before 1.0, a breaking change bumps the minor.</li><li>The first release is **0.2.0**.</li><li>`CHANGELOG.md` uses the Keep a Changelog format, and the release notes are the version's section.</li></ul> |
| D82 | **crates.io:** publish core, io, engine, vcf, ops and `gloci`, so `cargo install gloci` works. <ul><li>The internal crates say "no stable API" (§2).</li><li>The genome tables move into `crates/gloci-io/data/genomes/`. Today `builtin.rs` includes a path outside the package, which breaks `cargo publish` and the sdist.</li></ul> |
| D83 | **Python reuses the CLI wiring** through `gloci::embed`. Defaults, validation and messages are then identical, instead of being re-implemented and drifting. gloci-py holds no command logic. |
| D84 | **In-memory endpoints** use names that start with `"\0gloci:"`. No path contains NUL, so they can't collide with real files. <ul><li>They are read or written once, like stdin, and are never indexed.</li><li>Messages use their display name: `<list>`, `<generator>`, `<pyarrow.Table>`, `f.name`, or `<intersect output>`.</li></ul> |
| D85 | **Pipelines are immutable and lazy.** <ul><li>Each step validates its kwargs when it is called (kwargs → argv → clap).</li><li>Terminals run each step on a thread, joined to the next by 4 × 1 MiB pipes.</li><li>The output equals the shell pipeline that `cli()` prints.</li><li>The earliest step whose error isn't a broken pipe or an "upstream step failed" is the one reported.</li><li>A consumer that stops early (`break`) ends the steps quietly.</li></ul> |
| D86 | **Kwargs.** <ul><li>A kwarg is the long flag name with `-`→`_`, or any one-letter short flag (for example `f`, `F`, `s`, `S`, `b`, `g`).</li><li>`True` means the flag is present; `False`/`None` mean it is absent; a list repeats the flag; any other value is `str()`.</li><li>Non-input positionals are kwargs named by their clap id.</li><li>Reserved on steps: output/format/index/compression/threads (these go to terminals), quiet/verbose (logger level), and the bound input.</li><li>Messages say `flag=` rather than `--flag`.</li></ul> |
| D87 | **API** (amended in Task 4's follow-up, 2026-11-12: every P3 command is placed, and `UNBOUND` is empty). <ul><li>Sources: `gloci.open(src, *, format=None)`, `gloci.windows(**kw)` and `gloci.random(**kw)`.</li><li>One step method per stream command (Task 4 `BINDINGS`), named by its path joined with `_` (`vcf_filter`, `pair_intersect`, `stats_jaccard`). The pipeline is the input:<ul><li>as `-a`: `intersect`, `subtract`, `closest`, `map`, `coverage`, `pair_intersect`, `pair_pair`, `stats_jaccard`, `stats_fisher` and `stats_reldist` (the summary table comes back through `to_bytes()` and the table terminals);</li><li>as the inputs of an n-ary command, which takes extra positional sources: `merge`, `partition`, `complement`, `sort`, `vcf_concat`, `vcf_merge`, `cat`, `depth` and `samples`;</li><li>as the one positional: `view`, `convert`, `slop`, `flank`, `shift`, `cluster`, `spacing`, `groupby`, `genomecov`, `vcf_filter`, `vcf_extract`, `vcf_split`, `vcf_norm`, `vcf_annotate`, `getfasta` and `nuc` (with `fasta=`), `sample`, `shuffle`, `faidx` (with `regions=`), `dict`, `flagstat`, `idxstats`, `head`, and `bedcov` (with `files=` for the BAMs);</li><li>as a flag: `maskfasta` (`--fasta`, with `b=` for the masks) and `windows` (`-i`);</li><li>after another positional: `reheader`, whose input follows the required `header=`.</li></ul>`b=` takes any source or a list of them, including Pipelines. Visible aliases (`diff_indv_map=`) are kwargs too.</li><li>Terminals are the spec's `save`, `to_arrow`, `to_polars`, `to_pandas`, `count` and iteration, plus `to_bytes()` (the CLI's exact stdout) and `cli()`.</li><li>Module functions: `vcf_stats`, `vcf_extract_matrix`, `vcf_compare` (prefix outputs, like `vcf_stats`) and `index`.</li><li>`genomes`, `check` and `flags` (a lookup with no input) stay CLI-only. `idxstats` and `head` gained `-o` so they can be steps.</li></ul> |
| D88 | **Sources.** <ul><li>Paths.</li><li>File objects (anything with `.read`; text is encoded as UTF-8).</li><li>Iterables of `str` lines, `gloci.Record`s, or tuples (tab-joined `str()` values, with `None` written as `.`).</li><li>Objects with `__arrow_c_stream__`: integer `chrom start end` columns come first, other nulls become `.`, and a null coordinate or a tab is a `DataError` that names the column and row.</li></ul>`.line` is the 1-based item or row number. |
| D89 | **Outputs.** <ul><li>Iteration yields `gloci.Record`s whose coordinates come from re-reading the output with gloci's readers. Unrecognised output gives `None` coordinates.</li><li>`to_arrow()` builds the whole table in memory:<ul><li>`start`/`end` are int64;</li><li>VCF columns are `pos id ref alt qual filter info format <samples>`; GFF columns are `source feature score strand frame attributes`;</li><li>other columns are named by a `#` header line, else `col<N>`; `names=` overrides.</li></ul></li><li>Column types are inferred: int64, then float64, then string.</li><li>Arrow is handed over through the PyCapsule interface. pyarrow, polars and pandas are optional extras, and a missing one raises an `ImportError` that names the extra.</li></ul> |
| D90 | **Exceptions** follow §8.5. <ul><li>`UnsortedError` adds `.records`, `.lines` and `.hint`.</li><li>`str(exc)` is the CLI's message without `error: `.</li><li>A panic raises `GlociError` with the internal-error text.</li><li>Exceptions from user sources propagate unchanged.</li></ul> |
| D91 | **Runtime.** <ul><li>A driver loop holds the GIL only to pump sources (≤ 1 MiB each), forward logs to the logger `gloci` and check signals. It releases the GIL for waits of ≤ 50 ms.</li><li>Cancellation is a per-run `gloci_io::cancel::CancelToken` (amended in Task 3, fix round 1). The host makes one token per run and enters it (`RunGuard::enter`) on the thread that calls Task 4 `execute`, which installs it on every thread of the run (its rayon pool's workers included); the run's error is judged under it (`classify`, or `classify_under` on another thread). Task 5 `StageSet` owns it, so `StageSet::cancel()` is `token.cancel()`. Readers, outputs, sorters and pools capture the thread's token when they are built, and threads a run spawns enter it at their start. `KeyboardInterrupt` calls it, which makes that run's reads and writes fail (reported as `Cancelled`), removes partial output files and re-raises; other runs are unaffected, and there is no global flag to clear.</li><li>Each run gets its own rayon pool.</li><li>The CLI's signal and panic handlers are never installed.</li><li>`[profile.release-py]` sets `panic = "unwind"`.</li></ul> |
| D92 | **Wheels:** <ul><li>abi3 for CPython ≥ 3.10;</li><li>manylinux2014 (glibc 2.17) x86_64 and aarch64, and macOS arm64 (11.0);</li><li>an sdist elsewhere (needs Rust ≥ 1.91);</li><li>on Windows, `compile_error!` names WSL.</li></ul> |
| D93 | **Bioconda has two recipes.** <ul><li>`gloci` is the CLI. It builds from the tag tarball, has no run requirements, installs man pages and completions, and targets linux-64, linux-aarch64 and osx-arm64.</li><li>`py-gloci` is the library. It builds from the PyPI sdist and ships no `bin/gloci`.</li></ul>The `gloci` name stays with the CLI (§3). |
| D94 | **The book** is in `docs/book/`. <ul><li>It includes `migration.md`, `differences.md`, `CHANGELOG.md` and the gloci-py `README.md` (which is also the PyPI page).</li><li>The command reference and `SUMMARY.md` are generated and not committed.</li><li>PRs build the book; tags deploy it to Pages.</li></ul> |
| D95 | **Paper datasets** are pinned in `xtask/bench/datasets.toml` with provider checksums and `paper = true`. <ul><li>`[[derived]]` inputs are made by a command in the bench image and checked by sha256.</li><li>`paper/fetch-datasets.sh` is generated and committed.</li><li>CI never fetches any of these.</li></ul> |
| D96 | **The paper benchmark.** <ul><li>Real-data cases are in `xtask/paper/bench.toml`.</li><li>They use the existing runner and refer to nightly cases by name.</li><li>Every bench run records peak RSS with GNU `time`.</li><li>The Python benchmark has one binding target, `gloci-py < pybedtools`. polars-bio, PyRanges v1 and bioframe are yardsticks.</li></ul> |
| D97 | **Concordance.** <ul><li>`n_matching` is the multiset intersection of record lines (positional for `numeric`).</li><li>`match_rate = n_matching / max(n_gloci, n_reference)`, or 1.0 when both sides are empty.</li><li>Runs: a *default* run uses gloci's defaults, and a *compat* run adds `compat` flags.</li><li>A rate below 1.0 needs anchors: `expected_difference` for default runs, `compat_difference` for compat runs. Each anchor must exist in `docs/differences.md` and be pinned by a golden.</li><li>An unanchored mismatch fails, and is resolved through a minimal diff case.</li></ul> |
| D98 | **`cargo xtask paper`.** <ul><li>It regenerates `paper/tables/` and `paper/numbers.{json,tex}` from the newest `paper/results/<v>/`, purely and deterministically.</li><li>Subcommands: `fetch`, `run`, `bundle`.</li><li>`run` snapshots repository facts into `repo-snapshot.json`.</li><li>Zenodo gets the software from GitHub releases through its integration, which reads `CITATION.cff`, and the evidence through a deterministic bundle the owner uploads.</li></ul> |

## Review Focus

Each item has a test in the task named.

1. **Old glibc, no Rust.** A wheel must install with `--no-index` on glibc 2.17 with CPython 3.10, then import and run. `auditwheel show` must report `manylinux_2_17`. Task 11 `wheel-smoke`.
2. **Long operations.**
   - Other threads keep running during `save()`.
   - Ctrl-C raises within 2 s, leaves no `*.tmp-*` file, and the next call works.
   - A panic raises `GlociError` and the interpreter survives.

   Task 6 `test_gil_released`, `test_keyboard_interrupt_cancels` and `test_panic_is_glocierror`.
3. **Huge results.** Iterating 5M records grows RSS by less than 200 MB. A `break` stops every worker within 1 s and leaves no temp file. Task 8 `test_iteration_memory_flat` and `test_break_stops_workers`.
4. **Misbehaving inputs.** Task 7 `test_generator_exception_propagates`, `test_text_and_binary_file_objects`, `test_tuple_values` and `test_frame_sources` cover:
   - a generator that raises mid-stream (the same object propagates);
   - text and binary file objects;
   - tuples that hold `None`, floats or tabs;
   - frames with null or int32/uint64 coordinates.
5. **conda channel priority.** `gloci` and `py-gloci` must co-install with `--strict-channel-priority -c conda-forge -c bioconda`, report the same version, and clash on no file. Task 12: the nightly `conda-recipe` job.

## Lanes

| Lane | Tasks | Cross-lane waits |
|---|---|---|
| A: Python library | 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10 | — |
| B: Packaging & docs | 1 → 2 → 11 → 12 → 13 | 11 and 13 start after Task 6 merges |
| C: Paper evidence | 14 → 15 → 16 → 17 → 18 | 18 starts after Task 8 merges |

**Shared files that will conflict:**
- `Cargo.lock` (Tasks 1, 6, 8, 17): never hand-merge it. Take one side, then run `cargo update --workspace`.
- Root `Cargo.toml`: Task 1 changes the metadata; Task 6 changes the profile.
- `xtask/src/main.rs` and `lib.rs`: Tasks 2, 9 and 11–17 append subcommands and modules.
- `xtask/tests/workflows.rs`: each task appends its own test, and only Task 11 edits the release test.
- `release.yml`: Tasks 2, 11 and 12 (lane B).
- `nightly.yml`: Tasks 12 and 18.
- `CHANGELOG.md`: every task after Task 2 adds one bullet under `[Unreleased]`.
- `.gitignore`: Tasks 6, 13 and 14.
- `gloci-io/src/lib.rs` and `builtin.rs`: Tasks 1, 3 and 4.
- `crates/gloci/src/lib.rs`: Task 4. Any late P2a command arm moves into `embed`.
- `SKILL.md` and `cli.rs`: Task 1 only.
- `docs/differences.md`, `xtask/diff/cases.toml` and the goldens: Task 16 only.
- `xtask/bench/*` and `Dockerfile.bench`: Tasks 14, 15 and 18 (lane C).

## File structure

```
LICENSE README.md CHANGELOG.md CITATION.cff                         T1–T2
crates/gloci-io/{data/genomes/, src/memio.rs, src/cancel.rs}        T1, T3
crates/gloci/src/embed/{mod,spec,args,bindings,pipe,stages}.rs      T4–T5
crates/gloci-py/{Cargo.toml,pyproject.toml,README.md,src/,python/gloci/,tests/}  T6–T10
packaging/bioconda/{gloci,py-gloci}/   docs/book/                   T12, T13
xtask/src/{release,pystub,conda,book}.rs, paper/; xtask/paper/; xtask/bench/py/
paper/{README.md,fetch-datasets.sh,results/,tables/,numbers.*}      T14–T17
.github/workflows/{python,docs}.yml (new); release, nightly, ci
```

---

### Task 1: Licence, crate metadata, crates.io readiness (D79, D82)

**Files:**
- Create: `LICENSE` (MIT, "Copyright (c) 2026 Duarte Molha") and `README.md` (install lines, one example, links to the book and migration guide).
- Modify:
  - root `Cargo.toml`: `[workspace.package]` gets `license = "MIT"`, `repository = "https://github.com/duartemolha/gloci"`, `homepage = "https://duartemolha.github.io/gloci/"` and `publish = true`; every internal path dependency gets `version = "0.1.0"`;
  - crate manifests: `description`, `license.workspace` and `repository.workspace`; `gloci` also gets `readme = "../../README.md"`; testkit and xtask get `publish = false`;
  - `builtin.rs`: `concat!("../data/genomes/", $id, ".genome")`;
  - `crates/gloci-io/tests/builtin_genomes.rs` and `xtask/src/genomes.rs` (`DEFAULT_DIR`);
  - the help text at `cli.rs:77`, then regenerate SKILL.md;
  - `ci.yml`: add a `publish-dry-run` job;
  - `xtask/tests/workflows.rs`: the SHA-pin check runs over all workflows.
- Move: `git mv data/genomes crates/gloci-io/data/genomes`.
- Test: `xtask/tests/packaging.rs`.

- [ ] **Step 1: Write the failing tests:**
  - `published_crates_have_metadata` (from `cargo metadata`): the six crates have a description, license `MIT` and a repository; testkit and xtask have `publish == []`.
  - `internal_dep_versions_match_workspace`.
  - `licence_file_matches_metadata`.
  - `no_include_escapes_its_crate`: every `include_str!`/`include_bytes!` path, canonicalised, stays inside its crate.
  - `every_workflow_action_is_pinned`.
- [ ] **Step 2: Run** `cargo test -p xtask --test packaging`. Expected: FAIL.
- [ ] **Step 3: Implement.** The CI job runs `cargo publish --dry-run --locked -p gloci-core -p gloci-io -p gloci-engine -p gloci-vcf -p gloci-ops -p gloci`. Cargo ≥ 1.90 publishes several packages at once and resolves workspace dependencies.
- [ ] **Step 4: Run:**
  - `cargo test -p xtask --test packaging --test docs`;
  - `cargo test -p gloci-io --test builtin_genomes`;
  - `cargo package -p gloci-io --list | grep data/genomes/hg38.genome`.

  Expected: PASS.
- [ ] **Step 5: Commit** `build: MIT licence, crate metadata, genomes inside gloci-io (crates.io ready)`

### Task 2: Versioning, changelog, release checks, CITATION.cff (D81, D98)

**Files:**
- Create: `CHANGELOG.md` (`## [Unreleased]`, then `## [0.1.0] - 2026-11-08` with five bullets on P0–P1), `CITATION.cff` and `xtask/src/release.rs`.
- Modify:
  - `xtask/src/main.rs`: `ReleaseCheck { --tag }` and `ReleaseNotes { --version }`;
  - `release.yml`: the `version` job runs `cargo xtask release-check --tag "$GITHUB_REF_NAME"`, and the release step's `body_path` comes from `release-notes`;
  - `xtask/tests/workflows.rs`.
- Test: `xtask/tests/release.rs`.

**Interfaces:**
- `pub struct Section { pub version: Option<String>, pub date: Option<String>, pub body: String }`
- `pub fn changelog_sections(text: &str) -> Result<Vec<Section>>`
- `pub fn release_notes(changelog: &str, version: &str) -> Result<String>`
- `pub fn release_check(root: &Path, tag: Option<&str>) -> Result<Vec<String>>` returns problems, each with its fix. Tasks 11 and 12 extend it.

`CITATION.cff`:
- `cff-version: 1.2.0`, `type: software`, `title: gloci`;
- author Duarte Molha;
- `version`, `date-released`, `license: MIT`;
- `repository-code`, `url` (the book);
- `abstract`: the spec's first sentence.

- [ ] **Step 1: Write the failing tests:**
  - `changelog_parses_keep_a_changelog`: Unreleased comes first, versions decrease, and dates are ISO.
  - `release_notes_extracts_section`: the result is the exact body. An unknown version is an error that names the missing heading.
  - `release_check_reports_each_mismatch`: on a fixture repo, each of these gives one problem:
    - the tag doesn't match the version;
    - the version has no section, or a section without a bullet;
    - `CITATION.cff` has the wrong `version` or `date-released` (the date must be the section's);
    - there is no `LICENSE`;
    - the internal dependency versions don't match.
  - `repo_metadata_is_consistent`: the real repo has no problems.
- [ ] **Step 2: Run** `cargo test -p xtask --test release`. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p xtask --test release --test workflows`. Expected: PASS.
- [ ] **Step 5: Commit** `build: one version (D81), CHANGELOG, CITATION.cff, xtask release-check/release-notes`

### Task 3: In-memory endpoints and cancellation in gloci-io (D84, D91)

**Files:**
- Create: `crates/gloci-io/src/{memio.rs, cancel.rs}`.
- Modify:
  - `input.rs`: `open_input` resolves mem paths and wraps every source in `CancelRead`;
  - `output.rs`: a mem path gives `direct: true` and `Sink::Memory`/`MemoryBgzf`, and `Output::write` checks for cancellation;
  - `error.rs`: add `IoError::Cancelled` ("interrupted");
  - `lib.rs`.
- Test: `crates/gloci-io/tests/memio.rs`.

**Interfaces:**
- `pub const MEM_PREFIX: &str = "\0gloci:";`
- `pub fn register_input(display: &str, src: Box<dyn Read + Send>) -> String`
- `pub fn register_output(display: &str, sink: Box<dyn Write + Send>) -> String`
- `pub fn is_mem_path(&str) -> bool`
- `pub fn display_name(&str) -> String`
- `pub fn discard(&str)`
- `pub struct CancelToken` (`new()`, `cancel()`, `is_cancelled()`; clones share the flag), one per run.
- `pub struct RunGuard`, created with `RunGuard::enter(&CancelToken)`: installs the token on the thread and restores the previous one on drop.
- `pub struct Watch`: a token captured at construction (`Watch::current()`), polled with `check()`, installed on a spawned thread with `enter()`.
- `pub fn is_cancelled() -> bool` (this thread's token) and `pub fn is_cancellation(&dyn Error) -> bool`.

- [ ] **Step 1: Write the failing tests:**
  - `mem_input_reads_like_stdin`:
    - the records are read, and `name()` is `<list>`;
    - the path is not a regular file;
    - a second open fails with `can be read only once`.
  - `mem_output_receives_exact_bytes`:
    - the bytes arrive exactly;
    - with `-O bed.gz` they arrive BGZF-compressed;
    - `--write-index` is a `UsageError` containing `needs -o FILE`.
  - `mem_names_never_collide_with_files`: a real file named `gloci:1` is still opened as a file.
  - `cancel_fails_reads_and_writes`:
    - `read_batch` returns `IoError::Cancelled`;
    - a dropped atomic output leaves no `PATH.tmp-*`;
    - once the guard drops, reads work again.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test memio`. Expected: FAIL.
- [ ] **Step 3: Implement.** The registry is a `Mutex<HashMap>` keyed by an `AtomicU64` counter. Cancellation is one relaxed load per `read`, and the CLI never cancels.
- [ ] **Step 4: Run** `cargo test -p gloci-io && cargo test -p gloci --test golden --test behaviour`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(io): in-memory endpoints and cooperative cancellation`

### Task 4: `gloci::embed`: execute, specs, kwargs, bindings (D83, D86, D87)

**Files:**
- Create: `crates/gloci/src/embed/{mod.rs, spec.rs, args.rs, bindings.rs}`.
- Modify:
  - `crates/gloci/src/lib.rs`: `run` becomes `parse` → `logging::init` → `install_cleanup_handlers` → `execute` → `classify`. The per-command thread match moves to `embed::thread_and_tmp`.
  - `builtin.rs`: make `edit_distance` `pub` and re-export it.
- Test: `crates/gloci/tests/embed_args.rs`.

**Interfaces:**
- `pub fn parse(args: &[OsString]) -> Result<cli::Cli, UsageError>` (no ANSI, no "For more information").
- `pub struct ExecOptions { pub command_line: String, pub tmp_dir: Option<PathBuf> }`
- `pub fn execute(cli: &cli::Cli, o: &ExecOptions) -> anyhow::Result<bool>`. It runs in a per-call rayon pool sized by `-t`, and never prints or installs handlers.
- `pub enum Outcome { Usage(String), BrokenPipe, Failure(String) }` and `pub fn classify(&anyhow::Error) -> Outcome`.
- `pub struct ArgSpec { id, long: Option<String>, short: Option<char>, takes_value, multiple, positional, required, possible_values: Vec<String>, help }`
- `pub struct CommandSpec { path: Vec<String>, about: String, args: Vec<ArgSpec> }`
- `pub fn command_specs() -> Vec<CommandSpec>`, one per leaf command.
- `pub enum KwValue { Flag(bool), Value(String), Values(Vec<String>), Absent }`
- `pub fn kwargs_to_argv(spec: &CommandSpec, kwargs: &[(String, KwValue)]) -> Result<Vec<OsString>, UsageError>`
- `pub fn rewrite_flags_as_kwargs(msg: &str, spec: &CommandSpec) -> String`
- `pub enum Bind { Flag(&'static str), Positional { nary: bool } }`
- `pub struct Binding { path: &'static [&'static str], method: &'static str, input: Bind }`
- `pub const BINDINGS: &[Binding]`:
  - `-a`: `intersect`, `subtract`, `closest`, `map`, `coverage`;
  - n-ary positional: `merge`, `partition`, `complement`, `sort`, `vcf concat`;
  - single positional: `view`, `convert`, `slop`, `flank`, `shift`, `cluster`, `spacing`, `groupby`, `genomecov`, `vcf filter`, `vcf extract`;
  - `-i`: `windows`;
  - the method name is the path joined with `_`.
- `pub const MODULE_FUNCTIONS: &[(&[&str], &str)]`: `vcf stats`, `vcf extract` (`--genotype-matrix`), `index`, and `windows` (`-g`).
- `pub const CLI_ONLY`: `genomes`, `check`.
- (Follow-up: every P3 command is placed as D87, amended, lists; `flags` joins `CLI_ONLY`, `reheader` uses `Bind::PositionalAfter("header")`, and `UNBOUND` is empty.)

- [ ] **Step 1: Write the failing tests:**
  - `parse_errors_are_usage_errors`: `--bogus` gives a `UsageError` that names the flag and contains no `\x1b`.
  - `execute_equals_golden`: for the exact goldens `intersect-pairs-vs-bedtools`, `merge-default-vs-bedtools`, `closest-d-vs-bedtools` and `partition-all-vs-bedops-intersect`, `execute` with `-o <mem>` gives the bytes in `expected`.
  - `kwargs_rules`:
    - `f=0.5` → `-f 0.5`;
    - `min_frac_a="50%"` → `--min-frac-a 50%`;
    - `s=True` → `-s`;
    - `split=False` → nothing;
    - `c=[a,b]` → `-c a -c b`;
    - `mode="bogus"` → clap's list of `any, none, pairs, clip, count`;
    - `output=` → `output is set by save(path)`;
    - `threads=` → `threads belongs to the terminal method`;
    - `fooo=1` → `unexpected keyword argument 'fooo'`, with `did you mean` when the edit distance is ≤ 2.
  - `rewrite_flags`: `'--max-gap <N>'` → `'max_gap='`.
  - `every_cli_command_is_bound_or_excluded`:
    - every leaf is in `BINDINGS` or `MODULE_FUNCTIONS`, or else only in `CLI_ONLY` (`windows` and `vcf extract` are in both of the first two);
    - every `Flag` names an arg that exists.
- [ ] **Step 2: Run** `cargo test -p gloci --test embed_args`. Expected: FAIL.
- [ ] **Step 3: Implement** the specs by walking `gloci::command()`.
- [ ] **Step 4: Run** `cargo test -p gloci`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci::embed (parse/execute split, command specs, kwargs→argv, bindings)`

### Task 5: Stage runtime (D85)

**Files:** Create `crates/gloci/src/embed/{pipe.rs, stages.rs}`. Test: `crates/gloci/tests/embed_stages.rs`.

**Interfaces:**
- `pub const PIPE_CHUNKS: usize = 4; pub const PIPE_CHUNK: usize = 1 << 20;`
- `pub fn pipe(chunks, chunk_bytes) -> (PipeWriter, PipeReader)`. Dropping the writer is EOF. `PipeWriter::fail(msg)` gives the reader `upstream step failed: {msg}`. Dropping the reader gives the writer `BrokenPipe`.
- `pub struct StageSpec { argv: Vec<OsString>, command_line: String, tmp_dir: Option<PathBuf> }`
- `pub struct StageError { stage: usize, error: anyhow::Error }`
- `StageSet`:
  - `spawn(Vec<StageSpec>)`
  - `#[doc(hidden)] spawn_fns(Vec<Box<dyn FnOnce() -> anyhow::Result<bool> + Send>>)`
  - `wait(&mut self, Duration) -> Option<Result<bool, StageError>>`
  - `cancel(&self)`
  - Drop joins the threads.
- Each stage runs under a `RunGuard` and `catch_unwind`; a panic becomes `IoError::Internal`.

- [ ] **Step 1: Write the failing tests:**
  - `chained_stages_equal_shell_pipeline`: `intersect --mode any` → `merge --max-gap 50` equals the same `sh -c` pipeline using the `cargo_bin` gloci.
  - `upstream_data_error_wins`: a bad line in stage 0 gives `stage == 0` and the CLI's message.
  - `downstream_early_close_is_quiet`: every stage returns `Ok`, and Linux `Threads:` is back to baseline within 1 s.
  - `pipe_is_bounded`: the writer blocks at 4 MiB.
  - `deterministic_across_threads`.
  - `panic_becomes_internal_error`: the message contains `this is a bug`.
- [ ] **Step 2: Run** `cargo test -p gloci --test embed_stages`. Expected: FAIL.
- [ ] **Step 3: Implement** the pipes as bounded crossbeam channels of `Vec<u8>`.
- [ ] **Step 4: Run** `cargo test -p gloci --test embed_stages --test embed_args`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: embed stage runtime (bounded pipes, error priority, cancellation)`

### Task 6: gloci-py: crate, path pipelines, errors, runtime (D87, D90–D92)

**Files:**
- Create:
  - `crates/gloci-py/Cargo.toml`: `publish = false`; `[lib] name = "_native"`, `crate-type = ["cdylib"]`, `test = false`; `pyo3` with `abi3-py310` and `extension-module` (so libpython is never linked); feature `test-hooks`.
  - `pyproject.toml`:
    - `name = "gloci"`, `requires-python = ">=3.10"`, `dynamic = ["version"]`, `license = "MIT"`, `dependencies = []`;
    - extras `arrow` (`pyarrow>=14`), `polars` (`polars>=1.0`), `pandas` (`pandas>=2.0`, `pyarrow>=14`) and `test`;
    - `[tool.maturin]`: `python-source = "python"`, `module-name = "gloci._native"`, `profile = "release-py"`.
  - `README.md`, `requirements-ci.txt`, `src/{lib,run,errors,logs}.rs`, `python/gloci/{__init__.py, errors.py, _pipeline.py, py.typed}`, a placeholder `__init__.pyi`, `tests/{conftest.py, test_pipeline_paths.py, test_errors.py, test_runtime.py}` and `.github/workflows/python.yml`.
- Modify: root `Cargo.toml` (`[profile.release-py] inherits = "release"`, `panic = "unwind"`), `.gitignore` and `xtask/tests/workflows.rs`.

**Interfaces:**
- Python:
  - `gloci.open`;
  - `Pipeline.save(path, *, output_format=None, write_index=False, compression_level=None, threads=None)`, `.to_bytes(*, threads=None)`, `.count(*, threads=None)` and `.cli()`;
  - step methods generated from `_native.bindings()`;
  - `__version__`;
  - the D90 classes.
- `_native`:
  - `command_specs()`
  - `bindings()`
  - `kwargs_to_argv(path, kwargs) -> list[str]`
  - `run_plan(plan: dict) -> object`. The plan is `{"stages": [{"argv", "command_line"}], "sources": {token: obj}, "terminal": "save"|"bytes"|"count"|"iter"|"arrow", "names"}`, and the tokens `"\0src:<k>"` and `"\0out"` become mem paths.
- Each step's `command_line` is its CLI equivalent:
  - argv0 is `gloci`;
  - paths are kept as given;
  - the stream input is shown as `-`, and other in-memory inputs as `<python object N>`;
  - `-o` appears only for `save`;
  - `-t` never appears.

  `cli()` joins the steps with ` | `.

- [ ] **Step 1: Write the failing tests.** `conftest.py` provides `gloci_bin` (`$GLOCI_BIN`, which is required when `CI=true`), `cli(args, cwd)` (argv0 `gloci`) and `pair(n, seed)`.
  - `test_save_equals_cli`: `open("a.bed").intersect(b="b.bed", mode="pairs").save("o.tsv")` gives the same bytes as the CLI with `-o o.tsv`.
  - `test_chain_equals_shell_pipeline`: `.intersect(b=…).merge(max_gap=50).to_bytes()` equals `sh -c` of `cli()`, which is exactly `gloci intersect -a a.bed -b b.bed | gloci merge --max-gap 50 -`.
  - `test_vcf_header_records_command`: the output contains `##glociCommand=gloci intersect -a calls.vcf -b t.bed -o out.vcf; Version=<__version__>`.
  - `test_lazy_and_immutable`, `test_count` and `test_threads_deterministic`.
  - `test_usage_error_is_eager`:
    - `fooo=1` raises a `UsageError` that is also a `ValueError`, with "did you mean";
    - `merge(max_gap="x")` names `max_gap=`.
  - `test_data_error_attrs`: line 3 of `bad.bed` gives `.path == "bad.bed"` and `.line == 3`, and `str(e)` matches the CLI's stderr.
  - `test_unsorted_error_attrs`: `.records`, `.lines == (2, 3)`, and `.hint` contains `gloci sort`.
  - `test_open_error`.
  - `test_gil_released`: during a `save()` over `pair(500_000)`, a background counter advances by at least 100.
  - `test_keyboard_interrupt_cancels`: `Timer(0.3, _thread.interrupt_main)` raises within 2 s, no `*.tmp-*` is left, and the next call works.
  - `test_panic_is_glocierror`: `_native._test_panic()`, with `test-hooks`.
  - `test_warnings_go_to_logging`: a chr-prefix mismatch logs a WARNING on `gloci`.
- [ ] **Step 2: Run** `maturin develop -m crates/gloci-py/Cargo.toml --features test-hooks && pytest crates/gloci-py/tests -m "not large"`. Expected: FAIL.
- [ ] **Step 3: Implement:**
  - the driver loop (D91);
  - the error mapping (D90);
  - log forwarding, with the max level set from the logger at each run;
  - the Windows `compile_error!`.

  `python.yml` runs on PRs and pushes to main:
  - matrix: ubuntu-24.04 with 3.10 and 3.14, ubuntu-24.04-arm with 3.12, macos-latest with 3.12;
  - steps: `cargo build -p gloci` (export `GLOCI_BIN`), `maturin develop --features test-hooks`, `pytest -m "not large"`;
  - actions pinned by SHA.
- [ ] **Step 4: Run** the Step 2 command and `cargo test -p xtask --test workflows`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(py): gloci-py crate, path pipelines, typed exceptions, GIL-free runtime`

### Task 7: Python-object sources (D88)

**Files:**
- Create: `crates/gloci-py/src/sources.rs`.
- Modify: `_pipeline.py`, adding `classify_source(obj)`. It returns Path, FileObject, Lines, Frame or Nested, or raises a `TypeError` that lists the accepted kinds.
- Test: `tests/test_sources.py`.

A nested `Pipeline` runs as its own chain and writes into its registered pipe.

- [ ] **Step 1: Write the failing tests:**
  - `test_iterable_of_tuples`: `("chr1",15,25,None)` is fed as `chr1\t15\t25\t.`, and the result equals the file version.
  - `test_tuple_values`: floats are written with `str()`, and a tab raises `DataError` with `.line == 2` naming `<list>`.
  - `test_generator_exception_propagates`: `RuntimeError("boom")` at item 1,000 is raised as the same object, and no threads are left after 1 s.
  - `test_text_and_binary_file_objects`: `StringIO`, `BytesIO`, `open(p,"rb")`, `open(gz,"rb")` and `gzip.open` give identical bytes, and errors name `f.name` or `<BytesIO>`.
  - `test_records_and_lines_roundtrip`: `xfail(strict=True)` until Task 8.
  - `test_frame_sources` (`importorskip` pyarrow, polars and pandas):
    - `name, chrom, start, end` is fed as `chrom start end name`;
    - int32 and uint64 are accepted;
    - a null `start` at row 4 gives `.line == 4`;
    - a missing `end` raises `UsageError`.
  - `test_pipeline_as_b`: `intersect(b=open("b.bed").slop(b=10, g="hg38"))` equals the two-command shell version.
  - `test_two_inline_inputs`.
- [ ] **Step 2: Run** `pytest crates/gloci-py/tests/test_sources.py`. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `pytest crates/gloci-py/tests -m "not large"`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(py): iterables, file objects, Arrow frames and pipelines as sources`

### Task 8: Records, iteration, Arrow output (D89)

**Files:**
- Create: `crates/gloci-py/src/{records.rs, arrow.rs}`.
- Modify:
  - `Cargo.toml` (arrow-rs);
  - `_pipeline.py`;
  - `python.yml`: a `polars-only` job, and a nightly `pytest -m large` step built with `--profile release-py`.
- Test: `tests/{test_outputs.py, test_large.py}`.

**Interfaces:**
- `gloci.Record`: a frozen pyclass with `chrom`, `start` and `end` (each optional), `fields: tuple[str, ...]` and `line`. Sequence access is over `fields`; `str()` returns `line`.
- `Pipeline.__iter__`.
- `to_arrow(*, names=None, threads=None)`, `to_polars`, `to_pandas`.
- `ArrowStream.__arrow_c_stream__`, which `polars.from_arrow` consumes without pyarrow.
- The output is re-read through `gloci_io::open_reader_from`.

- [ ] **Step 1: Write the failing tests:**
  - `test_iteration_records`:
    - BED: `start=10, end=20`;
    - a VCF insertion: `start == end`, and `vcf_span="ref"` gives the REF span;
    - GFF column 4 = 101 gives `start == 100`.
  - `test_iteration_generic_tsv`: `vcf_extract(info="DP")` rows have `chrom is None` and no `#` line.
  - `test_count_matches_iteration`.
  - `test_to_arrow_schema`:
    - BED6 `intersect(mode="count")` gives `chrom:string start:int64 end:int64 col4:string col5:int64 col6:string col7:int64`;
    - VCF columns per D89, with `qual` as float64 where `.` becomes null;
    - `names=` with the wrong length raises `UsageError`.
  - `test_arrow_capsule_protocol`: `pyarrow.RecordBatchReader.from_stream`.
  - `test_missing_extra_message`: `pip install gloci[pandas]`.
  - `test_polars_without_pyarrow`.
  - `test_large.py`, marked `large`:
    - `test_iteration_memory_flat`: 5M records raise `ru_maxrss` by less than 200 MB;
    - `test_break_stops_workers`: `Threads:` returns to baseline within 1 s, and the step's `tmp_dir=` is empty.
- [ ] **Step 2: Run** `pytest crates/gloci-py/tests/test_outputs.py`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - Batches are about 4 MiB.
  - Records are built under the GIL.
  - Arrow buffers the text, infers the types, then builds the columns.
  - Dropping the iterator joins the steps with the GIL released.
- [ ] **Step 4: Run** `pytest crates/gloci-py/tests -m "not large"`, with Task 7's `xfail` removed. Expected: PASS.
- [ ] **Step 5: Commit** `feat(py): Record iteration, count, Arrow/polars/pandas output`

### Task 9: Typed stubs, signatures, Python guide

**Files:**
- Create: `xtask/src/pystub.rs`.
- Modify:
  - `xtask/src/main.rs`: `Pystub { --out }`;
  - `__init__.pyi` (now generated);
  - `_pipeline.py`: runtime `__signature__` and `__doc__`;
  - `crates/gloci-py/README.md`: sources, steps, terminals, frames, errors, threads and the GIL, and a ten-row "From pybedtools" table;
  - `python.yml`: add `mypy --strict`.
- Test: `xtask/tests/pystub.rs` and `tests/{test_api_surface.py, test_readme_examples.py, typing_example.py}`.

**Interfaces:** `pub fn generate_pyi(out: &Path) -> Result<()>`. Type mapping:

| Argument | Type |
|---|---|
| bool | `bool = False` |
| value | `str \| int \| float \| os.PathLike[str] \| None` |
| possible values | `Literal[...] \| None` |
| multiple | `str \| Sequence[str] \| None` |
| `b` | `Source \| Sequence[Source]` |

Docstrings hold the about text, `Args:` with the clap help, and `CLI: gloci <path>`.

- [ ] **Step 1: Write the failing tests:**
  - `pyi_is_current`, whose message is `run: cargo xtask pystub`;
  - `test_runtime_methods_match_stub`: the method and kwarg names from `ast` equal those from `inspect.signature`;
  - `test_help_works`: `Minimum fraction of A covered` appears;
  - `test_readme_examples`: every `python` block runs, unless it is tagged `no-run`;
  - `typing_example.py` passes `mypy --strict`.
- [ ] **Step 2: Run** `cargo test -p xtask --test pystub && pytest crates/gloci-py/tests/test_api_surface.py`. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the Step 2 commands and `mypy --strict crates/gloci-py/tests/typing_example.py`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(py): generated typed stubs, runtime signatures, Python guide`

### Task 10: Golden parity and the Python behaviour matrix (§8.5 testing)

**Files:** Create `crates/gloci-py/tests/parity/translate.py`, `tests/test_parity.py` and `tests/test_matrix.py`.

**Interfaces:** `translate(args) -> Callable[[], bytes]` splits an argv into the bound input, the other inputs and the kwargs, using `_native.command_specs()` and `_native.bindings()`. It raises `UntranslatableError` for `CLI_ONLY` commands.

- [ ] **Step 1: Write the tests:**
  - `test_golden_parity[<case>]`, over every `crates/gloci/tests/golden/*`, each in a temp copy:
    - with `GLOCI_BIN` set, the result equals the CLI's stdout byte for byte;
    - otherwise it equals `expected_gloci`, or `expected` for exact cases with no difference, and other cases are skipped;
    - for exit-status cases, exit 1 must raise `DataError` (with `str` equal to the stderr) and exit 2 must raise `UsageError`.
  - `test_parity_coverage`: at least 95% of cases translate, and the untranslated ones are all `CLI_ONLY`.
  - `test_ci_requires_gloci_bin`.
  - `test_matrix[<method>]`, over every binding:
    - `threads` 1, 2 and 8 give identical bytes;
    - empty input gives empty or header-only output;
    - unsorted input raises `UnsortedError`;
    - truncated BGZF raises `DataError`.
- [ ] **Step 2: Run** `pytest crates/gloci-py/tests/test_parity.py crates/gloci-py/tests/test_matrix.py`. Fix each real failure in its owning module, with a regression test.
- [ ] **Step 3: Run** `pytest crates/gloci-py/tests -m "not large"`. Expected: PASS.
- [ ] **Step 4: Commit** `test(py): golden parity against the CLI, Python behaviour matrix`

### Task 11: Release workflow: man/completions, wheels, PyPI, crates.io (D80, D92)

**Files:**
- Modify:
  - `release.yml`;
  - `xtask/src/main.rs`: `Man { --out }` and `Completions { --out }`;
  - `release.rs`: the pyproject must use `dynamic = ["version"]`;
  - `xtask/tests/workflows.rs`.
- Create: `crates/gloci-py/tests/test_smoke.py`. With only the wheel installed it checks:
  - the import works;
  - `__version__` equals `$GLOCI_EXPECT_VERSION`;
  - `open([("chr1",1,5)]).slop(b=1, g="hg38").to_bytes() == b"chr1\t0\t6\n"`, which proves the genome data is in the wheel;
  - `DataError` is raised for bad input.

**Workflow:**
- Triggers: tag `v*`, and `workflow_dispatch` with `dry_run`.
- `version`: `release-check`.
- `build`: each archive gains `man/`, `completions/` (from `cargo run --locked -p xtask -- man --out` and `… completions --out`) and `CHANGELOG.md`.
- `wheels` (maturin-action, `--profile release-py`): manylinux2014 x86_64 and aarch64, and aarch64-apple-darwin with `MACOSX_DEPLOYMENT_TARGET=11.0`.
- `sdist`.
- `wheel-smoke`:
  - Linux runs in `quay.io/pypa/manylinux2014_<arch>`: install with `/opt/python/cp310-cp310/bin/pip install --no-index --find-links dist gloci`, then `test_smoke.py`, then again with the newest `cp3*`; `auditwheel show` must print `manylinux_2_17`.
  - macOS uses Python 3.10 and 3.14.
- `sdist-smoke`: `python:3.12` with Rust 1.91 installed by rustup.
- `publish-pypi`: environment `pypi`, `id-token: write`, `skip-existing: true`; runs only `if: vars.GLOCI_PUBLISH_PYPI == 'true'` (D80, off by default).
- `publish-crates`: for each crate in order, skip it if `cargo info <crate>@<v>` succeeds, else run `cargo publish -p … --locked` with `CARGO_REGISTRY_TOKEN`; runs only `if: vars.GLOCI_PUBLISH_CRATES == 'true'` (off by default).
- The GHCR image push runs only `if: vars.GLOCI_PUBLISH_GHCR == 'true'` (off by default).
- Publish jobs run only when `github.event_name == 'push'` and `LICENSE` exists.
- The GitHub Release comes last and attaches the archives, wheels, sdist and `SHA256SUMS`. It `needs` the publish jobs but must still run when they are skipped: use `if: ${{ !failure() && !cancelled() && github.event_name == 'push' }}`.

- [ ] **Step 1: Write the failing assertions:**
  - The release test requires: `maturin`, `manylinux: 2014`, both manylinux2014 images, `auditwheel show`, `manylinux_2_17`, `gh-action-pypi-publish`, `skip-existing: true`, `id-token: write`, `cargo publish`, `cargo info`, `workflow_dispatch`, `dry_run`, `man --out`, `completions --out` and `release-check`.
  - The Release step must come after both publish jobs.
  - The old ban on `cargo publish`/`bioconda` is removed.
  - `publish_steps_are_gated`: every PyPI, crates.io and GHCR publish step carries its `vars.GLOCI_PUBLISH_*` condition, and the Release job runs when they are skipped.
- [ ] **Step 2: Run** `cargo test -p xtask --test workflows`. Expected: FAIL.
- [ ] **Step 3: Implement.** Pin each new action by SHA, found with `gh api repos/OWNER/REPO/commits/TAG`.
- [ ] **Step 4: Run** `cargo test -p xtask --test workflows --test docs` and actionlint. Expected: PASS. After merge, dispatch a `dry_run` and paste the job list and the `auditwheel show` output.
- [ ] **Step 5: Commit** `ci(release): wheels (abi3, manylinux2014, macOS arm64), PyPI, crates.io, man/completions, dry run`

### Task 12: Bioconda recipes (D93)

**Files:**
- Create: `packaging/bioconda/gloci/{meta.yaml, build.sh, test/{a,b,expected}.bed}`, `packaging/bioconda/py-gloci/meta.yaml` and `xtask/src/conda.rs`.
- Modify:
  - `xtask/src/main.rs`: `Conda { Render { --version, --sha256-cli, --sha256-py, --local, --out } }`;
  - `release.rs`: check the recipe versions;
  - `nightly.yml`: add `conda-recipe`;
  - `release.yml`: render the recipes with the real sha256 and attach them;
  - `xtask/tests/workflows.rs`.
- Test: `xtask/tests/conda.rs`.

**Interfaces:** `pub fn render(recipe: &Path, version: &str, sha256: &str, local_source: Option<&Path>) -> Result<String>`. `--local` swaps the `url` for a `path`.

**`gloci` recipe:**
- source `https://github.com/duartemolha/gloci/archive/refs/tags/v{{ version }}.tar.gz`;
- `run_exports pin_subpackage('gloci', max_pin="x.x")`;
- build requirements: rust and c compilers, `cargo-bundle-licenses`; no run requirements;
- tests:
  - `gloci --version`;
  - `gloci genomes show hg38`;
  - `gloci intersect -a a.bed -b b.bed | cmp - expected.bed`;
  - the man page and the three completions exist;
- `license_file: [LICENSE, THIRDPARTY.yml]`;
- `recipe-maintainers: [duartemolha]`;
- `additional-platforms: [linux-aarch64, osx-arm64]`.

**`build.sh`:**
1. `cargo-bundle-licenses --format yaml --output THIRDPARTY.yml`
2. `cargo install --no-track --locked --root "$PREFIX" --path crates/gloci`
3. `cargo run --locked -p xtask -- man --out "$PREFIX/share/man/man1"`
4. install the completions into bash-completion, zsh `site-functions` and fish `vendor_completions.d`

**`py-gloci` recipe:**
- source: the sdist attached to the GitHub Release (PyPI is off by default, D80);
- host: `python >=3.10`, pip, `maturin >=1.9,<2`;
- tests:
  - `import gloci`;
  - an `open([...]).count()` call;
  - `test ! -e $PREFIX/bin/gloci`.

- [ ] **Step 1: Write the failing tests:**
  - `recipes_render_with_versions`: no placeholders are left.
  - `cli_recipe_has_no_run_requirements`.
  - `py_recipe_installs_no_cli`.
  - `recipe_test_files_exist`.
  - `release_check` reports a recipe version mismatch.
  - The nightly workflow test asserts `--strict-channel-priority`, `-c conda-forge -c bioconda` and `--override-channels`.
- [ ] **Step 2: Run** `cargo test -p xtask --test conda`. Expected: FAIL.
- [ ] **Step 3: Implement the `conda-recipe` job**, in a pinned `condaforge/miniforge3`:
  1. set strict priority;
  2. render with `--local`;
  3. run `conda build -c conda-forge -c bioconda --override-channels` for each recipe;
  4. create a fresh env with `--strict-channel-priority -c local -c conda-forge -c bioconda gloci py-gloci`;
  5. assert that `gloci --version` matches `gloci.__version__` and that no file is installed twice.
- [ ] **Step 4: Run** `cargo test -p xtask --test conda --test release --test workflows`. Expected: PASS. Dispatch the nightly job once and paste the result.
- [ ] **Step 5: Commit** `build: bioconda recipes for gloci (CLI) and py-gloci (library)`

### Task 13: The mdBook site (D94)

**Files:**
- Create:
  - `docs/book/book.toml` (`build-dir = "../../target/book"`, `site-url = "/gloci/"`);
  - `docs/book/SUMMARY.in.md`, which holds a `{{reference}}` placeholder;
  - `docs/book/src/` pages:
    - `introduction.md`;
    - `install.md`;
    - `quickstart.md` (the §3 examples);
    - `semantics.md` (§4 for users, with the §4.2 and §4.3 tables verbatim);
    - include-only pages: `migration.md`, `differences.md`, `changelog.md`, `python.md` (the gloci-py README) and `python-api.md` (the stub);
    - `maintainers/releasing.md`: the release checklist, plus the first-release steps (PyPI trusted publisher, `CARGO_REGISTRY_TOKEN`, bioconda PR, Zenodo, Pages).
  - `xtask/src/book.rs` and `.github/workflows/docs.yml`.
- Modify: `xtask/src/main.rs` (`Book { --src }`), `.gitignore` (`docs/book/src/reference/` and `docs/book/src/SUMMARY.md`) and `xtask/tests/workflows.rs`.
- Test: `xtask/tests/book.rs`.

**Interfaces:**
- `pub fn generate_book(src: &Path) -> Result<()>` writes:
  - one `reference/<path-with-dashes>.md` per leaf command, containing the about text, `## Usage` (`render_long_help()`) and `## Examples` (with the `docs.rs` parser);
  - `SUMMARY.md`, in clap order.
- `pub fn check_links(src: &Path) -> Result<Vec<String>>`

- [ ] **Step 1: Write the failing tests:**
  - `book_generates_every_command`.
  - `book_includes_resolve`.
  - `book_links_resolve`: every relative link and `#anchor` resolves, either to an explicit `{#id}` or an mdBook heading slug.
  - `semantics_tables_match_spec`: byte-equal to §4.2 and the §4.3 worked examples.
  - `install_page_names_every_channel`: installing from the GitHub Release archives and wheels (`pip install <wheel URL>`), `cargo install --git`, and a clearly marked "once published" section for `cargo install gloci`, `pip install gloci`, `conda install -c conda-forge -c bioconda gloci`/`py-gloci` and `ghcr.io/`.
  - `docs_workflow_is_valid`: builds on PRs, and deploys only on `v*` tags, with `pages: write` and `id-token: write`.
- [ ] **Step 2: Run** `cargo test -p xtask --test book`. Expected: FAIL.
- [ ] **Step 3: Implement.** CI installs mdbook with `taiki-e/install-action` (SHA-pinned, `tool: mdbook@<pinned>`), builds with `cargo xtask book && mdbook build docs/book`, and deploys with `deploy-pages`.
- [ ] **Step 4: Run** `cargo test -p xtask --test book --test workflows && cargo xtask book && mdbook build docs/book`, which must produce no warnings. Expected: PASS.
- [ ] **Step 5: Commit** `docs: mdBook site (semantics, migration, differences, reference, Python)`

### Task 14: Paper datasets, derived inputs, fetch script (D95)

**Files:**
- Modify:
  - `xtask/bench/datasets.toml`;
  - `xtask/src/bench.rs`: add `Dataset.paper`, `pub struct Derived { name, from: Vec<String>, command, sha256 }` and `pub fn ensure_derived(d, cache, image, record) -> Result<PathBuf>`;
  - `Dockerfile.bench`: add GNU `time`, `gawk` and `ucsc-bigwigtobedgraph`, pinned;
  - `xtask/src/main.rs`: `Paper { Fetch { --record, --only }, FetchScript }`.
- Create: `xtask/src/paper/{mod,fetch}.rs`, `paper/README.md` and `paper/fetch-datasets.sh` (generated).
- Test: `xtask/tests/paper_datasets.rs`.

Datasets, checked on the providers' sites on 2026-11-10:

| name | url | checksum |
|---|---|---|
| `gencode-v46-gtf` | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz` | md5 from `release_46/MD5SUMS` |
| `encode-ccre-v4` | `https://www.encodeproject.org/files/ENCFF420VPZ/@@download/ENCFF420VPZ.bed.gz` | md5 `ea1820bd71581ba01c94fc600fbf794d` |
| `encode-ctcf-k562-peaks` | `https://www.encodeproject.org/files/ENCFF769AUF/@@download/ENCFF769AUF.bed.gz` | md5 `7d086cac19c5311a77b7e21e3d931435` |
| `encode-ctcf-k562-fc` | `https://www.encodeproject.org/files/ENCFF979PWH/@@download/ENCFF979PWH.bigWig` | md5 `4a1da9d7819bc622c27cd60950d28abb` |
| `gnomad-v4.1-exomes-chr22` | `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz` | md5 from the `x-goog-hash` header (base64 → hex); `curl -sI` first |
| `1kg-30x-chr20` | existing | existing |

Derived inputs are built with `LC_ALL=C` in the bench image. `from` names a dataset, a derived input, or `repo:<path>`.

| name | command |
|---|---|
| `exons-unsorted.bed` | `zcat {in} \| gawk -F'\t' -v OFS='\t' '$3=="exon"{match($9,/gene_name "[^"]+"/);print $1,$4-1,$5,substr($9,RSTART+11,RLENGTH-12),0,$7}'` |
| `exons.bed`, `ccre.bed`, `ctcf-peaks.bed` | `[zcat \| cut -f1-6 \|] sort -k1,1 -k2,2n -s` |
| `ccre.bed.gz`, `exons.bed.gz` | `bgzip -c {in} > {out} && tabix -p bed {out}` |
| `ctcf-fc.bedgraph` | `bigWigToBedGraph {in} /dev/stdout \| sort -k1,1 -k2,2n -s` |
| `1kg-chr20.bcf`; `1kg-chr20.p1/.p2.vcf.gz` | `bcftools view -Ob`; `-Oz -r chr20:1-30000000` / `-r chr20:30000001-` |
| `hg38.lex.genome` | from `repo:crates/gloci-io/data/genomes/hg38.genome`: `sort -k1,1` |

Resolve the bioRxiv 2025.02.19.639167 datasets (§9.2; deferred by D53) through the paper's data statement. If they are not public, record that in a comment, as was done for coitrees.

- [ ] **Step 1: Write the failing tests:**
  - `paper_datasets_have_provider_checksums`;
  - `derived_commands_reference_known_inputs`;
  - `nightly_never_fetches_paper_data`;
  - `no_workflow_runs_paper_downloads`: no workflow contains `xtask paper` or `fetch-datasets.sh`;
  - `bench_and_diff_images_pin_same_reference_tools`;
  - `fetch_script_is_current`.
- [ ] **Step 2: Run** `cargo test -p xtask --test paper_datasets`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - `paper fetch` reuses `ensure_dataset`.
  - An empty derived sha256 fails unless `--record` is given, which writes it.
  - The script is POSIX sh: curl, checksum, derive, with the tool versions in comments.
- [ ] **Step 4: Run** `cargo test -p xtask --test paper_datasets`. Expected: PASS. Then, on a workstation, run `cargo xtask paper fetch --record`, commit the digests, and paste the sizes and record counts.
- [ ] **Step 5: Commit** `feat(xtask): paper datasets (1kg, gnomAD, GENCODE, ENCODE) with derived inputs`

### Task 15: Real-data benchmarks and peak RSS (D96)

**Files:**
- Modify `xtask/src/bench.rs`:
  - `Case` gains `inputs: BTreeMap` (`a`, `b`, `c`), `label`, `table` and `threads`;
  - new placeholders `{in_a}`, `{in_b}`, `{in_c}`, `{genome}` and `{python}`;
  - after hyperfine, one run per command under `/usr/bin/time -f %M` writes `<case>.rss.json`.
- Create: `xtask/paper/bench.toml`.
- Test: `xtask/tests/bench_paper.rs`.

`from_nightly` lists the nightly 1kg `vcf filter` case and Task 20's `vcf stats`, `partition` and `genomecov` cases, by their committed names.

BED cases use `--chrom-order lex`. Yardsticks are marked (y).

| case | inputs | gloci | baselines | target |
|---|---|---|---|---|
| `sort-exons` | exons-unsorted | `sort` | sort-bed, bedtools sort; GNU sort (y) | `gloci < sort-bed` |
| `intersect-ccre-exons` | ccre, exons | `intersect --mode any -t 1` | `bedtools -u -sorted`, `bedops -e 1` | `≤ 0.5 × bedtools`, `≤ bedops` |
| `intersect-ccre-exons-bgzf` | `.bed.gz` copies | `-t 8` | same | `≤ 0.2 × bedtools` (§9.1) |
| `subtract-exons-ccre`, `merge-exons` | exons, ccre | `-t 1` | bedtools; `bedops -d` / `-m` | as intersect |
| `closest-peaks-exons` | peaks, exons | `closest --distance -t 1` | `bedtools closest -d`, `closest-features --dist` | as intersect |
| `map-ccre-signal` | ccre, ctcf-fc | `map -c 4:mean -t 1` | `bedtools map -o mean`, `bedmap --mean` | as intersect |
| `coverage-ccre-exons`, `genomecov-exons` | ccre, exons | `coverage`; `genomecov --bga -g {genome}` | bedtools | `≤ 0.5 × bedtools` |
| `partition-3` | ccre, exons, peaks | `partition --no-membership` | `bedops --partition`, `bedtools multiinter` | `≤ bedops` |
| `intersect-1kg-exons` | 1kg, exons | `-t 8` | `bedtools -u -header`; `bcftools view -R` (y) | `≤ 0.2 × bedtools` |
| `vcf-filter-gnomad` | gnomad | `vcf filter --pass-only -t 8` | `vcftools --remove-filtered-all --recode`; `bcftools view -f PASS` (y) | `≤ 0.2 × vcftools` |
| `vcf-filter-1kg-bcf` | 1kg BCF | as nightly | vcftools `--bcf` | `≤ 0.2 × vcftools` |

- [ ] **Step 1: Write the failing tests:** `paper_cases_parse_and_reference_known_inputs`, `every_paper_gloci_command_parses`, `from_nightly_names_exist`, `rss_file_parses` and `placeholders_substitute`.
- [ ] **Step 2: Run** `cargo test -p xtask --test bench_paper`. Expected: FAIL.
- [ ] **Step 3: Implement** `paper run --bench-only [--only CASE] --gloci BIN`, which writes to `paper/results/<v>/bench/`.
- [ ] **Step 4: Run** `cargo test -p xtask`. Expected: PASS. On a workstation, run `--only intersect-ccre-exons` and paste the medians and RSS.
- [ ] **Step 5: Commit** `feat(xtask): real-data paper benchmarks and peak RSS in every bench run`

### Task 16: The real-data concordance run (D97)

**Files:**
- Modify: `crates/gloci-testkit/src/differential.rs`.
- Create: `xtask/src/paper/concord.rs` and `xtask/paper/concordance.toml`.
- Modify as findings require: diff cases, goldens, `docs/differences.md` and `docs/migration.md`.
- Test: `xtask/tests/concordance.rs`.

**Interfaces:**
- `pub struct MatchStats { n_gloci, n_reference, n_matching: u64 }`, with `rate()`.
- `pub fn match_stats(mode, gloci: &mut dyn BufRead, reference: &mut dyn BufRead) -> io::Result<MatchStats>`. It streams, using a multiset of line hashes; `vcf-records` drops `#` lines.
- `concord::run` writes `concordance.csv` (`case,command,reference_tool,dataset,mode,n_gloci,n_reference,n_matching,match_rate,verdict,anchors`) and `concordance/<case>.<mode>.txt`, holding the first 10 mismatches from each side.
- The case schema is: `name`, `command`, `reference_tool`, `inputs`, `gloci`, `reference`, `compare`, `compat`, `expected_difference` and `compat_difference`.
  - Reference commands and compat flags copy the committed diff case for the same command.
  - Default runs add `--chrom-order lex`.
  - An anchor that never fires prints `anchor unused`.

| cases | reference | inputs | compare | anchors / compat |
|---|---|---|---|---|
| sort | sort-bed; bedtools sort | exons-unsorted | exact | compat `--full-key`; `#sort-tie-order` |
| intersect: any, none, pairs, clip, count, `-f 0.5 -r` | bedtools `-u/-v/-wo/(none)/-c/-f 0.5 -r`; bedops `-e/-n` | ccre, exons | exact | — |
| subtract, `--whole` | bedtools `subtract/-A`; `bedops -d` | exons, ccre | exact | — |
| closest, `--each-side` | `bedtools closest -d -t all`; `closest-features` | peaks, exons | exact | `#closest-ties-k`, `#closest-each-side` |
| map; bedmap | `bedtools map -o mean,max,count`; `bedmap --mean` | ccre, fc | numeric | `#map-empty-value` |
| coverage, `--hist`, `--mean` | bedtools coverage | ccre, exons | exact/numeric | `#coverage-point` |
| merge, `-c 4:count_distinct` | bedtools merge; `bedops -m` | exons | exact | — |
| complement, slop, flank, shift, windows, genomecov `--bga` | bedtools | exons, genome | exact | `#complement-past-end`, `#slop-flank-past-end`, `#shift-past-edge`, `#genomecov-chrom-order` |
| partition | `bedops --partition`; reshaped multiinter | ccre, exons, peaks | exact | `#partition-boundaries` |
| cluster, spacing, groupby | bedtools | exons(-unsorted) | exact/numeric | `#cluster-strand` |
| convert GTF | gtf2bed | gencode | exact | `#gff2bed-name` |
| intersect VCF vs BED | `bedtools -u -header` | 1kg, exons | vcf-records | `#vcf-allele-span`; compat `--vcf-span ref` |
| view (10 regions) | tabix | 1kg | exact | `#view-region-footprints`; compat `--vcf-span ref --inclusive-points` |
| vcf filter MAF+missing; `--pass-only`; `--include` | vcftools `--recode --recode-INFO-all`; `bcftools view -f PASS/-i` | 1kg; gnomad | vcf-records | `#maf-multiallelic`, `#expr-subset` |
| vcf stats: freq, counts, missing-site/-indv, het, hardy, tstv, singletons, snp-density, indel-hist | vcftools | 1kg | numeric | `#maf-multiallelic`, `#hwe-multiallelic`, `#all-missing-site-stats` |
| vcf extract `--info`, 012 (100 samples); vcf concat | vcftools `--get-INFO/--012`; bcftools concat | 1kg (p1, p2) | exact/vcf-records | `#vcf-extract-info`, `#vcf-extract-012` |

- [ ] **Step 1: Write the failing tests:**
  - `match_stats_modes`:
    - identical inputs give 1.0;
    - one changed line in 10 gives 9;
    - reordered lines give 1.0;
    - `vcf-records` ignores headers;
    - `numeric` treats `3` and `3.0` as equal;
    - two empty inputs give 1.0.
  - `concordance_cases_parse`.
  - `concordance_anchors_exist_and_are_pinned`: each anchor appears in some golden `expected_difference`.
  - `every_comparable_command_has_a_concordance_case`: every available command in `docs/migration.md` has a case. Exemptions, each with its reason, live in a constant.
  - `verdict_rules`.
- [ ] **Step 2: Run** `cargo test -p gloci-testkit && cargo test -p xtask --test concordance`. Expected: FAIL.
- [ ] **Step 3: Implement** `paper run --concordance-only [--only CASE]`, in the bench image with the cache mounted read-only.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Run the full concordance on a workstation. For each unanchored mismatch:
  1. cut a minimal repro and add it as a diff case;
  2. fix gloci, or add an anchored section and its migration link;
  3. run `cargo xtask diff --update-golden`;
  4. re-run.

  Paste the case count, the number of compat cases at 1.0, and the anchors hit.
- [ ] **Step 5: Commit** `feat(xtask): real-data concordance against bedtools/bedops/vcftools/bcftools/htslib`

### Task 17: `cargo xtask paper`: tables, numbers, Zenodo bundle (D98)

**Files:**
- Create: `xtask/src/paper/{tables,numbers,bundle}.rs`, `xtask/tests/fixtures/paper-results/0.0.0-fixture/` (two bench JSONs with RSS, a 4-row `concordance.csv`, `environment.json` and `repo-snapshot.json`) and `xtask/tests/fixtures/paper-expected/`.
- Modify:
  - `xtask/src/main.rs`: bare `paper` runs the tables; add `Run { --gloci, --wheel, --bench-only, --concordance-only, --only }` and `Bundle { --results, -o }`;
  - `xtask/Cargo.toml`: add `tar`;
  - `paper/README.md`.
- Test: `xtask/tests/paper_tables.rs`.

**Interfaces:**
- `paper run` writes two files:
  - `environment.json`: CPU, cores, RAM, kernel, image digests, every tool's `--version`, the gloci version and commit, and `started_at`/`finished_at`;
  - `repo-snapshot.json`: golden counts per tool (from the `-vs-<tool>` suffix), anchors with titles and counts, migration rows, test counts, and dataset records and bytes.
- `pub fn tables(results, out) -> Result<()>`
- `pub fn numbers(results) -> Result<BTreeMap<String, Value>>`
- `pub fn bundle(results, out) -> Result<()>`
- The CSVs:

  | File | Columns |
  |---|---|
  | `T1-datasets` | `dataset,description,release,url,checksum_algorithm,checksum,bytes,records` |
  | `T2-benchmarks` | `case,label,dataset,tool,role,threads,median_s,mean_s,stddev_s,min_s,max_s,runs,max_rss_mb` |
  | `T3-speedups` | `case,label,baseline,speedup,target,target_met` |
  | `T4-concordance` | the concordance columns, minus `verdict` |
  | `T5-differences` | `anchor,title,golden_cases,concordance_cases` |
  | `T6-command-coverage` | `tool,old_command,gloci_command,status` |
  | `T7-tests` | `category,count` |
  | `T8-software` | `tool,version,role` |

  Formats: seconds `%.3f`, speedup `%.2f`, rate `%.6f`, MB `%.1f`. Rows are sorted by the first two columns.
- The number keys are `bench.<case>.<tool>.median_s|max_rss_mb`, `bench.<case>.speedup.<tool>`, `concordance.{cases,compat.at_100,default.mean_match_rate}`, `concordance.<case>.<mode>.match_rate`, `tests.golden_cases`, `tests.golden.<tool>`, `differences.anchors`, `datasets.<name>.records` and `migration.<tool>.{available,total}`.
- `numbers.tex` starts with `\providecommand{\glnum}[1]{\csname glnum@#1\endcsname}`, then one sorted line per key: `\expandafter\def\csname glnum@<key>\endcsname{<display>}`.

- [ ] **Step 1: Write the failing tests:**
  - `tables_from_fixture_match_expected`;
  - `tables_are_deterministic`;
  - `numbers_tex_defines_every_key`;
  - `bundle_is_reproducible`:
    - two runs are byte-identical;
    - the files are exactly the results, tables, numbers, `datasets.toml`, the fetch script, both Dockerfiles, `xtask/paper/*.toml`, `CITATION.cff` and `README.txt`;
    - entries are sorted, mtime is `finished_at`, uid/gid are 0, and gzip mtime is 0;
  - `committed_tables_are_current`, only when `paper/results/` is non-empty.
- [ ] **Step 2: Run** `cargo test -p xtask --test paper_tables`. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p xtask`. Expected: PASS. Then, on a machine near the §9.1 reference (8-core x86_64, NVMe, 32 GB) and not the 4-core development host:
  1. `paper fetch`;
  2. `paper run --gloci <release musl> --wheel <x86_64 wheel>`;
  3. `cargo xtask paper`;
  4. commit `paper/results/<v>/`, `paper/tables/` and `paper/numbers.*`;
  5. `paper bundle`.

  Paste T3 and the T4 summary.
- [ ] **Step 5: Commit** `feat(xtask): cargo xtask paper (tables, numbers, Zenodo evidence bundle)`

### Task 18: Python library benchmarks (D96, §9.2)

**Files:**
- Create: `xtask/bench/py/<tool>_<op>.py`, for the tools gloci, pybedtools, polars_bio, pyranges1 and bioframe, and the ops intersect and merge. Each takes `IN_A [IN_B] OUT` and writes BED.
- Modify:
  - `xtask/bench/cases.toml`: nightly `py-intersect-10M` and `py-merge-10M`;
  - `xtask/paper/bench.toml`: `py-intersect-ccre-exons` and `py-merge-exons`;
  - `bench.rs`: `--wheel PATH`, which runs `pip install --no-deps` in the container;
  - `Dockerfile.bench`: Python 3.12 and the pinned Python tools;
  - `nightly.yml`: build an x86_64 wheel and pass `--wheel`.
- Test: `crates/gloci-py/tests/test_bench_scripts.py` and `xtask/tests/bench_paper.rs`.

Targets: `gloci-py < pybedtools` is binding (the §2 anti-goal); the others are yardsticks.

- [ ] **Step 1: Write the failing tests:**
  - `test_gloci_bench_scripts_match_cli`: the two gloci scripts write the exact bytes of `gloci intersect --chrom-order lex --mode any` and `gloci merge --chrom-order lex`.
  - `py_cases_parse`.
- [ ] **Step 2: Run** `pytest crates/gloci-py/tests/test_bench_scripts.py && cargo test -p xtask --test bench_paper`. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the Step 2 commands. Expected: PASS. Dispatch the nightly bench once, and paste the py medians.
- [ ] **Step 5: Commit** `bench: Python library vs pybedtools (yardsticks polars-bio, PyRanges v1, bioframe)`

---

## Execution notes for the controller

- **Owner first.** The licence is MIT (confirmed). Before the first public release, the owner confirms the copyright holder (default: Duarte Molha). The first public release also needs these owner-only steps:
  - nothing is published by default (D80). When the owner decides to publish, they create the PyPI project and its trusted publisher and set `GLOCI_PUBLISH_PYPI=true`;
  - for crates.io: add `CARGO_REGISTRY_TOKEN` and set `GLOCI_PUBLISH_CRATES=true`; for GHCR: set `GLOCI_PUBLISH_GHCR=true`;
  - enable the Zenodo GitHub integration;
  - enable GitHub Pages;
  - open the first bioconda PR, when the owner decides to.
- **Lanes.** Run the three lanes in separate worktrees, each with its own `target/`.
  - On this 4-core host, set `CARGO_BUILD_JOBS=1` and use scoped test commands. Run `cargo test --workspace` once per lane, before merging.
  - Tasks 11 and 13 wait for Task 6; Task 18 waits for Task 8.
- **Paper runs** (Step 4 of Tasks 14–17) need a workstation with Docker and about 60 GB free. They never run in CI.
- **Spec amendments.**
  - §8.5: one line on `to_bytes`/`cli` and the `py-gloci` name.
  - §13: the licence, once confirmed.
  - D59 needs no spec change.
