# Changelog

All notable changes to gloci are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/) (before 1.0, a breaking change bumps
the minor version).

## [Unreleased]

- P2a: the VCF commands `vcf filter`, `vcf extract`, `vcf stats` and `vcf concat`.
- P3a: BAM, CRAM, SAM and FASTA input; the samtools-style commands `sort`, `merge`, `cat`, `index` and `faidx`; and `getfasta`, `nuc` and `maskfasta`.
- P3b: `stats jaccard`, `fisher` and `reldist`; `random`, `sample` and `shuffle`; the BEDPE pair commands and `convert --bedpe`; parallel VCF selection; `vcf compare` (genotype concordance); and `vcf split`, `vcf norm`, `vcf annotate` and `vcf merge`.
- P3c: the samtools-style commands `flags`, `idxstats`, `flagstat`, `head`, `samples`, `depth` (with `--skip-mate-overlap`), `bedcov`, `dict` and `reheader`; `check` of alignment files and `check --quick` (quickcheck); and `view --count` and `view --regions-file` on alignments.
- P2b: `gloci::embed`, the in-process API the Python library is built on: the CLI's own parsing and dispatch with in-memory inputs and outputs, per-run cancellation, command specs and keyword arguments.
- P2b: the Python library `gloci` (crate `gloci-py`, built with maturin): lazy, immutable pipelines over files with one step method per stream command, run in-process with the CLI's output (`save()`, `to_bytes()`, `count()`, `cli()`); typed exceptions (`UsageError`, `DataError`, `UnsortedError`, `GlociError`); the GIL released while steps run, Ctrl-C cancelling the run, and warnings sent to the `logging` logger `gloci`.
- P2b: the documentation site (mdBook): install, quickstart, semantics, migration, differences, the generated per-command reference, the Python guide and API, and the maintainers' release checklist; `cargo xtask book` generates the reference and checks every link.
- P2b: Python-object sources for the Python library: file objects, iterables of lines or tuples, Arrow tables and streams (pyarrow, polars, pandas; read through the Arrow C stream interface) and nested pipelines, read as the run goes; the module sources `gloci.windows()` and `gloci.random()` and the functions `vcf_stats()`, `vcf_extract_matrix()`, `vcf_compare()` and `index()`. `sort(check=True)` is refused on the step: the check stays on the command line.
- P2b: records and tables from the Python library: iteration yields `gloci.Record`s (coordinates read back with gloci's readers, fields and line), handed over in 4 MiB batches with the steps stopped by `break`; `count()` counts what iteration yields, in Rust; `to_arrow()`, `to_polars()`, `to_pandas()` and `to_arrow_stream()` (a `gloci.ArrowStream`, through the Arrow PyCapsule interface; polars needs no pyarrow) build a table with `chrom start end` first and inferred column types. A step whose output is a table (`stats jaccard`, `vcf extract`, `idxstats`, ...) is marked in `gloci::embed`, so its header row names the columns rather than being a record.
- P2b: `idxstats` and `head` take `-o FILE` (atomic; a `.gz` name compresses).
- P2b: `flagstat` on a CRAM decodes only the data series it counts, and reads the same counts as before: a 258 MB CRAM of 7.5 million reads takes about 0.27 s where it took 13.4 s, against 0.45 s for `samtools flagstat` (output byte-identical, in text, JSON and TSV, at every `--threads`).
- P2b: the release archives carry the man pages, the shell completions and `CHANGELOG.md` next to the binary, `SKILL.md` and the licence.
- P2b: Python wheels (abi3, for CPython 3.10 and later; manylinux2014 on x86_64 and aarch64, macOS arm64), an sdist, and bioconda recipes for `gloci` and `py-gloci`; each release carries them with `SHA256SUMS`. Publishing to PyPI, crates.io, `ghcr.io` and GitHub Pages stays off until the owner enables it.
- P2b: the Python package is typed: `py.typed` and a generated stub (checked by `mypy --strict`), runtime signatures and docstrings from the CLI definition, and a Python guide in the book.
- P2b: the Python extras need pyarrow 15, polars 1.5 and pandas 2.2 or later (the lowest versions whose Arrow stream interface gloci uses); CI runs the tests on exactly those.
- P2b: a usage message that names the long option now says `--output-format` where it said `-O`.

- Stability tiers: `genomes`, `index`, `sort`, `view`, `intersect`, `merge`, `subtract`, `closest` and `coverage` are the stable core; every other command is experimental and its help starts with a line saying so. `contrib/gloci-shadow.sh` runs a pipeline step with the tool it uses today and with gloci, and logs whether the outputs are identical.
- `SKILL.md` is an index of 28 KB where it was 173 KB: the semantics, one entry per command and the pitfalls. The rules of a command stay in its `--help`, and the release archives now carry `docs/migration.md` and `docs/differences.md`, which it points to.
- Help text no longer cites sections of the internal design spec. `docs/migration.md` no longer lists `jaccard`, `fisher`, `reldist`, `shuffle`, `random`, `sample`, `pairtobed` and `pairtopair` as not yet available: they have been since P3b.
- Fuzzing: `cram_reader` and `fasta_index` targets; the nightly job runs every target (it left out `sam_reader` and `bedpe_reader`) for 180 s each.

## [0.1.0] - 2026-11-08

- Interval engine: half-open u64 coordinates, sorted-stream sweeps and an index for random access, with output byte-identical at every `--threads`.
- Interval commands with bedtools parity: intersect, merge, subtract, closest, coverage, genomecov, map, groupby, cluster, spacing, window, slop, flank and more, checked against pinned bedtools output through committed golden files.
- BED, GFF, VCF and BCF reading and writing, with bgzip and tabix (TBI and CSI) indexes, and built-in genome files.
- A single static binary (musl on Linux, plus macOS arm64) with man pages, shell completions and an agent usage guide (`SKILL.md`).
- A differential test harness (`cargo xtask diff`), a benchmark suite (`cargo xtask bench`) and synthetic data generators (`cargo xtask gen`).
