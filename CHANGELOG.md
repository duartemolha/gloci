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

## [0.1.0] - 2026-11-08

- Interval engine: half-open u64 coordinates, sorted-stream sweeps and an index for random access, with output byte-identical at every `--threads`.
- Interval commands with bedtools parity: intersect, merge, subtract, closest, coverage, genomecov, map, groupby, cluster, spacing, window, slop, flank and more, checked against pinned bedtools output through committed golden files.
- BED, GFF, VCF and BCF reading and writing, with bgzip and tabix (TBI and CSI) indexes, and built-in genome files.
- A single static binary (musl on Linux, plus macOS arm64) with man pages, shell completions and an agent usage guide (`SKILL.md`).
- A differential test harness (`cargo xtask diff`), a benchmark suite (`cargo xtask bench`) and synthetic data generators (`cargo xtask gen`).
