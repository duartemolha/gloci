# gloci P0 Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build spec Phase 0: the Section 4 semantics as tested library code, BED/bedGraph/GFF/VCF/BCF reading and writing with parallel BGZF and atomic output, the `sort`, `index`, `view` and `check` commands, and the benchmark, differential, property and fuzz harnesses running in CI.

**Architecture:** A Cargo workspace. `gloci-core` holds pure rules (coordinates, chromosomes, overlap, distance, footprints, order). `gloci-io` holds formats and files. `gloci-ops` implements commands as library functions. A thin clap binary `gloci` wires arguments. Every reader yields batches of `RecordMeta` (chromosome id, sort key, footprint, location, and a slice of the raw line), so records pass through unchanged. `xtask` holds the data generator and the benchmark and differential runners. `gloci-testkit` holds proptest strategies, naive reference implementations and fixture builders.

**Tech Stack:** Rust stable (edition 2024); clap 4.6 (+ clap_complete 4.6, clap_mangen 0.3); noodles-bgzf 0.52 (feature `libdeflate`), noodles-vcf 0.96, noodles-bcf 0.94, noodles-csi 0.62, noodles-tabix 0.68, noodles-core 0.21; flate2 1.1; memchr 2.8; itoa 1; rayon 1.12; crossbeam-channel 0.5; lz4_flex 0.14; tempfile 3.27; ctrlc 3.5 (feature `termination`); mimalloc 0.1; thiserror 2; anyhow 1; log 0.4; libc 0.2; proptest 1.12; insta 1.49; assert_cmd 2.2; predicates 3.1; rand 0.10 + rand_chacha 0.10; serde 1 + toml 1 + serde_json 1; cargo-fuzz; hyperfine.

**Spec:** `docs/superpowers/specs/2026-11-06-gloci-design.md`. This plan implements the P0 row of §10, using §4 (semantics), §5 (formats), §6 (CLI), §7.4 (commands), §8 (architecture) and §9 (testing).

**Scope:** P0 only. P1 (relations, `gloci-engine`, `gloci-vcf`), P2 and P3 get their own plans after P0 lands. Out of this plan: `--unsorted`, `--chrom-map`, per-chromosome parallelism, the mdBook, release binaries and bioconda.

**P0 is done when:** `cargo xtask bench --check-targets` shows `gloci sort` faster than `sort-bed`, and CI plus the nightly workflow run the property, differential, fuzz, benchmark and large-file harnesses (spec §10).

## Global Constraints

- Internally every feature is a **0-based, half-open** interval `[start, end)`; positions are **64-bit unsigned**.
- Rust stable, edition 2024. Release profile: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`.
- **No htslib.** Pure Rust plus vendored C; the binary must build as one static `musl` binary.
- Exit codes: **0** success, **1** data or runtime error, **2** usage error. Broken pipe (e.g. `| head`) exits quietly with status 0: no panic, no "Broken pipe" message.
- Errors say what, where and how to fix, e.g. `error: calls.vcf.gz line 1,204,331: POS "12a4" is not an integer` (line numbers with thousands separators).
- Output is byte-identical for every `--threads` value.
- `.gz` output is always **BGZF**. `--compression-level N` default **6**.
- Atomic file output: write `PATH.tmp-XXXX`, rename on success; on failure no partial file is left behind.
- Format detection uses content and magic bytes, never file extensions. `--format` overrides. `-` reads stdin.
- VCF output gets a `##glociCommand=<full command line>; Version=<x.y.z>` header line.
- `--vcf-span allele|ref`, default **`allele`**.
- `-t, --threads N` default `available_parallelism()`, env `GLOCI_THREADS`. `--tmp-dir DIR`, env `GLOCI_TMPDIR`. `sort --max-mem` default **2G**.
- Every command's `--help` has 2–4 real examples.
- Sortedness is checked on the fly in every streaming command and is never turned off. The first violation exits 1 and names the file, line number and both offending records, with a hint such as `run: gloci sort FILE -g <genome>`.
- The same flag means the same thing everywhere; flag names come from spec §6.1.
- Streaming memory stays flat as input grows (≤ ~200 MB RSS); `sort` obeys `--max-mem`.
- CI platforms: Linux x86_64 and aarch64, macOS aarch64.

## Decisions this plan adds

The spec is silent or ambiguous on these points. Each is the plan's call; change any of them before execution starts.

| # | Decision |
|---|---|
| D1 | `RecordMeta` carries both `key` (the position the file is sorted by: BED start, GFF start−1, VCF POS−1) and `span` (the footprint). Sortedness is checked on `(chrom, key)`. |
| D2 | `sort` orders VCF/BCF by the `ref` span (`POS`, then ref end), whatever `--vcf-span` says, so sorted VCF is always tabix-compatible. |
| D3 | The atomic temp file always sits next to `PATH` so the rename stays on one filesystem. `--tmp-dir` controls only `sort` spill files. |
| D4 | `-f/-F` fractions are parsed as exact decimals and compared with integer arithmetic, so `0.3` of 10 bases is exactly 3. |
| D5 | `--min-bases` counts as met when either side is a point, the same as the point rule for fractions. |
| D6 | Under `--signed a` or `--signed b`, an unknown strand is treated as `+`. |
| D7 | A chromosome missing from `--genome` is an error. When a VCF header declares contigs, a record on an undeclared contig is an error. A header with no contigs falls through to lexicographic order. |
| D8 | Region strings accept `chr`, `chr:s`, `chr:s-`, `chr:s-e` and `{name}:s-e`, with commas allowed in numbers. If the whole string is a known chromosome name, it means that chromosome. Otherwise it is split at the last `:`. If both readings are valid names, it is an error (samtools behaviour). |
| D9 | `view` region queries select records whose footprint (per `--vcf-span`) overlaps the region under §4.3, so an insertion on a region edge is excluded by default. `--vcf-span ref --inclusive-points` reproduces tabix. This goes in `docs/differences.md`. |
| D10 | `view` sorts and merges the requested regions, so each record is printed once, in file order. |
| D11 | For uncompressed BED without an index, `view` binary-searches to the first line of the region's chromosome, then scans that chromosome, so long features that start before the region are found. If the search can't find the chromosome, it falls back to a full scan with a warning, so a file sorted in another chromosome order still gives correct results. |
| D12 | `sort --check` checks `(chrom, start)` order, which is what streaming commands need. |
| D13 | `--unique` and `--duplicates` imply `--full-key` and behave like `uniq` and `uniq -d` on whole lines. |
| D14 | BED header lines (`#`, `track`, `browser`) before the first record are kept only with `--header`, and such lines after the first record are dropped. Blank lines are skipped. In GFF, leading `#`/`##` lines form the header, `###` lines are dropped, and `##FASTA` is an error in P0. |
| D15 | CRLF input is accepted, and output always uses LF. |
| D16 | bedGraph is detected only from a `track type=bedGraph` line or `--format bedgraph`; anything else BED-shaped is BED. |
| D17 | Output format comes from the `-o` extension (case-insensitive), or from `-O` for stdout. If both are given they must agree. Default: the input's format, uncompressed, on stdout. BCF or BGZF output to a terminal is a usage error. |
| D18 | `sort` doesn't convert formats (the `-o` family must match the input). VCF↔BCF conversion is `view`'s job. Other conversions are `convert` (P1). |
| D19 | When gloci replaces `PATH`, it deletes any existing `PATH.tbi`/`PATH.csi` unless it writes a fresh index. Stale indexes are never left behind. |
| D20 | Indexes are built in a second pass over the finished temp file. TBI is used unless a record's end exceeds 2^29, in which case CSI is used (min_shift 14, depth chosen to cover the largest end). BCF always gets CSI. Index spans use the ref span, extended by INFO/END when present (a superset of htslib's span, so htslib queries stay correct). A zero-length BED record is indexed as `[s, s+1)`. |
| D21 | `--compression-level` accepts 0–12 (libdeflate's range). |
| D22 | BCF fast path: hand-written record framing and a minimal typed-value decoder (CHROM, POS, rlen, alleles, INFO/END, INFO/SVLEN), verified against noodles. BCF↔VCF conversion uses noodles. VCF headers are parsed by hand for `##contig` and `#CHROM` only, so real-world headers that noodles rejects still work. |
| D23 | `--max-mem` takes binary units (`K`, `M`, `G`, `T` = powers of 1024), with a minimum of 16M. |
| D24 | Spill files are LZ4 frames (`lz4_flex`) in a temp directory under `--tmp-dir`, then `GLOCI_TMPDIR`, then the system temp dir. Merge fan-in is 128. Chunk boundaries depend only on `--max-mem` and the input, never on threads. |
| D25 | In `check`, duplicates are warnings (exit 0). Parse, order and bounds problems are errors (exit 1). |
| D26 | An index older than its data file triggers a warning (htslib behaviour). |
| D27 | `--genome` accepts a BAM header now, since it is a small binary parse. |
| D28 | VCF `POS 0` is invalid. REF/ALT comparison is ASCII case-insensitive. `END` is used only for records with a symbolic ALT. |
| D29 | Multiple VCF inputs to `sort` must have identical sample columns; the header comes from the first file. |
| D30 | A seekable BGZF file without the 28-byte EOF marker is an error ("probably truncated"). stdin can't be checked. |
| D31 | Being killed by SIGINT, SIGTERM or SIGHUP removes registered temp files and spill dirs, then exits 130. A panic hook does the same before abort. |

## Review Focus

The spec implies these inputs but no main-path test exercises them. They are the ones most likely to bite a user. Each has a test in the owning task.

1. **Empty inputs.** A zero-record BED, a header-only VCF/BCF, empty stdin, or an empty BGZF file through `sort`, `view`, `check` and `index` should give valid output (headers only, BGZF EOF block, a loadable index) and exit 0. Tests in Tasks 12, 15, 16 and 18.
2. **Job killed mid-write** (SLURM time limit, Ctrl-C). Neither `PATH` nor `PATH.tmp-*` should remain, and the sort spill dir should be gone. Test in Task 22.
3. **Rewriting outputs in place.** `gloci sort x.bed -o x.bed` should work, and replacing an indexed `x.bed.gz` without `--write-index` should leave no stale `.tbi`/`.csi`. Tests in Tasks 7 and 15.
4. **Unusual chromosome names.** Names containing `:` (GRCh38 `HLA-A*01:01:01:01`) need correct region parsing, and assemblies with 100k+ contigs should sort correctly in reasonable time. Tests in Tasks 2, 15 and 17.
5. **Very long records.** Multi-MB VCF lines (100k+ samples) under a small `--max-mem`, and records that straddle BGZF block boundaries, should sort and read intact. Tests in Tasks 8 and 14.

---

## File Structure

```
Cargo.toml                         workspace, [workspace.dependencies], release profile
rust-toolchain.toml                channel = "stable", components rustfmt, clippy
.gitignore
.github/workflows/ci.yml           fmt, clippy, test (3 platforms), static musl build
.github/workflows/nightly.yml      large-file tests, bench, differential, fuzz
crates/gloci-core/src/
  lib.rs        re-exports
  error.rs      CoordError, RegionError, GenomeError, UnknownChrom, FootprintError, FractionError
  fmt.rs        thousands()
  interval.rs   Interval, Strand, parse_coord, one_based_closed
  chrom.rs      ChromId, ChromTable, OrderSource
  genome.rs     Genome (+ .fai/.genome/.dict text parsers)
  region.rs     Region parsing (D8)
  overlap.rs    overlaps, overlap_len, Fraction, OverlapCriteria, StrandRule, Window
  distance.rs   distance, signed_distance, SignFrame
  footprint.rs  VcfSpan, Variant, footprint, alt_footprints
  order.rs      OrderChecker, OrderViolation
crates/gloci-core/tests/           coords.rs, overlap.rs, distance.rs, footprint.rs, order.rs
crates/gloci-io/src/
  lib.rs, error.rs (IoError, Loc, OrderError, UsageError)
  detect.rs     Compression, Format, detect_compression, detect_format
  input.rs      open_input (stdin, gzip, parallel BGZF, EOF check)
  atomic.rs     AtomicFile, cleanup registry, signal + panic handlers
  output.rs     OutFormat, OutputSpec, parse_output_spec, Output, open_output
  record.rs     RecordMeta, RecordBatch, Header, RecordReader, ReaderOptions, open_reader(_from)
  diag.rs       InvalidPolicy, Diagnostics
  bed.rs        BED/bedGraph reader, bed12_blocks
  gff.rs        GFF3/GTF reader
  vcf.rs        VcfHeader, VCF text reader
  bcf.rs        BCF framing reader, BcfDicts, typed-value decoder
  genome.rs     load_genome, chrom_table_for
  sorted.rs     OrderGuard
  writer.rs     RecordWriter, open_record_writer, WriterOptions
  convert.rs    convert_variants (VCF<->BCF via noodles)
  index.rs      IndexKind, IndexInfo, build_index, find_index, open_region_reader
crates/gloci-ops/src/
  lib.rs
  sort/mod.rs, sort/chunk.rs, sort/spill.rs, sort/merge.rs
  view.rs, bedsearch.rs, check.rs
crates/gloci-testkit/src/          lib.rs, fixtures.rs, strategies.rs, naive.rs
crates/gloci/src/
  main.rs       allocator, run(), exit-code mapping
  lib.rs        pub fn command() -> clap::Command; pub fn run(...)
  cli.rs        clap structs, shared arg groups, parse_size
  logging.rs    stderr logger for -q/-v
  commands/     mod.rs, sort.rs, view.rs, check.rs, index.rs
crates/gloci/tests/                common/mod.rs, cli_basics.rs, index.rs, sort.rs, view.rs,
                                   view_regions.rs, check.rs, behaviour.rs, large.rs
xtask/src/                         main.rs, gen.rs, bench.rs, diff.rs, docs.rs
xtask/bench/cases.toml, xtask/diff/cases.toml, xtask/docker/Dockerfile
fuzz/                              Cargo.toml, fuzz_targets/*.rs
docs/differences.md
```

---

### Task 1: Workspace, binary skeleton, CI

**Files:**
- Create: `Cargo.toml`, `rust-toolchain.toml`, `.gitignore`, `.cargo/config.toml`, `crates/{gloci-core,gloci-io,gloci-ops,gloci-testkit,gloci}/Cargo.toml` + `src/lib.rs`, `crates/gloci/src/main.rs`, `xtask/Cargo.toml`, `xtask/src/main.rs`, `.github/workflows/ci.yml`
- Test: `crates/gloci/tests/cli_basics.rs`

**Interfaces:**
- Produces: workspace members `crates/*` and `xtask` (`fuzz/` excluded). `[workspace.package] version = "0.1.0"`, `edition = "2024"`, `rust-version = "1.85"`, `publish = false` (the license is an open question, spec §13). All third-party versions from the Tech Stack line go in `[workspace.dependencies]`. Crate deps: core ← io ← ops ← gloci; testkit is dev-only. `cargo xtask` alias in `.cargo/config.toml`.

- [ ] **Step 1: `git init`, then write the failing test**

```rust
// crates/gloci/tests/cli_basics.rs
use assert_cmd::Command;
#[test] fn version_prints_name_and_version() {
    Command::cargo_bin("gloci").unwrap().arg("--version").assert().success().stdout("gloci 0.1.0\n");
}
#[test] fn no_args_is_usage_error() {
    Command::cargo_bin("gloci").unwrap().assert().code(2);
}
```

- [ ] **Step 2: Run** `cargo test -p gloci --test cli_basics`. Expected: FAIL (no binary).
- [ ] **Step 3: Implement.** Create the workspace manifests. `main.rs` sets `#[global_allocator] static GLOBAL: mimalloc::MiMalloc` and parses a clap `Cli` with `#[command(name = "gloci", version, arg_required_else_help = true)]`. Add the release profile from Global Constraints. In `ci.yml`, a matrix over `ubuntu-24.04`, `ubuntu-24.04-arm` and `macos-14` runs `cargo fmt --check`, `cargo clippy --workspace --all-targets -- -D warnings` and `cargo test --workspace`. A `musl` job on ubuntu installs `musl-tools`, builds `--release --target x86_64-unknown-linux-musl -p gloci`, and asserts that `ldd` reports "not a dynamic executable".
- [ ] **Step 4: Run** `cargo test --workspace && cargo clippy --workspace --all-targets -- -D warnings`. Expected: PASS, no warnings.
- [ ] **Step 5: Commit** `git add -A && git commit -m "chore: workspace skeleton, gloci binary, CI"`

---

### Task 2: Core coordinates, chromosomes, genome, regions

**Files:**
- Create: `crates/gloci-core/src/{error.rs,fmt.rs,interval.rs,chrom.rs,genome.rs,region.rs}`
- Test: `crates/gloci-core/tests/coords.rs`

**Interfaces:**
- Produces:
  - `pub struct Interval { pub start: u64, pub end: u64 }` (Copy, Eq, Ord, Hash, Debug), with `Interval::new(start, end) -> Result<Interval, CoordError>`, `Interval::point(p) -> Interval`, `.is_point() -> bool` and `.len() -> u64`.
  - `pub enum Strand { Plus, Minus, Unknown }` with `Strand::parse(&[u8]) -> Option<Strand>` (`+`, `-`, `.`; GFF `?` maps to Unknown).
  - `pub fn parse_coord(field: &[u8]) -> Result<u64, CoordError>` and `pub fn one_based_closed(s: u64, e: u64) -> Result<Interval, CoordError>` (gives `[s-1, e)` and requires `s >= 1`, `s <= e`).
  - `CoordError` Display strings: `"{v}" is not an integer`, `"{v}" is negative`, `"{v}" is too large`, `start {s} is greater than end {e}`, `start must be >= 1 (1-based coordinates)`.
  - `pub fn thousands(n: u64) -> String`.
  - `pub struct ChromId(pub u32)`; `pub enum OrderSource { GenomeFile(PathBuf), VcfHeader(String), Lexicographic, FirstAppearance }`.
  - `ChromTable` constructors: `::lexicographic()`, `::first_appearance()` and `::closed(&Genome, OrderSource)`. Methods: `.intern(&[u8]) -> Result<ChromId, UnknownChrom>`, `.get(&[u8]) -> Option<ChromId>`, `.name(ChromId) -> &[u8]`, `.len()`, `.cmp(a, b) -> Ordering`, `.ranks() -> Vec<u32>` (dense rank per id in the table's order) and `.order_source() -> &OrderSource`.
  - `Genome::new(Vec<(Vec<u8>, Option<u64>)>) -> Result<Genome, GenomeError>` (duplicate names are an error), `Genome::parse_sizes(&[u8])` (`.fai`/`.genome`: name TAB length …), `Genome::parse_dict(&[u8])` (`@SQ SN: LN:`), plus `.len()`, `.name(i)` and `.length(i) -> Option<u64>`.
  - `pub struct Region { pub chrom: Vec<u8>, pub start: u64, pub end: Option<u64> }` (0-based; `end: None` means to the chromosome end), with `Region::parse(s: &str, is_known: impl Fn(&[u8]) -> bool) -> Result<Region, RegionError>` (D8) and `.interval(chrom_len: Option<u64>) -> Interval`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-core/tests/coords.rs (excerpt; one test fn per group)
assert_eq!(parse_coord(b"007").unwrap(), 7);
assert_eq!(parse_coord(b"12a4").unwrap_err().to_string(), r#""12a4" is not an integer"#);
assert_eq!(parse_coord(b"-1").unwrap_err().to_string(), r#""-1" is negative"#);
assert!(parse_coord(b"18446744073709551616").is_err() && parse_coord(b"").is_err() && parse_coord(b"+5").is_err());
assert!(Interval::new(10, 5).is_err());
assert!(Interval::new(10, 10).unwrap().is_point());
assert_eq!(one_based_closed(100, 200).unwrap(), Interval::new(99, 200).unwrap());   // GFF, region strings
assert!(one_based_closed(0, 5).is_err());
assert_eq!(thousands(1204331), "1,204,331");
assert_eq!(Strand::parse(b"?"), Some(Strand::Unknown));
// regions (D8)
let none = |_: &[u8]| false;
assert_eq!(Region::parse("chr1:1,000-2,000", none).unwrap(), Region { chrom: b"chr1".to_vec(), start: 999, end: Some(2000) });
assert_eq!(Region::parse("chr1:100-100", none).unwrap().interval(None), Interval::new(99, 100).unwrap());
assert_eq!(Region::parse("chr1", none).unwrap().end, None);
assert_eq!(Region::parse("chr1:500", none).unwrap(), Region { chrom: b"chr1".to_vec(), start: 499, end: None });
assert_eq!(Region::parse("chr1:500-", none).unwrap().start, 499);
assert!(Region::parse("chr1:0-10", none).is_err() && Region::parse("chr1:20-10", none).is_err());
let hla = |n: &[u8]| n == b"HLA-A*01:01:01:01";
assert_eq!(Region::parse("HLA-A*01:01:01:01", hla).unwrap().chrom, b"HLA-A*01:01:01:01");
assert_eq!(Region::parse("{HLA-A*01:01:01:01}:5-10", none).unwrap().start, 4);
let both = |n: &[u8]| n == b"x:1" || n == b"x";
assert!(Region::parse("x:1", both).is_err());                       // ambiguous
// chromosome tables
let mut t = ChromTable::lexicographic();
let (c2, c10, c1) = (t.intern(b"chr2").unwrap(), t.intern(b"chr10").unwrap(), t.intern(b"chr1").unwrap());
assert_eq!(t.ranks(), vec![2, 1, 0]);
assert_eq!(t.cmp(c10, c2), Ordering::Less);
let g = Genome::parse_sizes(b"chr1\t100\nchr2\t200\nchr10\t50\n").unwrap();
let mut t = ChromTable::closed(&g, OrderSource::Lexicographic);
assert!(t.intern(b"chrUn").is_err());
assert_eq!(t.cmp(t.get(b"chr2").unwrap(), t.get(b"chr10").unwrap()), Ordering::Less);
assert!(Genome::parse_sizes(b"chr1\t1\nchr1\t2\n").is_err());       // duplicate
assert_eq!(Genome::parse_dict(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n").unwrap().length(0), Some(248956422));
// 100k contigs intern + ranks stays fast (Review Focus 4)
let mut t = ChromTable::lexicographic();
for i in 0..100_000 { t.intern(format!("scaffold_{i}").as_bytes()).unwrap(); }
assert_eq!(t.ranks().len(), 100_000);
```

- [ ] **Step 2: Run** `cargo test -p gloci-core --test coords`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** the files above. Use thiserror for the error enums. `ChromTable` uses a `HashMap<Box<[u8]>, ChromId>` plus a `Vec` of names. In `closed`, ids `0..n` follow genome order, so `cmp` compares ids. `lexicographic` compares name bytes. `first_appearance` compares ids, which are assigned in order of first appearance.
- [ ] **Step 4: Run** `cargo test -p gloci-core`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(core): intervals, coordinates, chromosome tables, genome, regions"`

---

### Task 3: Core overlap, fractions, windows, distance (§4.3, §4.4)

**Files:**
- Create: `crates/gloci-core/src/{overlap.rs,distance.rs}`
- Test: `crates/gloci-core/tests/{overlap.rs,distance.rs}`

**Interfaces:**
- Consumes: `Interval`, `Strand` (Task 2).
- Produces:
  - `pub fn overlaps(a: Interval, b: Interval, inclusive_points: bool) -> bool` and `pub fn overlap_len(a, b) -> u64` (0 if either is a point).
  - `Fraction::parse(&str) -> Result<Fraction, FractionError>`: a plain decimal in (0, 1], at most 18 fractional digits, stored as `n / 10^k`. `Fraction::met(&self, part: u64, whole: u64) -> bool` is exact, using u128 (D4).
  - `pub enum StrandRule { Any, Same, Opposite }`.
  - `pub struct OverlapCriteria { pub min_frac_a: Option<Fraction>, pub min_frac_b: Option<Fraction>, pub reciprocal: bool, pub either: bool, pub min_bases: Option<u64>, pub strand: StrandRule, pub inclusive_points: bool }` with `Default`, and `.accepts(a, a_strand, b, b_strand) -> bool`.
  - `pub struct Window { pub left: u64, pub right: u64, pub strand_aware: bool }` with `.apply(b: Interval, a_strand: Strand) -> Interval`. It pads B, saturating at 0, and swaps left/right when `strand_aware` and A is on `Minus`.
  - `pub fn distance(a, b, inclusive_points) -> u64` and `pub enum SignFrame { Ref, A, B }`.
  - `pub fn signed_distance(a, a_strand, b, b_strand, frame, inclusive_points) -> i64`.

`accepts` semantics: `overlaps && frac_ok && bases_ok && strand_ok`, where:
- `fa = min_frac_a.map_or(true, |f| a.is_point() || f.met(ov, a.len()))`.
- `fb` is the same on B; when `reciprocal` and `min_frac_b` is None, B uses `min_frac_a`.
- `frac_ok = if either { fa || fb } else { fa && fb }` (with `either` and only one threshold, it reduces to that threshold).
- `bases_ok = min_bases.map_or(true, |n| a.is_point() || b.is_point() || ov >= n)` (D5).
- `Same` and `Opposite` require both strands to be known.

`signed_distance` signs (bedtools `-D`):
- `Ref`: negative when B lies at lower coordinates than A.
- `A`: same as `Ref`, flipped when A is `Minus`.
- `B`: negative when A lies at lower coordinates than B, flipped when B is `Minus`.
- Unknown strand counts as `+` (D6). Overlap gives 0 in every frame.

- [ ] **Step 1: Write the failing tests.** Include the spec §4.3 worked-example table verbatim.

```rust
// crates/gloci-core/tests/overlap.rs
fn iv(s: u64, e: u64) -> Interval { Interval::new(s, e).unwrap() }
#[test] fn spec_4_3_worked_examples() {
    // (a, b, default, with --inclusive-points)
    let cases = [
        (iv(10,10), iv(5,15),  true,  true),   // insertion site inside a feature IS found
        (iv(5,5),   iv(5,15),  false, true),
        (iv(15,15), iv(5,15),  false, true),
        (iv(10,10), iv(10,10), true,  true),
        (iv(10,10), iv(11,11), false, false),
        (iv(5,10),  iv(10,15), false, false),  // bookended
        (iv(5,11),  iv(10,15), true,  true),
    ];
    for (a, b, dflt, incl) in cases {
        assert_eq!(overlaps(a, b, false), dflt, "{a:?} {b:?}");
        assert_eq!(overlaps(b, a, false), dflt, "symmetric {a:?} {b:?}");
        assert_eq!(overlaps(a, b, true), incl, "inclusive {a:?} {b:?}");
    }
}
#[test] fn fractions_are_exact() {
    assert!(Fraction::parse("0.3").unwrap().met(3, 10));
    assert!(!Fraction::parse("0.3").unwrap().met(2, 10));
    assert!(Fraction::parse("0.1").unwrap().met(1, 10));
    assert!(Fraction::parse("1").unwrap().met(10, 10) && Fraction::parse("1.0").is_ok());
    assert!(Fraction::parse("0.5").unwrap().met(u64::MAX / 2 + 1, u64::MAX));
    for bad in ["0", "1.5", "-0.1", "abc", "5e-1", ""] { assert!(Fraction::parse(bad).is_err(), "{bad}"); }
}
#[test] fn criteria() {
    let (a, b) = (iv(0, 10), iv(5, 20));                       // overlap 5: 0.5 of A, 1/3 of B
    let f = |s| Some(Fraction::parse(s).unwrap());
    let c = |m: OverlapCriteria| m.accepts(a, Strand::Plus, b, Strand::Plus);
    assert!(c(OverlapCriteria { min_frac_a: f("0.5"), ..Default::default() }));
    assert!(!c(OverlapCriteria { min_frac_a: f("0.5"), reciprocal: true, ..Default::default() }));
    assert!(c(OverlapCriteria { min_frac_a: f("0.9"), min_frac_b: f("0.3"), either: true, ..Default::default() }));
    assert!(!c(OverlapCriteria { min_bases: Some(6), ..Default::default() }));
    let point_ok = OverlapCriteria { min_frac_a: f("0.9"), min_bases: Some(6), ..Default::default() };
    assert!(point_ok.accepts(Interval::point(10), Strand::Plus, b, Strand::Plus));
    let same = OverlapCriteria { strand: StrandRule::Same, ..Default::default() };
    assert!(!same.accepts(a, Strand::Plus, b, Strand::Minus) && !same.accepts(a, Strand::Unknown, b, Strand::Unknown));
    let opp = OverlapCriteria { strand: StrandRule::Opposite, ..Default::default() };
    assert!(opp.accepts(a, Strand::Plus, b, Strand::Minus) && !opp.accepts(a, Strand::Unknown, b, Strand::Minus));
}
#[test] fn windows() {
    assert_eq!(Window { left: 5, right: 5, strand_aware: false }.apply(Interval::point(10), Strand::Plus), iv(5, 15));
    assert_eq!(Window { left: 50, right: 0, strand_aware: false }.apply(iv(10, 20), Strand::Plus), iv(0, 20));
    assert_eq!(Window { left: 5, right: 0, strand_aware: true }.apply(iv(10, 20), Strand::Minus), iv(10, 25));
}
```

```rust
// crates/gloci-core/tests/distance.rs
assert_eq!(distance(iv(5,15), iv(10,20), false), 0);
assert_eq!(distance(iv(5,10), iv(10,15), false), 1);          // bookended
assert_eq!(distance(iv(0,10), iv(15,20), false), 6);          // gap 5
assert_eq!(distance(Interval::point(10), iv(5,15), false), 0);
assert_eq!(distance(Interval::point(5), iv(5,15), false), 1);
assert_eq!(distance(Interval::point(5), iv(5,15), true), 0);
assert_eq!(distance(Interval::point(10), Interval::point(12), false), 3);
let (a, up, down) = (iv(100,200), iv(50,60), iv(250,260));
use SignFrame::*; use Strand::*;
assert_eq!(signed_distance(a, Plus, up, Plus, Ref, false), -41);
assert_eq!(signed_distance(a, Minus, up, Plus, A, false), 41);
assert_eq!(signed_distance(a, Plus, up, Plus, B, false), 41);
assert_eq!(signed_distance(a, Plus, up, Minus, B, false), -41);
assert_eq!(signed_distance(a, Plus, down, Plus, Ref, false), 51);
assert_eq!(signed_distance(a, Minus, down, Plus, A, false), -51);
assert_eq!(signed_distance(a, Plus, down, Plus, B, false), -51);
assert_eq!(signed_distance(a, Unknown, up, Plus, A, false), -41);   // D6
assert_eq!(signed_distance(a, Plus, iv(150,160), Minus, B, false), 0);
```

- [ ] **Step 2: Run** `cargo test -p gloci-core --test overlap --test distance`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `overlap.rs` and `distance.rs` per the Interfaces block. A point is the interval `[p,p)`, so for non-overlapping inputs `gap = max(b.start - a.end, a.start - b.end)` with saturating subtraction, and the distance is `gap + 1`.
- [ ] **Step 4: Run** `cargo test -p gloci-core`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(core): overlap, exact fractions, windows, distance (spec 4.3, 4.4)"`

---

### Task 4: Core variant footprints (§4.2)

**Files:**
- Create: `crates/gloci-core/src/footprint.rs`
- Test: `crates/gloci-core/tests/footprint.rs`

**Interfaces:**
- Consumes: `Interval` (Task 2).
- Produces:
  - `pub enum VcfSpan { Allele, Ref }` (Default = Allele).
  - `pub struct Variant<'a> { pub pos: u64, pub reference: &'a [u8], pub alts: &'a [u8], pub end: Option<u64>, pub svlen: Option<&'a [u8]> }`. `alts` is the raw ALT column; `svlen` is the raw INFO/SVLEN value, which may be comma-separated, one per ALT, or a single value for all ALTs.
  - `pub fn footprint(v: &Variant, span: VcfSpan) -> Result<Interval, FootprintError>`.
  - `pub fn alt_footprints(v: &Variant) -> Result<Vec<Option<Interval>>, FootprintError>` returns the allele footprint per ALT, with None for `*`, `<*>` and `<NON_REF>`.
  - `FootprintError` variants: `PosZero`, `EmptyRef`, `MissingSvEnd` (Display: `symbolic ALT {alt} needs INFO/END or INFO/SVLEN`), `EndBeforePos`, `BadSvlen`, `SvlenCount`.

Rules beyond the spec text:
- Trim the shared prefix first, then the shared suffix, comparing ASCII case-insensitively.
- A symbolic ALT is `<…>`, classified by its text before the first `:`. `INS` gives a point at `POS`. `*` and `NON_REF` give a gVCF block `[POS-1, END or ref end)`. Anything else (`DEL`, `DUP`, `INV`, `CNV`, `CN0`, …) is an SV ALT and gives `[POS, E)` in allele mode, where `E` is END, else `POS + |SVLEN|`. An SV ALT with neither gives `MissingSvEnd` in both modes.
- A breakend is any ALT containing `[`, `]`, or a leading/trailing `.` with bases. It gives `[POS-1, POS-1+len(REF))`.
- END is used only when a symbolic ALT is present (D28). `END < POS` gives `EndBeforePos`.
- `ref` mode is `[POS-1, max(POS-1+len(REF), E))`, where `E` is END, or the SVLEN-derived end for SV ALTs, when the record has a symbolic ALT. Otherwise it is `[POS-1, POS-1+len(REF))`.

- [ ] **Step 1: Write the failing tests.** Include the spec §4.2 table verbatim and the further rules.

```rust
// crates/gloci-core/tests/footprint.rs
fn iv(s: u64, e: u64) -> Interval { Interval::new(s, e).unwrap() }
fn fp(pos: u64, r: &str, a: &str, end: Option<u64>, svlen: Option<&str>, span: VcfSpan) -> Result<Interval, FootprintError> {
    footprint(&Variant { pos, reference: r.as_bytes(), alts: a.as_bytes(), end, svlen: svlen.map(str::as_bytes) }, span)
}
#[test] fn spec_4_2_table() {
    let p = Interval::point;
    let cases = [
        (100, "A", "G", None, iv(99,100), iv(99,100)),
        (100, "AC", "GT", None, iv(99,101), iv(99,101)),
        (100, "ACGT", "A", None, iv(99,103), iv(100,103)),
        (100, "A", "ACGT", None, iv(99,100), p(100)),
        (100, "ATG", "AC", None, iv(99,102), iv(100,102)),
        (100, "ACGT", "A,AC", None, iv(99,103), iv(100,103)),
        (100, "N", "<DEL>", Some(500), iv(99,500), iv(100,500)),
        (100, "N", "<INS>", None, iv(99,100), p(100)),
        (100, "N", "<DUP>", Some(700), iv(99,700), iv(100,700)),
        (100, "N", "<INV>", Some(700), iv(99,700), iv(100,700)),
        (100, "N", "<CNV>", Some(700), iv(99,700), iv(100,700)),
        (100, "G", "G]17:198982]", None, iv(99,100), iv(99,100)),
        (100, "A", "<*>", Some(250), iv(99,250), iv(99,250)),
        (100, "A", "<NON_REF>", Some(250), iv(99,250), iv(99,250)),
        (100, "A", ".", None, iv(99,100), iv(99,100)),
    ];
    for (pos, r, a, end, want_ref, want_allele) in cases {
        assert_eq!(fp(pos, r, a, end, None, VcfSpan::Ref).unwrap(), want_ref, "ref {r}>{a}");
        assert_eq!(fp(pos, r, a, end, None, VcfSpan::Allele).unwrap(), want_allele, "allele {r}>{a}");
    }
}
#[test] fn further_rules() {
    use VcfSpan::*;
    assert_eq!(fp(100, "N", "<DEL>", None, Some("-400"), Allele).unwrap(), iv(100, 500));     // SVLEN fallback
    assert_eq!(fp(100, "N", "<DEL>", Some(300), Some("-400"), Allele).unwrap(), iv(100, 300)); // END wins
    assert_eq!(fp(100, "N", "<DEL>,<DEL>", None, Some("-100,-300"), Allele).unwrap(), iv(100, 400));
    assert!(matches!(fp(100, "N", "<DEL>", None, None, Allele), Err(FootprintError::MissingSvEnd { .. })));
    assert!(matches!(fp(100, "N", "<DEL>", Some(50), None, Allele), Err(FootprintError::EndBeforePos { .. })));
    assert_eq!(fp(100, "ACGT", "A,*", None, None, Allele).unwrap(), iv(100, 103));          // * ignored
    assert_eq!(fp(100, "A", "G,<NON_REF>", Some(900), None, Allele).unwrap(), iv(99, 100)); // NON_REF ignored
    assert_eq!(fp(100, "A", "*", None, None, Allele).unwrap(), iv(99, 100));                // * alone
    assert_eq!(fp(100, "acgt", "A", None, None, Allele).unwrap(), iv(100, 103));            // case-insensitive
    assert_eq!(fp(100, "AA", "A", None, None, Allele).unwrap(), iv(100, 101));              // prefix first
    assert_eq!(fp(100, "N", "<INS:ME:ALU>", Some(100), None, Allele).unwrap(), Interval::point(100));
    assert_eq!(fp(100, "N", "<DEL:ME>", Some(400), None, Allele).unwrap(), iv(100, 400));
    assert_eq!(fp(100, "G", "G.", None, None, Allele).unwrap(), iv(99, 100));               // single breakend
    assert_eq!(fp(100, "ACGTT", "A", Some(999), None, Ref).unwrap(), iv(99, 104));          // END ignored, not symbolic
    assert!(matches!(fp(0, "A", "G", None, None, Allele), Err(FootprintError::PosZero)));
    let per_alt = alt_footprints(&Variant { pos: 100, reference: b"ACGT", alts: b"A,AC,*", end: None, svlen: None }).unwrap();
    assert_eq!(per_alt, vec![Some(iv(100, 103)), Some(iv(101, 103)), None]);
}
```

- [ ] **Step 2: Run** `cargo test -p gloci-core --test footprint`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `footprint.rs`. Split ALTs on `,` with `memchr`. The record footprint in allele mode is the smallest interval covering all `Some` per-ALT footprints. If there are none (only `*`, `.` or gVCF ALTs), fall back to the ref/gVCF span.
- [ ] **Step 4: Run** `cargo test -p gloci-core`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(core): allele-aware VCF footprints (spec 4.2)"`

---

### Task 5: Core chromosome-order checking (§4.5)

**Files:**
- Create: `crates/gloci-core/src/order.rs`
- Test: `crates/gloci-core/tests/order.rs`

**Interfaces:**
- Consumes: `ChromTable`, `ChromId` (Task 2).
- Produces:
  - `pub struct OrderChecker` with `::new()` and `.check(&mut self, table: &ChromTable, chrom: ChromId, key: u64) -> Result<(), OrderViolation>`. After a violation it adopts the new record as "previous" so callers like `check` can keep scanning.
  - `pub enum OrderViolation { ChromBackwards { prev: ChromId, next: ChromId }, StartBackwards { chrom: ChromId, prev: u64, next: u64 } }`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-core/tests/order.rs
#[test] fn lexicographic() {
    let mut t = ChromTable::lexicographic(); let mut c = OrderChecker::new();
    for (n, k) in [("chr1", 5), ("chr1", 5), ("chr10", 0), ("chr2", 3)] { let id = t.intern(n.as_bytes()).unwrap(); c.check(&t, id, k).unwrap(); }
    let mut c = OrderChecker::new();
    c.check(&t, t.get(b"chr2").unwrap(), 0).unwrap();
    assert!(matches!(c.check(&t, t.get(b"chr10").unwrap(), 0), Err(OrderViolation::ChromBackwards { .. })));
}
#[test] fn genome_order_and_reappearance() {
    let g = Genome::parse_sizes(b"chr1\t1000\nchr2\t1000\nchr10\t1000\n").unwrap();
    let t = ChromTable::closed(&g, OrderSource::Lexicographic);
    let id = |n: &[u8]| t.get(n).unwrap();
    let mut c = OrderChecker::new();
    c.check(&t, id(b"chr2"), 0).unwrap(); c.check(&t, id(b"chr10"), 0).unwrap();
    assert!(c.check(&t, id(b"chr2"), 50).is_err());                     // chr2 reappears
    let mut c = OrderChecker::new();
    c.check(&t, id(b"chr1"), 100).unwrap();
    assert!(matches!(c.check(&t, id(b"chr1"), 99), Err(OrderViolation::StartBackwards { prev: 100, next: 99, .. })));
    c.check(&t, id(b"chr1"), 120).unwrap();                             // continues after violation
}
#[test] fn first_appearance() {
    let mut t = ChromTable::first_appearance(); let mut c = OrderChecker::new();
    for n in ["chr2", "chr1", "chr10"] { let id = t.intern(n.as_bytes()).unwrap(); c.check(&t, id, 0).unwrap(); }
    assert!(c.check(&t, t.get(b"chr1").unwrap(), 0).is_err());
}
```

- [ ] **Step 2: Run** `cargo test -p gloci-core --test order`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `order.rs`. Compare chromosomes with `table.cmp`. End order is ignored (spec §4.5).
- [ ] **Step 4: Run** `cargo test -p gloci-core`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(core): on-the-fly sort order checker (spec 4.5)"`

---

### Task 6: IO detection and input opening

**Files:**
- Create: `crates/gloci-io/src/{error.rs,detect.rs,input.rs}`, `crates/gloci-testkit/src/fixtures.rs`
- Test: `crates/gloci-io/tests/input.rs`

**Interfaces:**
- Produces:
  - `pub enum Compression { None, Gzip, Bgzf }` and `pub enum Format { Bed, BedGraph, Vcf, Bcf, Gff3, Gtf }` with `FromStr` (`bed`, `bedgraph`, `vcf`, `bcf`, `gff`, `gff3`, `gtf`).
  - `pub fn detect_compression(head: &[u8]) -> Compression`. BGZF means gzip magic plus FEXTRA with a `BC` subfield.
  - `pub fn detect_format(head: &[u8]) -> Result<Format, String>`, applied to decompressed bytes.
  - `pub struct InputOptions { pub threads: NonZeroUsize, pub format: Option<Format> }` and `pub struct Input { pub name: String, pub format: Format, pub compression: Compression, pub reader: Box<dyn BufRead + Send> }`.
  - `pub fn open_input(path: &str, opts: &InputOptions) -> Result<Input, IoError>`, plus `pub fn open_input_from(name: &str, src: Box<dyn Read + Send>, opts) -> Result<Input, IoError>`. `"-"` opens stdin with the display name `stdin`.
  - `IoError` (thiserror) variants: `Open { path, source }`, `Read { path, source }`, `Parse { path, loc: Loc, msg }`, `Detect { path, msg }`, `Truncated { path }`, `Order(OrderError)` (Task 11) and `Usage(UsageError)`.
  - `pub enum Loc { Line(u64), Record(u64) }`. Display: `"{path} line {thousands(n)}: {msg}"` or `"{path} record {thousands(n)}: {msg}"`.
  - testkit: `fixtures::bgzip(data: &[u8]) -> Vec<u8>` (via noodles-bgzf, with EOF block), `fixtures::gzip(data) -> Vec<u8>` and `fixtures::write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf`.

Detection rules:
- `BCF\x02` → Bcf.
- `##fileformat=VCF` → Vcf.
- `##gff-version 3` → Gff3; `##gff-version 2` → Gtf.
- Otherwise, take the first line that isn't `#`, `track`, `browser` or blank. With ≥9 tab columns where columns 4 and 5 are integers and column 2 isn't, it is Gtf if column 9 contains `"`, else Gff3. With ≥3 columns where columns 2 and 3 are integers, it is Bed, or BedGraph if a preceding `track` line has `type=bedGraph`.
- An empty input is Bed.
- Anything else is `Detect` with the message `could not recognise the format of {path} (first data line: {line:.80}); pass --format`.

BGZF uses `noodles_bgzf::io::MultithreadedReader::with_worker_count(threads, file)`. Plain gzip uses `flate2::read::MultiGzDecoder`. Seekable BGZF files are checked for the 28-byte EOF block before reading (D30).

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/input.rs: one #[test] per line group
// detection is by content: each fixture is written under a misleading name ("x.txt")
for (bytes, fmt) in [(b"chr1\t1\t2\n".to_vec(), Format::Bed), (b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec(), Format::Vcf),
                     (b"##gff-version 3\nchr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=a\n".to_vec(), Format::Gff3),
                     (b"chr1\tsrc\texon\t1\t9\t.\t+\t.\tgene_id \"g\";\n".to_vec(), Format::Gtf),
                     (b"track type=bedGraph\nchr1\t0\t5\t1.5\n".to_vec(), Format::BedGraph), (Vec::new(), Format::Bed)] {
    for (wrap, comp) in [(plain as fn(&[u8]) -> Vec<u8>, Compression::None), (fixtures::gzip, Compression::Gzip), (fixtures::bgzip, Compression::Bgzf)] {
        let p = fixtures::write(dir.path(), "x.txt", &wrap(&bytes));
        let i = open_input(p.to_str().unwrap(), &opts(4)).unwrap();
        assert_eq!((i.format, i.compression), (fmt, comp));
        assert_eq!(read_all(i.reader), bytes);                       // detection peek is not lost
    }
}
assert!(open_input(write("y", b"hello world\n"), &opts(1)).unwrap_err().to_string().contains("pass --format"));
assert_eq!(open_input(write("z", b"hello\n"), &InputOptions { format: Some(Format::Bed), ..opts(1) }).unwrap().format, Format::Bed);
// truncated BGZF: drop the EOF block
let mut b = fixtures::bgzip(b"chr1\t1\t2\n"); b.truncate(b.len() - 28);
assert!(matches!(open_input(write("t.gz", &b), &opts(1)), Err(IoError::Truncated { .. })));
// parallel decompression gives identical bytes on a 5 MB input
assert_eq!(read_all(open_input(big, &opts(1)).unwrap().reader), read_all(open_input(big, &opts(8)).unwrap().reader));
```

- [ ] **Step 2: Run** `cargo test -p gloci-io --test input`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `detect.rs` and `input.rs`. Peek up to 64 KiB raw (enough for the BGZF header), then up to 64 KiB decompressed, and re-chain each peek with `Cursor::new(peek).chain(rest)` so no bytes are lost.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): content-based format detection, gzip/parallel BGZF input"`

---

### Task 7: Atomic output, output spec, BGZF writing, cleanup handlers

**Files:**
- Create: `crates/gloci-io/src/{atomic.rs,output.rs}`
- Test: `crates/gloci-io/tests/output.rs`

**Interfaces:**
- Consumes: `Format`, `IoError` (Task 6).
- Produces:
  - `pub struct UsageError(pub String)`.
  - `pub enum OutFormat { Bed, BedGraph, Vcf, Bcf, Gff3, Gtf, Tsv }` with `OutFormat::family(Format) -> OutFormat`.
  - `pub struct OutputSpec { pub path: Option<PathBuf>, pub format: Option<OutFormat>, pub bgzf: bool }` and `pub fn parse_output_spec(output: Option<&Path>, output_format: Option<&str>) -> Result<OutputSpec, UsageError>`. Accepted extensions and `-O` values: `bed`, `bed.gz`, `bedgraph`, `bedgraph.gz`, `vcf`, `vcf.gz`, `bcf` (always BGZF), `tsv`, `tsv.gz`, `gff`, `gff.gz`, `gff3`, `gff3.gz`, `gtf`, `gtf.gz` (D17).
  - `pub struct OutputOptions { pub threads: NonZeroUsize, pub level: u8, pub write_index: bool }` and `pub fn check_level(n: u8) -> Result<u8, UsageError>` (0–12, D21).
  - `pub fn refuse_binary_to_tty(is_tty: bool, f: OutFormat, bgzf: bool) -> Result<(), UsageError>`.
  - `pub struct Output` (impl `Write`) and `pub fn open_output(spec: &OutputSpec, f: OutFormat, opts: &OutputOptions) -> Result<Output, IoError>`. `Output::finish(self) -> Result<(), IoError>` finishes BGZF (EOF block), builds the index when requested (calls `index::write_index`, Task 12; until then `write_index` returns `Usage("not yet supported")`), renames, and deletes stale `PATH.tbi`/`PATH.csi` (D19).
  - `AtomicFile::create(dest: &Path) -> io::Result<AtomicFile>` (temp `dest.tmp-XXXXXX` next to dest, registered for cleanup), `.commit(self) -> io::Result<()>`, and `Drop` removes an uncommitted temp file.
  - `pub fn register_temp(path: &Path) -> TempGuard` (the guard unregisters on drop; used for spill dirs).
  - `pub fn install_cleanup_handlers()`: `ctrlc` (termination) and a panic hook remove every registered path; signals exit 130 (D31).
  - `pub fn cleanup_registered()`, used by the handlers and by tests.

BGZF writing always uses `noodles_bgzf::io::MultithreadedWriter` (builder: `set_worker_count`, `set_compression_level`), even with one thread. This makes block boundaries independent of the thread count. Stdout goes through a 1 MiB `BufWriter` over `stdout().lock()`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/output.rs
let s = parse_output_spec(Some(Path::new("o.bed.gz")), None).unwrap(); assert_eq!((s.format, s.bgzf), (Some(OutFormat::Bed), true));
assert_eq!(parse_output_spec(Some(Path::new("O.VCF")), None).unwrap().format, Some(OutFormat::Vcf));
assert!(parse_output_spec(Some(Path::new("o.bcf")), None).unwrap().bgzf);
assert!(parse_output_spec(Some(Path::new("o.txt")), None).unwrap_err().0.contains(".bed.gz"));   // lists accepted extensions
assert_eq!(parse_output_spec(None, Some("vcf.gz")).unwrap().bgzf, true);
assert!(parse_output_spec(Some(Path::new("a.bed")), Some("vcf")).is_err());
assert!(check_level(12).is_ok() && check_level(13).is_err());
assert!(refuse_binary_to_tty(true, OutFormat::Bcf, true).is_err() && refuse_binary_to_tty(false, OutFormat::Bcf, true).is_ok());
assert!(refuse_binary_to_tty(true, OutFormat::Bed, false).is_ok());
// atomic: commit publishes, drop leaves nothing, existing dest untouched on failure
{ let mut f = AtomicFile::create(&dest).unwrap(); f.write_all(b"x").unwrap(); }          // dropped
assert!(!dest.exists() && no_tmp_files(dir.path()));
fs::write(&dest, b"old").unwrap();
{ let mut f = AtomicFile::create(&dest).unwrap(); f.write_all(b"new").unwrap(); }
assert_eq!(fs::read(&dest).unwrap(), b"old");
let mut f = AtomicFile::create(&dest).unwrap(); f.write_all(b"new").unwrap(); f.commit().unwrap();
assert_eq!(fs::read(&dest).unwrap(), b"new"); assert!(no_tmp_files(dir.path()));
// stale index removed when PATH is replaced without --write-index (Review Focus 3)
fs::write(dir.path().join("o.bed.gz.tbi"), b"stale").unwrap();
write_via_output("o.bed.gz", b"chr1\t1\t2\n", /*write_index*/ false);
assert!(!dir.path().join("o.bed.gz.tbi").exists());
// BGZF determinism across worker counts, with EOF block
assert_eq!(bgzf_bytes(&five_mb, 1), bgzf_bytes(&five_mb, 8));
assert!(bgzf_bytes(&five_mb, 8).ends_with(&BGZF_EOF));
// cleanup registry
let p = dir.path().join("spill"); fs::create_dir(&p).unwrap(); let _g = register_temp(&p);
cleanup_registered(); assert!(!p.exists());
```

- [ ] **Step 2: Run** `cargo test -p gloci-io --test output`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `atomic.rs` and `output.rs` per the Interfaces block. Generate the 6-character temp suffix from the PID and a counter mixed with the time, so no extra RNG dependency is needed. Keep the registry in a `Mutex<Vec<PathBuf>>` behind a `static`.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): atomic outputs, output spec, deterministic BGZF writer, cleanup on signals"`

---

### Task 8: Record batches, diagnostics, BED/bedGraph/GFF readers

**Files:**
- Create: `crates/gloci-io/src/{record.rs,diag.rs,bed.rs,gff.rs}`
- Test: `crates/gloci-io/tests/text_readers.rs`

**Interfaces:**
- Consumes: `open_input_from` (Task 6); `ChromTable`, `Interval`, `parse_coord`, `one_based_closed`, `VcfSpan` (core).
- Produces:
  - `pub struct RecordMeta { pub chrom: ChromId, pub key: u64, pub span: Interval, pub loc: u64, pub off: usize, pub len: u32 }` (D1).
  - `pub struct RecordBatch { pub buf: Vec<u8>, pub recs: Vec<RecordMeta> }` with `.raw(i) -> &[u8]`, `.clear()`, `.len()`.
  - `pub enum Header { None, Lines(Vec<u8>), Vcf(VcfHeader) }` (`VcfHeader` comes in Task 9; add the variant then).
  - `pub trait RecordReader: Send` with:
    - `fn name(&self) -> &str`
    - `fn format(&self) -> Format`
    - `fn compression(&self) -> Compression`
    - `fn header(&self) -> &Header`
    - `fn read_batch(&mut self, chroms: &mut ChromTable, batch: &mut RecordBatch, max: usize) -> Result<bool, IoError>`, which appends up to `max` records and returns `Ok(false)` once the input is exhausted
    - `fn diagnostics(&self) -> &Diagnostics`
  - `pub struct ReaderOptions { pub threads: NonZeroUsize, pub format: Option<Format>, pub vcf_span: VcfSpan, pub invalid: InvalidPolicy }`.
  - `pub fn open_reader(path: &str, o: &ReaderOptions) -> Result<Box<dyn RecordReader>, IoError>` and `pub fn open_reader_from(name: &str, src: Box<dyn Read + Send>, o) -> …`. The header is read during open.
  - `pub enum InvalidPolicy { Fail, Skip }` and `pub struct Diagnostics { pub skipped: u64, pub examples: Vec<String> }` (at most 10 examples, each `"line N: msg"`), with `.summary(name) -> Option<String>` giving `skipped {thousands} invalid records in {name}; first {k}:\n  …`.
  - `pub fn bed12_blocks(line: &[u8]) -> Result<Option<Vec<Interval>>, String>` returns None below 12 columns. It errors on a count mismatch, unsorted or overlapping blocks, the first block not at offset 0, or the last block not ending at chromEnd.

Reader rules:
- Raw bytes exclude the line terminator (D15).
- `loc` is the 1-based physical line number.
- Field errors read `{field} "{v}" is not an integer` (field is `start`, `end` or `POS`). Start after end reads `start 10 is greater than end 5`.
- Fewer than 3 columns gives `expected at least 3 tab-separated columns, found {n}`, plus ` (the line contains spaces; BED must be tab-separated)` when the line has spaces.
- An unknown chromosome in a closed table is an invalid record: `chromosome "{name}" is not in {order source}`.
- Under `Skip`, invalid records are counted, the first 10 are kept, and reading continues. Under `Fail`, the first one is returned as `IoError::Parse`.
- GFF uses `key = start-1` and `span = one_based_closed(start, end)`. `##FASTA` gives `IoError::Parse` with the message `##FASTA sections are not supported yet; remove the FASTA section first`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/text_readers.rs (helper `read_all(bytes, opts) -> (Header, Vec<(String, u64, Interval, u64, Vec<u8>)>, Diagnostics)`)
let (h, recs, _) = read_all(b"#c\ntrack name=x\nchr1\t5\t9\tn\r\n\nchr1\t7\t7\n#mid\nchr2\t0\t3\n", bed());
assert_eq!(h, Header::Lines(b"#c\ntrack name=x\n".to_vec()));
assert_eq!(recs, vec![("chr1".into(), 5, iv(5,9), 3, b"chr1\t5\t9\tn".to_vec()),
                      ("chr1".into(), 7, Interval::point(7), 5, b"chr1\t7\t7".to_vec()),
                      ("chr2".into(), 0, iv(0,3), 7, b"chr2\t0\t3".to_vec())]);
assert_eq!(err(b"chr1\t1\t2\nchr1\t10\t5\n"), "x.bed line 2: start 10 is greater than end 5");
assert_eq!(err(b"chr1\t-1\t5\n"), r#"x.bed line 1: start "-1" is negative"#);
assert!(err(b"chr1 1 5\n").ends_with("(the line contains spaces; BED must be tab-separated)"));
let (_, recs, d) = read_all(b"chr1\t1\t2\nbad\nchr1\t3\t1\nchr1\t4\t5\n", skip());
assert_eq!((recs.len(), d.skipped), (2, 2)); assert!(d.examples[0].starts_with("line 2:"));
// batch size does not change results
assert_eq!(read_all_max(input, 1), read_all_max(input, 1000));
// GFF
let (h, recs, _) = read_all(b"##gff-version 3\nchr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g\n###\n", gff());
assert_eq!((recs[0].1, recs[0].2), (99, iv(99, 200))); assert_eq!(recs.len(), 1);
assert!(err_gff(b"##gff-version 3\n##FASTA\n>x\nACGT\n").contains("##FASTA"));
// BED12
assert_eq!(bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t2\t10,20,\t0,80,").unwrap().unwrap(), vec![iv(100,110), iv(180,200)]);
assert!(bed12_blocks(b"c\t100\t200\tn\t0\t+\t100\t200\t0\t2\t10,20,\t0,70,").is_err());    // last block ends at 190
// Review Focus 5: 40 KB lines straddling 64 KiB BGZF blocks read intact
let long = (0..20).map(|i| format!("chr1\t{i}\t{}\t{}\n", i + 1, "x".repeat(40_000))).collect::<String>();
let (_, recs, _) = read_all(&fixtures::bgzip(long.as_bytes()), bed());
assert!(recs.iter().all(|r| r.4.len() > 40_000) && recs.len() == 20);
```

- [ ] **Step 2: Run** `cargo test -p gloci-io --test text_readers`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `record.rs`, `diag.rs`, `bed.rs` and `gff.rs`. Lines are split with `memchr` over a `BufRead::fill_buf` window that grows for long lines. Records are appended into `batch.buf` with no per-record allocation. `open_reader_from` dispatches on the detected format; Vcf and Bcf return `Detect` "not yet supported" until Tasks 9 and 10.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): record batches, skip-invalid diagnostics, BED/bedGraph/GFF readers"`

---

### Task 9: VCF text reader and header

**Files:**
- Create: `crates/gloci-io/src/vcf.rs`
- Modify: `crates/gloci-io/src/record.rs` (add `Header::Vcf`, dispatch)
- Test: `crates/gloci-io/tests/vcf.rs`

**Interfaces:**
- Consumes: `footprint`, `Variant`, `VcfSpan` (Task 4); `RecordReader` (Task 8).
- Produces:
  - `pub struct VcfHeader { pub raw: Vec<u8>, pub contigs: Vec<(Vec<u8>, Option<u64>)>, pub samples: Vec<Vec<u8>>, pub bcf: Option<BcfDicts> }`. `raw` holds every header line, including `#CHROM`, with newlines. `bcf` is filled by Task 10; it is `None` for text input.
  - `VcfHeader::parse(raw: Vec<u8>) -> Result<VcfHeader, String>` parses `##contig=<ID=…,length=…>` and `#CHROM` by hand (D22).
  - `VcfHeader::with_command_line(&self, cmd: &str, version: &str) -> Vec<u8>` returns `raw` with `##glociCommand={cmd}; Version={version}` inserted immediately before `#CHROM`.
  - `pub fn shell_quote_args(args: &[OsString]) -> String` single-quotes args that contain whitespace or shell metacharacters, and escapes tabs and newlines as `\t` and `\n` so the header stays one line.

Record rules:
- Parse the first 8 columns only.
- `key = POS-1`; `span = footprint(.., vcf_span)`.
- INFO is scanned for `END=` and `SVLEN=` only when ALT contains `<`.
- Errors:
  - `POS "{v}" is not an integer`
  - `POS 0 (telomere) is not supported`
  - `expected at least 8 tab-separated columns, found {n}`
  - The footprint error text.
- A missing `#CHROM` line gives `Detect`: `VCF header has no #CHROM line`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/vcf.rs
const H: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr2,length=500>\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";
let (h, recs, _) = read_all(format!("{H}chr1\t100\t.\tACGT\tA\t.\tPASS\t.\tGT\t0/1\t0/0\n").as_bytes(), vcf(VcfSpan::Allele));
let Header::Vcf(h) = h else { panic!() };
assert_eq!(h.contigs, vec![(b"chr2".to_vec(), Some(500)), (b"chr1".to_vec(), None)]);
assert_eq!(h.samples, vec![b"s1".to_vec(), b"s2".to_vec()]);
assert_eq!((recs[0].1, recs[0].2, recs[0].3), (99, iv(100, 103), 5));          // key POS-1, allele span, line 5
assert_eq!(read_all(same_input, vcf(VcfSpan::Ref)).1[0].2, iv(99, 103));
let out = String::from_utf8(h.with_command_line("gloci sort in.vcf", "0.1.0")).unwrap();
assert!(out.contains("##contig=<ID=chr1>\n##glociCommand=gloci sort in.vcf; Version=0.1.0\n#CHROM"));
assert_eq!(err(format!("{H}chr1\t12a4\t.\tA\tG\t.\t.\t.\n")), r#"x.vcf line 5: POS "12a4" is not an integer"#);
assert!(err(format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\t.\n")).contains("needs INFO/END or INFO/SVLEN"));
assert_eq!(read_all(format!("{H}chr1\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=500\n"), vcf(VcfSpan::Allele)).1[0].2, iv(100, 500));
assert!(err(format!("{H}chr1\t0\t.\tA\tG\t.\t.\t.\n")).contains("POS 0"));
assert_eq!(shell_quote_args(&["gloci".into(), "view".into(), "a b.vcf".into()]), "gloci view 'a b.vcf'");
// Review Focus 1: header-only VCF
let (h, recs, _) = read_all(H.as_bytes(), vcf(VcfSpan::Allele)); assert!(recs.is_empty() && matches!(h, Header::Vcf(_)));
```

- [ ] **Step 2: Run** `cargo test -p gloci-io --test vcf`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `vcf.rs` (header reader and `VcfReader: RecordReader`) and the dispatch in `open_reader_from`.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): VCF reader with allele-aware spans and provenance header line"`

---

### Task 10: BCF raw reader, decoder and writer framing

**Files:**
- Create: `crates/gloci-io/src/bcf.rs`
- Modify: `crates/gloci-testkit/src/fixtures.rs` (add `vcf_to_bcf`)
- Test: `crates/gloci-io/tests/bcf.rs`

**Interfaces:**
- Consumes: `VcfHeader` (Task 9), `footprint` (Task 4).
- Produces:
  - `pub struct BcfDicts { pub contigs: Vec<Vec<u8>>, pub strings: Vec<Vec<u8>> }`. These are the BCF contig and string dictionaries, honouring `IDX=`, with `PASS` at string 0.
  - `pub fn read_bcf_header(r: &mut dyn Read) -> Result<VcfHeader, String>`. It accepts magic `BCF\2\1` and `BCF\2\2` and fills `VcfHeader.bcf`.
  - `pub fn write_bcf_header(w: &mut dyn Write, header_text: &[u8]) -> io::Result<()>` writes `BCF\2\2`, `l_text` and the NUL-terminated text.
  - `BcfReader: RecordReader`. Each raw record is the whole framed record (`l_shared`, `l_indiv`, shared, indiv), so writing is a byte copy. `loc` is the 1-based record number. `key = POS` (0-based on disk). `span` comes from the decoded alleles, INFO/END and INFO/SVLEN via `footprint`.
  - `pub(crate) mod typed` decodes BCF typed values: int8/16/32 with missing and end-of-vector sentinels, float, char strings, and the overflow length encoding.
  - testkit: `fixtures::vcf_to_bcf(vcf_text: &[u8]) -> Vec<u8>`, written with noodles-bcf so the fixtures come from an independent implementation.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/bcf.rs
// one shared VCF fixture with SNV, deletion, insertion, multi-allelic, <DEL> with END, <DEL> with SVLEN only,
// <*> gVCF block, breakend, ALT ".", contigs declared with IDX= out of order, two samples
for span in [VcfSpan::Allele, VcfSpan::Ref] {
    let from_vcf = read_all(FIXTURE_VCF, vcf(span)).1;                      // (chrom, key, span) per record
    let from_bcf = read_all(&fixtures::vcf_to_bcf(FIXTURE_VCF), bcf(span)).1;
    assert_eq!(strip_raw(from_vcf), strip_raw(from_bcf), "{span:?}");     // identical keys and footprints
}
// dictionaries match noodles' string maps for the same header
assert_eq!(read_header(bcf_bytes).bcf.unwrap().contigs, noodles_contig_names(FIXTURE_VCF));
// raw passthrough: header + records copied byte-for-byte re-read by noodles with identical records
let copied = copy_via_raw(&bcf_bytes);
assert_eq!(noodles_records(&copied), noodles_records(&bcf_bytes));
// truncated record → error naming the record number
let mut t = bgzf_decompress(&bcf_bytes); t.truncate(t.len() - 5);
assert!(err_bcf(&fixtures::bgzip(&t)).contains("record 9"));
// typed-value decoder unit cases live in bcf.rs #[cfg(test)]: int8 [1, missing], int16 END, int32 SVLEN vector with end-of-vector, overflow length
```

- [ ] **Step 2: Run** `cargo test -p gloci-io --test bcf`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `bcf.rs`. Shared-section layout: CHROM i32, POS i32, rlen i32, QUAL f32, n_allele_info u32, n_fmt_sample u32, then ID, the alleles, FILTER and INFO pairs. Look up INFO keys `END`/`SVLEN` through `strings`. Map contig index to name to `ChromTable` id through a cached `Vec<Option<ChromId>>`.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): BCF framing reader and typed decoder verified against noodles"`

---

### Task 11: Genome loading, chromosome-table resolution, OrderGuard

**Files:**
- Create: `crates/gloci-io/src/{genome.rs,sorted.rs}`
- Test: `crates/gloci-io/tests/genome_order.rs`

**Interfaces:**
- Consumes: `Genome`, `ChromTable`, `OrderChecker` (core); `open_input`, `read_bcf_header`, `VcfHeader` (Tasks 6, 9 and 10).
- Produces:
  - `pub fn load_genome(path: &Path) -> Result<Genome, IoError>`. It detects by content: `@HD`/`@SQ` means `.dict`, `##fileformat=VCF` means a VCF header (plain or bgz), `BCF\2` means BCF, `BAM\1` (after BGZF) means the BAM binary reference list (D27), and anything else is parsed as `.fai`/`.genome`. Errors carry the file and line.
  - `pub fn chrom_table_for(genome: Option<&Path>, inputs: &[&dyn RecordReader]) -> Result<ChromTable, IoError>`. Precedence: `-g` gives `closed(.., GenomeFile)`, then the first VCF/BCF header with ≥1 contig gives `closed(.., VcfHeader(name))`, then `lexicographic()` (spec §4.5, D7).
  - `pub struct OrderError { pub path: String, pub prev_loc: Loc, pub prev: String, pub loc: Loc, pub rec: String, pub hint: String }`. Its Display is:
    ```
    {path} is not sorted: {loc} ({rec}) comes after {prev_loc} ({prev})
    hint: {hint}
    ```
    Records are shown as the raw line cut to 120 chars, or as `CHROM POS` for BCF.
  - `pub struct OrderGuard` with `::new(name: &str, format: Format)` and `.check(&mut self, chroms: &ChromTable, batch: &RecordBatch) -> Result<(), IoError>`.
  - `pub fn sort_hint(name: &str, source: &OrderSource, chrom_backwards: bool) -> String`:
    - `GenomeFile(p)` → `run: gloci sort {name} -g {p}`
    - otherwise → `run: gloci sort {name}`
    - `Lexicographic` with `chrom_backwards` additionally appends ` (or pass -g <genome> if the file uses another chromosome order)`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/genome_order.rs
assert_eq!(names(load_genome(&write("g.fai", b"chr1\t248956422\t52\t60\t61\nchr2\t242193529\t253404903\t60\t61\n")).unwrap()), ["chr1", "chr2"]);
assert_eq!(names(load_genome(&write("g.genome", b"chrA\t10\n")).unwrap()), ["chrA"]);
assert_eq!(load_genome(&write("g.dict", DICT)).unwrap().length(0), Some(248956422));
assert_eq!(names(load_genome(&write("h.vcf.gz", &fixtures::bgzip(VCF_WITH_CONTIGS))).unwrap()), ["chr2", "chr1"]);
assert_eq!(names(load_genome(&write("h.bcf", &fixtures::vcf_to_bcf(VCF_WITH_CONTIGS))).unwrap()), ["chr2", "chr1"]);
assert_eq!(names(load_genome(&write("h.bam", &fixtures::bam_header(&[("chrM", 16569)]))).unwrap()), ["chrM"]);
assert!(load_genome(&write("bad.genome", b"chr1\tlong\n")).unwrap_err().to_string().contains("bad.genome line 1"));
// precedence
assert_eq!(*chrom_table_for(Some(&g), &[&vcf_reader]).unwrap().order_source(), OrderSource::GenomeFile(g.clone()));
assert!(matches!(chrom_table_for(None, &[&vcf_reader]).unwrap().order_source(), OrderSource::VcfHeader(_)));
assert_eq!(*chrom_table_for(None, &[&bed_reader]).unwrap().order_source(), OrderSource::Lexicographic);
assert_eq!(*chrom_table_for(None, &[&vcf_without_contigs]).unwrap().order_source(), OrderSource::Lexicographic);
// guard message names file, both lines (thousands), both records, hint
let e = guard_err(b"chr2\t5\t9\nchr10\t1\t2\n").to_string();
assert_eq!(e, "x.bed is not sorted: line 2 (chr10\t1\t2) comes after line 1 (chr2\t5\t9)\n\
               hint: run: gloci sort x.bed (or pass -g <genome> if the file uses another chromosome order)");
```
(Add `fixtures::bam_header(refs: &[(&str, u32)]) -> Vec<u8>` to testkit: a BGZF-wrapped `BAM\1`, `l_text=0`, `n_ref` and the refs.)

- [ ] **Step 2: Run** `cargo test -p gloci-io --test genome_order`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `genome.rs` and `sorted.rs`. `OrderGuard` keeps a copy of the previous record's raw bytes in a reused `Vec<u8>`, so the error can quote both records.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): --genome loading, chromosome order resolution, sortedness guard"`

---

### Task 12: Index build and indexed region reader

**Files:**
- Create: `crates/gloci-io/src/index.rs`
- Modify: `crates/gloci-io/src/output.rs` (call `write_index` from `Output::finish`)
- Create: `crates/gloci-testkit/src/{strategies.rs,naive.rs}` (first consumers here)
- Test: `crates/gloci-io/tests/index.rs`

**Interfaces:**
- Consumes: readers (Tasks 8–10), `OrderGuard` (Task 11), `overlaps` (Task 3), noodles-tabix `index::Indexer` / noodles-csi `binning_index::Indexer` and header presets `Builder::{bed,gff,vcf}()` plus `set_line_skip_count`.
- Produces:
  - `pub enum IndexKind { Tbi, Csi }` and `pub struct IndexInfo { pub kind: IndexKind, pub path: PathBuf, pub stale: bool }`.
  - `pub fn find_index(data: &Path) -> Option<IndexInfo>` looks for `PATH.tbi` and `PATH.csi`, prefers the newer, and sets `stale` when the index is older than the data (D26).
  - `pub fn write_index(data: &Path, final_data_path: &Path, format: Format, threads: NonZeroUsize) -> Result<PathBuf, IoError>` builds from the finished BGZF file. It checks order in first-appearance mode, chooses TBI or CSI (D20), writes atomically, and returns the index path. Non-BGZF input gives `IoError::Usage` with `{name} is not BGZF-compressed; recompress with: gloci view {name} -o {name}.gz`.
  - `pub fn open_region_reader(path: &Path, regions: &[(Vec<u8>, Interval)], o: &ReaderOptions, inclusive_points: bool) -> Result<Box<dyn RecordReader>, IoError>` yields only the records whose `span` overlaps some region under §4.3. Each record is yielded once, in file order. Regions are ordered by the index's reference order and merged (D10).
  - testkit:
    - `strategies::bed_records() -> impl Strategy<Value = Vec<TestRec>>`. It covers points, bookended features, duplicates, nesting, huge intervals (up to 2^40), mixed strands, and genomes that include chromosomes with no records.
    - `TestRec { chrom: String, start: u64, end: u64, strand: char, name: String }` with `.bed_line()`.
    - `strategies::regions(max: usize) -> impl Strategy<Value = Vec<(Vec<u8>, Interval)>>` (overlapping, duplicate and out-of-order regions included).
    - `naive::matches(rec: &TestRec, chrom: &[u8], iv: Interval, inclusive: bool) -> bool` (the §4.3 rule on the record's span).
    - `naive::region_query(recs, chrom, iv, inclusive) -> Vec<TestRec>`.
    - `naive::sort(recs, chrom_rank: impl Fn(&str) -> u32) -> Vec<TestRec>` (stable, by `(rank, start, end)`).

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-io/tests/index.rs
proptest! {
  #[test] fn region_reader_matches_naive(recs in strategies::bed_records(), regions in strategies::regions(4)) {
    let sorted = naive::sort(recs, lexicographic_rank);
    let p = write_bgzf_bed(&sorted); write_index(&p, &p, Format::Bed, nz(2)).unwrap();
    for incl in [false, true] {
      let got = lines(open_region_reader(&p, &regions, &bed_opts(), incl).unwrap());
      // filtering the sorted list keeps each record once, in file order
      let want: Vec<_> = sorted.iter().filter(|r| regions.iter().any(|(c, iv)| naive::matches(r, c, *iv, incl))).map(TestRec::bed_line).collect();
      prop_assert_eq!(got, want);
    }
  }
}
#[test] fn kinds() {
    assert_eq!(write_and_find(BED_SMALL, Format::Bed).kind, IndexKind::Tbi);
    assert_eq!(write_and_find(b"chr1\t600000000\t600000010\n", Format::Bed).kind, IndexKind::Csi);  // > 2^29
    assert_eq!(write_and_find_bcf(VCF_FIXTURE).kind, IndexKind::Csi);
    noodles_tabix::fs::read(write_and_find(BED_SMALL, Format::Bed).path).unwrap();               // interop
}
#[test] fn zero_length_bed_record() {             // the user's insertion-site case through an index
    let p = indexed(b"chr1\t10\t10\tins\n");
    assert_eq!(count(&p, "chr1", iv(5, 15), false), 1);
    assert_eq!(count(&p, "chr1", iv(10, 15), false), 0);
    assert_eq!(count(&p, "chr1", iv(10, 15), true), 1);
}
#[test] fn vcf_insertion_at_region_edge() {         // D9
    let p = indexed_vcf("chr1\t200\t.\tA\tACGT\t.\t.\t.\n");   // point [200,200)
    assert_eq!(count_vcf(&p, iv(99, 200), VcfSpan::Allele, false), 0);
    assert_eq!(count_vcf(&p, iv(99, 200), VcfSpan::Allele, true), 1);
    assert_eq!(count_vcf(&p, iv(99, 200), VcfSpan::Ref, false), 1);
}
#[test] fn unsorted_and_non_bgzf_are_refused() {
    assert!(matches!(write_index(&bgzf_bed(b"chr1\t5\t6\nchr1\t1\t2\n"), ..), Err(IoError::Order(_))));
    assert!(write_index(&plain_bed(), ..).unwrap_err().to_string().contains("not BGZF-compressed"));
}
#[test] fn empty_inputs_index_fine() {               // Review Focus 1
    let p = bgzf_bed(b""); write_index(&p, &p, Format::Bed, nz(1)).unwrap();
    assert_eq!(count(&p, "chr1", iv(0, 100), false), 0);
    let v = bgzf_vcf(HEADER_ONLY); write_index(&v, &v, Format::Vcf, nz(1)).unwrap();
}
#[test] fn stale_index_flagged() { /* touch data after index → find_index(..).unwrap().stale */ }
```

- [ ] **Step 2: Run** `cargo test -p gloci-io --test index`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `index.rs`. Decompress with `MultithreadedReader` and take chunk virtual positions from `virtual_position()`. If TBI construction fails with `EndExceedsMaxPosition`, rebuild as CSI. The BED preset uses `set_line_skip_count(n)` for leading `track`/`browser` lines. Index spans follow D20. Wire `Output::finish` to `write_index` when `write_index` is set.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(io): tbi/csi index build and index-backed region reader"`

---

### Task 13: CLI framework and `gloci index`

**Files:**
- Create: `crates/gloci/src/{lib.rs,cli.rs,logging.rs,commands/mod.rs,commands/index.rs}`, `xtask/src/docs.rs`, `crates/gloci/tests/{common/mod.rs,index.rs}`
- Modify: `crates/gloci/src/main.rs`, `crates/gloci/tests/cli_basics.rs`, `xtask/src/main.rs`

**Interfaces:**
- Consumes: `write_index`, `install_cleanup_handlers`, `IoError`, `UsageError` (io).
- Produces:
  - `pub fn command() -> clap::Command`, used by xtask.
  - `pub fn run(args: Vec<OsString>) -> std::process::ExitCode`.
  - `pub struct Ctx { pub threads: NonZeroUsize, pub tmp_dir: PathBuf, pub command_line: String, pub version: &'static str }`.
  - Shared flatten groups in `cli.rs` (names and flags from spec §6.1):
    - `OutputArgs { -o/--output, -O/--output-format, --write-index, --compression-level (default 6) }`
    - `ThreadArgs { -t/--threads, env GLOCI_THREADS }`
    - `GenomeArg { -g/--genome }`
    - `VcfSpanArg { --vcf-span allele|ref (default allele) }`
    - `InvalidArg { --skip-invalid }`
    - `TmpDirArg { --tmp-dir, env GLOCI_TMPDIR }`
    - `RegionArgs { --region (repeatable), --regions-file }`
    - `HeaderArg { --header }`
    - `FormatArg { --format }`
  - Global `-q/--quiet` and `-v/--verbose` (count).
  - `pub fn parse_size(s: &str) -> Result<u64, String>` (D23). It takes an integer with a required `K`/`M`/`G`/`T` suffix. Each command enforces its own minimum; `sort` requires `--max-mem >= 16M`.
  - The xtask subcommands `cargo xtask man` and `cargo xtask completions` write `target/assets/man/*.1` and `target/assets/completions/{gloci.bash,_gloci,gloci.fish}`.

Exit mapping in `run`:
- clap usage errors → 2.
- An error chain containing `UsageError` → print `error: {msg}`, exit 2.
- An error chain containing `io::ErrorKind::BrokenPipe` → exit 0 with nothing printed.
- Anything else → print `error: {err:#}`, exit 1.
- Before dispatch, call `install_cleanup_handlers()`, then build the global rayon pool with `threads`.

Every subcommand carries an `after_long_help` block that starts with `Examples:` and has 2–4 lines of the form `  gloci …`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci/tests/cli_basics.rs (add)
#[test] fn every_command_has_2_to_4_examples() {
    for sub in gloci::command().get_subcommands_mut() {
        let help = sub.render_long_help().to_string();
        let n = help.split("Examples:").nth(1).expect(sub.get_name()).lines().filter(|l| l.trim_start().starts_with("gloci ")).count();
        assert!((2..=4).contains(&n), "{}: {n} examples", sub.get_name());
    }
}
#[test] fn bad_threads_is_usage_error() { gloci().args(["index", "x.gz", "-t", "0"]).assert().code(2);
                                         gloci().env("GLOCI_THREADS", "0").args(["index", "x.gz"]).assert().code(2); }
#[test] fn parse_sizes() { assert_eq!(gloci::cli::parse_size("2G"), Ok(2 << 30)); assert_eq!(gloci::cli::parse_size("512M"), Ok(512 << 20));
                           assert!(gloci::cli::parse_size("10").is_err()); assert!(gloci::cli::parse_size("1.5G").is_err()); }
// crates/gloci/tests/index.rs
#[test] fn index_writes_tbi() { let p = bgzf_bed(dir, SORTED_BED); gloci().args(["index", p]).assert().success().stdout("");
                                assert!(exists(format!("{p}.tbi"))); }
#[test] fn index_plain_file_explains() { gloci().args(["index", plain]).assert().code(1)
                                         .stderr(contains("error: ").and(contains("not BGZF-compressed; recompress with: gloci view"))); }
#[test] fn index_unsorted_names_lines() { gloci().args(["index", unsorted_bgzf]).assert().code(1).stderr(contains("is not sorted: line 2")); }
#[test] fn man_and_completions_generate() { /* run `cargo xtask man` and `completions` via std::process; assert files exist */ }
```

- [ ] **Step 2: Run** `cargo test -p gloci`. Expected: FAIL.
- [ ] **Step 3: Implement** `cli.rs`, `logging.rs` (a `log` backend printing `warning: …` / `info: …` to stderr, with level from `-q`/`-v`), `commands/index.rs`, `lib.rs` and `main.rs`, plus `xtask/src/docs.rs` using `clap_mangen` and `clap_complete`.
- [ ] **Step 4: Run** `cargo test -p gloci && cargo xtask man && cargo xtask completions`. Expected: PASS; the files exist.
- [ ] **Step 5: Commit** `git commit -am "feat(cli): shared options, exit codes, help examples, gloci index, man/completions"`

---

### Task 14: Sort engine (in-memory + external merge)

**Files:**
- Create: `crates/gloci-ops/src/sort/{mod.rs,chunk.rs,spill.rs,merge.rs}`
- Test: `crates/gloci-ops/tests/sort_engine.rs`

**Interfaces:**
- Consumes: `RecordBatch`, `RecordReader`, `OrderGuard`, `register_temp` (io); `ChromTable::ranks` (core); testkit `strategies`, `naive`.
- Produces:
  - `pub enum Tie { InputOrder, FullLine }` and `pub enum Dedup { All, Unique, Duplicates }`.
  - `pub struct SortConfig { pub max_mem: u64, pub tmp_dir: PathBuf, pub threads: NonZeroUsize, pub tie: Tie, pub dedup: Dedup }`.
  - `pub struct Sorter` with:
    - `::new(cfg) -> Result<Sorter, IoError>`
    - `.push(&mut self, chroms: &ChromTable, batch: &RecordBatch) -> Result<(), IoError>`, which may spill
    - `.finish(self, chroms: &ChromTable, out: &mut dyn FnMut(&[u8]) -> io::Result<()>) -> Result<SortStats, IoError>`
  - `pub struct SortStats { pub records: u64, pub spills: u32 }`.
  - `pub fn merge_sorted(readers: &mut [Box<dyn RecordReader>], chroms: &mut ChromTable, tie: Tie, out: &mut dyn FnMut(&[u8]) -> io::Result<()>) -> Result<u64, IoError>` (`-m`; every input goes through an `OrderGuard`).

Order is `(rank[chrom], key, span.end, tiebreak)`, where the tiebreak is the global input sequence number, or the raw bytes then the sequence number under `FullLine`. Because the sequence number is unique, sorting is deterministic under any thread count (rayon `par_sort_unstable_by`).

Spilling:
- A chunk's arena and entry vector are pre-sized from `max_mem` and never grow past it. The exception is a single record larger than the budget, which gets its own chunk.
- Each run is written as an LZ4 frame of `(chrom u32, key u64, end u64, seq u64, len u32, bytes)` in a `tempfile::TempDir` under `tmp_dir`, registered for cleanup.
- Merge uses a `BinaryHeap` with fan-in 128 and cascading passes when there are more runs (D24).
- Ranks are recomputed from `ChromTable` at each chunk sort and at the merge; spill files store `ChromId`.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci-ops/tests/sort_engine.rs
proptest! {
  #[test] fn matches_naive_in_memory_and_spilled(recs in strategies::bed_records()) {
    let want = naive::sort(recs.clone(), lexicographic_rank).iter().map(TestRec::bed_line).collect::<Vec<_>>();
    prop_assert_eq!(sort_lines(&recs, cfg(64 << 20, Tie::InputOrder, Dedup::All)), want.clone());
    prop_assert_eq!(sort_lines(&recs, cfg(4 << 10, Tie::InputOrder, Dedup::All)), want);   // forces many spills
  }
  #[test] fn unique_and_duplicates_like_uniq(recs in strategies::bed_records()) {
    let full = naive::sort_full_line(recs.clone(), lexicographic_rank);
    prop_assert_eq!(sort_lines(&recs, cfg(4 << 10, Tie::FullLine, Dedup::Unique)), naive::uniq(&full));
    prop_assert_eq!(sort_lines(&recs, cfg(4 << 10, Tie::FullLine, Dedup::Duplicates)), naive::uniq_d(&full));
  }
}
#[test] fn many_runs_cascade() { let s = sort_stats(&gen_records(200_000), cfg(64 << 10, ..)); assert!(s.spills > 128); /* + equals naive */ }
#[test] fn thread_count_does_not_change_output() { assert_eq!(sort_bytes(&input, threads(1)), sort_bytes(&input, threads(8))); }
#[test] fn record_larger_than_budget() {      // Review Focus 5
    let big = format!("chr1\t5\t6\t{}", "x".repeat(3 << 20));
    assert_eq!(sort_lines_raw(&[&big, "chr1\t1\t2"], cfg(1 << 20, ..)), vec!["chr1\t1\t2".to_string(), big]);
}
#[test] fn spill_dir_removed_after_finish_and_on_drop() { /* tmp_dir empty after finish(); and after dropping a Sorter mid-push */ }
#[test] fn merge_sorted_rejects_unsorted_input() { assert!(matches!(merge(&[sorted_a, unsorted_b]), Err(IoError::Order(_)))); }
```

- [ ] **Step 2: Run** `cargo test -p gloci-ops --test sort_engine`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `sort/` per the Interfaces block, using testkit's `naive::{sort_full_line, uniq, uniq_d}` (add them to testkit).
- [ ] **Step 4: Run** `cargo test -p gloci-ops`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat(ops): external merge sort with bounded memory and deterministic order"`

---

### Task 15: `gloci sort` command

**Files:**
- Create: `crates/gloci/src/commands/sort.rs`, `crates/gloci-io/src/writer.rs`
- Test: `crates/gloci/tests/sort.rs`

**Interfaces:**
- Consumes: `Sorter`, `merge_sorted` (Task 14); `open_reader`, `chrom_table_for`, `parse_output_spec`, `open_output`, `refuse_binary_to_tty` (io); `Ctx`, shared args (Task 13).
- Produces:
  - `pub trait RecordWriter { fn write_record(&mut self, raw: &[u8]) -> io::Result<()>; fn finish(self: Box<Self>) -> Result<(), IoError>; }`.
  - `pub struct WriterOptions { pub write_header: bool, pub command_line: Option<String>, pub version: &'static str }`.
  - `pub fn open_record_writer(out: Output, records: Format, to: OutFormat, header: &Header, o: &WriterOptions) -> Result<Box<dyn RecordWriter>, IoError>`. Text formats get `raw + '\n'`; BCF records are copied as framed bytes. The header is written first: `Lines` only when `write_header`, VCF/BCF always with `with_command_line`. A family mismatch is a `UsageError` (D18).
- CLI: `gloci sort FILES... [-g] [--max-mem SIZE=2G] [--tmp-dir] [-t] [--full-key] [--unique|--duplicates] [-m/--merge-sorted] [--check] [--header] [--skip-invalid] [OutputArgs]`.
  - Readers open with `VcfSpan::Ref` (D2).
  - Multiple VCF inputs need identical sample lists (D29). Mixed formats are a usage error.
  - `--check` writes nothing and exits 0 or 1 (D12).
  - Skip-invalid diagnostics are printed as warnings at the end.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci/tests/sort.rs
#[test] fn bed_lexicographic_default() {
    gloci().args(["sort", &w("u.bed", "chr2\t5\t10\nchr10\t1\t2\nchr1\t7\t9\nchr1\t3\t4\n")]).assert().success()
      .stdout("chr1\t3\t4\nchr1\t7\t9\nchr10\t1\t2\nchr2\t5\t10\n");
}
#[test] fn genome_order() { /* -g g.genome (chr1,chr2,chr10) → chr1,chr1,chr2,chr10 */ }
#[test] fn ties_keep_input_order_and_full_key() { /* "chr1\t1\t5\tb" before "chr1\t1\t5\ta" kept; with --full-key a before b */ }
#[test] fn vcf_contig_order_and_provenance() {
    let out = stdout_of(["sort", &vcf_with_contigs_chr2_chr1]);
    assert!(out.contains("##glociCommand=") && out.contains("; Version=0.1.0\n#CHROM"));
    assert!(out.find("\nchr2\t").unwrap() < out.find("\nchr1\t").unwrap());
}
#[test] fn vcf_undeclared_contig_is_error() { gloci().args(["sort", v]).assert().code(1).stderr(contains("chromosome \"chrUn\" is not in")); }
#[test] fn bcf_roundtrip() { /* -o s.bcf; noodles reads sorted records equal to naive-sorted fixture */ }
#[test] fn gff_and_header() { /* GFF keyed by start-1, leading ## kept; BED --header keeps "#h" line, default drops it */ }
#[test] fn bgzf_with_index() { /* -o s.bed.gz --write-index → BGZF magic, s.bed.gz.tbi loads with noodles */ }
#[test] fn check_mode() {
    gloci().args(["sort", "--check", sorted]).assert().success().stdout("");
    gloci().args(["sort", "--check", unsorted]).assert().code(1).stderr(contains("is not sorted: line 2").and(contains("hint: run: gloci sort")));
}
#[test] fn merge_sorted_inputs() { /* -m a.bed b.bed → merged; -m with unsorted input → exit 1 */ }
#[test] fn unique_duplicates() { /* dup lines → --unique one copy each; --duplicates only the repeated ones */ }
#[test] fn tiny_max_mem_matches_default() { assert_eq!(out(["sort", big, "--max-mem", "16M"]), out(["sort", big])); } // ~50 MB generated input
#[test] fn usage_errors() {
    gloci().args(["sort", bed, "-o", "x.vcf.gz"]).assert().code(2);
    gloci().args(["sort", bed, "--max-mem", "1M"]).assert().code(2);
    gloci().args(["sort", vcf_samples_ab, vcf_samples_ac]).assert().code(1).stderr(contains("sample"));
}
#[test] fn in_place_rewrite() { gloci().args(["sort", p, "-o", p]).assert().success(); assert_eq!(read(p), SORTED); }  // Review Focus 3
#[test] fn empty_inputs() {                  // Review Focus 1
    gloci().args(["sort", empty_bed]).assert().success().stdout("");
    assert!(stdout_of(["sort", header_only_vcf]).ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    gloci().args(["sort", empty_bed, "-o", "e.bed.gz", "--write-index"]).assert().success(); assert!(exists("e.bed.gz.tbi"));
    gloci().args(["sort", "-"]).write_stdin("").assert().success().stdout("");
}
#[test] fn many_contigs() { /* 100k scaffolds × 3 records, shuffled → sorted lexicographically, completes < 10 s */ }   // Review Focus 4
```

- [ ] **Step 2: Run** `cargo test -p gloci --test sort`. Expected: FAIL.
- [ ] **Step 3: Implement** `writer.rs` and `commands/sort.rs`. Order of operations: open the readers, resolve the `ChromTable`, then stream batches into the `Sorter` (or call `merge_sorted` for `-m`, or run the `OrderGuard` only for `--check`), then `finish` into the `RecordWriter`.
- [ ] **Step 4: Run** `cargo test -p gloci`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat: gloci sort (external, -g/VCF contig order, --check, -m, --unique/--duplicates)"`

---

### Task 16: `gloci view` streaming and conversion

**Files:**
- Create: `crates/gloci-ops/src/view.rs`, `crates/gloci-io/src/convert.rs`, `crates/gloci/src/commands/view.rs`
- Test: `crates/gloci/tests/view.rs`

**Interfaces:**
- Consumes: `open_reader`, `open_record_writer` (Task 15), `open_output`, noodles-vcf/bcf readers and writers.
- Produces:
  - `pub fn convert_variants(input: &str, out: Output, to: OutFormat, o: &WriterOptions, threads) -> Result<u64, IoError>` handles VCF↔BCF via noodles `variant::io::Write`, keeping the original header lines plus `##glociCommand`.
  - `pub struct ViewOptions { pub header_only: bool, pub no_header: bool, pub regions: Vec<Region>, pub inclusive_points: bool, pub genome: Option<PathBuf> }`.
  - `pub fn run_view(path: &str, opts: &ViewOptions, ro: &ReaderOptions, out_spec: &OutputSpec, oo: &OutputOptions, wo: &WriterOptions) -> Result<(), IoError>`. Regions arrive in Task 17; with none given it streams.
- CLI: `gloci view FILE [REGION...] [--regions-file BED] [--header-only] [--no-header] [--vcf-span] [--inclusive-points] [-g] [--format] [--skip-invalid] [OutputArgs]`.
  - A cross-family `-O` (e.g. VCF to BED) is a usage error: `use gloci convert for format conversion (coming in a later release)`.
  - `--no-header` with BCF output is a usage error.
  - BCF or BGZF output to a TTY is a usage error (D17).

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci/tests/view.rs
#[test] fn vcf_passthrough_unchanged() {
    let out = stdout_of(["view", vcf_gz]);
    assert_eq!(records_of(&out), records_of(FIXTURE_VCF));                        // byte-identical record lines
    assert!(out.contains("##glociCommand=gloci view"));
}
#[test] fn bcf_to_vcf_and_back() {
    gloci().args(["view", bcf, "-O", "vcf"]).assert().success();                   // noodles-parsed records equal fixture
    gloci().args(["view", vcf, "-o", "o.bcf"]).assert().success();                 // noodles reads o.bcf, records equal
}
#[test] fn header_flags() {
    assert!(stdout_of(["view", vcf, "--header-only"]).lines().all(|l| l.starts_with('#')));
    assert!(!stdout_of(["view", vcf, "--no-header"]).contains("##"));
    gloci().args(["view", vcf, "--no-header", "-O", "bcf"]).assert().code(2);
}
#[test] fn bed_recompress_with_index() { gloci().args(["view", bed, "-o", "o.bed.gz", "--write-index"]).assert().success(); }
#[test] fn cross_family_refused() { gloci().args(["view", vcf, "-O", "bed"]).assert().code(2).stderr(contains("gloci convert")); }
#[test] fn bcf_bytes_to_pipe() { assert!(stdout_bytes(["view", bcf]).starts_with(&[0x1f, 0x8b])); }   // not a TTY under assert_cmd
#[test] fn empty_inputs() { gloci().args(["view", empty_bed]).assert().success().stdout("");
                            gloci().args(["view", header_only_bcf, "-O", "vcf"]).assert().success(); }  // Review Focus 1
#[test] fn stdin_input() { gloci().args(["view", "-"]).write_stdin(FIXTURE_VCF).assert().success(); }
```

- [ ] **Step 2: Run** `cargo test -p gloci --test view`. Expected: FAIL.
- [ ] **Step 3: Implement** `convert.rs`, `view.rs` (the streaming path) and `commands/view.rs`.
- [ ] **Step 4: Run** `cargo test -p gloci`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat: gloci view streaming, VCF<->BCF conversion, header options"`

---

### Task 17: `gloci view` region queries

**Files:**
- Create: `crates/gloci-ops/src/bedsearch.rs`
- Modify: `crates/gloci-ops/src/view.rs`, `crates/gloci/src/commands/view.rs`
- Test: `crates/gloci/tests/view_regions.rs`, `crates/gloci-ops/tests/bedsearch.rs`

**Interfaces:**
- Consumes: `find_index`, `open_region_reader` (Task 12); `Region::parse` (Task 2); `OrderGuard`; testkit `naive::region_query`.
- Produces: `pub fn find_chrom_block(file: &mut (impl Read + Seek), len: u64, chrom: &[u8], table: &ChromTable) -> Result<ChromBlock, IoError>`, where `pub enum ChromBlock { At(u64), NotFoundByBisection }`.

Region path in `run_view`:
1. Parse regions with `is_known` taken from the index's reference names, the VCF contigs, or `-g`. Unknown chromosomes print `warning: chromosome "1" not found in x.vcf.gz`, with ` (did you mean "chr1"?)` added when the `chr` prefix differs either way.
2. If `find_index` finds an index, use `open_region_reader`. A stale index triggers a warning.
3. Else, for uncompressed seekable BED (D11), run `find_chrom_block` per region chromosome. It bisects byte offsets for the first line whose chromosome is ≥ the target in `table` order, and accepts only if that line's chromosome equals the target. From there it scans the chromosome under an `OrderGuard`, stopping at the first record with `key >= region.end`. `NotFoundByBisection` falls back to a full scan with `warning: could not locate chr10 by binary search; scanning the whole file (sort it with gloci sort, or index it)`.
4. Otherwise (gzip, BGZF without an index, VCF without an index), print `warning: no index for {name}; scanning the whole file (create one with: gloci index {name})` and stream-filter under an `OrderGuard` in first-appearance mode.

In every path, records are selected with `overlaps(span, region, inclusive_points)`, each once, in file order (D9, D10).

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci/tests/view_regions.rs
#[test] fn insertion_site_inside_feature() {          // the user's §4.3 example, every access path
    for f in [plain_bed("chr1\t10\t10\tins\n"), bgzf_bed_indexed(..), bgzf_bed_unindexed(..)] {
        assert_eq!(stdout_of(["view", &f, "chr1:6-15"]), "chr1\t10\t10\tins\n");
        assert_eq!(stdout_of(["view", &f, "chr1:11-15"]), "");                       // point on left edge of [10,15)
        assert_eq!(stdout_of(["view", &f, "chr1:11-15", "--inclusive-points"]), "chr1\t10\t10\tins\n");
    }
}
#[test] fn vcf_span_controls_region_hits() { /* insertion POS 200 vs chr1:100-200: allele → none; --vcf-span ref → hit */ }
#[test] fn overlapping_regions_once_in_order() {
    assert_eq!(stdout_of(["view", idx, "chr1:150-300", "chr1:100-200"]), naive_once_sorted(&[(150,300),(100,200)]));
}
#[test] fn long_feature_starting_before_region_found_unindexed() {
    let f = plain_bed("chr1\t0\t1000000\tbig\nchr1\t500\t600\tsmall\nchr1\t999990\t999999\tx\n");
    assert_eq!(stdout_of(["view", &f, "chr1:999991-999995"]), "chr1\t0\t1000000\tbig\nchr1\t999990\t999999\tx\n");
}
#[test] fn karyotype_order_without_genome_still_correct() {
    // chr1, chr2, ..., chr10 order; lexicographic bisection may miss → fallback full scan, correct answer + warning
    gloci().args(["view", &kary, "chr10:1-100"]).assert().success().stdout(CHR10_RECS).stderr(contains("scanning the whole file"));
}
#[test] fn unsorted_plain_bed_errors_not_wrong() { gloci().args(["view", &unsorted, "chr1:1-100"]).assert().code(1).stderr(contains("is not sorted")); }
#[test] fn no_index_warns_and_scans() { gloci().args(["view", &vcf_gz_no_index, "chr1:1-1000"]).assert().success().stderr(contains("no index for")); }
#[test] fn chr_prefix_hint() { gloci().args(["view", &idx_chr, "1:1-100"]).assert().success().stderr(contains("did you mean \"chr1\"")); }
#[test] fn colon_contig_names() {                      // Review Focus 4
    assert_eq!(count(["view", &hla_vcf_gz, "HLA-A*01:01:01:01"]), 3);
    assert_eq!(count(["view", &hla_vcf_gz, "{HLA-A*01:01:01:01}:1-100"]), 1);
}
#[test] fn regions_file() { /* --regions-file r.bed equals positional regions */ }
// crates/gloci-ops/tests/bedsearch.rs: proptest — for sorted random BED + random region, view via bisection == naive::region_query
```

- [ ] **Step 2: Run** `cargo test -p gloci --test view_regions && cargo test -p gloci-ops --test bedsearch`. Expected: FAIL.
- [ ] **Step 3: Implement** `bedsearch.rs` and the region path in `view.rs` and `commands/view.rs`.
- [ ] **Step 4: Run** `cargo test -p gloci -p gloci-ops`. Expected: PASS.
- [ ] **Step 5: Commit** `git commit -am "feat: gloci view region queries (index, bisection on plain BED, verified fallback)"`

---

### Task 18: `gloci check`

**Files:**
- Create: `crates/gloci-ops/src/check.rs`, `crates/gloci/src/commands/check.rs`
- Test: `crates/gloci/tests/check.rs` (insta snapshots in `crates/gloci/tests/snapshots/`)

**Interfaces:**
- Consumes: readers with `InvalidPolicy::Skip`, `OrderChecker`, `chrom_table_for`, `bed12_blocks`, genome lengths.
- Produces:
  - `pub struct CheckReport { pub name: String, pub format: Format, pub compression: Compression, pub records: u64, pub chroms: Vec<(Vec<u8>, u64)>, pub order: String, pub errors: Problems, pub warnings: Problems }`.
  - `pub struct Problems { pub count: u64, pub first: Vec<String> }` (at most 10).
  - `pub fn check_file(path: &str, genome: Option<&Path>, ro: &ReaderOptions) -> Result<CheckReport, IoError>`.
  - `CheckReport::render(&self) -> String`.

Errors cover parse problems (from diagnostics), order violations (the checker continues after each), `end > chromosome length` and a chromosome missing from `-g` (when `-g` is given), BED12 block problems, and VCF data lines whose column count isn't `8`, or `9 + samples` when samples are declared. Warnings (D25) cover exact duplicate lines among records with the same `(chrom, key)` and blank lines.

The summary layout is fixed:

```
{name}: {FORMAT}, {compression}
  records:      {thousands}
  chromosomes:  {n} ({first three names}, … {n-3} more)
  order:        sorted ({order source}) | NOT sorted
  errors:       {count}
    {line N: msg} ×≤10
  warnings:     {count}
    …
result: OK | FAILED ({k} of {n} files have errors)
```

Exit 1 when any file has errors.

- [ ] **Step 1: Write the failing tests**

```rust
// crates/gloci/tests/check.rs
#[test] fn clean_file() { insta::assert_snapshot!(stdout_of(["check", clean_bed])); gloci().args(["check", clean_bed]).assert().success(); }
#[test] fn problems_listed_and_counted() {
    // fixture with exactly 14 problems: 1 order violation, start>end on line 3, 12 non-integer starts
    let o = gloci().args(["check", messy]).assert().code(1).get_output().stdout.clone();
    let s = String::from_utf8(o).unwrap();
    assert!(s.contains("errors:       14") && s.contains("line 3: start 10 is greater than end 5") && s.contains("order:        NOT sorted"));
    assert_eq!(s.lines().filter(|l| l.starts_with("    line ")).count(), 10);
}
#[test] fn genome_bounds() { gloci().args(["check", bed_past_end, "-g", g]).assert().code(1).stdout(contains("end 300 is past the end of chr1 (length 200)")); }
#[test] fn duplicates_warn_only() { gloci().args(["check", dup_bed]).assert().success().stdout(contains("warnings:     1")); }
#[test] fn vcf_column_count() { gloci().args(["check", vcf_missing_sample_col]).assert().code(1).stdout(contains("expected 11 columns, found 10")); }
#[test] fn bed12_blocks_checked() { gloci().args(["check", bad_bed12]).assert().code(1); }
#[test] fn multiple_files() { gloci().args(["check", clean_bed, messy]).assert().code(1).stdout(contains("result: FAILED (1 of 2 files have errors)")); }
#[test] fn empty_file() { gloci().args(["check", empty_bed]).assert().success().stdout(contains("records:      0")); }   // Review Focus 1
```

- [ ] **Step 2: Run** `cargo test -p gloci --test check`. Expected: FAIL.
- [ ] **Step 3: Implement** `check.rs` and `commands/check.rs`. Detect duplicates by hashing the raw lines of the current `(chrom, key)` group into a `HashSet<u64>`, cleared when the key changes.
- [ ] **Step 4: Run** `cargo test -p gloci --test check && cargo insta review`. Expected: PASS, snapshots accepted after reading them.
- [ ] **Step 5: Commit** `git commit -am "feat: gloci check with one-screen summary"`

---

### Task 19: Synthetic generator, benchmark runner, large-file tests

**Files:**
- Create: `xtask/src/{gen.rs,bench.rs}`, `xtask/bench/cases.toml`, `crates/gloci/tests/large.rs`, `.github/workflows/nightly.yml`
- Test: `xtask/src/bench.rs` (`#[cfg(test)]`), `xtask/tests/gen.rs`

**Interfaces:**
- Produces:
  - `cargo xtask gen bed --records N --chroms K --seed S --profile uniform|deep-nesting|whole-chrom|points [--len-min A --len-max B] -o FILE`. Lengths are log-uniform; `ChaCha8Rng::seed_from_u64(S)`; output is written streaming, unsorted unless `--sorted`.
  - `cargo xtask gen vcf --records N --samples M --seed S -o FILE` mixes SNV, MNV, ins, del, multi-allelic and `<DEL>` with END.
  - `cargo xtask bench [--cases xtask/bench/cases.toml] [--out bench-results] [--check-targets]` runs `hyperfine --warmup 1 --runs 5 --export-json` per case and stores `bench-results/<gloci version>/<case>.json`. `--check-targets` evaluates each case's `targets` and exits 1 on a miss.
  - `cases.toml` schema: `[[case]] name, setup (gen args), commands = { gloci = "...", sort-bed = "...", gnu-sort = "LC_ALL=C sort -k1,1 -k2,2n --parallel=8 -S 2G {in}", bedtools = "..." }, targets = ["gloci < sort-bed", "gloci <= 1.5 * gnu-sort"]`.
  - P0 cases: `sort-plain-10M` (BED, 10M records) and `sort-bgzf-10M` (`.bed.gz` input and output, `-t 8`).

Large-file tests (`#[ignore]`, nightly):
- `sort` on a 2 GiB generated BED with `--max-mem 512M`: the output passes `gloci sort --check`, and child peak RSS is ≤ 512 MiB + 150 MiB.
- `view` of a 2 GiB `.bed.gz`: peak RSS ≤ 200 MB.
- Peak RSS is measured with `libc::getrusage(RUSAGE_CHILDREN)` → `ru_maxrss` after the child exits; it is KiB on Linux and bytes on macOS.

`nightly.yml` runs `cargo test --release -- --ignored` and then `cargo xtask bench --check-targets` inside the differential container image (Task 20). It uploads `bench-results/` as an artifact.

- [ ] **Step 1: Write the failing tests**

```rust
// xtask/tests/gen.rs
#[test] fn same_seed_same_bytes() { assert_eq!(gen(["bed", "--records", "1000", "--chroms", "5", "--seed", "7"]), gen([.. same ..])); }
#[test] fn profiles_are_valid_bed() { for p in ["uniform", "deep-nesting", "whole-chrom", "points"] { /* gloci sort | gloci check → exit 0 */ } }
// xtask/src/bench.rs #[cfg(test)]
#[test] fn targets_evaluate() {
    let m = medians(&[("gloci", 1.0), ("sort-bed", 1.5), ("gnu-sort", 0.8)]);
    assert!(eval("gloci < sort-bed", &m).unwrap());
    assert!(eval("gloci <= 1.5 * gnu-sort", &m).unwrap());
    assert!(!eval("gloci < gnu-sort", &m).unwrap());
    assert!(eval("gloci < nosuch", &m).is_err());
}
```

- [ ] **Step 2: Run** `cargo test -p xtask`. Expected: FAIL.
- [ ] **Step 3: Implement** `gen.rs`, `bench.rs`, `cases.toml`, `large.rs` and `nightly.yml`. Target expressions are `NAME (<|<=) [FACTOR *] NAME`, nothing more.
- [ ] **Step 4: Run** `cargo test -p xtask && cargo test -p gloci --release --test large -- --ignored`. Expected: PASS (the large tests take minutes). Then run `cargo xtask bench --check-targets` locally if `sort-bed` is installed. If the `gloci < sort-bed` target misses, profile with `perf` before tuning, and record the result in the commit message.
- [ ] **Step 5: Commit** `git commit -am "feat(xtask): synthetic data generator, benchmark runner with targets, nightly large-file tests"`

---

### Task 20: Differential harness and `docs/differences.md`

**Files:**
- Create: `xtask/src/diff.rs`, `xtask/diff/cases.toml`, `xtask/docker/Dockerfile`, `docs/differences.md`, `crates/gloci/tests/golden.rs`, `crates/gloci/tests/golden/<case>/{input.*,args,expected}` (committed reference outputs)
- Modify: `.github/workflows/nightly.yml`
- Test: `xtask/src/diff.rs` (`#[cfg(test)]`), `crates/gloci/tests/golden.rs`

**Interfaces:**
- Produces:
  - `xtask/docker/Dockerfile`: `FROM mambaorg/micromamba:1.5.10`, then `micromamba install -y -c conda-forge -c bioconda bedtools=2.31.1 bedops=2.4.41 vcftools=0.1.16 htslib=1.21 bcftools=1.21 hyperfine coreutils`. The image tag is `gloci-diff:p0`.
  - `cargo xtask diff [--image gloci-diff:p0] [--cases xtask/diff/cases.toml] [--only NAME]`. It builds `gloci` for `x86_64-unknown-linux-musl`, mounts it with the generated inputs into the container, and runs each case's gloci and reference commands.
  - `cases.toml` schema: `[[case]] name, input = { gen = [...] } | { file = "..." }, gloci = "...", reference = "...", compare = "exact" | "lines-unordered" | "vcf-records", expected_difference = "<anchor in docs/differences.md>"?`.

Semantics:
- A case without `expected_difference` must match.
- A case with `expected_difference` must differ, and the anchor must exist in `docs/differences.md`. A stale entry (outputs now equal) fails, so the differences list stays honest.

Golden files (bedtools agreement in every CI run, no Docker needed):
- `cargo xtask diff --update-golden` runs each case's reference command inside the container on a small committed input and writes `crates/gloci/tests/golden/<case>/expected`, plus `args` (the gloci arguments) and the input file. Golden inputs stay small (≤ 2,000 records) so they can be committed.
- `crates/gloci/tests/golden.rs` runs on every `cargo test`. For each case directory it runs the gloci binary with `args` on the input and compares the output with `expected` using the case's compare mode. The rules match the live harness: a case without `expected_difference` must match, and a case with one must differ (a stale entry fails).
- A golden case cannot drift silently. When a reference tool version changes, `--update-golden` rewrites `expected`, and the diff shows up in review.

P0 cases:
- `sort-bed-vs-sort-bed`
- `sort-bed-vs-bedtools-sort` (bedtools 2.31.1 `sort -i`)
- `sort-bed-vs-gnu-sort` (`LC_ALL=C sort -k1,1 -k2,2n -s`)
- `sort-vcf-vs-bcftools-sort` (vcf-records)
- `view-region-ref-inclusive-vs-tabix` (`--vcf-span ref --inclusive-points`; must match)
- `view-region-default-vs-tabix` (expected difference `#view-region-footprints`)
- `index-interop-tabix-queries` (tabix queries using the gloci-built index equal queries using the tabix-built index)

`docs/differences.md` starts with these sections, each with an anchor, the gloci behaviour, the old behaviour and the reason:
- default `--vcf-span allele`
- `view` region footprints and points (D9)
- `##glociCommand` header line
- `sort --check` checks `(chrom, start)` only (D12)
- a missing BGZF EOF marker is an error (D30)
- sort tie order versus `sort-bed`, if the harness shows a difference

- [ ] **Step 1: Write the failing tests**

```rust
// xtask/src/diff.rs #[cfg(test)]
#[test] fn compare_modes() {
    assert!(compare("exact", b"a\nb\n", b"a\nb\n"));
    assert!(compare("lines-unordered", b"a\nb\n", b"b\na\n") && !compare("exact", b"a\nb\n", b"b\na\n"));
    assert!(compare("vcf-records", b"##x=1\n#CHROM\nr1\n", b"##y=2\n#CHROM\nr1\n"));
}
#[test] fn expected_difference_rules() {
    let doc = "## View region footprints {#view-region-footprints}\n";
    assert!(validate(&case(Some("view-region-footprints")), doc).is_ok());
    assert!(validate(&case(Some("nope")), doc).is_err());
    assert_eq!(verdict(&case(Some("view-region-footprints")), /*equal*/ true), Verdict::StaleDifference);
    assert_eq!(verdict(&case(None), false), Verdict::Mismatch);
}
#[test] fn manifest_parses_and_all_anchors_exist() { /* load xtask/diff/cases.toml + docs/differences.md */ }
// crates/gloci/tests/golden.rs
#[test] fn every_golden_case_agrees_with_its_reference_tool() {
    let cases = golden_cases();                                 // one per crates/gloci/tests/golden/<case>/
    assert!(cases.iter().any(|c| c.name == "sort-bed-vs-bedtools-sort"));
    for c in cases { assert_eq!(c.verdict(run_gloci(&c.args, &c.input)), Verdict::Pass, "{}", c.name); }
}
```

- [ ] **Step 2: Run** `cargo test -p xtask`. Expected: FAIL.
- [ ] **Step 3: Implement** `diff.rs`, the Dockerfile, `cases.toml` and `differences.md`, and add a nightly job: `docker build -t gloci-diff:p0 xtask/docker && cargo xtask diff`.
- [ ] **Step 4: Run** `cargo test -p xtask`, then `docker build -t gloci-diff:p0 xtask/docker && cargo xtask diff --update-golden && cargo xtask diff`, then `cargo test -p gloci --test golden`. Expected: every case passes; any mismatch is fixed or documented. If the container can't run on this host, generate the golden files in CI and note it in the report.
- [ ] **Step 5: Commit** `git commit -am "feat(xtask): differential harness against pinned bedtools/bedops/vcftools/htslib"`

---

### Task 21: Fuzz targets

**Files:**
- Create: `fuzz/Cargo.toml`, `fuzz/fuzz_targets/{bed_reader,vcf_reader,bcf_reader,gff_reader,bgzf_input,region_parse}.rs`
- Modify: `.github/workflows/nightly.yml`

**Interfaces:**
- Consumes: `open_reader_from` (with `ReaderOptions { invalid: Skip, .. }` and with `Fail`), `open_input_from`, `Region::parse`.

Each target reads every batch and must never panic. `region_parse` additionally asserts that a successful parse gives `start <= end` whenever `end` is `Some`. `bcf_reader` and `bgzf_input` wrap the fuzz bytes with `fixtures::bgzip` half the time, so the fuzzer gets past the BGZF layer.

- [ ] **Step 1: Write the targets.** Each is ≤ 20 lines using `libfuzzer_sys::fuzz_target!`.
- [ ] **Step 2: Run** `cargo +nightly fuzz run <target> -- -max_total_time=60` for each target. Expected: no crashes. Fix any crash in the owning module, adding a regression unit test with the crashing input.
- [ ] **Step 3: Add a nightly job** that runs each target for 60 s on `ubuntu-24.04` with `cargo-fuzz` installed, and uploads crash artifacts.
- [ ] **Step 4: Commit** `git commit -am "test: fuzz targets for all readers, BGZF input and region parsing"`

---

### Task 22: Cross-command behaviour suite

**Files:**
- Test: `crates/gloci/tests/behaviour.rs`

**Interfaces:**
- Consumes: the binary only.

- [ ] **Step 1: Write the tests**

```rust
// crates/gloci/tests/behaviour.rs
#[test] fn thread_count_never_changes_bytes() {
    let cases: &[&[&str]] = &[&["sort", BED_100K], &["sort", VCF_GZ, "-O", "vcf.gz"], &["sort", BED_100K, "--max-mem", "16M", "-O", "bed.gz"],
                              &["view", BCF, "-O", "vcf"], &["view", VCF, "-O", "bcf"], &["view", IDX_VCF_GZ, "chr1:1-50000"], &["check", BED_100K]];
    for args in cases { let r: Vec<_> = ["1", "2", "8"].iter().map(|t| stdout_bytes_with(args, &["-t", t])).collect();
                        assert!(r.windows(2).all(|w| w[0] == w[1]), "{args:?}"); }
    // index bytes too
    assert_eq!(index_bytes_with("-t", "1"), index_bytes_with("-t", "8"));
}
#[test] fn broken_pipe_is_quiet_success() {
    let mut child = std::process::Command::new(bin()).args(["view", BIG_BED]).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn().unwrap();
    let mut first = [0u8; 16]; child.stdout.as_mut().unwrap().read_exact(&mut first).unwrap(); drop(child.stdout.take());
    let out = child.wait_with_output().unwrap();
    assert_eq!(out.status.code(), Some(0)); assert!(out.stderr.is_empty());
}
#[test] fn truncated_bgzf_errors() { gloci().args(["view", &truncated_mid_block]).assert().code(1);
                                    gloci().args(["view", &missing_eof_block]).assert().code(1).stderr(contains("truncated")); }
#[test] fn no_partial_outputs_on_error() {
    gloci().args(["sort", &bed_with_bad_line_at_90pct, "-o", &out]).assert().code(1);
    assert!(!Path::new(&out).exists() && no_tmp_files(dir));
}
#[test] fn killed_mid_write_leaves_nothing() {      // Review Focus 2
    let tmp = dir.join("spill"); fs::create_dir(&tmp).unwrap();
    let mut c = Command::new(bin()).args(["sort", BIG_BED, "--max-mem", "16M", "--tmp-dir", tmp, "-o", &out]).spawn().unwrap();
    wait_until(|| fs::read_dir(&tmp).unwrap().next().is_some());          // spilling has started
    unsafe { libc::kill(c.id() as i32, libc::SIGTERM) }; let st = c.wait().unwrap();
    assert_eq!(st.code(), Some(130));
    assert!(!Path::new(&out).exists() && no_tmp_files(dir) && fs::read_dir(&tmp).unwrap().next().is_none());
}
#[test] fn missing_index_region_query_warns() { /* covered in Task 17; assert exit 0 + warning here for .vcf.gz too */ }
#[test] fn error_format() { gloci().args(["sort", &bad_pos_vcf]).assert().code(1)
                            .stderr(r#"error: bad.vcf line 3: POS "12a4" is not an integer"# .to_string() + "\n"); }
```

- [ ] **Step 2: Run** `cargo test -p gloci --test behaviour`. Expected: PASS. Any failure is a defect in the owning task's module. Fix it there with a unit test, then re-run.
- [ ] **Step 3: Run the full suite** `cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings && cargo test --workspace`. Expected: PASS.
- [ ] **Step 4: Commit** `git commit -am "test: cross-command determinism, broken pipe, truncation, atomicity and signal behaviour"`
