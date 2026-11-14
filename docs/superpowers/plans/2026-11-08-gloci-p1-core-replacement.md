# gloci P1 Core Replacement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build spec Phase 1. That means:
- one sorted sweep engine, plus an unsorted mode;
- the relation commands built on that engine (`intersect`, `subtract`, `closest`, `map`, `coverage`);
- the region and transform commands (`merge`, `complement`, `slop`, `flank`, `windows`);
- `convert` and `vcf filter`;
- TSV indexing, per-chromosome parallelism for indexed inputs, P1 benchmarks, `SKILL.md`, and release binaries.

**Architecture:**
- **`gloci-engine`** (new crate) holds the core intersection primitive (spec §8.2). For each A record, it yields every B record whose footprint overlaps A, tagged by its source file. On request, it also yields the nearest B records on each side. The sorted path is a sweep with an active set. The unsorted path loads B into per-chromosome interval indexes.
- **`gloci-ops`** commands are reporting layers over that one primitive: a piece model (whole A, overlap piece, inverse piece, or nothing) plus appended columns.
- **`gloci-vcf`** (new crate) holds genotype decoding, site statistics and the `vcf filter` selection pipeline.
- Every operation takes `RecordReader` streams and writes to a `RecordWriter`. The CLI only wires arguments.

**Tech Stack:**
- Everything P0 uses.
- coitrees 0.4 (unsorted index).
- rust-lapper 1 and superintervals (benchmark comparison only, dev-dependencies of `gloci-engine`).
- criterion (latest stable; engine micro-benchmarks).
- GitHub Actions: `cargo-zigbuild` or `cross` for musl cross-compilation, and `docker/build-push-action` for the release image.

Use the newest release of each new crate at implementation time, and let `Cargo.lock` pin it.

**Spec:** `docs/design.md` (amended 2026-11-08: §4.3 point thresholds, §7.1 `--mode` has no short flag). This plan implements the P1 row of §10, using §4, §5, §6, §7.1–7.4, §7.8, §8 and §9.

**Builds on:** P0 (`docs/superpowers/plans/2026-11-06-gloci-p0-foundations.md`, merged as `cb30575`). P0's decisions D1–D33 stay in force, except D5, which D34 supersedes. The P0 follow-up list is `docs/superpowers/plans/2026-11-06-gloci-p0-followups.md`.

**Out of scope (P2/P3):**
- `--chrom-map`, `--region` on relation commands (D52), filter expressions, and the full `map` op set (D46);
- `partition`, `cluster`, `shift`, `spacing`, `groupby`, `genomecov`;
- `vcf stats`, `vcf extract`, `vcf concat`;
- the Python library, WIG input, BAM/CRAM/FASTA, the mdBook and bioconda.

**P1 is done when:**
- every command in the P1 row runs with its differential cases passing against pinned bedtools, bedops and vcftools;
- `docs/migration.md` maps the old commands to them;
- the nightly x86_64 run of `cargo xtask bench --check-targets` meets the §9.1 targets for interval operations and `vcf filter`.

## Global Constraints

- Internally every feature is a **0-based, half-open** interval `[start, end)`; positions are **64-bit unsigned**.
- Rust stable, edition 2024, `rust-version = "1.91"`. Release profile: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`.
- **No htslib.** Pure Rust plus vendored C; the binary must build as one static `musl` binary.
- Exit codes: **0** success, **1** data or runtime error, **2** usage error. Broken pipe exits quietly with status 0.
- Errors say what, where and how to fix, e.g. `error: calls.vcf.gz line 1,204,331: POS "12a4" is not an integer` (line numbers with thousands separators).
- Output is byte-identical for every `--threads` value (spec §6.4).
- `.gz` output is always **BGZF**. `--compression-level N` default **6**. Atomic file output: write `PATH.tmp-XXXX` and rename on success, leaving no partial file on failure.
- Format detection uses content, never extensions. `--format`, `--a-format` and `--b-format` override it. `-` reads stdin.
- VCF output gets a `##glociCommand=<full command line>; Version=<x.y.z>` header line.
- `--vcf-span allele|ref`, default **`allele`**.
- `-t, --threads N` default `available_parallelism()`, env `GLOCI_THREADS`. `--tmp-dir DIR`, env `GLOCI_TMPDIR`.
- **Overlap (spec §4.3).** Non-empty intervals overlap when they share a base. A point overlaps an interval when strictly inside it (`--inclusive-points` allows the edges), and equal points overlap. Every threshold counts as met when either side is a point. **Distance (spec §4.4):** overlapping features have distance 0; otherwise it is the gap plus 1.
- Sortedness is checked on the fly in every streaming command and is never turned off. The first violation exits 1, names the file, the line and both records, and gives a hint (`run: gloci sort FILE`, plus the lex note when the pair looks lexicographic).
- **Record streams (spec §8.1).** Every operation in `gloci-ops` takes `Box<dyn RecordReader>` inputs and writes to `&mut dyn RecordWriter`. It never takes only file paths.
- **Output formats (spec §4.6).** Commands that output input records unchanged write A's format: VCF A gives VCF with its header. Commands that append columns write TSV. A VCF A then gets a single `#`-prefixed column-name line and no meta-header.
- The same flag means the same thing everywhere. Flag names come from spec §6.1 and §7. Every command's `--help` has 2–4 real examples.
- Streaming memory stays flat as input grows (≤ ~200 MB RSS, excluding active-set depth). `--unsorted` memory is proportional to the loaded B side, and the help text says so.
- **Testing (spec §9.3).** Every two-input and n-ary operation has a property test against a naive O(n·m) implementation. Every command has a thread-determinism test (`-t 1`, `-t 2`, `-t 8`). Every command has differential cases against the pinned reference tools in `xtask/diff/cases.toml`; single-invocation cases become golden files checked by `cargo test`. Every difference from a reference tool is listed in `docs/differences.md` under an anchor and pinned by a case with `expected_difference`.
- CI platforms: Linux x86_64 and aarch64, macOS aarch64.

## Decisions this plan adds

P0's D1–D33 remain. These continue the numbering. Change any of them before execution starts.

| # | Decision |
|---|---|
| D34 | Every threshold (`-f -F -r -e --min-bases`) counts as met when either side is a point (spec §4.3, amended). Supersedes D5. `-r` with only one of `-f`/`-F` applies that one fraction to both sides. |
| D35 | `intersect --mode` has no short flag (spec §7.1, amended). |
| D36 | The engine yields **candidate hits**: every B whose window-padded span overlaps A's span under the exact §4.3 rule (with `--inclusive-points` when given). Thresholds and strand rules are applied by the reporting layer through `OverlapCriteria`, measured against the **unpadded** B span. |
| D37 | `-w`/`--left`/`--right` cannot be combined with `-f`/`-F`/`-r`/`-e`/`--min-bases` (usage error, exit 2). bedtools `window` has no fractions either, and a fraction of a padded interval has no clear meaning. |
| D38 | Hits for one A are ordered by B file (command-line order), then B input order. A records are reported in A's input order, in every mode, sorted or unsorted. |
| D39 | Two flags, one meaning each (revised 2026-11-08 with the user). **`--unsorted`** = "my input is in some order; work with it as it is": on `intersect`, `subtract`, `closest`, `map` and `coverage` only, all `-b` files are loaded into memory, A is streamed in any order, and results follow A's input order. **`--sort`** = "sort my inputs first": on `merge`, `complement` and the same five relation commands, all inputs go through the P0 external sorter (bounded by `--max-mem`, default 2G, spilling to `--tmp-dir`) and the normal streaming pass runs on the sorted streams, so results are in sorted order and equal `gloci sort …` piped into the command. `--sort` and `--unsorted` together are a usage error. An order error's hint names the flags the command offers. `windows -i` and transforms need neither. |
| D40 | Unsorted index: coitrees for a chromosome whose padded B coordinates all fit in `i32`. Otherwise an own `u64` implicit augmented interval tree (cgranges-style), `ITree64`. Queries are widened by one base on each side, and every candidate is re-tested exactly (spec §8.2). Task 27's benchmark may swap the default with evidence. |
| D41 | Per-chromosome parallelism is automatic when A and every B are regular BGZF files with a `.tbi`/`.csi` and `--unsorted` is not given. It applies to `intersect`, `subtract`, `closest`, `map` and `coverage`. Each chromosome's output goes to a chunk that keeps its first 64 MiB in memory and spills the rest to an LZ4 file in `--tmp-dir`. Chunks are written in chromosome-table order, and at most `threads` are in flight. The output equals the streaming path byte for byte. |
| D42 | Appended-column output is TSV. `.tsv`/`.tsv.gz` (or `-O tsv`) always work. When A is BED-family, `.bed`/`.bed.gz`/`.bedgraph` also work, since the result is BED-like. Any other extension is a usage error. Stdout defaults to TSV. VCF A gets the line `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO[\tFORMAT\t<samples>]\t<appended names>`. BED/GFF A get their header lines only with `--header`. |
| D43 | `--write-index` on TSV output is allowed when A is BED-family; the index uses columns seq 1, begin 2, end 3, zero-based. For VCF or GFF A it is a usage error whose hint gives the matching `gloci index -s … -b … -e …` command. |
| D44 | Coordinate rewriting (`--mode clip`, `subtract` pieces, `slop`, `flank`): BED columns 2–3; GFF/GTF columns 4–5 written as 1-based closed. BED12 thick and block columns are left unchanged (bedtools behaviour). `--mode clip` with VCF/BCF A is a usage error ("a variant cannot be clipped; use --mode any"). |
| D45 | Placeholders for a missing B (`--mode pairs --keep-unmatched`, `closest` with no candidate), byte-for-byte as bedtools 2.31.1 writes them (revised after Task 9's differential cases): BED-family `. -1 -1` then `.` per further column, with `-1` in the score column (BED5+); GFF `. . . -1 -1 . . . .`; VCF `. -1 . . . . . .`. A missing distance is `-1`. |
| D46 | P1 aggregation ops: `count count_distinct sum mean median min max stdev var stdev_pop var_pop collapse distinct distinct_sorted first last bases indicator` (`stdev_pop`/`var_pop` added during Task 14 so bedtools' population `stdev` can be reproduced). Any other spec §6.3 op (P2) is a usage error ("… is not available yet"). Value rules: `.` or empty values are missing. A value that isn't a number in a numeric op is a data error naming the file, line, column and op. `sum`/`min`/`max` print as integers when every input value was an integer literal. Other floats print in Rust's shortest round-trip form. With no values, `count`/`count_distinct` print 0 and the others print `--empty` (default `.`). `stdev`/`var` are sample statistics, and with one value they print `--empty`. `distinct` keeps first-seen order; `distinct_sorted` sorts bytewise. COL `0` is allowed only with `count`, `bases` and `indicator`. |
| D47 | `merge -s` keeps clusters separate per strand but adds no strand column (bedtools). Output is sorted by chromosome, start, end. |
| D48 | `complement` ignores points, since they cover no bases. Without `-g` it reports only gaps between features on the same chromosome. |
| D49 | `slop`/`flank` accept BED and GFF/GTF; VCF/BCF input is a usage error. `-s` treats unknown strand as `+`. A transform output that breaks sort order prints one warning naming the first out-of-order line, with the hint `pipe the output through gloci sort`. |
| D50 | `vcf filter` rules: <ul><li>A genotype is **missing** when any of its alleles is missing or GT is absent.</li><li>Masking sets every allele missing, keeping the ploidy and separators.</li><li>Missing DP or GQ never triggers masking.</li><li>FT passes when it is `PASS`, `.` or absent.</li><li>A missing QUAL fails `--min-qual`.</li><li>`--pass-only` keeps FILTER `PASS` and `.`.</li><li>A record's **type set** is the union of its ALT kinds. `--types` keeps a record whose set meets the list, and `--exclude-types` drops one whose set meets it. `indel` is `ins` or `del`.</li><li>`--min-alleles`/`--max-alleles` count listed alleles (REF plus ALTs; ALT `.` gives 1).</li><li>`--thin` runs after every other site filter.</li><li>Dropped INFO keys also lose their `##INFO` header lines.</li><li>`--recalc-ac` writes AF like C `%g` (6 significant digits).</li></ul> |
| D51 | Relation chromosome table: `chrom_table_for(-g, [A, B…], --chrom-order)`, so the first VCF/BCF header among A then the Bs wins. Every input is order-checked against that one table. |
| D52 | `--region`/`--regions-file` exist on `view` (P0) and `vcf filter`. On the relation commands they are deferred to P2. |
| D53 | The mdBook and the bioRxiv external benchmark datasets are deferred to P2. |
| D54 | Strand comes from BED column 6 when the line has at least 6 columns and the value is `+`, `-` or `.`, and from GFF/GTF column 7. Anything else, and all VCF/BCF records, is unknown. `-s`/`-S` warn once per input that has unknown-strand records. |
| D55 | `closest` also takes `-s`/`-S` (bedtools users rely on "nearest gene on the same strand"). The engine keeps neighbour candidates per (file, strand) when a strand rule is active. |
| D56 | `coverage` of a point A: count of B overlapping it (spec §4.3), bases 0, length 0, fraction `.`; `--mean` prints `.`; `--hist`/`--per-base` print no lines for it. |
| D57 | With several `-b` files: `--mode count` reports the total. `--mode pairs` and `closest` insert a label column (from `--names`, else the basename) before B's columns. `coverage` reports one covered fraction per file (bedtools `annotate`) and allows only the default report. `map` needs all `-b` files in one format family. |
| D58 | Chromosome naming. After a relation (or a multi-input `merge`/`complement`) in which two inputs with records shared no chromosome names, gloci warns. When adding or removing `chr` would make names match, the warning says so. A record on a chromosome missing from a built-in `-g` genome adds a hint naming the matching genome in the other naming style (`hg38`↔`GRCh38` and so on) when that genome has it. |
| D59 | The release workflow runs only on a pushed `v*` tag. It builds static musl x86_64/aarch64 and macOS aarch64 archives with SHA-256 sums, attaches them to the GitHub release and pushes a `ghcr.io` image. It publishes nothing to crates.io or bioconda while the licence is open (spec §13). |
| D60 | `subtract --split` removes B's blocks (exons) from A's span. `subtract --sum-overlaps` removes A when the summed overlap lengths of the B records passing the strand rule reach `-f` of A (bedtools `-N`). It takes only `-f`; `-F`/`-r`/`-e` with it are usage errors. |

## Review Focus

These are the inputs most likely to hurt a user that no single command's happy-path tests would catch. Each has a pinned test in the task named.

1. **Deep overlap.** One whole-chromosome B against 100,000 A, and B nested 1,000 deep: results equal the naive answer, and the peak active set is the depth, not the file size. Task 5 `deep_overlap_whole_chromosome`.
2. **A and B sorted or named differently.** A sorted in karyotype order and B in lex order must exit 1 naming B with the lex hint, never give a partial answer. `chr1` against `1` must warn. Task 5 `b_lex_sorted_is_an_error`, Task 25 `chr_prefix_mismatch_warns`.
3. **Points everywhere.** Every row of the §4.3 table through `intersect`; edge points in `closest` (distance 1); points in `subtract`, `merge`, `complement`, `coverage` and `map`. Task 9 `worked_examples_table` and each command's `points` test.
4. **Empty inputs.** Empty A, empty B, a header-only VCF, and B with no record on A's chromosomes. Each command exits 0 with header-only or A-only output. Task 9 `empty_inputs`, Task 28 `empty_inputs_matrix`.
5. **VCF A through every relation.** Header plus `##glociCommand` in passthrough modes, the `#CHROM…` line in TSV modes, and a trimmed deletion whose span starts after its sort key. Task 5 `eviction_uses_a_key_not_span`, Task 28 `vcf_a_matrix`.

## File structure

```
crates/gloci-core/src/        overlap.rs (D34, blocks), variant_kind.rs (new: kinds, trimming)
crates/gloci-io/src/          record.rs (+strand), fields.rs (new), tsv.rs (new),
                              index/{mod,build,query,tsv}.rs (split of index.rs + TSV index)
crates/gloci-engine/          new: hit.rs, sweep.rs, active.rs, nearest.rs, unsorted.rs,
                              itree64.rs, kmerge.rs, parallel.rs, chunk.rs; benches/index.rs
crates/gloci-ops/src/         relate/{mod,criteria,out,rewrite}.rs, agg/{mod,spec,ops,num}.rs,
                              intersect.rs, subtract.rs, closest.rs, map.rs, coverage.rs,
                              merge.rs, complement.rs, transform.rs, windows.rs,
                              convert/{mod,vcf,gff,bed12}.rs, vcf_filter.rs
crates/gloci-vcf/             new: gt.rs, bcf_fmt.rs, stats.rs, hwe.rs,
                              filter/{mod,site,lists,samples,rewrite}.rs
crates/gloci/src/             cli.rs (+relation arg groups), commands/{intersect,subtract,closest,
                              map,coverage,merge,complement,slop,flank,windows,convert,vcf}.rs
crates/gloci-testkit/src/     naive.rs (+relations), strategies.rs (+relation inputs), mem.rs (new)
xtask/                        diff cases + inputs, bench cases + datasets, docs.rs (SKILL.md), gen.rs (VCF)
.github/workflows/            release.yml (new), nightly.yml, ci.yml
docs/                         migration.md, differences.md; SKILL.md at the repo root
```

---

### Task 1: Core rules for P1

**Files:**
- Modify: `crates/gloci-core/src/overlap.rs`, `crates/gloci-core/src/lib.rs`
- Create: `crates/gloci-core/src/variant_kind.rs`
- Test: `crates/gloci-core/tests/overlap_p1.rs`, `crates/gloci-core/tests/variant_kind.rs`; update the P0 point-threshold tests in `crates/gloci-core/tests/` that encode the old D5 rule

**Interfaces:**
- Produces:
  - `OverlapCriteria::accepts` (same signature, D34 behaviour);
  - `Fraction`, with manual `PartialEq` by value and a stored scale;
  - `pub fn blocks_overlap(a: &[Interval], b: &[Interval], inclusive_points: bool) -> bool`;
  - `pub fn blocks_overlap_len(a: &[Interval], b: &[Interval]) -> u64` (block lists sorted and disjoint);
  - `pub enum VariantKind { Snv, Mnv, Ins, Del, Sv, Ref, Other }`;
  - `pub fn alt_kinds(v: &Variant) -> Vec<VariantKind>`;
  - `pub fn trimmed_alleles<'a>(pos0: u64, ref_: &'a [u8], alt: &'a [u8]) -> (u64, &'a [u8], &'a [u8])`.

  All are re-exported from the crate root.

- [ ] **Step 1: Write the failing tests**

```rust
// overlap_p1.rs
let half = Fraction::parse("0.5").unwrap();
let c = OverlapCriteria { min_frac_b: Some(half), ..Default::default() };
assert!(c.accepts(Interval::point(10), Strand::Unknown, iv(5, 15), Strand::Unknown)); // point A, -F 0.5
let c = OverlapCriteria { min_bases: Some(3), ..Default::default() };
assert!(c.accepts(Interval::point(10), Strand::Unknown, iv(5, 15), Strand::Unknown));
let c = OverlapCriteria { min_frac_a: Some(Fraction::parse("1").unwrap()), ..Default::default() };
assert!(c.accepts(iv(5, 15), Strand::Unknown, Interval::point(10), Strand::Unknown)); // point B, -f 1
assert!(!c.accepts(Interval::point(5), Strand::Unknown, iv(5, 15), Strand::Unknown)); // edge: no overlap at all
// D34: -r with only -F
let c = OverlapCriteria { min_frac_b: Some(half), reciprocal: true, ..Default::default() };
assert!(!c.accepts(iv(0, 10), Strand::Unknown, iv(5, 100), Strand::Unknown)); // 5/10 ok, 5/95 not
assert_eq!(Fraction::parse("1").unwrap(), Fraction::parse("1.0").unwrap());
assert_eq!(Fraction::parse("50%").unwrap(), Fraction::parse("0.5").unwrap());
assert_eq!(blocks_overlap_len(&[iv(0, 10), iv(20, 30)], &[iv(8, 22)]), 4);
assert!(!blocks_overlap(&[iv(0, 10)], &[iv(10, 20)], false));
```

```rust
// variant_kind.rs: §4.2 table rows (POS 100)
kinds("A", "G") == [Snv]; kinds("AC", "GT") == [Mnv]; kinds("ACGT", "A") == [Del];
kinds("A", "ACGT") == [Ins]; kinds("ATG", "AC") == [Other]; kinds("N", "<DEL>") == [Sv];
kinds("G", "G]17:198982]") == [Sv]; kinds("A", "<*>") == [Ref]; kinds("A", ".") == [Ref];
kinds("ACGT", "A,AC") == [Del, Del];
assert_eq!(trimmed_alleles(99, b"ATG", b"AC"), (100, &b"TG"[..], &b"C"[..]));
assert_eq!(trimmed_alleles(99, b"A", b"ACGT"), (100, &b""[..], &b"CGT"[..]));
assert_eq!(trimmed_alleles(99, b"ACGT", b"A"), (100, &b"CGT"[..], &b""[..]));
```

- [ ] **Step 2: Run them:** `cargo test -p gloci-core --test overlap_p1 --test variant_kind`. They fail to compile (missing items), or the point rows fail.
- [ ] **Step 3: Implement.**
  - In `accepts`, let `point = a.is_point() || b.is_point()`; every fraction and `min_bases` check is `point || …`. Make the reciprocal fallback symmetric.
  - `Fraction` stores `scale: u128` (10^digits) at parse time; `met` uses it.
  - Write the block helpers.
  - Write `variant_kind.rs`. Kind of a trimmed pair: lengths 1/1 → Snv; equal length >1 → Mnv; empty REF → Ins; empty ALT → Del; anything else → Other. Symbolic and breakend ALTs → Sv. `.`, `<*>` and `<NON_REF>` → Ref.
  - Update the P0 tests that assert the old point rule, and say why in a comment that cites spec §4.3 (amended).
- [ ] **Step 4: Run** `cargo test -p gloci-core`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(core): thresholds met for point overlaps; blocks; variant kinds and trimming`

### Task 2: Strand and column access on records

**Files:**
- Modify: `crates/gloci-io/src/record.rs` (`RecordMeta` gains `pub strand: Strand`), `crates/gloci-io/src/bed.rs`, `crates/gloci-io/src/gff.rs`, `crates/gloci-io/src/vcf.rs`, `crates/gloci-io/src/bcf.rs`, and every other place that builds a `RecordMeta` (`rg 'RecordMeta \{'`)
- Create: `crates/gloci-io/src/fields.rs`
- Test: `crates/gloci-io/tests/strand.rs`

**Interfaces:**
- Produces:
  - `RecordMeta { chrom, key, span, strand, loc, off, len }`;
  - `pub fn field(raw: &[u8], n: usize) -> Option<&[u8]>` (0-based tab-separated field);
  - `pub fn column_count(raw: &[u8]) -> usize`;
  - `pub fn split_fields(raw: &[u8]) -> impl Iterator<Item = &[u8]>` (memchr-based).

- [ ] **Step 1: Write the failing tests** (`strand.rs`):
  - BED6 `chr1\t0\t5\tn\t0\t+` → `Plus`, and `-` → `Minus`;
  - BED6 with `.` → `Unknown`; BED4 → `Unknown`; BED6 with column 6 `x` → `Unknown`;
  - bedGraph → `Unknown`;
  - GFF column 7 `+` → `Plus`, and `?` → `Unknown`;
  - VCF and BCF → `Unknown`;
  - `field(b"a\tb\tc", 1) == Some(&b"b"[..])`, `field(b"a\tb\tc", 3) == None`, `column_count(b"a\tb\tc") == 3`.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test strand`. Fails to compile.
- [ ] **Step 3: Implement** the D54 rules in each reader's parse path, with no extra allocation (the column scan already happens there).
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: all pass, and the P0 sort/view goldens are unchanged.
- [ ] **Step 5: Commit** `feat(io): record strand (BED col 6, GFF col 7) and field helpers`

### Task 3: TSV indexing and the index module split

**Files:**
- Move `crates/gloci-io/src/index.rs` to `crates/gloci-io/src/index/{mod.rs, build.rs, query.rs}`. This is a pure move with no behaviour change; it was a P0 follow-up.
- Create: `crates/gloci-io/src/index/tsv.rs`
- Modify: `crates/gloci/src/cli.rs` (`IndexArgs`), `crates/gloci/src/commands/index.rs`
- Test: `crates/gloci-io/tests/index_tsv.rs`, `crates/gloci/tests/index_tsv.rs`; cases in `xtask/diff/cases.toml`

**Interfaces:**
- Produces:
  - `pub struct TsvColumns { pub seq: usize, pub begin: usize, pub end: Option<usize>, pub zero_based: bool, pub skip_lines: usize, pub comment: u8 }` (column numbers are 1-based, as the user gives them);
  - `pub fn write_tsv_index(data: &Path, final_data_path: &Path, cols: &TsvColumns, threads: NonZeroUsize) -> Result<PathBuf, IoError>`.
- CLI: `gloci index FILE [-s/--seq-col N -b/--begin-col N [-e/--end-col N] [--zero-based] [--skip-lines N] [--comment CHAR]]` (spec §7.4).

Semantics, as tabix uses them:
- Without `--zero-based`, the interval is `[begin-1, end)`, with `end` = the end column, or `begin` when there is no end column.
- With `--zero-based`, the interval is `[begin, end)`, or `[begin, begin+1)` when there is no end column.
- The first `skip_lines` lines are skipped, as are lines starting with `comment` (default `#`). Both are recorded in the index header.
- The index header uses the generic format code with the column numbers, plus the UCSC zero-based flag (`0x10000`) when zero-based.
- The output is `.tbi`, or `.csi` when an end exceeds 2^29.
- Chromosomes must be contiguous and begins sorted within each one. A violation exits 1 naming the line.
- Any of `-s`/`-b`/`-e` selects TSV mode. `-s` and `-b` must come together. `--zero-based`, `--skip-lines` and `--comment` outside TSV mode are usage errors.

- [ ] **Step 1: Write the failing tests.**
  - `tsv_index_query_matches_lines`: index a BGZF BED-like TSV with extra columns using `-s 1 -b 2 -e 3 --zero-based`. Query `chr1:100-200` directly with noodles-tabix/csi. Expect exactly the lines overlapping `[99,200)`.
  - `one_based_without_end_column`: `-s 1 -b 2`, where a line with begin 150 covers `[149,150)`.
  - `skip_lines_and_comment`: `--skip-lines 1 --comment '%'`.
  - `csi_for_large_end`: an end of 600,000,000 produces `.csi`.
  - CLI: an unsorted TSV exits 1 naming line 3; `-s 1` alone exits 2; `--zero-based` without `-s` exits 2.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test index_tsv` and `cargo test -p gloci --test index_tsv`. They fail.
- [ ] **Step 3: Implement** the move first and run the full suite to prove it changed nothing. Then implement `tsv.rs` on the same builder machinery as P0 (BGZF virtual offsets, file-order plans) and wire up the CLI.
- [ ] **Step 4: Add differential cases.**
  - `index-tsv-bed-like-vs-tabix`: the gloci command is `gloci index -s 1 -b 2 -e 3 --zero-based input.tsv.gz && tabix input.tsv.gz chr1:100-200`; the reference is `tabix -s1 -b2 -e3 -0 input.tsv.gz && tabix input.tsv.gz chr1:100-200`; compare `exact`.
  - The same for a 1-based `-s1 -b2 -e2` file, for `--skip-lines`, and for the large-end CSI case.
  - Run `cargo xtask diff --only index-tsv` in the pinned image and commit the results.
- [ ] **Step 5: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 6: Commit** `feat(index): tabix-style indexing of any bgzipped TSV; split index module`

### Task 4: TSV results writer

**Files:**
- Create: `crates/gloci-io/src/tsv.rs`
- Modify: `crates/gloci-io/src/output.rs`, `crates/gloci-io/src/writer.rs`, `crates/gloci-io/src/lib.rs`
- Test: `crates/gloci-io/tests/tsv_writer.rs`

**Interfaces:**
- Consumes: `TsvColumns` and `write_tsv_index` (Task 3).
- Produces:
  - `pub enum TsvHeader { None, Lines(Vec<u8>), VcfColumns { samples: Vec<Vec<u8>> } }`;
  - `pub fn tsv_output_ok(a_format: Format, to: OutFormat) -> Result<(), UsageError>` (D42);
  - `pub fn open_tsv_writer(out: Output, a_format: Format, appended: &[String], header: TsvHeader, write_index: bool) -> Result<Box<dyn RecordWriter>, IoError>`. It writes the header and returns a writer. `write_record` takes a complete row (A columns plus appended columns, tab-joined, no newline).
  - `Output::set_tsv_index(cols: TsvColumns)`, used by D43.

- [ ] **Step 1: Write the failing tests.**
  - `vcf_a_column_line`: with samples `s1 s2` and appended `["count"]`, the first line is exactly `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\tcount\n`.
  - `vcf_a_without_samples` has no FORMAT column.
  - `bed_a_without_header` writes rows only.
  - `bed_a_with_header` copies the `Lines` bytes first.
  - `tsv_output_ok` accepts `(Bed, Tsv)`, `(Bed, Bed)` and `(Vcf, Tsv)`, and refuses `(Vcf, Vcf)` and `(Gff3, Bed)` with the D42 message.
  - `tsv_index_for_bed_a` writes a `.tbi` that noodles can query.
  - `tsv_index_refused_for_vcf_a` is a usage error containing `gloci index -s 1 -b 2 -e 2`.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test tsv_writer`. It fails.
- [ ] **Step 3: Implement.** `Output::finish` builds the TSV index from the temp file via Task 3, using the same publish order as P0: index built, stale index removed, data renamed, index renamed.
- [ ] **Step 4: Run** `cargo test -p gloci-io`. Expected: pass.
- [ ] **Step 5: Commit** `feat(io): TSV results writer with column-name line and BED-preset index`

### Task 5: `gloci-engine`, a sorted sweep with exact overlaps

**Files:**
- Create: `crates/gloci-engine/Cargo.toml`, and `src/{lib.rs, hit.rs, sweep.rs, active.rs}`. Add it to the workspace and to `[workspace.dependencies]`.
- Create: `crates/gloci-testkit/src/mem.rs` (in-memory readers).
- Modify: `crates/gloci-testkit/src/naive.rs` (`naive_hits`) and `crates/gloci-testkit/src/strategies.rs` (`relation_inputs()`: A and 1–3 B record sets over `CHROMS`, mixing points, bookended, duplicate, nested and huge features with mixed strands).
- Test: `crates/gloci-engine/tests/sweep.rs`

**Interfaces:**
- Consumes: `RecordReader`, `RecordMeta` (with `strand`, Task 2), `OrderGuard`, `ChromTable`, `Window`, `overlaps`.
- Produces:

```rust
pub struct ARec<'a> { pub meta: RecordMeta, pub raw: &'a [u8] }
pub struct Hit<'a> { pub file: u16, pub span: Interval, pub strand: Strand, pub raw: &'a [u8], pub loc: u64 }
pub struct Neighbours<'a> { /* filled by Task 6; empty when nearest is None */ }
impl<'a> Neighbours<'a> { pub fn left(&self, file: u16) -> &[Hit<'a>]; pub fn right(&self, file: u16) -> &[Hit<'a>]; }
pub struct NearestConfig { pub k: usize, pub max_distance: Option<u64>, pub strand: StrandRule } // used from Task 6
pub struct SweepConfig { pub window: Option<Window>, pub inclusive_points: bool, pub nearest: Option<NearestConfig> }
pub trait Visitor { fn visit(&mut self, a: &ARec<'_>, hits: &[Hit<'_>], near: &Neighbours<'_>) -> Result<(), IoError>; }
pub struct Inputs { pub a: Box<dyn RecordReader>, pub b: Vec<Box<dyn RecordReader>> }
pub struct SweepSummary { pub a_records: u64, pub b_records: Vec<u64>, pub shared_chroms: bool }
pub fn sweep_sorted(inputs: Inputs, chroms: &mut ChromTable, cfg: &SweepConfig, v: &mut dyn Visitor) -> Result<SweepSummary, IoError>
// testkit::mem
pub fn mem_reader(name: &str, text: &str, o: &ReaderOptions) -> Box<dyn RecordReader>
```

Algorithm (the decisions the tests can't determine):
- **Order checking.** Every stream gets its own `OrderGuard::new(name, format)` checked against the shared table, so a B violation names B's file.
- **Loading.** B records are read in key order. A B is loaded while its `key <= a.span.end + pad`, where `pad = max(left, right)`. Since `span.start >= key`, this loads every possible candidate.
- **Eviction.** A B leaves the active set when `b.span.end + pad < a.key`. Use A's **key**, not its span start: a later A's span can start before the current A's span (a trimmed VCF deletion), but never before the current key.
- **Matching.** Hits are the active B where `overlaps(a.span, window.apply(b.span, a.meta.strand), inclusive_points)` holds, ordered per D38.
- **Chromosome change.** When A moves to a new chromosome, the active set is cleared, B records on chromosomes ranked earlier are read past (still order-checked), and a B record ahead stays pending.
- **End of input.** Reading stops when A is exhausted; any B order after that point doesn't affect the result.
- **Storage.** B raw bytes live in an arena that is compacted when more than half of it is dead. No allocation per record.
- `SweepSummary.shared_chroms` is true when at least one chromosome had both A and B records (Task 25 uses it).

- [ ] **Step 1: Write the failing tests.**
  - `sweep_matches_naive`: proptest, 256 cases, using `relation_inputs()`. With and without a window (strand-aware or not), and with `inclusive_points` on and off: for each A in input order, `(file, B line)` of the hits equals `naive_hits`.
  - `deep_overlap_whole_chromosome`: one B `[0, 10^9)` against 100,000 A, then 1,000 nested B against 1,000 A. Results equal the naive answer, and a debug counter of the peak active-set size is at most depth + 1.
  - `eviction_uses_a_key_not_span`: A1 is VCF `chr1 101 . ATTTT A` (key 100, span `[101,105)`), A2 is `chr1 102 . A G` (key 101, span `[101,102)`), and B is BED `chr1 100 103`. Both A1 and A2 hit B.
  - `b_lex_sorted_is_an_error`: A is karyotype-sorted and B lex-sorted. The run fails with an `OrderError` naming B's file, with the lex note.
  - `edges_and_empties`:
    - A on a chromosome B lacks, and B on a chromosome A lacks;
    - empty A, empty B, and a header-only VCF A;
    - an A point on B's edge, without and with `inclusive_points`.
  - `in_memory_streams`: the inputs come from `mem_reader`. This is the spec §8.1 requirement.
- [ ] **Step 2: Run** `cargo test -p gloci-engine`. It fails.
- [ ] **Step 3: Implement** `hit.rs`, `active.rs` and `sweep.rs`.
- [ ] **Step 4: Run** `cargo test -p gloci-engine` and `cargo clippy -p gloci-engine --all-targets -- -D warnings`. Expected: pass.
- [ ] **Step 5: Commit** `feat(engine): sorted sweep yielding exact candidate hits per A`

### Task 6: Nearest neighbours in the sweep

**Files:**
- Create: `crates/gloci-engine/src/nearest.rs`
- Modify: `crates/gloci-engine/src/sweep.rs`, `crates/gloci-testkit/src/naive.rs` (`naive_neighbours`)
- Test: `crates/gloci-engine/tests/nearest.rs`

**Interfaces:**
- Produces: `Neighbours::left(file)` and `Neighbours::right(file)`, filled when `SweepConfig.nearest` is `Some`.
  - **Left** holds non-overlapping B with `b.span.end <= a.span.start`. **Right** holds non-overlapping B with `b.span.start >= a.span.end`. Each is sorted nearest first by core `distance()`, with ties in input order.
  - **Guarantee.** For each file, each side contains every B that could be among the `k` nearest on that side, including all ties at the k-th distance, within `max_distance`.
  - When `strand` is not `Any`, the guarantee holds for the B whose strand meets the rule for this A, and B with unknown strand are excluded.

Retention (the guarantee follows from these rules; the test checks it):
- **Left pool** (per file, and per strand when a strand rule is active): drop a B once `k` pool members have a strictly larger end, or once `b.span.end + max_distance < a.key`.
- **Right look-ahead:** read B ahead, per pool, until one of these holds:
  - more than `k` candidates exist and the last has a start strictly beyond the k-th candidate's start;
  - a start exceeds `a.span.end + max_distance`;
  - the chromosome or the input ends.
- Every B with an identical start is retained, since they are ties.

- [ ] **Step 1: Write the failing tests.**
  - `neighbours_match_naive`: proptest with `k` in 1..=3, every `StrandRule`, and `max_distance` on/off. For each A and file, the top k with ties from `left`/`right` equals `naive_neighbours`.
  - `edge_cases`: no B left of A; only B on another chromosome (both sides empty); 5 B at equal distance with `k=1` (all 5 present); `k` larger than the number of B; an edge point at distance 1.
- [ ] **Step 2: Run** `cargo test -p gloci-engine --test nearest`. It fails.
- [ ] **Step 3: Implement** `nearest.rs`, wired into the sweep's load and evict steps.
- [ ] **Step 4: Run** `cargo test -p gloci-engine`. Expected: pass.
- [ ] **Step 5: Commit** `feat(engine): nearest-neighbour candidates for closest`

### Task 7: Unsorted mode

**Files:**
- Create: `crates/gloci-engine/src/unsorted.rs`, `crates/gloci-engine/src/itree64.rs`
- Modify: `crates/gloci-engine/Cargo.toml` (add `coitrees`)
- Test: `crates/gloci-engine/tests/unsorted.rs`

**Interfaces:**
- Produces: `pub fn sweep_unsorted(inputs: Inputs, chroms: &mut ChromTable, cfg: &SweepConfig, v: &mut dyn Visitor) -> Result<SweepSummary, IoError>`, with the same `Visitor` contract as `sweep_sorted` (hit order D38, the `Neighbours` guarantee).
  - Every B is loaded first. A is streamed in input order with no order check, and B needs no order.
- Produces: `pub struct ITree64` with `fn new(items: Vec<(Interval, u32)>) -> ITree64` and `fn query(&self, lo: u64, hi: u64, out: &mut Vec<u32>)`. The query is closed, `[lo, hi]`, returning candidates.

Per D40:
- Use coitrees when every padded B coordinate on a chromosome is below `i32::MAX - 2`, and `ITree64` otherwise.
- Query `[a.start - pad - 1, a.end + pad + 1]` (saturating), re-test every candidate with the Task 5 rule, and sort the hits per D38.
- While A arrives sorted (keys non-decreasing so far), use coitrees' sorted querent. At the first out-of-order A, switch to plain queries.
- **Neighbours:** each chromosome keeps B sorted by start and by end, and binary search fills `left`/`right` under the Task 6 guarantee.

- [ ] **Step 1: Write the failing tests.**
  - `unsorted_equals_sorted`: proptest. Shuffle A, run `sweep_unsorted`, and compare per-A hits and neighbours, keyed by A's line, with `sweep_sorted` on the sorted inputs.
  - `itree64_fallback`: B at 3,000,000,000 on a long chromosome; `ITree64` is used (debug flag) and the results equal the naive answer.
  - `sorted_querent_switch`: sorted A, then one out-of-order A. The results are unchanged.
- [ ] **Step 2: Run** `cargo test -p gloci-engine --test unsorted`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p gloci-engine`. Expected: pass.
- [ ] **Step 5: Commit** `feat(engine): --unsorted via per-chromosome coitrees with a 64-bit fallback`

### Task 8: Relation reporting scaffold and shared CLI arguments

**Files:**
- Create: `crates/gloci-ops/src/relate/{mod.rs, criteria.rs, out.rs, rewrite.rs}`
- Modify: `crates/gloci-ops/src/lib.rs`, `crates/gloci-ops/Cargo.toml` (depends on `gloci-engine`), `crates/gloci/src/cli.rs` (shared groups)
- Test: `crates/gloci-ops/tests/relate.rs`

**Interfaces:**
- Consumes: `sweep_sorted` and `sweep_unsorted` (Tasks 5–7), `open_tsv_writer` (Task 4), `open_record_writer` (P0), `chrom_table_for` (P0), `blocks_overlap` and `blocks_overlap_len` (Task 1), `bed12_blocks` (P0).
- Produces:

```rust
pub struct RelateOptions {
    pub criteria: OverlapCriteria, pub window: Option<Window>, pub split: bool, pub unsorted: bool,
    pub labels: Vec<String>, pub genome: Option<String>, pub chrom_order: ChromOrder, pub header: bool,
}
pub struct RelateInputs { pub a: Box<dyn RecordReader>, pub b: Vec<Box<dyn RecordReader>> }
pub fn run_relation(inputs: RelateInputs, o: &RelateOptions, nearest: Option<NearestConfig>, v: &mut dyn Visitor) -> Result<SweepSummary, IoError>
pub fn accepted(a: &ARec, hit: &Hit, o: &RelateOptions) -> bool         // criteria (unpadded B), --split blocks, D37
pub fn rewrite_coords(format: Format, raw: &[u8], iv: Interval, out: &mut Vec<u8>)   // D44
pub fn b_placeholder(b_format: Format, b_columns: usize, out: &mut Vec<u8>)          // D45
pub enum ReportOutput { SameAsA, Tsv { appended: Vec<String> } }
pub fn open_report(a: &dyn RecordReader, kind: &ReportOutput, spec: &OutputSpec, oo: &OutputOptions, wo: &WriterOptions, header: bool) -> Result<Box<dyn RecordWriter>, IoError>
```

These clap groups go in `cli.rs` for Tasks 9–15 to flatten:
- `AArg { -a FILE, --a-format FMT }`;
- `BArgs { -b FILE... (required, 1+), --names A,B,..., --b-format FMT }`. `--names` must have one label per `-b`, or it is a usage error;
- `OverlapArgs { -f FRAC, -F FRAC, -r, -e, --min-bases N, -s, -S (conflicting), --split, --inclusive-points }`. FRAC is a decimal or a percentage, via `Fraction::parse`;
- `WindowArgs { -w N, --left N, --right N, --strand-window }`, which conflicts with every threshold flag (D37);
- `UnsortedArg { --unsorted }`.

Each flag's help text cites its §4.3 rule in one line.

Behaviour:
- `run_relation` builds the table per D51 and picks `sweep_unsorted` when `o.unsorted`.
- When `-s` or `-S` is set and an input has records with unknown strand, it warns once per input (D54).
- `accepted` handles `--split`: when either side is BED12, blocks come from `bed12_blocks`, and any other record is one block. Overlap is `blocks_overlap`, and thresholds use block lengths and `blocks_overlap_len`.

- [ ] **Step 1: Write the failing tests.**
  - `rewrite_coords`:
    - BED: `chr1\t5\t15\tx` with `[7,9)` → `chr1\t7\t9\tx`;
    - GFF: `chr1\tsrc\tgene\t6\t15\t.\t+\t.\tID=g` with `[6,9)` → columns 4–5 become `7` and `9`.
  - `b_placeholder`: BED3 → `.\t-1\t-1`; BED5 → `.\t-1\t-1\t.\t.`; VCF with 8 columns → eight `.` fields.
  - `accepted`: `-f 0.5` on an unpadded B; `--split` with two BED12 records whose blocks miss each other but whose spans overlap gives false.
  - `relation_from_memory`: `run_relation` on `mem_reader` inputs.
- [ ] **Step 2: Run** `cargo test -p gloci-ops --test relate`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p gloci-ops`. Expected: pass.
- [ ] **Step 5: Commit** `feat(ops): relation reporting scaffold and shared overlap arguments`

### Task 9: `gloci intersect`

**Files:**
- Create: `crates/gloci-ops/src/intersect.rs`, `crates/gloci/src/commands/intersect.rs`, `crates/gloci-ops/tests/intersect_props.rs`
- Modify: `crates/gloci/src/cli.rs` (`Intersect(IntersectArgs)`), `crates/gloci/src/commands/mod.rs`, `crates/gloci/src/lib.rs` (dispatch), `crates/gloci-testkit/src/naive.rs` (`naive_intersect`)
- Test: `crates/gloci/tests/intersect.rs`; cases and inputs under `xtask/diff/`; `docs/differences.md`

**Interfaces:**
- Consumes: Task 8 (`run_relation`, `accepted`, `rewrite_coords`, `b_placeholder`, `open_report`), Task 1 (`trimmed_alleles`).
- Produces:

```rust
pub enum Mode { Any, None, Pairs { overlap_bases: bool, keep_unmatched: bool }, Clip, Count }
pub enum Match { Overlap, Alleles }
pub struct IntersectOptions { pub mode: Mode, pub matching: Match, pub relate: RelateOptions }
pub fn intersect(inputs: RelateInputs, o: &IntersectOptions, w: &mut dyn RecordWriter) -> Result<SweepSummary, IoError>
```

- CLI: `gloci intersect -a A -b B... [--mode any|none|pairs|clip|count] [--overlap-bases] [--keep-unmatched] [--match overlap|alleles]`, plus the flattened `AArg`, `BArgs`, `OverlapArgs`, `WindowArgs`, `UnsortedArg`, genome, chrom-order, vcf-span, header, skip-invalid, threads and output groups. `--overlap-bases` and `--keep-unmatched` need `--mode pairs`, or they are a usage error.

Output (D42, D44, D45, D57):
- `any`/`none`: A's format, passed through.
- `clip`: A's format with rewritten coordinates, once per accepted hit. VCF/BCF A is a usage error.
- `pairs`: TSV. Columns: A raw, then a label when there is more than one `-b`, then B raw, then the overlap bases if requested. `--keep-unmatched` adds rows with the D45 placeholder and overlap 0.
- `count`: TSV of A raw plus the total of accepted hits.
- `--match alleles` needs every input to be VCF/BCF (usage error otherwise). A hit is accepted when some trimmed ALT of A equals some trimmed ALT of B: same chromosome, trimmed position, trimmed REF and trimmed ALT.

- [ ] **Step 1: Write the failing tests** in `intersect.rs`, using `assert_cmd` and fixtures:
  - `worked_examples_table`: each spec §4.3 row as `-a` and `-b` files run with `--mode any`, without and with `--inclusive-points`. The output matches the "Overlap?" column.
  - `modes`: one fixture with exact expected output for any, none, pairs, pairs with `--overlap-bases`, pairs with `--keep-unmatched`, clip and count.
  - `multiple_b_with_names`: pairs rows carry the labels, and count is the total.
  - `thresholds`: `-f 0.5`, `-F 50%`, `-r -f 0.5`, `-e -f 0.9 -F 0.1`, `--min-bases 5`.
  - `strand`: `-s`, `-S`, and an unknown-strand warning printed once.
  - `windows`: `-w 10`, `--left 5 --right 20`, `--strand-window` on a minus-strand A. `-w 10 -f 0.5` exits 2.
  - `vcf_a`: any mode keeps the header plus `##glociCommand`; pairs writes the `#CHROM…` line; clip exits 2.
  - `match_alleles`:
    - A `chr1 100 . A ACGT` matches B `chr1 100 . A ACGT,G` but not B `chr1 100 . A AC`;
    - A `chr1 100 . ATG AC` matches B `chr1 101 . TG C`;
    - a BED input exits 2.
  - `empty_inputs`: empty A gives no records (VCF: header only); empty B with `--mode none` gives all of A.
  - `order_errors`: an unsorted B exits 1 naming B; `--unsorted` on a shuffled A gives the same records in A's order.
  - `determinism`: every mode under `-t 1`, `-t 2` and `-t 8` gives identical bytes.

  Also `intersect_props.rs`: proptest of any/none/pairs/count against `naive_intersect`, sorted and `--unsorted`.
- [ ] **Step 2: Run** `cargo test -p gloci --test intersect` and `cargo test -p gloci-ops --test intersect_props`. They fail.
- [ ] **Step 3: Implement** the op and the CLI. The help text has 4 examples, says that the default mode is `any`, and says that `--unsorted` loads every `-b` file into memory.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases** (bedtools 2.31.1 `-sorted`, bedops 2.4.41; gloci runs with `--chrom-order lex`). Inputs: `rel-a.bed` and `rel-b.bed`, BED6 with strands, no points, lex-sorted.

  | gloci | reference |
  |---|---|
  | `--mode any` | `bedtools intersect -sorted -u` |
  | `--mode none` | `-v` |
  | `--mode pairs` | `-wa -wb` |
  | `--mode pairs --overlap-bases` | `-wo` |
  | `--mode pairs --keep-unmatched --overlap-bases` | `-wao` |
  | `--mode count` | `-c` |
  | `--mode clip` | the bedtools default |
  | `-f 0.5`, `-F 0.5`, `-r -f 0.5`, `-e -f 0.5 -F 0.5`, `-s`, `-S` | the same flags |
  | `-w 10`, `--left 5 --right 20 --strand-window` | `bedtools window -w 10`, `-l 5 -r 20 -sw` |
  | `-f 50%` | `bedops --element-of 50%` |
  | `--mode none` | `bedops --not-element-of 1` |
  | `--match alleles` (`rel-a.vcf`, `rel-b.vcf`) | `bcftools isec -n=2 -w1 -c none` |

  Pinned differences, each with its own `docs/differences.md` anchor:
  - **`#zero-length-overlap`:** points in A and B (`rel-points-*.bed`) against bedtools.
  - **`#vcf-allele-span`:** VCF A `-u` against bedtools.
  - Any other difference the harness shows.

  Run `cargo xtask diff --only intersect` in the image and commit the results.
- [ ] **Step 6: Commit** `feat: gloci intersect (any/none/pairs/clip/count, --match alleles, --unsorted)`

### Task 10: Aggregation engine (`-c COL:OP`)

**Files:**
- Create: `crates/gloci-ops/src/agg/{mod.rs, spec.rs, ops.rs, num.rs}`
- Modify: `crates/gloci-testkit/src/differential.rs` and `xtask/src/diff.rs` (new compare mode `numeric`)
- Test: `crates/gloci-ops/tests/agg.rs`, `crates/gloci-testkit` unit tests for `numeric`

**Interfaces:**
- Produces:

```rust
pub enum Op { Count, CountDistinct, Sum, Mean, Median, Min, Max, Stdev, Var, Collapse, Distinct, DistinctSorted, First, Last, Bases, Indicator }
pub enum Col { Whole, Index(usize), Name(String) }      // "0" → Whole; 1-based index; header/VCF name
pub struct AggSpec { pub col: Col, pub op: Op, pub label: String }   // label = text as typed, e.g. "5:mean"
pub fn parse_aggs(args: &[String]) -> Result<Vec<AggSpec>, UsageError>   // repeatable, comma-separated
pub struct ResolvedAgg { /* column index or VCF field, op */ }
pub fn resolve(specs: &[AggSpec], format: Format, header: &Header) -> Result<Vec<ResolvedAgg>, UsageError>
pub struct Accumulators { /* one per agg, reused across A records */ }
impl Accumulators {
    pub fn new(aggs: Vec<ResolvedAgg>, empty: String, delim: String) -> Self;
    pub fn add(&mut self, raw: &[u8], overlap_bases: u64, file: &str, line: u64) -> Result<(), IoError>;
    pub fn finish_into(&mut self, out: &mut Vec<u8>);   // appends "\tv1\tv2…" and resets
}
```

Rules are D46. Name resolution:
- `QUAL`, `ID` and `FILTER` are VCF columns 6, 3 and 7.
- `INFO/KEY` is KEY's value; a flag gives `1`; an absent key is missing.
- For BED/TSV, other names are matched against the fields of the last header line, with the leading `#` stripped. An unknown name is a usage error that lists the available names.

The `numeric` compare mode compares tab fields: two numbers are equal within a relative tolerance of 1e-5 (absolute when both are below 1e-9), and any other fields must match exactly.

- [ ] **Step 1: Write the failing tests.**
  - `parse`: `["5:mean,5:max", "4:distinct"]` gives three specs with labels `5:mean`, `5:max` and `4:distinct`.
  - `parse_errors`:
    - `5` (no op) and `5:avg` (unknown op) are usage errors;
    - `5:mode` is a usage error saying `mode` is not available yet;
    - `0:sum` is a usage error.
  - `values`: for values `1 2 3 4`:

    | op | output |
    |---|---|
    | sum | `10` |
    | mean | `2.5` |
    | median | `2.5` |
    | min | `1` |
    | max | `4` |
    | stdev | `1.2909944487358056` |
    | var | `1.6666666666666667` |
    | count | `4` |
    | count_distinct | `4` |
    | first | `1` |
    | last | `4` |

  - `missing_and_empty`:
    - `[1.5, 2]` sums to `3.5`;
    - `[".", 3]` has mean `3`;
    - no values: count `0`, sum `.`, and `--empty NA` gives `NA`;
    - one value: stdev prints `.`.
  - `text_ops`: values `b a b` give `distinct` → `b,a`, `distinct_sorted` → `a,b`, and `collapse` → `b,a,b`; `--delim ';'` changes the separator.
  - `non_numeric`: `abc` in `5:sum` gives an error containing `b.bed line 7: column 5 value "abc" is not a number (5:sum)`.
  - `vcf_names`: `INFO/DP:mean` over records with `DP=10` and `DP=20` gives `15`; `QUAL:max` gives the larger QUAL.
- [ ] **Step 2: Run** `cargo test -p gloci-ops --test agg`. It fails.
- [ ] **Step 3: Implement.** Numeric values parse with a fast path for integer literals (kept exact as `i128` sums) and fall back to `f64`.
- [ ] **Step 4: Run** `cargo test -p gloci-ops` and `cargo test -p gloci-testkit`. Expected: pass.
- [ ] **Step 5: Commit** `feat(ops): -c COL:OP aggregation engine (P1 core ops); numeric diff compare mode`

### Task 11: `gloci merge`

**Files:**
- Create: `crates/gloci-engine/src/kmerge.rs`, `crates/gloci-ops/src/merge.rs`, `crates/gloci/src/commands/merge.rs`
- Modify: `cli.rs`, the dispatch, `crates/gloci-testkit/src/naive.rs` (`naive_merge`)
- Test: `crates/gloci/tests/merge.rs`, `crates/gloci-ops/tests/merge_props.rs`; cases

**Interfaces:**
- Consumes: Task 10 aggregation.
- Produces:
  - `pub fn kmerge(inputs: Vec<Box<dyn RecordReader>>, chroms: &mut ChromTable, f: &mut dyn FnMut(usize, &RecordMeta, &[u8]) -> Result<(), IoError>) -> Result<(), IoError>`. It yields records in (chromosome rank, key, input index, input order), with every input order-checked. Task 16 reuses it.
  - `pub struct MergeOptions { pub max_gap: i64, pub per_strand: bool, pub count: bool, pub aggs: Vec<AggSpec>, pub empty: String, pub delim: String }`
  - `pub fn merge(inputs: Vec<Box<dyn RecordReader>>, chroms: &mut ChromTable, o: &MergeOptions, w: &mut dyn RecordWriter) -> Result<(), IoError>`
- CLI: `gloci merge FILES... [-d/--max-gap N] [-s] [--count] [-c COL:OP]... [--empty S] [--delim S]`, plus genome, chrom-order, vcf-span, format, threads and output groups. The output is BED; `.tsv` is also accepted.

Semantics:
- Clusters are built from footprints (`span`). A record joins the open cluster when `start <= cluster_end + max_gap`, in signed arithmetic. So `-d 0` joins overlapping and bookended features, and `-d -N` needs at least N bases of overlap.
- Spans are taken in **start order** through a min-heap fed in key order. A span leaves the heap once the next incoming key exceeds its start, because a span never starts before its key.
- `-s` keeps separate open clusters per strand (`+`, `-`, unknown) and emits them in (start, end, strand) order (D47).
- Output columns: chrom, start, end, then `[count]`, then the aggregation columns (a `bases` aggregation is a usage error here).
- A cluster of points at `p` is `[p,p)`.

- [ ] **Step 1: Write the failing tests.**
  - `gap_table` (BED3 input):

    | input | flags | output |
    |---|---|---|
    | `[1,5) [5,8)` | `-d 0` | one cluster `[1,8)` |
    | `[1,5) [6,8)` | `-d 0` | two clusters |
    | `[1,5) [6,8)` | `-d 1` | one cluster |
    | `[1,5) [4,8)` | `-d -2` | two clusters |
    | `[1,5) [3,8)` | `-d -2` | one cluster |

  - `points`: a point at 5 inside `[1,8)` joins it; two points at 9 alone give `[9,9)`.
  - `per_strand`, `count_and_aggs` (`--count -c 5:max -c 4:distinct`), and `multiple_inputs_kmerged`.
  - `vcf_spans_heap`: the Task 5 A1/A2 records merge into `[101,105)`.
  - `unsorted_input`: exit 1 with the sort hint.
  - `determinism`, plus `merge_props.rs` (proptest against `naive_merge`, including negative gaps and `-s`).
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:** `bedtools merge` (default, `-d 10`, `-d -2`, `-s`, `-c 5 -o max`, `-c 4,4 -o distinct,count`); `bedops --merge a b` against `gloci merge a b`. Pin any differences. Commit the goldens.
- [ ] **Step 6: Commit** `feat: gloci merge (k-way, signed --max-gap, -s, --count, -c)`

### Task 12: `gloci subtract`

**Files:**
- Create: `crates/gloci-ops/src/subtract.rs`, `crates/gloci/src/commands/subtract.rs`
- Modify: `cli.rs`, the dispatch, `naive.rs` (`naive_subtract`)
- Test: `crates/gloci/tests/subtract.rs`, `crates/gloci-ops/tests/subtract_props.rs`; cases

**Interfaces:**
- Produces:
  - `pub struct SubtractOptions { pub whole: bool, pub sum_overlaps: bool, pub relate: RelateOptions }`
  - `pub fn subtract(inputs: RelateInputs, o: &SubtractOptions, w: &mut dyn RecordWriter) -> Result<SweepSummary, IoError>`
- CLI: `gloci subtract -a A -b B... [--whole] [--sum-overlaps]` with the same groups as intersect (spec §7.1).

Semantics (D44, D60, spec §4.6):
- **Default:** A's span minus the union of the accepted hits' windowed spans. Each remaining piece is written as A with rewritten coordinates, in ascending order. A point B removes nothing.
- A point A is dropped when it has any accepted hit, and kept otherwise.
- VCF/BCF A records are removed whole: exactly the records `intersect --mode none` would drop with the same criteria.
- `--whole`: drop A when it has any accepted hit.
- `--sum-overlaps`: drop A when the summed overlap lengths reach `-f` of A; with no `-f`, it acts like `--whole`.
- `--split`: removes B's blocks.

- [ ] **Step 1: Write the failing tests.**
  - `pieces`:
    - `[0,100)` minus `[10,20) [30,40)` → `[0,10) [20,30) [40,100)`;
    - a B covering A removes it;
    - a bookended B leaves A unchanged;
    - a B point inside A leaves A unchanged;
    - an A point inside B is removed.
  - `criteria`: `-f 0.5` per B; `-s`.
  - `whole_and_sum`: two B each covering 30% of A. `--sum-overlaps -f 0.5` removes A; `-f 0.5` alone removes nothing; `--whole` removes A.
  - `gff_pieces`: 1-based columns.
  - `vcf_a_whole_records`; `split_blocks`; `determinism`; and `subtract_props.rs` (proptest of piece arithmetic against `naive_subtract`).
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:** `bedtools subtract` (default, `-A`, `-N -f 0.5`, `-f 0.5`, `-s`); `bedops --difference a b` against `gloci merge a.bed | gloci subtract -a - -b b.bed` (the spec §11 mapping). Pin differences, and pin VCF A against bedtools under `#vcf-allele-span`.
- [ ] **Step 6: Commit** `feat: gloci subtract (pieces, --whole, --sum-overlaps, VCF whole records)`

### Task 13: `gloci closest`

**Files:**
- Create: `crates/gloci-ops/src/closest.rs`, `crates/gloci/src/commands/closest.rs`
- Modify: `cli.rs`, the dispatch, `naive.rs` (`naive_closest`)
- Test: `crates/gloci/tests/closest.rs`, `crates/gloci-ops/tests/closest_props.rs`; cases

**Interfaces:**
- Consumes: `Neighbours` (Task 6), `distance` and `signed_distance` (P0 core).
- Produces:

```rust
pub enum Ties { All, First, Last }
pub struct ClosestOptions {
    pub k: usize, pub ties: Ties, pub distance: bool, pub signed: Option<SignFrame>,
    pub ignore_overlaps: bool, pub ignore_upstream: bool, pub ignore_downstream: bool,
    pub max_distance: Option<u64>, pub each_side: bool, pub per_file: bool, pub relate: RelateOptions,
}
pub fn closest(inputs: RelateInputs, o: &ClosestOptions, w: &mut dyn RecordWriter) -> Result<SweepSummary, IoError>
```

- CLI: `gloci closest -a A -b B... [-k N] [--ties all|first|last] [--distance] [--signed ref|a|b] [--ignore-overlaps] [--ignore-upstream] [--ignore-downstream] [--max-distance N] [--each-side] [--per-file] [-s|-S]`, plus `AArg`, `BArgs`, `UnsortedArg` and the common groups. No fractions or windows.

Semantics:
- **Candidates:** accepted overlapping hits at distance 0 (dropped with `--ignore-overlaps`), plus `left`/`right`.
- **Up and downstream:** the sign of `signed_distance` in the `--signed` frame, `ref` when the flag is absent. `--ignore-upstream`/`--ignore-downstream` drop candidates by that sign.
- **Selection:** the `k` nearest. Ties at the k-th distance are all kept with `all`; with `first`/`last`, the first or last by (file, input order).
- **`--per-file`:** select per B file, and write each file's rows in file order.
- **Output (TSV):** A raw, then the label (D57), then B raw or the D45 placeholder, then the distance with `--distance` (signed when `--signed` is given, `-1` with no candidate).
- **`--each-side`:** one row per A: A, then the nearest upstream B (plus distance), then the nearest downstream B (plus distance). Overlapping B count on both sides unless `--ignore-overlaps`. `-k` applies per side.

- [ ] **Step 1: Write the failing tests.**
  - `distances`: overlap 0, bookended 1, a gap of 10 gives 11, and an edge point gives 1.
  - `ties` (all, first, last) and `k2`.
  - `signed_frames`: `--ignore-upstream --signed a` on a minus-strand A flips which side is dropped.
  - `no_candidate`: placeholder row with `-1`.
  - `each_side`, `per_file`, `max_distance`, and `strand_rules` (`-s` skips the nearer opposite-strand B).
  - `determinism`, plus `closest_props.rs` (proptest against `naive_closest`, sorted and `--unsorted`).
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:** `bedtools closest` (`-d`, `-D ref`, `-D a`, `-k 2 -d`, `-t first`, `-io -d`, `-iu -D ref`, `-id -D ref`, `-s -d`, and `-mdb each` for two files with `--per-file`); `closest-features --dist a b | tr '|' '\t'` against `gloci closest --each-side --distance`. Pin differences.
- [ ] **Step 6: Commit** `feat: gloci closest (-k, ties, signed distances, --each-side, --per-file)`

### Task 14: `gloci map`

**Files:**
- Create: `crates/gloci-ops/src/map.rs`, `crates/gloci/src/commands/map.rs`
- Modify: `cli.rs`, the dispatch, `naive.rs` (`naive_map`)
- Test: `crates/gloci/tests/map.rs`, `crates/gloci-ops/tests/map_props.rs`; cases

**Interfaces:**
- Produces:
  - `pub struct MapOptions { pub aggs: Vec<AggSpec>, pub skip_unmatched: bool, pub empty: String, pub delim: String, pub relate: RelateOptions }`
  - `pub fn map(inputs: RelateInputs, o: &MapOptions, w: &mut dyn RecordWriter) -> Result<SweepSummary, IoError>`
- CLI: `gloci map -a A -b B... -c COL:OP... [--skip-unmatched] [--empty S] [--delim S]`, plus the overlap, window and unsorted groups (spec §7.1).

Semantics:
- Output is TSV: A raw, then one column per aggregation over A's accepted hits, in hit order (D38).
- Columns are resolved against the first `-b` file's format and header. Every `-b` must be in that format family (D57).
- `--skip-unmatched` drops A with no accepted hit.
- The TSV column names are the aggregation labels.

- [ ] **Step 1: Write the failing tests.**
  - `sum_count`: an exact fixture.
  - `no_overlap`: `.` and `0`.
  - `window_range`: `-w 100`.
  - `fractions`: `-f 0.5`.
  - `vcf_b_info`: `-c INFO/DP:mean`.
  - `vcf_a_header_line`: the `#CHROM` line ends `\t4:sum`.
  - `bases_and_indicator`: `0:bases`, `0:indicator`.
  - `skip_unmatched`, `determinism`, plus `map_props.rs` (proptest of sum, count, min and max against `naive_map`).
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases.**
  - `bedtools map -c 4 -o sum|mean|median|min|max|count|collapse|distinct` against `-c 4:OP`, compared `numeric`.
  - `bedmap --echo --count|--sum|--mean a b | tr '|' '\t'`.
  - `bedmap --echo --range 100 --count` against `-w 100 -c 0:count`.
  - `bedmap --echo --bases` against `-c 0:bases`.

  Pin differences (for example, the stdev definition if bedtools uses the population form).
- [ ] **Step 6: Commit** `feat: gloci map (P1 core aggregation ops)`

### Task 15: `gloci coverage`

**Files:**
- Create: `crates/gloci-ops/src/coverage.rs`, `crates/gloci/src/commands/coverage.rs`
- Modify: `cli.rs`, the dispatch, `naive.rs` (`naive_depth`)
- Test: `crates/gloci/tests/coverage.rs`, `crates/gloci-ops/tests/coverage_props.rs`; cases

**Interfaces:**
- Produces:
  - `pub enum CoverageReport { Summary, Hist, PerBase, Mean }`
  - `pub struct CoverageOptions { pub report: CoverageReport, pub relate: RelateOptions }`
  - `pub struct CoverageTotals { /* depth -> bases, total length; summed across chunks by Task 24 */ }`
  - `pub fn coverage(inputs: RelateInputs, o: &CoverageOptions, w: &mut dyn RecordWriter) -> Result<(SweepSummary, CoverageTotals), IoError>`
  - `pub fn write_hist_totals(t: &CoverageTotals, w: &mut dyn RecordWriter) -> Result<(), IoError>`
- CLI: `gloci coverage -a A -b B... [--hist | --per-base | --mean]`, plus `OverlapArgs` (no windows), `UnsortedArg` and the common groups.

Semantics:
- Each accepted hit, clipped to A, adds 1 to depth over its bases.
- **Summary:** A raw, count of accepted hits, bases with depth > 0, A length, fraction (bases ÷ length).
- **`--hist`:** for each depth present (including 0): A raw, depth, bases at that depth, A length, fraction. After all A, `all` rows: depth, bases, total length, fraction.
- **`--per-base`:** A raw, 1-based offset within A, depth.
- **`--mean`:** A raw, total depth ÷ length.
- **More than one `-b`:** A raw plus one covered fraction per file; with `--hist`, `--per-base` or `--mean` this is a usage error.
- Point A follows D56. Numbers follow D46 formatting.

- [ ] **Step 1: Write the failing tests:** an exact `summary` fixture, `hist_with_all`, `per_base`, `mean`, `annotate_multi_b`, `points`, `determinism`, plus `coverage_props.rs` (proptest of depth profiles against `naive_depth`).
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases**, compared `numeric`: `bedtools coverage` (default, `-hist`, `-d`, `-mean`, `-s`, `-f 0.5`); `bedtools annotate -i a -files b1 b2` against `gloci coverage -a a -b b1 b2`. Pin differences.
- [ ] **Step 6: Commit** `feat: gloci coverage (--hist, --per-base, --mean, multi -b fractions)`

### Task 16: `gloci complement`

**Files:**
- Create: `crates/gloci-ops/src/complement.rs`, `crates/gloci/src/commands/complement.rs`
- Modify: `cli.rs`, the dispatch, `naive.rs` (`naive_complement`)
- Test: `crates/gloci/tests/complement.rs`; cases

**Interfaces:**
- Consumes: `kmerge` (Task 11).
- Produces: `pub fn complement(inputs: Vec<Box<dyn RecordReader>>, chroms: &mut ChromTable, genome: Option<&Genome>, w: &mut dyn RecordWriter) -> Result<(), IoError>`.
- CLI: `gloci complement FILES... [-g NAME|FILE]`, plus chrom-order, vcf-span, format, threads and output groups.

Semantics (spec §7.2, D48):
- **With `-g`:** every genome chromosome in genome order, with its uncovered gaps over `[0, length)`, including chromosomes that have no features.
  - A genome without lengths (for example a VCF header without `length=`) is a usage error.
  - A feature past the chromosome end exits 1 naming the line.
- **Without `-g`:** only gaps between merged features on the same chromosome.
- Points are ignored. The output is BED3.

**Also in this task (D39, revised):** add `--sort`, with `--max-mem SIZE` (default 2G) and `TmpDirArg`, to both `merge` (Task 11) and `complement`. Put the sorting step in one reusable helper, because Task 16b uses it for the relation commands. With `--sort`, all inputs go through the P0 external sorter (`gloci_ops::Sorter`, same chromosome table and tie rules as `gloci sort`). The sorted records are then fed to the merge or complement logic through the sorter's ordered output (a callback, or a temp file in `--tmp-dir`), so memory stays within `--max-mem`. The output must equal `gloci sort INPUTS | gloci merge -` (or `complement`) byte for byte. Without the flag, the order-error hint names `--sort`. Tests: `merge_sort_equals_sort_pipe`, `complement_sort_equals_sort_pipe` (shuffled inputs, several files), and `sort_spills_within_max_mem` (a small `--max-mem` forces a spill, and the output is unchanged).

- [ ] **Step 1: Write the failing tests:** `with_genome` (an empty chromosome becomes one full-length row), `without_genome` (no leading or trailing gap), `points_ignored`, `bookended_no_gap`, `builtin_genome` (`-g sacCer3` gives 17 rows for a file with one feature on chrI), `past_end_error`, `determinism`, and a proptest against `naive_complement`.
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:** `bedtools complement -i a -g g.genome` (with gloci `-g g.genome`); `bedops --complement a` against `gloci complement a`. Pin differences.
- [ ] **Step 6: Commit** `feat: gloci complement (with and without -g)`

### Task 16b: `--sort` on the relation commands

**Files:**
- Modify: `crates/gloci-ops/src/relate/mod.rs`, `crates/gloci/src/cli.rs` (a `SortArg` group: `--sort`, `--max-mem`, `TmpDirArg`; `--sort` conflicts with `--unsorted`), and the five command files
- Test: `crates/gloci/tests/relation_sort.rs`

**Interfaces:**
- Consumes: Task 16's sorting helper, and `run_relation` (Task 8).
- Produces: `RelateOptions.sort: Option<SortSettings { max_mem: u64, tmp_dir: PathBuf }>`. When it is set, `run_relation` sorts A and each B with the helper before the sorted sweep. The visitor and the output code are unchanged, and the order checks run on the sorted streams as usual.

Semantics (D39): results are in sorted order and equal `gloci sort A > a.s; gloci sort B > b.s; gloci <cmd> -a a.s -b b.s …` byte for byte. Memory stays within `--max-mem` per sort. Each command's `--help` says when to use `--sort` (inputs too large for memory, sorted output wanted) and when to use `--unsorted` (keep A's order, B fits in memory).

- [ ] **Step 1: Write the failing tests.**
  - `sort_equals_presorted_pipeline` for intersect (any, pairs), subtract, closest, map and coverage, on shuffled A and B;
  - `sort_and_unsorted_conflict` (exit 2);
  - `sort_spills_within_max_mem`;
  - `order_error_hint_names_both_flags`.
- [ ] **Step 2: Run** `cargo test -p gloci --test relation_sort`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Commit** `feat: --sort on relation commands (external sort, then the sorted sweep)`

### Task 17: `gloci slop` and `gloci flank`

**Files:**
- Create: `crates/gloci-ops/src/transform.rs`, `crates/gloci/src/commands/{slop,flank}.rs`
- Modify: `cli.rs`, the dispatch
- Test: `crates/gloci/tests/transform.rs`; cases

**Interfaces:**
- Produces:
  - `pub enum Amount { Bases(u64), Fraction(Fraction) }` (`--pct`: an exact decimal fraction of the feature length, floored)
  - `pub struct Pad { pub left: Amount, pub right: Amount, pub strand_aware: bool }`
  - `pub fn slop(input: Box<dyn RecordReader>, chroms: &mut ChromTable, pad: &Pad, genome: Option<&Genome>, w: &mut dyn RecordWriter) -> Result<(), IoError>`
  - `pub fn flank(input: Box<dyn RecordReader>, chroms: &mut ChromTable, pad: &Pad, genome: Option<&Genome>, w: &mut dyn RecordWriter) -> Result<(), IoError>`
- CLI (both commands): `-b/--both N`, `--left N`, `--right N`, `-s`, `--pct`, `-g NAME|FILE`, plus format, threads and output groups. `-b` conflicts with `--left`/`--right`, and at least one amount is required.

Semantics (D44, D49):
- **slop:** `[start - L, end + R)`, clipped at 0, and at the chromosome length when `-g` is given. A chromosome missing from `-g` is an error.
- **flank:** a left flank `[start - L, start)` and a right flank `[end, end + R)`, clipped the same way, left flank first.
- `-s` swaps L and R for `-` strand records.
- Output keeps the input format.
- Order: a warning when the output goes out of order. That can happen with `-s` or `--pct` in slop, and with any flank.
- Zero-length flanks follow whatever bedtools does; the differential case decides, and the outcome is documented.

- [ ] **Step 1: Write the failing tests.**
  - `slop_table`:
    - `-b 10` on `[100,200)` → `[90,210)`;
    - clip at 0;
    - clip at the chromosome end with `-g`;
    - `-s` on a minus-strand record with `--left 5 --right 20` → `[80,205)`;
    - `--pct -b 0.5` on `[1000,2000)` → `[500,2500)`.
  - `flank_table`: `-b 10`, `--left 10` only.
  - Also `gff_rewrite`, `vcf_usage_error`, `order_warning_once`, `determinism`.
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:** `bedtools slop -b 10 -g g`, `-l 5 -r 20 -s -g g`, `-pct -b 0.1 -g g`; `bedtools flank -b 10 -g g`, `-l 5 -s -g g`. Pin differences.
- [ ] **Step 6: Commit** `feat: gloci slop and gloci flank`

### Task 18: `gloci windows`

**Files:**
- Create: `crates/gloci-ops/src/windows.rs`, `crates/gloci/src/commands/windows.rs`
- Modify: `cli.rs`, the dispatch
- Test: `crates/gloci/tests/windows.rs`; cases

**Interfaces:**
- Produces:
  - `pub enum Split { Size { size: u64, step: u64 }, Count(u64) }`
  - `pub struct WindowsOptions { pub split: Split, pub reverse: bool, pub id_column: bool }`
  - `pub fn windows_genome(genome: &Genome, o: &WindowsOptions, w: &mut dyn RecordWriter) -> Result<(), IoError>`
  - `pub fn windows_intervals(input: Box<dyn RecordReader>, chroms: &mut ChromTable, o: &WindowsOptions, w: &mut dyn RecordWriter) -> Result<(), IoError>`
- CLI: `gloci windows (-g NAME|FILE | -i FILE) (-w SIZE [--step S] | -n COUNT) [--reverse] [--id-column]`. `-g` and `-i` are mutually exclusive and one is required; the same goes for `-w` and `-n`. `--step` defaults to SIZE.

Semantics (bedtools `makewindows`):
- Windows are cut from each source interval (each chromosome with `-g`, each footprint with `-i`). The last window is truncated at the interval end.
- `-n` splits each interval into COUNT windows, with the remainder distributed the way bedtools does it (the differential case decides).
- `--id-column` appends a 1-based window number per source interval; `--reverse` numbers from the end.
- Points in `-i` give no windows. The output is BED3, plus the id column.

- [ ] **Step 1: Write the failing tests.**
  - `size_step`: `-w 100 --step 50` on a 250 bp chromosome gives `[0,100) [50,150) [100,200) [150,250) [200,250)`.
  - Also `count`, `id_and_reverse`, `intervals_input`, `points_skipped`, `determinism`.
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:** `bedtools makewindows -g g -w 1000 -s 500`, `-g g -n 3`, `-b src.bed -w 100 -i winnum`, `-b src.bed -n 4 -i winnum -reverse`; `bedops --chop 100 a.bed` against `gloci merge a.bed | gloci windows -i - -w 100`. Pin differences.
- [ ] **Step 6: Commit** `feat: gloci windows`

### Task 19: `gloci convert`

**Files:**
- Create: `crates/gloci-ops/src/convert/{mod.rs, vcf.rs, gff.rs, bed12.rs}`, `crates/gloci/src/commands/convert.rs`
- Modify: `cli.rs`, the dispatch
- Test: `crates/gloci/tests/convert.rs`; cases

**Interfaces:**
- Consumes: `alt_kinds`, `alt_footprints` (core), `bed12_blocks` (io), and the P0 BCF→VCF text rendering in `gloci-io/src/convert.rs`.
- Produces:
  - `pub struct ConvertOptions { pub split_alleles: bool, pub types: Option<Vec<VariantKind>>, pub features: Vec<String>, pub name_attr: Option<String>, pub blocks: bool }`
  - `pub fn convert(input: Box<dyn RecordReader>, chroms: &mut ChromTable, o: &ConvertOptions, w: &mut dyn RecordWriter) -> Result<(), IoError>`
- CLI: `gloci convert FILE [--split-alleles] [--types LIST] [--feature NAME]... [--name-attr KEY] [--blocks]`, plus vcf-span, format, threads and output groups. The output is BED.

Semantics (spec §7.4):
- **VCF/BCF → BED:** chrom, footprint start, footprint end, ID, REF, ALT, QUAL, FILTER.
  - `--split-alleles` gives one row per ALT, with that ALT's footprint and ALT value. `*`, `<*>` and `<NON_REF>` rows are skipped when other ALTs exist.
  - `--types` follows the D50 type-set rule (per ALT with `--split-alleles`).
- **GFF3/GTF → BED**, 10 columns (bedops `gff2bed` layout): chrom, start−1, end, name, score, strand, source, type, phase, attributes.
  - The name comes from `--name-attr`, defaulting to `ID` for GFF3 and `gene_id` for GTF; it is `.` when absent.
  - `--feature` filters column 3 (repeatable or comma-separated).
- **BED12 → BED6** with `--blocks`: one row per block, keeping name, score and strand.
- Any other input/flag combination is a usage error saying what `convert` does for that input.
- If the output goes out of order (allele footprints), there is one warning (D49).

- [ ] **Step 1: Write the failing tests.**
  - `vcf_rows`: the §4.2 table rows in `allele` and `ref` modes.
  - `split_alleles`: `chr1 100 . ACGT A,AC` gives `chr1 100 103 . ACGT A …` and `chr1 101 103 . ACGT AC …`.
  - Also `types_del`, `gff3_id_and_name_attr`, `gtf_gene_id_default`, `feature_filter`, `bed12_blocks`, `usage_errors`, `determinism`.
- [ ] **Step 2: Run** the tests. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases:**
  - `vcf2bed < in.vcf | cut -f1-3` against `gloci convert in.vcf --vcf-span ref | cut -f1-3`;
  - `vcf2bed --deletions` against `--types del`;
  - `gff2bed < in.gff3` against `gloci convert in.gff3`;
  - `gtf2bed < in.gtf` against `gloci convert in.gtf`;
  - `bedtools bed12tobed6 -i in.bed` against `gloci convert in.bed --blocks`.

  Pin differences, including allele spans against `vcf2bed` under `#vcf-allele-span`.
- [ ] **Step 6: Commit** `feat: gloci convert (VCF/GFF/GTF to BED, BED12 blocks)`

### Task 20: `gloci-vcf`: genotypes and site statistics

**Files:**
- Create: `crates/gloci-vcf/Cargo.toml` and `src/{lib.rs, gt.rs, bcf_fmt.rs, stats.rs, hwe.rs}`. Add it to the workspace.
- Test: `crates/gloci-vcf/tests/{gt.rs, stats.rs, hwe.rs}`

**Interfaces:**
- Consumes: `BcfDicts` and the typed-value decoding in `gloci-io/src/bcf.rs`. Expose the decoding helpers `pub(crate)` → `pub` where needed, keeping the names.
- Produces:

```rust
pub struct Genotype { pub alleles: Vec<Option<u32>>, pub phased: bool }   // None = missing allele; reused across samples
pub struct SampleFields { pub gt: Option<Genotype>, pub dp: Option<u32>, pub gq: Option<u32>, pub ft_pass: bool }
pub trait GenotypeSource { fn n_samples(&self) -> usize; fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String>; }
pub fn text_genotypes(line: &[u8]) -> Result<impl GenotypeSource + '_, String>      // FORMAT keys located once per line
pub fn bcf_genotypes<'a>(raw: &'a [u8], dicts: &'a BcfDicts) -> Result<impl GenotypeSource + 'a, String>
pub struct MaskRules { pub min_gq: Option<u32>, pub min_dp: Option<u32>, pub max_dp: Option<u32>, pub mask_filtered: bool }
pub struct SiteStats {
    pub allele_counts: Vec<u64>, pub called_alleles: u64, pub samples: u64, pub missing: u64,
    pub dp_sum: u64, pub dp_n: u64, pub hwe_counts: Option<(u64, u64, u64)>, pub all_phased: bool,
}
pub fn site_stats(src: &dyn GenotypeSource, keep: Option<&[bool]>, mask: &MaskRules, n_alleles: usize, masked: &mut Vec<bool>) -> Result<SiteStats, String>
impl SiteStats { pub fn maf(&self) -> Option<f64>; pub fn mac(&self) -> Option<u64>; pub fn missing_fraction(&self) -> f64; pub fn mean_dp(&self) -> Option<f64>; }
pub fn hwe_p(het: u64, hom1: u64, hom2: u64) -> f64     // Wigginton 2005 exact test, two-sided, capped at 1
```

Rules (D50):
- MAF is the second-highest allele count ÷ called alleles; MAC is that count.
- Mean DP is taken over selected, unmasked samples that have a DP.
- `hwe_counts` covers biallelic sites, counting fully called diploid genotypes as (0/0, 0/1, 1/1). It is `None` otherwise, and then `--hwe` passes the site.
- `all_phased` is false if any called genotype with ploidy ≥ 2 is unphased.
- The text path allocates nothing per sample; it iterates fields in place.
- `masked[i]` reports which samples were masked, for Task 22.

- [ ] **Step 1: Write the failing tests.**
  - `gt_text`:
    - `0/1` and `1|0` (phased);
    - `./.` (missing), and `./1` (missing because one allele is missing; allele 1 still counts);
    - haploid `2`, triploid `0/1/2`;
    - FORMAT `DP:GT`, and a sample column shorter than FORMAT (the trailing fields are missing).
  - `gt_bcf_equals_text`: a 4-sample record encoded with `testkit::vcf_to_bcf` gives `SampleFields` identical to the text path.
  - `stats`: counts `[7,3]` give MAF `0.3` and MAC `3`; counts `[5,3,2]` give MAF `0.3`.
  - `masking`: `min_dp 10` masks DP 5; missing DP is not masked; `mask_filtered` masks FT `LowQual` and not `PASS` or `.`.
  - `hwe`:
    - `hwe_p(0, 1, 1)` is 1/3 within 1e-12;
    - `hwe_p(2, 0, 0)` is `1.0`;
    - `hwe_p(0, 0, 0)` is `1.0`;
    - `hwe_p(57, 14, 50)` equals, to 1e-10, a brute-force enumeration written in the test: exact genotype-configuration probabilities from log-factorials over every possible het count, summing those no larger than the observed one.
- [ ] **Step 2: Run** `cargo test -p gloci-vcf`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p gloci-vcf`. Expected: pass.
- [ ] **Step 5: Commit** `feat(vcf): genotype decoding (VCF text and BCF), site statistics, HWE exact test`

### Task 21: `vcf filter` site filters

**Files:**
- Create: `crates/gloci-vcf/src/filter/{mod.rs, site.rs, lists.rs}`
- Test: `crates/gloci-vcf/tests/site_filters.rs`

**Interfaces:**
- Consumes: Task 20 and core `alt_kinds`.
- Produces:

```rust
pub enum VcfType { Snv, Mnv, Ins, Del, Indel, Sv, Ref, Other }
pub struct SiteFilters {
    pub chroms: Option<HashSet<Vec<u8>>>, pub exclude_chroms: HashSet<Vec<u8>>,
    pub positions: Option<PositionSet>, pub exclude_positions: Option<PositionSet>,
    pub ids: Option<HashSet<Vec<u8>>>, pub exclude_ids: Option<HashSet<Vec<u8>>>,
    pub types: Option<Vec<VcfType>>, pub exclude_types: Option<Vec<VcfType>>,
    pub min_alleles: Option<usize>, pub max_alleles: Option<usize>,
    pub min_qual: Option<f64>, pub pass_only: bool, pub keep_filter: Option<Vec<Vec<u8>>>, pub remove_filter: Option<Vec<Vec<u8>>>,
    pub require_info: Vec<Vec<u8>>, pub exclude_info: Vec<Vec<u8>>,
    pub min_maf: Option<f64>, pub max_maf: Option<f64>, pub min_mac: Option<u64>, pub max_mac: Option<u64>,
    pub max_missing: Option<f64>, pub max_missing_count: Option<u64>,
    pub min_mean_dp: Option<f64>, pub max_mean_dp: Option<f64>,
    pub hwe: Option<f64>, pub thin: Option<u64>, pub phased_only: bool,
    pub exclude_regions: Option<HashMap<Vec<u8>, RegionSet>>, pub inclusive_points: bool,
}
pub type PositionSet = HashSet<(Vec<u8>, u64)>;   // (chrom, 1-based POS)
pub struct SiteView<'a> { pub chrom: &'a [u8], pub pos: u64, pub id: &'a [u8], pub ref_: &'a [u8], pub alts: &'a [u8], pub qual: Option<f64>, pub filter: &'a [u8], pub info: &'a [u8], pub span: Interval }
impl SiteFilters { pub fn needs_genotypes(&self) -> bool; pub fn static_pass(&self, s: &SiteView) -> bool; pub fn stats_pass(&self, st: &SiteStats) -> bool; }
pub struct Thinner { /* last kept POS per chromosome */ } impl Thinner { pub fn keep(&mut self, chrom: ChromId, pos: u64) -> bool }
pub fn read_positions(path: &str) -> Result<PositionSet, IoError>      // "chrom<TAB>POS" (1-based), '#' comments
pub fn read_ids(path: &str) -> Result<HashSet<Vec<u8>>, IoError>      // one ID per line
```

The rules are D50 and spec §7.8, plus these:
- A record's IDs are split on `;`, and it matches when any ID is listed.
- `--require-info`/`--exclude-info` test whether the key is present, as a flag or as `KEY=…`.
- `--max-missing` is the missing fraction among selected samples after masking, and it passes when ≤ X.
- `--thin` keeps a site when it is the first on its chromosome, or when `POS - last_kept_POS >= N`.

- [ ] **Step 1: Write the failing tests**, one per filter on a small fixture with exact kept/dropped sets:
  - `pass_only_keeps_dot`, `missing_qual_fails_min_qual`;
  - `keep_and_remove_filter`, `types_type_set`, `alleles_counts`;
  - `info_flags`, `positions_and_ids`;
  - `maf_mac_missing_meandp_hwe_phased` (stats-based);
  - `thin_first_wins`, `exclude_regions_points`.
- [ ] **Step 2: Run** `cargo test -p gloci-vcf --test site_filters`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p gloci-vcf`. Expected: pass.
- [ ] **Step 5: Commit** `feat(vcf): site filters for vcf filter`

### Task 22: `vcf filter` samples, masking and record output

**Files:**
- Create: `crates/gloci-vcf/src/filter/{samples.rs, rewrite.rs}`
- Test: `crates/gloci-vcf/tests/rewrite.rs`

**Interfaces:**
- Consumes: Task 20 (`MaskRules`, `SiteStats`, the masked flags), the P0 VCF↔BCF conversion in `gloci-io/src/convert.rs`, and `VcfHeader`.
- Produces:

```rust
pub struct SampleSelection { pub keep: Vec<bool> }
pub fn select_samples(header_samples: &[Vec<u8>], include: Option<&[String]>, exclude: &[String]) -> Result<SampleSelection, UsageError>
pub struct RewritePlan { pub selection: Option<SampleSelection>, pub mask: MaskRules, pub keep_info: Option<Vec<Vec<u8>>>, pub drop_info: Vec<Vec<u8>>, pub recalc_ac: bool }
impl RewritePlan { pub fn is_passthrough(&self) -> bool }
pub fn rewrite_text(line: &[u8], plan: &RewritePlan, masked: &[bool], stats: Option<&SiteStats>, out: &mut Vec<u8>) -> Result<(), String>
pub fn rewrite_header(h: &VcfHeader, plan: &RewritePlan) -> Vec<u8>
pub fn fmt_g6(x: f64) -> String                      // C "%g": 6 significant digits, trailing zeros trimmed
```

Rules:
- An unknown sample name is a usage error naming the first unknown one, together with the header's first five samples. `--samples-file` and `--exclude-samples-file` are read one name per line.
- **Passthrough** (no selection change, no masking, no INFO pruning, no recalc): the raw record bytes are written unchanged. BCF raw bytes go to BCF output unchanged.
- **Text rewrite:**
  - Unselected sample columns are dropped.
  - A masked genotype's alleles become `.`, keeping the ploidy and separators (`./.`, `.|.`, `.`). Its other FORMAT fields are kept.
  - INFO pruning removes the keys, and `rewrite_header` removes their `##INFO` lines.
  - `--recalc-ac` writes AC (one per ALT), AN and AF (`fmt_g6`), replacing them where they exist and otherwise appending them in that order. An INFO of `.` is replaced.
- **BCF with a non-passthrough plan:** render the record to VCF text (P0 converter), rewrite it, then re-encode it for BCF output. This slower path applies only to rewrites.

- [ ] **Step 1: Write the failing tests.**
  - `subset_text`, `unknown_sample_is_usage_error`.
  - `mask_text`: `0/1` → `./.`, `1|1` → `.|.`, haploid `1` → `.`.
  - `drop_info_and_header`, `keep_info`.
  - `recalc_ac`: samples `0/1` and `1/1` give `AC=3;AN=4;AF=0.75`.
  - `fmt_g6`: `0.333333333` → `0.333333`, `0.5` → `0.5`, `1e-7` → `1e-07`.
  - `bcf_rewrite_equals_text`: BCF in and BCF out with a sample subset, decoded back to text, equals the text path.
  - `bcf_passthrough_bytes`: the records are byte-identical.
- [ ] **Step 2: Run** `cargo test -p gloci-vcf --test rewrite`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p gloci-vcf`. Expected: pass.
- [ ] **Step 5: Commit** `feat(vcf): sample selection, genotype masking, INFO pruning and --recalc-ac`

### Task 23: `gloci vcf filter`

**Files:**
- Create: `crates/gloci-ops/src/vcf_filter.rs`, `crates/gloci/src/commands/vcf.rs`
- Modify: `cli.rs` (`Vcf(VcfArgs)` with `VcfCommand::Filter(VcfFilterArgs)`), the dispatch, `crates/gloci-testkit/src/differential.rs` and `xtask/src/diff.rs` (compare mode `vcf-records`), `docs/differences.md`
- Test: `crates/gloci/tests/vcf_filter.rs`; cases

**Interfaces:**
- Consumes: Tasks 20–22, and the P0 view region selection (`open_view`/the `Selection` reader) for `--region`/`--regions-file`.
- Produces:
  - `pub struct VcfFilterOptions { pub include_samples: Option<Vec<String>>, pub exclude_samples: Vec<String>, pub mask: MaskRules, pub sites: SiteFilters, pub plan: RewritePlan, pub regions: Vec<String>, pub regions_file: Option<String> }`
  - `pub struct FilterSummary { pub records_in: u64, pub records_out: u64 }`
  - `pub fn vcf_filter(input: Box<dyn RecordReader>, chroms: &mut ChromTable, o: &VcfFilterOptions, w: &mut dyn RecordWriter) -> Result<FilterSummary, IoError>`
- CLI: every flag in spec §7.8 steps 1–4 except `--include`/`--exclude` (P2), plus `--keep-info`, `--drop-info`, `--recalc-ac`, `--inclusive-points`, and the vcf-span, format, threads, skip-invalid and output groups. The output is VCF or BCF from `-o`/`-O` (converting between them is allowed, as in `view`), and `--write-index` works. Help examples include `gloci vcf filter cohort.bcf --pass-only --min-qual 30 --min-maf 0.01 --max-missing 0.1 -o filtered.bcf`.

Pipeline order (spec §7.8):
1. sample selection;
2. masking;
3. statistics recomputed (only when `needs_genotypes()` or `--recalc-ac`);
4. static site filters, then stats filters, ANDed, with `--thin` last.

Regions:
- `--region`/`--regions-file` use the index when one exists and scan otherwise (P0 view behaviour).
- `--exclude-regions-file`, `--chrom` and `--exclude-chrom` filter as records stream past.
- `--thin` requires sorted input; it is order-checked.

The `vcf-records` compare mode compares the data lines and the `#CHROM` line, ignoring `##` lines.

Validation (usage errors): `--biallelic` with `--min-alleles`/`--max-alleles`; MAF and missing fractions outside `[0, 1]`; `--samples` with `--samples-file`.

- [ ] **Step 1: Write the failing tests.**
  - `end_to_end_categories`: one case per §7.8 group.
  - `pipeline_order`: a site whose ALT carriers all have DP < 10, run with `--min-dp 10 --min-mac 1`, is dropped.
  - `subset_then_maf`: MAF is recomputed on `--samples`.
  - `region_with_index_and_scan`.
  - `bcf_passthrough`: the output records are byte-identical to the input records.
  - `stdin_input`, `determinism`, `usage_errors`.
- [ ] **Step 2: Run** `cargo test -p gloci --test vcf_filter`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Add differential cases** against vcftools 0.1.16, compared `vcf-records`. The reference is `vcftools --vcf input.vcf --recode --recode-INFO-all --stdout FLAGS`.

  | vcftools | gloci |
  |---|---|
  | `--minQ 30` | `--min-qual 30` |
  | `--remove-filtered-all` | `--pass-only` |
  | `--maf 0.1` | `--min-maf 0.1` |
  | `--max-maf 0.4` | `--max-maf 0.4` |
  | `--mac 2` | `--min-mac 2` |
  | `--min-alleles 2 --max-alleles 2` | `--biallelic` |
  | `--max-missing 0.9` | `--max-missing 0.1` |
  | `--min-meanDP 10` | `--min-mean-dp 10` |
  | `--hwe 0.01` | `--hwe 0.01` |
  | `--thin 1000` | `--thin 1000` |
  | `--keep-only-indels` | `--types indel` |
  | `--remove-indels` | `--exclude-types indel` |
  | `--chr chr2 --from-bp 100 --to-bp 5000` | `--region chr2:100-5000` |
  | `--bed r.bed` | `--regions-file r.bed` |
  | `--positions p.txt` | `--positions p.txt` |
  | `--snps ids.txt` | `--ids ids.txt` |
  | `--indv s1 --indv s3` | `--samples s1,s3` |
  | `--minDP 10` | `--min-dp 10` |
  | `--phased` | `--phased-only` |

  Pin each difference, at least these:
  - `#info-kept-by-default`;
  - `#max-missing-inverted` (documentation; the case shows the equivalence);
  - `#maf-multiallelic`;
  - `#missing-genotype-definition`, if vcftools differs;
  - masking output, for example how vcftools writes masked genotypes.

  Add a `bcftools view` agreement case wherever vcftools differs.
- [ ] **Step 6: Commit** `feat: gloci vcf filter (selection pipeline without expressions)`

### Task 24: Per-chromosome parallelism for indexed inputs

**Files:**
- Create: `crates/gloci-engine/src/{parallel.rs, chunk.rs}`
- Modify:
  - `crates/gloci-ops/src/relate/mod.rs` (`run_relation` chooses the mode);
  - the five relation ops, which already write through `&mut dyn RecordWriter`, so they run unchanged per chunk;
  - `cli.rs`: add `TmpDirArg` to intersect, subtract, closest, map and coverage.
- Test: `crates/gloci-engine/tests/parallel.rs`, `crates/gloci/tests/parallel_relations.rs`

**Interfaces:**
- Consumes: `open_region_reader`, `find_index` and `index_reference_names` (P0), and the P0 spill (LZ4) and cleanup registry.
- Produces:

```rust
pub struct ChunkSink { /* RecordWriter: first 64 MiB in memory, the rest LZ4-spilled to tmp_dir */ }
impl RecordWriter for ChunkSink { … }
pub fn all_indexed(paths: &[&str]) -> bool        // every path a regular BGZF file with .tbi/.csi; "-" → false
pub fn run_per_chrom<R: Send>(
    chroms: Vec<Vec<u8>>, threads: NonZeroUsize, tmp_dir: &Path,
    work: &(dyn Fn(&[u8], &mut ChunkSink) -> Result<R, IoError> + Sync),
    out: &mut dyn RecordWriter,
) -> Result<Vec<R>, IoError>                       // results in chromosome order
```

Per D41:
- **Chromosome list:** A's index reference names, ordered by the relation chromosome table. Names missing from a closed table are an error, as in streaming.
- **Each chunk:** open whole-chromosome region readers for A and every B, then run `sweep_sorted` with the command's visitor into the chunk sink.
- The main thread writes the header first, then appends chunks in order as each completes. At most `threads` chunks are in flight.
- The `coverage --hist` `all` rows come from the summed `CoverageTotals` and are written after the last chunk.
- Spill files are registered for cleanup. The first chunk error stops all work, removes the spills, and publishes no output.
- A missing index on any input falls back to the streaming path silently. That is not an error, because the results are identical.

- [ ] **Step 1: Write the failing tests.**
  - `parallel_equals_streaming`: for each of the 5 commands on indexed fixtures with 6 chromosomes, the output bytes match `-t 1` streaming (inputs copied without an index), and `-t 1`, `-t 2` and `-t 8` all match.
  - `chunk_spill_path`: a test-only constructor sets the in-memory limit to 1 KiB; the output is identical.
  - `chunk_error_cleans_up`: a malformed record in chr2 exits 1, leaves no output file, and leaves the tmp dir empty.
  - `missing_index_falls_back`.
- [ ] **Step 2: Run** `cargo test -p gloci-engine --test parallel` and `cargo test -p gloci --test parallel_relations`. They fail.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Commit** `feat(engine): per-chromosome parallel relations for indexed inputs`

### Task 25: Chromosome naming checks

**Files:**
- Create: `crates/gloci-io/src/naming.rs`
- Modify: `crates/gloci-ops/src/relate/mod.rs`, `crates/gloci-ops/src/merge.rs`, `crates/gloci-ops/src/complement.rs` (the end-of-run check); `crates/gloci-io/src/builtin.rs` (alias table); `crates/gloci-io/src/record.rs` (the `intern_chrom` hint)
- Test: `crates/gloci/tests/naming.rs`

**Interfaces:**
- Produces:
  - `pub fn chr_prefix_hint(a_names: &[&[u8]], b_names: &[&[u8]]) -> Option<String>`
  - `pub fn builtin_alias(id: &str) -> Option<&'static str>`, using the pairs `hg19/GRCh37`, `hg38/GRCh38`, `mm10/GRCm38`, `mm39/GRCm39`, `danRer11/GRCz11`, `dm6/BDGP6`, `ce11/WBcel235`, `sacCer3/R64`, matched case-insensitively, in both directions.
- `SweepSummary.shared_chroms` (Task 5) drives the relation warning. `merge` and `complement` track names per input.

Messages (D58):
- The no-shared-names warning: `warning: a.bed and b.bed share no chromosome names (a.bed: chr1, chr2, …; b.bed: 1, 2, …)`. When the `chr` prefix explains the mismatch, it adds `; the names differ by a "chr" prefix — rename one file's chromosomes so they match`.
- The `-g` hint, appended to the P0 unknown-chromosome error when the built-in alias genome lists the name: `; did you mean -g GRCh38?`.

- [ ] **Step 1: Write the failing tests.**
  - `chr_prefix_mismatch_warns`: an intersect of `chr1`-named A against `1`-named B exits 0 with empty output and the warning on stderr.
  - `disjoint_names_without_prefix`: the warning without the prefix clause.
  - `g_alias_hint`: `gloci sort x.bed -g GRCh38` on `chr1` records fails with `did you mean -g hg38?`.
  - `no_warning_when_shared`.
- [ ] **Step 2: Run** `cargo test -p gloci --test naming`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Commit** `feat: warn on disjoint chromosome names; suggest the other naming style for -g`

### Task 26: P0 follow-ups

**Files:** as each item names. Tests sit next to the code each item touches.

Items, chosen from `docs/superpowers/plans/2026-11-06-gloci-p0-followups.md`. Each one gets a test:
1. **`sort` with many inputs.** With more inputs than `MAX_FAN_IN` (128), sort merges in rounds through spill runs, so 300 inputs never hold 300 open files. Merge input buffers count against `--max-mem`. Test: `sort_300_inputs` under `ulimit -n 200` (Unix only), with output equal to the naive answer.
2. **Region strings.** `chr1:1,,0-5` (malformed comma groups) is an error with a hint.
3. **Header prefixes.** `@SQX` is not a `@SQ` line, and `#CHROMX` is not the `#CHROM` line. `.genome` files skip `#` comment lines and reject empty names.
4. **Direct-write outputs.** These open with `O_TRUNC` when the target behind `/dev/fd/N` is a regular file. Test: `3<>file` holding old longer content gives exactly the new bytes.
5. **Readers.** A final line ending in `\r` with no `\n` loses the `\r`. `read_batch(0)` panics under `debug_assert!` (document it as a contract).
6. **BCF reader.** A zero-count typed value (a Flag-typed END) counts as absent. BCF dictionaries are stored in the reader, which removes the per-record `unreachable!()`.
7. **CI.** The nightly corpus copy is `cp -n … || true`, and the fuzz corpus is cached between nights (actions/cache). The `ci.yml` static check writes `readelf` output to a file before grepping.

- [ ] **Step 1: Write the failing tests** for items 1–6 (item 7 is CI configuration).
- [ ] **Step 2: Run** them. They fail.
- [ ] **Step 3: Implement** every item.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: pass.
- [ ] **Step 5: Commit** `fix: P0 follow-ups (sort fan-in rounds, parsing edges, O_TRUNC, BCF flags, CI hygiene)`

### Task 27: P1 benchmarks

**Files:**
- Modify: `xtask/bench/cases.toml`, `xtask/src/bench.rs`, `xtask/src/gen.rs`, `.github/workflows/nightly.yml`
- Create: `xtask/bench/datasets.toml`, `xtask/docker/Dockerfile.bench` (yardsticks), `crates/gloci-engine/benches/index.rs`
- Test: `xtask/tests/gen.rs` (the VCF generator), and the bench `cases_file_parses_and_targets_name_commands` test

**Interfaces:**
- `gen`:
  - `vcf --records N --samples S --seed X` writes a sorted VCF with S samples, a realistic mix of SNVs and indels, multi-allelic sites, missing genotypes, DP/GQ/FT, and QUAL/FILTER.
  - `pair --records N --depth D` writes A/B BED pairs with controlled overlap depth.
- `bench`:
  - `datasets.toml` entries `{ name, url, sha256, unpack }`, fetched to `$GLOCI_BENCH_CACHE` and verified.
  - A case may name `dataset = "…"`.
  - Commands whose tool is missing from the image are skipped with a note (yardsticks only). Binding baselines must exist.

Cases and targets (spec §9.1; all gloci commands use `--chrom-order lex`):
- **Plain BED, `-t 1`, 10M A against 10M B:**
  - `intersect --mode any`, `merge`, `subtract`, `closest -d`, `map -c 5:sum` and `coverage`;
  - each with targets `gloci <= 0.5 * bedtools`, and `gloci <= bedops` where bedops has the operation (intersect `-e`, merge `-m`, subtract `-d` via `merge | subtract`, map via `bedmap`, closest via `closest-features`);
  - yardsticks gia, bedder, granges, bedtk and grit-genomics are recorded with no target.
- **BGZF with indexes, `-t 8`:** the same six, with target `gloci <= 0.2 * bedtools`.
- **`vcf filter`** on dataset `1kg-30x-chr20` (URL and sha256 pinned by the implementer), with `--min-maf 0.01 --max-missing 0.1 --min-qual 30` against vcftools with the equivalent flags: target `gloci <= 0.2 * vcftools`. Yardsticks: `bcftools view`, vcfkit, vcfexpress. A synthetic 2,500-sample VCF case covers local runs.
- **Public workload:** the coitrees benchmark set (A ≈ 2.76M against B ≈ 62M, coverage counts) as `coverage` and `intersect --mode count`, recorded with no target.
- **Index micro-benchmark** (criterion, `cargo bench -p gloci-engine`): coitrees, superintervals, `ITree64` and rust-lapper on a synthetic set and on the coitrees set. If another index is at least 1.2× faster than coitrees on both, switch the D40 default and ledger a ruling with the numbers.

The nightly workflow builds `Dockerfile.bench`, caches the datasets with `actions/cache` keyed on the hash of `datasets.toml`, and runs `cargo xtask bench --check-targets`.

- [ ] **Step 1: Write the failing tests:** the generator output shape (sorted, S sample columns, deterministic for a seed), and the cases file parses with every target naming a defined command.
- [ ] **Step 2: Run** `cargo test -p xtask`. It fails.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test -p xtask`, then locally `cargo xtask bench --only intersect-plain-1M --records 1000000` (best effort on aarch64, R7). Expected: it runs and prints medians.
- [ ] **Step 5: Commit** `feat(xtask): P1 benchmark workloads, datasets, yardsticks and index micro-benchmarks`

### Task 28: Cross-cutting checks: large files, determinism, failure behaviour

**Files:**
- Modify: `crates/gloci/tests/large.rs` (nightly-ignored, as in P0) and `crates/gloci/tests/behaviour.rs`
- Test: those two files

Tests:
- **`large.rs`:**
  - RSS ceilings (≤ 200 MB, the P0 assertion helper) for streaming `intersect`, `merge`, `subtract`, `closest`, `map`, `coverage`, `complement` and `vcf filter` on generated inputs of several GB (`GLOCI_LARGE_BYTES`);
  - RSS assertions for `check` and `index` (a P0 follow-up).
- **`behaviour.rs`:**
  - `determinism_matrix`: every P1 command and mode under `-t 1`, `-t 2` and `-t 8`, plain and BGZF, gives identical bytes.
  - `broken_pipe_matrix`: each command piped to `head -1` exits 0 with nothing on stderr.
  - `empty_inputs_matrix`: empty A, empty B, a header-only VCF A, and B with no records on A's chromosomes; every command exits 0 with the expected header-only or A-only output.
  - `vcf_a_matrix`: VCF A through intersect (any, pairs, count), subtract, closest, map and coverage gives the header plus `##glociCommand` in passthrough modes and the `#CHROM…` line in TSV modes.
  - `failure_matrix`: unsorted A, unsorted B, a malformed B line, truncated BGZF, and output to a missing directory. Each exits with the documented code and message and leaves no partial output.

- [ ] **Step 1: Write the tests.** These are regression matrices over commands that already exist, so most rows pass when written. A row that fails is a real bug.
- [ ] **Step 2: Run** `cargo test -p gloci --test behaviour`, and `cargo test -p gloci --test large -- --ignored` with `GLOCI_LARGE_BYTES=200000000` locally. Expected: pass. Fix any real failure in the owning module with a regression test.
- [ ] **Step 3: Commit** `test: P1 determinism, broken-pipe, empty-input, VCF and failure matrices; RSS ceilings`

### Task 29: Migration guide and `SKILL.md`

**Files:**
- Modify: `docs/migration.md`, `xtask/src/docs.rs`, `xtask/src/main.rs`
- Create: `SKILL.md` (generated and committed)
- Test: `xtask/tests/docs.rs` (`skill_md_is_current`)

**Interfaces:**
- `pub fn generate_skill(out: &Path) -> anyhow::Result<()>` in `xtask/src/docs.rs`, run by `cargo xtask docs --skill SKILL.md`.

Content:
- **`docs/migration.md`:** add the spec §11 bedtools, bedops and vcftools rows for every command that exists after P1, each checked against `--help`. Keep "not yet available (P2/P3)" rows for the rest. Extend the "Other behaviour changes" list with links to every new `docs/differences.md` anchor from Tasks 9–23.
- **`SKILL.md`:**
  - YAML front matter: `name: gloci`, plus a one-sentence `description` saying when an agent should use gloci;
  - core semantics in under 40 lines: 0-based half-open coordinates, points, allele footprints, karyotype vs lex order, sortedness errors and their fixes, VCF passthrough vs TSV output;
  - a command reference generated from clap (each command's about text and its Examples block);
  - the migration tables, included from `docs/migration.md`;
  - pitfalls: `--max-missing` is inverted, `--vcf-span ref` reproduces bedtools, use `--chrom-order lex` for old sorted files, and `--mode` has no `-m`.

  The test regenerates the file into a temp dir and fails if it differs from the committed `SKILL.md`, with the hint `run: cargo xtask docs --skill SKILL.md`.

- [ ] **Step 1: Write the failing test** `skill_md_is_current`.
- [ ] **Step 2: Run** `cargo test -p xtask --test docs`. It fails.
- [ ] **Step 3: Implement** the generator, write the migration rows, and generate `SKILL.md`.
- [ ] **Step 4: Run** `cargo test -p xtask` and the P0 help-examples test (it covers every new subcommand). Expected: pass.
- [ ] **Step 5: Commit** `docs: P1 migration tables and generated SKILL.md`

### Task 30: Release workflow

**Files:**
- Create: `.github/workflows/release.yml`, `xtask/docker/Dockerfile.release`
- Modify: `.github/workflows/ci.yml` (add an `actionlint` step)
- Test: `actionlint` in CI; a local `docker build -f xtask/docker/Dockerfile.release` against a locally built binary

Per D59, on a pushed tag `v*`:
- **Linux:** build `x86_64-unknown-linux-musl` and `aarch64-unknown-linux-musl` release binaries (`cargo zigbuild`, or `cross` with pinned versions). Assert that neither has an ELF interpreter (the P0 readelf check).
- **macOS:** build `aarch64-apple-darwin` on the current macOS runner.
- Run `gloci --version` on each native build, package each as `gloci-<tag>-<target>.tar.gz` with `SHA256SUMS`, and create the GitHub release with the assets (a pinned release action).
- Build `Dockerfile.release` (`FROM scratch`, the static x86_64 binary, `ENTRYPOINT ["/gloci"]`) and push `ghcr.io/<owner>/gloci:<tag>` and `:latest` with `docker/build-push-action`, using `permissions: packages: write`.
- Every third-party action is pinned to a commit SHA.
- Nothing is published to crates.io or bioconda.

- [ ] **Step 1: Write** `release.yml` and `Dockerfile.release`, and add `actionlint` to `ci.yml`.
- [ ] **Step 2: Verify** that `actionlint` (or `docker run rhysd/actionlint`) is clean, and that `docker build -f xtask/docker/Dockerfile.release --build-arg BIN=target/release/gloci .` builds and `docker run --rm IMAGE --version` prints the version, using the local glibc build for the smoke test only.
- [ ] **Step 3: Commit** `ci: release workflow (static musl and macOS binaries, checksums, ghcr image)`

---

## Execution notes for the controller

- **Order:** Tasks 1–8 are the foundation. Tasks 9–19 depend on 8 (and 11 depends on 10), 20–23 form the VCF track, and 24–30 come last. Run them in numeric order.
- **Spec amendments:** made in the same commit as this plan, §4.3 and §7.1. If a differential case shows bedtools behaviour that conflicts with a decision here, rule on it (spec first, then the plan), ledger it, and pin the difference.
- **Docker:** every task that adds differential cases runs them in the pinned `gloci-diff:p0` image, extended only by Task 27's separate bench image. A task adding a reference tool to `xtask/docker/Dockerfile` pins its version.
