# gloci P2a Breadth Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the command-line half of spec Phase 2:
- fix the three P1 residual issues;
- add the remaining interval commands: `shift`, `cluster`, `spacing`, `groupby`, `partition`, `genomecov`;
- complete the `map` op set;
- add `--chrom-map`, WIG input, and `--region` on relation commands;
- add the VCF QC trio: `vcf stats`, `vcf extract`, `vcf concat`;
- add bcftools-style filter expressions;
- add native BCF genotype encoding;
- switch the unsorted-index library;
- add per-chromosome parallelism for plain seekable files;
- add P2 benchmarks and docs.

**Architecture:** No new crates.
- Interval commands go in gloci-ops on the existing engine. `partition`, `cluster` and `genomecov` use `kmerge` and a boundary sweep; `groupby`, `spacing` and `shift` are streaming transforms.
- VCF work goes in gloci-vcf: an expression parser and evaluator, and a one-pass statistics engine that shares the `vcf filter` selection pipeline. gloci-ops wraps them.
- `--chrom-map` is applied in gloci-io's readers, so every command gets it.
- The Python library, packaging and the mdBook are **P2b**, a separate plan.

**Tech Stack:** Everything P1 uses, plus rust-lapper (replacing coitrees, per D61) and a hand-written recursive-descent parser for expressions (no parser-generator dependency).

**Spec:** `docs/superpowers/specs/2026-11-06-gloci-design.md`, covering:
- §4.5 (`--chrom-map`);
- §5.1 (WIG);
- §6.3 (full op set);
- §7.2 (`partition`);
- §7.3 (`shift`, `cluster`, `spacing`, `groupby`, `genomecov`);
- §7.8 (`vcf stats`/`extract`/`concat`, expressions);
- §8.3 (seekable-file parallelism);
- §9 (benchmarks, testing).

**Builds on:** P1 (merged as 8b18147; follow-ups listed in `docs/superpowers/plans/2026-11-08-gloci-p1-followups.md`).
- P0's D1–D33 and P1's D34–D60 stay in force, as amended by P1's final rulings.
- In particular, `--sort` and `--unsorted` (D39) carry over unchanged.
- So do the TSV rules (D42), the placeholder rule (D45) and the aggregation rules (D46).

**Out of scope (P2b or later):**
- the Python library (spec §8.5) and its Arrow conversion;
- bioconda and PyPI wheels;
- the mdBook;
- the bioRxiv external benchmark datasets;
- the spec §8.3 stage-2 parse/format worker pools (D77);
- BAM/CRAM and the other P3 items.

**P2a is done when:**
- every listed command and option exists, with differential cases against pinned bedtools, bedops, vcftools and bcftools, and with every difference pinned;
- every vcftools row of `docs/migration.md` that is marked P2 maps to a working gloci command;
- the nightly x86_64 benchmark shows `vcf stats` with 5 statistics at least 10× faster than 5 separate vcftools runs.

## Global Constraints

These are P1's Global Constraints. Each still applies:
- 0-based half-open u64 coordinates.
- Rust stable, edition 2024, `rust-version = "1.91"`.
- No htslib; one static musl binary.
- Exit codes 0, 1 (data or runtime), 2 (usage); broken pipe exits 0.
- Errors say what, where and how to fix.
- Output is byte-identical for every `--threads`.
- BGZF `.gz` output, compression level 6 by default.
- Atomic output with registered cleanup (P1 Task 26a API).
- Content detection; `-` reads stdin.
- VCF output carries a `##glociCommand` line.
- `--vcf-span allele` is the default.
- The §4.3 overlap and §4.4 distance rules.
- Sortedness is checked on the fly.
- Record streams in gloci-ops (§8.1). No logic lives only in the CLI.
- The §4.6 output-format rules.
- One flag means one thing everywhere.
- Streaming memory is flat, ≤ about 200 MB excluding active-set depth.
- Property tests against naive references; thread-determinism tests.
- Differential cases with goldens, every difference anchored and pinned.
- 2–4 `--help` examples per command.
- CI on Linux x86_64/aarch64 and macOS aarch64.

New in P2a:
- `vcf stats` computes any number of statistics in **one pass**. Its output columns use vcftools' names wherever an equivalent exists (spec §7.8).
- Filter expressions follow bcftools' semantics for every construct gloci supports. Anything unsupported is a usage error naming the construct, never a silent misreading.

## Decisions this plan adds

| # | Decision |
|---|---|
| D61 | The unsorted index switches from coitrees to **rust-lapper**, by the D40 rule. P1's real-crate benchmark measured rust-lapper 2.02× and 1.80× faster on both synthetic sets. rust-lapper is generic over `u64`, so `ITree64` and the i32 fallback are removed. Queries are still widened by one base on each side and every candidate is re-tested exactly. `unsorted_equals_sorted` must keep passing unchanged. |
| D62 | `--chrom-map FILE`, a global option on every command that reads intervals or variants. <ul><li>The file has two tab-separated columns, `old new`, with `#` comment lines allowed.</li><li>A duplicate `old` value is a usage error. An unmapped name passes through unchanged.</li><li>Renaming happens in the readers before interning, so order checks, `-g` lookups and indexes all see the new name. That includes VCF/BCF header contig lines.</li><li>Output records carry the new name: BED/GFF/VCF text records get their chromosome column rewritten, and BCF output gets a renamed header dictionary.</li><li>Index queries translate a requested region name back to the file's own name.</li><li>`-g` names are never renamed.</li></ul> |
| D63 | WIG (`fixedStep`/`variableStep`, with optional `span`) is detected from its first step line and read as bedGraph-like records (chrom, start, end, value). `track` and `browser` lines are header lines. WIG works wherever an interval input is accepted, with value column 4, and `gloci convert` turns WIG into bedGraph (one row per WIG data line, no merging of equal values). |
| D64 | The remaining spec §6.3 ops: <ul><li>`cv` = sample stdev ÷ mean (`--empty` when the mean is 0).</li><li>`mad` = median of absolute deviations from the median, unscaled.</li><li>`quantile:Q`, Q ∈ [0,1], uses linear interpolation (type 7).</li><li>`tmean:LO-HI`, LO < HI ∈ [0,1], is the mean of the values whose rank fraction lies in [LO, HI] (as `quantile`).</li><li>`wmean` is the mean weighted by overlapping bases; it is valid only in two-input commands.</li><li>`mode` and `antimode` are the most and least frequent value (compared as text, missing values skipped); ties go to the byte-wise smallest value, as in bedtools [revised in Task 6 review].</li><li>`bases_uniq` is the union of overlapping bases.</li><li>`frac_a` = `bases_uniq` ÷ A's length (`.` for a point A). `frac_b` = summed overlap ÷ summed B length.</li><li>`min_record`, `max_record` and `records` write whole B records into one cell. Inside a record, tabs become `\|`; records are separated by `--delim`. They are valid only in two-input commands.</li></ul> |
| D65 | `cluster` appends a cluster id. Ids are global and increase through the file, as bedtools does. `-d/--max-gap` and `-s` mean the same as in `merge`. |
| D66 | `spacing` appends the distance to the previous record on the same chromosome, following the bedtools convention; differential cases decide the exact edges. That is `.` for the first record and 0 for bookended records; for overlapping records the bedtools value is pinned, or documented as a difference. |
| D67 | `groupby -g COLS -c COL:OP…` groups consecutive lines whose key columns are equal, without requiring sort order (bedtools semantics). <ul><li>Output is the group's key columns followed by the aggregations.</li><li>It works on any tab-separated text, including BED. VCF input is a usage error that suggests `vcf extract`.</li><li>`--header` takes key names from the header line.</li></ul> |
| D68 | `partition` follows spec §7.2. <ul><li>Membership is written as two columns: a count, then the labels joined with commas. Labels come from `--names`, defaulting to basenames.</li><li>With `--bedgraph`, each input contributes one value column instead of membership. Gaps are filled from `--empty`, and `-g` together with `--empty` covers the whole genome.</li><li>Points contribute no segments, since they cover no bases.</li><li>Inputs are k-way merged and order-checked; `--sort` is available.</li></ul> |
| D69 | `genomecov` (BED in P2a; BAM is P3) follows bedtools `genomecov`. <ul><li>`-g` is required.</li><li>Modes: `--hist` (the default) gives per-chromosome rows plus `genome` rows. `--bedgraph` and `--bga` (the latter including zero-depth runs) give runs. `--per-base` is 1-based. `--scale X`, `--strand +\|-`, `--split`, `--5p`, `--3p` and `--max-depth N` (capped into the last bin) are also supported.</li><li>Output columns follow bedtools.</li></ul> |
| D70 | `shift --by N` takes a signed N. <ul><li>`--pct` makes N a fraction of the feature's length.</li><li>`-s` negates the shift on `-` strand records, treating unknown strand as `+`.</li><li>Output is clipped as `slop` clips it: at 0, and at the chromosome end with `-g`.</li><li>The transform order warning (P1 rule) applies.</li></ul> |
| D71 | Filter expressions for `vcf filter --include EXPR` / `--exclude EXPR` (and `vcf stats`/`vcf extract`) use the grammar in Task 13: a supported subset of bcftools `-i`/`-e`. <ul><li>Semantics follow bcftools for every supported construct. A `FMT/` comparison holds when **any** sample satisfies it.</li><li>A comparison against a missing value is false. Only `!=` with a missing value is true, matching bcftools.</li><li>Unsupported functions or fields are usage errors naming the token, with a hint.</li><li>Expressions are evaluated in pipeline step 4, ANDed with the other site filters and before `--thin`.</li></ul> |
| D72 | `vcf stats IN --out-prefix P [stat flags]` writes `P.<ext>` for each statistic requested, with vcftools' file extensions and column names (Task 15 table). `--summary` writes `P.summary.txt`, and `--json` writes `P.summary.json`. Every `vcf filter` selection option applies first, through the same pipeline code, in one pass. |
| D73 | `vcf extract`: <ul><li>`--info KEYS` writes TSV with columns `CHROM POS ID REF ALT KEY…`; a missing value is `.` and a flag is `1`/`0`.</li><li>`--format KEYS --layout wide` writes one row per site with `<sample>:<KEY>` columns. `--layout long` writes one row per site and sample, as `CHROM POS ID REF ALT SAMPLE KEY…`.</li><li>`--genotype-matrix 012` writes vcftools' `P.012`, `P.012.indv` and `P.012.pos`. A cell holds the count of non-reference alleles, and `-1` when missing.</li></ul> |
| D74 | `vcf concat FILES…` writes the inputs one after another. <ul><li>Samples must be identical, in the same order.</li><li>Headers must be compatible, by P1's multi-VCF rule: identical INFO, FILTER, FORMAT, ALT and contig definitions, with PASS implicit.</li><li>Within the chromosome table, the last record of file *i* must sort at or before the first record of file *i+1*. Otherwise it is an error with a `gloci sort` hint.</li><li>`--naive` needs BGZF inputs with byte-identical headers. It copies the BGZF blocks after each file's header blocks without recompressing them, and drops every EOF marker except the last.</li></ul> |
| D75 | `--region`/`--regions-file` on the relation commands restrict **A** to records overlapping the regions, with view semantics (`--inclusive-points` honoured). An indexed A is read through its index. B is read normally, or through its index for the same regions widened by the window or `--max-distance`. |
| D76 | Per-chromosome parallelism for plain seekable files. An uncompressed BED/TSV/VCF that is a regular file gets one memchr scan that records each chromosome's byte range and checks that chromosomes are contiguous. Workers then read those byte ranges through the P1 per-chromosome driver. It applies when every input is either indexed or a plain seekable file, and output stays byte-identical to streaming. At `-t 1` gloci streams. |
| D77 | The spec §8.3 stage-2 parse/format worker pools are deferred to P3, unless the P2a benchmarks miss a target. |
| D78 | Native BCF genotype encoding: gloci writes the BCF FORMAT `GT` field itself, with the int8/int16 widths and end-of-vector padding for mixed ploidy. That lifts P1's `#bcf-mixed-ploidy-missing` limitation, so that docs/differences.md section and its pin are removed. |

## Review Focus

1. **`--chrom-map` round trips.** A renamed chromosome must:
   - pass order checks against the renamed `-g`;
   - query correctly through an index built on the original names;
   - appear renamed in BED, VCF and BCF output.

   Task 4 `chrom_map_everywhere`.
2. **Expression semantics against bcftools** on missing values, multi-valued INFO indexes, and `FMT/` any-sample. Task 13 `expr_matches_bcftools` (differential).
3. **The `vcf stats` one-pass promise.** Requesting 5 statistics must equal 5 single-statistic runs byte for byte, under sample selection and masking. Task 16 `one_pass_equals_separate_runs`.
4. **`partition`/`genomecov` with points, nested features and empty chromosomes.** Equal to the naive references. Tasks 10–11 property tests.
5. **`vcf concat --naive` at BGZF block edges.** Headers spanning several blocks, files without an EOF marker, and the output must be readable by bcftools. Task 18 `naive_concat_block_edges`.

## File structure

```
crates/gloci-engine/src/unsorted.rs (rust-lapper), itree64.rs (deleted)
crates/gloci-io/src/chrom_map.rs (new), wig.rs (new), bcf_encode.rs (new; native GT)
crates/gloci-ops/src/{shift, cluster, spacing, groupby, partition, genomecov}.rs (new);
                    agg/ops.rs (full op set); relate (regions; seekable chunks)
crates/gloci-vcf/src/expr/{lexer, parser, eval}.rs (new); stats/{engine, site, sample, report}.rs (new);
                    extract.rs, concat.rs (new)
crates/gloci-ops/src/{vcf_stats, vcf_extract, vcf_concat}.rs (new)
crates/gloci/src/commands/{shift, cluster, spacing, groupby, partition, genomecov}.rs; vcf.rs (stats/extract/concat)
xtask/bench/cases.toml (P2 cases); docs/migration.md, docs/differences.md, SKILL.md
```

---

### Task 1: P1 residual fixes

**Files:** crates/gloci-ops/src/relate/ (the indexed driver), crates/gloci/src/commands/sort.rs, crates/gloci-io/src/input.rs, xtask/bench/cases.toml, crates/gloci-engine/src/parallel.rs, crates/gloci/src/commands/vcf.rs, and the five relation help texts in cli.rs.

Items. The detail for each is in `docs/superpowers/plans/2026-11-08-gloci-p1-followups.md`, items 1–3 and the parked minors.
1. **Output changes with `--threads` on indexed inputs out of order.**
   - `run_indexed` checks each input's index reference order against the chromosome table. If the order differs, it raises the same order error streaming gives, with the lex-note hint where it applies.
   - Remove the false sentence "the output is the same as with one thread". Say instead that indexed and streaming give identical output for inputs in table order.
   - Test `indexed_and_streaming_refuse_alike`: a lex-sorted indexed BED and a VCF whose header contig order differs from its record order. Both `-t 1` and `-t 4` exit 1 with the same message.
2. **`sort` holds every regular-file reader open.**
   - Keep readers open, or spool them, only for non-regular inputs. Reopen regular files on each pass.
   - Treat `/dev/fd/*` as non-reopenable.
   - Test `sort_many_bgzf_inputs_bounded`: 100 BGZF inputs at `-t 8`, with fewer than 120 peak threads read from `/proc/self/status` on Linux. Peak RSS stays under 150 MB, measured as a `#[ignore]` large test.
3. **The degenerate 1kg benchmark case.** Remove `--min-qual` (and the vcftools/bcftools equivalents) from `vcf-filter-1kg-chr20`, so the case times the MAF and missingness pipeline. Add a bench test that fails if any binding `vcf filter` case's gloci command keeps zero records on its dataset's first 10,000 lines (a sampled dry run).
4. **Parked minors.**
   - `run_chunks` returns `IoError::Internal` when it gets fewer results than jobs.
   - `vcf filter` runs `validate()` before reading any list file.
   - Indexed-mode skip and parse messages name the line (`loc` is now numbered), not a BGZF offset, where the line is known.

- [ ] Step 1: write the failing tests above.
- [ ] Step 2: run `cargo test -p gloci --test sort --test parallel_relations` and confirm they fail.
- [ ] Step 3: implement.
- [ ] Step 4: run `cargo test --workspace`. Expected: pass.
- [ ] Step 5: commit `fix: P1 residuals (indexed order check, sort reopens regular files, bench case, minors)`

### Task 2: Unsorted index on rust-lapper (D61)

**Files:** crates/gloci-engine/{Cargo.toml, src/unsorted.rs}. Delete src/itree64.rs and its tests. Update benches/index.rs to keep the comparison.

**Interfaces:** `sweep_unsorted` keeps the same signature and `Visitor` contract. Stats fields that reported index counts are renamed to match the new index (or removed). Grep the tests that read them.

- [ ] Step 1: keep `unsorted_equals_sorted` and every unsorted test. Add `lapper_handles_u64_coordinates` (B at 3,000,000,000 and at `u64::MAX - 10`) and `lapper_points_and_bookends` (widen-and-retest).
- [ ] Step 2: implement the switch. Use a `Lapper<u64, u32>` per chromosome and file, holding indexes into the B table. Use lapper's sorted `seek` when A arrives sorted, and plain `find` otherwise.
- [ ] Step 3: run `cargo test -p gloci-engine` and `cargo bench -p gloci-engine --bench index -- --quick`. Paste the numbers.
- [ ] Step 4: commit `perf(engine): unsorted mode on rust-lapper (u64, no i32 fallback)`

### Task 3: Native BCF genotype encoding (D78)

**Files:** crates/gloci-io/src/bcf_encode.rs (new) and crates/gloci-io/src/convert.rs, which uses it for `GT` when re-encoding text records to BCF.

**Interfaces:** `pub fn encode_gt(samples: &[&[u8]], out: &mut Vec<u8>) -> Result<(), String>`. It takes text GT values (`0/1`, `1|0`, `.`, `./.`, haploid, triploid) and writes the BCF typed FORMAT GT block. The width is int8 when every allele index is ≤ 63, otherwise int16. Padding uses the end-of-vector sentinel, and missing alleles are 0.

- [ ] Step 1: tests:
  - `encode_gt_roundtrip`: text, then `encode_gt`, then gloci-vcf's `bcf_genotypes` decoder gives the same `SampleFields`. Cover mixed ploidy with a bare `.`, alleles ≥ 63, phased and unphased, and a leading phase prefix.
  - `vcf_to_bcf_mixed_ploidy`: the input that used to fail now converts, and `bcftools view` reads it identically. Do this as a differential case that replaces the `#bcf-mixed-ploidy-missing` pin.
- [ ] Step 2: implement, routing the converter's GT field through it.
- [ ] Step 3: remove the docs/differences.md section and update the pinned case to a matching case.
- [ ] Step 4: run `cargo test --workspace` and `cargo xtask diff --only vcf-to-bcf-mixed-ploidy-vs-bcftools`.
- [ ] Step 5: commit `feat(io): native BCF GT encoding (mixed ploidy with missing haploid)`

### Task 4: `--chrom-map` (D62)

**Files:**
- crates/gloci-io/src/chrom_map.rs (new);
- record.rs and every reader (BED, GFF, VCF, BCF, and WIG later), which apply the map before interning;
- vcf.rs (header contig lines) and bcf.rs (dictionary);
- writer.rs (rewriting output records);
- index/query.rs (translating region names back);
- crates/gloci/src/cli.rs (a global `ChromMapArg`, flattened into every command that reads records).

**Interfaces:**
- `pub struct ChromMap`, with `pub fn load(path: &str) -> Result<ChromMap, IoError>`, `pub fn rename<'a>(&self, name: &'a [u8]) -> &'a [u8] or Cow` and `pub fn reverse(&self, new: &[u8]) -> Option<&[u8]>`.
- `ReaderOptions.chrom_map: Option<Arc<ChromMap>>`.

- [ ] Step 1: tests:
  - `chrom_map_everywhere`, run for sort, intersect, merge, view region and vcf filter. A `1→chr1` map on `1`-named BED, VCF and BCF inputs against `-g hg38` gives renamed output and passes order checks. A region query `chr1:…` on an index built with `1` returns the right records.
  - `chrom_map_errors`: a duplicate old name exits 2; a malformed line exits 2 naming the line.
  - The P1 naming warning now suggests `--chrom-map`.
- [ ] Step 2: implement.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: --chrom-map renames chromosomes on input`

### Task 5: WIG input (D63)

**Files:** crates/gloci-io/src/wig.rs (new), detect.rs (detection), record.rs (dispatch), crates/gloci-ops/src/convert (WIG → bedGraph).

- [ ] Step 1: tests:
  - `wig_fixed_and_variable`: fixedStep with step and span, and variableStep with span. The records and their 0-based conversion are exact.
  - `wig_as_b_in_map`: `map -b x.wig -c 4:mean`.
  - `convert_wig_to_bedgraph`, as a differential case against UCSC `wigToBigWig` if present in the image. Otherwise use a hand-checked golden.
  - `wig_errors`: a data line before any step line exits 1 naming the line.
- [ ] Step 2: implement.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat(io): WIG input (fixedStep/variableStep) and WIG→bedGraph`

### Task 6: The full `map` op set (D64)

**Files:** crates/gloci-ops/src/agg/{spec.rs, ops.rs}, plus tests.

- [ ] Step 1: value tests with exact strings for each new op.
  - Fixed values for `1 2 3 4 10`: `mad` is `1` and `quantile:0.25` is `2`.
  - For `cv` and `tmean:0.2-0.8`, compute the expected string in the test with an independent helper written from the D64 formula. Don't hard-code it.
  - Also test:
    - `mode` and `antimode` ties;
    - `wmean` with overlap weights;
    - `bases_uniq` against summed bases on overlapping B;
    - `frac_a` on a point A gives `.`;
    - `records` with tabs gives `|` inside the record.
  - Usage errors:
    - `wmean`, `records`, `min_record`, `max_record`, `frac_*` and `bases_uniq` in `merge` and `groupby`;
    - Q out of range;
    - LO ≥ HI.
- [ ] Step 2: differential cases against `bedtools map -o OP` for every new op that bedtools 2.31.1 also has. List them from `bedtools map -h` in the image: at least `mode`, `antimode` and `median`, and `stdev` against `stdev_pop`. Compare numerically and pin each difference.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat(ops): full -c COL:OP set`

### Task 7: `gloci shift` (D70)

**Files:** crates/gloci-ops/src/shift.rs (reuse transform.rs machinery), crates/gloci/src/commands/shift.rs, cli.rs.

- [ ] Step 1: tests:
  - `shift_table`: `+10`, `-10` clipped at 0, `-g` clipping at the end, `-s` on a minus-strand record, and `--pct 0.5`;
  - the order warning (only when the shift itself breaks order);
  - GFF rewrite;
  - VCF is a usage error;
  - determinism.
- [ ] Step 2: differential cases against `bedtools shift -s 10 -g g`, `-m`/`-p` strand forms, and `-pct`.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: gloci shift`

### Task 8: `gloci cluster` and `gloci spacing` (D65, D66)

**Files:** crates/gloci-ops/src/{cluster.rs, spacing.rs}, commands, cli.

- [ ] Step 1: tests:
  - cluster: gap table (`-d 0`, `-d 10`, negative `-d`), `-s`, global ids across chromosomes, points, a VCF footprint input, and a property test against a naive version;
  - spacing: first record `.`, bookended, overlapping, chromosome change;
  - determinism for both.
- [ ] Step 2: differential cases against `bedtools cluster` (default, `-d 10`, `-s`) and `bedtools spacing`. Pin differences.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: gloci cluster and gloci spacing`

### Task 9: `gloci groupby` (D67)

**Files:** crates/gloci-ops/src/groupby.rs, commands, cli.

- [ ] Step 1: tests:
  - consecutive grouping, where unsorted input produces repeated groups, as in bedtools;
  - multiple key columns;
  - `--header` names;
  - every scalar op;
  - the usage errors from D64 for two-input-only ops;
  - VCF input is a usage error;
  - determinism.
- [ ] Step 2: differential cases against `bedtools groupby -g 1,2 -c 5 -o sum,mean` and similar, compared numerically.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: gloci groupby`

### Task 10: `gloci partition` (D68)

**Files:** crates/gloci-ops/src/partition.rs (`kmerge` plus a boundary sweep), commands, cli, `naive_partition` in the testkit.

- [ ] Step 1: tests:
  - membership columns for 3 inputs;
  - `--min-files`, `--max-files` and `--all`;
  - `--merge-adjacent`;
  - `--no-membership`;
  - `--bedgraph` with `--empty 0`, and `-g` with whole-genome fill;
  - points contribute nothing;
  - `--sort`;
  - a property test against `naive_partition`, with nested features, points and empty chromosomes.
- [ ] Step 2: differential cases, each comparing against `partition` with the flags given:
  - `bedops --partition` against `--no-membership`;
  - `bedops --intersect` against `--all --merge-adjacent --no-membership`;
  - `bedops --symmdiff` against `--max-files 1 --merge-adjacent --no-membership`;
  - `bedtools multiinter` against default membership, after column reshaping by a shell pipeline (live) or a committed golden twin;
  - `bedtools unionbedg` against `--bedgraph`.

  Pin differences.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: gloci partition (multiinter/unionbedg/bedops -p -i -s)`

### Task 11: `gloci genomecov` for BED (D69)

**Files:** crates/gloci-ops/src/genomecov.rs, commands, cli, `naive_genomecov` in the testkit.

- [ ] Step 1: tests:
  - every mode on a fixture with exact output: `--hist` with `genome` rows, `--bedgraph`, `--bga`, `--per-base`;
  - `--scale`, `--strand`, `--split`, `--5p`/`--3p`, `--max-depth`;
  - an empty chromosome;
  - a property test of depth runs against `naive_genomecov`;
  - determinism.
- [ ] Step 2: differential cases against `bedtools genomecov -i x.bed -g g` (default), `-bg`, `-bga`, `-d`, `-strand +`, `-5`, `-scale 2`, `-split` and `-max 3`. Pin differences.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: gloci genomecov (BED)`

### Task 12: `--region` on relation commands (D75)

**Files:** crates/gloci-ops/src/relate/ (stream and indexed entry points), the five relation commands in cli.rs (`RegionArgs`).

- [ ] Step 1: tests:
  - `region_restricts_a`: intersect, closest, map and coverage with `--region chr2:100-500`, against the full run filtered to A records overlapping the region;
  - indexed A uses its index (a debug log line);
  - `closest` still finds B outside the region;
  - `--inclusive-points` at region edges;
  - determinism.
- [ ] Step 2: implement.
- [ ] Step 3: run `cargo test --workspace`.
- [ ] Step 4: commit `feat: --region/--regions-file on relation commands`
### Task 13: Filter expression language (D71)

**Files:**
- crates/gloci-vcf/src/expr/{mod.rs, lexer.rs, parser.rs, eval.rs} (new)
- Test: crates/gloci-vcf/tests/expr.rs

**Interfaces:**
- `pub struct Expr` (the parsed AST).
- `pub fn parse(src: &str, header: &VcfHeader) -> Result<Expr, UsageError>`, which resolves field names and types against the header.
- `impl Expr { pub fn eval(&self, site: &SiteView, gt: Option<&dyn GenotypeSource>, stats: Option<&SiteStats>, scratch: &mut EvalScratch) -> Result<bool, String>; pub fn needs_genotypes(&self) -> bool }`

Grammar (a supported subset of bcftools `-i`/`-e`):
```
expr   := or
or     := and ( "||" and )*
and    := unary ( "&&" unary )*
unary  := "!" unary | cmp
cmp    := term ( OP term )?          OP := "==" | "=" | "!=" | "<" | "<=" | ">" | ">=" | "~" | "!~"
term   := NUMBER | STRING | "." | FIELD | FUNC "(" expr ")" | "(" expr ")"
FIELD  := QUAL | POS | CHROM | ID | REF | ALT | FILTER | TYPE | N_ALT | N_SAMPLES
        | INFO/KEY [ "[" INT "]" ] | KEY (bare, means INFO/KEY) | FMT/KEY [ "[" INT "]" ] | FORMAT/KEY | GT
        | MAF | AC | AN | F_MISSING | N_MISSING
FUNC   := N_PASS
```
Semantics follow bcftools:
- **Multi-valued fields.** A multi-valued INFO without an index matches if any value matches.
- **`FMT/` and `GT`.** A comparison holds if any sample holds. `N_PASS(e)` counts the samples where `e` holds.
- **`GT` keywords.** `GT` compares with `"het"`, `"hom"`, `"ref"`, `"alt"`, `"mis"` and `"hap"`.
- **`TYPE` keywords.** `TYPE` takes `"snp"`, `"mnp"`, `"indel"`, `"other"` and `"ref"`, mapped from gloci's `VariantKind`.
- **`~`.** It is a regular-expression match. Use the `regex` crate (pure Rust).
- **Missing values.** A comparison involving a missing value is false, except `!=`.
- **Computed fields.** `MAF`, `AC`, `AN`, `F_MISSING` and `N_MISSING` are computed from genotypes after masking. `INFO/AC` reads the INFO field instead.
- **Unsupported input.** The sample-wise operators `&` and `|`, and any function not listed, are usage errors that name the token and say "not supported; use && / ||" (or list the supported functions).

- [ ] **Step 1: Write the failing tests.**
  - The parser: precedence, parentheses, errors with the column of the bad token, and unknown INFO keys (a usage error, as in bcftools).
  - Evaluation tables:
    - `QUAL>30` with QUAL `.` is false;
    - `INFO/DP>10` and `DP>10`;
    - `INFO/AF[1]>0.1`;
    - `FILTER="PASS"`;
    - `FMT/DP>20` matches when any sample does;
    - `N_PASS(FMT/GQ>30)>=2`;
    - `GT="het"`;
    - `TYPE="indel"`;
    - `MAF<0.05`;
    - `ID~"^rs"`;
    - `!(QUAL<20 || DP<5)`.
  - BCF and text records give the same result.
- [ ] **Step 2: Run** `cargo test -p gloci-vcf --test expr`. It fails.
- [ ] **Step 3: Implement.** Evaluation must not allocate per record; regexes are compiled once.
- [ ] **Step 4: Run** `cargo test -p gloci-vcf`. Expected: pass.
- [ ] **Step 5: Commit** `feat(vcf): bcftools-style filter expressions`

### Task 14: `vcf filter --include/--exclude`

**Files:**
- crates/gloci-ops/src/vcf_filter.rs (expression in step 4, before `--thin`)
- crates/gloci/src/cli.rs (`--include EXPR` and `--exclude EXPR`; the two are mutually exclusive)

- [ ] **Step 1: Write the tests.**
  - The pipeline order: an expression on `MAF` sees masked and selected samples.
  - `--exclude` negates.
  - A bad expression exits 2 before reading input.
  - Determinism.
- [ ] **Step 2: Add differential cases.** Name them `expr_matches_bcftools`. For each expression in a set of 15 (covering every grammar construct), compare `gloci vcf filter in.vcf --include E` with `bcftools view -i E in.vcf`. Use `vcf-records` compare mode. Pin every difference, for example MAF on multi-allelic sites, referencing the existing anchor.
- [ ] **Step 3: Run** `cargo test --workspace` and the diff cases.
- [ ] **Step 4: Commit** `feat: vcf filter --include/--exclude expressions`

### Task 15: `vcf stats` engine and site statistics (D72)

**Files:**
- crates/gloci-vcf/src/stats/{mod.rs, engine.rs, site.rs} (new)
- Test: crates/gloci-vcf/tests/stats_site.rs

**Interfaces:**
- `pub trait Stat { fn site(&mut self, s: &SiteView, gt: Option<&dyn GenotypeSource>, st: Option<&SiteStats>) -> Result<(), String>; fn finish(&mut self, out: &mut dyn StatSink) -> io::Result<()>; fn needs_genotypes(&self) -> bool; fn ext(&self) -> &'static str }`
- `pub struct StatsEngine { stats: Vec<Box<dyn Stat>> }`, which runs every requested stat in one pass.

Site statistics. Extensions and columns follow vcftools 0.1.16, and the differential cases decide the remaining formatting details such as number precision:

| Flag | File | Columns |
|---|---|---|
| `--freq` | `.frq` | `CHROM POS N_ALLELES N_CHR {ALLELE:FREQ}` |
| `--counts` | `.frq.count` | `CHROM POS N_ALLELES N_CHR {ALLELE:COUNT}` |
| `--site-depth` | `.ldepth` | `CHROM POS SUM_DEPTH SUMSQ_DEPTH` |
| `--site-mean-depth` | `.ldepth.mean` | `CHROM POS MEAN_DEPTH VAR_DEPTH` |
| `--geno-depth` | `.gdepth` | `CHROM POS <sample>…` |
| `--site-quality` | `.lqual` | `CHROM POS QUAL` |
| `--missing-site` | `.lmiss` | `CHR POS N_DATA N_GENOTYPE_FILTERED N_MISS F_MISS` |
| `--hardy` | `.hwe` | `CHR POS OBS(HOM1/HET/HOM2) E(HOM1/HET/HOM2) ChiSq_HWE P_HWE P_HET_DEFICIT P_HET_EXCESS` |
| `--singletons` | `.singletons` | `CHROM POS SINGLETON/DOUBLETON ALLELE INDV` |
| `--filter-summary` | `.FILTER.summary` | `FILTER N_VARIANTS N_Ts N_Tv Ts/Tv` |
| `--snp-density BIN` | `.snpden` | `CHROM BIN_START SNP_COUNT VARIANTS/KB` |
| `--indel-hist` | `.indel.hist` | `LENGTH COUNT PCT` |
| `--tstv` | `.TsTv.summary` | `MODEL COUNT` |

`--hardy` needs the HWE heterozygote-deficit and heterozygote-excess p-values. Add them to gloci-vcf's `hwe` module (Wigginton's one-sided tails), with tests against a brute-force enumeration.

- [ ] **Step 1: Write the failing tests.** Exact output for each statistic on a small fixture, hand-checked. HWE one-sided tails against the brute force.
- [ ] **Step 2: Run** `cargo test -p gloci-vcf --test stats_site`. It fails.
- [ ] **Step 3: Implement.** Per-site work must not allocate. Statistics that need whole-file state (`--filter-summary`, `--snp-density`, `--indel-hist`, `--tstv`) accumulate in bounded maps.
- [ ] **Step 4: Run** `cargo test -p gloci-vcf`.
- [ ] **Step 5: Commit** `feat(vcf): vcf stats engine and site statistics`

### Task 16: `vcf stats` per-sample statistics, the summary, and the command

**Files:**
- crates/gloci-vcf/src/stats/{sample.rs, report.rs}
- crates/gloci-ops/src/vcf_stats.rs
- crates/gloci/src/commands/vcf.rs (`stats` subcommand)
- crates/gloci/src/cli.rs
- Test: crates/gloci/tests/vcf_stats.rs, plus cases

Per-sample statistics:

| Flag | File | Columns |
|---|---|---|
| `--depth` | `.idepth` | `INDV N_SITES MEAN_DEPTH` |
| `--missing-indv` | `.imiss` | `INDV N_DATA N_GENOTYPES_FILTERED N_MISS F_MISS` |
| `--het` | `.het` | `INDV O(HOM) E(HOM) N_SITES F` |

`--summary` writes `P.summary.txt`, a short human overview: counts by type, Ts/Tv, samples, mean depth and missingness. `--json` writes the same summary as `P.summary.json`.

The command:
- takes every `vcf filter` selection flag through a shared clap group, including `--include` and `--exclude`;
- requires `--out-prefix P`;
- writes its outputs atomically, so a failure publishes no file.

- [ ] **Step 1: Write the tests.**
  - `one_pass_equals_separate_runs`: 5 statistics in one run equal 5 single-statistic runs byte for byte, with `--samples` and `--min-dp` masking.
  - Per-sample output is exact on a fixture.
  - The summary text and JSON follow a schema.
  - A usage error when no statistic is requested.
  - Determinism.
- [ ] **Step 2: Add differential cases** against vcftools for each statistic: `vcftools --vcf in.vcf --freq --out ref`, then compare `ref.frq` with gloci's `P.frq`, and the same for every flag in both tables. Compare numerically where the output is floating point. Pin differences, such as MAF and allele handling on multi-allelic sites and F_MISS definitions, with anchors.
- [ ] **Step 3: Run** `cargo test --workspace` and the diff cases.
- [ ] **Step 4: Commit** `feat: gloci vcf stats (one pass, vcftools-compatible outputs)`

### Task 17: `vcf extract` (D73)

**Files:**
- crates/gloci-vcf/src/extract.rs
- crates/gloci-ops/src/vcf_extract.rs
- commands/vcf.rs, cli.rs
- Test: crates/gloci/tests/vcf_extract.rs, plus cases

- [ ] **Step 1: Write the tests.**
  - `--info DP,AF` with missing values and flags.
  - `--format GT,DP` in `--layout wide` and in `--layout long`.
  - `--genotype-matrix 012` files.
  - Selection flags apply.
  - BCF and VCF give the same output.
  - Determinism.
- [ ] **Step 2: Add differential cases:**
  - `vcftools --get-INFO DP` against `--info DP` (the column header is vcftools' layout `CHROM POS REF ALT DP`, so name the shared columns to match, or pin the difference);
  - `vcftools --extract-FORMAT-info DP` against `--format DP --layout wide`;
  - `vcftools --012` against `--genotype-matrix 012`.
- [ ] **Step 3: Run** `cargo test --workspace`.
- [ ] **Step 4: Commit** `feat: gloci vcf extract (INFO/FORMAT tables, 012 matrix)`

### Task 18: `vcf concat` (D74)

**Files:**
- crates/gloci-vcf/src/concat.rs
- crates/gloci-ops/src/vcf_concat.rs
- commands/vcf.rs, cli.rs
- Test: crates/gloci/tests/vcf_concat.rs, plus cases

- [ ] **Step 1: Write the tests.**
  - Concatenation of two per-chromosome files.
  - Sample mismatch exits 1 naming the first difference.
  - Header-definition mismatch exits 1, using P1's message.
  - Order across files violated exits 1 with a hint.
  - BCF inputs.
  - `naive_concat_block_edges`:
    - headers spanning more than one BGZF block;
    - an input without an EOF marker (accepted, with a warning);
    - output readable by `bcftools view`, and byte-equal after decompression to the non-naive output.
  - `--naive` with non-identical headers exits 2 with a hint.
- [ ] **Step 2: Add differential cases:** `bcftools concat a.vcf.gz b.vcf.gz` against `gloci vcf concat`, compared with `vcf-records`. Add `vcf-concat` from vcftools' perl scripts if it is in the image.
- [ ] **Step 3: Run** `cargo test --workspace`.
- [ ] **Step 4: Commit** `feat: gloci vcf concat (checked, and --naive block copy)`

### Task 19: Per-chromosome parallelism for plain seekable files (D76)

**Files:**
- crates/gloci-engine/src/seekable.rs (new: chromosome byte-range scan and range reader)
- crates/gloci-ops/src/relate/ (eligibility: every input indexed or plain seekable)
- Test: crates/gloci/tests/parallel_relations.rs

**Interfaces:**
- `pub fn chrom_ranges(path: &Path, format: Format, header_len: u64) -> Result<Vec<(Vec<u8>, Range<u64>)>, IoError>`. It is one memchr pass. Non-contiguous chromosomes give the streaming order error.
- `open_range_reader(path, range, &ReaderOptions)`, which yields the records in a byte range, with correct line numbers. The scan records each range's starting line number.

- [ ] **Step 1: Write the tests.**
  - `seekable_equals_streaming`, for the 5 relation commands on plain BED inputs, at `-t 2` and `-t 8`, compared byte for byte with `-t 1` streaming.
  - Mixed indexed and plain inputs.
  - Non-contiguous chromosomes exit 1, the same as streaming.
  - Error line numbers equal streaming's.
  - stdin is never treated as seekable.
- [ ] **Step 2: Implement** by reusing P1's per-chromosome driver (shared budget, ≤ 8 chunks in flight, grouping).
- [ ] **Step 3: Run** `cargo test --workspace`, and time a 2 GB plain BED intersect at `-t 1` and `-t 8` in a release build. Paste the timings.
- [ ] **Step 4: Commit** `feat(engine): per-chromosome parallelism for plain seekable files`

### Task 20: P2 benchmarks

**Files:** xtask/bench/cases.toml, xtask/src/gen.rs, xtask/docker/Dockerfile.bench (plink2), .github/workflows/nightly.yml.

Cases and targets:
- **`vcf stats` with 5 statistics** (`--freq --missing-site --missing-indv --het --tstv`) on the 1kg chr20 dataset and on the synthetic 2,500-sample VCF. The target is `gloci <= 0.1 * vcftools-5runs`, where the baseline command runs vcftools 5 times in sequence. Yardsticks: `bcftools stats` and `plink2 --freq --missing --het`.
- **New interval commands, `-t 1`, plain BED, 10M records:**
  - `partition` (against `bedops --partition`), `genomecov` (against bedtools), `cluster`, `groupby`, `shift`;
  - targets `gloci <= 0.5 * bedtools`, and `gloci <= bedops` where bedops has the operation.
- **`vcf filter --include`** (an expression) against `bcftools view -i`, as a yardstick only.

Every gloci command line must parse (P1's `every_gloci_bench_command_parses` covers new cases automatically). The sampled dry-run guard from Task 1 applies to the new VCF cases.

- [ ] **Step 1: Write the tests.** Bench cases parse, and the generator makes the 2,500-sample VCF with a missingness and depth mix that every requested statistic exercises.
- [ ] **Step 2: Implement.**
- [ ] **Step 3: Run** `cargo test -p xtask`, and run the synthetic `vcf stats` case locally (best effort on aarch64). Paste the medians.
- [ ] **Step 4: Commit** `feat(xtask): P2 benchmarks (vcf stats, partition, genomecov…)`

### Task 21: Behaviour matrices for P2 commands

**Files:** crates/gloci/tests/behaviour.rs and crates/gloci/tests/large.rs.

Extend P1's matrices to `shift`, `cluster`, `spacing`, `groupby`, `partition`, `genomecov`, `vcf stats`, `vcf extract` and `vcf concat`:
- determinism at -t 1/2/8;
- broken pipe;
- empty inputs;
- VCF input where it applies;
- the failure matrix: unsorted input, malformed lines, truncated BGZF, missing output directory;
- `--chrom-map` round trip.

Add RSS ceilings for streaming `partition`, `genomecov`, `cluster` and `vcf stats` on generated multi-GB inputs (`#[ignore]`).

- [ ] **Step 1: Write the matrices.**
- [ ] **Step 2: Run** `cargo test -p gloci --test behaviour`, and the large tests with `GLOCI_LARGE_BYTES=200000000`. Fix any real failure in the module that owns it, with a regression test.
- [ ] **Step 3: Commit** `test: P2 determinism, broken-pipe, empty-input and failure matrices; RSS ceilings`

### Task 22: Documentation

**Files:** docs/migration.md, docs/differences.md, SKILL.md (regenerated), xtask/src/docs.rs if new sections need it.

- [ ] **Step 1: Update `docs/migration.md`.** Turn every P2 row of the spec §11 tables into a working mapping, checked against `--help`:
  - bedtools: `multiinter`, `unionbedg`, `genomecov`, `cluster`, `spacing`, `groupby`, `shift`, and the full `map -o` ops;
  - bedops: `-i`, `-s`, `-p`, and `bedmap` ops;
  - vcftools: every `vcf stats` flag, `--get-INFO`/`--extract-FORMAT-info`/`--012`, and `vcf-concat`;
  - filter expressions: a short section mapping common bcftools `-i`/`-e` patterns, and stating the unsupported constructs.
- [ ] **Step 2: Regenerate SKILL.md** (`cargo xtask docs --skill SKILL.md`). The freshness test must pass. Add pitfalls for the expression subset and for `--chrom-map`.
- [ ] **Step 3: Run** `cargo test -p xtask` and the P0 help-examples test.
- [ ] **Step 4: Commit** `docs: P2a migration tables, expression guide and SKILL.md`

---

## Execution notes for the controller

- **Order.** Run Task 1 first, because it fixes P1 correctness. Tasks 2–5 are infrastructure. Tasks 6–12 are interval work. Tasks 13–18 are VCF work: 13 comes before 14 and 16, and 15 comes before 16. Tasks 19–22 come last.
- **Docker.** New reference tools (plink2, and UCSC `wigToBigWig` if it is used) go into the pinned images with pinned versions.
- **Spec amendments.** None are needed for D61–D78 except D78, which removes a documented limitation. The controller edits docs/differences.md through Task 3.
