# Coverage summaries on `gloci depth` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `gloci depth --summary` prints samtools coverage 1.21's per-reference table, and `gloci depth --out-prefix P [--by BED|SIZE] [--thresholds N,…]` writes mosdepth 0.3.10's files (names, formats, defaults), both streaming with flat memory and the same bytes at every `-t`.

**Architecture:** Both modes run on the existing depth core (`crates/gloci-ops/src/depth/`). Task 1 turns the core's text emitter into one implementation of a `DepthSink` trait, with hooks for per-read statistics. `--summary` is a sink that sums positions and reads per reference. The mosdepth mode adds mosdepth's mate-overlap rule and read filter to the core, plus a sink that keeps streaming accumulators per region, per threshold and per histogram. Its BGZF outputs get a CSI index through a new staged-output call in gloci-io, and all its files are published together.

**Tech Stack:** Rust (workspace crates gloci-io, gloci-ops, gloci, gloci-py, xtask). No new runtime dependency. Reference tools: samtools 1.21 (in `gloci-diff:p0` already) and mosdepth 0.3.10 (MIT; added to `gloci-diff:p0` by Task 6, pinned as in `xtask/docker/Dockerfile.bench`).

**Spec:** `docs/superpowers/specs/2026-11-17-tracks-query-coverage-liftover-design.md`: "Feature 3: coverage summaries on `gloci depth`", plus "Owner constraints", "Common to all four" and "Order of the work". Read with docs/decisions.md D102, D116, D126, D137–D139 and D87.

**Builds on:**
- the depth core: `DepthRing`, `CountRules` (ring.rs); `MateTable` (mates.rs, D139); `Plan`, `Input`, `Emitter`, `run_group` (rows.rs); `depth_sources`, `open_depth_inputs`, `depth_regions`, `check_same_references` (mod.rs);
- the command: `crates/gloci/src/commands/depth.rs`, `DepthArgs` in `crates/gloci/src/cli/alignments.rs`;
- gloci-io: `Output`, `open_output`, `OutputSpec`, `StagedSet`, `index/build.rs`;
- number formatting: `gloci_core::cfmt::{write_g, write_fixed}` (C `%g`/`%.Nf`, D145).

**Lane:** this plan is lane C's first half (liftover follows, in its own plan). Lanes A (bigWig/bigBed) and B (`vcf query`) run beside it in other worktrees. Tasks run in order 1 → 7; Task 4 has no dependency and may run anywhere before Task 5.

## Global Constraints

- **Owner constraints.**
  - gloci stays free for commercial use: no UCSC code, binary, library or licence anywhere. That covers gloci itself, its tests, `xtask/docker/Dockerfile` and `xtask/docker/Dockerfile.bench`.
  - A new runtime dependency must be MIT, Apache-2.0, BSD, zlib or similar, checked when it is added. This plan adds none.
  - mosdepth (MIT) is a test-only reference tool in the comparison image.
  - Never write the name of UCSC's code base (the four-letter word the owner banned) in code, comments, docs, tests or commit messages. Cite published format descriptions and tool behaviour instead.
- **Toolchain.** Rust 1.100.0 (`rust-toolchain.toml`), MSRV 1.91, so no std API newer than 1.91. `cargo clippy --workspace --all-targets -- -D warnings` and `cargo fmt --all -- --check` pass before every commit.
- **Help rules** (`crates/gloci/tests/cli_basics.rs`):
  - the summary is at most 64 characters and names the tools it replaces in brackets;
  - the long help's first paragraph names every tool the summary names;
  - 2 to 4 `gloci …` lines in an Examples block that ends at the first blank line;
  - the summary paragraph is at most 12 lines at 100 columns;
  - `depth` keeps its `EXPERIMENTAL` marker and its group in `crates/gloci/src/cli/mod.rs`. Neither file changes: this plan adds no command.
- **Output.**
  - D17 (format from `-o`'s extension).
  - Atomic outputs: `--summary` through `open_output`; the mosdepth files staged in one `StagedSet` and published together, so a failed run leaves the previous set untouched.
  - Output never changes with `--threads`.
  - Streaming with flat memory.
  - No silent wrong answers: errors name the line or the construct; exit 2 for usage errors, 1 for data errors.
- **Parity.**
  - Committed goldens:
    - `--summary` (stdout) through `xtask/diff/cases.toml` → `crates/gloci/tests/golden/` → `crates/gloci/tests/golden.rs`;
    - the mosdepth files through committed mosdepth outputs in `crates/gloci/tests/mosdepth/` (made by its `regen.sh` in `gloci-diff:p0`) → `crates/gloci/tests/depth_files.rs`. This is the same pattern as `crates/gloci-vcf/tests/golden/stats` → `crates/gloci/tests/vcf_stats.rs`.
  - Live cases run nightly (`cargo xtask diff`).
  - A case without `expected_difference` must match. One with it must differ, and names an anchor in `docs/differences.md` that has gloci/old/reason/pinned-by.
  - A cases.toml `gloci` line that is one plain `gloci ARGS` becomes a golden that `cargo test` runs *inside the committed golden directory*. So a case whose gloci side writes `--out-prefix` files must not be a plain invocation: use `… && cat …`, which makes it live-only. The one exception is a case that exits 2 before writing anything.
  - Regenerate a golden with `cargo xtask diff --only NAME --update-golden`; check it with `cargo xtask diff --only NAME`. `--only` takes one exact name.
- **Docs.**
  - D173 in `docs/decisions.md` (text below), appended after D170.
  - `docs/design.md` §7.4 and §11; `docs/migration.md` rows; `docs/differences.md` anchors; a `CHANGELOG.md` bullet under `[Unreleased]` → `### Added`.
  - After any help change: `cargo xtask pystub` and `cargo xtask docs --skill SKILL.md`. Check with `cargo test -p xtask --test pystub --test docs --test references --test book`.
  - `xtask/tests/references.rs` requires every cited `Dnnn` and spec `§` to exist. It forbids plan labels ("Task 3", "lane C") in code comments.
- **Shared files** (lanes A and B edit them too). Keep each edit small and self-contained: append entries, never reflow neighbours.
  - `CHANGELOG.md`, `docs/*.md`, `xtask/docker/Dockerfile` (one added line), `crates/gloci/src/embed/bindings.rs`, `xtask/diff/cases.toml` (one appended section).
  - Generated files (`SKILL.md`, `crates/gloci-py/python/gloci/__init__.pyi`, the book reference) are regenerated after a rebase, never merged by hand.
- **Tests.**
  - Run with `cargo test -p <crate> --test <file>`, or `cargo test -p <crate> --lib <filter>` for unit tests.
  - Large tests are `#[ignore]` and run nightly. Here, compile them only (`--no-run`).
  - Never run the full `--ignored` suite or build an image for another architecture on the laptop.
- **Disk.** Delete scratch BAMs, `target/diff-work/*` and timing data when a task ends. After rebuilding `gloci-diff:p0`, run `docker image prune -f`.

## Review Focus

1. **Overlapping proper pairs with multi-operation CIGARs and same-start mates.** Covers a deletion or intron inside the overlap, an insertion, soft clips, a supplementary of a pair, a non-proper pair and a mate on another reference. The expected depth is mosdepth's, base for base. Task 3 `pair_rule_matches_mosdepth_per_base`; Task 6 case `pairs-window-thresholds`.
2. **Means at a rounding boundary.** A window whose exact mean ends in 5 at the third decimal, such as 11/200 = 0.055, prints what mosdepth's base-by-base sum prints (`0.05`), not the exact quotient (`0.06`). A window mean of exactly .5 counts in region-dist bin 1. Task 5 `region_mean_accumulates_base_by_base`, `window_dist_rounds_half_away_from_zero`; Task 6 case `means-base-by-base`.
3. **References with no record, or only filtered records.**
   - mosdepth: a reference with only filtered records keeps its dist and summary lines; a reference with no record loses them.
   - samtools coverage: references no counted read reaches come last.
   - Tests: Task 5 `reference_without_records_is_left_out`; Task 2 `summary_lists_uncovered_references_last`.
4. **Huge headers.** 200,000 `@SQ` lines through `--summary` and `--by SIZE`, with nothing quadratic in the reference count. Task 2 `summary_huge_header`; Task 6 `by_windows_huge_header`.
5. **Failure and threads.**
   - A BED region past `LN`, or an unwritable prefix, publishes nothing and leaves the previous files as they were.
   - `-t 1` and `-t 4` write byte-identical files, `.bed.gz` and `.csi` included.
   - Tests: Task 6 `bed_region_past_end_publishes_nothing`, `threads_do_not_change_the_files`.

## Decision D173 (the text for `docs/decisions.md`)

Task 2 appends part A as a new paragraph after D170. Task 6 inserts part B's `<li>` items before part A's closing `</ul>`. Task 7 inserts part C the same way.

**Part A (Task 2):**

```markdown
<a id="d173"></a>**D173** **Coverage summaries on `depth`: `--summary` (samtools coverage) and mosdepth's files (`--out-prefix`)** (owner, 2026-11-17). <ul><li>**Modes.** `--summary` and `--out-prefix P` each replace the per-position rows and read through the depth core (D137, D138): one merged stream at every `-t`, threads decode, so no table or file changes with `-t`. Each conflicts with the other and with `--all-positions`, `--all-references`, `--count-deletions`, `--skip-mate-overlap`, `--column-names` and `--regions-file`.</li><li>**`--summary`** prints samtools coverage 1.21's table, `#rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq`: one row per reference, or one row for the one `--region`. A second `--region`, or a region on a reference the file lacks, is exit 2 (samtools exits 1: `#samtools-exit-codes`).<ul><li>Window: the reference, or the region with its end as given (it may pass `LN`; no end means `LN`). `startpos` is its 1-based start.</li><li>numreads: the records that pass the filters, mapped or not (in a region, those that overlap it, D116), summed over the inputs. meanmapq is their mean MAPQ.</li><li>covbases: the positions of the window whose depth, summed over the inputs, is above 0. meandepth: the sum of that depth over the window, divided by the window's length. Depth counts M, = and X bases whose quality is at least `--min-baseq`; a missing quality passes, and so does every base of a SEQ `*` read. D and N bases never count, there is no mate clipping, and reads are cut at `LN` (D137).</li><li>meanbaseq: the mean quality of the counted bases that have one. SEQ `*` reads have none; QUAL `*` counts as 255, as in samtools.</li><li>`--min-read-length` compares the CIGAR's query length (M, I, S, = and X), as samtools coverage `-l` does. `depth`'s own rule is D137's.</li><li>`--require-flags` keeps D102's meaning (all of the flags); samtools coverage `--rf` keeps reads with any of them: `#coverage-require-flags`.</li><li>Order: first the references a counted, mapped read reaches, in `@SQ` order, then the others in `@SQ` order. samtools prints the references its pileup never reached last.</li><li>Numbers: counts as integers, coverage and meandepth as C `%g`, meanbaseq and meanmapq as `%.3g`. A zero divisor prints 0.</li><li>Several FILEs give one table, as samtools coverage does with several inputs. Their `@SQ` lists must be identical (D126).</li><li>CRAM qualities are decoded. samtools coverage skips them unless `-Q` is given, and then prints 255: `#coverage-cram-baseq`.</li><li>Not provided: `-m`, `-A`, `-w`, `-D`, `-H`, `-b`, and `-d` (gloci never caps the depth; samtools' pileup stops at 1,000,000 reads).</li></ul></li></ul>
```

**Part B (Task 6), inserted before the final `</ul>`:**

```markdown
<li>**`--out-prefix P`** writes mosdepth 0.3.10's files, under its names and in its formats, for one FILE (exit 2 for more). It always writes `P.mosdepth.global.dist.txt` and `P.mosdepth.summary.txt`; with `--by`, also `P.mosdepth.region.dist.txt` and `P.regions.bed.gz`; with `--thresholds`, also `P.thresholds.bed.gz`. Each `.bed.gz` is BGZF with a `.csi`. There is no per-base file: that is `gloci depth`. The files are staged and published together. `--out-prefix` also conflicts with `--region`, `--min-baseq`, `--min-read-length`, `--require-flags`, `--summary` and `-o`. `--by` needs `--out-prefix`, and `--thresholds` needs `--by`.<ul><li>**Filters**: `--exclude-flags` (default `UNMAP,SECONDARY,QCFAIL,DUP`, mosdepth's 1796) and `--min-mapq` (default 0, mosdepth's). The depth core applies them after it notes each record's reference, and the reader filters nothing. The reason: mosdepth leaves a reference with no record at all out of the dist and summary files, but keeps one whose records all fail the filters, with depth 0.</li><li>**Depth**: M, = and X bases, with no quality gate; D and N bases never count. A read is cut at `LN` with one warning, where mosdepth stops: `#mosdepth-past-reference-end`.</li><li>**Mates** (mosdepth's rule): this applies to a read with PROPER_PAIR set and SUPPLEMENTARY clear. The read is held, with its aligned segments, when its mate is on the same reference with `end > mpos` and either `pos < mpos`, or `pos == mpos` with no read of its name held. Segments are runs of M, = and X; D and N end a run, while I, S, H and P do not. A later such read of the same name that is not held takes the entry, and its bases at positions the held segments cover are not counted. mosdepth's shortcut for two one-operation CIGARs subtracts up to the first mate's end even when the second mate ends before it: `#mosdepth-contained-mate`. An entry is dropped once a read starts past its mate position, and at a reference change, so the table grows with the depth only. mosdepth keeps entries to the end of the chromosome: `#mosdepth-mate-table`.</li><li>**Regions**: `--by SIZE` (digits, at least 1) makes windows `[k·SIZE, min((k+1)·SIZE, LN))` on every reference. `--by FILE` reads a BED with gloci's reader:<ul><li>comment, `track`, `browser` and blank lines are skipped; mosdepth stops reading at a blank line: `#mosdepth-bed-blank-line`;</li><li>a non-empty fourth column is the region's name;</li><li>each reference's regions are sorted by start, keeping file order for ties, and references follow `@SQ` order;</li><li>a BED chromosome the file lacks is skipped with one warning;</li><li>a region that ends past `LN` is exit 1 naming the line, where mosdepth crashes.</li></ul></li><li>**Means**: a region's mean is the sum over its bases, in order, of depth / length in f64, as mosdepth adds it; the printed rounding depends on this. It prints as `%.2f`; a zero-length region prints `0.00`.</li><li>**Dists**: the histograms are per base; a depth above 400,000 counts at 399,990, as in mosdepth. The region dist counts each base of the BED regions, a base in two regions twice, or each window by its mean rounded half away from zero, capped at 511. A dist is written from the highest depth down to 0 as the cumulative fraction, `%.2f`, skipping depths above 300 with no count and lines whose fraction is below 8e-5. The `total` lines come from the sum over references.</li><li>**Summary**: the header is `chrom length bases mean min max`. Each reference gets a `name` line and, with `--by`, a `name_region` line over the region bases (a base in two regions counted twice). Then come `total` and, with `--by`, `total_region`. A minimum over no bases prints 0. The mean is bases / length as `%.2f`, or 0 when the length is 0.</li><li>**A reference with no record** gets region rows of `0.00` and threshold counts of 0, and no dist or summary lines, as in mosdepth.</li><li>**Thresholds**: sorted ascending, duplicates kept, under the header `#chrom start end region NX…`. A region without a name is `unknown`. A threshold of 0 is exit 2: mosdepth counts 0X but writes 0 on references without records (`#mosdepth-thresholds-zero`).</li><li>**Index**: a CSI with min_shift 14 and the BED preset (0-based), as `gloci index` writes it. mosdepth's `.csi` has the generic 1-based preset, so tabix returns one more region at a boundary: `#mosdepth-csi-preset`.</li><li>**Memory**: the depth ring and the pair table (both grow with the depth only), the BED's regions with one accumulator per region of the current reference, and the histograms. There are no per-base arrays.</li><li>**Not provided**: mosdepth's per-base and quantized files, `--use-median`, `--fast-mode`, `--chrom`, `--read-groups`, `--include-flag`, `--min-frag-len`/`--max-frag-len` and `MOSDEPTH_PRECISION`. mosdepth needs an indexed BAM or CRAM; gloci reads any BAM, CRAM or SAM stream.</li></ul></li>
```

**Part C (Task 7), inserted before the final `</ul>`:**

```markdown
<li>**Python** (D87): `depth(summary=True)` is a step whose table is headed by its first line. `depth_files(src, out=P, by=…, thresholds=…)` is a module function like `vcf_stats`, and the `depth` step refuses `out_prefix`, `by` and `thresholds`.</li>
```

---

### Task 1: The depth core writes to a sink, with read hooks

Pure refactor: every existing `depth` and `bedcov` output stays byte-identical.

**Files:**
- Modify: `crates/gloci-ops/src/depth/rows.rs`: the `DepthSink` trait; `Emitter` implements it; `Printer`/`run_group` become generic over it; the `Plan` knobs; the hooks in `Input::advance`.
- Modify: `crates/gloci-ops/src/depth/mod.rs`: `run_sources`, `ModeRules`; `depth_sources` keeps its signature and calls `run_sources`; `#[cfg(test)] pub(crate) mod test_sink;`.
- Create: `crates/gloci-ops/src/depth/test_sink.rs` (`#[cfg(test)]`). It holds:
  - `Recorder`, a `DepthSink` that stores positions and zeros as merged runs `(rid, start, end, depth)`, and the `reference_done`, `read` (input, rid, read name, span) and `placed` events in order;
  - `fn sam_source(text: &str, ro: &ReaderOptions) -> Box<dyn AlnSource>`, built on `gloci_io::open_aln_source_from`.

  The unit tests of Tasks 1, 3 and 5 use both.
- Test: a `#[cfg(test)] mod tests` in `crates/gloci-ops/src/depth/rows.rs`. The trait is `pub(crate)`, so the tests live in the crate. In every test SAM, RNEXT is `*`, PNEXT and TLEN are 0, and SEQ and QUAL are `*` unless stated.

**Interfaces:**
- Consumes: the existing core (see "Builds on").
- Produces (Tasks 2, 3, 5 use these; all `pub(crate)`):

```rust
/// Where a depth run's positions go: `depth`'s rows (`Emitter`), the `--summary` table, or
/// mosdepth's files. Positions arrive in `@SQ` order, then position order, each at most once.
pub(crate) trait DepthSink {
    /// The depths at `pos` (0-based) of reference `rid`, one per input in input order.
    fn position(&mut self, rid: usize, pos: u64, depths: &[u32]) -> io::Result<()>;
    /// Positions `[a, b)` of `rid`, where every input's depth is 0.
    fn zeros(&mut self, rid: usize, a: u64, b: u64) -> io::Result<()>;
    /// Every position of `rid` that this run gives has been given.
    fn reference_done(&mut self, _rid: usize) -> io::Result<()> { Ok(()) }
    /// A record of input `input` on `rid` that passed every filter (reader, plan filter, read
    /// length), whether or not it covers a base. `span` is its D101 span cut at `LN`, or `None`
    /// when it is unmapped or the cut span is empty.
    fn read(&mut self, _input: usize, _rid: usize, _a: &AlnRef<'_>, _span: Option<Interval>) {}
    /// A record placed on `rid` (reference and position set), before the plan's filter.
    fn placed(&mut self, _rid: usize) {}
    fn flush(&mut self) -> io::Result<()> { Ok(()) }
}

pub(crate) enum ReadLength { SeqLessClips /* D137, depth -l */, CigarQuery /* M I S = X: coverage -l */ }
pub(crate) enum MateRule { Off, ClipAtMateEnd /* D139 */ }       // Task 3 adds SharedBases
pub(crate) struct ReadFilter { pub exclude_flags: u16, pub min_mapq: u8 }
impl ReadFilter { pub fn passes(&self, a: &AlnRef<'_>) -> bool }  // no excluded bit set and MAPQ >= min
// Plan gains: read_length: ReadLength, mates: MateRule (replaces skip_mates), filter: Option<ReadFilter>

/// What a mode sets on the plan beyond `DepthOptions`. `ModeRules::depth()` = SeqLessClips, no filter.
#[derive(Clone, Copy)]
pub(crate) struct ModeRules { pub read_length: ReadLength, pub filter: Option<ReadFilter> }

/// Check the inputs' headers (D126), reopen them through the regions when `o.regions` is set (as
/// `depth_sources` does today), build the sink from the shared header, run every reference
/// (`run_group`) and log the diagnostics. Returns the summary and the sink.
pub(crate) fn run_sources<S: DepthSink>(
    opened: Vec<Box<dyn AlnSource>>,
    paths: &[String],
    o: &DepthOptions,
    mode: ModeRules,
    ro: &ReaderOptions,
    make_sink: impl FnOnce(&Arc<AlnHeader>) -> Result<S, IoError>,
) -> Result<(DepthSummary, S), IoError>
```

- `run_group` calls `reference_done(r)` once per reference it prints, after that reference's last position. That covers each `pending` flush to `LN` at a reference change, the flush at the end, and each `whole(r)` under `AllPositions::AllRefs`.
- `Emitter` keeps today's formatting code (`RowHead`, `FLUSH_AT`). It takes the reference names and the column count at construction, and implements `position` as `row_start` + one `column` per depth + `row_end`.
- `Printer` is generic (`<S: DepthSink>`, monomorphised): no `dyn` call on the text path.
- Hook order inside `Input::advance`, for each record with a reference and position:
  1. the order check (unchanged) and `last` update;
  2. `sink.placed(rid)`;
  3. `plan.filter` (skip when it fails);
  4. the read-length check per `plan.read_length` (skip when short);
  5. compute the span as today, then `sink.read(k, rid, &a, span)`, *before* the existing `continue`s for a missing or empty span;
  6. the existing `cur` logic.

- [ ] **Step 1: Write the failing tests** in rows.rs `mod tests`. Each opens SAM text with `sam_source` and runs `run_sources` with a `Recorder`.
  - `sink_gets_every_position_once_in_order`. The header has `c1` LN 20, `c2` LN 5, `c3` LN 8; the reads are `r1 0 c1 3 60 3M2D2M`, `r2 0 c1 4 60 4M` and `r3 0 c3 2 60 2M`; the mode is `AllPositions::AllRefs`.
    - The runs cover `c1 [0,20)`, `c2 [0,5)` and `c3 [0,8)` exactly once, in order.
    - The c1 depths at 2..9 are `1,2,2,1,1,1,1`.
    - `reference_done` comes as `[0, 1, 2]`, each after its reference's last run.
  - `read_hook_sees_each_passing_record`. Header `c1` LN 100, `c2` LN 50; the reader's `exclude_flags` is 0.
    - Reads `m 0 c1 3 60 45M` and `u 4 c1 7 0 *` (placed unmapped), with `min_read_length` 0: the `read` events (name, span) are `m Some([2,47))` and `u None`.
    - Reads `long 0 c1 1 60 30S30M30S` and `noseq 0 c1 1 60 50M`, with `min_read_length` 60:
      - under `CigarQuery`, the only event is `long Some([0,30))` (query lengths 90 and 50);
      - under `SeqLessClips`, there is no event (lengths 30 and 50).
  - `placed_hook_runs_before_the_plan_filter`. `ReadFilter { exclude_flags: 0x400, min_mapq: 20 }`, reads `d 1024 c2 1 60 5M` and `q 0 c2 3 5 5M`: `placed(1)` comes twice, and `read` never comes.
- [ ] **Step 2: Run** `cargo test -p gloci-ops --lib depth::rows`. Expected: FAIL to compile (`DepthSink` and `run_sources` do not exist).
- [ ] **Step 3: Implement** the interfaces above. `depth_sources` writes the `--column-names` line, then calls `run_sources` with an `Emitter` over `out`, and returns `DepthSummary` as before.
- [ ] **Step 4: Run** `cargo test -p gloci-ops --lib depth::rows`, `cargo test -p gloci-ops --test depth_core --test depth_mates --test depth_threads` and `cargo test -p gloci --test depth --test bedcov --test golden`. Expected: PASS, with no golden changed.
- [ ] **Step 5: Check the text path's speed.**
  - Before changing anything, build the base commit with `cargo build --release -p gloci` and copy the binary to `$SCRATCH/gloci-before`.
  - Generate the input: `cargo xtask gen sam --records 2000000 --chroms 4 --seed 1 -o $SCRATCH/p.sam`, then `$SCRATCH/gloci-before view $SCRATCH/p.sam -o $SCRATCH/p.bam`.
  - Time `depth $SCRATCH/p.bam > /dev/null` 5 times with each binary.
  - Expected: the new binary's median is within 3% of the old one's, and the two outputs are byte-identical.
  
  Delete the files afterwards.
- [ ] **Step 6: Commit** `refactor(depth): positions go to a sink, with read hooks for summaries`

### Task 2: `gloci depth --summary` (samtools coverage)

**Files:**
- Create: `crates/gloci-ops/src/depth/coverage.rs`
- Modify: `crates/gloci-ops/src/depth/mod.rs` (`pub mod coverage;` and its re-exports), `crates/gloci/src/cli/alignments.rs` (`DepthArgs`, help), `crates/gloci/src/commands/depth.rs` (`run_summary`)
- Fixture: `xtask/diff/inputs/samtools/mosdepth.{sam,bam,bam.bai}`, built by a block appended to `make-samtools.sh` (below)
- Test: `crates/gloci/tests/depth_summary.rs`; goldens appended to `xtask/diff/cases.toml`
- Docs: D173 part A; `docs/differences.md` sections `#coverage-cram-baseq` and `#coverage-require-flags`, plus the new case appended to the "Pinned by" list of `#samtools-exit-codes`; `docs/migration.md` samtools table row; `SKILL.md` and the pystub regenerated

**Interfaces:**
- Consumes: `run_sources`, `ModeRules`, `ReadLength::CigarQuery`, `DepthSink` (Task 1).
- Produces:

```rust
pub struct CoverageOptions {
    pub min_baseq: u8,
    /// Skip reads whose CIGAR query length (M I S = X) is below this (samtools coverage -l).
    pub min_read_length: u64,
    /// The one region: reference name and window (0-based start; end as given, or LN when open).
    pub region: Option<(Vec<u8>, Interval)>,
}
pub struct CoverageSummary { pub rows: u64, pub max_ring: usize }
pub fn coverage_summary(opened: Vec<Box<dyn AlnSource>>, paths: &[String], o: &CoverageOptions,
    ro: &ReaderOptions, out: &mut dyn Write) -> Result<CoverageSummary, IoError>
```

- `coverage_summary` runs the core with:
  - `DepthOptions { rules: CountRules { deletions: false, ref_skips: false, min_baseq }, min_read_length, all: Covered, regions: region.map(|r| vec![r]), skip_mate_overlap: false, column_names: false }`;
  - `ModeRules { read_length: CigarQuery, filter: None }`.
- Its sink holds one record per reference (`Vec`, sized from the header): numreads, MAPQ sum, quality sum, quality bases, covered bases, depth sum, and "covered" (a `read` hook with `Some(span)` was seen).
  - `position`: add the depths' sum; count the base when the sum is above 0.
  - `read`: numreads += 1, MAPQ sum += MAPQ; when `span` is `Some`, mark the reference covered and add `quality_in(a, span, window, min_baseq)`.
- `fn quality_in(a: &AlnRef<'_>, span: Interval, window: Interval, min_baseq: u8) -> (u64, u64)` returns the quality sum and count of the M/=/X bases at positions inside `span ∩ window` that pass the gate (`q >= min_baseq` or `q == 0xFF`, the same gate as `DepthRing`). A read with `seq_len() == 0` gives `(0, 0)`.
- The rows are written at the end:
  - without a region: covered references in `@SQ` order, then uncovered ones;
  - with a region: the one row, whatever was counted.
  - Each row has the fields of D173 part A. Use `cfmt::write_g(out, 100.0 * covbases as f64 / len, 6)`, `write_g(depth_sum as f64 / len, 6)`, and `write_g(…, 3)` for the two means, where `len = (end - start) as f64`.
  - The header line is always first.
- `DepthArgs` gains `pub summary: bool`:
  - `#[arg(long, conflicts_with_all = ["all_positions", "all_references", "count_deletions", "skip_mate_overlap", "column_names", "regions_file"])]`;
  - doc: "Print one row per reference instead of one per position, as samtools coverage: reads, covered bases, mean depth, mean base and mapping quality. With --region, one row for that region."
- The help:
  - summary: `Depth of alignments (samtools depth, coverage)`;
  - the long help's first line: "Depth of coverage of alignments at each position (samtools depth), or as a table per reference (samtools coverage)."
  - Examples: the current first three, then `gloci depth sample.bam --summary` in place of the `--regions-file` example;
  - a new "Summary:" section after "Positions and regions:" that says: one table for several FILEs; references no counted read reaches come last; `--min-read-length` counts soft clips here; D and N never count.
- `run_summary` (commands/depth.rs):
  - more than one `--region` is `UsageError("--summary takes one --region, as samtools coverage -r")`;
  - parse the region with `Region::parse` against the first input's header; an unknown reference is a `UsageError` naming the region and the file;
  - window = `Interval { start: r.start, end: r.end.unwrap_or(LN) }`;
  - output through `own_tsv_output_spec` and `open_output` as `run` does.

**Fixture block** (append to `make-samtools.sh`; the SAM is written by the script, then converted like `mates.bam`):

```sh
# mosdepth.sam: overlapping proper pairs for mosdepth's mate rule, and references with no reads or
# only filtered reads (the depth --summary and --out-prefix cases). chr1: p partial overlap; x MAPQ 0;
# q both mates at 201, the second longer; r a deletion inside the overlap; s an intron inside it;
# t soft clips and an insertion; u an overlapping pair without PROPER_PAIR; v a pair and a
# supplementary (2147); w a proper pair with its mate on chr2; z the second mate MAPQ 5; y1-y4 dup,
# secondary, qcfail, placed unmapped. chr2: a duplicate and a placed unmapped read only. chr3: none.
# chr4: 320 one-base reads at 5001 and two reads overlapping by 10 (distribution skips above 300).
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:300\n@SQ\tSN:chr3\tLN:200\n@SQ\tSN:chr4\tLN:20000\n'
  sed 's/ \{1,\}/\t/g' <<'EOF'
p 99 chr1 101 60 50M = 131 80 * *
x 0 chr1 110 0 30M * 0 0 * *
p 147 chr1 131 60 50M = 101 -80 * *
q 99 chr1 201 60 40M = 201 60 * *
q 147 chr1 201 60 60M = 201 -60 * *
r 99 chr1 301 60 20M10D20M = 321 60 * *
r 147 chr1 321 60 40M = 301 -60 * *
s 99 chr1 401 60 10M50N10M = 411 70 * *
s 147 chr1 411 60 60M = 401 -70 * *
t 99 chr1 501 60 5S20M2I20M = 521 50 * *
t 147 chr1 521 60 30M5S = 501 -50 * *
u 97 chr1 601 60 30M = 611 40 * *
u 145 chr1 611 60 30M = 601 -40 * *
v 99 chr1 701 60 50M = 721 70 * *
v 147 chr1 721 60 50M = 701 -70 * *
v 2147 chr1 731 60 20M = 701 0 * *
w 99 chr1 801 60 40M chr2 101 0 * *
z 99 chr1 851 60 40M = 871 60 * *
z 147 chr1 871 5 40M = 851 -60 * *
y1 1024 chr1 901 60 20M * 0 0 * *
y2 256 chr1 905 60 20M * 0 0 * *
y3 512 chr1 910 60 20M * 0 0 * *
y4 4 chr1 920 0 * * 0 0 * *
d1 1024 chr2 50 60 20M * 0 0 * *
d2 4 chr2 80 0 * * 0 0 * *
EOF
  for i in $(seq 1 320); do printf 'deep%03d\t0\tchr4\t5001\t60\t1M\t*\t0\t0\t*\t*\n' "$i"; done
  printf 'w1\t0\tchr4\t6001\t60\t20M\t*\t0\t0\t*\t*\nw2\t0\tchr4\t6011\t60\t20M\t*\t0\t0\t*\t*\n'
} > mosdepth.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o mosdepth.bam mosdepth.sam
samtools index -@ 0 mosdepth.bam
```

- [ ] **Step 1: Build the fixture.** Run the block in the image (`docker run --rm -u "$(id -u):$(id -g)" -v "$PWD/xtask/diff/inputs:/work" -w /work/samtools gloci-diff:p0 sh make-samtools.sh`), then `git status`. Expected: only `mosdepth.{sam,bam,bam.bai}` are new, and no other fixture changed.
- [ ] **Step 2: Write the failing tests** in `crates/gloci/tests/depth_summary.rs`. Fixtures come from `xtask/diff/inputs/samtools`; `H` is the header line `#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq`. The expected rows are samtools 1.21's.
  - `summary_of_depth_bam`: `depth depth.bam --summary` gives `H`, then `chr1\t1\t500\t15\t190\t38\t0.546\t47.4\t55.9`, `chr2\t1\t300\t3\t60\t20\t0.2\t40\t60` and `chr3\t1\t200\t0\t0\t0\t0\t0\t0`.
  - `summary_sums_several_inputs`: `depth.bam depth2.bam --summary` gives `chr1\t1\t500\t19\t247\t49.4\t0.736\t45.5\t56.8`, `chr2\t1\t300\t4\t90\t30\t0.3\t40\t60` and `chr3\t1\t200\t1\t10\t5\t0.05\t40\t60`.
  - `summary_lists_uncovered_references_last`:
    - `mosdepth.bam --summary` gives rows in the order chr1, chr4, chr2, chr3: `chr1\t1\t1000\t19\t530\t53\t0.76\t0\t53.9`, `chr4\t1\t20000\t322\t31\t0.155\t0.018\t0\t60`, `chr2\t1\t300\t0\t0\t0\t0\t0\t0`, `chr3\t1\t200\t0\t0\t0\t0\t0\t0`.
    - With `--exclude-flags 0`: `chr1\t1\t1000\t23\t549\t54.9\t0.82\t0\t52.4`, `chr2\t1\t300\t2\t20\t6.66667\t0.0666667\t0\t30`, `chr4\t1\t20000\t322\t31\t0.155\t0.018\t0\t60`, `chr3\t1\t200\t0\t0\t0\t0\t0\t0`.
  - `summary_region_window_as_given`, on `depth.bam`:
    - `--region chr1:100-200` gives the one row `chr1\t100\t200\t1\t3\t2.9703\t0.029703\t40\t60`;
    - `--region chr1:600-700` gives `chr1\t600\t700\t0\t0\t0\t0\t0\t0`;
    - `--region chr1:450` gives `chr1\t450\t500\t0\t0\t0\t0\t0\t0`;
    - `--region chr3` gives `chr3\t1\t200\t0\t0\t0\t0\t0\t0`.
  - `summary_min_read_length_counts_soft_clips`: `depth.bam --summary --min-read-length 50`; the chr1 row is `chr1\t1\t500\t2\t80\t16\t0.16\t40\t60`.
  - `summary_seq_star_reads_count_depth_not_quality`: on `depth-noseq.bam`, the chr1 row is `chr1\t1\t500\t3\t89\t17.8\t0.238\t40\t60`; with `--min-read-length 50` it is `chr1\t1\t500\t1\t50\t10\t0.1\t0\t60`.
  - `summary_qual_star_counts_255`: a temp SAM with depth.sam's header and its `p1`, `p1` and `nq` lines; the chr1 row is `chr1\t1\t500\t3\t35\t7\t0.1\t83\t60`.
  - `summary_filters`:
    - `--min-mapq 30 --min-baseq 30` gives chr1 `chr1\t1\t500\t14\t189\t37.8\t0.524\t47.9\t57.9`;
    - `--exclude-flags 0` gives chr1 `chr1\t1\t500\t20\t220\t44\t0.606\t46.7\t51.2`;
    - on `depth-pastend.bam`: chr1 `chr1\t1\t500\t1\t20\t4\t0.04\t40\t60` and chr2 `chr2\t1\t300\t3\t31\t10.3333\t0.106667\t40\t60`.
  - `summary_cram_reads_qualities`: `depth.cram --summary` gives the same rows as `summary_of_depth_bam` (meanbaseq `47.4`, not samtools' `255`).
  - `summary_usage_errors`, each exit 2:
    - two `--region`s (stderr contains `one --region`);
    - `--region chrZ:1-10` (stderr contains `chrZ`);
    - `--summary` with `--all-positions`, `--all-references`, `--count-deletions`, `--skip-mate-overlap`, `--column-names` or `--regions-file x.bed`.
  - `summary_threads_identical`: `depth.bam depth2.bam --summary` gives equal stdout at `-t 1` and `-t 4`.
  - `summary_writes_o_file`: `-o t.tsv` writes the same bytes as stdout; `-o t.tsv.gz` is BGZF.
  - `summary_huge_header`:
    - input: a temp SAM with 200,000 `@SQ` lines (`c0`…`c199999`, LN 10) and one read `r 0 c199999 1 60 5M`;
    - output: 200,001 lines; the first row is `c199999\t1\t10\t1\t5\t50\t0.5\t0\t60`, followed by `c0\t1\t10\t0\t0\t0\t0\t0\t0` and the rest in order;
    - it finishes in under 20 s in a debug build.
- [ ] **Step 3: Run** `cargo test -p gloci --test depth_summary`. Expected: FAIL (`--summary` is an unknown argument).
- [ ] **Step 4: Implement** coverage.rs, the CLI flag and help, and `run_summary`.
- [ ] **Step 5: Run** `cargo test -p gloci --test depth_summary --test depth --test cli_basics`. Expected: PASS.
- [ ] **Step 6: Add goldens.** Append a `# depth --summary (D173) vs samtools coverage 1.21` section to `xtask/diff/cases.toml`.
  - Each case below has `input` = the named fixture under `xtask/diff/inputs/samtools/`, `gloci = "gloci depth input.<ext> --summary …"` and `reference = "samtools coverage … input.<ext>"`, with the flags mapped as listed.
  - Exact cases (`compare = "exact"`):
    - `depth-summary-vs-samtools` (depth.bam);
    - `depth-summary-mates-vs-samtools` (mates.bam);
    - `depth-summary-order-vs-samtools` (mosdepth.bam);
    - `depth-summary-exclude-flags-zero-vs-samtools` (mosdepth.bam, `--exclude-flags 0` / `--ff 0`);
    - `depth-summary-two-files-vs-samtools` (aux `depth2.bam`; `input.bam depth2.bam` on both sides);
    - `depth-summary-sam-vs-samtools` (depth.sam);
    - `depth-summary-region-vs-samtools` (aux `depth-idx/input.bai`, `--region chr1:100-200` / `-r chr1:100-200`), and likewise `-region-whole` (`chr3`), `-region-open-end` (`chr1:450`) and `-region-past-end` (`chr1:600-700`);
    - `depth-summary-min-read-length-vs-samtools` (`--min-read-length 50` / `-l 50`);
    - `depth-summary-noseq-vs-samtools` and `depth-summary-noseq-min-read-length-vs-samtools` (depth-noseq.bam);
    - `depth-summary-mapq-baseq-vs-samtools` (`--min-mapq 30 --min-baseq 30` / `-q 30 -Q 30`);
    - `depth-summary-require-flags-vs-samtools` (`--require-flags PROPER_PAIR` / `--rf PROPER_PAIR`);
    - `depth-summary-pastend-vs-samtools` (depth-pastend.bam);
    - `depth-summary-empty-vs-samtools` (empty.bam);
    - `depth-summary-cram-min-baseq-vs-samtools` (depth.cram, `--min-baseq 1` / `-Q 1`).
  - Pinned cases:
    - `depth-summary-cram-vs-samtools` (depth.cram, no flags): `expected_difference = "#coverage-cram-baseq"`;
    - `depth-summary-require-any-vs-samtools` (`--require-flags PAIRED,SUPPLEMENTARY` / `--rf PAIRED,SUPPLEMENTARY`): `#coverage-require-flags`;
    - `depth-summary-unknown-region-vs-samtools` (`--region chrZ:1-10` / `-r chrZ:1-10`, `compare = "exit-status"`): `#samtools-exit-codes`.

  Write the two new sections in `docs/differences.md`:
  - `#coverage-cram-baseq`: gloci decodes CRAM qualities and prints their mean; samtools coverage asks htslib to skip QUAL unless `-Q` is given, and prints 255. Pinned by the case and `summary_cram_reads_qualities`.
  - `#coverage-require-flags`: gloci's `--require-flags` needs all of the flags (D102, as everywhere in gloci); samtools coverage `--rf` keeps a read with any of them. Pinned by the case.

  Then:
  - `for c in $(grep -o '^name = "depth-summary-[^"]*"' xtask/diff/cases.toml | cut -d'"' -f2); do cargo xtask diff --only "$c" --update-golden || break; done`;
  - the same loop without `--update-golden`. Expected: every case `ok`; the three pinned ones say `differs as documented`.
- [ ] **Step 7: Docs.**
  - Append D173 part A to `docs/decisions.md`.
  - In `docs/migration.md`'s samtools table, add a row after the `samtools depth` row: `samtools coverage [-r REG] [-l N] [-q N] [-Q N] [--rf F] [--ff F] [-o OUT] FILE...` → `gloci depth FILE... --summary [--region REG] [--min-read-length N] [--min-mapq N] [--min-baseq N] [--require-flags F] [--exclude-flags F] [-o OUT]`. Its notes: one table for several files, as samtools; `-q` is `--min-mapq` and `-Q` is `--min-baseq`; `-m`, `-A`, `-w`, `-D`, `-H`, `-b` and `-d` have no counterpart; links to the two new anchors.
  - Run `cargo xtask pystub` and `cargo xtask docs --skill SKILL.md`.
- [ ] **Step 8: Run** `cargo test -p gloci --test golden --test depth_summary --test cli_basics` and `cargo test -p xtask --test pystub --test docs --test references --test book`. Expected: PASS.
- [ ] **Step 9: Commit** `feat(depth): --summary prints samtools coverage's table (D173)`

### Task 3: mosdepth's mate rule and read filter in the core

**Files:**
- Create: `crates/gloci-ops/src/depth/pairs.rs` (`PairTable`, `aligned_segments`)
- Modify: `crates/gloci-ops/src/depth/ring.rs` (`add_except`), `crates/gloci-ops/src/depth/rows.rs` (`MateRule::SharedBases`, `Input` holds a `PairTable`), `crates/gloci-ops/src/depth/mod.rs` (`ModeRules` gains `shared_bases: bool`; `DepthSummary` reports pending pairs through the existing `max_pending_mates`)
- Fixture: `xtask/diff/inputs/samtools/mosdepth-stale.{sam,bam,bam.bai}`, built by a block appended to `make-samtools.sh`
- Test: `#[cfg(test)] mod tests` in `crates/gloci-ops/src/depth/pairs.rs`

**Interfaces:**
- Consumes: `run_sources`, `ModeRules`, `ReadFilter`, the `Recorder` test sink (Task 1); `mosdepth.bam` and `mates.bam`.
- Produces:

```rust
/// The aligned segments of `r` from `pos`: maximal runs of M, = and X positions (I, S, H and P do
/// not end a run, D and N do), cut at `end`. Written to `out`, which is cleared first.
pub fn aligned_segments(r: &AlnRef<'_>, pos: u64, end: u64, out: &mut Vec<Interval>);

/// mosdepth 0.3.10's pairs (D173): reads waiting for their mate, with their aligned segments.
#[derive(Default)]
pub struct PairTable { /* name -> (mpos, segments); min-heap on mpos, as MateTable */ }
impl PairTable {
    pub fn new() -> PairTable;
    /// The segments `r` must not count again, written to `out` (cleared first; empty when none).
    /// Only a read with PROPER_PAIR (0x2) set and SUPPLEMENTARY (0x800) clear takes part. It is
    /// held (replacing an entry of its name) when its mate is on its reference, `end > mpos`, and
    /// `pos < mpos` or (`pos == mpos` and no entry of its name is held). Otherwise it takes the
    /// entry of its name, if any, and `out` gets that entry's segments.
    pub fn shared_with_mate(&mut self, r: &AlnRef<'_>, pos: u64, end: u64, out: &mut Vec<Interval>);
    /// Drop the entries whose mate position is below `pos`.
    pub fn advance(&mut self, pos: u64);
    pub fn reset(&mut self);
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn peak(&self) -> usize;
}

// DepthRing:
/// As `add`, but bases at positions inside any of `except` (sorted, disjoint) are not counted.
pub fn add_except(&mut self, r: &AlnRef<'_>, span: Interval, rules: &CountRules, except: &[Interval]);
```

- In `Input::add_current` under `MateRule::SharedBases`:
  1. `pairs.shared_with_mate(&a, c.pos, c.span.end, &mut buf)`;
  2. then `pairs.advance(c.pos)`;
  3. then `ring.add_except(&a, c.span, &plan.rules, &buf)`.

  This is the same consume-or-store-then-drop order as D139. `run_group` resets the table at each reference change. `ModeRules { shared_bases: true }` sets `MateRule::SharedBases`. `ModeRules::depth()` and the `--summary` call site (Task 2) set `shared_bases: false`.
- Re-export `PairTable` and `aligned_segments` from `depth/mod.rs` beside `MateTable`.

**Fixture block** (append to `make-samtools.sh`):

```sh
# mosdepth-stale.sam: s waits for its mate (mpos 121), whose record is a duplicate; o starts past 121;
# a secondary alignment of s (355) at 301. With -F 1540 mosdepth lets the secondary take the stale
# entry and adds depth over [150,300); gloci dropped the entry at o (#mosdepth-mate-table).
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
  printf 's\t99\tchr1\t101\t60\t50M\t=\t121\t70\t*\t*\n'
  printf 's\t1171\tchr1\t121\t60\t50M\t=\t101\t-70\t*\t*\n'
  printf 'o\t0\tchr1\t131\t60\t10M\t*\t0\t0\t*\t*\n'
  printf 's\t355\tchr1\t301\t60\t50M\t=\t101\t0\t*\t*\n'
} > mosdepth-stale.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o mosdepth-stale.bam mosdepth-stale.sam
samtools index -@ 0 mosdepth-stale.bam
```

- [ ] **Step 1: Build the fixture** (as in Task 2, Step 1). Expected: only the three `mosdepth-stale.*` files are new.
- [ ] **Step 2: Write the failing tests.** Runs are `(ref, start, end, depth)` from the `Recorder`, under `AllPositions::AllRefs`, `CountRules::default()` and `ModeRules { read_length: SeqLessClips, filter: Some(ReadFilter { exclude_flags: 1796, min_mapq: 0 }), shared_bases: true }`. The reader filters nothing (`exclude_flags: 0`).
  - `aligned_segments_merge_across_insertions`:
    - `5S20M2I20M` at 500 → `[500,540)`;
    - `20M10D20M` at 300 → `[300,320) [330,350)`;
    - `10M50N10M` at 400 → `[400,410) [460,470)`;
    - `50M` at 0 cut at end 30 → `[0,30)`.
  - `add_except_skips_the_mate_bases`: `40M` at 320 with except `[300,320) [330,350)` gives depths 1 at 320..330, 0 at 330..350 and 1 at 350..360.
  - `pair_rule_matches_mosdepth_per_base`: `mosdepth.bam`'s runs equal mosdepth 0.3.10's per-base output:
    - chr1: `0-100:0 100-109:1 109-139:2 139-180:1 180-200:0 200-260:1 260-300:0 300-360:1 360-400:0 400-470:1 470-500:0 500-550:1 550-600:0 600-610:1 610-630:2 630-640:1 640-700:0 700-730:1 730-750:2 750-770:1 770-800:0 800-840:1 840-850:0 850-910:1 910-1000:0`;
    - chr2: `0-300:0`; chr3: `0-200:0`;
    - chr4: `0-5000:0 5000-5001:320 5001-6000:0 6000-6010:1 6010-6020:2 6020-6030:1 6030-20000:0`.
  - `pair_rule_mates_bam`: `mates.bam`'s chr1 runs equal mosdepth's except that the contained mate (template c) leaves `400-600:1 600-700:0`, where mosdepth has `400-500:1 500-700:0`. The full expected chr1 list is `0-100:0 100-250:1 250-300:0 300-380:1 380-400:0 400-600:1 600-700:0 700-750:1 750-900:0 900-950:1 950-1000:0 1000-1040:1 1040-1100:0 1100-1130:1 1130-1140:0 1140-1145:1 1145-1200:0 1200-1250:1 1250-1400:0 1400-1500:1 1500-1519:2 1519-1530:3 1530-1539:2 1539-1550:1 1550-1600:0 1600-1670:1 1670-1800:0 1800-1840:1 1840-1900:0 1900-1990:1 1990-2000:0`, and chr2 is `0-4:0 4-44:1 44-100:0 100-170:1 170-300:0 300-350:1 350-500:0`.
  - `stale_pair_entry_is_dropped`: `mosdepth-stale.bam` with `exclude_flags: 1540` gives `0-100:0 100-130:1 130-140:2 140-150:1 150-300:0 300-350:1 350-1000:0`. mosdepth has `140-350:1`.
  - `pair_table_bounded`: an in-memory SAM on one reference of LN 30,000,000.
    - 100,000 proper pairs (`99`/`147`, `150M`): pair i at 1-based positions `200·i+1` and `200·i+51`. `DepthSummary::max_pending_mates` is at most 2.
    - Then 1,000 more pairs on the same pattern whose second mate is a duplicate (`1171`), filtered by the plan: the peak stays at most 2.
- [ ] **Step 3: Run** `cargo test -p gloci-ops --lib depth::pairs`. Expected: FAIL to compile.
- [ ] **Step 4: Implement** the interfaces above.
- [ ] **Step 5: Run** `cargo test -p gloci-ops --lib depth`, `cargo test -p gloci-ops --test depth_mates --test depth_core`, `cargo test -p gloci --test depth --test depth_summary`. Expected: PASS.
- [ ] **Step 6: Commit** `feat(depth): mosdepth's mate rule and read filter in the depth core`

### Task 4: gloci-io stages a BGZF output with a CSI index

**Files:**
- Modify: `crates/gloci-io/src/index/build.rs` (`build_csi_index`), `crates/gloci-io/src/output.rs` (`Output::finish_staged_with_csi`), `crates/gloci-io/src/staged.rs` (`StagedSet::remove_on_publish`)
- Test: `crates/gloci-io/tests/output.rs`

**Interfaces:**
- Produces (Task 5 uses these):

```rust
// index/build.rs
/// The CSI index of the finished BGZF file `data` (min_shift 14, the smallest depth of at least 5
/// that covers its largest end, the tabix meta of `format`), whatever its coordinates; messages
/// name `final_data_path`. BCF's branch of `build_index_with`, for any text format.
pub(crate) fn build_csi_index(data: &Path, final_data_path: &Path, format: Format,
    threads: NonZeroUsize) -> Result<Vec<u8>, IoError>;

// output.rs
impl Output {
    /// Finish a BGZF file output and stage it, and its CSI index as `PATH.csi`, in `set`; both are
    /// published by `set.publish()`, which also removes a `PATH.tbi` left from earlier data. An
    /// output that is not a regular BGZF file is an internal error. An unsorted file is the error
    /// `unindexable_output` gives, and nothing is staged.
    pub fn finish_staged_with_csi(self, format: Format, set: &mut StagedSet) -> Result<(), IoError>;
}

// staged.rs
impl StagedSet {
    /// Remove `path`, if it exists, once every staged file is in place.
    pub fn remove_on_publish(&mut self, path: PathBuf);
}
```

- [ ] **Step 1: Write the failing tests** in `crates/gloci-io/tests/output.rs`. Each writes `#chrom\tstart\tend\tregion\t1X\nchr1\t0\t100\tunknown\t5\nchr1\t100\t200\tunknown\t0\nchr2\t0\t50\tunknown\t1\n` through `open_output(&OutputSpec { path: Some(dir/"o.thresholds.bed.gz"), format: Some(OutFormat::Bed), bgzf: true, direct: false }, OutFormat::Bed, &opts(t, false))`.
  - `staged_bgzf_gets_a_csi`:
    - after `finish_staged_with_csi(Format::Bed, &mut set)`, neither file exists;
    - after `set.publish()`, both exist;
    - the decompressed `.csi` starts with `b"CSI\x01"`, then `14u32` (min_shift), a depth `>= 5` and, at bytes 16..20, `0x10000` (the BED preset), followed by columns 1, 2, 3 and meta `#`;
    - `open_region_reader(dest, &[(b"chr1".to_vec(), Interval { start: 0, end: 100 })], &ro, false)` yields only the `chr1\t0\t100` line.
  - `staged_csi_removes_a_stale_tbi`: an existing `o.thresholds.bed.gz.tbi` is gone after `publish()`.
  - `staged_csi_unsorted_stages_nothing`: lines `chr1 100 200` then `chr1 0 100` → `finish_staged_with_csi` errs, and the message names `o.thresholds.bed.gz`. Publishing the (empty) set creates no file, and no temp file remains (`no_tmp_files`).
  - `staged_csi_is_thread_independent`: at threads 1 and 4, the data bytes and the `.csi` bytes are equal.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test output`. Expected: FAIL to compile.
- [ ] **Step 3: Implement.**
  - `build_csi_index`: a `pass::<BinnedIndex>` at `TBI_DEPTH`, redone at `csi_depth(max_end)` when it overflows, then `with_htslib_loffsets` and `csi_bytes`, as the BCF branch does.
  - `finish_staged_with_csi`: finish the `FileBgzf` sink, build from `af.temp_path()`, `stage_index(dest + ".csi", &bytes)`, push both onto the set, and `remove_on_publish(dest + ".tbi")`.
- [ ] **Step 4: Run** `cargo test -p gloci-io --test output --test index`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(io): stage a BGZF output with a CSI index`

### Task 5: mosdepth's tables and files

**Files:**
- Create: `crates/gloci-ops/src/depth/mosdepth.rs`
- Modify: `crates/gloci-ops/src/depth/mod.rs` (`pub mod mosdepth;` and re-exports)
- Fixtures: `xtask/diff/inputs/samtools/mosdepth.bed` and `mosdepth-blank.bed`, hand-written (exact bytes below), with a comment line appended to `make-samtools.sh` saying so
- Test: `#[cfg(test)] mod tests` in `mosdepth.rs`

**Interfaces:**
- Consumes: `run_sources`, `ModeRules { shared_bases: true, filter }`, `DepthSink` (Tasks 1, 3); `Output::finish_staged_with_csi`, `StagedSet` (Task 4); gloci's BED reader (`open_reader` with `Format::Bed`, `RecordBatch::raw`, `RecordMeta::loc`, `gloci_io::field`), as `view::for_each_region` uses it.
- Produces (Task 6 uses these):

```rust
pub enum By { Nothing, Windows(u64), Bed(String) }
pub struct MosdepthOptions {
    /// Output path stem: files are `<prefix>.mosdepth.global.dist.txt` and so on.
    pub prefix: PathBuf,
    pub by: By,
    /// Ascending, duplicates kept, each at least 1 (the caller checks).
    pub thresholds: Vec<u32>,
    pub exclude_flags: u16,
    pub min_mapq: u8,
}
pub struct MosdepthSummary { pub max_ring: usize, pub max_pending_pairs: usize, pub regions: u64, pub files: Vec<PathBuf> }
/// mosdepth's files for `src` (opened from `path` with no reader filters), as D173 part B says.
/// Every file is staged and published together; on an error none is.
pub fn mosdepth_files(src: Box<dyn AlnSource>, path: &str, o: &MosdepthOptions,
    ro: &ReaderOptions, oo: &OutputOptions) -> Result<MosdepthSummary, IoError>;

// crate-internal pieces with fixed behaviour:
pub(crate) fn write_distribution(name: &[u8], hist: &[u64], out: &mut Vec<u8>);
pub(crate) struct DepthStat { cum_length: u64, cum_depth: u64, min: u32 /* u32::MAX = none */, max: u32 }
// DepthStat: add(v, n), merge(&other), write_line(name, out) -> "name\tlen\tbases\t%.2f\tmin\tmax\n"
pub(crate) fn hist_slot(v: u32) -> usize;   // v > 400_000 -> 399_990, else v
pub(crate) struct BedRegion { pub start: u64, pub end: u64, pub name: Option<Box<[u8]>> }
pub(crate) fn load_bed_regions(bed: &str, hdr: &AlnHeader, input: &str, ro: &ReaderOptions)
    -> Result<Vec<Vec<BedRegion>>, IoError>;   // per reference id, stable-sorted by start
```

`write_distribution` is mosdepth's algorithm; the accumulation order fixes the rounding:

```rust
let sum: u64 = hist.iter().sum();
if sum < 1 { return; }
let mut cum = 0f64;
for d in (0..hist.len()).rev() {
    let v = hist[d];
    if d > 300 && v == 0 { continue; }
    cum += v as f64 / sum as f64;
    if cum < 8e-5 { continue; }
    // name \t d \t cfmt::write_fixed(cum, 2) \n
}
```

The sink, `Tables`, is built over five `&mut dyn Write` targets: global dist, summary, and the optional region dist, regions and thresholds. Unit tests drive it with `Vec<u8>`s; `mosdepth_files` builds it over the staged outputs. `mosdepth_files` runs the core with `AllPositions::AllRefs`, so every position of every reference arrives exactly once.
- Per reference it keeps:
  - "placed" (from the `placed` hook);
  - a per-base histogram and a `DepthStat`;
  - for the current reference's windows or BED regions: a mean accumulator, threshold counts and a `DepthStat`;
  - in BED mode, a per-base region histogram; in window mode, a per-window histogram of `min(mean.round() as usize, 511)` (`f64::round` rounds half away from zero, as Nim's `toInt`).
- Means: for a run of n bases at depth v, add `v as f64 / len as f64` n times; skip runs at depth 0, where the addition changes nothing.
- At `reference_done(rid)`:
  - write the reference's region rows (`chrom\tstart\tend\t[name\t]mean`) and threshold rows (`chrom\tstart\tend\tname-or-unknown\tc1…`);
  - when it was placed: write its dist lines and summary lines, and merge its histograms and stats into the totals;
  - drop the reference's accumulators.
- At the end:
  - write the `total` dist lines, then the summary `total` and `total_region` lines;
  - finish the text files with `finish_staged` and the `.bed.gz` files with `finish_staged_with_csi(Format::Bed, …)`;
  - `publish()`.
- The summary header line is written before the first summary line.
- Text files open through `plain_text_output_spec` + `open_output(OutFormat::Tsv)`. BGZF files open through `OutputSpec { format: Some(OutFormat::Bed), bgzf: true, … }`.

**`load_bed_regions` rules:**
- Lines are read through gloci's BED reader, which skips comment, `track`, `browser` and blank lines.
- A non-empty fourth field is the name.
- A chromosome the header lacks gets one warning (`{bed}: N region(s) on chrZ skipped: {input} has no such reference`).
- A region ending past `LN` is `IoError::Parse` naming the BED file and line, e.g. `line 4: chr1 480 1001 ends past the end of chr1 (LN 1000)`.

**Fixture bytes:**
- `mosdepth.bed`: `# targets for the mosdepth cases\ntrack name=targets\nbrowser position chr1:1-100\nchr4\t5000\t5001\tdeep\t0\t+\nchr1\t700\t800\tv\nchr1\t100\t180\tp\textra\nchr1\t500\t550\nchr1\t200\t260\tq\nchr1\t200\t210\tq-head\nchr1\t650\t650\tpoint\nchrZ\t1\t10\tnowhere\nchr3\t0\t200\tempty3\nchr2\t0\t300\tdups2\r\n`
- `mosdepth-blank.bed`: `chr1\t100\t180\ta\n\nchr1\t200\t260\tb\n`

- [ ] **Step 1: Write the failing tests** (unit tests that feed the sink positions directly; `refs` and writers are in-memory `Vec<u8>`):
  - `distribution_lines_follow_mosdepth`:
    - hist `{0: 19969, 1: 20, 2: 10, 320: 1}` named `chr4` → exactly `chr4\t2\t0.00\nchr4\t1\t0.00\nchr4\t0\t1.00\n`;
    - hist `{320: 1}` → 302 lines: `chr4\t320\t1.00`, then `chr4\t300\t1.00` down to `chr4\t0\t1.00`;
    - an all-zero or empty hist → no line.
  - `depth_stat_lines`: depths `[0, 0, 3]` → `x\t3\t3\t1.00\t0\t3\n`; no bases → `x_region\t0\t0\t0.00\t0\t0\n`.
  - `hist_caps_at_400000`: `hist_slot(400_001) == 399_990` and `hist_slot(400_000) == 400_000`.
  - `region_mean_accumulates_base_by_base`:
    - setup: one reference `c` LN 400, `By::Windows(200)`, depth 1 at `[0,11)` and `[200,225)`, 0 elsewhere;
    - the regions text is `c\t0\t200\t0.05\nc\t200\t400\t0.13\n` (the exact quotients would print `0.06` and `0.12`);
    - mosdepth 0.3.10 printed these two lines for this input.
  - `window_dist_rounds_half_away_from_zero`: windows of 100. A window with 50 bases at depth 1 counts in bin 1, one with 49 bases in bin 0, and one with mean 600 in bin 511.
  - `reference_without_records_is_left_out`: refs `c1` (placed, all depth 0) and `c2` (never placed), windows of 10.
    - global dist = `c1\t0\t1.00\ntotal\t0\t1.00\n`;
    - the summary has `c1` and `c1_region` lines and no `c2` line;
    - the regions text has `c2` rows with `0.00`.
  - `thresholds_header_and_rows`: thresholds `[1, 1, 2]` → header `#chrom\tstart\tend\tregion\t1X\t1X\t2X\n`; a window row ends `unknown\t…`.
  - `bed_regions_sorted_stably_and_named`: `mosdepth.bed` against `mosdepth.bam`'s header.
    - chr1 = `(100,180,p) (200,260,q) (200,210,q-head) (500,550,-) (650,650,point) (700,800,v)`;
    - chr2 = `(0,300,dups2)`; chr3 = `(0,200,empty3)`; chr4 = `(5000,5001,deep)`;
    - one warning names `chrZ`.
  - `bed_region_past_the_end_is_an_error`: `chr1\t480\t1001\n` with chr1 LN 1000 → `Err`, and the message contains `line 1` and `LN 1000`.
  - `blank_bed_lines_are_skipped`: `mosdepth-blank.bed` → regions `a` and `b`.
- [ ] **Step 2: Run** `cargo test -p gloci-ops --lib depth::mosdepth`. Expected: FAIL to compile.
- [ ] **Step 3: Implement** `mosdepth.rs`.
- [ ] **Step 4: Run** `cargo test -p gloci-ops --lib depth`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(depth): mosdepth's region, threshold, distribution and summary tables`

### Task 6: `depth --out-prefix/--by/--thresholds`, mosdepth in the comparison image, parity

**Files:**
- Modify:
  - `xtask/docker/Dockerfile`: add one line, `        mosdepth=0.3.10 \`, after the line that holds `samtools=1.21`;
  - `crates/gloci/src/cli/alignments.rs`: `DepthArgs` and help;
  - `crates/gloci/src/commands/depth.rs`: `run_files`;
  - `crates/gloci/tests/cli_basics.rs`: add `"mosdepth"` to `TOOLS` in `long_summary_names_every_tool_of_the_short_one`.
- Fixture: `xtask/diff/inputs/samtools/mosdepth-means.{sam,bam,bam.bai}` (block below)
- Create: `crates/gloci/tests/mosdepth/regen.sh`, the committed outputs under `crates/gloci/tests/mosdepth/<case>/`, and `crates/gloci/tests/depth_files.rs`
- Modify: `xtask/diff/cases.toml` (live cases plus one golden), `docs/differences.md` (six sections), `docs/decisions.md` (D173 part B), `docs/migration.md` (a `### mosdepth` subsection after `### samtools`), `SKILL.md` and the pystub regenerated

**Interfaces:**
- Consumes: `mosdepth_files`, `MosdepthOptions`, `By` (Task 5).
- Produces:
  - `DepthArgs`:
    - `pub out_prefix: Option<PathBuf>` (`--out-prefix P`), with `conflicts_with_all` = summary, output, all_positions, all_references, count_deletions, skip_mate_overlap, column_names, region, regions_file, min_baseq, min_read_length, require_flags;
    - `pub by: Option<String>` (`--by BED|SIZE`, `requires = "out_prefix"`);
    - `pub thresholds: Vec<u32>` (`--thresholds N,…`, `requires = "by"`, `value_delimiter = ','`).
  - `run_files` usage errors (exit 2): more than one FILE (`--out-prefix writes the files of one FILE, as mosdepth`); `--by 0`; a threshold of 0.
  - It opens the input with `AlignOptions { exclude_flags: 0, require_flags: 0, min_mapq: 0, … }`, passes `--exclude-flags` and `--min-mapq` in `MosdepthOptions`, and sorts the thresholds.
- Help:
  - summary: `Depth of alignments (samtools depth, coverage, mosdepth)` (56 characters);
  - the long help's first line adds "…, or as mosdepth's coverage files (mosdepth).";
  - Examples: replace `gloci depth a.bam b.bam --column-names …` with `gloci depth sample.bam --by targets.bed --thresholds 1,10,20 --out-prefix sample`;
  - a new section "Coverage files (--out-prefix P):" that lists the files and says:
    - one FILE;
    - `--by BED` keeps a fourth column as the name; `--by SIZE` makes windows;
    - the defaults are mosdepth's, the same `--exclude-flags` and `--min-mapq` as here;
    - the bases both mates of a proper pair cover count once, by mosdepth's rule, without `--skip-mate-overlap`;
    - memory follows the longest read and the BED, not the reference; `-t` never changes the files.

**Fixture block** (append to `make-samtools.sh`):

```sh
# mosdepth-means.sam: windows of 200 whose exact means end in 5 at the third decimal (11/200, 25/200);
# mosdepth's base-by-base sum prints 0.05 and 0.13.
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:400\n'
  printf 'r1\t0\tc\t1\t60\t11M\t*\t0\t0\t*\t*\nr2\t0\tc\t201\t60\t25M\t*\t0\t0\t*\t*\n'
} > mosdepth-means.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o mosdepth-means.bam mosdepth-means.sam
samtools index -@ 0 mosdepth-means.bam
```

**Committed mosdepth outputs.** One directory per case under `crates/gloci/tests/mosdepth/`. Each holds mosdepth's files without the prefix: `mosdepth.global.dist.txt`, `mosdepth.summary.txt`, and where written `mosdepth.region.dist.txt`, `regions.bed` (`zcat` of `regions.bed.gz`) and `thresholds.bed`. The case table is shared by `regen.sh` and `depth_files.rs`:

| case | input | mosdepth (`-n` always) | gloci (`--out-prefix P` always) | anchor |
|---|---|---|---|---|
| `depth-window-100` | depth.bam | `--by 100` | `--by 100` | |
| `depth-targets-thresholds` | depth.bam | `--by targets.bed -T 1,2,10` | `--by targets.bed --thresholds 1,2,10` | |
| `pairs-window-thresholds` | mosdepth.bam | `--by 100 -T 1,2,3` | `--by 100 --thresholds 1,2,3` | |
| `pairs-bed-mapq` | mosdepth.bam | `-Q 20 --by mosdepth.bed -T 2,1` | `--min-mapq 20 --by mosdepth.bed --thresholds 2,1` | |
| `pairs-no-by` | mosdepth.bam | (none) | (none) | |
| `pairs-flags-zero` | mosdepth.bam | `-F 0 --by 500` | `--exclude-flags 0 --by 500` | |
| `means-base-by-base` | mosdepth-means.bam | `--by 200` | `--by 200` | |
| `empty` | empty.bam | `--by 100` | `--by 100` | |
| `mates-contained` | mates.bam | `--by 100` | `--by 100` | `#mosdepth-contained-mate` |
| `stale-pair-entry` | mosdepth-stale.bam | `-F 1540 --by 100` | `--exclude-flags UNMAP,QCFAIL,DUP --by 100` | `#mosdepth-mate-table` |
| `bed-blank-line` | mosdepth.bam | `--by mosdepth-blank.bed` | `--by mosdepth-blank.bed` | `#mosdepth-bed-blank-line` |

Two extra gloci-only rows compare with `depth-window-100`'s files: `depth.cram` (`--by 100`), and `depth.sam` on stdin (`-` `--by 100`). For an anchored case, `regen.sh` writes mosdepth's files, and the test writes gloci's own pinned files to `<case>/gloci/` (generated once by the test under `GLOCI_BLESS=1`, then reviewed and committed).

`regen.sh` (run from its directory) is a `docker run --rm -u "$(id -u):$(id -g)"` in `gloci-diff:p0`. It mounts this directory at `/work` and `../../../../xtask/diff/inputs/samtools` read-only at `/in`. For each row it runs `mosdepth -n <args> <case>/o /in/<input>` (BED paths as `/in/<file>`), renames `o.mosdepth.*` to `mosdepth.*`, `zcat`s `o.*.bed.gz` into `*.bed`, and removes the `.gz` and `.csi` files.

- [ ] **Step 1: The comparison image gains mosdepth.**
  1. Before the change, save `docker run --rm gloci-diff:p0 micromamba list > $SCRATCH/before.txt`.
  2. Add the Dockerfile line, then `docker build -t gloci-diff:p0 xtask/docker`.
  3. Save the same list to `$SCRATCH/after.txt` and `diff` the two. Expected: only `mosdepth` and its new dependencies are added, and every pinned tool (bedtools, bedops, vcftools, htslib, bcftools, samtools, hyperfine) keeps its version and build.
  4. Check `docker run --rm gloci-diff:p0 mosdepth --version`. Expected: `mosdepth 0.3.10`.
  5. Spot-check `cargo xtask diff --only sort-bed-vs-gnu-sort`, `--only depth-default-vs-samtools` and `--only vcf-stats-freq-vs-vcftools`. Expected: `ok`.
  6. `docker image prune -f`.
- [ ] **Step 2: Build the fixture** (as in Task 2, Step 1), then run `sh regen.sh` and commit the outputs. Expected values, from mosdepth 0.3.10:
  - `pairs-window-thresholds/mosdepth.summary.txt` = `chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t1000\t600\t0.60\t0\t2\nchr1_region\t1000\t600\t0.60\t0\t2\nchr2\t300\t0\t0.00\t0\t0\nchr2_region\t300\t0\t0.00\t0\t0\nchr4\t20000\t360\t0.02\t0\t320\nchr4_region\t20000\t360\t0.02\t0\t320\ntotal\t21300\t960\t0.05\t0\t320\ntotal_region\t21300\t960\t0.05\t0\t320\n`;
  - `means-base-by-base/regions.bed` = `c\t0\t200\t0.05\nc\t200\t400\t0.13\n`;
  - `pairs-bed-mapq/regions.bed` starts `chr1\t100\t180\tp\t1.00\nchr1\t200\t260\tq\t1.00\nchr1\t200\t210\tq-head\t1.00\nchr1\t500\t550\t1.00\n`;
  - `empty/mosdepth.summary.txt` = `chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t0\t0\t0.00\t0\t0\ntotal_region\t0\t0\t0.00\t0\t0\n`, and its two dist files are empty.
- [ ] **Step 3: Write the failing tests** in `crates/gloci/tests/depth_files.rs`:
  - `files_match_mosdepth`, for every row of the table (plus the two gloci-only rows):
    - run `gloci depth <input> <gloci args> --out-prefix <tmp>/o`;
    - the set of files written equals the case directory's set plus `.gz` and `.csi`;
    - each text file is byte-equal; each `.bed.gz` is BGZF, and its decompressed text is byte-equal.
  - For an anchored row: gloci's files differ from mosdepth's and equal `<case>/gloci/`, and `check_anchor` finds the anchor in `docs/differences.md`. The expected differing lines are:
    - `mates-contained`: regions `chr1\t500\t600\t1.00` (mosdepth `0.00`); summary `chr1\t2000\t1055\t0.53\t0\t3` (mosdepth `955 0.48`) and `total\t2500\t1215\t0.49\t0\t3`;
    - `stale-pair-entry`: regions `chr1\t100\t200\t0.60` and `chr1\t200\t300\t0.00` (mosdepth `1.10`, `1.00`); summary `chr1\t1000\t110\t0.11\t0\t2` (mosdepth `260 0.26`);
    - `bed-blank-line`: regions `chr1\t100\t180\ta\t1.37\nchr1\t200\t260\tb\t1.00\n` (mosdepth stops after `a`).
  - `files_have_a_csi_with_the_bed_preset`: `depth.bam --by 100` → `o.regions.bed.gz.csi` exists, and `gloci view o.regions.bed.gz --region chr1:1-100` prints exactly `chr1\t0\t100\t0.81\n`. With mosdepth's `.csi`, tabix also prints `chr1\t100\t200\t0.03`.
  - `past_reference_end_is_cut_with_one_warning`: `depth-pastend.bam --by 100` exits 0; stderr has exactly one `past the end` warning; regions contain `chr2\t0\t100\t0.20` and `chr2\t200\t300\t0.12`.
  - `usage_errors`, each exit 2 with no `o.*` file: two FILEs; `--by 0`; `--by 100 --thresholds 0`; `--thresholds 1` without `--by`; `--by 100` without `--out-prefix`; `--out-prefix` with `-o`, `--region`, `--min-baseq 1`, `--min-read-length 1`, `--require-flags 1`, `--summary`, `--skip-mate-overlap` or `--all-positions`.
  - `bed_region_past_end_publishes_nothing`:
    - first a good run at prefix `o`;
    - then `--by bad.bed`, where `bad.bed` is `chr1\t480\t501\n` against depth.bam (chr1 LN 500); it exits 1, and stderr contains `bad.bed`, `line 1` and `LN 500`;
    - every `o.*` file from the good run is byte-identical afterwards, and no temp file remains.
  - `unwritable_prefix_fails_cleanly`: `--out-prefix <tmp>/missing-dir/o --by 100` exits 1, stderr names `o.mosdepth.global.dist.txt`, and `<tmp>` holds nothing new.
  - `threads_do_not_change_the_files`: `mosdepth.bam --by mosdepth.bed --thresholds 1,2` at `-t 1` and `-t 4` → every file byte-identical, `.bed.gz` and `.csi` included.
  - `by_windows_huge_header`: a temp SAM with 200,000 `@SQ` (LN 10) and one read on the last reference, `--by 5` → `o.regions.bed.gz` holds 400,000 rows; the run takes under 30 s in a debug build.
- [ ] **Step 4: Run** `cargo test -p gloci --test depth_files`. Expected: FAIL (unknown argument `--out-prefix`).
- [ ] **Step 5: Implement** the CLI, the help and `run_files`. Then run `GLOCI_BLESS=1 cargo test -p gloci --test depth_files files_match_mosdepth`, review the three `gloci/` directories against the lines above, and run without `GLOCI_BLESS`. Expected: PASS.
- [ ] **Step 6: Live cases and one golden.** Append a `# depth --out-prefix (D173) vs mosdepth 0.3.10` section to `xtask/diff/cases.toml`. The gloci side is `gloci depth input.bam <args> --out-prefix o && cat o.mosdepth.global.dist.txt o.mosdepth.region.dist.txt o.mosdepth.summary.txt && zcat o.regions.bed.gz`. The reference side is `samtools index input.bam && mosdepth -n <args> o input.bam && cat … && zcat o.regions.bed.gz`. BED files are `aux`.
  - Exact live cases: `depth-files-window-vs-mosdepth` (depth.bam, `--by 100`); `depth-files-targets-thresholds-vs-mosdepth` (aux targets.bed, adds `&& zcat o.thresholds.bed.gz` on both sides); `depth-files-pairs-bed-vs-mosdepth` (mosdepth.bam, aux mosdepth.bed, `-Q 20` / `--min-mapq 20`).
  - Pinned live cases:
    - `depth-files-contained-mate-vs-mosdepth` (mates.bam): `#mosdepth-contained-mate`;
    - `depth-files-past-end-vs-mosdepth` (depth-pastend.bam, `compare = "exit-status"`): `#mosdepth-past-reference-end`;
    - `depth-files-csi-preset-vs-mosdepth` (`… && tabix o.regions.bed.gz chr1:1-100` on both sides): `#mosdepth-csi-preset`.
  - The golden: `depth-files-threshold-zero-vs-mosdepth`, with `gloci = "gloci depth input.bam --by 100 --thresholds 0 --out-prefix o"`. This is a plain invocation that exits 2 before writing, so it is safe as a golden. `reference = "samtools index input.bam && mosdepth -n --by 100 -T 0 o input.bam"`, `compare = "exit-status"`, `expected_difference = "#mosdepth-thresholds-zero"`.
  - Then run `cargo xtask diff --only depth-files-threshold-zero-vs-mosdepth --update-golden`, and `cargo xtask diff --only NAME` for each new case. Expected: `ok` for every case; the pinned ones say `differs as documented`.
- [ ] **Step 7: Docs.**
  - D173 part B.
  - Six `docs/differences.md` sections, each with gloci / old behaviour / reason / pinned by (the test names above and the live case):
    - `#mosdepth-contained-mate`;
    - `#mosdepth-mate-table`;
    - `#mosdepth-past-reference-end`;
    - `#mosdepth-bed-blank-line`;
    - `#mosdepth-csi-preset`;
    - `#mosdepth-thresholds-zero`.
  - `docs/migration.md`: a new `### mosdepth` subsection with a table:
    - `mosdepth -n [--by BED|N] [-T N,...] [-Q N] [-F N] [-t N] P in.bam` → `gloci depth in.bam --out-prefix P [--by BED|N] [--thresholds N,...] [--min-mapq N] [--exclude-flags F] [-t N]`, noting that `-F 1796` is the default and naming the anchors;
    - mosdepth's per-base file → `gloci depth in.bam --all-references` (one row per position), or runs with `gloci genomecov in.bam --bga --split-deletions --exclude-flags UNMAP,SECONDARY,QCFAIL,DUP` (without mosdepth's mate correction);
    - not provided: `--quantize`, `--use-median`, `--fast-mode`, `--chrom`, `--read-groups`, `--include-flag`, `--min-frag-len`, `--max-frag-len`, `MOSDEPTH_PRECISION`.
  - `cargo xtask pystub`, `cargo xtask docs --skill SKILL.md`.
- [ ] **Step 8: Run** `cargo test -p gloci --test depth_files --test depth_summary --test depth --test golden --test cli_basics` and `cargo test -p xtask --test pystub --test docs --test references --test book`. Expected: PASS.
- [ ] **Step 9: Commit** `feat(depth): --out-prefix writes mosdepth's files; mosdepth joins the comparison image (D173)`

### Task 7: Python, design, changelog, benchmarks, large tests

**Files:**
- Modify:
  - `crates/gloci/src/embed/bindings.rs`: `nary(&["depth"], "depth").shaped(Shape::Table, &[("summary", Shape::HeaderLine)])`, and `(&["depth"], "depth_files")` in `MODULE_FUNCTIONS`;
  - `crates/gloci/src/embed/args.rs`: three `RESERVED_ON_STEP_OF` entries for `depth`: `out_prefix`, `by` and `thresholds`, each with the message `{kw} is set by depth_files(src, out=...): mosdepth's files are written to disk, not to a pipeline`;
  - `crates/gloci/src/embed/signature.rs`: `"depth_files"` joins the `vcf_stats` arm of the source help;
  - `crates/gloci-py/python/gloci/_functions.py` (`depth_files`, as `vcf_stats`) and `crates/gloci-py/python/gloci/__init__.py` (export and `__all__`).
- Test:
  - `crates/gloci/tests/embed_args.rs` (`the_depth_step_leaves_the_files_to_its_function`);
  - `crates/gloci-py/tests/test_api_surface.py` (add `"depth_files"` to the function set);
  - `crates/gloci-py/tests/test_sources.py` (`test_depth_files_match_the_cli`);
  - `crates/gloci-py/tests/typing_example.py` (one `gloci.depth_files("a.bam", out="d", by="500")` line).
- Modify:
  - `docs/design.md`: in §7.4, a bullet after the `depth, bedcov` bullet; in §11, a samtools-table row for `coverage` and a row for mosdepth; in §9.2, a benchmark-table row "`depth --summary`, `depth --by` | samtools coverage; mosdepth (binding)";
  - `docs/decisions.md`: D173 part C;
  - `CHANGELOG.md`: one bullet appended to `[Unreleased]` → `### Added`;
  - `xtask/bench/cases.toml` and `crates/gloci/tests/large.rs`: below.

**Interfaces:**
- Consumes: the CLI of Tasks 2 and 6.
- Produces: `gloci.depth_files(src, *, out, **kwargs) -> list[pathlib.Path]`.

- [ ] **Step 1: Write the failing tests.**
  - `the_depth_step_leaves_the_files_to_its_function`:
    - `argv_of(&["depth"], &[("out_prefix", val("x"))])` errs with exactly `out_prefix is set by depth_files(src, out=...): mosdepth's files are written to disk, not to a pipeline`;
    - `kwargs_to_argv_as(spec(&["depth"]), &kw(&[("by", val("100"))]), Role::Function)` gives `["--by", "100"]`;
    - `argv_of(&["depth"], &[("summary", KwValue::Bool(true))])` is `Ok`.
  - `test_depth_files_match_the_cli`:
    - write the SAM text of the `mosdepth-means.sam` block (Task 6) to `here/"in"/"m.sam"`;
    - `want = cli_files(cli, here, ["depth", "--out-prefix", "m", "--by", "200", sam], "c9")`;
    - `got = gloci.depth_files(sam, out=here/"py"/"m", by="200")`;
    - `[p.name for p in got] == ["m.mosdepth.global.dist.txt", "m.mosdepth.region.dist.txt", "m.mosdepth.summary.txt", "m.regions.bed.gz", "m.regions.bed.gz.csi"]` and `written(got, here/"py") == want`;
    - `gloci.open(sam).depth(summary=True).to_bytes()` equals the CLI's `depth --summary` stdout.
- [ ] **Step 2: Run** `cargo test -p gloci --test embed_args the_depth_step`. Expected: FAIL.
- [ ] **Step 3: Implement** the bindings, the reservations and the Python function. Run `cargo xtask pystub`.
- [ ] **Step 4: Run**:
  - `cargo test -p gloci --test embed_args`;
  - `maturin develop -m crates/gloci-py/Cargo.toml && pytest crates/gloci-py/tests/test_api_surface.py crates/gloci-py/tests/test_sources.py -k "depth or runtime_functions"`;
  - `cargo test -p xtask --test pystub`.
  
  Expected: PASS.
- [ ] **Step 5: Benchmarks and large tests.**
  - In `xtask/bench/cases.toml`, after `depth-bam-50M-t8`, add two cases. Both use the `sam --records 50000000 --read-len 100 --depth 30 --chroms 24 --spliced-frac 0.1 --seed 1` setup.
    - `depth-summary-bam-50M`, input `bam`, target `gloci <= samtools`: `gloci = "{gloci} depth -t 1 --summary {in} > {out}"`, `samtools = "samtools coverage {in} > {out}"`, with `same` = both.
    - `depth-by-bam-50M`, input `bam-indexed`, target `gloci-t4 <= mosdepth-t4`:
      - `gloci-t4 = "{gloci} depth -t 4 --by 500 --thresholds 1,10,20 --out-prefix {out}.g {in} && cat {out}.g.mosdepth.summary.txt {out}.g.mosdepth.global.dist.txt > {out}"`;
      - `mosdepth-t4 = "mosdepth -t 4 -n --by 500 -T 1,10,20 {out}.m {in} && cat {out}.m.mosdepth.summary.txt {out}.m.mosdepth.global.dist.txt > {out}"`;
      - `same` = both.
  - In `crates/gloci/tests/large.rs`, add two `#[ignore = "nightly: generates a 2 GB BAM"]` tests that call `p3c_rss_ceiling`:
    - `depth_summary_2gb_bam_streams`, on `big_sorted_bam()` with `--summary`;
    - `depth_by_2gb_pairs_streams`, on `big_pairs_bam()` with `--by 500 --thresholds 1,10,20 --out-prefix <tmp>`, and with `--by big_targets()`.
  - Run `cargo test -p xtask` and `cargo test -p gloci --test large --no-run`. Expected: PASS and compiles. Do not run the ignored tests.
  - The README's performance table gains a `depth --by` row only from the first nightly bench run that has these cases. That run is not part of this task.
- [ ] **Step 6: Docs.**
  - D173 part C.
  - `docs/design.md` §7.4 bullet: "**`depth --summary`, `depth --out-prefix`** **[D173, owner 2026-11-17]**: samtools coverage's table, and mosdepth 0.3.10's files under their names; both stream on the depth core, and the files are published together with CSI indexes."
  - `CHANGELOG.md` bullet: "`gloci depth --summary` prints samtools coverage's per-reference table, and `gloci depth --out-prefix P --by BED|SIZE [--thresholds 1,10,20]` writes mosdepth's global and region distributions, summary, `regions.bed.gz` and `thresholds.bed.gz` (with `.csi`) under mosdepth's names, so MultiQC reads them unchanged."
  - `cargo xtask docs --skill SKILL.md`.
- [ ] **Step 7: Run** `cargo test --workspace`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo fmt --all -- --check` and `cargo test -p xtask --test pystub --test docs --test references --test book`. Expected: PASS.
- [ ] **Step 8: Commit** `feat(py): depth_files() and depth(summary=True); depth summary docs and benchmarks`
