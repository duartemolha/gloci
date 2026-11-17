# bigWig and bigBed Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Read bigWig and bigBed wherever gloci reads intervals (whole files in the command's chromosome order, regions through the file's own index); write them from `gloci convert` and `gloci genomecov --bedgraph/--bga`; add `map`'s `wsum`; check all of it against pyBigWig; and take the last UCSC binary out of the benchmark image.

**Architecture:** A new module `gloci-io/src/bbi/` wraps the `bigtools` crate.
- **Reading** is a `RecordReader` (`BbiReader`) that yields bedGraph rows (bigWig) or BED rows (bigBed) as text, as the WIG reader does (D63). It fetches one chromosome at a time through the file's R-tree, in the order of the `ChromTable` its `read_batch` is handed; a region read is the same reader with a region list.
- **Writing** is a `BbiWriter`: it parses each bedGraph/BED line a command produces, validates it, and sends it through a bounded channel to a thread that runs bigtools' writer on a tokio runtime, into an `AtomicFile`. It is also a `RecordWriter`, so `genomecov` writes bigWig without touching its core.
- Only `convert` and `genomecov` accept a bigWig/bigBed output (`parse_track_output_spec`); every other command's `parse_output_spec` refuses one with a message naming those two.

**Tech Stack:** Rust 1.100.0 (MSRV 1.91); `bigtools` 0.5.8 (MIT; features `read` and `write` only); `tokio` 1 (`rt`, `rt-multi-thread`); pyBigWig 0.3.29 (MIT, on libBigWig) in the comparison image, for checking only.

**Spec:** `docs/superpowers/specs/2026-11-17-tracks-query-coverage-liftover-design.md`: "Feature 1: bigWig and bigBed", bound by "Owner constraints", "Common to all four" and "Order of the work".

**Builds on:** D9/D10 (region reads), D17 (output format), D62 (`--chrom-map`), D63 (WIG as bedGraph rows), D87 (Python placement), D91 (cancellation), D100 (input gate). This plan adds **D171**.

**Lane A** of three parallel lanes (B `vcf query`, C `depth` summaries then liftover). Task 2 can run beside Task 1; the rest run in order. Tasks 2 and 9 are image, script and benchmark-data work: light review (controller verifies). The others get the full review.

## Global Constraints

- **Owner (2026-11-17):** gloci stays free for commercial use, with no dependence on UCSC's licence. No UCSC code, binary, library or licence in gloci, its tests, `xtask/docker/Dockerfile` (comparison image) or `xtask/docker/Dockerfile.bench` (benchmark image); no `ucsc-*` package anywhere. Formats are implemented from their published description, never from UCSC source.
- **Never write the name of UCSC's source tree, or the surname of the bigWig paper's first author, anywhere:** code, comments, docs, tests, fixtures, commit messages. Cite the paper as "the bigWig/bigBed format paper (doi:10.1093/bioinformatics/btq351)". Tool names (`bedGraphToBigWig`, `bigWigAverageOverBed`, …) may be named, as the behaviour gloci replaces.
- **Licences.** New runtime dependencies must be permissive (MIT, Apache-2.0, BSD, zlib or similar). Verified with `cargo tree` on 2026-11-17 for `bigtools` 0.5.8 with `default-features = false, features = ["read", "write"]`:

  | Crate (version resolved) | Licence |
  |---|---|
  | bigtools 0.5.8, tokio 1.54, bytes 1.12, bincode 1.3.3, slab 0.4 | MIT |
  | byteordered 0.6, crossbeam-channel 0.5, crossbeam-utils 0.8, futures 0.3.34 (channel, core, executor, io, macro, sink, task, util), serde 1.0.229, smallvec 1.16, tempfile 3.27, thiserror 1.0.69, either 1.19, itertools 0.10.5, pin-project-lite 0.2 | MIT OR Apache-2.0 |
  | byteorder 1.5, memchr 2.8 | Unlicense OR MIT |
  | libdeflater / libdeflate-sys 1.26.1 (libdeflate itself is MIT) | Apache-2.0 |
  | unicode-ident (build-time, proc-macros only; already in the tree) | (MIT OR Apache-2.0) AND Unicode-3.0 |
  | **index_list 0.2.17** (a non-optional bigtools dependency) | **MPL-2.0: not permissive** |

  `index_list` is replaced by a workspace crate of gloci's own (MIT) through `[patch.crates-io]` until the owner decides (D171, Task 1). The check, run after every dependency change, must print nothing: `cargo tree --workspace -e normal --prefix none --format '{l} {p}' | grep -E 'MPL|GPL'`.
- **Test-only tools** may have other open licences if they never ship. Here: pyBigWig 0.3.29 (MIT) in the comparison image.
- **Toolchain:** `rust-toolchain.toml` pins 1.100.0; MSRV 1.91, so no std API newer than 1.91. `cargo clippy --workspace --all-targets -- -D warnings` and `cargo fmt --all -- --check` pass after every task.
- **Help rules** (`crates/gloci/tests/cli_basics.rs`):
  - a summary has at most 64 characters and names the tools it replaces in brackets;
  - the long help's first paragraph names every tool the summary names;
  - 2 to 4 `gloci ...` example lines, in an `Examples:` block that ends at the first blank line;
  - summary paragraphs have at most 12 lines;
  - `convert`, `genomecov` and `map` keep their experimental marker. `view` is core: bigWig/bigBed input is an experimental input format, said in its help and in design.md §5.1.
- **Output (D17):** the format comes from the `-o` extension (`.bw`, `.bigwig`, `.bb`, `.bigbed`, any letter case) or `-O bigwig|bigbed`. Outputs are atomic (temp file + rename: `gloci-io/src/atomic.rs`). The bytes never change with `--threads`. Memory stays flat in the input size. No silent wrong answers: an error names the line, record or construct. Exit 2 for usage errors, 1 for data errors.
- **Byte order:** bigtools writes the host's byte order. Every platform gloci builds for is little-endian, so gloci's files are little-endian. The reader accepts both byte orders.
- **Parity testing:**
  - Cases go in `xtask/diff/cases.toml` (`input = { file = ... }`, `aux = [...]`, `gloci`, `reference`, `compare`, optional `expected_difference = "#anchor"`).
  - Run one case with `cargo xtask diff --image gloci-diff:lane-a --only NAME`, then store its golden with the same command plus `--update-golden`. Goldens live in `crates/gloci/tests/golden/NAME/` and `cargo test -p gloci --test golden` replays them.
  - Only a single plain `gloci ARGS` invocation becomes a golden; a shell pipeline stays live-only.
  - An intentional difference gets a section `## Title {#anchor}` in `docs/differences.md` (gloci / Old behaviour / Reason / Pinned by), and its case pins `expected_gloci`.
- **Lanes:**
  - This lane builds and uses the image tag `gloci-diff:lane-a`. Never rebuild `gloci-diff:p0` here: the controller rebuilds it after merging the three lanes.
  - Shared files (`CHANGELOG.md`, `docs/*.md`, `SKILL.md`, `crates/gloci-py/python/gloci/__init__.pyi`, `xtask/docker/Dockerfile`, `crates/gloci/src/cli/mod.rs`, embed bindings): append self-contained entries and do not reflow neighbours.
  - Generated files are regenerated after a rebase, never merged by hand: `cargo xtask pystub && cargo xtask docs --skill SKILL.md`.
- **References:** code comments and Markdown outside `docs/superpowers/` never name a plan task ("Task 3") or lane; they cite D171 and design sections. `cargo test -p xtask --test references` enforces this.
- **Python (D87):** no new command, so no binding changes; `crates/gloci/tests/embed_args.rs` stays green. bigWig/bigBed are files. A path source reads one; `save()` on a `convert()` or `genomecov()` step writes one. An in-memory source or sink is refused with the usage error of Task 3/Task 5. `convert()` gains the keyword `autosql=` (regenerated stub).
- **Tests:** run with `cargo test -p <crate> --test <file>`. Large tests are `#[ignore]` and run nightly; never run `--ignored` suites locally.
- **Disk:** remove downloaded benchmark data, scratch outputs and the `gloci-diff:lane-a` image when the lane ends (Task 9).

## Review Focus

These are the five inputs most likely to bite a user that the spec implies but does not spell out. Each one's test is in the task named.

1. **Zero records** (an empty sample in a pipeline: `genomecov --bedgraph` of a BAM with no reads). bigtools cannot write an empty file. The run exits 1 saying there are no records, and leaves no file and no temp file (Task 5 `empty_input_is_an_error_and_leaves_no_file`; Task 7 `genomecov_empty_bedgraph_to_bigwig_is_an_error`).
2. **A truncated or corrupt bigWig/bigBed** (a cut-short download). Exit 1 naming the file and saying it is truncated or not valid bigWig, never "internal error" or a crash (Task 3 `truncated_or_corrupt_bigwig_is_a_data_error`).
3. **A bigBed record past its chromosome's end, or on a chromosome missing from `-g`.** bigtools checks only the start of a bigBed record. Exit 1 naming the line and the size source (Task 5 `end_past_chromosome_is_an_error`, `unknown_chromosome_names_the_sizes_source`).
4. **More than 65,535 chromosomes with records** (`genomecov --bga` on a draft assembly writes every scaffold). bigtools stores the chromosome count in 16 bits and would write a corrupt file. Exit 1 naming the limit, before any corruption (Task 5 `too_many_chromosomes_is_an_error`).
5. **Values a 32-bit float cannot hold or that are not finite** (`1e39`, `nan`, `inf`). Exit 1 naming the line; merely imprecise values are rounded, as documented in `#bigwig-float32` (Task 5 `value_outside_float32_is_an_error`; Task 6 golden `bigwig-float32-vs-input`).

## Decision this plan adds: D171

Task 1 writes this into `docs/decisions.md` as one entry in the house format: `<a id="d171"></a>**D171** **bigWig and bigBed (owner, 2026-11-17).** <ul><li>…</li></ul>`, one `<li>` per bullet below. Task 9 adds the measured memory.

- **Library.** `bigtools` 0.5.8 (MIT), an independent implementation of the bigWig/bigBed format paper (doi:10.1093/bioinformatics/btq351), with `read` and `write` only (no `remote`, no `cli`), and `tokio` for its runtime. Every dependency is MIT, Apache-2.0 or Unlicense, except `index_list` 0.2.17, which is MPL-2.0. That crate is replaced through `[patch.crates-io]` by `crates/index-list`, gloci's own MIT implementation of the part of its API bigtools uses, until the owner rules on MPL-2.0. A crates.io release of gloci must settle this first, because `[patch]` does not travel with a published crate.
- **Detection and input.** Formats `bigwig` and `bigbed`, from the first 4 bytes: `26 FC 8F 88` or `88 8F FC 26` for bigWig, `EB F2 89 87` or `87 89 F2 EB` for bigBed. `--format bigwig|bigbed` overrides. They are read only from a regular file, through its index. Standard input, a pipe, a `/dev/fd` path or an in-memory source is exit 2 saying so.
- **Records.** A bigWig record is a bedGraph row `chrom start end value`, with the value printed by Rust's shortest round-trip formatting of the 32-bit float: no exponent, `5` not `5.0`, `-0` for negative zero. A bigBed record is `chrom start end` plus the stored rest of the line, field for field. Both are BED family wherever WIG is (column 4 is the value). Column 6 of a bigBed gives the strand (D54). A reader does not consume to the end.
- **Order.** The reader sorts the file's chromosomes by the command's `ChromTable` (`cmp_names`: lexicographic by default, `--chrom-order`, or `-g`/a header's order). Chromosomes the table cannot place (a closed or first-appearance table that lacks them) follow, in the order of the file's chromosome list. So `convert` (first-appearance table) writes the file's own order. A chromosome is interned at its first record, so declared but empty chromosomes cost nothing. Within a chromosome, records come as stored (starts ascending).
- **Regions** (`view` positional regions, `--regions-file`; relations' region selection). They read through the file's index, never by scanning. Regions are merged per chromosome (D10) and visited in the command's order. bigWig records are **clipped to each region**, as bigtools and `bigWigToBedGraph -chrom -start -end` return them, so a record spanning two regions gives one piece in each. bigBed records are **whole**, kept when they overlap a region under spec §4.3 (with `--inclusive-points` for points), each once. Region names go through `--chrom-map` back to the file's names (D62). This differs from `view` of a BED, which never clips: a bigWig is a signal, and its region read is the signal over the region.
- **Writing.** bigWig is written by `convert` (from bedGraph, a 4-column BED, WIG or bigWig) and by `genomecov --bedgraph`/`--bga`. bigBed is written by `convert` (from BED or bigBed). Any other command, mode or input is exit 2 naming those that take it. The output must be a regular file: standard output, a device or an in-memory sink is exit 2. `--write-index` is exit 2 (the file carries its own index). Chromosome sizes come from `-g`, else from the input's own header: an alignment header for `genomecov`, the bigWig or bigBed being converted. Without them the run is exit 2 naming `-g`.
- **Input rules** (exit 1, naming the line, or the record number for genomecov output):
  - each chromosome's records are contiguous with starts ascending, in any chromosome order, because bigtools numbers chromosomes by first appearance; another order names `gloci sort`;
  - the chromosome is in the sizes and the end is within its length;
  - at least one record;
  - at most 65,535 chromosomes with records;
  - positions and lengths fit 32 bits;
  - text is UTF-8.
- **bigWig rows:** exactly 4 columns; no overlaps; no zero-length records; a value that parses as a finite number whose nearest 32-bit float is finite. That float is stored, so an inexact value is rounded (`#bigwig-float32`).
- **bigBed layout.** The first record decides it. The standard columns are the longest prefix of size 3, 4, 5, 6, 8, 9 or 12 whose fields are valid: name non-empty; score, thickStart and thickEnd 32-bit unsigned integers; strand `+`, `-` or `.`; itemRgb an unsigned integer or `R,G,B` of 0–255; blockCount/blockSizes/blockStarts valid BED12 (`bed12_blocks`). The other columns are extra. Every record must have the same column count, its standard columns valid. A record at `[0, 0)` is refused, because bigtools' reader rejects one.
- **autoSql** names the standard columns by their BED names, and the extra ones `field1`, `field2`, … as `string`. `--autosql FILE` supplies the text instead; its field count must equal the column count. bigtools writes the header's `definedFieldCount` equal to the field count, so gloci rewrites it (bytes 34–35, in the file's byte order) to the number of standard columns before the rename.
- **Runtime and determinism.** At `-t 1`, a current-thread runtime; otherwise a multi-thread runtime with `--threads` workers, as bigtools' own tools do. Sections reach the file in input order through per-chromosome channels of join handles, so timing cannot change the bytes; a test compares `-t 1`, `2` and `8`.
- **Defaults kept from bigtools.** zlib deflate through libdeflater at level 6; 1024 items per section; R-tree block size 256; zoom levels of 160 bases, ×4 each, at most 10, no larger than the longest chromosome; a level is kept only while its data is at most half the full data and has fewer sections than the level before. bigtools spools zoom data and sections to unnamed temporary files in the system temporary directory (`$TMPDIR`), not `--tmp-dir`.
- **`map` `wsum`** is the sum of value × overlapping bases (two-input commands only). With `-c 0:bases_uniq`, `-c 4:wmean` and `--empty 0` it gives `bigWigAverageOverBed`'s `sum`, `covered` and `mean`; `mean0` is `sum / (end - start)`.
- **Checking.** No UCSC tool is used. Files gloci writes are read back by pyBigWig (libBigWig, an independent C implementation) in the comparison image, value by value. gloci reads a pyBigWig-written bigWig (all three section types) and ENCODE's CTCF bigWig (the benchmark's derived input). The benchmark image's `ucsc-bigwigtobedgraph` is replaced by `gloci convert`.

## Files

| File | Responsibility |
|---|---|
| `crates/index-list/` (new crate `index_list`) | MIT implementation of the `IndexList`/`ListIndex` API bigtools calls |
| `crates/gloci-io/src/bbi/mod.rs` (new) | magic numbers, `bbi_format`, `bbi_chrom_sizes`, shared error wording |
| `crates/gloci-io/src/bbi/read.rs` (new) | `BbiReader`: whole-file and region reads |
| `crates/gloci-io/src/bbi/write.rs` (new) | `BbiWriter`: validation, writer thread, runtime, atomic commit |
| `crates/gloci-io/src/bbi/autosql.rs` (new) | bigBed layout detection, field checks, autoSql text |
| `crates/gloci-io/src/{detect,input,record,output,writer,atomic,lib}.rs` | formats, detection, dispatch, output specs, `AtomicFile::try_clone_file` |
| `crates/gloci-ops/src/view.rs` | region selection through the bbi index |
| `crates/gloci-ops/src/convert/{mod,track}.rs` | bbi → text (dispatch), text → bbi (`convert_to_track`) |
| `crates/gloci-ops/src/agg/{spec,ops}.rs` | `wsum` |
| `crates/gloci/src/cli/{files,intervals,relate,shared}.rs`, `crates/gloci/src/commands/{convert,genomecov}.rs` | options, help, command wiring |
| `xtask/diff/inputs/bigwig/` (new) | fixtures and the pyBigWig scripts |

---

### Task 1: Dependencies, the licence check and the `index_list` replacement

**Files:**
- Modify: `Cargo.toml` (workspace dependencies, `[patch.crates-io]`), `crates/gloci-io/Cargo.toml`, `Cargo.lock`
- Create: `crates/index-list/Cargo.toml`, `crates/index-list/src/lib.rs`, `crates/index-list/tests/model.rs`
- Create: `crates/gloci-io/tests/bbi_deps.rs`
- Modify: `docs/decisions.md` (D171), `docs/design.md` §8.4 (Libraries table)

**Interfaces:**
- Produces: workspace dependencies `bigtools = { version = "=0.5.8", default-features = false, features = ["read", "write"] }` and `tokio = { version = "1", default-features = false, features = ["rt", "rt-multi-thread"] }`; `gloci-io` depends on both (`bigtools.workspace = true`, `tokio.workspace = true`).
- Produces: crate `index_list` (path `crates/index-list`, version `0.2.17`, `publish = false`, `license = "MIT"`):

```rust
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct ListIndex(/* private; Default is "no element" */);
impl ListIndex { pub fn is_some(&self) -> bool; pub fn is_none(&self) -> bool; }
pub struct IndexList<T> { /* private */ }
impl<T> IndexList<T> {
    pub fn new() -> Self;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn insert_first(&mut self, elem: T) -> ListIndex;
    pub fn insert_last(&mut self, elem: T) -> ListIndex;
    /// An `index` not in the list appends, as `insert_last`.
    pub fn insert_after(&mut self, index: ListIndex, elem: T) -> ListIndex;
    pub fn remove_first(&mut self) -> Option<T>;
    pub fn get(&self, index: ListIndex) -> Option<&T>;
    pub fn get_mut(&mut self, index: ListIndex) -> Option<&mut T>;
    pub fn get_first(&self) -> Option<&T>;
    pub fn get_last(&self) -> Option<&T>;
    pub fn first_index(&self) -> ListIndex;
    pub fn next_index(&self, index: ListIndex) -> ListIndex;
    pub fn drain_iter(&mut self) -> impl Iterator<Item = T> + '_;
}
// plus Default and (for T: Debug) Debug
```

- [ ] **Step 1: The owner gate.** Ask the controller whether the owner has accepted MPL-2.0 for `index_list`. If so, skip Steps 4–7 (no `[patch]`, no crate), the licence check of Step 8 expects that one `MPL-2.0` line, and write D171's first bullet as "`index_list` 0.2.17 (MPL-2.0) is used unmodified, accepted by the owner on DATE". Otherwise do every step.
- [ ] **Step 2: Pin bigtools' bigBed behaviour against the original crate.**
  - Add the two workspace dependencies and `gloci-io`'s dependencies, but no `[patch]` yet.
  - Write `crates/gloci-io/tests/bbi_deps.rs::bigbed_zoom_counts_overlaps`. It exercises the code in bigtools' bigBed writer that uses `IndexList`.
    - Write with bigtools directly: `BigBedWrite::new(File, sizes {"chr1": 1000})`, `options.input_sort_type = InputSortType::START`, `options.manual_zoom_sizes = Some(vec![1000])`, entries `[0,100) "a"`, `[50,150) "b"`, `[120,200) "c"` through `BedParserStreamingIterator::wrap_iter(items.into_iter(), true)` (items `Result<(String, BedEntry), BedValueError>`), on a current-thread runtime.
    - Read with `BigBedRead::open(File)`. `get_zoom_interval("chr1", 0, 1000, 1000)`, summed over its records, gives coverage depth: `bases_covered == 200`, `sum == 280.0`, `sum_squares == 440.0`, max of `max_val == 2.0`, min of `min_val == 1.0`.
    - `get_interval("chr1", 0, u32::MAX)` returns the three entries, rests `a`, `b`, `c`.
- [ ] **Step 3: Run** `cargo test -p gloci-io --test bbi_deps`. Expected: PASS against crates.io's `index_list`. These numbers follow from coverage depth over the three entries; if bigtools gives others, stop and report them, as the shim must reproduce what bigtools does with the original.
- [ ] **Step 4: Write the failing tests** in `crates/index-list/tests/model.rs`:
  - `matches_a_vec_model` (proptest, up to 200 ops): ops `InsertFirst(u16)`, `InsertLast(u16)`, `InsertAfter(k, u16)` (k picks a live handle), `RemoveFirst`, `SetVia(k, u16)` (`*get_mut(handle) = v`), checked against a `Vec<(id, u16)>` model. After every op:
    - walking `first_index()`/`next_index()` yields the model's values in order;
    - `get_first`, `get_last` and `len` agree with the model;
    - `remove_first` returns the model's first.
  - `ends_are_none`: `IndexList::<u8>::new().first_index().is_none()`; on `[7]`, `next_index(first_index()).is_none()`; `get(ListIndex::default())` is `None`.
  - `drain_iter_empties_in_order`: `[1, 2, 3]` drains to `vec![1, 2, 3]`, then `len() == 0`.
  - `insert_after_unknown_index_appends`: on `[1]`, `insert_after(ListIndex::default(), 2)` gives `[1, 2]`.
- [ ] **Step 5: Run** `cargo test -p index_list --test model`. Expected: FAIL (no such package).
- [ ] **Step 6: Implement.**
  - `IndexList`: a `Vec` of slots with prev/next links, a free-slot list and head/tail. Write it from the method list above and bigtools' call sites (`src/bbi/bigbedwrite.rs`, `src/bed/indexer.rs` under `~/.cargo/registry/src/*/bigtools-0.5.8`). Do not open or copy `index_list`'s own source.
  - Crate manifest: description "gloci's own MIT implementation of the index_list API that bigtools uses (D171); replaces the MPL-2.0 crate through [patch.crates-io]"; `edition.workspace`, `rust-version.workspace`; dev-dependency `proptest.workspace = true`.
  - Root `Cargo.toml`: `[patch.crates-io]` with `index_list = { path = "crates/index-list" }`, under a comment citing D171. `crates/*` already makes the crate a member.
- [ ] **Step 7: Run** `cargo test -p index_list --test model && cargo test -p gloci-io --test bbi_deps`. Expected: PASS; `bbi_deps` now runs on the shim and must give Step 3's numbers. `cargo tree -i index_list` prints `index_list v0.2.17 (…/crates/index-list)`.
- [ ] **Step 8: Licence check.** `cargo tree --workspace -e normal --prefix none --format '{l} {p}' | grep -E 'MPL|GPL'` prints nothing. Every crate the change added is in the Global Constraints table. A crate that is not, or a licence that differs: stop and report.
- [ ] **Step 9: Build checks CI runs.**
  - If `rustup toolchain list` shows 1.91, run `cargo +1.91 check --workspace --locked`; otherwise CI's `msrv` job covers it.
  - Run CI's `cargo publish --dry-run --locked -p gloci-core -p gloci-io -p gloci-engine -p gloci-vcf -p gloci-ops -p gloci` (slow). If it fails over `index_list` or the lock file, stop and report the error: that job must stay green, and the fix is the controller's call.
- [ ] **Step 10: Docs.**
  - `docs/decisions.md`: append D171 (this plan's D171 section, every bullet).
  - `docs/design.md` §8.4: add a row `| bigWig, bigBed | bigtools 0.5.8 (read, write) on a tokio runtime | MIT, an independent implementation of the published format; index_list replaced by an MIT crate (D171) |`.
- [ ] **Step 11: Run** `cargo clippy --workspace --all-targets -- -D warnings && cargo fmt --all -- --check && cargo test -p xtask --test references`. Expected: PASS.
- [ ] **Step 12: Commit**

```bash
git add Cargo.toml Cargo.lock crates/index-list crates/gloci-io/Cargo.toml crates/gloci-io/tests/bbi_deps.rs docs/decisions.md docs/design.md
git commit -m "build: bigtools 0.5.8 (read, write) and tokio for bigWig and bigBed; index_list replaced by an MIT crate (D171)"
```

### Task 2: pyBigWig in the comparison image, and its scripts and fixture (light review)

**Files:**
- Modify: `xtask/docker/Dockerfile` (append only)
- Create: `xtask/diff/inputs/bigwig/bbi_dump.py`, `xtask/diff/inputs/bigwig/make_pybigwig.py`, `xtask/diff/inputs/bigwig/pybigwig.bw` (generated), `xtask/diff/inputs/bigwig/pybigwig.lex.bedgraph` (written by hand)

**Interfaces:**
- Produces: image `gloci-diff:lane-a`, with `python3 -I bbi_dump.py` usable as a case's `reference` (list `xtask/diff/inputs/bigwig/bbi_dump.py` in the case's `aux`).
- Produces: `bbi_dump.py FILE [--lex] [--region CHROM:START-END]...`:
  - opens FILE with pyBigWig and prints one line per record to stdout;
  - chromosomes come in `chroms()` order, or sorted by name with `--lex`;
  - `--region` (1-based, closed, like gloci's region strings) prints only those regions, in the order given: bigWig intervals clipped to `[START-1, END)`, bigBed entries whole when `start < END` and `end > START-1`;
  - a bigWig line is `chrom start end value`, with the value from `f32()` below; a bigBed line is `chrom start end`, then a tab and the rest when the rest is non-empty.
- Produces: fixture `pybigwig.bw`, header order `chr1 5000, chr2 3000, chr10 2000, chrM 100`, with one section of each type and an empty chromosome.

- [ ] **Step 1: Append to the Dockerfile.** Use its own `RUN` line, so the other lanes' edits do not touch it:

```dockerfile
# pyBigWig (MIT, on libBigWig) reads gloci's bigWig and bigBed files back, value by value (D171).
RUN micromamba install -y -n base -c conda-forge -c bioconda pybigwig=0.3.29 \
    && micromamba clean --all --yes
```

- [ ] **Step 2: Write `bbi_dump.py`.** Standard library plus `pyBigWig` only. Its value formatter, the shortest decimal that reads back as the same 32-bit float, without an exponent (what Rust's `{}` prints for an `f32`):

```python
def f32(v: float) -> str:
    if v == 0.0:
        return "-0" if math.copysign(1.0, v) < 0 else "0"
    for digits in range(1, 10):
        s = "%.*e" % (digits - 1, v)
        if struct.unpack("<f", struct.pack("<f", float(s)))[0] == v:
            break
    t = format(decimal.Decimal(s), "f")
    return t.rstrip("0").rstrip(".") if "." in t else t
```

- [ ] **Step 3: Write `make_pybigwig.py`.** It writes `pybigwig.bw` in the current directory:

```python
import pyBigWig
bw = pyBigWig.open("pybigwig.bw", "w")
bw.addHeader([("chr1", 5000), ("chr2", 3000), ("chr10", 2000), ("chrM", 100)])
bw.addEntries(["chr1", "chr1", "chr1"], [0, 100, 250], ends=[100, 200, 300], values=[1.5, 0.25, -2.0])
bw.addEntries("chr2", [10, 20, 30], values=[3.0, 4.5, 0.1], span=5)
bw.addEntries("chr10", 100, values=[7.0, 8.0, 9.0], span=10, step=20)
bw.close()
```

- [ ] **Step 4: Write `pybigwig.lex.bedgraph` by hand** (tab-separated, `\n` endings):

```
chr1	0	100	1.5
chr1	100	200	0.25
chr1	250	300	-2
chr10	100	110	7
chr10	120	130	8
chr10	140	150	9
chr2	10	15	3
chr2	20	25	4.5
chr2	30	35	0.1
```

- [ ] **Step 5: Build and check the pins.** Run `docker build -t gloci-diff:lane-a xtask/docker`, then `docker run --rm gloci-diff:lane-a sh -c 'samtools --version | head -1; bcftools --version | head -1; bedtools --version; python3 -I -c "import pyBigWig; print(pyBigWig.__version__)"'`. Expected: `samtools 1.21`, `bcftools 1.21`, `bedtools v2.31.1`, `0.3.29`. If the solver changed a pinned tool, add `htslib=1.21 bcftools=1.21 samtools=1.21` to the new `RUN` line and rebuild.
- [ ] **Step 6: Make the fixture and check the script.** Run `docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD/xtask/diff/inputs/bigwig:/w" -w /w gloci-diff:lane-a python3 -I make_pybigwig.py`. Then run the same with `python3 -I bbi_dump.py --lex pybigwig.bw | diff - pybigwig.lex.bedgraph`. Expected: no output.
- [ ] **Step 7: Commit**

```bash
git add xtask/docker/Dockerfile xtask/diff/inputs/bigwig
git commit -m "test(diff): pyBigWig 0.3.29 in the comparison image, a dump script and a pyBigWig-written bigWig fixture"
```

### Task 3: Reading bigWig and bigBed

**Files:**
- Create: `crates/gloci-io/src/bbi/mod.rs`, `crates/gloci-io/src/bbi/read.rs`
- Modify:
  - `crates/gloci-io/src/detect.rs`: variants, `name`, `FromStr` and its message, `is_bbi`;
  - `crates/gloci-io/src/input.rs`: detection, refusal of non-files;
  - `crates/gloci-io/src/record.rs`: the `dispatch_raw` arm;
  - `crates/gloci-io/src/output.rs`: `OutFormat::family`;
  - `crates/gloci-io/src/lib.rs`;
  - every `match` on `Format` the compiler names (in `tsv.rs`, `index/`, `gloci-ops/src/{check,map}.rs`, `relate/{files,rewrite}.rs`, `convert/mod.rs`);
  - `crates/gloci/src/commands/index.rs`: refusal;
  - `crates/gloci/src/cli/shared.rs`: `--format`/`--a-format` help;
  - `crates/gloci/Cargo.toml`: dev-dependencies `bigtools.workspace = true`, `tokio.workspace = true`.
- Test: `crates/gloci-io/tests/bbi_read.rs`, `crates/gloci/tests/bigwig.rs`
- Modify: `xtask/diff/cases.toml`, `crates/gloci/tests/golden/bigwig-view-pybigwig-file-vs-pybigwig/` (generated), `docs/design.md` §5.1, `docs/migration.md`, `CHANGELOG.md`, `SKILL.md`, `crates/gloci-py/python/gloci/__init__.pyi`

**Interfaces:**
- Produces in `gloci_io` (re-exported from `lib.rs`):

```rust
pub enum Format { /* … */ BigWig, BigBed }      // name(): "bigWig", "bigBed"; FromStr: "bigwig", "bigbed"
impl Format { pub fn is_bbi(self) -> bool }
pub const BIGWIG_MAGIC: u32 = 0x888F_FC26;
pub const BIGBED_MAGIC: u32 = 0x8789_F2EB;
/// The bbi format whose magic number, in either byte order, starts `head`.
pub fn bbi_format(head: &[u8]) -> Option<Format>;
/// The chromosomes a bigWig or bigBed declares, as (file's own name, length), in its chromosome list's order.
pub fn bbi_chrom_sizes(path: &Path) -> Result<Vec<(Vec<u8>, u64)>, IoError>;
pub(crate) fn open_bbi(path: &str, format: Format, o: &ReaderOptions) -> Result<Box<dyn RecordReader>, IoError>;
```

- `OutFormat::family(Format::BigWig) == OutFormat::BedGraph`, `family(Format::BigBed) == OutFormat::Bed`.
- Rule for the compiler's list of `match`es: wherever `Format::Wig` sits in a BED-family arm, add `Format::BigWig | Format::BigBed` to that arm. `check.rs` names them `BIGWIG` and `BIGBED`.

- [ ] **Step 1: Write the failing tests.**
  - `crates/gloci-io/tests/bbi_read.rs::bbi_magic_both_byte_orders`: `bbi_format` maps `[0x26,0xFC,0x8F,0x88]` and `[0x88,0x8F,0xFC,0x26]` to `BigWig`, `[0xEB,0xF2,0x89,0x87]` and `[0x87,0x89,0xF2,0xEB]` to `BigBed`, and `b"chr1\t0\t1\n"` and a 3-byte slice to `None`.
  - The fixtures in `crates/gloci/tests/bigwig.rs` come from two helpers that call bigtools directly (not gloci), records in the order given: `bigtools_bigwig(path, sizes: &[(&str, u32)], recs: &[(&str, u32, u32, f32)])` and `bigtools_bigbed(path, sizes, recs: &[(&str, u32, u32, &str)])`, with `InputSortType::START` and `wrap_iter(.., true)`.
    - **K** (`k.bw`): sizes `chr1 1000, chr2 1000, chr10 1000`; records in stored order `chr2 0 10 0.1`, `chr2 10 20 5`, `chr10 5 6 (1.0f32/3.0)`, `chr1 100 200 -2.5`, `chr1 300 301 1e-7`.
    - **B** (`b.bb`): sizes `chr1 1000`; `10 20 "a\t0\t+"`, `15 40 "b\t5\t-"`, `500 900 "c\t9\t."`.
  - `reads_bigwig_in_lex_order`: `gloci view k.bw` prints exactly:

```
chr1	100	200	-2.5
chr1	300	301	0.0000001
chr10	5	6	0.33333334
chr2	0	10	0.1
chr2	10	20	5
```

  - `reads_in_the_commands_order`:
    - `view --chrom-order karyotype k.bw` gives the chromosome blocks `chr1, chr2, chr10`;
    - `view -g order.genome k.bw`, where `order.genome` lists `chr10 1000`, `chr2 1000`, `chr1 1000`, gives `chr10, chr2, chr1`;
    - `convert k.bw` gives the stored order `chr2, chr10, chr1`, rows as above.
  - `map_reads_a_bigwig_b_in_lex_order`: `map -a a.bed -b k.bw -c 4:count`, with `a.bed` = `chr1 0 1000`, `chr10 0 1000`, `chr2 0 1000`, exits 0 and its last columns are `2`, `1`, `2`.
  - `empty_chromosomes_and_g_without_them`: `xtask/diff/inputs/bigwig/pybigwig.bw` (Task 2) declares `chrM` with no records (bigtools never writes an empty chromosome). `view -g g3.genome pybigwig.bw`, with `g3.genome` = `chr1 5000`, `chr2 3000`, `chr10 2000`, exits 0 with the 9 rows.
  - `bigbed_records_are_bed_lines`:
    - `view b.bb` prints `chr1\t10\t20\ta\t0\t+`, `chr1\t15\t40\tb\t5\t-`, `chr1\t500\t900\tc\t9\t.`;
    - `intersect -a b.bb -b q.bed -s`, with `q.bed` = `chr1 0 1000 q 0 +`, prints only record `a`;
    - `convert b.bb` prints the same three rows.
  - `chrom_map_renames_bigwig_chromosomes`: with `map.tsv` mapping `chr1→1`, `chr2→2`, `chr10→10`, `view --chrom-map map.tsv k.bw` prints chromosomes `1`, `10`, `2` in that order.
  - `bigwig_must_come_from_a_file`:
    - `gloci view -` with `k.bw` on stdin exits 2, and stderr contains `bigWig` and `file`;
    - `gloci view /dev/stdin < k.bw` exits 2 the same way.
  - `format_override`:
    - `view --format bigwig a.bed` fails, and stderr contains `a.bed` and `not a bigWig`;
    - `view --format bigbed k.bw` fails, and stderr contains `not a bigBed`;
    - `view --format bigwig k.bw` succeeds.
  - **Review Focus 2**, `truncated_or_corrupt_bigwig_is_a_data_error`: write K with 3,000 records on `chr1` (`[i*10, i*10+5)`, value `i`).
    - Truncate copies to half their length and to 40 bytes: `view` of each exits 1, and stderr names the file, contains `truncated` and does not contain `internal error`.
    - A copy with 64 bytes at offset 1024 set to `0xFF`: `view` exits 1 without `internal error`.
  - `index_refuses_bigwig`: `gloci index k.bw` exits 2, and stderr contains `own index`.
  - `reading_does_not_depend_on_threads`: `view k.bw` at `-t 1` and `-t 8` prints the same bytes.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test bbi_read && cargo test -p gloci --test bigwig`. Expected: FAIL (no `bbi_format`; bigWig files not recognised).
- [ ] **Step 3: Implement.**
  - **`input.rs`.** In `open_input`, a path that `is_regular_file` and whose first 4 bytes are a bbi magic, or whose `opts.format` is `BigWig`/`BigBed`, returns `Input { name: path, format, compression: None, reader: empty }` without decompressing. In `open_input_from` (stdin, pipes, fd paths, in-memory and spooled sources), a raw peek with a bbi magic, or such an `opts.format`, is `UsageError`: "NAME is bigWig, which gloci reads through its own index and so only from a file, not from standard input, a pipe or memory; save it to a file first".
  - **`record.rs`.** `Format::BigWig | Format::BigBed => open_bbi(&input.name, input.format, o)`.
  - **`read.rs`.** `BbiReader` opens `std::fs::File` with `BigWigRead::open`/`BigBedRead::open`, not `open_file`, which prints to stderr.
    - On the first `read_batch` it sorts the file's chromosomes (renamed by `--chrom-map`; `ChromMap::check_file_names` first) by D171's order rule, using `chroms.cmp_names`.
    - It walks each chromosome with `get_interval_move(file_name, 0, u32::MAX)`, getting the reader back through `Into<BigWigRead<_>>` when a chromosome ends.
    - Each row is formatted into the batch (bigWig: `write!(buf, "{name}\t{start}\t{end}\t{value}")` with `value: f32`; bigBed: `chrom\tstart\tend` then `\t` and the rest when non-empty). The `RecordMeta` gets the chromosome from `intern_chrom` at its first record, `key = start`, `span = [start, end)`, the strand from column 6 for bigBed (D54) and `Unknown` for bigWig, and `loc = Loc::Record(n)` (n counts records yielded).
    - `header()` is `Header::None`, `compression()` is `Compression::None`, `consume_to_end()` is false.
    - Every bigtools call runs under `crate::atomic::panic_note(format!("{path} is not a valid bigWig file (it may be truncated)"))`, or bigBed, so a bigtools panic is a data error. bigtools' I/O errors map to `IoError::Detect` with the same wording plus bigtools' message.
  - **`index.rs`.** `gloci index` on a bbi file: `UsageError("NAME is bigWig, which carries its own index; gloci reads it directly")`.
  - **`convert`.** `gloci-ops/src/convert/mod.rs` dispatches `Format::BigWig | Format::BigBed` to a pass-through like `wig::convert`. `what()` gains "bigWig"/"bigBed", and `SUMMARY` gains "bigWig to bedGraph, bigBed to BED". The existing option checks then refuse conversion options for them.
  - **Help.** The `--format` help in `shared.rs` adds "bigwig or bigbed (from a file)"; so does `--a-format`.
- [ ] **Step 4: Run** the commands of Step 2. Expected: PASS.
- [ ] **Step 5: Golden case** in `xtask/diff/cases.toml`:

```toml
[[case]]
name = "bigwig-view-pybigwig-file-vs-pybigwig"
input = { file = "xtask/diff/inputs/bigwig/pybigwig.bw" }
aux = ["xtask/diff/inputs/bigwig/bbi_dump.py"]
# A file pyBigWig wrote, with bedGraph, variableStep and fixedStep sections, in header (karyotype) order.
gloci = "gloci view input.bw"
reference = "python3 -I bbi_dump.py --lex input.bw"
compare = "exact"
```

Run `cargo xtask diff --image gloci-diff:lane-a --only bigwig-view-pybigwig-file-vs-pybigwig`; expected `ok (matches)`. Then add `--update-golden`. The test `golden_expected_is_the_hand_written_listing` in `crates/gloci/tests/bigwig.rs` asserts that the golden's `expected` equals `xtask/diff/inputs/bigwig/pybigwig.lex.bedgraph`.
- [ ] **Step 6: Docs.**
  - `docs/design.md` §5.1: add a table row `| bigWig, bigBed (experimental) | read from a file through its own index, as bedGraph and BED records, in the command's chromosome order (D171) | — |`, and add the two magic numbers to the detection bullet.
  - `docs/migration.md`: a new subsection `### bigWig and bigBed tools` immediately before `### Filter expressions (bcftools \`-i\` / \`-e\`)`, with rows `bigWigToBedGraph in.bw out.bedGraph` → `gloci convert in.bw -o out.bedgraph` (the file's own order; `gloci view in.bw` gives lexicographic order) and `bigBedToBed in.bb out.bed` → `gloci convert in.bb -o out.bed`.
  - `CHANGELOG.md` `### Added`: one bullet, "bigWig and bigBed input (experimental): read from a file through its index, wherever intervals are read, …".
- [ ] **Step 7: Regenerate and check.** `cargo xtask pystub && cargo xtask docs --skill SKILL.md`, then `cargo test -p gloci-io && cargo test -p gloci --test cli_basics --test golden --test input_gate --test embed_args --test wig --test convert && cargo test -p xtask --test pystub --test docs --test references --test book`, then clippy and fmt. Expected: PASS.
- [ ] **Step 8: Commit** `feat(io): read bigWig and bigBed through their index, in the command's chromosome order (D171)`

### Task 4: Region reads through the file's index

**Files:**
- Modify: `crates/gloci-io/src/bbi/read.rs`, `crates/gloci-io/src/lib.rs`, `crates/gloci-ops/src/view.rs` (`select`, `resolve_regions`), `crates/gloci/src/cli/files.rs` (`view` help, "Regions:" section)
- Test: `crates/gloci/tests/bigwig.rs`; case and golden `bigwig-region-vs-pybigwig`
- Modify: `docs/migration.md` (one row), `docs/design.md` §7.4 (`view` paragraph), `SKILL.md`, `__init__.pyi`

**Interfaces:**
- Consumes: `BbiReader` (Task 3).
- Produces:

```rust
/// The records of the bigWig or bigBed at `path` within `regions` (named as records are after
/// --chrom-map), read through its index: bigWig rows clipped to each region, bigBed rows whole and
/// once, by spec §4.3 with `inclusive_points`; chromosomes in the order of the table read_batch gets.
pub fn open_bbi_regions(path: &Path, regions: &[(Vec<u8>, Interval)], o: &ReaderOptions,
                        inclusive_points: bool) -> Result<Box<dyn RecordReader>, IoError>;
```

- `view::select`: for a bbi format, `Src::Indexed(open_bbi_regions(..))`, never `Src::Scan`. `resolve_regions` takes the file's chromosome names (renamed) from `bbi_chrom_sizes` as complete `known` names, so a region on a missing chromosome warns as it does for an indexed file.

- [ ] **Step 1: Write the failing tests** in `crates/gloci/tests/bigwig.rs`:
  - `region_read_clips_bigwig_records`: `bigtools_bigwig` with `chr1 0 100 1.5`, `chr1 100 200 0.25` (size 1000). `view r.bw chr1:51-150` prints exactly `chr1\t50\t100\t1.5\nchr1\t100\t150\t0.25\n`.
  - `bigwig_record_spanning_two_regions_gives_two_pieces`: `view r.bw chr1:11-20 chr1:31-40` prints `chr1\t10\t20\t1.5\nchr1\t30\t40\t1.5\n`.
  - `bigbed_region_returns_whole_records_once`: a bigBed record `chr1 10 2000 "big"` (size 5000) and regions `chr1:21-30 chr1:1001-1010` print `chr1\t10\t2000\tbig` once.
  - `bigbed_region_excludes_bookended`: a record `chr1 0 10 "x"` and region `chr1:11-20` print nothing (bigtools' own query also returns records that only touch the region; gloci applies §4.3).
  - `regions_file_reads_bigwig`: `view --regions-file reg.bed r.bw`, with `reg.bed` = `chr1\t50\t150`, prints the same as the first test.
  - `region_on_a_missing_chromosome_warns`: `view r.bw chrZ:1-10` exits 0 with empty stdout, and stderr contains `chrZ` and `r.bw`.
  - `region_names_go_through_the_chrom_map`: with `map.tsv` mapping `chr1→1`, `view --chrom-map map.tsv r.bw 1:51-150` prints the first test's rows with chromosome `1`.
  - `bigbed_region_points`: a bigBed point `chr1 50 50 "p"` and region `chr1:1-50` print it only with `--inclusive-points` (§4.3).
- [ ] **Step 2: Run** `cargo test -p gloci --test bigwig region`. Expected: FAIL (regions scan instead: wrong rows).
- [ ] **Step 3: Implement.**
  - **bigWig.** For each merged region, call `get_interval_move(file_name, start, end)`; the rows come back clipped.
  - **bigBed, the once-only rule.** Regions on a chromosome are sorted and disjoint, so a record that overlaps region *i* and an earlier one also overlaps region *i−1*. Emit a record from region *i*'s query unless it overlaps region *i−1*, then filter with `gloci_core::overlaps(span, region, inclusive_points)`.
  - **Region names.** Translate each to the file's name with `ChromMap::file_name`; a name the file lacks yields nothing.
  - **View help.** In the "Regions:" paragraph, add: "A bigWig or bigBed (experimental) is read through its own index: bigWig records are clipped to each region, bigBed records are whole (D171)."
- [ ] **Step 4: Run** the tests of Step 1. Expected: PASS.
- [ ] **Step 5: Golden case**

```toml
[[case]]
name = "bigwig-region-vs-pybigwig"
input = { file = "xtask/diff/inputs/bigwig/pybigwig.bw" }
aux = ["xtask/diff/inputs/bigwig/bbi_dump.py"]
# chr1 [50,150) cuts two bedGraph records; chr2 [11,31) cuts variableStep records at both ends.
gloci = "gloci view input.bw chr1:51-150 chr2:12-31"
reference = "python3 -I bbi_dump.py --region chr1:51-150 --region chr2:12-31 input.bw"
compare = "exact"
```

Run it with `cargo xtask diff --image gloci-diff:lane-a --only bigwig-region-vs-pybigwig`, then with `--update-golden`. Its expected output is `chr1 50 100 1.5`, `chr1 100 150 0.25`, `chr2 11 15 3`, `chr2 20 25 4.5`, `chr2 30 31 0.1`.
- [ ] **Step 6: Docs.**
  - `docs/migration.md` (bigWig subsection): add a row `bigWigToBedGraph -chrom=chr1 -start=0 -end=1000 in.bw out` → `gloci view in.bw chr1:1-1000` (records clipped to the region; several regions and `--regions-file` work too).
  - `docs/design.md` §7.4, `view` paragraph: add one sentence with the same rule, citing D171.
  - Regenerate SKILL.md and the stub.
- [ ] **Step 7: Run** `cargo test -p gloci --test bigwig --test golden --test cli_basics --test relation_regions && cargo test -p gloci-ops`, plus the xtask checks, clippy and fmt. Expected: PASS.
- [ ] **Step 8: Commit** `feat(view): bigWig and bigBed regions through the file's index; bigWig rows clipped to the region (D171)`

### Task 5: The bigWig/bigBed writer

**Files:**
- Create: `crates/gloci-io/src/bbi/write.rs`, `crates/gloci-io/src/bbi/autosql.rs`
- Modify:
  - `crates/gloci-io/src/output.rs`: `OutFormat::{BigWig, BigBed}`, the lookup table, `parse_track_output_spec`, the refusal;
  - `crates/gloci-io/src/writer.rs`: the `family_name`, `family_mismatch` and `check_family` arms;
  - `crates/gloci-io/src/atomic.rs`: `try_clone_file`;
  - `crates/gloci-io/src/lib.rs`.
- Test: `crates/gloci-io/tests/bbi_write.rs`

**Interfaces:**
- Consumes: `bbi_format`, `BIGBED_MAGIC` (Task 3); `gloci_core::Genome`; `crate::bed12_blocks`.
- Produces:

```rust
pub enum OutFormat { /* … */ BigWig, BigBed }   // name(): "bigWig", "bigBed"; index_format(): None
/// As parse_output_spec, accepting .bw/.bigwig/.bb/.bigbed and -O bigwig|bigbed: for convert and genomecov.
pub fn parse_track_output_spec(output: Option<&Path>, output_format: Option<&str>) -> Result<OutputSpec, UsageError>;
// parse_output_spec: a bigWig or bigBed extension or -O is UsageError(bbi_output_refusal(f)).
pub fn bbi_output_refusal(f: OutFormat) -> String;
impl AtomicFile { pub fn try_clone_file(&self) -> io::Result<std::fs::File>; }

#[derive(Clone, Copy, Debug, PartialEq, Eq)] pub enum BbiKind { BigWig, BigBed }
/// Chromosome sizes and how to name their source in errors ("-g hg38", "the header of reads.bam").
pub struct BbiSizes { pub genome: Genome, pub source: String }
pub struct BbiWriteOptions { pub kind: BbiKind, pub sizes: BbiSizes, pub threads: NonZeroUsize, pub autosql: Option<String> }
pub enum PushError { /// The record is refused; the caller names where it came from.
                     Record(String), /// The writer stopped; its error.
                     Writer(IoError) }
pub struct BbiWriter { /* AtomicFile, sender, writer thread, checks */ }
impl BbiWriter {
    /// Usage error for stdout, a device or pipe (`spec.direct`), an in-memory sink, or a spec that is not bigWig/bigBed.
    pub fn create(spec: &OutputSpec, o: BbiWriteOptions) -> Result<BbiWriter, IoError>;
    /// One bedGraph (BigWig) or BED (BigBed) line, no newline.
    pub fn push(&mut self, raw: &[u8]) -> Result<(), PushError>;
    /// Ends the input, joins the writer, fixes definedFieldCount (bigBed), renames into place.
    pub fn finish(self) -> Result<(), IoError>;
}
impl RecordWriter for BbiWriter {}  // write_record = push (a Record error says "record N of the output: …"); finish = finish
// Dropping a BbiWriter without finish aborts the write; the temp file is removed.

// autosql.rs (pub, for tests)
pub struct BedLayout { pub columns: usize, pub standard: usize }
pub fn bed_layout(first: &[u8]) -> Result<BedLayout, String>;
pub fn check_bed_line(line: &[u8], layout: &BedLayout) -> Result<(), String>;
pub fn bed_autosql(layout: &BedLayout) -> String;
pub fn autosql_field_count(text: &str) -> Result<usize, String>;  // through bigtools::bed::autosql::parse::parse_autosql
```

- [ ] **Step 1: Write the failing tests** in `crates/gloci-io/tests/bbi_write.rs`. Sizes are built with `Genome::new`; files are read back with bigtools; the output spec comes from `parse_track_output_spec`.
  - `output_spec_rules`:
    - `parse_output_spec(Some("x.bw"), None)` is `Err`, with a message containing `gloci convert` and `genomecov --bedgraph`;
    - so is `parse_output_spec(Some("x.bb"), None)`, its message containing `gloci convert`;
    - `parse_track_output_spec` gives `BigWig` for `x.bw`, `x.BigWig` and (`None`, `-O bigwig`), and `BigBed` for `x.bb` and `x.bigbed`;
    - `-o x.bw -O bigbed` disagree (Err).
  - `refuses_non_file_outputs`: `create` with `path: None`, with `/dev/null` (direct), and with an in-memory name: each is `IoError::Usage`, its message containing `regular file`.
  - `bigwig_round_trip`: push `chr2 0 10 0.5`, `chr10 0 5 3`, `chr1 0 1 -1`. `BigWigRead` gives `chroms()` names `chr2, chr10, chr1` (any chromosome order, kept) and the three values.
  - `bytes_do_not_depend_on_threads`: 3,000 records on 3 chromosomes, written at `threads` 1, 2 and 8: identical bytes.
  - `split_chromosome_is_an_error`: `chr1`, `chr2`, `chr1` gives `Record` containing `chr1` and `gloci sort`.
  - `starts_must_ascend`: `chr1 50 60 1`, then `chr1 10 20 1`, gives `Record` containing `gloci sort`.
  - `bigwig_overlap_is_an_error`: `chr1 0 10 1`, then `chr1 5 20 2`, gives `Record` containing `overlap` and `chr1:0-10`.
  - `bigwig_row_shape`:
    - `chr1 0 10` (3 columns) and `chr1 0 10 1 x` (5 columns) give `Record` containing `4 columns`;
    - `chr1 0 10 geneA` gives `Record` containing `geneA` and `not a number`;
    - `chr1 5 5 1` gives `Record` containing `zero-length`.
  - **Review Focus 5**, `value_outside_float32_is_an_error`: `1e39`, `-1e39`, `nan`, `inf` and `-inf` each give `Record` naming the value. `0.123456789` is accepted and reads back as `0.12345679_f32`.
  - **Review Focus 3**:
    - `end_past_chromosome_is_an_error`, for both kinds: size `chr1 1000`, record `chr1 900 1001 …`, gives `Record` containing `1000`;
    - `unknown_chromosome_names_the_sizes_source`: a record on `chrZ` gives `Record` containing `chrZ` and the `source` string;
    - a chromosome whose length the genome does not give (`None`) gives `Record` containing `length`.
  - `chromosome_longer_than_32_bits_is_an_error`: size `chrBig 5_000_000_000`; a record on it gives `Record` containing `32-bit`.
  - **Review Focus 4**, `too_many_chromosomes_is_an_error`: 65,536 chromosomes of length 10, one record each. The 65,536th push gives `Record` containing `65,535`.
  - **Review Focus 1**, `empty_input_is_an_error_and_leaves_no_file`: `finish` with no pushes is `Err`, containing `no records`. The destination does not exist, and its directory holds no `*.tmp-*` file.
  - `drop_without_finish_removes_the_temp_file`: push 10 records and drop: no destination, no temp file.
  - `bigbed_layout_rules` (`bed_layout`):

    | First record | `standard` | `columns` |
    |---|---|---|
    | `c 0 1` | 3 | 3 |
    | `c 0 1 n` | 4 | 4 |
    | `c 0 1 n x` (score not an integer) | 4 | 5 |
    | `c 0 1 n 5 x` (strand invalid) | 5 | 6 |
    | `c 0 1 n 5 + 12.5 x y z` (narrowPeak-like) | 6 | 10 |
    | `c 0 1 n 5 + 0` (BED7) | 6 | 7 |
    | `c 0 10 n 5 + 0 10 255,0,0` | 9 | 9 |
    | valid BED12 `c 0 10 n 5 + 0 10 0 2 3,4, 0,6,` | 12 | 12 |
    | the same with blockCount 3 | 9 | 12 |

  - `bigbed_autosql_and_header`: BED6+2 records written (`chr1 0 10 n 5 + 1.5 x`). The file's autoSql (`BigBedRead::autosql`) equals exactly:

```
table bed
"BED records written by gloci"
(
    string chrom; "Chromosome"
    uint chromStart; "Start, 0-based"
    uint chromEnd; "End, exclusive"
    string name; "Name"
    uint score; "Score"
    char[1] strand; "Strand: +, - or ."
    string field1; "Extra column 1"
    string field2; "Extra column 2"
)
```

    with a final newline. Header bytes `[32..34]` read as u16 LE are 8 (fieldCount) and `[34..36]` are 6 (definedFieldCount). The standard rows above it in `bed_autosql` are `uint thickStart; "Start of the thick part"`, `uint thickEnd; "End of the thick part"`, `uint itemRgb; "Colour, R,G,B"`, `int blockCount; "Number of blocks"`, `int[blockCount] blockSizes; "Block sizes"` and `int[blockCount] blockStarts; "Block starts, relative to chromStart"`. `autosql_field_count` of each generated text equals `columns`.
  - `bigbed_records_must_keep_the_layout`:
    - a second record with 7 columns after an 8-column first gives `Record` containing `8 columns`;
    - a second record whose score is `x`, under a layout with `standard >= 5`, gives `Record` containing `score`.
  - `bigbed_given_autosql_must_match`: `autosql: Some(text with 5 fields)` and an 8-column first record give `Record` containing `5` and `8`.
  - `bigbed_bytes_and_zero_length`:
    - a line with byte `0xE9` gives `Record` containing `UTF-8`;
    - `chr1 0 0 x` gives `Record` containing `0`;
    - `chr1 5 5 x` is accepted.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test bbi_write`. Expected: FAIL to compile (no `BbiWriter`).
- [ ] **Step 3: Implement.**
  - **`create`** makes the `AtomicFile` and hands bigtools `try_clone_file()`, because bigtools drops the writer it is given. It builds the runtime per D171: `Builder::new_current_thread()` at 1, `new_multi_thread().worker_threads(n)` otherwise. It spawns a named thread (`gloci-bbi-writer`) running `BigWigWrite::new(file, sizes).write(..)` or `BigBedWrite` (with `options.input_sort_type = InputSortType::START`, `autosql`), over `BedParserStreamingIterator::wrap_iter(iter, true)`. The iterator flattens batches of 4,096 records from a `crossbeam_channel::bounded(4)`; an abort message yields one `Err(BedValueError::InvalidInput(..))`.
  - **Sizes.** bigtools gets a `HashMap<String, u32>` of the UTF-8 names with lengths at most `u32::MAX`; the others are kept aside for the `Record` errors above.
  - **`push`** validates in this order: UTF-8; shape and value (bigWig) or layout (bigBed: `bed_layout` at the first record, then `check_bed_line`); coordinates parsed as u64; chromosome known with a length, at most 32 bits, end within it; the new-chromosome count; contiguity, start order and overlap (bigWig). A refused record sends nothing. A send that fails joins the thread and returns `Writer`.
  - **`finish`** with zero records returns its own "no records to write; a bigWig or bigBed needs at least one record" error without waiting on bigtools. It maps `BBIProcessError::IoError` to `IoError::Write { path: dest }` and anything else to `IoError::Internal`. For bigBed it rewrites bytes 34–35 to `standard`, in the byte order the magic at offset 0 shows. Then it commits.
  - **`bbi_output_refusal`:**
    - bigWig: "bigWig output is written only by gloci convert (from bedGraph, WIG or bigWig) and gloci genomecov --bedgraph or --bga; write bedGraph here (-o x.bedgraph) and convert it";
    - bigBed: "bigBed output is written only by gloci convert (from BED or bigBed); write BED here (-o x.bed) and convert it with gloci convert x.bed -g GENOME -o x.bb".
  - **`check_family`** returns the same refusal for these two outputs, before the family check.
- [ ] **Step 4: Run** `cargo test -p gloci-io --test bbi_write`. Expected: PASS.
  - If `too_many_chromosomes_is_an_error` fails with "Too many open files": bigtools spools each chromosome still waiting for the file to an unnamed temporary file. Raise the soft `RLIMIT_NOFILE` to the hard limit in `create` (`libc::setrlimit`, Unix only; `libc` is a workspace dependency) and note it in D171. If that is not enough, stop and report.
- [ ] **Step 5: Run** `cargo test --workspace --no-fail-fast`, which catches every `match` on `OutFormat` and the commands that now refuse `.bw`; then clippy and fmt. Expected: PASS.
- [ ] **Step 6: Commit** `feat(io): a bigWig and bigBed writer on bigtools: checked input, atomic output, the same bytes at any thread count (D171)`

### Task 6: `gloci convert` to and from bigWig and bigBed

**Files:**
- Create: `crates/gloci-ops/src/convert/track.rs`
- Modify:
  - `crates/gloci-ops/src/convert/mod.rs`;
  - `crates/gloci/src/commands/convert.rs`;
  - `crates/gloci/src/cli/files.rs`: `ConvertArgs` gets `--autosql` and a `-g` that no longer `requires = "bam"`, plus the summary, long help, example and a new section.
- Create fixtures in `xtask/diff/inputs/bigwig/`: `signal.bedgraph`, `signal.genome`, `signal.bw`, `peaks.bed`, `peaks.genome`, `peaks.bb`, `genes12.bed`, `genes12.bb`, `float32.bedgraph`, `float32.genome`, `float32.bw`, `regions.bed`
- Test: `crates/gloci/tests/bigwig.rs`, `crates/gloci-py/tests/test_outputs.py`; cases and goldens below
- Modify: `docs/differences.md` (`#bigwig-float32`), `docs/migration.md`, `docs/design.md` §5.2 and §7.4 (`convert`), `CHANGELOG.md`, `SKILL.md`, `__init__.pyi`

**Interfaces:**
- Consumes: `BbiWriter`, `BbiWriteOptions`, `BbiSizes`, `PushError`, `parse_track_output_spec` (Task 5); `bbi_chrom_sizes` (Task 3); `OrderGuard::contiguous` (`gloci-io/src/sorted.rs`).
- Produces:

```rust
/// The input formats each output takes: bigWig from BedGraph, Bed (4 numeric columns), Wig, BigWig;
/// bigBed from Bed, BigBed. Anything else: UsageError naming them.
pub fn check_track_input(format: Format, kind: BbiKind) -> Result<(), UsageError>;
/// Every record of `input` to `w`, in input order: OrderGuard::contiguous on each batch, then
/// push; a refused record is IoError::Parse { path: input name, loc: the record's (Line for text,
/// Record for binary input), msg }.
pub fn convert_to_track(input: Box<dyn RecordReader>, chroms: &mut ChromTable, w: &mut BbiWriter) -> Result<(), IoError>;
```

- CLI:
  - `ConvertArgs.autosql: Option<PathBuf>` (`--autosql FILE`): "bigBed output: the autoSql describing the columns, stored as given (bedToBigBed -as)".
  - `-g` help: "With --bam, or for a bigWig or bigBed output: …".
- `convert::run` parses `-o` with `parse_track_output_spec`. A bigWig/bigBed output:
  - with `--bam`, `--fastq`, `--bedpe`, any conversion option (`--blocks`, `--split-alleles`, …) or `--write-index`: exit 2;
  - `--autosql` without a bigBed output: exit 2;
  - `-g` without `--bam` or a bigWig/bigBed output: exit 2, saying where `-g` applies.
- Sizes: `-g` (`resolve_genome`, source `-g ARG`), else `bbi_chrom_sizes` of a bigWig/bigBed input (renamed by `--chrom-map`; source `the header of FILE`), else exit 2 naming `-g`.

- [ ] **Step 1: Make the fixtures.**
  - `signal.bedgraph`, by:

```sh
awk 'BEGIN{OFS="\t"; for(i=0;i<3000;i++) print "chr1", i*10, i*10+7, (i%40-20)/8;
  for(i=0;i<50;i++) print "chr10", i*100, i*100+100, i%7;
  for(i=0;i<40;i++) print "chr2", i*25, i*25+5, (i+1)/4}' > signal.bedgraph
```

    Its values are exact eighths and quarters, so awk's text is the shortest float text. `signal.genome` is `chr1 40000`, `chr10 6000`, `chr2 2000`.
  - `peaks.bed` (BED6+2, lexicographic): `chr1 100 200 p1 500 + 12.5 x`, `chr1 150 400 p2 0 - 3 y`, `chr1 1000 1001 p3 1000 . 0.001 z`, `chr2 5 50 p4 7 + 99 w`; `peaks.genome` is `chr1 5000`, `chr2 1000`.
  - `genes12.bed`: two valid BED12 records on `chr1` (genome `peaks.genome`).
  - `float32.bedgraph`: `chr1 0 10 0.123456789`, `chr1 10 20 16777217`, `chr1 20 30 1e-50`, `chr1 30 40 3.4028235e38`; `float32.genome` is `chr1 100`.
  - `regions.bed` (BED4, sorted): `chr1 0 100 r1`, `chr1 95 1000 r2`, `chr10 450 1250 r3`, `chr2 0 30 r4`, `chr2 1000 1500 r5` (no signal).
  - The `.bw`/`.bb` files come from this task's `gloci convert X -g X.genome -o X.bw|.bb` once Step 4 passes, and are committed.
- [ ] **Step 2: Write the failing tests** in `crates/gloci/tests/bigwig.rs`:
  - `convert_bedgraph_to_bigwig_and_back`: `convert signal.bedgraph -g signal.genome -o t.bw`, then `view t.bw`, is byte-equal to `signal.bedgraph`. So is `convert t.bw -o t.bedgraph`.
  - `convert_wig_to_bigwig`: `crates/gloci/tests/wig-golden/input.wig` to `.bw` with a genome (`chr1 1000`, `chr2 1000`), then `view`, equals `expected.bedgraph` sorted lexicographically.
  - `convert_bedgraph_detected_as_bed`: a 4-column file without a track line converts. A BED4 with a name in column 4 (`chr1 0 10 geneA`) exits 1, its stderr containing `line 1` and `not a number`. A BED6 exits 1, its stderr containing `4 columns`.
  - `convert_bigwig_to_bigwig_keeps_its_sizes`: `convert t.bw -o u.bw` without `-g` succeeds, and `bbi_chrom_sizes(u.bw)` lengths equal the genome's.
  - `convert_bed_to_bigbed_and_back`: `peaks.bed` and `genes12.bed` round-trip byte-equal through `.bb` and `view`.
  - `convert_autosql_option`: `--autosql my.as` (8 fields) is stored verbatim (`BigBedRead::autosql`). A 5-field file exits 1, its stderr containing `5` and `8`. `--autosql` with `-o x.bw` exits 2.
  - `convert_to_bigwig_needs_sizes`: `convert signal.bedgraph -o x.bw` exits 2, its stderr containing `-g`.
  - `convert_errors_name_the_input_line`:
    - a bedGraph whose third line overlaps the second exits 1, stderr containing `line 3` and `overlap`;
    - a bedGraph `chr1 …, chr2 …, chr1 …` exits 1, stderr containing `gloci sort` (OrderGuard);
    - a record on a chromosome missing from `-g` exits 1, stderr containing that chromosome and `-g`.
  - `convert_track_option_conflicts`: each exits 2:
    - `--blocks -o x.bb`;
    - `--bam -g hg38 -o x.bw`;
    - `--write-index -o x.bw`;
    - `-g hg38 -o x.bed` (`-g` without `--bam` or a track output);
    - a VCF to `-o x.bw` (stderr names bedGraph, WIG and bigWig);
    - a bedGraph with `type=bedGraph` to `-o x.bb`.
  - `other_commands_refuse_bigwig_output`: `sort signal.bedgraph -o y.bw` and `view peaks.bed -o y.bb` exit 2, stderr containing `gloci convert`.
  - `bigwig_output_must_be_a_file`: `convert signal.bedgraph -g signal.genome -O bigwig` (stdout) exits 2.
  - `convert_threads_do_not_change_the_file`: the same conversion at `-t 1`, `2` and `8` gives byte-identical files.
  - `failed_convert_leaves_nothing`: an error at line 3,000 of a long bedGraph leaves no `x.bw` and no `x.bw.tmp-*`.
  - Python `crates/gloci-py/tests/test_outputs.py::test_bigwig_save_and_open`:
    - `gloci.open(signal).convert(genome=str(genome)).save(tmp / "x.bw")`;
    - `gloci.open(str(tmp / "x.bw")).to_bytes()` equals the bedGraph's bytes;
    - `gloci.open(open(tmp / "x.bw", "rb")).to_bytes()` raises `gloci.UsageError`, with `bigWig` in the message.
- [ ] **Step 3: Run** `cargo test -p gloci --test bigwig convert`. Expected: FAIL.
- [ ] **Step 4: Implement** `check_track_input`, `convert_to_track` and `run_track` in `commands/convert.rs`.
  - `run_track` opens the input, checks it, builds the sizes, creates the writer (threads `ctx.threads`), converts, then calls `finish` with `with_output_name`.
  - A `PushError::Writer(e)` returns `e`.
  - Help:
    - summary `Convert to BED/FASTQ/BAM/bigWig/bigBed (convert2bed, bedtobam)` (62 characters);
    - the first paragraph adds "bedGraph, WIG or bigWig to bigWig and BED to bigBed (bedGraphToBigWig, wigToBigWig, bedToBigBed), bigWig and bigBed to bedGraph and BED (bigWigToBedGraph, bigBedToBed)";
    - the example `gloci convert reads.byname.bam --bedpe -o pairs.bedpe` becomes `gloci convert cov.bedgraph -g hg38 -o cov.bw`;
    - a new section `bigWig and bigBed:` before `Order and output:`, saying D171's write rules in at most 8 lines.
- [ ] **Step 5: Run** the tests of Step 2. Expected: PASS. Run the Python test with `maturin develop -m crates/gloci-py/Cargo.toml --features test-hooks && pytest crates/gloci-py/tests -k bigwig` where a virtualenv has maturin; otherwise CI's python job runs it.
- [ ] **Step 6: Golden cases**, then `cargo xtask diff --image gloci-diff:lane-a --only NAME` and `--update-golden` for each:

```toml
[[case]]
name = "bigwig-written-by-gloci-vs-pybigwig"
input = { file = "xtask/diff/inputs/bigwig/signal.bw" }
aux = ["xtask/diff/inputs/bigwig/bbi_dump.py"]
# gloci's bigWig (3,000 records on chr1: several sections and zoom levels) read back by pyBigWig.
gloci = "gloci view input.bw"
reference = "python3 -I bbi_dump.py --lex input.bw"
compare = "exact"

[[case]]
name = "bigbed-written-by-gloci-vs-pybigwig"
input = { file = "xtask/diff/inputs/bigwig/peaks.bb" }
aux = ["xtask/diff/inputs/bigwig/bbi_dump.py"]
gloci = "gloci view input.bb"
reference = "python3 -I bbi_dump.py --lex input.bb"
compare = "exact"

[[case]]
name = "bigbed12-written-by-gloci-vs-pybigwig"
input = { file = "xtask/diff/inputs/bigwig/genes12.bb" }
aux = ["xtask/diff/inputs/bigwig/bbi_dump.py"]
gloci = "gloci view input.bb"
reference = "python3 -I bbi_dump.py --lex input.bb"
compare = "exact"

[[case]]
name = "bigwig-float32-vs-input"
input = { file = "xtask/diff/inputs/bigwig/float32.bw" }
aux = ["xtask/diff/inputs/bigwig/float32.bedgraph"]
# bigWig holds 32-bit floats: the values read back are the nearest ones, not the text written.
gloci = "gloci view input.bw"
reference = "cat float32.bedgraph"
compare = "exact"
expected_difference = "#bigwig-float32"
```

  - Round-trip tests in `bigwig.rs`: `golden_round_trips` asserts that each golden's `expected` equals its source text (`signal.bedgraph`, `peaks.bed`, `genes12.bed`).
  - The float32 golden's `expected_gloci` is `0.12345679`, `16777216`, `0`, `340282350000000000000000000000000000000` in column 4.
- [ ] **Step 7: Docs.**
  - `docs/differences.md`: append `## bigWig stores 32-bit floats {#bigwig-float32}`.
    - **gloci:** a bigWig value is the nearest 32-bit float to the number written, printed back as the shortest decimal that reads as that float (`0.123456789` → `0.12345679`; `16777217` → `16777216`; `1e-50`, too small for a 32-bit float, → `0`). A value whose nearest float is infinite (`1e39`) or that is not finite is an error naming the line.
    - **Old behaviour:** the bedGraph keeps the text as written.
    - **Reason:** the bigWig/bigBed format paper (doi:10.1093/bioinformatics/btq351) stores values as 32-bit floats, so every writer rounds. Refusing inexact values would refuse most real signal files.
    - **Pinned by:** `bigwig-float32-vs-input` and `value_outside_float32_is_an_error`.
  - `docs/migration.md` (bigWig subsection), rows:
    - `bedGraphToBigWig in.bedGraph chrom.sizes out.bw` → `gloci convert in.bedgraph -g chrom.sizes -o out.bw` (any chromosome order with each chromosome contiguous; a 4-column BED without a track line is a bedGraph);
    - `wigToBigWig in.wig chrom.sizes out.bw` → `gloci convert in.wig -g chrom.sizes -o out.bw`;
    - `bedToBigBed [-as=x.as] in.bed chrom.sizes out.bb` → `gloci convert in.bed -g chrom.sizes -o out.bb [--autosql x.as]` (bedN+M is found from the first record; no `-type`).
  - `docs/design.md`:
    - §5.2: add `.bw`, `.bigwig`, `.bb`, `.bigbed` to the format bullet, "convert and genomecov only; a regular file only (D171)";
    - §7.4: one `convert` bullet for bigWig/bigBed.
  - `CHANGELOG.md`: one `### Added` bullet.
  - Regenerate SKILL.md and the stub (`convert()` gains `autosql=`).
- [ ] **Step 8: Run** `cargo test -p gloci --test bigwig --test convert --test convert_bam --test golden --test cli_basics --test embed_args && cargo test -p gloci-ops && cargo test -p xtask --test pystub --test docs --test references --test book`, clippy and fmt. Expected: PASS.
- [ ] **Step 9: Commit** `feat(convert): bedGraph, WIG and bigWig to bigWig, BED to bigBed, and back; checked against pyBigWig (D171)`

### Task 7: `gloci genomecov --bedgraph/--bga -o cov.bw`

**Files:**
- Modify: `crates/gloci/src/commands/genomecov.rs`, `crates/gloci/src/cli/intervals.rs` (`genomecov` help)
- Test: `crates/gloci/tests/bigwig.rs`; live case `genomecov-bga-bigwig-vs-bedtools`
- Modify: `docs/migration.md`, `CHANGELOG.md`, `SKILL.md`, `__init__.pyi`

**Interfaces:**
- Consumes: `BbiWriter` as `RecordWriter`, `BbiSizes`, `parse_track_output_spec` (Task 5); `genome_of` (already in `genomecov.rs`), whose source string is `-g ARG` or `the header of FILE`.
- Produces: `genomecov -o x.bw` (or `-O bigwig` with `-o`) in `GenomecovMode::BedGraph | BedGraphAll` writes a bigWig.
  - Any other mode with a bigWig output: exit 2, naming `--bedgraph` and `--bga`.
  - A bigBed output: exit 2 (`bbi_output_refusal(BigBed)`).
  - `--write-index`: exit 2.

- [ ] **Step 1: Write the failing tests** in `crates/gloci/tests/bigwig.rs`:
  - `genomecov_bga_to_bigwig_matches_bedgraph`: for `xtask/diff/inputs/gc-a.bed` with `-g xtask/diff/inputs/gc.genome`, `view` of the `-o g.bw` output equals the `--bga` stdout sorted lexicographically (`LC_ALL=C sort -k1,1 -k2,2n`). The same holds for `xtask/diff/inputs/aln/cigars-bases.bam` (sizes from its header, no `-g`).
  - `genomecov_bedgraph_to_bigwig`: the same for `--bedgraph` and for `--bedgraph --scale 0.5`.
  - `genomecov_bigwig_needs_a_bedgraph_mode`: each exits 2:
    - the default `--hist` with `-o x.bw` (stderr contains `--bedgraph` and `--bga`);
    - `--per-base -o x.bw`;
    - `--bga -o x.bb`;
    - `--bga -o x.bw --write-index`.
  - `genomecov_bigwig_threads`: an indexed BAM (an existing fixture with its `.bai`) at `-t 1` and `-t 8` gives byte-identical `x.bw`.
  - **Review Focus 1**, `genomecov_empty_bedgraph_to_bigwig_is_an_error`: an empty BED with `-g gc.genome --bedgraph -o x.bw` exits 1, its stderr containing `no records`, and leaves no `x.bw`.
- [ ] **Step 2: Run** `cargo test -p gloci --test bigwig genomecov`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - `genomecov::run` parses `-o` with `parse_track_output_spec`. For a bigWig output it checks the mode, then creates `BbiWriter::create(&spec, BbiWriteOptions { kind: BigWig, sizes: BbiSizes { genome: genome.clone(), source }, threads: ctx.threads, autosql: None })` and passes it as `w`. Everything else in the function stays.
  - Help:
    - the example `gloci genomecov reads.bam --bga --split` becomes `gloci genomecov reads.bam --bga --split -o cov.bw`;
    - add to the "Output:" paragraph: "-o FILE.bw writes --bedgraph or --bga as bigWig, with the chromosome sizes of -g or of the BAM header (D171)."
- [ ] **Step 4: Run** the tests of Step 1. Expected: PASS.
- [ ] **Step 5: Live case** (a pipeline, so not a golden):

```toml
[[case]]
name = "genomecov-bga-bigwig-vs-bedtools"
input = { file = "xtask/diff/inputs/gc-a.bed" }
aux = ["xtask/diff/inputs/gc.genome", "xtask/diff/inputs/bigwig/bbi_dump.py"]
# gloci writes bigWig; pyBigWig reads it back in the file's order, which is the genome's.
gloci = "gloci genomecov input.bed -g gc.genome --bga -o g.bw && python3 -I bbi_dump.py g.bw"
reference = "bedtools genomecov -i input.bed -g gc.genome -bga"
compare = "exact"
```

Run `cargo xtask diff --image gloci-diff:lane-a --only genomecov-bga-bigwig-vs-bedtools`. Expected: `ok (matches)`.
- [ ] **Step 6: Docs.**
  - `docs/migration.md` (bigWig subsection): add a row `bedtools genomecov -ibam x.bam -bga | bedGraphToBigWig - chrom.sizes cov.bw` → `gloci genomecov x.bam --bga -o cov.bw`.
  - `CHANGELOG.md`: one bullet.
  - Regenerate SKILL.md and the stub.
- [ ] **Step 7: Run** `cargo test -p gloci --test bigwig --test genomecov --test genomecov_alignments --test cli_basics`, the xtask checks, clippy and fmt. Expected: PASS.
- [ ] **Step 8: Commit** `feat(genomecov): --bedgraph and --bga write bigWig to -o FILE.bw (D171)`

### Task 8: `map`'s `wsum`, and `bigWigAverageOverBed`'s columns

**Files:**
- Modify:
  - `crates/gloci-ops/src/agg/spec.rs`: `Op::WSum`, `"wsum"`, `needs_second_input`, `needs_overlap`, the three op lists in messages;
  - `crates/gloci-ops/src/agg/ops.rs`: state and finish;
  - `crates/gloci/src/cli/relate.rs`: `map` help;
  - `crates/gloci/src/cli/intervals.rs`: the `groupby` doc's list of two-input ops.
- Test: `crates/gloci-ops/tests/agg.rs`, `crates/gloci/tests/bigwig.rs`; case and golden `map-bigwig-wsum-vs-bedtools-per-base`
- Modify: `docs/design.md` §6.3, `docs/migration.md` (`map` row ops list; bigWig subsection), `CHANGELOG.md`, `SKILL.md`, `__init__.pyi`

**Interfaces:**
- Produces: `Op::WSum` (`"wsum"`), a `State::Numeric` op whose output is `n.wxsum` through `fmt_f64`, or `--empty` when no B record counted (`n.n == 0`). It needs a second input and the overlap. Column 0 is refused as for `wmean`.

- [ ] **Step 1: Write the failing tests.**
  - `agg.rs::wsum_multiplies_values_by_overlapping_bases`, with `run_hits` as `wmean`'s test uses it:
    - hits `("10", 10 bases)` and `("20", 30 bases)` give `"700"`;
    - `("0.5", 3 bases)` gives `"1.5"`;
    - two hits with 0 shared bases give `"0"`;
    - no hits give `"."`.
  - Add `"wsum"` to `two_input_ops_are_refused_without_a_second_input`.
  - `parse_aggs(["0:wsum"])` errors, mentioning `column 0`.
  - `bigwig.rs::map_over_bigwig_gives_bigwigaverageoverbed_columns`: B is `r.bw` (`chr1 0 100 1.5`, `chr1 100 200 0.25`, size 1000), A is `a.bed` (`chr1 50 150 x` and `chr1 500 600 y`). `map -a a.bed -b r.bw -c 4:wsum -c 0:bases_uniq -c 4:wmean --empty 0` prints exactly `chr1\t50\t150\tx\t87.5\t100\t0.875\nchr1\t500\t600\ty\t0\t0\t0\n` (1.5×50 + 0.25×50 = 87.5 over 100 bases).
- [ ] **Step 2: Run** `cargo test -p gloci-ops --test agg wsum && cargo test -p gloci --test bigwig map_over`. Expected: FAIL (unknown op `wsum`).
- [ ] **Step 3: Implement.**
  - Add the op and its finish arm.
  - In `map`'s help "Operations:" paragraph, insert `wsum` after `wmean` in the op list, and add: "wsum is the sum of value times overlapping bases (over a bigWig, -c 4:wsum -c 0:bases_uniq -c 4:wmean --empty 0 gives bigWigAverageOverBed's sum, covered and mean)".
  - Add `wsum` to the op lists in `spec.rs`'s three messages and to the `groupby` doc.
- [ ] **Step 4: Run** the tests of Step 1. Expected: PASS.
- [ ] **Step 5: Golden case** (bedtools over a per-base expansion of the same signal):

```toml
[[case]]
name = "map-bigwig-wsum-vs-bedtools-per-base"
input = { file = "xtask/diff/inputs/bigwig/signal.bw" }
aux = ["xtask/diff/inputs/bigwig/signal.bedgraph", "xtask/diff/inputs/bigwig/regions.bed"]
# One base per row makes bedtools' sum, count and mean the weighted sum, covered bases and weighted mean.
# bedtools map has no -null, so its `.` for a region without signal becomes 0, as --empty 0 writes.
gloci = "gloci map -a regions.bed -b input.bw -c 4:wsum -c 0:bases_uniq -c 4:wmean --empty 0"
reference = '''awk 'BEGIN{OFS="\t"}{for(i=$2;i<$3;i++) print $1,i,i+1,$4}' signal.bedgraph > pb.bed && bedtools map -a regions.bed -b pb.bed -c 4,4,4 -o sum,count,mean -prec 17 | awk 'BEGIN{OFS="\t"}{for(i=5;i<=NF;i++) if($i==".") $i=0; print}' '''
compare = "numeric"
```

Run it with `cargo xtask diff --image gloci-diff:lane-a --only map-bigwig-wsum-vs-bedtools-per-base`, then with `--update-golden`.
- [ ] **Step 6: Docs.**
  - `docs/design.md` §6.3: add `wsum` (value × overlapping bases) beside `wmean`.
  - `docs/migration.md`:
    - the `map` row's op list gains `wsum`;
    - the bigWig subsection gains a row `bigWigAverageOverBed in.bw in.bed out.tab` → `gloci map -a in.bed -b in.bw -c 4:wsum -c 0:bases_uniq -c 4:wmean --empty 0` (columns sum, covered and mean after A's; size is `end - start`, mean0 is `sum / size`; A keeps its order and all its columns).
  - `CHANGELOG.md`: one bullet.
  - Regenerate SKILL.md and the stub.
- [ ] **Step 7: Run** `cargo test -p gloci-ops && cargo test -p gloci --test map --test groupby --test bigwig --test golden --test cli_basics`, the xtask checks, clippy and fmt. Expected: PASS.
- [ ] **Step 8: Commit** `feat(map): wsum, value times overlapping bases; bigWigAverageOverBed's columns from map (D171)`

### Task 9: The benchmark image, real-data checks, the nightly guard (light review)

**Files:**
- Modify: `xtask/docker/Dockerfile.bench`, `xtask/bench/datasets.toml`, `paper/fetch-datasets.sh` (regenerated), `xtask/tests/paper_datasets.rs`, `crates/gloci/tests/large.rs`
- Modify: `docs/decisions.md` (D171: a memory bullet), `CHANGELOG.md` (`### Changed`)

**Interfaces:**
- Consumes: `gloci convert` of bigWig, bedGraph and BED (Tasks 3, 6).
- Produces: the derived input `ctcf-fc.bedgraph`, built by gloci; no `ucsc-` text in either Dockerfile.

- [ ] **Step 1: Write the failing tests.**
  - In `xtask/tests/paper_datasets.rs::bench_and_diff_images_pin_same_reference_tools`, the tool list becomes `["gawk=", "time="]`.
  - New test `no_image_installs_ucsc_tools`: neither `xtask/docker/Dockerfile` nor `xtask/docker/Dockerfile.bench` contains `ucsc-` or `bigWigToBedGraph`.
  - `datasets_build_the_ctcf_bedgraph_with_gloci`: the derived `ctcf-fc.bedgraph` has `from` containing `repo:target/release/gloci` and a command starting `{in2} convert {in}`.
- [ ] **Step 2: Run** `cargo test -p xtask --test paper_datasets`. Expected: FAIL.
- [ ] **Step 3: Implement.**
  - `Dockerfile.bench`: drop the `ucsc-bigwigtobedgraph=469 \` line, and drop "bigWigToBedGraph prints no version;" from the pin-check comment.
  - `datasets.toml` `ctcf-fc.bedgraph`: `from = ["encode-ctcf-k562-fc", "repo:target/release/gloci"]`, `command = "{in2} convert {in} | sort -k1,1 -k2,2n -s > {out}"`, and a comment "gloci reads the ENCODE bigWig itself (D171); build it first with cargo build --release". Empty the `sha256` of `ctcf-fc.bedgraph` and `ctcf-fc.bed5` (both change if the printed values do).
  - Regenerate the script with `cargo xtask paper fetch-script`.
  - The bench image is not rebuilt here; it is rebuilt the next time benchmarks run.
- [ ] **Step 4: Record the digests and check against the real data.** Run `cargo build --release -p gloci && cargo xtask paper fetch --only ctcf-fc.bed5 --record --no-docker`. This downloads the 503 MB ENCODE bigWig and the GENCODE GTF and records both sha256 values in `datasets.toml`. Then, with `C=target/bench-data/cache` and `T` a scratch directory:
  - `/usr/bin/time -v target/release/gloci convert $C/ENCFF979PWH.bigWig -o $T/ctcf.bw` (reads ENCODE's file, writes a bigWig with its sizes);
  - `/usr/bin/time -v target/release/gloci convert $C/ctcf-fc.bedgraph -g hg38 -o $T/ctcf2.bw`; if `hg38` lacks a chromosome the error names it, so record that and skip this run;
  - `/usr/bin/time -v target/release/gloci convert $C/exons.bed -g hg38 -o $T/exons.bb` (GENCODE exons, BED6) (re-run `paper fetch --only exons.bed --no-docker` if absent);
  - `target/release/gloci view $T/ctcf.bw | cmp - $C/ctcf-fc.bedgraph` (expected: no output) and `target/release/gloci view $T/exons.bb | cmp - $C/exons.bed` (expected: no output).

  Note each run's "Maximum resident set size" and elapsed time.
- [ ] **Step 5: Nightly guard.** In `crates/gloci/tests/large.rs`, `#[ignore] fn track_writes_stay_flat`:
  - `convert big_a.bed -g <chr1..chr22 in the order lex_chroms(22) gives, each 4,000,000,000 long> -o x.bb` (the shared 2 GiB BED6);
  - `genomecov big_a.bed -g <same> --bga -o x.bw`;
  - each peak RSS (`wait_rss`) at most the larger of 256 MiB and twice the largest RSS of Step 4, printed as the other large tests print theirs.

  Check that it compiles (`cargo test -p gloci --test large --no-run`); it runs nightly, not here.
- [ ] **Step 6: Docs.**
  - D171 gains a `<li>` with the three peak RSS figures and times from Step 4, and the round-trip results. The figures must be flat: below 256 MiB on 0.5 GB of input.
  - `CHANGELOG.md` `### Changed`: "The benchmark image no longer installs UCSC's bigWigToBedGraph; gloci converts the ENCODE bigWig itself."
- [ ] **Step 7: Run** `cargo test -p xtask --test paper_datasets --test references`, `cargo test --workspace --no-fail-fast`, clippy, fmt, and the licence check of the Global Constraints. Expected: PASS, and the licence check prints nothing.
- [ ] **Step 8: Clean up** with `rm -rf $T target/bench-data/cache/{ENCFF979PWH.bigWig,ctcf-fc.bedgraph,ctcf-fc.bed5,gencode.v46.annotation.gtf.gz,exons*.bed*}`. Remove the `bbi-deps` scratch project if one was made, and run `docker image rm gloci-diff:lane-a` once no later lane step needs it.
- [ ] **Step 9: Commit** `chore(bench): gloci builds the CTCF bedGraph from ENCODE's bigWig; no UCSC tool in any image (D171)`
