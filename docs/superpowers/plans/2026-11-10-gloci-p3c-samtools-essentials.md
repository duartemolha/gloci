# gloci P3c samtools Essentials Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cover the samtools commands a genomics team runs every day, on top of P3a's BAM/CRAM/FASTA support:
- `idxstats`, `flagstat`, `flags`, `head`, `samples`;
- `quickcheck` (as `check --quick`) and full `check` of alignment files;
- `reheader` (BAM block copy; CRAM header container);
- `faidx` region extraction and `dict`;
- `view --count`, and `view --regions-file` (samtools `-L`) on alignments;
- `depth` and `bedcov`;
- a samtools section in `docs/migration.md`.

**Architecture:** No new crates.
- Per-record commands read through P3a's `AlnSource` (D120): every record, in file order, unmapped included.
- `depth` and `bedcov` share one per-base counting core (`gloci-ops/src/depth/`), a ring buffer bounded by the longest read span. genomecov's event sweep does not fit: samtools counts each base against its own quality, prints zero rows inside read spans, and clips overlapping mates per base.
- `faidx` and `getfasta` share P3a's sequence path: `IndexedFasta::fetch`/`fetch_rev` and `seq::complement`. `faidx` adds only samtools' region strings, headers and line wrapping.
- `reheader` reuses P3a's `bgzf_copy::copy_from` (Task 15) for BAM. For CRAM it rewrites the file definition and the header container, then copies every data container byte for byte.

**Tech Stack:** Everything P3a uses, with no new dependency: `md-5` 0.10 is already a workspace dependency and is added to gloci-ops for `dict`. CRC32 comes from `flate2::Crc`, which gloci-io already depends on. Reference tool: samtools 1.21 (and bcftools 1.21 for `head` on VCF), in `gloci-diff:p0` since P3a Task 1.

**Spec:** `docs/superpowers/specs/2026-11-06-gloci-design.md`, covering:
- §10, the P3c row and the Order note;
- §4.1 (region strings), §4.3 (overlap), §4.5 (order and sortedness);
- §5.2 (atomic output);
- §6.1 (global flags), §6.4 (exit codes, errors, determinism);
- §7.4 (`check`, `view`, `index`);
- §8.1–§8.3;
- §9 (performance, differential tests, threads, large files);
- §11 (migration).

**Builds on:** P3a, fully merged: Tasks 1–20 and D99–D124, as amended on 2026-11-10. In particular:
- D101 and D116: interval commands skip unmapped reads, and record-preserving commands keep them.
- D120: `AlnSource`, `open_aln_source`, `sam_line_to_bam`.
- D102: the read filters. D107: the alignment indexes. D108: FASTA. D109: the CRAM reference.
- D118: the `bam` compare mode and `@PG` normalisation.
- D123 and Task 15: `bgzf_copy::copy_from`.

**Out of scope:**
- The spec's list: `markdup`, `fixmate`, `calmd`, `mpileup`, `consensus`, `stats`, `ampliconclip`.
- Everything else D125 lists as not planned.
- P3b; Python.

**P3c is done when:**
- every command and flag in D125–D140 has golden cases against samtools 1.21;
- every intentional difference is anchored in `docs/differences.md`, pinned by a case and linked from the new samtools table in `docs/migration.md`;
- the Task 13 benchmarks show gloci no slower than samtools on each command it replaces (§2 anti-goal).

## Global Constraints

P3a's Global Constraints all apply:
- u64 half-open coordinates; `rust-version = "1.91"`; no htslib;
- exit codes 0, 1 (data or runtime) and 2 (usage), with errors that say what, where and how to fix;
- atomic output; content detection; record streams in gloci-ops;
- 2–4 `--help` examples per command;
- CI on three platforms;
- goldens committed under `crates/gloci/tests/golden/` and replayed by plain `cargo test`;
- output byte-identical at `-t 1/2/8`, with a test for every new path;
- nothing quadratic in the contig count.

New in P3c, binding on every task:
- **Parity.** Output is byte-identical to samtools 1.21 unless the difference is anchored in `docs/differences.md` and pinned by a case (`expected_difference`). The task that pins a difference writes its section.
- **Goldens discriminate.** Each task's fixtures (D141) are chosen so that a plausible wrong rule gives different output; the task lists the trap each fixture sets.
- **Memory bounds,** each asserted in its task:

  | Path | Bound |
  |---|---|
  | `idxstats` through an index | the header + the index's per-reference counts |
  | `idxstats`/`flagstat`/`view --count` scans | one record + P3a's BGZF/CRAM pool |
  | `check --quick` | the header + the first record + 28 bytes (BGZF), or the header + the first container + 38 bytes (CRAM) |
  | `reheader` | the old and new headers + one BGZF block (64 KiB) |
  | `head` | the header + one batch of up to 16384 records (text formats; a fixed record batch) |
  | `faidx` | ≤ 1 MiB of sequence + one output line |
  | `dict` | one FASTA line + the MD5 state |
  | `depth` | per input, `4 B × next_pow2(longest span + 1)`, plus pending mates (∝ depth) |
  | `bedcov` | per worker, `depth`'s bound for one region + one batch of ≤ 4 MiB of BED lines |

- **Huge headers.** 200,000-`@SQ` tests guard `idxstats`, `reheader`, `samples --match-reference` and `depth --all-references` (Tasks 2, 5, 6, 11). Name lookups use `HashMap`s built once.
- **Exit codes** follow §6.4 even where samtools differs (D126, `#samtools-exit-codes`).

## Decisions this plan adds

| # | Decision |
|---|---|
| D125 | **Scope and names.** <ul><li>`idxstats`, `flagstat`, `reheader`, `faidx`, `dict`, `head`, `samples`, `flags`, `depth` and `bedcov` keep samtools' names, as `cat` does (D123): gloci has no command with those meanings.</li><li>`quickcheck` is `check --quick` (D129), because `check` already means "validate files" (§7.4).</li><li>`view -c` is `view --count`, because `-c` is `COL:OP` (§6.3).</li><li>`view -L` is the existing `--regions-file` (§6.1), which already reads through an index (D133).</li><li>samtools' short flags are not copied, because several letters already mean something else in gloci (`-c -n -i -s -a -b -r -f`). Every samtools flag maps to a long flag named in its D-row and in the migration table.</li><li>The D100 gate adds `idxstats, flagstat, check, reheader, head, samples, depth, bedcov` to the alignment readers, and `faidx, dict` to the FASTA readers.</li><li>**Not planned** (each a migration row):<ul><li>`fqidx` and `faidx --fastq`;</li><li>custom index paths (`-X`, `--fai-idx`, `--gzi-idx`);</li><li>file lists (`depth -f`, path lists on stdin);</li><li>`reheader -c` (pipe `view --header-only` instead) and `--in-place` (`-o` set to the input replaces it atomically);</li><li>`depth --incl-flags` and `-d`/`-m` (ignored by samtools);</li><li>`bedcov --max-depth` (samtools' default is no cap);</li><li>`samples -F` (repeat `--match-reference`).</li></ul></li></ul> |
| D126 | **Shared rules.** <ul><li>Per-record commands read through `AlnSource`. Only `view`, `depth` and `bedcov` take the D102 filters (and `--chrom-map`). `idxstats`, `flagstat`, `check`, `head`, `samples` and `reheader` describe the whole file and take neither.</li><li>`--column-names` writes one `#`-prefixed column-name line first (samtools `depth -H`, `bedcov -H`, `samples -h`).</li><li>Exit status is §6.4's: 1 for data errors and failed checks, 2 for usage errors. samtools exits 1 for usage errors, 2 for `bedcov` errors, and a bit mask from `quickcheck`: `#samtools-exit-codes`.</li><li>**CRAM.** Commands that decode records need `--reference` whenever a container needs reference bases (D109): the `idxstats` and `flagstat` scans, `head`, `view`, `depth` and `bedcov`. samtools decodes only the fields it needs, so its `flagstat`, `idxstats` and `view -c` read such a file without one: `#cram-counts-need-reference`. `samples`, `reheader` and `check --quick` never need a reference.</li><li>`depth` and `bedcov` inputs must have identical `@SQ` lists (names, `LN`, order), else exit 1 naming the first difference. samtools compares reference ids across files unchecked: `#depth-bedcov-same-sq`.</li></ul> |
| D127 | **`gloci idxstats FILE`.** <ul><li>Per `@SQ`: `name TAB LN TAB mapped TAB unmapped`, then `* TAB 0 TAB 0 TAB unplaced`.</li><li>A BAM with a `.bai`/`.csi` is answered from the index's metadata pseudo-bins and `n_no_coor`, without reading records.</li><li>Otherwise one warning, and every record is counted through `AlnSource`. This covers a BAM without an index (samtools also scans), a SAM, and a CRAM (a `.crai` holds no counts). The warning reads: `x.bam has no index; counting every record (index it for speed: gloci index x.bam)`.</li><li>A placed unmapped read counts as unmapped on its reference (D107).</li><li>The scan does not need grouped references. samtools' scan stops with `file is not position sorted`: `#idxstats-unsorted`.</li><li>An index lacking `n_no_coor`, or with a reference that has bins but no metadata pseudo-bin, counts as "no counts": one warning and a scan. htslib prints `0 0` for it: `#idxstats-index-without-counts`.</li><li>Stale indexes warn (D26).</li></ul> |
| D128 | **`gloci flagstat FILE [-O text\|tsv\|json] [-o PATH]`.** <ul><li>The counts follow samtools' `flagstat_loop` exactly, with a QC-passed column and a QC-failed column (`0x200`) (Task 3).</li><li>`-O` keeps its §6.1 meaning, output format; `-o x.json`/`x.tsv` infers it; the default is `text`.</li><li>The text, TSV and JSON layouts are samtools' byte for byte.</li><li>Percentages are computed as samtools computes them: `(n as f32 / total as f32) as f64 * 100.0`, printed `{:.2}`. With a total of 0 they print `N/A` (`null` in JSON).</li><li>Every record counts, unmapped and unplaced included.</li></ul> |
| D129 | **`check` reads alignments; `check --quick` is quickcheck.** <ul><li>**Full `check`** of a BAM, SAM or CRAM reads every record through `AlnSource`.<ul><li>Errors: records P3a rejects (Task 2's list); a missing BGZF EOF marker; positions out of order when `@HD SO:coordinate` is declared (a decrease in `(ref id, pos)`, or a placed read after an unplaced one).</li><li>Warnings: a span past `LN`; identical records at one position (the P0 duplicate rule).</li><li>Under any other `SO`, order is reported but is not an error.</li><li>The summary adds mapped, unmapped-placed and unplaced counts, and the `SO`.</li></ul></li><li>**`--quick`** works on any format `check` reads. It checks:<ul><li>that the header parses;</li><li>for alignments, at least one `@SQ` unless `--unaligned` (samtools `-u`);</li><li>that the first record parses, or, for CRAM, the first container's header (so no reference is needed);</li><li>the BGZF EOF marker, or the CRAM EOF container, when the input is a regular file. stdin is not checked, as htslib.</li></ul></li><li>`--list-failed` prints only the failing paths, one per line, on stdout (samtools `-v`; gloci's `-v` is verbose logging). Without it, the usual summary is printed. `-q` silences warnings.</li><li>`--chrom-map` is a usage error (exit 2) for BAM, CRAM and SAM; `-g` and `--chrom-order` warn once and are ignored, as the header gives the references and the order.</li><li>The `--quick` read is bounded: header + first record + 28 bytes (BGZF), or header + first container + 38 bytes (CRAM).</li><li>CRAM 2.x is refused by name and version wherever the prelude is read: `#cram-limits`.</li><li>Any failing file → exit 1.</li><li>Pinned under `#quickcheck`: exit 1, not samtools' bit mask; VCF, BCF and BGZF BED pass (samtools: not sequence data); a corrupt first record fails (samtools never reads records).</li></ul> |
| D130 | **`gloci reheader HEADER IN [-o OUT]`.** <ul><li>HEADER is SAM header text, or a SAM/BAM/CRAM whose header is taken (as samtools); `-` reads stdin.</li><li>IN is a BAM or CRAM regular file. stdin is a usage error suggesting `gloci view - -o x.bam` first. `-o IN` replaces IN atomically.</li><li>**`@SQ` check.**<ul><li>The new list must keep the input's references in order with equal `LN`. Names may change, and new `@SQ` lines may follow the last old one, because records point at references by position.</li><li>For CRAM, an `M5` present in both lists must be equal.</li><li>Otherwise exit 1, naming both lines: `@SQ 3 of new.sam (SN:chr3 LN:1000) does not match @SQ 3 of x.bam (SN:chr3 LN:198295559); reheader keeps each read's reference by its position in the @SQ list, so change only SN: to rename, or use gloci view --chrom-map`.</li><li>samtools checks nothing: `#reheader-sq-check`.</li></ul></li><li>`@PG` follows D104 (`with_program`).</li><li>**BAM.**<ul><li>The input's EOF marker is checked first (`open_input`'s rule).</li><li>The new header is written as BGZF and flushed to a block edge.</li><li>Then `copy_from(input, out, extent)` runs, where extent is the old header's byte length in the decompressed stream.</li><li>One EOF marker ends the output.</li></ul></li><li>**CRAM.** The same file definition (version, file id), then a new header container: one raw `FILE_HEADER` block holding `i32 length + text`, with container and block CRC32 on 3.x. Then every byte from the old first data container to the end, verbatim. noodles' header writer is not used: it re-serialises the text and `expect`s a reference for every `@SQ` without `M5`.</li></ul> |
| D131 | **`gloci faidx REF REGION… [--regions-file BED] [-o OUT]`.** <ul><li>REF needs `.fai` (+ `.gzi`) (D108). Without one: exit 2 with `gloci index REF`. samtools builds it: `#fasta-implicit-index`.</li><li>With no regions: exit 2 naming `gloci index REF` (samtools then builds the index).</li><li>Regions use `Region::parse` (D8, htslib's rules). D8 gains htslib's `chr:-E` = `chr:1-E`. A region that does not parse is a usage error.</li><li>Each region prints `>` + the region text as typed + the strand suffix, then the bases wrapped at `--line-width N`. The default is the sequence's `.fai` line width; `0` means one line. A BED region prints as `chrom:start+1-end`; a zero-length BED region is skipped with one count warning.</li><li>An end past the sequence is cut there with a warning (samtools `Truncated sequence`). A start past it prints the header and no bases (samtools `Zero length sequence`).</li><li>An unknown sequence is a data error before any output. `--skip-invalid` skips the region entirely; samtools `-c` still prints its `>` line: `#faidx-missing-sequence`.</li><li>`--reverse-complement` (samtools `-i`) uses `fetch_rev` + `seq::complement`, the table samtools' `comp_base` uses.</li><li>`--mark-strand rc\|no\|sign\|custom,POS,NEG`: default `rc`, meaning `/rc` on reverse-complemented regions only. `sign` gives `(+)`/`(-)`.</li><li>`-o x.fa.gz` is BGZF, and `--write-index` writes `.fai`/`.gzi` (P3a Task 11's FASTA output path).</li></ul> |
| D132 | **`gloci dict REF [-o OUT]`** streams the FASTA (plain, gzip or BGZF; no index). <ul><li>It writes `@HD VN:1.0 SO:unsorted` unless `--no-header`.</li><li>Then, per sequence: `@SQ SN:<first word> LN:<n> M5:<md5>`, followed by `AH:*` (`--alt FILE`), `AN:` (`--alias`), `UR:`, `AS:` (`--assembly`) and `SP:` (`--species`), in samtools' order.</li><li>`LN` and `M5` cover only bytes `!`–`~`, with `a`–`z` uppercased. IUPAC codes, `*` and `-` count; spaces, tabs and `\r` do not.</li><li>`UR` defaults to `file://` + the canonical absolute path (none for stdin); `--uri` overrides it.</li><li>`--alias` follows samtools: strip or add `chr`; `M` adds `,chrMT,MT` and `MT` adds `,chrM,M`.</li><li>A duplicate name is a data error naming both lines. samtools writes both: `#dict-duplicate-names`.</li><li>Lengths are u64.</li></ul> |
| D133 | **`view --count` and alignment regions.** <ul><li>`--count` prints the number of records `view` would write (regions and D102 filters applied), as `N\n`, on stdout or `-o`. It works for every format `view` reads. It conflicts with `--header-only`, `-O` and `--write-index`.</li><li>`--regions-file BED` is samtools `-L`. Reads match a region under §4.3 on the D101 span (D116). A zero-length region matches reads under §4.3 (a point inside the span). It reads through the index when there is one; samtools `-L` scans.</li><li>Each record is written or counted once, however many regions it overlaps. samtools repeats it per positional region: `#view-overlapping-regions`.</li><li>Positional regions and `--regions-file` combine as a union (P0), and overlapping regions write each record once. samtools intersects them: `#view-regions-union`.</li></ul> |
| D134 | **`gloci head FILE [--headers N] [--records N]`** (samtools and bcftools `head`): all header lines and no records by default. <ul><li>Alignments print SAM text: the first N lines of the header text, then records through `open_alignment_writer` in its no-header SAM mode.</li><li>VCF/BCF prints VCF text with the header as stored, without the PASS line bcftools adds (`#head-vcf-pass-line`); BED, GFF and TSV print their header and data lines verbatim.</li><li>No `@PG` or `##glociCommand` line is added, because `head` shows the file as stored.</li><li>A BAM whose text has no `@SQ` for its binary references gets synthesised `@SQ SN LN` lines. samtools 1.21 prints nothing for it: `#head-notext-header`.</li><li>`-n`/`-h` are not copied (`--records`, `--headers`).</li></ul> |
| D135 | **`gloci samples FILE…`** reads only headers. <ul><li>Per file, one line per distinct `--tag` value (default `SM`) of its `@RG` lines, in first-seen order, or `.` when there are none.</li><li>Columns: `SAMPLE PATH`; then `Y`/`N` with `--show-index` (D107's index search); then, with `--match-reference FASTA` (repeatable), the first matching FASTA in samtools' order, or `.`. samtools checks the last-given FASTA first. A FASTA matches when its `.fai` has the same names and lengths in the same order as the `@SQ` lines.</li><li>`--column-names` prints `#<TAG> PATH [INDEX] [REFERENCE]`. `--tag` must be two characters.</li><li>samtools prints samples in hash order: `#samples-order`. A FASTA without `.fai` is a usage error (D108; samtools builds one): `#fasta-implicit-index`.</li></ul> |
| D136 | **`gloci flags [FLAGS…]`, and FLAGS parsing** (amends D102 for every FLAGS option). <ul><li>Each argument prints `0x<hex> TAB <dec> TAB <names>`. Names are comma-joined in bit order: `PAIRED PROPER_PAIR UNMAP MUNMAP REVERSE MREVERSE READ1 READ2 SECONDARY QCFAIL DUP SUPPLEMENTARY`; none gives an empty field.</li><li>With no arguments, it prints samtools' 12-line table (`%#6x %5d  %-15s%s`), without samtools' usage text.</li><li>FLAGS is a decimal number, a `0x` hex number, or comma-separated names, case-insensitive, up to `0xFFFF`.</li><li>A decimal with a leading zero (`010`) is a usage error, because samtools reads it as octal 8: `#flags-leading-zero`.</li><li>Trailing text (`12abc`) is a usage error, where samtools reads 12.</li></ul> |
| D137 | **`gloci depth FILE…`: what counts.** <ul><li>**Records.** `AlnSource`, in file order. Records with no reference are skipped.</li><li>**Filters.** `--exclude-flags` defaults to `UNMAP,SECONDARY,QCFAIL,DUP`, samtools' set; `--exclude-flags 0` turns it off. The flag keeps its D102 meaning; only this command's default differs, and its help says so. Then `--require-flags` (all), `--min-mapq` (`-Q`), and `--min-read-length` (`-l`). The read length is `l_seq` minus its end soft clips, or, with SEQ `*`, the sum of `M I = X`.</li><li>**CIGAR walk** from `pos`:<ul><li>`M = X`: each base counts when `--min-baseq` (`-q`) is 0 or its quality is at least that. A missing quality (`0xFF`) passes. With SEQ `*`, every base passes; samtools reads past the record there.</li><li>`D`: counts only with `--count-deletions` (`-J`), each base gated by the quality of the *next* query base, as samtools does.</li><li>`N`: never counts.</li><li>`I S`: advance the query. `H P`: nothing.</li></ul></li><li>**The span** is the D101 span (`bam_endpos`; a placed unmapped read, kept only when `UNMAP` is not excluded, spans `[pos, pos+1)`), cut at `LN` with one warning. samtools prints rows past `LN`: `#depth-past-reference-end`.</li><li>**Order.** Any decrease in `(ref id, pos)` within an input is exit 1 with the `gloci sort` hint (D103). samtools notices only a read that ends before its last printed row: `#depth-order-check`.</li></ul> |
| D138 | **`depth`: which rows.** <ul><li>Each row is `chrom TAB pos(1-based) TAB depth…`, one depth column per input in argument order, in `@SQ` order and then position order.</li><li>A position is printed when it lies in a passing read's span in any input. Zero-depth positions inside spans (deletions, introns, filtered bases) print 0, as samtools.</li><li>`--all-positions` (`-a`) adds every position of each reference that has a passing read. `--all-references` (`-aa`) adds every position of every reference.</li><li>`--region` (repeatable, samtools `-r`) and `--regions-file` (`-b`) limit rows to their union, each position once. Reads come through the index when there is one; otherwise the input is scanned with a warning, where samtools `-r` stops with an error.</li><li>With regions, `--all-positions` prints every position in the regions. That equals samtools `-a -r` and `-aa -b`. samtools `-a -b` skips regions on references with no passing read: `#depth-regions`.</li><li>`--column-names`: `#CHROM POS <path>…` (paths as given).</li><li>With indexes and `-t` > 1, references run as `run_chunks` jobs (consecutive small references grouped) and are concatenated in `@SQ` order.</li></ul> |
| D139 | **`depth --skip-mate-overlap`** (samtools `-s`), per input. <ul><li>It applies to a read with `PAIRED` set and `MUNMAP` clear.</li><li>If the read's name is pending, the read takes that entry's end as its clip, and the entry is removed. Bases (and counted deletions) before the clip are not counted.</li><li>Otherwise the read stores `name → end` when its mate is on the same reference with `pos ≤ mpos ≤ end`.</li><li>An entry is dropped once a read starts after its `mpos`, and at each reference change. So the table holds only reads whose mates are still ahead, which is ∝ depth.</li><li>samtools also stores reads whose `mpos` is −1 or behind them, and never drops entries. A later alignment with the same name, anywhere in the file, is then clipped by a stale end, even on another reference. Its memory grows with the number of pairs. Results differ only for such templates (a supplementary alignment after its resolved pair, an unknown mate position): `#depth-mate-overlap`.</li></ul> |
| D140 | **`gloci bedcov BED FILE…`** <ul><li>Each FILE needs a `.bai`/`.csi`/`.crai`. Without one: exit 2 with `gloci index FILE`, as samtools exits 2.</li><li>BED lines may come in any order, as in `getfasta`: each region is a random-access query and the output follows the BED. Comment, `track`, `browser` and blank lines are skipped and not echoed.</li><li>Each data line is echoed verbatim. Then come one column per FILE for each of the sum, the `--bases-at-depth N` count (`-d`), and the `--read-count` (`-c`), in that order, as samtools.</li><li>**Sum.** The total, over `[start, end)`, of the per-base depth of passing reads. It uses depth's core with `M = X` always, `D` unless `--split-deletions`, and `N` unless `--split`; these are the D105 flags with their meaning, so samtools `-j` = `--split --split-deletions`.</li><li>**Bases at depth.** Uses the same depth. samtools drops `D` and `N` from both columns whenever `-d` is given: `#bedcov-depth-threshold`.</li><li>**Read count.** Passing reads whose D101 span overlaps the region under §4.3. A zero-length region has no bases and counts no reads, as samtools 1.21 does.</li><li>**Filters.** `--exclude-flags` defaults to samtools' `UNMAP,SECONDARY,QCFAIL,DUP`; `--min-mapq` applies. Unmapped reads never count, as htslib's pileup skips them.</li><li>**Errors.** A BED chromosome not in the `@SQ` list is exit 1 naming the line, or skipped and counted with `--skip-invalid`. samtools prints an error, goes on, and exits 2: `#bedcov-bad-line`.</li><li>`--column-names`: `#chrom start end <path>_cov… <path>_depth… <path>_count…`.</li><li>Batches of 256 lines run on ≤ threads workers, each with its own indexed readers. Output stays in BED order.</li><li>Use it where bedtools `multicov` was used (§2 non-goal).</li></ul> |
| D141 | **Fixtures and compare modes.** <ul><li>P3c fixtures live in `xtask/diff/inputs/samtools/`, each SAM next to its BAM and `.bai`, rebuilt by `make-samtools.sh` in `gloci-diff:p0`, as D119.</li><li>Index-path goldens name the index as an `aux` file (`gloci idxstats flags.bam`), because the input is always renamed `input.<ext>`.</li><li>New compare mode `cram`: both outputs are decoded through `open_aln_source` (the fixtures need no reference); the header text minus `@PG` lines, the references and the records must be equal.</li><li>CRAM goldens decode with `samtools view -b -T ref.fa` first (as P3a Task 14), except the `#cram-counts-need-reference` pins, which run samtools on the CRAM with `REF_PATH=/nonexistent`.</li></ul> |

## Review Focus

1. **An unindexed BAM given to `idxstats`** prints exactly what samtools' scan prints, with one warning. So do a BAM whose index lacks metadata pseudo-bins, and a BAM whose references are not grouped. Task 2 `idxstats_scan_equals_index`, `idxstats_index_without_counts`, `idxstats_unsorted_scan`; goldens `idxstats-scan-vs-samtools`, `idxstats-index-without-counts-vs-samtools`.
2. **`depth` with overlapping mate pairs** (same start, partial overlap, containment, a supplementary after the pair, quality-filtered overlap) equals `samtools depth -s` except the pinned stale-entry template, with a bounded mate table. Task 7 `mate_overlap_table`, `mate_table_bounded`; goldens `depth-s-*`.
3. **Huge headers through `reheader`.** 200,000 `@SQ` lines, a header over many BGZF blocks ending mid-block, and renamed `SN`s copy every later block verbatim in linear time. Task 11 `reheader_huge_header`.
4. **`flagstat` on supplementary, secondary, duplicate-secondary and QC-failed reads,** and on mates on another reference at MAPQ 4 and 5, gives samtools' counts in both columns. Task 3 `flagstat_loop_table`; golden `flagstat-categories-vs-samtools`.
5. **Lowercase, IUPAC, `*`, `-`, spaces and CRLF in `dict` input** give samtools' `LN` and `M5` (the CRAM normalised digest). Task 10 `dict_md5_normalisation`; golden `dict-edge-vs-samtools`.

## Lanes

Task 1 runs first, alone, and merges before the lanes start. It adds the shared accessors, the gate list, the FLAGS rules, the CRAM prelude reader, `--column-names`, and the only manifest change. Then three lanes run (4-core host):

| Lane | Tasks | Cross-lane waits |
|---|---|---|
| A: file summaries | 2 → 3 → 4 → 5 | — |
| B: depth | 6 → 7 → 8 | — |
| C: FASTA, reheader, view | 9 → 10 → 11 → 12 | — |
| Join | 13 → 14 → 15 | after all lanes |

**Shared files that will conflict.** Never hand-merge a generated file.
- `SKILL.md`: each task regenerates it as its last step. On a conflict, take either side and rerun `cargo xtask docs --skill SKILL.md`.
- `Cargo.toml`, `crates/gloci-ops/Cargo.toml`, `Cargo.lock`: Task 1 only (`md-5`). A later conflict means an unplanned dependency; stop and ask.
- `crates/gloci/src/cli.rs`: Tasks 1–12. Append new structs at the end under `// P3c Task N` and new variants at the end of `Commands`; keep both sides. `CheckArgs` (Task 4), `ViewArgs` (Task 12).
- `crates/gloci/src/commands/mod.rs`, `crates/gloci/src/lib.rs`, `crates/gloci-ops/src/lib.rs`: one line per new command or module; keep both sides.
- `xtask/diff/cases.toml`: one appended `# P3c Task N` section per task. `docs/differences.md`: one section per anchor, written by its task. `docs/migration.md`: Task 15 only.
- `xtask/diff/inputs/samtools/make-samtools.sh`: one appended `# Task N` block per task, and new fixture paths with one owner each.
- Single owners: `gloci-io/src/record.rs` gate (Task 1), `gloci-ops/src/check.rs` (Task 4), `gloci-ops/src/view.rs` (Task 12), `gloci-io/src/index/query.rs` (Task 2), `gloci-testkit/src/differential.rs` (Task 11), `xtask/bench/*` and `Dockerfile.bench` (Task 13).

## File structure

```
crates/gloci-io/src/align/{record.rs (+accessors), flags.rs, source.rs (+regions), cram_prelude.rs}   T1
crates/gloci-io/src/index/query.rs (+alignment_counts)                                               T2
crates/gloci-io/src/align/reheader.rs                                                                T11
crates/gloci-ops/src/{aln_flags, idxstats, flagstat, head, samples, faidx, dict, reheader, bedcov}.rs T1–T5, T8–T11
crates/gloci-ops/src/depth/{mod, ring, rows, mates}.rs                                               T6, T7
crates/gloci-ops/src/{check, view}.rs                                                                T4, T12
crates/gloci/src/commands/{flags, idxstats, flagstat, head, samples, depth, bedcov, faidx, dict, reheader}.rs; check.rs, view.rs
crates/gloci-testkit/src/differential.rs (`cram` mode)                                               T11
xtask/diff/inputs/samtools/; xtask/bench/cases.toml                                                  T1–T13
```

---

### Task 1: Foundation and `gloci flags` (D125, D126, D136, D141)

**Files:**
- Modify:
  - `gloci-io/src/align/record.rs`: accessors on `Aln`, `SamAln` and `AlnRef`;
  - `gloci-io/src/align/source.rs`: `open_aln_source_regions`;
  - the module holding P3a's `parse_flags`: the D136 rules;
  - `gloci-io/src/record.rs`: the D100 gate lists (D125);
  - `crates/gloci-ops/Cargo.toml` (`md-5.workspace = true`), `Cargo.lock`;
  - `cli.rs`: `ColumnNamesArg`, and `Commands::Flags(FlagsArgs)` (`FLAGS…`; 2 examples);
  - `commands/mod.rs`, `gloci-ops/src/lib.rs`;
  - SKILL.md, cases.toml, differences.md: `#samtools-exit-codes`, `#flags-leading-zero`, and `#fasta-implicit-index`. The last is pinned by Tasks 5 and 9, which run in different lanes, so this task writes the section for both.
- Create:
  - `gloci-io/src/align/flags.rs`, `gloci-io/src/align/cram_prelude.rs`;
  - `gloci-ops/src/aln_flags.rs`, `gloci/src/commands/flags.rs`;
  - `xtask/diff/inputs/samtools/make-samtools.sh`, and `samtools/empty.{sam,bam}` (header only; later tasks use it).
- Test: `gloci-io/tests/p3c_access.rs`, `gloci/tests/flags.rs`.

**Interfaces (produces):**
- On `Aln<'a>`, `SamAln<'a>` and `AlnRef<'a>`:
  - `fn mate_ref_id(&self) -> Option<usize>`, `fn mate_pos(&self) -> Option<u64>`, `fn seq_len(&self) -> usize`;
  - `fn quals(&self) -> Cow<'a, [u8]>`: raw Phred values; empty for SEQ `*`; `0xFF` × `seq_len` for a missing QUAL (SAM `*`).
- `pub const SAM_FLAGS: [(u16, &str, &str); 12]` (bit, name, samtools' description) and `pub fn flag_names(flag: u16) -> String`.
- `pub fn open_aln_source_regions(path: &str, regions: &[(Vec<u8>, Interval)], o: &ReaderOptions) -> Result<Box<dyn AlnSource>, IoError>`.
  - It reads through `.bai`/`.csi`/`.crai` when present (the reader P3a Task 6 builds for `view`, and `CramFile` from Task 14), else scans with one warning.
  - It yields the records whose D101 span (`[pos, pos+1)` for a placed unmapped read) overlaps a region under §4.3, each once, in file order.
- `pub struct CramPrelude { pub major: u8, pub minor: u8, pub file_id: [u8; 20], pub header_text: Vec<u8>, pub data_start: u64 }`
  - `pub fn read_cram_prelude(path: &Path) -> Result<CramPrelude, IoError>`: noodles' `header_reader`, `raw_sam_header_reader` and `discard_to_end`.
  - `pub fn cram_eof_present(path: &Path) -> io::Result<bool>`: the 3.x 38-byte EOF container; for 2.1, its EOF container.
- `pub struct ColumnNamesArg { #[arg(long)] pub column_names: bool }`.
- `pub fn flags_line(arg: &str) -> Result<String, UsageError>` and `pub fn flags_table() -> String`, in `aln_flags.rs`.

- [ ] **Step 1: Write the failing tests.**
  - `parse_flags_samtools_rules`:
    - `0x704`→0x704, `1796`→1796, `secondary,dup`→0x500, `Supplementary`→0x800, `0`→0, `4096`→0x1000, `65535`→0xFFFF;
    - errors: `65536`, `12abc`, `-1`, `0x`, `""`;
    - `010` → Err containing `leading zero` and `0x8`;
    - `bogus` → Err listing `PAIRED`.
  - `flag_names_table`: 0x704 → `UNMAP,SECONDARY,QCFAIL,DUP`; 0 → `""`; 0xFFF → all 12; 0x1000 → `""`.
  - `aln_accessors_bam_equals_sam`, on P3a's `flags.bam`/`flags.sam`, per record:
    - equal `mate_ref_id`, `mate_pos`, `seq_len` and `quals`;
    - a SAM QUAL `*` gives `0xFF`s, and SEQ `*` gives empty `quals`.
  - `regions_source_equals_filtered_scan`: on indexed `order.bam` and `flags.bam`, overlapping and bookended regions, a point, and a placed unmapped read. It gives the same records as `open_aln_source` filtered by §4.3, each once; unindexed, the same with one warning.
  - `cram_prelude_reads_header`: P3a's 3.0 and 3.1 fixtures give their header text, and `data_start` equals the first data container's offset; `cram_eof_present` is false after cutting the last 38 bytes.
  - `gate_lists_p3c_commands`: `gloci merge x.bam` names `depth` and `flagstat`; `gloci sort ref.fa` names `faidx` and `dict`.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test p3c_access && cargo test -p gloci --test flags`. Expected: FAIL.
- [ ] **Step 3: Add golden cases** (`input` is `samtools/empty.sam`, which neither side reads):
  - `flags-vs-samtools`: `gloci flags 0x704 1796 secondary,dup UNMAP,MUNMAP 0 4096` against `samtools flags …`;
  - `flags-table-vs-samtools`: `gloci flags` against `samtools flags | grep -E '^ +0x'`;
  - `flags-bad-name-vs-samtools` (compare `exit-status`, pinned `#samtools-exit-codes`);
  - `flags-octal-vs-samtools`: `010` (compare `exit-status`, pinned `#flags-leading-zero`).
- [ ] **Step 4: Implement.** Write `make-samtools.sh` with the D119 header (`--no-PG`, one thread, level 6) and a `# Task 1` block that builds `empty.bam`. Run `cargo xtask diff --only flags-`, then `--update-golden`.
- [ ] **Step 5: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 6: Commit** `feat: P3c foundation (alignment accessors, regions source, CRAM prelude, FLAGS rules) and gloci flags`

### Task 2: `gloci idxstats` (D127)

**Files:**
- Create: `gloci-ops/src/idxstats.rs`, `gloci/src/commands/idxstats.rs`.
- Fixtures:
  - `samtools/unsorted-chroms.{sam,bam}`: chr2 reads between two chr1 runs;
  - `samtools/nometa/flags.bam` (a copy of P3a's) and `samtools/nometa/flags.bam.bai`: its index without pseudo-bins, written by the `#[ignore]` test `make_nometa_index`. The directory keeps the basenames that `aux` files need.
- Modify:
  - `gloci-io/src/index/query.rs`: `alignment_counts`;
  - `cli.rs`: `Commands::Idxstats(IdxstatsArgs)` (`FILE`, `--reference`, threads); 2 examples;
  - `commands/mod.rs`, `gloci-ops/src/lib.rs`;
  - SKILL.md, cases.toml, differences.md (`#idxstats-unsorted`, `#idxstats-index-without-counts`, `#cram-counts-need-reference`).
- Test: `gloci/tests/idxstats.rs`.

**Interfaces:**
- Consumes: `IndexedFile` (P3a Task 6), `open_aln_source` (P3a Tasks 3 and 13).
- Produces:
  - `pub struct RefCounts { pub mapped: u64, pub unmapped: u64 }`
  - `pub struct AlnIndexCounts { pub per_ref: Vec<Option<RefCounts>>, pub unplaced: u64 }`. `None` marks a reference with bins but no pseudo-bin; a reference with no bins is `Some(0, 0)`.
  - `impl IndexedFile { pub fn alignment_counts(&self) -> Option<AlnIndexCounts> }` (BAI/CSI over a BAM only)
  - `pub enum IdxstatsPath { Index, Scan }`
  - `pub fn idxstats(path: &str, o: &ReaderOptions, out: &mut dyn Write) -> Result<IdxstatsPath, IoError>`

- [ ] **Step 1: Write the failing tests.**
  - `idxstats_index_matches_samtools`: `flags.bam` + `flags.bam.bai` → `Index`, byte-equal to the committed samtools output.
  - `idxstats_scan_equals_index`: the same BAM without its index → `Scan`, the same bytes, and the D127 warning.
  - `idxstats_index_without_counts`: `nometa/flags.bam.bai` → one warning, `Scan`, and the true counts.
  - `idxstats_unsorted_scan`: `unsorted-chroms.bam` → exit 0 with correct counts.
  - `idxstats_cram`: P3a's CRAM with `--reference ref.fa` gives the BAM's output; the reference-needing CRAM without `--reference` exits 1 with D109's message.
  - `idxstats_huge_header`: 200,000 `@SQ`, reads on every 1,000th, a gloci-built index. Index path < 5 s (debug), 200,001 lines.
  - Determinism of the scan at -t 1/2/8.
- [ ] **Step 2: Add golden cases** (`samtools idxstats`):
  - `idxstats-index-vs-samtools`: aux `flags.bam`, `flags.bam.bai`; `gloci idxstats flags.bam`;
  - `idxstats-scan-vs-samtools`: input `order.bam`, no index;
  - `idxstats-sam-vs-samtools`;
  - `idxstats-cram-vs-samtools`: gloci with `--reference ref.fa`, against `samtools view -b -T ref.fa input.cram -o t.bam && samtools idxstats t.bam`;
  - `idxstats-unsorted-vs-samtools` (compare `exit-status`, pinned `#idxstats-unsorted`);
  - `idxstats-index-without-counts-vs-samtools`: aux `samtools/nometa/flags.bam` and `samtools/nometa/flags.bam.bai` (pinned `#idxstats-index-without-counts`);
  - `idxstats-cram-no-reference-vs-samtools`: `REF_PATH=/nonexistent samtools idxstats input.cram` (compare `exit-status`, pinned `#cram-counts-need-reference`). If samtools also fails, the case goes stale; then drop the pin and the anchor.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci idxstats (index metadata, scan fallback)`

### Task 3: `gloci flagstat` (D128)

**Files:**
- Create: `gloci-ops/src/flagstat.rs`, `gloci/src/commands/flagstat.rs`.
- Fixtures:
  - `samtools/flagstat.{sam,bam}`: 160 QC-passed records, 153 of them mapped, so samtools prints `95.63%` where f64 gives `95.62%`. 160 are paired with 7 singletons (`4.37%` against `4.38%`). Also: secondary, supplementary, duplicate primaries and duplicate secondaries, read1/read2, proper pairs including an unmapped one with `0x2`, mates on another reference at MAPQ 4 and 5, a mapped read whose RNEXT is `*`, placed and unplaced unmapped reads, and 12 QC-failed records with no paired reads (`N/A`).
- Modify: `cli.rs` (`Commands::Flagstat(FlagstatArgs)`: `FILE`, `-o`, `-O`, `--reference`, threads; 3 examples), `commands/mod.rs`, `gloci-ops/src/lib.rs`, SKILL.md, cases.toml.
- Test: `gloci-ops/tests/flagstat.rs`, `gloci/tests/flagstat.rs`.

**Interfaces:**
- Consumes: `open_aln_source`, `AlnRef::{flag, ref_id, mate_ref_id, mapq}` (Task 1).
- Produces:
  - `#[derive(Default)] pub struct FlagStats`, with one `[u64; 2]` per counter (index 1 = QC-failed): `total primary secondary supplementary duplicates primary_duplicates mapped primary_mapped paired read1 read2 proper_pair both_mapped singletons mate_other_ref mate_other_ref_q5`;
  - `fn add(&mut self, flag: u16, ref_id: Option<usize>, mate_ref_id: Option<usize>, mapq: u8)`;
  - `pub enum FlagstatFormat { Text, Tsv, Json }`;
  - `pub fn percent(n: u64, total: u64, json: bool) -> String`;
  - `pub fn render_flagstat(s: &FlagStats, f: FlagstatFormat) -> String`;
  - `pub fn flagstat(src: Box<dyn AlnSource>) -> Result<FlagStats, IoError>`.

`add` is samtools' `flagstat_loop`:

```rust
let w = (flag & 0x200 != 0) as usize;
total[w] += 1;
if flag & 0x100 != 0 { secondary[w] += 1 }
else if flag & 0x800 != 0 { supplementary[w] += 1 }
else {
    primary[w] += 1;
    if flag & 0x1 != 0 {
        paired[w] += 1;
        if flag & 0x2 != 0 && flag & 0x4 == 0 { proper_pair[w] += 1 }
        if flag & 0x40 != 0 { read1[w] += 1 }
        if flag & 0x80 != 0 { read2[w] += 1 }
        if flag & 0x8 != 0 && flag & 0x4 == 0 { singletons[w] += 1 }
        if flag & 0xC == 0 {
            both_mapped[w] += 1;
            if mate_ref_id != ref_id { mate_other_ref[w] += 1; if mapq >= 5 { mate_other_ref_q5[w] += 1 } }
        }
    }
    if flag & 0x4 == 0 { primary_mapped[w] += 1 }
    if flag & 0x400 != 0 { primary_duplicates[w] += 1 }
}
if flag & 0x4 == 0 { mapped[w] += 1 }
if flag & 0x400 != 0 { duplicates[w] += 1 }
```

The percentage bases are: mapped/total, primary mapped/primary, proper/paired, and singletons/paired.

- [ ] **Step 1: Write the failing tests.**
  - `flagstat_loop_table`, (flag, ref, mate ref, MAPQ) → the counters that move:
    - `0x41|0x2|0x20` (primary paired proper read1) moves `total primary paired proper_pair read1 both_mapped mapped primary_mapped`;
    - `0x100|0x400` moves `total secondary mapped duplicates` (not `primary_duplicates`);
    - `0x800` moves `total supplementary mapped`;
    - `0x1|0x8` moves `singletons`;
    - `0x1|0x4|0x2` moves `paired` and not `proper_pair`;
    - `0x1`, ref 0, mate ref 1, MAPQ 4 moves `mate_other_ref` only; MAPQ 5 also moves `_q5`;
    - `0x1`, mate ref `None` counts as another reference;
    - `0x200|…` moves column 1 only.
  - `percent_f32`: (153,160) → `95.63%`; (7,160) → `4.37%`; (1,160) → `0.63%`; (0,0) → `N/A`; JSON (0,0) → `null`, (153,160) → `95.63`.
  - `flagstat_format_inference`: `-o s.json` gives JSON; `-O tsv -o s.json` gives TSV; `-O csv` → exit 2 listing `text, tsv, json`.
  - `flagstat_cram_equals_bam`; determinism at -t 1/2/8.
- [ ] **Step 2: Add golden cases** (`samtools flagstat`):
  - `flagstat-categories-vs-samtools` on `flagstat.bam`, with `-O tsv` and `-O json`;
  - SAM input; `empty.bam`;
  - CRAM: decode with samtools first;
  - `flagstat-cram-no-reference-vs-samtools` (pinned `#cram-counts-need-reference`, as Task 2).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci flagstat (text, tsv, json; samtools' f32 percentages)`

### Task 4: `check` for alignments, and `check --quick` (D129)

**Files:**
- Modify: `gloci-ops/src/check.rs`, `gloci/src/commands/check.rs`, `cli.rs` (`CheckArgs`: `--quick`, `--list-failed`, `--unaligned`, `--reference`; example `gloci check --quick --list-failed *.bam > bad.txt`), SKILL.md, cases.toml, differences.md (`#quickcheck`).
- Fixtures, from `make-samtools.sh` plus the `#[ignore]` helper `make_broken_bams`:
  - `samtools/noeof.bam`;
  - `samtools/truncated.bam`, cut mid-block;
  - `samtools/nosq.{sam,bam}`, unaligned with `@RG`;
  - `samtools/badfirst.bam`: header and EOF fine, first `block_size` 20;
  - `samtools/coord-unsorted.{sam,bam}`: `SO:coordinate`, one read out of place.
- Test: `gloci/tests/check_alignments.rs`.

**Interfaces:**
- Consumes: `open_aln_source`, `read_cram_prelude`, `cram_eof_present` (Task 1), `open_input`'s EOF rule.
- Produces:
  - `pub struct QuickOptions { pub unaligned: bool }`
  - `pub fn quick_check(path: &str, o: &ReaderOptions, q: &QuickOptions) -> CheckReport` (problems recorded; `has_errors()` as today);
  - `check_file` handles alignments (D129).

- [ ] **Step 1: Write the failing tests.**
  - `check_alignment_full`:
    - `coord-unsorted.bam` → exit 1, naming `record 7 of chr1` and both positions;
    - the same records with `SO:queryname` → exit 0 with `sorted: no`;
    - a read past `LN` → a warning;
    - the summary shows mapped, unmapped-placed and unplaced counts.
  - `quick_check_matrix` (exit status, `--list-failed` stdout):
    - good BAM, SAM, CRAM, VCF.gz → 0 and empty;
    - `noeof.bam`, `truncated.bam`, `badfirst.bam`, a missing path → 1 and the path;
    - `nosq.bam` → 1, or 0 with `--unaligned`;
    - CRAM with its EOF container cut → 1;
    - a reference-needing CRAM → 0 with no `--reference`;
    - stdin → no EOF check.
  - `quick_reads_little`: on a 100 MB BAM, at most the header, one BGZF block and the last 28 bytes are read (a counting reader).
  - `list_failed_only_paths`: with three files, stdout is exactly the failing paths.
- [ ] **Step 2: Add golden cases** (`samtools quickcheck -v` stdout; compare `exact` unless noted):
  - good BAM; `noeof.bam`; `nosq.bam` without and with `--unaligned` (`-u`);
  - `badfirst.bam` (pinned `#quickcheck`);
  - `input.vcf.gz` (pinned `#quickcheck`);
  - `noeof.bam` with compare `exit-status` (pinned `#samtools-exit-codes`).
- [ ] **Step 3: Implement, then update the goldens.** For CRAM, `--quick` reads the first data container at `data_start` with noodles' `read_container`, which parses its header and blocks without decoding slices, so no reference is needed.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci check reads alignments; check --quick replaces samtools quickcheck`

### Task 5: `gloci head` and `gloci samples` (D134, D135)

**Files:**
- Create: `gloci-ops/src/{head.rs, samples.rs}`, `gloci/src/commands/{head.rs, samples.rs}`.
- Fixtures:
  - `samtools/rg.{sam,bam}`: four `@RG` lines with samples `NA12878`, `NA12891`, `NA12878` and `HG002`, a `LB` tag on two, and one without `SM`. Choose the names so samtools' khash order differs from header order (checked when the golden is made);
  - `samtools/notext.bam`: binary references with empty header text, from `sam_to_bam` with the text dropped.
- Modify: `cli.rs` (`Commands::Head(HeadArgs)`: `FILE`, `--headers N`, `--records N`, `--reference`; `Commands::Samples(SamplesArgs)`: `FILES…`, `--tag`, `--show-index`, `--match-reference`, `ColumnNamesArg`, `-o`; 2–3 examples each), `commands/mod.rs`, `gloci-ops/src/lib.rs`, SKILL.md, cases.toml, differences.md (`#samples-order`; a case pinning `#fasta-implicit-index`).
- Test: `gloci/tests/head.rs`, `gloci/tests/samples.rs`.

**Interfaces:**
- Consumes: `open_aln_source`, `open_alignment_writer` (P3a Task 3), `find_index` (D107), `IndexedFasta::entries` (P3a Task 8), and the VCF/BED readers.
- Produces:
  - `pub struct HeadOptions { pub headers: Option<u64>, pub records: u64 }`
  - `pub fn head(path: &str, o: &HeadOptions, ro: &ReaderOptions, out: &mut dyn Write) -> Result<(), IoError>`
  - `pub struct SamplesOptions { pub tag: [u8; 2], pub show_index: bool, pub match_reference: Vec<PathBuf>, pub column_names: bool }`
  - `pub fn samples(paths: &[String], o: &SamplesOptions, out: &mut dyn Write) -> Result<(), IoError>`

- [ ] **Step 1: Write the failing tests.**
  - `head_defaults_and_limits`:
    - the default prints the whole header and no record;
    - `--headers 0 --records 3` prints 3 SAM lines;
    - `--headers 2` cuts after 2 lines;
    - more records than the file holds → all, exit 0;
    - `notext.bam` prints synthesised `@SQ\tSN:…\tLN:…` lines;
    - no `@PG` is added;
    - a truncated BAM prints what it read and exits 1.
  - `head_vcf_and_bed`: VCF.gz and BCF print VCF text with no `##glociCommand`; BED prints `#` lines and data lines verbatim.
  - `samples_rows`:
    - `rg.bam` → `NA12878`, `NA12891`, `HG002` in that order;
    - `--tag LB`;
    - `--tag X` → exit 2;
    - `nosq.bam` without `@RG` → `.`;
    - `--show-index` Y/N;
    - `--match-reference` with two matching FASTAs → the last given;
    - a FASTA without `.fai` → exit 2 with `gloci index`.
  - `samples_huge_header_match`: 200,000 `@SQ` against a 200,000-entry `.fai`: < 5 s (debug).
- [ ] **Step 2: Add golden cases:**
  - head (`samtools head`):
    - `head-vs-samtools`: default;
    - `head-records-vs-samtools`: `--headers 2 --records 3` against `-h 2 -n 3`;
    - `--headers 0 --records 5`;
    - SAM input; `notext.bam`;
    - `head-vcf-vs-bcftools`: `gloci head input.vcf.gz --records 2` against `bcftools head -n 2`.
  - samples:
    - `samples-vs-samtools`: compare `lines-unordered` on `rg.bam`;
    - `samples-order-vs-samtools`: compare `exact`, pinned `#samples-order`;
    - `--column-names --tag LB` against `-h -T LB`;
    - `--show-index` against `-i`;
    - `--match-reference ref.fa` against `-f ref.fa` (aux `ref.fa`, `ref.fa.fai`);
    - two files;
    - `samples-fasta-no-fai-vs-samtools` (compare `exit-status`, pinned `#fasta-implicit-index`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci head (any format) and gloci samples`

### Task 6: The depth core and `gloci depth` (D137, D138)

**Files:**
- Create: `gloci-ops/src/depth/{mod.rs, ring.rs, rows.rs}`, `gloci/src/commands/depth.rs`.
- Fixtures:
  - `samtools/depth.{sam,bam,bam.bai}`: chr1, chr2 (`LN` 300) and an empty chr3. It has reads with `D`, `N` (a 100 bp intron), `I`, `S`, `H`, `=`/`X` and a mapped `*` CIGAR; a deletion whose next base has quality 10; QUAL `*`; a read past chr2's `LN`; duplicate, secondary, QC-fail and supplementary reads; and a placed unmapped read.
  - `samtools/depth2.{sam,bam,bam.bai}`: the same `@SQ`, other reads;
  - `samtools/depth-reordered.{sam,bam}`: the same names in another order;
  - `samtools/depth-unsorted.{sam,bam}`: a read 3 bp before its predecessor but ending after samtools' last row;
  - `samtools/targets.bed`: chr1 regions, one on chr3, and overlapping ones.
- Modify: `cli.rs` (`Commands::Depth(DepthArgs)`: `FILES…`, `--all-positions`, `--all-references`, `--count-deletions`, `--min-baseq`, `--min-read-length`; its own `--exclude-flags` (default `UNMAP,SECONDARY,QCFAIL,DUP`) beside `--min-mapq`, `--require-flags` and `--reference`, from which it builds `AlignOptions`; `RegionArgs`, `ChromMapArg`, `ColumnNamesArg`, `-o`, threads; 4 examples), `commands/mod.rs`, `gloci-ops/src/lib.rs`, SKILL.md, cases.toml, differences.md (`#depth-past-reference-end`, `#depth-order-check`, `#depth-regions`, `#depth-bedcov-same-sq`).
- Test: `gloci-ops/tests/depth_core.rs`, `gloci/tests/depth.rs`.

**Interfaces:**
- Consumes: `open_aln_source`, `open_aln_source_regions`, `AlnRef::{for_each_op, quals, seq_len, span}` (Task 1), `run_chunks` (P2a/P3a), `IndexedFile::weight` (P3a Task 6).
- Produces (Tasks 7 and 8 use them):
  - `pub struct CountRules { pub deletions: bool, pub ref_skips: bool, pub min_baseq: u8 }`
  - `pub struct DepthRing`: one input's counts, keyed by absolute position.
    - `fn add(&mut self, r: &AlnRef<'_>, span: Interval, rules: &CountRules, clip: u64)`: zero the positions beyond the current covered end, then walk the CIGAR. `clip` is Task 7's (0 = none).
    - `fn covered_end(&self) -> u64`
    - `fn take(&mut self, pos: u64) -> u32`: the depth at `pos`; positions before it are released.
    - `fn reset(&mut self)`
    - `fn capacity(&self) -> usize`
  - `pub enum AllPositions { Covered, RefsWithReads, AllRefs }`
  - `pub struct DepthOptions { pub rules: CountRules, pub min_read_length: u64, pub all: AllPositions, pub regions: Option<Vec<(Vec<u8>, Interval)>>, pub skip_mate_overlap: bool, pub column_names: bool }`
  - `pub struct DepthSummary { pub rows: u64, pub max_ring: usize, pub max_pending_mates: usize }`
  - `pub fn depth(paths: &[String], o: &DepthOptions, ro: &ReaderOptions, threads: NonZeroUsize, out: &mut dyn Write) -> Result<DepthSummary, IoError>`

Rows (samtools' `add_depth` order):
- Within a reference, keep `next`, the first position not yet printed.
- When a record at `pos` arrives from any input (records merged by `(ref id, pos)`, ties to the lower input index), print each position in `[next, pos)` that is below some input's `covered_end`, or every one under `--all-*`.
- At a reference change, and at the end, flush up to the largest `covered_end`, or to `LN` under `--all-*`.
- `--all-references` also prints references with no record, in `@SQ` order.
- Every printed position is filtered by the regions.

- [ ] **Step 1: Write the failing tests.**
  - `depth_counting_table` (one `DepthRing`): (CIGAR, quals, rules) → per-position depth:
    - `5M`, quals 30, `-q 20` → 1s;
    - `2M3D2M` → `1 1 0 0 0 1 1`, or all 1 with deletions;
    - deletions with `-q 20`, next base quality 10 → the `D` positions are 0;
    - `3M100N3M` → 100 zeros;
    - `2S3M2I3M` → 6 positions;
    - `=`/`X` count;
    - QUAL `0xFF` passes `-q 50`;
    - SEQ `*` with `-q 20` → counts.
  - `depth_spliced_rows`: `depth.bam` prints every intron position of the spliced read, with 0.
  - `depth_all_modes`:
    - `--all-positions` prints chr1 and chr2 whole but not chr3;
    - `--all-references` prints chr3 too;
    - `--region chr3:1-10 --all-positions` prints 10 zero rows.
  - `depth_inputs`:
    - two files give two columns;
    - `depth-reordered.bam` with `depth.bam` → exit 1 naming `@SQ 1`;
    - `depth-unsorted.bam` → exit 1 with `gloci sort`.
  - `depth_past_ln`: one warning; no row past 300.
  - `depth_memory_bounded`: 100,000 reads of 150 bp plus one read with a 500 kb `N` → `max_ring` ≤ 2^20.
  - `depth_per_chromosome_equals_streaming`: indexed, -t 1/2/8, default, `--all-references`, and with regions.
  - `depth_all_refs_huge_header`: 200,000 `@SQ` of `LN` 10 with `--all-references` → 2,000,000 rows in < 20 s (debug).
- [ ] **Step 2: Add golden cases** (`samtools depth`, on `depth.bam` unless noted):
  - default, `-a`, `-aa`, `-J`, `-q 20`, `-Q 30`, `-J -q 20`, `-l 50`;
  - `--exclude-flags UNMAP,SECONDARY,QCFAIL` against `-g DUP`;
  - `--exclude-flags UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` against `-G SUPPLEMENTARY`;
  - `--exclude-flags 0` against `-g 0x704` (the placed unmapped read prints a 0 row in both);
  - `--require-flags PROPER_PAIR`;
  - `--column-names` against `-H`;
  - `--region chr1:100-200` against `samtools index input.bam && samtools depth -r chr1:100-200 input.bam`;
  - `--regions-file targets.bed` against `-b`;
  - `-aa -b`;
  - two files (aux `depth2.bam`) with `-H`;
  - SAM input; CRAM, decoded with samtools first;
  - pinned:
    - `-a -b` (`#depth-regions`);
    - the past-`LN` read (`#depth-past-reference-end`);
    - `depth-reordered.bam` with `depth.bam` (compare `exit-status`, `#depth-bedcov-same-sq`);
    - `depth-unsorted.bam` (compare `exit-status`, `#depth-order-check`).
- [ ] **Step 3: Implement, then update the goldens.** Format rows with `itoa` into a reused buffer.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci depth (samtools counting rules, all-positions modes, per-chromosome)`

### Task 7: `depth --skip-mate-overlap` (D139)

**Files:**
- Create: `gloci-ops/src/depth/mates.rs`.
- Fixture: `samtools/mates.{sam,bam}`, with templates a: [100,200) + [150,250); b: both at 300; c: [400,600) containing [450,500); d: a non-overlapping pair; e: a deletion before the clip; f: `MPOS` −1 with `MUNMAP` clear; g: a supplementary (PAIRED, MUNMAP clear) starting inside its second mate after the pair resolved; h: overlap bases of quality 5 in the second mate.
- Modify: `depth/mod.rs`, `cli.rs` (`DepthArgs --skip-mate-overlap`; example), SKILL.md, cases.toml, differences.md (`#depth-mate-overlap`).
- Test: `gloci-ops/tests/depth_mates.rs`.

**Interfaces:**
- Consumes: `DepthRing::add`'s `clip` (Task 6), `AlnRef::{name, flag, mate_ref_id, mate_pos}`.
- Produces: `pub struct MateTable`, with:
  - `fn clip_for(&mut self, r: &AlnRef<'_>, pos: u64, end: u64) -> u64`: consume, store or neither, per D139;
  - `fn advance(&mut self, pos: u64)`: drop entries with `mpos < pos`;
  - `fn reset(&mut self)`;
  - `fn len(&self) -> usize`.
  
  It holds a `HashMap<Box<[u8]>, (u64 end, u64 mpos)>` and a min-heap on `mpos` for eviction, with stale heap entries skipped.

- [ ] **Step 1: Write the failing tests.**
  - `mate_overlap_table` (depth per position):
    - a gives 1 over [100,250);
    - b gives 1 over the longer read;
    - c gives 1 over [400,600);
    - d gives both counted, with `len() == 0` afterwards;
    - e with `--count-deletions`: the deletion before the clip is not counted;
    - f is not stored;
    - h with `-q 10`: the second mate's bases after the clip use its own qualities.
  - `mate_table_bounded`: 1,000,000 non-overlapping pairs at depth 30 → `max_pending_mates` ≤ 60.
  - `mates_per_chromosome_equals_streaming`: -t 1/2/8.
- [ ] **Step 2: Add golden cases** (`samtools depth -s` on `mates.bam` minus template g): `-s`, `-s -J`, `-s -q 10`, `-s -a`. Template g alone, pinned `#depth-mate-overlap`.
- [ ] **Step 3: Implement, then update the goldens.** For each record, call `clip_for` (consume or store) first, then `advance(pos)`. A mate at exactly `mpos` then finds its entry before the entry can be dropped.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: depth --skip-mate-overlap (samtools -s with a bounded mate table)`

### Task 8: `gloci bedcov` (D140)

**Files:**
- Create: `gloci-ops/src/bedcov.rs`, `gloci/src/commands/bedcov.rs`.
- Fixtures: `samtools/bedcov.bed` (lines out of order with 3–6 columns, one CRLF; `#`, `track` and blank lines; a zero-length region inside a read; regions only over the intron and only over a deletion; overlapping regions; one past chr2's `LN`), and `bedcov-badchrom.bed`, which adds a `chrZ` line.
- Modify: `cli.rs` (`Commands::Bedcov(BedcovArgs)`: `BED`, `FILES…`, `--split`, `--split-deletions`, `--bases-at-depth N`, `--read-count`, the filter flags as `DepthArgs` declares them (D140 default), `ChromMapArg`, `InvalidArg`, `ColumnNamesArg`, `-o`, threads; 3 examples), `commands/mod.rs`, `gloci-ops/src/lib.rs`, SKILL.md, cases.toml, differences.md (`#bedcov-depth-threshold`, `#bedcov-bad-line`; `#depth-past-reference-end` gains its bedcov line).
- Test: `gloci-ops/tests/bedcov.rs`, `gloci/tests/bedcov.rs`.

**Interfaces:**
- Consumes: `DepthRing`, `CountRules` (Task 6); `open_aln_source_regions` (Task 1); the BED reader without order checks (as `getfasta`).
- Produces:
  - `pub struct RegionTotals { pub sum: u64, pub at_depth: u64, pub reads: u64 }`
  - `pub fn region_totals(src: &mut dyn AlnSource, ref_len: u64, iv: Interval, rules: &CountRules, at_depth: Option<u32>) -> Result<RegionTotals, IoError>`
  - `pub struct BedcovOptions { pub rules: CountRules, pub bases_at_depth: Option<u32>, pub read_count: bool, pub column_names: bool }`
  - `pub fn bedcov(bed: Box<dyn RecordReader>, files: &[String], o: &BedcovOptions, ro: &ReaderOptions, threads: NonZeroUsize, out: &mut dyn Write) -> Result<u64, IoError>` (returns the lines written)

- [ ] **Step 1: Write the failing tests.**
  - `bedcov_equals_naive`: 500 seeded reads with `D`/`N` and 300 regions (points, bookends, past `LN`). For each of the four `--split`/`--split-deletions` combinations, `sum`, `at_depth` (N = 2) and `reads` equal a naive per-base count.
  - `bedcov_output_order`: unordered BED → output in BED order, at -t 1/2/8 with 10,000 regions.
  - `bedcov_needs_index`: exit 2 naming `gloci index depth2.bam`.
  - `bedcov_bad_chrom`: exit 1 naming line 9; `--skip-invalid` → one count warning, and the other lines are written.
  - `bedcov_unmapped_never_counted`: `--exclude-flags 0` leaves the placed unmapped read uncounted in `--read-count`.
  - `bedcov_whole_chromosome_region`: a 250 Mb region over 1,000,000 reads keeps `capacity()` ≤ 2^12.
- [ ] **Step 2: Add golden cases** (`samtools bedcov`; aux `depth.bam`, `depth.bam.bai`, `depth2.bam`, `depth2.bam.bai`):
  - default, `-c`;
  - `-j` against `--split --split-deletions`;
  - `--split` alone and `--split-deletions` alone (each against itself; no samtools flag);
  - `-d 2 -j` against `--bases-at-depth 2 --split --split-deletions` (match);
  - `-Q 30`;
  - `-g DUP` (as Task 6);
  - two files with `-H -c -d 1 -j`;
  - CRAM, decoded first;
  - pinned:
    - `-d 2` alone against `--bases-at-depth 2` (`#bedcov-depth-threshold`);
    - `bedcov-badchrom.bed` (compare `exit-status`, `#bedcov-bad-line`);
    - the past-`LN` region (`#depth-past-reference-end`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci bedcov (per-region depth sums, thresholds, read counts)`

### Task 9: `gloci faidx` (D131)

**Files:**
- Create: `gloci-ops/src/faidx.rs`, `gloci/src/commands/faidx.rs`.
- Modify:
  - `gloci-core/src/region.rs`: `chr:-E`;
  - `cli.rs`: `Commands::Faidx(FaidxArgs)`: `REF`, `REGION…`, `RegionArgs`, `--line-width`, `--reverse-complement`, `--mark-strand`, `InvalidArg`, `OutputArgs`; 3 examples;
  - `commands/mod.rs`, `gloci-ops/src/lib.rs`;
  - SKILL.md, cases.toml, differences.md (`#faidx-missing-sequence`; a case pinning `#fasta-implicit-index`).
- Test: `gloci-core` region tests, `gloci-ops/tests/faidx.rs`, `gloci/tests/faidx.rs`.

**Interfaces:**
- Consumes: `IndexedFasta::{open, len_of, entries, fetch, fetch_rev}` (P3a Task 8), `seq::complement` (P3a Task 9), the FASTA output path with `--write-index` (P3a Task 11).
- Produces:
  - `pub enum MarkStrand { Rc, No, Sign, Custom(String, String) }`, with `FromStr`;
  - `pub struct FaidxOptions { pub line_width: Option<u64>, pub reverse_complement: bool, pub mark: MarkStrand }`;
  - `pub struct FaidxRegion { pub label: String, pub chrom: Vec<u8>, pub iv: Interval }`;
  - `pub fn faidx(fasta: &mut IndexedFasta, regions: &[FaidxRegion], o: &FaidxOptions, out: &mut dyn Write) -> Result<(), IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `region_open_start`: `chr1:-100` = `chr1:1-100`; `chr1:0-5` stays an error.
  - `faidx_wrapping`: the default width equals the `.fai` line width (60, and 70 for a second FASTA); `--line-width 7`; `--line-width 0` gives one line; a 3 MiB region wraps across 1 MiB chunks with no short line except the last.
  - `faidx_complement_table`: `seq::complement(b, false)` equals samtools' `comp_base` for bytes 0–255 (a constant copied from `faidx.c`).
  - `faidx_edges`:
    - an end past the sequence → cut, with a warning;
    - a start past it → the header only;
    - an unknown name → exit 1 before any output; with `--skip-invalid`, skipped;
    - no regions → exit 2 naming `gloci index`;
    - no `.fai` → exit 2.
  - `faidx_regions_file`: a BED gives `>chr1:5-70`-style headers; a zero-length line is skipped with a count warning.
  - `faidx_memory`: a whole 300 Mb sequence holds ≤ 1 MiB plus one line.
- [ ] **Step 2: Add golden cases** (`samtools faidx`; aux `ref.fa` and `ref.fa.fai` from P3a's `fa/`):
  - `faidx-vs-samtools`: `gloci faidx ref.fa chr1:5-70 chr2:60000-60100 chrM`;
  - `-n 7` and `-n 0`;
  - `-i` and `-i --mark-strand sign|no|custom,+,-`;
  - lowercase and IUPAC under `-i`;
  - `chr1:1,000-1,050`, `{chr1}:1-10`, `chr1:-100`;
  - a region past the end; a start past the end;
  - BGZF `ref.fa.gz` with `.fai`/`.gzi`;
  - pinned:
    - an unknown name with `--skip-invalid` against `-c` (`#faidx-missing-sequence`);
    - `input.fa` without `.fai` (compare `exit-status`, `#fasta-implicit-index`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci faidx (samtools regions, wrapping, reverse complement)`

### Task 10: `gloci dict` (D132)

**Files:**
- Create: `gloci-ops/src/dict.rs`, `gloci/src/commands/dict.rs`.
- Fixtures:
  - `fa/dict-edge.fa` (owned here): lowercase, IUPAC, `*` and `-`, a space and a tab inside sequence lines, CRLF, an empty sequence, `chrM`/`MT`/`M` names, and a description after the name;
  - `samtools/dup.fa`;
  - `samtools/alt.txt` (a `.alt`-style list).
- Modify: `cli.rs` (`Commands::Dict(DictArgs)`: `REF`, `--no-header`, `--assembly`, `--species`, `--uri`, `--alias`, `--alt`, `-o`; 2 examples), `commands/mod.rs`, `gloci-ops/src/lib.rs`, SKILL.md, cases.toml, differences.md (`#dict-duplicate-names`).
- Test: `gloci-ops/tests/dict.rs`, `gloci/tests/dict.rs`.

**Interfaces:**
- Consumes: `FastaStream` (P3a Task 8), `md-5`.
- Produces:
  - `pub struct DictOptions { pub hd_line: bool, pub assembly: Option<String>, pub species: Option<String>, pub uri: Option<String>, pub alias: bool, pub alt: Option<HashSet<Vec<u8>>> }`
  - `pub fn normalized_md5(state: &mut md5::Md5, line: &[u8]) -> u64` (feeds the kept bytes and returns how many were kept)
  - `pub fn dict(fasta: FastaStream, path_for_uri: Option<&Path>, o: &DictOptions, out: &mut dyn Write) -> Result<(), IoError>`

- [ ] **Step 1: Write the failing tests.**
  - `dict_md5_normalisation`:
    - `ACGT` and ` AC\tgt\r` → `f1f8f4bf413b16ad135722aa4591043e` with `LN:4`;
    - `acgtNNRYKM*-` counts 12 bytes, and its MD5 equals that of the uppercased bytes;
    - an empty sequence → `LN:0 M5:d41d8cd98f00b204e9800998ecf8427e`.
  - `dict_tags`: `--alias` on `chr1`, `1`, `chrM`, `M` and `MT` → `AN:1`, `AN:chr1`, `AN:M,chrMT,MT`, `AN:chrM,chrMT,MT` and `AN:chrMT,chrM,M`. Tags come in samtools' order. `UR:file://<canonical path>` by default; none for stdin.
  - `dict_duplicate_names`: exit 1 naming both lines.
  - `dict_round_trip`: `gloci dict ref.fa -o ref.dict`, then `-g ref.dict` equals `-g ref.fa.fai` in `gloci sort`.
  - `dict_streams`: a 1 GB generated FASTA holds ≤ one line.
- [ ] **Step 2: Add golden cases** (`samtools dict`; both sides pass `--uri x`/`-u x` so `UR` is not machine-specific):
  - `dict-vs-samtools` on `ref.fa`;
  - `dict-edge-vs-samtools`;
  - `--no-header` against `-H`;
  - `--alias --assembly GRCh38 --species human` against `-A -a … -s …`;
  - `--alt alt.txt` against `-l`;
  - `ref.fa.gz` (BGZF); plain gzip;
  - `dup.fa` (compare `exit-status`, pinned `#dict-duplicate-names`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci dict (samtools-identical LN and M5)`

### Task 11: `gloci reheader` (D130, D141 `cram` mode)

**Files:**
- Create: `gloci-io/src/align/reheader.rs`, `gloci-ops/src/reheader.rs`, `gloci/src/commands/reheader.rs`.
- Fixtures: `samtools/new-header.sam` (renamed `SN`s, a changed `@RG`, one appended `@SQ`), `samtools/bad-ln.sam`, `samtools/bad-order.sam`.
- Modify: `gloci-testkit/src/differential.rs` (`cram` mode), `cli.rs` (`Commands::Reheader(ReheaderArgs)`: `HEADER`, `IN`, `-o`, `--compression-level`; 3 examples, one being `gloci view --header-only x.bam -O sam | sed 's/SN:chr/SN:/' | gloci reheader - x.bam -o y.bam`), `commands/mod.rs`, `gloci-ops/src/lib.rs`, SKILL.md, cases.toml, differences.md (`#reheader-sq-check`).
- Test: `gloci-io/tests/reheader.rs`, `gloci/tests/reheader.rs`.

**Interfaces:**
- Consumes: `AlnHeader::{read_bam, parse_sam, write_bam, with_program}` (P3a Tasks 2 and 3), `bgzf_copy::copy_from` (P3a Task 15), `read_cram_prelude` (Task 1).
- Produces:
  - `pub fn read_bam_header_with_extent(r: &mut dyn Read) -> Result<(AlnHeader, u64), String>`: wraps the decompressed stream in an unbuffered counting reader, so the extent (the header's byte length) is exact even when `l_text` includes NUL padding.
  - `pub fn check_sq_compatible(old: &AlnHeader, new: &AlnHeader, cram: bool) -> Result<(), String>` (D130 message)
  - `pub fn write_cram_header_container(major: u8, text: &[u8], out: &mut dyn Write) -> io::Result<()>` (ITF8/LTF8, CRC32 via `flate2::Crc` on 3.x)
  - `pub fn reheader(header_src: &str, input: &Path, out: Output, o: &WriterOptions, command_line: &str) -> Result<(), IoError>`

- [ ] **Step 1: Write the failing tests.**
  - `sq_compatibility_table`:
    - renamed names → Ok;
    - an appended `@SQ` → Ok;
    - a changed `LN`, a swapped order, or a dropped `@SQ` → Err naming `@SQ N` of both files;
    - CRAM with a different `M5` → Err; `M5` missing on one side → Ok.
  - `reheader_bam_blocks`: after the new header's blocks, the output's blocks equal the input's from the header-end block onward (that block recompressed, the rest byte for byte). It decodes to the input records under the new header, plus `@PG`.
  - `reheader_huge_header`:
    - 200,000 `@SQ` (a header of about 6 MB over many blocks, ending mid-block), renamed;
    - < 10 s (debug);
    - the largest buffer ≤ both headers + 64 KiB;
    - the records are unchanged.
  - `reheader_cram`: the 3.0 and 3.1 fixtures; the data containers are copied byte for byte from `data_start`; P3a's `CramReader` decodes the output with the new names (and with `--reference` renamed to match).
  - `reheader_inputs`:
    - stdin as IN → exit 2;
    - `-` as HEADER works;
    - a SAM IN → exit 2;
    - a missing EOF → exit 1;
    - `-o IN` replaces atomically; a failure leaves IN untouched.
  - Determinism at -t 1/2/8.
- [ ] **Step 2: Add the `cram` compare mode** with a unit test: equal header text minus `@PG`, references and records → equal; a changed record → not equal.
- [ ] **Step 3: Add golden cases:**
  - `reheader-bam-vs-samtools`: `gloci reheader new-header.sam input.bam` against `samtools reheader -P new-header.sam input.bam` (compare `bam`);
  - a BAM header file as HEADER;
  - `reheader-cram-vs-samtools` on P3a's `embed_ref` CRAM (compare `cram`);
  - pinned: `bad-ln.sam` and `bad-order.sam` (compare `exit-status`, `#reheader-sq-check`).
- [ ] **Step 4: Implement, then update the goldens.**
- [ ] **Step 5: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 6: Commit** `feat: gloci reheader (BAM block copy, CRAM header container, @SQ checks)`

### Task 12: `view --count`, and `--regions-file` on alignments (D133)

**Files:**
- Modify: `gloci-ops/src/view.rs` (`ViewOptions.count`), `gloci/src/commands/view.rs`, `cli.rs` (`ViewArgs --count`; examples `gloci view --count reads.bam chr1` and `gloci view reads.bam --regions-file targets.bed -o on_target.bam`), SKILL.md, cases.toml, differences.md (`#view-overlapping-regions`, `#view-regions-union`).
- Test: `gloci/tests/view_count.rs`.

**Interfaces:**
- Consumes: `open_aln_source_regions` (Task 1), P3a's alignment `view` (Task 3).
- Produces: `ViewOptions.count: bool`. `run_view` writes `N\n` instead of records.

- [ ] **Step 1: Write the failing tests.**
  - `view_count_formats`: BAM, SAM, CRAM (with `--reference`), VCF, BCF and BED, with and without regions, equal the record count of the same `view` without `--count`.
  - `view_count_conflicts`: with `--header-only`, `-O bam` or `--write-index` → exit 2.
  - `view_count_filters`: `--exclude-flags UNMAP` and `--min-mapq 30` on `flags.bam` equal samtools' `-F 4` and `-q 30` counts (committed numbers).
  - `regions_file_uses_index`: an indexed BAM with a 3-line BED reads at most the index's chunks (a counting reader); unindexed, it gives the same records with one warning.
  - Determinism at -t 1/2/8.
- [ ] **Step 2: Add golden cases** (`samtools view`; gloci reads `input.bam` unindexed, and the reference indexes first where samtools needs it):
  - `view-count-vs-samtools`: `--count input.bam` against `-c`;
  - `--count input.bam chr1:1-150 chr2:95-105`;
  - `--count --exclude-flags UNMAP,SECONDARY` against `-c -F 0x104`;
  - `view-regions-file-vs-samtools`: `-O sam --no-header --regions-file targets.bed` against `-L targets.bed`, with a zero-length region and a placed unmapped read;
  - `--count --regions-file targets.bed` against `-c -L`;
  - pinned:
    - `--count input.bam chr1:1-150 chr1:100-200` (`#view-overlapping-regions`);
    - `--count input.bam chr1 --regions-file targets.bed` against `-c -L targets.bed input.bam chr1`, where `targets.bed` also covers chr2 (`#view-regions-union`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: view --count, and samtools -L parity for --regions-file on alignments`

### Task 13: P3c benchmarks

**Files:** `xtask/bench/cases.toml` and `xtask/src/bench.rs`.
- The input kind `bam-indexed` is P3a Task 18's `bam`, then `gloci index`.
- The `bedcov` targets are every 30th window of `gloci windows -g x.bam -w 100000`, so their names match the BAM.

samtools 1.21 is the binding baseline here, because it is the tool these commands replace (§2 anti-goal). mosdepth is a yardstick for `depth`.

| Case | Target |
|---|---|
| `flagstat-bam-50M` (`-t 1` and `-t 8`) | `gloci <= samtools flagstat` (`-@ 0` / `-@ 7`) |
| `idxstats-bam-200k-sq` (indexed) | `gloci <= samtools idxstats` |
| `view-count-bam-50M` | `gloci <= samtools view -c` |
| `depth-bam-50M` (`-t 1`; and `--all-positions`) | `gloci <= samtools depth`; yardstick `mosdepth -t 8 --fast-mode` |
| `bedcov-bam-10k-targets` | `gloci <= 0.5 * samtools bedcov` (samtools uses pileup) |
| `faidx-10k-regions`, `dict-3G`, `reheader-bam-50M` | `gloci <= samtools` |

- [ ] **Step 1: Write the failing tests.** `every_gloci_bench_command_parses` covers the new cases; the dry-run guard keeps zero-output `depth`/`bedcov` cases out.
- [ ] **Step 2: Implement.**
- [ ] **Step 3: Run** `cargo test -p xtask` and `cargo xtask bench --only depth-bam-50M --records 1000000`. Paste the medians.
- [ ] **Step 4: Commit** `bench: P3c cases against samtools 1.21`

### Task 14: Behaviour matrices and large-file tests

**Files:** `crates/gloci/tests/behaviour.rs`, `crates/gloci/tests/large.rs`.

Extend the matrices to `idxstats`, `flagstat`, `check --quick`, `reheader`, `head`, `samples`, `flags`, `faidx`, `dict`, `view --count`, `depth` and `bedcov`:
- determinism at -t 1/2/8;
- broken pipe (`depth --all-references | head`, `faidx | head`);
- empty inputs:
  - a header-only BAM: `flagstat` zeros with `N/A`; `idxstats` zero rows; `depth` no rows, or all-zero rows with `--all-references`; `bedcov` zeros; `head` the header; `samples` `.`;
  - a FASTA with no sequences;
- failures:
  - a BAM missing its EOF block (every reader exits 1);
  - a truncated CRAM;
  - a missing `.bai` for `bedcov` and a missing `.fai` for `faidx` (exit 2);
  - an unsorted BAM for `depth` (hint);
  - a missing output directory (no partial file);
- a `--chrom-map` round trip for `depth` and `bedcov` (BED names and output rows).

RSS ceilings (`#[ignore]`):

| Run | Ceiling |
|---|---|
| `flagstat`, `idxstats` (scan), `view --count`, `depth --all-positions`, `bedcov` with 10,000 targets, `reheader`, each on a 2 GB BAM | < 200 MB |
| `faidx` of a whole 300 Mb sequence; `dict` of 3 GB | < 200 MB |
| `depth --skip-mate-overlap` on 2 GB of pairs | < 200 MB |

- [ ] **Step 1: Write the matrices.**
- [ ] **Step 2: Run** `cargo test -p gloci --test behaviour`, then the large tests with `GLOCI_LARGE_BYTES=200000000`. Fix each real failure where it lives, with a regression test.
- [ ] **Step 3: Commit** `test: P3c determinism, broken-pipe, empty-input and failure matrices; RSS ceilings`

### Task 15: Documentation

**Files:** `docs/migration.md`, `docs/differences.md` (links), SKILL.md.

- [ ] **Step 1: Add a `### samtools` table to `docs/migration.md`,** checked against `--help`.
  - Move P3a's `index`, `faidx` (indexing), `sort`, `merge`, `cat` and `view` rows into it.
  - Add one row per P3c command, mapping every samtools flag its D-row names (D125–D140). Include the pipes: `xargs gloci faidx ref.fa < list` for `faidx -r`, and `view --header-only | sed | reheader -` for `reheader -c`.
  - Then a "not planned" row: `markdup`, `fixmate`, `calmd`, `mpileup`, `consensus`, `stats`, `ampliconclip`, `fqidx`, and the D125 options.
  - Add the bedtools row `multicov` → `bedcov --read-count`.
  - Add a "samtools behaviour" section:
    - `depth` and `bedcov` exclude `UNMAP,SECONDARY,QCFAIL,DUP` by default, unlike `coverage`/`genomecov` (D102);
    - `depth` prints zero rows inside spliced reads, as samtools;
    - each read is counted once by `view --count`;
    - CRAM needs `--reference` for record-reading commands;
    - exit codes;
    - no implicit FASTA indexes.
  - Link every P3c anchor from its row.
- [ ] **Step 2: Regenerate SKILL.md,** adding pitfalls for:
  - the `depth`/`bedcov` filter default;
  - `check --quick --list-failed` (not `-v`);
  - `view --count` (not `-c`);
  - `--split` and `--split-deletions` in `bedcov` (samtools `-j`);
  - `faidx` needing `gloci index` first.
- [ ] **Step 3: Run** `cargo test -p xtask` and the help-examples test. Expected: PASS.
- [ ] **Step 4: Commit** `docs: samtools migration table, behaviour section and SKILL.md for P3c`

---

## Execution notes for the controller

- **Order.** Task 1 merges first. Then lanes A, B and C run in separate worktrees, each with its own `target/`, with `CARGO_BUILD_JOBS=1`, scoped test commands, and one `cargo test --workspace` per task before merging. Tasks 13–15 run last, in order.
- **Docker.** No image change: samtools and bcftools 1.21 are in `gloci-diff:p0` since P3a. Each task runs its `make-samtools.sh` block in that image and commits the fixtures it creates.
- **Spec amendments** (controller):
  - §7.4: `check` reads alignments, and `--quick` (D129); `view --count` and the regions union (D133); new entries for `idxstats`, `flagstat`, `reheader`, `head`, `samples`, `flags`, `faidx`, `dict`, `depth` and `bedcov` (D127–D140);
  - §6.1: FLAGS parsing (D136); `--column-names` is per command;
  - §4.1 / D8: `chr:-E` (D131);
  - §10: the P3c row's "done when";
  - §11: the samtools table.
- **P3a interfaces this plan assumes.** The named P3c task adds each one where P3a left it out:
  - `Aln` mate, quality and length accessors (Task 1);
  - a public indexed `AlnSource` for regions; P3a Task 6 builds one inside `view` (Task 1);
  - `AlnHeader::read_bam` reading without its own buffering, so the header extent can be counted (Task 11);
  - htslib's synthesised `@SQ` lines for a BAM with empty text, in the SAM rendering (Task 5);
  - `IndexedFile::alignment_counts` over the BAI/CSI pseudo-bins (Task 2).
