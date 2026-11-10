# gloci P3a Sequence and Alignment I/O Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the sequence and alignment half of spec Phase 3:
- BAM, CRAM and SAM as interval input to `intersect`, `coverage`, `genomecov`, `convert` (bedtools `bamtobed`) and `view`;
- BAM and SAM output;
- `.bai`, `.csi` and `.crai`, with region queries and per-chromosome parallelism;
- FASTA with `.fai` and `.gzi`;
- `getfasta`, `maskfasta` and `nuc`;
- PSL and RepeatMasker in `convert`;
- alignment sort (coordinate, `-n`, `-N`), merge and cat, with samtools 1.21's ordering and headers: `gloci sort`, `gloci sort -m`, `gloci cat`.

**Architecture:** No new crates.
- **Alignments share one layout.** BAM records stay raw bytes in the batch buffer, decoded by hand like BCF (D22). CRAM records are decoded by noodles-cram and transcoded to BAM bytes, so spans, blocks, filters, BED rendering, BAM output and indexing handle one layout. SAM keeps its line for unchanged passthrough, behind the same accessors.
- **Commands opt in.** A new `ReaderOptions.accept` lets eight commands read alignments; the rest refuse them with a hint.
- **Two read paths.** Interval commands read alignments as `RecordReader` intervals, with unmapped reads skipped. Record-preserving commands (`view`, `sort`, `sort -m`, `cat`) read every record through `AlnSource` (D120). Alignment sort reuses the external sorter with a pluggable `Order`.
- **FASTA** is its own reader family (`gloci-io::fasta`): a chunked random-access fetcher over `.fai`/`.gzi`, and a sequential line stream.
- **Indexes** extend the binning-index code. BAI is noodles' `Index<LinearIndex>`, the TBI type.

**Tech Stack:** Everything P2a uses, plus:
- noodles-sam 0.93, noodles-bam 0.98, noodles-cram 0.102 (feature `libdeflate`) and noodles-fasta 0.68;
- samtools 1.21 in both images, and mosdepth 0.3.10 in the bench image.

**Spec:** `docs/superpowers/specs/2026-11-06-gloci-design.md`: §2, §4.1–§4.6, §5.1–§5.2, §6.1, §6.4, §7.1, §7.3 (`genomecov`), §7.4 (`index` rows and `convert`), §7.5, §8.3, §8.4, §9, §10 (the P3 row and the Order note, as amended in d9a1635), and §11 (bedtools rows).

**Builds on:**
- P2a, fully merged, plus P2b Task 1. D1–D79 and D82 stay in force; D80–D98 are untouched.
- P2a Task 19's validation contract applies unchanged: plain inputs are drained, indexed ones are not, and errors are the same at every `-t` except the documented double fault. Alignments are BGZF-like inputs: never drained, and read per chromosome when indexed.

**Out of scope:**
- P3c, the next phase, has its own plan: `idxstats`, `flagstat`, `quickcheck`, `reheader`, `faidx` regions, `dict`, `view -c`/`-L`, `head`, `samples`, `flags`, `depth` and `bedcov`. It builds on `AlnSource` (Task 3), the writers (Task 3), the indexes (Tasks 6, 13) and `IndexedFasta` (Task 8);
- P3b: `stats`, `sample`/`random`/`shuffle`, BEDPE (with bamtobed `-bedpe`), `vcf popgen` and `vcf compare`;
- writing CRAM, and samtools features beyond indexing;
- the §2 BAM utilities (`multicov`, `tag`, …);
- `REF_PATH`/`REF_CACHE`;
- the D77 pools other than CRAM's (D110);
- Python.

**P3a is done when:**
- every item has differential cases against bedtools 2.31.1, samtools 1.21 or bedops 2.4.41, with committed goldens;
- every intentional difference is anchored in `docs/differences.md` and linked from `docs/migration.md`;
- these §11 rows map to working commands: `bamtobed`, `getfasta`, `maskfasta`, `nuc`, `genomecov -ibam`, `intersect -abam`, `coverage` with BAM, and `--cram-ref`;
- samtools `sort [-n|-N]`, `merge [-n|-N]`, `cat` and `index` map to `gloci sort`, `sort -m`, `cat` and `index`, with records and headers equal to samtools 1.21's;
- the nightly bench meets the Task 18 targets.

## Global Constraints

P2a's Global Constraints all apply (u64 half-open coordinates, `rust-version = "1.91"`, no htslib, exit codes 0/1/2, errors that say what/where/how, atomic BGZF output, content detection, record streams, 2–4 help examples, CI on three platforms). New in P3a, binding on every task:
- **Goldens.** Every bedtools- or samtools-facing behaviour is a `xtask/diff/cases.toml` case whose gloci side is one plain invocation, so it is committed under `crates/gloci/tests/golden/` and replayed by plain `cargo test`. Index interop needs samtools, so it is a live case, also pinned by a committed-fixture test.
- **Differences.** The task that pins a difference writes its anchored section in `docs/differences.md`. Task 20 links each from `docs/migration.md`.
- **Goldens discriminate.** Fixtures have contig orders where lexicographic, karyotype and header order all differ (`chrM, chr1, chr2, chr10, chrX`), tied reads whose names sort differently as text and as numbers (`read10`, `read9`), and a `nuc` interval whose `%` value differs between f32 and f64.
- **Threads.** Output is byte-identical at `-t 1/2/8`, with a test for every new path.
- **Memory bounds**, each asserted in its task:

  | Path | Bound |
  |---|---|
  | BAM/SAM streaming | ≤ 4 MiB batch + the BGZF pool; slim B reads (D115); active set = reads overlapping the current A |
  | CRAM | ≤ 2 × threads containers + one reference sequence per decoder |
  | alignment `genomecov` | pending events ∝ depth, never ∝ chromosome length |
  | `getfasta`/`nuc` | ≤ 1 MiB of sequence per record |
  | `maskfasta` | one line + 16 B per merged mask |
  | alignment `sort`/`sort -m` | `--max-mem` (record bytes + 40 B per record); a record over the budget gets a chunk of its own. Temp disk is about the input BAM's size (LZ4 runs of the decoded records), up to 2× during a cascaded merge of more than 128 runs |
  | `cat` | one BGZF block |
- **Contig count.** Nothing may be quadratic in the number of contigs; 200,000-`@SQ` tests guard this (Tasks 2 and 6).
- **No panics.** noodles-cram's `expect` on a missing reference, and its `todo!` in `fs::index`, are never reached (D109).

## Decisions this plan adds

| # | Decision |
|---|---|
| D99 | **Libraries.** noodles-sam 0.93, -bam 0.98, -cram 0.102 and -fasta 0.68 are the latest releases. They pin the locked noodles-bgzf 0.52, -csi 0.62 and -core 0.21, and rust 1.91, and they are pure Rust, CRAM 3.1 codecs included. rust-htslib is excluded by §8.4. <ul><li>BAM is decoded by hand, like BCF.</li><li>noodles handles the SAM header, BAM↔SAM, CRAM, and BAI/CSI/CRAI/FAI/GZI I/O.</li><li>noodles-fasta's `IndexedReader` is not used, because it chooses BGZF by extension (against §5.1).</li></ul> |
| D100 | **Formats.** <ul><li>Detection by content:<ul><li>BAM: a BGZF stream that decompresses to `BAM\1`;</li><li>CRAM: `CRAM` with major version 2 or 3; 4 is an error;</li><li>SAM: an `@HD/@SQ/@RG/@PG/@CO` + TAB line, or ≥ 11 tab fields with integer fields 2, 4, 5, 8 and 9;</li><li>FASTA: a leading `>`;</li><li>PSL and RepeatMasker: their header, or their column shape.</li></ul></li><li>`--format` adds `bam`, `cram`, `sam`, `psl` and `rmsk`.</li><li>Only `intersect`, `coverage`, `genomecov`, `convert`, `view`, `index`, `sort` and `cat` read alignments. Elsewhere: `x.bam is BAM; this command does not read alignments (intersect, coverage, genomecov, convert, view, index, sort and cat do); convert it to BED first: gloci convert x.bam -o x.bed` (exit 2). In `gloci merge`, the message adds `to merge sorted alignment files as samtools merge does, use: gloci sort -m` (D122).</li><li>PSL and RepeatMasker are read by `convert` only. FASTA is never an interval input.</li><li>SAM is never split by D76.</li><li>Output: `.bam`, `.sam`, and `.fa`/`.fasta`/`.fna` (`.gz` = BGZF), or `-O bam\|sam\|fasta\|fasta.gz`. `.cram` is a usage error suggesting `samtools view -C`.</li></ul> |
| D101 | **An alignment as an interval.** <ul><li>The span is `[pos, pos + reflen)`, with reflen = the `M D N = X` lengths, and at least 1 for a mapped read with no reference ops (htslib's `bam_endpos`).</li><li>Clips never widen it. The strand is `-` for flag `0x10`, else `+`.</li><li>Unmapped reads (`0x4`, placed or not) are skipped by every interval command, counted at `-v`. bedtools `intersect -v` writes them: `#unmapped-reads`. The record-preserving commands keep every record (D120).</li><li>More than 65,535 ops come from `CG:B:I`.</li><li>A read past its reference's `LN` (chrM wrap) is cut at `LN`, with one warning per input. bedtools `genomecov` clips the same way.</li><li>Errors name `record N of CHROM`, counted per chromosome as the indexed reader counts, so every `-t` reports one place. SAM errors name the line.</li></ul> |
| D102 | **Read filters** (the spec had none): `--min-mapq N`, `--require-flags FLAGS` and `--exclude-flags FLAGS`. <ul><li>FLAGS is a number or comma-separated samtools names.</li><li>The filters apply to every alignment input first.</li><li>The defaults filter nothing, as bedtools. `--exclude-flags 0x704` gives samtools `depth` and mosdepth semantics.</li><li>Counts are reported at `-v`.</li></ul> |
| D103 | **Header order.** <ul><li>An alignment header's `@SQ` order is a declared order, like a VCF header in §4.5 rule 2 (a spec gap). The first header in input order wins, closed. `OrderSource::AlignmentHeader(path)` hints `gloci sort b.bed -g x.bam`.</li><li>An out-of-order alignment hints `gloci sort x.bam -o x.sorted.bam` (naming `SO:queryname` when declared), or `--sort` (D124). `--unsorted` holds B alignments in memory.</li><li>Reference id → chromosome id is a vector built once.</li><li>`--chrom-map` renames the binary reference list, `@SQ SN:`, and SAM `RNAME`/`RNEXT`.</li><li>`-g` also takes a CRAM.</li></ul> |
| D104 | **Alignments in output.** <ul><li>`intersect --mode any\|none` and `view` write alignments: BAM for a BAM or CRAM A, SAM for a SAM A, or as `-o`/`-O` choose. Asking for BED is a usage error suggesting `… -O bam \| gloci convert - -o x.bed` (D18).</li><li>Where a command prints an alignment's columns (pairs, count, `coverage`, B columns), it prints the BED6 line `convert` writes. bedtools prints BED12 with `0,0,0` and trailing commas: `#bam-as-bed-columns`.</li><li>`--mode clip` of an alignment is a usage error.</li><li>Output gets `@PG ID:gloci[.N] PN:gloci PP:… VN:… CL:…` (samtools' id rule).</li></ul> |
| D105 | **Blocks.** <ul><li>`--split` uses the reference segments between `N` ops; `D` stays inside (bedtools `intersect`/`coverage`/`bamtobed -split`).</li><li>`--split-deletions` also cuts at `D` (bamtobed `-splitD`), with or without `--split`. It is a usage error with no alignment input.</li><li>Empty blocks are dropped; with none left, the span is the block.</li><li>The rule is one for every command. So `genomecov` counts deleted bases, unlike bedtools `-ibam` (which cuts at `D` unless `-ignoreD`): `#bam-genomecov-deletions`. Parity is `--split-deletions`.</li></ul> |
| D106 | **`convert` of alignments** (bamtobed). <ul><li>Default BED6: `RNAME start end QNAME[/1\|/2] MAPQ strand`.</li><li>`--bed12`: thick = span; `--color R,G,B`, default `255,0,0`.</li><li>`--blocks` is `-split`; `--blocks --split-deletions` is `-splitD`.</li><li>`--score-tag TAG` (`-tag`; `-ed` = `NM`). A read missing the tag is a data error.</li><li>`--cigar` (`-cigar`) on BED6 only.</li><li>Input order is kept and unmapped reads are skipped.</li><li>bamtobed's extra column under `-split -tag` is not copied: `#bamtobed-split-tag`.</li></ul> |
| D107 | **Alignment indexes.** <ul><li>`gloci index x.bam` writes `.bai`, or `.csi` with `--csi` or when an `LN` exceeds 2^29.</li><li>Spans are htslib's. A placed unmapped read is indexed at `[pos, pos+1)` in the metadata counts, and unplaced reads go to `n_no_coor`.</li><li>gloci finds `x.bam.bai`, `x.bai`, `x.bam.csi`, `x.cram.crai` and `x.crai`. Stale ones warn (D26).</li><li>`--write-index` on BAM output works. On SAM output, and `gloci index` of a SAM, it is a usage error suggesting `gloci view x.sam -o x.bam`.</li><li>samtools and gloci must read each other's indexes with identical query results.</li></ul> |
| D108 | **FASTA.** <ul><li>Random access needs `.fai`, plus `.gzi` for BGZF. A missing index is a usage error with the exact `gloci index FILE` command; gloci never builds one implicitly (bedtools does). Plain gzip is an error with a `bgzip` recipe.</li><li>`gloci index` writes `.fai`/`.gzi` byte-identical to samtools/htslib: first-word names, equal line lengths (else an error naming the line), CRLF allowed.</li><li>Each fetch checks the line ends against `.fai` (else `does not match; rebuild it: gloci index ref.fa`) and streams chunks of ≤ 1 MiB, forwards or backwards.</li><li>`-g ref.fa` reads its `.fai`.</li></ul> |
| D109 | **The CRAM reference** is `--reference FASTA` only (bedtools `--cram-ref`): no `REF_PATH`, no cache directory, no network. <ul><li>Given one, each `@SQ` in the `.fai` must match its `LN`, else an error before any record.</li><li>A container that needs reference bases the FASTA lacks, or with no `--reference`: `x.cram needs the reference sequence of chr7 (M5 …, UR …) to decode its reads; pass it with --reference FASTA`.</li><li>A slice MD5 mismatch: `the sequence chr1 in ref.fa does not match the one x.cram was made with (chr1:10,001-20,000)`.</li><li>gloci's repository adapter never returns `None` (a missing name is `Some(Err)`), and slice reference ids are checked first. That keeps noodles' `expect`s unreachable.</li><li>gloci builds `.crai` itself: noodles' `fs::index` decodes multi-reference slices with an empty repository and can `todo!()`.</li><li>Each decoder clears its repository when the reference id changes, so it holds ≤ 1 sequence; per-chromosome runs hold ≤ threads.</li><li>Records are transcoded to BAM bytes.</li><li>`--reference` without a CRAM warns.</li></ul> |
| D110 | **CRAM decode pool.** Containers are read in order and decoded on ≤ threads workers, ≤ 2 × threads in flight, and released in order: the CRAM analogue of the §8.3 BGZF pool. Other D77 pools stay deferred unless Task 18 misses a target. |
| D111 | **`getfasta --fasta REF FILE`** <ul><li>Output is in input order, one line per sequence.</li><li>The header is `chrom:start-end` (0-based start, as bedtools). `--name` gives `name::chrom:start-end`; `--name-only` gives `name`.</li><li>`-s` appends `(+)`/`(-)`/`(.)` and reverse-complements `-` with bedtools' IUPAC-aware, case-keeping table. `--rna` writes U.</li><li>`--split` joins BED12 blocks. `--tab` writes `header TAB seq`; `--bed-out` writes the record, then the sequence.</li><li>A record outside the FASTA is a data error (bedtools skips it: `#getfasta-outside-fasta`); `--skip-invalid` skips and counts. Points are skipped with one count warning.</li><li>`--chrom-map` renames the intervals only. VCF follows `--vcf-span`. There is no `-fullHeader`.</li></ul> |
| D112 | **`nuc --fasta REF FILE`** <ul><li>TSV: the record's columns, then `pct_at pct_gc num_A num_C num_G num_T num_N num_oth seq_len`, then `seq` (`--seq`) and the pattern count (`--pattern P`).</li><li>bedtools' header line comes first: `#1_usercol … K_usercol K+1_pct_at …`, with K from the first record. An empty input has no header.</li><li>Counts are case-insensitive, `U` counts as T, and anything else is other.</li><li>`%` is computed in f32 and printed with 6 decimals, as bedtools prints it.</li><li>`-s` reverse-complements `-` records. Pattern matches overlap and are case-sensitive unless `--ignore-case` (`-C`).</li><li>Errors as `getfasta`.</li></ul> |
| D113 | **`maskfasta --fasta IN -b MASK… [-o OUT]`** <ul><li>The FASTA (any compression, no index) is streamed. Masked bases become `--mask-char C` (default `N`) or lowercase (`--soft`).</li><li>Headers and line breaks are kept exactly. bedtools re-wraps to the first line's width and drops descriptions: `#maskfasta-layout`, `#maskfasta-header`.</li><li>Masks (any order) are merged per chromosome. Absent chromosomes get one warning. A mask past its sequence's end is a data error (`#maskfasta-past-end`).</li><li>`.fa.gz` output is BGZF, and `--write-index` writes `.fai`/`.gzi`.</li></ul> |
| D114 | **PSL and RepeatMasker → BED** (convert only; this plan counts them as alignment and sequence I/O, so they are in P3a). PSL spans `[tStart, tEnd)` on `tName`; RepeatMasker spans `[begin-1, end)`, with strand `C` = `-`. Columns follow bedops 2.4.41 `psl2bed`/`rmsk2bed`, pinned against `convert2bed --do-not-sort`. Input order is kept. |
| D115 | **Slim records.** When a command never writes an input's alignments (every B; A outside `--mode any\|none`/`view`), the reader keeps only the fixed fields, the name, the resolved CIGAR and the `--score-tag` tag. `AlignOptions.full_records` disables this. A slim record never reaches a writer (debug assertion). |
| D116 | **`view` on alignments:** BAM↔SAM and CRAM→BAM/SAM, writing every record, unmapped included (D120); region queries through the index or by scan, with §4.3 on the D101 span (a placed unmapped read spans `[pos, pos+1)`, as samtools finds it); `--header-only`; `--no-header` for SAM; the D102 filters. It is not `samtools view`. |
| D117 | **Alignment `genomecov`.** <ul><li>`-g` defaults to the header's lengths, with empty chromosomes in header order (`#genomecov-chrom-order`).</li><li>`--5p`/`--3p` take the first/last aligned base by strand, as bedtools.</li><li>An indexed input at `-t` > 1 runs per chromosome, with `genome` rows summed and output identical to `-t 1`.</li><li>`-pc`, `-fs` and `-du` are not planned: no §7.3 flag, and a §2 non-goal.</li></ul> |
| D118 | **Tools and normalisation.** <ul><li>`gloci-diff` adds `samtools=1.21`; the bench image adds `samtools=1.21` and `mosdepth=0.3.10` (both on htslib 1.21).</li><li>`observe()` drops `@PG\tID:gloci` lines, and references use `samtools view --no-PG`.</li><li>New compare mode `bam`: decompress both, then compare the header text minus `@PG`, the references, and the record bytes.</li></ul> |
| D119 | **Binary fixtures** live in `xtask/diff/inputs/{aln,fa}/`, each next to its SAM/FASTA source. `make-<set>.sh` rebuilds them in `gloci-diff:p0` (`--no-PG`, one thread, level 6). A test checks that each BAM decodes to its source. |
| D120 | **Every-record alignment source.** <ul><li>`AlnSource` yields every record of a BAM, SAM or CRAM in file order, unmapped included, with the D102 filters when given.</li><li>Records come in the reader's raw layout: BAM bytes for BAM and CRAM, the line for SAM.</li><li>`view`, `sort`, `sort -m` and `cat` use it, and so will P3c's per-record commands. Interval commands keep `RecordReader` (D101).</li><li>Sort and merge encode SAM records to BAM bytes once, at read.</li></ul> |
| D121 | **`gloci sort` for alignments** follows samtools 1.21 `bam_sort.c` exactly. <ul><li>Coordinate (the default): `(tid, with -1 as n_ref so unplaced reads go last; pos+1; the reverse-strand bit; input order)`. samtools' LSD radix sort is stable and its run merge breaks ties by run, then index, so this is samtools' total order.</li><li>`-n`: `(strnum_cmp(QNAME), flag rank, input order)`, with samtools' `strnum_cmp` ported byte for byte: unsigned bytes, digit runs compared as numbers after leading zeros are skipped, so `r01` equals `r1`. The flag rank is `((f&0xC0)<<8)\|((f&0x100)<<3)\|((f&0x800)>>3)`: neither mate, then READ1, then READ2; within each, primary, then supplementary, then secondary.</li><li>`-N`: as `-n`, but QNAMEs compare as bytes (`strcmp`). `-n`/`-N` are the owner's §10 flags; `windows -n COUNT` keeps its meaning (an owner exception, like `groupby -g`). The long forms are `--by-name` and `--by-name-bytes`.</li><li>`@HD`, as `sam_hdr_update_hd`: `SO` replaced in place, or appended. Coordinate removes `SS`; `-n` sets `SS:queryname:natural`; `-N` sets `SS:queryname:lexicographical`. `GO` is removed. With no `@HD`, `@HD VN:1.6 SO:…` comes first. `@PG` follows D104.</li><li>One input (BAM, SAM or CRAM). CRAM writes BAM. `-o x.sam`/`-O sam` writes SAM.</li><li>`--write-index` works on coordinate output only (`.bai`, or `.csi` per D107).</li><li>`--check` checks the chosen order.</li><li>`-g`, `--full-key`, `--unique` and `--duplicates` with alignments are usage errors. The order is the header's.</li><li>The external sorter takes an `Order` (interval, coordinate, name natural, name bytes). Spill runs hold raw records with their `block_size`.</li></ul> |
| D122 | **Alignment merge is `gloci sort -m`.** `-m` already means "merge sorted inputs" in gloci (bedops `-u`, `sort -m`). `gloci merge` keeps its interval meaning, one name for one meaning (§2), and its alignment message points to `sort -m`. <ul><li>The keys are D121's (`-n`/`-N` select name order). Ties go to input file order, then record order (samtools' `heap_lt`).</li><li>Each input is order-checked as it is read. An input whose `@HD SO`/`SS` contradicts the mode is an error up front. An absent `SS` is accepted.</li><li>`@SQ` lists must be identical (names, `LN`, order), else an error naming the first difference. samtools appends unknown `@SQ` lines and only warns when coordinate order is lost: `#alignment-merge-sq`.</li><li>`@HD` comes from the first input, with `SO`/`SS` set as D121.</li><li>`@RG`/`@PG` with the same ID are combined: an identical line is kept once (samtools `-c -p`). A same-ID `@RG` with a different line is an error naming both lines (samtools `-c` keeps the first): `#alignment-merge-rg`. A same-ID `@PG` keeps the first. samtools' default suffixes colliding IDs with `lrand48` seeded by the clock, which cannot be reproduced. Goldens use `samtools merge -c -p`.</li><li>`@CO` lines come from every input, in order.</li><li>`-r` (an RG from the file name) is not supported.</li></ul> |
| D123 | **`gloci cat`** (samtools `cat`) joins BAMs byte-level. Its meaning is `vcf concat --naive`'s, not `vcf concat`'s checked join, so it takes samtools' name. <ul><li>The header is the first input's, plus each later input's `@RG` lines with new IDs (samtools).</li><li>Every input's `@SQ` list must equal the first's. samtools does not check, so records silently point at the wrong reference: `#alignment-cat-sq`.</li><li>Each later file's block holding its header end is recompressed from there; every later block is copied verbatim; EOF markers are dropped and one is written at the end. The format-neutral copy loop moves from `gloci-vcf/src/concat.rs` into `gloci-io` and serves both.</li><li>No order check. `@HD` is kept as samtools keeps it.</li><li>BAM only. CRAM is refused, suggesting `samtools cat`; SAM is refused, suggesting `gloci view x.sam -o x.bam`.</li></ul> |
| D124 | **`--sort` with alignments** (amends D103): <ul><li>relation commands and `genomecov` sort alignment inputs with the interval sorter, like BED, into a level-1 temporary BAM (slim records allowed there), then stream;</li><li>the table is the header's, so BED inputs must follow `@SQ` order or be sorted too;</li><li>memory stays within `--max-mem`.</li></ul> |

## Review Focus

1. **Huge alignment headers** (200,000 `@SQ`) stay linear. The known risks are `IndexedFile::rid_of` (a linear `position()` today) and the per-chromosome planner's disjoint-names warning (`Vec::contains` inside `any`). Task 2 `bam_huge_header_linear`, Task 6 `huge_header_per_chromosome`.
2. **CRAM with a missing, mismatched or `.gzi`-less reference, or truncated,** exits 1 with the contig and the fix, never an abort. A CRAM needing no reference decodes without one. Task 13 `cram_missing_reference`, `cram_mismatched_reference`, `cram_reference_bgzf`, `cram_truncated`, `cram_no_ref_needed`.
3. **Unmapped, secondary, supplementary, duplicate, MAPQ-0 and chrM-overhang reads** follow one rule in `convert`, `intersect`, `coverage` and `genomecov`, matching bedtools by default. Task 2 `bam_unmapped_and_filters`; the `flags.bam` goldens in Tasks 4, 5 and 7; Task 7 `genomecov-ibam-overhang-vs-bedtools`.
4. **CIGAR shapes** (clips, `I D N = X`, `*`, all-clip, > 65,535 ops) give bedtools' spans, blocks and BED12. Task 2 `bam_spans_table`, `bam_long_cigar_cg_tag`; Task 4 goldens on `cigars.bam`.
5. **Unusual FASTA:** lowercase, IUPAC, `U`, `N`, CRLF, a short last line, descriptions, empty sequences, BGZF without `.gzi`. Task 8 `fai_matches_samtools`, `fai_errors`; Task 9 `getfasta-strand-iupac-vs-bedtools`; Task 10 `nuc-lowercase-iupac-vs-bedtools`; Task 11 goldens.
6. **Sorting alignments** (owner-requested addition):
   - **tie order:** the same position on both strands, same-position reads kept in input order, unplaced reads last by strand bit; `r2`/`r10`, `a1b2`/`a1b10`, `r01`/`r1` and READ1/READ2/supplementary/secondary under `-n` and `-N`;
   - **huge records:** a long-read record larger than `--max-mem`;
   - **a name sort over millions of reads,** which stays within `--max-mem` and spills.

   Task 16 `aln-sort-ties-vs-samtools`, `aln-sort-name-vs-samtools`, `aln_sort_huge_record`, `strnum_cmp_matches_samtools`; Task 19 `sort -n` RSS ceiling.

## Lanes

Task 1 runs first, alone, and merges before the lanes start: it adds every shared enum variant and the `ReaderOptions` field. There are three lanes (4-core host).

| Lane | Tasks | Cross-lane waits |
|---|---|---|
| A: alignments | 2 → 3 → 4 → 5 → 6 → 7 | — |
| B: FASTA, then alignment sort | 8 → 9 → 10 → 11 → 16 → 17 | 16 after Task 6; 17 after Task 13 (CRAM inputs) |
| C: PSL, cat, CRAM, bench | 12 → 15 → 13 → 14 → 18 | 15 after Task 3; 13 after 3 and 8; 14 after 6 and 7; 18 after 11, 14 and 17 |
| Join | 19 → 20 | after all lanes |

Tasks 15–17 were added on the owner's request (alignment cat, sort and merge). They are numbered after the CRAM tasks so that Tasks 1–14, which P3c's plan cites, keep their numbers.

**Shared files that will conflict.** Never hand-merge a generated file.
- `SKILL.md`: each task that changes help regenerates it as its last step. On a conflict, take either side and rerun `cargo xtask docs --skill SKILL.md`.
- `Cargo.toml` and `Cargo.lock`: Task 1 only. A later conflict means an unplanned dependency; stop and ask.
- `crates/gloci/src/cli.rs`: Tasks 3–7 (A), 9–11 and 16–17 (B: new `Commands` variants, `SortArgs`), and 12, 13 and 15 (C). Append new structs at the end under `// P3a Task N` and new variants at the end of `Commands`; keep both sides.
- `crates/gloci/src/lib.rs`, `commands/mod.rs`: Tasks 9–11, 15 (`cat`) and Task 5 (`mod.rs`).
- `crates/gloci-ops/src/sort/*` and `commands/sort.rs`: Task 5 (`--sort` writes a temporary BAM), then Tasks 16–17, which start after it merges.
- `crates/gloci-vcf/src/concat.rs`: Task 15 only (the copy loop moves to gloci-io).
- In gloci-io: `lib.rs` and `align/mod.rs` (Tasks 2, 3, 4, 8, 12, 13); `record.rs` `dispatch` (2, 3, 12, 13); `genome.rs` (2, 8, 13). Elsewhere: `commands/index.rs` (6, 8, 13) and `convert/mod.rs` (4, 12). Keep both sides.
- `xtask/diff/cases.toml`: one appended `# P3a Task N` section per task. `docs/differences.md`: per-task sections (Tasks 3, 4, 5, 7, 9, 11, 15 and 17, and 12–14 if they pin one). `docs/migration.md`: Task 20 only.
- `crates/gloci/tests/golden/` and `xtask/diff/inputs/`: new paths, one owner each. `xtask/docker/Dockerfile`: Task 1. `Dockerfile.bench`, `xtask/bench/*` and `gen.rs`: Task 18.

## File structure

```
crates/gloci-io/src/align/{mod,header,record,bam,sam,writer,bed,cram}.rs   T1–T4, T13
crates/gloci-io/src/fasta/{mod,index,fetch,stream}.rs; {psl,rmsk}.rs; index/crai.rs   T8, T12, T13
crates/gloci-ops/src/{seq,getfasta,nuc,maskfasta}.rs; convert/{align,psl,rmsk}.rs   T4, T9–T12
crates/gloci-ops/src/{intersect,coverage,genomecov,view}.rs, relate/*   T5–T7, T14
crates/gloci-io/src/{align/source.rs, bgzf_copy.rs}; gloci-ops/src/{aln_cat,aln_sort}.rs, sort/order.rs   T3, T15–T17
crates/gloci/src/commands/{getfasta,nuc,maskfasta}.rs; cli.rs
xtask/diff/inputs/{aln,fa,convert}/; xtask/docker/*; xtask/bench/cases.toml
```

---

### Task 1: Foundation: dependencies, formats, detection, the input gate (D99, D100, D118)

**Files:**
- Modify:
  - root `Cargo.toml` (`noodles-sam = "0.93"`, `noodles-bam = "0.98"`, `noodles-cram = { version = "0.102", features = ["libdeflate"] }`, `noodles-fasta = "0.68"`), and the gloci-io and gloci-testkit manifests;
  - `gloci-io/src/{detect.rs, record.rs, output.rs}`;
  - every exhaustive `match` on `Format`/`OutFormat` (≈16 files: `grep -rn 'Format::Wig\|OutFormat::Tsv'`). New arms handle the variant or return a `UsageError`, never panic;
  - every `ReaderOptions { … }` literal (82: `grep -rn 'ReaderOptions {'`), each adding `accept: Accept::default()`;
  - `crates/gloci/src/cli.rs` (`AlignArgs`, not yet flattened);
  - `gloci-testkit/src/differential.rs`;
  - `xtask/docker/Dockerfile` (`samtools=1.21`).
- Create: `gloci-io/src/align/mod.rs`.
- Test: `gloci-io/tests/detect_p3.rs`, `gloci/tests/input_gate.rs`.

**Interfaces (produces):**
- `Format` gains `Bam, Cram, Sam, Fasta, Psl, Rmsk`, with `fn is_alignment(self) -> bool` and `fn is_convert_only(self) -> bool`.
- `OutFormat` gains `Bam, Sam, Fasta`.
- `#[derive(Debug, Clone, Default)] pub struct AlignOptions { pub min_mapq: u8, pub require_flags: u16, pub exclude_flags: u16, pub full_records: bool, pub reference: Option<PathBuf> }`
- `#[derive(Debug, Clone, Default)] pub struct Accept { pub alignments: Option<AlignOptions>, pub convert_only: bool }`, as `ReaderOptions.accept`.
- `pub fn parse_flags(s: &str) -> Result<u16, String>`
- `pub struct AlignArgs` (`--min-mapq`, `--require-flags`, `--exclude-flags`, `--reference`), with `fn options(&self, full_records: bool) -> Result<AlignOptions, UsageError>`.

Until a reader lands (Tasks 2, 3, 12, 13), an accepted format is an `IoError::Detect` naming it.

- [ ] **Step 1: Write the failing tests.**
  - `detect_p3_formats`, bytes → format:
    - `bgzip(b"BAM\x01…")` → Bam;
    - `b"CRAM\x03\x00"` and `b"CRAM\x02\x01"` → Cram;
    - `b"CRAM\x04\x00"` → Err containing `CRAM 4`;
    - `b"@HD\tVN:1.6\n"` → Sam;
    - `b"r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"` → Sam;
    - `b">chr1 d\nACGT\n"` and its gzip → Fasta;
    - a `psLayout` header and a bare PSL line → Psl;
    - a RepeatMasker header and a bare data line → Rmsk;
    - BED and VCF unchanged.
  - `input_gate`:
    - `gloci merge x.bam` → exit 2, with `does not read alignments`, `gloci convert x.bam -o x.bed` and `gloci sort -m`;
    - `gloci sort ref.fa` → exit 2, naming `getfasta, maskfasta or nuc`;
    - `gloci merge x.psl` → exit 2, with `gloci convert x.psl`;
    - `gloci sort a.bed -o o.cram` → exit 2, with `samtools view -C`.
  - `parse_flags_table`: `0x704`→0x704, `1796`→1796, `secondary,dup`→0x500, `UNMAP,MUNMAP`→0xC, `Supplementary`→0x800; `""` → Err; `bogus` → Err listing `PAIRED`.
  - `compare_bam_ignores_pg_and_compression`: one BAM at levels 1 and 9 with different `@PG` lines compares equal; one changed record byte does not.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test detect_p3 && cargo test -p gloci --test input_gate`. Expected: FAIL.
- [ ] **Step 3: Implement.** CRAM is checked on the raw bytes, and SAM before BED.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat(io): P3a formats, detection and the alignment/FASTA input gate`

### Task 2: BAM reader (D101–D103, D115, D119)

**Files:**
- Create: `gloci-io/src/align/{header.rs, record.rs, bam.rs}`, and the fixtures `xtask/diff/inputs/aln/{cigars,flags,order}.sam` with `make-alignments.sh`, the committed `.bam`/`.bam.bai`, and `targets.bed` (in `@SQ` order).
- Modify:
  - `record.rs`: `Header::Alignment(Arc<AlnHeader>)`, and `dispatch` → `BamReader`;
  - `genome.rs`: the table comes from the first header in input order, VCF or alignment;
  - `gloci-core/src/chrom.rs`: `OrderSource::AlignmentHeader(String)`;
  - `sorted.rs`: `sort_hint` (`-g {path}`, and `gloci sort x.bam -o x.sorted.bam` or `--sort` for alignment inputs);
  - `gloci-testkit/src/fixtures.rs`: `pub fn sam_to_bam(sam: &str) -> Vec<u8>`.
- Test: `gloci-io/tests/bam.rs`, and a fuzz body in `fuzz_bodies.rs` and `tests/fuzz_smoke.rs`.

**Interfaces (produces):**
- `pub struct AlnHeader { pub text: Vec<u8>, pub refs: Vec<(Vec<u8>, u64)>, pub sort_order: SortOrder }`, with `pub enum SortOrder { Coordinate, Queryname, Unsorted, Unknown }`.
- `impl AlnHeader`: `fn read_bam(r: &mut dyn Read) -> Result<Self, String>`, `fn parse_sam(text: Vec<u8>) -> Result<Self, String>`, `fn write_bam(&self, w: &mut dyn Write) -> io::Result<()>` and `fn renamed(&self, map: &ChromMap) -> Self`.
- `pub struct Aln<'a>`, a view of one record without `block_size`:
  - `fn parse(raw: &'a [u8], h: &AlnHeader) -> Result<Self, String>`;
  - `ref_id() -> Option<usize>`, `pos() -> Option<u64>`, `flag() -> u16`, `mapq() -> u8`, `name() -> &'a [u8]`;
  - `for_each_op(&self, f: &mut dyn FnMut(u8, u32))` (SAM letters, `CG`-resolved);
  - `int_tag(&self, tag: [u8; 2]) -> Option<i64>`;
  - `span(&self, ref_len: u64) -> Option<(Interval, bool /*clipped*/)>`;
  - `strand(&self) -> Strand`.
- `pub enum AlnRef<'a> { Bam(Aln<'a>) }` (Task 3 adds `Sam`), with the same methods, and `pub fn aln_ref<'a>(f: Format, raw: &'a [u8], h: &AlnHeader) -> Result<AlnRef<'a>, String>`.
- `pub fn blocks(r: &AlnRef<'_>, split_deletions: bool, out: &mut Vec<Interval>)` (D105)
- `pub struct BamReader: RecordReader` (format Bam, compression Bgzf, not drained; `loc` = record number in the file).

The data errors, each naming `record N of CHROM` and skipped and counted under `--skip-invalid`:
- `block_size` < 32 or past the data;
- an empty or unterminated name;
- `refID` outside `-1..n_ref`;
- `pos` < −1;
- a mapped read with `refID` or `pos` −1;
- a CIGAR op > 8;
- `l_seq` > 0 and different from the CIGAR's query length.

- [ ] **Step 1: Write the failing tests.**
  - `bam_spans_table`, (CIGAR, flag, pos, LN) → (span, strand, clipped):

    | CIGAR, flag | pos, LN | span | strand | clipped |
    |---|---|---|---|---|
    | `10M` | 100 | `[100,110)` | `+` | |
    | `5S10M5S`, `3H10M2H`, `5M3I5M`, `4=1X5=` | 100 | `[100,110)` | | |
    | `5M100N5M` | 100 | `[100,210)` | | |
    | `5M2D5M` | 100 | `[100,112)` | | |
    | `*` and `10S5I`, mapped | 100 | `[100,101)` | | |
    | `10M`, flag 16 | 100 | | `−` | |
    | `20M` | 16,560, LN 16,569 | `[16560,16569)` | | true |
  - `bam_long_cigar_cg_tag`: 70,000 ops of `1M1I` → span 35,000, and 70,000 ops seen.
  - `bam_unmapped_and_filters` on `flags.bam`:
    - default: every mapped read is kept (secondary, supplementary, duplicate and QC-fail included), and the 2 unmapped reads are skipped and counted;
    - `min_mapq 20` drops MAPQ 0 and 10;
    - `exclude_flags 0x900` drops secondary and supplementary;
    - `require_flags 0x2` keeps proper pairs.
  - `bam_header_is_an_order_source`: `order.bam` gives `chrM, chr1, chr2, chr10, chrX` from `AlignmentHeader`. A lex-sorted BED hints `gloci sort b.bed -g order.bam`. A queryname-sorted BAM hints `gloci sort` and names `SO:queryname`.
  - `bam_malformed_records`: each listed shape → an error naming `record N of chr1`.
  - `bam_huge_header_linear`: 200,000 `@SQ` with a read on every 1,000th. Open, read and build the table in < 10 s (debug), giving 200 records.
  - `bam_chrom_map_renames_header`: `chr1→1` renames the records and `@SQ\tSN:1`.
  - `bam_slim_records`: slim raw is shorter, with the same name, flag, MAPQ, ops, span and `NM`.
  - `fixtures_match_sources`.
  - `fuzz_bam`: 2,000 seeded flips and truncations never panic.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test bam`. Expected: FAIL.
- [ ] **Step 3: Implement.** Read through the decoded `Input`, so the BGZF pool applies. Build the reference-id → `ChromId` vector once.
- [ ] **Step 4: Run** `cargo test -p gloci-io -p gloci-core`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(io): BAM reader (CIGAR spans, unmapped and filters, header order)`

### Task 3: SAM input, BAM/SAM output, `AlnSource`, `view` of alignments (D104, D116, D120)

**Files:**
- Create: `align/{sam.rs, writer.rs, source.rs}`, and the fixture `aln/tags.sam` (tags `i Z A f B`, with `f` values `0.1` and `1e-07`).
- Modify:
  - `align/record.rs` (`AlnRef::Sam`), `record.rs`;
  - `writer.rs`: `check_family` for alignments, with `open_record_writer` delegating;
  - `gloci-ops/src/view.rs` and `gloci/src/commands/view.rs`: read through `AlnSource`, so every record is written; regions as D116; `--no-header` with BAM is a usage error;
  - `cli.rs`: `ViewArgs` flattens `AlignArgs`; example `gloci view reads.bam -O sam | head`;
  - SKILL.md, cases.toml, and differences.md if `#sam-float-tags` is pinned.
- Test: `gloci-io/tests/sam.rs`, `gloci/tests/view_alignments.rs`.

**Interfaces:**
- Consumes: Task 2.
- Produces:
  - `pub struct SamAln<'a>`, with `Aln`'s methods;
  - `pub struct SamReader` (raw = the line; an `RNAME` not in `@SQ` is a data error);
  - `impl AlnHeader { pub fn with_program(&self, command_line: &str, version: &str) -> AlnHeader }`;
  - `pub fn open_alignment_writer(out: Output, records: Format, to: OutFormat, h: &AlnHeader, o: &WriterOptions) -> Result<Box<dyn RecordWriter>, IoError>`. It writes:
    - BAM/CRAM → BAM: the raw bytes;
    - BAM/CRAM → SAM: each record read back with `noodles_bam::io::Reader::from` over a reused `block_size`-prefixed cursor, then `noodles_sam::io::Writer::write_alignment_record` with a `sam::Header` built from `refs` (gloci writes the header text itself);
    - SAM → SAM: the line;
    - SAM → BAM: a `noodles_sam::Record` written by `noodles_bam::io::Writer::from(Vec<u8>)`.
  - `pub fn sam_line_to_bam(line: &[u8], h: &AlnHeader, out: &mut Vec<u8>) -> Result<(), String>`, the same encoder the writer uses (Tasks 16 and 17 use it);
  - `pub trait AlnSource: Send` (D120). Tasks 13 and 15–17 and P3c consume it:
    - `fn name(&self) -> &str`;
    - `fn layout(&self) -> Format` (Bam for BAM and CRAM, Sam for SAM);
    - `fn header(&self) -> &Arc<AlnHeader>`;
    - `fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError>`: appends the next record (no `block_size`) and returns its record number;
    - `fn diagnostics(&self) -> &Diagnostics`.
  - `pub fn open_aln_source(path: &str, o: &ReaderOptions) -> Result<Box<dyn AlnSource>, IoError>` (BAM and SAM here; CRAM in Task 13).

- [ ] **Step 1: Write the failing tests.**
  - `aln_source_keeps_every_record`: on `flags.bam` and `flags.sam`, `AlnSource` yields every record, unmapped included, in file order; with `--exclude-flags UNMAP` the unmapped ones are dropped.
  - `sam_reader_equals_bam_reader`: `flags.sam` and `flags.bam` give identical metas and skip counts.
  - `sam_passthrough_bytes`: `view -O sam x.sam` copies each record line exactly.
  - `pg_line_ids`: an existing `@PG ID:gloci` gives `ID:gloci.1 … PP:gloci`.
  - `view x.bam -o y.bam` is identical at -t 1/2/8 (compare `bam`).
- [ ] **Step 2: Add golden cases:**
  - `view-bam-to-sam-vs-samtools`: `gloci view -O sam input.bam` against `samtools view -h --no-PG`;
  - `view-sam-to-bam-vs-samtools` (compare `bam`);
  - `view-bam-tags-vs-samtools`: match, or pin `#sam-float-tags`.
- [ ] **Step 3: Implement.** Run `cargo xtask diff --only view-`, then `--update-golden`.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat(io): SAM input, BAM/SAM output with @PG, view for alignments`

### Task 4: `convert` for alignments (bamtobed) (D106)

**Files:**
- Create: `align/bed.rs`, `gloci-ops/src/convert/align.rs`.
- Modify:
  - `convert/mod.rs`;
  - `cli.rs`: `ConvertArgs` gets `--bed12`, `--color`, `--score-tag`, `--cigar`, `--split-deletions` and `AlignArgs`; example `gloci convert reads.bam --bed12 -o reads.bed`;
  - `commands/convert.rs`: slim records keep the score tag;
  - SKILL.md, cases.toml, differences.md.
- Test: `gloci/tests/convert_alignments.rs`.

**Interfaces:**
- Consumes: `aln_ref`, `blocks` (Task 2), `SamAln` (Task 3).
- Produces:
  - `pub enum AlnBedLayout { Bed6, Bed12, Blocks }`
  - `pub struct AlnBedOptions { pub layout: AlnBedLayout, pub score_tag: Option<[u8; 2]>, pub cigar: bool, pub split_deletions: bool, pub color: String }` (Default: Bed6, MAPQ, no CIGAR, N-only, `"255,0,0"`)
  - `pub struct AlnBed`, with:
    - `fn for_reader(r: &dyn RecordReader, o: AlnBedOptions) -> Option<AlnBed>`;
    - `fn write(&mut self, raw: &[u8], loc: u64, out: &mut Vec<u8>) -> Result<(), IoError>`;
    - `fn bed6<'a>(&'a mut self, raw: &[u8], loc: u64) -> Result<&'a [u8], IoError>`.

- [ ] **Step 1: Write the failing tests:**
  - usage errors: `--color` without `--bed12`, and `--cigar` with `--bed12`/`--blocks`;
  - a non-integer `--score-tag` type is a data error;
  - `--write-index` on BED from a sorted BAM;
  - determinism.
- [ ] **Step 2: Add golden cases** (`bedtools bamtobed -i input.bam …` against `gloci convert input.bam …`):
  - default on `flags.bam`;
  - `-bed12` and `-split` (`--blocks`) on `cigars.bam`;
  - `-splitD`, and `-bed12 -splitD`;
  - `-ed` (`--score-tag NM`), `-tag AS`, `-cigar`;
  - `-bed12 -color 0,128,255`;
  - `order.bam` (ties keep file order);
  - SAM input (reference `samtools view -b input.sam | bedtools bamtobed -i stdin`);
  - `--min-mapq 30 --exclude-flags SECONDARY` (reference `samtools view -b -q 30 -F 0x100 … | bedtools bamtobed -i stdin`);
  - an absent `--score-tag XY` (compare `exit-status`);
  - `-split -tag NM`, pinned `#bamtobed-split-tag`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: convert BAM/CRAM/SAM to BED (bamtobed, -bed12, -split, -splitD, -tag)`

### Task 5: Alignments in `intersect` and `coverage` (D104, D105, D115, D124)

**Files:**
- Modify:
  - `gloci-ops/src/relate/criteria.rs`: block source per input (BED12, alignment, none);
  - `relate/mod.rs`: `RelateOptions.split_deletions: bool`;
  - `relate/files.rs`: `b_width` 6 for alignments;
  - `intersect.rs`: `Text` renders alignments via `AlnBed::bed6`; `validate` refuses clip of an alignment and `--split-deletions` without one; passthrough through `open_record_writer`;
  - `coverage.rs`: alignment A blocks;
  - `cli.rs`: `IntersectArgs` and `CoverageArgs` get `--split-deletions` (not `OverlapArgs`, whose other commands read no alignments) and `AlignArgs`; examples `gloci intersect -a reads.bam -b targets.bed -o on_target.bam` and `gloci coverage -a targets.bed -b reads.bam --split`;
  - `commands/{intersect,coverage,mod}.rs`: full A only in `--mode any|none`, slim B;
  - `gloci-ops/src/sort/inputs.rs` (`sort_readers`): `--sort` of an alignment input writes a level-1 temporary BAM through `open_alignment_writer` and reopens it with `BamReader` (D124);
  - SKILL.md, cases.toml, differences.md.
- Test: `gloci/tests/relate_alignments.rs`.

**Interfaces:**
- Consumes: `AlnBed::bed6` (Task 4), `blocks` (Task 2), `open_alignment_writer` (Task 3).
- Produces: `RelateOptions.split_deletions`.

- [ ] **Step 1: Write the failing tests.**
  - `alignment_relations_deterministic`: intersect (any/none/count/pairs) and coverage (summary/hist/per-base), with BAM as A and as B, at -t 1/2/8.
  - `clip_refused_for_alignments`, `split_deletions_needs_alignment`: each exits 2.
  - `sort_flag_for_alignments` (D124): `intersect --sort -a shuffled.bam -b targets.bed` and `coverage --sort -a targets.bed -b shuffled.bam` equal the runs on the sorted BAM, within `--max-mem 16M`; the temp directory is gone afterwards. Without `--sort`, the order error hints `gloci sort shuffled.bam -o shuffled.sorted.bam` or `--sort`.
  - `bed_output_for_alignment_a_refused`: exits 2 naming the `gloci convert -` pipe.
  - `slim_b_records`: 10,000 B reads of 150 bp at one locus; peak active-set raw bytes ≤ ⅓ of full records.
  - `unsorted_alignment_b`: `--unsorted` with a shuffled SAM B equals sorted streaming, in A's order.
  - A hand-checked `--split-deletions` coverage case.
- [ ] **Step 2: Add golden cases** (`order.bam`; the targets are bookended to a read, inside an intron, over only a deletion, and a point):
  - `intersect -abam input.bam -b targets.bed -u` against `gloci intersect -a input.bam -b targets.bed -O bam` (compare `bam`); and `-u -split`, `-u -s`, `-u -S`, `-u -f 0.5`;
  - `-v` without unmapped reads (match), and on `flags.bam` (pinned `#unmapped-reads`);
  - `-a targets.bed -b input.bam -c` and `-c -split`;
  - `-abam … -bed -wa -wb`, and `-a targets.bed -b input.bam -wa -wb` (pinned `#bam-as-bed-columns`);
  - `coverage -a targets.bed -b input.bam`, with `-split`, `-hist`, `-mean` and `-d`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: BAM/SAM in intersect and coverage (--split on CIGAR N, --split-deletions)`

### Task 6: Alignment indexes, regions, per-chromosome parallelism (D107)

**Files:**
- Modify:
  - `gloci-io/src/index/mod.rs`: `IndexKind::Bai`, and `find_index` for alignments;
  - `index/build.rs`: `Parser::Bam`, BAI via `noodles_bam::bai`, CSI via `--csi` or LN > 2^29, metadata and `n_no_coor`;
  - `index/query.rs`: `Loaded::Bai`, names from `AlnHeader`, `rid_of` through a `HashMap`, `RegionReader` BAM frames;
  - `gloci-engine/src/parallel.rs`: `all_indexed`;
  - `gloci-ops/src/relate/indexed.rs`: the disjoint-names warning and the `--region` chunk filter switch to `HashSet`s;
  - `relate/files.rs`;
  - `view.rs`: alignment regions are read through the index as an `AlnSource` that keeps placed unmapped reads (D116). The `RecordReader` region path keeps skipping them;
  - the output path (`--write-index` for BAM);
  - `cli.rs`/`commands/index.rs`: `--csi`;
  - SKILL.md, cases.toml.
- Test: `gloci-io/tests/bam_index.rs`, `gloci/tests/parallel_alignments.rs`.

**Interfaces:**
- Consumes: Tasks 2 and 3.
- Produces: `IndexedFile` works for BAM (`Header::Alignment`), with `chromosomes_reader`, `regions_reader` and `weight`. Tasks 7 and 14 use them.

- [ ] **Step 1: Write the failing tests.**
  - `bai_queries_match_samtools_bai`: on `cigars`, `flags` and `order.bam` with their samtools `.bai`, 300 seeded regions (random; read ends ±1; bin edges 16,384 and 131,072; zero-reference-length reads). Reads through samtools' `.bai` = through gloci's = a §4.3-filtered scan.
  - `bai_metadata`: the mapped/unmapped counts and `n_no_coor` equal samtools'.
  - `csi_for_long_references`: LN 600,000,000 → `.csi`; `--csi` forces it.
  - `bam_per_chromosome_equals_streaming`: intersect (any/count) and coverage (summary/hist), with an indexed BAM as A and as B, at -t 1/2/8, and with `--region chr2:100-500`.
  - `indexed_bam_b_is_not_drained`, mirroring `indexed_b_is_not_drained`.
  - `bam_error_same_at_every_t`: a corrupt record on an A chromosome of an indexed B.
  - `huge_header_per_chromosome`: 200,000 `@SQ`, with reads and a BED A on 20,000; -t 8 finishes in < 20 s (debug) and equals -t 1.
  - `write_index_on_bam_output`.
- [ ] **Step 2: Add cases.** Golden `view-bam-region-vs-samtools`: `gloci view -O sam --no-header input.bam chr1:1-150 chr2:95-105` against `samtools index input.bam && samtools view input.bam chr1:1-150 chr2:95-105` (one region per chromosome, so no read repeats). Live: `index-interop-bam-samtools-view` and `index-interop-bam-csi-samtools`.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace` and `cargo xtask diff --only index-interop-bam`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: BAM .bai/.csi indexes, region queries and per-chromosome parallelism`

### Task 7: `genomecov` on alignments (D105, D117, D124)

**Files:**
- Modify:
  - `gloci-ops/src/genomecov.rs`;
  - `cli.rs`: `GenomecovArgs`: `-g` required unless the input is an alignment; `--split-deletions`; `AlignArgs`; example `gloci genomecov reads.bam --bga --split`;
  - `commands/genomecov.rs`;
  - SKILL.md, cases.toml, differences.md.
- Test: `gloci/tests/genomecov_alignments.rs`.

**Interfaces:**
- Consumes: `blocks` (Task 2), the BAM `IndexedFile` (Task 6), `run_chunks`.
- Produces: `GenomecovOptions.split_deletions: bool`.

- [ ] **Step 1: Write the failing tests.**
  - `genomecov_alignment_parallel_equals_streaming`: hist, `--bga` and `--per-base` at -t 1/2/8.
  - `genomecov_bam_events_bounded`: 100,000 reads of 100 bp at depth 30; peak pending events ≤ 120.
  - `genomecov_bam_no_genome_needed`: without `-g`, header lengths; with `-g`, the given genome.
  - `genomecov_sort_flag_bam`: `--sort` on a shuffled BAM equals the run on the sorted BAM (D124).
- [ ] **Step 2: Add golden cases** against `bedtools genomecov -ibam input.bam`:
  - default against gloci's default (pinned `#bam-genomecov-deletions`), and against `--split-deletions` (match);
  - `-bg`, `-bga`, `-d` and `-split`, each with `--split-deletions`;
  - `-ignoreD` (gloci's default), and `-split -ignoreD` (`--split`);
  - `-strand -`, `-5`, `-3`, `-max 3`, `-scale 0.5 -bg`;
  - `genomecov-ibam-overhang-vs-bedtools`;
  - `order.bam` (header order, empty chromosomes).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: genomecov on BAM/SAM (header genome, --split-deletions, per-chromosome)`

### Task 8: FASTA indexes and random access (D108)

**Files:**
- Create:
  - `gloci-io/src/fasta/{mod.rs, index.rs, fetch.rs, stream.rs}`;
  - fixtures `xtask/diff/inputs/fa/`:
    - `ref.fa`: chrM, chr1, chr2, chr10, chrX; 60-column lines and a short last line; lowercase, `N` and IUPAC runs; chr1:100-129 holds 27 A/T of 29;
    - `ref.fa.gz`, multi-block (a 300 kb chr2);
    - `crlf.fa`, `desc.fa`, `empty-seq.fa`, `irregular.fa`;
    - `make-fasta.sh`, and the samtools `.fai`/`.gzi`.
- Modify:
  - `commands/index.rs`: FASTA → `.fai` (+ `.gzi`);
  - `cli.rs`: `IndexArgs` help; example `gloci index ref.fa.gz`;
  - `genome.rs`: `-g ref.fa`;
  - `record.rs`: the FASTA gate message now names `-g ref.fa` as well;
  - SKILL.md, cases.toml.
- Test: `gloci-io/tests/fasta.rs`.

**Interfaces (produces):**
- `pub struct FaiEntry { pub name: Vec<u8>, pub len: u64, pub offset: u64, pub line_bases: u64, pub line_width: u64 }`
- `pub fn build_fai(src: &mut dyn BufRead, name: &str) -> Result<Vec<FaiEntry>, IoError>`
- `pub fn write_fasta_index(path: &Path, threads: NonZeroUsize) -> Result<Vec<PathBuf>, IoError>` (atomic)
- `pub struct IndexedFasta`, with:
  - `open(path: &Path) -> Result<Self, IoError>`, `name()`, `entries() -> &[FaiEntry]`, `len_of(chrom: &[u8]) -> Option<u64>`, `reopen(&self) -> Result<Self, IoError>` (shares the index);
  - `fetch(&mut self, chrom: &[u8], iv: Interval, f: &mut dyn FnMut(&[u8]) -> Result<(), IoError>) -> Result<(), IoError>`: chunks ≤ 1 MiB, in order, newline-free, case kept;
  - `fetch_rev(...)`: the same chunks, last first.
- `pub enum FastaLine<'a> { Header(&'a [u8]), Seq(&'a [u8]) }`
- `pub struct FastaStream` (plain, gzip or BGZF via `open_input`), with `open(path: &str, threads: NonZeroUsize)`, `next(&mut self) -> Result<Option<FastaLine<'_>>, IoError>`, `line() -> u64` and `crlf() -> bool`.

- [ ] **Step 1: Write the failing tests.**
  - `fai_matches_samtools`: `.fai` for `ref`, `crlf`, `desc`, `empty-seq` and `ref.fa.gz`, and `.gzi` for `ref.fa.gz`, each byte-equal to the committed samtools/htslib files.
  - `fai_errors`:
    - `irregular.fa` → exit 1 naming the line;
    - a duplicate name → exit 1;
    - no `.fai` → exit 2 with `gloci index ref.fa`;
    - BGZF without `.gzi` → exit 2 with `gloci index ref.fa.gz`;
    - plain gzip → exit 2 with `bgzip`;
    - an edited FASTA → `does not match` on fetch;
    - an old `.fai` → a stale warning.
  - `fetch_equals_naive`: 1,000 seeded intervals (line and BGZF block boundaries, whole sequences, zero-length). The chunks, and `fetch_rev`'s chunks last first, equal naive slicing. A 3 MiB sequence gives ≥ 3 chunks, none over 1 MiB.
  - `genome_from_fasta`: `-g ref.fa` equals `-g ref.fa.fai`; without a `.fai`, it names `gloci index`.
  - `fuzz_fasta`: mutated FASTA (through `build_fai`) and mutated `.fai` lines (through `open` and `fetch`) never panic.
- [ ] **Step 2: Add live cases** `index-interop-faidx` and `index-interop-faidx-bgzf`: gloci's index, then `samtools faidx ref.fa chr1:5-70 chr2:60000-60100`, against samtools' own index.
- [ ] **Step 3: Implement** (`.gzi` lookups via `noodles_bgzf::gzi`).
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat(io): FASTA .fai/.gzi (samtools-identical) and chunked random access`

### Task 9: `getfasta` (D111)

**Files:**
- Create: `gloci-ops/src/{seq.rs, getfasta.rs}`, `gloci/src/commands/getfasta.rs`, and fixtures `fa/{gf.bed, gf12.bed, gf.gff3, gf.vcf}` (strands `+`/`-`/`.`; lowercase and IUPAC under `-`; a 3-block BED12).
- Modify:
  - `cli.rs`: `Commands::Getfasta(GetfastaArgs)`: `--fasta`, `FILE`, `-s`, `--name`/`--name-only`, `--split`, `--tab`/`--bed-out`, `--rna`, plus output, chrom-map, vcf-span, format, invalid and threads args; 3 examples;
  - `lib.rs`, `commands/mod.rs`;
  - SKILL.md, cases.toml, differences.md.
- Test: `gloci-ops/tests/seq.rs`, `gloci/tests/getfasta.rs`.

**Interfaces:**
- Consumes: `IndexedFasta` (Task 8).
- Produces (Task 10 reuses `seq.rs` and `SeqSummary`):
  - `pub fn complement(b: u8, rna: bool) -> u8`; `pub fn revcomp_in_place(seq: &mut [u8], rna: bool)`
  - `pub enum FastaName { Coords, NameCoords, NameOnly }`; `pub enum SeqLayout { Fasta, Tab, BedOut }`
  - `pub struct GetfastaOptions { pub strand: bool, pub name: FastaName, pub split: bool, pub layout: SeqLayout, pub rna: bool }`
  - `pub struct SeqSummary { pub written: u64, pub points: u64, pub skipped: u64 }`
  - `pub fn getfasta(input: Box<dyn RecordReader>, fasta: &mut IndexedFasta, chroms: &mut ChromTable, o: &GetfastaOptions, out: &mut dyn Write) -> Result<SeqSummary, IoError>`

- [ ] **Step 1: Write the failing tests.**
  - `revcomp_table`, over bytes 0–255:
    - `A↔T`, `C↔G`, `U→A`, `R↔Y`, `K↔M`, `B↔V`, `D↔H`;
    - `S`, `W`, `N`, `X` unchanged;
    - lowercase alike; every other byte unchanged;
    - `rna`: `A→U`, `a→u`.
  - `getfasta_streams_large_records`: a 5 MiB record on `+` and `−` equals the naive answer, holding ≤ 1 MiB plus one header.
  - Points warn with a count.
  - `-o out.fa.gz` is BGZF.
  - `--chrom-map` renames the intervals.
  - Determinism.
- [ ] **Step 2: Add golden cases** (`bedtools getfasta -fi ref.fa -bed input.bed`):
  - default, `-name`, `-nameOnly`;
  - `getfasta-strand-iupac-vs-bedtools` (`-s`), and `-s -name`;
  - `-tab`, `-bedOut`;
  - `-split` and `-split -s`;
  - `-rna -s`;
  - GFF3;
  - VCF with `--vcf-span ref`;
  - BGZF FASTA;
  - out-of-FASTA records: match with `--skip-invalid`; by default, compare `exit-status`, pinned `#getfasta-outside-fasta`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci getfasta`

### Task 10: `nuc` (D112)

**Files:**
- Create: `gloci-ops/src/nuc.rs`, `gloci/src/commands/nuc.rs`.
- Modify: `cli.rs` (`Commands::Nuc(NucArgs)`: `--fasta`, `FILE`, `-s`, `--seq`, `--pattern`, `--ignore-case`, and the common args; 3 examples), `lib.rs`, `commands/mod.rs`, SKILL.md, cases.toml.
- Test: `gloci-ops/tests/nuc.rs`, `gloci/tests/nuc.rs`.

**Interfaces:**
- Consumes: Tasks 8 and 9.
- Produces:
  - `pub struct NucOptions { pub strand: bool, pub seq: bool, pub pattern: Option<Vec<u8>>, pub ignore_case: bool }`
  - `pub struct BaseCounts { pub a: u64, pub c: u64, pub g: u64, pub t: u64, pub n: u64, pub other: u64 }`, with `fn add(&mut self, chunk: &[u8])`
  - `pub fn pct_f32(count: u64, len: u64) -> String`
  - `pub fn nuc(input: Box<dyn RecordReader>, fasta: &mut IndexedFasta, chroms: &mut ChromTable, o: &NucOptions, out: &mut dyn Write) -> Result<SeqSummary, IoError>`

  Pattern counting streams with a (len−1)-byte carry. BCF columns render through `Text`.

- [ ] **Step 1: Write the failing tests.**
  - `pct_f32_table`: (27,29)→`0.931035`, (39,41)→`0.951219`, (16,51)→`0.313726`, (1,3)→`0.333333`, (0,10)→`0.000000`, (10,10)→`1.000000`.
  - `base_counts_table`: `ACGTUNacgtunRYKM-*` → a 2, c 2, g 2, t 4, n 2, other 6.
  - `pattern_across_chunks`: counted once across a 1 MiB boundary.
  - Points warn; determinism.
- [ ] **Step 2: Add golden cases** (`bedtools nuc -fi ref.fa -bed input.bed`):
  - BED3 and BED6;
  - `-s`, `-seq`, `-pattern GC`, `-pattern gc -C`, `-s -seq -pattern ACG`;
  - `nuc-lowercase-iupac-vs-bedtools`;
  - `nuc-f32-rounding-vs-bedtools` (chr1:100-129);
  - VCF with `--vcf-span ref`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci nuc`

### Task 11: `maskfasta` (D113)

**Files:**
- Create: `gloci-ops/src/maskfasta.rs`, `gloci/src/commands/maskfasta.rs`, fixture `fa/masks.bed` (overlapping masks, one across line breaks, one on an absent chromosome).
- Modify:
  - `cli.rs`: `Commands::Maskfasta(MaskfastaArgs)`: `--fasta`, `-b FILE…` (required), `--soft`/`--mask-char`, output (FASTA, `--write-index`), chrom-map, vcf-span and invalid args; 3 examples;
  - the output path: `OutFormat::Fasta`, with `--write-index` → `write_fasta_index`;
  - `lib.rs`, `commands/mod.rs`;
  - SKILL.md, cases.toml, differences.md.
- Test: `gloci/tests/maskfasta.rs`.

**Interfaces:**
- Consumes: `FastaStream`, `write_fasta_index` (Task 8).
- Produces:
  - `pub enum Mask { Hard(u8), Soft }`
  - `pub struct MaskSet`, with `fn load(readers: Vec<Box<dyn RecordReader>>, chroms: &mut ChromTable) -> Result<MaskSet, IoError>`
  - `pub struct MaskSummary { pub masked_bases: u64, pub missing_chroms: Vec<Vec<u8>> }`
  - `pub fn maskfasta(fasta: FastaStream, masks: &MaskSet, mask: Mask, out: &mut dyn Write) -> Result<MaskSummary, IoError>`

- [ ] **Step 1: Write the failing tests.**
  - `maskfasta_streams`: a 50 MiB FASTA holds ≤ one line.
  - `mask_char_validation`: `NN` or a non-printable character → exit 2.
  - `soft_mask_keeps_non_letters`.
  - `-o out.fa.gz --write-index` → a `.fai`/`.gzi` equal to `gloci index` of the output.
  - Determinism.
- [ ] **Step 2: Add golden cases** (`bedtools maskfasta -fi ref.fa -bed input.bed -fo /dev/stdout`):
  - default, `-soft`, `-mc X`;
  - masks on an absent chromosome (match);
  - `desc.fa` by default (pinned `#maskfasta-header`) and with `-fullHeader` (match);
  - `irregular.fa` and `crlf.fa` (match, or pinned `#maskfasta-layout`);
  - a mask past the end (compare `exit-status`, pinned `#maskfasta-past-end`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci maskfasta (streamed, layout kept)`

### Task 12: PSL and RepeatMasker in `convert` (D114)

**Files:**
- Create: `gloci-io/src/{psl.rs, rmsk.rs}` (`TextKind`), `gloci-ops/src/convert/{psl.rs, rmsk.rs}`, fixtures `xtask/diff/inputs/convert/{blat.psl, blat-header.psl, repeats.rmsk.out}`.
- Modify: `record.rs` (dispatch, `convert_only`), `convert/mod.rs`, `commands/convert.rs`, the convert help (example `gloci convert hits.psl -o hits.bed`), SKILL.md, cases.toml.
- Test: `gloci/tests/convert_psl_rmsk.rs`.

- [ ] **Step 1: Write the failing tests.**
  - Spans and strands on both formats.
  - The headers (5 and 3 lines) are skipped.
  - A malformed line → exit 1 naming the line.
  - Out-of-order input warns.
  - `merge x.psl` is refused (Task 1's gate).
- [ ] **Step 2: Add golden cases** against `convert2bed --input=psl --do-not-sort` (headered and bare) and `convert2bed --input=rmsk --do-not-sort`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: convert PSL and RepeatMasker .out to BED`

### Task 13: CRAM input, the reference, and `.crai` (D109, D110)

**Files:**
- Create:
  - `gloci-io/src/align/cram.rs`: `CramReader`, the pool, transcoding, `RefAdapter`;
  - `gloci-io/src/index/crai.rs`;
  - `aln/make-cram.sh`, and CRAMs made by samtools 1.21 from the three SAMs with `fa/ref.fa`: 3.0, 3.1, multi-reference slices (`seqs_per_slice=2,multi_seq_per_slice=1`), `embed_ref=1`, `no_ref=1`, plus their `.crai`.
- Modify: `record.rs`, `genome.rs` (`-g x.cram`), `index/mod.rs`, `commands/index.rs` and `cli.rs` (`IndexArgs --reference`), SKILL.md (reference and memory pitfalls).
- Test: `gloci-io/tests/cram.rs`, `gloci/tests/cram_reference.rs`.

**Interfaces:**
- Consumes: `IndexedFasta` (Task 8); the `BamReader` parse/filter/slim path (Task 2), which every transcoded record takes; `view` (Task 3) for CLI tests.
- Produces:
  - `pub struct CramReader: RecordReader` (format Cram, compression None; raw records are BAM records). CRAM also implements `AlnSource` (Task 3) for every record, and `open_aln_source` dispatches to it;
  - `pub(crate) struct RefAdapter`, which implements `noodles_fasta::repository::Adapter` and returns `Some(Err)`, never `None`;
  - `pub fn write_crai(path: &Path, reference: Option<&Path>, threads: NonZeroUsize) -> Result<PathBuf, IoError>`;
  - `pub struct CramIndex`: entries of reference, start, span, container offset, slice offset and size. Task 14 uses it.

- [ ] **Step 1: Write the failing tests.**
  - `cram_equals_bam`: each CRAM yields the BAM fixture's metas and records (name, flag, ref, pos, MAPQ, CIGAR, tags, SEQ, QUAL). Through `AlnSource`, it yields every record, unmapped included, equal to the BAM's.
  - `cram_missing_reference`: exit 1 (not 101 or 134), with `needs the reference sequence of chr1` and `--reference`.
  - `cram_mismatched_reference`: a base changed in a used region → `does not match`; a changed `LN` → `was not made with`, before any record.
  - `cram_reference_bgzf`: `ref.fa.gz` with `.fai`/`.gzi` works; without `.gzi` it gives Task 8's error.
  - `cram_no_ref_needed`: the `no_ref` and `embed_ref` fixtures decode with no reference.
  - `cram_reference_cache_bounded`: a counting adapter sees ≤ 1 live sequence per decoder.
  - `cram_truncated`: cuts at every container boundary and mid-container → exit 1 naming the file.
  - `cram_pool_deterministic`: -t 1/2/8.
  - `crai_matches_samtools`: entry sets equal.
  - `index_cram_cli`: the multi-reference fixture without `--reference` → exit 2 with the fix; with it, `.crai` is written.
  - `-g x.cram` needs no reference; `--reference` with no CRAM input warns.
- [ ] **Step 2: Add the live case** `index-interop-cram-samtools-view` (gloci's `.crai`, then `samtools view -T ref.fa` regions).
- [ ] **Step 3: Implement.** Check slice reference ids against the header before decoding. Map noodles' checksum `InvalidData` and the adapter's errors to the D109 messages. Clear each decoder's repository when the reference id changes.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat(io): CRAM input (--reference, bounded cache, decode pool) and .crai`

### Task 14: CRAM in every alignment command, with regions and per-chromosome runs

**Files:** `align/cram.rs` (`CramFile`: names, `weight` from `.crai`, `chromosomes_reader`, `regions_reader`), `relate/indexed.rs` (`Src::Cram`), `relate/files.rs`, `gloci-engine/src/parallel.rs`, `view.rs`, `genomecov.rs`, cases.toml. Test: `gloci/tests/parallel_cram.rs`.

**Interfaces:**
- Consumes: `CramIndex` (Task 13), the chunk driver (Tasks 6 and 7), `RelateOptions.split_deletions` (Task 5).

- [ ] **Step 1: Write the failing tests.**
  - `cram_per_chromosome_equals_streaming`: intersect, coverage, genomecov and view with an indexed CRAM, at -t 1/2/8.
  - `cram_worker_reference_memory`: ≤ threads live sequences.
  - `cram_error_same_at_every_t`.
- [ ] **Step 2: Add golden cases.** gloci reads `input.cram` with `--reference ref.fa`; the reference command decodes it with `samtools view -b -T ref.fa` first:
  - `convert` (against bamtobed);
  - `intersect -a targets.bed -b input.cram --mode count`;
  - `coverage`;
  - `genomecov --split-deletions`;
  - `view-cram-region-vs-samtools`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: CRAM regions and per-chromosome parallelism in every alignment command`

### Task 15: `gloci cat` for BAM (D123)

**Files:**
- Create:
  - `gloci-io/src/bgzf_copy.rs`, `gloci-ops/src/aln_cat.rs`, `gloci/src/commands/cat.rs`;
  - fixtures `aln/cat-{a,b,c}.sam` with their BAMs (`make-alignments.sh`): `cat-b` adds an `@RG`, and `cat-c` has a different `@SQ` list.
- Modify:
  - `gloci-vcf/src/concat.rs`: its block-copy loop moves to `bgzf_copy`. `vcf concat` behaviour and tests are unchanged;
  - `cli.rs`: `Commands::Cat(CatArgs)`, with `FILES…`, `-o`, `--write-index` and threads; 2 examples;
  - `lib.rs`, `commands/mod.rs`;
  - SKILL.md, cases.toml, differences.md (`#alignment-cat-sq`).
- Test: `gloci-io/tests/bgzf_copy.rs`, `gloci/tests/cat.rs`.

**Interfaces:**
- Consumes: `AlnHeader::{read_bam, write_bam, with_program}` (Tasks 2, 3).
- Produces:
  - `pub fn copy_from(src: &mut dyn Read, out: &mut dyn Write, start: u64) -> io::Result<CopyStats>`: reads the BGZF stream sequentially (stdin works). The block holding uncompressed offset `start` is recompressed from there, later blocks are copied verbatim, and EOF markers are dropped. `vcf concat --naive` calls it;
  - `pub fn cat_header(inputs: &[(&str, &AlnHeader)]) -> Result<AlnHeader, IoError>`;
  - `pub fn cat_bams(paths: &[String], out: Output, o: &WriterOptions) -> Result<(), IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `copy_from_mid_block_and_edge`: a header ending mid-block recompresses only that block; on a block edge, everything after it is copied verbatim. The two existing `vcf concat` block-edge tests still pass.
  - `cat_header_rules`:
    - the extra `@RG` is appended;
    - different `@SQ` → exit 1, naming the first difference;
    - CRAM → exit 2 suggesting `samtools cat`;
    - SAM → exit 2 with the `gloci view` hint.
  - `cat_memory_one_block`: the largest buffer held is ≤ one BGZF block plus the header.
  - `--write-index` on joined sorted pieces gives a `.bai`.
  - Determinism.
- [ ] **Step 2: Add golden cases** against `samtools cat --no-PG` (compare `bam`):
  - `gloci cat input.bam cat-b.bam` (with the extra `@RG`);
  - three inputs;
  - `cat-c` (compare `exit-status`, pinned `#alignment-cat-sq`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci cat (BAM block copy, as samtools cat)`

### Task 16: `gloci sort` for alignments (D121)

**Files:**
- Create:
  - `gloci-ops/src/sort/order.rs`, `gloci-ops/src/aln_sort.rs`;
  - fixtures `aln/sort.sam` (with its BAM), `aln/sort-nohd.sam` (no `@HD`) and `aln/sort-hd.sam` (`@HD … SO:unsorted GO:query SS:x`). `sort.sam` has: reads at chr1:100 on `-` then `+`; `+` reads there named `read10` then `read9`; one at chr1:99; chrM reads after chr1's (`@SQ` has chrM first); a placed unmapped read at chr1:100; interleaved unplaced reads with and without `0x10`; the names `r2 r10 r1 r01 a1b2 a1b10 R1 r1a`; and one name with READ1, READ2, secondary and supplementary records.
- Modify:
  - `gloci-ops/src/sort/{mod.rs, chunk.rs, merge.rs, spill.rs}`: `SortConfig.order: Order` replaces `tie`; `Sorter::push_raw`; the comparators dispatch on `Order`; interval sorting is unchanged;
  - `align/header.rs`: `AlnSortOrder`, `sub_sort`, `with_sort_order`;
  - `gloci/src/commands/sort.rs`: dispatch for alignments;
  - `cli.rs`: `SortArgs` gets `-n/--by-name` and `-N/--by-name-bytes`; help and the example `gloci sort reads.bam -o reads.sorted.bam --write-index -t 8`;
  - SKILL.md, cases.toml.
- Test: `gloci-ops/tests/aln_sort.rs`, `gloci/tests/sort_alignments.rs`.

**Interfaces:**
- Consumes: `AlnSource`, `open_aln_source`, `sam_line_to_bam` and `open_alignment_writer` (Task 3); BAM `--write-index` (Task 6).
- Produces (Task 17 uses all of them):
  - `pub enum AlnSortOrder { Coordinate, Name, NameBytes }`; `AlnHeader.sub_sort: Option<Vec<u8>>`; `impl AlnHeader { pub fn with_sort_order(&self, o: AlnSortOrder) -> AlnHeader }`;
  - `pub enum Order { Interval(Tie), Aln(AlnSortOrder) }`, as `SortConfig.order`;
  - `pub fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering`; `pub fn name_flag_rank(flag: u16) -> u32`; `pub fn coordinate_key(raw: &[u8], n_ref: u32) -> (u32, u64)` (`(tid, or n_ref for -1; ((pos+1) << 1) | reverse)`);
  - `impl Sorter { pub fn push_raw(&mut self, raw: &[u8]) -> Result<(), IoError> }` (the input sequence number is the last tie);
  - `pub struct AlnSortOptions { pub order: AlnSortOrder, pub max_mem: u64, pub tmp_dir: PathBuf, pub threads: NonZeroUsize }`;
  - `pub fn sort_alignments(src: Box<dyn AlnSource>, o: &AlnSortOptions, out: Output, to: OutFormat, w: &WriterOptions) -> Result<SortStats, IoError>`;
  - `pub fn check_alignments(src: Box<dyn AlnSource>, order: AlnSortOrder) -> Result<bool, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `strnum_cmp_matches_samtools`:
    - natural: `r2` < `r10`; `a1b2` < `a1b10`; `r01` = `r1`; `007` = `7`; `r1` < `r1a`; `R1` < `r1`; `a10` > `a9b`; `x` < `x0`;
    - bytes (`-N`): `r10` < `r2`.
  - `name_flag_rank_table`: none < supplementary < secondary < READ1 < READ1+supplementary < READ1+secondary < READ2.
  - `aln_sort_huge_record`: a 40 MiB record under `--max-mem 16M` sorts in a chunk of its own, equal to the in-memory result.
  - `aln_sort_spills_equal_in_memory`: 200,000 records at `--max-mem 16M` (≥ 4 spills) equal `--max-mem 2G`, for all three orders, at -t 1/2/8.
  - `aln_sort_header`: the D121 `@HD` rules, with and without `@HD`.
  - `aln_sort_cli_errors`, each exit 2:
    - two inputs (suggesting `sort -m` or `cat`);
    - `-n` on a BED;
    - `--write-index` with `-n`;
    - `-g` with a BAM;
    - `-o x.cram`.
  - `aln_sort_check`: `--check` exits 0 when sorted, and 1 naming both records when not.
- [ ] **Step 2: Add golden cases** against `samtools sort --no-PG` (compare `bam`):
  - `aln-sort-ties-vs-samtools` (coordinate, on `sort.bam`);
  - `aln-sort-name-vs-samtools` (`-n`) and `aln-sort-name-bytes-vs-samtools` (`-N`);
  - `aln-sort-sam-input-vs-samtools`;
  - `aln-sort-no-hd-vs-samtools` and `aln-sort-hd-vs-samtools`;
  - `aln-sort-flags-vs-samtools` (`flags.bam`, unmapped placed and unplaced).
- [ ] **Step 3: Implement.** Port `strnum_cmp` line by line from samtools 1.21 `bam_sort.c`. Coordinate keys go in the existing `RecHead` (`chrom` = tid, `key` = pos/strand, `end` = 0). Name orders compare the QNAME and flag from the raw bytes. Then update the goldens.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci sort for BAM/SAM/CRAM (samtools coordinate, -n and -N order)`

### Task 17: `gloci sort -m` for alignments (D122)

**Files:**
- Create:
  - `gloci-ops/src/aln_merge.rs`;
  - fixtures `aln/merge-{a,b,c}.sam` (coordinate-sorted, with name-sorted twins): records tying across files at one position and strand; one `@RG` shared identically; `@PG ID:bwa` with different `CL`s; a twin whose `@RG` reuses the ID with another `SM`; `merge-c` with a different `@SQ` list.
- Modify: `gloci-ops/src/sort/merge.rs` (`merge_inputs` over an `Order` and `AlnSource`s), `commands/sort.rs` (`-m` with alignments), SKILL.md, cases.toml, differences.md (`#alignment-merge-sq`, `#alignment-merge-rg`).
- Test: `gloci/tests/merge_alignments.rs`.

**Interfaces:**
- Consumes: `Order`, `AlnSortOrder` and `with_sort_order` (Task 16); `AlnSource` (Task 3) and its CRAM implementation (Task 13).
- Produces:
  - `pub fn merge_alignment_headers(inputs: &[(&str, &AlnHeader)], order: AlnSortOrder) -> Result<AlnHeader, IoError>`;
  - `pub fn merge_alignments(srcs: Vec<Box<dyn AlnSource>>, order: AlnSortOrder, tmp_dir: &Path, out: Output, to: OutFormat, w: &WriterOptions) -> Result<u64, IoError>` (cascades above 128 inputs, as `merge_sorted` does).

- [ ] **Step 1: Write the failing tests.**
  - `merge_header_rules`:
    - an identical `@RG` is kept once;
    - a same-ID `@RG` with a different `SM` → exit 1 naming both lines;
    - for `@PG ID:bwa`, the first is kept;
    - `@CO` lines from every input;
    - different `@SQ` → exit 1 naming the first difference;
    - an `SO:queryname` input in a coordinate merge → exit 1 with the hint `pass -n/-N, or sort it: gloci sort`.
  - `merge_checks_order`: an unsorted input → exit 1 naming the file and record.
  - `merge_ties_by_input`: equal keys across three inputs come out in input order.
  - `merge_many_inputs`: 200 inputs (a cascade) equal the sort of their concatenation.
  - `merge_and_sort_cram_sam_inputs`: CRAM (`--reference`) and SAM inputs give the same output as their BAM twins, for `sort` and for `sort -m`.
  - Determinism.
- [ ] **Step 2: Add golden cases** against `samtools merge -c -p --no-PG -o - …` (compare `bam`):
  - coordinate, `-n` and `-N`;
  - plain `samtools merge` (no `-c -p`) on inputs with no colliding IDs (match);
  - the `@RG` conflict (compare `exit-status`, pinned `#alignment-merge-rg`);
  - the `@SQ` mismatch (compare `exit-status`, pinned `#alignment-merge-sq`).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS. Regenerate SKILL.md.
- [ ] **Step 5: Commit** `feat: gloci sort -m merges sorted BAM/SAM/CRAM (samtools merge)`

### Task 18: P3a benchmarks

**Files:**
- `xtask/src/gen.rs`: `gen sam` (`--records --read-len --depth --chroms --spliced-frac --shuffled --seed --ref-out FASTA`) and `gen fasta`.
- `xtask/src/bench.rs`, input kinds:
  - `bam`: `gloci view x.sam -o x.bam`, then `gloci index`;
  - `cram`: `samtools view -C -T`, then `samtools index`;
  - `fa`: `gloci index`.
- `xtask/bench/cases.toml`.
- `xtask/docker/Dockerfile.bench`: `samtools=1.21`, `mosdepth=0.3.10`.

bedtools is binding for its commands. samtools is binding for `sort`, `merge` and `cat`, which gloci now replaces (the §2 anti-goal). mosdepth and `samtools index` are yardsticks. samtools `-m` is per thread, so `-@ 8 -m 512M` matches gloci's `--max-mem 4G`.

| Case | Target |
|---|---|
| `convert-bam-50M` | `gloci <= 0.5 * bedtools` (bamtobed) |
| `genomecov-bam-50M` (`-t 1 --bga --split-deletions`) | `gloci <= 0.5 * bedtools`; yardsticks mosdepth `-t 8` and gloci `-t 8` |
| `coverage-bam-targets` | `gloci <= 0.5 * bedtools` (`-sorted`) |
| `intersect-bam-any` (`-o out.bam`) | `gloci <= 0.5 * bedtools` (`-abam -u`) |
| `coverage-cram-targets` (`-t 8`) | `gloci <= bedtools --cram-ref` |
| `getfasta-1M`, `nuc-1M`, `maskfasta-3G` | `gloci <= 0.5 * bedtools` |
| `index-bam-50M` | yardstick `samtools index` |
| `sort-bam-50M` (`-t 8 --max-mem 4G`, shuffled input) | `gloci <= samtools sort -@ 8 -m 512M` |
| `sort-name-bam-50M` (`-n`, coordinate-sorted input) | `gloci <= samtools sort -n -@ 8 -m 512M` |
| `merge-bam-8x` (`sort -m -t 8`) | `gloci <= samtools merge -@ 8` |
| `cat-bam-8x` | `gloci <= samtools cat` |

- [ ] **Step 1: Write the failing tests.**
  - `every_gloci_bench_command_parses` covers the new cases.
  - `gen_sam_is_valid`: the record count, depth and spliced fraction are as requested, `--shuffled` output is not sorted, and the reads match `--ref-out` at their positions.
- [ ] **Step 2: Implement.**
- [ ] **Step 3: Run** `cargo test -p xtask` and `cargo xtask bench --only genomecov-bam-50M --records 1000000`. Paste the medians.
- [ ] **Step 4: Commit** `bench: P3a BAM/CRAM/FASTA/sort cases (bedtools and samtools binding; mosdepth yardstick)`

### Task 19: Behaviour matrices and large-file tests

**Files:** `crates/gloci/tests/behaviour.rs`, `crates/gloci/tests/large.rs`.

Extend the matrices to `getfasta`, `nuc` and `maskfasta`; to `convert`, `view`, `intersect`, `coverage` and `genomecov` on BAM, CRAM and SAM; and to `sort`, `sort -m` and `cat` on alignments:
- determinism at -t 1/2/8;
- broken pipe;
- empty inputs: a header-only BAM/SAM, a FASTA with no sequences, an empty BED;
- failures:
  - a BAM missing its EOF block;
  - a truncated CRAM;
  - a missing `.bai`/`.crai`/`.fai`;
  - an unsorted BAM (the `gloci sort`/`--sort` hint);
  - a missing output directory;
  - a full `--tmp-dir` during `sort` (no spill files are left behind);
- a `--chrom-map` round trip (BAM header; `getfasta` intervals).

RSS ceilings (`#[ignore]`):

| Run | Ceiling |
|---|---|
| streaming `genomecov`, `convert`, and `coverage` with B, each on a 2 GB BAM | < 200 MB |
| `getfasta` of whole 300 Mb chromosomes | < 200 MB |
| `maskfasta` of 1 GB | < 200 MB |
| CRAM `convert` | < 200 MB + one reference sequence |
| `sort` and `sort -n` of a 2 GB BAM at `--max-mem 512M` | < 512 MB + 200 MB; temp disk ≤ 2× the input; spill directory removed |

- [ ] **Step 1: Write the matrices.**
- [ ] **Step 2: Run** `cargo test -p gloci --test behaviour`, then the large tests with `GLOCI_LARGE_BYTES=200000000`. Fix each real failure where it lives, with a regression test.
- [ ] **Step 3: Commit** `test: P3a determinism, broken-pipe, empty-input and failure matrices; RSS ceilings`

### Task 20: Documentation

**Files:** `docs/migration.md`, `docs/differences.md` (links), SKILL.md.

- [ ] **Step 1: Update `docs/migration.md`.** Replace each P3a "not yet available" row, checked against `--help`:
  - `bamtobed [-bed12 -split -splitD -ed -tag -cigar -color]`;
  - `getfasta`, `maskfasta`, `nuc` with their flags;
  - `genomecov -ibam` (`--split-deletions` parity; `-pc`/`-fs`/`-du` not planned);
  - `intersect -abam [-bed]` (the `-O bam | gloci convert -` pipe);
  - `coverage -b x.bam`;
  - `--cram-ref` → `--reference`;
  - `merge -i x.bam` → `gloci convert x.bam | gloci merge -`;
  - `samtools index` and `samtools faidx` in the indexing table;
  - `samtools sort [-n -N -m -@ -T -o --write-index]` → `gloci sort [-n -N --max-mem -t --tmp-dir -o --write-index]` (`-m` is per thread in samtools, a total in gloci; temp disk about the input size);
  - `samtools merge [-n -N -c -p]` → `gloci sort -m` (IDs are combined by default, and a conflicting `@RG` is an error);
  - `samtools cat` → `gloci cat`.

  Add a "BAM, CRAM and FASTA behaviour" section:
  - interval commands skip unmapped reads; `view`, `sort`, `sort -m` and `cat` keep them;
  - `--sort` works with BAM inputs;
  - the filter defaults, and the `--exclude-flags 0x704` recipe;
  - B BEDs follow the BAM's `@SQ` order (`gloci sort b.bed -g x.bam`);
  - no implicit indexes;
  - CRAM reference memory per thread.

  Link every P3a anchor from its row.
- [ ] **Step 2: Regenerate SKILL.md,** adding pitfalls for:
  - `--split` versus `--split-deletions`;
  - unmapped reads;
  - `--reference`;
  - the eight alignment commands;
  - `gloci merge` not being samtools merge (`sort -m` is);
  - `sort -n` versus `windows -n`.
- [ ] **Step 3: Run** `cargo test -p xtask` and the help-examples test. Expected: PASS.
- [ ] **Step 4: Commit** `docs: P3a migration rows, BAM/CRAM/FASTA behaviour section and SKILL.md`

---

## Execution notes for the controller

- **Order.** Task 1 merges first, then the three lanes run in separate worktrees, each with its own `target/`. Use `CARGO_BUILD_JOBS=1`, scoped test commands, and one `cargo test --workspace` per task before merging. Tasks 19 and 20 run last.
- **Docker.** Rebuild `gloci-diff:p0` after Task 1, for samtools. Fixture scripts (D119) run in that image, and each task commits the fixtures it creates. The bench image changes in Task 18.
- **Spec amendments** (controller):
  - §4.5 rule 2 (D103);
  - §5.1 SAM, PSL and RepeatMasker detection (D100, D114);
  - §6.1: `--min-mapq`, `--require-flags`, `--exclude-flags`, `--reference`, `--split-deletions` (D102, D105, D109);
  - §7.3 `genomecov` (D117);
  - §7.4 `view` (D116); `sort` for alignments with `-n`/`-N` (an owner exception to one flag, one meaning, D121); `sort -m` as samtools merge (D122); a new `cat` entry (D123); `--sort` with alignments (D124);
  - §11: a samtools table (`sort`, `merge`, `cat`, `index`, `faidx`);
  - §7.5 flag names (D111–D113);
  - §8.3 the CRAM pool (D110).
