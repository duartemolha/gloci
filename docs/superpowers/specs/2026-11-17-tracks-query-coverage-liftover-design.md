# bigWig/bigBed, `vcf query`, coverage summaries and liftover: design

Date: 2026-11-17. Status: approved by the owner, 2026-11-17 (bigtools on its MIT licence).

## Intent

gloci replaces the everyday parts of bedtools, bedops, vcftools, tabix, bgzip and
samtools. Pipelines still call other tools beside it for four jobs that come up
constantly in a genomics team's work on multi-GB files:

1. coverage and signal tracks as **bigWig**, and annotation tracks as **bigBed**
   (today UCSC's `bedGraphToBigWig`, `bigWigToBedGraph`, `bigWigAverageOverBed`,
   `bedToBigBed`, `bigBedToBed`);
2. getting fields out of a VCF with a **format string** (`bcftools query -f`);
3. **coverage QC summaries**: per-chromosome tables (`samtools coverage`) and
   per-region means and depth thresholds (mosdepth `--by`, `--thresholds`);
4. moving intervals between assemblies with **chain files** (UCSC `liftOver`,
   CrossMap).

The owner asked for all four (2026-11-17). Each lands as an experimental command or
option (design.md stability tiers), checked against the tool it replaces where a
reference tool may be used (below).

## Owner constraints (2026-11-17)

- **gloci must stay free for commercial use, with no dependence on UCSC's licence.**
  UCSC's programs and source are free only for academic, non-profit and personal
  use. So:
  - no UCSC code, binary or library in gloci, its tests, its comparison image or its
    benchmark image. The benchmark image's `ucsc-bigwigtobedgraph` (used to prepare
    the map case's input) is replaced by gloci's own bigWig reader once it exists;
  - formats and behaviour are implemented from their published descriptions (the
    bigWig/bigBed format paper and spec; the chain format description), never from
    UCSC source;
  - every new runtime dependency is under a permissive licence (MIT, Apache-2.0,
    BSD, zlib or similar), checked when it is added;
  - reference tools used only in tests may have other open licences (CrossMap) as
    long as they never ship with gloci;
  - gloci's code, comments, docs, tests and commit messages never name UCSC's code
    base; they cite the published format descriptions and the behaviour of tools only
    (owner, 2026-11-17).
- If the bigWig format turned out to be hard to implement without UCSC, the owner
  would rather drop it. It is not hard: an independent MIT library exists (below).
- Standing rules apply: streaming with flat memory on multi-GB inputs; output that
  never changes with `--threads`; one static binary; no silent wrong answers; golden
  files committed for parity; intentional differences documented in
  `docs/differences.md`.

## Approaches considered

| Feature | Picked | Rejected and why |
|---|---|---|
| bigWig/bigBed | The `bigtools` crate (0.5.8, MIT; the owner accepted it on its licence, 2026-11-17; deps tokio, futures, serde, bincode, crossbeam, libdeflater: all MIT or Apache-2.0), with its `read` and `write` features only (no `remote`, no `cli`) | Writing our own reader and writer: several thousand lines for a fiddly format (B+ tree, R-tree index, zoom levels), no gain in behaviour. Skipping the format: the owner would only skip it if it needed UCSC, which it does not |
| `vcf query` | A new `gloci vcf query` subcommand on the shared selection stage (D157) | Extending `vcf extract` with a format string: `extract` writes fixed tables and matrices; two output models in one command would confuse both |
| Coverage QC | New options on `gloci depth` (`--summary`; `--by`, `--thresholds`, `--out-prefix`) | A new command: these are summaries of the depth `depth` already computes, and `coverage` is taken by bedtools coverage (stable core) |
| Liftover | Our own chain mapper with whole-record rules (UCSC `liftOver`'s documented behaviour: a record maps whole or goes to the unmapped file) | CrossMap-style splitting of a record across chain gaps: less familiar, and it silently changes record counts. Wrapping CrossMap: a Python dependency, GPL |

## Feature 1: bigWig and bigBed

### Reading

- Two new formats, `bigwig` and `bigbed`, detected from the magic number in the first
  4 bytes (bigWig `0x888FFC26`, bigBed `0x8789F2EB`, little-endian; also accept the
  byte-swapped forms). `--format bigwig|bigbed` overrides.
- Read only from files: the format needs random access. Standard input, a pipe or a
  `/dev/fd` path is a usage error that says so.
- **Records.** bigWig records are bedGraph-like (`chrom start end value`), as WIG is
  read (D63), and work wherever an interval input is accepted, with value column 4.
  The value is a 32-bit float, printed as the shortest decimal that reads back as the
  same float. bigBed records are BED records: `chrom start end` plus the stored rest
  of the line, field for field.
- **Order (builder's question, settled as my call).** The reader fetches each
  chromosome through the file's index, in the order the reading command needs (its
  `--chrom-order` or `-g`, lexicographic by default), whatever order the file stores
  them in: files written through pyBigWig (deepTools `bamCoverage`, for one) store the
  BAM header's order, and reading them in stored order would make `map -b cov.bw`
  fail its order check. So the reader is told the command's order when it opens the
  file; within a chromosome, starts ascend as stored.
- **Regions.** `gloci view x.bw chr1:1-1000` (and `--regions-file`) read through the
  file's own index, no scan. `--chrom-map` applies (D62), with region names
  translated back to the file's own names.
- **`bigWigAverageOverBed`**: `gloci map -a regions.bed -b signal.bw -c 4:wmean -c
  0:bases_uniq` gives its `mean` and `covered` columns; `map` gains a `wsum` operation
  (value × overlapping bases), which gives `sum`. `mean0` is `sum / length`. The
  migration guide says so.
- `gloci convert x.bw` writes bedGraph; `gloci convert x.bb` writes BED.

### Writing

- **bigWig** when `-o` ends in `.bw` or `.bigwig`, or with `-O bigwig` (D17), from
  output that is bedGraph-shaped: `genomecov --bedgraph` and `--bga`, and `convert`
  of bedGraph, WIG or bigWig. Any other command or mode given a bigWig output is a
  usage error naming the ones that take it.
- **bigBed** when `-o` ends in `.bb` or `.bigbed`, or with `-O bigbed`, from `convert`
  of BED (3 to 12 standard columns plus any extra). The autoSql describes the standard
  columns by their BED names; extra columns get `field1`, `field2`, … as strings,
  unless `--autosql FILE` gives one.
- **Chromosome sizes** come from `-g` (a built-in genome or a file), or from the
  input's own header where it has one (an alignment header for `genomecov`, a bigWig
  being converted). Without them the run is a usage error naming `-g`.
- Written to a file only: the format is written with seeks. Standard output is a
  usage error.
- **Input rules.** Records must be sorted in the order the format requires (check
  what `bigtools` needs; gloci's lexicographic default is expected to satisfy it);
  an input in another order is an error naming the fix (`gloci sort`). bigWig records
  must not overlap and must have a numeric value; an overlap or a non-number is an
  error naming the line. A value that is not exactly representable as a 32-bit float
  is stored rounded: a documented difference (`#bigwig-float32`), not an error.
- **Determinism.** The bytes written are the same for every `--threads`. `bigtools`
  encodes on a tokio runtime; if a multi-threaded runtime makes the output depend on
  timing, use a current-thread runtime and say so in D171.
- Zoom levels and compression: `bigtools`' defaults (deflate through libdeflater),
  recorded in D171.
- Peak memory stays flat in the input size apart from what the format itself needs
  (the chromosome list and per-section index entries); measure it on the GENCODE
  exons and the CTCF bedGraph and record it.

### Checking

No UCSC tool is used. bigWig and bigBed written by gloci are read back by pyBigWig
(MIT, libBigWig; an independent C implementation) in the comparison image, and
compared value by value with the bedGraph or BED they came from; gloci reads files
written by pyBigWig, and the ENCODE bigWig the benchmark already downloads. Round
trips (bedGraph → bigWig → bedGraph) are byte-for-byte golden cases.

## Feature 2: `gloci vcf query`

- `gloci vcf query FILE -f FORMAT [selection options] [-H] [-o OUT]`. Every selection
  option of `vcf filter`/`vcf extract` applies (samples, regions, `--include`/
  `--exclude`, MAF, …), through the same stage (D157), so `-t` parallelises it and
  the output never changes with `-t`.
- **Format language: bcftools 1.21 `query`'s**, for the common subset:
  - fixed fields `%CHROM %POS %POS0 %END %END0 %ID %REF %ALT %FIRST_ALT %QUAL
    %FILTER %TYPE`, `%INFO` (the whole column), `%INFO/TAG` and `%TAG` (an INFO tag
    outside brackets), subscripts `%TAG{N}`, and `%LINE`;
  - a per-sample block `[ … ]` with `%SAMPLE %GT %TGT` and FORMAT tags (`%DP %AD
    %GQ …`), subscripts included;
  - escapes `\t`, `\n`, `\\`, and literal text;
  - `-H` writes bcftools' header line (`# [1]CHROM [2]POS …`).
- Missing values print as bcftools prints them (`.`), vectors joined with `,`, flags
  as bcftools prints them. Check every case against bcftools 1.21.
- **Sample-level filters (builder's question, settled as my call).** When
  `--include`/`--exclude` tests FORMAT or GT values, bcftools `query` prints a site if
  any sample passes and, inside the `[ … ]` block, only the samples that pass
  (`-i 'GT="alt"' -f '[%SAMPLE\n]'` lists the carriers). `vcf query` does the same.
  The shared stage keeps its site rule (D71) for every other command; the expression
  engine gains a per-sample evaluation that returns a pass mask, carried with each
  batch through the D157 workers, used only by `vcf query`. Verify the exact rule
  (sites with no passing sample, `-e`, mixed site and sample terms) against bcftools
  1.21 and record it in D172. Keeping D71's rule here would change ported scripts'
  answers without a word, which the anti-goals forbid.
- **Anything outside the subset is a usage error naming the construct** (exit 2),
  never a blank or a guess. The list of supported tags is in `--help` and in D172.
- Output: text, to standard output or `-o` (`.gz` compresses, D17).

### Checking

Golden cases generated with bcftools 1.21 `query` on the existing VCF/BCF fixtures
(haploid, missing, multi-allelic, phased, Number=A/R/G, flags, symbolic ALTs), and
live cases in `xtask/diff/cases.toml`.

## Feature 3: coverage summaries on `gloci depth`

`depth` already streams per-position depth from BAM, CRAM and SAM with samtools'
read filters. Two new output modes; each replaces the per-position rows.

### `--summary`: the `samtools coverage` table

- One row per reference, in `@SQ` order: `#rname startpos endpos numreads covbases
  coverage meandepth meanbaseq meanmapq`, with samtools coverage 1.21's definitions,
  default read filters and number formatting. `--region` limits it to one region, as
  samtools `-r` does.
- With several inputs, do what samtools coverage 1.21 does with several inputs
  (check it), not one table per input.
- samtools coverage's histogram mode (`-m`) is not provided.

### `--by` and `--thresholds`: mosdepth's region outputs

- `gloci depth reads.bam --by regions.bed|SIZE [--thresholds 1,10,20] --out-prefix P`
  writes mosdepth 0.3.10's files, with mosdepth's names and formats, so tools that
  parse them (MultiQC among them) read gloci's unchanged:
  `P.mosdepth.global.dist.txt`, `P.mosdepth.region.dist.txt`,
  `P.mosdepth.summary.txt`, `P.regions.bed.gz` and, with `--thresholds`,
  `P.thresholds.bed.gz` (each `.bed.gz` BGZF with its `.csi`, as mosdepth writes).
- `--by SIZE` makes windows of SIZE bases over every reference; `--by FILE` uses the
  BED's regions (a fourth column is kept as the region name). `--thresholds` needs
  `--by`; `--by` needs `--out-prefix`. No per-base file: that is `gloci depth` itself.
- Defaults follow mosdepth 0.3.10 where they decide the numbers: its read filters
  (flag exclusion 1796, MAPQ 0) and its handling of overlapping mates. Where gloci's
  existing `depth` defaults differ (samtools'), this mode uses mosdepth's and the help
  says so.
- Memory stays flat: per-reference depth arrays as mosdepth keeps them are not
  required; a streaming accumulator per region and per threshold is.
- `-t` never changes the files.

### Checking

Golden cases against samtools coverage 1.21 and mosdepth 0.3.10 on the existing BAM
and CRAM fixtures; mosdepth moves from the benchmark image into the comparison image
too. The README benchmark gains a `depth --by` row beside mosdepth if it is cheap to
add.

## Feature 4: `gloci liftover`

- `gloci liftover FILE --chain CHAIN [-o OUT] [--unmapped FILE] [--min-match FRAC]
  [--multiple]`. `CHAIN` is a chain file, plain or gzip/BGZF; the user supplies it
  (UCSC and Ensembl both publish them as free data). Nothing is downloaded and no
  chain is built in.
- Inputs in this release: BED (3 to 12 columns), bedGraph and GFF/GTF. Every other
  column is carried unchanged. **VCF liftover is out of scope** for this release
  (owner default, 2026-11-17): it needs REF checks against the target FASTA, REF/ALT
  swaps and INFO decisions, and comes later as its own design.
- **Mapping rules: whole records.** For each chain whose source span overlaps the
  record, count the record's bases that fall in the chain's aligned blocks. A chain
  whose share of the record's bases reaches `--min-match` (default 0.95) is a
  candidate; the record's target interval runs from the target position of its first
  aligned base to that of its last, plus one. A record with no candidate goes to the
  unmapped file with the reason as a comment line before it: `#Deleted in new` (no
  chain overlaps it), `#Partially deleted in new` (none reaches `--min-match`),
  `#Duplicated in new` (more than one candidate, without `--multiple`; with it, every
  candidate is written). These are `liftOver`'s documented reasons, so scripts that
  read its unmapped files read gloci's.
- **Strand.** A chain onto the target's reverse strand maps the interval onto the
  reverse coordinates and flips the strand column (BED column 6, GFF column 7; `.`
  stays `.`).
- **BED12.** Each block is mapped; the record maps when its span maps by the rules
  above and every block lands in the same chain; `thickStart`/`thickEnd` are mapped
  the same way. The exact block rules are the builder's to settle and record in D174,
  with fixtures.
- **Zero-length records** (insertion points) map when their position lies inside an
  aligned block.
- **Output order** is input order, with target chromosome names from the chain. The
  output is often unsorted; the help says to pipe it through `gloci sort`. The
  command never claims sorted output.
- Without `--unmapped` the unmapped records are counted in a warning (how many, and
  how many per reason), never dropped silently.
- Memory: the chain's blocks, indexed by source chromosome (hg19 to hg38 is a few MB);
  records stream.

### Checking

No UCSC `liftOver`. Hand-built chain fixtures, each case's expected output written by
hand: gaps, `--min-match` on both sides of the threshold, duplicates with and without
`--multiple`, reverse-strand chains, BED12, points. Real-data cases with Ensembl's
GRCh37-to-GRCh38 chain (free data) against CrossMap, which runs only in the
comparison image: the two must agree on every record that falls whole inside one
aligned block; where CrossMap splits a record and gloci does not, the difference is
documented (`#liftover-whole-records`).

## Common to all four

- Each feature gets a decision: D171 (bigWig/bigBed), D172 (`vcf query`), D173
  (`depth` summaries), D174 (liftover), in `docs/decisions.md`; design.md's format,
  command and non-goal sections are updated; the migration guide gains rows for every
  replaced command; differences get anchors and golden cases; CHANGELOG entries.
- Help follows the enforced rules (`crates/gloci/tests/cli_basics.rs`): summaries of
  at most 64 characters naming the tools replaced, 2 to 4 examples, the experimental
  marker. SKILL.md, the Python stub and the book reference are regenerated.
- Python: new stream commands are placed by the existing rules (D87): `vcf query`
  and `liftover` as steps, the `depth` prefix outputs like `vcf_stats`'. bigWig and
  bigBed are files to `save()`, not records to iterate, unless that is cheap.
- The comparison image (`xtask/docker/Dockerfile`) gains pyBigWig, mosdepth and
  CrossMap from bioconda/conda-forge; no `ucsc-*` package anywhere. The benchmark
  image drops `ucsc-bigwigtobedgraph` once gloci converts bigWig itself.
- Builder's calls (not the owner's): internal module layout, the exact text of
  messages, how bigWig reading plugs into the existing reader traits, zoom-level
  parameters, and the order of the work within each feature.

## Order of the work

Three lanes in parallel (worktrees, as before): **A** bigWig/bigBed; **B** `vcf
query`; **C** `depth` summaries, then liftover. Each task reviewed, each lane merged
into local main when green; nothing pushed or published.
