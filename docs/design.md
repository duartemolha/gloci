# gloci — Design Spec

This is the design specification gloci was built from. The source cites it by section (`spec §4.3`, or a bare
`§4.3`) and cites the numbered design decisions by number (`D37`); the decisions are in
[`decisions.md`](decisions.md), grouped by the phase that made them. The section numbers and headings below
are cited and stay as written.

Two markers in the text say where a choice came from:

- `[owner, DATE]` (also written `[you, DATE]` or `[you]`, in the spec's early drafts, where "you" is the
  project owner) is a decision the project owner made, on that date when one is given.
- `[my call]` is a default the spec's author chose instead of asking.

Text with no marker follows from those decisions. Later amendments are marked where they apply.

**Status:** Draft for review · **Date:** 2026-11-06

A single, fast, correct Rust command-line tool that replaces **bedtools**, **bedops** and **vcftools** for genomic interval and variant work on files of many gigabytes.

> Marking convention: **[you]** = decided by the project owner in conversation. **[my call]** = chosen by the spec author as a sensible default; change freely. Unmarked text is a direct consequence of those decisions.

---

## 1. Why

- bedtools appears to be no longer actively maintained. **[you]**
- The team wants to remove bedtools, bedops and vcftools as dependencies and run one maintained tool instead. **[you]**
- The current tools give quiet wrong answers at the edges. The motivating example: `bedtools subtract` with a VCF builds each variant's interval from position and REF length and ignores the alleles, so a subtraction can remove the wrong number of bases. **[you]**
- Real inputs are large: several GB per file, sometimes bigger. **[you]**
- The old tools are inconsistent with each other: different flags for the same idea, different coordinate assumptions, and one statistic per vcftools run.

## 2. Goals, non-goals, anti-goals

### Goals

1. **Replace the three tools** for the work a genomics team actually does: interval set operations, mapping and aggregation, coverage, sorting and indexing, format conversion, VCF filtering and VCF statistics. **[you]**
2. **Handle multi-GB inputs** with streaming and flat memory use, multiple cores, and direct reading and writing of bgzipped and indexed files. **[you, confirmed playback]**
3. **Be correct before being fast.** Coordinate handling is precisely defined (Section 4), variant intervals come from alleles, and nothing is silently skipped or assumed. **[you, confirmed playback]**
4. **Provide one coherent interface.** One binary with subcommands, the same flag for the same idea everywhere, and output that is safe to pipe. **[you, confirmed playback]**
5. **Ship as one static binary** with no runtime dependencies such as Perl, Python or htslib, so it is easy to install on HPC clusters. **[my call]**
6. **Be usable as a Python library** (P2), filling the role pybedtools fills for bedtools. It calls the same Rust engine in-process, with no subprocess or temp files (Section 8.5). **[you]**

### Non-goals

- **Command or flag compatibility with the old tools.** Pipelines will be rewritten. **[you]** A migration guide (Section 11) maps old commands to new ones.
- Byte-identical output with the old tools. Differences are allowed but must be documented (Section 9.3).
- Variant calling, annotation (VEP-like), imputation and phasing. Left-alignment and normalization against a reference were a non-goal too, until the owner's VCF additions to P3b (Section 10): `vcf norm -f` (bcftools `norm -f`) and `vcf split` (bcftools `norm -m-`) are built. **[owner, 2026-11-11]**
- Becoming a general replacement for bcftools or samtools. Merging samples across VCFs (`vcf-merge`) was not planned for v1; it is built as `vcf merge` (bcftools `merge`), one of the owner's VCF additions to P3b (Section 10). **[owner, 2026-11-11]**
- The bedops **Starch** format, and legacy genotype export formats (IMPUTE, BEAGLE, LDhat, PLINK ped/tped). **[my call]**
- bedtools' BAM-specific utilities were a non-goal **[my call]**, reversed on the owner's request of 2026-11-15: `igv`, `links`, `bamtofastq`, `bedtobam`, `bedpetobam` and `tag` are built (D165 is their scope; D169, D170; `bamtofastq` as `gloci convert --fastq`, D166; `bedtobam` and `bedpetobam` as `gloci convert --bam`, D167; `tag`, D168). **[owner, 2026-11-15]** `multicov` is `gloci bedcov X.bed A.bam B.bam --read-count`.
- VCF liftover (CrossMap `vcf`) is not in this release; it needs REF checks against the target FASTA and comes later as its own design **[owner default, 2026-11-17]**.
- Native Windows. WSL works. **[my call]**
- A stable public Rust library API in v1. The crates are internal; the stable, versioned API is the Python one from P2 (Section 8.5). **[my call, revised: Python bindings moved to goals]**

### Anti-goals (failure even if it "works")

- **Silent wrong answers.** That includes accepting unsorted input without saying so, skipping malformed lines without saying so, guessing at coordinates, and building VCF intervals that ignore alleles when the user didn't ask for that.
- **Memory that grows with file size** in streaming operations.
- **Output that changes with `--threads`.**
- **Being slower than the tool it replaces** on any operation it claims to replace, as measured on the benchmark suite.
- **A grab-bag of inconsistent flags**, which is the experience being replaced.
- **Partial output files left behind** after an error.

## 3. Shape of the tool

```
gloci <command> [options] [inputs]
gloci vcf <command> [options] <input>
```

```sh
gloci sort big.bed -g hg38.fa.fai --max-mem 8G -o big.sorted.bed.gz --write-index
gloci intersect -a calls.vcf.gz -b targets.bed -o on_target.vcf.gz --write-index
gloci subtract -a targets.bed -b calls.vcf.gz                  # allele-aware by default
gloci map -a genes.bed -b signal.bedgraph -c 4:mean -c 4:max
gloci merge peaks_*.bed.gz --max-gap 50 -c 5:max
gloci vcf filter cohort.bcf --pass-only --min-qual 30 --min-maf 0.01 --max-missing 0.1 -o filtered.bcf
gloci vcf stats cohort.vcf.gz --freq --missing-site --missing-indv --het --tstv-summary --out-prefix qc/cohort
```

**Name:** `gloci`, for "genomic loci", pronounced **GLOH-sigh** (/ˈɡloʊsaɪ/), with the ending said the way English speakers say "loci". **[you]** It is used for the binary, the crate, the repository and the bioconda package. On 2026-11-06 the name was free on crates.io, bioconda, conda-forge, PyPI, CRAN and Debian. The nearest existing name is `cloci`, an unrelated Python gene-cluster tool.

---

## 4. Core semantics

This section is the contract. Every rule here has a table-driven unit test, and a rule never changes without this section changing first.

### 4.1 Coordinates

- Internally every feature is a **0-based, half-open** interval `[start, end)` on a named chromosome.
- Positions are **64-bit unsigned**. Some plant and lungfish chromosomes are longer than 2³² bp. **[my call]**
- Conversion on input:

| Format | On disk | Internal |
|---|---|---|
| BED, bedGraph | 0-based half-open | unchanged |
| VCF, BCF | 1-based `POS` | see 4.2 |
| GFF3, GTF | 1-based closed `[s, e]` | `[s-1, e)` |
| Region strings `chr:s-e` | 1-based closed (samtools convention) | `[s-1, e)` |

- A **point** is a zero-length interval `[p, p)`. It sits between 0-based bases `p-1` and `p`. Points are valid input: BED allows `start == end`, and insertions are points.
- `start > end` is an error. Negative or non-integer coordinates are errors.

### 4.2 Variant footprints (VCF/BCF to interval)

Selected with the global `--vcf-span <allele|ref>`. **The default is `allele`.** **[my call, driven by the motivating bug]**

- **`ref`:** `[POS-1, POS-1+len(REF))`. This is the bases REF covers, and it matches what bedtools does.
- **`allele`:** the reference bases each ALT actually changes. For each ALT, trim the **shared prefix first, then the shared suffix**, from REF and ALT. The remaining REF bases are the footprint. If no REF bases remain, the footprint is a **point** at the trim position.
- A record's footprint in `allele` mode is the smallest interval covering the footprints of all its ALTs. Use `gloci convert --split-alleles` to get one interval per ALT.

| POS | REF | ALT | Kind | `ref` | `allele` |
|---|---|---|---|---|---|
| 100 | A | G | SNV | [99,100) | [99,100) |
| 100 | AC | GT | MNV | [99,101) | [99,101) |
| 100 | ACGT | A | deletion of CGT | [99,103) | [100,103) |
| 100 | A | ACGT | insertion after base 100 | [99,100) | point [100,100) |
| 100 | ATG | AC | complex, TG→C | [99,102) | [100,102) |
| 100 | ACGT | A,AC | multi-allelic deletions | [99,103) | [100,103) |
| 100 | N | `<DEL>` with END=500 | symbolic deletion | [99,500) | [100,500) |
| 100 | N | `<INS>` | symbolic insertion | [99,100) | point [100,100) |
| 100 | N | `<DUP>`/`<INV>`/`<CNV>` with END=E | symbolic SV | [99,E) | [100,E) |
| 100 | G | `G]17:198982]` | breakend | [99,100) | [99,100) |
| 100 | A | `<*>` or `<NON_REF>` only, END=E | gVCF reference block | [99,E) | [99,E) |
| 100 | A | `.` | no ALT | [99,100) | [99,100) |

Further rules:

- **Symbolic SVs:** use `INFO/END` when present. Otherwise use `POS + |SVLEN|` as the 1-based last affected base. If neither is present, it is an error, or a warning with `--skip-invalid`.
- **`*`, `<*>` and `<NON_REF>` alternative alleles** contribute no footprint when other ALTs exist.
- **Breakend mates** are ignored in v1.
- **gloci does not left-align footprints.** Inside repeats the footprint depends on how the caller wrote the variant. The docs recommend normalized input (`gloci vcf norm -f REF`, or `bcftools norm`).

### 4.3 Overlap

- Two non-empty intervals overlap when `a.start < b.end && b.start < a.end`. They must share at least one base. **Bookended intervals do not overlap.**
- **A point `p` vs a non-empty interval `[s,e)`:** they overlap when `s < p < e`, meaning the point is strictly inside with interval bases on both sides. `--inclusive-points` relaxes this to `s <= p <= e`. **[my call]**
- **Point vs point:** they overlap when the points are equal.
- Worked examples (BED coordinates). A zero-length BED interval is a point, so an insertion site inside a feature is found. bedtools has historically missed this case.

| A | B | Overlap? | Why |
|---|---|---|---|
| `chr1 10 10` | `chr1 5 15` | **yes** | point 10 is strictly inside `[5,15)` |
| `chr1 5 5` | `chr1 5 15` | no (yes with `--inclusive-points`) | point sits on B's left edge, before base 5 |
| `chr1 15 15` | `chr1 5 15` | no (yes with `--inclusive-points`) | point sits on B's right edge, after base 14 |
| `chr1 10 10` | `chr1 10 10` | yes | equal points |
| `chr1 10 10` | `chr1 11 11` | no | different points |
| `chr1 5 10` | `chr1 10 15` | no | bookended, no shared base |
| `chr1 5 11` | `chr1 10 15` | yes | share base 10 |
- **Fraction thresholds:**
  - `-f` is the minimum fraction of A covered, and `-F` is the minimum fraction of B covered. Both accept a decimal (`0.5`) or a percentage (`50%`), as bedops `-e 50%` does. **[you, 2026-11-07]**
  - `-r` requires both, and `-e` requires either.
  - `--min-bases N` sets an absolute minimum overlap.
  - Every threshold (`-f`, `-F`, `-r`, `-e`, `--min-bases`) counts as met when either side of the overlap is a point. A point covers no bases, so no threshold can be measured, and an insertion inside a feature is never dropped by a threshold. **[my call, revised 2026-11-08: the earlier rule applied only to the point's own side, so `-F` dropped insertions while `--min-bases` kept them]**
- **Strand:**
  - `-s` keeps only overlaps on the same strand, and `-S` only overlaps on opposite strands.
  - Features with unknown strand (`.`, or any VCF record) never match under `-s` or `-S`, and gloci warns once.
- **BED12 blocks:** `--split` uses blocks (exons) instead of the full span.
- **Windows:**
  - `-w N` pads B by N on both sides before testing.
  - `--left` and `--right` set the padding asymmetrically.
  - `--strand-window` interprets left and right relative to A's strand.

### 4.4 Distance (closest, spacing)

- Overlapping features have distance **0**.
- Non-overlapping features have distance equal to the gap in bases **+ 1**, so bookended features have distance 1. This matches bedtools and keeps "overlapping" distinct from "adjacent".
- Signed distances (`--signed ref|a|b`) are negative for upstream features, with the reference frame chosen by the flag.

### 4.5 Chromosome order and sortedness

- Streaming operations need inputs sorted by **(chromosome, start)**. The sort order of `end` doesn't matter for correctness.
- **Chromosome order is chosen in this order** **[my call]**:
  1. `--genome NAME|FILE` order, if given. NAME is a built-in genome such as `hg38` (Section 5.3). FILE can be `.fai`, a 2-column `.genome`, `.dict`, or a VCF/BCF/BAM/SAM/CRAM header.
  2. Otherwise, the contig order declared in the header of an input, if any: a VCF/BCF header's `##contig` lines, or an alignment header's `@SQ` lines (BAM, SAM, CRAM), taken from the first such header in input order, and closed (a name outside it is an error). A VCF/BCF header that lists no contigs declares no order, so the input falls through to rule 3 like any other **[owner, 2026-11-15]**. **[D103, as amended by D124]** The order of an alignment header is a declared order like a VCF header's. An alignment out of that order is an error whose hint is `gloci sort x.bam -o x.sorted.bam` (naming `SO:queryname` when the header declares it) or `--sort`; a B BED must follow the same order (`gloci sort b.bed -g x.bam`), or be sorted by `--sort` too. `--unsorted` holds B alignments in memory. `--chrom-map` renames the binary reference list and the `@SQ SN:` values, so BAM and CRAM records, which name references by number, carry the new names in any output, SAM included. SAM text is written as it is read, so it keeps its `RNAME`/`RNEXT`: `--chrom-map` is refused whenever SAM input lines would be written as records (`view`, `sort` and `sort -m` of a SAM input, whatever the output format, and `intersect --mode any|none` of a SAM A; convert to BAM first). The commands that only compare SAM reads as intervals (`intersect --mode pairs|count`, `coverage`, `genomecov`, `convert`) rename their chromosomes.
  3. Otherwise, **lexicographic order** (byte order of the names) **[owner, 2026-11-15, reversing D32]**, the order of `sort -k1,1` (with `LC_ALL=C`), bedops `sort-bed` and `bedtools sort`, so files from the tools gloci replaces just work. `--chrom-order karyotype` selects **karyotype order** instead, which matches how people read chromosomes and how reference genomes (and therefore most VCFs) are laid out. Karyotype order compares names with any leading `chr` prefix ignored:
     1. names that are a number, optionally followed by letters (`1`, `2`, `22`, `2L`, `2A`), by number and then suffix;
     2. then `X`, `Y`, and `M`/`MT`, in that order;
     3. then everything else (`chr1_KI270706v1_random`, `chrUn_…`, `chrEBV`, `HLA-…`) in natural order, with digit runs compared as numbers;
     4. remaining ties (`1` vs `chr1`, `chr01` vs `chr1`) are broken by raw bytes, so the order is total.

     This gives `chr1 … chr22, chrX, chrY, chrM` for human, the GRCh38 layout. Names that don't follow these patterns, such as yeast's Roman numerals, need `-g` in either order to get the reference order.
- **Validation:** sortedness is **checked on the fly** in every streaming command.
  - At the first violation, the command stops with exit code 1.
  - The error names the file, line number and both offending records, and gives a hint such as: `run: gloci sort FILE -g <genome>`.
  - This check is cheap and is never turned off.
- **`--unsorted`** and **`--sort`** are two different choices for unsorted input, each with one meaning everywhere. **[you, 2026-11-08]**
  - `--unsorted` means "work with my input in the order it is". It exists in the relation commands, where it loads the B side into an interval index and streams A in any order. Results follow A's order, and memory is proportional to the loaded file; the help text says so.
  - `--sort` means "sort my inputs first". It exists in `merge`, `complement`, `genomecov` and the relation commands. The inputs go through the external sorter, so memory stays within `--max-mem`, and results come out in sorted order. **[D124]** With alignment inputs (BAM, SAM, CRAM) the relation commands and `genomecov` sort them with the interval sorter into a temporary file and then stream: a SAM input spills a plain SAM temp, which preserves its layout, and BAM and CRAM spill a level-1 BAM. Memory stays within `--max-mem`. The alignments are keyed by `gloci sort`'s coordinate order (D121: reference, position, strand, then input order; the table is the header's) and the temporary header says `SO:coordinate`, so `intersect --sort` writes what `gloci sort` followed by `intersect` writes, `@HD` included. A `-g` table that orders the references unlike the header is followed instead, and the header says `SO:unsorted` (P3c final fix wave I5: keeping a header's `SO:coordinate` there was false).
- **Indexed inputs** (bgzip with `.tbi`/`.csi`) can be read chromosome by chromosome through the index. Chromosome order differences between inputs therefore don't matter in that mode.
- **Chromosome naming:**
  - Names must match exactly.
  - If two inputs share no chromosome names, gloci warns. It detects the common `chr1` vs `1` case and suggests `--chrom-map FILE`, a 2-column rename table applied on the fly (Phase 2).

### 4.6 Records and passthrough

- Input lines are carried through **unchanged** wherever the command outputs input records. Extra BED columns and the full VCF line are kept, and numbers are not reformatted.
- Commands that output input records unchanged write **A's format**. Examples are `intersect --mode any|none`, VCF `subtract`, and `vcf filter`. A VCF A gives VCF out with its header.
- Commands that **append columns**, such as `--mode pairs`, `count`, `map`, `closest` and `coverage`, write **TSV**.
  - BED-in gives BED-like TSV out.
  - VCF-in gives VCF columns plus appended columns, with no VCF meta-header, preceded by a single `#`-prefixed column-name line.
- VCF output gets a `##glociCommand=<full command line>; Version=<x.y.z>` header line for provenance. **[my call]**
- **Subtracting from VCF records:** a variant can't be partly removed. When A is VCF, `subtract` removes whole records, the same as `intersect --mode none` with the same criteria.

---

## 5. Formats

### 5.1 Input

| Format | Notes | Phase |
|---|---|---|
| BED3–BED12, BED-like (BED3 + any columns) | tab-separated; skips `#`, `track`, `browser` lines; CRLF tolerated | 1 |
| bedGraph | 4th column numeric | 1 |
| VCF 4.1–4.5 | plain or bgzip; header required | 1 |
| BCF 2.1/2.2 | | 1 |
| GFF3, GTF | for `convert` and as interval input | 1 |
| WIG (fixed/variableStep) | convert to bedGraph | 2 |
| BEDPE | for pair ops | 3 |
| BAM, CRAM, SAM | as intervals (alignment spans or blocks) for `intersect`, `coverage`, `genomecov` and `convert`; every record for `view`, `index`, `sort` and `cat` (`cat`: BAM only) | 3 |
| FASTA + `.fai` | for `getfasta`, `maskfasta`, `nuc` | 3 |
| PSL, RepeatMasker `.out` | convert only | 3 |
| Chain files | for liftover --chain only: plain, gzip or BGZF; whitespace and ids as published (D174) | after P3 |
| bigWig, bigBed (experimental) | read from a file through its own index, as bedGraph and BED records, in the command's chromosome order (D171) | — |

- **Format detection** uses content and magic bytes: the BGZF/gzip header, `BCF\2`, the bigWig magic number (`26 FC 8F 88`, or `88 8F FC 26` in the other byte order) and the bigBed one (`EB F2 89 87` or `87 89 F2 EB`), `##fileformat=VCF`, `##gff-version`, and BED column shape. File extensions are not used. `--format` / `--a-format` / `--b-format` override detection. **[D100, D114]** The alignment and sequence formats are detected by content too: BAM is a BGZF stream that decompresses to `BAM\1`; CRAM is `CRAM` with major version 2 or 3 (4 is an error); SAM is an `@HD`/`@SQ`/`@RG`/`@PG`/`@CO` line followed by a TAB, or at least 11 tab fields whose fields 2, 4, 5, 8 and 9 are integers; FASTA starts with `>`, after any blank lines (as samtools skips them); PSL and RepeatMasker are recognised by their header or their column shape. `--format` adds `bam`, `cram`, `sam`, `psl` and `rmsk`. Only `intersect`, `coverage`, `genomecov`, `convert`, `view`, `index`, `sort`, `cat` and, **[D125]**, `idxstats`, `flagstat`, `check`, `reheader`, `head`, `samples`, `depth` and `bedcov` read alignments; any other command exits 2 naming them and suggesting `gloci convert x.bam -o x.bed` (and, in `merge`, `gloci sort -m`). PSL and RepeatMasker are read by `convert` only, and FASTA is never an interval input.
- **Compression:** gzip and BGZF are read transparently. Plain gzip is read single-threaded, because it can't be split. The docs recommend bgzip.
- **stdin:** `-` reads from stdin wherever a file is accepted, for at most one input of a run: `-` given for two inputs (`-a - -b -`) is a usage error that names both, since the second would read nothing.

### 5.2 Output

- The default output is stdout, or `-o PATH`.
- **Format and compression** come from the `-o` extension (`.bed`, `.bed.gz`, `.vcf`, `.vcf.gz`, `.bcf`, `.tsv`, `.tsv.gz`, `.bedgraph.gz`, and `.bw`, `.bigwig`, `.bb`, `.bigbed`: `convert` and `genomecov` only; a regular file only (D171)), or from `-O FMT` for stdout. `.gz` output is always **BGZF**, so it can be indexed.
- `--write-index` creates `.tbi`. It switches to `.csi` automatically for chromosomes over 512 Mbp and always uses `.csi` for BCF.
- **Atomic file output:** gloci writes to `PATH.tmp-XXXX` and renames on success. On failure the temp file is removed and no partial file is left behind. **[my call]** The temp file is not forced to disk before the rename (no `fsync`) **[owner request, 2026-11-17]**: the rename is what keeps a failed or interrupted run from leaving a partial file, and surviving an OS crash or power loss just after the run is left to the file system, as samtools, bcftools and bedtools leave it. The flush added up to 0.7 s to `sort -o` of the 58 MB GENCODE exons on a laptop disk (0.56 to 1.18 s with it, 0.35 to 0.60 s without, alternating runs) and costs a server round trip on NFS or Lustre.

---

### 5.3 Built-in genomes **[you]**

gloci ships chromosome names and sizes for common assemblies, so `-g hg38` works without a file. Passing a file still works.

- **How `-g` resolves.** A value containing `/`, or naming an existing file, is a file. Otherwise it is looked up as a built-in name (case-insensitive). An unknown name is an error listing the built-ins, with a "did you mean" suggestion, that points at `gloci genomes list` for their species and assemblies; `-g`'s error also names the file forms, `gloci genomes show`'s does not, since `show` takes only a built-in id. **[owner request, 2026-11-16]**
- **Naming styles are separate entries**, because names must match exactly (Section 4.5):
  - UCSC ids use `chr`-prefixed names.
  - GRC/Ensembl ids use the Ensembl names (`1 … 22, X, Y, MT`).
  - Each entry is built from its own project's published tables. They are not renamed copies: hg19 `chrM` (16,571 bp) is not GRCh37 `MT` (16,569 bp).
- **Order.** Built-ins list chromosomes in karyotype order (Section 4.5, `--chrom-order karyotype`), so `-g NAME` adds lengths and the exact chromosome set but never reorders. The exception is genomes with Roman-numeral chromosome names (yeast `sacCer3`/`R64`, worm `ce11`/`WBcel235`). These list their Roman-numeral chromosomes in numeric order (I, II, … XVI), followed by the rest.
- **Initial list:**
  - Human: `hg19`, `hg38`, `hs1` (T2T-CHM13v2.0), `GRCh37`, `GRCh38`.
  - Mouse: `mm10`, `mm39`, `GRCm38`, `GRCm39`.
  - Rat: `rn7`.
  - Zebrafish: `danRer11`, `GRCz11`.
  - Fruit fly: `dm6`, `BDGP6`.
  - *C. elegans*: `ce11`, `WBcel235`.
  - Budding yeast: `sacCer3`, `R64`.
  - *Arabidopsis*: `TAIR10`.
- **What each entry contains.** The source URL, retrieval date and checksum of each entry are recorded in the repository.
  - UCSC entries contain every sequence in the UCSC table, including alts, randoms and unplaced contigs, each with its real length.
  - Ensembl-named entries use Ensembl's primary assembly: chromosomes plus unlocalized and unplaced contigs, with LRG regions excluded. They leave out haplotype and patch regions, because Ensembl stores those as N-padded copies the length of a whole chromosome. Those copies would distort genome-wide operations such as `complement -g` or `windows -g`. A file that uses such contigs can pass a `.fai` instead. **[revised 2026-11-07]**
- **The data is compiled into the binary** (a few hundred KB), so installing and building need no network. A maintainer script refreshes it.
- **Built-in chains for `gloci liftover`** (D174, R13) follow the same rule for `--chain`: a value with `/` or naming an existing file is a file, otherwise a built-in id (case-insensitive), and an unknown id is exit 2 with "did you mean". There are six, all human: `GRCh37-to-GRCh38` and `GRCh38-to-GRCh37` (Ensembl, no restrictions), and `GRCh38-to-CHM13v2`, `CHM13v2-to-GRCh38`, `GRCh37-to-CHM13v2` and `CHM13v2-to-GRCh37` (T2T consortium, CC0). They are the published files, gzip -9, in the `gloci-chains` and `gloci-chains-t2t` crates (about 1 MB and 8.5 MB), with their source, hash and licence in `data/SOURCES.toml`. No other species has one: no freely licensed GRCm38-to-GRCm39 chain exists. Mitochondrial names are looked up and written as for a chain file, so `GRCh38-to-GRCh37` writes GRCh37's `MT` as `chrMT` for a `chr` input, never as `chrM`; `GRCh37-to-GRCh38` leaves `chrM` and `M` unmapped with a warning, because hg19 `chrM` is not GRCh37 `MT` (see above); the T2T chains have no mitochondrial chain. **[owner, 2026-11-17; amended 2026-11-17]**
- **`gloci genomes list`** lists the built-ins: id, species, assembly, naming style, sequence count and total length, tab-separated with no header line. `gloci genomes` alone prints the same list; `list` is there so the listing shows under Commands in `gloci genomes --help` and has a page in the command reference **[owner request, 2026-11-16]**. `gloci genomes show NAME` prints one as a 2-column `.genome` file.

## 6. CLI conventions

### 6.1 Global options (where relevant to a command)

| Option | Meaning |
|---|---|
| `-o, --output PATH` | output file (default stdout) |
| `-O, --output-format FMT` | `bed`, `bed.gz`, `vcf`, `vcf.gz`, `bcf`, `tsv`, `tsv.gz`, `bedgraph`, `bedgraph.gz`; with alignments `bam`, `sam`; for FASTA `fasta`, `fasta.gz` (a `.cram` output is a usage error suggesting `samtools view -C`) |
| `--write-index` | write `.tbi`/`.csi` next to the output |
| `--compression-level N` | BGZF level (default 6) |
| `-t, --threads N` | default: `available_parallelism()`, which honours cgroups and CPU affinity (so SLURM allocations); env `GLOCI_THREADS` |
| `-g, --genome NAME\|FILE` | chromosome sizes and order: a built-in genome name (`hg38`, `GRCh38`, `mm39` …, see `gloci genomes list`) or a file (`.fai`, `.genome`, `.dict`, VCF/BCF/BAM header) The one exception is `groupby`, where `-g COLS` names the group columns, as in bedtools groupby; `groupby` takes no genome. **[owner, 2026-11-10: keep this exception]** |
| `--chrom-order karyotype\|lex` | chromosome order when neither `-g` nor the header of a VCF/BCF or alignment input sets one (Section 4.5), default `lex` |
| `--vcf-span allele\|ref` | variant footprint rule (Section 4.2), default `allele` |
| `--unsorted` | process unsorted input in its own order via the in-memory index; results follow A's order (Section 4.5) |
| `--sort` | sort the inputs first with bounded memory (`--max-mem`), then run the streaming pass; results in sorted order (Section 4.5) |
| `--region chr:s-e` (repeatable), `--regions-file BED` | restrict work to regions; uses the index when present |
| `--chrom-map FILE` | rename chromosomes on input (Phase 2) |
| `--header` | carry BED/TSV header lines from A to the output |
| `--min-mapq N`, `--require-flags FLAGS`, `--exclude-flags FLAGS` | read filters, applied first to every alignment input of `intersect`, `coverage`, `genomecov`, `convert` and `view`. FLAGS is a number or comma-separated samtools names. The defaults filter nothing, as in bedtools; `--exclude-flags 0x704` gives `samtools depth` and mosdepth semantics. Counts are reported at `-v`. **[D102]** |
| `--reference FASTA` | the FASTA a CRAM was made with (bedtools `--cram-ref`), with its `.fai` (and `.gzi` when BGZF). There is no `REF_PATH`, cache directory or download; a CRAM written without a reference needs none, and neither does `gloci index x.cram`, where `--reference` only checks lengths. It warns when given with no CRAM input. **[D109]** |
| `--split-deletions` | for alignments, also cut reads at CIGAR `D` (`--split` cuts at `N` only); a usage error with no alignment input. **[D105]** |
| `--skip-invalid` | skip malformed records instead of failing; prints a count and the first 10 examples to stderr |
| `--tmp-dir DIR` | temp dir for `sort` and atomic outputs; env `GLOCI_TMPDIR` |
| `-q, --quiet` / `-v, --verbose` | logging to stderr |

### 6.2 Input naming

- **Two-input relations** take `-a FILE` (the records being reported) and `-b FILE...` (one or more).
  - With several B files, `--names` labels them in the output. Otherwise gloci uses the file basenames.
- **N-ary commands** (`merge`, `partition`, `sort`, `complement`) take positional files.

### 6.3 Aggregation syntax (merge, map, groupby, cluster)

- Use `-c COL:OP`, repeatable or comma-separated, for example `-c 5:mean,5:max -c 4:distinct`.
  - `COL` is a 1-based column number, or a name: a header column, or for VCF `QUAL`, `ID`, `INFO/<KEY>`, `FILTER`.
  - This replaces bedtools' two parallel lists `-c 5,5 -o mean,max`, which are easy to misalign, and frees `-o` for output.
- **Ops:**
  - `count`, `count_distinct`
  - `sum`, `mean`, `median`, `min`, `max`
  - `stdev`, `var` (sample), `stdev_pop`, `var_pop`, `cv`, `mad`
  - `quantile:Q`, `tmean:LO-HI`, `wmean` (weighted by overlapping bases), `wsum` (value × overlapping bases)
  - `mode`, `antimode`, `first`, `last`
  - `collapse`, `distinct`, `distinct_sorted`
  - `bases` (overlapping bases), `bases_uniq`, `frac_a`, `frac_b`
  - `min_record`, `max_record` (the whole B record with the min or max value)
  - `records` (all B records), `indicator` (0/1)
- `--empty STR` sets the value printed when nothing overlaps (default `.`). `--delim` sets the list separator (default `,`).

### 6.4 Behaviour

- **Exit codes:** 0 success, 1 data or runtime error, 2 usage error.
- **Broken pipe** (e.g. `| head`) exits quietly with status 0. There is no panic and no "Broken pipe" message.
- **Errors** say what, where and how to fix it. Example: `error: calls.vcf.gz line 1,204,331: POS "12a4" is not an integer`.
- **Determinism:** output is byte-identical for every `--threads` value.
- **Help:** every command's `--help` has 2–4 real examples. Shell completions and man pages are generated from the CLI definition.
- **Help layout [owner, 2026-11-14].** `gloci --help` lists the commands in groups, the stable core first and then the experimental groups (`Stable core`, then headings such as `Intervals (experimental)`; `command_groups` in `crates/gloci/src/cli/mod.rs`), each with a one-line summary of at most 64 characters that names the tool it replaces (`Report A records that overlap B (bedtools intersect)`), and a line under the list saying what the two tiers mean. A command's `--help` has the layout: a summary paragraph (at most 12 lines at 100 columns, naming every tool the one-line summary names), the usage, the options, an `Examples:` block, and headed sections for the detail. Tests in `crates/gloci/tests/cli_basics.rs` enforce the grouping (every command is in exactly one group), the 64 characters, the 80-column list, the summary length and the examples.
- **Stability tiers** (added in commit bcbf8672, 2026-11-13; no decision in `docs/decisions.md` records them). The stable core is `genomes`, `index`, `sort`, `view`, `intersect`, `merge`, `subtract`, `closest` and `coverage` (`CORE_COMMANDS` in `crates/gloci/src/lib.rs`): their options and output change only with a version bump and a changelog entry. Every other command is experimental and its help starts with a line saying so; a test fails when a command is neither in the core nor marked. A command is promoted when it has been run beside the tool it replaces on real data (each difference fixed or listed in `docs/differences.md`), has a large-input test with a memory ceiling in the nightly run, and has had its option names reviewed once more (`docs/book/src/stability.md`). `contrib/gloci-shadow.sh` runs a step with the old tool and with gloci and logs whether the outputs are identical.

---

## 7. Commands

Phase markers are explained in Section 10.

### 7.1 Interval relations (A vs B)

**`gloci intersect`** — P1. Replaces bedtools `intersect` and `window`, bedops `-e`/`-n`, vcftools `--bed`/`--exclude-bed`, and `vcf-isec`.

- **`--mode`** sets what is reported. The default is `any` **[my call]**, because it's the most common intent and never duplicates records. It has no short flag: `-m` is `sort --merge-sorted`, and one flag keeps one meaning. **[my call, revised 2026-11-08]**
  - `any`: each A once if it overlaps ≥1 B.
  - `none`: A with no overlap.
  - `pairs`: one line per A–B pair, A columns then B columns. `--overlap-bases` adds the overlap length. `--keep-unmatched` keeps A rows with no partner, padded with `.`.
  - `clip`: the overlapping part of A, once per overlap. This is the bedtools default.
  - `count`: A plus the number of overlapping B.
- **Overlap criteria:** all of Section 4.3, meaning `-f -F -r -e --min-bases -s -S --split -w --left --right --strand-window --inclusive-points`.
- **VCF vs VCF:** `--match overlap|alleles`.
  - `overlap` (default) uses footprints.
  - `alleles` requires the same chromosome, POS, REF and at least one shared ALT, compared after the trimming in Section 4.2.

**`gloci subtract`** — P1. Replaces bedtools `subtract` and bedops `-d`.

- **Default:** removes the bases B covers from each A. An A may split into several pieces, each keeping A's other columns.
- `--whole`: removes an A entirely if any overlap meets the criteria. This is bedtools `-A`.
- `--sum-overlaps`: with `-f`, adds up coverage from all B instead of testing each B alone. This is bedtools `-N`.
- If A is VCF, records are removed whole (Section 4.6).
- `--match overlap|alleles` on VCF/BCF inputs, as `intersect`. With `alleles`, an A record is dropped only by an overlapping B that shares an allele (trimmed POS, REF and ALT, Section 4.2). A multi-allelic A is dropped whole when any ALT is shared, with no allele trimming of the record: an owner decision. The result equals `intersect --mode none --match alleles`. **[owner addition, Task 21]**
- Supports the same overlap criteria as `intersect`.

**`gloci closest`** — P1. Replaces bedtools `closest` and bedops `closest-features`.

- `-k N`: the N nearest.
- `--ties all|first|last`: default `all`.
- `--distance`: appends a distance column. `--signed ref|a|b` sets the sign convention.
- `--ignore-overlaps`, `--ignore-upstream`, `--ignore-downstream`.
- `--max-distance N`.
- `--each-side`: reports the nearest upstream and nearest downstream separately, like closest-features.
- `--per-file`: with several B files, reports the nearest feature in each file.

**`gloci map`** — P1 core ops, P2 full op set. Replaces bedtools `map` and bedops `bedmap`.

- Echoes each A, then appends one column per `-c COL:OP` computed over the overlapping B (Section 6.3).
- Supports the same overlap criteria as `intersect`, including `-w`, which is bedops `--range`.
- `--skip-unmatched` drops A with no overlap.
- All requested aggregations are computed in one pass.

**`gloci coverage`** — P1. Replaces bedtools `coverage` and `annotate`.

- For each A it reports the number of overlapping B, the bases covered, A's length, and the fraction covered.
- `--hist`: depth histogram per A plus an `all` summary.
- `--per-base`: depth at each base.
- `--mean`: mean depth.
- With several `-b` files, appends one fraction column per file, which replaces `annotate`.

### 7.2 Region algebra (computed regions out)

**`gloci merge FILES...`** — P1. Replaces bedtools `merge` and bedops `-m`.

- Inputs are k-way merged as one stream.
- `-d, --max-gap N`: the default 0 merges overlapping and bookended intervals. A negative N requires at least |N| bases of overlap.
- `-s`: merge per strand.
- `-c COL:OP`: aggregations.
- `--count`: adds the number of merged records.

**`gloci complement FILES...`** — P1. Replaces bedtools `complement` and bedops `-c`.

- With `-g`, reports uncovered regions across whole chromosomes, including chromosomes with no features. This is bedtools behaviour.
- Without `-g`, reports only the gaps between features. This is bedops behaviour.

**`gloci partition FILES...`** — P2. Replaces bedops `-p`, `-i` and `-s`, plus bedtools `multiinter` and `unionbedg`.

- Splits all inputs at every boundary into disjoint segments.
- Each segment line carries the number of files covering it and their names.
- `--min-files K` / `--max-files K` / `--all`: filter segments by how many files cover them.
- `--merge-adjacent`: joins neighbouring segments with the same membership.
- `--no-membership`: coordinates only.
- `--bedgraph`: inputs are bedGraphs and the output has one value column per file. `--empty 0` fills gaps, and `-g` with `--empty` covers the whole genome.

How the old commands map onto `partition`:

| Old | New |
|---|---|
| bedops `-i` | `partition --all --merge-adjacent --no-membership` |
| bedops `-s` | `partition --max-files 1 --merge-adjacent --no-membership` |
| bedops `-p` | `partition --no-membership` |

### 7.3 Single-input transforms

| Command | Replaces | Phase | Key options |
|---|---|---|---|
| `slop` | bedtools slop, bedops `--range` | P1 | `-b/--both N`, `--left N`, `--right N`, `-s` (strand-aware), `--pct`; clips at 0 and, with `-g`, at the chromosome end |
| `flank` | bedtools flank | P1 | same as slop |
| `windows` | bedtools makewindows, bedops `--chop` | P1 | `-g` or `-i BED`; `-w SIZE`, `--step S`, `-n COUNT`, `--reverse`, `--id-column` |
| `shift` | bedtools shift | P2 | `--by N`, `-s`, `--pct`, clipping as slop |
| `cluster` | bedtools cluster | P2 | `-d/--max-gap`, `-s`; appends a cluster id |
| `spacing` | bedtools spacing | P2 | appends the distance to the previous feature |
| `groupby` | bedtools groupby | P2 | `-g COLS -c COL:OP` over consecutive equal keys |
| `genomecov` | bedtools genomecov | P2 (BED), P3 (BAM, CRAM, SAM) | `--bedgraph`, `--bga`, `--hist` (default), `--per-base`, `--scale X`, `--strand +/-`, `--split`, `--split-deletions`, `--5p/--3p`, `--max-depth N`, `--sort` |
| `sample`, `random`, `shuffle` | bedtools | P3 | `--seed` required for reproducibility; `shuffle --exclude/--include BED` |
| `liftover` | CrossMap `bed`, `gff`, `region`; UCSC `liftOver` | after P3 (D174) | `--chain`, `--unmapped`, `--chrom-style`, `--ignore-blocks`, `--merge-segments`, `--min-ratio`; input order out |

**`liftover`** **[D174]**

- A record is lifted through the aligned blocks of a chain, by CrossMap's rules, not UCSC `liftOver`'s: a record that crosses a chain gap is split into one record per piece (R5), and a record that loses bases to a gap is written shortened. The command warns once about both, as CrossMap does neither aloud.
- Records that do not lift go to `--unmapped` as read, or are counted by reason in one warning (R6). The output is in input order and may be unsorted, so pipe it through `gloci sort`.
- A BED record of 12 or more columns is lifted block by block (R8): each block needs one piece, the blocks are rebuilt in input order and the thick ends follow their blocks; `--ignore-blocks` lifts columns 1 to 3 only (CrossMap `--naive-bed-parsing`).
- `--merge-segments [--min-ratio F]` (R10, CrossMap `region -r`) writes a record cut at chain gaps as one record from the first start to the last end, when at least F of it maps and its pieces share a chromosome and a strand; it is not for GFF, nor for BED12 without `--ignore-blocks`.
- The chain is held in memory as an interval index per source chromosome (R1); records stream.

**`genomecov` on alignments (BAM, CRAM, SAM)** **[D117]**

- `-g` is optional: the default is the header's chromosomes and lengths, in header order, empty chromosomes included. `-g` overrides it, where bedtools ignores `-g` with `-ibam`.
- A read counts its whole reference span, deleted bases included. `--split` cuts at `N` only; `--split-deletions` also cuts at `D` and is bedtools `-ibam -split`; plain `-ibam` (cuts at `D`, not `N`) has no gloci counterpart (docs/differences.md#bam-genomecov-deletions).
- `--5p`/`--3p` take the first or last aligned base by strand, as bedtools. Unmapped reads are left out, and the read filters (Section 6.1) apply first.
- An indexed BAM, or a CRAM with a `.crai`, at `-t` above 1 runs per chromosome, with `genome` rows summed and output identical to `-t 1`.
- `-pc`, `-fs` and `-du` are not planned (no flag in this table, and a Section 2 non-goal).

**Sort order of transform output:**

- Transforms that keep order intact, such as uniform `slop` or `flank` without `-s`/`--pct`, produce sorted output.
- Transforms that can break order say so in `--help` and warn when that actually happens.
- How to fix order cheaply (a bounded reorder buffer or a recommendation to run `sort`) is left to the builder.

### 7.4 Sorting, indexing, viewing, conversion

**`gloci sort FILES...`** — P0. Replaces bedops `sort-bed`, `vcf-sort`, and `sort -k1,1 -k2,2n`.

- **Order:** chromosome order (Section 4.5), then start, then end.
  - Ties keep input order (a stable sort).
  - `--full-key` breaks ties on the whole line instead.
- **Input:** BED, VCF/BCF (keeping the header; chromosome order comes from `##contig` unless `-g` is given), GFF, and alignments (below).
- **External merge sort:**
  - Memory is limited by `--max-mem` (default 2G), with chunks sorted in parallel and spill files in `--tmp-dir` (compressed with a fast codec).
  - Inputs larger than RAM must work.
- `-m, --merge-sorted`: k-way merge of already sorted inputs (bedops `-u` / `--everything`, `sort -m`).
- `--check`: checks order without writing anything. Exits 1 and prints the first out-of-order pair if unsorted.
- `--unique` / `--duplicates`: keep only unique lines, or only duplicated ones.
- Output can be BGZF with `--write-index`.

**`gloci sort` for alignments** (one BAM, SAM or CRAM, or stdin) follows samtools 1.21 `bam_sort.c`. **[D121]**

- **Coordinate order** (the default): reference in header order with unplaced reads last, then position, then the reverse-strand bit, then input order.
- **`-n`** sorts by read name in samtools' natural order (digit runs compare as numbers after leading zeros, so `r2` precedes `r10` and `r01` equals `r1`), and **`-N`** by plain byte order (`--by-name`, `--by-name-bytes`). In both, ties go READ1 before READ2 and primary before supplementary before secondary, then input order. `-n` and `-N` are the owner's flags, an exception to one flag, one meaning: `windows -n COUNT` keeps its own meaning, as `groupby -g` does. **[owner exception]**
- **`@HD`:** `SO` is set (`coordinate`, or `queryname` with `SS:queryname:natural` for `-n` and `SS:queryname:lexicographical` for `-N`); a coordinate sort removes `SS`; `GO` is removed; with no `@HD`, `@HD VN:1.6 SO:…` is written first. gloci adds an `@PG` line (`ID:gloci`, with a suffix when taken).
- **Output format** follows the input's: SAM for SAM, BAM for BAM and CRAM (CRAM is never written); `-o` and `-O` override. samtools always writes BAM. **[D121]**
- Every record is kept, unmapped reads included. `--write-index` works on coordinate-ordered BAM output only (`.bai`, or `.csi` for long references). `--check` checks the chosen order. `-g`, `--full-key`, `--unique` and `--duplicates` are usage errors with alignments: the order is the header's.
- Memory is bounded by `--max-mem`, a total for the whole sort (samtools' `-m` is per thread); spill runs hold raw records.

**`gloci sort -m` for alignments is `samtools merge`.** **[D122]** `-m` already means "merge sorted inputs" in gloci, and `gloci merge` keeps its interval meaning (one name, one meaning); its alignment error points to `sort -m`.

- The keys are the sort's (`-n`/`-N` select name order). Ties go to input file order, then record order (samtools' `heap_lt`). Each input is order-checked as it is read, and an input whose `@HD SO`/`SS` contradicts the mode is an error up front.
- `@SQ` lists must be identical (names, lengths, order), else an error naming the first difference. samtools appends unknown `@SQ` lines instead. `@HD` comes from the first input with `SO`/`SS` set; `@CO` lines come from every input; `@RG` and `@PG` lines with the same ID are combined (an identical line is kept once, as `samtools merge -c -p`), and a same-ID `@RG` with different fields is an error. samtools' default random-suffix renaming is not reproduced; `-r` is not supported.
- The output follows the inputs' format: SAM when every input is SAM, otherwise BAM; `-o`/`-O` choose.
- Many inputs are merged in cascade rounds: inputs are opened per round, at most fan-in at a time (128, fewer under a small `--max-mem`), with the decode threads shared within `-t`, and the rounds' runs go through `--tmp-dir`. A pipe input is spooled once. Open files and memory stay bounded.

**`gloci cat FILE...`** — P3a. Replaces `samtools cat` for BAM. Its meaning is `vcf concat --naive`'s (a byte-level join), not `vcf concat`'s checked join, so it takes samtools' name. **[D123]**

- The header is the first input's plus the `@RG` lines of later inputs with new IDs, and an `@PG` line. Every input's `@SQ` list must equal the first's, else an error (samtools does not check, and the records silently point at the wrong reference).
- The block holding each later file's header end is recompressed from there; every later block is copied verbatim; end-of-file markers are dropped and one is written at the end. The order of the records is not checked.
- BAM only: CRAM is refused (suggesting `samtools cat`) and SAM is refused (suggesting `gloci view x.sam -o x.bam`). `-` reads one input from stdin. `--write-index` writes the `.bai` for coordinate-sorted output.

**`gloci index FILE`** — P0, extended in P1 and P3. Replaces `tabix`, `bcftools index`, `samtools index` and `samtools faidx` for the formats gloci reads. **[you]**

- `gloci index FILE` chooses the index type from the file's content, with no format flags. The one exception is arbitrary TSV (below).
- Index types by phase:

| Input | Index | Replaces | Phase |
|---|---|---|---|
| bgzipped VCF, BED, bedGraph, GFF/GTF | `.tbi`, or `.csi` when an end exceeds 2^29 | `tabix -p vcf\|bed\|gff` | P0 |
| BCF | `.csi` | `bcftools index` | P0 |
| any bgzipped tab-separated file | `.tbi`/`.csi`, using `-s/--seq-col`, `-b/--begin-col`, `-e/--end-col`, `--zero-based`, `--skip-lines N`, `--comment CHAR` | `tabix -s -b -e` | P1, so TSV outputs of `map`, `coverage` and `intersect --mode pairs` can be indexed |
| BAM | `.bai`, or `.csi` with `--csi` or when a reference exceeds 2^29 (`--csi` is for BAM) | `samtools index` | P3a |
| CRAM | `.crai`, with samtools' entries; needs no reference, even for multi-reference slices (`--reference`, if given, only checks lengths) **[D107, D109]** | `samtools index` | P3a |
| FASTA, bgzipped FASTA | `.fai`, plus `.gzi` for BGZF | `samtools faidx` | P3, with `getfasta`/`maskfasta`/`nuc` |

- Each new index type ships with differential tests showing that the reference tool (htslib, samtools) reads gloci's index and returns identical query results (Section 9.3).
- Indexing formats gloci reads is in scope; becoming a general samtools replacement is not (Section 2).

**`gloci view FILE [REGION...]`** — P0. Replaces bedops `bedextract`, `tabix` region queries, and vcftools `--recode` as a pure format conversion.

- Region queries use the index if one exists. A bigWig or bigBed is always read through its own index: bigWig records are clipped to each region, bigBed records are whole and kept once when they overlap a region under Section 4.3. A zero-length region (a BED region with start equal to end; `chr1:5-5` is the one base 5) on a bigWig yields nothing, because a record clipped to no bases is empty. **[D171]**
- On **uncompressed** BED without an index, gloci scans the file and checks its order, so an unsorted file is an error and never a wrong answer. With `--assume-sorted`, regions are found by binary search on the file instead, like bedextract. This is fast, but it checks only the chromosome block it reads, so it trusts the user that the file is sorted. Indexing (`gloci sort … -o x.bed.gz --write-index`) gives both speed and safety. **[you: correctness first; revised 2026-11-07]**
- Also handles format conversion (`-O vcf.gz|bcf|...`), `--header-only` and `--no-header`.
- **`view` on alignments (BAM, SAM, CRAM)** **[D116]** converts BAM↔SAM and CRAM→BAM/SAM (as `-O` or the `-o` extension say; the input's format by default, BAM for CRAM) and writes every record, unmapped reads included (D120), with an `@PG` line added. `--header-only`, `--no-header` (SAM output only) and the read filters apply. Regions are read through the index (`.bai`, `.csi`, `.crai`) or by scanning, with a warning, and select by Section 4.3 on the D101 span; a placed unmapped read spans `[pos, pos+1)`, as samtools finds it. The `.bai` and `.crai` answers for a region past a reference's `LN` differ as samtools' do. It is not `samtools view`: there is no `-b`, `-h`, `-c` or `-L`. **[D116]**

**`gloci genomes [list | show NAME]`** — P0. Lists the built-in genomes, or prints one as a `.genome` file (Section 5.3). `list` was added on 2026-11-16; `gloci genomes` alone still lists.

**`gloci check FILE...`** — P0.

- Validates format, coordinates (start ≤ end, and within chromosome lengths when `-g` is given), sortedness and duplicates.
- Prints a one-screen summary: record count, chromosomes, and the first problems.
- Exits 1 if there are errors.

**`gloci convert`** — P1. Replaces bedops `convert2bed`, `vcf2bed`, `gff2bed`, `gtf2bed`, and bedtools `bed12tobed6`.

- VCF to BED (default columns `chrom start end ID REF ALT QUAL FILTER`):
  - The footprint follows `--vcf-span`.
  - `--split-alleles` writes one row per ALT.
  - `--types snv,mnv,ins,del,sv,ref,other` filters by type. This replaces vcf2bed `--snvs`/`--insertions`/`--deletions`.
- GFF/GTF to BED: `--feature exon`, and `--name-attr gene_name` for the name column.
- BED12 to BED6: `--blocks`.
- bigWig and bigBed **[D171]**: `-o x.bw` writes a bigWig from bedGraph, a 4-column BED, WIG or bigWig, and `-o x.bb` a bigBed from BED or bigBed (`bedGraphToBigWig`, `wigToBigWig`, `bedToBigBed`), with the chromosome sizes of `-g` or of the bigWig or bigBed converted; `--autosql FILE` supplies a bigBed's autoSql. Each chromosome's records must be together and sorted by start, a bigWig takes no overlaps and stores 32-bit floats (`#bigwig-float32`), and the output is a regular file, written whole or not at all. A bigWig or bigBed input converts to bedGraph or BED (`bigWigToBedGraph`, `bigBedToBed`).
- WIG to bedGraph (P2). BAM/CRAM/SAM to BED (bedtools `bamtobed`: BED6 by default, `--bed12`, `--blocks`, `--split-deletions`, `--score-tag`, `--cigar`, `--color`; **[D106]**), PSL and RepeatMasker to BED (**[D114]**, columns of bedops 2.4.41 `psl2bed` and `rmsk2bed`) (P3).
- **An alignment as an interval [D101]:** the span is `[pos, pos + reflen)`, reflen being the `M D N = X` lengths, with a read that has no reference operations counting one base everywhere (convert, intersect, genomecov); clips never widen it; the strand is `-` for flag `0x10`. A read at or past its reference's `LN` is skipped and counted, and one that runs past it is cut at `LN`, with one warning per input. Unmapped reads are skipped by every interval command. Errors name the place: SAM `line N`; BAM `record N of CHROM (record M of the file)`; CRAM `record N of CHROM`, without the file-wide number because a `.crai` has no counts. A read through an index for a region, which starts mid-file, is `record at BGZF offset X:Y` (its virtual offset). `genomecov -g` names a read past the end of its `-g` chromosome as `record N of CHROM (counting only the reads genomecov counts)`, so the place is the same streamed and per chromosome.
- **Alignments in output [D104]:** `intersect --mode any|none` writes BAM for a BAM or CRAM A and SAM for a SAM A (`-o`/`-O` choose); asking for BED is a usage error suggesting `-O bam | gloci convert -`. Where a command prints an alignment's columns it prints BED6, and an unmatched B placeholder for an alignment is BED6 wide. `--mode clip` of reads is a usage error. `--chrom-map` is refused when `intersect` writes a SAM A as alignments (`--mode any|none`), because the reads keep their names.

**P3c samtools essentials** **[D125–D141, owner 2026-11-10]**: `idxstats`, `flagstat`, `reheader`, `faidx`, `dict`, `head`, `samples`, `flags`, `depth` and `bedcov` keep samtools' names (no gloci command has those meanings); `quickcheck` is `check --quick`, `view -c` is `view --count` (`-c` is `COL:OP`), `view -L` is `--regions-file`. samtools' short flags are not copied; each maps to a long flag in the migration table. Exit status is §6.4's, not samtools' (D126). CRAM record decoding needs `--reference` (D109); `samples`, `reheader` and `check --quick` never do. `depth` and `bedcov` inputs need identical `@SQ` lists.

- **`idxstats`** **[D127]**: from the index's pseudo-bins when a BAM has one, else a scan with one warning. An index lacking `n_no_coor`, or a reference's pseudo-bin, counts as "no counts": warn and scan.
- **`flagstat`** **[D128]**: samtools' counts and layouts (text, TSV, JSON) byte for byte; `-O` is the output format.
- **`check`** **[D129]**: reads every alignment record; `--quick` is `quickcheck`: header parses, first record parses (CRAM: first container), EOF marker. `--list-failed` prints failing paths; any failure exits 1. `--chrom-map` is refused for BAM/CRAM/SAM (exit 2); `-g` and `--chrom-order` warn once. The `--quick` read is header + first record + 28 bytes (BGZF), or header + first container + 38 bytes (CRAM).
- **`reheader`** **[D130]**: BAM or CRAM, records are copied; the new `@SQ` list must keep every old reference in place with the same `LN`, else exit 1. Memory: old + new headers + one 64 KiB block.
- **`faidx`, `dict`** **[D131, D132]**: need `gloci index` first (`faidx`); `dict` streams.
- **`view --count`, `--regions-file`** **[D133]**: a record is written or counted once however many regions it overlaps. `--regions-file` combined with positional regions is a union.
- **`head`** **[D134]**: all header lines, no records by default. A BAM with no `@SQ` text gets synthesised `@SQ` lines (samtools 1.21 prints nothing). A VCF header is printed as stored, with no PASS line. Text formats are read in batches of up to 16384 records (a fixed record batch).
- **`samples`, `flags`** **[D135, D136]**: FLAGS is a decimal, `0x` hex or comma-separated names; a leading-zero decimal and trailing text are usage errors. §6.1's `--column-names` is per command.
- **`depth`, `bedcov`** **[D137–D140]**: exclude `UNMAP,SECONDARY,QCFAIL,DUP` by default (D102). `depth` reads its inputs as one merged stream at every `-t`; threads decode BGZF/CRAM (D138). `--skip-mate-overlap` is samtools `-s`. A `bedcov` region of zero length counts no reads, as samtools 1.21 does. `chr:-E` is `chr:1-E` (D131, §4.1).
- **`depth --summary`, `depth --out-prefix`** **[D173, owner 2026-11-17]**: samtools coverage's table, and mosdepth 0.3.10's files under their names; both stream on the depth core, and the files are published together with CSI indexes.
- **CRAM 2.x** is refused by name and version wherever the prelude is read (`#cram-limits`).

### 7.5 Sequence (P3)

`getfasta`, `maskfasta` and `nuc` work as in bedtools, with strand-aware reverse-complement, `--split` for BED12, and soft or hard masking. `getfasta` and `nuc` need FASTA plus `.fai`; bgzipped FASTA plus `.gzi` is also supported. gloci never builds an index implicitly: a missing one is a usage error naming `gloci index FILE` **[D108]**. Their flag names are:

- **`getfasta FILE --fasta REF`** **[D111]**: output in input order; the header is `chrom:start-end` (0-based start, as bedtools), `--name` gives `name::chrom:start-end` and `--name-only` gives `name`; `-s` appends `(+)`/`(-)`/`(.)` and reverse-complements `-` records with bedtools' IUPAC-aware, case-keeping table; `--rna` writes U; `--split` joins BED12 blocks; `--tab` writes `header TAB seq`; `--bed-out` writes the record, then the sequence. A record outside the FASTA is a data error (bedtools skips it), and `--skip-invalid` skips and counts. Zero-length records are skipped with one counted warning. `--chrom-map` renames the intervals only, and VCF follows `--vcf-span`. There is no `-fullHeader`.
- **`nuc FILE --fasta REF`** **[D112]**: TSV of the record's columns, then `pct_at pct_gc num_A num_C num_G num_T num_N num_oth seq_len`, then the sequence (`--seq`) and the pattern count (`--pattern P`, overlapping matches, case-sensitive unless `--ignore-case`/`-C`), after bedtools' `#1_usercol …` header line. Counts are case-insensitive and `U` counts as T; percentages are computed in single precision and printed with 6 decimals, as bedtools. `-s` reverse-complements `-` records. Errors as `getfasta`.
- **`maskfasta --fasta IN -b MASK... [-o OUT]`** **[D113]**: streams the FASTA (any compression, no index needed). Masked bases become `--mask-char C` (default `N`; a space is refused) or lowercase (`--soft`). Headers and line breaks are kept exactly, where bedtools re-wraps and drops descriptions. Masks may come in any order and are merged per chromosome; a chromosome the FASTA lacks gets one warning, and a mask past its sequence's end is a data error. `.fa.gz` output is BGZF, and `--write-index` writes `.fai`/`.gzi`.

### 7.6 Interval statistics (P3)

`gloci stats jaccard|fisher|reldist` work as in bedtools, with clearly documented formulas.

### 7.7 Pairs (P3)

`gloci pair-intersect` covers bedtools `pairtobed`, and `gloci pair-pair` covers `pairtopair`. Both take BEDPE input.

### 7.8 VCF commands

All VCF commands share one **selection pipeline**, applied in this fixed, documented order **[my call; it matches vcftools semantics]**:

1. Sample selection: `--samples a,b`, `--samples-file F`, `--exclude-samples`, `--exclude-samples-file`.
2. Genotype masking. These set the genotype to missing; the record is not removed.
   - `--min-gq N`, `--min-dp N`, `--max-dp N`.
   - `--mask-filtered-genotypes` masks genotypes whose FORMAT/FT is not PASS.
3. Per-site statistics are recomputed from the remaining genotypes: allele counts, frequencies, missingness, mean depth.
4. Site filters, combined with AND:
   - **Regions:** `--region`, `--regions-file`, `--exclude-regions-file`, `--chrom`, `--exclude-chrom`.
   - **Positions and IDs:** `--positions F` / `--exclude-positions F` (chromosome and 1-based POS), `--ids F` / `--exclude-ids F`.
   - **Type:** `--types` / `--exclude-types` with `snv,mnv,ins,del,indel,sv,ref,other`. `--min-alleles N`, `--max-alleles N`, `--biallelic`.
   - **Quality and FILTER:** `--min-qual X`, `--pass-only`, `--keep-filter TAGS`, `--remove-filter TAGS`.
   - **INFO flags:** `--require-info FLAG`, `--exclude-info FLAG`.
   - **Frequency:** `--min-maf X`, `--max-maf X`, `--min-mac N`, `--max-mac N`.
   - **Missingness:** `--max-missing FRAC` (fraction of missing genotypes; see migration note), `--max-missing-count N`.
   - **Depth:** `--min-mean-dp X`, `--max-mean-dp X`.
   - **Other:** `--hwe P`, `--thin N` (keep sites at least N bp apart, first wins), `--phased-only`.
   - **Expressions (P2):** `--include EXPR` / `--exclude EXPR`, in an expression language modelled on bcftools so users already know it. Example: `QUAL>30 && INFO/DP>10 && GT="het"`.

Definitions:

- **MAF** is the frequency of the second most common allele among called alleles. For biallelic sites this equals `min(p, 1-p)`. **[my call]** Any difference from vcftools on multi-allelic sites goes in the differences list.
- **HWE** uses the exact test (Wigginton et al. 2005). It is only run for biallelic diploid sites; other sites pass.
- Mixed ploidy is supported. Missing alleles inside a genotype (`./1`) count as missing for that allele.

**`gloci vcf filter IN`** — P1, without expressions.

- Writes the selected records as VCF/BCF.
- **INFO is kept by default**. vcftools drops it unless you pass `--recode-INFO-all`. `--keep-info KEYS` / `--drop-info KEYS` prune it.
- `--recalc-ac` updates INFO/AC, AN and AF after subsetting samples. This is off by default, so INFO stays exactly as it was in the input. **[my call]**

**`gloci vcf stats IN --out-prefix P [stat flags]`** — P2.

- **Computes any number of statistics in one pass** and writes `P.<stat>.tsv` for each. vcftools needs one run per statistic, so this is a large saving on big files.
- Accepts every selection option above, so filtering and stats happen in one pass without a VCF round-trip.
- Flags and the vcftools options they replace:

| Flag | vcftools option |
|---|---|
| `--freq` | `--freq/--freq2` |
| `--counts` | `--counts/--counts2` |
| `--depth` | `--depth` |
| `--site-depth` | `--site-depth` |
| `--site-mean-depth` | `--site-mean-depth` |
| `--geno-depth` | `--geno-depth` |
| `--site-quality` | `--site-quality` |
| `--missing-indv` | `--missing-indv` |
| `--missing-site` | `--missing-site` |
| `--het` | `--het` |
| `--hardy` | `--hardy` |
| `--tstv` | `--TsTv-summary`, plus `--tstv-by-count`, `--tstv-by-qual` |
| `--filter-summary` | `--FILTER-summary` |
| `--snp-density BIN` | `--SNPdensity` |
| `--singletons` | `--singletons` |
| `--indel-hist` | `--hist-indel-len` |
| `--summary` | (new) a short human-readable overview, plus `--json` |

- Output columns use vcftools' column names wherever an equivalent exists, to ease migration. **[my call]**

**`gloci vcf extract IN`** — P2. Replaces `--get-INFO`, `--extract-FORMAT-info` and `--012`.

- `--info KEYS` writes TSV with one row per site.
- `--format KEYS --layout wide|long` writes per-sample values.
- `--genotype-matrix 012` writes `.012`, `.indv` and `.pos` files.

**`gloci vcf query IN -f FORMAT`** — built on the owner's request of 2026-11-17 (D172). Replaces `bcftools query`.

- bcftools 1.21's format language for the common fields; anything else is a usage error naming it.
- Every selection option runs through the shared stage, where `--include`/`--exclude` that read genotypes pick samples as well as records, as bcftools query does, while the other `vcf` commands keep or drop whole records (D71).
- Differences are in `docs/differences.md` (`#vcf-query-values` to `#vcf-query-format-strict`).

**`gloci vcf concat FILES...`** — P2. Replaces `vcf-concat`.

- Checks that samples, headers and order are compatible.
- `--naive` concatenates BGZF blocks without decompressing when headers are identical. This is much faster.

**`gloci vcf popgen IN`** — built on the owner's request of 2026-11-14. Replaces `--site-pi`, `--window-pi`, `--TajimaD`, `--weir-fst-pop`/`--fst-window-*`, `--geno-r2`, `--hap-r2`, `--interchrom-*`, `--ld-window*`, `--relatedness` and `--relatedness2`.

- Multiple statistics are computed in one pass where possible.
- Formulas are documented, and each is cross-checked against vcftools 0.1.16 by differential cases and committed goldens.
- **Statistics [D158–D162, owner 2026-11-14].** Nucleotide diversity per site (`--site-pi`) and in windows (`--window-pi`, `--window-pi-step`), Tajima's D in bins (`--tajima-d`) (D159); Weir and Cockerham's Fst for two or more populations (`--fst-pop FILE`, once per population), per site, in windows (`--fst-window-size`, `--fst-window-step`) and as a summary (D160); linkage disequilibrium (`--geno-r2`, `--hap-r2`, `--interchrom-geno-r2`, `--interchrom-hap-r2`, with `--ld-window*` and `--min-r2`), its sites held within `--max-mem` (D161); and relatedness of every pair of samples (`--relatedness` Yang's Ajk, `--relatedness2` Manichaikul's kinship) (D162). The usual `vcf` selection flags apply first (D158). Python has `vcf_popgen()`.
- **Differences from vcftools** are in `docs/differences.md` (the `vcf popgen` sections, `#popgen-sorted-input` to `#relatedness-rounding`): sorted input for windows, an error for a window step larger than the window, 64-bit counts, populations that vcftools mishandles refused, and LD and relatedness stopping where vcftools fills a file or reads masked genotypes as called.

**`gloci vcf compare A B`** — P3b. Genotype concordance only **[owner, 2026-11-11]**: replaces vcftools `--diff-site-discordance`, `--diff-indv-discordance`, `--diff-discordance-matrix` and `--diff-indv-map`. `--diff-site` and `--diff-indv` are `gloci intersect`/`subtract --match alleles` and the `#CHROM` line.

- Reports site overlap (shared, A-only and B-only), matched by alleles.
- Reports per-sample genotype concordance and a discordance matrix.

---

## 8. Architecture

### 8.1 Workspace **[my call]**

```
crates/
  gloci-core     Interval, Strand, ChromId interning, genome/order, overlap & footprint rules (Section 4)
  gloci-chains   the built-in Ensembl GRCh37/GRCh38 liftover chains, as data (D174, R13)
  gloci-chains-t2t  the built-in T2T-CHM13v2 liftover chains, as data (split off for crates.io's size limit)
  gloci-io       format detection; BGZF (parallel, libdeflate); BED/bedGraph/GFF parsers;
               VCF/BCF reader+writer; tbi/csi read+write; atomic writer; chain files
  gloci-engine   sorted k-way sweep, active set, unsorted-mode interval index,
               per-chromosome parallel driver, ordered output
  gloci-ops      one module per command, as library functions over record streams
  gloci-vcf      genotype decoding, selection pipeline, stats, popgen
  gloci        clap-based `gloci` binary, published as `gloci` so `cargo install gloci` works;
               no logic beyond argument wiring
  gloci-py     PyO3 extension module for the Python library (P2); no logic beyond wiring,
               like the CLI crate
xtask/         benchmark runner, differential-test runner, synthetic data generator (not shipped)
```

**Requirement from P1 on [you, enabling Section 8.5]:** every operation in `gloci-ops` takes its inputs as record streams (the `RecordReader` interface) and writes to a record sink. It never takes only file paths. That way the CLI, tests and the Python library can feed files, in-memory data or Arrow tables to the same code.

### 8.2 Data flow

```
 reader ─► BGZF decompress (pool) ─► line split + parse (pool, batches) ─► engine (sweep)
        ─► format (pool, ordered) ─► BGZF compress (pool, ordered) ─► atomic writer
```

- Records move in **batches** of thousands, not one at a time. Parsing doesn't allocate per record. Each record keeps a slice of its raw line for unchanged passthrough, and other fields are parsed lazily.
- Chromosome names are interned to `u32` IDs per run.
- **The sweep engine** uses the classic sorted "chromsweep":
  - It advances the B streams alongside A.
  - It keeps an active set of B ordered by end.
  - It evicts B whose `end + window <= a.start`.
  - Memory is limited by overlap depth, not file size.
  - Pathological depth, such as a whole-chromosome feature against everything, is documented, and the active set must stay correct when it happens.
- **One core intersection primitive (P1) [you, 2026-11-07, after reviewing bedder-rs].** The engine emits, for each A record, every overlapping B record tagged by its source file. Optionally it also emits the n nearest B within a maximum distance. `intersect`, `subtract`, `closest`, `map` and `coverage` are reporting layers over this one sweep, not separate algorithms.
  - Internally, what to report is described by a "piece" model, after bedder's `--a-piece`/`--b-piece`: whole A, the overlapping piece, the inverse piece (which is `subtract`), or none. The user-facing `--mode` flags in Section 7 map onto it.
- **Unsorted mode** loads one side into a cache-friendly interval index per chromosome. The index is `rust-lapper` (pure Rust, generic over 64-bit coordinates). It was chosen by benchmark in P2a over `coitrees`, superintervals and a cgranges-style tree. **[revised 2026-11-09, D61]** Intervals are stored end-inclusive (`[start, end+1)`), so the index is used only to generate candidates:
  - Queries are widened by one base on each side, and every candidate is re-tested with Section 4.3's exact rules. Points and bookended intervals therefore behave exactly as in sorted mode.
  - When the streamed side happens to be sorted, the index's sorted-query mode is used.

### 8.3 Parallelism

1. BGZF decompression and compression across threads, the biggest win on `.gz` inputs.
2. Parsing and formatting of batches across threads, with output reassembled in order.
3. **Per-chromosome (or per-region-chunk) parallelism** when every input is indexed, or is a plain seekable file that can be split. Chunks run independently and output is concatenated in chromosome order. This is P1 for indexed inputs where cheap, and P2 at the latest.
4. Streaming from stdin is limited to stages 1–2. The docs say so.

5. **CRAM decode pool (P3a) [D110].** CRAM containers are read in order and decoded on at most `-t` workers, with at most 2 × `-t` in flight, and released in order: the CRAM analogue of stage 1. Workers have 16 MiB stacks (noodles-cram's order-1 rANS and name tokeniser need them). The workers share one lease-based cache of reference sequences, each read from the FASTA once; with a `.crai` and `-t` above 1 a CRAM is read per chromosome, one chromosome per worker. Other pools stay deferred unless a benchmark misses a target.
6. **Waiting is sleeping.** Threads that wait for each other sleep on a mutex and condition variables (`crates/gloci-io/src/handoff.rs`): the hand-off between the selection stage of the `vcf` commands and its workers, and the workers of the BGZF reader and writer (`crates/gloci-io/src/bgzf.rs`), which replace `noodles_bgzf`'s multithreaded pair. A crossbeam channel's blocking calls and an idle rayon worker retry around `sched_yield`, and on Linux a yield on a CPU another thread wants gives up a whole scheduler slice (milliseconds); with more runnable threads than cores, which is a busy machine or `-t` above the core count, every block or batch then costs a slice. The compressed bytes are the same at every thread count, and a bad block fails with noodles' own error.

**CRAM limits [D109].** A reference is `--reference` only. Every `@SQ` in the `.fai` must match its CRAM `LN`, else an error before any record. gloci builds the `.crai` itself and never reaches a noodles code path that panics on a missing reference (the repository adapter returns an error, not `None`, and slice reference ids are checked first). Decoding is bounded: every codec's clear-text decoded size is checked against the block size, a block is capped at 1 GiB, and a read may overhang the end of its reference (`LN`) by at most max(1 Mb, 4 × the read length). The residual risk on crafted files with valid checksums is documented in `docs/differences.md#cram-limits`.

Output must be byte-identical for every thread count (Section 6.4).

### 8.4 Libraries **[my call; the builder may swap any of these with evidence]**

| Need | Choice | Why |
|---|---|---|
| CLI | `clap` (derive) + `clap_complete` + `clap_mangen` | completions and man pages for free |
| VCF/BCF/tabix/CSI/FASTA/BAM | `noodles` | pure Rust, spec-faithful; keeps the static build simple. Hot paths (BED parse, VCF fixed columns, GT decode) may get hand-written fast parsers if profiling says so, verified against noodles |
| Deflate | `libdeflater` | much faster than zlib for BGZF blocks |
| Threads | a mutex and condition variables (`gloci-io`'s `Handoff` and `Pool`) for the BGZF pools and the `vcf` selection stage; `crossbeam-channel` and `rayon` elsewhere (amended 2026-11-14) | bounded queues give back-pressure and bounded memory; waiters that sleep instead of yielding the CPU (item 6 of Section 8.3) |
| Bytes and numbers | `memchr`, `atoi`/hand-rolled, `itoa`, `ryu` | |
| Allocator | `mimalloc` | Built with `no_thp` **[owner request, 2026-11-16]**: by default mimalloc asks Linux for transparent huge pages for its whole arena, so every 2 MiB it touches is resident; without them peak memory falls by half or more (the merge of 1.67 million exons 26 → 11 MiB, an 8-thread BGZF intersect 258 → 119 MiB, an 8-thread `vcf filter` of the 1000 Genomes chr20 BCF 83 → 32 MiB) at the same speed (median of 11 interleaved `sort` runs 0.34 s either way). They also cost stalls while the kernel compacts memory to find 2 MiB pages (`sort` once spent 0.76 s of system time where it needs 0.07). The binary also sets mimalloc's purge delay to 10 ms (`crates/gloci/src/main.rs`; `MIMALLOC_PURGE_DELAY` in the environment still wins): at mimalloc's 1000 ms, longer than most runs, every buffer freed during a run stayed resident to its end (8-thread BGZF intersect 118 → 80 MiB, a CRAM `flagstat` 20 → 17 MiB, at the same speed; 0 ms made the CRAM decoder's frees and allocations 57 % slower). |
| Errors | `thiserror` (libs), `anyhow` (CLI) | |
| bigWig, bigBed | bigtools 0.5.8 (read, write) on a tokio runtime | MIT, an independent implementation of the published format; its one MPL-2.0 dependency, `index_list`, is accepted unmodified (owner, 2026-11-17; D171) |
| Testing | `proptest`, `cargo-fuzz`, `insta` (snapshots), `criterion` + `hyperfine` | |

- **No htslib.** A pure-Rust plus vendored-C build gives one static `musl` binary.
- Rust edition 2024, built with the toolchain pinned in `rust-toolchain.toml` (1.100.0) instead of the newest stable, so a new Rust release does not change a build by itself **[owner, 2026-11-14]**. The minimum supported Rust (`rust-version`, 1.91) is checked by its own CI job, `msrv`, with `cargo +1.91 check --workspace --locked`.
- Release profile: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`.
- Portable x86-64 builds, plus an optional `x86-64-v3` build.

### 8.5 Python library (P2) **[you; design my call]**

- **Packaging.** `pip install gloci` installs prebuilt wheels built with maturin: abi3, Linux x86_64 and aarch64 (manylinux), and macOS arm64. There are no runtime dependencies beyond Python. The bioconda package for the library is the separate `py-gloci` recipe (the `gloci` recipe is the CLI). The name was free on PyPI on 2026-11-06.
- **Terminal methods.** `save()`, `to_bytes()` (the output as `bytes`), `count()`, iteration and the Arrow/polars/pandas converters run the chain; `cli()` returns the equivalent `gloci` command line without running it.
- **In-process, not a CLI wrapper.** Calls go straight into `gloci-ops` through PyO3 and release the GIL during work. There is no subprocess and no temp-file shuttling, unlike pybedtools. Results and semantics are identical to the CLI because it is the same code.
- **Fluent, lazy pipelines.** Each step is recorded, and the chain runs as one streaming pass when a terminal method is called. Keyword arguments use the CLI flag names with `-` replaced by `_`.

  ```python
  import gloci
  (gloci.open("calls.vcf.gz")                 # path, file object, DataFrame/Arrow table, or iterable
        .intersect(b="targets.bed")           # same semantics and defaults as `gloci intersect`
        .merge(max_gap=50)
        .save("on_target.bed.gz", write_index=True))
  ```

  - **Terminal methods:** `.save(path, ...)`, `.to_arrow()`, `.to_polars()`, `.to_pandas()`, `.count()`, and iteration over records (the slow, convenient path).
- **Data frames.** Conversion goes through Apache Arrow, with zero-copy hand-off to `pyarrow` and `polars`. `pandas` goes through `pyarrow`. Coordinates follow Section 4 (0-based half-open), and the columns are named `chrom`, `start`, `end`, then format-specific columns.
- **Errors are typed exceptions** carrying the same information as the CLI messages:
  - a base `gloci.GlociError`;
  - `gloci.UsageError` (also a `ValueError`);
  - `gloci.DataError`, with `.path` and `.line`;
  - `gloci.UnsortedError(DataError)`, with both offending records and the hint.
- **Stability.** The Python API is versioned with semantic versioning. Rust crate APIs stay internal (Section 2).
- **Testing.** Python tests (pytest) run each fluent pipeline against the equivalent CLI command and assert identical output. That reuses the bedtools golden files (Section 9.3).

---

## 9. Correctness, performance, testing

### 9.1 Performance targets **[my call; firm numbers set after the Phase 0 baseline]**

Reference machine: 8-core x86_64, NVMe, 32 GB RAM. Baselines: bedtools 2.31.x, bedops 2.4.4x, vcftools 0.1.16.

| Area | Target |
|---|---|
| Sorted interval ops, plain BED, `-t 1` | ≥ 2× faster than `bedtools -sorted`; never slower than bedops |
| Same, bgzipped inputs, `-t 8` | ≥ 5× faster than bedtools |
| `vcf filter`, 1000 Genomes 30x chr20 (≈2.5k samples) | ≥ 5× faster than vcftools `--recode` |
| `vcf stats` with 5 stats | ≥ 10× faster than 5 separate vcftools runs |
| `sort` | faster than `sort-bed`; within 1.5× of GNU `sort --parallel`; obeys `--max-mem` |
| Memory, streaming ops on sorted input | flat as input grows; ≤ ~200 MB RSS excluding active-set depth |

### 9.2 Benchmark suite (lives in `xtask`, run in CI on every release and nightly)

- **Public data:**
  - 1000 Genomes 30x chr20 VCF/BCF.
  - gnomAD v4 sites chr22, which has heavy INFO.
  - GENCODE GTF exons.
  - ENCODE cCREs and peaks.
  - A large bedGraph signal track.
- **Synthetic generator:** produces 10M–1B intervals with controlled length distributions, overlap depth and chromosome counts, including pathological cases (deep nesting, a whole-chromosome feature, many points).
- **Public workload:** the coitrees benchmark set (A ≈ 2.76M intervals against B ≈ 62M, coverage counts), a realistic large case next to the synthetic data.
- **Baselines.** The binding comparisons (Section 2) are against the tools gloci replaces: bedtools, bedops and vcftools. The rest are yardsticks, so gloci's numbers stand up next to current alternatives. All are pinned in the benchmark container. **[you, 2026-11-07]**

| Area | Phase | Baselines |
|---|---|---|
| Interval operations (sort, intersect, merge, subtract, closest, map, coverage) | P0–P1 | bedtools, bedops; yardsticks gia, bedder, granges, bedtk, grit-genomics |
| `vcf filter` | P1 | vcftools `--recode`; yardsticks bcftools `view`, vcfkit, vcfexpress |
| `vcf stats` (several statistics in one pass) | P2 | vcftools, one run per statistic; yardsticks bcftools `stats`, PLINK 2, VariantFlow |
| Python library | P2 | pybedtools; yardsticks polars-bio, PyRanges v1, bioframe |
| BAM/CRAM `coverage`/`genomecov` | P3 | bedtools; yardstick mosdepth |
| `depth --summary`, `depth --by` | P3c | samtools coverage; mosdepth (binding) |
| Unsorted-mode interval index (internal) | P1 | coitrees (preferred) against superintervals, a cgranges-style tree and rust-lapper |

- **External datasets:** the datasets of the published interval-query benchmark (bioRxiv 2025.02.19.639167, which found AIList, bedtk and granges consistently fastest). This lets gloci's results be compared with published numbers.
- Results are stored per release so regressions show up.

### 9.3 Testing

1. **Semantics tables:** every table and rule in Section 4 has a table-driven unit test.
2. **Property tests:** every two-input and n-ary operation is compared with a **naive O(n·m) reference implementation** on random inputs. Inputs include points, bookended features, duplicates, nesting, huge intervals, empty chromosomes and strand mixes.
3. **Differential tests:**
   - Pinned bedtools, bedops and vcftools run in a container on a test corpus.
   - Every difference must be either fixed as a bug or listed in `docs/differences.md` with the reason. Examples of listed differences: default `--vcf-span allele`, MAF on multi-allelic sites, and INFO kept by default.
   - A difference that comes from the platform of the reference tool (its C library or CPU, such as vcftools rounding sums differently on aarch64) is pinned with `expected_difference_arch` in the case, naming the architectures where it shows; on any other the outputs must be equal.
4. **Fuzzing:** BED, VCF, BCF, GFF and BGZF readers, plus region-string parsing.
5. **Thread determinism:** every command runs with `-t 1`, `-t 2` and `-t 8`, and the outputs must be byte-identical. The behaviour tests (matrices of every command over its input formats, thread counts, failures and standard input) are one test binary, `crates/gloci/tests/behaviour/main.rs`, with a module for each matrix and for its fixtures.
6. **Large-file tests (nightly):**
   - Generated multi-GB inputs, with an RSS ceiling assertion for streaming commands.
   - `sort` on an input 4× larger than `--max-mem`.
7. **Failure behaviour:** unsorted input, malformed lines, truncated BGZF, a missing index, chromosome-name mismatch, and broken pipe. Each must give the documented message or exit code, and no partial output files.
8. **CI platforms:** Linux x86_64 and aarch64, macOS aarch64.

---

## 10. Phases **[my call — you didn't name the most-used commands, so this follows typical lab usage]**

| Phase | Contents | Done when |
|---|---|---|
| **P0 Foundations** | workspace; Section 4 core types and rules; BED/VCF/BCF/GFF I/O; parallel BGZF; sort validation; atomic writer; built-in genomes (`-g hg38`); `sort`, `index`, `view`, `check`, `genomes`; benchmark + differential + property-test harnesses | `sort` beats `sort-bed`; harnesses run in CI |
| **P1 Core replacement** | `intersect`, `subtract`, `closest`, `map` (core ops), `coverage`, `merge`, `complement`, `slop`, `flank`, `windows`, `convert` (VCF/GFF/BED12), `vcf filter` (no expressions); unsorted mode; per-chromosome parallelism for indexed inputs; every op callable on in-memory record streams (Section 8.1); tabix-style indexing of any bgzipped TSV | the team's most common bedtools/bedops/vcftools pipeline steps can be rewritten on gloci; P1 performance targets met |
| **P2 Breadth** | `partition`, `cluster`, `spacing`, `shift`, `groupby`, `genomecov` (BED), full `map` ops, `vcf stats`, `vcf extract`, `vcf concat`, filter expressions, `--chrom-map`, WIG input; bioconda package; **Python library** (`pip install gloci`, Section 8.5) | vcftools no longer needed for QC/filtering |
| **P3 Extended** | BAM/CRAM input (intersect, coverage, genomecov, convert) and indexing (`.bai`/`.csi`, `.crai`); FASTA indexing (`.fai`/`.gzi`); `getfasta`, `maskfasta`, `nuc` alignment `sort` (coordinate, `-n`, `-N`), `merge` and `cat` **[owner, 2026-11-10]**; `stats jaccard/fisher/reldist`; `sample`/`random`/`shuffle`; BEDPE pair ops; `vcf compare` (genotype concordance only) **[owner, 2026-11-11]**; `vcf popgen` (left out **[owner, 2026-11-11]**, then built on the owner's request of 2026-11-14) | all three tools can be uninstalled |
| **P3c samtools essentials [owner, 2026-11-10]** | on top of P3a's BAM/CRAM/FASTA support: `idxstats`, `flagstat`, `quickcheck`, `reheader`, `faidx` region extraction, `dict`, `view -c`/`-L`, `head`, `samples`, `flags`, `depth`, `bedcov`. Out of scope: `markdup`, `fixmate`, `calmd`, `mpileup`, `consensus`, `stats`, `ampliconclip` | samtools is needed only for what gloci does not cover |

**Order [revised 2026-11-10, owner]:** after P2a (the P2 command breadth), the build order is:
P2b Task 1 only (MIT licence; genome tables moved inside `gloci-io` so `cargo install` and source packages work),
then **P3a** (BAM/CRAM, FASTA, `getfasta`/`maskfasta`/`nuc`, BAM in `coverage`/`genomecov`/`intersect`/`convert`, alignment `sort`/`merge`/`cat`),
then **P3c** (samtools essentials, which build directly on P3a),
then **P3b** (`stats`, `sample`/`random`/`shuffle`, BEDPE, `vcf compare`; `vcf popgen` came later),
then the rest of P2b (Python library, packaging, mdBook, paper evidence suite) over the complete command set.
Nothing is published to PyPI, crates.io, GHCR or bioconda until the owner enables it.

**P3b scope and amendments [owner, 2026-11-11]:**

- **Note.** `vcf popgen` was left out of P3b, and built later on the owner's request of 2026-11-14. `vcf compare` is narrowed to genotype concordance: `--site-discordance`, `--sample-discordance`, `--discordance-matrix` and `--sample-map`; site and sample set differences are `intersect`/`subtract --match alleles` and the headers. P3b also covers the owner's VCF additions, `vcf split` (bcftools `norm -m-`), `vcf norm -f` (bcftools `norm -f`), `vcf annotate` and `vcf merge` (bcftools `annotate`, `merge`).
- **D157 amendment (parallel selection).** The batch cap (at most max(1, 524,288 / samples) records) applies at `-t 1` too. The earliest failing record wins at every `-t`. At `-t` above 1 a read-ahead warning can precede an earlier record's error, because workers run ahead of the sink.
- **D163 amendment (`vcf compare`).** The discordance walk follows vcftools. The matrix scope in `--match alleles` is the trimmed allele (`AT>A` and `ATG>AG` count together); in `--match position` it is the same REF and ALT text. `--thin` thins each file on its own, as `vcf filter --thin` does, whatever the two headers list and with `-g` **[owner, 2026-11-15, lifting the refusal this amendment first made]**: the selection checks only that each chromosome is one block. Threads are at most max(`-t`, 3): A and B are read on threads of their own beside the walking thread.
- **Memory bounds.** `shuffle` holds about 52 bytes per allowed interval; `random` holds one batch of 262,144 rows; `pair-intersect` and `pair-pair` hold B at about 5 to 6 times its text size.
- **Randomness.** `--seed` reproduces gloci's own output at any `-t`, never bedtools'; without it a fresh seed is printed as a warning.

**Distribution:**

- Static Linux binaries (x86_64 and aarch64 musl) and macOS binaries on releases from P1.
- `cargo install` and a container image.
- Bioconda from P2.
- Python wheels on PyPI (`pip install gloci`) from P2.

**Docs:**

- An mdBook with a semantics chapter (Section 4), the migration guide (Section 11) and the differences list.
- A `SKILL.md` agent usage guide, generated from the CLI help and the migration guide, so AI coding agents use gloci correctly (P1). This follows bedder-rs. **[you, 2026-11-07]**

---

## 11. Migration guide (old → new)

The user-facing version of this guide lives in `docs/migration.md` (later part of the mdBook). It starts with the behaviour changes below, because these make pipelines fail or give different output even after the flags are rewritten.

### Behaviour changes to know about

**Chromosome order.** Without `-g` or a VCF/BCF header, gloci sorts and expects chromosomes in lexicographic byte order (`chr1, chr10, chr11 … chr2`), as bedtools, bedops `sort-bed` and `sort -k1,1` do (`LC_ALL=C`). `--chrom-order karyotype` selects `chr1, chr2 … chr22, chrX, chrY, chrM` instead. In practice:

- Files sorted by the old tools work as they are, and `gloci sort` writes the same bytes as `sort-bed` and `bedtools sort`.
- A file in karyotype order (the output of `gloci convert` from a VCF or BAM, which keeps the header's order, or of `--chrom-order karyotype`) is reported as **not sorted** by the streaming commands. The error names both records. When the pair is in karyotype order, the hint says so and gives the fixes: pass `-g GENOME` or `--chrom-order karyotype`, or re-sort with `gloci sort FILE`. The other direction is hinted the same way for runs with `--chrom-order karyotype`.
- `-g genome.fai` and VCF/BCF header order still take precedence over the default.
- Mixing the two orders within one pipeline is reported as an error, never as a silent wrong answer.

**BAM, CRAM and FASTA.** Interval commands skip unmapped reads (bedtools `intersect -v` writes them); `view`, `sort`, `sort -m` and `cat` keep every record. The read filters default to nothing. CRAM needs `--reference FASTA` (bedtools `--cram-ref`), and gloci never writes a CRAM. `genomecov` of a BAM counts deleted bases, so bedtools' `-ibam -split` is `--split-deletions`. gloci never builds an index implicitly (bedtools builds a `.fai`). `gloci merge` is the interval merge; `samtools merge` is `gloci sort -m`. See `docs/migration.md`.

### bedtools

| bedtools | gloci |
|---|---|
| `intersect -u` | `intersect` (default `--mode any`) |
| `intersect -v` | `intersect --mode none` |
| `intersect` (no mode flag) | `intersect --mode clip` |
| `intersect -wa -wb` | `intersect --mode pairs` |
| `intersect -wo` | `intersect --mode pairs --overlap-bases` |
| `intersect -wao` / `-loj` | `intersect --mode pairs --keep-unmatched [--overlap-bases]` |
| `intersect -c` | `intersect --mode count` |
| `-f -F -r -e -s -S -split` | same short flags |
| `-sorted` | default (sortedness is checked); unsorted input needs `--unsorted` |
| `window -w/-l/-r/-sw` | `intersect -w / --left / --right / --strand-window` |
| `subtract`, `-A`, `-N` | `subtract`, `--whole`, `--sum-overlaps` |
| `merge -d -s -c -o` | `merge --max-gap -s -c COL:OP` |
| `complement -g` | `complement -g` |
| `closest -d -D -k -t -io -iu -id -mdb` | `closest --distance --signed -k --ties --ignore-overlaps --ignore-upstream --ignore-downstream --per-file` |
| `map -c -o` | `map -c COL:OP` |
| `coverage [-hist -d -mean]`, `annotate` | `coverage [--hist --per-base --mean]`, `coverage -b f1 f2 ...` |
| `genomecov` | `genomecov` |
| `slop`, `flank`, `shift` | same names |
| `makewindows` | `windows` |
| `multiinter`, `unionbedg` | `partition`, `partition --bedgraph` |
| `cluster`, `spacing`, `groupby`, `sort` | same names |
| `bed12tobed6`, `bamtobed [-bed12 -split -splitD -ed -tag -cigar -color]` | `convert --blocks`, `convert [--bed12 --blocks --split-deletions --score-tag --cigar --color]` (P3a) |
| `getfasta -fi`, `maskfasta -fi`, `nuc -fi` | same names, with `--fasta` for `-fi` (P3a) |
| `intersect -abam`, `coverage -b x.bam`, `genomecov -ibam`, `--cram-ref` | `intersect -a x.bam`, `coverage -b x.bam`, `genomecov x.bam`, `--reference` (P3a) |
| `jaccard`, `fisher`, `reldist` | `stats jaccard / fisher / reldist` (P3) |
| `shuffle`, `random`, `sample` | same names (P3) |
| `pairtobed`, `pairtopair` | `pair-intersect`, `pair-pair` (P3) |
| `igv -i F [-path D -sess S -sort K -clps -name -slop N -img T]` | `gloci igv F [--path D --session S --sort K --collapse --name --slop N --image T]` (D169 **[owner, 2026-11-15]**) |
| `links -i F [-base U -org O -db B]` | `gloci links F [--base U --org O --db B]` (D170 **[owner, 2026-11-15]**) |
| `bamtofastq -i FILE -fq R1.fq [-fq2 R2.fq] [-tags]` | `convert FILE --fastq [-o R1.fq] [--fastq2 R2.fq] [--tags]` (D166 **[owner, 2026-11-15]**) |
| `bedtobam [-mapq -bed12 -ubam]`, `bedpetobam [-mapq -ubam]` | `convert --bam -g G [--mapq --bed12 --compression-level 0]`, `convert --format bedpe --bam -g G` (D167 **[owner, 2026-11-15]**) |
| `multicov -bams A.bam B.bam -bed F` | `bedcov F A.bam B.bam --read-count` (D165 **[owner, 2026-11-15]**) |
| `tag -i IN.bam -files F… [-labels L…] [-tag XX] [-names -scores -intervals -s -S -f F]` | `gloci tag IN.bam -b F… [--labels L,… --tag XX --names --scores --intervals -s -S -f F]` (D168 **[owner, 2026-11-15]**) |
| `expand`, `split`, `summary`, `overlap` | not planned (see non-goals; `overlap` → `--mode pairs --overlap-bases`) |

### samtools **[P3a, owner addition]**

| samtools | gloci |
|---|---|
| `sort [-n -N -m -@ -T -o --write-index]` | `sort [-n -N --max-mem -t --tmp-dir -o --write-index]` (`-m` is per thread in samtools and a total in gloci; the output follows the input's format, where samtools writes BAM; temporary disk is about the input size) |
| `merge [-n -N -c -p]` | `sort -m [-n -N]` (IDs are combined by default; a conflicting `@RG` and differing `@SQ` lists are errors) |
| `cat` | `cat` (BAM only; `@SQ` lists must match) |
| `index` | `index` (`.bai`, `.csi` or `.crai`) |
| `faidx REF.fa` | `index REF.fa` |
| `faidx REF.fa REGION...`, `-r`, `-c`, `-i`, `-n`, `--mark-strand` | `faidx REF REGION... [--regions-file --skip-invalid --reverse-complement --line-width --mark-strand]` (needs the index first) |
| `dict`, `flags`, `idxstats`, `flagstat`, `head`, `samples`, `reheader` | the same names |
| `quickcheck [-u -v]` | `check --quick [--unaligned --list-failed]` |
| `view -c`, `view -L BED` | `view --count`, `view --regions-file BED` |
| `depth [-a -aa -J -q -Q -l -g -G -s -r -b -H]`, `bedcov [-Q -g -G -j -d -c -H]` | `depth`, `bedcov` (long flags, see the migration guide); bedtools `multicov` is `bedcov --read-count` |
| `coverage [-r -q -Q -l --rf --ff]` | `depth --summary [--region --min-mapq --min-baseq --min-read-length --require-flags --exclude-flags]` (D173 **[owner, 2026-11-17]**) |
| mosdepth `[-n] [-b BED\|SIZE] [-T 1,10,20] [-Q -F] PREFIX FILE` (a binding comparison) | `depth --out-prefix PREFIX [--by BED\|SIZE] [--thresholds 1,10,20] [--min-mapq --exclude-flags] FILE` (D173 **[owner, 2026-11-17]**) |

### bedops

| bedops | gloci |
|---|---|
| `bedops -u` / `--everything` | `sort -m` |
| `bedops -m` | `merge` |
| `bedops -c` | `complement` (no `-g`) |
| `bedops -d A B...` | `merge A \| subtract -a - -b B...` |
| `bedops -e [N%]` / `-n` | `intersect [-f 0.N]` / `intersect --mode none` |
| `bedops -i` / `-s` / `-p` | `partition` (see 7.2) |
| `bedops --chop N` | `merge ... \| windows -i - -w N` |
| `--range N` | `-w N` (relations) or `slop` |
| `bedmap --echo --count/--sum/--mean/...` | `map -c COL:count/sum/mean/...` |
| `bedmap --echo-map`, `--echo-map-id` | `map -c 0:records`, `map -c 4:collapse` |
| `bedmap --bases`, `--fraction-*` | `map -c 0:bases`, `-f/-F/-e` |
| `closest-features` | `closest --each-side [--distance]` |
| `sort-bed [--max-mem --tmpdir --check-sort --unique --duplicates]` | `sort [--max-mem --tmp-dir --check --unique --duplicates]`; `sort-bed`'s chromosome order is the default |
| `bedextract` | `view FILE REGION` |
| `convert2bed`, `vcf2bed [--snvs/--insertions/--deletions]`, `gff2bed`, `gtf2bed` | `convert [--types snv/ins/del]` |
| `starch`, `unstarch`, `starchcat` | not planned |

### vcftools

| vcftools | gloci |
|---|---|
| `--vcf/--gzvcf/--bcf` | input is auto-detected |
| `--recode`, `--recode-bcf` | `vcf filter ... -o out.vcf.gz` / `out.bcf` |
| `--recode-INFO-all` | default (INFO is kept); prune with `--keep-info/--drop-info` |
| `--chr/--not-chr/--from-bp/--to-bp` | `--chrom/--exclude-chrom/--region` |
| `--bed/--exclude-bed` | `--regions-file/--exclude-regions-file` |
| `--positions/--exclude-positions`, `--snp/--exclude` | `--positions/--exclude-positions`, `--ids/--exclude-ids` |
| `--keep-only-indels/--remove-indels` | `--types indel` / `--exclude-types indel` |
| `--remove-filtered-all`, `--keep-filtered/--remove-filtered` | `--pass-only`, `--keep-filter/--remove-filter` |
| `--keep-INFO/--remove-INFO` | `--require-info/--exclude-info` |
| `--maf/--max-maf/--mac/--max-mac/--min-alleles/--max-alleles/--minQ` | `--min-maf/--max-maf/--min-mac/--max-mac/--min-alleles/--max-alleles/--min-qual` |
| `--min-meanDP/--max-meanDP`, `--hwe`, `--thin`, `--phased` | `--min-mean-dp/--max-mean-dp`, `--hwe`, `--thin`, `--phased-only` |
| **`--max-missing X`** (X = minimum fraction *called*) | **`--max-missing (1−X)`** (fraction *missing*). Inverted on purpose; vcftools' meaning is a known trap. |
| `--indv/--keep/--remove-indv/--remove` | `--samples/--samples-file/--exclude-samples/--exclude-samples-file` |
| `--minGQ/--minDP/--maxDP`, `--remove-filtered-geno-all` | `--min-gq/--min-dp/--max-dp`, `--mask-filtered-genotypes` |
| `--freq`, `--depth`, `--missing-*`, `--het`, `--hardy`, `--TsTv*`, `--FILTER-summary`, `--SNPdensity`, `--singletons`, `--hist-indel-len` … | `vcf stats` flags, all in one run (P2) |
| `--get-INFO`, `--extract-FORMAT-info`, `--012` | `vcf extract` (P2) |
| `--site-pi`, `--window-pi`, `--TajimaD`, `--weir-fst-pop`, `--*-r2`, `--relatedness*` | `vcf popgen` |
| `--diff-site-discordance`, `--diff-indv-discordance`, `--diff-discordance-matrix`, `--diff-indv-map` | `vcf compare` (P3b); `--diff-site` and `--diff-indv` are `intersect`/`subtract --match alleles` and the `#CHROM` line |
| `vcf-concat`, `vcf-sort`, `vcf-isec`, `vcf-subset`, `vcf-stats` | `vcf concat` (P2), `sort`, `intersect --match alleles`, `vcf filter --samples`, `vcf stats` |
| `vcf-merge` | `vcf merge` (bcftools `merge`; an owner addition to P3b) |
| `--plink*`, `--IMPUTE`, `--BEAGLE-*`, `--ldhat*` | not planned |

---

## 12. Left to the builder

- Batch sizes, channel depths, and how threads are split between stages.
- The interval index used in unsorted mode, chosen by benchmark.
- Where hand-written parsers replace `noodles`, driven by profiling and verified by the differential and fuzz tests.
- How transforms that break sort order are reported or repaired (Section 7.3).
- Exact TSV column names for new outputs that have no vcftools equivalent.
- The concrete grammar of the P2 expression language, modelled on bcftools `-i/-e`.

## 13. Open questions (non-blocking)

- **License, and public open source or internal.** The licence is MIT (confirmed). The copyright holder (the LICENSE names Duarte Molha) is still the owner's call and should be settled before the first public release.
- **Which commands the team uses most.** The answer may reorder P1 and P2, but it doesn't change the design.
