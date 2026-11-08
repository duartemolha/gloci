# gloci — Design Spec

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
- Variant calling, left-alignment or normalization against a reference, annotation (VEP-like), imputation and phasing.
- Becoming a general replacement for bcftools or samtools. Merging samples across VCFs (`vcf-merge`) is not planned for v1.
- The bedops **Starch** format, and legacy genotype export formats (IMPUTE, BEAGLE, LDhat, PLINK ped/tped). **[my call]**
- bedtools' BAM-specific utilities: `bamtofastq`, `bedtobam`, `bedpetobam`, `tag`, `multicov`, `igv`, `links`. **[my call]**
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
gloci vcf stats cohort.vcf.gz --freq --missing-site --missing-indv --het --tstv --out-prefix qc/cohort
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
- **gloci does not left-align.** Inside repeats the footprint depends on how the caller wrote the variant. The docs recommend normalized input (e.g. `bcftools norm`).

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
  1. `--genome NAME|FILE` order, if given. NAME is a built-in genome such as `hg38` (Section 5.3). FILE can be `.fai`, a 2-column `.genome`, `.dict`, or a VCF/BCF/BAM header.
  2. Otherwise, the contig order declared in a VCF/BCF header among the inputs, if any.
  3. Otherwise, **karyotype order** **[you]**, which matches how people read chromosomes and how reference genomes (and therefore most VCFs) are laid out. Names are compared with any leading `chr` prefix ignored:
     1. names that are a number, optionally followed by letters (`1`, `2`, `22`, `2L`, `2A`), by number and then suffix;
     2. then `X`, `Y`, and `M`/`MT`, in that order;
     3. then everything else (`chr1_KI270706v1_random`, `chrUn_…`, `chrEBV`, `HLA-…`) in natural order, with digit runs compared as numbers;
     4. remaining ties (`1` vs `chr1`, `chr01` vs `chr1`) are broken by raw bytes, so the order is total.

     This gives `chr1 … chr22, chrX, chrY, chrM` for human, the GRCh38 layout. `--chrom-order lex` selects lexicographic byte order instead, the same as `sort -k1,1`, bedops `sort-bed` and `bedtools sort`. Names that don't follow these patterns, such as yeast's Roman numerals, need `-g`.
- **Validation:** sortedness is **checked on the fly** in every streaming command.
  - At the first violation, the command stops with exit code 1.
  - The error names the file, line number and both offending records, and gives a hint such as: `run: gloci sort FILE -g <genome>`.
  - This check is cheap and is never turned off.
- **`--unsorted`:** opts into an in-memory mode. The smaller or B-side input is loaded into an interval index and the other input is streamed in any order. Memory is then proportional to the loaded file, and the help text says so.
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
| BAM, CRAM | as intervals (alignment spans or blocks) | 3 |
| FASTA + `.fai` | for `getfasta`, `maskfasta`, `nuc` | 3 |
| PSL, RepeatMasker `.out` | convert only | 3 |

- **Format detection** uses content and magic bytes: the BGZF/gzip header, `BCF\2`, `##fileformat=VCF`, `##gff-version`, and BED column shape. File extensions are not used. `--format` / `--a-format` / `--b-format` override detection.
- **Compression:** gzip and BGZF are read transparently. Plain gzip is read single-threaded, because it can't be split. The docs recommend bgzip.
- **stdin:** `-` reads from stdin wherever a file is accepted.

### 5.2 Output

- The default output is stdout, or `-o PATH`.
- **Format and compression** come from the `-o` extension (`.bed`, `.bed.gz`, `.vcf`, `.vcf.gz`, `.bcf`, `.tsv`, `.tsv.gz`, `.bedgraph.gz`), or from `-O FMT` for stdout. `.gz` output is always **BGZF**, so it can be indexed.
- `--write-index` creates `.tbi`. It switches to `.csi` automatically for chromosomes over 512 Mbp and always uses `.csi` for BCF.
- **Atomic file output:** gloci writes to `PATH.tmp-XXXX` and renames on success. On failure the temp file is removed and no partial file is left behind. **[my call]**

---

### 5.3 Built-in genomes **[you]**

gloci ships chromosome names and sizes for common assemblies, so `-g hg38` works without a file. Passing a file still works.

- **How `-g` resolves.** A value containing `/`, or naming an existing file, is a file. Otherwise it is looked up as a built-in name (case-insensitive). An unknown name is an error listing the built-ins, with a "did you mean" suggestion.
- **Naming styles are separate entries**, because names must match exactly (Section 4.5):
  - UCSC ids use `chr`-prefixed names.
  - GRC/Ensembl ids use the Ensembl names (`1 … 22, X, Y, MT`).
  - Each entry is built from its own project's published tables. They are not renamed copies: hg19 `chrM` (16,571 bp) is not GRCh37 `MT` (16,569 bp).
- **Order.** Built-ins list chromosomes in the same karyotype order gloci uses without `-g` (Section 4.5), so `-g NAME` adds lengths and the exact chromosome set but never reorders. The exception is genomes with Roman-numeral chromosome names (yeast `sacCer3`/`R64`, worm `ce11`/`WBcel235`). These list their Roman-numeral chromosomes in numeric order (I, II, … XVI), followed by the rest.
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
- **`gloci genomes`** lists the built-ins: id, species, assembly, naming style, sequence count and total length. `gloci genomes show NAME` prints one as a 2-column `.genome` file.

## 6. CLI conventions

### 6.1 Global options (where relevant to a command)

| Option | Meaning |
|---|---|
| `-o, --output PATH` | output file (default stdout) |
| `-O, --output-format FMT` | `bed`, `bed.gz`, `vcf`, `vcf.gz`, `bcf`, `tsv`, `tsv.gz`, `bedgraph`, `bedgraph.gz` |
| `--write-index` | write `.tbi`/`.csi` next to the output |
| `--compression-level N` | BGZF level (default 6) |
| `-t, --threads N` | default: `available_parallelism()`, which honours cgroups and CPU affinity (so SLURM allocations); env `GLOCI_THREADS` |
| `-g, --genome NAME\|FILE` | chromosome sizes and order: a built-in genome name (`hg38`, `GRCh38`, `mm39` …, see `gloci genomes`) or a file (`.fai`, `.genome`, `.dict`, VCF/BCF/BAM header) |
| `--chrom-order karyotype\|lex` | chromosome order when neither `-g` nor a VCF/BCF header sets one (Section 4.5), default `karyotype` |
| `--vcf-span allele\|ref` | variant footprint rule (Section 4.2), default `allele` |
| `--unsorted` | allow unsorted input via the in-memory index (Section 4.5) |
| `--region chr:s-e` (repeatable), `--regions-file BED` | restrict work to regions; uses the index when present |
| `--chrom-map FILE` | rename chromosomes on input (Phase 2) |
| `--header` | carry BED/TSV header lines from A to the output |
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
  - `quantile:Q`, `tmean:LO-HI`, `wmean` (weighted by overlapping bases)
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
| `genomecov` | bedtools genomecov | P2 (BED), P3 (BAM) | `--bedgraph`, `--bga`, `--hist` (default), `--per-base`, `--scale X`, `--strand +/-`, `--split`, `--5p/--3p`, `--max-depth N` |
| `sample`, `random`, `shuffle` | bedtools | P3 | `--seed` required for reproducibility; `shuffle --exclude/--include BED` |

**Sort order of transform output:**

- Transforms that keep order intact, such as uniform `slop` or `flank` without `-s`/`--pct`, produce sorted output.
- Transforms that can break order say so in `--help` and warn when that actually happens.
- How to fix order cheaply (a bounded reorder buffer or a recommendation to run `sort`) is left to the builder.

### 7.4 Sorting, indexing, viewing, conversion

**`gloci sort FILES...`** — P0. Replaces bedops `sort-bed`, `vcf-sort`, and `sort -k1,1 -k2,2n`.

- **Order:** chromosome order (Section 4.5), then start, then end.
  - Ties keep input order (a stable sort).
  - `--full-key` breaks ties on the whole line instead.
- **Input:** BED, VCF/BCF (keeping the header; chromosome order comes from `##contig` unless `-g` is given) and GFF.
- **External merge sort:**
  - Memory is limited by `--max-mem` (default 2G), with chunks sorted in parallel and spill files in `--tmp-dir` (compressed with a fast codec).
  - Inputs larger than RAM must work.
- `-m, --merge-sorted`: k-way merge of already sorted inputs (bedops `-u` / `--everything`, `sort -m`).
- `--check`: checks order without writing anything. Exits 1 and prints the first out-of-order pair if unsorted.
- `--unique` / `--duplicates`: keep only unique lines, or only duplicated ones.
- Output can be BGZF with `--write-index`.

**`gloci index FILE`** — P0, extended in P1 and P3. Replaces `tabix`, `bcftools index`, `samtools index` and `samtools faidx` for the formats gloci reads. **[you]**

- `gloci index FILE` chooses the index type from the file's content, with no format flags. The one exception is arbitrary TSV (below).
- Index types by phase:

| Input | Index | Replaces | Phase |
|---|---|---|---|
| bgzipped VCF, BED, bedGraph, GFF/GTF | `.tbi`, or `.csi` when an end exceeds 2^29 | `tabix -p vcf\|bed\|gff` | P0 |
| BCF | `.csi` | `bcftools index` | P0 |
| any bgzipped tab-separated file | `.tbi`/`.csi`, using `-s/--seq-col`, `-b/--begin-col`, `-e/--end-col`, `--zero-based`, `--skip-lines N`, `--comment CHAR` | `tabix -s -b -e` | P1, so TSV outputs of `map`, `coverage` and `intersect --mode pairs` can be indexed |
| BAM | `.bai`, or `.csi` with `--csi` or for long references | `samtools index` | P3, with BAM input |
| CRAM | `.crai` | `samtools index` | P3, with CRAM input and reference handling |
| FASTA, bgzipped FASTA | `.fai`, plus `.gzi` for BGZF | `samtools faidx` | P3, with `getfasta`/`maskfasta`/`nuc` |

- Each new index type ships with differential tests showing that the reference tool (htslib, samtools) reads gloci's index and returns identical query results (Section 9.3).
- Indexing formats gloci reads is in scope; becoming a general samtools replacement is not (Section 2).

**`gloci view FILE [REGION...]`** — P0. Replaces bedops `bedextract`, `tabix` region queries, and vcftools `--recode` as a pure format conversion.

- Region queries use the index if one exists.
- On **uncompressed** BED without an index, gloci scans the file and checks its order, so an unsorted file is an error and never a wrong answer. With `--assume-sorted`, regions are found by binary search on the file instead, like bedextract. This is fast, but it checks only the chromosome block it reads, so it trusts the user that the file is sorted. Indexing (`gloci sort … -o x.bed.gz --write-index`) gives both speed and safety. **[you: correctness first; revised 2026-11-07]**
- Also handles format conversion (`-O vcf.gz|bcf|...`), `--header-only` and `--no-header`.

**`gloci genomes [show NAME]`** — P0. Lists the built-in genomes, or prints one as a `.genome` file (Section 5.3).

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
- WIG to bedGraph (P2). BAM/SAM to BED, PSL and RepeatMasker (P3).

### 7.5 Sequence (P3)

`getfasta`, `maskfasta` and `nuc` work as in bedtools, with strand-aware reverse-complement, `--split` for BED12, and soft or hard masking. They need FASTA plus `.fai`; bgzipped FASTA plus `.gzi` is also supported.

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

**`gloci vcf concat FILES...`** — P2. Replaces `vcf-concat`.

- Checks that samples, headers and order are compatible.
- `--naive` concatenates BGZF blocks without decompressing when headers are identical. This is much faster.

**`gloci vcf popgen IN`** — P3. Replaces `--site-pi`, `--window-pi`, `--TajimaD`, `--weir-fst-pop`/`--fst-window-*`, `--geno-r2`, `--hap-r2`, `--interchrom-*`, `--ld-window*`, `--relatedness` and `--relatedness2`.

- Multiple statistics are computed in one pass where possible.
- Formulas are documented, and each is cross-checked against vcftools.

**`gloci vcf compare A B`** — P3. Replaces vcftools `--diff*` and `vcf-compare`.

- Reports site overlap (shared, A-only and B-only), matched by alleles.
- Reports per-sample genotype concordance and a discordance matrix.

---

## 8. Architecture

### 8.1 Workspace **[my call]**

```
crates/
  gloci-core     Interval, Strand, ChromId interning, genome/order, overlap & footprint rules (Section 4)
  gloci-io       format detection; BGZF (parallel, libdeflate); BED/bedGraph/GFF parsers;
               VCF/BCF reader+writer; tbi/csi read+write; atomic writer
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
- **Unsorted mode** loads one side into a cache-friendly interval index per chromosome. The preferred index is `coitrees`: pure Rust, with SIMD on AVX2 and on ARM NEON. The builder may swap it with benchmark evidence (a cgranges-style implicit tree, `rust-lapper`, superintervals). Because coitrees intervals are end-inclusive integers with no zero-length intervals, it is used only to generate candidates:
  - Queries are widened by one base on each side, and every candidate is re-tested with Section 4.3's exact rules. Points and bookended intervals therefore behave exactly as in sorted mode.
  - Chromosomes with coordinates beyond the index's integer range fall back to a 64-bit implementation.
  - When the streamed side happens to be sorted, the index's sorted-query mode is used.

### 8.3 Parallelism

1. BGZF decompression and compression across threads, the biggest win on `.gz` inputs.
2. Parsing and formatting of batches across threads, with output reassembled in order.
3. **Per-chromosome (or per-region-chunk) parallelism** when every input is indexed, or is a plain seekable file that can be split. Chunks run independently and output is concatenated in chromosome order. This is P1 for indexed inputs where cheap, and P2 at the latest.
4. Streaming from stdin is limited to stages 1–2. The docs say so.

Output must be byte-identical for every thread count (Section 6.4).

### 8.4 Libraries **[my call; the builder may swap any of these with evidence]**

| Need | Choice | Why |
|---|---|---|
| CLI | `clap` (derive) + `clap_complete` + `clap_mangen` | completions and man pages for free |
| VCF/BCF/tabix/CSI/FASTA/BAM | `noodles` | pure Rust, spec-faithful; keeps the static build simple. Hot paths (BED parse, VCF fixed columns, GT decode) may get hand-written fast parsers if profiling says so, verified against noodles |
| Deflate | `libdeflater` | much faster than zlib for BGZF blocks |
| Threads | `crossbeam-channel` + a fixed worker pool (or `rayon` where it fits) | bounded queues give back-pressure and bounded memory |
| Bytes and numbers | `memchr`, `atoi`/hand-rolled, `itoa`, `ryu` | |
| Allocator | `mimalloc` | |
| Errors | `thiserror` (libs), `anyhow` (CLI) | |
| Testing | `proptest`, `cargo-fuzz`, `insta` (snapshots), `criterion` + `hyperfine` | |

- **No htslib.** A pure-Rust plus vendored-C build gives one static `musl` binary.
- Rust stable, edition 2024.
- Release profile: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`.
- Portable x86-64 builds, plus an optional `x86-64-v3` build.

### 8.5 Python library (P2) **[you; design my call]**

- **Packaging.** `pip install gloci` installs prebuilt wheels built with maturin: abi3, Linux x86_64 and aarch64 (manylinux), and macOS arm64. There are no runtime dependencies beyond Python. The same package goes to bioconda. The name was free on PyPI on 2026-11-06.
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
| Unsorted-mode interval index (internal) | P1 | coitrees (preferred) against superintervals, a cgranges-style tree and rust-lapper |

- **External datasets:** the datasets of the published interval-query benchmark (bioRxiv 2025.02.19.639167, which found AIList, bedtk and granges consistently fastest). This lets gloci's results be compared with published numbers.
- Results are stored per release so regressions show up.

### 9.3 Testing

1. **Semantics tables:** every table and rule in Section 4 has a table-driven unit test.
2. **Property tests:** every two-input and n-ary operation is compared with a **naive O(n·m) reference implementation** on random inputs. Inputs include points, bookended features, duplicates, nesting, huge intervals, empty chromosomes and strand mixes.
3. **Differential tests:**
   - Pinned bedtools, bedops and vcftools run in a container on a test corpus.
   - Every difference must be either fixed as a bug or listed in `docs/differences.md` with the reason. Examples of listed differences: default `--vcf-span allele`, default karyotype chromosome order, MAF on multi-allelic sites, and INFO kept by default.
4. **Fuzzing:** BED, VCF, BCF, GFF and BGZF readers, plus region-string parsing.
5. **Thread determinism:** every command runs with `-t 1`, `-t 2` and `-t 8`, and the outputs must be byte-identical.
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
| **P3 Extended** | BAM/CRAM input (intersect, coverage, genomecov, convert) and indexing (`.bai`/`.csi`, `.crai`); FASTA indexing (`.fai`/`.gzi`); `getfasta`, `maskfasta`, `nuc`; `stats jaccard/fisher/reldist`; `sample`/`random`/`shuffle`; BEDPE pair ops; `vcf popgen`; `vcf compare` | all three tools can be uninstalled |

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

**Chromosome order.** Without `-g` or a VCF/BCF header, gloci sorts and expects chromosomes in karyotype order (`chr1, chr2 … chr22, chrX, chrY, chrM`). bedtools, bedops `sort-bed` and `sort -k1,1` use lexicographic order (`chr1, chr10, chr11 …`). In practice:

- Files sorted by the old tools are reported as **not sorted** by gloci's streaming commands. The error names both records. When the pair is in lexicographic order, the hint says the file looks sorted by `sort`/`sort-bed`/`bedtools sort`, and gives both fixes: re-sort with `gloci sort FILE`, or pass `--chrom-order lex`.
- `gloci sort` output differs from `sort-bed`/`bedtools sort` output in chromosome order only. Within a chromosome, records are ordered the same way. `gloci sort --chrom-order lex` reproduces the old order.
- `-g genome.fai` and VCF/BCF header order still take precedence, so pipelines that already pass a genome file behave as before.
- Re-sort stored intermediate files once with `gloci sort`, or keep `--chrom-order lex` on every gloci command in a pipeline that must stay byte-compatible with old outputs. Mixing the two orders within one pipeline is reported as an error, never as a silent wrong answer.

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
| `bed12tobed6`, `bamtobed` | `convert --blocks`, `convert` (P3) |
| `getfasta`, `maskfasta`, `nuc` | same names (P3) |
| `jaccard`, `fisher`, `reldist` | `stats jaccard / fisher / reldist` (P3) |
| `shuffle`, `random`, `sample` | same names (P3) |
| `pairtobed`, `pairtopair` | `pair-intersect`, `pair-pair` (P3) |
| `multicov`, `tag`, `bamtofastq`, `bedtobam`, `bedpetobam`, `igv`, `links`, `expand`, `split`, `summary`, `overlap` | not planned (see non-goals; `overlap` → `--mode pairs --overlap-bases`) |

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
| `sort-bed [--max-mem --tmpdir --check-sort --unique --duplicates]` | `sort [--max-mem --tmp-dir --check --unique --duplicates]`; add `--chrom-order lex` to reproduce sort-bed's chromosome order |
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
| `--site-pi`, `--window-pi`, `--TajimaD`, `--weir-fst-pop`, `--*-r2`, `--relatedness*` | `vcf popgen` (P3) |
| `--diff*`, `vcf-compare` | `vcf compare` (P3) |
| `vcf-concat`, `vcf-sort`, `vcf-isec`, `vcf-subset`, `vcf-stats` | `vcf concat` (P2), `sort`, `intersect --match alleles`, `vcf filter --samples`, `vcf stats` |
| `vcf-merge`, `--plink*`, `--IMPUTE`, `--BEAGLE-*`, `--ldhat*` | not planned |

---

## 12. Left to the builder

- Batch sizes, channel depths, and how threads are split between stages.
- The interval index used in unsorted mode, chosen by benchmark.
- Where hand-written parsers replace `noodles`, driven by profiling and verified by the differential and fuzz tests.
- How transforms that break sort order are reported or repaired (Section 7.3).
- Exact TSV column names for new outputs that have no vcftools equivalent.
- The concrete grammar of the P2 expression language, modelled on bcftools `-i/-e`.

## 13. Open questions (non-blocking)

- **License, and public open source or internal.** This affects nothing in P0–P1 but should be settled before the first release.
- **Which commands the team uses most.** The answer may reorder P1 and P2, but it doesn't change the design.
