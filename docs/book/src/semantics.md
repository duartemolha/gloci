# Semantics

This page is the contract. Every rule here has a table-driven test, and a rule never changes without this
page changing first. Where gloci differs from bedtools, bedops or vcftools on purpose, see
[Differences](differences.md).

## Coordinates

- Internally every feature is a **0-based, half-open** interval `[start, end)` on a named chromosome.
- Positions are **64-bit unsigned**, so chromosomes longer than 2³² bp work.
- Conversion on input:

| Format | On disk | Internal |
|---|---|---|
| BED, bedGraph | 0-based half-open | unchanged |
| VCF, BCF | 1-based `POS` | see 4.2 |
| GFF3, GTF | 1-based closed `[s, e]` | `[s-1, e)` |
| Region strings `chr:s-e` | 1-based closed (samtools convention) | `[s-1, e)` |

- A **point** is a zero-length interval `[p, p)`. It sits between 0-based bases `p-1` and `p`. Points are
  valid input: BED allows `start == end`, and insertions are points.
- `start > end` is an error. Negative or non-integer coordinates are errors.

## Variant footprints (VCF/BCF to interval) {#variant-footprints}

Selected with the global `--vcf-span <allele|ref>`. **The default is `allele`.**

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
- **Breakend mates** are ignored.
- **gloci does not left-align.** Inside repeats the footprint depends on how the caller wrote the variant. Normalise the input first (for example with `gloci vcf norm`, or `bcftools norm`).

## Overlap {#overlap}

- Two non-empty intervals overlap when `a.start < b.end && b.start < a.end`. They must share at least one base. **Bookended intervals do not overlap.**
- **A point `p` vs a non-empty interval `[s,e)`:** they overlap when `s < p < e`, meaning the point is strictly inside with interval bases on both sides. `--inclusive-points` relaxes this to `s <= p <= e`.
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

- **Fraction thresholds:** `-f` is the minimum fraction of A covered, and `-F` is the minimum fraction of B covered. Both accept a decimal (`0.5`) or a percentage (`50%`). `-r` requires both, and `-e` requires either. `--min-bases N` sets an absolute minimum overlap. Every threshold counts as met when either side of the overlap is a point, because a point covers no bases; an insertion inside a feature is never dropped by a threshold.
- **Strand:** `-s` keeps only overlaps on the same strand, and `-S` only overlaps on opposite strands. Features with unknown strand (`.`, or any VCF record) never match under `-s` or `-S`, and gloci warns once.
- **BED12 blocks:** `--split` uses blocks (exons) instead of the full span.
- **Windows:** `-w N` pads B by N on both sides before testing; `--left` and `--right` set the padding asymmetrically, and `--strand-window` reads them relative to A's strand.

## Distance

- Overlapping features have distance **0**.
- Non-overlapping features have distance equal to the gap in bases **+ 1**, so bookended features have distance 1. This matches bedtools.
- Signed distances (`--signed ref|a|b`) are negative for upstream features, with the reference frame chosen by the flag.

## Chromosome order and sortedness {#chromosome-order}

Streaming commands need inputs sorted by **(chromosome, start)**. The chromosome order is the first of:

1. `--genome NAME|FILE` (or `-g`): a built-in genome such as `hg38`, or a `.fai`, a two-column `.genome`, a `.dict`, or a VCF/BCF/BAM/SAM/CRAM header.
2. The contig order declared in the header of the first input that has one (`##contig` lines of a VCF/BCF, `@SQ` lines of an alignment file).
3. **Karyotype order**: `chr1 ... chr22, chrX, chrY, chrM`, with a leading `chr` ignored and everything else in natural order. `--chrom-order lex` selects byte order instead, the same as `sort -k1,1`, bedops `sort-bed` and `bedtools sort`.

Sortedness is checked on the fly in every streaming command. At the first violation the command exits 1, naming the file, the line, both records and a hint such as `run: gloci sort FILE -g <genome>`. Two switches handle unsorted input, each with one meaning everywhere:

- `--sort` sorts the inputs first, within `--max-mem` (results come out sorted).
- `--unsorted` works with the input in the order it is: B is loaded into memory and A is streamed in its own order.

Chromosome names must match exactly. When two inputs share no names, gloci warns and, for the common `chr1` vs `1` case, suggests `--chrom-map FILE`.

## Records and passthrough

- Input lines are carried through **unchanged** wherever a command outputs input records: extra BED columns and the full VCF line are kept, and numbers are not reformatted.
- Commands that output input records unchanged write **A's format** (`intersect --mode any|none`, `subtract`, `vcf filter`). A VCF A gives VCF out with its header.
- Commands that **append columns** (`--mode pairs`, `count`, `map`, `closest`, `coverage`) write **TSV**. VCF in gives the VCF columns plus the appended ones, with a single `#`-prefixed column-name line and no meta-header.
- VCF output gets a `##glociCommand=<full command line>; Version=<x.y.z>` header line for provenance.
- `subtract` on a VCF A removes whole records, the same as `intersect --mode none` with the same criteria.
