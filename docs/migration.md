# Migrating from bedtools, bedops and vcftools

gloci (pronounced GLOH-sigh, like "loci") is one tool for genomic interval and
variant work.

This guide is for people moving existing pipelines to gloci. It starts with the
behaviour changes that break or alter a pipeline even after the flags are
rewritten.

## Chromosome order (behaviour change)

Without `-g` and without a VCF/BCF header that lists contigs, gloci uses
**karyotype order**: `chr1, chr2 ... chr22, chrX, chrY, chrM`. `sort -k1,1`,
bedops `sort-bed` and `bedtools sort` use **lexicographic** order (for `sort`,
only with `LC_ALL=C`; other locales order names differently again):
`chr1, chr10, chr11 ... chr19, chr2, chr20 ...`.

Same records, `chr1, chr2, chr10, chrX`:

```
# sort-bed / bedtools sort / sort -k1,1 -k2,2n      # gloci sort
chr1    100  200                                    chr1    100  200
chr10   100  200                                    chr2    100  200
chr2    100  200                                    chr10   100  200
chrX    100  200                                    chrX    100  200
```

Within a chromosome the order is the same. Only the chromosome order differs.

Records that share a start are ordered by end by `gloci sort`. `sort -k1,1 -k2,2n -s`
leaves them in input order; without `-s`, GNU `sort` falls back to comparing the
whole line, so ties are not in input order either. `gloci sort` keeps input order
for records that tie on chromosome, start and end, and `--full-key` orders such
ties by the whole line (see [differences.md](differences.md#sort-tie-order)).

### What you will see

A file sorted by an old tool is rejected by gloci's streaming commands:

```
x.bed is not sorted: line 3 (chr2	0	3) comes after line 2 (chr10	1	2)
hint: run: gloci sort x.bed (this file looks sorted by sort -k1,1, sort-bed or bedtools sort; re-sort it with gloci sort, or pass --chrom-order lex; see docs/migration.md)
```

gloci reports this as an error. It never silently gives a wrong answer from
files in a mixed order.

### Two fixes

1. Re-sort once, then keep the files: `gloci sort x.bed > x.sorted.bed`.
2. Keep the old order everywhere: add `--chrom-order lex` to every gloci command
   in the pipeline. Use this when outputs must stay byte-compatible with old
   results. `gloci sort x.bed --chrom-order lex` writes the old order.

### What still takes precedence

`-g genome.fai` (or `.genome`, `.dict`, a BAM header) and the contig order in a
VCF/BCF header win over `--chrom-order`. A pipeline that already passes `-g`
behaves as before. `--chrom-order` only applies when neither is present.

Names that follow no common pattern, such as yeast's Roman numerals, need `-g`.

### Built-in genomes

bedtools needs a chrom-sizes file (`-g human.hg38.genome`). gloci has the common
assemblies built in, so `-g hg38` replaces it: `gloci sort x.bed -g hg38`. Files
still work, and a file named like a built-in wins (with a warning). `gloci genomes`
lists the built-ins (UCSC names such as `hg38` and `mm39`, Ensembl names such as
`GRCh38` and `GRCm39`) and `gloci genomes show hg38` prints one as a `.genome`
file. Yeast (`-g sacCer3`) comes in chromosome order, I to XVI.

### Worked example

Before:

```sh
LC_ALL=C sort -k1,1 -k2,2n peaks.bed > peaks.sorted.bed
LC_ALL=C sort -k1,1 -k2,2n genes.bed > genes.sorted.bed
bedtools intersect -sorted -a peaks.sorted.bed -b genes.sorted.bed -u > hits.bed
```

After, with karyotype order (`gloci sort` is available now; `gloci intersect`
arrives in a later release, so its lines show the intended form):

```sh
gloci sort peaks.bed > peaks.sorted.bed
gloci sort genes.bed > genes.sorted.bed
gloci intersect -a peaks.sorted.bed -b genes.sorted.bed > hits.bed
```

After, keeping the old order so `hits.bed` is unchanged (again, `intersect` is
the intended form):

```sh
gloci sort peaks.bed --chrom-order lex > peaks.sorted.bed
gloci sort genes.bed --chrom-order lex > genes.sorted.bed
gloci intersect -a peaks.sorted.bed -b genes.sorted.bed --chrom-order lex > hits.bed
```

## Command mapping

These are the commands P0 ships. Every row was checked against `gloci --help`.

| Old command | gloci |
|---|---|
| `sort-bed FILE`, `bedtools sort -i FILE`, `LC_ALL=C sort -k1,1 -k2,2n FILE` | `gloci sort FILE --chrom-order lex` (drop `--chrom-order lex` to get karyotype order) |
| `bedtools sort -i FILE -g GENOME` | `gloci sort FILE -g GENOME` (`-g` also takes a built-in such as `hg38`) |
| `vcf-sort FILE`, `bcftools sort FILE` | `gloci sort FILE` (the order comes from the `##contig` lines; `-g` overrides) |
| `sort-bed --max-mem 8G --tmpdir DIR FILE` | `gloci sort FILE --max-mem 8G --tmp-dir DIR` |
| `sort-bed --check-sort FILE` | `gloci sort --check FILE` (checks chromosome and start, which is what the streaming commands need) |
| `bedops -u A B`, `sort -m A B` | `gloci sort -m A B` (inputs must already be sorted in the same order) |
| `sort-bed FILE \| uniq`, `sort -u` | `gloci sort --unique FILE` (`--duplicates` is `uniq -d`) |
| `bgzip -c FILE > FILE.gz` | `gloci view FILE -o FILE.gz` (always BGZF; `--compression-level N`, `-t` threads) |
| `bgzip -dc FILE.gz > FILE` | `gloci view FILE.gz -o FILE` (the name needs a known extension such as `x.bed`; otherwise add `-O bed`) |
| `tabix -p bed\|vcf\|gff FILE.gz`, `tabix -C ...` | `gloci index FILE.gz` (the type comes from the content; TBI, or CSI when a record ends past 2^29) |
| `bcftools index FILE.bcf`, `bcftools index -c\|-t FILE.vcf.gz` | `gloci index FILE.bcf` (a BCF always gets a CSI; `.vcf.gz` gets a TBI, or a CSI when needed) |
| `bgzip FILE && tabix -p vcf FILE.gz` | `gloci sort FILE -o FILE.gz --write-index`, or `gloci view FILE -o FILE.gz --write-index` for a file that is already sorted |
| `tabix FILE.gz REGION...`, `bcftools view -r REGION FILE`, `bcftools view -H -r ...` | `gloci view FILE.gz REGION...` (add `--no-header` for `-H`) |
| `tabix -h FILE.gz REGION`, `tabix -H FILE.gz` | `gloci view FILE.gz REGION` keeps the header; `gloci view --header-only FILE.gz` prints only it |
| `tabix -R TARGETS FILE.gz`, `bcftools view -R TARGETS FILE` | `gloci view FILE.gz --regions-file TARGETS.bed` (0-based BED, like `bedextract`) |
| `bedextract CHR FILE`, `bedextract FILE REGIONS.bed` | `gloci view FILE CHR`, or `gloci view FILE --regions-file REGIONS.bed`; add `--assume-sorted` for the binary search on an uncompressed BED file |
| `bcftools view -Ob FILE.vcf`, `bcftools view -Ov FILE.bcf` | `gloci view FILE.vcf -o FILE.bcf`, `gloci view FILE.bcf -O vcf` |

The full old-to-new tables for the P1 to P3 commands (bedtools, bedops and
vcftools) are in section 11 of the design spec
(`docs/superpowers/specs/2026-11-06-gloci-design.md`). They will move into this
guide as those commands arrive.

## Other behaviour changes

- Region and point semantics differ from tabix: a zero-length BED record or an
  insertion on a region's edge is not returned by default.
  [differences.md](differences.md#view-region-footprints); `--vcf-span ref
  --inclusive-points` reproduces tabix.
- A BGZF file without its end-of-file marker is an error, not a warning.
  [differences.md](differences.md#bgzf-eof-marker)
- Records on a contig the VCF header does not declare: `view` and `index`
  accept them with a warning, as htslib does, but `sort` and BCF output refuse
  them. [differences.md](differences.md#undeclared-contig-sort)
- `gloci index` refuses a BCF whose contigs are not in dictionary order, which
  `bcftools index` accepts. Run `gloci sort` first.
  [differences.md](differences.md#bcf-index-contig-order)
- VCF records are footprinted by their alleles by default, not by REF.
  [differences.md](differences.md#default-vcf-span-allele)
