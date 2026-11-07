# Migrating from bedtools, bedops and vcftools

This guide is for people moving existing pipelines to gloci. It starts with the
behaviour changes that break or alter a pipeline even after the flags are
rewritten.

## Chromosome order (behaviour change)

Without `-g` and without a VCF/BCF header that lists contigs, gloci uses
**karyotype order**: `chr1, chr2 ... chr22, chrX, chrY, chrM`. `sort -k1,1`,
bedops `sort-bed` and `bedtools sort` use **lexicographic** order:
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

### Worked example

Before:

```sh
sort -k1,1 -k2,2n peaks.bed > peaks.sorted.bed
sort -k1,1 -k2,2n genes.bed > genes.sorted.bed
bedtools intersect -sorted -a peaks.sorted.bed -b genes.sorted.bed -u > hits.bed
```

After, with karyotype order (`sort` and `intersect` land in later releases; the
command lines below are the intended forms):

```sh
gloci sort peaks.bed > peaks.sorted.bed
gloci sort genes.bed > genes.sorted.bed
gloci intersect -a peaks.sorted.bed -b genes.sorted.bed > hits.bed
```

After, keeping the old order so `hits.bed` is unchanged:

```sh
gloci sort peaks.bed --chrom-order lex > peaks.sorted.bed
gloci sort genes.bed --chrom-order lex > genes.sorted.bed
gloci intersect -a peaks.sorted.bed -b genes.sorted.bed --chrom-order lex > hits.bed
```

## Command mapping

The full old-to-new tables for bedtools, bedops and vcftools are in
section 11 of the design spec
(`docs/superpowers/specs/2026-11-06-gloci-design.md`). They will move into this
guide in P1.
