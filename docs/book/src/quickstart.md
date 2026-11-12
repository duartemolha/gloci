# Quickstart

gloci has two command shapes:

```
gloci <command> [options] [inputs]
gloci vcf <command> [options] <input>
```

Input formats are detected from content, never from file names (`-` reads stdin, and `.gz` output is
always BGZF). `-o FILE` writes atomically, and `--write-index` builds the `.tbi` or `.csi` beside it.

```sh
gloci sort big.bed -g hg38.fa.fai --max-mem 8G -o big.sorted.bed.gz --write-index
gloci intersect -a calls.vcf.gz -b targets.bed -o on_target.vcf.gz --write-index
gloci subtract -a targets.bed -b calls.vcf.gz                  # allele-aware by default
gloci map -a genes.bed -b signal.bedgraph -c 4:mean -c 4:max
gloci merge peaks_*.bed.gz --max-gap 50 -c 5:max
gloci vcf filter cohort.bcf --pass-only --min-qual 30 --min-maf 0.01 --max-missing 0.1 -o filtered.bcf
gloci vcf stats cohort.vcf.gz --freq --missing-site --missing-indv --het --tstv --out-prefix qc/cohort
```

What each line does:

- `sort` orders a large BED within an 8 GiB memory budget and indexes the result.
- `intersect` keeps the VCF records that fall in the targets; the output is a VCF with its header.
- `subtract` removes the variants from the targets, aware of the alleles by default.
- `map` adds the mean and the maximum of column 4 of the bedGraph to each gene.
- `merge` merges the peaks of many files, with a gap of 50 allowed, keeping the highest score.
- `vcf filter` and `vcf stats` are the vcftools-style filtering and QC reports.

Unsorted input stops the command with exit 1 and a hint; `--sort` sorts first. The
[command reference](reference/sort.md) has every command's flags, and
[Migrating](migration.md) maps bedtools, bedops, vcftools and samtools invocations to gloci ones.

The same steps from Python:

```python
import gloci

(gloci.open("calls.vcf.gz")
      .intersect(b="targets.bed")
      .merge(max_gap=50)
      .save("on_target.bed.gz", write_index=True))
```
