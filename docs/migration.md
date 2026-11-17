# Migrating from bedtools, bedops and vcftools

gloci (pronounced GLOH-sigh, like "loci") is one tool for genomic interval and
variant work.

This guide is for people moving existing pipelines to gloci. It starts with the
behaviour changes that break or alter a pipeline even after the flags are
rewritten.

## Chromosome order

Without `-g`, without a VCF/BCF header that lists contigs and without an
alignment (BAM, CRAM, SAM) input, whose header lists its references, gloci uses
**lexicographic order**, the byte order of the names: `chr1, chr10, chr11 ...
chr19, chr2, chr20 ...`. That is the order of `sort-bed`, `bedtools sort` and
`LC_ALL=C sort -k1,1` (other locales order names differently), so files sorted
by those tools work as they are, and `gloci sort x.bed` writes the same bytes as
`sort-bed` and `bedtools sort`. Within a chromosome the order is the same too.
Beside a header that sets the order, a BED must follow that order instead: for
`gloci intersect -a reads.bam -b x.bed`, sort it with `gloci sort x.bed -g
reads.bam`, or pass `--sort`.

Records that share a start are ordered by end by `gloci sort`. `sort -k1,1 -k2,2n -s`
leaves them in input order; without `-s`, GNU `sort` falls back to comparing the
whole line, so ties are not in input order either. `gloci sort` keeps input order
for records that tie on chromosome, start and end, and `--full-key` orders such
ties by the whole line (see [differences.md](differences.md#sort-tie-order)).

### Karyotype order

`--chrom-order karyotype` selects `chr1, chr2 ... chr22, chrX, chrY, chrM`
instead (a leading `chr` is ignored, everything else is in natural order). A VCF,
BCF or BAM header usually lists its contigs in that order, and `gloci convert`
keeps the header's order, so a BED made from one is in karyotype order.

### What you will see

A file in karyotype order, read without `-g` or `--chrom-order karyotype`, is
rejected by gloci's streaming commands:

```
x.bed is not sorted: line 3 (chr2	0	3) comes after line 2 (chr10	1	2)
hint: run: gloci sort x.bed (this file looks sorted in karyotype order, as a VCF or BAM header lists it; pass -g GENOME or --chrom-order karyotype, or re-sort it with gloci sort)
```

The reverse also holds: with `--chrom-order karyotype`, a file sorted by
`sort-bed` or `sort -k1,1` is rejected, and the hint says so. gloci reports this
as an error. It never silently gives a wrong answer from files in a mixed order.

### Fixes

1. Pass the genome the file follows: `-g hg38` or `-g genome.fai` (or any
   `.genome`, `.dict`, BAM header).
2. Or pass `--chrom-order karyotype` to every gloci command in the pipeline.
3. Or re-sort once, then keep the files: `gloci sort x.bed > x.sorted.bed` (or add
   `--sort` to a streaming relation, `merge` or `complement` command to sort on
   the fly).

### What still takes precedence

`-g genome.fai` (or `.genome`, `.dict`, a BAM header, a built-in such as `hg38`),
the contig order in a VCF/BCF header and the reference order in the header of
an alignment input (BAM, CRAM, SAM) win over `--chrom-order`, in that order.
`--chrom-order` only applies when none of them is present.

Names that follow no common pattern, such as yeast's Roman numerals, need `-g`
to get the reference order.

### Built-in genomes

bedtools needs a chrom-sizes file (`-g human.hg38.genome`). gloci has the common
assemblies built in, so `-g hg38` replaces it: `gloci sort x.bed -g hg38`. Files
still work, and a file named like a built-in wins (with a warning). `gloci genomes list`
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

After, with the same order and the same `hits.bed`:

```sh
gloci sort peaks.bed > peaks.sorted.bed
gloci sort genes.bed > genes.sorted.bed
gloci intersect -a peaks.sorted.bed -b genes.sorted.bed > hits.bed
```

Files that are already sorted by the old tools need no `gloci sort` step.

## Command mapping

Every row below was checked against `gloci <cmd> --help`.

Rows marked "not planned" name tools gloci does not intend to replace.
BAM, CRAM and FASTA have their own section, [BAM, CRAM and FASTA
behaviour](#bam-cram-and-fasta-behaviour), after the tables. The samtools
commands have their own table, [samtools](#samtools), with [samtools
behaviour](#samtools-behaviour) after the tables.

### Sorting, indexing and viewing (htslib, bcftools, bedops, bedtools)

| Old command | gloci |
|---|---|
| `sort-bed FILE`, `bedtools sort -i FILE`, `LC_ALL=C sort -k1,1 -k2,2n FILE` | `gloci sort FILE` (add `--chrom-order karyotype` for karyotype order) |
| `bedtools sort -i FILE -g GENOME` | `gloci sort FILE -g GENOME` (`-g` also takes a built-in such as `hg38`) |
| `vcf-sort FILE`, `bcftools sort FILE` | `gloci sort FILE` (the order comes from the `##contig` lines; `-g` overrides) |
| `sort-bed --max-mem 8G --tmpdir DIR FILE` | `gloci sort FILE --max-mem 8G --tmp-dir DIR` |
| `sort-bed --check-sort FILE` | `gloci sort --check FILE` (checks chromosome and start, which is what the streaming commands need) |
| `bedops -u A B`, `sort -m A B` | `gloci sort -m A B` (inputs must already be sorted in the same order) |
| `sort-bed FILE \| uniq`, `sort -u` | `gloci sort --unique FILE` (`--duplicates` is `uniq -d`) |
| `bcftools annotate --rename-chrs MAP FILE`, `sed 's/^1\t/chr1\t/'` | any command with `--chrom-map MAP` (a tab-separated table of old and new name; unlisted names are kept; order checks, `-g` and the output all see the new names). With `view --assume-sorted` the file must be sorted by the *new* names, see Other behaviour changes below |
| `bgzip -c FILE > FILE.gz` | `gloci view FILE -o FILE.gz` (always BGZF; `--compression-level N`, `-t` threads) |
| `bgzip -dc FILE.gz > FILE` | `gloci view FILE.gz -o FILE` (the name needs a known extension such as `x.bed`; otherwise add `-O bed`) |
| `tabix -p bed\|vcf\|gff FILE.gz`, `tabix -C ...` | `gloci index FILE.gz` (the type comes from the content; TBI, or CSI when a record ends past 2^29) |
| `tabix -s S -b B [-e E] [-0] [-S N] [-c C] FILE.gz` | `gloci index FILE.gz -s S -b B [-e E] [--zero-based] [--skip-lines N] [--comment C]` |
| `bcftools concat A B`, `bcftools concat -n A.vcf.gz B.vcf.gz` | `gloci vcf concat A B` (the output format follows `-o`/`-O`; the header is the first file's, and gloci adds no command-line header line to it). `-n`/`--naive` is `gloci vcf concat --naive A.vcf.gz B.vcf.gz`: BGZF VCF only (BCF exits 2), byte-identical headers (exit 2 otherwise), BGZF blocks copied without recompression, and the order across files is still checked. `bcftools concat -a`, `-r`, `-d` and `-D` (overlapping files, regions, deduplication) are not available; use `gloci sort` and `gloci view` |
| `bcftools norm -m- --no-version FILE`, `bcftools norm -m-snps\|-indels\|-both\|-any FILE` (no `-f`) | `gloci vcf split FILE [--types snps\|indels\|both\|any]`: one record per ALT; GT keeps the ALT as 1 and recodes the others to 0; Number=A/R/G INFO and FORMAT fields are subset (haploid G too). Alleles are not trimmed. Values keep their text, and a record that is not split is written as it is: [differences.md#vcf-split-values](differences.md#vcf-split-values). A wrong number of values for a field is an error (bcftools: `--force`). `norm -m+` (merging) is not available; with `-f`, `-d` or `-c` see `gloci vcf norm` below |
| `bcftools norm -f REF.fa --no-version FILE`, `bcftools norm -m- -f REF.fa FILE`, `-c e\|w\|x\|s`, `-d snps\|indels\|both\|all\|exact`, `-w N` | `gloci vcf norm FILE -f REF.fa [--split[=snps\|indels\|both\|any]] [--check-ref e\|w\|x\|s] [--rm-dup snps\|indels\|both\|all\|exact] [--window N]`: left-aligns and trims every record, re-sorting within `--window` bases (default 1000), and checks REF against the FASTA. The FASTA needs a `.fai` (and a `.gzi` for BGZF): `gloci index REF.fa`; gloci never builds one. `--split` splits multi-allelic records first, as `-m-`. A few edge cases differ (symbolic ALTs such as `<DEL>`, `-c x`, `-c s` on an ALT equal to REF, float text): [differences.md#vcf-norm-edge-cases](differences.md#vcf-norm-edge-cases). `-a`/`--atomize`, `-m+`, `--old-rec-tag`, `-N`, `-i`/`-e`, `-r`/`-t` and `-g` are not available |
| `bcftools annotate --no-version -a SRC -c LIST FILE`, `-x LIST`, `-h FILE`, `-I FORMAT`, `-m [+-]TAG`, `--min-overlap` | `gloci vcf annotate FILE [--annotations SRC --columns LIST] [--remove LIST] [--header-lines FILE] [--set-id FORMAT] [--mark-sites [+-]TAG] [--min-overlap ANN[:VCF]]`: transfers ID, QUAL, FILTER and INFO from a VCF/BCF, tab or BED source; removes fields and their header lines; sets IDs; marks sites. The rules (which records pair, `+TAG`, `NEW:=OLD`, Number=A/R by allele, the first matching row, `^` lists, `-x ^FORMAT/TAG`) are bcftools' 1.21. Both files must be sorted in one chromosome order and are streamed, no index needed ([differences.md#vcf-annotate-sorted-source](differences.md#vcf-annotate-sorted-source), [differences.md#vcf-annotate-lenient](differences.md#vcf-annotate-lenient)); values keep their text ([differences.md#vcf-annotate-values](differences.md#vcf-annotate-values)). `-h` is `--header-lines` (long form only; `-h` is help). `-c -TAG` and `-c =TAG` (bcftools: not implemented / `Number=.`), FORMAT transfers, `ALT`/`REF` replacement, `~ID`, `-s`/`-S`, `-i`/`-e` (pipe through `gloci vcf filter`), `--rename-chrs` (`gloci view FILE --chrom-map MAP`), `-l`/`--merge-logic`, `--single-overlaps`, `--rename-annots`, `--pair-logic`, `-C`, `-H`, `-k`, `-r`/`-R` (select regions first: `gloci view FILE REGION \| gloci vcf annotate - ...`) and `-I %ALT` are not available. `-c ...,ALT` without `REF` is accepted and ignored, as in bcftools; a source chromosome missing from FILE's `##contig` lines is an error. A piped or process-substituted FILE or source works |
| `bcftools index FILE.bcf`, `bcftools index -c\|-t FILE.vcf.gz` | `gloci index FILE.bcf` (a BCF always gets a CSI; `.vcf.gz` gets a TBI, or a CSI when needed) |
| `bgzip FILE && tabix -p vcf FILE.gz` | `gloci sort FILE -o FILE.gz --write-index`, or `gloci view FILE -o FILE.gz --write-index` for a file that is already sorted |
| `tabix FILE.gz REGION...`, `bcftools view -r REGION FILE`, `bcftools view -H -r ...` | `gloci view FILE.gz REGION...` (or `--region REGION`; add `--no-header` for `-H`) |
| `tabix -h FILE.gz REGION`, `tabix -H FILE.gz` | `gloci view FILE.gz REGION` keeps the header; `gloci view --header-only FILE.gz` prints only it |
| `tabix -R TARGETS FILE.gz`, `bcftools view -R TARGETS FILE` | `gloci view FILE.gz --regions-file TARGETS.bed` (0-based BED, like `bedextract`) |
| `bedextract CHR FILE`, `bedextract FILE REGIONS.bed` | `gloci view FILE CHR`, or `gloci view FILE --regions-file REGIONS.bed`; add `--assume-sorted` for the binary search on an uncompressed BED file |
| `bcftools view -Ob FILE.vcf`, `bcftools view -Ov FILE.bcf` | `gloci view FILE.vcf -o FILE.bcf`, `gloci view FILE.bcf -O vcf` |
| `bcftools query -f FMT [-H] [-i/-e EXPR] [-s S\|-S FILE] [-r R\|-R BED] [-t T\|-T BED] FILE` | `gloci vcf query FILE -f FMT [-H] [--include/--exclude EXPR] [--samples S\|--samples-file FILE] [--region R\|--regions-file BED]`: bcftools 1.21's fields, `[ ]`, `{N}` and `\t \n \\` (anything else exits 2); `-i/-e` on GT or FORMAT print only the passing samples inside `[ ]`, as bcftools does, and `&&` there lists the samples that pass either test (bcftools' per-sample `&` and `\|` are not available); the expression subset is `vcf filter`'s; `-t/-T` are `--region`/`--regions-file`; `-HH`, `-u`, `-a` and `--force-samples` are not available; `-l` is `gloci view --header-only FILE \| tail -n1 \| cut -f10- \| tr '\t' '\n'`; differences: [values](differences.md#vcf-query-values), [subscript past the end](differences.md#vcf-query-index-past-end), [`-H` names](differences.md#vcf-query-header-info), [sample subscripts](differences.md#vcf-query-samples), [strict formats](differences.md#vcf-query-format-strict) |

### bedtools

| bedtools | gloci |
|---|---|
| `intersect -u` | `gloci intersect -a A -b B` (default `--mode any`; note `--mode` has no `-m`) |
| `intersect -abam A.bam -b B.bed -u` (`-v`, `-split`, `-s`, `-S`, `-f`, `-c`) | `gloci intersect -a A.bam -b B.bed [--mode none] [--split] [-s] [-S] [-f FRAC]` (BAM, CRAM or SAM as A or B, which writes BAM (also for a CRAM A, whose reads need `--reference FASTA`), or SAM for a SAM A, and `-o x.bam` or `-O bam\|sam` set it; output is the same reads as bedtools 2.31.1; `--split` cuts at CIGAR `N`, and `--split-deletions` also at `D`; unmapped reads are never written, even by `--mode none`, see [differences.md#unmapped-reads](differences.md#unmapped-reads); `-abam -bed -wa -wb` and `--mode pairs` with an alignment A print BED6, not BED12, see [differences.md#bam-as-bed-columns](differences.md#bam-as-bed-columns); `--mode clip` of reads is an error; for the BED that `-abam -bed` writes, pipe the reads: `gloci intersect -a A.bam -b B.bed -O bam \| gloci convert - -o hits.bed`) |
| `intersect -v` | `gloci intersect -a A -b B --mode none` |
| `intersect` (no mode flag) | `gloci intersect -a A -b B --mode clip` |
| `intersect -wa -wb` | `gloci intersect -a A -b B --mode pairs` |
| `intersect -wo` | `gloci intersect -a A -b B --mode pairs --overlap-bases` |
| `intersect -wao` / `-loj` | `gloci intersect -a A -b B --mode pairs --keep-unmatched [--overlap-bases]` |
| `intersect -c` | `gloci intersect -a A -b B --mode count` |
| `-f -F -r -e -s -S -split` | `-f` (`--min-frac-a`), `-F` (`--min-frac-b`), `-r` (`--reciprocal`), `-e` (`--either`), `-s` (`--same-strand`), `-S` (`--opposite-strand`), `--split` (plus `--min-bases N`, `--inclusive-points`) |
| `-sorted` | default (sortedness is checked); unsorted input needs `--unsorted` (loads `-b` into memory, keeps A's order) or `--sort` (external sort first) |
| `window -w/-l/-r/-sw` | `gloci intersect -a A -b B -w/--window N` or `--left L --right R` with `--strand-window` (also on `subtract` and `map`) |
| `subtract`, `-A`, `-N` | `gloci subtract -a A -b B`, `--whole`, `--sum-overlaps [-f FRAC]` (plus `--split` to cut by B's BED12 blocks); `--match alleles` on VCF/BCF drops an A only for a B with the same allele (`bcftools isec -C -c none`, trimmed: [differences.md#trimmed-allele-match](differences.md#trimmed-allele-match)) |
| `merge -d -s -c -o` | `gloci merge FILE... -d/--max-gap N -s --count -c/--column COL:OP` (`--delim`, `--empty`, `--sort`) |
| `merge -i x.bam` (BAM input to `merge`) | `gloci convert x.bam \| gloci merge - -g x.bam` (the BED keeps the BAM header's chromosome order, which `-g x.bam` names; the BAM must be coordinate-sorted, or use `--sort` in place of `-g x.bam`; `gloci merge` reads intervals, not alignments, and says so, naming `convert`; to merge sorted BAM files as `samtools merge` does, use `gloci sort -m`) |
| `--cram-ref REF.fa` (`intersect`, `coverage`, `genomecov`, `bamtobed`, any bedtools command reading a CRAM) | `--reference REF.fa` on `intersect`, `coverage`, `genomecov`, `convert`, `view`, `sort`, `index` (optional there). There is no `REF_PATH`, cache or download. See [BAM, CRAM and FASTA behaviour](#bam-cram-and-fasta-behaviour) |
| `complement -g` | `gloci complement FILE... [-g GENOME]` (`-g` takes a built-in such as `hg38` or a `.genome`/`.fai`/`.dict` file; `--sort`) |
| `closest -d -D -k -t -io -iu -id -mdb` | `gloci closest -a A -b B --distance --signed ref\|a\|b -k N --ties all\|first\|last --ignore-overlaps --ignore-upstream --ignore-downstream --per-file` (plus `--max-distance N`, `-s`/`--same-strand`, `-S`/`--opposite-strand`) |
| `map -c -o` | `gloci map -a A -b B -c/--column COL:OP` (`--delim`, `--empty`, `--skip-unmatched`; ops: `count`, `count_distinct`, `sum`, `mean`, `median`, `min`, `max`, `stdev`, `var`, `stdev_pop`, `var_pop`, `collapse`, `distinct`, `distinct_sorted`, `first`, `last`, `bases`, `indicator`, `cv`, `mad`, `quantile:Q`, `tmean:LO-HI`, `wmean`, `mode`, `antimode`, `bases_uniq`, `frac_a`, `frac_b`, `min_record`, `max_record`, `records`) |
| `bcftools merge --no-version [-m MODE] [--force-samples] [-i RULES] [-0] [-f LIST] [-F x\|+] [-l FILE] A.vcf.gz B.vcf.gz ...` | `gloci vcf merge A B ... [--merge none\|snps\|indels\|both\|all\|id] [--force-samples] [--info-rules RULES] [--missing-to-ref] [--apply-filters LIST] [--filter-logic x\|+] [--file-list FILE]`: combines the samples of sorted VCF/BCF files; the records at one position are paired, their ALTs united, GT and Number=A/R/G values re-indexed and INFO, QUAL, FILTER and ID merged by bcftools' 1.21 rules (`-m both` default; AC and AN recomputed from the merged genotypes, `DP`/`DP4` summed). The inputs are sorted and streamed, no index is needed ([differences.md#vcf-merge-no-index](differences.md#vcf-merge-no-index)); contigs are written in the header's order ([differences.md#vcf-merge-contig-order](differences.md#vcf-merge-contig-order)); conflicting tag types stop the merge ([differences.md#vcf-merge-conflicting-types](differences.md#vcf-merge-conflicting-types)) and an undeclared INFO tag is accepted ([differences.md#vcf-merge-undeclared-info](differences.md#vcf-merge-undeclared-info)); which records of one position share a row and in what order is bcftools' (`can_merge`/`stage_line`, ported). The output keeps the first input's format (`-O vcf` for text from BCF inputs). `-l` is `--file-list`. `-g`/`--gvcf`, `-L`/`--local-alleles`, `-M`/`--missing-rules`, `--print-header`, `--use-header`, `-r`/`-R`/`--regions-overlap` (select regions after merging: `gloci vcf merge ... \| gloci view - REGION`) and `-m snp-ins-del`, `*`, `**` are not available |
| `coverage [-hist -d -mean]`, `annotate` | `gloci coverage -a A -b B [--hist \| --per-base \| --mean]`, `gloci coverage -a A -b B1 B2... [--names ...]` |
| `coverage -a A.bed -b reads.bam [-split -hist -d -mean]` | `gloci coverage -a A.bed -b reads.bam [--split] [--hist \| --per-base \| --mean]` (BAM, CRAM or SAM as B, a CRAM with `--reference FASTA`; `--split-deletions` also cuts reads at `D`; `--min-mapq`, `--require-flags` and `--exclude-flags` filter reads; an alignment A prints as BED6, see [differences.md#bam-as-bed-columns](differences.md#bam-as-bed-columns)) |
| `slop -b -l -r -s -pct -g` | `gloci slop FILE -b/--both N` or `--left L --right R` with `-s`, `--pct`, `-g GENOME` |
| `flank -b -l -r -s -pct -g` | `gloci flank FILE -b/--both N` or `--left L --right R` with `-s`, `--pct`, `-g GENOME` |
| `shift -s N`, `shift -p X -m Y`, `-pct`, `-g` | `gloci shift FILE --by N [-g GENOME]` (signed; `--pct` reads N as an exact decimal fraction of the feature length); `-p X -m -X` is `--by X -s`. bedtools' `-p`/`-m` with unrelated magnitudes has no gloci equivalent: run `shift` per strand after filtering |
| `makewindows -g/-b -w -s -n -i` | `gloci windows -g GENOME \| -i/--intervals FILE` with `-w/--size SIZE [--step S]` or `-n/--count COUNT`, plus `--id-column`, `--reverse` |
| `sort -i FILE [-g GENOME]` | `gloci sort FILE [-g GENOME]` (byte order is the default; `--chrom-order karyotype` gives karyotype order) |
| `bed12tobed6 -i FILE` | `gloci convert FILE --blocks` |
| `genomecov -i x.bed -g G` (`-bg`, `-bga`, `-d`, `-scale`, `-strand`, `-split`, `-5`, `-3`, `-max`) | `gloci genomecov x.bed -g G [--bedgraph \| --bga \| --per-base] [--scale X] [--strand +\|-] [--split] [--5p \| --3p] [--max-depth N]` (`-g` is required, as in bedtools; the input must be sorted, `--sort` sorts it; `-dz` and `-trackline` have no counterpart; the `--hist` fractions divide in single precision and print six significant digits, as bedtools does; a `genome` row that counts more than 2^32 bases at one depth is correct in gloci, see the `genomecov` note under Other behaviour changes; BAM input is `-ibam`, see the next row) |
| `genomecov -ibam x.bam` (`-bg`, `-bga`, `-d`, `-scale`, `-strand`, `-split`, `-5`, `-3`, `-max`, `-ignoreD`) | `gloci genomecov x.bam [--bedgraph \| --bga \| --per-base] [--scale X] [--strand +\|-] [--split \| --split-deletions] [--5p \| --3p] [--max-depth N]` (the chromosomes and sizes come from the header, in header order with empty chromosomes included, where bedtools writes a contig without reads after the others, see [differences.md#genomecov-chrom-order](differences.md#genomecov-chrom-order); `-g` overrides them where bedtools ignores it, see [differences.md#genomecov-bam-genome-arg](differences.md#genomecov-bam-genome-arg); a read with no reference operations counts one base, see [differences.md#bamtobed-no-reference-ops](differences.md#bamtobed-no-reference-ops); reads are cut at the end of their reference and unmapped reads are left out, counted at `-v`; `--min-mapq`, `--require-flags` and `--exclude-flags` filter first; bedtools cuts a read at deletions unless `-ignoreD`, gloci counts deleted bases by default, so `-split` is `--split-deletions`, `-split -ignoreD` is `--split`, and `-ignoreD` is the default, see [differences.md#bam-genomecov-deletions](differences.md#bam-genomecov-deletions); an indexed BAM, or a CRAM with a `.crai` (and `--reference FASTA`), at `-t` above 1 runs per chromosome with the same output; `-pc`, `-fs` and `-du` are not planned) |
| `multiinter -i A B C [-names]` | `gloci partition A B C --merge-adjacent [--names X,Y,Z]` (the count and label columns; the per-file 0/1 columns are not written; without `-names` bedtools labels the files `1,2,3` and gloci uses the file names (pass `--names 1,2,3` to match); `-cluster` and `-empty -g` have no counterpart (for the latter use `gloci complement -g`)) |
| `unionbedg -i A B [-filler -empty -g]` | `gloci partition A B --bedgraph [--empty VALUE [-g GENOME]]` (`-filler X` is `--empty X`; `--empty` also writes the uncovered segments, as it implies `--min-files 0`; an explicit `--min-files` overrides that) |
| `groupby -g -c -o [-header]` | `gloci groupby FILE -g COLS -c/--column COL:OP` (bedtools `-o OP` becomes the `:OP` of `-c COL:OP`, and gloci's `-o` is the output file; `--header`, `--delim`, `--empty`; same ops as `map` except those that need a second input; `distinct` keeps first-seen order where bedtools sorts it, use `distinct_sorted`, see [differences.md#map-distinct-order](differences.md#map-distinct-order); VCF input is refused, use `gloci vcf extract`) |
| `cluster -d -s` | `gloci cluster FILE -d/--max-gap N -s` (`--sort`; the same gap and strand rules as `gloci merge`; ids are global, and `-s` keeps input order, see [differences.md#cluster-strand](differences.md#cluster-strand); a VCF input is clustered by its footprints and the output is TSV with a `#CHROM ... cluster` line, see [differences.md#append-vcf-header](differences.md#append-vcf-header); BCF is refused, use `gloci view -O vcf`; bedtools accepts chromosomes in any grouped order, gloci needs a sorted chromosome order, see [Chromosome order](#chromosome-order)) |
| `spacing` | `gloci spacing FILE` (`.` for the first record of a chromosome, `0` bookended, `-1` overlapping, only the previous record counts, as in bedtools; a zero-length record is a point, see [differences.md#spacing-zero-length](differences.md#spacing-zero-length); VCF input and chromosome order as for `cluster`) |
| `bamtobed -i FILE` (`-bed12`, `-split`, `-splitD`, `-ed`, `-tag TAG`, `-cigar`, `-color R,G,B`, `-bedpe`, `-mate1`) | `gloci convert FILE [--bed12] [--blocks] [--blocks --split-deletions] [--score-tag NM\|TAG] [--cigar] [--color R,G,B] [--bedpe [--mate1-first]]` (BAM, CRAM or SAM input, a CRAM with `--reference FASTA`, which bedtools reads as BAM only; output matches bedtools 2.31.1 in input order except for the pinned differences below, with `/1` and `/2` on mates, MAPQ as the score and unmapped reads skipped; `-ed` is `--score-tag NM`; `--min-mapq`, `--require-flags` and `--exclude-flags` filter; a read with no reference bases is one base, not zero-length, see [differences.md#bamtobed-no-reference-ops](differences.md#bamtobed-no-reference-ops); reads are cut at the end of their reference, and a BED12 row opens at its first kept block and closes at its last, with no size-0 blocks, see [differences.md#bamtobed-past-reference-end](differences.md#bamtobed-past-reference-end); a negative tag score is printed signed, where bedtools wraps it to unsigned 32-bit, see [differences.md#bamtobed-negative-score](differences.md#bamtobed-negative-score); `--blocks --score-tag` has no extra column and `--cigar` with `--blocks` is an error, see [differences.md#bamtobed-split-tag](differences.md#bamtobed-split-tag); `-bedpe` is `gloci convert FILE --bedpe` on a name-sorted file (`gloci sort -n`; a coordinate-sorted one is exit 2), `-mate1` is `--mate1-first` and `-ed` is `--score-tag NM`; the output is BEDPE (`.bedpe`, `.bedpe.gz`) and feeds `pair-intersect` and `pair-pair`; it matches bedtools, except that secondary and supplementary alignments are not mates, see [differences.md#bamtobed-bedpe](differences.md#bamtobed-bedpe)) |
| `getfasta -fi REF -bed FILE` (`-name`, `-nameOnly`, `-s`, `-split`, `-tab`, `-bedOut`, `-rna`) | `gloci getfasta FILE --fasta REF [--name \| --name-only] [-s] [--split] [--tab \| --bed-out] [--rna]` (output is byte-identical to bedtools 2.31.1 for BED, GFF and VCF input, including the IUPAC-aware reverse complement, except that `--split` refuses BED12 records whose blocks are out of order, overlap, or do not span the record where bedtools writes them unchecked, see [differences.md#getfasta-split-blocks](differences.md#getfasta-split-blocks); the FASTA needs a `.fai` (and `.gzi` for BGZF) made with `gloci index`, which gloci never builds itself; VCF is bedtools' REF span with `--vcf-span ref`; zero-length intervals are skipped with a counted warning, as bedtools does; an interval on a chromosome the FASTA lacks, or past its end, is an error where bedtools skips it, `--skip-invalid` gives bedtools' behaviour, see [differences.md#getfasta-outside-fasta](differences.md#getfasta-outside-fasta); there is no `-fullHeader`) |
| `maskfasta -fi F -bed B [-soft -mc C -fullHeader]` | `gloci maskfasta --fasta F -b B [--soft \| --mask-char C] [-o OUT.fa[.gz] --write-index]` (streamed, any compression, no index needed; headers are kept whole and line layout is kept, so `-fullHeader` is the default, see [differences.md#maskfasta-header](differences.md#maskfasta-header) and [differences.md#maskfasta-layout](differences.md#maskfasta-layout); a mask past its sequence's end is an error, see [differences.md#maskfasta-past-end](differences.md#maskfasta-past-end)) |
| `nuc -fi REF -bed FILE` (`-s`, `-seq`, `-pattern P`, `-C`) | `gloci nuc FILE --fasta REF [-s] [--seq] [--pattern P [-C]]` (output is byte-identical to bedtools 2.31.1, header line and f32 percentages included, for BED, GFF and VCF input; the counts and the pattern describe the reverse complement under `-s`; the FASTA needs a `.fai` made with `gloci index`; VCF is bedtools' REF span with `--vcf-span ref`; zero-length intervals are skipped with a counted warning, as bedtools does; an interval outside the FASTA is an error where bedtools skips it, `--skip-invalid` gives bedtools' behaviour, see [differences.md#getfasta-outside-fasta](differences.md#getfasta-outside-fasta)) |
| `jaccard -a A -b B [-f -F -r -e -s -S]` | `gloci stats jaccard -a A -b B` (same filters, plus `--min-bases`, `--inclusive-points`, `--sort`; header line and one row, `%g` numbers; the point rule is [differences.md#jaccard-points](differences.md#jaccard-points)) |
| `fisher -a A -b B -g G [-m -s -S -f -F -r -e]` | `gloci stats fisher -a A -b B -g G [-m/--merge-overlaps] [--exclude BED]` (`-g` is required and must have lengths; bedtools' text report; `-S` is a flag here; the empty-input estimate is [differences.md#fisher-empty-estimate](differences.md#fisher-empty-estimate), the `--exclude` estimate [differences.md#fisher-exclude-estimate](differences.md#fisher-exclude-estimate), tiny tail p-values [differences.md#fisher-lgamma](differences.md#fisher-lgamma)) |
| `reldist -a A -b B [-detail]` | `gloci stats reldist -a A -b B [--detail]` (midpoints past 2^31: [differences.md#reldist-midpoint-overflow](differences.md#reldist-midpoint-overflow)) |
| `pairtobed -a PAIRS -b FILE [-type T] [-f F] [-s \| -S]` | `gloci pair-intersect -a PAIRS.bedpe -b FILE [--type T] [-f F] [-s \| -S]` (`-type` becomes `--type`, same nine values; the output matches bedtools 2.31.1, including `-f`, strands and the BEDPE column counts 6, 7, 8 and 10 or more, except for the pinned differences below; B is BED, GFF or VCF, held in memory (about 5 to 6 times its size), in any order; BEDPE is never detected, so A is read as BEDPE whatever it looks like, and a VCF, BCF, GFF or alignment A is an error, with `gloci convert x.bam --bedpe` for alignments; hits come in B's input order, not bedtools' bin order, see [differences.md#pair-hit-order](differences.md#pair-hit-order); an unknown strand never matches `-s`/`-S`, see [differences.md#pair-unknown-strand](differences.md#pair-unknown-strand); `notboth` pads a pair with no hit, see [differences.md#pair-notboth-padding](differences.md#pair-notboth-padding); points follow section 4.3, see [differences.md#pair-points](differences.md#pair-points); the spans of overlapping or nested ends, see [differences.md#pair-span-overlapping-ends](differences.md#pair-span-overlapping-ends); bedtools accepts some lines gloci refuses, see [differences.md#bedpe-strictness](differences.md#bedpe-strictness); span types leave out a pair with an absent end, see [differences.md#pair-span-absent-ends](differences.md#pair-span-absent-ends); also `-F`, `-r`, `-e`, `--min-bases`, `--inclusive-points`, `--vcf-span`, `--chrom-map`, `--skip-invalid` and `-t`, which bedtools lacks; no `-abam`, `-ubam`, `-bedpe` or `-ed`: convert alignments with `gloci convert x.bam --bedpe`) |
| `pairtopair -a PAIRS -b PAIRS [-type T] [-slop N] [-ss] [-is] [-rdn] [-f F]` | `gloci pair-pair -a PAIRS.bedpe -b PAIRS.bedpe [--type T] [-w N \| --left N --right N [--strand-window]] [-s \| -S] [--require-different-names]` (`-type` becomes `--type`, values `both` (default), `either`, `neither`, `notboth`; bedtools enforces strand unless `-is`, gloci ignores it unless `-s`/`-S`, and an unknown strand never matches, see [differences.md#pair-unknown-strand](differences.md#pair-unknown-strand); `-slop N` pads A in bedtools and `-w N` pads B in gloci, which finds the same pairs, and `-ss -slop N` is `--left N --strand-window`; `-rdn` becomes `--require-different-names`; a B that hits in both orientations is written once, see [differences.md#pair-pair-once](differences.md#pair-pair-once); `either` writes B as it reads, see [differences.md#pair-pair-b-columns](differences.md#pair-pair-b-columns); hits in B's input order, see [differences.md#pair-hit-order](differences.md#pair-hit-order); bedtools accepts some lines gloci refuses, see [differences.md#bedpe-strictness](differences.md#bedpe-strictness); also `-F`, `-r`, `-e`, `--min-bases`, `--inclusive-points`, `--chrom-map`, `--skip-invalid` and `-t`; both files are BEDPE and B is held in memory, about 5 to 6 times its size) |
| `shuffle -i FILE -g G [-seed S -chrom -chromFirst -excl BED -incl BED -noOverlapping -maxTries N]` | `gloci shuffle FILE -g G [--seed S] [--same-chrom] [--chrom-first] [--exclude BED] [--include BED] [--no-overlapping] [--max-tries N] [--skip-unplaceable] [--header]` (other places than bedtools picks for the same seed, [differences.md#rng-stream](differences.md#rng-stream); the whole record stays inside `--include`, and a record that fits nowhere is exit 1 where bedtools drops it, `--skip-unplaceable` drops it, see [differences.md#shuffle-placement](differences.md#shuffle-placement); BEDPE is opt-in with `--format bedpe` and gloci moves both ends where bedtools moves only end 1, see [differences.md#shuffle-bedpe-ends](differences.md#shuffle-bedpe-ends); output is in input order and unsorted; `-f`, `-allowBeyondChromEnd` and `-incl` strand copying are not planned; VCF and BAM input are refused, use `gloci convert`) |
| `random -g G -l L -n N -seed S` | `gloci random -g G -l L -n N --seed S` (BED6 in id order; same seed gives other intervals than bedtools, [differences.md#rng-stream](differences.md#rng-stream); `-l` longer than every chromosome is exit 2) |
| `sample -i FILE -n N -seed S [-s forward\|reverse] [-header]` | `gloci sample FILE -n N --seed S [--strand +\|-] [--header]` (output in input order, [differences.md#sample-order](differences.md#sample-order); other records than bedtools picks, [differences.md#rng-stream](differences.md#rng-stream); `-ubam`/BAM: `samtools view -s`) |
| `multicov -bams A.bam B.bam -bed X.bed` | `gloci bedcov X.bed A.bam B.bam --read-count` (one count column per file, after the base-depth sum; each file needs an index; the filters are samtools' defaults, so add `--exclude-flags 0` for every read, see [the samtools behaviour section](#samtools-behaviour)) |
| `igv -i F [-path D] [-sess S] [-sort K] [-clps] [-name] [-slop N] [-img T]` | `gloci igv F [--path D] [--session S] [--sort K] [--collapse] [--name] [--slop N] [--image T]` (the script is byte for byte bedtools'; `--image` takes `png`, `eps` or `svg` and `--slop` a whole number from 0, see [differences.md#igv-option-values](differences.md#igv-option-values); nothing is written when the input cannot be recognised, see [differences.md#browser-unrecognised-input](differences.md#browser-unrecognised-input); every BED or GFF record needs the same number of columns, as in bedtools; `-i -` is `gloci igv -`) |
| `links -i F [-base U] [-org O] [-db B]` | `gloci links F [--base U] [--org O] [--db B]` (the page is byte for byte bedtools', hg18 links by default; the title is the input's name, `stdin` for `-`; nothing is written when the input cannot be recognised, see [differences.md#browser-unrecognised-input](differences.md#browser-unrecognised-input)) |
| `bamtofastq -i FILE -fq R1.fq [-fq2 R2.fq] [-tags]` | `gloci convert FILE --fastq [-o R1.fq] [--fastq2 R2.fq] [--tags]` (BAM, CRAM with `--reference FASTA`, or SAM input; `.gz` output names compress; without `-fq2` bedtools writes every read twice and gloci once, [differences.md#bamtofastq-single-twice](differences.md#bamtofastq-single-twice); `--fastq2` needs input sorted by name, `gloci sort -n`, [differences.md#bamtofastq-coordinate-sorted](differences.md#bamtofastq-coordinate-sorted), refuses a name with three records, [differences.md#bamtofastq-three-records](differences.md#bamtofastq-three-records), writes no empty read for an unpaired first record, [differences.md#bamtofastq-first-unpaired](differences.md#bamtofastq-first-unpaired), and `--tags` is [differences.md#bamtofastq-tags-missing](differences.md#bamtofastq-tags-missing) and [differences.md#bamtofastq-tags-orphan](differences.md#bamtofastq-tags-orphan); every read is written, so add `--exclude-flags 0x900` to leave out secondary and supplementary alignments) |
| `bedtobam -i FILE -g G [-mapq N] [-bed12] [-ubam]` | `gloci convert FILE --bam -g G [--mapq N] [--bed12] [--compression-level 0] [-o OUT.bam]` (BED, GFF or VCF input; BAM to stdout or `-o`, the same header and records as bedtools but for the `@PG` line, see the introduction of [differences.md](differences.md); `-g` also takes a built-in genome, a `.fai` or a `.dict`; `-ubam` is `--compression-level 0`, which bedtools 2.31.1 compresses anyway; a chromosome the genome lacks is exit 1 where bedtools puts the record on the first reference, see [differences.md#bedtobam-missing-chromosome](differences.md#bedtobam-missing-chromosome), and `--chrom-map` renames the input's chromosomes first; a zero-length record at 0 is the base after it, see [differences.md#bedtobam-zero-length-at-start](differences.md#bedtobam-zero-length-at-start); `--bed12` refuses a record that is not valid BED12, see [differences.md#bedtobam-bed12-blocks](differences.md#bedtobam-bed12-blocks); a name of more than 254 characters is exit 1, see [differences.md#bedtobam-unrepresentable](differences.md#bedtobam-unrepresentable)) |
| `bedpetobam -i FILE -g G [-mapq N] [-ubam]` | `gloci convert FILE --format bedpe --bam -g G [--mapq N] [--compression-level 0] [-o OUT.bam]` (`--format bedpe` is needed, as BEDPE is never detected; the `@SQ` lines keep the genome's order where bedtools sorts them by name, see [differences.md#bedpetobam-header-order](differences.md#bedpetobam-header-order); an absent end (`. -1 -1`) is an unmapped read placed at its mate, see [differences.md#bedpetobam-absent-end](differences.md#bedpetobam-absent-end); 7 and 8 columns are read, see [differences.md#bedpetobam-columns](differences.md#bedpetobam-columns); a strand other than `+`, `-` or `.` is exit 1, see [differences.md#bedpe-strictness](differences.md#bedpe-strictness); a chromosome the genome lacks is exit 1, see [differences.md#bedtobam-missing-chromosome](differences.md#bedtobam-missing-chromosome)) |
| `tag -i IN.bam -files A.bed B.bed -labels LA LB [-tag XX] [-names \| -scores \| -intervals] [-s \| -S] [-f F]` | `gloci tag IN.bam -b A.bed B.bed --labels LA,LB [--tag XX] [--names \| --scores \| --intervals] [-s \| -S] [-f F]` (`-files` becomes `-b` and `-labels` takes one comma-separated list; the tags are bedtools 2.31.1's: the same reads tagged, the same values in the same order, bedtools' bin order for several features of one file, a second tag beside an existing one; under `-s` and `-S` a feature of unknown strand (`.`, none, or any text but `+` or `-`) matches neither, as in `gloci intersect`, where bedtools counts it as the other strand under `-S`, see [differences.md#tag-unknown-strand](differences.md#tag-unknown-strand); the records are written in input order and otherwise unchanged, as BAM for BAM and CRAM input and SAM for SAM, `-o x.sam` or `-O sam` for SAM text; the features are BED or GFF, held in memory (about 140 bytes each), and a VCF is converted first (`gloci convert x.vcf -o x.bed`), see [differences.md#tag-feature-formats](differences.md#tag-feature-formats); `--tag` must be two characters, `--labels` one per `-b` and `-f` above 0 and at most 1, see [differences.md#tag-option-values](differences.md#tag-option-values); also CRAM and SAM input, `--reference`, `--write-index` and `-t`, which bedtools lacks) |
| `expand`, `split`, `summary`, `overlap` | not planned (`overlap` → `gloci intersect --mode pairs --overlap-bases`) |

### samtools

| samtools | gloci |
|---|---|
| `samtools index FILE.bam`, `samtools index -c FILE.bam` | `gloci index FILE.bam` (a `.bai`, or a `.csi` when a reference is longer than 2^29); `gloci index FILE.bam --csi` always writes a `.csi`. The BAM must be coordinate-sorted. SAM text has no index: `gloci view FILE.sam -o FILE.bam --write-index`. samtools and gloci read each other's indexes. [differences.md#bam-index-long-reference](differences.md#bam-index-long-reference) |
| `samtools view FILE.bam REGION...` | `gloci view -O sam FILE.bam REGION...` (through `FILE.bam.bai`, `FILE.bai` or `FILE.bam.csi`, or for a CRAM `FILE.cram.crai` or `FILE.crai`, which also needs `--reference FASTA`; without an index gloci scans the file and warns, where samtools exits 1; [differences.md#bam-region-without-index](differences.md#bam-region-without-index)). `gloci view -O sam` writes the header, as `samtools view -h` does; add `--no-header` for `samtools view` without `-h`. Regions are 1-based and closed, a read is selected by where its CIGAR puts it, and a placed unmapped read counts as one base, as samtools finds them |
| `samtools sort IN.bam -o OUT.bam`, `samtools sort -n`, `samtools sort -N` | `gloci sort IN.bam -o OUT.bam`, `gloci sort -n IN.bam -o OUT.bam`, `gloci sort -N IN.bam -o OUT.bam`. The order is samtools 1.21's: reference in header order with unplaced reads last, position, strand, then input order; `-n` is natural name order (`r2` before `r10`, `r01` equal to `r1`) and `-N` plain byte order, each with READ1 before READ2 and primary before supplementary before secondary. The `@HD` line gets `SO` (and `SS:queryname:natural` or `SS:queryname:lexicographical`) as samtools sets it. gloci keeps the input's format (BAM for BAM and CRAM, since CRAM is never written, and a CRAM needs `--reference FASTA`; SAM for SAM), where samtools always writes BAM: pass `-O bam` for samtools' default (`-o x.bam` and `-o x.sam` also choose); `--write-index` writes the `.bai` (or `.csi`) for coordinate order. The flag names differ: samtools' `-@ N` is `-t N`, `-m SIZE` (memory per thread) is `--max-mem SIZE` for the whole sort (gloci's `-m` merges sorted inputs), `-T PREFIX` is `--tmp-dir DIR`, `-l N` is `--compression-level N`, and `-u` is `--compression-level 0`. Not available: `-t TAG`, `-M` and `--template-coordinate`. `gloci sort --check FILE` checks the order without writing |
| `samtools merge -c -p -o OUT.bam A.bam B.bam`, `samtools merge -n`, `samtools merge -N` | `gloci sort -m A.bam B.bam -o OUT.bam`, `gloci sort -m -n A.bam B.bam -o OUT.bam`, `gloci sort -m -N A.bam B.bam -o OUT.bam`. The inputs (BAM, SAM, or CRAM with `--reference`) must already be in the order; each is checked as it is read, and one whose `@HD SO`/`SS` says another order is refused up front. The keys and ties are `gloci sort`'s (ties go to the earlier input, then record order, as samtools). The header follows `samtools merge -c -p`: `@RG` and `@PG` lines with the same ID are kept once, `@CO` lines come from every input, `@HD` is the first input's with `SO`/`SS` set. Two differences: `@SQ` lists must be identical ([differences.md#alignment-merge-sq](differences.md#alignment-merge-sq)) and a same-ID `@RG` with different fields is an error ([differences.md#alignment-merge-rg](differences.md#alignment-merge-rg)); samtools' default (no `-c`) renames a colliding ID with a random suffix, which gloci does not. The output follows the inputs' format: SAM when every input is SAM, otherwise BAM (a CRAM input, or a mix of formats, writes BAM); `-o`/`-O` choose. Very large input counts are cascaded: at most 128 inputs (fewer under a small `--max-mem`) are open at a time, with the `-t` threads shared among them, and the rounds' runs go through `--tmp-dir`, so open files and memory stay bounded (1,100 inputs merge under `ulimit -n 1024`); a pipe input is spooled once to `--tmp-dir` first. Not available: `-r`, `-R`, `-h`, `-b`, `-f` |
| `samtools view -b -o X.bam --write-index FILE`, then `samtools index X.bam` | `gloci view FILE -o X.bam --write-index` (writes `X.bam.bai`; the data must be coordinate-sorted). An old `X.bam.bai`, `X.bai` or `X.bam.csi` beside an output written without `--write-index` is removed |
| `samtools cat -o X.bam A.bam B.bam` | `gloci cat A.bam B.bam -o X.bam [--write-index]` (BAM only: CRAM and SAM are refused with a hint; blocks are copied without recompression, the order of the records is not checked, and `@SQ` lists that differ are an error, where samtools joins them and the later records point at the wrong contigs, [differences.md#alignment-cat-sq](differences.md#alignment-cat-sq); `-` reads one input from stdin; `samtools cat -h` and `-b` are not available) |
| `samtools faidx REF.fa`, `samtools faidx REF.fa.gz` (BGZF) | `gloci index REF.fa` (`.fai`, and `.gzi` for BGZF, byte-identical to samtools; stricter about uneven lines and duplicate names, [differences.md#fasta-index-strictness](differences.md#fasta-index-strictness); plain gzip is refused with a `bgzip` recipe). gloci never builds the index implicitly: every command that reads a FASTA at random (`getfasta`, `nuc`, a CRAM `--reference`) needs it first |
| `samtools faidx REF.fa chr1:1-100`, `samtools faidx -r REGIONS REF.fa` | `gloci faidx REF.fa chr1:1-100` (byte-identical output, see the `samtools faidx REF REGION...` row below), or `gloci faidx REF.fa --regions-file BED`; `gloci getfasta BED --fasta REF.fa` over a BED of the regions (0-based, so `chr1:1-100` is `chr1 0 100`; the header is `chr1:0-100`, bedtools' form). `-r FILE` of `chr:START-END` text: `xargs gloci faidx ref.fa < list` ([differences.md#faidx-regions-file](differences.md#faidx-regions-file)) |
| `samtools faidx REF REGION...` (`-n N`, `-i`, `--mark-strand T`, `-c`, `-r FILE`, `-o OUT`) | `gloci faidx REF REGION... [--line-width N] [--reverse-complement] [--mark-strand rc\|no\|sign\|custom,POS,NEG] [--skip-invalid] [--regions-file BED] [-o OUT[.gz] [--write-index]]` (output is byte-identical to samtools 1.21, headers as typed, wrapping at each sequence's `.fai` width by default, `chr:-END` accepted; the FASTA needs a `.fai` (and `.gzi` for BGZF) made with `gloci index`, which gloci never builds, so `samtools faidx REF` with no region, and a missing index, exit 2 naming `gloci index REF`, see [differences.md#fasta-implicit-index](differences.md#fasta-implicit-index); `-c` is `--skip-invalid` and skips the region without writing its header, see [differences.md#faidx-missing-sequence](differences.md#faidx-missing-sequence); a region htslib guesses at (`chr1:0-5`, `chr1:0`, `chr1:+5`, an ambiguous `ctgA:5`) is a usage error, see [differences.md#faidx-region-syntax](differences.md#faidx-region-syntax); `-r` read `chr:START-END` text, `--regions-file` reads BED, see [differences.md#faidx-regions-file](differences.md#faidx-regions-file); `-X`/`--fai-idx`, `--gzi-idx` and `-f` (FASTQ) are not available) |
| `samtools dict [-H] [-A] [-l FILE] [-a STR] [-s STR] [-u STR] [-o OUT] REF` | `gloci dict REF [--no-header] [--alias] [--alt FILE] [--assembly STR] [--species STR] [--uri STR] [-o OUT]` (output is byte-identical to samtools 1.21; a repeated sequence name, [differences.md#dict-duplicate-names](differences.md#dict-duplicate-names), and an empty file, [differences.md#dict-empty-input](differences.md#dict-empty-input), exit 1) |
| `samtools flags FLAGS...` | `gloci flags [FLAGS...]` (output is byte-identical to samtools 1.21, and with no argument the table of the twelve flags, without the usage text; a decimal with a leading zero such as `010`, which samtools reads as octal 8 (`00` to `07` mean the same either way and are read), and trailing text such as `12abc` exit 2, see [differences.md#flags-leading-zero](differences.md#flags-leading-zero); an unknown name exits 2 where samtools exits 1, see [differences.md#samtools-exit-codes](differences.md#samtools-exit-codes)) |
| `samtools idxstats FILE` | `gloci idxstats FILE [--reference FA] [-o FILE]` (output is byte-identical to samtools 1.21; a BAM with a `.bai`/`.csi` is answered from the index, anything else is counted record by record with one warning; a file whose references are not grouped is counted, where samtools stops, see [differences.md#idxstats-unsorted](differences.md#idxstats-unsorted); an index without a reference's counts is scanned, where samtools prints zeros, see [differences.md#idxstats-index-without-counts](differences.md#idxstats-index-without-counts); a CRAM that needs a reference needs `--reference`, see [differences.md#cram-counts-need-reference](differences.md#cram-counts-need-reference)) |
| `samtools quickcheck [-u] [-v] FILE...` | `gloci check --quick [--unaligned] [--list-failed] FILE...` (`-u` is `--unaligned`; `--list-failed` prints the failing paths as `-v` does; exit status 1 for any failing file, where samtools exits with a bit mask, see [differences.md#samtools-exit-codes](differences.md#samtools-exit-codes); the first record must also parse, and a VCF, BCF or BGZF BED passes, see [differences.md#quickcheck](differences.md#quickcheck)) |
| `samtools head [-h N] [-n N] FILE`, `bcftools head` | `gloci head FILE [--headers N] [--records N] [-o FILE]` (any format; a BAM with no `@SQ` text gets them written, see [differences.md#head-notext-header](differences.md#head-notext-header); a VCF header is printed as stored, without the PASS line bcftools adds, see [differences.md#head-vcf-pass-line](differences.md#head-vcf-pass-line)) |
| `samtools samples [-T TAG] [-i] [-h] [-f FA]... FILE...` | `gloci samples FILE... [--tag TAG] [--show-index] [--column-names] [--match-reference FA]...` (samples in header order, see [differences.md#samples-order](differences.md#samples-order); a FASTA with no `.fai` exits 2, see [differences.md#fasta-implicit-index](differences.md#fasta-implicit-index)) |
| `samtools depth [-a\|-aa] [-J] [-q N] [-Q N] [-l N] [-g F] [-G F] [-s] [--require-flags F] [-r REG] [-b BED] [-H] [-o OUT] FILE...` | `gloci depth FILE... [--all-positions\|--all-references] [--count-deletions] [--min-baseq N] [--min-mapq N] [--min-read-length N] [--exclude-flags F] [--require-flags F] [--skip-mate-overlap] [--region REG] [--regions-file BED] [--column-names] [-o OUT]` (output is byte-identical to samtools 1.21 except where [differences.md#depth-past-reference-end](differences.md#depth-past-reference-end), [differences.md#depth-order-check](differences.md#depth-order-check), [differences.md#depth-regions](differences.md#depth-regions) and [differences.md#depth-bedcov-same-sq](differences.md#depth-bedcov-same-sq) say; `--exclude-flags` defaults to `UNMAP,SECONDARY,QCFAIL,DUP`, so `-g DUP` is `--exclude-flags UNMAP,SECONDARY,QCFAIL` and `-G X` adds X to that list; `-s` is `--skip-mate-overlap`, which counts the bases two mates share once and differs only for templates with a stale or unknown mate position, see [differences.md#depth-mate-overlap](differences.md#depth-mate-overlap)) |
| `samtools coverage [-r REG] [-l N] [-q N] [-Q N] [--rf F] [--ff F] [-o OUT] FILE...` | `gloci depth FILE... --summary [--region REG] [--min-read-length N] [--min-mapq N] [--min-baseq N] [--require-flags F] [--exclude-flags F] [-o OUT]` (one table for several files, as samtools; `-q` is `--min-mapq` and `-Q` is `--min-baseq`; `-m`, `-A`, `-w`, `-D`, `-H`, `-b` and `-d` have no counterpart; [differences.md#coverage-cram-baseq](differences.md#coverage-cram-baseq) and [differences.md#coverage-require-flags](differences.md#coverage-require-flags) say where the output differs) |
| `samtools bedcov [-Q N] [-g F] [-G F] [-j] [-d N] [-c] [-H] BED FILE...` | `gloci bedcov BED FILE... [--min-mapq N] [--exclude-flags F] [--require-flags F] [--split --split-deletions] [--bases-at-depth N] [--read-count] [--column-names] [-o OUT]` (each file needs an index; `-j` is `--split --split-deletions`; `-d N` alone drops D and N bases from both columns in samtools, see [differences.md#bedcov-depth-threshold](differences.md#bedcov-depth-threshold); a BED line on a missing chromosome exits 1, see [differences.md#bedcov-bad-line](differences.md#bedcov-bad-line); reads past the end of their reference are cut, see [differences.md#depth-past-reference-end](differences.md#depth-past-reference-end)) |
| `samtools reheader [-P] HEADER IN.bam`, `samtools reheader [-P] HEADER IN.cram` | `gloci reheader HEADER IN [-o OUT]` (HEADER is SAM header text, or a SAM, BAM or CRAM file whose header is taken, `-` is stdin; IN is a BAM or CRAM regular file, never stdin or SAM; the new `@SQ` list must keep every old reference in place with the same `LN` (names may change, new `@SQ` lines may follow), else exit 1, where samtools writes a file whose reads point at the wrong contigs, [differences.md#reheader-sq-check](differences.md#reheader-sq-check); a gloci `@PG` line is always added, so `-P` has no equivalent; `-i` is `-o IN`, which replaces the file atomically; `-c CMD` is replaced by the `view --header-only | sed | reheader -` pipe in the samtools table) |
| `samtools flagstat FILE`, `samtools flagstat -O json\|tsv FILE` | `gloci flagstat FILE [-O text\|tsv\|json] [-o OUT] [--reference FA]` (the layouts and counts are samtools 1.21's byte for byte; every record counts, unmapped and unplaced included; a CRAM that needs a reference needs `--reference`, see [differences.md#cram-counts-need-reference](differences.md#cram-counts-need-reference); `-@ N` is `-t N`) |
| `samtools view -c FILE [REGION...]`, `samtools view -c -L BED FILE` | `gloci view --count FILE [REGION...] [--regions-file BED]` (`-c` is `COL:OP` in gloci, so the flag is long; prints the number of records `view` would write, with `--min-mapq`, `--require-flags` and `--exclude-flags` applied, and counts a read once however many regions it overlaps, see [differences.md#view-overlapping-regions](differences.md#view-overlapping-regions); `-L` is `--regions-file BED` (0-based BED), read through the index when there is one, and combines with positional regions as a union, where samtools intersects them, see [differences.md#view-regions-union](differences.md#view-regions-union); a CRAM that needs a reference needs `--reference`, see [differences.md#cram-counts-need-reference](differences.md#cram-counts-need-reference)) |
| `samtools view -H FILE \| sed ... \| samtools reheader -c`, `samtools reheader -c CMD IN.bam` | `gloci view --header-only IN.bam \| sed ... \| gloci reheader - IN.bam -o OUT.bam` (`-c CMD` is not available; the pipe is its replacement) |
| `samtools markdup`, `fixmate`, `calmd`, `mpileup`, `consensus`, `stats`, `ampliconclip`, `fqidx`, `faidx --fastq` / `-f`, `-X` / `--fai-idx` / `--gzi-idx`, `depth -f` and `--incl-flags` / `-d` / `-m`, `bedcov --max-depth`, `samples -F` | not planned (for `samples -F`, repeat `--match-reference`; the other options have no counterpart) |

### bedops

| bedops | gloci |
|---|---|
| `bedops -u` / `--everything` | `gloci sort -m A B...` (for sorted inputs) or `gloci sort A B...` |
| `bedops -m` / `--merge` | `gloci merge A B...` |
| `bedops -c` / `--complement` | `gloci complement A B...` (no `-g`: gaps between features on each chromosome) |
| `bedops -d A B...` / `--difference` | `gloci merge A \| gloci subtract -a - -b B...` |
| `bedops -e [N%]` / `--element-of`, `-n` / `--not-element-of` | `gloci intersect -a A -b B... [-f 0.N]` / `gloci intersect -a A -b B... --mode none` |
| `bedops --chop N` | `gloci merge A \| gloci windows -i - -w N` |
| `--range N`, `--range L:R` | `-w/--window N` or `--left L --right R` (on `intersect`, `subtract`, `map`) or `gloci slop FILE -b N` / `--left L --right R` |
| `bedmap --echo --count/--sum/--mean/--min/--max/--median/--stdev/--variance/...` | `gloci map -a A -b B -c COL:count/sum/mean/min/max/median/stdev_pop/var_pop/...` (`--delim`, `--skip-unmatched` for `--skip-unmapped`; `gloci map` always writes A's record first, which is `--echo`; `--stdev` and `--variance` are sample statistics (n - 1): use `stdev` and `var`, not `stdev_pop` and `var_pop`; for one value bedmap prints `NAN` and gloci prints `--empty`, default `.`) |
| `bedmap --mad`, `--cv`, `--wmean` | `-c COL:mad`, `-c COL:cv`, `-c COL:wmean` (`mad` is the unscaled median absolute deviation) |
| `bedmap --tmean F1 F2` | `-c COL:tmean:F1-(1-F2)`, with the arithmetic done by you: `--tmean 0.1 0.3` trims the bottom 10% and the top 30%, so it is `tmean:0.1-0.7`. The two can differ for small sets of values, and for asymmetric trims even at eight values: for 2, 9, 4, 20, 6, 11, 3, 8, bedmap `--tmean 0.1 0.3` gives 6 and `tmean:0.1-0.7` gives 5.25; for the values 4 and 8, bedmap `--tmean 0.2 0.2` gives 6 and `tmean:0.2-0.8` gives `.` |
| `bedmap --kth F` | no exact equivalent. bedmap takes the value with fraction F of the scores below it, without interpolation; `-c COL:quantile:Q` interpolates (type 7). The two agree only at some F: for 1, 3, 5, 7, 10 at F = 0.1 bedmap gives 2 and `quantile:0.1` gives 1.8; for 4 and 8 at F = 0.25 bedmap gives 6 and gloci gives 5 |
| `bedmap --min-element`, `--max-element` | `-c COL:min_record`, `-c COL:max_record` (the whole B record with the smallest or largest value; gloci joins the record's fields with `\|`, while bedmap writes them tab-separated and prints the score with 6 decimals, so `3` is `3.000000`; see [differences.md#map-record-format](differences.md#map-record-format)) |
| `bedmap --indicator`, `--bases-uniq`, `--bases-uniq-f` | `-c 0:indicator`, `-c 0:bases_uniq`, `-c 0:frac_a` (for an element with no mapped records `bases_uniq` and `indicator` print 0 like bedmap, but `frac_a` prints `--empty` (`.`) where `--bases-uniq-f` prints `0.000000`, see [differences.md#map-empty-value](differences.md#map-empty-value)) |
| bedmap numeric output (`--prec`) | bedmap prints every number with 6 fixed decimals (`--prec 6`), such as `19.000000` for a sum, `NAN` for none; gloci prints the shortest exact form (`19`) and `--empty` for none. The values agree; compare them as numbers, or format them in a pipeline (`awk '{ printf "%.6f\n", $4 }'`) |
| `bedmap --echo-map-score` | `-c 5:collapse` (column 5 of B) |
| `bedmap --echo-ref-name`, `--echo-ref-size`, `--echo-map-range` | no counterpart (A's columns are already in the output; use `cut` or `awk`) |
| `bedmap --echo-map-id`, `--echo-map-id-uniq` | `gloci map -a A -b B -c 4:collapse`, `-c 4:distinct_sorted` (`--echo-map-id-uniq` sorts the IDs bytewise, as `distinct_sorted` does; plain `distinct` keeps first-seen order, see [differences.md#map-distinct-order](differences.md#map-distinct-order); `--echo-map` whole-record list is `0:records`, tabs shown as `\|`) |
| `bedmap --bases`, `--fraction-ref/--fraction-map/--fraction-both/--fraction-either` | `gloci map -a A -b B -c 0:bases`, with `-f` / `-F` / `-r` / `-e`; `--bp-ovr N` is `--min-bases N` |
| `closest-features` | `gloci closest -a A -b B --each-side [--distance]` (add `--ignore-overlaps` when an overlapping B should not count on both sides) |
| `sort-bed [--max-mem --tmpdir --check-sort --unique --duplicates]` | `gloci sort [--max-mem --tmp-dir --check --unique --duplicates]`; `sort-bed`'s chromosome order is the default (`--full-key` reproduces whole-line tie breaking) |
| `bedextract CHR FILE`, `bedextract FILE REGIONS.bed` | `gloci view FILE CHR`, `gloci view FILE --regions-file REGIONS.bed` (add `--assume-sorted` for binary search on an uncompressed BED file) |
| `convert2bed`, `vcf2bed [--snvs/--insertions/--deletions]`, `gff2bed`, `gtf2bed` | `gloci convert FILE [--types snv/ins/del] [--split-alleles] [--feature TYPE] [--name-attr KEY]` (rows in input order, where bedops sorts them unless `--do-not-sort`: pipe to `gloci sort --full-key` for its order, see [differences.md#convert-input-order](differences.md#convert-input-order); a one-base GFF/GTF feature keeps its attributes, without bedops' added `zero_length_insertion=True`, see [differences.md#gff2bed-zero-length](differences.md#gff2bed-zero-length)) |
| `bedops -p` / `--partition` | `gloci partition FILES --no-membership` |
| `bedops -i` / `--intersect` | `gloci partition FILES --all --merge-adjacent --no-membership` |
| `bedops -s` / `--symmdiff` | `gloci partition FILES --max-files 1 --merge-adjacent --no-membership` |
| `starch`, `unstarch`, `starchcat` | not planned (use BGZF `.bed.gz` + `gloci index`) |

### vcftools

| vcftools | gloci |
|---|---|
| `--vcf` / `--gzvcf` / `--bcf` | `gloci vcf filter FILE` (input format is auto-detected from content; `--format` overrides) |
| `--recode`, `--recode-bcf` | `gloci vcf filter FILE ... -o out.vcf.gz` / `-o out.bcf` (or stdout with `-O vcf` / `-O bcf`) |
| `--recode-INFO-all` | default in `gloci vcf filter` (INFO is kept); prune with `--keep-info KEYS` or `--drop-info KEYS` |
| `--chr` / `--not-chr` / `--from-bp` / `--to-bp` | `--chrom CHROM` / `--exclude-chrom CHROM` / `--region CHR:START-END` |
| `--bed` / `--exclude-bed` | `--regions-file FILE.bed` / `--exclude-regions-file FILE.bed` |
| `--positions` / `--exclude-positions`, `--snp` / `--exclude` | `--positions FILE` / `--exclude-positions FILE`, `--ids FILE` / `--exclude-ids FILE` |
| `--keep-only-indels` / `--remove-indels` | `--types indel` / `--exclude-types indel` (also `snv,mnv,ins,del,sv,ref,other`) |
| `--remove-filtered-all`, `--keep-filtered` / `--remove-filtered` | `--pass-only`, `--keep-filter TAGS` / `--remove-filter TAGS` |
| `--keep-INFO` / `--remove-INFO` | `--require-info FLAG` / `--exclude-info FLAG` |
| `--maf` / `--max-maf` / `--mac` / `--max-mac` / `--min-alleles` / `--max-alleles` / `--minQ` | `--min-maf` / `--max-maf` / `--min-mac` / `--max-mac` / `--min-alleles` / `--max-alleles` / `--min-qual` (plus `--biallelic`) |
| `--min-meanDP` / `--max-meanDP`, `--hwe`, `--thin`, `--phased` | `--min-mean-dp` / `--max-mean-dp`, `--hwe`, `--thin`, `--phased-only` |
| **`--max-missing X`** (X = minimum fraction *called*) | **`--max-missing (1 - X)`** (maximum fraction *missing*; e.g. vcftools `0.9` becomes gloci `0.1`; `--max-missing-count N` bounds missing genotype count). Inverted on purpose; vcftools' meaning is a known trap. |
| `--indv` / `--keep` / `--remove-indv` / `--remove` | `--samples` / `--samples-file` / `--exclude-samples` / `--exclude-samples-file` (add `--recalc-ac` to recompute `INFO/AC`, `AN`, `AF`) |
| `--minGQ` / `--minDP` / `--maxDP`, `--remove-filtered-geno-all` | `--min-gq` / `--min-dp` / `--max-dp`, `--mask-filtered-genotypes` |
| `bcftools view -i EXPR` / `-e EXPR` (subset: FORMAT fields GT, DP, GQ only; see [differences.md#expr-format-fields](differences.md#expr-format-fields)) | `gloci vcf filter FILE --include EXPR` / `--exclude EXPR` |
| `vcf-sort`, `vcf-isec`, `vcf-subset` | `gloci sort`, `gloci intersect --match alleles`, `gloci vcf filter --samples ... --recalc-ac` |
| `--freq` / `--counts` | `gloci vcf stats FILE --out-prefix P --freq` / `--counts` (write `P.frq` / `P.frq.count`, same columns) |
| `--depth` / `--site-depth` / `--site-mean-depth` / `--geno-depth` / `--site-quality` | `--depth` / `--site-depth` / `--site-mean-depth` / `--geno-depth` / `--site-quality` (`P.idepth`, `P.ldepth`, `P.ldepth.mean`, `P.gdepth`, `P.lqual`) |
| `--missing-indv` / `--missing-site` | `--missing-indv` / `--missing-site` (`P.imiss`, `P.lmiss`) |
| `--het` / `--hardy` / `--singletons` | `--het` / `--hardy` / `--singletons` (`P.het`, `P.hwe`, `P.singletons`; the singleton carrier after masking differs, [differences.md#vcf-stats-singletons-masked](differences.md#vcf-stats-singletons-masked)) |
| `--FILTER-summary` / `--TsTv-summary` / `--hist-indel-len` | `--filter-summary` / `--tstv-summary` / `--indel-hist` (`P.FILTER.summary`, `P.TsTv.summary`, `P.indel.hist`; the vcftools spellings are accepted as aliases) |
| `--SNPdensity BIN` | `--snp-density BIN` (`P.snpden`; needs sorted input, [differences.md#vcf-stats-snp-density-unsorted](differences.md#vcf-stats-snp-density-unsorted)) |
| `--minGQ` / `--minDP` / `--maxDP` with a statistic | the same `vcf stats` run takes `--min-gq` / `--min-dp` / `--max-dp` (a `.` DP or GQ is not masked, [differences.md#vcf-stats-mask-missing](differences.md#vcf-stats-mask-missing)) |
| several statistics (vcftools allows one per run), `vcf-stats` | one `gloci vcf stats` run with every flag, plus `--summary` / `--json` for an overview (`P.summary.txt`, `P.summary.json`) |
| `--get-INFO KEY` | `gloci vcf extract FILE --info KEY` (adds an ID column; `.` for missing; [differences.md#vcf-extract-info](differences.md#vcf-extract-info)) |
| `--extract-FORMAT-info KEY` | `gloci vcf extract FILE --format KEY --layout wide` (columns are `SAMPLE:KEY`; `--layout long` gives a row per sample; [differences.md#vcf-extract-format-wide](differences.md#vcf-extract-format-wide)) |
| `--012` (writes `out.012`, `.012.indv`, `.012.pos`) | `gloci vcf extract FILE --genotype-matrix 012 --out-prefix out` (same files; haploid and half-missing reference calls count 0; [differences.md#vcf-extract-012](differences.md#vcf-extract-012)) |
| `vcf-concat A.vcf.gz B.vcf.gz` | `gloci vcf concat A.vcf.gz B.vcf.gz -o OUT.vcf.gz` (identical samples; compatible header definitions; each file must start at or after the end of the one before, or it exits 1 with a `gloci sort` hint; inputs may be VCF, VCF.gz or BCF) |
| `--site-pi` / `--window-pi SIZE` / `--window-pi-step STEP` / `--TajimaD BIN` | `gloci vcf popgen FILE --out-prefix P --site-pi` / `--window-pi SIZE` / `--window-pi-step STEP` / `--tajima-d BIN` (write `P.sites.pi`, `P.windowed.pi`, `P.Tajima.D`, same columns; any of them in one run; the windowed ones need sorted input, [differences.md#popgen-sorted-input](differences.md#popgen-sorted-input); a step larger than the window is an error, [differences.md#popgen-window-step](differences.md#popgen-window-step); counts never wrap, [differences.md#popgen-integer-overflow](differences.md#popgen-integer-overflow)) |
| `--weir-fst-pop FILE` (twice or more) / `--fst-window-size SIZE` / `--fst-window-step STEP` | `gloci vcf popgen FILE --out-prefix P --fst-pop FILE` (once per population) / `--fst-window-size SIZE` / `--fst-window-step STEP` (write `P.weir.fst` and, with a window, `P.windowed.weir.fst`, same columns, plus `P.weir.fst.summary` with the two estimates vcftools logs; the populations select no samples, [differences.md#popgen-independent-options](differences.md#popgen-independent-options); one population, an empty one or a sample in two is an error, [differences.md#fst-populations](differences.md#fst-populations)) |
| `--geno-r2` / `--hap-r2` / `--interchrom-geno-r2` / `--interchrom-hap-r2` with `--ld-window N` / `--ld-window-bp BP` / `--ld-window-min N` / `--ld-window-bp-min BP` / `--min-r2 R2` | `gloci vcf popgen FILE --out-prefix P --geno-r2` / `--hap-r2` / `--interchrom-geno-r2` / `--interchrom-hap-r2` with the same window flags and `--min-r2` (write `P.geno.ld`, `P.hap.ld`, `P.interchrom.geno.ld`, `P.interchrom.hap.ld`, same columns; any of them in one run; `--hap-r2` uses only fully phased sites but filters none from the other statistics, add `--phased-only` for vcftools' `--phased`, [differences.md#popgen-independent-options](differences.md#popgen-independent-options); `--geno-r2` and `--hap-r2` need sorted input, [differences.md#popgen-sorted-input](differences.md#popgen-sorted-input); a window that admits no pair is an error, [differences.md#popgen-window-step](differences.md#popgen-window-step); the sites held stay within `--max-mem SIZE` (default 2G), with fewer than two sites the files hold their headers, and an excluded sample's phase is not read, [differences.md#ld-rules](differences.md#ld-rules)) |
| `--relatedness` / `--relatedness2` | `gloci vcf popgen FILE --out-prefix P --relatedness` / `--relatedness2` (write `P.relatedness`, `P.relatedness2`, same columns; both in one run, beside the other statistics; a masked or half-called genotype is missing in both, [differences.md#relatedness-genotypes](differences.md#relatedness-genotypes)) |
| `--diff B`, `--gzdiff B`, `--diff-bcf B` with `--diff-site-discordance` / `--diff-indv-discordance` / `--diff-discordance-matrix`, `--diff-indv-map FILE` | `gloci vcf compare A B --out-prefix P` with `--site-discordance` / `--sample-discordance` / `--discordance-matrix`, `--sample-map FILE` (B is the second positional file, any of VCF, VCF.gz or BCF, so `--gzdiff` and `--diff-bcf` need no flag; `P.diff.sites`, `P.diff.indv`, `P.diff.discordance_matrix`, same columns; genotype concordance only. Sites are paired by `--match alleles` by default (same CHROM and POS and a shared ALT after trimming); add `--match position` for vcftools' site pairing, [differences.md#compare-allele-matching](differences.md#compare-allele-matching). The matrix counts biallelic pairs of the same trimmed allele under `--match alleles` (so `AT>A` and `ATG>AG` count together) and of the same REF and ALT text under `--match position`. A half-missing `./1` is missing, [differences.md#compare-half-missing](differences.md#compare-half-missing). `--thin` thins each file on its own, whatever the headers list and with `-g`) |
| `--diff B --diff-site` | sites in both files: `gloci intersect -a A -b B --match alleles`; in A only: `gloci subtract -a A -b B --match alleles` (swap A and B for B only). `--match alleles` is the default; `--match position` pairs by POS and REF like vcftools. `vcf compare` does not write these site sets |
| `--diff B --diff-indv` | the samples of each file are the `#CHROM` line: `gloci view --header-only FILE \| tail -n1 \| cut -f10-`; `vcf compare` writes every sample of the union in `P.diff.indv` with a discordance row, not the in-one-file lists |
| `vcf-compare` (text report) | `gloci vcf compare A B --out-prefix P` gives the discordance tables; the text report itself is not planned |
| `vcf-merge` | `gloci vcf merge A B ...` (see the bcftools `merge` row above) |
| `--plink*`, `--IMPUTE`, `--BEAGLE-*`, `--ldhat*` | not planned |

### Filter expressions (bcftools `-i` / `-e`)

`bcftools view -i EXPR` is `gloci vcf filter FILE --include EXPR`, and
`bcftools view -e EXPR` is `--exclude EXPR`. The same flags work on
`vcf extract` and `vcf stats`, so a selection can be applied while extracting or
counting. The grammar is a subset of bcftools 1.21: what it accepts means what
bcftools means (the differences are in
[differences.md#expr-subset](differences.md#expr-subset)), and everything else is a usage error (exit 2) that
names the token and its column. Quote the expression in the shell.

| bcftools expression | gloci expression (same text unless shown) |
|---|---|
| `QUAL>30`, `POS>=1000`, `CHROM="1"`, `REF="A"`, `ALT="C"` | unchanged |
| `FILTER="PASS"`, `FILTER!="PASS"`, `FILTER~"Low"` | unchanged (an `=` or `!=` tag must be declared in the header, and one tag per comparison: `FILTER="A,B"` is an error, use `FILTER="A" \|\| FILTER="B"`. For `~`, bcftools also errors on a pattern that is not a declared tag (`The filter "Low" not present in the VCF header`); gloci accepts any regular expression) |
| `ID~"^rs"`, `ID="."`, `ID!="."` | unchanged (regular expressions are Rust `regex` syntax, not POSIX) |
| `TYPE="snp"`, `TYPE="indel"` | unchanged (`snp`, `mnp`, `indel`, `other`, `ref`) |
| `INFO/DP>10 && INFO/AF<0.3`, `INFO/DP>=10 \|\| QUAL<20` | unchanged (`&&`, `\|\|`, parentheses, and `!` before a group or comparison, which bcftools rejects) |
| `INFO/AF[0]>0.2`, `INFO/CSQ~"missense"` | `[i]` is a 0-based element index, as in bcftools. An index past the end is missing in gloci, so the comparison is false, where bcftools 1.21 can select the record: [differences.md#expr-index-past-end](differences.md#expr-index-past-end) |
| `INFO/DB` (a flag), `INFO/DB=1` | unchanged (a bare flag is its presence; bcftools selects nothing for the bare form) |
| `QUAL=.`, `INFO/AF=.`, `ID=.` | unchanged (`.` is the missing value; in bcftools only `ID="."` works) |
| `FMT/DP>6`, `FMT/GQ>=40`, `GT="het"`, `GT="alt"`, `GT="mis"` | unchanged: a bare or `FMT/` comparison holds when any sample satisfies it (the keywords are `ref`, `alt`, `het`, `hom`, `mis`, and the literal genotypes such as `"1/1"`) |
| `FMT/GQ[1]>=40`, `GT[0]="het"` | unchanged (sample subscripts are 0-based over the whole header, even with `--samples`) |
| `N_PASS(GT="het")>=2`, `N_PASS(FMT/DP>5 && FMT/GQ>30)>=1` | unchanged (`N_PASS` is the only function; its test must be a sample-wise comparison joined with `&&` or `\|\|`) |
| `F_MISSING<0.1`, `N_MISSING=0`, `MAF>0.05`, `AC>1`, `AN=6`, `N_ALT>=1`, `N_SAMPLES>2` | unchanged (`MAF` is bcftools' first-ALT MAF, not `--min-maf`; with `--samples` these see only the selected samples) |
| `DP>10` when `DP` is both INFO and FORMAT | `INFO/DP>10` or `FMT/DP>10` (bcftools refuses the bare form too) |
| `bcftools view -f PASS` | `--keep-filter PASS`. Not `--pass-only`: that keeps FILTER `PASS` and `.` (vcftools `--remove-filtered-all`), where `bcftools view -f PASS` drops `.` |
| `-v snps`, `-m2 -M2`, `-s S`, `-r R`, `-R BED` | `--types snv`, `--biallelic`, `--samples S`, `--region R`, `--regions-file BED` |
| `bcftools filter -i EXPR -s TAG` (soft filter), `-g`, `-G` (SnpGap, IndelGap), `-S .` | not available: gloci removes records, it never sets FILTER; `--pass-only`, `--keep-filter` and `--remove-filter` select on an existing FILTER |

Not supported, each a usage error (exit 2) naming the construct:

- Functions other than `N_PASS`: `ABS`, `STRLEN`, `MIN`, `MAX`, `AVG`, `SUM`, `COUNT`, `SMPL_*`, `BINOM`, `PHRED` and the rest. Compute the value outside, or filter with bcftools.
- FORMAT fields other than `GT`, `DP` and `GQ`: `FMT/AD`, `FMT/FT`, `FMT/PL` and so on. `vcf filter` has `--mask-filtered-genotypes` for FT.
- Comparing two fields (`QUAL>POS`, `INFO/DP>FMT/DP`) or arithmetic (`QUAL*2>10`): a field is compared with a number, a string or `.`.
- The sample-wise `&` and `|`: write `&&` and `||`, which combine whole-record conditions; use `N_PASS` for several conditions on one sample.
- Index ranges and wildcards (`INFO/AF[0-1]`, `AF[*]`), `@file` lists, `GT` subscripts past the sample count, and regular expression features Rust's `regex` lacks (back-references, look-around), which are usage errors; POSIX classes such as `[[:alpha:]]` work.
- Words for the operators (`and`, `or`) and tags the header does not define (`MAC`, `ILEN`, any undeclared INFO key).

Before relying on a new expression, run it once on a small file with
`gloci vcf filter small.vcf --include EXPR` and compare with
`bcftools view -i EXPR small.vcf`.

## BAM, CRAM and FASTA behaviour

What to expect when a pipeline moves from bedtools and samtools to gloci for
alignments and sequence files. Each point links to the pinned difference where
there is one.

- **Unmapped reads.** `intersect`, `coverage`, `genomecov` and `convert` treat
  an alignment as an interval, and an unmapped read (flag `0x4`, placed or not)
  has none: they skip it and count it at `-v`. bedtools `intersect -abam -v`
  writes them ([differences.md#unmapped-reads](differences.md#unmapped-reads)).
  The record-preserving commands, `view`, `sort`, `sort -m` and `cat`, keep
  every record, unmapped reads included. A read with no reference operations
  is one base, and a read that runs past its reference's end is cut there
  ([differences.md#bamtobed-no-reference-ops](differences.md#bamtobed-no-reference-ops),
  [differences.md#bamtobed-past-reference-end](differences.md#bamtobed-past-reference-end)).
- **Filters.** `--min-mapq N`, `--require-flags FLAGS` and `--exclude-flags
  FLAGS` (a number such as `0x704`, or comma-separated samtools names) apply to
  every alignment input of `intersect`, `coverage`, `genomecov`, `convert` and
  `view` before anything else. The defaults filter nothing, as in bedtools. For
  the semantics of `samtools depth` and mosdepth (no unmapped, secondary,
  QC-fail or duplicate reads) pass `--exclude-flags 0x704`. Given to a run
  with no BAM, CRAM or SAM input, they are not used, and gloci says so once.
- **`--split` and `--split-deletions`.** `--split` cuts a read at CIGAR `N`
  only, and `D` stays inside a block (bedtools `-split`). `--split-deletions`
  also cuts at `D` (bedtools `bamtobed -splitD`). `genomecov` counts deleted
  bases by default, where bedtools `-ibam` cuts at them unless `-ignoreD`, so
  `bedtools genomecov -ibam -split` is `--split-deletions` and
  `-split -ignoreD` is `--split`
  ([differences.md#bam-genomecov-deletions](differences.md#bam-genomecov-deletions)).
- **Sorted input, and `--sort`.** Streaming commands need the alignments in
  the order of the header's `@SQ` lines, which is checked as they are read;
  the error names `gloci sort x.bam -o x.sorted.bam`. `--sort` works with BAM,
  CRAM and SAM inputs as with BED: the reads are sorted first, bounded by
  `--max-mem` (a SAM input is spilled as SAM text and keeps its layout; BAM and
  CRAM spill a level-1 BAM), and `--unsorted` holds a B alignment file in
  memory.
- **B BEDs follow the BAM's `@SQ` order.** The header's order is the
  chromosome order, and it beats `--chrom-order`. Sort a BED to it with
  `gloci sort b.bed -g x.bam` (a BAM, SAM or CRAM header is a valid `-g`), or
  hand the BED to the command with `--sort`.
- **No implicit indexes.** gloci never builds an index to answer a query. A
  region query or a per-chromosome run (`-t` above 1 on `intersect`,
  `coverage` or `genomecov`) uses `x.bam.bai`, `x.bai`, `x.bam.csi`,
  `x.cram.crai` or `x.crai` when it exists, and otherwise scans the file (`view`
  warns, and samtools exits 1
  there, [differences.md#bam-region-without-index](differences.md#bam-region-without-index)).
  Build one with `gloci index x.bam` or `--write-index`. An index older than its
  file is warned about. FASTA is the same: `.fai` (and `.gzi` for BGZF) come
  from `gloci index ref.fa`, where bedtools builds a `.fai` beside the FASTA on
  its own. SAM text has no index.
- **CRAM.** The reference is `--reference FASTA` only (bedtools `--cram-ref`),
  with no `REF_PATH`, cache or download, and gloci never writes a CRAM (`sort`
  and `view` write BAM or SAM). `gloci index x.cram` needs no reference. See the
  CRAM entry under Other behaviour changes, and
  [differences.md#cram-limits](differences.md#cram-limits) for what gloci
  refuses to decode.
- **CRAM reference memory.** The decoder threads share one cache of reference
  sequences: each sequence is read from the FASTA once and held once, however
  many threads use it. Sorted input holds at most about two chromosomes at a
  time (chr1 is about 250 MB, and about twice that while it loads). With a
  `.crai` and `-t` above 1, a CRAM is read per chromosome, one chromosome per
  worker, so at most `-t` chromosomes at once (workers are capped at 8, so
  8 x 250 MB for chr1-sized chromosomes). Unsorted or name-sorted input keeps
  every chromosome its slices touch, as samtools does. A CRAM written without a
  reference holds a zero-filled stand-in as long as each reference it touches.
  Up to 2 x `-t` containers are in flight, each decoded on a worker with a 16
  MiB stack.
- **Alignment output.** A command that writes reads (`intersect --mode any` or
  `none`, `view`, `sort`, `sort -m`) writes the input's format: BAM for BAM
  and CRAM, SAM for SAM, and `-o`/`-O` override. samtools always writes BAM
  from `sort` and `merge`. BED output of reads is a usage error that suggests
  `-O bam | gloci convert -`. Every output gets an `@PG` line.
  Where an alignment is printed as columns (`--mode pairs`, `count`, B columns,
  `coverage`) it is BED6, not bedtools' BED12 with `0,0,0`
  ([differences.md#bam-as-bed-columns](differences.md#bam-as-bed-columns)).
- **Indexing, merging and joining.** `gloci sort -m` is `samtools merge`, and
  `gloci merge` stays the interval merge of bedtools and bedops. `-n` is
  `sort -n` (by read name) on `sort` and `windows -n` (a window count) on
  `windows`: the same letter, different commands. `gloci cat` joins BAM files
  that are already in order.
- **FASTA.** `getfasta` and `nuc` fail on an interval outside the FASTA
  (`--skip-invalid` skips and counts, as bedtools does), `maskfasta` streams
  the FASTA and keeps its header and line layout, and a mask past the end of its
  sequence is an error (a mask on a chromosome the FASTA lacks is a warning)
  ([differences.md#getfasta-outside-fasta](differences.md#getfasta-outside-fasta),
  [differences.md#maskfasta-layout](differences.md#maskfasta-layout),
  [differences.md#maskfasta-past-end](differences.md#maskfasta-past-end)).

## samtools behaviour

What to expect when `samtools` commands move to gloci. Each point links to the
pinned difference where there is one.

- **Filter defaults.** `depth` and `bedcov` exclude `UNMAP,SECONDARY,QCFAIL,DUP`
  by default, as samtools does, where `coverage` and `genomecov` filter nothing
  by default, as bedtools does (D102). `--exclude-flags 0` turns the default off
  (unmapped reads still never count in pileups); samtools `-g DUP` is
  `--exclude-flags UNMAP,SECONDARY,QCFAIL`.
- **Zero rows inside spliced reads.** `depth` prints a 0 row at each position of
  an intron, a deletion or a filtered base inside a read's span, as samtools
  does; a position outside every span prints no row without `--all-positions`.
  [differences.md#depth-regions](differences.md#depth-regions)
- **`view --count` counts a read once.** However many regions it overlaps;
  samtools repeats it per positional region.
  [differences.md#view-overlapping-regions](differences.md#view-overlapping-regions)
- **CRAM needs `--reference`** for every command that decodes records (`view`,
  `depth`, `bedcov`, `head`, and the `flagstat` and `idxstats` scans) when a
  container needs reference bases; `samples`, `reheader` and `check --quick`
  never do. samtools reads `flagstat`, `idxstats` and `view -c` without one.
  [differences.md#cram-counts-need-reference](differences.md#cram-counts-need-reference)
- **Exit codes.** 0 for success, 1 for data errors and failed checks, 2 for usage
  errors. samtools exits 1 for usage errors, 2 for `bedcov` errors and a bit mask
  from `quickcheck`.
  [differences.md#samtools-exit-codes](differences.md#samtools-exit-codes)
- **No implicit FASTA indexes.** `faidx`, `samples --match-reference` and a CRAM
  `--reference` need the `.fai` (and `.gzi`) from `gloci index ref.fa` first,
  where samtools builds it.
  [differences.md#fasta-implicit-index](differences.md#fasta-implicit-index)
- **Long flags.** samtools' short flags are not copied, since `-c -n -i -s -a -b
  -r -f` mean other things in gloci; every samtools flag maps to a long flag in
  its row, and `-@ N` is `-t N`.

## Other behaviour changes

- **VCF allele footprints:** VCF/BCF records are footprinted by the reference
  bases their ALT alleles change by default (`--vcf-span allele`), not by REF.
  This applies to `view`, `intersect`, `subtract` and `convert`. Pass
  `--vcf-span ref` for `[POS-1, POS-1+len(REF))`.
  [differences.md#default-vcf-span-allele](differences.md#default-vcf-span-allele),
  [differences.md#vcf-allele-span](differences.md#vcf-allele-span)
- **Points and region/interval edges:** A zero-length BED record or a pure VCF
  insertion is a point between bases and overlaps a region or interval only when
  strictly inside it. Pass `--inclusive-points` to include points on either edge
  (`--vcf-span ref --inclusive-points` reproduces tabix on VCF).
  [differences.md#view-region-footprints](differences.md#view-region-footprints),
  [differences.md#zero-length-overlap](differences.md#zero-length-overlap)
- **Zero-length BED records in `sort`:** `gloci sort` accepts zero-length BED
  records as points, whereas `sort-bed` rejects them.
  [differences.md#zero-length-bed-records](differences.md#zero-length-bed-records)
- **`sort --check` checks `(chrom, start)` only:** Records with the same start
  may have any end order (`sort-bed --check-sort` also requires sorted ends).
  [differences.md#sort-check-chrom-start](differences.md#sort-check-chrom-start)
- **`##glociCommand` header line:** VCF and BCF output from `sort`, `view`,
  `intersect`, `subtract` and `vcf filter` appends one `##glociCommand`
  provenance header line.
  [differences.md#glocicommand-header](differences.md#glocicommand-header)
- **BGZF end-of-file marker:** A seekable BGZF file without its 28-byte EOF
  block is an error, not a warning.
  [differences.md#bgzf-eof-marker](differences.md#bgzf-eof-marker)
- **`view --assume-sorted`:** Uncompressed BED region queries without an index
  scan the file and check sort order by default; `--assume-sorted` binary-searches
  like `bedextract` while also finding nested or long spanning features.
  [differences.md#view-assume-sorted](differences.md#view-assume-sorted)
- **Symbolic SVs with only `SVLEN`:** gloci spans symbolic SVs without `INFO/END`
  over `POS + |SVLEN|` when indexing and querying; rebuild htslib-built indexes
  with `gloci index` to query inside them.
  [differences.md#htslib-index-svlen](differences.md#htslib-index-svlen)
- **Symbolic insertions with `END` under `--vcf-span ref`:** an `<INS>` (or
  `<INS:ME:...>`) with INFO/END spans to END, as tabix has it; `bedtools
  intersect` spans it by its REF base alone.
  [differences.md#vcf-span-ins-end](differences.md#vcf-span-ins-end)
- **VCF header validation on BCF conversion:** Converting between VCF and BCF
  refuses malformed headers that `noodles` rejects (such as duplicate INFO IDs).
  [differences.md#vcf-header-conversion](differences.md#vcf-header-conversion)
- **Clipping GFF/GTF by a point:** `intersect --mode clip` skips the empty
  piece a point leaves in a GFF/GTF feature (with one warning) where bedtools
  writes one base around the point.
  [differences.md#clip-gff-points](differences.md#clip-gff-points)
- **All-missing sites:** `vcf filter --min-maf` drops a site with no called
  genotype (vcftools `--maf` keeps it; bcftools drops it).
  [differences.md#all-missing-site-stats](differences.md#all-missing-site-stats)
- **Undeclared VCF contigs:** Records on a contig missing from the VCF header
  are accepted with a warning by `view` and `index`, as htslib does, but `sort`
  and BCF output refuse them (add `##contig` or pass `-g`).
  [differences.md#undeclared-contig-sort](differences.md#undeclared-contig-sort)
- **BCF contig dictionary order in `index`:** `gloci index` refuses a BCF whose
  contigs are not in header dictionary order. Run `gloci sort` first.
  [differences.md#bcf-index-contig-order](differences.md#bcf-index-contig-order)
- **FASTA indexing is stricter than `samtools faidx`:** `gloci index ref.fa` exits 1
  for a duplicate sequence name, an empty name, or a space or non-ASCII byte in a
  sequence line, where samtools indexes the file. Valid FASTA gets a
  byte-identical `.fai`/`.gzi`.
  [differences.md#fasta-index-strictness](differences.md#fasta-index-strictness)
- **`getfasta --split` checks the BED12 blocks:** a record with out-of-order or
  overlapping blocks, or blocks that do not span it, makes `gloci getfasta`
  exit 1 where bedtools writes the blocks unchecked; `--skip-invalid` skips them.
  [differences.md#getfasta-split-blocks](differences.md#getfasta-split-blocks)
- **`getfasta` and `nuc` outside the FASTA:** an interval on a chromosome the FASTA lacks,
  or ending past its end, makes `gloci getfasta` and `gloci nuc` exit 1 naming the interval,
  where bedtools prints a message, skips it and exits 0. `--skip-invalid`
  skips and counts them.
  [differences.md#getfasta-outside-fasta](differences.md#getfasta-outside-fasta)
- **`faidx` does not build a missing index:** `gloci faidx ref.fa REGION` exits 2
  naming `gloci index ref.fa` when `ref.fa.fai` is absent, and `samtools faidx
  ref.fa` with no region (which only builds the index) is `gloci index ref.fa`.
  [differences.md#fasta-implicit-index](differences.md#fasta-implicit-index)
- **`faidx` regions:** an unknown sequence exits 1 before any output, and
  `--skip-invalid` skips the region without the header `samtools faidx -c`
  still writes; a region htslib guesses at (`chr1:0`, `chr1:+5`, `chr1:0-5`, an ambiguous contig name) is exit 2 before any output; the
  regions file is BED, not `chr:START-END` text.
  [differences.md#faidx-missing-sequence](differences.md#faidx-missing-sequence)
  [differences.md#faidx-region-syntax](differences.md#faidx-region-syntax)
  [differences.md#faidx-regions-file](differences.md#faidx-regions-file)
- **`dict` (replaces `samtools dict`):** `gloci dict REF [--no-header] [--alias]
  [--alt FILE] [--assembly STR] [--species STR] [--uri STR] [-o OUT]` writes the
  same dictionary byte for byte (`-H` is `--no-header`, `-A` is `--alias`, `-l`
  is `--alt`, `-a` is `--assembly`, `-s` is `--species`, `-u` is `--uri`; `-o` is
  the same). Two inputs differ: a
  repeated sequence name exits 1 where samtools writes both, and an empty file
  exits 1 where samtools writes only `@HD`.
  [differences.md#dict-duplicate-names](differences.md#dict-duplicate-names)
  [differences.md#dict-empty-input](differences.md#dict-empty-input)
- **`flags` (replaces `samtools flags`):** `gloci flags [FLAGS...]` prints the same
  lines. A decimal with a leading zero (`010`) and trailing text (`12abc`) exit 2
  where samtools reads octal 8 and 12 (`00` to `07` mean the same either way and
  are read), and a usage error exits 2 where samtools exits 1.
  [differences.md#flags-leading-zero](differences.md#flags-leading-zero)
  [differences.md#samtools-exit-codes](differences.md#samtools-exit-codes)
- **`depth` (replaces `samtools depth`):** `gloci depth FILE...` prints the same
  rows for BAM, CRAM and SAM (`-a` is `--all-positions`, `-aa` is
  `--all-references`, `-J` is `--count-deletions`, `-q` is `--min-baseq`, `-Q` is
  `--min-mapq`, `-l` is `--min-read-length`, `-r` is `--region`, `-b` is
  `--regions-file`, `-H` is `--column-names`). Four things differ: a read past
  the end of its reference is cut there where samtools prints rows beyond it; an
  input out of order is an error where samtools counts it; `-a` with regions
  prints the regions on references with no read, as `-r` does; and inputs with
  different `@SQ` lists are refused where samtools reads their ids as equal.
  [differences.md#depth-past-reference-end](differences.md#depth-past-reference-end)
  [differences.md#depth-order-check](differences.md#depth-order-check)
  [differences.md#depth-regions](differences.md#depth-regions)
  [differences.md#depth-bedcov-same-sq](differences.md#depth-bedcov-same-sq)
  A read with SEQ `*` passes every `--min-baseq`; samtools reads past the record there.
  [differences.md#depth-seq-star-baseq](differences.md#depth-seq-star-baseq)
  `depth --skip-mate-overlap` (samtools `-s`) does not clip by a stale entry.
  [differences.md#depth-mate-overlap](differences.md#depth-mate-overlap)
- **`reheader` (replaces `samtools reheader`):** `gloci reheader HEADER IN [-o OUT]`
  writes the same BAM or CRAM as samtools, except that a new `@SQ` list that
  moves, drops or resizes a reference exits 1 where samtools writes the file.
  [differences.md#reheader-sq-check](differences.md#reheader-sq-check)
- **`view --count` and `--regions-file` (replace `samtools view -c` and `-L`):** a read that
  overlaps several regions counts and is written once, and positional regions and
  `--regions-file` combine as a union.
  [differences.md#view-overlapping-regions](differences.md#view-overlapping-regions)
  [differences.md#view-regions-union](differences.md#view-regions-union)
- **TSV indexing without `-e`:** `gloci index -s -b` without `-e` sets the end
  column equal to the begin column (1-base feature) instead of tabix's GFF
  default of column 5.
  [differences.md#tsv-index-no-end-column](differences.md#tsv-index-no-end-column)
- **Per-record overlap fractions:** `-f`, `-F`, `-r`, `-e` and `--min-bases` are
  evaluated per B record (like bedtools), not summed across B records (like
  `bedops --element-of`; run `gloci merge` on B first to sum).
  [differences.md#per-record-overlap-fraction](differences.md#per-record-overlap-fraction)
- **Trimmed allele matching (`--match alleles`):** `intersect --match alleles` and `subtract --match alleles`
  compare ALT alleles after trimming shared REF/ALT prefix and suffix and
  matches when any ALT is shared.
  [differences.md#trimmed-allele-match](differences.md#trimmed-allele-match)
- **`merge -s` keeps unknown-strand records:** Records with `.` or no strand are
  clustered together as a third group (with a warning) instead of being silently
  dropped as in `bedtools merge -s`.
  [differences.md#merge-unknown-strand](differences.md#merge-unknown-strand)
- **`subtract --sum-overlaps` threshold:** `--sum-overlaps -f F` drops A when
  the summed overlap reaches at least `F` (`>= F`, whereas `bedtools subtract -N`
  requires `> F`), and works without `-f`.
  [differences.md#sum-overlaps-threshold](differences.md#sum-overlaps-threshold)
- **`closest --each-side` overlaps and placeholders:** `--each-side` counts an
  overlapping B on both sides at distance 0 unless `--ignore-overlaps` is
  passed, and writes standard placeholders (`. -1 -1 ...`) instead of `NA`.
  [differences.md#closest-each-side](differences.md#closest-each-side)
- **`closest -k N --ties first|last`:** Cuts ties only at the N-th distance so
  at most N rows are written, rather than collapsing ties at every closer
  distance.
  [differences.md#closest-ties-k](differences.md#closest-ties-k)
- **Sample `stdev` and `var` in `map` and `merge`:** `stdev` and `var` divide by
  `n - 1` (returning `--empty` for one value, matching `bedtools sstdev`); use
  `stdev_pop` and `var_pop` for population statistics dividing by `n` (matching
  `bedtools map -o stdev`). `bedmap --stdev` and `--variance` are sample
  statistics too, so they are `stdev` and `var`.
  [differences.md#map-stdev-sample](differences.md#map-stdev-sample)
- **First-seen order for `distinct`:** `map` and `merge` `distinct` preserves
  first-seen order; use `distinct_sorted` for `bedtools`' bytewise-sorted list.
  [differences.md#map-distinct-order](differences.md#map-distinct-order)
- **Missing values in `mode` and `antimode`:** `map` skips `.` and empty
  values; `bedtools map` counts `.` as a value, so it can print `.` while real
  values exist.
  [differences.md#map-mode-missing](differences.md#map-mode-missing)
- **`groupby --header`:** writes the key names and the `-c` labels as its first
  line and lets `-g` and `-c` use names, where `bedtools groupby -header` copies
  the input header unchanged.
  [differences.md#groupby-header](differences.md#groupby-header)
- **`cluster -s` order and ids:** `gloci cluster -s` keeps the input order and numbers
  clusters as their first record is read; `bedtools cluster -s` regroups by strand and
  drops `.`-strand records.
  [differences.md#cluster-strand](differences.md#cluster-strand)
- **`spacing` on points:** a zero-length record is a point, where `bedtools spacing`
  widens it by one base on each side.
  [differences.md#spacing-zero-length](differences.md#spacing-zero-length)
- **A missing column in `groupby`, `map` and `merge`:** a line without a `-g` or `-c`
  column is an error naming the line (exit 1), also for `N:count`; `bedtools groupby`,
  `bedtools map` and `bedtools merge -c` stop with a message and exit 0 when every line
  lacks the column (they fail on ragged input, as gloci does).
  [differences.md#groupby-missing-column](differences.md#groupby-missing-column)
- **Non-numeric values in `map` and `groupby`:** a value that is not a number in
  a numeric op is an error naming the line; bedtools warns, skips it and exits 0.
  [differences.md#agg-non-numeric](differences.md#agg-non-numeric)
- **Empty placeholder in `map`:** An aggregation with no values prints `--empty`
  (default `.`), not `bedmap`'s `NAN`.
  [differences.md#map-empty-value](differences.md#map-empty-value)
- **`coverage` of a point A:** A zero-length A covers 0 bases (covered bases `0`,
  length `0`, fraction `.`, `--mean` `.`, and no `--hist` or `--per-base` rows).
  [differences.md#coverage-point](differences.md#coverage-point)
- **`coverage --split` measures A by its blocks:** With `--split`, A's length and
  `--per-base` offsets run over A's BED12 blocks, excluding introns.
  [differences.md#coverage-split-length](differences.md#coverage-split-length)
- **`coverage --split` counts a spliced B once:** A BED12 or BAM/SAM B that meets A through
  several blocks is one hit, not one per block.
  [differences.md#coverage-split-count](differences.md#coverage-split-count)
- **`-f`/`-F` under `--split` judge each hit alone:** bedtools decides them for all hits of
  an A together when B is a BAM or BED12 file; gloci uses each B record's own blocks.
  [differences.md#split-fraction](differences.md#split-fraction)
- **`complement` ignores points:** Zero-length records cover no bases and do not
  split gaps in `complement`.
  [differences.md#complement-points](differences.md#complement-points)
- **`complement -g` rejects features past chromosome end:** A feature extending
  past its chromosome length in `-g` exits 1 instead of warning.
  [differences.md#complement-past-end](differences.md#complement-past-end)
- **`partition` cuts at every boundary:** without `--merge-adjacent`,
  neighbouring segments with equal membership stay apart (as `bedops -p`, unlike
  `multiinter`).
  [differences.md#partition-boundaries](differences.md#partition-boundaries)
- **`partition` ignores points:** a zero-length record adds no segment and no
  boundary (`multiinter` widens it, `bedops -p` writes a zero-length segment).
  [differences.md#partition-points](differences.md#partition-points)
- **`partition --bedgraph --empty -g` covers every chromosome:** also those
  without features (`unionbedg -empty -g` omits them).
  [differences.md#partition-empty-chromosomes](differences.md#partition-empty-chromosomes)
- **`genomecov` prints chromosomes in genome order:** a chromosome without
  features comes in its place, not after the others.
  [differences.md#genomecov-chrom-order](differences.md#genomecov-chrom-order)
- **`genomecov` of a BAM counts deleted bases:** `bedtools genomecov -ibam` cuts reads at
  `D` unless `-ignoreD`; gloci covers the whole reference span. `--split-deletions` is
  bedtools `-split`.
  [differences.md#bam-genomecov-deletions](differences.md#bam-genomecov-deletions)
- **`genomecov` of a BAM uses `-g`, and counts a read with no reference operations as one
  base:** bedtools ignores `-g` with `-ibam` and counts no base for such a read.
  [differences.md#genomecov-bam-genome-arg](differences.md#genomecov-bam-genome-arg),
  [differences.md#bamtobed-no-reference-ops](differences.md#bamtobed-no-reference-ops)
- **`genomecov` ignores points, treats an unknown strand as plus and refuses
  features past the chromosome end:**
  [points](differences.md#genomecov-points),
  [strand](differences.md#genomecov-unknown-strand),
  [past the end](differences.md#genomecov-past-end)
- **`genomecov` counts a genome of more than 2^32 bases:** the `genome` rows
  of `bedtools genomecov` wrap at 2^32 bases; gloci's do not.
  [differences.md#genomecov-genome-total](differences.md#genomecov-genome-total)
- **`genomecov --max-depth` is for the histogram only:** `bedtools -max` is
  silently ignored with `-bg`, `-bga` and `-d`; gloci refuses it.
  [differences.md#genomecov-max-depth](differences.md#genomecov-max-depth)
- **One-sided `slop` and `flank`:** `--left N` or `--right N` may be given alone
  without requiring both.
  [differences.md#slop-flank-one-side](differences.md#slop-flank-one-side)
- **`flank` of a point:** Flanking `[p, p)` writes `[p - L, p)` and `[p, p + R)`
  touching the point, without a 1-base gap.
  [differences.md#flank-points](differences.md#flank-points)
- **`slop` and `flank` past chromosome end:** A feature starting at or beyond its
  chromosome length in `-g` exits 1.
  [differences.md#slop-flank-past-end](differences.md#slop-flank-past-end)
- **`shift` of a point past an edge:** A zero-length feature shifted wholly
  beyond 0 or the chromosome end stays zero-length at the edge; bedtools writes
  an inverted interval (such as `1 0`).
  [differences.md#shift-past-edge](differences.md#shift-past-edge)
- **`shift --pct`:** Exact decimal arithmetic instead of single-precision float
  (so `-0.1` of 10 bases is -1, not bedtools' -2); a zero-length feature does not
  move.
  [differences.md#pct-exact-decimal](differences.md#pct-exact-decimal)
- **`flank -s` on unknown strand:** Treats `.` or missing strand as `+` (matching
  `slop -s`) instead of dropping the record.
  [differences.md#flank-unknown-strand](differences.md#flank-unknown-strand)
- **`windows -i` on a point:** A zero-length interval produces no windows.
  [differences.md#windows-points](differences.md#windows-points)
- **`windows -g` on a zero-length chromosome:** Skips a 0-length contig with a
  warning instead of exiting 1.
  [differences.md#windows-zero-length-chromosome](differences.md#windows-zero-length-chromosome)
- **`convert` writes one row per VCF record by default:** Multi-allelic records
  stay one row unless `--split-alleles` is passed (unlike `vcf2bed`).
  [differences.md#convert-multiallelic-rows](differences.md#convert-multiallelic-rows)
- **`convert` variant spans and `--types`:** Variant spans and `--types` follow
  trimmed allele footprints (or full REF/END/SVLEN spans under `--vcf-span ref`)
  rather than `vcf2bed`'s string-length heuristics.
  [differences.md#vcf2bed-variant-spans](differences.md#vcf2bed-variant-spans)
- **`convert` names GFF3 rows by `ID`:** Column 4 of GFF3 output defaults to the
  `ID` attribute (`--name-attr` overrides), whereas `gff2bed` writes `.`.
  [differences.md#gff2bed-name](differences.md#gff2bed-name)
- **`convert --blocks` validates BED12:** Records with fewer than 12 columns or
  invalid block coordinates exit 1 instead of passing through unchanged.
  [differences.md#bed12-short-line](differences.md#bed12-short-line)
- **`convert` reads BLAT PSL headers and unterminated last lines:** `gloci
  convert x.psl` skips the five-line header and keeps a last line without a
  newline; `convert2bed --input=psl` fails on the header and drops that line.
  [differences.md#psl-rmsk-bedops-quirks](differences.md#psl-rmsk-bedops-quirks)
- **`maskfasta` keeps the FASTA as it is:** `gloci maskfasta` keeps header
  descriptions, line widths, CRLF and a missing final newline; `bedtools
  maskfasta` drops descriptions (unless `-fullHeader`), re-wraps to the first
  line's width, and masks nothing in a CRLF file; it also adds a final line
  break that was missing and drops sequences with no bases.
  [differences.md#maskfasta-header](differences.md#maskfasta-header)
  [differences.md#maskfasta-layout](differences.md#maskfasta-layout)
- **`maskfasta` refuses a mask past the sequence end:** exit 1, where `bedtools
  maskfasta` writes extra bases past the end (or aborts with exit 134 when the
  mask starts past it). With `--skip-invalid` gloci clips such masks to the
  sequence end, so the sequence keeps its length.
  [differences.md#maskfasta-past-end](differences.md#maskfasta-past-end)
- **`maskfasta --mask-char` refuses a space:** exit 2, where `bedtools
  maskfasta -mc ' '` writes spaces that make the FASTA unindexable.
  [differences.md#maskfasta-mask-char](differences.md#maskfasta-mask-char)
- **`cat` refuses BAM files with different `@SQ` lists:** exit 1, where
  `samtools cat` joins them and the records of the later files point at the
  wrong contigs. gloci's `cat` is BAM only (CRAM and SAM are refused with a hint).
  [differences.md#alignment-cat-sq](differences.md#alignment-cat-sq)
  It also refuses an input without the BGZF end-of-file marker (exit 1; samtools
  warns and exits 0), because the output gets a new marker and would hide the truncation:
  [differences.md#bgzf-eof-marker](differences.md#bgzf-eof-marker)
- **`sort -m` of alignments refuses different `@SQ` lists and conflicting `@RG` lines:**
  exit 1, where `samtools merge` appends the unknown `@SQ` lines (warning only when
  the order is lost) and `-c` keeps the first of two different `@RG` lines silently.
  [differences.md#alignment-merge-sq](differences.md#alignment-merge-sq),
  [differences.md#alignment-merge-rg](differences.md#alignment-merge-rg)
- **`vcf filter` keeps INFO by default:** All INFO fields are preserved unless
  pruned with `--keep-info` or `--drop-info` (`vcftools --recode` drops INFO
  unless `--recode-INFO-all` is given).
  [differences.md#info-kept-by-default](differences.md#info-kept-by-default)
- **REF and ALT as read:** `vcf filter` output, and the allele columns of
  `vcf stats` and `vcf extract`, keep the case of REF and ALT, where vcftools
  upper-cases them (a symbolic `<DEL:ME:Alu>` becomes `<DEL:ME:ALU>`).
  [differences.md#vcf-allele-case](differences.md#vcf-allele-case)
- **`vcf filter` keeps values as written:** where `bcftools view -i/-e/-f`
  prints every record again (`0.600000` becomes `0.6`), gloci writes the kept
  records unchanged.
  [differences.md#vcf-filter-values](differences.md#vcf-filter-values)
- **Inverted `vcf filter --max-missing`:** `--max-missing FRAC` is the maximum
  fraction of *missing* genotypes (`0.1` = at most 10% missing), the inverse of
  `vcftools --max-missing` (`0.9` = at least 90% called).
  [differences.md#max-missing-inverted](differences.md#max-missing-inverted)
- **`vcf filter` MAF and MAC on multi-allelic sites:** `--min-maf`/`--max-maf`
  and `--min-mac`/`--max-mac` test the second most common allele, not the rarest
  allele.
  [differences.md#maf-multiallelic](differences.md#maf-multiallelic)
- **`vcf filter --hwe` on multi-allelic sites:** `--hwe` tests diploid biallelic
  sites and lets multi-allelic sites pass (add `--biallelic` to drop them).
  [differences.md#hwe-multiallelic](differences.md#hwe-multiallelic)
- **`vcf filter --regions-file` on headerless BED:** Every non-header line of a
  BED file is used (`vcftools --bed` unconditionally skips line 1).
  [differences.md#bed-header-skipping](differences.md#bed-header-skipping)
- **`vcf extract` differs from `--get-INFO`, `--extract-FORMAT-info` and `--012`:**
  `--info` adds an ID column and writes `.` (not `?`) for missing values,
  `--format --layout wide` names columns `SAMPLE:KEY` and keeps every site, and
  `--genotype-matrix 012` counts a haploid or half-missing reference call as 0
  (vcftools counts 1).
  [differences.md#vcf-extract-info](differences.md#vcf-extract-info),
  [differences.md#vcf-extract-format-wide](differences.md#vcf-extract-format-wide),
  [differences.md#vcf-extract-012](differences.md#vcf-extract-012)
- **`vcf stats` replaces the vcftools statistic flags:** every flag is one
  `gloci vcf stats FILE --out-prefix P` run that writes `P.<ext>` in vcftools'
  formats, after the same selection flags as `vcf filter`, and publishes files
  only when the whole run succeeds. Differences: a singleton's carrier is
  searched among the unmasked genotypes, `--snp-density` rejects unsorted
  input, a `.` DP or GQ is never masked, and a polyploid genotype stops the run
  only in a sample that is selected and not masked.
  [differences.md#vcf-stats-singletons-masked](differences.md#vcf-stats-singletons-masked),
  [differences.md#vcf-stats-snp-density-unsorted](differences.md#vcf-stats-snp-density-unsorted),
  [differences.md#vcf-stats-mask-missing](differences.md#vcf-stats-mask-missing),
  [differences.md#polyploid-unused-samples](differences.md#polyploid-unused-samples).
  The layout of `--json` and the note on `vcf filter ... | vcf stats` (which
  keeps the original DP and GQ of masked genotypes) are in
  [differences.md](differences.md), under "`vcf stats --json`: the summary.json
  schema" and "`vcf filter ... | vcf stats` is `vcf stats` without masking".
- **`vcf popgen` replaces `--site-pi`, `--window-pi`, `--TajimaD`,
  `--weir-fst-pop`, the LD options, `--relatedness` and `--relatedness2`:** one
  `gloci vcf popgen FILE --out-prefix P` run writes any of `P.sites.pi`,
  `P.windowed.pi`, `P.Tajima.D`, `P.weir.fst`, `P.windowed.weir.fst`,
  `P.geno.ld`, `P.hap.ld`, `P.interchrom.{geno,hap}.ld`, `P.relatedness` and
  `P.relatedness2` in vcftools' formats, after the same selection
  flags as `vcf filter` (`--minDP 5` is `--min-dp 5`, `--indv s1` is `--samples
  s1`). `--weir-fst-pop` is `--fst-pop`, and `P.weir.fst.summary` holds the mean
  and weighted Fst vcftools prints in its log. Differences: the windowed
  statistics (`--geno-r2` and `--hap-r2` among them) reject unsorted input, a
  window step larger than the window or an LD window that admits no pair is an
  error rather than reset or empty, no count wraps past 32-bit limits, the
  population files select no samples (vcftools keeps only their samples, so its
  site filters count only them) and `--hap-r2` filters no sites from the other
  statistics, and one population, a population with no selected sample or a
  sample in two populations is an error. Fst prints the same digits on every
  platform, where vcftools on aarch64 rounds a few steps differently, and so do
  LD, where it prints rounding noise for an r² of 0, and relatedness. LD holds
  its sites in memory within `--max-mem`, where vcftools uses a temporary file,
  writes its headers for fewer than two sites, where vcftools stops, and does
  not read the phase of an excluded sample. Relatedness reads a masked or
  half-called genotype as missing, where vcftools counts a masked one in
  `--relatedness2` and reads `./1` as `0/0` in `--relatedness`. A polyploid
  genotype stops the run only in a sample that is used, where vcftools stops on
  any.
  [differences.md#popgen-sorted-input](differences.md#popgen-sorted-input),
  [differences.md#popgen-window-step](differences.md#popgen-window-step),
  [differences.md#popgen-integer-overflow](differences.md#popgen-integer-overflow),
  [differences.md#popgen-independent-options](differences.md#popgen-independent-options),
  [differences.md#fst-populations](differences.md#fst-populations),
  [differences.md#fst-rounding](differences.md#fst-rounding),
  [differences.md#ld-rules](differences.md#ld-rules),
  [differences.md#ld-rounding](differences.md#ld-rounding),
  [differences.md#relatedness-genotypes](differences.md#relatedness-genotypes),
  [differences.md#relatedness-rounding](differences.md#relatedness-rounding),
  [differences.md#polyploid-unused-samples](differences.md#polyploid-unused-samples)
- **`vcf filter --include/--exclude` is a subset of `bcftools view -i/-e`:**
  fields against numbers, strings or `.`, `N_PASS`, and the FORMAT fields GT,
  DP and GQ only (other FORMAT keys such as AD and FT are rejected, exit 2).
  Differences from bcftools 1.21: `!~` is a real negation, a bare INFO flag is
  its presence, `!(...)` and `.` literals parse, an index past the end is
  missing. The expression `MAF` is bcftools' first-ALT MAF, unlike `--min-maf`.
  [differences.md#expr-subset](differences.md#expr-subset)
  [differences.md#expr-format-fields](differences.md#expr-format-fields)
  [differences.md#expr-not-match](differences.md#expr-not-match)
  [differences.md#expr-bare-flag](differences.md#expr-bare-flag)
  [differences.md#expr-grammar-extensions](differences.md#expr-grammar-extensions)
  [differences.md#expr-index-past-end](differences.md#expr-index-past-end)
  [differences.md#expr-maf-first-alt](differences.md#expr-maf-first-alt)
- **`--samples` with `--include/--exclude`:** `bcftools view -s S -i E` filters on
  all samples; gloci filters on the selected samples. The bcftools equivalent
  of gloci's behaviour is the two-step pipeline
  `bcftools view -s S | bcftools view -i E`.
  [differences.md#expr-samples-selection](differences.md#expr-samples-selection)
- **Symbolic ALT without END or SVLEN:** a record such as `<DEL>` with neither
  INFO/END nor SVLEN is an error under the default `--vcf-span allele` (it has
  no allele span); add the field, pass `--vcf-span ref` (the REF span, as in
  bcftools and tabix) or `--skip-invalid`.
  [differences.md#vcf-symbolic-no-end](differences.md#vcf-symbolic-no-end)

- **`--chrom-map` with `--assume-sorted` means "sorted by the new names":**
  `view --assume-sorted` binary-searches a plain BED file and trusts the claim
  that it is sorted. With `--chrom-map` the order it trusts is that of the
  renamed chromosomes. Two cases:
  - A map that reverses the order (`a`, `b` renamed to `z`, `y`) is safe: the
    search does not find the chromosome, warns "could not locate ... scanning
    the whole file" and falls back to a scan, so no record is lost.
  - A map that renames onto a name already in the file is not. With the lines
    `a`, `a`, `b`, `c`, the map `a` to `c` and `-g c,b`, `view --assume-sorted
    c:1-1000` silently returns only the renamed `a` records, while the plain
    scan fails with "not sorted: line 4 (c ...) comes after ...". The claim
    that the file is sorted in new-name space is false there, and gloci does
    not check it.
  To query such a file, rename and re-sort first:
  `gloci view FILE --chrom-map MAP | gloci sort - -o renamed.bed.gz --write-index`,
  and query that. (A note, not a pinned difference: the old tools have no
  rename table to compare with. Both cases are tested in
  `crates/gloci/tests/chrom_map.rs`.)
- **`genomecov` counts past 32 bits:** a `genome` row counts bases in 64 bits.
  `bedtools genomecov` wraps that count at 32 bits when one depth covers more
  than about 4.29 Gb in total across the genome (for instance a 6 Gb genome
  with 5 Gb at depth 1), so its `genome` row is wrong there; it also needs
  about 24 GB of memory to get that far, so no differential case can pin it.
  gloci's side is pinned by `genome_row_counts_past_32_bits` in
  `crates/gloci/tests/genomecov.rs` (`genome 0 1000000000 6000000000 0.166667`
  and `genome 1 5000000000 6000000000 0.833333`). The fractions are single
  precision as in bedtools, including the `genome` rows
  (`genomecov-hist-genome-row-float-vs-bedtools`).
- **`vcf stats` validates FORMAT fields lazily:** only the fields a requested
  statistic or a mask rule reads are decoded, so a malformed DP, GQ or FT value
  that nothing reads is not an error (`--freq` ignores a bad DP; `--depth` or
  `--min-dp` reports it, exit 1, naming the record and sample).
- **SAM input is checked more strictly than `samtools view` checks it:** a
  RNAME or RNEXT that is not in the header is an error (samtools warns and makes
  the read unmapped), and a few values samtools wraps or truncates are refused.
  [unknown reference](differences.md#sam-unknown-rname),
  [parse strictness](differences.md#sam-parse-strictness)
- **CRAM needs `--reference FASTA` (bedtools `--cram-ref`):** there is no
  `REF_PATH`, cache or download; the FASTA needs its `.fai` (and `.gzi` when BGZF)
  and must have the lengths of the CRAM's `@SQ` lines. Every command that reads
  alignments takes it (`intersect`, `coverage`, `genomecov`, `convert`, `view`,
  `sort`), and warns when it is given with no CRAM input. A CRAM written without a
  reference (`no_ref`, `embed_ref`) needs none, and neither does `gloci index
  x.cram`. With its `.crai` a CRAM is read by region, and per chromosome at `-t`
  above 1, like an indexed BAM (the reference is read once per chromosome, and
  at most one chromosome per worker is held); an error in a record names "record N
  of CHROM" without the file-wide number a BAM adds. A region that starts past a chromosome's end finds the reads that run past it in a BAM but not in a CRAM, as in samtools ([differences.md#region-past-reference-end](differences.md#region-past-reference-end)). At `-t` above 1 the container of reads with no place is never decoded, so damage there is seen only when streaming ([differences.md#cram-unplaced-tail](differences.md#cram-unplaced-tail)). MD and NM come out as stored, where samtools recomputes them (compare
  with `samtools view --input-fmt-option decode_md=0`):
  [differences.md#cram-md-nm](differences.md#cram-md-nm). A CRAM without its
  end-of-file container is refused (exit 1; samtools warns and exits 0):
  [differences.md#cram-eof-container](differences.md#cram-eof-container). A read
  more than 1 Mb past its chromosome's end, and a few shapes only a corrupt or crafted
  file has, are refused where samtools decodes them:
  [differences.md#cram-limits](differences.md#cram-limits)

### Randomness

`--seed N` makes `random`, `sample` and `shuffle` reproduce gloci's own output
at any `-t`, on any platform, and never bedtools': the same seed gives other
intervals and records than `bedtools random|sample|shuffle -seed N`, which use
the C++ standard library's generator. Without `--seed` a fresh seed is chosen
and printed as a warning, so a run can be repeated from its log. Samples are
as uniform as bedtools' (checked by chi-squared tests at fixed seeds), not the
same samples. [differences.md#rng-stream](differences.md#rng-stream)

### Memory of the statistics commands

- `shuffle` holds its allowed intervals: about 52 bytes per interval of
  `--include` (or of the genome minus `--exclude`), plus the records it moves
  as it streams them.
- `random` holds one batch of 262,144 rows at a time, whatever `-n`.
- `pair-intersect` and `pair-pair` hold B in memory, at about 5 to 6 times
  its text size; A streams.
- `stats jaccard`, `fisher` and `reldist` stream sorted inputs in flat memory
  (`reldist --unsorted` holds every B midpoint).
- `vcf compare` reads A and B on threads of their own; it uses at most
  max(`-t`, 3) threads.


### Threads and errors

`-t N` never changes a successful result: output is byte-identical at every
thread count. After an error the exit code is the same too, and so is the
message for a single malformed input. When more than one input is malformed,
the error reported may name a different input at different thread counts; the
exit code is the same. How much output was written to standard output before an
error also varies with `-t` (`-o FILE` is written atomically).
