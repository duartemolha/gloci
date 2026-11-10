# Migrating from bedtools, bedops and vcftools

gloci (pronounced GLOH-sigh, like "loci") is one tool for genomic interval and
variant work.

This guide is for people moving existing pipelines to gloci. It starts with the
behaviour changes that break or alter a pipeline even after the flags are
rewritten.

## Chromosome order (behaviour change)

Without `-g` and without a VCF/BCF header that lists contigs, gloci uses
**karyotype order**: `chr1, chr2 ... chr22, chrX, chrY, chrM` (see
[differences.md](differences.md#karyotype-chromosome-order)). `sort -k1,1`,
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

1. Re-sort once, then keep the files: `gloci sort x.bed > x.sorted.bed` (or add
   `--sort` to a streaming relation, `merge` or `complement` command to sort on
   the fly).
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

After, with karyotype order:

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

Every row below was checked against `gloci <cmd> --help`.

Rows marked "not yet available" name the command that will replace the tool
and the phase that builds it. After the P2 command set (this page), the build
order is: the packaging prerequisites of P2b (licence, genome tables inside
the crate), then P3a (BAM/CRAM, FASTA, `getfasta`, `maskfasta`, `nuc`, and BAM
input to `coverage`, `genomecov`, `intersect` and `convert`), then P3b
(`stats`, `sample`, `random`, `shuffle`, BEDPE, `vcf popgen`, `vcf compare`),
then the rest of P2b (Python library, packaging, mdBook, paper evidence suite).

### Sorting, indexing and viewing (htslib, bcftools, bedops, bedtools)

| Old command | gloci |
|---|---|
| `sort-bed FILE`, `bedtools sort -i FILE`, `LC_ALL=C sort -k1,1 -k2,2n FILE` | `gloci sort FILE --chrom-order lex` (drop `--chrom-order lex` to get karyotype order) |
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
| `bcftools norm -m- --no-version FILE`, `bcftools norm -m-snps\|-indels\|-both\|-any FILE` (no `-f`) | `gloci vcf split FILE [--types snps\|indels\|both\|any]`: one record per ALT; GT keeps the ALT as 1 and recodes the others to 0; Number=A/R/G INFO and FORMAT fields are subset (haploid G too). Alleles are not trimmed. Values keep their text, and a record that is not split is written as it is: [differences.md#vcf-split-values](differences.md#vcf-split-values). A wrong number of values for a field is an error (bcftools: `--force`). `norm -m+` (merging), `-f`, `-d` and `-c` are not available |
| `bcftools index FILE.bcf`, `bcftools index -c\|-t FILE.vcf.gz` | `gloci index FILE.bcf` (a BCF always gets a CSI; `.vcf.gz` gets a TBI, or a CSI when needed) |
| `bgzip FILE && tabix -p vcf FILE.gz` | `gloci sort FILE -o FILE.gz --write-index`, or `gloci view FILE -o FILE.gz --write-index` for a file that is already sorted |
| `tabix FILE.gz REGION...`, `bcftools view -r REGION FILE`, `bcftools view -H -r ...` | `gloci view FILE.gz REGION...` (or `--region REGION`; add `--no-header` for `-H`) |
| `tabix -h FILE.gz REGION`, `tabix -H FILE.gz` | `gloci view FILE.gz REGION` keeps the header; `gloci view --header-only FILE.gz` prints only it |
| `tabix -R TARGETS FILE.gz`, `bcftools view -R TARGETS FILE` | `gloci view FILE.gz --regions-file TARGETS.bed` (0-based BED, like `bedextract`) |
| `bedextract CHR FILE`, `bedextract FILE REGIONS.bed` | `gloci view FILE CHR`, or `gloci view FILE --regions-file REGIONS.bed`; add `--assume-sorted` for the binary search on an uncompressed BED file |
| `samtools index FILE.bam`, `samtools index -c FILE.bam` | `gloci index FILE.bam` (a `.bai`, or a `.csi` when a reference is longer than 2^29); `gloci index FILE.bam --csi` always writes a `.csi`. The BAM must be coordinate-sorted. SAM text has no index: `gloci view FILE.sam -o FILE.bam --write-index`. samtools and gloci read each other's indexes. [differences.md#bam-index-long-reference](differences.md#bam-index-long-reference) |
| `samtools view FILE.bam REGION...` | `gloci view -O sam FILE.bam REGION...` (through `FILE.bam.bai`, `FILE.bai` or `FILE.bam.csi`; without an index gloci scans the file and warns, where samtools exits 1; [differences.md#bam-region-without-index](differences.md#bam-region-without-index)). `gloci view -O sam` writes the header, as `samtools view -h` does; add `--no-header` for `samtools view` without `-h`. Regions are 1-based and closed, a read is selected by where its CIGAR puts it, and a placed unmapped read counts as one base, as samtools finds them |
| `samtools view -b -o X.bam --write-index FILE`, `samtools cat -o X.bam A.bam B.bam`, then `samtools index X.bam` | `gloci view FILE -o X.bam --write-index`, `gloci cat A.bam B.bam -o X.bam --write-index` (writes `X.bam.bai`; the data must be coordinate-sorted). An old `X.bam.bai`, `X.bai` or `X.bam.csi` beside an output written without `--write-index` is removed |
| `bcftools view -Ob FILE.vcf`, `bcftools view -Ov FILE.bcf` | `gloci view FILE.vcf -o FILE.bcf`, `gloci view FILE.bcf -O vcf` |

### bedtools

| bedtools | gloci |
|---|---|
| `intersect -u` | `gloci intersect -a A -b B` (default `--mode any`; note `--mode` has no `-m`) |
| `intersect -abam A.bam -b B.bed -u` (`-v`, `-split`, `-s`, `-S`, `-f`, `-c`) | `gloci intersect -a A.bam -b B.bed [--mode none] [--split] [-s] [-S] [-f FRAC]` (BAM or SAM as A or B, which writes BAM, or SAM for a SAM A, and `-o x.bam` or `-O bam\|sam` set it; output is the same reads as bedtools 2.31.1; `--split` cuts at CIGAR `N`, and `--split-deletions` also at `D`; unmapped reads are never written, even by `--mode none`, see [differences.md#unmapped-reads](differences.md#unmapped-reads); `-abam -bed -wa -wb` and `--mode pairs` with an alignment A print BED6, not BED12, see [differences.md#bam-as-bed-columns](differences.md#bam-as-bed-columns); `--mode clip` of reads is an error; for BED write BAM and pipe to `gloci convert -`) |
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
| `complement -g` | `gloci complement FILE... [-g GENOME]` (`-g` takes a built-in such as `hg38` or a `.genome`/`.fai`/`.dict` file; `--sort`) |
| `closest -d -D -k -t -io -iu -id -mdb` | `gloci closest -a A -b B --distance --signed ref\|a\|b -k N --ties all\|first\|last --ignore-overlaps --ignore-upstream --ignore-downstream --per-file` (plus `--max-distance N`, `-s`/`--same-strand`, `-S`/`--opposite-strand`) |
| `map -c -o` | `gloci map -a A -b B -c/--column COL:OP` (`--delim`, `--empty`, `--skip-unmatched`; ops: `count`, `count_distinct`, `sum`, `mean`, `median`, `min`, `max`, `stdev`, `var`, `stdev_pop`, `var_pop`, `collapse`, `distinct`, `distinct_sorted`, `first`, `last`, `bases`, `indicator`, `cv`, `mad`, `quantile:Q`, `tmean:LO-HI`, `wmean`, `mode`, `antimode`, `bases_uniq`, `frac_a`, `frac_b`, `min_record`, `max_record`, `records`) |
| `coverage [-hist -d -mean]`, `annotate` | `gloci coverage -a A -b B [--hist \| --per-base \| --mean]`, `gloci coverage -a A -b B1 B2... [--names ...]` |
| `coverage -a A.bed -b reads.bam [-split -hist -d -mean]` | `gloci coverage -a A.bed -b reads.bam [--split] [--hist \| --per-base \| --mean]` (BAM or SAM as B; `--split-deletions` also cuts reads at `D`; `--min-mapq`, `--require-flags` and `--exclude-flags` filter reads; an alignment A prints as BED6, see [differences.md#bam-as-bed-columns](differences.md#bam-as-bed-columns)) |
| `slop -b -l -r -s -pct -g` | `gloci slop FILE -b/--both N` or `--left L --right R` with `-s`, `--pct`, `-g GENOME` |
| `flank -b -l -r -s -pct -g` | `gloci flank FILE -b/--both N` or `--left L --right R` with `-s`, `--pct`, `-g GENOME` |
| `shift -s N`, `shift -p X -m Y`, `-pct`, `-g` | `gloci shift FILE --by N [-g GENOME]` (signed; `--pct` reads N as an exact decimal fraction of the feature length); `-p X -m -X` is `--by X -s`. bedtools' `-p`/`-m` with unrelated magnitudes has no gloci equivalent: run `shift` per strand after filtering |
| `makewindows -g/-b -w -s -n -i` | `gloci windows -g GENOME \| -i/--intervals FILE` with `-w/--size SIZE [--step S]` or `-n/--count COUNT`, plus `--id-column`, `--reverse` |
| `sort -i FILE [-g GENOME]` | `gloci sort FILE [-g GENOME]` (`--chrom-order lex` reproduces byte order) |
| `bed12tobed6 -i FILE` | `gloci convert FILE --blocks` |
| `genomecov -i x.bed -g G` (`-bg`, `-bga`, `-d`, `-scale`, `-strand`, `-split`, `-5`, `-3`, `-max`) | `gloci genomecov x.bed -g G [--bedgraph \| --bga \| --per-base] [--scale X] [--strand +\|-] [--split] [--5p \| --3p] [--max-depth N]` (`-g` is required, as in bedtools; the input must be sorted, `--sort` sorts it; `-dz` and `-trackline` have no counterpart; the `--hist` fractions divide in single precision and print six significant digits, as bedtools does; a `genome` row that counts more than 2^32 bases at one depth is correct in gloci, see the `genomecov` note under Other behaviour changes; BAM input (`-ibam`) is not available until P3a) |
| `multiinter -i A B C [-names]` | `gloci partition A B C --merge-adjacent [--names X,Y,Z]` (the count and label columns; the per-file 0/1 columns are not written; without `-names` bedtools labels the files `1,2,3` and gloci uses the file names (pass `--names 1,2,3` to match); `-cluster` and `-empty -g` have no counterpart (for the latter use `gloci complement -g`)) |
| `unionbedg -i A B [-filler -empty -g]` | `gloci partition A B --bedgraph [--empty VALUE [-g GENOME]]` (`-filler X` is `--empty X`; `--empty` also writes the uncovered segments, as it implies `--min-files 0`; an explicit `--min-files` overrides that) |
| `groupby -g -c -o [-header]` | `gloci groupby FILE -g COLS -c/--column COL:OP` (bedtools `-o OP` becomes the `:OP` of `-c COL:OP`, and gloci's `-o` is the output file; `--header`, `--delim`, `--empty`; same ops as `map` except those that need a second input; `distinct` keeps first-seen order where bedtools sorts it, use `distinct_sorted`, see [differences.md#map-distinct-order](differences.md#map-distinct-order); VCF input is refused, use `gloci vcf extract`) |
| `cluster -d -s` | `gloci cluster FILE -d/--max-gap N -s` (`--sort`; the same gap and strand rules as `gloci merge`; ids are global, and `-s` keeps input order, see [differences.md#cluster-strand](differences.md#cluster-strand); a VCF input is clustered by its footprints and the output is TSV with a `#CHROM ... cluster` line, see [differences.md#append-vcf-header](differences.md#append-vcf-header); BCF is refused, use `gloci view -O vcf`; bedtools accepts chromosomes in any grouped order, gloci needs a sorted chromosome order, see [Chromosome order](#chromosome-order-behaviour-change)) |
| `spacing` | `gloci spacing FILE` (`.` for the first record of a chromosome, `0` bookended, `-1` overlapping, only the previous record counts, as in bedtools; a zero-length record is a point, see [differences.md#spacing-zero-length](differences.md#spacing-zero-length); VCF input and chromosome order as for `cluster`) |
| `bamtobed -i FILE` (`-bed12`, `-split`, `-splitD`, `-ed`, `-tag TAG`, `-cigar`, `-color R,G,B`) | `gloci convert FILE [--bed12] [--blocks] [--blocks --split-deletions] [--score-tag NM\|TAG] [--cigar] [--color R,G,B]` (BAM, CRAM or SAM input, which bedtools reads as BAM only; output matches bedtools 2.31.1 in input order except for the pinned differences below, with `/1` and `/2` on mates, MAPQ as the score and unmapped reads skipped; `-ed` is `--score-tag NM`; `--min-mapq`, `--require-flags` and `--exclude-flags` filter; a read with no reference bases is one base, not zero-length, see [differences.md#bamtobed-no-reference-ops](differences.md#bamtobed-no-reference-ops); reads are cut at the end of their reference, see [differences.md#bamtobed-past-reference-end](differences.md#bamtobed-past-reference-end); a negative tag score is printed signed, where bedtools wraps it to unsigned 32-bit, see [differences.md#bamtobed-negative-score](differences.md#bamtobed-negative-score); `--blocks --score-tag` has no extra column and `--cigar` with `--blocks` is an error, see [differences.md#bamtobed-split-tag](differences.md#bamtobed-split-tag)) |
| `getfasta -fi REF -bed FILE` (`-name`, `-nameOnly`, `-s`, `-split`, `-tab`, `-bedOut`, `-rna`) | `gloci getfasta FILE --fasta REF [--name \| --name-only] [-s] [--split] [--tab \| --bed-out] [--rna]` (output is byte-identical to bedtools 2.31.1 for BED, GFF and VCF input, including the IUPAC-aware reverse complement, except that `--split` refuses BED12 records whose blocks are out of order, overlap, or do not span the record where bedtools writes them unchecked, see [differences.md#getfasta-split-blocks](differences.md#getfasta-split-blocks); the FASTA needs a `.fai` (and `.gzi` for BGZF) made with `gloci index`, which gloci never builds itself; VCF is bedtools' REF span with `--vcf-span ref`; zero-length intervals are skipped with a counted warning, as bedtools does; an interval on a chromosome the FASTA lacks, or past its end, is an error where bedtools skips it, `--skip-invalid` gives bedtools' behaviour, see [differences.md#getfasta-outside-fasta](differences.md#getfasta-outside-fasta); there is no `-fullHeader`) |
| `maskfasta -fi F -bed B [-soft -mc C -fullHeader]` | `gloci maskfasta --fasta F -b B [--soft \| --mask-char C] [-o OUT.fa[.gz] --write-index]` (streamed, any compression, no index needed; headers are kept whole and line layout is kept, so `-fullHeader` is the default, see [differences.md#maskfasta-header](differences.md#maskfasta-header) and [differences.md#maskfasta-layout](differences.md#maskfasta-layout); a mask past its sequence's end is an error, see [differences.md#maskfasta-past-end](differences.md#maskfasta-past-end)) |
| `nuc -fi REF -bed FILE` (`-s`, `-seq`, `-pattern P`, `-C`) | `gloci nuc FILE --fasta REF [-s] [--seq] [--pattern P [-C]]` (output is byte-identical to bedtools 2.31.1, header line and f32 percentages included, for BED, GFF and VCF input; the counts and the pattern describe the reverse complement under `-s`; the FASTA needs a `.fai` made with `gloci index`; VCF is bedtools' REF span with `--vcf-span ref`; zero-length intervals are skipped with a counted warning, as bedtools does; an interval outside the FASTA is an error where bedtools skips it, `--skip-invalid` gives bedtools' behaviour, see [differences.md#getfasta-outside-fasta](differences.md#getfasta-outside-fasta)) |
| `jaccard`, `fisher`, `reldist` | not yet available (`stats jaccard / fisher / reldist`, P3b) |
| `shuffle`, `random`, `sample` | not yet available (`shuffle`, `random`, `sample`, P3b) |
| `pairtobed`, `pairtopair` | not yet available (`pair-intersect`, `pair-pair`, P3b) |
| `multicov`, `tag`, `bamtofastq`, `bedtobam`, `bedpetobam`, `igv`, `links`, `expand`, `split`, `summary`, `overlap` | not planned (`overlap` → `gloci intersect --mode pairs --overlap-bases`) |

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
| `sort-bed [--max-mem --tmpdir --check-sort --unique --duplicates]` | `gloci sort [--max-mem --tmp-dir --check --unique --duplicates]`; add `--chrom-order lex` to reproduce `sort-bed`'s chromosome order (`--full-key` reproduces whole-line tie breaking) |
| `bedextract CHR FILE`, `bedextract FILE REGIONS.bed` | `gloci view FILE CHR`, `gloci view FILE --regions-file REGIONS.bed` (add `--assume-sorted` for binary search on an uncompressed BED file) |
| `convert2bed`, `vcf2bed [--snvs/--insertions/--deletions]`, `gff2bed`, `gtf2bed` | `gloci convert FILE [--types snv/ins/del] [--split-alleles] [--feature TYPE] [--name-attr KEY]` |
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
| `--site-pi`, `--window-pi`, `--TajimaD`, `--weir-fst-pop`, `--*-r2`, `--relatedness*` | not yet available (`vcf popgen`, P3b) |
| `--diff*`, `vcf-compare` | not yet available (`vcf compare`, P3b) |
| `vcf-merge`, `--plink*`, `--IMPUTE`, `--BEAGLE-*`, `--ldhat*` | not planned |

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
- **`genomecov` ignores points, treats an unknown strand as plus and refuses
  features past the chromosome end:**
  [points](differences.md#genomecov-points),
  [strand](differences.md#genomecov-unknown-strand),
  [past the end](differences.md#genomecov-past-end)
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
- **`vcf filter` keeps INFO by default:** All INFO fields are preserved unless
  pruned with `--keep-info` or `--drop-info` (`vcftools --recode` drops INFO
  unless `--recode-INFO-all` is given).
  [differences.md#info-kept-by-default](differences.md#info-kept-by-default)
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
  input, and a `.` DP or GQ is never masked.
  [differences.md#vcf-stats-singletons-masked](differences.md#vcf-stats-singletons-masked),
  [differences.md#vcf-stats-snp-density-unsorted](differences.md#vcf-stats-snp-density-unsorted),
  [differences.md#vcf-stats-mask-missing](differences.md#vcf-stats-mask-missing).
  The layout of `--json` and the note on `vcf filter ... | vcf stats` (which
  keeps the original DP and GQ of masked genotypes) are in
  [differences.md](differences.md), under "`vcf stats --json`: the summary.json
  schema" and "`vcf filter ... | vcf stats` is `vcf stats` without masking".
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

### Threads and errors

`-t N` never changes a successful result: output is byte-identical at every
thread count. After an error the exit code is the same too, and so is the
message for a single malformed input. When more than one input is malformed,
the error reported may name a different input at different thread counts; the
exit code is the same. How much output was written to standard output before an
error also varies with `-t` (`-o FILE` is written atomically).
