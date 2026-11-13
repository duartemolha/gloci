# gloci

gloci is a single fast binary for genomic interval, VCF, alignment and
sequence work. It replaces the everyday parts of bedtools, bedops, vcftools,
tabix and bgzip, and common samtools commands, and it is built to handle
multi-GB VCF, BED and BAM files.

It reads BAM, CRAM and SAM in the interval commands (`intersect`, `coverage`,
`genomecov`, `convert`), and covers samtools `view`, `sort`, `merge`
(`gloci sort -m`), `cat` and `index` (`gloci index ref.fa` for FASTA). FASTA
tools are `getfasta`, `maskfasta`, `nuc`, `faidx` and `dict`. The samtools
essentials are there too: `flags`, `idxstats`, `flagstat`, `head`, `samples`,
`depth`, `bedcov`, `reheader`, `view --count`, and `check` (`check --quick` is
`samtools quickcheck`).

## Install

From a GitHub Release: download the archive for your platform from the
[Releases page](https://github.com/duartemolha/gloci/releases), unpack it and
put `gloci` on your `PATH`.

From source, with a Rust toolchain (1.91 or newer):

```sh
cargo install --git https://github.com/duartemolha/gloci gloci
```

Not yet published: gloci is not on crates.io, PyPI or bioconda yet, so
`cargo install gloci`, `pip install gloci` and `conda install gloci` do not
work today.

## Example

```sh
gloci genomes show hg38 > hg38.genome
gloci sort peaks.bed -g hg38 > peaks.sorted.bed
gloci sort reads.bam -o reads.sorted.bam --write-index
gloci genomecov reads.sorted.bam --bedgraph -o coverage.bedgraph
```

`gloci --help` lists every subcommand, and `gloci <subcommand> --help` gives
its options and examples.

## Performance

Seconds, median of 5 runs (3 for the VCF row), on real data: GENCODE v46
exons (1.67 million intervals), ENCODE cCREs (2.35 million), an ENCODE CTCF
ChIP-seq BAM and CRAM (7.6 million reads) and the 1000 Genomes chr20 panel
(BCF). One thread unless noted. The other tools are bedtools 2.31.1, bedops
2.4.41 (sort-bed in the sort row), samtools and bcftools 1.21 and vcftools
0.1.16.

The commands of a row give the same result, which the benchmark runner now
checks before it times them. Where the tools answer different questions there
is a row for each form, with gloci doing that same work.

| Task | gloci | bedtools | bedops | other |
|---|---:|---:|---:|---|
| sort exons | 0.36 | 1.64 | 0.92 | GNU sort 1.69 |
| merge exons | 0.20 | 0.39 | 0.36 | |
| intersect cCREs with exons | 0.40 | 1.15 | 0.92 | |
| the same, BGZF inputs, 8 threads | 0.27 | 1.42 | 1.16 | |
| subtract cCREs from exons, each exon cut | 0.63 | 2.21 | | |
| the same on merged exons (`bedops -d`) | 0.48 | | 0.87 | |
| closest exons to each CTCF peak, ties kept | 0.23 | 0.58 | | |
| one nearest exon on each side (`closest-features`) | 0.21 | | 0.41 | |
| coverage of cCREs by exons | 1.10 | 2.61 | | |
| genomecov of exons | 0.29 | 14.40 | | |
| partition of three files at every boundary | 0.62 | | 1.38 | |
| the same with equal neighbours joined (`multiinter`) | 0.95 | 4.55 | | |
| map CTCF signal onto cCREs | 9.77 | 19.76 | 39.87 | |
| flagstat, BAM | 1.48 | | | samtools 1.53 |
| flagstat, CRAM | 0.24 | | | samtools 0.46 |
| bedcov of CTCF peaks, BAM | 4.08 | | | samtools 4.46 |
| depth, BAM | 20.6 | | | samtools 24.3 |
| depth, CRAM | 27.8 | | | samtools 23.2 |
| vcf filter, BCF to compressed VCF, 8 threads | 209.6 | | | bcftools 63.7, vcftools 504 |

Where gloci is slower. Filtering a BCF on genotypes and writing VCF, bcftools
was more than three times as fast: gloci counted genotypes and wrote VCF text
several times slower than bcftools does, and with more threads than cores
(this machine has four) it lost more: with four threads it took 123 to 141
seconds. That row was measured before gloci's BCF-to-VCF writer was replaced
by one that writes the text as htslib does; the change took the command's CPU
time at one thread from 215 to 110 seconds, and the row has not been timed
again, so read its gloci figure as an upper bound. Counting genotypes is
still several times slower than in bcftools. `depth` on a CRAM costs gloci
16.9 seconds of CPU against 9.5 for samtools, which decodes only the fields
`depth` reads; both `depth` rows are mostly the time to write 3.3 GB of
output.

In the BGZF row bedtools has one thread, and bedops, which cannot read BGZF, is
timed with the decompression of both inputs to files.

Memory: gloci holds 18 to 67 MiB in the single-thread rows, 260 MiB for
`sort` (which keeps the file), 288 MiB for the 8-thread intersect and 81 MiB
for the VCF filter. bedops and samtools stream with 3 to 14 MiB;
`bedtools genomecov` holds 1.9 GiB.

These numbers are indicative. They were measured on 13 November 2026 on a
4-core aarch64 machine (WSL2, 11 GiB), where long runs vary: the `depth` rows
by 2 to 4 seconds between runs, and vcftools between 262 and 504 seconds. The
commands were checked on the same day by comparing what each tool writes; the
first version of this table had errors in both directions (see
[paper/README.md](paper/README.md)). To measure on your own machine:

```sh
docker build -t gloci-bench:p1 -f xtask/docker/Dockerfile.bench xtask/docker
cargo xtask paper fetch              # about 6.7 GB of public data
cargo xtask paper run --bench-only
```

## Stability

gloci is young. `genomes`, `index`, `sort`, `view`, `intersect`, `merge`,
`subtract`, `closest` and `coverage` are the stable core: their options and
output change only with a version bump and a changelog entry. Every other
command is experimental and says so at the top of its `--help`: it is checked
against the tool it replaces, but its options and output may still change.

Before moving a pipeline step to gloci, run it beside the tool it replaces on
your own data. [contrib/gloci-shadow.sh](contrib/gloci-shadow.sh) does that
for one step and logs whether the outputs were identical; see
[Stability](docs/book/src/stability.md).

## Coming from bedtools, bedops, vcftools or samtools

See the [migration guide](docs/migration.md) for the command mapping, and
[docs/differences.md](docs/differences.md) for the places where gloci
intentionally differs.

## Licence

MIT. See [LICENSE](LICENSE).
