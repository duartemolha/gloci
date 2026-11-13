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

Seconds, median of 5 runs, on real data: GENCODE v46 exons (1.67 million
intervals), ENCODE cCREs (2.35 million), an ENCODE CTCF ChIP-seq BAM and CRAM
(7.6 million reads) and the 1000 Genomes chr20 panel (BCF). One thread unless
noted. The other tools are bedtools 2.31.1, bedops 2.4.41, samtools and
bcftools 1.21, vcftools 0.1.16 and mosdepth 0.3.10.

| Task | gloci | bedtools | bedops | other |
|---|---:|---:|---:|---|
| sort exons | 0.57 | 1.74 | 1.03 | GNU sort 0.68 |
| merge exons | 0.33 | 0.49 | 0.45 | |
| intersect cCREs with exons | 0.53 | 1.15 | 0.91 | |
| the same, BGZF inputs, `-t 8` | 0.35 | 1.43 | 1.12 | |
| subtract cCREs from exons | 1.08 | 2.26 | 0.86 | |
| closest exon to each CTCF peak | 0.33 | 0.62 | 0.44 | |
| coverage of cCREs by exons | 1.24 | 3.28 | | |
| genomecov of exons | 0.29 | 15.51 | | |
| partition of three files | 0.59 | 4.68 | 1.46 | |
| map CTCF signal onto cCREs | 10.28 | 20.42 | 1.36 | |
| flagstat, BAM | 1.47 | | | samtools 1.51 |
| flagstat, CRAM | 0.25 | | | samtools 0.48 |
| bedcov of CTCF peaks, BAM | 4.17 | | | samtools 4.53 |
| depth, BAM | 23.5 | | | samtools 30.3, mosdepth 31.5 |
| depth, CRAM | 33.9 | | | samtools 27.6 |
| vcf filter, BCF | 123.2 | | | bcftools 107.3, vcftools 262.4 |

gloci is not the fastest everywhere: bedops is quicker at `map` (by a wide
margin) and `subtract`, samtools at `depth` on CRAM, and bcftools at filtering
a BCF. gloci also holds more memory than bedops and samtools, which stream
with a few MiB: 14 to 67 MiB in the single-thread runs, 126 MiB for the
8-thread intersect and 257 MiB for `sort`, which keeps the file. It holds far
less than `bedtools genomecov` (1.9 GiB) and mosdepth (4.7 GiB).

These numbers are indicative: they were measured on 13 November 2026 on a
4-core aarch64 machine (WSL2, 11 GiB), where the two `depth` rows varied by
several seconds between runs. To measure on your own machine:

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
