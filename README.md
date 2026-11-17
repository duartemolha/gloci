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
gloci genomes list
gloci genomes show hg38 > hg38.genome
gloci sort peaks.bed -g hg38 > peaks.sorted.bed
gloci sort reads.bam -o reads.sorted.bam --write-index
gloci genomecov reads.sorted.bam --bedgraph -o coverage.bedgraph
```

`gloci --help` lists every subcommand, and `gloci <subcommand> --help` gives
its options and examples.

## Performance

Each cell is the wall time in seconds, then the peak memory in MiB: `0.32 / 118`
is 0.32 seconds and 118 MiB. The time is the fastest of 3 runs: on the machine
below, writing a large output to disk sometimes stalls for a second or more, and
the fastest run leaves those stalls out. The memory is the peak resident set
size of one more run, from GNU `time`; for a pipeline it is the largest process,
not the sum. The data is real: GENCODE v46 exons (1.67 million intervals),
ENCODE cCREs (2.35 million), an ENCODE CTCF ChIP-seq BAM and CRAM (7.6 million
reads) and the 1000 Genomes chr20 panel (BCF). One thread unless noted. The
other tools are bedtools 2.31.1, bedops 2.4.41 (sort-bed in the sort row),
samtools and bcftools 1.21 and vcftools 0.1.16.

The commands of a row give the same result, which the benchmark runner checks
before it times them. Where the tools answer different questions there is a row
for each form, with gloci doing that same work.

| Task | gloci | bedtools | bedops | other |
|---|---:|---:|---:|---|
| sort exons | 0.32 / 118 | 1.43 / 669 | 0.85 / 92 | GNU sort 1.09 / 133 |
| merge exons | 0.24 / 11 | 0.42 / 4.9 | 0.36 / 5.5 | |
| intersect cCREs with exons | 0.48 / 16 | 1.18 / 37 | 0.90 / 5.5 | |
| the same, BGZF inputs, 8 threads | 0.41 / 81 | 1.53 / 39 | 2.56 / 5.5 | |
| subtract cCREs from exons, each exon cut | 1.12 / 17 | 2.32 / 21 | | |
| the same on merged exons (`bedops -d`) | 0.51 / 13 | | 0.85 / 5.5 | |
| closest exons to each CTCF peak, ties kept | 0.25 / 14 | 0.53 / 6.4 | | |
| one nearest exon on each side (`closest-features`) | 0.24 / 14 | | 0.36 / 5.3 | |
| coverage of cCREs by exons | 1.93 / 16 | 3.41 / 101 | | |
| genomecov of exons | 0.37 / 11 | 15.34 / 1,903 | | |
| partition of three files at every boundary | 0.93 / 12 | | 1.48 / 5.5 | |
| the same with equal neighbours joined (`multiinter`) | 0.94 / 12 | 4.85 / 4.1 | | |
| map CTCF signal onto cCREs | 11.07 / 14 | 20.72 / 37 | 40.58 / 11 | |
| flagstat, BAM | 1.44 / 8.8 | | | samtools 1.50 / 3.0 |
| flagstat, CRAM | 0.27 / 16 | | | samtools 0.45 / 7.8 |
| bedcov of CTCF peaks, BAM | 4.14 / 28 | | | samtools 4.42 / 14 |
| depth, BAM | 31.6 / 10 | | | samtools 34.6 / 2.8 |
| depth, CRAM | 37.9 / 44 | | | samtools 33.4 / 7.8 |
| vcf filter, BCF to compressed VCF, 8 threads | 22.9 / 30 | | | bcftools 34.9 / 15, vcftools 257 / 5.6 |

Where gloci is slower. `depth` on a CRAM costs gloci 19.9 seconds of CPU against
9.7 for samtools, which decodes only the fields `depth` reads; both `depth` rows
are mostly the time to write 3.3 GB of output. In the VCF row gloci used 67
seconds of CPU across its 8 threads and bcftools 118.

In the BGZF row bedtools has one thread, and bedops, which cannot read BGZF, is
timed with the decompression of both inputs to files.

Memory: bedops and samtools stream, holding 2.8 to 14 MiB (sort-bed, which sorts
in memory, holds 92). gloci holds more: 8.8 to 44 MiB in the single-thread rows
other than `sort`, 118 MiB for `sort` (the file plus 32 bytes a record), 81 MiB
for the 8-thread intersect and 30 MiB for the 8-thread VCF filter. bedtools
holds the most where it loads a whole file: 669 MiB for `sort`, 101 MiB for
`coverage` and 1.9 GiB for `genomecov`.

These numbers are indicative. They were measured on 17 November 2026 on a 4-core
aarch64 laptop (WSL2, 11 GiB) on its fastest power setting, each case started
only once a short gloci run showed the machine at full speed. Runs that write
large outputs vary most there: `coverage` took gloci 1.93 to 4.06 seconds and
bedtools 3.41 to 4.87, and the `depth` rows varied by up to 13 seconds. The
commands were audited on 13 November by comparing what each tool writes; the
first version of this table had errors in both directions (see
[paper/README.md](paper/README.md)). To measure on your own machine:

```sh
docker build -t gloci-bench:p1 -f xtask/docker/Dockerfile.bench xtask/docker
cargo xtask paper fetch              # about 6.7 GB of public data
cargo xtask paper run --bench-only --runs 3
```

## Stability

gloci is young. `genomes`, `index`, `sort`, `view`, `intersect`, `merge`,
`subtract`, `closest` and `coverage` are the stable core: their options and
output change only with a version bump and a changelog entry. Every other
command is experimental and says so at the top of its `--help`: it is checked
against the tool it replaces, but its options and output may still change.
`gloci --help` lists the commands in groups and marks every group but the
stable core as experimental.

Before moving a pipeline step to gloci, run it beside the tool it replaces on
your own data. [contrib/gloci-shadow.sh](contrib/gloci-shadow.sh) does that
for one step and logs whether the outputs were identical; see
[Stability](docs/book/src/stability.md).

## Coming from bedtools, bedops, vcftools or samtools

See the [migration guide](docs/migration.md) for the command mapping, and
[docs/differences.md](docs/differences.md) for the places where gloci
intentionally differs.

## Design

The design specification gloci was built from is [`docs/design.md`](docs/design.md), and the numbered design
decisions the source cites (`D37`) are in [`docs/decisions.md`](docs/decisions.md).

## Licence

MIT. See [LICENSE](LICENSE).
