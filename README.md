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

## Coming from bedtools, bedops, vcftools or samtools

See the [migration guide](docs/migration.md) for the command mapping, and
[docs/differences.md](docs/differences.md) for the places where gloci
intentionally differs.

## Licence

MIT. See [LICENSE](LICENSE).
