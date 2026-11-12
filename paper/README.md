# Paper datasets

The real inputs of the paper's benchmarks and concordance runs (D95). CI never
fetches them: no workflow runs `cargo xtask paper` or `fetch-datasets.sh`, and
no nightly case uses them (apart from `1kg-30x-chr20`, which nightly already
used).

```sh
docker build -t gloci-bench:p1 -f xtask/docker/Dockerfile.bench xtask/docker
cargo xtask paper fetch --record   # first time: also writes the derived sha256s
cargo xtask paper fetch            # afterwards: verifies everything
sh paper/fetch-datasets.sh [NAME]  # the same without cargo
```

The script needs bash (derived commands run under `bash -o pipefail`), curl and
the tools of `Dockerfile.bench` at the pinned versions. `NAME` builds one item
and what it needs; an unknown name is an error. While a derived sha256 is
empty the script exits 1 and lists the inputs as UNVERIFIED; set
`GLOCI_ALLOW_UNPINNED=1` only when recording.

Downloads total about 6.7 GB (gnomAD chr22 is 5.1 GB, the 1000 Genomes chr20
panel 0.7 GB, the K562 BAM and signal 0.4 and 0.5 GB); the derived inputs take
about as much again. Each download is checked against the provider's own md5
(a bad file is moved to `.bad`). A derived input is checked against the sha256
in `xtask/bench/datasets.toml`; an empty one is TO-PIN and `--record` fills it.

**Pinned (Task 15):** every derived input was built twice into different
cache directories and the sha256 of every file (inputs, indexes and the CRAM
`.crai`) was identical. That check found one non-determinism, the CRAM written
with `-@ 4` (htslib picks codecs per slice by trial, so the choice depends on
slice arrival order); it is now written single-threaded. To re-check after
editing a command, derive into two directories and compare.

`paper/fetch-datasets.sh` is generated: edit `xtask/bench/datasets.toml`, then
run `cargo xtask paper fetch-script`. A test fails when it is stale.

## Concordance (D97)

```sh
cargo xtask paper run --concordance-only [--only CASE] --data-dir DIR --out DIR
```

runs each case of `xtask/paper/concordance.toml` in the bench image (cache
mounted read-only): the reference command once, gloci in the default mode
(`--chrom-order lex` added) and, for cases with `compat` flags, in the compat
mode. Outputs are matched line by line as multisets (`exact`; `numeric`,
where numbers agree within the decimals the less precise side printed or a
relative 1e-5; or `vcf-records` without `#` lines), and the order is checked
too. `<out>/<version>/concordance.csv` gets a row per run (rows of
cases not run are kept), and `<out>/<version>/concordance/<case>.<mode>.txt`
the commands, counts and the first 10 unmatched lines of each side. A run
that differs passes only when the case names the documented difference
(`docs/differences.md`); `anchor unused` means it did not show on this data.
Reference outputs are cached in `target/concordance-work` while their command
and the image stay the same; delete that directory to run them again.

## Providers: licence and citation

- ENCODE (cCREs, CTCF peaks, signal, BAM): open access, no restrictions;
  cite the ENCODE Project Consortium and the accessions in `datasets.toml`.
- GENCODE v46: open; cite Frankish et al., GENCODE.
- gnomAD v4.1: released under CC0 1.0; cite Chen et al. 2024 (gnomAD v4).
- 1000 Genomes 30x (NYGC): open; cite Byrska-Bishop et al. 2022, Cell.
