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

**Before pinning (Task 15):** derive every input twice into different cache
directories (`--data-dir A`, then `--data-dir B`) and compare the digests; pin
only inputs that agree. The commands are written to be byte-deterministic
(`--no-version`, `--no-PG`, fixed CRAM file name), but only a second run proves it.

`paper/fetch-datasets.sh` is generated: edit `xtask/bench/datasets.toml`, then
run `cargo xtask paper fetch-script`. A test fails when it is stale.

## Providers: licence and citation

- ENCODE (cCREs, CTCF peaks, signal, BAM): open access, no restrictions;
  cite the ENCODE Project Consortium and the accessions in `datasets.toml`.
- GENCODE v46: open; cite Frankish et al., GENCODE.
- gnomAD v4.1: released under CC0 1.0; cite Chen et al. 2024 (gnomAD v4).
- 1000 Genomes 30x (NYGC): open; cite Byrska-Bishop et al. 2022, Cell.
