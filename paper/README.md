# Paper datasets

The real inputs of the paper's benchmarks and concordance runs (D95). CI never
fetches them: no workflow runs `cargo xtask paper` or `fetch-datasets.sh`, and
no nightly case uses them (apart from `1kg-30x-chr20`, which nightly already
used).

```sh
docker build -t gloci-bench:p1 -f xtask/docker/Dockerfile.bench xtask/docker
cargo xtask paper fetch --record   # first time: also writes the derived sha256s
cargo xtask paper fetch            # afterwards: verifies everything
sh paper/fetch-datasets.sh         # the same without cargo (needs the tools of Dockerfile.bench)
```

Downloads total about 6.2 GB (gnomAD chr22 is 5.1 GB, the 1000 Genomes chr20
panel 0.7 GB); the derived inputs take about as much again. Each download is
checked against the provider's own md5. A derived input is checked against the
sha256 in `xtask/bench/datasets.toml`; an empty one is TO-PIN and `--record`
fills it.

`paper/fetch-datasets.sh` is generated: edit `xtask/bench/datasets.toml`, then
run `cargo xtask paper fetch-script`. A test fails when it is stale.
