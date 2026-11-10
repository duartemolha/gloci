# gloci P2a follow-ups

P2a ("breadth") merged into `main` on 2026-11-10 at 4fe2c31. These items were found during P2a and left open on purpose. Fold each one into the plan named.

**Merged state:** `cargo test --workspace` gave 1081 passed, 0 failed, 21 ignored. fmt and clippy `-D warnings` were clean.

## Residuals with a planned home

| # | Item | Evidence | Proposed home |
|---|---|---|---|
| 1 | `vcf stats` with 5 statistics runs at 0.22× of five vcftools runs. The target is 0.1×. | Release build, synthetic 2,500 samples × 50k records: 9.1 s against 41 s. Shared decode and selective field decode are already in place. | P3: make the stats stage parallel. Do it together with D77's deferred stage-2 pools. |
| 2 | The per-chromosome parallel speedup for plain seekable files is not yet shown. | 2 GB BED intersect on a loaded 4-core WSL host: `-t 1`/`4`/`8` = 22.8/24.2/22.5 s. One cleaner run gave 18.3 s at `-t 1` against 12.1 s at `-t 4`. | Re-time on a quiet machine before quoting numbers. P2b's paper benchmarks run on such a machine. |
| 3 | Single-threaded `view`/`intersect`/`closest` are 3–4% slower in CPU time than at P1. | Min of 10 interleaved runs, 316 MB BED. The drain wrapper, the per-record renamer and the deferred-error path are ruled out. The cause is probably code layout or LTO drift. | Profile with `perf`, or bisect, if single-thread throughput becomes a goal. |
| 4 | noodles' BCF decoder panics (`todo!`) when a GT allele is 63 or higher (int16 GT). | Any bcftools BCF holding such a call. The panic is caught and reported as an internal error, exit 1. | P3a: it touches the BCF reader next to BAM/CRAM decoding. Decode GT natively, as P2a now encodes it. |

## Minor gaps (documented or harmless)

- `END=.`: the BCF reader accepts it and prints it, but gloci's VCF text reader rejects it, while bcftools accepts both.
- Malformed INFO number tokens differ from bcftools: a bare sign, a leading space or hex in a Float, `nanx`/`infinityx`, and a Flag given a value.
- BCF→VCF percent-encodes `%` in string values. bcftools prints the raw bytes.
- A BCF record whose GT is not the first FORMAT key can't be printed as VCF (a limit of the noodles writer). The VCF spec forbids this order anyway.
- A reserved FORMAT key with a non-spec type (e.g. `HQ` as Float) is rejected with a message that names no line.
- The nesting-limit message says "more than 128 levels", but the effective limit is 127.
- bedmap joins multiple values with `;` and gloci joins them with `,`. The `--echo-map-id` migration rows don't say so.
- The run-level "share no chromosome names" check misses two warnings that `-t 1` gives. One case is an indexed A whose `--region` is on another chromosome. The other is a chromosome whose A lines are all invalid under `--skip-invalid`.
- Task 1 leftovers:
  - "record N of CHROM" is a record number, and the doc comment could say so.
  - `/proc/<pid>/task/<tid>/fd/N` is still treated as reopenable.
  - `/dev/stdin` redirected from a regular file is spooled, an extra copy.
  - The non-identical failing allele pair from the P1 bug has no test of its own.
- `quantile:-0` is accepted; it is harmless.

## Owner decisions recorded

- MIT licence. The copyright line defaults to Duarte Molha. Confirm the holder before the first public release.
- Nothing is published to PyPI, crates.io, GHCR or bioconda until the owner enables it (P2b D80).
- Build order after P2a: P2b Task 1, then P3a, then P3b, then the rest of P2b (spec §10).
- `groupby -g` names the group columns, as in bedtools (spec §5, global flags).
