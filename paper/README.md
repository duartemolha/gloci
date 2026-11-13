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
mode. `exact` and `vcf-records` (without `#` lines) outputs are matched line
by line as multisets, and the order is checked too. `numeric` outputs are
compared by position (D97): line i of gloci with line i of the reference, the
same text fields and every number within the diff harness's tolerance
(relative 1e-5, absolute below 1e-9, integer literals exactly), or within half
a unit of the last decimal for a reference that prints a fixed number of them
(`ref_decimals`). `<out>/<version>/concordance.csv` gets a row per run (rows of
cases not run are kept), and `<out>/<version>/concordance/<case>.<mode>.txt`
the commands, counts and the first 10 unmatched lines of each side. A run
that differs passes only when the case names the documented difference
(`docs/differences.md`); `anchor unused` means it did not show on this data.
Reference outputs are cached in `target/concordance-work`, keyed by the image
Id (or, with `--no-docker`, the host tool's `--version`), the pinned digest of
every input in `datasets.toml` and the command; a change to any of them runs
the reference again. Delete that directory to run them all again.

## Audit of the benchmark commands (2026-11-13)

The runner times the commands; it does not compare what they write. On
2026-11-13 every command of the sixteen cases the README reports was run once
in the bench image, its output kept, and the outputs of each case compared.
All exited 0 with nothing on standard error (vcftools prints its banner).

| Case | Outputs | Finding |
|---|---|---|
| sort-exons | the same lines from all four, each in (chromosome, start) order; ties differ (`differences.md#sort-tie-order`) | gloci had no `-t` and used every core in a case labelled one thread, and GNU sort had `--parallel=8`; both now ask for one |
| merge-exons | identical: gloci, bedtools, bedops | none |
| intersect-ccre-exons | identical: gloci, bedtools, bedops | none |
| intersect-ccre-exons-bgzf | identical: gloci, bedtools, bedops | by design: gloci has 8 threads, bedtools has one, and bedops is timed with the decompression of both inputs to files |
| subtract-exons-ccre | gloci = bedtools; bedops writes the merged difference | `gloci-merged` (`merge` piped into `subtract`) gives bedops' bytes and is what the bedops target now names |
| closest-peaks-exons | gloci = bedtools; closest-features writes one neighbour on each side | `gloci-each-side` does that work (its rows differ where an exon overlaps the peak, `differences.md#closest-each-side`); the bedops target now names it |
| map-ccre-signal | the same 2,348,854 means from all three | bedmap reads column 5: given the four-column bedGraph it printed `NAN` for every row but the first, in about a second, without reading the file. It now reads `ctcf-fc.bed5` |
| coverage-ccre-exons | the same counts, covered bases and lengths; gloci prints the fraction in full and bedtools to seven decimals | none |
| genomecov-exons | the same rows; the chromosome order differs (`differences.md#genomecov-chrom-order`) | none |
| partition-3 | gloci = bedops; multiinter joins segments and names the files | bedops needs sort-bed order and `exons.bed` was sorted by start alone: bedops wrote 42,593 rows with start past end. `exons.bed` is now sorted by start and end. `gloci-merge-adjacent` gives multiinter's first five columns byte for byte and is what the bedtools target names |
| flagstat-ctcf-bam, flagstat-ctcf-cram | identical: gloci, samtools | none |
| bedcov-ctcf-peaks | identical: gloci, samtools | none |
| depth-ctcf-bam, depth-ctcf-cram | identical: gloci, samtools (199,427,781 lines) | mosdepth uses 8 threads and writes a compressed per-base file: a yardstick, not the same work |
| vcf-filter-1kg-bcf | the same 370,398 records from all three; vcftools upper-cases five symbolic ALTs, and gloci writes INFO numbers and one string differently from bcftools | gloci had no `-t` (every core) and bcftools one thread; both now ask for 8 |

The same sorting fault changed `closest-features`' answer for 4,832 peaks;
`bedops -e`, `-m` and `-d` gave the same output on either order.

Not run in the audit, only read: `vcf-filter-gnomad`, `vcf-compare-1kg`, the
`py-*` cases and the cases of `xtask/bench/cases.toml`. Reading those found:

- `intersect-plain-10M` and `intersect-bgzf-10M`: bedtools had no `-u`, so it
  wrote every overlap where gloci, bedops and the other tools write each A
  once. Fixed.
- Thread counts that differ within a case and are not labelled:
  `sort-plain-10M` (gloci every core, GNU sort 8, sort-bed one), the
  `vcf-filter-*` and `vcf-stats-1kg-chr20` cases (gloci every core, bcftools
  one), `intersect-bam-any` and `coverage-cram-targets` (gloci 8, bedtools
  one). Not changed: they are the nightly targets.

What would keep this true is a step in the runner that runs each command once
into its own file and compares the files by a rule the case states; until then
a new or edited case needs the same check by hand.

## Providers: licence and citation

- ENCODE (cCREs, CTCF peaks, signal, BAM): open access, no restrictions;
  cite the ENCODE Project Consortium and the accessions in `datasets.toml`.
- GENCODE v46: open; cite Frankish et al., GENCODE.
- gnomAD v4.1: released under CC0 1.0; cite Chen et al. 2024 (gnomAD v4).
- 1000 Genomes 30x (NYGC): open; cite Byrska-Bishop et al. 2022, Cell.

## Results, tables and the Zenodo bundle (D98)

On a machine near the reference one (§9.1: 8-core x86_64, NVMe, 32 GB):

```sh
cargo xtask paper fetch
cargo xtask paper run --gloci <release musl gloci> --wheel <x86_64 wheel>
cargo xtask paper                  # tables and numbers from paper/results/<version>
cargo xtask paper bundle -o gloci-paper-evidence.tar.gz
```

`paper run` writes `paper/results/<version>/`: `bench/` (hyperfine JSON and
`<case>.rss.json` per case), `concordance.csv` (with `concordance/`),
`environment.json` and `repo-snapshot.json`.

- `environment.json` records the host (CPU model, cores, RAM, kernel,
  architecture), the bench image id and repo digests, every tool's
  `--version` (run in the image), the gloci version and commit, the wheel's
  file name and sha256 (`--wheel` only records it; nothing is run with it) and
  `started_at` / `finished_at`. A value that cannot be read is `null`, never a
  guess; so is the version of a tool whose `--version` fails.
- `repo-snapshot.json` records the repository at that commit: golden cases per
  tool (the `-vs-<tool>` suffix), the differences anchors with their titles
  and golden counts, the migration tables of `docs/migration.md` (status is
  `unavailable` when the gloci cell starts with "not available", `partial`
  when it says so later), test counts, and bytes and records of every dataset
  and derived input.

`cargo xtask paper` writes `paper/tables/T1..T8-*.csv`, `paper/numbers.tex`
(`\glnum{bench.merge-exons.gloci.median_s}`) and `paper/numbers.json`:

| Table | Content |
|---|---|
| T1 | datasets: release, url, provider checksum, bytes, records |
| T2 | every benchmark command: role (subject, baseline = named by a target, else yardstick), median/mean/stddev/min/max seconds, runs, peak RSS in MiB |
| T3 | speedup of each tool's median over gloci's, with the target that names it and whether it is met |
| T4 | concordance per case and mode, without the verdict |
| T5 | differences anchors: golden cases and concordance cases that exercise each |
| T6 | the migration tables, one row per old command |
| T7 | test counts |
| T8 | software versions and roles |

Seconds are `%.3f`, speedups `%.2f`, rates `%.6f`, MiB `%.1f`; rows are sorted.
A median from the native timing runner (no hyperfine) has no mean, stddev, min
or max: those cells are empty. `committed_tables_are_current` fails when the
committed tables differ from what the committed results produce (it is skipped
while `paper/results/` is empty).

`cargo xtask paper bundle` writes a `.tar.gz` of the results, the tables and
numbers, `datasets.toml`, the fetch script, both Dockerfiles,
`xtask/paper/*.toml`, `CITATION.cff` and a `README.txt`. It is reproducible
(sorted entries, mtime = `finished_at`, uid/gid 0, modes 0644/0755, gzip mtime
0), so two runs give identical bytes; the deflate stream is the same for one
build of `flate2`, not promised across versions of it.
