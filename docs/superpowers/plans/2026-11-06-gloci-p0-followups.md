# gloci P0 follow-ups

These items came out of the P0 task reviews and the final whole-branch review of 6 November. On 2026-11-14 the list was checked against the code: items that are fixed or no longer apply were dropped (git history keeps them), and what is still true is listed below with its current place.

## Settle before P1

- **Point vs interval thresholds.** A point A inside interval B passes `--min-bases` but fails `-F`/`-r`. The Task 3 ruling evaluates a threshold on the non-point side normally. Decide this together with D5 before P1 intersect. The code change is one line in `OverlapCriteria::accepts`. Settled: D34 counts every threshold as met when either side is a point (`crates/gloci-core/src/overlap.rs:160`).
- **`-m` flag collision.** In the spec, `-m` means `--merge-sorted` for `sort` (§7.4) but `--mode` for `intersect` (§7.1). That breaks "same flag, same meaning". Settled: `intersect --mode` has no short flag.
- **Tie order vs bedtools.** gloci breaks ties on end, then input order (spec §7.4). bedtools and sort-bed keep input order. The difference is documented in `docs/differences.md#sort-tie-order`. Confirm that this default is wanted. Settled: kept, with `--full-key` to reproduce `sort-bed`, and pinned by golden cases.
- **Run the nightly workflow once** with `workflow_dispatch` right after merge. It covers the x86_64 bench (the P0 performance gate), the differential cases, fuzzing and the large-file tests. Scheduled and dispatched workflows only exist on the default branch. Call P0 done only once this run is green.

## Open on 2026-11-14

### Parsing

- A symbolic ALT without its closing `>` (`<DEL`) is read as a sequence (`crates/gloci-core/src/footprint.rs:68`); bcftools and tabix accept it as is, so refusing it needs a `docs/differences.md` section.
- A gzip header with the BGZF `BC` subfield not first, or other subfields beside it, is called BGZF and then fails to decode, where htslib reads it as gzip (`crates/gloci-io/src/detect.rs:97`).
- `detect_compression` does not check the CM byte; flate2 rejects such a header itself (`crates/gloci-io/src/detect.rs:93`).
- `detect_format` returns messages with a literal `{path}` for the caller to fill; every caller does so today (`crates/gloci-io/src/detect.rs:240`, `:312`).
- A line is buffered whole before the 4 GiB check, so memory is bounded only by the line length (`crates/gloci-io/src/record.rs:563`, `:641`).
- The BCF `IDX` cap (`raw.len() + 65536`) refuses headers whose lines were pruned but kept their high `IDX` values, which htslib reads (`crates/gloci-io/src/bcf.rs:38`).
- Chromosome ids are cast with `as u32` unchecked; reachable only past 2^32 distinct names (`crates/gloci-core/src/chrom.rs:118`).
- A chromosome name containing `}` cannot be written in the `{name}:start-end` form, as in htslib (`crates/gloci-core/src/region.rs:25`).
- `sv_end` runs twice per SV ALT under `--vcf-span allele` (`crates/gloci-core/src/footprint.rs:236`, `:188`); `footprint` and `alt_footprints` repeat their preamble (`:212`, `:224`).
- The VCF reader collects its first eight columns into a `Vec` per record (`crates/gloci-io/src/vcf.rs:295`).
- `header_entries` runs twice for a VCF that sets the chromosome order (`crates/gloci-io/src/genome.rs:264`).
- `BuiltinGenome::len` and `total_length` re-parse the table on each call (`crates/gloci-io/src/builtin.rs:55`).

### Errors and messages

- The GFF column-count error lacks the "line contains spaces" hint BED's has (`crates/gloci-io/src/gff.rs:61`).
- A malformed first data line fails format detection ("pass --format") before the line parser can name the problem (`crates/gloci-io/src/detect.rs:310`).
- `run: gloci ...` hints do not shell-quote the paths they name (`crates/gloci-io/src/sorted.rs:15` and the other hints).
- Index-backed region reads report a bad record without a how-to-fix clause (`crates/gloci-io/src/index/query.rs:483`).
- An unknown chromosome under `-g FILE` gets no how-to-fix clause (`crates/gloci-io/src/naming.rs:116`).
- `check` says `value ""` for a bedGraph line with no value column (`crates/gloci-ops/src/check.rs:478`).
- The logger's `eprintln!` panics when stderr is a pipe whose reader has exited; the panic hook then exits 1 after cleanup (`crates/gloci/src/logging.rs:21`).
- `install_cleanup_handlers` ignores a failure to set the signal handler (`crates/gloci-io/src/atomic.rs:195`).

### Options and documentation

- `--tmp-dir ""` is taken as unset (`crates/gloci/src/embed/mod.rs:214`).
- `sort --check` ignores `--unique` without a word (`crates/gloci/src/commands/sort.rs:155`); its help does not say each input is checked on its own (`crates/gloci/src/cli.rs:3290`).
- `view --no-header --write-index` on VCF writes a headerless `.vcf.gz` and its `.tbi`; only BCF output is refused (`crates/gloci-ops/src/view.rs:666`).
- That `view` keeps BED header lines by default while `sort` needs `--header` is not documented (`crates/gloci/src/cli.rs:3220`).
- `check`'s help does not say that "records" counts the records that parsed (`crates/gloci/src/cli.rs:222`).
- The migration example header `sort -k1,1 -k2,2n` lacks `LC_ALL=C`; the prose above it says it (`docs/migration.md:22`).

### Indexing and region queries

- `--assume-sorted` reads a whole batch past a chromosome's block, so a bad line or a chromosome missing from `-g` just after it fails the query (`crates/gloci-ops/src/bedsearch.rs:327`).
- An order error on the index path cites the emitted-record count as a line number; it can fire only when an index does not match its file, which D26 allows with a warning (`crates/gloci-io/src/index/query.rs:610`). The fix needs a location form other than `Loc::Line` for the order error, and `Loc` is matched exhaustively.
- The `--assume-sorted` fallback warning advises sorting also when the chromosome is simply absent (`crates/gloci-ops/src/view.rs:374`).
- One chromosome bisection cannot find turns a multi-region query into a full scan, and `--assume-sorted` is ignored silently where it does not apply; neither is in `docs/differences.md#view-assume-sorted` (`crates/gloci-ops/src/view.rs:363`).
- The index is read twice for a region query (`crates/gloci-ops/src/view.rs:253`, `:362`); each bisection probe makes a fresh 8 KiB `BufReader` (`crates/gloci-ops/src/bedsearch.rs:41`); selected records are copied twice (`crates/gloci-ops/src/view.rs:508`, `:572`).
- The TBI pass keeps a sparse linear index it does not use (`crates/gloci-io/src/index/build.rs:425`).
- An indexed `.vcf.gz` converted to BCF finds an undeclared contig only at its first record, not from the index names before the header (`crates/gloci-io/src/convert.rs:215`).

### Sort

- A corrupt spill record length is trusted for an allocation of up to 4 GiB (`crates/gloci-ops/src/sort/spill.rs:178`).
- After the 1/16 jump the arena/entries split is fixed from the early mix of record sizes (`crates/gloci-ops/src/sort/chunk.rs:139`); `raw.len() as u32` there is unchecked, though readers refuse records over 4 GiB (`:90`).
- A BCF header mismatch between inputs is reported as `IoError::Detect` (`crates/gloci-ops/src/sort/inputs.rs:450`).

### Output and internals

- No fsync of the parent directory after the rename that publishes an output (`crates/gloci-io/src/atomic.rs:282`).
- Writers report errors under the placeholder path `output` and rely on callers to rename it (`crates/gloci-io/src/writer.rs:166`, `crates/gloci-io/src/error.rs:175`).
- `Diagnostics` keeps its bookkeeping in public fields (`crates/gloci-io/src/diag.rs:20`).
- `check` repeats OrderGuard's hint and lex-note logic (`crates/gloci-ops/src/check.rs:571`) and finds duplicates by a 64-bit hash, which can warn falsely (`:312`).
- `strip_thread_options` reads an option with an optional value (`vcf norm --split`) as taking the next argument; no recorded command line shows it today (`crates/gloci/src/lib.rs:55`).

### Tests

- Region strings: `chr1:`, `chr1:5-10-20` and a ranks permutation at scale (`crates/gloci-core/tests/coords.rs`).
- Overlap: `either` with one threshold, `Fraction::parse(".5")` and `TooPrecise`, 18-digit `met` near `u64::MAX`, a strand-aware window on Plus, inclusive points in `accepts`, `signed_distance` with points (`crates/gloci-core/tests/overlap.rs`, `distance.rs`).
- Footprints: ref-mode multi-SV max end, `.G`, a hull of SEQ and SV ALTs, `alt_footprints` with an SV, a huge END (`crates/gloci-core/tests/footprint.rs`).
- `OrderViolation::ChromBackwards` ids are never asserted (`crates/gloci-core/tests/order.rs:13`).
- Input: gzip magic alone, BGZF shorter than 28 bytes, a CRLF first line (`crates/gloci-io/tests/input.rs`).
- Output: dropping an unfinished `Output` leaves no temp file; compression levels 0 and 12 (`crates/gloci-io/tests/output.rs`).
- Text readers: skip-policy locations after a header, GFF CRLF, BED12 block overflow, a `\r`-only line (`crates/gloci-io/tests/text_readers.rs`).
- VCF header: `with_command_line` over an existing `##glociCommand` line, contig lines with extra fields or `ID` after `length`, quoting an empty argument (`crates/gloci-io/tests/vcf.rs`).
- Genome loading: malformed BAM bounds, a VCF without contigs, a gzipped `.fai`, CRLF; an OrderError over 120 characters or non-ASCII (`crates/gloci-io/tests/genome_order.rs`).
- Index: GFF/GTF spans, a header-only VCF queried, a header-only BCF, the "has no index" error, `--skip-invalid` on the offsets path (`crates/gloci-io/tests/index.rs`).
- Spill runs truncated at an LZ4 block boundary (`crates/gloci-ops/src/sort/spill.rs`); the sort engine tests use only lexicographic order (`crates/gloci-ops/tests/sort_engine.rs`); the karyotype proptest generator rarely makes chr/X/Y/MT or digit-letter names (`crates/gloci-core/tests/coords.rs:327`).
- `-q` with `-v` is untested (`crates/gloci/src/cli.rs:68`); the bad `--tmp-dir` test asserts only failure (`crates/gloci/tests/sort.rs:478`).
- Thread determinism of index bytes covers only `.tbi` on VCF (`crates/gloci/tests/behaviour.rs:133`); the kill tests fire during the spill, never during the final write (`crates/gloci/tests/behaviour.rs:281`).
- Differential cases: exit-status cases pin no error text (`xtask/diff/cases.toml:200`); the PASS-insertion half of the `bcftools sort` header claim is unpinned (`cases.toml:272`); `--update-golden` can bless a regression into `expected_gloci` (`xtask/src/diff.rs:227`).
- `xtask gen`: no golden digest guards `log_uniform` against last-ulp differences across platforms (`xtask/src/gen.rs:62`); the deep-nesting test does not check containment (`xtask/tests/gen.rs:76`).

### CI and tooling

- The 50 MB sort test runs in normal CI (~25 s in debug) (`crates/gloci/tests/sort.rs:401`), and the kill test writes 100 MB (`crates/gloci/tests/behaviour.rs:281`).
- The nightly large-file job sets no `GLOCI_LARGE_BYTES`, so its disk peak is not bounded to the 3 GB ruling (`.github/workflows/nightly.yml:11`).
- cargo-fuzz is installed every night instead of cached (`.github/workflows/nightly.yml:155`).
- coreutils is unpinned in the reference-tools image (`xtask/docker/Dockerfile:5`).
- `xtask gen`'s argument parser repeats the clap defaults (`xtask/src/gen.rs:1133`, `xtask/src/main.rs:264`); `genomes::by_id` is public and unused (`xtask/src/genomes.rs:313`).
