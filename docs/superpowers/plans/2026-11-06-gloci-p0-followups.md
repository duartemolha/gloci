# gloci P0 follow-ups

These items came out of the P0 task reviews and the final whole-branch review. The final reviewer triaged all of them as safe to defer past the P0 merge. The P1 plan should start from this list.

## Settle before P1

- **Point vs interval thresholds.** A point A inside interval B passes `--min-bases` but fails `-F`/`-r`. The Task 3 ruling evaluates a threshold on the non-point side normally. Decide this together with D5 before P1 intersect. The code change is one line in `OverlapCriteria::accepts`.
- **`-m` flag collision.** In the spec, `-m` means `--merge-sorted` for `sort` (§7.4) but `--mode` for `intersect` (§7.1). That breaks "same flag, same meaning".
- **Tie order vs bedtools.** gloci breaks ties on end, then input order (spec §7.4). bedtools and sort-bed keep input order. The difference is documented in `docs/differences.md#sort-tie-order`. Confirm that this default is wanted.
- **Run the nightly workflow once** with `workflow_dispatch` right after merge. It covers the x86_64 bench (the P0 performance gate), the differential cases, fuzzing and the large-file tests. Scheduled and dispatched workflows only exist on the default branch. Call P0 done only once this run is green.
- **Deferred from the final review:**
  - `sort` opens every input at once, so more than about 1,000 files hits the open-file limit;
  - merge buffers are left out of the `--max-mem` accounting;
  - no RSS assertions for `check` and `index`;
  - no "did you mean `-g hg38`" suggestion when `-g GRCh38` is used on chr-named files.

## Deferred minors (from review ledger)

The final fix wave (20b3182..b1cb2c5) already fixed these: the `rust-version` bump, the xtask `genomes.rs` doc line, the `--write-index` order hint, the multi-VCF header check, the BCF encoder tests and the `Sorter::push` assert. Everything else below is still open.

- Task 1: minor (deferred): rust-version 1.85 not enforced by any MSRV job; lockfile deps may need newer rustc.
- Task 1: minor (deferred): io→core, ops→io and gloci→testkit edges unused today (plan Interfaces mandate the chain; kept).
- Task 1: minor (deferred): empty [dependencies] sections with trailing blank lines in gloci-core/testkit Cargo.toml.
- Task 1: minor (deferred): ci.yml readelf|grep -q relies on no pipefail; safer to write readelf output to a file first.
- Task 2: minor (deferred): region.rs:601 comma stripping accepts malformed groups (chr1:1,,0-5).
- Task 2: minor (deferred): region.rs:555 `{a}b` reports Unterminated (misleading); names containing `}` inexpressible.
- Task 2: minor (deferred): genome.rs:366 `@SQ` prefix check matches `@SQX`; parse_sizes accepts empty names, no `#` comment skip.
- Task 2: minor (deferred): chrom.rs:177 unchecked `as u32` cast for ids.
- Task 2: minor (deferred): coords.rs test gaps (brace errors, chr1:, chr1:5-10-20, ranks permutation at scale, closed+GenomeFile label).
- Task 2: minor (deferred): closed() ignores OrderSource label for ordering — later tasks must pass GenomeFile/VcfHeader only.
- Task 3: minor (deferred): reciprocal with only min_frac_b set gives A no threshold (untested/undocumented).
- Task 3: minor (deferred): Fraction derived PartialEq distinguishes "1" vs "1.0".
- Task 3: minor (deferred): Fraction::met recomputes 10^digits per call on hot path; store scale.
- Task 3: minor (deferred): test gaps (either with one threshold, reciprocal fallback, parse .5/TooPrecise, 18-digit met near u64::MAX, strand window on Plus, inclusive in accepts, signed_distance with points).
- Task 4: minor (deferred): footprint.rs:235-291 unchecked POS arithmetic (wraps near u64::MAX in release); SVLEN overflow reported as BadSvlen (misleading).
- Task 4: minor (deferred): sv_end computed twice per SV ALT; footprint/alt_footprints duplicate preamble.
- Task 4: minor (deferred): unterminated `<DEL` ALT treated as a sequence instead of an error.
- Task 4: minor (deferred): test gaps (ref-mode multi-SV max end, `.G`, lowercase ALT, mixed hull, alt_footprints with SV, huge END).
- Task 5: minor (deferred): order.rs tests — no isolated equal-key test, ChromBackwards payload ids not asserted.
- Task 6: minor (deferred): open_input pre-classifies compression from 18 bytes; BGZF with non-first BC subfield / XLEN>6 skips the D30 EOF check.
- Task 6: minor (deferred): detect_compression ignores CM byte (head[2]==8).
- Task 6: minor (deferred): detect_format Err string carries literal "{path}" placeholder (trap for other callers).
- Task 6: minor (deferred): test gaps (short inputs, gzip magic only, BGZF < 28 bytes, BCF detection, CRLF first line, `-` stdin).
- Task 6: minor (deferred): process — implementer did not capture a true RED run.
- Task 7: minor (deferred): finish() returns Err if stale-index removal fails after publish (consider warning).
- Task 7: minor (deferred): -o a.vcf.gz -O vcf rejected (compares format+bgzf); message doesn't mention compression. Strict reading kept.
- Task 7: minor (deferred): install_cleanup_handlers swallows ctrlc set_handler errors; no parent-dir fsync.
- Task 7: minor (deferred): test gaps (drop unfinished Output leaves no temp; levels 0/12; stdout BGZF path).
- Task 8: minor (deferred): record.rs:568 final line ending in `\r` without `\n` keeps `\r` (D15 edge).
- Task 8: minor (deferred): read_batch(max=0) returns Ok(true) forever — clamp or debug_assert.
- Task 8: minor (deferred): ">4 GiB line" message lacks where/how-to-fix; checked only after buffering.
- Task 8: minor (deferred): GFF column-count error lacks the spaces hint BED has.
- Task 8: minor (deferred): format detection rejects malformed first data line before line-level parse errors (user sees Detect error, not line error).
- Task 8: minor (deferred): test gaps (skip-policy loc after header scan, GFF CRLF, bed12 overflow, `\r`-only line).
- Task 9: minor (deferred): no tests for with_command_line with existing glociCommand lines, contig extra attrs / ID-after-length, empty-arg quoting.
- Task 9: minor (deferred): vcf.rs:491 allocates Vec<&[u8]> per record (use fixed array).
- Task 9: minor (deferred): `#CHROM` prefix check matches `#CHROMX`; quoted contig values with commas unsupported.
- Task 10: minor (deferred): bcf.rs:355 render treats zero-count typed values as present (Flag-typed END).
- Task 10: minor (deferred): IDX cap raw.len()+65536 could reject pruned headers.
- Task 10: minor (deferred): dicts() unreachable!() per record; store dicts in reader.
- Task 11: minor (deferred): genome.rs line numbers off by leading blank lines; VCF header errors in load_genome always "line 1"; BAM name NUL not verified; hint paths not shell-quoted; header entries computed twice; test gaps (malformed BAM bounds, no-contig load_genome error, >120-char / non-ASCII OrderError, gzipped .fai, CRLF).
- Task 12: minor (deferred): --write-index order-violation hint names the unpublished output path.
- Task 12: minor (deferred): coverage gaps (GFF/GTF index spans, header-only VCF index load+query, header-only BCF, "has no index" error, skip-invalid offsets path).
- Task 12: minor (deferred): region-reader parse errors (index.rs:1053) lack a how-to-fix suffix.
- Task 12: minor (deferred): index.rs is 1344 lines — split into index/build.rs + index/query.rs.
- Task 12: minor (deferred): sparse Linear maintained during TBI pass unnecessarily.
- Task 13: minor (deferred): help-examples test should skip the auto `help` subcommand and count only lines in the Examples block; -q/-v conflict untested and usage text odd; logging eprintln! panics if stderr closed.
- Task 14: minor (deferred): merge_sorted memory grows with input count (256 KiB/input); corrupt spill length trusted for allocation; truncated LZ4 at block boundary untested; no Dedup+large-cascade test; tests lexicographic only.
- Task 14b: minor (deferred): sort_hint unused `_chrom_backwards` param (drop when Task 15 touches sorted.rs); karyotype proptest generator rarely hits chr/X/Y/MT/digit-letter classes; gloci crate depends on gloci-core only for ChromOrder (re-export from gloci-io instead); stray blank line in gloci-core Cargo.toml.
- Task 14b: minor (deferred): migration.md:18 table header `sort -k1,1 -k2,2n` lacks LC_ALL=C (prose covers it).
- Task 15: minor (deferred): `--tmp-dir ""` treated as unset; --check ignores --unique silently; multi-VCF later-header definitions dropped without warning; --check per-input not documented; writer.rs "output" placeholder relies on callers' name_output (pass dest name in WriterOptions); usage_errors bad --tmp-dir asserts only .failure(); 50 MB test ~25 s in debug (consider release/nightly).
- Task 15: minor (deferred): BCF header mismatch reported via IoError::Detect (repurposed variant).
- Task 16: minor (deferred): duplicate wo.write_header assignment; --no-header with --write-index on variant formats should be rejected; --vcf-span exposed on view without effect until Task 17; noodles-rejected headers can't convert to BCF (add to differences.md, Task 20); document view keeps BED headers by default vs sort's opt-in --header.
- Task 17: minor (deferred): block scan parses up to a whole batch past the block (errors from unrelated lines; -g unknown later chrom) — bound with a line-aware adapter; fallback hint misleading when chromosome simply absent/empty file; no-fix-command hint for non-regular files (/dev/fd); index-path order errors cite "line 0"; index read twice; fresh 8 KiB BufReader per probe; Selection copies records twice; Diagnostics private bookkeeping in public fields; append_matching belongs in gloci-io/view.rs; one missing chromosome turns a multi-region query into a full scan (document).
- Task 17: minor (deferred): --assume-sorted silently ignored where it doesn't apply (indexed/compressed/non-BED/stdin) — optional warning.
- Task 18: minor (deferred): `value ""` wording for missing bedGraph value; help should say records = parsed records; hint/lex-note logic duplicated with OrderGuard (shared helper); 64-bit hash duplicate false positives (warning only); no BCF/stdin check tests.
- Task 19: minor (deferred): gen.rs run_args duplicates clap defaults; log_uniform uses libm (cross-platform last-ulp risk — add golden checksum); large test disk peak ~4 GiB > 3 GB ruling (GLOCI_LARGE_BYTES not used in nightly), "2 GiB .bed.gz" read as uncompressed; musl static binary in CI bench may bias gate vs glibc tools; weak deep-nesting test.
- Task 19: minor (deferred): after the 1/16 jump the arena/entries split is fixed from the early mix — skewed inputs may close chunks early (deterministic; untested); coreutils unpinned in Dockerfile.
- Task 19: minor (deferred): chunk.rs:90 `raw.len() as u32` (unreachable: readers reject >4 GiB records).
- Task 20: minor (deferred): exit-status pins only exit codes (not error text); PASS-insertion half of bcftools header claim unpinned; "any fully nested feature" bedextract claim unpinned (only long-feature miss shown); no right-edge point or zero-length query region case for BED; --update-golden can bless gloci regressions into expected_gloci (visible in review diffs); x86_64 nightly unverified.
- Task 21: minor (deferred): persist/cache fuzz corpus between nights; cache cargo-fuzz install.
- Task 21: minor (deferred): nightly `cp -n` under bash -e could fail if corpus pre-exists (add `|| true`).
- Task 22: minor (deferred): check determinism case trivially satisfied (empty stdout); index determinism only .tbi on VCF; kill test fires during spill not final write; kill test heavy (~10 s, 100 MB); broken pipe only for view plain BED.
- Task 22: minor (deferred): optional-value options (num_args 0..=1) would be treated as always taking a value in strip_thread_options (none known); -q/-v remain in provenance (VCF bytes differ by -q/-v, allowed).
- Task 23: minor (deferred): xtask/src/genomes.rs:1 first doc-comment line accidentally deleted (dangling fragment); BuiltinGenome::len/total_length re-parse per call; xtask by_id pub but unused.
- Final: minor (deferred): direct-write open lacks O_TRUNC (stale tail possible with `3<>file`); indexed .vcf.gz→BCF could detect undeclared contigs before the header using index names; -g order errors carry no "add ##contig" hint (names the genome instead).
