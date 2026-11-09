# gloci P1 follow-ups

P1 was merged to `main` at 8b18147. The final whole-branch review and its single fix wave fixed every finding. The re-review of that fix wave then found three new Important issues, listed first below. Under the process they go to the next plan, not a second fix wave. The P2 plan should start with them.

## Open Important issues (fix first in P2)

1. **Output changes with `--threads` on indexed inputs out of table order.**
   - `-t 1` streams, and so refuses an indexed input whose file order differs from the chromosome table (a lex-sorted `.bed.gz`, or a VCF header that lists contigs in a different order from its records).
   - `-t ≥2` uses the indexed path, which silently reorders the output and succeeds.
   - Fix: `run_indexed` checks each input's index reference order against the table and raises the same order error as streaming. Also correct the "same as with one thread" sentence in the five relation `--help` texts.
2. **`sort` keeps every regular-file reader open with its own BGZF pool.**
   - The rule is `keep = paths.len() <= fan_in` (commands/sort.rs).
   - Measured on 100 BGZF inputs: 298 MB and 1,618 threads at `-t 16`, against 95 MB and 49 threads before.
   - Fix: keep or spool only non-regular inputs, and reopen regular files on each pass. Also treat `/dev/fd/*` as non-reopenable, because on macOS reopening it shares the file offset.
3. **The `vcf-filter-1kg-chr20` benchmark is degenerate.**
   - Every record in the dataset has QUAL `.`, so `--min-qual 30` drops everything for every tool.
   - Fix: drop `--min-qual` from that case for all tools, so it times the MAF and missingness pipeline.

## Deliberate P2 work carried from P1

- **Unsorted-mode index (D40).** The real-crate benchmark shows rust-lapper about 2× faster than coitrees, and superintervals about 1.2–1.3×, on both synthetic sets. Decide whether to switch.
- **BCF output with a missing haploid genotype next to diploid calls** is currently a clear error (`#bcf-mixed-ploidy-missing`). Add native GT encoding.
- **`--region`/`--regions-file` on the relation commands** (D52).
- **The mdBook and the bioRxiv external benchmark datasets** (D53).
- **Spec §8.3 parallelism for plain seekable files** ("P2 at the latest"), and the stage-2 parse/format worker pools.
- **The full `map` op set**, beyond the P1 core plus `stdev_pop`/`var_pop`.
- Merge and complement `--sort` sorted-record numbering is now applied in ops. Check that no CLI-only behaviour is left (spec §8.1) before the Python library work.

## Parked minors from the final review

- `run_chunks` can return Ok with fewer results than jobs; it should return Internal. This is unreachable today under `catch_unwind`.
- `vcf filter` `validate()` now runs after the list files are read, so a usage error can be reported as a file error.
- Indexed-mode stderr (skip examples, parse errors) names BGZF offsets, while `-t 1` names line numbers.

## Deferred minors (from the P1 review ledger)

- Task 1: minor (deferred): variant_kind.rs ~L52 `v.pos - 1` underflows for POS 0 (telomere breakends); `(0,0) => Ref` for identical REF/ALT undocumented; D34 test gaps (-r with only -f, -r/either with a point, -F-only point-on-B); redundant `min_frac_a.or(one)`; blocks_overlap_len lacks debug_assert on sorted input; blocks_overlap is O(n·m).
- Task 2: minor (deferred): no region-reader (open_region_reader) strand test; bed.rs/gff.rs still split with slice::split instead of the new split_fields; no empty-column-6 BED test; behaviour::killed_mid_write_leaves_nothing flaked once under load (fold into Task 26 if it recurs).
- Task 3: minor (deferred): query test helper encodes the 1-based no-end rule itself (differential case is the independent check).
- Task 3: minor (deferred): tsv_sort_hint interpolates path/comment char unquoted into a bash-specific ($'\t') shell line; hint test checks substrings only.
- Task 4: minor (deferred): open_tsv_writer's write_index param vs Output::wants_index() two sources of truth (debug_assert or drop); IoError::Write path "output" hides destination (expose Output::dest_name); TsvPending reported as IoError::Detect; tests lack failed-index-no-partial, stale .tbi removal on TSV path, line_skip_count assert, '#' header + index, direct/stdout refusal; skip count assumes every header line ends in \n.
- Task 5: minor (deferred): hits compared as strings (loc not compared; duplicate-line order unchecked); SweepSummary.peak_active is a doc(hidden) pub test hook; testkit↔io dev-dependency cycle; ChromTable::cmp lacks a==b shortcut (~12% read time, karyotype) — fold into Task 26.
- Task 6: minor (deferred): look-ahead max_distance cap is conservative by one base (counts start == end+max); flag ignored under a strand rule (documented).
- Task 7: minor (deferred): padded-end i32 rule sends large-window chromosomes to ITree64 needlessly (speed only); O(n) edge-point walk under --inclusive-points; shuffled-A closest ~4× sorted (Task 27); "sorted A" querent condition comment.
- Task 8: minor (deferred): StrandWatch warns only for records read and lacks a how-to-fix; GFF point rewrite writes start+1 > end; short GFF rows pass through without rewrite; extra long flag names (--window, --min-frac-a/-b, --reciprocal, --either, --same-strand, --opposite-strand) not in spec — record in docs (Task 29).
- Task 8: minor (deferred): bedGraph with ≥12 columns parsed as BED12 under --split; global test logger via set_logger in tests/relate.rs (one per binary).
- Task 9: minor (deferred): validate() runs twice; ColumnWatch Relaxed atomic; placeholder width only from B file 1.
- Task 9: minor (deferred): cases.toml comment claims a "same deletion" pair that is not the same variant; #trimmed-allele-match prose claims (partial ALT lists, bcftools norm agreement) unpinned.
- Task 10: minor (deferred): i128 overflow forces min/max float printing; text ops allocate per value (to_vec) and SipHash sets; COL names containing ':' unaddressable; compare-mode names duplicated in xtask.
- Task 10: minor (deferred): compare() doc comment doesn't mention integer literals must match exactly.
- Task 11: minor (deferred): property test aggregates only 4:collapse; kmerge_with rank refresh skipped on empty refill (edge, no failing case).
- Task 11: minor (deferred): differences.md line over the file's wrap width.
- Task 12: minor (deferred): --split piece arithmetic not in the property test (split:false); naive shares Window::apply/Fraction::met/spans_overlap and D60 semantics.
- Task 13: minor (deferred): validate() twice; ColumnWatch widened to pub(crate); over-long help line; Reporter put/start_row duplicate intersect helpers; reporter allocates a candidate Vec per A (Task 27).
- Task 14: minor (deferred): BED header-name column path untested; no BCF -b test.
- Task 14: re-review: all 4 addressed; minor (deferred): map help says "as bedtools stdev and var do" (bedtools has no var).
- Task 26a: minor (deferred): stress test exempts an empty spill/ dir; SIGTERM-under-lock test's pre-signal read loop has no timeout; report lacks pasted raw output for the 600-run and broken-variant claims; test hook protected by doc only.
- Task 16: minor (deferred): no RED run; spill-test size near but not proven minimal; --sort on mixed formats is a usage error.
- Task 17: minor (deferred): stale "output line" comment in tests/transform.rs:209; no -l/-r short flags; huge-pad CLI test.
- Task 19: minor (deferred): GTF backslash skip also applies outside quotes; #vcf-allele-span convert bullet should say "by default" for --split-alleles; ref-split has no differential pin (unit test only); doc line wrapping.
- Task 20: minor (deferred): setting gt=None directly skips refreshing spare (realloc only); partial alleles left in gt on parse error (pre-existing, caller gets Err).
