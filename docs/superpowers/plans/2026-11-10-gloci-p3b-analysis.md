# gloci P3b Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the analysis half of spec Phase 3:
- `stats jaccard`, `stats fisher` and `stats reldist`;
- `random`, `sample` and `shuffle` on a pinned, portable random stream;
- BEDPE input, `pair-intersect` (bedtools `pairtobed`), `pair-pair` (`pairtopair`) and `convert --bedpe` (bedtools `bamtobed -bedpe`/`-mate1`, moved here from P3a);
- `vcf popgen`: nucleotide diversity, Tajima's D, Weir and Cockerham Fst, LD and relatedness;
- `vcf compare`: site overlap, sample overlap and genotype discordance;
- a parallel selection stage for every `vcf` command, which closes P2a follow-up 1.

**Architecture:** No new crates.
- Interval statistics are reporting layers over the existing sweep (`run_relation`), fed by merged-stream adapters. `pair-*` query a new in-memory `IntervalIndex` in gloci-engine (the D61 lapper rule) per BEDPE end.
- Randomness is one module, `gloci-ops/src/rng.rs`: ChaCha8 keyed by `--seed`, one stream per unit of work, so output never depends on threads.
- `vcf popgen` statistics implement P2a's `Stat` trait and run through the `vcf stats` engine in one pass. `vcf compare` walks two selected streams in the shared chromosome table's order.
- `run_selection` decodes and filters batches on worker threads and hands records to its sink in input order.

**Tech Stack:** Everything P3c uses, plus `libm` 0.2 (pure-Rust `lgamma`/`exp`, identical on every platform) and `rand_chacha` 0.10 as a normal dependency of gloci-ops (already locked for xtask). No reference tool or Docker change: every bedtools 2.31.1 and vcftools 0.1.16 option used here exists in `gloci-diff:p0` (checked against their sources; Task 1 re-checks in the image).

**Spec:** `docs/superpowers/specs/2026-11-06-gloci-design.md`: §2, §4.1–§4.6, §5.1–§5.2, §6.1–§6.4, §7.3, §7.6, §7.7, §7.8 (`vcf popgen`, `vcf compare`), §8.1–§8.4, §9, §10 (P3 row, Order note) and §11.

**Builds on:** P3a (D99–D124) and P3c (D125–D141), fully merged; D1–D79 and D82 stay in force; P2b's D80–D98 are untouched. In particular: D50/D71/D72 (selection, expressions, the `vcf stats` engine: `Stat`, `StatsEngine`, `SelectStage`, `run_selection`); D36–D39, D42, D45 (sweep, criteria, `--sort`/`--unsorted`, TSV, placeholders); D61 (store `[start, end+1)`, widen queries by one, re-test exactly); D62 (`--chrom-map`); D101, D102, D120 (alignment spans, read filters, `AlnSource`).

**Out of scope:** P2b Tasks 2+ (they run after this plan); D77's BED parse pools; everything D142 lists as not planned.

**P3b is done when:**
- every bedtools- and vcftools-facing behaviour in D146–D163 is pinned by committed goldens replayed by plain `cargo test`;
- every intentional difference is anchored in `docs/differences.md`, pinned, and linked from `docs/migration.md`;
- the §11 rows for `jaccard`, `fisher`, `reldist`, `shuffle`, `random`, `sample`, `pairtobed`, `pairtopair`, `bamtobed -bedpe`, `--site-pi`, `--window-pi`, `--TajimaD`, `--weir-fst-pop`, `--*-r2`, `--relatedness*`, `--diff*` and `vcf-compare` map to working commands;
- the nightly bench meets Task 17's targets, including `vcf stats` with 5 statistics at ≤ 0.1× five vcftools runs.

## Global Constraints

P3c's Global Constraints all apply: u64 half-open coordinates; `rust-version = "1.91"`; no htslib; exit codes 0/1/2 with errors that say what, where and how to fix; atomic output; content detection; record streams in gloci-ops; 2–4 `--help` examples per command; CI on three platforms; output byte-identical at `-t 1/2/8`, tested for every new path; nothing quadratic in the contig count.

New in P3b, binding on every task:
- **Goldens.** Single-output commands get `xtask/diff/cases.toml` cases whose gloci side is one plain invocation, committed under `crates/gloci/tests/golden/`. Multi-file outputs (`vcf popgen`, `vcf compare`) follow `vcf stats`: vcftools' files are committed under `crates/gloci-vcf/tests/golden/{popgen,compare}/` by an append-only `regen.sh` in each, compared byte for byte by `crates/gloci/tests/vcf_popgen.rs` and `vcf_compare.rs`, and mirrored by live `&& cat P.ext` cases.
- **Goldens discriminate.** Each task names the trap its fixtures set: an input where a plausible wrong rule gives different output.
- **Differences.** The task that pins a difference writes its anchored section; Task 19 links each from `docs/migration.md`.
- **Randomness** follows D143–D144: no platform RNG, no `rand` distribution, nothing that depends on scheduling.
- **Statistics** follow D145: the reference formula in its operation order and precision, printed as the reference prints it, with an oracle test written from the published formula.
- **Memory bounds,** each asserted in its task:

  | Path | Bound |
  |---|---|
  | `stats jaccard`/`fisher` | the sweep's active set + one merged interval per input; no per-record vectors |
  | `stats reldist` | 8 B × the B records of one chromosome (`--unsorted`: of all chromosomes) |
  | `random` | ≤ 2 × threads output batches |
  | `sample -n N` | N × (record bytes + 24 B) |
  | `shuffle` | 32 B × allowed intervals; `--no-overlapping` adds ≈ 48 B per placed record |
  | `pair-intersect`/`pair-pair` | B's raw bytes + 32 B per indexed end; A in ≤ 2 × threads batches |
  | `convert --bedpe` | the records of one read name |
  | parallel selection (D157) | ≤ 2 × workers batches of ≤ 524,288 decoded genotypes (≈ 24 B each) |
  | windowed popgen statistics | ⌈size/step⌉ + 1 open windows |
  | LD | held sites × samples × 4 bits ≤ `--max-mem` (default 2G), else exit 1 |
  | `--relatedness2` / `--relatedness` | 8 B / 12 B × N(N+1)/2 + a 64-site block, N = selected samples |
  | `vcf compare` | the records of one position in both files + 24 B per sample of the union |
- **Contig count.** Per-chromosome state lives in vectors indexed by `ChromId`; genome projections use prefix sums and binary search. 200,000-contig tests guard `random`, `shuffle`, `pair-intersect` and `vcf compare`.

## Decisions this plan adds

| # | Decision |
|---|---|
| D142 | **Scope and names.** <ul><li>New: `stats jaccard\|fisher\|reldist`, `random`, `sample`, `shuffle`, `pair-intersect`, `pair-pair`, `convert --bedpe`, `vcf popgen`, `vcf compare`.</li><li>Renamed where a bedtools letter already means something in gloci: fisher `-m` → `--merge-overlaps`; shuffle `-chrom` → `--same-chrom` (`--chrom` is a `vcf` filter), `-chromFirst` → `--chrom-first`, `-excl`/`-incl` → `--exclude`/`--include` (§7.3), `-maxTries` → `--max-tries`, `-noOverlapping` → `--no-overlapping`; random `-l`/`-n` → `-l/--length`, `-n/--count` (as `windows`); sample `-s forward\|reverse` → `--strand +\|-` (as `genomecov`); `-type` → `--type`; pairtopair `-slop N` → `-w N`, `-ss -slop N` → `--left N --strand-window`, its default strand enforcement → `-s`, `-rdn` → `--require-different-names`; bamtobed `-bedpe`/`-mate1`/`-ed` → `convert --bedpe`/`--mate1-first`/`--score-tag NM`. vcftools spellings are visible aliases (`TajimaD`, `weir-fst-pop`, `diff-site`, …).</li><li>`--include`/`--exclude` take a BED in `shuffle`/`fisher` (§7.3) and an expression in `vcf` commands (§7.8); no command has both. Raised with the owner as a third exception beside `groupby -g` and alignment `sort -n`.</li><li>`sample` (bedtools) and P3c's `samples` (samtools) are both spec names; each help text names the other.</li><li>**Not planned** (migration rows): jaccard/fisher `-split` (run `convert --blocks` first), jaccard `-S +\|-`; shuffle `-f`, `-allowBeyondChromEnd`, `-incl` strand copying; `sample` of BAM (`samtools view -s`), `-ubam`; pairtobed `-abam`/`-bedpe`/`-ed`/`-ubam` (pipe `convert --bedpe`); vcftools `--geno-chisq`, `--*-r2-positions`, `--hapcount`, `--LROH`, `--pca`, `--mendel`, `--diff-switch-error`; vcf-compare's text report.</li></ul> |
| D143 | **The random stream.** `rand_chacha` 0.10 `ChaCha8Rng`, used only through `next_u64`, `set_stream` and `set_word_pos`. Key = `seed.to_le_bytes()` then 24 zero bytes; stream = the unit of work (D144). Draws are gloci's own: `below(n)` = Lemire's multiply-shift with rejection (`m = x·n` as u128; if `(m as u64) < n`, `t = n.wrapping_neg() % n` and redraw while `(m as u64) < t`; return `m >> 64`); `unit()` = `(x >> 11) × 2⁻⁵³`; `coin()` = top bit. A known-answer test pins the first outputs; changing them is a breaking change for the release notes. `--seed N` is a u64; without it gloci seeds from `RandomState` and warns `no --seed given; using --seed N (pass it to reproduce this output)`, which is how §7.3's "`--seed` required for reproducibility" is read. |
| D144 | **What randomness guarantees.** Same seed, input and options → the same bytes at every `--threads`, on every platform, in every release unless its notes say the stream changed. Units: `random` interval *i* uses stream *i*; `shuffle` data record *j* uses stream *j* (retries included); `sample` uses stream 0 with record *j*'s key at word 2(*j*−1). Distributions are exact (no modulo or cast bias) and defined in D150–D153. Not guaranteed: bedtools' output for the same seed (`#rng-stream`), or gloci's under other options. Tested by the known-answer test, χ² uniformity at fixed seeds (deterministic), constraint property tests and `-t` determinism. |
| D145 | **Numbers.** <ul><li>The reference's operation order and precision: f32 where bedtools divides in `float` (jaccard); f64 for `double`; exact u128 rationals where bedtools uses x87 `long double` (fisher's estimate; a 1-unit golden difference pins `#fisher-estimate-exact`). vcftools' 32-bit products and positions are computed in u64 (`#popgen-integer-overflow`).</li><li>`lgamma`/`exp` come from `libm`; `sqrt` is IEEE-exact.</li><li>Printing goes through `gloci_core::cfmt`: `write_g(x, prec)` (`%.Pg`; C++ `<<` is precision 6), `write_e6`, `write_fixed(x, prec)` (`%.Nf`; Rust rounds ties to even, as glibc). NaN prints `nan`, infinities `inf`/`-inf`; fixed reference strings (fisher's `-nan`) are printed as such. A golden where the reference prints `-nan` for an arithmetic NaN pins `#nan-sign` (a NaN's sign differs between x86 and Arm; gloci never prints it).</li><li>Each statistic has an oracle test from the published formula (hypergeometric Fisher; Nei–Li π; Tajima 1989; Weir–Cockerham 1984; Pearson r², Lewontin D′; Yang 2010; Manichaikul 2010), never from the port.</li></ul> |
| D146 | **BEDPE** (§5.1). <ul><li>Never detected, since `intersect --mode pairs` of two BED3 files has its shape. `pair-intersect -a`, `pair-pair -a/-b` read BEDPE; `shuffle --format bedpe` opts in; `--format bedpe` elsewhere is exit 2 naming those commands. An input detected as VCF, BCF, GFF or alignments is exit 2 (alignments: `gloci convert x.bam --bedpe`).</li><li>Columns `chrom1 start1 end1 chrom2 start2 end2 [name score strand1 strand2 extra…]`: 6, 7, 8 or ≥ 10 per line (9 is an error, as bedtools); lines may differ. Leading `#`/`track`/`browser` lines are headers; CRLF tolerated.</li><li>chrom `.` = an absent end (`. -1 -1`): integer coordinates, ignored, overlaps nothing. Otherwise non-negative, `start ≤ end`, points allowed. Strands `+ - .`; none or `.` is unknown and never matches `-s`/`-S`.</li><li>Records are written unchanged; `--chrom-map` rewrites both chrom columns; `--skip-invalid` applies.</li><li>Output `.bedpe[.gz]` / `-O bedpe\|bedpe.gz` (BGZF for `.gz`); `--write-index` is exit 2 suggesting `gloci index -s 1 -b 2 -e 3 --zero-based x.bedpe.gz`.</li></ul> |
| D147 | **`stats jaccard -a A -b B`.** Each input is merged as `gloci merge` does (overlapping and bookended; per strand with `-s`). Intersection = Σ overlap bases of merged A × merged B pairs accepted by §4.3 (`-f -F -r -e --min-bases -s -S --inclusive-points`, measured on merged intervals); n_intersections = those pairs; union = Σ merged A + Σ merged B − intersection, both inputs read to the end. Output `intersection\tunion\tjaccard\tn_intersections` and one row, jaccard = `I as f32 / U as f32` printed `%g`. Sorted inputs, order-checked; `--sort`; VCF per `--vcf-span`; points cover no bases but count as intersections. |
| D148 | **`stats fisher -a A -b B -g G`.** <ul><li>Counts, unmerged (merged with `--merge-overlaps`): records nA, nB; summed lengths LA, LB; n11 = accepted pairs; both inputs read to the end. G = Σ `-g` lengths (required) minus the merged length of `--exclude BED`, which only shrinks G, as bedtools.</li><li>Table as bedtools: n12 = max(0, nA−n11), n21 = max(0, nB−n11); `n22_full = max(n11+n12+n21, ⌊G·nA·nB / (2·nA·nB + LB·nA + LA·nB)⌋)` (= `⌊G/(2 + LA/nA + LB/nB)⌋` exactly; u128, f64 on overflow); an empty A or B gives `n11+n12+n21` (C's `(long long)NaN`); n22 = max(0, n22_full − n12 − n21 − n11).</li><li>p-values: klib's `kt_fisher_exact`, operation for operation, `lbinom` on `libm::lgamma`; ratio = `(n11/n12)/(n21/n22)`.</li><li>Output is bedtools' text byte for byte: the `#` lines (`%-12lld` cells), `left\tright\ttwo-tail\tratio`, `%.5g` p-values and `%.3f` ratio or the literal `-nan`.</li></ul> |
| D149 | **`stats reldist -a A -b B [--detail]`.** Midpoints `⌊(start+end)/2⌋` in u64 (bedtools casts the sum to `int`: `#reldist-midpoint-overflow`). Per A, bedtools' `lower_bound` neighbours among B midpoints on its chromosome; A before the first, after the last, or without B is skipped; rel = 0 when a distance is 0, else `min(dl, dr)/(right − left)` in f64. Summary: bins `⌊rel·100⌋`, rows `%.2f\t%lu\t%lu\t%.3f` under `reldist\tcount\ttotal\tfraction`, ascending. `--detail`: A's record, a tab, `%.3f`, in A's order (§4.6). Sorted streaming reads each chromosome's B once; `--unsorted` loads all B midpoints. |
| D150 | **`random -g G [-l 100] [-n 1000000]`.** Valid starts on chromosome c are `0..=len_c − l`. Interval *i* draws uniformly over all valid starts (prefix sums, binary search), then `coin()` → `+`: bedtools' rejection distribution without rejection. Rows `chrom start end i l strand` in *i* order, not sorted (help says so). `l` above every length is exit 2 naming the longest chromosome (bedtools loops forever: `#random-too-long`); `-g` without lengths is exit 2. |
| D151 | **`sample FILE -n N`.** Complying records (with `--strand`, that strand; unknown never complies) are numbered *j*; key *j* = stream 0 at word 2(*j*−1); the N smallest (key, *j*) are kept, a uniform N-subset, written in input order (bedtools: reservoir order, `#sample-order`). Headers per `--header`; VCF keeps its header. Fewer than N complying records: exit 1 with the count, as bedtools, and nothing written. |
| D152 | **`shuffle FILE -g G`.** <ul><li>Allowed set = `--include` regions (else the genome) minus `--exclude`, merged per chromosome. A record of length L may go to `(c, s)` with `a ≤ s ≤ b−L` inside one allowed `[a, b)`: the whole record stays inside (bedtools only places the start in an `-incl` interval: `#shuffle-placement`).</li><li>Uniform over all placements; `--same-chrom` over the record's chromosome; `--chrom-first` picks uniformly among chromosomes with a placement, then within it. Exact: allowed intervals sorted by length (ties: chromosome rank, start) with prefix sums make the count and the *r*-th placement binary searches.</li><li>No placement: exit 1 naming file, line and length, with the fix (shrink `--exclude`, widen `--include`, or `--skip-unplaceable`, which skips, counts and lists 10); nothing written. bedtools drops the record with a message.</li><li>`--no-overlapping`: input order; a placement overlapping an earlier one is redrawn from the same stream up to `--max-tries` (1000), then is unplaceable.</li><li>Other columns and strand kept (GFF via `rewrite_coords`); VCF/BCF is exit 2 suggesting `gloci convert`; the order warning applies.</li></ul> |
| D153 | **`shuffle --format bedpe`.** A one-chromosome pair is placed as its outer span, keeping its layout; a two-chromosome pair places each end independently; absent ends stay absent. With `--exclude` both ends must avoid excluded regions (the gap may cross them, as bedtools), redrawn up to `--max-tries`. `--include`, `--same-chrom`, `--chrom-first`, `--no-overlapping` are exit 2 (bedtools ignores them). |
| D154 | **`pair-intersect -a PAIRS -b FILE [--type T]`.** <ul><li>B (BED, GFF, VCF footprints) goes into an `IntervalIndex`; A streams in any order in ≤ 2 × threads batches, queried in parallel, written in A's order. Each present end is tested under §4.3 with the end as A.</li><li>`--type`: `either` (default) A+B per hit of end 1, then end 2; `both` the same when both ends hit; `xor` one end only; `neither` A alone, no hit; `notboth` A alone with no hit, A+B with one; `ispan`/`ospan` A+B per B hitting the inner `[end1, start2)` (`[end2, start1)` when end 1 is right of end 2) or outer span; `notispan`/`notospan` A alone when nothing hits it.</li><li>Span types need both ends on one chromosome (others skipped, counted in one warning); an empty inner span hits nothing; `-s`/`-S` with them is exit 2, as bedtools.</li><li>Hits in B's input order (bedtools: bin order, `#pair-hit-order`). A+B rows are TSV (D42); `notboth`'s A-alone rows are padded with `.` (D45, `#pair-notboth-padding`); A-alone types write BEDPE.</li></ul> |
| D155 | **`pair-pair -a PAIRS -b PAIRS [--type both]`.** B's end 1 and end 2 each get an `IntervalIndex`; orientation 1 tests A1–B1 and A2–B2, orientation 2 A1–B2 and A2–B1. `both`: B with both ends hit in either orientation; `either`: any end; `neither`: A alone with no B end hit; `notboth`: A alone with no B fully hit. Each (A, B) is written once with B's line as written (bedtools prints one row per orientation, `#pair-pair-once`, and in `either` swaps B's ends when only its end 2 hit, `#pair-pair-b-columns`), in B's input order. Strand is ignored unless `-s`/`-S` (bedtools enforces strand unless `-is`, matching `.` with `.`: `#pair-unknown-strand`). `-w`/`--left`/`--right`/`--strand-window` pad B per §4.3 (D37 applies). `--require-different-names` drops B named like A (column 7); a line without a name is then a data error. |
| D156 | **`convert --bedpe`** (alignments). Reads `AlnSource` with the D102 filters in file order; `@HD SO:coordinate` is exit 2 suggesting `gloci sort -n x.bam -o x.byname.bam`. Consecutive records sharing a QNAME form a group whose primary READ1 and READ2 (`0x40`/`0x80`, not `0x900`) make the pair; secondary and supplementary records are ignored (bedtools pairs whatever is adjacent: `#bamtobed-bedpe`); a group with one primary mate warns once with a count. Ends: `RNAME start end strand` (D101) or `. -1 -1 .` when unmapped; swapped so `(chrom1, start1) ≤ (chrom2, start2)` comparing names as bytes (`.` first), unless `--mate1-first`. Name = QNAME; score = min MAPQ when both mapped, else 0, or with `--score-tag TAG` the sum over mapped ends (`NM` = `-ed`; missing on a mapped end is a data error). Conflicts with `--bed12`, `--blocks`, `--cigar`, `--color`. |
| D157 | **Parallel selection** (P2a follow-up 1 lands here: `vcf popgen`/`compare` are new consumers of the stage, and P2b's evidence suite needs §9.1's `vcf stats` target). With `--threads` > 1, `run_selection` hands batches to `threads − 1` workers, each with its own copy of the stage (expression, scratch, BCF dictionaries, masks), returning per record pass/fail, `SiteStats`, mask flags and, under shared decoding, the samples in a `PackedGenotypes` arena. The driver keeps the order guard and `--thin`, takes batches back in order (≤ 2 × workers in flight) and calls the sink in record order. Batches hold ≤ max(1, 524,288 / samples) records. The earliest failing record is reported, so errors match at every `-t`. `-t 1` keeps today's path. Every `vcf` command uses it; the other D77 pools stay deferred. |
| D158 | **`vcf popgen IN --out-prefix P`.** <ul><li>Selection first (D50), then every requested statistic in one pass through the `vcf stats` engine; files publish atomically (D72).</li><li>Files (vcftools' names and columns): `--site-pi` `.sites.pi`; `--window-pi SIZE [--window-pi-step S]` `.windowed.pi`; `--tajima-d BIN` `.Tajima.D`; `--fst-pop FILE`×≥2 `.weir.fst` + `.weir.fst.summary`, with `--fst-window-size N [--fst-window-step S]` `.windowed.weir.fst`; `--geno-r2` `.geno.ld`; `--hap-r2` `.hap.ld`; `--interchrom-geno-r2`/`--interchrom-hap-r2` `.interchrom.{geno,hap}.ld`; `--relatedness[2]` `.relatedness[2]`.</li><li>Windowed statistics and window LD need input sorted by chromosome and position (the `--snp-density` rule: `#popgen-sorted-input`). A step larger than its window is exit 2 (vcftools resets it: `#popgen-window-step`).</li><li>Options never leak between statistics: populations select no samples, and hap-r2's phase rule removes no sites from others (vcftools does both: `#popgen-independent-options`). vcftools `--het` is `vcf stats --het`.</li></ul> |
| D159 | **π and Tajima's D**, as vcftools 0.1.16. <ul><li>A site counts only when every unmasked selected sample is diploid (`SiteStats` neither `haploid` nor `polyploid`); allele counts include `./1`'s called allele.</li><li>Site π = `Σ c_a(n−c_a) / (n(n−1))` over called alleles; NaN when n < 2.</li><li>Window π: windows `[k·step+1, k·step+size]`, k from `⌈(pos−size)/step⌉` (≥ 0) to `⌈pos/step⌉−1`. Polymorphic sites add `n(n−1)` comparisons and their mismatches; the rest of the window is monomorphic at `N(N−1)`, N = 2 × selected samples. N_VARIANTS = sites with `c_ref < n`; windows without one are not printed.</li><li>Tajima's D: biallelic diploid sites; bin = `⌊POS × (1/BIN)⌋` in f64 (POS 49, BIN 49 → bin 0); n = 2 × selected samples throughout; polymorphic sites add `p(1−p)`. Rows run from the first bin with S > 0 to the chromosome's last bin, empty ones `0 nan`; BIN_START = `bin × BIN`. Constants a1…e2 in f64 (vcftools' 32-bit `n²` wraps above 32,768 samples: `#popgen-integer-overflow`).</li></ul> |
| D160 | **Weir and Cockerham Fst** (vcftools' multi-allele form). Population files: one name per line (first word, blank lines skipped), after selection. Fewer than 2 files, a file with no selected sample, or a sample in two files is exit 2 naming it (vcftools skips, divides by zero or double-counts: `#fst-populations`); unknown names warn once with a count. Per diploid site: n_i, p_ij, h_j from per-population genotype counts, then a, b, c per allele in vcftools' operation order; fst = Σa/Σ(a+b+c) over alleles where none is NaN; every diploid site gets a row (an all-missing population gives NaN). `.weir.fst.summary` = `N_SITES\tMEAN_FST\tWEIGHTED_FST` over non-NaN sites (vcftools only logs them). Windows (D159's) sum a, a+b+c, fst and count non-NaN sites; printed when Σ(a+b+c) ≠ 0 and the count > 0; written as positions pass them. |
| D161 | **LD** (biallelic sites, 4 bits per sample). `--geno-r2`: dosage 2/1/0 for `0/0`/het/`1/1`; missing, haploid, polyploid or masked genotypes skipped; r² = cov²/(var₁var₂) over samples called at both, N_INDV. `--hap-r2`: haplotypes of phased diploid genotypes, D, D′, r², N_CHR as `calc_hap_r2`, using only sites where every selected genotype is phased (vcftools' implied `--phased`). Pairs (i, j>i) on one chromosome with `j−i ≤ --ld-window`, `j−i ≥ max(1, --ld-window-min)`, `--ld-window-bp-min ≤ Δpos ≤ --ld-window-bp`; `--min-r2 X > 0` drops r² < X and NaN; order (i, j); a site's pairs are written when its window closes. `--interchrom-*`: every pair on different chromosomes, (i, j) order, at the end. Held sites past `--max-mem` (2G): exit 1 naming the fix (`--ld-window-bp`, `--ld-window`, `--max-mem`); vcftools spills to a temporary file (`#ld-rules`). |
| D162 | **Relatedness** (biallelic diploid sites, selected samples). `--relatedness2` (Manichaikul): per pair N_AaAa, N_AAaa, N1_Aa, N2_Aa and φ = (N_AaAa − 2·N_AAaa)/(N1_Aa + N2_Aa), every ordered pair with self, header order. `--relatedness` (Yang Ajk, unadjusted): sites with ε < f < 1−ε; per pair i < j, Σ(x_i−2f)(x_j−2f)/(2f(1−f)) over samples called at both, divided by that count; self = 1 + Σ(x²−(1+2f)x+2f²)/(2f(1−f))/count. Masked and half-missing genotypes are missing (vcftools counts masked ones in `relatedness2` and `./1` as 0 in `relatedness`: `#relatedness-genotypes`). 64-site blocks become per-sample bit masks (pair counts add popcounts); matrix rows split over scoped threads, each pair summing sites in input order, so every `-t` gives the same output. |
| D163 | **`vcf compare A B --out-prefix P`.** <ul><li>Both inputs selected (D50; `--samples` applies per file by name, a name in neither is exit 2), walked in the shared table's order (§4.5: `-g`, else A's header, else karyotype/lex), each order-checked; a chromosome in one file only is fine (vcftools stops).</li><li>`--match alleles` (default, §7.8): same CHROM and POS and a shared ALT after §4.2 trimming (`intersect --match alleles`); each A record takes the first unmatched matching B record; leftovers at the position pair as `O`, the rest are `1`/`2`. No normalisation: the same change at another POS is another site. `--match position` reproduces vcftools' walk (same POS → `B`, differing REF → `O`, `N`/`.` REF taking the other's; overlapping REF spans → `O`).</li><li>Files (all unless some are requested): `--site-overlap` `.diff.sites_in_files`, `--sample-overlap` `.diff.indv_in_files`, `--site-discordance` `.diff.sites`, `--sample-discordance` `.diff.indv`, `--discordance-matrix` `.diff.discordance_matrix`.</li><li>Samples: the union by name after `--sample-map FILE` renames B's (`--diff-indv-map`), byte order. Genotypes compare as unordered allele strings; any missing allele makes it missing (vcftools counts `./1` as called: `#compare-half-missing`). The matrix counts biallelic pairs with identical REF/ALT, skipping haploid genotypes, as vcftools.</li></ul> |
| D164 | **Fixtures and goldens.** Interval fixtures in `xtask/diff/inputs/p3b/{stats,random,pairs}/`; alignment fixtures extend P3a's `aln/` through `make-alignments.sh`; VCF fixtures beside their goldens. Random commands cannot match bedtools' stream, so their goldens are forced cases with one possible placement (which must match) and `#rng-stream` cases pinning gloci's own output (`expected_gloci`). No new compare mode. |

## Review Focus

1. **`shuffle` with an `--exclude` (or `--include`) that leaves no room for a record:** exit 1 naming the record, its length and the fix, no output file; `--skip-unplaceable` writes the rest and counts the skipped. Task 6 `shuffle_exclude_leaves_no_room`, `shuffle_include_shorter_than_feature`.
2. **`fisher` with zero overlaps, every A overlapping, or an empty A:** bedtools' table, p-values and ratio text exactly (`0.000`, `inf`, `-nan`). Task 3 `fisher_table_heuristic`; goldens `fisher-no-overlap-vs-bedtools`, `fisher-all-overlap-vs-bedtools`, `fisher-empty-a-vs-bedtools`.
3. **Fst with a population whose genotypes are all missing at a site, or that names no sample:** vcftools' NaN rows, windows and summary skip them; no sample is exit 2. Task 12 `fst_population_all_missing`, `fst_population_rules`; golden `missingpop.weir.fst`.
4. **BEDPE with `.` chromosomes (`-1 -1` ends) and `.` strands:** an absent end overlaps nothing, `neither`/`notboth` still report the pair, span types skip it with a count, `.` never matches `-s`/`-S`. Task 7 `bedpe_absent_ends`, Task 8 `pair_pair_absent_ends`; goldens `pairtobed-absent-ends-vs-bedtools`, `pairtopair-unknown-strand-vs-bedtools` (pinned).
5. **`vcf compare` on files with different sample sets** (disjoint, partly shared, reordered, renamed): the union in byte order, one-file samples `0 0 nan`, discordance over shared samples only, `--sample-map` joining renamed ones. Task 15 `compare_sample_union`, Task 16 `compare_different_samples`; golden `diffsamples.diff.indv`.

## Lanes

Task 1 runs first, alone, and merges before the lanes start: it owns every manifest change, the format variants, the BEDPE reader, the random stream, `cfmt` and `SeedArg`. Then three lanes (4-core host):

| Lane | Tasks | Cross-lane waits |
|---|---|---|
| A: interval statistics, randomness | 2 → 3 → 4 → 5 → 6 | — |
| B: pairs, then compare | 7 → 8 → 9 → 15 → 16 | 15 after Task 10 (`Hit::genotypes`, `set_threads`) |
| C: selection stage, popgen | 10 → 11 → 12 → 13 → 14 | — |
| Join | 17 → 18 → 19 | after all lanes |

**Shared files that will conflict.** Never hand-merge a generated file.
- `SKILL.md`: each task regenerates it last (`cargo xtask docs --skill SKILL.md`); on a conflict take either side and rerun.
- `Cargo.toml`, crate manifests, `Cargo.lock`: Task 1 only; a later conflict means an unplanned dependency: stop and ask.
- `crates/gloci/src/cli.rs`: every task. Append structs at the end under `// P3b Task N` and variants at the end of `Commands`/`VcfCommand`/`StatsCommand`; keep both sides. `ConvertArgs` is Task 9's.
- `crates/gloci/src/{lib.rs, commands/mod.rs}`, `gloci-ops/src/lib.rs`, `gloci-vcf/src/lib.rs`: one line per module or command; keep both sides.
- `xtask/diff/cases.toml`: one appended `# P3b Task N` section per task. `docs/differences.md`: one section per anchor, by its task. `docs/migration.md`: Task 19 only.
- `crates/gloci-vcf/tests/golden/popgen/regen.sh` (Tasks 11–14) and `compare/regen.sh` (15–16): one appended `# Task N` block each; `make-alignments.sh`: Task 9's block.
- Single owners: `vcf_select.rs` and `gloci-vcf/src/gt.rs` (Task 10); `vcf_stats.rs` (Task 10, then Task 11's `run_stats`); `gloci-engine/src/index.rs` (7); `merge.rs` (2); `xtask/bench/*`, `gen.rs` (17); `behaviour.rs`, `large.rs` (18).

## File structure

```
crates/gloci-core/src/cfmt.rs (moved from gloci-vcf stats/fmt.rs)                    T1
crates/gloci-io/src/bedpe.rs; detect.rs, output.rs (Bedpe); gloci-ops/src/rng.rs     T1
crates/gloci-ops/src/stats/{mod,merged,jaccard,fisher,reldist}.rs                     T2–T4
crates/gloci-ops/src/{random,sample}.rs; shuffle/{mod,allowed,bedpe}.rs               T5–T6
crates/gloci-engine/src/index.rs; gloci-ops/src/pairs/{mod,intersect,pair_pair}.rs    T7–T8
crates/gloci-ops/src/convert/align_bedpe.rs                                           T9
crates/gloci-ops/src/vcf_select.rs; gloci-vcf/src/gt.rs (PackedGenotypes)             T10
crates/gloci-vcf/src/popgen/{mod,diversity,fst,ld,relatedness}.rs; gloci-ops/src/vcf_popgen.rs  T11–T14
crates/gloci-vcf/src/compare/{mod,walk,genotypes}.rs; gloci-ops/src/vcf_compare.rs    T15–T16
crates/gloci/src/commands/{stats,random,sample,shuffle,pairs}.rs; vcf.rs; cli.rs
```

Every task's CLI step adds the flags its D-row names, with 2–4 help examples, and its dispatch line; `stats` and `vcf popgen` help states each formula (§7.6, §7.8). Its last step regenerates SKILL.md. Task 18 checks every new command at `-t 1/2/8`; tasks test threads only where they add their own.

---

### Task 1: Foundation (D143, D145, D146)

**Files:**
- Modify: root `Cargo.toml` (`libm = "0.2"`); gloci-ops manifest (`rand_chacha`, `libm`); gloci-core manifest (`itoa`); `gloci-io/src/{detect.rs, output.rs, record.rs, lib.rs}`; every exhaustive `match` on `Format`/`OutFormat` (`grep -rn 'Format::Wig\|OutFormat::Tsv'`), each new arm a `UsageError`; `gloci-vcf/src/stats/*.rs` (now `gloci_core::cfmt`); `cli.rs` (`SeedArg`).
- Create: `gloci-core/src/cfmt.rs` (git-move of `gloci-vcf/src/stats/fmt.rs`), `gloci-io/src/bedpe.rs`, `gloci-ops/src/rng.rs`.
- Test: `gloci-io/tests/bedpe.rs`, a body in `gloci-io/src/fuzz_bodies.rs` run by `tests/fuzz_smoke.rs`, `gloci-ops/tests/rng.rs`, `gloci/tests/bedpe_gate.rs`, unit tests in `cfmt.rs`.

**Interfaces (produces):**
- `gloci_core::cfmt`: `write_g(out: &mut Vec<u8>, x: f64, prec: usize)`, `write_g6`, `write_e6`, `write_fixed(out: &mut Vec<u8>, x: f64, prec: usize)` (replacing `write_f1/f2/f5`), `write_u64`, `write_i64`.
- `Format::Bedpe` (`--format bedpe`, never detected); `OutFormat::Bedpe` (`.bedpe[.gz]`, `-O bedpe|bedpe.gz`), written like TSV.
- `gloci_io::bedpe`: `pub struct BedpeEnd { pub chrom: ChromId, pub span: Interval, pub strand: Strand }`; `pub struct BedpeMeta { pub ends: [Option<BedpeEnd>; 2], pub loc: u64, pub off: usize, pub len: u32 }`; `#[derive(Default)] pub struct BedpeBatch { pub buf: Vec<u8>, pub recs: Vec<BedpeMeta> }` (`raw(i)`, `len()`, `clear()`); `pub struct BedpeReader` (`open(path: &str, o: &ReaderOptions) -> Result<Self, IoError>`, `name()`, `header() -> &[Vec<u8>]`, `read_batch(&mut self, chroms: &mut ChromTable, batch: &mut BedpeBatch, max: usize) -> Result<bool, IoError>`, `diagnostics()`); `pub fn bedpe_name(raw: &[u8]) -> Option<&[u8]>` (column 7).
- `gloci_ops::rng`: `pub struct Stream` (`new(seed: u64, stream: u64)`, `at_word(seed: u64, stream: u64, word: u128)`, `next_u64()`, `below(n: u64) -> u64`, `unit() -> f64`, `coin() -> bool`); `pub fn seed_or_entropy(seed: Option<u64>) -> u64` (warns per D143 on `None`).
- `cli.rs`: `pub struct SeedArg { #[arg(long, value_name = "N")] pub seed: Option<u64> }`.

- [ ] **Step 1: Write the failing tests.**
  - `bedpe_reader_table`: `chr1\t10\t20\tchr2\t30\t40` → both ends, unknown strands; 10 columns with `+`/`-` → strands set; `.\t-1\t-1\tchr1\t5\t5\tr1\t0\t.\t+` → end 1 absent, end 2 the point 5 on `+`; 9 columns → error with `6, 7, 8 or at least 10 columns` and `line 1`; `start1 > end1` → error naming `start1`; `-5` on a named chromosome → `negative`; strand `x` → error naming column 9; leading `#`/`track` lines → `header()`; `--chrom-map 1→chr1` rewrites columns 1 and 4; `--skip-invalid` skips and counts.
  - `fuzz_bedpe` (a body in `gloci-io/src/fuzz_bodies.rs`, run by `tests/fuzz_smoke.rs`, §9.3): 2,000 seeded byte flips and truncations never panic.
  - `bedpe_gate`: a BEDPE-shaped file detects as BED; `gloci merge --format bedpe x.bedpe` → exit 2 naming `pair-intersect`, `pair-pair`, `shuffle --format bedpe`; a VCF opened as BEDPE → exit 2 naming VCF; `-o x.bedpe.gz --write-index` → exit 2 with `gloci index -s 1 -b 2 -e 3 --zero-based`.
  - `rng_known_answers`: the first four `Stream::new(42, 0).next_u64()` equal constants committed at implementation (comment: changing them is breaking); `new(42, 1)` differs; `at_word(42, 0, 2).next_u64()` = the second value of `new(42, 0)`.
  - `rng_below`: `below(1) == 0`; `below(u64::MAX)` returns; 300,000 draws of `below(3)` at seed 7 → χ² < 13.82 (2 dof, p = 0.001); `unit()` ∈ [0, 1); `coin()` frequency 0.5 ± 0.005 over 100,000.
  - `cfmt_table`: `write_g(1.234567e-4, 5)` → `0.00012346`; `(123456789.0, 5)` → `1.2346e+08`; `(1.0, 5)` → `1`; `(2.5e-5, 6)` → `2.5e-05`; NaN → `nan`; `write_fixed(0.125, 2)` → `0.12`; `(0.375, 2)` → `0.38`; `(∞, 3)` → `inf`.
- [ ] **Step 2: Run** `cargo test -p gloci-io --test bedpe && cargo test -p gloci-ops --test rng && cargo test -p gloci-core cfmt && cargo test -p gloci --test bedpe_gate`. Expected: FAIL.
- [ ] **Step 3: Implement.** The BEDPE reader uses the shared line machinery (`Class`), splits with `memchr`, interns both chromosomes. `below` is written as D143 states.
- [ ] **Step 4: Run** `cargo test --workspace` (all existing `vcf stats` tests pass on `cfmt`). Run `docker run --rm gloci-diff:p0 sh -c 'bedtools fisher 2>&1; bedtools shuffle 2>&1; vcftools --help'` and confirm `-exclude`, `-noOverlapping`, `--weir-fst-pop`.
- [ ] **Step 5: Commit** `feat: P3b foundation (BEDPE reader, seeded ChaCha8 stream, C number formatting)`

### Task 2: `stats jaccard` (D147)

**Files:** create `gloci-ops/src/stats/{mod.rs, merged.rs, jaccard.rs}`, `commands/stats.rs`, `xtask/diff/inputs/p3b/stats/`; modify `merge.rs` (factor its clustering; `gloci merge` output unchanged), testkit `naive.rs` (`naive_jaccard`), `cli.rs` (`Commands::Stats(StatsArgs)`, `StatsCommand::Jaccard(JaccardArgs)`: `OverlapArgs` minus `--split`, `RelationSortArg`, the relation input flags); test `gloci/tests/stats_jaccard.rs`.

**Interfaces:**
- Consumes: `run_relation`, `RelateInputs`, `OverlapCriteria::accepts`, `cfmt::write_g6`.
- Produces:
  - `pub struct MergedReader` (`RecordReader`): `new(inner: Box<dyn RecordReader>, by_strand: bool) -> Self`; raw `chrom\tstart\tend`, meta strand = the cluster's.
  - `pub struct Tally` (`RecordReader`, `consume_to_end() == true`): `new(inner: Box<dyn RecordReader>) -> (Tally, Arc<TallyCounts>)`; `pub struct TallyCounts { pub records: AtomicU64, pub bases: AtomicU64 }`.
  - `pub struct JaccardOptions { pub criteria: OverlapCriteria, pub sort: Option<SortSettings>, pub genome: Option<String>, pub chrom_order: ChromOrder }`; `pub struct JaccardResult { pub intersection: u64, pub union: u64, pub n_intersections: u64 }`; `pub fn jaccard(inputs: RelateInputs, o: &JaccardOptions) -> Result<JaccardResult, IoError>`; `pub fn write_jaccard(r: &JaccardResult, out: &mut Vec<u8>)`.

- [ ] **Step 1: Write the failing tests.**
  - `jaccard_table`: A `chr1 10 20`, `chr1 20 30`, `chr1 50 60`; B `chr1 15 25`, `chr1 55 70` → `15\t40\t0.375\t2` (unmerged A would give 3).
  - `jaccard_strand`: `-s` merges strands apart and drops opposite pairs; unknown strands never match, one warning. `jaccard_fraction`: `-f 0.5` drops a 5-base overlap of a 20-base merged A.
  - `jaccard_matches_naive`: proptest against `naive_jaccard` with points, bookends, nesting.
  - Empty inputs → `0\t0\tnan\t0`; VCF A uses footprints; unsorted input exits 1 with the `--sort` hint and `--sort` equals the sorted run.
- [ ] **Step 2: Add golden cases** (`bedtools jaccard -a input.bed -b b.bed …`): the table fixture (trap: bookended A must merge); `-s`; `-f 0.5`; `-f 0.5 -r`; a generated pair (`gen bed`, 2,000 records, 6 chromosomes); VCF A with `--vcf-span ref`; points inside B (match, or pin `#zero-length-overlap`); both empty (match, or `#nan-sign`).
- [ ] **Step 3: Implement**, run `cargo xtask diff --only stats-jaccard`, then `--update-golden`.
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci stats jaccard`

### Task 3: `stats fisher` (D148)

**Files:** create `gloci-ops/src/stats/fisher.rs`; modify `cli.rs` (`StatsCommand::Fisher(FisherArgs)`: Task 2's flags, `-g` required, `--merge-overlaps`, `--exclude FILE`); tests `gloci-ops/tests/fisher_exact.rs`, `gloci/tests/stats_fisher.rs`.

**Interfaces:**
- Consumes: `MergedReader`, `Tally`.
- Produces: `pub struct FisherOptions { pub criteria: OverlapCriteria, pub merge: bool, pub exclude: Option<String>, pub sort: Option<SortSettings> }`; `pub struct FisherCounts { pub n_a: u64, pub n_b: u64, pub overlaps: u64, pub len_a: u64, pub len_b: u64, pub genome_size: u64 }`; `pub struct FisherTable { pub n11: i64, pub n12: i64, pub n21: i64, pub n22: i64, pub n22_full: i64 }`; `pub fn fisher_counts(inputs: RelateInputs, genome: &Genome, o: &FisherOptions) -> Result<FisherCounts, IoError>`; `pub fn fisher_table(c: &FisherCounts) -> FisherTable`; `pub fn kt_fisher_exact(n11: i64, n12: i64, n21: i64, n22: i64) -> (f64, f64, f64)` (left, right, two-tail); `pub fn write_fisher(c: &FisherCounts, t: &FisherTable, out: &mut Vec<u8>)`.

- [ ] **Step 1: Write the failing tests.**
  - `kt_fisher_exact_matches_enumeration`: every table with n ≤ 40 (cells in steps of 3) equals a direct hypergeometric sum written in the test within 1e-9 relative; tables hitting klib's `q == 0` shortcut are checked for its fixed answers instead.
  - `fisher_table_heuristic`: nA 10, LA 1000, nB 20, LB 1000, n11 4, G 100,000 → (4, 6, 16, 631), `n22_full` 657; nA 3, LA 4, nB 3, LB 5, G 1000 → `n22_full` exactly 200 (f80 may give 199); nA 0 → `n22_full` = nB, n22 = 0.
  - `fisher_genome_needs_lengths`: no `-g`, or a `-g` without lengths → exit 2. `--exclude` subtracts the merged length of overlapping exclude records once.
- [ ] **Step 2: Add golden cases** (`bedtools fisher -a input.bed -b b.bed -g genome.txt …`): default; `-m`; `-f 0.5`; `-exclude gaps.bed`; `fisher-no-overlap-vs-bedtools`; `fisher-all-overlap-vs-bedtools` (ratio `inf`); `fisher-empty-a-vs-bedtools` (`1 1 1 -nan`); `fisher-estimate-integer-vs-bedtools` on the thirds counts (match, or pin `#fisher-estimate-exact`); a generated pair with p-values below 1e-5 (trap: `%.5g` exponent form).
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci stats fisher (bedtools table, klib exact test)`

### Task 4: `stats reldist` (D149)

**Files:** create `gloci-ops/src/stats/reldist.rs`; modify `cli.rs` (`StatsCommand::Reldist(ReldistArgs)`: `--detail`, `UnsortedArg`, `RelationSortArg`); test `gloci/tests/stats_reldist.rs`.

**Interfaces (produces):** `pub struct ReldistOptions { pub detail: bool, pub unsorted: bool, pub sort: Option<SortSettings> }`; `pub struct ReldistSummary { pub bins: BTreeMap<u64, u64>, pub total: u64 }` (key `⌊rel·100⌋`); `pub fn reldist(inputs: RelateInputs, o: &ReldistOptions, detail: Option<&mut dyn RecordWriter>) -> Result<ReldistSummary, IoError>`; `pub fn write_reldist(s: &ReldistSummary, out: &mut Vec<u8>)`.

- [ ] **Step 1: Write the failing tests.**
  - `reldist_table`: B midpoints 100, 200, 300; A midpoints 150, 129, 100, 300, 50, 400 → `0.00\t2\t4\t0.500`, `0.28\t1\t4\t0.250`, `0.50\t1\t4\t0.250` (trap: rel 0.29 lands in bin `0.28` because `0.29 × 100` is 28.999…; 50 and 400 are skipped).
  - `reldist_detail`: A's line, a tab and `0.290` for the 129 record, in A's order.
  - `reldist_unsorted_equals_sorted`; A on a chromosome without B is skipped; header only when nothing qualifies.
- [ ] **Step 2: Add golden cases** (`bedtools reldist -a input.bed -b b.bed [-detail]`): the table fixture; `-detail`; a generated pair; a fixture at 1,500,000,000 bp pinned `#reldist-midpoint-overflow`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci stats reldist`

### Task 5: `random` and `sample` (D144, D150, D151)

**Files:** create `gloci-ops/src/{random.rs, sample.rs}`, `commands/{random.rs, sample.rs}`, `xtask/diff/inputs/p3b/random/`; modify `cli.rs` (`Commands::Random(RandomArgs)`, `Commands::Sample(SampleArgs)` per D150–D151 with `SeedArg`); tests `gloci/tests/{random.rs, sample.rs}`.

**Interfaces:**
- Consumes: `Stream`, `seed_or_entropy`.
- Produces: `pub struct RandomOptions { pub length: u64, pub count: u64, pub seed: u64 }`; `pub fn random_intervals(genome: &Genome, o: &RandomOptions, threads: NonZeroUsize, w: &mut dyn Write) -> Result<(), IoError>`; `pub struct SampleOptions { pub count: u64, pub seed: u64, pub strand: Option<Strand> }`; `pub fn sample(reader: Box<dyn RecordReader>, chroms: &mut ChromTable, o: &SampleOptions, w: &mut dyn RecordWriter) -> Result<u64, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `random_pinned`: chromosomes of 100, 200, 700 bp, `-l 10 -n 5 --seed 42` → committed exact output.
  - `random_uniform`: 100,000 intervals of `-l 1` → chromosome counts χ² < 18.42 (2 dof, p = 1e-4) against 100:200:700; starts on chr3 in 10 bins χ² < 33.72; every interval inside its chromosome; ids 1..n.
  - `random_too_long`: `-l 800` → exit 2 naming `chr3 (700 bp)`. `random_huge_contig_count`: 200,000 contigs, `-n 10000`, < 2 s (debug).
  - `seed_warning_reproduces`: without `--seed` the warning names a seed, and rerunning with it gives the same bytes.
  - `sample_is_uniform`: 1,000 records, `-n 100`, seeds 1..=2000 → inclusion counts χ² < 1174 (999 dof, p = 1e-4) against 200.
  - `sample_input_order`; `sample_all_records` (`-n` = count returns the input); `sample_too_few` (exit 1 naming `37 records`, no output file); `sample_strand` (`--strand +` keeps `+` only, unknown never).
  - VCF keeps its header.
- [ ] **Step 2: Add golden cases:** `sample-all-vs-bedtools` (`-n` = count, `lines-unordered`; trap: extra columns and header); `sample-strand-vs-bedtools` (`-s forward` against `--strand +`, `lines-unordered`); `sample-too-few-vs-bedtools` (`exit-status`); `sample-subset-vs-bedtools` and `random-vs-bedtools`, pinned `#sample-order` and `#rng-stream`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci random and gloci sample (seeded, thread-independent)`

### Task 6: `shuffle`, including BEDPE (D152, D153)

**Files:** create `gloci-ops/src/shuffle/{mod.rs, allowed.rs, bedpe.rs}`, `commands/shuffle.rs`; modify `cli.rs` (`Commands::Shuffle(ShuffleArgs)` per D152 with `SeedArg`, `FormatArg`); tests `gloci-ops/tests/shuffle_allowed.rs`, `gloci/tests/shuffle.rs`.

**Interfaces:**
- Consumes: `Stream`, `BedpeReader`, `rewrite_coords`.
- Produces:
  - `pub struct Allowed`: `build(genome: &Genome, chroms: &ChromTable, include: Option<&[(ChromId, Interval)]>, exclude: &[(ChromId, Interval)]) -> Allowed`; `placements(&self, len: u64, chrom: Option<ChromId>) -> u128`; `place(&self, len: u64, chrom: Option<ChromId>, rng: &mut Stream) -> Option<(ChromId, u64)>`; `place_chrom_first(&self, len: u64, rng: &mut Stream) -> Option<(ChromId, u64)>`.
  - `pub struct ShuffleOptions { pub seed: u64, pub include: Option<String>, pub exclude: Option<String>, pub same_chrom: bool, pub chrom_first: bool, pub no_overlapping: bool, pub max_tries: u32, pub skip_unplaceable: bool }`; `pub struct ShuffleSummary { pub records: u64, pub unplaceable: u64 }`.
  - `pub fn shuffle(reader: Box<dyn RecordReader>, chroms: &mut ChromTable, genome: &Genome, o: &ShuffleOptions, threads: NonZeroUsize, w: &mut dyn RecordWriter) -> Result<ShuffleSummary, IoError>`; `pub fn shuffle_bedpe(reader: BedpeReader, chroms: &mut ChromTable, genome: &Genome, o: &ShuffleOptions, w: &mut dyn Write) -> Result<ShuffleSummary, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `allowed_counts`: allowed intervals of 10 and 30 bases → `placements(5)` 32 (6 + 26), `placements(0)` 42, `placements(31)` 0; a `chrom` restriction counts that chromosome only.
  - `shuffle_constraints`: proptest over genomes, include/exclude sets and features: every output keeps its length and other columns, lies in one include interval, overlaps no exclude interval (a point on an exclude edge is allowed), stays on its chromosome under `--same-chrom`, and under `--no-overlapping` overlaps no other output.
  - `shuffle_uniform`: 100,000 length-1 records on 100/200/700 bp → χ² < 18.42; the 10 + 30 allowed set with L = 5 → interval counts 6:26 (χ² < 10.83, 1 dof); `--chrom-first` → equal chromosome counts (χ² < 18.42).
  - `shuffle_exclude_leaves_no_room` (Review Focus 1): chr1 of 100 bp, exclude `chr1 0 95`, a 10-bp record → exit 1 naming `input.bed line 1`, `length 10`, `--exclude`, `--skip-unplaceable`; no output file. With `--skip-unplaceable`: exit 0, the other records written, warning `1 record could not be placed`.
  - `shuffle_include_shorter_than_feature` (include intervals of 5 bp); `shuffle_no_overlapping_max_tries` (records filling the genome → exit 1 after `--max-tries 5`); `shuffle_vcf_refused` (exit 2 with `gloci convert`).
  - `shuffle_bedpe`: a one-chromosome pair keeps its end lengths and gap; two-chromosome ends placed independently; a `.` end stays `. -1 -1`; with `--exclude` neither end overlaps; `--include`/`--same-chrom`/`--chrom-first`/`--no-overlapping` → exit 2.
  - `shuffle_huge_contig_count` (200,000 contigs); GFF rewritten 1-based.
- [ ] **Step 2: Add golden cases:** `shuffle-forced-vs-bedtools` (chromosomes `chrA` 10 bp and `chrB` 25 bp, BED6 records of those lengths, `-chrom` against `--same-chrom`; trap: start 0 is the only placement, so columns and passthrough must match); `shuffle-forced-bedpe-vs-bedtools` (a pair spanning a whole chromosome); `shuffle-include-containment` and `shuffle-vs-bedtools`, pinned `#shuffle-placement` and `#rng-stream`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci shuffle (exact uniform placement, include/exclude, BEDPE)`

### Task 7: `IntervalIndex` and `pair-intersect` (D154)

**Files:** create `gloci-engine/src/index.rs`, `gloci-ops/src/pairs/{mod.rs, intersect.rs}`, `commands/pairs.rs`, `xtask/diff/inputs/p3b/pairs/`; modify `gloci-engine/src/lib.rs`, testkit (`naive_pair_intersect`), `cli.rs` (`Commands::PairIntersect(PairIntersectArgs)`: `-a`, `-b`, `--type`, `OverlapArgs` minus `--split`); tests `gloci-engine/tests/interval_index.rs`, `gloci/tests/pair_intersect.rs`.

**Interfaces:**
- Consumes: `BedpeReader`, `OverlapCriteria::accepts`, D42 B rendering, `b_placeholder`.
- Produces:
  - `pub struct IndexEntry { pub chrom: ChromId, pub span: Interval, pub id: u32 }`; `pub struct IntervalIndex` with `build(entries: Vec<IndexEntry>) -> IntervalIndex` and `query(&self, chrom: ChromId, a: Interval, window: Option<Window>, a_strand: Strand, inclusive_points: bool, out: &mut Vec<u32>)` (ids ascending; the D61 rule; lappers in a `Vec` indexed by `ChromId`).
  - `pub enum PairType { Either, Neither, Both, Xor, NotBoth, InnerSpan, OuterSpan, NotInnerSpan, NotOuterSpan }` (values `either neither both xor notboth ispan ospan notispan notospan`); `pub struct PairIntersectOptions { pub kind: PairType, pub criteria: OverlapCriteria }`; `pub struct PairSummary { pub a_pairs: u64, pub span_skipped: u64 }`; `pub fn pair_intersect(a: BedpeReader, b: Box<dyn RecordReader>, chroms: &mut ChromTable, o: &PairIntersectOptions, threads: NonZeroUsize, out: &mut dyn Write) -> Result<PairSummary, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `interval_index_rules`: §4.3's point and bookend table with and without `--inclusive-points`; B at 3,000,000,000 and `u64::MAX − 10` found; a window pads B; ids ascending.
  - `pair_type_table`: A pairs with an end-1 hit only, end-2 only, both, neither → exact rows for all nine types.
  - `bedpe_absent_ends` (Review Focus 4): `.\t-1\t-1\tchr1\t100\t200` → end 1 never hits; `either` reports end-2 hits; `neither`/`notboth` report it when end 2 misses; span types skip it, warning count 1.
  - `pair_unknown_strand`: `-s` with `.` strands → no hits, one warning.
  - `pair_span_types`: end 1 left of end 2, right of it, touching (no inner hits; `notispan` reports A), two chromosomes (skipped, counted); `-s` with `ispan` → exit 2.
  - `pair_intersect_matches_naive` (proptest); `notboth_padding`; VCF B footprints; unsorted B accepted; `pair_huge_contig_count` (200,000 contigs).
- [ ] **Step 2: Add golden cases** (`bedtools pairtobed -a input.bedpe -b b.bed …`): all nine `-type` values on a fixture whose hits share one 16 kb bin; `-f 0.5`; `-s`, `-S` on `+`/`-` strands; `pairtobed-absent-ends-vs-bedtools`; 6-, 7-, 8- and 11-column A; hits in different bins (`lines-unordered` match plus an exact case pinned `#pair-hit-order`); `.` strands under `-s` pinned `#pair-unknown-strand`; `notboth` pinned `#pair-notboth-padding`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: BEDPE pair-intersect (bedtools pairtobed) on a reusable interval index`

### Task 8: `pair-pair` (D155)

**Files:** create `gloci-ops/src/pairs/pair_pair.rs`; modify `cli.rs` (`Commands::PairPair(PairPairArgs)`: `--type` default `both`, `OverlapArgs` minus `--split`, `WindowArgs`, `--require-different-names`); test `gloci/tests/pair_pair.rs`.

**Interfaces:**
- Consumes: `IntervalIndex`, `PairSummary`.
- Produces: `pub enum PairPairType { Neither, Either, Both, NotBoth }`; `pub struct PairPairOptions { pub kind: PairPairType, pub criteria: OverlapCriteria, pub window: Option<Window>, pub different_names: bool }`; `pub fn pair_pair(a: BedpeReader, b: BedpeReader, chroms: &mut ChromTable, o: &PairPairOptions, threads: NonZeroUsize, out: &mut dyn Write) -> Result<PairSummary, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `pair_pair_orientations`: a B matching only in orientation 2 is found by `both`; a B matching in both is written once.
  - `pair_pair_types`: exact rows for the four types on one fixture.
  - `pair_pair_window`: `-w 10` finds a B 10 bases away; `--left 10 --strand-window` extends downstream of each A end by its strand; `-w` with `-f` → exit 2 (D37).
  - `pair_pair_names`: `--require-different-names` drops a self hit; a 6-column line under it → exit 1 naming the line.
  - `pair_pair_absent_ends`; `-s` with `.` strands → no hits; proptest against a naive four-way comparison.
- [ ] **Step 2: Add golden cases** (`bedtools pairtopair -a input.bedpe -b b.bedpe …`): `-is -type both|either|neither|notboth`; default (strand enforced) against `-s` on `+`/`-` strands; `-is -slop 10` against `-w 10`; `-is -ss -slop 10` against `--left 10 --strand-window`; `-is -rdn` against `--require-different-names`; `pairtopair-unknown-strand-vs-bedtools` (`.` strands, default against `-s`) pinned `#pair-unknown-strand`; a both-orientation B pinned `#pair-pair-once`; an `either` hit on B's end 2 only pinned `#pair-pair-b-columns`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci pair-pair (bedtools pairtopair)`

### Task 9: `convert --bedpe` (D156)

**Files:** create `gloci-ops/src/convert/align_bedpe.rs`, fixtures `xtask/diff/inputs/aln/{pairs.sam, supp.sam}` (+ BAMs via an appended `make-alignments.sh` block); modify `convert/mod.rs`, `commands/convert.rs`, `cli.rs` (`ConvertArgs`: `--bedpe`, `--mate1-first` requiring it, the D156 conflicts); test `gloci/tests/convert_bedpe.rs`.

**Interfaces:**
- Consumes: P3a `open_aln_source`, `AlnSource`, `aln_ref` (`name`, `flag`, `mapq`, `ref_id`, `pos`, `span`, `strand`, `int_tag`), `AlnHeader.sort_order`.
- Produces: `pub struct BedpeOptions { pub mate1_first: bool, pub score_tag: Option<[u8; 2]> }`; `pub struct BedpeSummary { pub pairs: u64, pub missing_mate: u64, pub unpaired: u64 }`; `pub fn alignments_to_bedpe(src: Box<dyn AlnSource>, o: &BedpeOptions, out: &mut dyn Write) -> Result<BedpeSummary, IoError>`.

- [ ] **Step 1: Write the failing tests** (`bedpe_pairs_table` on `pairs.sam`): a proper pair → `chr1 100 150 chr1 300 350 r1 60 + -`; mate 2 unmapped → `. -1 -1 chr1 100 150 r2 0 . +`; both unmapped → `. -1 -1 . -1 -1 r3 0 . .`; a `chr10`/`chr2` pair → the `chr10` end first; `--mate1-first` keeps READ1 first; `--score-tag NM` sums NM, and a mapped end without NM → exit 1 naming the read; `supp.sam` (R1, R1 supplementary, R2) → one pair from the primaries; an orphan mate warns `1 read name has one mate`; `SO:coordinate` → exit 2 with `gloci sort -n`; SAM, BAM and CRAM agree.
- [ ] **Step 2: Add golden cases** (`bedtools bamtobed -bedpe -i input.bam`): `pairs.bam`; `-mate1`; `-ed` against `--score-tag NM`; SAM via `samtools view -b`; `supp.bam` pinned `#bamtobed-bedpe`.
- [ ] **Step 3: Implement, then update the goldens.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: convert --bedpe (bamtobed -bedpe, -mate1)`

### Task 10: Parallel selection stage (D157)

**Files:** modify `gloci-ops/src/vcf_select.rs`, `gloci-vcf/src/{gt.rs, lib.rs}`, `gloci-ops/src/{vcf_stats.rs, vcf_filter.rs, vcf_extract.rs}` (threads in; `Hit::genotypes` replaces `hit.decoded.reader()` + `MaskedGenotypes`), `commands/vcf.rs`; test `gloci/tests/vcf_parallel.rs`, unit tests in `gt.rs`.

**Interfaces (produces):**
- `impl SelectStage { pub fn set_threads(&mut self, n: NonZeroUsize) }`.
- `impl Hit<'_> { pub fn genotypes(&self) -> Option<&dyn GenotypeSource> }`: decoded samples, masked ones reading as missing; `None` unless `share_decoded` was called. The `decoded` field goes.
- `pub struct PackedGenotypes` (`clear()`, `push_record(&mut self, src: &dyn GenotypeSource, fields: GtFields, masked: &[bool]) -> Result<RecordSlot, String>`, `record(&self, slot: RecordSlot) -> PackedRecord<'_>`), with `PackedRecord: GenotypeSource`.
- `run_selection`'s signature is unchanged.

- [ ] **Step 1: Write the failing tests.**
  - `packed_roundtrip`: mixed ploidy, phased and unphased, alleles ≥ 63, DP/GQ/FT and masked flags read back equal.
  - `parallel_selection_equals_sequential`: `vcf filter`, `extract` and `stats` (all 15 statistics) byte-identical at `-t 1/2/8` on VCF and BCF with `--samples`, `--min-dp`, an expression, `--thin 50` and `--region`.
  - `parallel_selection_first_error`: a malformed record in the third batch → the same message at `-t 1` and `-t 8`, no output file.
  - `parallel_selection_bounded`: 20,000 samples at `-t 8` never hold more than 2 × workers × 524,288 decoded genotypes in flight (a test-only counter).
  - `one_pass_equals_separate_runs` passes at `-t 8`.
- [ ] **Step 2: Run** `cargo test -p gloci --test vcf_parallel`. Expected: FAIL.
- [ ] **Step 3: Implement** with scoped workers and bounded `crossbeam-channel`s.
- [ ] **Step 4: Run** `cargo test --workspace`. Time the synthetic 2,500-sample × 50k-record 5-statistic case at `-t 1/4/8` (release) and paste the numbers. If `-t 8` misses 0.1× of five vcftools runs, accumulate the integer per-sample counters of `--depth`, `--missing-indv` and `--het` on the workers per batch, summed in batch order (floating sums stay on the driver), and re-time.
- [ ] **Step 5: Commit** `perf(vcf): parallel selection stage for every vcf command`

### Task 11: `vcf popgen`: command, π, Tajima's D (D158, D159)

**Files:** create `gloci-vcf/src/popgen/{mod.rs, diversity.rs}`, `gloci-ops/src/vcf_popgen.rs`, `crates/gloci-vcf/tests/golden/popgen/{regen.sh, popgen.vcf}` + vcftools' outputs; modify `vcf_stats.rs` (extract `run_stats`), `cli.rs` (`VcfCommand::Popgen(Box<VcfPopgenArgs>)`: FILE, `VcfSelectArgs`, `--out-prefix`, this task's D158 flags, `--max-mem`, `--vcf-span` default `ref`, `--chrom-map`, threads, invalid), `commands/vcf.rs`; tests `gloci-vcf/tests/popgen_diversity.rs`, `gloci/tests/vcf_popgen.rs`.

**Interfaces:**
- Consumes: `Stat`, `StatsEngine`, `SiteStats`, `Hit::genotypes`, `set_threads`.
- Produces:
  - `pub(crate) fn run_stats(input: SelectedInput, select: &VcfSelect, prefix: &Path, fields: GtFields, need_genotypes: bool, threads: NonZeroUsize, command: &str, build: &dyn Fn(&[(usize, String)], &VcfHeader) -> Result<Vec<Box<dyn Stat>>, IoError>) -> Result<StatsSummary, IoError>`; `vcf_stats` calls it.
  - `pub struct SitePi`; `pub struct WindowPi` (`new(size: u64, step: u64)`); `pub struct TajimaD` (`new(bin: u64, n_samples: u64)`); each `impl Stat`.
  - `#[derive(Default)] pub struct PopgenRequest { pub site_pi: bool, pub window_pi: Option<(u64, u64)>, pub tajima_d: Option<u64> }` (Task 12 adds `fst: Option<FstRequest>`, Task 13 `ld: LdRequest`, Task 14 `relatedness`/`relatedness2: bool`).
  - `pub struct PopgenOptions { pub select: VcfSelect, pub request: PopgenRequest, pub out_prefix: PathBuf, pub max_mem: u64 }`; `pub fn vcf_popgen(input: SelectedInput, o: &PopgenOptions, threads: NonZeroUsize) -> Result<StatsSummary, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `site_pi_table`: counts (3, 1) → `0.5`; monomorphic → `0`; one called allele → `nan`; a site with a haploid sample skipped; a masked sample excluded.
  - `window_pi_table`: 2 samples, size 10, step 10, a site at POS 5 with counts (3, 1) → `chr1\t1\t10\t1\t0.05`; with step 5, a site at POS 7 is in windows 1–10 and 6–15 only.
  - `tajima_bin_reciprocal`: BIN 49, a polymorphic site at POS 49 → row `chr1\t0\t1\t…` (bin 0); POS 98 → bin 1. `tajima_trailing_bins`: bins after the first S > 0 bin print `0\tnan`, none before it.
  - `tajima_oracle`: D from Tajima (1989) with constant n, written in the test, within 1e-12.
  - `popgen_unsorted_window_stats`: unsorted input → exit 1 with `gloci sort`, no file; `--site-pi` alone accepts it. `--window-pi 10 --window-pi-step 20` → exit 2; no statistic → exit 2 listing the flags.
  - `popgen_one_pass_equals_separate_runs`.
- [ ] **Step 2: Add goldens:** `regen.sh` runs vcftools on `popgen.vcf` (diploid with missing genotypes, a haploid site, a multiallelic site, sites at POS 49 and 98, DP on every genotype) for `--site-pi`, `--window-pi 10 [--window-pi-step 5]`, `--TajimaD 49|100`, plus `--minDP 5` and `--indv` variants; `vcf_popgen.rs` compares each file; live cases mirror them.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci vcf popgen (site and window pi, Tajima's D)`

### Task 12: Weir and Cockerham Fst (D160)

**Files:** create `gloci-vcf/src/popgen/fst.rs`, fixtures `fst.vcf`, `missingpop.vcf`, `pop1.txt`–`pop3.txt`; modify `vcf_popgen.rs`, `cli.rs` (`--fst-pop` repeatable, `--fst-window-size`, `--fst-window-step`), `regen.sh`; tests `gloci-vcf/tests/popgen_fst.rs`, `vcf_popgen.rs`.

**Interfaces (produces):** `pub struct FstRequest { pub pop_files: Vec<PathBuf>, pub window: Option<(u64, u64)> }`; `pub struct Populations { pub names: Vec<String>, pub members: Vec<Vec<usize>> }` with `load(files: &[PathBuf], selected: &[(usize, String)]) -> Result<Populations, IoError>`; `pub struct WeirFst` (`new(pops: Populations, window: Option<(u64, u64)>)`, `impl Stat`).

- [ ] **Step 1: Write the failing tests.**
  - `fst_oracle`: per-allele Weir–Cockerham (1984) written in the test, biallelic and triallelic, 2 and 3 populations, within 1e-12.
  - `fst_population_all_missing` (Review Focus 3): population 2 all missing at site 2 → that row `nan`, its window skips it, the summary counts the others.
  - `fst_population_rules`: one file → exit 2; a file with no selected sample → exit 2 naming it; a sample in two files → exit 2 naming both; unknown names → one warning with the count.
  - `fst_windows_progressive`: equal to vcftools' end-of-file output; unsorted input with a window → exit 1. `--min-dp` changes n per population.
- [ ] **Step 2: Add goldens:** `--weir-fst-pop` with 2 and 3 populations, `--fst-window-size 100 --fst-window-step 50`, `--minDP 5`, and `missingpop.weir.fst`. The summary file is checked against vcftools' two logged estimates printed with `write_g6`.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: vcf popgen Weir and Cockerham Fst (per site, windowed, summary)`

### Task 13: Linkage disequilibrium (D161)

**Files:** create `gloci-vcf/src/popgen/ld.rs`, fixtures `ld.vcf` (unphased) and `phased.vcf` (one site with an unphased genotype); modify `vcf_popgen.rs`, `cli.rs` (the D161 flags), `regen.sh`; tests `gloci-vcf/tests/popgen_ld.rs`, `vcf_popgen.rs`.

**Interfaces (produces):** `pub struct LdWindow { pub snps: u64, pub bp: u64, pub min_snps: u64, pub min_bp: u64, pub min_r2: Option<f64> }` (defaults `u64::MAX`, `u64::MAX`, 1, 0, `None`); `#[derive(Default)] pub struct LdRequest { pub geno: bool, pub hap: bool, pub interchrom_geno: bool, pub interchrom_hap: bool, pub window: LdWindow }`; `pub struct GenoR2`, `pub struct HapR2` (each `new(window: LdWindow, budget: u64)`), `pub struct InterchromR2` (`new(haplotypes: bool, budget: u64)`), each `impl Stat`.

- [ ] **Step 1: Write the failing tests.**
  - `geno_r2_oracle` (Pearson r² of dosages over samples called at both) and `hap_r2_oracle` (Lewontin D, D′, r² from haplotype counts), both written in the test.
  - `ld_window_table`: `--ld-window 2`, `--ld-window-bp 100`, `--ld-window-min 2`, `--ld-window-bp-min 10` → exact pair lists; no pair crosses a chromosome; `--min-r2 0.2` drops low and NaN pairs.
  - `ld_genotype_rules`: haploid, polyploid, masked and half-missing genotypes skipped; `phased.vcf`'s unphased site is absent from `hap.ld` and present in `geno.ld` in the same run.
  - `ld_memory_guard`: `--max-mem 1K` with an unbounded window → exit 1 naming `--ld-window-bp`. `interchrom_order`: (i, j) order, different chromosomes only. Unsorted input with a window → exit 1.
- [ ] **Step 2: Add goldens:** `--geno-r2`; `--hap-r2` on `phased.vcf`; `--ld-window 2`; `--ld-window-bp 100`; `--ld-window-min 2`; `--min-r2 0.2`; `--interchrom-geno-r2`; `--interchrom-hap-r2`.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: vcf popgen LD (geno/hap r2, windows, interchromosomal)`

### Task 14: Relatedness (D162)

**Files:** create `gloci-vcf/src/popgen/relatedness.rs`, fixture `related.vcf` (a parent–child pair, unrelated samples, a half-missing genotype, DP for masking); modify `vcf_popgen.rs`, `cli.rs` (`--relatedness`, `--relatedness2`), `regen.sh`; tests `gloci-vcf/tests/popgen_relatedness.rs`, `vcf_popgen.rs`.

**Interfaces (produces):** `pub struct RelatednessYang` and `pub struct RelatednessManichaikul`, each `new(samples: Vec<(usize, String)>, threads: NonZeroUsize)` and `impl Stat`.

- [ ] **Step 1: Write the failing tests.**
  - `relatedness2_oracle` (φ from per-pair counts) and `relatedness_oracle` (Yang's equation 6), written in the test, within 1e-12.
  - `relatedness_genotype_rules`: a masked genotype is missing in both; a half-missing genotype is missing in `--relatedness`.
  - `relatedness_blocks`: 130 sites (two full blocks and a partial one) equal the site-by-site oracle.
  - Rows: every ordered pair in `relatedness2`, i ≤ j in `relatedness`, header order; `--samples` restricts them. `-t 1/2/8` identical; a 5,000-sample run stays under the D162 bound (`#[ignore]`, RSS).
- [ ] **Step 2: Add goldens:** `--relatedness`, `--relatedness2`, both with `--indv`; masking against `--minDP 5` (match, or `#relatedness-genotypes`); the half-missing fixture pinned `#relatedness-genotypes`.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: vcf popgen relatedness (Yang Ajk, Manichaikul phi)`

### Task 15: `vcf compare`: walk, site and sample overlap (D163)

**Files:** create `gloci-vcf/src/compare/{mod.rs, walk.rs}`, `gloci-ops/src/vcf_compare.rs`, `crates/gloci-vcf/tests/golden/compare/{regen.sh, a.vcf, b.vcf, diffsamples-a.vcf, diffsamples-b.vcf, map.txt}` + vcftools' outputs; modify `cli.rs` (`VcfCommand::Compare(Box<VcfCompareArgs>)`: `A B`, `VcfSelectArgs`, `--out-prefix`, `--match`, `--sample-map`, `--site-overlap`, `--sample-overlap`, `-g`, chrom order, `--chrom-map`, threads), `commands/vcf.rs`; tests `gloci-vcf/tests/compare_walk.rs`, `gloci/tests/vcf_compare.rs`.

**Interfaces:**
- Consumes: `set_threads`, `Hit::genotypes`, `gloci_core::trimmed_alleles`.
- Produces: `pub enum MatchRule { Alleles, Position }`; `pub enum SiteClass { Both, First, Second, Overlap }` (`B 1 2 O`); `pub fn share_allele(a: &SiteView, b: &SiteView) -> bool`; `pub struct CompareRequest { pub site_overlap: bool, pub sample_overlap: bool, pub site_discordance: bool, pub sample_discordance: bool, pub discordance_matrix: bool }` (all false = all); `pub struct CompareOptions { pub select: VcfSelect, pub matching: MatchRule, pub sample_map: Option<PathBuf>, pub request: CompareRequest, pub out_prefix: PathBuf }`; `pub struct CompareSummary { pub both: u64, pub first: u64, pub second: u64, pub overlap: u64, pub files: Vec<PathBuf> }`; `pub fn vcf_compare(a: SelectedInput, b: SelectedInput, o: &CompareOptions, threads: NonZeroUsize) -> Result<CompareSummary, IoError>`.

- [ ] **Step 1: Write the failing tests.**
  - `compare_walk_table` (alleles): identical → `B`; same POS, different ALT → `O`; A `A>C,G` vs B `A>C` → `B`; A `A>C,G` vs B `A>C` and `A>G` → `B` and `2`; `AT>A` vs `ATG>AG` at one POS → `B`; a chromosome only in B → `2` rows.
  - `compare_walk_position`: vcftools' rules, including REF spans overlapping at different positions → `O` and a REF of `N`.
  - `compare_sample_union` (Review Focus 5): disjoint, partly shared and reordered sample lists → `.diff.indv_in_files` in byte order with `B`/`1`/`2`; `--sample-map` joins renamed samples; `--samples x` with x in neither file → exit 2.
  - `compare_chrom_tables`: header order `chr1, chr2, chr10` in both files works; an out-of-order record → exit 1 naming the file. `compare_huge_contig_count` (200,000 contigs).
- [ ] **Step 2: Add goldens:** vcftools `--vcf a.vcf --diff b.vcf --diff-site` and `--diff-indv` (with and without `--diff-indv-map map.txt`) against `--match position`; alleles mode on the same files pinned `#compare-allele-matching`.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: gloci vcf compare (site and sample overlap, allele or position matching)`

### Task 16: `vcf compare`: genotype discordance (D163)

**Files:** create `gloci-vcf/src/compare/genotypes.rs`; modify `vcf_compare.rs`, `cli.rs` (`--site-discordance`, `--sample-discordance`, `--discordance-matrix`), `regen.sh`; tests `gloci-vcf/tests/compare_genotypes.rs`, `vcf_compare.rs`.

**Interfaces (produces):** `pub enum GtCompare { BothMissing, FirstMissing, SecondMissing, Concordant, Discordant }`; `pub fn compare_genotypes(a: &SampleFields, a_site: &SiteView, b: &SampleFields, b_site: &SiteView) -> GtCompare`; `pub struct DiscordanceMatrix { pub cells: [[u64; 4]; 4] }` (file 1 × file 2: `0/0`, `0/1`, `1/1`, `./.`).

- [ ] **Step 1: Write the failing tests.**
  - `genotype_compare_table`: `0/1` vs `1/0` → concordant; ALT `C,G` vs `G,C` with `0/1` vs `0/2` → concordant; `0/1` vs `./.` → second missing; `./1` → missing; a masked sample → missing; haploid `1` vs `1` → concordant.
  - `discordance_matrix_rules`: only biallelic pairs with identical REF and ALT; haploid genotypes skipped.
  - `compare_different_samples` (Review Focus 5): `.diff.indv` lists the union with one-file samples `0\t0\tnan`, shared samples counting shared sites only; `.diff.sites` gives `1`/`2` rows `0\t0\t0\tnan`.
- [ ] **Step 2: Add goldens:** `--diff-site-discordance`, `--diff-indv-discordance`, `--diff-discordance-matrix` on `a.vcf`/`b.vcf` and on the different-sample pair (`diffsamples.diff.indv`), against `--match position`; the half-missing fixture pinned `#compare-half-missing`.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: vcf compare genotype discordance (per site, per sample, matrix)`

### Task 17: P3b benchmarks

**Files:** `xtask/bench/cases.toml`, `xtask/src/gen.rs` (`gen bedpe --pairs N --chroms C --inter-frac F --seed S`; `gen vcf --pops K --pop-prefix P`; `gen sam --paired --by-name`), `xtask/tests/gen.rs`.

| Case | Target |
|---|---|
| `jaccard-plain-10M`, `fisher-plain-10M`, `reldist-plain-10M`, `sample-plain-10M` (`-n 1000000`) | `gloci <= 0.5 * bedtools` |
| `shuffle-plain-10M` (with `--exclude`), `random-10M`, `pair-intersect-1M-vs-10M`, `pair-pair-1M` | `gloci <= bedtools` |
| `convert-bedpe-50M` (name-sorted BAM) | `gloci <= bedtools bamtobed -bedpe` |
| `vcf-stats-synthetic-2500` (existing, now `-t 8`) | `gloci <= 0.1 * vcftools-5runs` (D157) |
| `popgen-synthetic-2500` (`--site-pi --window-pi 10000 --tajima-d 10000 --fst-pop p1 --fst-pop p2`) | `gloci <= 0.2 * vcftools-4runs` |
| `relatedness2-synthetic-500`, `geno-r2-synthetic-500` (`--ld-window-bp 10000`) | `gloci <= vcftools` |
| `vcf-compare-synthetic-2500` (all outputs) | `gloci <= 0.2 * vcftools-3runs` |

- [ ] **Step 1: Write the tests.** Bench cases parse (P1's test). `gen_bedpe_is_valid`: the requested pairs, inter-chromosomal fraction ±2%, every line valid for `BedpeReader`. `gen_vcf_pops`: population files partition the samples.
- [ ] **Step 2: Implement.**
- [ ] **Step 3: Run** `cargo test -p xtask`, then the synthetic `vcf-stats`, `popgen` and `jaccard` cases locally (best effort). Paste the medians.
- [ ] **Step 4: Commit** `feat(xtask): P3b benchmarks (stats, randomness, pairs, popgen, compare)`

### Task 18: Behaviour matrices and large-file tests

**Files:** `crates/gloci/tests/behaviour.rs`, `crates/gloci/tests/large.rs`.

Extend P2a's matrices to every P3b command: determinism at `-t 1/2/8` (seeded commands with `--seed`); broken pipe; empty inputs; the failure matrix (unsorted input, malformed lines, truncated BGZF, missing output directory, no partial file); the `--chrom-map` round trip; and `#[ignore]` RSS ceilings on generated multi-GB inputs for `jaccard`, `fisher`, sorted `reldist`, `sample -n 1000`, `pair-intersect` (≤ 2 × B's size + 200 MB), windowed `vcf popgen` and `vcf compare`.

- [ ] **Step 1: Write the matrices.**
- [ ] **Step 2: Run** `cargo test -p gloci --test behaviour`, and the large tests with `GLOCI_LARGE_BYTES=200000000`. Fix real failures in their owning module, with a regression test.
- [ ] **Step 3: Commit** `test: P3b determinism, broken-pipe, empty-input and failure matrices; RSS ceilings`

### Task 19: Documentation

**Files:** `docs/migration.md`, `docs/differences.md` (links only), `SKILL.md`.

- [ ] **Step 1: Update `docs/migration.md`**, checked against `--help`, each row linking its anchors:
  - bedtools: `jaccard`, `fisher`, `reldist`, `random`, `sample`, `shuffle`, `pairtobed`, `pairtopair`, `bamtobed -bedpe`/`-mate1`/`-ed`, with every D142 rename and not-planned row;
  - vcftools: `--site-pi`, `--window-pi[-step]`, `--TajimaD`, `--weir-fst-pop`, `--fst-window-size/-step`, `--geno-r2`, `--hap-r2`, `--interchrom-*`, `--ld-window*`, `--min-r2`, `--relatedness[2]`, `--het` (`vcf stats --het`), `--diff`/`--gzdiff`/`--diff-bcf` (second file positional), every `--diff-*` output, `--diff-indv-map` → `--sample-map`, the not-planned rows; `vcf-compare` → `vcf compare`;
  - a "Randomness" note under Other behaviour changes: `--seed` reproduces gloci's output at any `-t`, never bedtools'.
- [ ] **Step 2: Regenerate SKILL.md** with pitfalls: BEDPE is never detected (`--format bedpe` for `shuffle`); popgen windows need sorted input; `vcf compare` matches by alleles unless `--match position`; `sample` is not `samples`.
- [ ] **Step 3: Run** `cargo test -p xtask` and the help-examples test.
- [ ] **Step 4: Commit** `docs: P3b migration rows, randomness note and SKILL.md`

---

## Spec gaps this plan rules on

1. BEDPE detection, absent ends, column counts and output (D146).
2. `--include`/`--exclude` as a BED in `shuffle`/`fisher` and an expression in `vcf` (D142; owner to confirm as a third exception).
3. `sample` beside P3c's `samples` (D142).
4. §7.3's "`--seed` required for reproducibility" read as optional with a logged seed; what randomness guarantees without bedtools parity (D143, D144).
5. `shuffle` placement (whole record inside, exact uniform, unplaceable records are errors) and BEDPE flags (D152, D153).
6. `sample` output order (D151).
7. `pair-pair`'s strand default, one row per pair, B's columns as written (D155); `pair-intersect`'s `notboth` padding and span edges (D154).
8. Where `bamtobed -bedpe` lives, and pairing only primary mates (D156).
9. P2a follow-up 1 (`vcf stats` parallelism) belongs in P3b (D157).
10. Popgen file names, sorted input for windows, options that never leak, the new `.weir.fst.summary` (D158, D160).
11. Fst population rules (D160); LD memory with unbounded windows (D161); relatedness with masked and half-missing genotypes (D162).
12. "Matched by alleles" in `vcf compare` made precise, plus `--match position` for vcftools parity (D163).
13. NaN text and platform-independent transcendental functions (D145).

## Execution notes for the controller

- **Order.** Task 1 merges first; then the three lanes run in separate worktrees with their own `target/`, `CARGO_BUILD_JOBS=1`, scoped test commands and one `cargo test --workspace` per task before merging. Task 15 starts after Task 10 merges. Tasks 17–19 run last, in order.
- **Docker.** No image change. Regenerate goldens in `gloci-diff:p0` with `cargo xtask diff --update-golden` and the two `regen.sh` scripts.
- **Follow-ups.** When Task 10 meets its target, mark P2a follow-up 1 closed in `docs/superpowers/plans/2026-11-10-gloci-p2a-followups.md`.
- **Raise with the owner:** D142's `--include`/`--exclude` exception and §7.3's `--seed` wording (D143).
