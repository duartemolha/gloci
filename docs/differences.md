# Differences from bedtools, bedops, vcftools and htslib

gloci aims to produce the same output as the tools it replaces. Where it
deliberately does not, the difference is listed here with the reason. Nothing
else may differ: the differential harness (`cargo xtask diff`) runs gloci and
pinned reference tools (bedtools 2.31.1, bedops 2.4.41, vcftools 0.1.16,
htslib 1.21, bcftools 1.21, GNU coreutils) on the same inputs, and
`crates/gloci/tests/golden/` holds the reference outputs, which every
`cargo test` checks without Docker.

Each section has an anchor and lists the cases that pin it. A case in
`xtask/diff/cases.toml` that names a section as its `expected_difference` must
produce different output (or a different exit status) from the reference tool,
and gloci's own output for it is committed (`expected_gloci`) and must stay
byte for byte the same, so a difference cannot quietly turn into a different
wrong answer. If the behaviour stops differing, the case fails as stale until
the section is removed. Cases without an `expected_difference` must match the
reference tool. Most cases run on every `cargo test` from committed reference
outputs; the few marked "live only" need Docker and run nightly (`cargo xtask
diff`). The commit history of `crates/gloci/tests/golden/` shows every change in
a reference tool's output.

Each section gives the gloci behaviour, the old behaviour and the reason.

## Default `--vcf-span allele` {#default-vcf-span-allele}

- **gloci:** a VCF record's footprint is the reference bases its ALT alleles
  change: the shared prefix and suffix of REF and ALT are trimmed, and a pure
  insertion is a point. `--vcf-span ref` gives the REF span.
- **Old behaviour:** bedtools, tabix and bcftools use `[POS-1, POS-1+len(REF))`
  (and INFO/END for symbolic alleles). That is `--vcf-span ref`.
- **Reason:** the REF span ignores the alleles, so for example
  `bedtools subtract` against a VCF removes the wrong number of bases for an
  indel or a padded SNV. Pass `--vcf-span ref` for the old result.
- **Pinned by:** `view-region-allele-span-vs-tabix` (a deletion whose anchor
  base is in the region: must differ) and `view-region-ref-inclusive-vs-tabix`
  (`--vcf-span ref`: must match).

## Default karyotype chromosome order {#karyotype-chromosome-order}

- **gloci:** without `-g` and without a VCF/BCF header that lists contigs,
  chromosomes are ordered `chr1, chr2 ... chr22, chrX, chrY, chrM`.
  `--chrom-order lex` selects byte order.
- **Old behaviour:** `sort-bed`, `bedtools sort`, and `LC_ALL=C sort -k1,1` order
  lexicographically: `chr1, chr10, chr11 ... chr2`.
- **Reason:** karyotype order is what people expect from a sorted genome file,
  and a file sorted either way is detected, so a mismatch is an error with a
  hint, never a silently wrong answer. See [migration.md](migration.md) for the
  two fixes (re-sort once, or pass `--chrom-order lex` everywhere).
- **Pinned by:** `sort-bed-karyotype-default-vs-bedtools-sort` (must differ).
  Every other sort case passes `--chrom-order lex` or `-g`, and matches.

## Sort tie order {#sort-tie-order}

- **gloci:** `sort` orders by chromosome, start and end. Records that tie on
  all three keep their input order (the sort is stable), unless `--full-key`
  breaks the tie on the whole line.
- **Old behaviour:**
  - `sort-bed` orders by chromosome, start, end and then the rest of the line,
    which is what `--full-key` does: `sort --full-key` matches `sort-bed`
    (`sort-bed-vs-sort-bed`), plain `sort` does not on records that differ only
    after the end column (`sort-bed-tie-order-vs-sort-bed`).
  - `bedtools sort` and `sort -k1,1 -k2,2n -s` order by chromosome and start
    only and keep input order for records with the same start, so a longer
    feature that came first stays first (`sort-bed-tie-order-vs-bedtools-sort`).
    The two agree on distinct starts (`sort-bed-vs-bedtools-sort`,
    `sort-bed-vs-gnu-sort`).
  - `bcftools sort` places VCF records at the same POS in an order that is not
    input order (`sort-vcf-tie-order-vs-bcftools-sort`). gloci keeps input order
    for equal POS and REF end. The files agree whenever POS is distinct
    (`sort-vcf-vs-bcftools-sort`).
- **Reason:** ordering by start and then end gives one defined order, and
  records that tie on all three keys keep their input order, so repeated runs
  give identical output. Use `--full-key` to reproduce `sort-bed` exactly.
- **Pinned by:** `sort-bed-tie-order-vs-sort-bed`,
  `sort-bed-tie-order-vs-bedtools-sort`, `sort-vcf-tie-order-vs-bcftools-sort`
  (must differ) and, on inputs without ties, `sort-bed-vs-sort-bed` (with
  `--full-key`), `sort-bed-vs-bedtools-sort`, `sort-bed-vs-gnu-sort`,
  `sort-bed-vs-bedtools-sort-genome`, `sort-vcf-vs-bcftools-sort` (must match).

## `view` region footprints and points {#view-region-footprints}

- **gloci:** `view` selects records whose footprint (per `--vcf-span`)
  overlaps the region (half-open intervals, so bookended features do not
  overlap). A point (a zero-length BED record, an insertion under `allele`)
  overlaps a region only when strictly inside it, so an insertion on a region's
  edge is excluded. `--inclusive-points` includes points on either edge.
- **Old behaviour:** for VCF, `tabix` returns every record whose REF span
  touches the 1-based region, so an insertion anchored on the first or last
  base of the region is returned. For BED, `tabix` treats a zero-length record
  as empty: it returns a point strictly inside the region and never one on the
  edge, which is gloci's default; `--inclusive-points` therefore adds edge
  points that tabix does not return.
- **Reason:** an insertion between two bases belongs to neither base, so it is
  in a region only when the region contains both neighbours.
  `--vcf-span ref --inclusive-points` reproduces tabix for VCF.
- **Pinned by:** `view-region-default-vs-tabix` (insertions on region edges:
  must differ), `view-region-allele-inclusive-vs-tabix` (same queries with
  `--inclusive-points`: must match), `view-bed-region-points-inclusive-vs-tabix`
  (must differ), `view-bed-region-points-vs-tabix` and `view-bed-region-vs-tabix`
  (must match).

## `##glociCommand` header line {#glocicommand-header}

- **gloci:** VCF and BCF output of `sort` and `view` carries the original header
  plus one `##glociCommand=<command line>; Version=<version>` line.
- **Old behaviour:** `bcftools sort` (htslib 1.21) adds no provenance line, but
  it does rewrite the header: it puts `##FILTER=<ID=PASS,Description="All
  filters passed">` second (adding it when the input lacks it), which gloci does
  not do. `tabix` prints no header for a region query; GNU sort and `sort-bed`
  never add header lines.
- **Reason:** provenance in the file, in one line a tool can recognise. The
  `##glociCommand` line is normalised in the committed gloci outputs because it
  holds the command line.
- **Pinned by:** `sort-vcf-header-vs-bcftools-sort` (the input already has the
  PASS line second, so the `##glociCommand` line is the only difference: must
  differ). The `vcf-records` cases ignore header lines.

## `sort --check` checks `(chrom, start)` only {#sort-check-chrom-start}

- **gloci:** `sort --check` accepts a file ordered by chromosome and start. Two
  records with the same start may be in any order of end (exit status 0).
- **Old behaviour:** `sort-bed --check-sort` also requires ends to be ordered
  when starts are equal ("not properly sorted by end coordinates when start
  coordinates are identical"), and `sort -c -k1,1 -k2,2n` rejects such a file
  through its last-resort whole-line comparison.
- **Reason:** `(chrom, start)` is the order every streaming command needs, and
  files sorted by `bedtools sort` or `sort -k1,1 -k2,2n`,
  which leave ends unordered, must be accepted. A file that `gloci sort` wrote
  passes as well.
- **Pinned by:** `sort-check-start-ties-vs-sort-bed` (exit status 0 against
  `sort-bed --check-sort` exit status 1).

## A missing BGZF end-of-file marker is an error {#bgzf-eof-marker}

- **gloci:** a seekable BGZF file without the 28-byte end-of-file block is
  rejected: `BGZF end-of-file marker missing; the file is probably truncated`
  (exit status 1). Standard input cannot be checked.
- **Old behaviour:** htslib (`tabix`, `bcftools`) prints the warning `EOF
  marker is absent. The input may be truncated` and carries on with exit
  status 0.
- **Pinned by:** `truncated-bgzf-vs-bgzip` (a file cut at a block boundary:
  gloci exit status 1, `bgzip -dc` exit status 0).
- **Reason:** a missing marker is the signature of a truncated download or an
  interrupted write, and gloci never gives a quiet partial answer.

## `view --assume-sorted` {#view-assume-sorted}

- **gloci:** for an uncompressed BED file with no index, `view` scans the
  whole file and checks the order of every record it reads, so an unsorted file is an
  error. `--assume-sorted` binary-searches instead (like `bedextract`) and trusts
  the file to be sorted in the order given by `-g` or `--chrom-order`; it checks
  only the block it reads, and on an unsorted file the answer can be incomplete.
- **Old behaviour:** `bedextract` always trusts the order, and it silently
  misses a long feature that starts before the region and any fully nested
  feature; gloci's search finds them. bedtools and tabix read an index or the
  whole file.
- **Pinned by:** `view-bed-unsorted-vs-bedextract` (unsorted file: gloci exit
  status 1, `bedextract` exit status 0), `view-bed-assume-sorted-nested-vs-bedextract`
  (must differ) and `view-bed-assume-sorted-vs-bedextract` (sorted, non-nested
  file: must match).
- **Reason:** correct by default, fast on request. Index the file
  (`gloci sort FILE -o FILE.gz --write-index`) to get speed and checking.

## Indexes built by htslib and SVs with only `SVLEN` {#htslib-index-svlen}

- **gloci:** the footprint of a symbolic SV without INFO/END is
  `POS + |SVLEN|`, and the index spans gloci builds include it. Queries on a
  gloci-built index find the record from any region it touches, and htslib reads
  that index correctly (`index-interop-tabix-queries`).
- **Old behaviour:** htslib ignores SVLEN. An index built by `tabix` or
  `bcftools index` spans such a record only over its REF bases, so a region that
  touches only the tail of an SVLEN-only SV (in another index bin) is missed by
  `tabix`, and also by gloci when it reads that htslib-built index.
- **Reason:** the SVLEN-derived end is part of the footprint rules. Rebuild the
  index with `gloci index` if a file has SVLEN-only SVs and you query regions
  inside them. Records with INFO/END are unaffected.
- **Pinned by:** `view-svlen-only-sv-vs-tabix` (scan: gloci finds the SV,
  tabix does not), and live only: `index-svlen-only-sv-vs-tabix` (gloci-built
  index against a tabix-built one), `index-htslib-built-svlen-gloci-view-vs-tabix`
  (must match: gloci misses it through an htslib-built index, like tabix).

## VCF headers noodles rejects cannot be converted {#vcf-header-conversion}

- **gloci:** `view` between VCF and BCF parses the header with the `noodles`
  library for record encoding. A header it rejects, for example a duplicated
  INFO ID, stops the conversion: `the VCF header cannot be used for
  conversion (parse: duplicate INFO ID: X); fix the header lines it names`.
  Plain VCF streaming and sorting of the same file do not parse the header
  that way and still work.
- **Old behaviour:** `bcftools view -Ob` 1.21 reports no error for a duplicated
  INFO ID.
- **Reason:** gloci does not convert records against a header it cannot
  interpret, because a wrong BCF dictionary gives wrong records. Fix the
  header (the error names the line) and convert again.
- **Pinned by:** `vcf-to-bcf-bad-header-vs-bcftools` (a duplicated INFO ID:
  gloci exit status 1, `bcftools view -Ob` exit status 0).

## Zero-length BED records {#zero-length-bed-records}

- **gloci:** accepted, and treated as points.
- **Old behaviour:** `sort-bed` rejects them (`Genomic end coordinate is less
  than (or equal to) start coordinate`); `bedtools sort` and GNU sort accept
  them. The other sort-bed cases therefore use an input without points
  (`ties-nopoints.bed`).
- **Pinned by:** `sort-bed-points-vs-sort-bed` (gloci exit status 0, `sort-bed`
  exit status 1). The bedtools and GNU sort cases include points and must match.

## `sort` refuses records on undeclared contigs {#undeclared-contig-sort}

- **gloci:** `sort`, `sort --check` and `sort -m` exit 1 on a VCF whose
  records use a contig that its `##contig` lines omit, naming the line and
  saying to add `##contig=<ID=...>` to the header or to pass `-g` with a genome
  that lists it. Writing such a record to BCF (`view -O bcf`) is refused the
  same way, because BCF needs every contig in the header. `view` (with or
  without regions, with `--write-index`) and `index` accept the records, with
  one warning per contig on stderr, as htslib does.
- **Old behaviour:** `bcftools sort` accepts the file and appends the contig
  to the order, in order of first appearance.
- **Reason:** a sort needs a chromosome order, and a contig the header does not
  declare has none, so appending it would silently pick an order the
  header does not state. Add the `##contig` line, or pass `-g`, and the sort
  is exact.
- **Pinned by:** `sort-undeclared-contig-vs-bcftools-sort` (exit status: must
  differ). `view-region-undeclared-contig-vs-tabix` and
  `view-region-undeclared-contig-vs-bcftools-view` pin the agreeing side.

## `index` refuses a BCF whose contigs are out of dictionary order {#bcf-index-contig-order}

- **gloci:** `gloci index x.bcf` exits 1 when the records of one contig are
  contiguous but the contigs are not in the order of the header's contig
  dictionary (for example `chr2` records, then `chr1` records, with `chr1`
  declared first). The error names the first out-of-order record and says to
  run `gloci sort x.bcf`.
- **Old behaviour:** `bcftools index` accepts such a file and indexes each
  contig by its dictionary id.
- **Reason:** the indexer gloci uses needs reference ids to increase through the
  file. A BCF in dictionary order, which is what `gloci sort` writes, indexes
  and queries the same as with htslib.
- **Pinned by:** `index-bcf-contigs-out-of-dictionary-order-vs-bcftools` (exit
  status: must differ).

## `index -s -b` without `-e` writes the begin column as the end column {#tsv-index-no-end-column}

- **gloci:** `gloci index -s 1 -b 2 FILE` writes `col_end` = the begin column
  (2) in the tabix header. A feature covers one base.
- **Old behaviour:** `tabix -s1 -b2 FILE` leaves `col_end` at its built-in GFF
  default, 5, and `tabix -e 0` is refused when the file is 1-based.
- **Reason:** 5 is an accident: a file with a fifth column would have it read
  as the end. An end column equal to the begin column is what htslib itself
  treats as "no end column", so queries return the same lines (the
  `index-tsv-skip-lines-vs-tabix` case pins this).
- **Pinned by:** `index-tsv-no-end-header-vs-tabix` (header words: must
  differ).

## Points overlap an interval only when inside it {#zero-length-overlap}

- **gloci:** a zero-length record (a point, such as an insertion site) overlaps
  an interval only when it lies strictly inside it, and two points overlap
  only when they are equal. `--inclusive-points` also lets a point on the
  edge of an interval overlap it. This applies to every relation command, so
  `intersect` keeps a point on B's edge out of `--mode any`, `pairs` and
  `count`.
- **Old behaviour:** `bedtools intersect` treats a point on either edge of an
  interval as overlapping it (`chr1 100 100` against `chr1 100 200` and
  `chr1 200 200` against the same interval are both reported). That is the
  result of `--inclusive-points`.
- **Reason:** a point on the edge sits between two bases, so it does not lie
  within the interval's bases; counting it makes an insertion next to a
  feature look like it hit the feature. Spec section 4.3 has the table.
- **Pinned by:** `intersect-points-vs-bedtools` (must differ: bedtools reports
  the two edge points, gloci does not). The `intersect` tests also check every
  row of the spec table, with and without `--inclusive-points`.

## `intersect` on a VCF uses the allele span {#vcf-allele-span}

- **gloci:** a VCF A or B record covers the bases its ALT alleles change
  (see [Default `--vcf-span allele`](#default-vcf-span-allele)). An insertion
  is a point and a deletion covers only the deleted bases, so a region that
  holds only the anchor base of an indel does not select it. `--vcf-span ref`
  gives the old result.
- **Old behaviour:** `bedtools intersect` uses `[POS-1, POS-1+len(REF))`, so a
  region over the anchor base alone selects the insertion and the deletion.
- **Reason:** the same as the default span: the anchor base is not changed by
  the variant.
- **Pinned by:** `intersect-vcf-allele-span-vs-bedtools` (must differ: bedtools
  also reports the insertion and the deletion whose anchor base is in a
  region) and `intersect-vcf-ref-span-vs-bedtools` (`--vcf-span ref`: must
  match). `subtract` on a VCF A removes the same whole records as
  `intersect --mode none`, so it has the same difference:
  `subtract-vcf-allele-span-vs-bedtools` (must differ, against
  `bedtools subtract -A`) and `subtract-vcf-ref-span-vs-bedtools` (must match).

## Overlap fractions are measured per B record {#per-record-overlap-fraction}

- **gloci:** `-f`, `-F`, `-r`, `-e` and `--min-bases` are measured for each B
  record separately, as bedtools does. A that is covered 80% by three B
  records of 12%, 28% and 40% does not meet `-f 0.5` for any of them.
- **Old behaviour:** `bedops --element-of 50%` adds up the coverage of all
  B elements, so that A is kept.
- **Reason:** per-record thresholds are what `bedtools intersect` users have,
  and `--mode pairs` can then report which B records met them. To count the
  summed coverage, run `gloci merge` on B first.
- **Pinned by:** `intersect-fraction-vs-bedops` (must differ: bedops keeps a
  record gloci drops). `intersect-f-vs-bedtools` and the other threshold cases
  match bedtools.

## `--match alleles` compares trimmed alleles {#trimmed-allele-match}

- **gloci:** `intersect --match alleles` keeps a record when some ALT of it
  equals some ALT of a B record after trimming the shared prefix and suffix
  of REF and ALT, at the same trimmed position. `A>AG` at 100 matches `C>GC`
  at 101 (an insertion of G), and records match when they share any one ALT of
  several. The overlap options (`-f`, `-s`, `-w`, ...) are refused with it.
- **Old behaviour:** `bcftools isec -c none` pairs records at the same
  position with identical REF and ALT lists, and does not normalise.
- **Reason:** the same variant is often written with different anchors, and a
  multi-allelic record shares a variant with another that lists more or fewer
  ALTs. Run `bcftools norm` on both inputs and the results then agree for
  records with identical ALT lists.
- **Pinned by:** `intersect-alleles-unnormalised-vs-bcftools` (must differ);
  `intersect-alleles-vs-bcftools` and `intersect-alleles-multiallelic-vs-bcftools`
  (must match).

## `merge -s` keeps records of unknown strand {#merge-unknown-strand}

- **gloci:** `merge -s` clusters `+` and `-` records separately and clusters
  records of unknown strand (`.`, no strand column, VCF records) together, as
  a third group, so they appear in the output (the pinned case has `.` in
  column 6). It warns once per input that has such records. Output is sorted by chromosome, start, end.
- **Old behaviour:** `bedtools merge -s` writes nothing for records whose
  column 6 is `.`; they vanish from the output. (Other inputs, such as BED3,
  were not compared.)
- **Reason:** silently dropping records is a data-loss trap, and a feature of
  unknown strand is still a feature. To drop them anyway, filter the input
  first (`awk '$6 != "."'`).
- **Pinned by:** `merge-strand-vs-bedtools` (must differ: bedtools drops the
  `.`-strand records `chr1 600 700` and `chr1 690 720`). The other `merge`
  cases match bedtools and bedops.

## `subtract --sum-overlaps` removes A when the sum reaches `-f` {#sum-overlaps-threshold}

- **gloci:** `--sum-overlaps -f F` removes an A record when the bases shared
  with the B records that pass the strand rule add up to at least `F` of A,
  the same "at least" every other threshold uses. Without `-f` it removes A
  when any such B overlaps it, like `--whole`.
- **Old behaviour:** `bedtools subtract -N -f F` removes A only when the sum
  is more than `F` of A: an A covered exactly 50% is kept by `-f 0.5`. Without
  `-f`, `-N` is refused with "-f must be in the range (0.0, 1.0]".
- **Reason:** one rule for "at least" across all thresholds (`-f 0.5` on
  `intersect` accepts a 50% overlap), and a command that works without `-f`.
  Raise `-f` slightly (for example to `0.51`) to get the bedtools result.
- **Pinned by:** `subtract-sum-overlaps-equal-vs-bedtools` (must differ: A
  `chr1 600 700` is covered exactly 50% and bedtools keeps it).
  `subtract-sum-overlaps-vs-bedtools` (`-f 0.6`) must match.

## `closest --each-side` counts an overlap on both sides {#closest-each-side}

- **gloci:** `closest --each-side` (the replacement for `bedops
  closest-features`) treats an overlapping B as the nearest B on both sides,
  at distance 0, unless `--ignore-overlaps` is given. A side with no B is the
  D45 placeholder (`. -1 -1 ...`, distance `-1`), as in plain `closest`. `-k`
  and `--ties` apply per side, and several equally near B on one side add
  rows instead of choosing one. With `--distance` the distances are signed
  in the `--signed` frame (default `ref`). With `--ignore-upstream` or
  `--ignore-downstream`, an overlapping B (sign 0) still shows in the ignored
  column at distance 0, since an overlap is neither upstream nor downstream;
  add `--ignore-overlaps` to drop it.
- **Old behaviour:** `closest-features` gives an overlapping element
  precedence on one side and reports the nearest non-overlapping element on
  the other ("elements reported as closest to the left and right edges are
  never the same"); it writes `NA` for a missing side, and always one row.
- **Reason:** counting an overlap on both sides is the same rule as
  `closest` without `--each-side`, and the placeholder keeps the columns of
  every row the same width.
- **Pinned by:** `closest-each-side-overlap-vs-bedops` and
  `closest-each-side-missing-vs-bedops` (must differ). With
  `--ignore-overlaps` and every side present, `closest-each-side-vs-bedops`
  and `closest-each-side-nodist-vs-bedops` match exactly.

## `closest -k N --ties first|last` cuts only the ties at the N-th distance {#closest-ties-k}

- **gloci:** `-k N` reports the N nearest B and every B tied with the N-th.
  `--ties first|last` keeps, among the B tied at the N-th distance, only as
  many as N still needs (the first or last by file and input order). With
  `-k 1` that is exactly one row, as in bedtools.
- **Old behaviour:** `bedtools closest -k N -t first|last` collapses the
  ties of every distance to one B and counts distinct distances, so two B at
  the same distance below the N-th use up one of the N places only once and a
  farther B follows. Its default `-t all` agrees with gloci.
- **Reason:** gloci's rule keeps the output to N rows without ties and is
  the one `-t all` already follows; bedtools' counting of distinct distances
  needs an unbounded look at farther B.
- **Pinned by:** `closest-k2-ties-first-vs-bedtools` (must differ: the
  equally near `chr1 450 500` pair gives `b6` and `b7` in gloci, `b6` and
  the farther `b5` in bedtools). `-k 1 --ties first|last` and `-k 2` with
  `--ties all` match exactly.
