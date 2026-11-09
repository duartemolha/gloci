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

## Symbolic ALTs with neither `END` nor `SVLEN` {#vcf-symbolic-no-end}

- **gloci:** under the default `--vcf-span allele`, a record whose ALT is a
  symbolic SV allele (`<DEL>`, `<DUP>`, `<CNV>`, ...) with neither INFO/END nor
  INFO/SVLEN has no span, so it is an invalid record: an error that says to
  "add INFO/END or INFO/SVLEN, or use --vcf-span ref, or --skip-invalid"
  (`--skip-invalid` skips it with a warning, like other invalid records). Under
  `--vcf-span ref` (and in `sort`, which orders by the REF end) the span is the
  REF span and the record is kept. `<INS>`, `<*>` and `<NON_REF>` never need
  either field.
- **Old behaviour:** bcftools, tabix and htslib use the REF span and accept
  such records.
- **Reason:** the allele span of a deletion or duplication of unknown length
  cannot be guessed, and a silent guess would change what overlaps it.
- **Pinned by:** `view-symbolic-no-end-default-vs-bcftools-view` (exit status
  differs: gloci 1, bcftools 0); `view-symbolic-no-end-ref-vs-bcftools-view`
  and `vcf-filter-symbolic-no-end-ref-vs-bcftools-view` (with
  `--vcf-span ref` gloci matches bcftools).

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
  match). In default or `--whole` mode `subtract` on a VCF A removes the same
  whole records as `intersect --mode none` (under `--sum-overlaps -f` it drops
  a record when the summed overlaps reach `-f` of its footprint), so it has the
  same difference: `subtract-vcf-allele-span-vs-bedtools` (must differ, against
  `bedtools subtract -A`) and `subtract-vcf-ref-span-vs-bedtools` (must match).
- **`convert`:** the same holds for `gloci convert`. Its default rows start a
  deletion after its anchor base where `vcf2bed` starts at `POS-1`
  (`--vcf-span ref` gives `vcf2bed`'s start), an insertion becomes a
  zero-width row (`400 400`) where `vcf2bed` writes its anchor base, and
  `--split-alleles` gives each ALT its own allele footprint by default (while
  honouring `--vcf-span ref`). Pinned by `convert-vcf-allele-span-vs-vcf2bed`
  and `convert-vcf-split-deletions-vs-vcf2bed` (must differ);
  `convert-vcf-ref-vs-vcf2bed` (`--vcf-span ref`: must match).

## `intersect --mode clip` of a GFF/GTF feature by a point {#clip-gff-points}

- **gloci:** a clip piece is the shared span of A and B. When B is a point
  (zero-length) inside a GFF/GTF feature, that span is empty, and GFF/GTF
  cannot hold an empty feature (its end would be less than its start), so the
  piece is not written; one warning gives the count of such pieces and
  suggests `--mode any` or `--mode pairs`. A BED A keeps its zero-length
  piece (`chr1 15 15`), since BED can hold it.
- **Old behaviour:** `bedtools intersect` 2.31.1 treats the point as one base
  and writes a two-base feature (`chr1 src gene 15 16 ...`) for a point at
  BED position 15.
- **Reason:** gloci's points have no bases (spec §4.3), so a feature clipped
  to one has no extent; writing a base the point does not cover would be a
  wrong answer, and writing `16 15` would be invalid GFF.
- **Pinned by:** `intersect-clip-gff-points-vs-bedtools`.

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

## `map` `stdev` and `var` are sample statistics {#map-stdev-sample}

- **gloci:** `stdev` and `var` divide by n - 1. With one value they print the
  `--empty` value (default `.`), since a sample of one has no spread.
  `stdev_pop` and `var_pop` divide by n and print 0 for one value.
- **Old behaviour:** `bedtools map -o stdev` is the population
  form, dividing by n, and gives 0 for one value. `sstdev` is the
  sample form. bedmap's `--stdev` is the population form as well.
- **Reason:** the sample form is what a user summarising a handful of
  overlapping records expects, and it agrees with bedtools' `sstdev`. Use
  `stdev_pop` / `var_pop` to reproduce bedtools' `stdev` / `var`.
- **Pinned by:** `map-stdev-vs-bedtools` (must differ). `map-sstdev-vs-bedtools`,
  and `map-stdev-pop-vs-bedtools` match (bedtools map has no `var` op).

## `map` `distinct` keeps first-seen order {#map-distinct-order}

- **gloci:** `distinct` lists each value once in the order the B records were
  seen (B file, then input order). `distinct_sorted` lists them in bytewise
  order.
- **Old behaviour:** `bedtools map -o distinct` lists them in bytewise order
  (`10,3,5`).
- **Reason:** first-seen order follows the records, which makes `distinct`
  line up with `collapse`; the sorted form is kept under its own name.
- **Pinned by:** `map-distinct-vs-bedtools` (must differ).
  `map-distinct-sorted-vs-bedtools` matches.

## `map` `mode` and `antimode` skip missing values {#map-mode-missing}

- **gloci:** `mode` and `antimode` skip `.` and empty values, like every other
  op, so `.` and `3` give `3`. Ties between equally frequent values go to the
  bytewise smallest value, as in bedtools.
- **Old behaviour:** `bedtools map -o mode` counts `.` as a value, so `.` and
  `3` give `.`.
- **Reason:** missing is missing. bedtools can report `.` as the mode while real
  values exist.
- **Pinned by:** `map-mode-missing-vs-bedtools` and
  `map-antimode-missing-vs-bedtools` (must differ). `map-mode-vs-bedtools`,
  `map-antimode-vs-bedtools`, `map-mode-ties-vs-bedtools` and
  `map-antimode-ties-vs-bedtools` match.

## `groupby --header` writes the column names {#groupby-header}

- **gloci:** with `--header` the first line names the columns, so `-g` and `-c`
  accept names, and the output starts with a line holding the key names (with the
  input's leading `#`, if it had one) and the `-c` labels as typed, such as
  `#gene` and `expr:sum`. Only line 1 is the header: later lines that start with
  `#` are comments and are skipped.
- **Old behaviour:** `bedtools groupby -header` copies the input header lines
  unchanged (every leading `#` line) and does not accept names, so the header has
  the input's columns over the output's.
- **Reason:** the copied header names columns that are no longer there.
- **Pinned by:** `groupby-header-vs-bedtools` (must differ).

## A missing column is an error, not a message and exit 0 {#groupby-missing-column}

- **gloci:** a line with fewer tab-separated columns than a `-g` key or a `-c`
  column needs is an error naming the file, the line, the column and the line's
  column count, with exit status 1. An empty value or `.` is still a missing
  value and is skipped. `map -c` does the same for a B line.
- **Old behaviour:** `bedtools groupby` prints a message that names the column
  but not the line, stops with the output so far, and exits 0.
- **Reason:** a pipeline must not take truncated output for a result.
- **Pinned by:** `groupby-missing-column-vs-bedtools` and
  `groupby-absent-column-vs-bedtools` (must differ, by exit status).
  `map-absent-column-vs-bedtools` matches: both tools fail.

## A value that is not a number is an error in `map` and `groupby` {#agg-non-numeric}

- **gloci:** a numeric op (`sum`, `mean`, `median`, ...) given a value that is
  not a number stops with an error naming the file, line and value.
- **Old behaviour:** `bedtools map` and `bedtools groupby` print a warning, skip
  the value and exit 0.
- **Reason:** skipping hides bad data, such as a header line read as data.
- **Pinned by:** `groupby-non-numeric-vs-bedtools` and
  `map-non-numeric-vs-bedtools` (must differ, by exit status).

## `map` prints the empty value, not `NAN`, for no values {#map-empty-value}

- **gloci:** an aggregation with no values prints `--empty` (default `.`);
  `count` and `count_distinct` print 0.
- **Old behaviour:** `bedmap --sum` and `--mean` print `NAN` for an element
  with no mapped records. bedtools map prints `.` like gloci.
- **Reason:** one placeholder for every operation, chosen by the user; the
  bedtools default is the same.
- **Pinned by:** `map-empty-sum-vs-bedops` (must differ). With every A
  overlapped, `map-sum-vs-bedops` and `map-mean-vs-bedops` match.

## `coverage` of a point A has no bases {#coverage-point}

- **gloci:** a point A (a zero-length feature) covers no bases. `coverage`
  writes the number of B that overlap it (spec 4.3), 0 covered bases,
  length 0 and fraction `.`. `--mean` writes `.`, and `--hist` and
  `--per-base` write no rows for it. Points add nothing to the `all` rows.
- **Old behaviour:** `bedtools coverage` treats the point as two bases long:
  `chr1 10 10` with an overlapping B reports 2 covered bases of length 2
  (fraction 1), two `-d` rows and a `-hist` row at depth 1.
- **Reason:** a feature with no bases has no covered fraction, and counting
  invented bases would skew the `all` totals.
- **Pinned by:** `coverage-point-default-vs-bedtools`,
  `coverage-point-hist-vs-bedtools`, `coverage-point-per-base-vs-bedtools` and
  `coverage-point-mean-vs-bedtools` (must differ).

## `coverage --split` measures A by its blocks {#coverage-split-length}

- **gloci:** with `--split`, depth counts the bases shared by a block of A
  and a block of B, and A's length is the sum of A's blocks, so the
  fraction is covered block bases over block bases. `--hist` has no depth-0
  row for the gap between blocks, and `--per-base` offsets run along the
  blocks.
- **Old behaviour:** `bedtools coverage -split` counts the same covered
  bases but divides by A's whole span, so introns count as uncovered, and
  `-hist` reports them at depth 0. `bedtools coverage -d -split` ignores
  `-split` altogether: it walks A's whole span and counts B's whole spans, so
  a base in the gap between A's blocks gets a row and can have depth.
- **Reason:** for a spliced feature the fraction of its own bases covered
  is the useful number, and it agrees with the block lengths that
  `-f`/`-F` use under `--split`.
- **Pinned by:** `coverage-split-vs-bedtools`,
  `coverage-split-hist-vs-bedtools` and `coverage-split-per-base-vs-bedtools`
  (must differ).

## `complement` ignores points {#complement-points}

- **gloci:** a point (a zero-length feature) covers no bases, so it is ignored:
  it neither splits a gap nor makes a chromosome appear (D48).
- **Old behaviour:** `bedtools complement` treats the point at `p` as a base
  and splits the gap into `[.., p-1)` and `[p+1, ..)` (`chr1 200 299`,
  `chr1 301 500` for a point at 300). A point at 0 gives a row of negative
  width (`chr10 0 -1`, then `chr10 1 300`). `bedops --complement` splits the
  gap at the point into two abutting rows (`chr1 200 300`, `chr1 300 500`).
- **Reason:** a feature with no bases leaves no base uncovered or covered, and
  inventing one would remove a base from the output.
- **Pinned by:** `complement-point-vs-bedtools` and
  `complement-point-vs-bedops` (must differ).

## `complement` rejects a feature past the chromosome end {#complement-past-end}

- **gloci:** with `-g`, a feature that ends beyond its chromosome's length is
  an error (exit 1) naming the file and line, with a hint to check the genome
  build.
- **Old behaviour:** `bedtools complement` prints `WARNING: ... exceeds the
  length of chromosome` and exits 0, writing the gaps it can.
- **Reason:** a feature past the end almost always means the genome build
  does not match the input, and the output would silently be wrong.
- **Pinned by:** `complement-past-end-vs-bedtools` (must differ; exit status).

## `partition` cuts at every boundary unless asked to merge {#partition-boundaries}

- **gloci:** `partition` cuts the inputs at every start and end, like
  `bedops --partition`, so neighbouring segments with the same membership stay
  separate (`chr1 0 10` and `chr1 10 20`, both covered by the same files).
  `--merge-adjacent` joins them.
- **Old behaviour:** `bedtools multiinter` joins neighbouring segments whose
  membership is equal, even where a feature ends and another begins;
  `bedtools unionbedg` and `bedops --partition` do not.
- **Reason:** one default for the whole family; the segment list of
  `bedops -p` is the finest one and `--merge-adjacent` derives the coarser
  ones (`-i` and `-s` need it, and so does `multiinter`).
- **Pinned by:** `partition-boundaries-vs-bedtools-multiinter` (must differ);
  `partition-merge-adjacent-vs-bedtools-multiinter` and
  `partition-min-files-vs-bedtools-multiinter` match with `--merge-adjacent`.

## `partition` ignores points {#partition-points}

- **gloci:** a point (a zero-length feature) covers no bases, so it adds no
  segment and no boundary (D68).
- **Old behaviour:** `bedtools multiinter` widens an isolated point at `p` to
  `[p-1, p+1)` (a point inside an interval is ignored, so it adds no row
  there), and
  `bedops --partition` writes it as a zero-length segment (`chr1 150 150`) and
  cuts the segments around it there.
- **Reason:** a feature with no bases covers nothing; inventing bases changes
  coverage counts.
- **Also differs:** default labels are file basenames in gloci; `multiinter`
  numbers the inputs (`1,2,...`) unless `-names` is given.
- **Pinned by:** `partition-points-vs-bedtools-multiinter` and
  `partition-points-vs-bedops-partition` (must differ).

## `partition --bedgraph --empty -g` covers chromosomes without features {#partition-empty-chromosomes}

- **gloci:** with `--empty VALUE` and `-g GENOME` (a genome with lengths), every
  chromosome of the genome is written over `[0, length)`, including those no
  input has a feature on.
- **Old behaviour:** `bedtools unionbedg -empty -g` fills the gaps of the
  chromosomes that have a feature in some input and omits the others.
- **Reason:** "the whole genome" should not depend on which chromosomes the
  inputs happen to touch; the same rule as `complement -g`.
- **Pinned by:** `partition-bedgraph-empty-genome-vs-bedtools-unionbedg` (must
  differ).

## `slop` and `flank` accept one side alone {#slop-flank-one-side}

- **gloci:** `--left N` or `--right N` alone pads only that side.
- **Old behaviour:** `bedtools slop` and `bedtools flank` exit 1 with "Need
  both -l and -r" unless `-b` is given.
- **Reason:** a one-sided promoter or flank is the common use, and writing
  `-r 0` to say "nothing" only invites mistakes.
- **Pinned by:** `flank-one-side-vs-bedtools` and `slop-one-side-vs-bedtools`
  (must differ; exit status).

## `flank` of a point {#flank-points}

- **gloci:** a point `[p, p)` has its strips `[p - L, p)` and `[p, p + R)`,
  as `slop` pads it to `[p - L, p + R)`.
- **Old behaviour:** `bedtools flank` writes `[p - L - 1, p - 1)` and
  `[p + 1, p + R + 1)`, leaving a base between the point and each strip.
  `bedtools slop` agrees with gloci.
- **Reason:** a point sits between bases (spec §4.3); the strips should touch
  it.
- **Pinned by:** `flank-point-vs-bedtools` (must differ). `slop-point-vs-bedtools`
  must match.

## `slop` and `flank` reject a feature beyond the chromosome end {#slop-flank-past-end}

- **gloci:** with `-g`, a feature that starts at or after the end of its
  chromosome (and so cannot be clipped to anything) is an error (exit 1)
  naming the file and line, with a hint to check the genome build. A feature
  that only straddles the end is clipped, as in bedtools.
- **Old behaviour:** `bedtools slop` writes a 1-base feature at the very end
  of the chromosome (`[len - 1, len)`) and `bedtools flank` writes a left
  strip that lies beyond the end, both exit 0.
- **Reason:** a feature past the end almost always means the genome build
  does not match the input; the output would be meaningless (and an empty
  GFF record is not representable).
- **Pinned by:** `slop-past-end-vs-bedtools` and `flank-past-end-vs-bedtools`
  (must differ; exit status).

## `shift` keeps a zero-length feature valid at a chromosome edge {#shift-past-edge}

- **gloci:** a zero-length feature (`[p, p)`) shifted wholly past 0 or (with
  `-g`) the chromosome end is written zero-length at the edge (`[0, 0)` or
  `[len, len)`). A non-empty feature pushed wholly over an edge becomes the
  1-base interval there (`[len - 1, len)` or `[0, 1)`), as in bedtools.
- **Old behaviour:** `bedtools shift` writes an inverted interval such as
  `1 0` (or `1000 999` at the end) for a zero-length feature pushed over an
  edge, and exits 0.
- **Reason:** an inverted interval is invalid and breaks every later tool.
- **Pinned by:** `shift-point-edge-vs-bedtools` (must differ).

## `shift --pct` uses exact decimals {#pct-exact-decimal}

- **gloci:** `shift --pct --by F` moves a feature of length L by exactly
  `floor(F * L)` bases (F read as an exact decimal, so `-0.1` on 10 bases is
  -1 and `0.7` is 7). A zero-length feature does not move (0 * F = 0). A
  non-empty feature whose shift leaves nothing of it inside the chromosome
  is written as the 1-base edge interval (`[0, 1)` for `[0,5)` shifted by
  `-0.9`: floor(-4.5) is -5).
- **Old behaviour:** `bedtools shift -pct` computes in single-precision
  float, so non-dyadic fractions can land one base away (10 bases: -0.1
  gives -2, -0.3 gives -4, 0.7 gives 6, 0.9 gives 8). It treats a
  zero-length feature as 2 bases long (`[10,10)` with -0.5 gives `9 9`),
  and adds the fractional shift to each coordinate before truncating, so
  `[0,5)` with -0.9 gives `0 0`.
- **Reason:** results must not depend on float rounding; `slop` and `flank`
  `--pct` already use exact decimals (and agree with bedtools on the
  fractions tried, 0.1 to 1.1).
- **Pinned by:** `shift-pct-decimal-neg-vs-bedtools`,
  `shift-pct-decimal-pos-vs-bedtools`, `shift-pct-point-vs-bedtools` and
  `shift-pct-edge-vs-bedtools` (must differ).

## `flank -s` on a record of unknown strand {#flank-unknown-strand}

- **gloci:** with `-s`, a record whose strand is `.` (or missing) counts as
  `+`, as in `slop -s` (D49).
- **Old behaviour:** `bedtools flank -s` writes no strip at all for such a
  record (`bedtools slop -s` treats it as `+`).
- **Reason:** silently dropping features is worse than a defined default, and
  the two commands should agree.
- **Pinned by:** `flank-unknown-strand-vs-bedtools` (must differ).

## `windows` of a point {#windows-points}

- **gloci:** a point `[p, p)` has no bases, so `windows -i` writes nothing
  for it.
- **Old behaviour:** `bedtools makewindows -b` reads a zero-length record as
  the single base at 1-based position `p` and writes the window `[p - 1, p + 1)`
  for it (with `-w 4`).
- **Reason:** a point sits between bases (spec §4.3) and covers nothing to
  cut.
- **Pinned by:** `windows-point-vs-bedtools` (must differ).

## `windows -g` skips a zero-length chromosome {#windows-zero-length-chromosome}

- **gloci:** a chromosome of length 0 in `-g` gets no windows; one
  `warning:` names it and the run goes on.
- **Old behaviour:** `bedtools makewindows -g` exits 1 with "Each chromosome
  must have non-zero length".
- **Reason:** zero-length contigs occur in `.fai` and `.dict` files, and one
  of them should not fail a genome-wide run.
- **Pinned by:** `windows-zero-length-chromosome-vs-bedtools` (must differ;
  exit status).

## `convert` writes one row per VCF record {#convert-multiallelic-rows}

- **gloci:** `gloci convert` writes one row per VCF/BCF record, with the ALT
  column as it is (`A,AC`) and the footprint of the whole record.
  `--split-alleles` gives one row per ALT instead.
- **Old behaviour:** `vcf2bed` always writes one row per ALT allele.
- **Reason:** one row per record keeps the record count and the ALT column
  intact for a plain conversion; spec §4.2 makes the per-ALT rows an explicit
  choice.
- **Pinned by:** `convert-vcf-multiallelic-vs-vcf2bed` (must differ). With
  `--split-alleles` and SNV alleles the output matches:
  `convert-vcf-split-alleles-vs-vcf2bed`.

## `convert` spans variants by their alleles, `vcf2bed` by their lengths {#vcf2bed-variant-spans}

- **gloci:** with `--vcf-span ref` an MNV covers its whole REF, a symbolic ALT
  runs to INFO/END (or POS plus SVLEN) and a breakend covers its REF base.
  `--types` follows the variant kinds of spec §4.2, decided on the alleles
  with their shared prefix and suffix trimmed: `ATG` to `AC` (`TG` to `C`) is
  `other`, a breakend is `sv`.
- **Old behaviour:** `vcf2bed` gives an MNV (`AC` to `GT`) a single base, a
  symbolic ALT such as `<DEL>` its single anchor base whatever END says, and
  a deletion `len(REF) - len(ALT) + 1` bases (so `ATG` to `AC` gets two) with
  `--deletions` but only its anchor base without it.
  `--deletions` and `--insertions` pick records by the lengths of REF and
  ALT: `ATG` to `AC` is a deletion and a breakend is an insertion.
- **Reason:** the footprint is what the variant changes on the reference
  (spec §4.2), and a type should say what the variant is, not how long its
  strings are.
- **Pinned by:** `convert-vcf-spans-vs-vcf2bed`,
  `convert-vcf-ref-deletions-vs-vcf2bed`,
  `convert-vcf-types-del-vs-vcf2bed` and `convert-vcf-types-ins-vs-vcf2bed`
  (must differ). For single-base SNVs, plain insertions and deletions the
  output matches: `convert-vcf-ref-vs-vcf2bed`,
  `convert-vcf-snvs-vs-vcf2bed`, `convert-vcf-insertions-vs-vcf2bed` and
  `convert-vcf-deletions-vs-vcf2bed`.

## `convert` names GFF3 features by their ID {#gff2bed-name}

- **gloci:** the name column of a GFF3 row is the `ID` attribute (`.` when the
  feature has none); `--name-attr KEY` picks another attribute. A GTF row is
  named by `gene_id` by default, as `gtf2bed` does.
- **Old behaviour:** `gff2bed` 2.4.41 writes `.` in the name column of every
  GFF3 row.
- **Reason:** a BED name of `.` throws away the one identifier GFF3 features
  have.
- **Pinned by:** `convert-gff3-name-vs-gff2bed` (must differ). With a
  `--name-attr` the feature does not have, every other column matches:
  `convert-gff3-vs-gff2bed`.

## `convert --blocks` rejects a record that is not valid BED12 {#bed12-short-line}

- **gloci:** a record with fewer than 12 columns, a first block that does not
  start at offset 0, blocks that overlap or are out of order, or a last block
  that does not end at chromEnd is a data error (exit 1) naming the file and
  the line. `--skip-invalid` does not skip it. A blockCount of 0 or one that
  disagrees with the block lists is an error too.
- **Old behaviour:** `bedtools bed12tobed6` writes a record with fewer than 12
  columns unchanged and exits 0, and for the other cases above (except a
  zero or mismatched blockCount, which it also rejects) writes whatever
  blocks the lists give.
- **Reason:** such a record has no well-formed blocks, and passing it through
  hides a wrong input file.
- **Pinned by:** `convert-bed12-short-line-vs-bedtools` and
  `convert-bed12-invalid-vs-bedtools` (must differ; exit status).
  `convert-bed12-zero-blocks-vs-bedtools` must match (both exit 1).

## INFO kept by default in `vcf filter` {#info-kept-by-default}

- **gloci:** `gloci vcf filter` keeps all INFO fields by default. `--keep-info`
  and `--drop-info` selectively prune them.
- **Old behaviour:** `vcftools --recode` removes all INFO fields unless
  `--recode-INFO-all` is explicitly passed.
- **Reason:** dropping data by default is unexpected in modern genomic pipelines;
  INFO annotations should be preserved unless explicitly dropped.
- **Pinned by:** `vcf-filter-info-default-vs-vcftools` (must differ). With
  `--recode-INFO-all` vcftools matches gloci: `vcf-filter-min-qual-vs-vcftools`.
  `bcftools view` also keeps INFO by default, agreeing with gloci:
  `vcf-filter-info-default-vs-bcftools` (must match).

## Inverted `--max-missing` parameter {#max-missing-inverted}

- **gloci:** `--max-missing FRAC` specifies the maximum allowed fraction of
  *missing* genotypes (between 0.0 and 1.0, where 0.1 means at most 10% missing).
- **Old behaviour:** `vcftools --max-missing FRAC` specifies the minimum allowed
  fraction of *present* (called) genotypes (where 0.9 means at least 90% called).
- **Reason:** In gloci, all `--max-*` flags specify upper bounds on what is
  being filtered; having `--max-missing` specify a minimum present rate was
  a historical vcftools naming inversion.
- **Pinned by:** `vcf-filter-max-missing-inverted-vs-vcftools` (must differ),
  `vcf-filter-max-missing-equiv-vs-vcftools` (gloci `--max-missing 0.1` matches
  vcftools `--max-missing 0.9`), and `vcf-filter-max-missing-vs-bcftools`
  (gloci `--max-missing 0.1` matches `bcftools view -i 'F_MISSING<=0.1'`).

## Minor allele frequency and count on multi-allelic sites {#maf-multiallelic}

- **gloci:** MAF is defined as the frequency of the second most common allele
  among called alleles (spec §7.8). On biallelic sites this is `min(p, 1-p)`.
  Similarly, `--min-mac` / `--max-mac` test the count of the second most common allele.
- **Old behaviour:** `vcftools --maf` and `--mac` evaluate multi-allelic sites
  against the least frequent alternate allele.
- **Reason:** The second most common allele is the standard mathematical
  definition of minor allele frequency for arbitrary ploidy and allele counts.
  Note that `bcftools view -q` filters on alternate allele frequency (frequency
  of the alternate allele or alleles via `-m`/`-M`), not on minor allele frequency
  in the second-most-common sense.
- **Pinned by:** `vcf-filter-maf-multiallelic-vs-vcftools` and
  `vcf-filter-min-mac-vs-vcftools` (must differ).

## `--min-maf` on a site with no called genotype {#all-missing-site-stats}

- **gloci:** a site where every genotype is missing has no called allele, so
  its MAF and MAC are undefined, and an undefined statistic fails every
  threshold on it: `--min-maf` and `--max-maf` both drop the site, as do
  `--min-mac`/`--max-mac`, and `--min-mean-dp`/`--max-mean-dp` drop a site with
  no FORMAT/DP. `--hwe` keeps a site it cannot test. `--phased-only` reads each
  genotype's separators, called or not, as vcftools does: a site of `./.`
  no-calls is unphased and dropped, one of `.|.` is kept.
- **Old behaviour:** `vcftools --maf 0.1` 0.1.16 keeps such a site, while its
  `--max-maf` drops it. `bcftools view -q 0.1:minor` 1.21 drops it, as gloci
  does.
- **Reason:** one rule for every threshold on a statistic the site does not
  have; a site with no data is not a common variant. vcftools' two MAF flags
  disagree with each other on it.
- **Pinned by:** `vcf-filter-all-missing-maf-vs-vcftools` (must differ) and
  `vcf-filter-all-missing-maf-vs-bcftools` (matches); the `--max-maf`,
  `--mac`, `--min-meanDP`, `--hwe` and `--phased` cases on the same input
  match vcftools.

## HWE filtering on multiallelic sites {#hwe-multiallelic}

- **gloci:** `--hwe P` evaluates Hardy-Weinberg equilibrium on diploid biallelic
  sites; multiallelic or non-diploid sites pass the filter (spec §7.8).
- **Old behaviour:** `vcftools --hwe` automatically removes all non-biallelic sites.
- **Reason:** Filtering non-biallelic sites should be an explicit user choice
  (via `--biallelic`) rather than an implicit side-effect of testing for HWE.
  Note that `bcftools view` has no standalone `--hwe` flag; `bcftools +fill-tags`
  computes an HWE tag only on biallelic diploid sites and leaves multiallelic sites
  unassigned.
- **Pinned by:** `vcf-filter-hwe-vs-vcftools` (must differ).

## BED files without headers in region filtering {#bed-header-skipping}

- **gloci:** `gloci vcf filter --regions-file` treats every line of a standard
  BED file as a region (unless it starts with `#`, `browser`, or `track`).
- **Old behaviour:** `vcftools --bed` unconditionally skips the first line of
  any BED file, assuming it is a header line. If the BED file has no header, the
  first interval is silently discarded.
- **Reason:** Standard BED files created by bedtools or pipelines rarely contain
  header lines; unconditionally discarding the first data row drops variants.
- **Pinned by:** `vcf-filter-regions-file-vs-vcftools` (must differ) and
  `vcf-filter-regions-file-vs-bcftools` (bcftools view -R does not skip line 1
  of headerless BED files, matching gloci).




## Filter expressions: `!~` negates the match {#expr-not-match}

- **gloci:** in `vcf filter --include/--exclude`, `FIELD !~ "re"` is true when
  the regular expression does not match, and (like `!=`) when the value is
  missing.
- **Old behaviour:** `bcftools view -i 'ID!~"^rs"'` (1.21) selects the same
  records as `ID~"^rs"`.
- **Reason:** the operator is documented as the negation of `~`; bcftools 1.21
  does not implement it.
- **Pinned by:** `expr-differs-from-bcftools-not-match` (must differ) and
  `difference_not_match_is_a_real_negation` in
  `crates/gloci-vcf/tests/expr.rs`.

## Filter expressions: an index past the end is missing {#expr-index-past-end}

- **gloci:** `INFO/AF[5]` on a record with fewer values is a missing value, so
  `INFO/AF[5]>0` is false.
- **Old behaviour:** bcftools 1.21 evaluates it on a record with a single `AF`
  value as if the index were in range and selects the record.
- **Reason:** an out-of-range index has no value to compare.
- **Pinned by:** `expr-differs-from-bcftools-index-past-end` (must differ).

## Filter expressions: a bare INFO flag is its presence {#expr-bare-flag}

- **gloci:** `INFO/DB` alone is true on records that carry the flag, so
  `--include INFO/DB` and `--exclude INFO/DB` work as written. `INFO/DB=1` is
  the same.
- **Old behaviour:** bcftools 1.21 selects no record for a bare flag
  (`-i 'INFO/DB'`), though `INFO/DB=1` works.
- **Reason:** a bare flag that matches nothing is never what the user means.
- **Pinned by:** `expr-differs-from-bcftools-bare-flag` (must differ).

## Filter expressions: `!` and the missing literal {#expr-grammar-extensions}

- **gloci:** `!` negates a parenthesised group (`!(QUAL<20 || DP<5)`) or a bare
  comparison (`!POS=10 && POS<30`). `.` is the missing value for text and
  numeric fields (`QUAL=.`, `QUAL!=.`, `ID=.`, `FILTER=.`, `INFO/AF=.`).
- **Old behaviour:** bcftools 1.21 fails with "Could not parse the expression"
  on `!(QUAL<20)`, `!POS=10` and `!(POS=10)`, and with `the tag "." is not
  defined` on `QUAL=.`, `INFO/AF=.` and `ID=.` (only `ID="."` works).
- **Reason:** these are the natural spellings; rejecting them gains nothing.
- **Pinned by:** `expr-differs-from-bcftools-negated-group`,
  `expr-differs-from-bcftools-missing-literal` and
  `expr-differs-from-bcftools-negated-comparison` (must differ).

## Filter expressions: a supported subset of bcftools {#expr-subset}

- **gloci:** the grammar of `vcf filter --include/--exclude` accepts one
  function (`N_PASS`), the FORMAT fields `GT`, `DP` and `GQ` (see
  [the FORMAT-field limitation](#expr-format-fields)), a field against a
  number, string or `.` (never two fields), and Rust `regex` syntax for `~`.
  Anything else is a usage error (exit 2) naming the token and its column. The
  sample-wise `&` and `|` are rejected with "use && / ||".
- **Old behaviour:** bcftools also has `ABS`, `STRLEN`, `MIN`, `MAX`, `SUM`,
  `COUNT`, other FORMAT fields, field-against-field comparisons and POSIX
  regular expressions.
- **Reason:** gloci implements the subset its users need, with bcftools
  semantics for everything it does accept.
- **Masked genotypes:** with `--min-gq`, `--min-dp`, `--max-dp` or
  `--mask-filtered-genotypes`, a masked sample reads as missing inside the
  expression (no GT, DP or GQ), while the output record still prints its
  original values. `--include 'FMT/GQ>=30'` after `--min-gq 30` therefore never
  sees a masked sample's GQ.
- **What matches bcftools 1.21 (not a difference):**
  - Sample subscripts: `FMT/GQ[1]`, `FMT/DP[1]` and `GT[1]` select sample 1
    (0-based), numbered over the whole header, so with `--samples` the
    subscript still means the original column. An index at or past the sample
    count is a usage error at parse time. Inside `N_PASS`, a subscripted field
    holds only for its own sample. For INFO fields `[i]` is an element index.
  - Bare `AC`, `AN` and `MAF` read INFO when the header defines the key
    (any-value for a multi-valued INFO such as `AC=3,2`) and are computed from
    the genotypes otherwise. `F_MISSING` and `N_MISSING` are always computed,
    even when the header defines them as INFO keys (`INFO/F_MISSING` reads the
    header key). An unindexed computed `AC` is the first ALT only, and missing
    when no allele is called.
  - A bare key defined as both INFO and FORMAT (such as `DP`) is an error that
    says to write `INFO/DP` or `FMT/DP`; bcftools also refuses it.
  - Any-sample FORMAT comparisons, `N_PASS`, the `GT` keywords, `TYPE`,
    `FILTER` and `ID` follow bcftools; each is pinned by an
    `expr-matches-bcftools-*` case.
- **Pinned by:** `expr-differs-from-bcftools-field-vs-field` and
  `expr-differs-from-bcftools-function-abs` (must differ), and
  `unsupported_constructs_are_usage_errors_with_hints` in
  `crates/gloci-vcf/tests/expr.rs`.

## Filter expressions: only the FORMAT fields GT, DP and GQ {#expr-format-fields}

- **gloci:** `FMT/GT`, `FMT/DP` and `FMT/GQ` (and bare `GT`, `DP`, `GQ`) are
  the only FORMAT fields an expression can read. `FMT/AD`, `FMT/FT`, `FMT/PL`
  and every other FORMAT key are a usage error, exit 2: "FMT/AD is not
  supported; supported FORMAT fields: GT, DP, GQ".
- **Old behaviour:** bcftools evaluates any FORMAT field, such as
  `bcftools view -i 'FMT/AD>1'`.
- **Reason:** gloci decodes only these three fields per sample. Selecting on
  another field is not possible yet; filter with bcftools for those.
- **Pinned by:** `expr-differs-from-bcftools-format-field-ad` (exit status
  differs: gloci 2, bcftools 0) and `unsupported_format_key_is_a_usage_error`.

## Filter expression `MAF` is not `--min-maf` {#expr-maf-first-alt}

- **gloci:** in an expression, `MAF` is bcftools' `min(f, 1-f)` with `f` the
  frequency of the first ALT among called alleles, and an unindexed `AC` is the
  first ALT's count (`AC[i]` picks another). `vcf filter --min-maf` keeps the
  second-most-common-allele definition (see
  [multi-allelic MAF](#maf-multiallelic)).
- **Old behaviour:** the two definitions agree on biallelic sites and differ
  on multi-allelic ones (counts 2,6,2 of 10: expression `MAF` 0.4, `--min-maf`
  0.2). `bcftools view -i 'MAF>=0.3'` keeps such a site; `--min-maf 0.3` drops
  it.
- **Reason:** `--include`/`--exclude` expressions must select what bcftools
  selects; the filter flags keep the documented definition.
- **Pinned by:**
  `expr-differs-from-bcftools-min-maf-is-not-expression-maf` (must differ) and
  `difference_maf_is_first_alt_not_second_most_common`.

## Filter expressions with `--samples` see the selected samples {#expr-samples-selection}

- **gloci:** with `--samples` / `--exclude-samples` (and their `-file` forms),
  `--include`/`--exclude` evaluate on the selected samples only, after masking
  (pipeline step 4): `MAF`, `AC`, `AN`, `N_MISSING`, `F_MISSING`, `N_PASS`, any-sample
  `FMT/` comparisons and `GT` keywords all use them. Sample subscripts still
  number the whole header.
- **Old behaviour:** `bcftools view -s s1,s2 -i 'MAF<0.3'` (1.21) evaluates the
  expression on ALL samples and then drops the unselected ones from the output.
  gloci's result equals the two-step pipeline
  `bcftools view -s s1,s2 | bcftools view -i 'MAF<0.3'`.
- **Reason:** the filters of one gloci command all describe the output cohort,
  as `--min-maf` already does with `--samples`.
- **Pinned by:** `expr-differs-from-bcftools-samples-selection` (must differ;
  its reference is the single-command bcftools form) and `maf_sees_the_selected_samples`
  in `crates/gloci/tests/vcf_filter_expr.rs`.
