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
- **`cat`:** `gloci cat` refuses a BAM input without the marker (exit 1, naming
  the file), where `samtools cat` warns (`EOF marker is absent`) and exits 0.
  This is stricter than `vcf concat --naive`, which warns and copies the blocks
  as they are: `cat` writes one new marker at the end of its output, and that
  would make a truncated input look complete. Pinned by
  `cat-no-eof-marker-vs-samtools` (exit status: must differ).

## `view --assume-sorted` {#view-assume-sorted}

- **gloci:** for an uncompressed BED file with no index, `view` scans the
  whole file and checks the order of every record it reads, so an unsorted file is an
  error. `--assume-sorted` binary-searches instead (like `bedextract`) and trusts
  the file to be sorted in the order given by `-g` or `--chrom-order`; it checks
  only the block it reads, and on an unsorted file the answer can be incomplete.
  With `--chrom-map` the trusted order is that of the *new* names: a map that
  renames onto a name already in the file can give a partial answer under
  `--assume-sorted` where the scan reports an error (see the note in
  [migration.md](migration.md#other-behaviour-changes)).
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
  "add INFO/END or INFO/SVLEN, use --vcf-span ref, or skip invalid records where the command supports it"
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

- **gloci:** `intersect --match alleles` (and `subtract --match alleles`, which
  drops an A record exactly when `intersect --mode none --match alleles` does,
  one shared ALT of a multi-allelic A dropping it whole) keeps a record when
  some ALT of it equals some ALT of a B record after trimming the shared
  prefix and suffix of REF and ALT, at the same trimmed position. `A>AG` at
  100 matches `C>GC` at 101 (an insertion of G), and records match when they
  share any one ALT of several. The overlap options (`-f`, `-s`, `-w`, ...) are
  refused with it.
- **Old behaviour:** `bcftools isec -c none` pairs records at the same
  position with identical REF and ALT lists, and does not normalise. So `isec
  -C` keeps a multi-allelic `C>A,G` that shares only `G` with a B `C>G`, where
  `subtract --match alleles` drops it.
- **Reason:** the same variant is often written with different anchors, and a
  multi-allelic record shares a variant with another that lists more or fewer
  ALTs. Run `bcftools norm` on both inputs and the results then agree for
  records with identical ALT lists.
- **Pinned by:** `intersect-alleles-unnormalised-vs-bcftools`,
  `subtract-alleles-unnormalised-vs-bcftools` and
  `subtract-alleles-shared-alt-vs-bcftools` (`bcftools isec -C`) (must differ);
  `intersect-alleles-vs-bcftools`, `intersect-alleles-multiallelic-vs-bcftools`,
  `subtract-alleles-vs-bcftools` and `subtract-alleles-multiallelic-vs-bcftools`
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

## `cluster -s` keeps input order and unknown-strand records {#cluster-strand}

- **gloci:** `gloci cluster -s` writes the records in input order, each with
  the id of its strand's cluster. Ids still start at 1 and increase through the
  file, and a cluster takes the next id when its first record is read, so a
  `+` and a `-` cluster that interleave get ids in file order. Records of
  unknown strand (`.`, no strand column, VCF records) form a third group, as in
  `merge -s`, and a warning says so.
- **Old behaviour:** `bedtools cluster -s` regroups the records by strand
  within each chromosome (all `+` records, then all `-` ones), numbers the
  clusters in that order, and writes nothing for records whose column 6 is `.`.
- **Reason:** a cluster id column should not reorder the file, and dropping
  records is a data-loss trap (as for `merge -s`, see
  [merge-unknown-strand](#merge-unknown-strand)). Without `-s`, and with `-s` on
  a file whose records all share one strand, gloci and bedtools agree. Sort by
  the id column to regroup (`sort -k7,7n`).
- **Pinned by:** `cluster-strand-vs-bedtools` and `cluster-strand-bed3-vs-bedtools`
  (a BED3 input, where bedtools writes nothing and exits 0) (must differ);
  `cluster-strand-single-strand-vs-bedtools` (must match). `cluster` without
  `-s`, with `-d 10`, with a negative `-d` and across chromosomes matches
  bedtools: `cluster-default-vs-bedtools`, `cluster-gap-vs-bedtools`,
  `cluster-negative-gap-vs-bedtools`, `cluster-global-ids-vs-bedtools`.

## `cluster` and `spacing` on VCF write a column-name line {#append-vcf-header}

- **gloci:** a VCF input is clustered (or measured) by its footprints, `--vcf-span
  allele` by default and `--vcf-span ref` as bedtools does. The output is TSV,
  as for every command that appends a column: a `#CHROM ... INFO cluster` (or
  `spacing`) line, then each VCF line with the value appended. The `##` meta
  lines are not copied, because the result is no longer a VCF.
- **Old behaviour:** `bedtools cluster` and `bedtools spacing` print the VCF
  records with the value appended and no header line at all.
- **Reason:** the same as `map` and `coverage` on VCF: the column the
  value sits in should be named.
- **Pinned by:** `cluster-vcf-vs-bedtools` and `spacing-vcf-vs-bedtools` (must
  differ by the header line only; the records are identical with
  `--vcf-span ref`). The goldens use `--vcf-span ref`: the default allele
  footprint can differ from bedtools, which uses `POS` to `POS + len(REF)`, so
  clusters and distances on indels may differ unless `--vcf-span ref` is given.

## `spacing` treats a zero-length record as a point {#spacing-zero-length}

- **gloci:** the distance to the previous record is the record's start minus the
  previous record's end (`0` when bookended), and `-1` when the record starts
  before the previous one ends. A zero-length record is a point on one boundary,
  not a base wide: a point at `5 5` after `0 5` is bookended (`0`), and a record
  starting at `5` after a point `5 5` is also `0`.
- **Old behaviour:** `bedtools spacing` widens a zero-length record by one base on
  each side. A point `10 10` after `0 5` is at distance 4 (gloci 5), a point
  `5 5` after `0 5` overlaps (`-1`, gloci `0`), and the distance between two points
  is two less than the gap.
- **Reason:** the footprint of a record is the same in every gloci command
  (a point does not cover neighbouring bases).
- **Pinned by:** `spacing-points-vs-bedtools` (must differ). `spacing` on
  records with a positive length (overlaps, nested records, bookended records,
  chromosome changes) matches: `spacing-vs-bedtools`. Both tools compare a record
  only with the one before it in the file, not with the one reaching furthest.

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
  sample form. bedmap's `--stdev` and `--variance` are sample statistics
  (n - 1) too, and print `NAN` for one value.
- **Reason:** the sample form is what a user summarising a handful of
  overlapping records expects, and it agrees with bedtools' `sstdev`. Use
  `stdev_pop` / `var_pop` to reproduce bedtools' `stdev` / `var`.
- **Pinned by:** `map-stdev-vs-bedtools` (must differ). `map-sstdev-vs-bedtools`,
  and `map-stdev-pop-vs-bedtools` match (bedtools map has no `var` op).

## `map` and `groupby` `distinct` keep first-seen order {#map-distinct-order}

- **gloci:** `distinct` lists each value once in the order the records were
  seen (in `map`: B file, then input order; in `groupby` and `merge`: input
  order). `distinct_sorted` lists them in bytewise order.
- **Old behaviour:** `bedtools map -o distinct` and `bedtools groupby -o
  distinct` list them in bytewise order (`10,3,5`); `bedmap --echo-map-id-uniq`
  does too.
- **Reason:** first-seen order follows the records, which makes `distinct`
  line up with `collapse`; the sorted form is kept under its own name.
- **Pinned by:** `map-distinct-vs-bedtools` and `groupby-distinct-vs-bedtools`
  (must differ). `map-distinct-sorted-vs-bedtools` matches.

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
  value and is skipped. `map -c` and `merge -c` do the same for a B or input
  line, and so does `N:count` with an explicit column (as `N:sum` does); column
  0 is the whole record and cannot be missing.
- **Old behaviour:** `bedtools groupby` prints a message that names the column
  but not the line, stops with the output so far, and exits 0. `bedtools map`
  (uniform B) and `bedtools merge -c` (uniform input) that lack the column on
  every line print the same kind of message and exit 0, for `count` as well as
  `sum`. A *ragged* file (lines of different widths) makes bedtools fail with
  exit 1, which agrees with gloci.
- **Reason:** a pipeline must not take truncated output for a result.
- **Pinned by:** `groupby-missing-column-vs-bedtools`,
  `groupby-absent-column-vs-bedtools`, `map-uniform-absent-column-vs-bedtools`
  and `merge-uniform-absent-column-vs-bedtools`, and their `-o count` twins
  `map-uniform-absent-column-count-vs-bedtools` and
  `merge-uniform-absent-column-count-vs-bedtools` (must differ, by exit status).
  `map-absent-column-vs-bedtools` (ragged B) matches: both tools fail.

## A value that is not a number is an error in `map` and `groupby` {#agg-non-numeric}

- **gloci:** a numeric op (`sum`, `mean`, `median`, ...) given a value that is
  not a number stops with an error naming the file, line and value.
- **Old behaviour:** `bedtools map` and `bedtools groupby` print a warning, skip
  the value and exit 0.
- **Reason:** skipping hides bad data, such as a header line read as data.
- **Pinned by:** `groupby-non-numeric-vs-bedtools` and
  `map-non-numeric-vs-bedtools` (must differ, by exit status).

## `min_record` and `max_record` write the record with `|` {#map-record-format}

- **gloci:** `-c COL:min_record` and `-c COL:max_record` (and `records`) write
  the B record's fields joined with `|`, each field as it is in the file.
- **Old behaviour:** `bedmap --min-element` and `--max-element` write the
  record with tabs between the fields and print the score column with 6
  decimals (`3.000000`).
- **Reason:** one output column per operation, so the result can be split on
  the column delimiter; the record's own text is not reformatted.
- **Pinned by:** `map-min-record-vs-bedops` (must differ).

## `map` prints the empty value, not `NAN`, for no values {#map-empty-value}

- **gloci:** an aggregation with no values prints `--empty` (default `.`);
  `count`, `count_distinct`, `bases`, `bases_uniq` and `indicator` print 0,
  because they always have a value. `frac_a` has none for an element with no
  mapped records and prints `--empty`.
- **Old behaviour:** `bedmap --sum` and `--mean` print `NAN` for an element
  with no mapped records. bedtools map prints `.` like gloci.
  `bedmap --bases-uniq-f` prints `0.000000` where gloci's `frac_a` prints `.`
  (while `bedmap --bases-uniq` and gloci's `bases_uniq` both print 0).
- **Reason:** one placeholder for every operation, chosen by the user; the
  bedtools default is the same.
- **Pinned by:** `map-empty-sum-vs-bedops` and `map-frac-a-empty-vs-bedops`
  (must differ). With every A
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

## `coverage --split` counts a spliced B once {#coverage-split-count}

- **gloci:** the hit count (column 7 of the default report) is the number of B
  features that overlap A. A spliced read, or a BED12 feature, that meets A through
  two of its blocks is one feature and counts once.
- **Old behaviour:** `bedtools coverage -split` counts the read once for each block that
  overlaps A, so a target spanning an intron reports 2. The covered bases and the
  fraction are the same as gloci's.
- **Reason:** the column is "number of B features overlapping A"; a count that grows
  with the number of blocks says nothing about how many reads there are. It applies
  to a BED12 B as much as to a BAM or SAM B.
- **Pinned by:** `coverage-bed-bam-split-count-vs-bedtools` (must differ).

## `-f`, `-F`, `-r` and `-e` under `--split` judge each hit on its own {#split-fraction}

- **gloci:** with `--split`, a hit counts when the bases A and B's blocks share, over the
  summed length of A's blocks (`-f`) or of B's blocks (`-F`), meet the thresholds. Each B
  record is judged alone, and a BED record without blocks is one block (D105).
  This holds for BED12 and for BAM/SAM reads (blocks cut at `N`).
- **Old behaviour:** `bedtools intersect -split -f` with a BAM or BED12 B does not decide
  hit by hit. The outcome for one A depends on how many other B records overlap it:
  three reads covering 90, 60 and 50 bases of a 100-base target give no hit at
  `-f 0.7`, although the first read alone covers 90%, while two of the same reads give
  both, though one covers only 60%. On a random set of spliced reads, `-f 0.5 -split`
  gave far fewer rows in bedtools than the per-hit rule (1,756 against 3,678 in one
  comparison). Without `-split`, bedtools judges each hit alone, as gloci does.
- **Reason:** a threshold that changes with the other reads cannot be reasoned about,
  and a single read gives the same answer in both tools (checked on reads and targets
  covering all of block, intron and bookended overlaps).
- **Pinned by:** `intersect-bed-bam-split-fraction-vs-bedtools` (must differ).

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

## `convert` reads PSL and RepeatMasker the way the files are written {#psl-rmsk-bedops-quirks}

- **gloci:** `gloci convert x.psl` skips the BLAT header (the five lines from
  `psLayout version 3` to the dashed line) and reads a final line that has no
  newline. RepeatMasker `.out` loses its three header lines the same way. A
  line that is not a PSL or RepeatMasker row (too few columns, a coordinate
  that is not a number, a strand that is not `+`/`-` or `C`) is a data error
  (exit 1) naming the file and the line; `--skip-invalid` skips it. That
  includes a PSL row whose blockCount differs from the number of entries in
  blockSizes, qStarts or tStarts, and (as `convert2bed` does) a row with more
  than 21 columns (exit 1; `convert2bed` exits 22).
- **Old behaviour:** `convert2bed --input=psl` (bedops 2.4.41) exits 22 with
  "Invalid field count (19)" on a header line, so a headered BLAT file must be
  stripped first, it silently drops a last line that has no newline, and it
  copies a blockCount that disagrees with the block lists unchecked.
- **Reason:** BLAT writes the header by default, and a file cut short or
  written without a final newline should not lose its last hit unseen.
- **Pinned by:** `convert-psl-header-unskipped-vs-convert2bed` (exit status;
  must differ) and `convert-psl-no-final-newline-vs-convert2bed` (must
  differ), and `convert-psl-block-count-vs-convert2bed` (exit status; must
  differ; both tools refuse a 22nd column, with exit 1 and 22:
  `convert-psl-22-columns-vs-convert2bed`, must differ).
  Every column otherwise matches: `convert-psl-vs-convert2bed`,
  `convert-psl-header-vs-convert2bed` (the header cut off for the reference)
  and `convert-rmsk-vs-convert2bed`.

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

## `genomecov` prints chromosomes in genome order {#genomecov-chrom-order}

- **gloci:** every chromosome of `-g` is written in genome order (the order the
  sorted input already has), a chromosome without features in its place with
  one zero-depth row or run.
- **Old behaviour:** `bedtools genomecov` writes the chromosomes that have
  features first, in input order, and appends the featureless ones at the end.
  The `-bga` and `-d` outputs are then not in genome order.
- **Reason:** a `--bga` bedGraph must stay sorted for bedGraphToBigWig and
  `gloci` tools that check order; the input is sorted, so streaming in genome
  order needs no buffering. When the featureless chromosomes are last in `-g`
  the outputs are identical.
- **BAM input:** `bedtools genomecov -ibam` also writes a header contig with no reads after
  the others; gloci keeps the header order.
- **Pinned by:** `genomecov-empty-chrom-middle-hist-vs-bedtools` and
  `genomecov-empty-chrom-middle-bga-vs-bedtools` (must differ);
  `genomecov-ibam-empty-mid-hist-vs-bedtools`, `genomecov-ibam-empty-mid-bga-vs-bedtools` and
  `genomecov-ibam-empty-mid-d-vs-bedtools` (must differ; `midempty.bam` has no reads on the
  middle contig of its header; `tags.bam` and `order.bam` have the empty contig last and
  match); the other
  `genomecov-*` cases use a genome that lists the empty chromosome last and
  match.


## `genomecov` of alignments counts deleted bases {#bam-genomecov-deletions}

- **gloci:** a read counts its whole reference span, so the bases under a `D` operation are
  covered (D105, the same rule as every command). `--split` counts the blocks between `N`
  operations, and `--split-deletions` counts the blocks cut at `N` and at `D`.
- **Old behaviour:** `bedtools genomecov -ibam` cuts every read at `D` (and not at `N`)
  unless `-ignoreD` is given, so deleted bases are not covered.
- **Reason:** one block rule for `intersect`, `coverage`, `convert` and `genomecov`, so the
  commands agree with each other. For parity use `--split-deletions`, which is bedtools
  `-split` (cuts at `N` and `D`). Plain bedtools `-ibam` (cut at `D`, not at `N`) has no gloci
  counterpart; `-ignoreD` is gloci's default, and `-split -ignoreD` is `--split`.
- **Pinned by:** `genomecov-ibam-default-vs-bedtools` (must differ; `cigars-bases.bam` has a
  read with `2D`). `genomecov-ibam-ignored-vs-bedtools`,
  `genomecov-ibam-bga-split-deletions-vs-bedtools` and
  `genomecov-ibam-split-ignored-vs-bedtools` match.

## `genomecov -g` with a BAM is used, not ignored {#genomecov-bam-genome-arg}

- **gloci:** `-g` gives the chromosomes and lengths even for a BAM or SAM input, in place of
  the header's. A read past the end of its `-g` chromosome is an error naming the read's
  place, as for BED input.
- **Old behaviour:** `bedtools genomecov -ibam x.bam -g G` ignores `-g` ("Genome (-g) files
  are ignored when BAM input is provided") and uses the header.
- **Reason:** a silently ignored option hides a mismatch between the genome and the reads;
  a different genome is a real request (a subset of the chromosomes, or the same build
  under other names).
- **Pinned by:** `genomecov-ibam-genome-arg` (must differ; the `-g` file gives a longer
  `chr1` than the header's). Without `-g` the header is used, as bedtools does.

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

## `genomecov` ignores points {#genomecov-points}

- **gloci:** a zero-length feature covers no bases and is ignored (D68).
- **Old behaviour:** `bedtools genomecov` widens a point at `p` to
  `[p-1, p+1)`, as `multiinter` does (`chr1 5 5` gives `chr1 4 6 1`). A file
  holding only `chr1 0 0` makes it abort with `corrupted size vs. prev_size`.
- **Reason:** inventing bases changes depths.
- **Pinned by:** `genomecov-points-vs-bedtools` (must differ).

## `genomecov` strand handling {#genomecov-unknown-strand}

- **gloci:** with `--5p`/`--3p` a feature whose strand is `.` or missing counts
  as plus. With `--strand + | -` a feature with unknown strand is not counted,
  and a BED record with fewer than 6 columns is an error naming the line.
- **Old behaviour:** `bedtools genomecov -strand` stops with "Invalid strand
  value (.)" on a `.` strand, and `-5`/`-3` use the end-1 base for an unknown
  strand for both ends.
- **Reason:** one rule for unknown strands across gloci (as `slop -s`); a
  feature on an unknown strand is neither + nor -.
- **Pinned by:** `genomecov-unknown-strand-5p-vs-bedtools` and
  `genomecov-unknown-strand-filter-vs-bedtools` (must differ).

## `genomecov` rejects features past the chromosome end {#genomecov-past-end}

- **gloci:** a feature that extends past the end of its chromosome in `-g`
  exits 1 naming the file, line and chromosome. A feature on a chromosome
  missing from `-g` is also an error.
- **Old behaviour:** `bedtools genomecov` silently clips the feature at the
  chromosome end (a chromosome missing from `-g` is an error).
- **Reason:** a mismatched genome build gives wrong depths without a sign.
- **Pinned by:** `genomecov-past-end-vs-bedtools` (exit status must differ).

## `genomecov --max-depth` is histogram-only {#genomecov-max-depth}

- **gloci:** `--max-depth N` with `--bedgraph`, `--bga` or `--per-base` is a
  usage error.
- **Old behaviour:** `bedtools genomecov -max N` is documented as irrelevant
  for `-d` and `-bg` and is silently ignored there.
- **Reason:** a silently ignored cap suggests the depths were capped.
- **Pinned by:** `genomecov-max-depth-bedgraph-vs-bedtools` (exit status must
  differ).

## `vcf extract --info`: an ID column, `.` for missing and `0` for an unset flag {#vcf-extract-info}

- **gloci:** `gloci vcf extract --info KEYS` writes `CHROM POS ID REF ALT KEY...`.
  A key a record lacks is `.`; an INFO flag (header `Type=Flag`) is `1` when
  set and `0` when not.
- **Old behaviour:** `vcftools --get-INFO KEY` writes `CHROM POS REF ALT KEY...`
  (no ID column), writes `?` for a key a record lacks, and `?` for an unset flag.
- **Reason:** decision D73 keeps ID so a row can be joined on it; `.` is the
  VCF missing value and `?` is not.
- **Pinned by:** `vcf-extract-info-vs-vcftools` (must differ) and the
  `info_table_*` tests in `crates/gloci/tests/vcf_extract.rs`.

## `vcf extract --format --layout wide`: `SAMPLE:KEY` columns, every site {#vcf-extract-format-wide}

- **gloci:** `--format KEYS --layout wide` writes `CHROM POS ID REF ALT` and one
  `<sample>:<KEY>` column per sample and key, for every site that passes the
  selection. A sample without the key (or a site whose FORMAT lacks it) reads
  `.`. Any FORMAT key works, as text.
- **Old behaviour:** `vcftools --extract-FORMAT-info KEY` writes `CHROM POS`
  and a column named by the sample alone, and leaves out a site whose FORMAT
  has no such key. It takes one key.
- **Reason:** decision D73 names columns `<sample>:<KEY>` so several keys can
  share a row, and writes one row per site so rows line up across runs.
- **Pinned by:** `vcf-extract-format-wide-vs-vcftools` (must differ) and the
  `format_*` tests in `crates/gloci/tests/vcf_extract.rs`.

## `vcf extract --genotype-matrix 012`: counts non-reference alleles {#vcf-extract-012}

- **gloci:** each cell is the number of non-reference alleles among the alleles
  that are called: `0/0` and a haploid `0` are 0, `0/1` and a haploid `1` are 1,
  `1/1` is 2, `0/.` is 0, `./1` is 1, and a genotype with no called allele is
  `-1`. The files, their layout (`.012` has a row per sample led by its 0-based
  index; `.012.indv`; `.012.pos`) and the omission of multi-allelic sites
  (more than one ALT) match vcftools. A genotype of ploidy above 2, or one that
  names an allele the ALT column lacks, is an error with the record's position
  (with `--skip-invalid` the latter is written as `-1` and counted in a warning).
- **Old behaviour:** `vcftools --012` writes 1 for a haploid reference call
  (`0`) and for `0/.` or `./0`, writes 1 for `2/2` on a one-ALT site, and stops
  at a ploidy above 2.
- **Reason:** decision D73 defines the cell as the count of non-reference
  alleles; vcftools counts a haploid or half-missing reference call as 1, which
  misreports male X and mitochondrial calls.
- **Pinned by:** `vcf-extract-012-haploid-vs-vcftools` (must differ) and
  `vcf-extract-012-vs-vcftools` (must match, diploid input with a multi-allelic
  site), plus `allele_beyond_alt_errors_or_is_skipped_and_vcftools_writes_1` and `genotype_matrix_matches_vcftools_files` in
  `crates/gloci/tests/vcf_extract.rs`.

## `vcf stats --singletons` names a carrier that masking left {#vcf-stats-singletons-masked}

- **gloci:** with genotype masking (`--min-gq`, `--min-dp`, `--max-dp`,
  `--mask-filtered-genotypes`), a singleton (`S`) or doubleton (`D`) allele is
  counted over the genotypes masking left, and its carrier is searched among
  those genotypes only. A doubleton needs an unmasked homozygote.
  `PREFIX.singletons` agrees with vcftools on every site where no masked
  genotype carries the allele; every other statistic agrees everywhere.
- **Old behaviour:** `vcftools --singletons --minGQ 20` counts the alleles of the
  unmasked genotypes but names the first sample whose raw genotype holds the
  allele, masked or not, and decides `D` from a raw homozygote. In
  `masked.vcf` at `1:30` the surviving `G` is `s2`'s, and vcftools prints
  `S G s1` (`s1`'s `1:4:10` is masked). In `ftmask.vcf` with
  `--remove-filtered-geno-all`, `1:50` has `D G s1` in vcftools and
  `D G s2` in gloci (`s1`'s `1/1` is filtered, `s2`'s `1/1` is not), and `1:60`
  has the row `D A s4` in vcftools (two `A` in unmasked heterozygotes, and a
  masked `0/0` of `s4` taken as the homozygote) where gloci writes none.
- **Reason:** naming a sample whose genotype was removed contradicts the
  counts on the same line; the masked genotype is missing data (D50).
- **Pinned by:** `vcf-stats-singletons-masked-vs-vcftools` and
  `vcf-stats-singletons-doubleton-masked-vs-vcftools` (must differ) and
  `singletons_name_an_unmasked_carrier_unlike_vcftools` in
  `crates/gloci/tests/vcf_stats.rs`, which holds the exact substitutions
  against the committed vcftools files.

## `vcf stats --snp-density` needs sorted input {#vcf-stats-snp-density-unsorted}

- **gloci:** `--snp-density BIN` streams its rows, so it needs the input sorted
  by chromosome and position. A chromosome that comes back after another, or a
  position in an earlier bin, stops the run with the record's position and
  `run: gloci sort FILE`; no file of the run is written.
- **Old behaviour:** `vcftools --SNPdensity` accepts unsorted input and writes
  the bins as it meets them, repeating a chromosome's rows when it returns.
- **Reason:** rows already written cannot be merged, and holding every bin of a
  multi-gigabyte file in memory to repair unsorted input is what the streaming
  design avoids. A silent duplicate row would be a wrong density.
- **Pinned by:** `vcf-stats-snp-density-unsorted-vs-vcftools` (must differ) and
  `unsorted_input_stops_snp_density_with_a_hint_and_publishes_nothing` in
  `crates/gloci/tests/vcf_stats.rs`.

## Masking never applies to a missing DP or GQ {#vcf-stats-mask-missing}

- **gloci:** `--min-gq`, `--min-dp` and `--max-dp` mask a genotype only when it
  has a value below the minimum (or above the maximum). A genotype whose DP or
  GQ is `.` is kept (D50). This holds in every command that masks, `vcf filter`
  and `vcf extract` as well as `vcf stats`.
- **Old behaviour:** vcftools reads a `.` as -1 when the FORMAT lists the field,
  so `--minGQ` and `--minDP` filter it (`GQ` or `DP` `.` is counted in
  `N_GENOTYPES_FILTERED` of `--missing-indv` and drops out of every statistic).
  When the FORMAT does not list the field, vcftools keeps the genotype, as
  gloci does. `--maxDP` does not filter a `.` in either tool.
  The same holds for a sample cell cut short before DP or GQ (`0/1:9` under
  `GT:DP:GQ`, or a bare `0/0`): vcftools masks it, gloci keeps it, because the
  value is missing just as a `.` is.
- **Reason:** a missing value is not evidence the genotype is bad. D50 keeps a
  call whose quality or depth was never recorded.
- **Pinned by:** `vcf-stats-mask-missing-vs-vcftools` (must differ),
  `truncated_sample_cells_are_not_masked_unlike_vcftools` (the cut-short cells,
  from the committed vcftools files `mqt.*`) and
  `missing_dp_gq_is_not_masked_unlike_vcftools` in
  `crates/gloci/tests/vcf_stats.rs` (both tools' counts, from the committed
  vcftools files `mq.*`).

## `vcf stats --json`: the summary.json schema

`--json` writes `PREFIX.summary.json`; `--summary` writes the same figures as
text in `PREFIX.summary.txt`. These are not vcftools files: vcftools has no
overview, and gloci's layout is its own. Keys are sorted at every level, the
file ends with a newline, and `schema_version` is bumped when a key is renamed
or removed (adding a key does not bump it).

| Key | Type | Meaning |
| --- | --- | --- |
| `counts_by_type` | object | `indel`, `mnp`, `other`, `ref`, `snp`: records of each type (integers) |
| `depth.mean` | number or null | mean FORMAT DP over the unmasked genotypes that have a DP; null when there are none |
| `depth.values` | integer | how many genotypes that mean is over |
| `genotypes.total` | integer | genotypes of the selected records and samples |
| `genotypes.missing` | integer | missing genotypes; masked ones count here too |
| `genotypes.filtered` | integer | the masked ones among them (`--min-dp`, `--max-dp`, `--min-gq`, `--remove-filtered-geno-all`) |
| `genotypes.missing_fraction` | number or null | `missing / total`; null when `total` is 0 |
| `multiallelic_sites` | integer | records whose ALT column lists more than one allele |
| `records` | integer | selected records |
| `samples` | integer | selected samples |
| `schema_version` | integer | 1 |
| `ts_tv.transitions`, `ts_tv.transversions` | integer | biallelic single-base SNPs of each kind |
| `ts_tv.ratio` | number or null | `transitions / transversions`; null when there are no transversions |

Numbers that are not integers are written with six significant digits, as
vcftools writes them (`%g`).

The record type is that of the ALT alleles when they agree: `snp`, `mnp` or
`indel` (insertions and deletions are both `indel`), so a multi-allelic record
whose ALTs are all SNPs is `snp`. Alleles of different types, such as a SNP
with an indel, make the record `other`, and `other` also takes structural
variants (there is no `sv` count; a symbolic ALT such as `<DEL>` is `other`). A
`<*>` ALT is ignored beside a real one, so `T,<*>` is `snp`; a record whose only
ALT is `<*>` or `.` is `ref`.

## `vcf filter ... | vcf stats` is `vcf stats` without masking

Running `vcf filter` with the selection flags (`--samples`, regions, site
filters) and piping it to `vcf stats` gives the same files as `vcf stats` with
those flags. The masking flags (`--min-dp`, `--max-dp`, `--min-gq`,
`--remove-filtered-geno-all`) are the exception: `vcf filter` writes a masked
genotype as `./.` but keeps the record's original DP and GQ, so a later
`vcf stats` sees a missing genotype with a depth. Allele and genotype counts
agree, but depth statistics still include the masked genotype's DP and
missingness counts it as missing rather than as filtered. To get the masked
statistics, give the masking flags to `vcf stats` itself.

## `vcf stats` and the vcftools behaviour it reproduces

These are not differences: gloci writes what vcftools writes. They are listed
because a user may take them for bugs. Each is pinned by the byte-for-byte
comparison of `crates/gloci/tests/vcf_stats.rs` with the committed vcftools
0.1.16 files in `crates/gloci-vcf/tests/golden/stats`.

- `--snp-density` counts one SNP per position (two records at one position count
  once), counts every record (indels, multi-allelic) in its bin, starts bins at
  multiples of BIN on the 1-based position, and writes zero rows from a
  chromosome's first bin to its last.
- `--hardy` writes no row for a site that is not biallelic or that has any
  one-allele genotype (`1`, `.`). A half-missing `./1` is left out of the genotype
  classes but its allele feeds the expected counts. An all-missing site still
  gets a row.
- `--het` uses only sites with exactly two alleles, both seen among the called
  alleles, and no one-allele genotype among the unmasked samples; a sample adds a
  site only with a fully called genotype. `E(HOM)` is the sum over those sites of
  `1 - 2pq * n / (n - 1)` (p the REF frequency, n the called alleles) printed
  with one decimal, and `F = (O - E) / (N - E)` with five. A sample with no such
  site has no row. A multi-allelic site is skipped, not split.
- `--missing-indv` counts a genotype as missing when its first allele is `.`
  (`./1` and a lone `.` are, `1/.` is not) or the sample has no GT;
  `N_DATA` counts the unmasked genotypes of every site and `F_MISS` is `nan`
  when there are none. `--missing-site` counts alleles instead, as vcftools does,
  so the two `F_MISS` columns do not describe the same thing.
- `--depth` averages the DP of the genotypes that have one, `nan` for none.
- `--FILTER-summary` sorts the tags inside a FILTER, orders rows by count
  descending then name descending, and takes Ts/Tv from the first ALT;
  `--TsTv-summary` uses biallelic single-base SNPs only.
- `--indel-hist` has the header `LENGTH COUNT PRCT`; a SNP-only site counts once
  at length 0.
- A genotype of ploidy above 2 stops `--freq`, `--counts`, `--missing-site`,
  `--missing-indv`, `--hardy`, `--het` and `--singletons` with the record's
  position, as vcftools does; depth, quality, FILTER, Ts/Tv, indel and the
  summary read through such records.
- Numbers use C `%g` (6 digits), `%e` for the HWE columns, `%.2f`, `%.1f` and
  `%.5f` where vcftools does, and `nan` in lower case.

`vcf stats` also differs from vcftools in ways that are not about its output:
it takes every `vcf filter` selection flag, runs any number of statistics in one
pass (vcftools allows one output function per run), writes each file only when
the whole run succeeds, and defaults `--vcf-span` to `ref` because no statistic
depends on the span.

## `index` on FASTA refuses what `samtools faidx` quietly mishandles {#fasta-index-strictness}

- **gloci:** `gloci index ref.fa` exits 1, naming the line, for three inputs:
  a sequence name used twice, a header with no name (`>` alone), and a sequence
  line holding a space, a control character or a non-ASCII byte. No index is
  written.
- **Old behaviour:** `samtools faidx` exits 0. It keeps the first of two equal
  names and warns, indexes an empty name, and does not count a space as a base
  (so the index offsets no longer match the file).
- **Reason:** each of these gives an index that returns the wrong bases or hides
  a sequence; an error that names the line is safer. Every FASTA both tools
  accept gets a byte-identical `.fai` (and `.gzi` for BGZF). Empty sequences
  are left out of the index as samtools does, with a warning.
- **Also:** `*` and `-` are accepted, because the check is "ASCII graphic". A
  trailing space or tab on a sequence line is the rejection users will most
  often meet: samtools indexes it wrongly and returns truncated sequence.
- **Stale indexes:** detection at fetch time checks only the bytes a fetch
  reads. A same-length edit, or a sequence appended after the last entry, is
  not detected; this is the same trust level as samtools/htslib. The "index
  older than FASTA" mtime warning covers most real edits.
- **Pinned by:** `index-fasta-duplicate-name-vs-samtools`,
  `index-fasta-space-in-sequence-vs-samtools` and
  `index-fasta-empty-name-vs-samtools` (exit status: must differ);
  `index-fasta-irregular-lines-vs-samtools` pins the agreeing side (both exit
  1), and `index-interop-faidx` and `index-interop-faidx-bgzf` pin identical
  indexes.

## `index` of a BAM with a long reference writes a `.csi` {#bam-index-long-reference}

- **gloci:** `gloci index x.bam` writes `x.bam.bai`, or `x.bam.csi` when a reference is
  longer than 2^29 (536,870,912) bases or a read ends past that. `--csi` asks for a `.csi`
  at any length. The `.csi` has min_shift 14 and the smallest depth (at least 5) that
  holds the longest reference. A `.csi` or `.bai` left from before is removed, so one
  index is left.
- **Old behaviour:** `samtools index x.bam` exits 1 with `Region ... cannot be stored in a
  bai index. Try using a csi index`, and `samtools index -c` is the way to a `.csi`.
- **Reason:** the format follows from the data, as `gloci index` of a BED, GFF or VCF
  switches from `.tbi` to `.csi` on its own, so one command line works for every genome.
  samtools reads the `.csi` gloci writes, and gloci reads samtools' `.bai` and `.csi`;
  queries through either return the same reads. The two tools write the same content
  (bins, chunks, linear index, metadata) for a file that fits a `.bai`; the order of the
  bins within a reference may differ (htslib writes them in hash order, gloci sorted), so
  the bytes can differ.
- **Pinned by:** `index-bam-long-reference-vs-samtools` (exit status: must differ);
  `index-interop-bam-csi-samtools` pins the agreeing side (samtools reads gloci's
  `.csi`), and `index-interop-bam-samtools-view` the `.bai`.

## `view` of a BAM region without an index scans the file {#bam-region-without-index}

- **gloci:** `gloci view x.bam chr1:100-200` reads through `x.bam.bai`, `x.bai` or
  `x.bam.csi` when there is one. With none, or with SAM or a BAM on standard input, it
  reads the whole file and keeps the reads that overlap a region, with a warning that
  says to run `gloci index x.bam`. A scan and an index give the same reads.
- **Old behaviour:** `samtools view x.bam chr1:100-200` exits 1 without an index
  (`Random alignment retrieval only works for indexed SAM.gz, BAM or CRAM files`).
- **Reason:** a slow answer with a warning serves a pipeline better than a failure; the
  warning names the fix. An index that is older than its BAM is used with a warning too.
- **Pinned by:** `view-bam-region-no-index-vs-samtools` (exit status: must differ);
  `view-bam-region-vs-samtools` pins the indexed side.

## `getfasta` and `nuc` fail on an interval outside the FASTA {#getfasta-outside-fasta}

- **gloci:** `gloci getfasta` and `gloci nuc` exit 1 at the first interval whose chromosome is
  not in the FASTA's `.fai`, or whose end is past the chromosome's length. The
  message names the input line, the interval and the length, and points to
  `--skip-invalid`. Nothing more is written; a file given with `-o` is not
  created. With `--skip-invalid` the interval is skipped and one warning gives
  the count and the first cause; the exit status is 0.
- **Old behaviour:** `bedtools getfasta` prints `chromosome (X) was not found in
  the FASTA file. Skipping.` or `Feature (...) beyond the length of ... size.
  Skipping.` to stderr, skips the record and exits 0; `bedtools nuc` prints
  `Feature (...) beyond the length of ... size (N bp).  Skipping.` and does the
  same (a chromosome the FASTA lacks is the same message with a size of 0). A script that checks only
  the exit status gets a short output and no failure.
- **Reason:** a chromosome-naming mismatch (`chr1` against `1`) or a FASTA of
  another assembly makes every interval fall outside; with the exit status 0 the
  output is empty or partial and the mistake goes unnoticed. `--chrom-map`
  fixes naming, and `--skip-invalid` gives bedtools' behaviour.
- **Unchanged:** a zero-length interval is skipped with one warning that counts
  them, as in bedtools (it has no sequence); it is not an error. `nuc` writes
  its header line for the first record even when that record is skipped, as
  bedtools does.
- **Pinned by:** `getfasta-outside-vs-bedtools` and
  `getfasta-missing-chrom-vs-bedtools` (exit status: must differ);
  `getfasta-outside-skip-vs-bedtools` and `getfasta-missing-chrom-skip-vs-bedtools`
  pin the matching output under `--skip-invalid`. The same for `nuc`:
  `nuc-outside-vs-bedtools`, `nuc-missing-chrom-vs-bedtools`,
  `nuc-outside-skip-vs-bedtools` and `nuc-missing-chrom-skip-vs-bedtools`.

## `getfasta --split` refuses BED12 records with inconsistent blocks {#getfasta-split-blocks}

- **gloci:** with `--split`, a BED12 record whose blocks are not self-consistent
  exits 1, naming the input line and the problem (`invalid BED12: ...`). A record
  is inconsistent when its blocks are out of order or overlap, its first block is
  not at offset 0, or its last block does not end at chromEnd. With
  `--skip-invalid` the record is skipped and one warning gives the count and the
  first cause (together with any intervals outside the FASTA).
- **Old behaviour:** `bedtools getfasta -split` does not check the blocks. It
  concatenates them in the order listed and writes the sequence, exit 0.
- **Reason:** gloci applies the same BED12 validation as `convert --blocks`; a
  record whose blocks disagree with its span has no well-defined exons, and the
  sequence bedtools writes for it is not the transcript the record describes.
- **Workaround:** fix the BED12 record, or pass `--skip-invalid` to skip it.
- **Pinned by:** `getfasta-split-bad-blocks-vs-bedtools` (exit status: must
  differ) and `getfasta-split-bad-blocks-skip-vs-bedtools`, which pins the
  matching output for the consistent records under `--skip-invalid`.

## `maskfasta` keeps the FASTA layout; bedtools re-wraps and cannot match CRLF names {#maskfasta-layout}

- **gloci:** `gloci maskfasta` writes every line break where the input had it:
  uneven line widths, blank lines, CRLF line ends and a missing final line break
  are all kept, and only the masked bases change. A CRLF FASTA is masked like
  any other.
- **Old behaviour:** `bedtools maskfasta` re-wraps each sequence to the width of
  the first line, so an irregular FASTA comes out with different line breaks.
  On a CRLF FASTA it keeps the CR in the header name (`chrM\r`), so no mask
  finds its sequence and nothing is masked, with exit 0.
- **Reason:** the output is the input with fewer bases, so `diff` shows only the
  masked bases, and the `.fai` of the output matches the `.fai` of the input.
  Silently doing nothing on CRLF input is a wrong answer.
- **Also:** bedtools ends the output with a line break even when the input's
  last line had none, and drops sequences that have no bases (for
  `>e`, `>chr1` with bases, `>f` it prints only `>chr1`). gloci keeps the last
  line as it was and keeps `>e` and `>f`.
- **Pinned by:** `maskfasta-irregular-lines-vs-bedtools`,
  `maskfasta-crlf-vs-bedtools`, `maskfasta-no-final-newline-vs-bedtools` and
  `maskfasta-empty-sequences-vs-bedtools` (outputs must differ).

## `maskfasta` keeps the whole header; bedtools drops the description {#maskfasta-header}

- **gloci:** the header line is copied as it is (`>chr1 Homo sapiens chromosome 1`).
  A mask names a sequence by the first word of its header.
- **Old behaviour:** `bedtools maskfasta` writes only the first word (`>chr1`)
  unless `-fullHeader` is given.
- **Reason:** dropping the description loses information the user did not ask to
  remove. `-fullHeader` output is identical in both tools
  (`maskfasta-full-header-vs-bedtools`).
- **Pinned by:** `maskfasta-description-vs-bedtools` (output must differ).

## `maskfasta` refuses a mask past the end of its sequence {#maskfasta-past-end}

- **gloci:** a mask whose end is beyond the length of its sequence is an error
  (exit 1) that names the sequence, the mask and the length. With `-o`, nothing
  is published. With `--skip-invalid` the mask is clipped to the end of the
  sequence instead (one that starts at or past the end masks nothing), and one
  warning gives the number of such masks and the first.
- **Old behaviour:** `bedtools maskfasta` accepts a mask that starts inside the
  sequence and ends past it (exit 0) and writes the whole mask length, so the
  sequence grows by the overhang (a mask 990-1010 on a 1,000 base chr1 leaves
  chr1 with 1,010 bases); it aborts with a C++ exception (exit 134) when the
  mask starts past the end.
- **Reason:** such a mask nearly always means the intervals belong to another
  assembly. Masking a different sequence silently gives a wrong reference.
- **Pinned by:** `maskfasta-past-end-vs-bedtools` (exit status: must differ);
  `maskfasta-past-end-skip-invalid-vs-bedtools` pins gloci's clipped output
  against bedtools' lengthened sequence (output must differ).

## `maskfasta --mask-char` refuses a space {#maskfasta-mask-char}

- **gloci:** `--mask-char` takes one printable, non-space ASCII character; a
  space exits 2 with the reason.
- **Old behaviour:** `bedtools maskfasta -mc ' '` writes spaces into the
  sequence, exit 0.
- **Reason:** a space inside a sequence line makes the FASTA impossible to
  index (`gloci index` and `samtools faidx` reject it), so the output would
  be unusable for most tools.
- **Pinned by:** `maskfasta-mask-char-space-vs-bedtools` (exit status: must
  differ).

## SAM input: a read whose reference is not in the header is an error {#sam-unknown-rname}

- **gloci:** a record whose RNAME or RNEXT is not an `@SQ` name of the header
  is a data error naming the file and line. The message says to add the `@SQ`
  line, fix the name, or pass `--skip-invalid`, which skips such records and
  counts them. A header with no `@SQ` lines refuses every placed read the same way.
- **Old behaviour:** `samtools view` warns `unrecognized reference name` (or
  `unrecognized mate reference name`), writes the read as unmapped with its
  position kept (RNAME `*`, flag |= 4, or RNEXT `*`), and exits 0. With no `@SQ`
  lines at all it fails.
- **Reason:** turning reads into unmapped ones quietly loses their placement, for
  a typo or a header from another build.
- **Pinned by:** `view-sam-unknown-rname-vs-samtools` and
  `view-sam-unknown-rnext-vs-samtools` (exit status must differ).

## SAM input: junk that samtools lets through is refused {#sam-parse-strictness}

- **gloci:** these are data errors naming the file and line, honoured by
  `--skip-invalid`:
  - a character (`A`) optional field of other than one character (`XA:A:ab`);
  - a `B` array element out of range for its type (`B:c,200`) or empty (`B:c,1,`);
  - a float that is not plain decimal text (`XF:f: 1.5`, `XF:f:1.5x`, hex floats
    such as `0x1p3`);
  - a QUAL character outside `!` to `~` (0x7f);
  - a CIGAR length of more than 9 digits (so leading zeros as in `0000000004M`).

  The checks run when the line is read, whatever the output format.
- **Old behaviour:** `samtools view` reads each of these, storing the bytes it
  parsed (`ab` as `a`, the out-of-range element wrapped, the trailing comma as a
  zero element, the number by `strtod`, QUAL 0x7f as 94).
- **Also:** a SAM line that gloci does not normalise is copied as written (D120), spelling
  included: lowercase bases, `0x10` for a FLAG, `I` and `c` tag types, `0.10` and `+7`.
  samtools rewrites all of these, because it parses every line into a BAM record and prints
  that. A line that htslib reads as a different record (a mapped read with CIGAR `*`, with no
  position or no reference, or whose RNEXT has PNEXT 0) is written as the record samtools
  writes, flag 4 and `*` included, so the output agrees with BAM output and with the filters.
  SAM to BAM encodes what samtools encodes.
- **A `@PG` line with no ID:** samtools cannot chain its own `@PG` line to it and fails the
  run; gloci adds its line and goes on.
- **Reason:** each of these is more likely a broken file than a meant value.
- **Pinned by:** `view-sam-strict-char-two-vs-samtools`,
  `view-sam-strict-b-range-vs-samtools`, `view-sam-strict-b-comma-vs-samtools`,
  `view-sam-strict-float-space-vs-samtools`,
  `view-sam-strict-float-trailing-vs-samtools`,
  `view-sam-strict-qual-del-vs-samtools` and
  `view-sam-strict-cigar-zeros-vs-samtools` and `view-sam-pg-without-id-vs-samtools` (exit
  status must differ).

## `cat` refuses BAM files whose `@SQ` lists differ {#alignment-cat-sq}

- **gloci:** `gloci cat` compares the reference lists (`@SQ` names and lengths, in
  order) of all its inputs with the first one's. If they differ it exits 1,
  names the file and the first difference, and writes nothing.
- **Old behaviour:** `samtools cat` does not check. It writes the first file's
  header and copies the records of every file as they are (exit 0).
- **Reason:** a BAM record names its reference by number, so records copied from
  a file with another reference list point at the wrong contigs, and the joined
  file looks fine until a variant caller or a genome browser uses it.
  The other merge rules are samtools': the header is the first file's, plus the
  `@RG` lines of later files whose ID it does not have (appended at the end;
  an ID already there is skipped without comparing the rest of the line);
  `@PG`, `@CO` and `@SQ` lines of later files are dropped. The order of the
  records is not checked. samtools also keeps empty interior (non-EOF) BGZF
  blocks of the inputs; gloci drops them, and the records are identical.
- **Pinned by:** `cat-different-sq-vs-samtools` (the lists differ in order) and
  `cat-different-sq-length-vs-samtools` (in a length) (exit status: must differ).


## `sort -m` refuses alignment files whose `@SQ` lists differ {#alignment-merge-sq}

- **gloci:** `gloci sort -m` compares the reference lists (`@SQ` names and lengths, in
  order) of all its alignment inputs with the first one's. If they differ it exits 1, names
  the file and the first difference, and writes nothing.
- **Old behaviour:** `samtools merge` appends the `@SQ` lines it does not know and
  only warns when that loses the coordinate order (exit 0); records of the later
  files are renumbered to the merged list.
- **Reason:** a merge of files aligned to different references, or to one
  reference in another order, is almost always a mistake, and the merged file
  looks fine until a variant caller or a genome browser uses it.
- **Pinned by:** `aln-merge-different-sq-vs-samtools` (a length differs; exit status:
  must differ).

## `sort -m` refuses conflicting `@RG` lines {#alignment-merge-rg}

- **gloci:** two inputs that have an `@RG` line with the same `ID` but different
  other fields are an error (exit 1) naming both lines. An identical line is kept
  once, and `@PG` lines with the same `ID` keep the first.
- **Old behaviour:** `samtools merge -c` keeps the first `@RG` line and drops the
  other without a word (exit 0); without `-c` it renames the colliding ID with a
  random suffix (`g1-58152F2E`), which gloci cannot reproduce.
- **Reason:** reads of the second file would be attributed to the first file's
  sample or library.
- **Pinned by:** `aln-merge-rg-conflict-vs-samtools` (exit status: must differ).


## `convert` of alignments: `--blocks` with a score tag or `--cigar` {#bamtobed-split-tag}

- **gloci:** `convert --blocks --score-tag NM` writes the BED6 rows of `--blocks` with
  the tag as the score, the six columns of every other `convert` row. `--cigar` with
  `--blocks`, `--bed12` or `--split-deletions` is a usage error (exit 2) that says to
  drop one of them.
- **Old behaviour:** `bedtools bamtobed -split -tag NM` (also `-splitD`) writes an extra
  column, the block start again, after the start: `chr1 600 600 605 c_skip 0 +`.
  `-splitD` also writes an empty trailing block for a CIGAR that ends in `D`, which
  gloci drops (D105: empty blocks are dropped). `-split -ed` is refused (exit 1, "Edit distances cannot be computed for each
  'chunk'"). `-cigar -split` drops the CIGAR without a word (exit 0).
- **Reason:** the extra column makes the file no longer BED, and a flag that is quietly
  ignored hides a mistake.
- **Pinned by:** `bamtobed-split-tag-vs-bedtools` (must differ),
  `bamtobed-cigar-split-vs-bedtools` (exit status: gloci 2, bedtools 0) and
  `bamtobed-split-ed-vs-bedtools` (exit status: gloci 0, bedtools 1).

## `convert` of alignments cuts reads at the end of their reference, and blocks at their bases {#bamtobed-past-reference-end}

- **gloci:** a read that runs past the `LN` of its reference (a chrM read wrapping
  round) is written up to `LN`, with one warning per input. A read that starts at or
  past `LN` covers no base and is skipped and counted at `-v`. Every interval command
  does the same (D101). The blocks of `--bed12` and `--blocks` are the ones every
  `--split` command compares (one shared rule): blocks with no base are left out. A
  BED12 row opens at the first kept block (chromStart and thickStart) and closes at
  the last (chromEnd and thickEnd), so its first block starts at offset 0 and its last
  ends at chromEnd, as BED12 requires.
- **Old behaviour:** `bedtools bamtobed` writes the interval as the CIGAR gives it,
  past the end of the reference (`chr1 994 1004` on a 1,000 base chr1), and writes a
  read that starts past the end as well. With `-bed12` a start beyond 2^31 gives a
  negative end. A spliced read `5M20N5M` at 991 on a 1,000 base chr1 is
  `990 995 ... 1 5 0` in gloci's BED12; BED6 keeps the cut span (`990 1000`). A read
  that starts in `N` (`5N10M` at 101), or in `D` under `-splitD`, gets a size-0 first
  block from bedtools (`100 115 ... 2 0,10 0,5`), and one that ends in `N` or `D` a
  size-0 last block (`10,0 0,15`); `-split` writes those blocks as zero-length rows
  (`100 100`). gloci writes `105 115 ... 1 10 0` and `400 410 ... 1 10 0`, and no
  zero-length rows. BED6 is the same in both (`100 115`).
- **Reason:** an interval past its chromosome fails in a genome browser and in
  `bedToBigBed`, and the other gloci commands refuse it. So would a BED12 row whose
  first block does not start at chromStart: `gloci convert x.bam --bed12 | gloci
  intersect -a - -b y.bed --split` must read what convert writes.
- **Pinned by:** `bamtobed-past-end-vs-bedtools`, `bamtobed-edge-vs-bedtools`,
  `bamtobed-past-end-bed12-vs-bedtools`, `bamtobed-ends-bed12-vs-bedtools`,
  `bamtobed-ends-bed12-splitd-vs-bedtools` and `bamtobed-ends-split-vs-bedtools` (must
  differ); `bamtobed-ends-vs-bedtools` (BED6, matches).

## `convert` and `genomecov` of a read with no reference bases: one base {#bamtobed-no-reference-ops}

- **gloci:** a mapped read whose CIGAR has no `M D N = X` operation (`10S5I`) is the one
  base at its position, `[pos, pos+1)`, in every layout: BED6, BED12 (one block of
  size 1) and `--blocks`. That is the span htslib's `bam_endpos` gives and every gloci
  interval command uses (D101), so `gloci convert x.bam | gloci intersect ...` and
  `gloci intersect x.bam ...` agree. `gloci genomecov` counts that base too.
- **Old behaviour:** `bedtools bamtobed` writes a zero-length interval (`chr1 800 800`),
  with `-bed12` a single block of size 0 (`1 0 0`), and `-split` the same zero-length row.
  `bedtools genomecov -ibam` counts no base for such a read.
- **Reason:** one span for a read in all commands. A zero-length row also overlaps
  nothing under the rules of [Points overlap an interval only when inside
  it](#zero-length-overlap), so the read would vanish from the next step of a pipeline.
- **Pinned by:** `bamtobed-no-reference-ops-vs-bedtools`,
  `bamtobed-no-reference-ops-bed12-vs-bedtools` and
  `bamtobed-no-reference-ops-split-vs-bedtools` (must differ; they run on `cigars.bam`,
  which has one such read, while the matching cases use `cigars-bases.bam`, the same
  file without it) and `genomecov-ibam-no-reference-ops` (must differ, on `cigars.bam`).

## `convert --score-tag` prints a negative score as a negative number {#bamtobed-negative-score}

- **gloci:** the score of `--score-tag AS` is the tag's integer, signed: `-5` is `-5`.
- **Old behaviour:** `bedtools bamtobed -tag AS` (and `-ed`) prints a negative integer as an
  unsigned 32-bit number: `-1` is `4294967295` and `-5` is `4294967291`.
- **Reason:** the unsigned form is a wraparound of the tag's own value, not a score
  anyone can sort or filter on. Non-negative tags are identical.
- **Pinned by:** `bamtobed-negative-score-vs-bedtools` (must differ).
## `intersect -v` leaves unmapped reads out {#unmapped-reads}

- **gloci:** every interval command skips unmapped reads (flag `0x4`, placed or not) and
  counts them at `-v` (D101). `gloci intersect -a x.bam -b t.bed --mode none` writes the
  mapped reads that have no hit, and no unmapped read.
- **Old behaviour:** `bedtools intersect -abam x.bam -b t.bed -v` writes the unmapped reads
  too, because a read with no position has no hit. With `-u` it never writes them, so
  `-u` and gloci agree.
- **Reason:** an unmapped read is not an interval, so it cannot be "off target". To keep
  every record, use `gloci view` (it writes every record) or `samtools view -f 4`.
- **Pinned by:** `intersect-bam-flags-v-vs-bedtools` (must differ; `flags.bam` has two
  unmapped reads). `intersect-bam-v-vs-bedtools` runs on `order.bam`, which has none, and
  matches.

## An alignment prints as BED6 where bedtools prints BED12 {#bam-as-bed-columns}

- **gloci:** where a command prints the columns of an alignment (`intersect --mode pairs`
  or `--mode count`, `coverage`, an alignment B), it prints the BED6 row `gloci convert`
  writes: `chr1 100 110 read10 60 +`.
- **Old behaviour:** `bedtools intersect -abam x.bam -b t.bed -bed -wa -wb` prints the
  reads as BED12 with `0,0,0` as the colour and trailing commas in the block lists
  (`1 10, 0,`). A BAM `-b` (with `-a x.bed -wa -wb`) is BED6 in bedtools too, and matches. For an A
  with no hit, `-loj` and `-wao` with a BAM `-b` write a 12-column placeholder
  (`. -1 -1 . -1 . . . . . . .`) although the rows with a hit have 6 columns, so
  the file has rows of two widths.
- **Reason:** one row layout for a read everywhere. `gloci convert x.bam --bed12` writes
  bedtools' BED12 for pipelines that need the blocks.
- **Unmatched A:** gloci writes the BED6 placeholder `. -1 -1 . -1 .` (plus `0` with
  `--overlap-bases`), the width of every row.
- **Pinned by:** `intersect-abam-wa-wb-vs-bedtools`, `intersect-bed-bam-loj-vs-bedtools`
  and `intersect-bed-bam-wao-vs-bedtools` (must differ);
  `intersect-bed-bam-wa-wb-vs-bedtools` (BAM as B) matches.

  `view-sam-strict-cigar-zeros-vs-samtools` (exit status must differ).

## CRAM: MD and NM come out as stored {#cram-md-nm}

- **gloci:** a CRAM record keeps the optional fields the file stores, in their
  order: `gloci view x.cram` writes what `samtools view --input-fmt-option
  decode_md=0` writes.
- **Old behaviour:** samtools (htslib's default `decode_md=1`) computes MD and NM
  from the reference for every mapped read that lacks them and appends them after
  its other fields. htslib's CRAM writer drops MD and NM when they match the
  reference, so a BAM converted to CRAM and back reads with them under samtools
  and without them under gloci; a read that never had them gets them from
  samtools.
- **Also:** CRAM keeps neither the MAPQ of an unmapped read (both read 0) nor `=`
  and `X` in a CIGAR (both read `M`); gloci reads these as samtools does.
- **Effect:** `convert --score-tag NM` (bamtobed `-ed`) finds no NM on reads whose
  NM the writer dropped. Recompute the tags with `samtools calmd` first.
- **Reason:** gloci writes what the file holds; computing tags is an edit, and
  needs the reference even where the reads do not.
- **RNEXT:** for an unmapped read placed on a contig whose mate reference was `=`, samtools
  decodes a CRAM with RNEXT `*` where the BAM and gloci have `=`.
- **Pinned by:** `view-cram-md-nm-vs-samtools` (md.cram: samtools adds MD and NM to
  the reads that lack them; must differ),
  `view-cram-md-nm-decode-md-0-vs-samtools` (equal to `decode_md=0`) and
  `view-cram-unmapped-rnext-vs-samtools` (rnext-eq.cram; must differ).

## A CRAM without its end-of-file container is an error {#cram-eof-container}

- **gloci:** a CRAM that does not end with the 38-byte end-of-file container is
  rejected: `the CRAM end-of-file container is missing (...); the file is probably
  truncated` (exit status 1). A file is checked before any record is read; standard
  input fails when it ends.
- **Old behaviour:** samtools prints `EOF marker is absent. The input is probably
  truncated` and decodes what is there, with exit status 0.
- **Pinned by:** `view-cram-no-eof-vs-samtools` (flags-noeof.cram, flags.cram
  without its end-of-file container: exit status must differ).
- **Reason:** as for BGZF ([#bgzf-eof-marker](#bgzf-eof-marker)): a missing end is
  the signature of a cut-off file, and gloci never gives a quiet partial answer.

## An indexed CRAM run at `-t` above 1 does not read the unplaced reads {#cram-unplaced-tail}

- **gloci:** a CRAM with a `.crai` is read chromosome by chromosome at `-t` above 1,
  and the container of reads with no place (unmapped reads at the end of a
  coordinate-sorted file) is in no chromosome's chunk: it is counted from its slice
  headers (the `-v` unmapped count) but never decoded, so it is never verified. A
  corrupt unplaced container is an error (exit 1) when the file is streamed (`-t 1`,
  no `.crai`, standard input) and goes unnoticed (exit 0) in an indexed run at `-t 4`.
  An indexed BAM behaves alike for the reads at the end of a BAM. The same holds for
  a container no chunk reaches for any other reason: only the containers the chunks
  read are checked.
- **Old behaviour:** samtools decodes sequentially and reports the first damage it
  finds (and tolerates some, as it does a checksum mismatch in a block it can read).
- **Pinned by:** the cargo test `cram_unplaced_tail_is_read_only_when_streaming`
  (it needs a damaged file, which a differential case cannot compare).
- **Reason:** reading every container just to verify one that cannot change a result
  would cost the parallel run its point; `gloci check` is the way to verify a file.

## `view` of a region past the end of a reference: BAM and CRAM differ, as samtools' indexes do {#region-past-reference-end}

- **gloci:** for a region that starts past the end of its reference (`chr1:1001-1010`
  with `LN:1000`), `gloci view x.bam` finds the reads that run past the end and
  overlap it (the BAM index bins them by their uncut span), and `gloci view x.cram`
  finds none (the `.crai` clips a slice's span at the end of the reference, and gloci
  selects slices by it exactly).
- **Old behaviour:** samtools does the same: its `.bai` finds the reads and its `.crai`
  does not.
- **Pinned by:** `view-bam-region-past-end-vs-samtools` and
  `view-cram-region-past-end-vs-samtools`, both equal to samtools.
- **Reason:** both follow what the index of that format can find; a region inside the
  reference finds the same reads in both.

## CRAM: limits on what gloci decodes {#cram-limits}

- **gloci:** a read may run past the end of its reference (htslib and gloci decode the bases
  there against `N`) by at most 1 Mb, or four times the longest read of its slice when that
  is more; a read further past the end is an error naming gloci's limit. A few other shapes
  are refused before noodles-cram decodes a record, because it would crash on them:
  - a header block that would expand past 256 MiB;
  - a slice claiming more records than 4,096 plus its bytes;
  - reads claiming more features than the slice stores feature codes;
  - in a slice of several references, a mapped read with no reference;
  - a mate given as a distance past the end of its slice;
  - a block that says it decodes to more than 1 GiB, or whose codec's clear-text header
    (rANS 4x8, rANS Nx16, the arithmetic coder, fqzcomp, the name tokeniser) gives a
    decoded size other than its block header's;
  - a name-tokeniser block that claims more names than its bytes or than the slice's
    records;
  - the Golomb, Subexp and Golomb-Rice encodings, which htslib never writes (rewrite the
    file with `samtools view -C`).
- **Old behaviour:** samtools decodes a read however far past the end it runs, and reads the
  other encodings.
- **Remaining risk:** gloci checks every container, block and slice header (each has a
  CRC32 checksum, so random corruption is always an error), the values above where a data
  series has a block of its own or a constant code (as htslib writes them), and the decoded
  size each block states, in its header and its codec's. A file crafted on purpose, with
  checksums recomputed, can still reach what noodles-cram 0.102 trusts beyond these: its own
  `expect`s and unchecked arithmetic (a record's reference, its mates' template length, its
  alignment span, a base it reconstructs without a reference), and the sizes inside
  compressed streams past their headers (rANS Nx16's RLE and PACK metadata, the tokeniser's
  token streams). A panic there ends gloci with `the CRAM decoder failed on malformed data in
  FILE: container N (at byte B); the file is corrupt or was crafted (see #cram-limits …)`
  (exit status 1, also in release builds: gloci's panic hook exits before the abort); a size
  that asks for more memory than the machine has aborts it (exit status 134). Do not feed
  gloci CRAMs from untrusted sources you would not give to a library that panics; samtools
  has its own share of such limits.
- **Reason:** a read a megabase past its chromosome's end is not real data, and decoding it
  costs that much memory; the rest keeps a corrupt count from asking for gigabytes.
- **Pinned by:** `view-cram-overhang-limit-vs-samtools` (overhang.cram: a deletion ending 1.5 Mb
  past chr1; exit status must differ), and the tests `cram_overhang_limit`,
  `prechecks_refuse_counts_the_slice_cannot_hold`, `tokeniser_header_is_checked`,
  `block_sizes_are_checked`, `crafted_cram_panic_points_to_the_limits` and
  `cram_corrupt_bytes_are_errors_not_panics`.
## A FASTA without an index is an error, never built on the fly {#fasta-implicit-index}

- **gloci:** every command that reads a FASTA at random (`faidx`, `getfasta`,
  `nuc`, a CRAM `--reference`) exits 2 when the `.fai` (or, for BGZF, the
  `.gzi`) is missing, and the message names the fix: `gloci index REF`.
  `gloci faidx REF` with no region exits 2 the same way.
- **Old behaviour:** `samtools faidx REF` builds a missing `.fai` itself, with
  or without a region, and so does any samtools command given a FASTA
  reference.
- **Reason:** an index written as a side effect of a read can be stale for the
  next run, can land in a read-only directory (and fail there), and is a
  file the user never asked for. Building it is one explicit step.
- **Pinned by:** `faidx-no-index-vs-samtools` (exit status: must differ).

## `faidx` skips a missing sequence without writing its header {#faidx-missing-sequence}

- **gloci:** an unknown sequence name is a data error (exit 1) before any
  output, naming the region, as `samtools faidx` also exits 1. With
  `--skip-invalid`, regions on unknown sequences are skipped entirely, with
  one warning that counts them, and the rest are written.
- **Old behaviour:** `samtools faidx -c` writes the `>` line of the region it
  could not read, with no bases, then goes on. Without `-c` it exits 1 after
  writing that header.
- **Reason:** a header with no sequence reads as an empty sequence in the
  output; skipping the region leaves a FASTA that holds only what exists.
- **Pinned by:** `faidx-missing-sequence-continue-vs-samtools` (output must
  differ); `faidx-missing-sequence-exit-vs-samtools` pins the agreeing exit 1.

## `faidx` is strict about region syntax, where samtools guesses {#faidx-region-syntax}

- **gloci:** a region is parsed under gloci's region rules (D8): a 1-based
  start of at least 1, well-formed comma groups, plain digits, and an
  unambiguous name. One that does not parse exits 2 before any output.
  `chr:-END` is accepted as `chr:1-END`, as htslib does, but `chr:-0` is not.
  A region whose end is before its start (`chr1:5-4`) is a header with no
  bases and exit 0, like a start past the end. A coordinate beyond 2^64 is an
  error ("too large").
- **Old behaviour:** `samtools faidx` guesses:
  - `chr1:0`, `chr1:`, `chr1:+5`, `chr1: 5`, `chr1:,5`, `chr1:1e3` and `chr1:-0`
    return real sequence and exit 0 (the whole sequence for `chr1:0`, `chr1:`
    and `chr1:-0`; from 5 for `chr1:+5`, `chr1: 5` and `chr1:,5`; position 1000
    for `chr1:1e3`), and `chr1:1,0` reads as 10;
  - `chr1:0-5` writes the header with no bases and exits 0 (a "Zero length
    sequence" message);
  - `chr1:5-1`, `chr1:5-4` and `chr1:abc` write the header and exit 1;
  - when a contig is literally named `ctgA:5` and another is `ctgA`, the region
    `ctgA:5` (or `ctgA:5-9`) writes the header and exits 1, where gloci exits 2
    and asks for `{ctgA}:5` or `{ctgA:5}`;
  - `chr1:1-99999999999999999999` is clamped to the end of the sequence.
- **Reason:** a typo such as `chr1:0-5`, `chr1:1,0` or `chr1:+5` quietly gives
  the wrong bases or an empty record; an error that shows the region is
  safer. The `{name}` form removes every ambiguity.
- **Pinned by:** `faidx-region-syntax-vs-samtools` (`chr1:0-5`) and
  `faidx-region-ambiguous-vs-samtools` (exit status: must differ);
  `faidx-region-braces-vs-samtools` pins the `{name}` forms, which agree.

## `faidx --regions-file` reads BED {#faidx-regions-file}

- **gloci:** `--regions-file FILE` reads BED lines (0-based, half-open), as in
  every gloci command that takes regions. A region prints as
  `chrom:START-END` with a 1-based start; a zero-length line is skipped with
  one warning that counts them.
- **Old behaviour:** `samtools faidx -r FILE` reads one `chr:START-END` text
  region per line. Given a BED file it treats each whole line as a sequence
  name and exits 1.
- **Reason:** one regions-file format across commands. Text regions are
  positional arguments, or `--region`.
- **Pinned by:** `faidx-regions-file-text-vs-samtools` (exit status: must
  differ); `faidx-regions-file-vs-samtools` pins identical output for the same
  regions written as BED and as text.

## `dict` refuses a repeated sequence name {#dict-duplicate-names}

- **gloci:** `gloci dict` exits 1 when two sequences share a name (the first
  word of the header), naming both lines, and writes no dictionary to `-o FILE`.
- **Old behaviour:** `samtools dict` writes an `@SQ` line for each and exits 0.
- **Reason:** a dictionary with two `SN:a` lines is rejected or misread by the
  tools that take it as `-g`, and the duplicate is nearly always a mistake worth
  catching where the dictionary is made. Rename one of the sequences.
- **Pinned by:** `dict-duplicate-names-vs-samtools` (exit status: must differ).

## `dict` of an empty file, or of text that does not start with `>`, is an error {#dict-empty-input}

- **gloci:** a zero-byte file exits 1 saying it is empty and has no records
  (`X is empty: it has no records`, the message every gloci command gives for
  empty input), and a file
  whose first line is not a `>` header exits 1 saying the format is not
  recognised, as every gloci command that reads a FASTA does.
- **Old behaviour:** `samtools dict` writes the `@HD` line for an empty file and
  exits 0; for a file with sequence lines before the first header it ignores
  them and writes the dictionary of the records that follow.
- **Reason:** both are nearly always a failed upstream step or the wrong file;
  an empty or partial dictionary would hide it.
- **Pinned by:** `dict-empty-file-vs-samtools` and
  `dict-text-before-header-vs-samtools` (exit status: must differ).

## Exit codes follow gloci's, not samtools' {#samtools-exit-codes}

- **gloci:** every command exits 0 on success, 1 for a data or runtime error
  (and for a failed check), and 2 for a usage error: an unknown flag name, a
  bad option value, a missing argument.
- **Old behaviour:** samtools exits 1 for a usage error (it prints its usage
  text and exits 1, as `samtools flags bogus` does), 2 for an error in
  `bedcov`, and from `quickcheck` a bit mask (2 for a file that cannot be
  opened and 16 for a missing end-of-file block, among others, summed when
  several apply).
- **Reason:** one meaning per code across every gloci command, so that a
  pipeline can tell a mistake in the command line from a bad file without
  reading the message. `check --quick` writes the problems it finds to stderr
  and exits 1.
- **Pinned by:** `flags-bad-name-vs-samtools` (exit status: must differ).

## `flags` refuses a leading zero and trailing text {#flags-leading-zero}

- **gloci:** a FLAGS argument written in decimal with a leading zero (`010`)
  exits 2, and the message says samtools reads it as octal 8 and what to write
  instead (`8` or `0x8`). Trailing text (`12abc`) exits 2 as well. This applies
  to every option that takes FLAGS (`--require-flags`, `--exclude-flags`), as
  well as to `gloci flags`.
- **Old behaviour:** `samtools flags 010` prints `0x8 8 MUNMAP`, because C's
  `strtol` with base 0 reads a leading zero as octal; `samtools flags 12abc`
  prints the line for 12 and ignores `abc`.
- **Reason:** `010` is written to mean ten far more often than eight, and a
  flag set that silently means something else filters the wrong reads. Write
  `0x8` or `8` for eight, and `10` for ten.
- **Pinned by:** `flags-octal-vs-samtools` (exit status: must differ), and the
  tests `parse_flags_samtools_rules` and `flags_usage_errors_exit_2_and_print_nothing`.

## `depth` cuts a read at the end of its reference {#depth-past-reference-end}

- **gloci:** a read that runs past the end of its reference (`LN` in `@SQ`) is
  counted only up to `LN`, and one warning says so. No row is written past `LN`,
  and a read that starts at or past `LN` writes nothing.
- **Old behaviour:** `samtools depth` writes rows for positions past `LN`, one
  for each base the read covers there (`chr2 301 2` for a reference of length
  300).
- **Reason:** a position past the end of a reference does not exist; the same
  cut applies to the spans of every other gloci command that reads alignments
  (D101). A pipeline that wants those rows needs the reads, not a depth.
- **Pinned by:** `depth-past-reference-end-vs-samtools` (must differ).

## `depth` stops at a read out of order, wherever it ends {#depth-order-check}

- **gloci:** every decrease of `(reference, position)` within an input is an
  error (exit 1) that names both records and says to run `gloci sort`.
- **Old behaviour:** `samtools depth` notices only a read that ends before its
  last printed row, so a read a few bases before its predecessor that ends
  further on is counted only from where the pileup has got to (its first bases
  have no rows), and the status is 0.
- **Reason:** an unsorted file gives depths that depend on how far out of order
  it is, and a pipeline gets no sign. Sorting is one command.
- **Pinned by:** `depth-unsorted-vs-samtools` (exit status: must differ).

## `depth --all-positions` with regions prints regions on every reference {#depth-regions}

- **gloci:** with `--region` or `--regions-file`, `--all-positions` prints every
  position of the regions, as `--all-references` does, whether or not a read is
  on that reference. That is what `samtools depth -a -r REGION` and
  `samtools depth -aa -b FILE` print.
- **Old behaviour:** `samtools depth -a -b FILE` skips the regions on a
  reference that has no passing read, so a region on an empty reference writes
  nothing, where one on a reference with reads writes its zero rows.
- **Reason:** the rows should not depend on whether some other read of the
  reference passed a filter, and `-r` and `-b` should agree.
- **Pinned by:** `depth-all-positions-regions-file-vs-samtools` (must differ),
  and the test `depth_all_modes`.

## `depth` and `bedcov` refuse inputs with different `@SQ` lists {#depth-bedcov-same-sq}

- **gloci:** every input must have the same `@SQ` lines (names, lengths and
  order). If not, the command exits 1 before writing anything and names the first
  line that differs (`@SQ 1 is chr2 (LN 300) in the first and chr1 (LN 500) in
  the second`).
- **Old behaviour:** `samtools depth` and `samtools bedcov` compare reference ids
  across files without checking, so files that list the same references in
  another order are read as if their ids meant the same references, and the rows
  carry the first file's names for the other file's depths. The status is 0.
- **Reason:** the depths of different references in one row are wrong, and
  nothing in the output shows it.
- **Pinned by:** `depth-reordered-references-vs-samtools` (exit status: must
  differ).

## `depth --min-baseq` counts every base of a read with SEQ `*` {#depth-seq-star-baseq}

- **gloci:** a read whose SEQ is `*` has no base qualities, so every base of it
  passes `--min-baseq`, as a missing quality (QUAL `*`) does.
- **Old behaviour:** `samtools depth -q N` looks up the quality of each base
  past the end of such a record, where the quality bytes would be. It reads
  whatever follows the CIGAR in the record: the record's own aux bytes, and
  past those stale heap memory. Some bases count and some do not, in a pattern
  nobody chose.
- **Reason:** an answer that depends on the memory layout cannot be matched;
  counting the read is the only reading that does not invent a quality.
- **Pinned by:** `depth-seq-star-min-baseq-vs-samtools` (must differ; the
  fixture's SEQ `*` reads carry an aux array as long as their CIGAR, so the
  bytes samtools reads are the file's own and its output is the same on every
  run) and the cargo test
  `depth_noseq_reads_use_their_cigar_for_the_length`, which checks gloci's side.
