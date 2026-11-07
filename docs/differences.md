# Differences from bedtools, bedops, vcftools and htslib

gloci aims to produce the same output as the tools it replaces. Where it
deliberately does not, the difference is listed here with the reason. Nothing
else may differ: the differential harness (`cargo xtask diff`) runs gloci and
pinned reference tools (bedtools 2.31.1, bedops 2.4.41, vcftools 0.1.16,
htslib 1.21, bcftools 1.21, GNU coreutils) on the same inputs, and
`crates/gloci/tests/golden/` holds the reference outputs, which every
`cargo test` checks without Docker.

Each section has an anchor. A case in `xtask/diff/cases.toml` that names a
section as its `expected_difference` must produce different output from the
reference tool, so a section whose behaviour stops differing fails the build
until the section is removed. Cases without an `expected_difference` must match
the reference output exactly (or under their compare mode); `cases.toml` and
the golden directory list them all.

Each section gives the gloci behaviour, the old behaviour and the reason.

## Default `--vcf-span allele` {#default-vcf-span-allele}

- **gloci:** a VCF record's footprint is the reference bases its ALT alleles
  change (spec section 4.2): the shared prefix and suffix of REF and ALT are
  trimmed, and a pure insertion is a point. `--vcf-span ref` gives the REF span.
- **Old behaviour:** bedtools, tabix and bcftools use `[POS-1, POS-1+len(REF))`
  (and INFO/END for symbolic alleles). That is `--vcf-span ref`.
- **Reason:** the REF span ignores the alleles, so for example
  `bedtools subtract` against a VCF removes the wrong number of bases for an
  indel or a padded SNV. Pass `--vcf-span ref` for the old result.
- **Pinned by:** `view-region-ref-inclusive-vs-tabix` (must match) and
  `view-region-default-vs-tabix` (must differ).

## Default karyotype chromosome order {#karyotype-chromosome-order}

- **gloci:** without `-g` and without a VCF/BCF header that lists contigs,
  chromosomes are ordered `chr1, chr2 ... chr22, chrX, chrY, chrM` (decision
  D32). `--chrom-order lex` selects byte order.
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
- **Reason:** a total, documented order (start, then end) makes the output
  independent of input order for the common case, and stability keeps the
  rest predictable. Use `--full-key` to reproduce `sort-bed` exactly.

## `view` region footprints and points {#view-region-footprints}

- **gloci:** `view` selects records whose footprint (per `--vcf-span`)
  overlaps the region under the half-open rules of spec section 4.3. A point
  (a zero-length BED record, an insertion under `allele`) overlaps a region
  only when strictly inside it, so an insertion on a region's edge is
  excluded. `--inclusive-points` includes points on the edge.
- **Old behaviour:** `tabix` returns every record whose REF span touches the
  1-based region, so an insertion anchored on the first or last base of the
  region is returned.
- **Reason:** an insertion between two bases belongs to neither base, so it is
  in a region only when the region contains both neighbours (decision D9).
  `--vcf-span ref --inclusive-points` reproduces tabix.
- **Pinned by:** `view-region-ref-inclusive-vs-tabix` (must match),
  `view-region-default-vs-tabix` (must differ).

## `##glociCommand` header line {#glocicommand-header}

- **gloci:** VCF and BCF output of `sort` and `view` carries the original header
  plus one `##glociCommand=<command line>; Version=<version>` line.
- **Old behaviour:** `bcftools sort` (1.21) copies the header unchanged;
  `tabix` prints no header for a region query; GNU sort and `sort-bed` never
  add header lines.
- **Reason:** provenance in the file, in one line a tool can recognise.
- **Pinned by:** the `vcf-records` compare mode, which ignores header lines
  (`sort-vcf-vs-bcftools-sort`). Header text is therefore not compared against
  the reference tools.

## `sort --check` checks `(chrom, start)` only {#sort-check-chrom-start}

- **gloci:** `sort --check` accepts a file ordered by chromosome and start. Two
  records with the same start may be in any order of end.
- **Old behaviour:** `sort-bed --check-sort` also requires ends to be ordered
  when starts are equal ("not properly sorted by end coordinates when start
  coordinates are identical"), and `sort -c -k1,1 -k2,2n` rejects such a file
  through its last-resort whole-line comparison.
- **Reason:** `(chrom, start)` is the order every streaming command needs
  (decision D12), and files sorted by `bedtools sort` or `sort -k1,1 -k2,2n`,
  which leave ends unordered, must be accepted. A file that `gloci sort` wrote
  passes as well.

## A missing BGZF end-of-file marker is an error {#bgzf-eof-marker}

- **gloci:** a seekable BGZF file without the 28-byte end-of-file block is
  rejected: `BGZF end-of-file marker missing; the file is probably truncated`
  (decision D30). Standard input cannot be checked.
- **Old behaviour:** htslib (`tabix`, `bcftools`) prints the warning `EOF
  marker is absent. The input may be truncated` and carries on with exit
  status 0.
- **Reason:** a missing marker is the signature of a truncated download or an
  interrupted write, and gloci never gives a quiet partial answer.

## `view --assume-sorted` {#view-assume-sorted}

- **gloci:** for an uncompressed BED file with no index, `view` scans the
  chromosome's records and checks their order, so an unsorted file is an
  error. `--assume-sorted` binary-searches instead (like `bedextract`) and trusts
  the file to be sorted in the order given by `-g` or `--chrom-order`; it checks
  only the block it reads, and on an unsorted file the answer can be incomplete
  (decision D11, revised).
- **Old behaviour:** `bedextract` always trusts the order. bedtools and tabix
  read an index or the whole file.
- **Reason:** correct by default, fast on request. Index the file
  (`gloci sort FILE -o FILE.gz --write-index`) to get speed and checking.

## Indexes built by htslib and SVs with only `SVLEN` {#htslib-index-svlen}

- **gloci:** the footprint of a symbolic SV without INFO/END is
  `POS + |SVLEN|` (spec section 4.2), and the index spans gloci builds include
  it. Queries on a gloci-built index are correct, and htslib reads that index
  correctly (`index-interop-tabix-queries`).
- **Old behaviour:** htslib ignores SVLEN. An index built by `tabix` or `bcftools
  index` spans such a record only over its REF bases, so a region that touches
  only the tail of an SVLEN-only SV is missed when gloci reads it through that
  index.
- **Reason:** the SVLEN-derived end is part of the footprint rules. Rebuild the
  index with `gloci index` if a file has SVLEN-only SVs and you query regions
  inside them. Records with INFO/END are unaffected.

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

## Zero-length BED records

- **gloci:** accepted as points (spec section 4.1).
- **Old behaviour:** `sort-bed` rejects them (`Genomic end coordinate is less
  than (or equal to) start coordinate`); `bedtools sort` and GNU sort accept
  them. The sort-bed cases therefore use an input without points
  (`ties-nopoints.bed`), and the points are covered by the bedtools and GNU
  cases.
