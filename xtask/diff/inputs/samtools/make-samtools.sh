#!/bin/sh
# Rebuilds the inputs of the samtools-essentials cases (P3a faidx/dict, P3c) that are not
# hand-written (D119, D141). Run it in the reference image, from the directory above this one:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work/samtools gloci-diff:p0 sh make-samtools.sh
# One thread, level 6, no @PG line, so the output is reproducible.
set -eu

# Task 10: the inputs of the `gloci dict` cases that are not hand-written FASTA: a plain (non-BGZF)
# gzip of fa/ref.fa and an empty file. dup.fa (a repeated name) and alt.txt (an alternate-locus
# list: first tab-separated field is the name, "MT\r" ends in CRLF, "M extra" and "1 " match
# nothing) and nohdr.fa (sequence text before the first header) are written by hand.
gzip -n -c ../fa/ref.fa > ref-gzip.fa.gz
: > empty.fa

# Task 1: empty.sam and empty.bam, a header and no records, for the cases whose input neither
# side reads (`flags`) and for the commands that must cope with a file of no reads. The header
# has two references and a read group; the BAM has its .bai.
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@RG\tID:g1\tSM:s1\n' > empty.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o empty.bam empty.sam
samtools index -@ 0 empty.bam

# Task 2: unsorted-chroms.{sam,bam}, the reads of chr2 between two runs of chr1 reads and an
# unmapped read with no place at the end: `idxstats` must count by reference, not by run, where
# samtools' scan stops at "file is not position sorted". No @HD SO: the file claims no order. The
# BAM has no index (it could not be sorted into one). nometa/flags.bam is a copy of the P3a BAM
# (../aln/flags.bam); its .bai without metadata pseudo-bins is written by the `#[ignore]` test
# `make_nometa_index` (crates/gloci/tests/idxstats.rs), not by samtools.
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr3\tLN:3000\n' > unsorted-chroms.sam
{
  printf 'a1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'a2\t16\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'b1\t0\tchr2\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'b2\t4\tchr2\t30\t0\t*\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'a3\t0\tchr1\t50\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'a4\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n'
} >> unsorted-chroms.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o unsorted-chroms.bam unsorted-chroms.sam
mkdir -p nometa
cp ../aln/flags.bam nometa/flags.bam

# Task 2: a CSI of flags.bam written by `samtools index -c`, in csi/ beside a copy of the BAM so
# the case can name both as aux files.
mkdir -p csi
cp ../aln/flags.bam csi/flags.bam
samtools index -c -@ 0 csi/flags.bam

# Task 3: flagstat.{sam,bam}, built by the awk below and sorted by samtools. QC-passed: 320 records,
# 306 mapped (95.63% in samtools' f32 arithmetic, 95.62% in f64), 160 paired primaries of which 7
# are singletons (4.37% against 4.38%), 14 unmapped (7 mates of singletons, 7 single-end; placed
# and unplaced, one mate carrying 0x2), 60+26+20+40 mapped pairs (proper; proper duplicates; not
# proper; mate on another reference at MAPQ 4, MAPQ 5, and two with RNEXT *), 40 secondary (10
# duplicates, 10 flagged paired), 40 supplementary (10 duplicates), 8 single-end duplicates. Then
# 12 QC-failed (0x200) records, none paired, so the paired percentages there print N/A.
{
  printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n'
  awk 'function e(n, f, rn, p, mq, cg, rx, px) {
         printf "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t0\tACGTA\tIIIII\n", n, f, rn, p, mq, cg, rx, px
       }
       function pair(n, f1, f2, r1, r2, q1, q2) {
         p += 10
         e(n, f1, r1, p, q1, "5M", (r1 == r2 ? "=" : r2), p + 100)
         e(n, f2, r2, p + 100, q2, "5M", (r1 == r2 ? "=" : r1), p)
       }
       BEGIN {
         p = 100
         for (i = 1; i <= 30; i++) pair("pa" i, 99, 147, "chr1", "chr1", 60, 60)
         for (i = 1; i <= 13; i++) pair("pf" i, 1123, 1171, "chr1", "chr1", 60, 60)
         for (i = 1; i <= 10; i++) pair("pb" i, 65, 145, "chr2", "chr2", 60, 60)
         for (i = 1; i <= 9; i++) pair("pq4" i, 65, 129, "chr1", "chr2", 4, 4)
         for (i = 1; i <= 10; i++) pair("pq5" i, 65, 129, "chr1", "chr2", 5, 60)
         for (i = 1; i <= 2; i++) { p += 10; e("pstar" i, 65, "chr2", p, 60, "5M", "*", 0) }
         for (i = 1; i <= 7; i++) {
           p += 10
           e("s" i, 73, "chr1", p, 60, "5M", "=", p)
           if (i <= 4) e("s" i, 133, "chr1", p, 0, "*", "=", p)
           else if (i <= 6) e("s" i, 133, "*", 0, 0, "*", "chr1", p)
           else e("s" i, 135, "chr1", p, 0, "*", "=", p)
         }
         for (i = 1; i <= 73; i++) {
           p += 10
           f = (i % 2) ? 0 : 16
           if (i <= 8) f += 1024
           e("e" i, f, "chr2", p, 30, "5M", "*", 0)
         }
         for (i = 1; i <= 4; i++) { p += 10; e("u" i, 4, "chr2", p, 0, "*", "*", 0) }
         for (i = 5; i <= 7; i++) e("u" i, 4, "*", 0, 0, "*", "*", 0)
         for (i = 1; i <= 40; i++) {
           p += 10
           f = (i <= 20) ? 256 : (i <= 30) ? 1280 : 257
           e("sec" i, f, "chr1", p, 0, "5M", "*", 0)
         }
         for (i = 1; i <= 40; i++) {
           p += 10
           f = (i <= 20) ? 2048 : (i <= 30) ? 2064 : 3072
           e("sup" i, f, "chr2", p, 60, "5M", "*", 0)
         }
         split("512 512 512 512 528 528 516 516 768 2560 1536 1536", qf, " ")
         for (i = 1; i <= 12; i++) {
           p += 10
           if (i == 8) e("q" i, qf[i], "*", 0, 0, "*", "*", 0)
           else e("q" i, qf[i], "chr1", p, (qf[i] == 516 ? 0 : 20), (qf[i] == 516 ? "*" : "5M"), "*", 0)
         }
       }'
} > flagstat-unsorted.sam
samtools sort -@ 0 --no-PG -O sam -o flagstat.sam flagstat-unsorted.sam
rm flagstat-unsorted.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o flagstat.bam flagstat.sam
samtools index -@ 0 flagstat.bam

# Task 4: the inputs of the `check` cases. noeof.bam is flagstat.bam without its 28-byte BGZF
# end-of-file block; truncated.bam is cut in the middle of a block. nosq.{sam,bam} are unaligned
# reads with a read group and no @SQ line. coord-unsorted.{sam,bam} claim SO:coordinate and have
# one read (record 7) before the one that precedes it; no case uses it (samtools has nothing to
# say about order), the tests of `gloci check` do. badfirst.bam (a header and an end-of-file
# block that are fine, a first record whose block_size is 20) is written by the `#[ignore]` test
# `make_broken_bams` (crates/gloci/tests/check_alignments.rs), not here.
head -c -28 flagstat.bam > noeof.bam
head -c 1500 flagstat.bam > truncated.bam
printf '@HD\tVN:1.6\n@RG\tID:g1\tSM:s1\n' > nosq.sam
{
  printf 'u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g1\n'
  printf 'u2\t4\t*\t0\t0\t*\t*\t0\t0\tGGGTA\tIIIII\tRG:Z:g1\n'
  printf 'u3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTA\tIIIII\tRG:Z:g1\n'
} >> nosq.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o nosq.bam nosq.sam
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n' > coord-unsorted.sam
{
  for p in 100 200 300 400 500 600; do
    printf 'r%s\t0\tchr1\t%s\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n' "$p" "$p"
  done
  printf 'rback\t0\tchr1\t150\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  for p in 700 800; do
    printf 'r%s\t0\tchr1\t%s\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n' "$p" "$p"
  done
  printf 'q1\t0\tchr2\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
} >> coord-unsorted.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o coord-unsorted.bam coord-unsorted.sam
# calls.vcf.gz: a small BGZF VCF (bgzip's own output, with its EOF block), for the check case that
# samtools quickcheck rejects as "not sequence data" and gloci accepts.
bgzip -c ../cluster-in.vcf > calls.vcf.gz
