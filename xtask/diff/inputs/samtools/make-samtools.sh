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
