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

# Task 11: the inputs of the `gloci reheader` cases. reheader-in.sam is hand-written (three
# references of different lengths, so a swapped order shows in the LNs, and one unplaced read);
# new-header.sam renames its references to 1, 2, 3, appends a fourth @SQ and changes the @RG;
# bad-ln.sam changes chr2's LN and bad-order.sam swaps chr2 and chr3 (both are hand-written).
# new-header.bam is new-header.sam as a BAM header file. new-header-cram.sam is the header
# for aln/flags.embed.cram (chr1 1000, chr2 1234, with its M5s), renamed, with a third @SQ added.
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o reheader-in.bam reheader-in.sam
samtools index -@ 0 reheader-in.bam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o new-header.bam new-header.sam
