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

# Task 6: the `depth` fixtures, all written by hand. depth.sam sets the traps: a deletion whose
# next base has quality 10, `N`, `S`, `I`, `=`/`X`, a CIGAR `*` (samtools reads a mapped read with
# one as unmapped), QUAL `*`, MAPQ 29 and 30, duplicate, secondary, QC-fail and supplementary
# reads, a placed unmapped read, a read ending at chr2's LN, reads of 30, 49 and 50 bases once soft
# clips are cut, and an empty chr3. depth2.sam has other reads on the same references;
# depth-reordered.sam lists them in another order; depth-unsorted.sam has a read 3 bp before its
# predecessor that ends after samtools' last row. depth-noseq.sam has reads with SEQ `*` (samtools
# reads past the record for their base qualities, so only -l is compared on it). depth-pastend.sam
# has reads that run past chr2's LN (samtools prints rows beyond it, so no other fixture has one).
# depth.cram is
# written without a reference (every base stored), so reading it needs none. depth-idx/input.bai
# is the index under the name samtools and gloci find beside a file called input.bam.
for s in depth depth2 depth-reordered depth-unsorted depth-noseq depth-pastend; do
    samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o "$s.bam" "$s.sam"
done
samtools index -@ 0 depth.bam
samtools index -@ 0 depth2.bam
samtools index -@ 0 depth-noseq.bam
samtools index -@ 0 depth-pastend.bam
mkdir -p depth-idx
samtools index -@ 0 depth.bam depth-idx/input.bai
samtools view -C -@ 0 --no-PG --output-fmt-option version=3.0 --output-fmt-option no_ref=1 -o depth.cram depth.sam
samtools index -@ 0 depth.cram
