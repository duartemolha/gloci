#!/bin/sh
# Rebuilds the BAM fixtures and their BAI indexes from the SAM sources (D119).
# Run it in the reference image, from this directory:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work gloci-diff:p0 sh make-alignments.sh
# One thread, level 6, no @PG line, so the output is reproducible.
set -eu
# cat-big: 1300 reads of 100 bases with pseudo-random sequence and qualities, so the BAM spans
# several BGZF blocks (cat tests copy whole blocks). Same references as cat-a.
awk 'BEGIN { srand(7); print "@HD\tVN:1.6\tSO:coordinate"; print "@SQ\tSN:chr1\tLN:100000";
    print "@SQ\tSN:chr2\tLN:100000"; print "@RG\tID:g4\tSM:s4"
    for (i = 0; i < 1300; i++) { seq = ""; qual = ""
        for (j = 0; j < 100; j++) { seq = seq substr("ACGT", int(rand() * 4) + 1, 1)
            qual = qual sprintf("%c", 33 + int(rand() * 40)) }
        printf "big%d\t0\tchr1\t%d\t60\t100M\t*\t0\t0\t%s\t%s\tRG:Z:g4\n", i, i + 1, seq, qual } }' > cat-big.sam
# cigars-bases: cigars without the all-clip read, which bedtools writes zero-length (D101).
grep -v '^c_allclip' cigars.sam > cigars-bases.sam
for set in cat-a cat-b cat-c cat-d cat-big cigars cigars-bases flags order tags edge pastend negtag splice-end; do
    samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o "$set.bam" "$set.sam"
    # cat-* are not indexed; edge.sam is not in coordinate order, so samtools cannot index it.
    case "$set" in cat-*|edge) ;; *) samtools index -@ 0 "$set.bam" ;; esac
done
# cat-noeof: cat-a without its BGZF end-of-file block (a truncated copy).
head -c -28 cat-a.bam > cat-noeof.bam
# cat-g.bam is not made here: it is `gloci view cat-b.sam -o cat-g.bam`, so its header ends in
# the middle of a block (samtools puts the header in a block of its own).
