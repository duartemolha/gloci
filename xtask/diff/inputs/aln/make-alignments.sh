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
# midempty: reads on c1 and c3 only, so c2 (in the middle of the header) has none.
# cigars-small: cigars-bases on a 1000-base reference, so a per-base (-d) golden stays small.
sed 's/LN:100000/LN:1000/' cigars-bases.sam > cigars-small.sam
for set in cat-a cat-b cat-c cat-d cat-big cigars cigars-bases cigars-small midempty flags order tags edge pastend negtag splice-end splitfrac longref; do
    samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o "$set.bam" "$set.sam"
    # cat-* are not indexed; edge.sam is not in coordinate order, so samtools cannot index it;
    # longref has a 600 Mb reference, which a .bai cannot hold (it needs `samtools index -c`).
    case "$set" in cat-*|edge|longref) ;; *) samtools index -@ 0 "$set.bam" ;; esac
done
# longref.sam is written by hand: a 600,000,000-base reference, which needs a .csi.
# index/input.bai: order.bam's index under the name gloci finds beside `input.bam` (x.bam -> x.bai),
# for the golden case view-bam-region-vs-samtools, whose input is always called input.bam.
mkdir -p index
cp order.bam.bai index/input.bai
# pos-minus-one.bam: reads with a reference and POS 0 (BAM pos -1), as the aligners write unmapped mates.
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o pos-minus-one.bam pos-minus-one.sam
samtools index -@ 0 pos-minus-one.bam
# cat-noeof: cat-a without its BGZF end-of-file block (a truncated copy).
head -c -28 cat-a.bam > cat-noeof.bam
# cat-g.bam is not made here: it is `gloci view cat-b.sam -o cat-g.bam`, so its header ends in
# the middle of a block (samtools puts the header in a block of its own).
