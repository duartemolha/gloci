#!/bin/sh
# Rebuilds the BAM fixtures and their BAI indexes from the SAM sources (D119).
# Run it in the reference image, from this directory:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work gloci-diff:p0 sh make-alignments.sh
# One thread, level 6, no @PG line, so the output is reproducible.
set -eu
for set in cat-a cat-b cat-c cigars flags order tags edge pastend; do
    samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o "$set.bam" "$set.sam"
    # cat-* are not indexed; edge.sam is not in coordinate order, so samtools cannot index it.
    case "$set" in cat-*|edge) ;; *) samtools index -@ 0 "$set.bam" ;; esac
done
