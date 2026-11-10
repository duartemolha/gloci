#!/bin/sh
# Rebuilds the BAM fixtures and their BAI indexes from the SAM sources (D119).
# Run it in the reference image, from this directory:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work gloci-diff:p0 sh make-alignments.sh
# One thread, level 6, no @PG line, so the output is reproducible.
set -eu
for set in cigars flags order tags; do
    samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o "$set.bam" "$set.sam"
    samtools index -@ 0 "$set.bam"
done
