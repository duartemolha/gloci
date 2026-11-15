#!/bin/sh
# Rebuilds the BAM and CRAM inputs of the `gloci tag` cases from reads.sam. Run it in the reference
# image, from this directory:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work gloci-diff:p0 sh make-tag.sh
# One thread, level 6, no @PG line, so the output is reproducible. reads.sam is written by hand:
# its reads are in no order (the output keeps it), and each one sets a trap that the comment
# above the tag cases in xtask/diff/cases.toml describes. The CRAM is written without a
# reference (no_ref), so the cases read it alone; it keeps `=` and `X` operations as `M`, so
# reads-cram.bam, its records decoded back to BAM, is what bedtools tags in the CRAM case.
set -eu
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o reads.bam reads.sam
samtools view -C -@ 0 --output-fmt-option no_ref=1 --no-PG -o reads.cram reads.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o reads-cram.bam reads.cram
