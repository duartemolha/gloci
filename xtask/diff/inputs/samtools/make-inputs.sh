#!/bin/sh
# Regenerates the inputs of the `gloci dict` cases that are not hand-written
# FASTA: a plain (non-BGZF) gzip of fa/ref.fa and an empty file. dup.fa
# (a repeated name) and alt.txt (an alternate-locus list: first tab-separated
# field is the name, "MT\r" ends in CRLF, "M extra" and "1 " match nothing)
# and nohdr.fa (sequence text before the first header) are written by hand. Run from this directory.
set -eu
gzip -n -c ../fa/ref.fa > ref-gzip.fa.gz
: > empty.fa
