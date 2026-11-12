"""gloci-py: A records overlapping any B record (bedtools intersect -u). Usage: IN_A IN_B OUT."""

import sys

import gloci

a, b, out = sys.argv[1:4]
gloci.open(a).intersect(b=b, mode="any", chrom_order="lex").save(out)
