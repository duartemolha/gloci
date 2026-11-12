"""gloci-py: merge overlapping and bookended records. Usage: IN_A OUT."""

import sys

import gloci

a, out = sys.argv[1:3]
gloci.open(a).merge(chrom_order="lex").save(out)
