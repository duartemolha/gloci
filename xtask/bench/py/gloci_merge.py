"""gloci-py: merge overlapping and bookended records.

Usage: IN_A OUT [--threads N]. Without --threads the library's default is used.
"""

import sys

import gloci

args = sys.argv[1:]
threads = int(args.pop(args.index("--threads") + 1)) if "--threads" in args else None
if threads is not None:
    args.remove("--threads")
a, out = args
gloci.open(a).merge(chrom_order="lex").save(out, threads=threads)
