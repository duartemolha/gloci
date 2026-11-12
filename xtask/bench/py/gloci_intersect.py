"""gloci-py: A records overlapping any B record (bedtools intersect -u).

Usage: IN_A IN_B OUT [--threads N]. Without --threads the library's default is used.
"""

import sys

import gloci

args = sys.argv[1:]
threads = int(args.pop(args.index("--threads") + 1)) if "--threads" in args else None
if threads is not None:
    args.remove("--threads")
a, b, out = args
gloci.open(a).intersect(b=b, mode="any", chrom_order="lex").save(out, threads=threads)
