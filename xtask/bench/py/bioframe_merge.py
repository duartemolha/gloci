"""bioframe: merge overlapping and bookended records. Usage: IN_A OUT."""

import sys

import bioframe

a, out = sys.argv[1:3]
df = bioframe.read_table(a, names=["chrom", "start", "end"], usecols=[0, 1, 2])
bioframe.merge(df, min_dist=0)[["chrom", "start", "end"]].to_csv(
    out, sep="\t", header=False, index=False
)
