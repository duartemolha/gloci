"""bioframe: A records overlapping any B record. Usage: IN_A IN_B OUT."""

import sys

import bioframe

COLS = ["chrom", "start", "end"]
a, b, out = sys.argv[1:4]


def read(path):
    return bioframe.read_table(path, names=COLS, usecols=[0, 1, 2])


df_a, df_b = read(a), read(b)
hit = bioframe.count_overlaps(df_a, df_b)["count"] > 0
df_a[hit].to_csv(out, sep="\t", header=False, index=False)
