"""Dump a bigWig or bigBed through pyBigWig, one line per record.

usage: bbi_dump.py FILE [--lex] [--region CHROM:START-END]...
Regions are 1-based and closed. bigWig intervals are clipped to the region;
bigBed entries are printed whole when they overlap it.
"""
import argparse
import decimal
import math
import struct
import sys

import pyBigWig


def f32(v: float) -> str:
    if v == 0.0:
        return "-0" if math.copysign(1.0, v) < 0 else "0"
    for digits in range(1, 10):
        s = "%.*e" % (digits - 1, v)
        if struct.unpack("<f", struct.pack("<f", float(s)))[0] == v:
            break
    t = format(decimal.Decimal(s), "f")
    return t.rstrip("0").rstrip(".") if "." in t else t


def parse_region(text):
    chrom, _, span = text.rpartition(":")
    start, _, end = span.partition("-")
    return chrom, int(start) - 1, int(end)


def dump_span(bw, is_bw, out, chrom, start, end):
    if start >= end:
        return
    if is_bw:
        for s, e, v in bw.intervals(chrom, start, end) or []:
            s, e = max(s, start), min(e, end)
            out.write("%s\t%d\t%d\t%s\n" % (chrom, s, e, f32(v)))
    else:
        for s, e, rest in bw.entries(chrom, start, end) or []:
            if s < end and e > start:
                line = "%s\t%d\t%d" % (chrom, s, e)
                if rest:
                    line += "\t" + rest
                out.write(line + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--lex", action="store_true")
    ap.add_argument("--region", action="append", default=[])
    args = ap.parse_args()
    bw = pyBigWig.open(args.file)
    if bw is None:
        sys.exit("cannot open " + args.file)
    is_bw = bw.isBigWig()
    chroms = bw.chroms()
    out = sys.stdout
    if args.region:
        for r in args.region:
            chrom, start, end = parse_region(r)
            dump_span(bw, is_bw, out, chrom, start, end)
    else:
        names = sorted(chroms) if args.lex else list(chroms)
        for chrom in names:
            dump_span(bw, is_bw, out, chrom, 0, chroms[chrom])
    bw.close()


if __name__ == "__main__":
    main()
