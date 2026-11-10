#!/usr/bin/env python3
"""Writes input.vcf: FORMAT and INFO fields of every type and Number, each
missing, cut short, mixed and present, at several FORMAT positions."""
import itertools

types = ["Integer", "Float", "String", "Character"]
numbers = [("1", "1"), ("A", "A"), ("R", "R"), ("G", "G"), (".", "D"), ("2", "2")]
short = {"Integer": "I", "Float": "F", "String": "S", "Character": "C"}

def fid(t, n):
    return short[t] + n[1]

# values of a field per sample ("present" and a second flavour), 1 ALT diploid
def vals(t, n, k):
    cnt = {"1": 1, "A": 1, "R": 2, "G": 3, ".": [1, 3, 2][k % 3], "2": 2}[n[0]]
    pool = {
        "Integer": ["1", "20", "300", "70000", "-5", "7"],
        "Float": ["0.5", "1.25", "300.5", "-2", "7", "0.125"],
        "String": ["ab", "c", "defgh", "x", "yy", "zzz"],
        "Character": ["a", "b", "c", "d", "e", "f"],
    }[t]
    return ",".join(pool[(k + i) % len(pool)] for i in range(cnt))

lines = ["##fileformat=VCFv4.3", "##contig=<ID=chr1,length=100000>"]
lines.append('##FILTER=<ID=q10,Description="low">')
lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">')
for t, n in itertools.product(types, numbers):
    lines.append(f'##FORMAT=<ID={fid(t, n)},Number={n[0]},Type={t},Description="{t} Number={n[0]}">')
for t, n in itertools.product(["Integer", "Float", "String"], numbers[:2] + numbers[4:5]):
    lines.append(f'##INFO=<ID=N{fid(t, n)},Number={n[0]},Type={t},Description="{t} Number={n[0]}">')
lines.append('##INFO=<ID=FL,Number=0,Type=Flag,Description="flag">')
lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
lines.append('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="svlen">')
lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3")

only = []
pos = 0
def rec(fmt, cells, info=".", alt="G", ref="A"):
    global pos
    pos += 1
    lines.append("\t".join(["chr1", str(pos), ".", ref, alt, ".", ".", info, fmt] + cells))

gts = ["0/1", "1/1", "./."]
layouts = [("F:DP", 0), ("DP:F", 1), ("GT:F:DP", 1), ("GT:DP:F", 2)]
for t, n in itertools.product(types, numbers):
    f = fid(t, n)
    for layout, fpos in layouts:
        keys = layout.replace("F", f).split(":")
        def cells_for(kind):
            out = []
            for k in range(3):
                parts = []
                for j, key in enumerate(keys):
                    if key == "GT":
                        parts.append(gts[k])
                    elif key == "DP":
                        parts.append(str(5 + k))
                    else:
                        v = {"present": vals(t, n, k)}.get(kind, ".")
                        if kind in ("mixed", "mixedcut"):
                            v = [vals(t, n, k), ".", vals(t, n, k + 1)][k]
                        if kind == "dots2":
                            v = [".", ".,.", vals(t, n, k)][k] if n[0] != "1" else [".", ".", vals(t, n, k)][k]
                        parts.append(v)
                if kind in ("cut", "mixedcut"):
                    # drop the field and everything after it, in every sample
                    # (or only in the first two for mixedcut)
                    if kind == "cut" or k < 2:
                        parts = parts[:fpos]
                        if not parts:
                            parts = ["."]
                out.append(":".join(parts))
            return out
        for kind in ["dot", "cut", "mixed", "mixedcut", "dots2", "present"]:
            rec(layout.replace("F", f), cells_for(kind))

# empty subfield and trailing-empty
rec("GT:DP:S1:C1", ["0/1:5::x", "1/1:6:ab:y", "./.:7:c:z"])
rec("GT:DP:SD:CD", ["0/1:5:a,b:x", "1/1:6::y", "./.:7:c:"])
rec("GT:DP:ID", ["0/1:5:", "1/1:6:1,2", "./.:7:"])
rec("GT:DP:FD", ["0/1:5:", "1/1:6:", "./.:7:"])

# multi-allelic records
for t, n in itertools.product(types, numbers):
    f = fid(t, n)
    def v(k):
        cnt = {"1": 1, "A": 2, "R": 3, ".": 2 + k % 2, "2": 2, "G": 6}[n[0]]
        pool = {"Integer": "1 2 3 4 5 6", "Float": "0.5 1.5 2.5 3.5 4.5 5.5",
                "String": "a bb ccc d ee fff", "Character": "a b c d e f"}[t].split()
        return ",".join(pool[(k + i) % 6] for i in range(cnt))
    rec(f"GT:{f}:DP", ["0/1:" + v(0) + ":5", "1/2:.:6", "2/2:" + v(1) + ":7"], alt="G,C")
    rec(f"GT:{f}:DP", ["0/1:.:5", "1/2:.:6", "2/2:.:7"], alt="G,C")

# INFO
infos = []
for t in ["Integer", "Float", "String"]:
    for n in [numbers[0], numbers[1], numbers[4]]:
        f = "N" + fid(t, n)
        scal = {"Integer": ["5", "70000"], "Float": ["0.5", "2"], "String": ["abc", "d"]}[t]
        infos += [f"{f}=.", f"{f}={scal[0]}", f"{f}={scal[1]}"]
        if n[0] != "1":
            infos += [f"{f}=.,.", f"{f}={scal[0]},.", f"{f}=.,{scal[1]}", f"{f}={scal[0]},{scal[1]}"]
for i in infos:
    rec("GT:DP", ["0/1:5", "1/1:6", "./.:7"], info=i)
rec("GT:DP", ["0/1:5", "1/1:6", "./.:7"], info="NI1=.;NF1=.;NS1=.;FL")
rec("GT:DP", ["0/1:5", "1/1:6", "./.:7"], info="NI1=3;NFD=.,1.5;NSA=.;FL")
rec("GT:DP", ["0/1:5", "1/1:6", "./.:7"], info="FL", alt="G")
# INFO/END sets the record length (an END below POS or not a whole number is refused by the reader) (END - POS + 1); SVLEN does not
gtdp = ["0/1:5", "1/1:6", "./.:7"]
for info, alt, ref in [
    ("END=5000", "<DEL>", "A"),
    ("END=%d" % (pos + 40), "<NON_REF>", "A"),
    ("END=%d" % (pos + 41), "<DEL>", "ACGT"),
    ("END=%d" % (pos + 42), "<DEL>", "ACGTACGT"),
    ("END=%d;SVLEN=-50" % (pos + 43), "<DEL>", "A"),
    ("SVLEN=-50", "<DEL>", "A"),
    ("SVLEN=-50;END=%d" % (pos + 45), "<DEL>", "A"),
    ("END=3000000000", "G", "A"),
]:
    rec("GT:DP", gtdp, info=info, alt=alt, ref=ref)
# INFO values htslib reads leniently or as missing
# a key with no `=` is `00` for bcftools, which prints it like a flag; gloci
# prints `NI1=.`, so these records are compared as bytes only
for info in ["NI1", "NF1", "NS1", "NS1="]:
    only.append(("GT:DP", gtdp, info, "G", "A"))
for info in ["NF1=", "NFD=1,,2", "NFD=", "NI1=1.5", "NI1=12abc", "NI1=abc", "NF1=2x", "NF1=abc"]:
    rec("GT:DP", gtdp, info=info)
# FORMAT `.`, and a column beyond the samples
rec(".", [".", ".", "."])
rec("GT:DP", ["0/1:5", "1/1:6", "./.:7", "0/0:8"])
# an empty String after the first key: the width counts only the key's own text
only.append(("GT:DP:S1", ["0/1:5:", "1/1:6:", "./.:7:"], ".", "G", "A"))
rec("GT:DP:S1", ["0/1:5", "1/1:6", "./.:7:"])
only.append(("GT:S1:DP", ["0/1::5", "1/1", "./.::7"], ".", "G", "A"))
rec("GT:DP:S1", ["0/1:5:", "1/1:6:x", "./.:7:"])
only.append(("GT:DP:S1:C1", ["0/1:5::", "1/1:6::", "./.:7::"], ".", "G", "A"))
rec("S1:DP", [":5", ":6", ":7"])
open("input.vcf", "w").write("\n".join(lines) + "\n")
print(pos, "records")

# Records compared as bytes only: bcftools prints them in ways gloci does not
# (a key with no value like a flag, an empty String as a space) and noodles'
# VCF writer wants GT first, so gloci cannot print the GT-later ones back.
lines = lines[:lines.index("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3") + 1]
pos = 0
for r in only:
    rec(*r)
rec("DP:GT", ["5", "6", "7"])
rec("DP:GT:S1", ["5", "6", "7"])
rec("DP:GT", ["5", "6", "7"], info="NI1=.")
rec("DP:GT:ID", ["5", "6", "7"])
rec("DP:S1:GT", ["5:a:0/1", "6:b:1/1", "7:c"])
open("bytes-only.vcf", "w").write("\n".join(lines) + "\n")
