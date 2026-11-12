#!/bin/sh
# Rebuilds the mates.*.cram fixtures (flagstat's partial CRAM decoding, Extra X1). make-cram.sh
# runs it last; it runs alone the same way, in the reference image, from the directory above:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work/aln gloci-diff:p0 sh make-mates.sh
# One thread, no @PG line and fixed seeds, so the output is reproducible.
#
# 80 templates on chr1, chr2 and chr10 of ../fa/ref.fa, 20-base reads: pairs whose mates lie close
# (htslib attaches them within a slice: the mate is a distance, NF), far apart or across slices,
# and on other references (detached: MF, NS, NP, TS); pairs with an unmapped mate (fixmate places
# it at its mate) and with both unmapped (no place); secondary and supplementary reads of paired
# templates, duplicates, QC-failed pairs, low-MAPQ mates on other references, and single reads.
# samtools fixmate fills in the mate fields, TLEN and the proper-pair bit, so htslib attaches every
# pair it can. Variants: mates.cram (CRAM 3.0, slices of 20 reads, two to a container),
# mates.v31.cram (3.1), mates.multi.cram (slices that hold several references: RI),
# mates.noref.cram (no reference), mates.namesorted.cram (name order: mates side by side, slices
# of several references, positions in the core block) and mates.raw.cram (3.0 with every block
# stored raw, level 0: the fallback tests rewrite its encodings).
set -eu
REF=../fa/ref.fa
awk 'BEGIN { srand(43); OFS = "\t"
    print "@HD", "VN:1.6", "SO:unsorted"
    print "@SQ", "SN:chr1", "LN:1000"; print "@SQ", "SN:chr2", "LN:1234"
    print "@SQ", "SN:chr10", "LN:600"
    split("chr1 chr2 chr10", ref, " "); split("1000 1234 600", len, " ")
    q = "IIIIIIIIIIIIIIIIIIII"
    for (t = 0; t < 80; t++) {
        kind = t % 8; c = 1 + int(rand() * 3); p = 1 + int(rand() * (len[c] - 300))
        mq1 = 20 + int(rand() * 40); mq2 = 20 + int(rand() * 40)
        extra = (t % 11 == 5) ? 512 : (t % 13 == 7) ? 1024 : 0
        proper = (kind < 2 || kind == 3 || kind == 7) ? 2 : 0
        f1 = 65 + proper + extra; f2 = 129 + 16 + proper + extra
        c2 = c; p2 = p + 30 + int(rand() * 100); u1 = 0; u2 = 0
        if (kind == 3) p2 = p + 200 + int(rand() * 60)
        if (kind == 4) { c2 = 1 + (c % 3); p2 = 1 + int(rand() * (len[c2] - 40)); mq1 = int(rand() * 10) }
        if (kind == 2) f2 -= 16
        if (kind == 5) { u2 = 1; f2 = 129 + 4 + extra }
        if (kind == 6) { u1 = 1; u2 = 1; f1 = 65 + 4 + extra; f2 = 129 + 4 + extra }
        read("t" t, f1, c, p, mq1, u1)
        read("t" t, f2, c2, p2, mq2, u2)
        if (kind == 7) {
            read("t" t, f1 + 256, c, 1 + int(rand() * (len[c] - 40)), 3, 0)
            read("t" t, f2 + 2048, 1 + (c % 3), 1 + int(rand() * 500), 15, 0)
        }
    }
    for (i = 0; i < 15; i++) { c = 1 + int(rand() * 3)
        read("s" i, (i % 2) * 16, c, 1 + int(rand() * (len[c] - 40)), int(rand() * 60), 0) }
    for (i = 0; i < 3; i++) read("u" i, 4, 0, 0, 0, 1)
}
function read(name, flag, c, pos, mq, unmapped,   s, j) {
    s = ""; for (j = 0; j < 20; j++) s = s substr("ACGT", int(rand() * 4) + 1, 1)
    if (unmapped) print name, flag, "*", 0, 0, "*", "*", 0, 0, s, q
    else print name, flag, ref[c], pos, mq, "20M", "*", 0, 0, s, q
}' > tmp-mates.sam
samtools sort -n -@ 0 --no-PG -O sam -o tmp-mates-n.sam tmp-mates.sam
samtools fixmate -@ 0 --no-PG -O sam tmp-mates-n.sam tmp-mates-f.sam
samtools sort -@ 0 --no-PG -O sam -o tmp-mates-c.sam tmp-mates-f.sam
cram() { out=$1; in=$2; shift 2; samtools view -C -@ 0 --no-PG -o "$out" "$@" "$in"; }
cram mates.cram tmp-mates-c.sam -T "$REF" --output-fmt-option version=3.0 \
    --output-fmt-option seqs_per_slice=20 --output-fmt-option slices_per_container=2
cram mates.v31.cram tmp-mates-c.sam -T "$REF" --output-fmt-option version=3.1 \
    --output-fmt-option seqs_per_slice=20 --output-fmt-option slices_per_container=2
cram mates.multi.cram tmp-mates-c.sam -T "$REF" --output-fmt-option version=3.0 \
    --output-fmt-option seqs_per_slice=20 --output-fmt-option multi_seq_per_slice=1
cram mates.noref.cram tmp-mates-c.sam --output-fmt-option version=3.0 \
    --output-fmt-option no_ref=1 --output-fmt-option seqs_per_slice=20
cram mates.namesorted.cram tmp-mates-f.sam -T "$REF" --output-fmt-option version=3.0 \
    --output-fmt-option seqs_per_slice=20 --output-fmt-option multi_seq_per_slice=1
cram mates.raw.cram tmp-mates-c.sam -T "$REF" --output-fmt-option version=3.0 \
    --output-fmt-option seqs_per_slice=20 --output-fmt-option level=0
rm tmp-mates.sam tmp-mates-n.sam tmp-mates-f.sam tmp-mates-c.sam
