#!/bin/sh
# Writes the real-data liftover inputs from the Ensembl GRCh37-to-GRCh38 slice
# (D174 R12). Deterministic: no randomness. Run from any directory:
#   sh xtask/diff/inputs/liftover/make-real-inputs.sh
set -eu
dir=$(dirname "$0")
chain="$dir/ensembl-grch37-to-grch38.21-22.chain"

# Walks chains whose source is 21 or 22; k numbers blocks across the file.
# Emits "kind<TAB>chrom<TAB>start<TAB>end<TAB>name<TAB>strand<TAB>k" lines for BED6
# and "B12<TAB>..." lines for BED12.
awk '
function flush(   i, n, kk, ok, st, ce, sizes, starts, cs) {
    if (nb == 0) return
    for (i = 1; i <= nb; i++) {
        kk = k0 + i - 1
        s = bs[i]; e = s + bz[i]
        if (bz[i] >= 300 && kk % 5 == 0)
            printf "BED\t%s\t%d\t%d\tin%d\t%s\t%d\n", chrom, s + 100, s + 200, kk, (kk % 10 == 0 ? "+" : "-"), kk
        if (kk % 13 == 0 && bz[i] >= 60)
            printf "BED\t%s\t%d\t%d\tpt%d\t+\t%d\n", chrom, s + 50, s + 50, kk, kk
        if (i < nb && kk % 7 == 0)
            printf "BED\t%s\t%d\t%d\tgap%d\t+\t%d\n", chrom, e - 40, bs[i + 1] + 40, kk, kk
        if (kk % 9 == 0 && i + 2 <= nb && bz[i] >= 120 && bz[i + 1] >= 120 && bz[i + 2] >= 120) {
            cs = bs[i] + 10
            ce = bs[i + 2] + 110
            printf "B12\t%s\t%d\t%d\ttx%d\t%s\t%d\t%d\t3\t100,100,100,\t0,%d,%d,\n", chrom, cs, ce, kk, ((kk / 9) % 2 == 0 ? "+" : "-"), cs + 5, ce - 5, bs[i + 1] - bs[i], bs[i + 2] - bs[i]
        }
    }
    printf "BED\t%s\t%d\t%d\toff%d\t+\t%d\n", chrom, tend, tend + 500, k0, k0
    k0 += nb
    nb = 0
}
/^#/ { next }
/^chain/ {
    flush()
    inchain = ($3 == "21" || $3 == "22")
    if (inchain) { chrom = $3; tend = $7; pos = $6 }
    next
}
NF == 0 { next }
inchain {
    nb++; bs[nb] = pos; bz[nb] = $1
    if (NF == 3) pos += $1 + $2
    else pos += $1
}
END { flush() }
BEGIN { k0 = 0; nb = 0 }
' "$chain" > "$dir/.real.tmp"

awk -F'\t' '$1=="BED"{print $2"\t"$3"\t"$4"\t"$5"\t0\t"$6}' "$dir/.real.tmp" > "$dir/grch37-21-22.bed"
awk -F'\t' '$1=="BED"{print "chr"$2"\t"$3"\t"$4"\t"$5"\t0\t"$6}' "$dir/.real.tmp" > "$dir/grch37-chr21-22.bed"
awk -F'\t' '$1=="B12"{print $2"\t"$3"\t"$4"\t"$5"\t0\t"$6"\t"$7"\t"$8"\t0\t"$9"\t"$10"\t"$11}' "$dir/.real.tmp" > "$dir/grch37-21-22-bed12.bed"
awk -F'\t' '$1=="BED" && $5 !~ /^pt/ {printf "%s\tensembl\texon\t%d\t%d\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s.1\";\n", $2, $3 + 1, $4, $6, $5, $5}' "$dir/.real.tmp" > "$dir/grch37-21-22.gtf"
awk -F'\t' '$1=="BED"{print $2"\t"$3"\t"$4"\t"$7}' "$dir/.real.tmp" > "$dir/grch37-21-22.bedGraph"
rm -f "$dir/.real.tmp"
