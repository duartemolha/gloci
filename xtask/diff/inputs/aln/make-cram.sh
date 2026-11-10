#!/bin/sh
# Rebuilds the CRAM fixtures and their .crai indexes from the SAM sources and ../fa/ref.fa (D109,
# D119). Run it in the reference image, from the directory above this one:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work/aln gloci-diff:p0 sh make-cram.sh
# One thread and no @PG line, so the output is reproducible (htslib writes the output name as the
# CRAM file id, and the reference path /work/fa/ref.fa as each @SQ line's UR).
#
# gloci checks each @SQ length against the reference's .fai (D109), and cigars, flags and tags were
# written for longer references: their @SQ lengths are replaced by those of ref.fa (every read lies
# within ref.fa, so records and spans are those of the BAM fixtures). pastend has reads that run
# past the end of chr1 and md has reads whose MD and NM tags samtools drops when they are right.
#
# Variants: NAME.cram (CRAM 3.0), NAME.v31.cram (3.1), NAME.multi.cram (slices of two reads that
# may hold several references), NAME.embed.cram (the reference bases embedded in each slice) and
# NAME.noref.cram (no reference: every base stored). The .crai files are samtools' (gloci's must
# equal them).
set -eu
REF=../fa/ref.fa
for set in cigars flags tags pastend md; do
    awk -F'\t' -v OFS='\t' 'NR == FNR { len[$1] = $2; next }
        /^@SQ/ { for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) n = substr($i, 4)
                 for (i = 2; i <= NF; i++) if ($i ~ /^LN:/) $i = "LN:" len[n] }
        { print }' "$REF.fai" "$set.sam" > "tmp-$set.sam"
    cram() { out=$1; shift; samtools view -C -@ 0 --no-PG -o "$out" "$@" "tmp-$set.sam"; }
    cram "$set.cram" -T "$REF" --output-fmt-option version=3.0
    case "$set" in
        md) ;;
        pastend) cram "$set.embed.cram" -T "$REF" --output-fmt-option version=3.0 \
                     --output-fmt-option embed_ref=1 ;;
        *)
            cram "$set.v31.cram" -T "$REF" --output-fmt-option version=3.1
            cram "$set.multi.cram" -T "$REF" --output-fmt-option version=3.0 \
                --output-fmt-option seqs_per_slice=2 --output-fmt-option multi_seq_per_slice=1
            cram "$set.embed.cram" -T "$REF" --output-fmt-option version=3.0 \
                --output-fmt-option embed_ref=1
            cram "$set.noref.cram" --output-fmt-option version=3.0 --output-fmt-option no_ref=1
            ;;
    esac
    rm "tmp-$set.sam"
done
# cat-big.v31.cram: cat-big's 1,300 reads with qualities that drift along each read (so htslib
# codes them with order-1 rANS) and names for the tokeniser, as CRAM 3.1 without a reference: the
# codecs whose tables noodles keeps on the stack.
awk -F'\t' -v OFS='\t' 'BEGIN { srand(11) } /^@/ { print; next } {
    q = ""; v = 38
    for (j = 1; j <= length($10); j++) {
        v += int(rand() * 3) - 1 - (rand() < 0.08); if (v < 2) v = 2; if (v > 40) v = 40
        q = q sprintf("%c", 33 + v) }
    $11 = q; print }' cat-big.sam > tmp-big31.sam
samtools view -C -@ 0 --no-PG --output-fmt-option version=3.1 --output-fmt-option no_ref=1 \
    -o cat-big.v31.cram tmp-big31.sam
rm tmp-big31.sam
# contigs.noref.cram: chr1 of 10 Mb with 30 reads near its end, then 200 contigs of 50 kb with two
# reads each, without a reference, in slices of 100 reads that hold several references: the
# placeholders noodles needs there must stay bounded.
awk 'BEGIN { srand(13); OFS = "\t"; print "@HD", "VN:1.6", "SO:coordinate"
    print "@SQ", "SN:chr1", "LN:10000000"
    for (c = 0; c < 200; c++) print "@SQ", "SN:ctg" c, "LN:50000"
    q = ""; for (j = 0; j < 100; j++) q = q "I"
    for (i = 0; i < 30; i++) print "p" i, 0, "chr1", 9990001 + i * 100, 60, "100M", "*", 0, 0, bases(), q
    for (c = 0; c < 200; c++) for (j = 0; j < 2; j++)
        print "c" c "_" j, 0, "ctg" c, 1001 + j * 20000, 60, "100M", "*", 0, 0, bases(), q }
    function bases(  s, j) { s = ""; for (j = 0; j < 100; j++) s = s substr("ACGT", int(rand() * 4) + 1, 1); return s }' \
    > tmp-contigs.sam
samtools view -C -@ 0 --no-PG --output-fmt-option version=3.0 --output-fmt-option no_ref=1 \
    --output-fmt-option seqs_per_slice=100 --output-fmt-option multi_seq_per_slice=1 \
    -o contigs.noref.cram tmp-contigs.sam
rm tmp-contigs.sam
# Name-sorted CRAMs, in slices that hold several references, so every slice needs the same
# references: flags.namesorted.cram (with ref.fa) and pairs.namesorted.noref.cram (chr1 and chr2
# of 10 Mb, 150 reads each, without a reference). The decoders must share one copy of each.
awk -F'\t' -v OFS='\t' 'NR == FNR { len[$1] = $2; next }
    /^@SQ/ { for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) n = substr($i, 4)
             for (i = 2; i <= NF; i++) if ($i ~ /^LN:/) $i = "LN:" len[n] }
    { print }' "$REF.fai" flags.sam > tmp-flags.sam
samtools sort -n -@ 0 --no-PG -O sam -o tmp-flags-n.sam tmp-flags.sam
samtools view -C -@ 0 --no-PG -T "$REF" --output-fmt-option version=3.0 \
    --output-fmt-option seqs_per_slice=4 --output-fmt-option multi_seq_per_slice=1 \
    -o flags.namesorted.cram tmp-flags-n.sam
awk 'BEGIN { srand(17); OFS = "\t"; print "@HD", "VN:1.6", "SO:queryname"
    print "@SQ", "SN:chr1", "LN:10000000"; print "@SQ", "SN:chr2", "LN:10000000"
    q = ""; for (j = 0; j < 100; j++) q = q "I"
    for (i = 0; i < 300; i++) {
        s = ""; for (j = 0; j < 100; j++) s = s substr("ACGT", int(rand() * 4) + 1, 1)
        printf "q%03d\t0\tchr%d\t%d\t60\t100M\t*\t0\t0\t%s\t%s\n", i, i % 2 + 1, int(rand() * 9999000) + 1, s, q } }' \
    > tmp-pairs.sam
samtools view -C -@ 0 --no-PG --output-fmt-option version=3.0 --output-fmt-option no_ref=1 \
    --output-fmt-option seqs_per_slice=50 --output-fmt-option multi_seq_per_slice=1 \
    -o pairs.namesorted.noref.cram tmp-pairs.sam
rm tmp-flags.sam tmp-flags-n.sam tmp-pairs.sam
for f in *.cram; do
    case "$f" in
        flags-noeof.cram|unsorted.cram|overhang.cram|*.namesorted.*) ;;
        *) samtools index "$f" ;;
    esac
done
# unsorted.cram: chr1 reads at 501 and 502, then at 101 and 102, in slices of two: samtools index
# refuses it (not sorted by position), and so must gloci index.
printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n' > tmp-unsorted.sam
for p in 501 502 101 102; do
    printf 'u%s\t0\tchr1\t%s\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n' "$p" "$p" >> tmp-unsorted.sam
done
samtools view -C -@ 0 --no-PG -T "$REF" --output-fmt-option version=3.0 \
    --output-fmt-option seqs_per_slice=2 -o unsorted.cram tmp-unsorted.sam
rm tmp-unsorted.sam
# overhang.cram: a read whose deletion ends 1.5 Mb past the end of chr1 (samtools decodes it;
# gloci decodes reads at most 1 Mb, or four read lengths, past their reference's end), and one
# 0.5 Mb past it, which both decode.
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n' > tmp-overhang.sam
printf 'far\t0\tchr1\t995\t60\t5M1500000D5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n' >> tmp-overhang.sam
printf 'near\t0\tchr1\t996\t60\t5M500000D5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n' >> tmp-overhang.sam
samtools view -C -@ 0 --no-PG -T "$REF" --output-fmt-option version=3.0 -o overhang.cram \
    tmp-overhang.sam
rm tmp-overhang.sam
# flags-noeof: flags.cram without its 38-byte end-of-file container (a truncated copy).
head -c -38 flags.cram > flags-noeof.cram
