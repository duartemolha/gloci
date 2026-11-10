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
for f in *.cram; do samtools index "$f"; done
# flags-noeof: flags.cram without its 38-byte end-of-file container (a truncated copy).
head -c -38 flags.cram > flags-noeof.cram
