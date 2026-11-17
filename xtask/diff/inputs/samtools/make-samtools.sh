#!/bin/sh
# Rebuilds the inputs of the samtools-essentials cases (faidx, dict and the rest) that are not
# hand-written (D119, D141). Run it in the reference image, from the directory above this one:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -w /work/samtools gloci-diff:p0 sh make-samtools.sh
# One thread, level 6, no @PG line, so the output is reproducible.
set -eu

# the inputs of the `gloci dict` cases that are not hand-written FASTA: a plain (non-BGZF)
# gzip of fa/ref.fa and an empty file. dup.fa (a repeated name) and alt.txt (an alternate-locus
# list: first tab-separated field is the name, "MT\r" ends in CRLF, "M extra" and "1 " match
# nothing) and nohdr.fa (sequence text before the first header) are written by hand.
gzip -n -c ../fa/ref.fa > ref-gzip.fa.gz
: > empty.fa

# empty.sam and empty.bam, a header and no records, for the cases whose input neither
# side reads (`flags`) and for the commands that must cope with a file of no reads. The header
# has two references and a read group; the BAM has its .bai.
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@RG\tID:g1\tSM:s1\n' > empty.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o empty.bam empty.sam
samtools index -@ 0 empty.bam

# unsorted-chroms.{sam,bam}, the reads of chr2 between two runs of chr1 reads and an
# unmapped read with no place at the end: `idxstats` must count by reference, not by run, where
# samtools' scan stops at "file is not position sorted". No @HD SO: the file claims no order. The
# BAM has no index (it could not be sorted into one). nometa/flags.bam is a copy of the alignment BAM
# (../aln/flags.bam); its .bai without metadata pseudo-bins is written by the `#[ignore]` test
# `make_nometa_index` (crates/gloci/tests/idxstats.rs), not by samtools.
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr3\tLN:3000\n' > unsorted-chroms.sam
{
  printf 'a1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'a2\t16\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'b1\t0\tchr2\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'b2\t4\tchr2\t30\t0\t*\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'a3\t0\tchr1\t50\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  printf 'a4\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n'
} >> unsorted-chroms.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o unsorted-chroms.bam unsorted-chroms.sam
mkdir -p nometa
cp ../aln/flags.bam nometa/flags.bam
# The copy is newer than its .bai, which would draw a stale-index warning; keep the index newer.
touch nometa/flags.bam.bai

# a CSI of flags.bam written by `samtools index -c`, in csi/ beside a copy of the BAM so
# the case can name both as aux files.
mkdir -p csi
cp ../aln/flags.bam csi/flags.bam
samtools index -c -@ 0 csi/flags.bam

# flagstat.{sam,bam}, built by the awk below and sorted by samtools. QC-passed: 320 records,
# 306 mapped (95.63% in samtools' f32 arithmetic, 95.62% in f64), 160 paired primaries of which 7
# are singletons (4.37% against 4.38%), 14 unmapped (7 mates of singletons, 7 single-end; placed
# and unplaced, one mate carrying 0x2), 60+26+20+40 mapped pairs (proper; proper duplicates; not
# proper; mate on another reference at MAPQ 4, MAPQ 5, and two with RNEXT *), 40 secondary (10
# duplicates, 10 flagged paired), 40 supplementary (10 duplicates), 8 single-end duplicates. Then
# 12 QC-failed (0x200) records, none paired, so the paired percentages there print N/A.
{
  printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n'
  awk 'function e(n, f, rn, p, mq, cg, rx, px) {
         printf "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t0\tACGTA\tIIIII\n", n, f, rn, p, mq, cg, rx, px
       }
       function pair(n, f1, f2, r1, r2, q1, q2) {
         p += 10
         e(n, f1, r1, p, q1, "5M", (r1 == r2 ? "=" : r2), p + 100)
         e(n, f2, r2, p + 100, q2, "5M", (r1 == r2 ? "=" : r1), p)
       }
       BEGIN {
         p = 100
         for (i = 1; i <= 30; i++) pair("pa" i, 99, 147, "chr1", "chr1", 60, 60)
         for (i = 1; i <= 13; i++) pair("pf" i, 1123, 1171, "chr1", "chr1", 60, 60)
         for (i = 1; i <= 10; i++) pair("pb" i, 65, 145, "chr2", "chr2", 60, 60)
         for (i = 1; i <= 9; i++) pair("pq4" i, 65, 129, "chr1", "chr2", 4, 4)
         for (i = 1; i <= 10; i++) pair("pq5" i, 65, 129, "chr1", "chr2", 5, 60)
         for (i = 1; i <= 2; i++) { p += 10; e("pstar" i, 65, "chr2", p, 60, "5M", "*", 0) }
         for (i = 1; i <= 7; i++) {
           p += 10
           e("s" i, 73, "chr1", p, 60, "5M", "=", p)
           if (i <= 4) e("s" i, 133, "chr1", p, 0, "*", "=", p)
           else if (i <= 6) e("s" i, 133, "*", 0, 0, "*", "chr1", p)
           else e("s" i, 135, "chr1", p, 0, "*", "=", p)
         }
         for (i = 1; i <= 73; i++) {
           p += 10
           f = (i % 2) ? 0 : 16
           if (i <= 8) f += 1024
           e("e" i, f, "chr2", p, 30, "5M", "*", 0)
         }
         for (i = 1; i <= 4; i++) { p += 10; e("u" i, 4, "chr2", p, 0, "*", "*", 0) }
         for (i = 5; i <= 7; i++) e("u" i, 4, "*", 0, 0, "*", "*", 0)
         for (i = 1; i <= 40; i++) {
           p += 10
           f = (i <= 20) ? 256 : (i <= 30) ? 1280 : 257
           e("sec" i, f, "chr1", p, 0, "5M", "*", 0)
         }
         for (i = 1; i <= 40; i++) {
           p += 10
           f = (i <= 20) ? 2048 : (i <= 30) ? 2064 : 3072
           e("sup" i, f, "chr2", p, 60, "5M", "*", 0)
         }
         split("512 512 512 512 528 528 516 516 768 2560 1536 1536", qf, " ")
         for (i = 1; i <= 12; i++) {
           p += 10
           if (i == 8) e("q" i, qf[i], "*", 0, 0, "*", "*", 0)
           else e("q" i, qf[i], "chr1", p, (qf[i] == 516 ? 0 : 20), (qf[i] == 516 ? "*" : "5M"), "*", 0)
         }
       }'
} > flagstat-unsorted.sam
samtools sort -@ 0 --no-PG -O sam -o flagstat.sam flagstat-unsorted.sam
rm flagstat-unsorted.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o flagstat.bam flagstat.sam
samtools index -@ 0 flagstat.bam

# the inputs of the `check` cases. noeof.bam is flagstat.bam without its 28-byte BGZF
# end-of-file block; truncated.bam is cut in the middle of a block. nosq.{sam,bam} are unaligned
# reads with a read group and no @SQ line. coord-unsorted.{sam,bam} claim SO:coordinate and have
# one read (record 7) before the one that precedes it; no case uses it (samtools has nothing to
# say about order), the tests of `gloci check` do. badfirst.bam (a header and an end-of-file
# block that are fine, a first record whose block_size is 20) is written by the `#[ignore]` test
# `make_broken_bams` (crates/gloci/tests/check_alignments.rs), not here.
head -c -28 flagstat.bam > noeof.bam
head -c 1500 flagstat.bam > truncated.bam
printf '@HD\tVN:1.6\n@RG\tID:g1\tSM:s1\n' > nosq.sam
{
  printf 'u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g1\n'
  printf 'u2\t4\t*\t0\t0\t*\t*\t0\t0\tGGGTA\tIIIII\tRG:Z:g1\n'
  printf 'u3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTA\tIIIII\tRG:Z:g1\n'
} >> nosq.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o nosq.bam nosq.sam
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n' > coord-unsorted.sam
{
  for p in 100 200 300 400 500 600; do
    printf 'r%s\t0\tchr1\t%s\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n' "$p" "$p"
  done
  printf 'rback\t0\tchr1\t150\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
  for p in 700 800; do
    printf 'r%s\t0\tchr1\t%s\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n' "$p" "$p"
  done
  printf 'q1\t0\tchr2\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n'
} >> coord-unsorted.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o coord-unsorted.bam coord-unsorted.sam
# calls.vcf.gz: a small BGZF VCF (bgzip's own output, with its EOF block), for the check case that
# samtools quickcheck rejects as "not sequence data" and gloci accepts.
bgzip -c ../cluster-in.vcf > calls.vcf.gz

# the inputs of the `head` and `samples` cases. rg.{sam,bam} have five @RG lines: samples
# NA12878, NA12891, NA12878 (again) and HG002, a LB on two (libA, libB) and a fifth with a LB (libC)
# and no SM, so samtools' hash order (HG002, NA12878, NA12891) is not the header order and the
# tag cases tell the samples apart. nosm.sam has @RG lines with no SM at all. ref.fa matches the
# @SQ lines of rg.sam; ref2.fa is a copy of it (two FASTAs that both match); other.fa has a
# chr2 one base short (it matches nothing). The notext.bam of the head case is aln/notext.bam.
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@RG\tID:g1\tSM:NA12878\tLB:libA\n@RG\tID:g2\tSM:NA12891\n@RG\tID:g3\tSM:NA12878\tLB:libB\n@RG\tID:g4\tSM:HG002\n@RG\tID:g5\tLB:libC\n' > rg.sam
printf 'r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g1\nr2\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g2\nr3\t0\tchr2\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g4\n' >> rg.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o rg.bam rg.sam
samtools index -@ 0 rg.bam
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:g5\tLB:libC\n@RG\tID:g6\tPL:x\n' > nosm.sam
awk 'BEGIN{srand(3);for(c=1;c<=2;c++){print ">chr" c; n=(c==1?1000:2000); s=""; for(i=0;i<n;i++){s=s substr("ACGT",int(rand()*4)+1,1)}; for(i=1;i<=n;i+=60) print substr(s,i,60)}}' > ref.fa
samtools faidx ref.fa
cp ref.fa ref2.fa
cp ref.fa.fai ref2.fa.fai
awk 'BEGIN{print ">chr1"; s=""; for(i=0;i<1000;i++)s=s "A"; print s; print ">chr2"; s=""; for(i=0;i<999;i++)s=s "A"; print s}' > other.fa
samtools faidx other.fa

# BGZF and BCF copies of ../hdr-distinct.vcf (its header has ##FILTER=PASS,
# which bcftools head would otherwise add), for the head cases on compressed input.
bgzip -c ../hdr-distinct.vcf > hdr-pass.vcf.gz
bcftools view -O b -o hdr-pass.bcf ../hdr-distinct.vcf
# the `depth` fixtures, all written by hand. depth.sam sets the traps: a deletion whose
# next base has quality 10, `N`, `S`, `I`, `=`/`X`, a CIGAR `*` (samtools reads a mapped read with
# one as unmapped), QUAL `*`, MAPQ 29 and 30, duplicate, secondary, QC-fail and supplementary
# reads, a placed unmapped read, a read ending at chr2's LN, reads of 30, 49 and 50 bases once soft
# clips are cut, and an empty chr3. depth2.sam has other reads on the same references;
# depth-reordered.sam lists them in another order; depth-unsorted.sam has a read 3 bp before its
# predecessor that ends after samtools' last row. depth-noseq.sam has reads with SEQ `*`, each with an
# XB:B:C aux array as long as its CIGAR: samtools -q reads the bytes after the CIGAR as qualities,
# which are then the array's own, so its output is stable (gloci passes every base). depth-pastend.sam
# has reads that run past chr2's LN (samtools prints rows beyond it, so no other fixture has one).
# depth.cram is
# written without a reference (every base stored), so reading it needs none. depth-idx/input.bai
# is the index under the name samtools and gloci find beside a file called input.bam.
for s in depth depth2 depth-reordered depth-unsorted depth-noseq depth-pastend; do
    samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o "$s.bam" "$s.sam"
done
samtools index -@ 0 depth.bam
samtools index -@ 0 depth2.bam
samtools index -@ 0 depth-noseq.bam
samtools index -@ 0 depth-pastend.bam
mkdir -p depth-idx
samtools index -@ 0 depth.bam depth-idx/input.bai
samtools view -C -@ 0 --no-PG --output-fmt-option version=3.0 --output-fmt-option no_ref=1 -o depth.cram depth.sam
samtools index -@ 0 depth.cram

# the `depth --skip-mate-overlap` fixture, written by hand (mates.sam). Templates: a overlaps
# its mate in part; b has both mates at one position, the longer second; c's second mate lies inside
# the first; d does not overlap; e's second mate has a deletion before the first mate's end; f has
# the mate-unmapped flag clear and no mate position; g is a resolved pair followed by two
# supplementary alignments of one name inside the second mate (samtools keeps a stale entry for
# them, gloci does not: #depth-mate-overlap); h's second mate has quality 5 on the bases it shares;
# j's mate is on chr2; k is a pair on chr2; n's first mate has MAPQ 5, so under -Q 10 it neither waits nor clips.
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o mates.bam mates.sam
samtools index -@ 0 mates.bam

# the `bedcov` fixtures bedcov.bed and bedcov-badchrom.bed are written by hand (see the
# comment of the bedcov cases). bedcov-badchrom.bed is bedcov.bed with a chrZ line inserted as
# line 9. Nothing to rebuild: the BAMs are depth.bam, depth2.bam and depth-pastend.bam above.
# the inputs of the `gloci reheader` cases. reheader-in.sam is hand-written (three
# references of different lengths, so a swapped order shows in the LNs, and one unplaced read);
# new-header.sam renames its references to 1, 2, 3, appends a fourth @SQ and changes the @RG;
# bad-ln.sam changes chr2's LN and bad-order.sam swaps chr2 and chr3 (both are hand-written).
# new-header.bam is new-header.sam as a BAM header file. new-header-cram.sam is the header
# for aln/flags.embed.cram (chr1 1000, chr2 1234, with its M5s), renamed, with a third @SQ added.
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o reheader-in.bam reheader-in.sam
samtools index -@ 0 reheader-in.bam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o new-header.bam new-header.sam

# view-targets.bed, the `view --regions-file` / `-L` regions, is hand-written; the reads are
# aln/flags.bam's. chr1 100 100 and chr1 110 110 sit on the edges of the 10-base reads at
# chr1:101 (so a zero-length region must match only reads that strictly contain it), chr1 105 105
# is inside them, chr1 250 260 and chr1 255 258 overlap each other and r_mapq10 (a read must
# count once), chr1 700 701 is the placed unmapped read, and chr2 55 56 is inside both chr2 reads.
# Nothing here needs building.

# mosdepth.sam: overlapping proper pairs for mosdepth's mate rule, and references with no reads or
# only filtered reads (the depth --summary and --out-prefix cases). chr1: p partial overlap; x MAPQ 0;
# q both mates at 201, the second longer; r a deletion inside the overlap; s an intron inside it;
# t soft clips and an insertion; u an overlapping pair without PROPER_PAIR; v a pair and a
# supplementary (2147); w a proper pair with its mate on chr2; z the second mate MAPQ 5; y1-y4 dup,
# secondary, qcfail, placed unmapped. chr2: a duplicate and a placed unmapped read only. chr3: none.
# chr4: 320 one-base reads at 5001 and two reads overlapping by 10 (distribution skips above 300).
{
  printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:300\n@SQ\tSN:chr3\tLN:200\n@SQ\tSN:chr4\tLN:20000\n'
  sed 's/ \{1,\}/\t/g' <<'EOF'
p 99 chr1 101 60 50M = 131 80 * *
x 0 chr1 110 0 30M * 0 0 * *
p 147 chr1 131 60 50M = 101 -80 * *
q 99 chr1 201 60 40M = 201 60 * *
q 147 chr1 201 60 60M = 201 -60 * *
r 99 chr1 301 60 20M10D20M = 321 60 * *
r 147 chr1 321 60 40M = 301 -60 * *
s 99 chr1 401 60 10M50N10M = 411 70 * *
s 147 chr1 411 60 60M = 401 -70 * *
t 99 chr1 501 60 5S20M2I20M = 521 50 * *
t 147 chr1 521 60 30M5S = 501 -50 * *
u 97 chr1 601 60 30M = 611 40 * *
u 145 chr1 611 60 30M = 601 -40 * *
v 99 chr1 701 60 50M = 721 70 * *
v 147 chr1 721 60 50M = 701 -70 * *
v 2147 chr1 731 60 20M = 701 0 * *
w 99 chr1 801 60 40M chr2 101 0 * *
z 99 chr1 851 60 40M = 871 60 * *
z 147 chr1 871 5 40M = 851 -60 * *
y1 1024 chr1 901 60 20M * 0 0 * *
y2 256 chr1 905 60 20M * 0 0 * *
y3 512 chr1 910 60 20M * 0 0 * *
y4 4 chr1 920 0 * * 0 0 * *
d1 1024 chr2 50 60 20M * 0 0 * *
d2 4 chr2 80 0 * * 0 0 * *
EOF
  for i in $(seq 1 320); do printf 'deep%03d\t0\tchr4\t5001\t60\t1M\t*\t0\t0\t*\t*\n' "$i"; done
  printf 'w1\t0\tchr4\t6001\t60\t20M\t*\t0\t0\t*\t*\nw2\t0\tchr4\t6011\t60\t20M\t*\t0\t0\t*\t*\n'
} > mosdepth.sam
samtools view -b -@ 0 --output-fmt-option level=6 --no-PG -o mosdepth.bam mosdepth.sam
samtools index -@ 0 mosdepth.bam
