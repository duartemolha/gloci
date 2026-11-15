#!/bin/sh
# Regenerate the vcftools 0.1.16 goldens here (needs the gloci-diff:p0 image).
# Usage: sh regen.sh   (from this directory)
# vcftools writes one statistic per run. A run's files are PREFIX.<ext>, with
# PREFIX = popgen[.MASK][.STAT]: MASK is dp5 (--minDP 5) or indv (--indv s2
# --indv s3), STAT is w10 (--window-pi 10), w10s5 and w10s3 (with
# --window-pi-step 5 or 3), t49 and t100 (--TajimaD 49 or 100), none for
# --site-pi. vcftools prints a 0/0 site pi as `nan` on Arm and `-nan` on x86;
# every file is rewritten with `nan` so the goldens do not depend on the CPU
# that made them.
set -e
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
run() { # PREFIX [vcftools options]
  p=$1; shift
  v() { vcftools --vcf popgen.vcf "$@" >/dev/null 2>&1; }
  v --site-pi "$@" --out $p
  v --window-pi 10 "$@" --out $p.w10
  v --window-pi 10 --window-pi-step 5 "$@" --out $p.w10s5
  v --window-pi 10 --window-pi-step 3 "$@" --out $p.w10s3
  v --TajimaD 49 "$@" --out $p.t49
  v --TajimaD 100 "$@" --out $p.t100
}
run popgen
run popgen.dp5 --minDP 5
run popgen.indv --indv s2 --indv s3
rm -f *.log
sed -i "s/-nan/nan/g" *.sites.pi *.windowed.pi *.Tajima.D
'

# Weir and Cockerham Fst (--weir-fst-pop) on fst.vcf and missingpop.vcf, with
# the populations of pop1.txt and pop2.txt (and pop3.txt for p3). A run's
# files are PREFIX.weir.fst, or PREFIX.windowed.weir.fst with a window, with
# PREFIX = fst[.p3][.dp5][.w100|.w100s50] (--minDP 5; --fst-window-size 100,
# with --fst-window-step 50) or missingpop[.w100]. vcftools writes no summary
# file: it prints the mean and the weighted estimate in its log (on stderr),
# kept here as PREFIX.weir.fst.log for the check of gloci's
# PREFIX.weir.fst.summary.
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
fst() { # VCF PREFIX [vcftools options]
  f=$1; p=$2; shift 2
  vcftools --vcf $f "$@" --out $p >$p.err 2>&1
  grep "Fst estimate:" $p.err > $p.weir.fst.log
  rm -f $p.err
}
P2="--weir-fst-pop pop1.txt --weir-fst-pop pop2.txt"
P3="$P2 --weir-fst-pop pop3.txt"
fst fst.vcf fst $P2
fst fst.vcf fst.p3 $P3
fst fst.vcf fst.dp5 $P2 --minDP 5
fst fst.vcf fst.w100 $P2 --fst-window-size 100
fst fst.vcf fst.w100s50 $P2 --fst-window-size 100 --fst-window-step 50
fst fst.vcf fst.p3.w100s50 $P3 --fst-window-size 100 --fst-window-step 50
fst fst.vcf fst.dp5.w100s50 $P2 --minDP 5 --fst-window-size 100 --fst-window-step 50
fst missingpop.vcf missingpop $P2
fst missingpop.vcf missingpop.w100 $P2 --fst-window-size 100
sed -i "s/-nan/nan/g" *.weir.fst *.weir.fst.log
'

# Linkage disequilibrium on ld.vcf (unphased: --geno-r2, --interchrom-geno-r2)
# and phased.vcf (phased but for chr1:40, chr1:60 and chr1:80: --hap-r2,
# --interchrom-hap-r2, and --geno-r2, which keeps those sites). A run's file
# is PREFIX.geno.ld, PREFIX.hap.ld or PREFIX.interchrom.{geno,hap}.ld, with
# PREFIX = ld|phased[.TAG]: w2 and w1 (--ld-window 2 or 1), bp100
# (--ld-window-bp 100), min2 (--ld-window-min 2), bpmin10 (--ld-window-bp-min
# 10), mix (--ld-window 3 --ld-window-bp 150 --ld-window-min 2
# --ld-window-bp-min 5), r2 (--min-r2 0.2), dp5 (--minDP 5), indv
# (--remove-indv s1 --remove-indv s5). vcftools prints 0/0 as -nan on x86;
# the files are rewritten with nan.
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
r2() { # VCF PREFIX [vcftools options]
  f=$1; p=$2; shift 2
  vcftools --vcf $f "$@" --out $p >/dev/null 2>&1
}
r2 ld.vcf ld --geno-r2
r2 ld.vcf ld.w2 --geno-r2 --ld-window 2
r2 ld.vcf ld.bp100 --geno-r2 --ld-window-bp 100
r2 ld.vcf ld.min2 --geno-r2 --ld-window-min 2
r2 ld.vcf ld.bpmin10 --geno-r2 --ld-window-bp-min 10
r2 ld.vcf ld.mix --geno-r2 --ld-window 3 --ld-window-bp 150 --ld-window-min 2 --ld-window-bp-min 5
r2 ld.vcf ld.r2 --geno-r2 --min-r2 0.2
r2 ld.vcf ld.dp5 --geno-r2 --minDP 5
r2 ld.vcf ld.indv --geno-r2 --remove-indv s1 --remove-indv s5
r2 ld.vcf ld --interchrom-geno-r2
r2 ld.vcf ld.r2 --interchrom-geno-r2 --min-r2 0.2
r2 phased.vcf phased --hap-r2
r2 phased.vcf phased --geno-r2
r2 phased.vcf phased.w1 --hap-r2 --ld-window 1
r2 phased.vcf phased.dp5 --hap-r2 --minDP 5
r2 phased.vcf phased.r2 --hap-r2 --min-r2 0.2
r2 phased.vcf phased --interchrom-hap-r2
sed -i "s/-nan/nan/g" *.geno.ld *.hap.ld
'

# A sample without GT on nogt.vcf (phased otherwise): vcftools reads it as
# ./., unphased, so --hap-r2 leaves out chr1:20, whose FORMAT has no GT, and
# chr1:40, where n2's column has none. PREFIX.hap.ld with PREFIX = nogt or
# nogt.w1 (--ld-window 1).
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
vcftools --vcf nogt.vcf --hap-r2 --out nogt >/dev/null 2>&1
vcftools --vcf nogt.vcf --hap-r2 --ld-window 1 --out nogt.w1 >/dev/null 2>&1
sed -i "s/-nan/nan/g" nogt.hap.ld nogt.w1.hap.ld
'

# Relatedness (--relatedness, --relatedness2) on related.vcf: dad, mom, their
# child kid, unrelated u1, u2 and u3, h1 and h2 (homozygous everywhere, so
# --relatedness2 divides by 0) and gap (called at four sites, none of them
# u3's); and on halfcall.vcf, which holds ./1 and 1/. genotypes. A run's file
# is PREFIX.relatedness or PREFIX.relatedness2, with PREFIX = related,
# related.dp5 (--minDP 5), related.indv (--indv u2 --indv kid --indv dad) or
# halfcall. related.dp5.relatedness2 and halfcall.relatedness differ from
# gloci (docs/differences.md#relatedness-genotypes): vcftools counts masked
# genotypes in --relatedness2 and reads ./1 and 1/. as 0 in --relatedness.
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
rel() { # VCF PREFIX [vcftools options]
  f=$1; p=$2; shift 2
  vcftools --vcf $f --relatedness "$@" --out $p >/dev/null 2>&1
  vcftools --vcf $f --relatedness2 "$@" --out $p >/dev/null 2>&1
}
rel related.vcf related
rel related.vcf related.dp5 --minDP 5
rel related.vcf related.indv --indv u2 --indv kid --indv dad
rel halfcall.vcf halfcall
rm -f related.log related.dp5.log related.indv.log halfcall.log
sed -i "s/-nan/nan/g" *.relatedness *.relatedness2
'
