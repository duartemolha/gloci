#!/bin/sh
# Regenerate the vcftools 0.1.16 goldens here (needs the gloci-diff:p0 image).
# Usage: sh regen.sh   (from this directory)
# vcftools --diff-site-discordance et al. compare B (a.vcf, --diff b.vcf) by
# position, so gloci is run with --match position. vcftools prints a 0/0 ratio
# as `nan` on Arm and `-nan` on x86; every file is rewritten with `nan` so the
# goldens do not depend on the CPU that made them.
set -e
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
run() { # PREFIX A B [vcftools options]
  p=$1; a=$2; b=$3; shift 3
  for o in site-discordance indv-discordance discordance-matrix; do
    vcftools --vcf $a --diff $b --diff-$o "$@" --out $p >/dev/null 2>&1 || true
  done
}
run main a.vcf b.vcf
run ds diffsamples-a.vcf diffsamples-b.vcf
run dsmap diffsamples-a.vcf diffsamples-b.vcf --diff-indv-map map.txt
# half-missing genotypes: vcftools counts ./1 as called and stops on it in the
# matrix, so no matrix file is written (#compare-half-missing)
run half half-a.vcf half-b.vcf
rm -f *.log
sed -i "s/-nan/nan/g" *.diff.sites *.diff.indv
'
