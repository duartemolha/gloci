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
