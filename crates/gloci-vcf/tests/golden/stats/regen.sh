#!/bin/sh
# Regenerate the vcftools 0.1.16 goldens here (needs the gloci-diff:p0 image).
# Usage: sh regen.sh   (from this directory)
set -e
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" gloci-diff:p0 sh -c '
cd /work
for f in main mixed quality ploidy filters density hwegrid hwebig masked; do
  for o in freq:frq counts:frq.count site-depth:ldepth site-mean-depth:ldepth.mean \
           geno-depth:gdepth site-quality:lqual missing-site:lmiss hardy:hwe singletons:singletons \
           FILTER-summary:FILTER.summary hist-indel-len:indel.hist TsTv-summary:TsTv.summary; do
    vcftools --vcf $f.vcf --${o%%:*} --out $f >/dev/null 2>&1 || true
  done
done
for f in main density; do for b in 1000 7 2500; do vcftools --vcf $f.vcf --SNPdensity $b --out $f.b$b >/dev/null 2>&1; done; done
vcftools --vcf masked.vcf --minGQ 20 --out gq20 --freq >/dev/null 2>&1
for o in freq:frq counts:frq.count site-depth:ldepth site-mean-depth:ldepth.mean geno-depth:gdepth \
         missing-site:lmiss hardy:hwe singletons:singletons; do
  vcftools --vcf masked.vcf --minGQ 20 --${o%%:*} --out gq20 >/dev/null 2>&1
done
for o in site-mean-depth:ldepth.mean geno-depth:gdepth missing-site:lmiss; do
  vcftools --vcf masked.vcf --minDP 5 --${o%%:*} --out dp5 >/dev/null 2>&1
done
# vcftools stops at the polyploid genotype of ploidy.vcf for these
rm -f ploidy.frq ploidy.frq.count ploidy.hwe ploidy.lmiss ploidy.singletons
rm -f *.log
'
