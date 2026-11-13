# BCF as VCF text, as bcftools writes it

`input.bcf` is `input.vcf` written by bcftools 1.21, and `expected.txt` is what
bcftools prints for it:

```sh
bcftools view --no-version -Ob input.vcf -o input.bcf
bcftools view --no-version -H input.bcf > expected.txt
```

`gloci view -O vcf --no-header input.bcf` must give the same bytes
(`tests/bcf_text_parity.rs`). The records are the awkward ones: a comma, `=`
and `%` inside strings, numbers on both sides of every change of notation,
short and missing vectors, every integer width, mixed ploidy and phasing, a
genotype that needs 16 bits, strings of different widths, a sample cut short
before its GT (which htslib prints as `-65`), and breakend and symbolic
alleles.

To regenerate after editing `input.vcf`, run the two commands above with
bcftools 1.21 (the `gloci-diff:p0` image has it).
