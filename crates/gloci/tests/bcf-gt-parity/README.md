# BCF GT byte parity fixture

`expected.bcf` is bcftools 1.21 (htslib 1.21) converting `input.vcf`; the test
`bcf_gt_parity.rs` requires gloci's BCF records to be byte-identical to it.
bcftools 1.21 rejects a leading phase separator (`|1`, `/0|1`), so the input
has none; gloci refuses them (see `bcf_output_refuses_a_leading_phase_separator`).

Regenerate (from the repository root):

    chmod 777 crates/gloci/tests/bcf-gt-parity
    docker run --rm -u $(id -u) -v $PWD/crates/gloci/tests/bcf-gt-parity:/w -w /w gloci-diff:p0 \
        bcftools view -Ob -l0 input.vcf -o expected.bcf
    chmod 755 crates/gloci/tests/bcf-gt-parity
