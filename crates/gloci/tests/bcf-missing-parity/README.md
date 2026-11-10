# BCF missing String FORMAT / FILTER byte parity fixture

`expected.bcf` is bcftools 1.21 (htslib 1.21) converting `input.vcf`; the test
`bcf_missing_parity.rs` requires gloci's BCF records to be byte-identical to it.
The records cover a String FORMAT field that is `.` in every sample (first or
later key, with and without GT) and declared FILTER ids.

htslib writes such a field with width 1 for the first FORMAT key and width 2
(`.` and NUL) for any later key. For later keys it also pads non-missing
String values by one byte; gloci does not (the decoded values are equal), so
this fixture holds no such record.

bcftools 1.21 refuses a FILTER that the header does not declare
(`Tag not defined in header`); gloci refuses it too, naming the line and the
`##FILTER` line to add (see `bcf_output_refuses_an_undeclared_filter`).

Regenerate (from the repository root):

    chmod 777 crates/gloci/tests/bcf-missing-parity
    docker run --rm -u $(id -u) -v $PWD/crates/gloci/tests/bcf-missing-parity:/w -w /w gloci-diff:p0 \
        bcftools view -Ob -l0 input.vcf -o expected.bcf
    chmod 755 crates/gloci/tests/bcf-missing-parity
