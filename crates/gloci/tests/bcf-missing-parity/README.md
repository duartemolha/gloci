# BCF String FORMAT / FILTER byte parity fixture

`expected.bcf` is bcftools 1.21 (htslib 1.21) converting `input.vcf`; the test
`bcf_missing_parity.rs` requires gloci's BCF records to be byte-identical to it.
The records cover a String FORMAT field that is `.` in every sample (first or
later key, with and without GT) and declared FILTER ids.

htslib writes a String FORMAT field with the width of its longest value for the
first FORMAT key, and one byte more (a NUL after each value) for every later
key; a field that is `.` everywhere is a lone `.` (first key) or `.` and NUL.
Widths of 15 and up use the escaped length form (the input has values of 14,
15, 16, 127 and 130 characters).

bcftools 1.21 refuses a FILTER that the header does not declare
(`Tag not defined in header`); gloci refuses it too, naming the line and the
`##FILTER` line to add (see `bcf_output_refuses_an_undeclared_filter`).

Regenerate (from the repository root):

    chmod 777 crates/gloci/tests/bcf-missing-parity
    docker run --rm -u $(id -u) -v $PWD/crates/gloci/tests/bcf-missing-parity:/w -w /w gloci-diff:p0 \
        bcftools view -Ob -l0 input.vcf -o expected.bcf
    chmod 755 crates/gloci/tests/bcf-missing-parity
