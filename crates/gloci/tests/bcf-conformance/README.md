# BCF conformance matrix

`input.vcf` holds FORMAT fields of every type (Integer, Float, String,
Character) and Number (1, A, R, G, `.`, 2), each missing in every sample, cut
short, mixed and present, at several FORMAT positions, plus empty subfields,
multi-allelic records and INFO values of every type that are `.`, partly `.`
or present. `bytes-only.vcf` has records compared as bytes only: GT after another key
(noodles' VCF writer wants GT first), and ones bcftools prints in ways gloci
does not (a key with no `=`, an empty String shown as a space).
`gen_input.py` writes both.

`expected.bcf` and `bytes-only.bcf` are bcftools 1.21 (htslib 1.21) converting
them; `expected.vcf` is `bcftools view -H expected.bcf`. The test
`bcf_conformance.rs` requires gloci's record bytes to equal them, gloci to
print `expected.bcf` as `expected.vcf`, and `gloci check` to pass on all of
them.

htslib rules the matrix pins (see `crates/gloci-io/src/bcf_fields.rs`): a
numeric field is as wide as the longest list; a field that is `.` or cut short
everywhere is one missing element; a String field is as wide as the longest
text any sample has for the key (a `.` standing for a missing or cut-short
value does not count) plus one NUL byte unless it is the first FORMAT key; a String or GT
field that every sample is cut short before is zero-width; an empty float is
0 (missing at the very end of the line) and an empty String is all NULs;
INFO/END sets the record length (END - POS + 1) when it is a single value of
at least POS, SVLEN never does; an INFO number is read by its leading digits
(`12abc` is 12, `abc` missing), an empty INFO float is missing, and an INFO key
with no `=` is a typed-missing value `00`; a FORMAT of `.` has no fields, and a
sample column beyond the header's samples is ignored.

Regenerate (from the repository root):

    python3 crates/gloci/tests/bcf-conformance/gen_input.py   # run inside the directory
    chmod 777 crates/gloci/tests/bcf-conformance
    docker run --rm -u $(id -u) -v $PWD/crates/gloci/tests/bcf-conformance:/w -w /w gloci-diff:p0 sh -c \
        'bcftools view -Ob -l0 input.vcf -o expected.bcf && bcftools view -H expected.bcf > expected.vcf \
         && bcftools view -Ob -l0 bytes-only.vcf -o bytes-only.bcf'
    chmod 755 crates/gloci/tests/bcf-conformance
