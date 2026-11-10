#!/bin/sh
# Regenerates the FASTA fixtures in this directory, deterministically, and the
# samtools 1.21 / htslib 1.21 indexes that gloci's .fai/.gzi must equal byte for
# byte. Run inside the reference image (gloci-diff:p0):
#
#   docker run --rm --user "$(id -u)" -v "$PWD":/w -w /w gloci-diff:p0 sh make-fasta.sh
#
# Contig order is chrM, chr1, chr2, chr10, chrX on purpose: lexicographic,
# karyotype and file order all differ.
set -eu

# gen NAME LEN DESC  -> a FASTA record body on stdout (60-column lines).
gen() {
    perl -e '
        my ($name, $len, $desc) = @ARGV;
        my $s = 12345 + length($name) * 7919;
        my $seq = "";
        for (1 .. $len) {
            $s = ($s * 1103515245 + 12345) % 2147483648;
            $seq .= substr("ACGT", ($s >> 16) & 3, 1);
        }
        if ($name eq "chr1") {
            # chr1:100-129 (0-based, half-open): 27 A/T and 2 G/C of 29 bases.
            substr($seq, 100, 29) = substr("ATATATATATATGATATATATATATCATAT", 0, 29);
            substr($seq, 200, 60) = lc(substr($seq, 200, 60));       # soft-masked run
            substr($seq, 300, 40) = "N" x 40;                        # N run
        }
        if ($name eq "chrX") { substr($seq, 50, 20) = "RYKMSWBDHVNNRYKMSWBD"; }
        if ($name eq "chr10") { substr($seq, 10, 12) = "ACGUUACGUUAC"; substr($seq, 30, 6) = "nnacgt"; }
        if ($name eq "chrM") { substr($seq, 0, 4) = "acgt"; }
        print ">$name", ($desc eq "-" ? "" : " $desc"), "\n";
        print substr($seq, $_, 60), "\n" for map { $_ * 60 } 0 .. int(($len - 1) / 60);
    ' "$@"
}

# ref.fa: five contigs; chr10 is exactly ten full lines, the rest end short.
{
    gen chrM 150 -
    gen chr1 1000 -
    gen chr2 1234 -
    gen chr10 600 -
    gen chrX 777 -
} > ref.fa

# ref.fa.gz: the same contigs with a 300 kb chr2, BGZF over several blocks.
{
    gen chrM 150 -
    gen chr1 1000 -
    gen chr2 300000 -
    gen chr10 600 -
    gen chrX 777 -
} > big.fa.tmp
bgzip -c big.fa.tmp > ref.fa.gz
rm big.fa.tmp

# crlf.fa: ref.fa with CRLF line ends.
sed 's/$/\r/' ref.fa > crlf.fa

# desc.fa: descriptions after the name (spaces and a tab); the index keeps the first word.
{
    gen chrM 150 "mitochondrion circular"
    gen chr1 1000 "Homo sapiens	chromosome 1"
    gen chr2 1234 "x"
} > desc.fa

# empty-seq.fa: empty sequences (first, middle) are left out of the index by samtools.
{
    printf '>empty1\n>chr1 first\n'
    gen chrM 150 - | tail -n +2
    printf '>empty2 nothing here\n>chr2\n'
    gen chr10 130 - | tail -n +2
} > empty-seq.fa

# irregular.fa: the second line of chr1 is shorter than the first.
printf '>chr1\nACGTACGTAC\nACGTAC\nACGTACGTAC\n' > irregular.fa

for f in ref.fa crlf.fa desc.fa empty-seq.fa; do samtools faidx "$f"; done
samtools faidx ref.fa.gz

# Inputs on which gloci is stricter than samtools (docs/differences.md#fasta-index-strictness):
# samtools keeps the first of two equal names (with a warning), counts a space as
# nothing, and accepts an empty name; gloci refuses each with exit 1.
printf '>chr1\nACGT\n>chr2\nAC\n>chr1\nGG\n' > dup-name.fa
printf '>chr1\nAC GT\nACGT\n' > space.fa
printf '>\nACGT\n' > noname.fa
