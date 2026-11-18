use super::shared::{
    AlignArgs, ChromMapArg, ChromOrderArg, FormatArg, GenomeArg, HeaderArg, InvalidArg, OutputArgs,
    RegionArgs, ThreadArgs, TmpDirArg, VcfSpanArg, parse_size,
};
use clap::{Args, Subcommand, ValueEnum};
use std::path::PathBuf;

fn parse_comment(s: &str) -> Result<u8, String> {
    match s.as_bytes() {
        [b] if b.is_ascii() => Ok(*b),
        _ => Err(format!(
            "invalid comment character {s:?}: use one ASCII character"
        )),
    }
}

/// A variant kind for `convert --types` (spec 4.2).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum KindArg {
    Snv,
    Mnv,
    Ins,
    Del,
    Sv,
    Ref,
    Other,
}

impl From<KindArg> for gloci_core::VariantKind {
    fn from(k: KindArg) -> Self {
        use gloci_core::VariantKind as K;
        match k {
            KindArg::Snv => K::Snv,
            KindArg::Mnv => K::Mnv,
            KindArg::Ins => K::Ins,
            KindArg::Del => K::Del,
            KindArg::Sv => K::Sv,
            KindArg::Ref => K::Ref,
            KindArg::Other => K::Other,
        }
    }
}

/// Convert VCF/BCF, GFF3/GTF, WIG, PSL, RepeatMasker .out, BED12 or BAM/CRAM/SAM to BED,
/// BAM/CRAM/SAM to FASTQ, BED, GFF, VCF and BEDPE to BAM, bedGraph, WIG or bigWig to bigWig, BED
/// to bigBed, and bigWig and bigBed to bedGraph and BED (bedops convert2bed, vcf2bed, gff2bed,
/// gtf2bed; bedtools bed12tobed6, bamtobed, bamtofastq, bedtobam, bedpetobam; bedGraphToBigWig,
/// wigToBigWig, bedToBigBed, bigWigToBedGraph, bigBedToBed).
#[derive(Args)]
#[command(
    about = "Convert to BED/FASTQ/BAM/bigWig/bigBed (convert2bed, bedtobam)",
    long_about
)]
#[command(after_long_help = "\
Examples:
  gloci convert calls.vcf.gz -o calls.bed
  gloci convert cov.bedgraph -g hg38 -o cov.bw
  gloci convert reads.byname.bam --fastq -o r1.fq.gz --fastq2 r2.fq.gz
  gloci convert peaks.bed --bam -g hg38 -o peaks.bam

VCF and BCF:
VCF/BCF gives one row per record: chrom, start, end, ID, REF, ALT, QUAL, FILTER, with the \
footprint chosen by --vcf-span. --split-alleles writes one row per ALT instead, each with \
that ALT's own footprint (with --vcf-span ref, the record's REF footprint on every row; \
`*`, `<*>` and `<NON_REF>` rows are skipped when other ALTs exist). --types keeps the \
records (with --split-alleles, the ALTs) of the listed kinds: snv, mnv, ins, del, sv, ref, \
other; a record is kept when any of its ALTs is of a listed kind.

GFF3 and GTF:
GFF3/GTF gives 10 columns (bedops gff2bed): chrom, start, end, name, score, strand, \
source, type, phase, attributes. The name is the --name-attr attribute (default ID for \
GFF3, gene_id for GTF; `.` when absent). --feature keeps the listed column-3 types.

WIG:
WIG (fixedStep or variableStep, with optional span) gives one bedGraph row per data line: \
chrom, start, end, value; equal values are not merged. fixedStep needs start= and step=; \
values must be finite. A span larger than the step gives overlapping rows.

PSL:
PSL (BLAT; the header is skipped) gives bedops psl2bed's 21 columns: tName, tStart, tEnd, \
qName, matches, strand, then the other PSL columns in order.

RepeatMasker:
RepeatMasker .out gives bedops rmsk2bed's columns: sequence, begin-1, end, repeat, score, \
strand (`C` is `-`), the divergence, deletion and insertion percentages, the `(left)` of \
the query, class/family, the repeat positions as written, ID and `*`.

Alignments:
BAM, CRAM and SAM give one BED6 row per mapped read, in input order: chrom, start, end, \
name (with /1 or /2 for the first or second mate), MAPQ, strand; unmapped reads are \
skipped and counted at -v, and --min-mapq, --require-flags and --exclude-flags filter. The \
span is the reference bases the CIGAR covers (M D N = X), cut at the end of the reference. \
--bed12 writes one BED12 row per read with the thick span equal to the read and --color as \
itemRgb; --blocks writes one BED6 row per block, cut at N; --split-deletions cuts at D as \
well. --score-tag NM takes the score from an integer tag (NM is bedtools -ed), and --cigar \
appends the CIGAR to BED6 rows.

BEDPE from alignments:
--bedpe (bedtools bamtobed -bedpe) writes BEDPE instead, one row per read name from \
alignments sorted by name (gloci sort -n; a coordinate sorted file is exit 2): chrom1 \
start1 end1 chrom2 start2 end2 name score strand1 strand2, from the primary READ1 and \
READ2 that are next to each other (secondary and supplementary reads are ignored). An \
unmapped mate is `. -1 -1` with strand `.`. The end that sorts first by (chromosome as \
bytes, start) comes first; --mate1-first keeps READ1 first (bedtools -mate1). The score is \
the smaller MAPQ of two mapped ends, else 0; with --score-tag it is the sum over the \
mapped ends (--score-tag NM is bedtools -ed). A name with one mate is not written and is \
counted in a warning. The output is .bedpe (or .bedpe.gz) and feeds pair-intersect and \
pair-pair; --bedpe cannot be used with --bed12, --blocks, --cigar, --color or \
--split-deletions.

FASTQ from alignments:
--fastq (bedtools bamtofastq) writes the reads as FASTQ instead of BED: every record, in \
input order, with secondary, supplementary, unmapped, QC-fail and duplicate reads \
included (--exclude-flags 0x900 drops the first two), under its own name with no /1 or /2. \
A read on the reverse strand is written reverse-complemented with its qualities reversed; \
a missing QUAL is `!` for every base and a missing SEQ is an empty read. The output is \
.fq or .fastq, with .gz for BGZF, or plain text on stdout. --fastq2 FILE writes the second \
mates there (bedtools -fq2) and needs the input sorted by name (gloci sort -n; a \
coordinate sorted file is exit 2): the records of a name that are next to each other are \
the pair, READ1 goes to -o as NAME/1 and READ2 to FILE as NAME/2, in either order in the \
input. A name with one record, or two that are not a READ1 and a READ2, is skipped and \
counted in a warning; a name with more than two records is an error (drop the secondary \
and supplementary alignments with --exclude-flags 0x900). With -o and --fastq2 both on \
stdout, or both the same file, the ends are interleaved. --tags (bedtools -tags) takes \
the second end from the R2 and Q2 tags of the first record of each pair, as they are \
(--fastq2 is needed). --fastq cannot be used with the options of the BED outputs.

BAM from features:
--bam (bedtools bedtobam) writes BED, GFF or VCF records as BAM instead of BED, one record \
per feature in input order, to -o (.bam) or stdout. The header is @HD SO:unsorted and one \
@SQ per chromosome of -g, in its order, with AS: the -g value as given; a feature on a \
chromosome -g lacks is exit 1. A record is the name (BED column 4: BED3 is exit 1; the GFF \
feature type; REF/ALT_ID for VCF), FLAG 16 on the - strand and 0 otherwise, the start, \
--mapq (default 255), the length as M, and no mate, bases or qualities. A zero-length BED \
record is one base wider on each side, as bedtools reads it. --bed12 makes each block of a \
BED12 record an M and each gap an N (adjacent blocks stay separate M operations); a record \
that is not valid BED12 is exit 1. With --format bedpe (bedtools bedpetobam) each BEDPE \
pair, which needs a name (column 7), is two records: READ1 (FLAG 65) and READ2 (FLAG 129), \
plus 16 when the end is on - and 32 when its mate is, and 2 when both ends are on one \
chromosome with the first on + and the second on -; the mate fields point at the other end \
and TLEN on both is start2 - start1 (0 across chromosomes). An absent end (. -1 -1) is \
unmapped (FLAG 4, and 8 on its mate) and placed at its mate. --compression-level 0 is \
bedtools -ubam.

BED12:
BED12 with --blocks gives one BED6 row per block; a line that is not valid BED12 (fewer \
than 12 columns, block counts that disagree, blocks out of order, a last block not ending \
at chromEnd) is an error, also with --skip-invalid.

bigWig and bigBed:
A .bw or .bigwig output is a bigWig from bedGraph, a 4-column BED, WIG or bigWig; a .bb or \
.bigbed output is a bigBed from BED or bigBed (-O bigwig|bigbed names it too). The output \
must be a regular file; the chromosome sizes come from -g, else from the bigWig or bigBed \
being converted. Each chromosome's records must be together and sorted by start (gloci \
sort); a bigWig takes no overlaps and stores 32-bit floats. The first record sets a \
bigBed's columns; --autosql FILE must agree, and its leading fields that are BED's by name \
and type are the standard columns; a bigBed input keeps its own autoSql and standard \
columns when the autoSql's field count equals the records' columns (else a warning, and a \
generated autoSql). A record past its chromosome or on one -g lacks, an empty input and --write-index \
are errors; --skip-invalid skips lines that do not parse, never a record the writer refuses \
(a value NA or nan, an overlap). The file is checked, then renamed into place: written whole \
or not at all.

Order and output:
The input is converted in its own order and needs no sorting. A conversion whose output \
goes out of order (allele footprints can start after their record) prints one warning; \
--write-index then fails.")]
pub struct ConvertArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// VCF/BCF: one row per ALT, with that ALT's footprint.
    #[arg(long)]
    pub split_alleles: bool,
    /// VCF/BCF: keep records (with --split-alleles, ALTs) of these kinds.
    #[arg(long, value_enum, value_delimiter = ',', value_name = "LIST")]
    pub types: Vec<KindArg>,
    /// GFF3/GTF: keep this column-3 type (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',', value_name = "NAME")]
    pub feature: Vec<String>,
    /// GFF3/GTF: the attribute for the name column (default ID for GFF3,
    /// gene_id for GTF).
    #[arg(long, value_name = "KEY", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    pub name_attr: Option<String>,
    /// BED12, BAM, CRAM, SAM: write one BED6 row per block (for alignments the
    /// segments between N operations, as bedtools bamtobed -split).
    #[arg(long)]
    pub blocks: bool,
    /// BAM, CRAM, SAM: write one BED12 row per read, its blocks cut at N
    /// (bedtools bamtobed -bed12); with --bam: read BED12 blocks as M and the gaps as N
    /// (bedtools bedtobam -bed12).
    #[arg(long)]
    pub bed12: bool,
    /// BAM, CRAM, SAM with --bed12: the itemRgb column, as R,G,B.
    #[arg(long, value_name = "R,G,B")]
    pub color: Option<String>,
    /// BAM, CRAM, SAM: take the score column from this integer tag (NM, AS, ...)
    /// instead of the MAPQ; every read needs it (bedtools bamtobed -tag, -ed is NM).
    #[arg(long, value_name = "TAG")]
    pub score_tag: Option<String>,
    /// BAM, CRAM, SAM: append the CIGAR as a seventh column (BED6 only).
    #[arg(long)]
    pub cigar: bool,
    /// BAM, CRAM, SAM: also cut blocks at deletions (bedtools bamtobed -splitD);
    /// implies --blocks unless --bed12 is given.
    #[arg(long)]
    pub split_deletions: bool,
    /// BAM, CRAM, SAM sorted by name: write BEDPE pairs of READ1 and READ2
    /// instead of BED (bedtools bamtobed -bedpe).
    #[arg(long, conflicts_with_all = ["bed12", "blocks", "cigar", "color", "split_deletions"])]
    pub bedpe: bool,
    /// With --bedpe: keep READ1 as the first end instead of putting the end
    /// that sorts first there (bedtools bamtobed -mate1).
    #[arg(long, requires = "bedpe")]
    pub mate1_first: bool,
    /// BAM, CRAM, SAM: write the reads as FASTQ instead of BED, to -o or stdout
    /// (bedtools bamtofastq).
    #[arg(long, conflicts_with_all = [
        "split_alleles", "types", "feature", "name_attr", "blocks", "bed12", "color",
        "score_tag", "cigar", "split_deletions", "bedpe", "mate1_first",
    ])]
    pub fastq: bool,
    /// With --fastq: write the second mate of each pair to FASTQ and the first to -o;
    /// the input must be sorted by name (bedtools -fq2).
    #[arg(long, value_name = "FASTQ", requires = "fastq")]
    pub fastq2: Option<PathBuf>,
    /// With --fastq2: take the second end from the R2 and Q2 tags of the first record of
    /// each pair (bedtools -tags).
    #[arg(long, requires = "fastq2")]
    pub tags: bool,
    /// BED, GFF, VCF, or BEDPE with --format bedpe: write the records as BAM instead of BED,
    /// to -o or stdout (bedtools bedtobam, bedpetobam); needs -g.
    #[arg(long, requires = "genome", conflicts_with_all = [
        "split_alleles", "types", "feature", "name_attr", "blocks", "color", "score_tag",
        "cigar", "split_deletions", "bedpe", "mate1_first", "fastq", "min_mapq",
        "require_flags", "exclude_flags", "reference", "vcf_span",
    ])]
    pub bam: bool,
    /// With --bam, or for a bigWig or bigBed output: the references of the BAM header, or the
    /// chromosome sizes, in their order: a built-in genome
    /// such as hg38 (see gloci genomes list), or a .fai, .genome, .dict, VCF/BCF or BAM header
    /// file.
    #[arg(short, long, value_name = "NAME|FILE")]
    pub genome: Option<String>,
    /// bigBed output: the autoSql to store, as given (bedToBigBed -as); its field count must equal
    /// the records' column count, which the first record sets.
    #[arg(long, value_name = "FILE")]
    pub autosql: Option<PathBuf>,
    /// With --bam: the MAPQ of every record, 0 to 255 (default 255; bedtools -mapq).
    #[arg(long, value_name = "N", requires = "bam")]
    pub mapq: Option<u8>,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Validate files and print a one-screen summary of each (samtools quickcheck).
///
/// check reads every file to the end, keeps going past problems and reports the record count, the
/// chromosomes, whether the records are sorted, and the first 10 errors and warnings with their line
/// numbers. --quick only checks the header, the first record and the end-of-file marker, as samtools
/// quickcheck does.
#[derive(Args)]
#[command(about = "Validate files and summarise each (samtools quickcheck)")]
#[command(after_long_help = "\
Examples:
  gloci check peaks.bed
  gloci check calls.vcf.gz -g hg38.fa.fai
  gloci check reads.bam
  gloci check --quick --list-failed *.bam > bad.txt

Errors and warnings:
Errors: records that do not parse, records out of order, an end past the chromosome length \
or a chromosome missing from -g, bad BED12 blocks, a bedGraph value that is not a number, \
a VCF line with the wrong number of columns. Warnings: exact duplicate lines at the same \
position, and blank lines. Exit status 1 when any file has errors. The order is checked \
against -g, else the ##contig lines of a VCF/BCF, else --chrom-order (lex by \
default).

Alignments:
A BAM, CRAM or SAM is read record by record (a CRAM needs --reference): the summary adds \
the mapped, unmapped-but-placed and unplaced counts and the @HD sort order. With \
SO:coordinate declared, a read before the previous one, or a placed read after an unplaced \
one, is an error; a read past the end of its reference and identical records at one \
position are warnings. A BAM without its BGZF end-of-file block (a CRAM without its \
end-of-file container) is an error. -g, --chrom-order and --chrom-map do not apply to them \
(-g and --chrom-order give one warning, --chrom-map is a usage error): the header says the \
references and the order. Duplicate warnings assume coordinate-sorted input, as they \
compare only records at the same position.

Quick check:
--quick is samtools quickcheck, for any format: the header must parse (alignments: with an \
@SQ line, unless --unaligned), the first record must parse (a CRAM: its first container, \
so no reference is needed) and a file (not stdin) must end with its end-of-file marker. \
Only those parts are read. --list-failed prints just the failing paths. Exit status 1 for \
any failing file, not samtools' bit mask.")]
pub struct CheckArgs {
    /// Input files.
    #[arg(value_name = "FILES", required = true)]
    pub files: Vec<String>,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Only check that each file is usable, reading just its header, its first record and
    /// its end-of-file marker: the last 28 bytes of a BGZF file, 38 of a CRAM (samtools
    /// quickcheck).
    #[arg(long)]
    pub quick: bool,
    /// With --quick, accept a BAM, CRAM or SAM header with no @SQ line (samtools quickcheck -u).
    #[arg(long, requires = "quick")]
    pub unaligned: bool,
    /// Print only the paths of the files that fail, one per line, and what is wrong with each
    /// on stderr (samtools quickcheck -v).
    #[arg(long)]
    pub list_failed: bool,
    /// The FASTA a CRAM input was made with, with its .fai (and .gzi when BGZF). Needed to
    /// read the records of a CRAM; --quick never needs it.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
}

/// Print the records in regions, stream or recompress a file, convert VCF to BCF and back, or print
/// only the header (samtools view).
///
/// With no regions the records are written as read (comment and track lines after the first record
/// and blank lines are dropped, and CRLF line ends become LF), so view can recompress a file (-o
/// FILE.gz is BGZF) or index it (--write-index). VCF and BCF convert to each other with -O or the -o
/// extension; the header is kept as it is, plus a ##glociCommand line. BED and GFF header lines are
/// kept by default. For conversions to BED from VCF/BCF, GFF3/GTF or BED12, use gloci convert.
#[derive(Args)]
#[command(about = "Print regions, recompress or convert a file (samtools view)")]
#[command(after_long_help = "\
Examples:
  gloci view calls.vcf.gz chr1:1,000,000-2,000,000 chrX
  gloci view reads.cram --reference ref.fa -o reads.bam
  gloci view --count reads.bam chr1
  gloci view reads.bam --regions-file targets.bed -o on_target.bam

Alignments:
BAM, CRAM and SAM input is written as BAM or SAM, as -O or the -o extension say (the same \
as the input by default; BAM for CRAM), every record, unmapped reads included, with an @PG \
line added to the header. --header-only, --no-header (SAM output only), --min-mapq, \
--require-flags and --exclude-flags apply. Regions select the reads whose footprint \
overlaps them; an unmapped read placed at a position counts as one base there. It is not \
samtools view: there is no -b, -h, -c or -L. CRAM needs the FASTA it was made with \
(--reference, with its .fai), unless it was written without one; MD and NM tags come out \
as stored, not regenerated as samtools does.

Regions:
With regions, view prints the records whose footprint (per --vcf-span) overlaps a region, \
each once, in file order. A point (a zero-length BED record, an insertion) overlaps a \
region only when strictly inside it, unless --inclusive-points is given. Regions are read \
through a .tbi/.csi index when the file has one, and otherwise by scanning the whole file, \
with a warning. Either way the order of every record read is checked, so an unsorted file \
is an error, never a wrong answer. With --assume-sorted, an uncompressed BED file without \
an index is binary-searched instead (like bedextract); that trusts the file to be sorted. \
Index the file (gloci sort FILE -o FILE.gz --write-index) for speed and checking both. \
A bigWig or bigBed (experimental) is read through its own index: bigWig records are \
clipped to each region, bigBed records are whole (D171).")]
pub struct ViewArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Regions such as chr1, chr1:1,000-2,000 or chr1:5,000- (1-based,
    /// inclusive); write {name}:start-end for a name containing ':'.
    #[arg(value_name = "REGION")]
    pub regions: Vec<String>,
    #[command(flatten)]
    pub region: RegionArgs,
    /// A point on the edge of a region overlaps it.
    #[arg(long)]
    pub inclusive_points: bool,
    /// Find regions in an uncompressed BED file without an index by binary
    /// search, like bedextract. Fast, but it trusts that the file is sorted
    /// (in -g or --chrom-order order) and checks only the chromosome block
    /// it reads; on an unsorted file the answer can be incomplete.
    #[arg(long)]
    pub assume_sorted: bool,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    /// Write only the header lines.
    #[arg(long, conflicts_with = "no_header")]
    pub header_only: bool,
    /// Write only the records, without the header.
    #[arg(long)]
    pub no_header: bool,
    /// Write the number of records that would be written (regions and read filters
    /// applied), as one line, instead of the records; as samtools view -c. Without
    /// regions, a CRAM decodes only the flags and MAPQ, as samtools does.
    #[arg(long, conflicts_with_all = ["header_only", "no_header", "output_format", "write_index"])]
    pub count: bool,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Sort BED/bedGraph/GFF/VCF/BCF files, or BAM/SAM/CRAM alignments, with bounded memory (bedtools sort, samtools sort).
///
/// Records are ordered by chromosome, then start, then end. The chromosome order comes from -g, else
/// from the ##contig lines of a VCF/BCF header, else from --chrom-order (lex by default).
/// Records that tie keep their input order unless --full-key is given.
#[derive(Args)]
#[command(about = "Sort BED, VCF or alignments (bedtools sort, samtools sort)")]
#[command(after_long_help = "\
Examples:
  gloci sort reads.bam -o reads.sorted.bam --write-index -t 8
  gloci sort peaks.bed > peaks.sorted.bed
  gloci sort huge.bed --max-mem 8G --tmp-dir /scratch -t 16 -o huge.sorted.bed.gz
  gloci sort --check peaks.sorted.bed

Records and output:
BED header lines are kept only with --header; the leading ## lines of a GFF and the header \
of a VCF/BCF are always kept. VCF/BCF records are ordered by POS and REF end whatever \
--vcf-span says, so the output is tabix-compatible. The output format is the input format: \
-o must match it. Use -o with the same path as an input to sort a file in place.

Memory:
Sorting spills to --tmp-dir when the input is larger than --max-mem. That directory must \
exist and be writable; it is checked before any input is read.

Alignments:
Alignments (one BAM, SAM or CRAM file, or stdin) are sorted as samtools sort does: by \
reference in header order (unplaced reads last), position and strand, ties in input order; \
-n by read name in natural order (r2 before r10) and -N by plain byte order, each with \
READ1 before READ2 and primary before supplementary before secondary. Every record is \
kept, unmapped reads included. The header's @HD line gets SO (and SS for -n and -N) and \
gloci adds an @PG line. The output keeps the input's format (BAM for BAM and CRAM, SAM for \
SAM; -o and -O override), where samtools sort always writes BAM. --write-index works for \
coordinate-sorted BAM output only. -g, --full-key, --unique and --duplicates do not apply: \
the order is the header's. To merge sorted BAM files use sort -m (or cat for files already \
in order); -n reads.bam sorts by name.")]
pub struct SortArgs {
    /// Input files (default: stdin).
    #[arg(value_name = "FILES")]
    pub files: Vec<String>,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub format: FormatArg,
    /// Memory budget for sorting, in binary units (K, M, G, T); at least 16M.
    #[arg(long, value_name = "SIZE", default_value = "2G", value_parser = parse_size)]
    pub max_mem: u64,
    #[command(flatten)]
    pub tmp_dir: TmpDirArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Break ties on the whole record instead of keeping input order.
    #[arg(long)]
    pub full_key: bool,
    /// Write one copy of each distinct record (implies --full-key).
    #[arg(long, conflicts_with = "duplicates")]
    pub unique: bool,
    /// Write one copy of each record that occurs more than once (implies --full-key).
    #[arg(long)]
    pub duplicates: bool,
    /// Merge inputs that are already sorted instead of sorting them.
    #[arg(short = 'm', long)]
    pub merge_sorted: bool,
    /// Only check that the inputs are sorted; write nothing, exit 1 if not.
    #[arg(long)]
    pub check: bool,
    /// Alignments: sort by read name, in natural order (r2 before r10), as samtools sort -n.
    #[arg(short = 'n', long, conflicts_with = "by_name_bytes")]
    pub by_name: bool,
    /// Alignments: sort by read name as bytes (r10 before r2), as samtools sort -N.
    #[arg(short = 'N', long)]
    pub by_name_bytes: bool,
    /// The FASTA a CRAM input was made with (bedtools --cram-ref), with its .fai (and .gzi
    /// when BGZF). Not needed for a CRAM written without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// List the built-in genomes that -g accepts, or print one.
///
/// gloci genomes list prints one line per built-in genome, and gloci
/// genomes with no subcommand prints the same list. gloci genomes show
/// NAME prints one genome as a .genome file.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci genomes list
  gloci genomes show hg38
  gloci sort peaks.bed -g hg38 > peaks.sorted.bed")]
#[command(args_conflicts_with_subcommands = true)]
pub struct GenomesArgs {
    #[command(subcommand)]
    pub action: Option<GenomesAction>,
}

#[derive(Subcommand)]
pub enum GenomesAction {
    /// List the built-in genomes, one tab-separated line each.
    ///
    /// The columns are id, species, assembly, naming (ucsc or ensembl),
    /// sequences and total_bp, with no header line. The ids are what -g and
    /// gloci genomes show take. The lengths come from UCSC chrom.sizes files
    /// and from the Ensembl core databases; the sources are recorded in
    /// crates/gloci-io/data/genomes/SOURCES.toml in the source tree.
    #[command(after_long_help = "\
Examples:
  gloci genomes list
  gloci genomes list | grep -i sapiens
  gloci genomes list | cut -f1")]
    List,
    /// Print one built-in genome as name and length, tab-separated.
    ///
    /// The chromosomes come in the order -g uses for sorting. The output is
    /// a valid .genome file, so it can also be saved and edited.
    #[command(after_long_help = "\
Examples:
  gloci genomes show hg38
  gloci genomes show GRCh38 > GRCh38.genome
  gloci genomes show mm39 | cut -f1")]
    Show {
        /// Genome id (case-insensitive), as listed by gloci genomes list.
        name: String,
    },
}

/// Build a .bai/.csi index for a BAM, a .crai for a CRAM, a .tbi/.csi for a sorted BGZF file, or a
/// .fai for a FASTA (samtools index, tabix).
///
/// A coordinate-sorted BAM gets PATH.bai, or PATH.csi with --csi or when a reference is longer than
/// 2^29; samtools and gloci read each other's indexes. An unsorted BAM is refused: gloci sort x.bam
/// -o sorted.bam --write-index sorts and indexes in one step. SAM text has no index: write BAM with
/// gloci view x.sam -o x.bam.
///
/// Plain gzip is refused: recompress with bgzip.
#[derive(Args)]
#[command(about = "Index a BAM, CRAM, FASTA or BGZF file (samtools index, tabix)")]
#[command(after_long_help = "\
Examples:
  gloci index sorted.bam
  gloci index calls.vcf.gz
  gloci index reads.cram
  gloci index scores.tsv.gz -s 1 -b 2 -e 3 --zero-based

CRAM:
A CRAM file gets PATH.crai, with the entries samtools index writes; it needs no reference, \
even for slices that hold several references.

FASTA:
A FASTA file (plain or BGZF) gets PATH.fai, and PATH.gzi when it is BGZF, byte-identical \
to samtools faidx. Every line of a sequence but the last must have the same length. \
Reading a FASTA at random (getfasta, nuc, a CRAM reference) needs these files and never \
builds them itself; maskfasta streams its FASTA and needs no index.

Text and VCF files:
For any other BGZF file, BED, GFF/GTF, VCF and BCF files \
are recognised by content. Any other tab-separated file is indexed like tabix does when \
you name its columns with -s (sequence), -b (begin) and -e (end, optional). The columns \
are numbered from 1. Without --zero-based the begin column is 1-based and a feature covers \
[begin-1, end]; with it, [begin, end). Without -e a feature covers one base. Lines \
starting with --comment (default #) and the first --skip-lines lines are not indexed; both \
are recorded in the index. Chromosomes must be contiguous and sorted by begin. The index \
is a .tbi, or a .csi when an end passes 2^29.")]
pub struct IndexArgs {
    /// File to index: a coordinate-sorted BAM (.bai or .csi), a CRAM (.crai), a FASTA (.fai, and
    /// .gzi when BGZF), or a sorted BGZF-compressed BED/GFF/VCF/BCF or tab-separated file (.tbi or
    /// .csi).
    pub file: PathBuf,
    /// Sequence (chromosome) column, numbered from 1. Selects TSV mode; needs -b.
    #[arg(short = 's', long, value_name = "N")]
    pub seq_col: Option<usize>,
    /// Begin column, numbered from 1. Selects TSV mode; needs -s.
    #[arg(short = 'b', long, value_name = "N")]
    pub begin_col: Option<usize>,
    /// End column, numbered from 1; without it a feature covers one base.
    /// Selects TSV mode.
    #[arg(short = 'e', long, value_name = "N")]
    pub end_col: Option<usize>,
    /// The begin column is 0-based (TSV mode).
    #[arg(long)]
    pub zero_based: bool,
    /// Skip the first N lines (TSV mode).
    #[arg(long, value_name = "N")]
    pub skip_lines: Option<usize>,
    /// Lines starting with this character are comments (TSV mode, default #).
    #[arg(long, value_name = "CHAR", value_parser = parse_comment)]
    pub comment: Option<u8>,
    /// Write a .csi instead of a .bai for a BAM file. A BAM with a reference longer than
    /// 2^29 always gets one.
    #[arg(long)]
    pub csi: bool,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// The FASTA a CRAM was made with. Indexing a CRAM reads no reference bases, so it is
    /// optional; when given, its sequence lengths are checked against the CRAM's header.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
}

/// Print the header and the first records of any input as text (samtools and bcftools head).
///
/// With no options it prints the whole header and no record. --headers N keeps the first N
/// header lines and --records N prints N records after them. BAM, CRAM and SAM print SAM
/// text; a BAM whose header text lacks the @SQ lines of its binary references gets them
/// written out (samtools prints nothing for such a file). VCF and BCF print VCF text, and
/// BED, GFF and the other text formats print their header and data lines as stored.
/// Nothing is added: no @PG line, no ##glociCommand line. Only the header and the records
/// asked for are read. A file that is damaged after the lines printed exits 1 with the
/// lines before the damage already written. A CRAM that needs a reference sequence needs
/// --reference. head describes the whole file, so it takes no read filters and no
/// --chrom-map.
#[derive(Args)]
#[command(about = "Print the header and first records (samtools, bcftools head)")]
#[command(after_long_help = "\
Examples:
  gloci head reads.bam
  gloci head reads.bam --headers 0 --records 5
  gloci head calls.vcf.gz --records 2")]
pub struct HeadArgs {
    /// The input: BAM, CRAM, SAM, VCF, BCF, BED, GFF or another text format (- is standard
    /// input).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Print the first N header lines (default: all of them).
    #[arg(long, value_name = "N")]
    pub headers: Option<u64>,
    /// Print N records after the header lines.
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub records: u64,
    /// The FASTA a CRAM input was made with, with its .fai (and .gzi when BGZF). Not needed
    /// for a CRAM written without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the lines to FILE instead of stdout, as text whatever its extension (.sam, .vcf,
    /// .txt, ...); .gz compresses them. A file is written only when head succeeds.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}
