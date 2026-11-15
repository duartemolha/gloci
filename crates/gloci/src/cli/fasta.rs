use super::shared::{
    ChromMapArg, FormatArg, InvalidArg, OutputArgs, RegionArgs, ThreadArgs, VcfSpanArg,
};
use clap::Args;
use std::path::PathBuf;

/// The mask character of `maskfasta`: one printable ASCII character.
fn parse_mask_char(s: &str) -> Result<u8, String> {
    match s.as_bytes() {
        [b] if b.is_ascii_graphic() => Ok(*b),
        [b' '] => Err(
            "a space is refused as the mask character: a space inside a sequence line makes the \
FASTA impossible to index (gloci index and samtools faidx reject it); use a printable character \
such as N or X"
                .into(),
        ),
        _ => Err(format!(
            "invalid mask character {s:?}: use exactly one printable ASCII character such as N or X (or --soft to lowercase instead)"
        )),
    }
}

/// Mask bases of a FASTA under intervals, streaming it (bedtools maskfasta).
///
/// Copies --fasta (plain, gzip or BGZF; no index needed) to the output with every base under a -b
/// interval replaced by --mask-char (default N) or, with --soft, lowercased (letters only; other
/// characters are kept). Headers and line layout are kept exactly: line widths, CRLF line ends, blank
/// lines and a missing final line break. (bedtools re-wraps lines to the first line's width, drops
/// header descriptions and cannot match names in a CRLF file; it also adds a final line break that
/// was missing and drops empty sequences, which gloci keeps.)
#[derive(Args)]
#[command(about = "Mask FASTA bases under intervals (bedtools maskfasta)")]
#[command(after_long_help = "\
Examples:
  gloci maskfasta --fasta hg38.fa -b repeats.bed -o hg38.masked.fa
  gloci maskfasta --fasta hg38.fa.gz -b rmsk.bed segdups.bed --soft -o soft.fa.gz --write-index
  gloci maskfasta --fasta ref.fa -b low-mappability.bed --mask-char X | head

Memory:
The FASTA is read in pieces of at most 1 MiB, so a chromosome of any size needs no memory \
beyond the masks.

Masks:
Masks come from any number of files in any order, in any interval format (BED, bedGraph, \
GFF/GTF, VCF/BCF); overlapping masks are merged so a base is masked once. A mask names a \
sequence by the first word of its header; use --chrom-map when the names differ. Masks on \
a chromosome the FASTA lacks are skipped with one warning. A mask that reaches past the \
end of its sequence is an error, unless --skip-invalid is given: then it is clipped to the \
sequence end (one that starts at or past the end masks nothing) with one warning that \
counts them. --mask-char takes one printable character; a space is refused because it \
makes the FASTA unindexable.

Output:
Output goes to stdout or -o FILE (.fa, .fasta, .fna; a .gz name is BGZF). --write-index \
also writes FILE.fai, and FILE.gzi for BGZF, the same files gloci index writes.")]
pub struct MaskfastaArgs {
    /// FASTA to mask (plain, gzip or BGZF; `-` for stdin).
    #[arg(long, value_name = "FILE", required = true)]
    pub fasta: String,
    /// Interval files whose bases are masked (-b x y or -b x -b y).
    #[arg(
        short = 'b',
        long = "b",
        value_name = "FILE",
        required = true,
        num_args = 1..,
        action = clap::ArgAction::Append
    )]
    pub masks: Vec<String>,
    /// Lowercase masked letters instead of replacing them.
    #[arg(long, conflicts_with = "mask_char")]
    pub soft: bool,
    /// Character that replaces masked bases.
    #[arg(long, value_name = "C", default_value = "N", value_parser = parse_mask_char)]
    pub mask_char: u8,
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

/// Write the sequence under each interval (bedtools getfasta).
///
/// --fasta REF.fa is read at random through REF.fa.fai (and REF.fa.gz.gzi for BGZF); build them with
/// gloci index. gloci never builds them itself. The intervals (BED, GFF/GTF, VCF/BCF) are written in
/// input order, one sequence each, on a single line with the case of the FASTA kept and memory that
/// does not grow with the interval length. The header is chrom:start-end with a 0-based start, as
/// bedtools; --name gives name::chrom:start-end and --name-only gives name (BED column 4, the GFF
/// feature type, or REF/ALT_ID for VCF).
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci getfasta --fasta ref.fa genes.bed
  gloci getfasta --fasta ref.fa genes.bed12 --split -s --name-only -o transcripts.fa
  gloci getfasta --fasta ref.fa.gz calls.vcf.gz --vcf-span ref --tab

Sequence:
-s reverse-complements minus-strand records with the IUPAC-aware, case-keeping table of \
bedtools and adds (+), (-) or (.) to the header; --rna writes U for the complement of A \
(as bedtools, it changes reverse-complemented sequences only). --split joins the blocks of \
BED12 records; with -s a minus-strand record is reverse-complemented as a whole. VCF \
follows --vcf-span (--vcf-span ref is bedtools' rule).

Output:
--tab writes the header, TAB, the sequence; --bed-out writes the input record, TAB, the \
sequence. -o FILE.fa.gz writes BGZF.

Invalid intervals:
A zero-length interval has no sequence: it is skipped with one warning that counts them, \
as bedtools does. An interval on a chromosome the FASTA lacks, or ending past its end, is \
an error; bedtools skips it with a message and exits 0. --skip-invalid skips such \
intervals with one warning carrying the count.")]
pub struct GetfastaArgs {
    /// Intervals to read (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The FASTA to read (plain or BGZF), with its .fai (and .gzi) beside it.
    #[arg(short = 'f', long, value_name = "FILE", required = true)]
    pub fasta: PathBuf,
    /// Reverse-complement minus-strand records and add the strand to the
    /// header.
    #[arg(short = 's')]
    pub strand: bool,
    /// Header `name::chrom:start-end`.
    #[arg(long, conflicts_with = "name_only")]
    pub name: bool,
    /// Header `name`.
    #[arg(long)]
    pub name_only: bool,
    /// BED12: join the blocks (exons) of each record.
    #[arg(long)]
    pub split: bool,
    /// Write `header TAB sequence` lines.
    #[arg(long, conflicts_with = "bed_out")]
    pub tab: bool,
    /// Write each input record, TAB, the sequence.
    #[arg(long)]
    pub bed_out: bool,
    /// Write U for the complement of A (reverse-complemented sequences only).
    #[arg(long)]
    pub rna: bool,
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

/// Base composition of the sequence under each interval (bedtools nuc).
///
/// --fasta REF.fa is read at random through REF.fa.fai (and REF.fa.gz.gzi for BGZF); build them with
/// gloci index. gloci never builds them itself. Each interval (BED, GFF/GTF, VCF/BCF) is written as
/// it was read, then pct_at pct_gc num_A num_C num_G num_T num_N num_oth seq_len, preceded by
/// bedtools' header line (#1_usercol ...). The counts ignore case, U counts as T and any other byte
/// (IUPAC codes, -, *) is other; the percentages divide by the interval length, are computed in f32
/// and print with six decimals, as bedtools prints them.
#[derive(Args)]
#[command(about = "Base composition under each interval (bedtools nuc)")]
#[command(after_long_help = "\
Examples:
  gloci nuc --fasta ref.fa genes.bed
  gloci nuc --fasta ref.fa genes.bed -s --seq -o genes.nuc.tsv
  gloci nuc --fasta ref.fa.gz calls.vcf.gz --vcf-span ref --pattern CG

Strand and extra columns:
-s reverse-complements minus-strand records first, so everything below describes the \
reverse complement. --seq adds the sequence; --pattern P adds the number of times P occurs \
(overlapping matches count, case-sensitive unless --ignore-case).

Memory and invalid intervals:
Memory does not grow with the interval length. A zero-length interval is skipped with one \
warning that counts them, as bedtools does. An interval on a chromosome the FASTA lacks, \
or ending past its end, is an error; bedtools skips it with a message and exits 0. \
--skip-invalid skips such intervals with one warning carrying the count.")]
pub struct NucArgs {
    /// Intervals to read (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The FASTA to read (plain or BGZF), with its .fai (and .gzi) beside it.
    #[arg(short = 'f', long, value_name = "FILE", required = true)]
    pub fasta: PathBuf,
    /// Reverse-complement minus-strand records before counting.
    #[arg(short = 's')]
    pub strand: bool,
    /// Add the sequence as a column.
    #[arg(long)]
    pub seq: bool,
    /// Add a column counting this pattern (overlapping matches).
    #[arg(long, value_name = "P")]
    pub pattern: Option<String>,
    /// Match --pattern without regard to case.
    #[arg(short = 'C', long, requires = "pattern")]
    pub ignore_case: bool,
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

const FAIDX_REGION_HELP: &str = "Another region, as the positional REGION (repeatable).";
const FAIDX_REGIONS_FILE_HELP: &str = "Read regions from a BED FILE (0-based, half-open); each prints as chrom:START-END with a 1-based start.";

/// Write regions of an indexed FASTA (samtools faidx).
///
/// REF.fa is read at random through REF.fa.fai (and REF.fa.gz.gzi for BGZF); build them with gloci
/// index. gloci never builds them itself (samtools faidx does): without them, or with no region, this
/// exits 2 naming gloci index REF. Each REGION is chr, chr:START, chr:START-, chr:START-END, chr:-END
/// (1-based, closed, commas allowed) or {name}:START-END for a name containing a colon, read as
/// htslib reads it; --region adds more, and --regions-file reads BED lines (0-based, half-open;
/// samtools -r reads chr:START-END text instead). Each region prints ">" and its text as typed, then
/// its bases wrapped at --line-width (default: that sequence's line width in the .fai; 0 is one
/// line), the case of the FASTA kept.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci faidx ref.fa chr1:1,000-2,000 chrM
  gloci faidx ref.fa chr2:500-900 --reverse-complement --mark-strand sign -o rc.fa
  gloci faidx ref.fa.gz --regions-file genes.bed --line-width 0 -o genes.fa.gz --write-index

Region details:
A BED region is named chrom:START-END with a 1-based start, and a zero-length one is \
skipped with one warning that counts them. An end past the sequence is cut there with a \
warning; a start past it writes the header alone. A sequence the index lacks is an error \
before anything is written, unless --skip-invalid skips that region entirely (samtools -c \
still writes its header).

Output:
--reverse-complement writes the reverse complement, with samtools' table (IUPAC-aware, \
case kept), and --mark-strand chooses the header suffix: rc (/rc on reverse-complemented \
regions, the default), no, sign ((+) or (-)) or custom,POS,NEG. Memory does not grow with \
the region length. -o FILE.fa.gz writes BGZF, and --write-index also writes FILE.fai (and \
FILE.gzi).")]
#[command(mut_arg("region", |a| a.help(FAIDX_REGION_HELP)))]
#[command(mut_arg("regions_file", |a| a.help(FAIDX_REGIONS_FILE_HELP)))]
pub struct FaidxArgs {
    /// The FASTA to read (plain or BGZF), with its .fai (and .gzi) beside it.
    #[arg(value_name = "REF")]
    pub fasta: PathBuf,
    /// Regions to write: chr, chr:START, chr:START-END, chr:-END or
    /// {name}:START-END (1-based, closed, commas allowed).
    #[arg(value_name = "REGION")]
    pub regions: Vec<String>,
    #[command(flatten)]
    pub region: RegionArgs,
    /// Bases per output line; 0 writes each sequence on one line (default:
    /// the line width of that sequence in the .fai).
    #[arg(long, value_name = "N")]
    pub line_width: Option<u64>,
    /// Write the reverse complement of each region.
    #[arg(long)]
    pub reverse_complement: bool,
    /// Strand suffix of the header: rc (/rc on reverse-complemented regions),
    /// no, sign ((+)/(-)) or custom,POS,NEG.
    #[arg(long, value_name = "TYPE", default_value = "rc")]
    pub mark_strand: String,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Write a SAM sequence dictionary of a FASTA (samtools dict).
///
/// Streams REF (plain, gzip or BGZF; standard input is -); no index is
/// needed. Writes "@HD VN:1.0 SO:unsorted" (unless --no-header), then one
/// "@SQ SN:<name> LN:<length> M5:<md5>" line per sequence, where the name
/// is the first word of the header. Length and MD5 count only the bytes
/// "!" to "~" with a-z uppercased, as samtools does: IUPAC codes, "*" and
/// "-" count; spaces, tabs and carriage returns do not. More tags follow in
/// samtools' order: AH:* (--alt), AN: (--alias), UR: and AS:, SP:. UR
/// defaults to file:// and the canonical path of REF (none for standard
/// input); --uri replaces it. --alias adds the name with "chr" removed or
/// added; M and MT also list chrMT,MT or chrM,M. --alt reads one name per
/// line (up to the first tab). A name used twice is an error naming both
/// lines (samtools writes both). The result works as -g for gloci
/// commands.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci dict ref.fa -o ref.dict
  gloci dict ref.fa.gz --alias --assembly GRCh38 --species human --alt ref.fa.alt -o ref.dict")]
pub struct DictArgs {
    /// The FASTA to read (plain, gzip or BGZF; - is standard input).
    #[arg(value_name = "REF")]
    pub fasta: String,
    /// Leave out the @HD line.
    #[arg(long)]
    pub no_header: bool,
    /// Assembly name for the AS: tag.
    #[arg(long, value_name = "STR")]
    pub assembly: Option<String>,
    /// Species for the SP: tag.
    #[arg(long, value_name = "STR")]
    pub species: Option<String>,
    /// URI for the UR: tag (default: file:// and the absolute path of REF).
    #[arg(long, value_name = "STR")]
    pub uri: Option<String>,
    /// Add an AN: tag by adding or removing chr.
    #[arg(long)]
    pub alias: bool,
    /// Add AH:* to the sequences named in FILE (one name per line, up to the
    /// first tab).
    #[arg(long, value_name = "FILE")]
    pub alt: Option<PathBuf>,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the dictionary to FILE instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}
