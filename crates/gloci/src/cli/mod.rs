//! clap definitions: global flags, shared option groups and the subcommands.
use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

mod files;
mod intervals;
mod relate;
mod sampling;
mod shared;
mod vcf;
mod vcf_edit;

pub use files::*;
pub use intervals::*;
pub use relate::*;
pub use sampling::*;
pub use shared::*;
pub use vcf::*;
pub use vcf_edit::*;

/// The first line of the help of every command outside the stable core
/// ([`crate::CORE_COMMANDS`]); `every_command_is_core_or_marked_experimental` keeps the two in step.
pub const EXPERIMENTAL: &str = "Experimental: checked against the tool it replaces, but its options and output may \
still change between releases.";

#[derive(Parser)]
#[command(
    name = "gloci",
    version,
    about = "Genomic interval toolkit",
    arg_required_else_help = true,
    after_help = "Stable core: genomes, index, sort, view, intersect, merge, subtract, closest, coverage.\n\
Every other command is experimental: checked against the tool it replaces, but its options and\n\
output may still change between releases."
)]
pub struct Cli {
    /// Only print errors.
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    pub quiet: bool,
    /// Print more progress information (repeat for debug output).
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    Genomes(GenomesArgs),
    Index(IndexArgs),
    Sort(SortArgs),
    View(ViewArgs),
    #[command(before_help = EXPERIMENTAL)]
    Cat(CatArgs),
    #[command(before_help = EXPERIMENTAL)]
    Check(CheckArgs),
    Intersect(IntersectArgs),
    Merge(MergeArgs),
    #[command(before_help = EXPERIMENTAL)]
    Cluster(ClusterArgs),
    #[command(before_help = EXPERIMENTAL)]
    Spacing(SpacingArgs),
    #[command(before_help = EXPERIMENTAL)]
    Complement(ComplementArgs),
    #[command(before_help = EXPERIMENTAL)]
    Partition(PartitionArgs),
    #[command(before_help = EXPERIMENTAL)]
    Genomecov(GenomecovArgs),
    #[command(before_help = EXPERIMENTAL)]
    Slop(SlopArgs),
    #[command(before_help = EXPERIMENTAL)]
    Flank(FlankArgs),
    #[command(before_help = EXPERIMENTAL)]
    Shift(ShiftArgs),
    #[command(before_help = EXPERIMENTAL)]
    Windows(WindowsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Convert(ConvertArgs),
    #[command(before_help = EXPERIMENTAL)]
    Getfasta(GetfastaArgs),
    #[command(before_help = EXPERIMENTAL)]
    Nuc(NucArgs),
    Subtract(SubtractArgs),
    Closest(ClosestArgs),
    #[command(before_help = EXPERIMENTAL)]
    Map(MapArgs),
    Coverage(CoverageArgs),
    #[command(before_help = EXPERIMENTAL)]
    Groupby(GroupbyArgs),
    #[command(before_help = EXPERIMENTAL)]
    Maskfasta(MaskfastaArgs),
    #[command(before_help = EXPERIMENTAL)]
    Vcf(Box<VcfArgs>),
    // P3b Task 2
    #[command(before_help = EXPERIMENTAL)]
    Stats(StatsArgs),
    // P3b Task 7
    #[command(before_help = EXPERIMENTAL)]
    PairIntersect(PairIntersectArgs),
    // P3b Task 8
    #[command(before_help = EXPERIMENTAL)]
    PairPair(PairPairArgs),
    // P3b Task 5
    #[command(before_help = EXPERIMENTAL)]
    Random(RandomArgs),
    #[command(before_help = EXPERIMENTAL)]
    Sample(SampleArgs),
    // P3b Task 6
    #[command(before_help = EXPERIMENTAL)]
    Shuffle(ShuffleArgs),
    #[command(before_help = EXPERIMENTAL)]
    Faidx(FaidxArgs),
    #[command(before_help = EXPERIMENTAL)]
    Dict(DictArgs),
    #[command(before_help = EXPERIMENTAL)]
    Flags(FlagsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Idxstats(IdxstatsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Flagstat(FlagstatArgs),
    #[command(before_help = EXPERIMENTAL)]
    Head(HeadArgs),
    #[command(before_help = EXPERIMENTAL)]
    Samples(SamplesArgs),
    #[command(before_help = EXPERIMENTAL)]
    Depth(DepthArgs),
    #[command(before_help = EXPERIMENTAL)]
    Bedcov(BedcovArgs),
    #[command(before_help = EXPERIMENTAL)]
    Reheader(ReheaderArgs),
}

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
/// Copies --fasta (plain, gzip or BGZF; no index needed) to the output
/// with every base under a -b interval replaced by --mask-char (default
/// N) or, with --soft, lowercased (letters only; other characters are
/// kept). Headers and line layout are kept exactly: line widths, CRLF
/// line ends, blank lines and a missing final line break. (bedtools
/// re-wraps lines to the first line's width, drops header descriptions
/// and cannot match names in a CRLF file; it also adds a final line break
/// that was missing and drops empty sequences, which gloci keeps.) The FASTA is read in pieces
/// of at most 1 MiB, so a chromosome of any size needs no memory beyond
/// the masks.
///
/// Masks come from any number of files in any order, in any interval
/// format (BED, bedGraph, GFF/GTF, VCF/BCF); overlapping masks are
/// merged so a base is masked once. A mask names a sequence by the
/// first word of its header; use --chrom-map when the names differ.
/// Masks on a chromosome the FASTA lacks are skipped with one warning.
/// A mask that reaches past the end of its sequence is an error, unless
/// --skip-invalid is given: then it is clipped to the sequence end (one
/// that starts at or past the end masks nothing) with one warning that
/// counts them. --mask-char takes one printable character; a space is
/// refused because it makes the FASTA unindexable.
///
/// Output goes to stdout or -o FILE (.fa, .fasta, .fna; a .gz name is
/// BGZF). --write-index also writes FILE.fai, and FILE.gzi for BGZF,
/// the same files gloci index writes.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci maskfasta --fasta hg38.fa -b repeats.bed -o hg38.masked.fa
  gloci maskfasta --fasta hg38.fa.gz -b rmsk.bed segdups.bed --soft -o soft.fa.gz --write-index
  gloci maskfasta --fasta ref.fa -b low-mappability.bed --mask-char X | head")]
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
/// --fasta REF.fa is read at random through REF.fa.fai (and REF.fa.gz.gzi
/// for BGZF); build them with gloci index. gloci never builds them itself.
/// The intervals (BED, GFF/GTF, VCF/BCF) are written in input order, one
/// sequence each, on a single line with the case of the FASTA kept and
/// memory that does not grow with the interval length. The header is
/// chrom:start-end with a 0-based start, as bedtools; --name gives
/// name::chrom:start-end and --name-only gives name (BED column 4, the
/// GFF feature type, or REF/ALT_ID for VCF). -s reverse-complements
/// minus-strand records with the IUPAC-aware, case-keeping table of
/// bedtools and adds (+), (-) or (.) to the header; --rna writes U for
/// the complement of A (as bedtools, it changes reverse-complemented
/// sequences only). --split joins the blocks of BED12 records; with -s a
/// minus-strand record is reverse-complemented as a whole. --tab writes
/// the header, TAB, the sequence; --bed-out writes the input record, TAB,
/// the sequence. VCF follows --vcf-span (--vcf-span ref is bedtools'
/// rule). A zero-length interval has no sequence: it is skipped with one
/// warning that counts them, as bedtools does. An interval on a
/// chromosome the FASTA lacks, or ending past its end, is an error; bedtools skips it
/// with a message and exits 0. --skip-invalid skips such intervals with one
/// warning carrying the count. -o FILE.fa.gz writes BGZF.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci getfasta --fasta ref.fa genes.bed
  gloci getfasta --fasta ref.fa genes.bed12 --split -s --name-only -o transcripts.fa
  gloci getfasta --fasta ref.fa.gz calls.vcf.gz --vcf-span ref --tab")]
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
/// --fasta REF.fa is read at random through REF.fa.fai (and REF.fa.gz.gzi
/// for BGZF); build them with gloci index. gloci never builds them itself.
/// Each interval (BED, GFF/GTF, VCF/BCF) is written as it was read, then
/// pct_at pct_gc num_A num_C num_G num_T num_N num_oth seq_len, preceded
/// by bedtools' header line (#1_usercol ...). The counts ignore case, U
/// counts as T and any other byte (IUPAC codes, -, *) is other; the
/// percentages divide by the interval length, are computed in f32 and
/// print with six decimals, as bedtools prints them. -s reverse-complements
/// minus-strand records first, so everything below describes the reverse
/// complement. --seq adds the sequence; --pattern P adds the number of
/// times P occurs (overlapping matches count, case-sensitive unless
/// --ignore-case). Memory does not grow with the interval length. A
/// zero-length interval is skipped with one warning that counts them, as
/// bedtools does. An interval on a chromosome the FASTA lacks, or ending
/// past its end, is an error; bedtools skips it with a message and exits
/// 0. --skip-invalid skips such intervals with one warning carrying the
/// count.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci nuc --fasta ref.fa genes.bed
  gloci nuc --fasta ref.fa genes.bed -s --seq -o genes.nuc.tsv
  gloci nuc --fasta ref.fa.gz calls.vcf.gz --vcf-span ref --pattern CG")]
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

/// Join BAM files by copying their compressed blocks, as samtools cat does.
///
/// The output header is the first file's, plus the @RG lines of the later
/// files whose ID it does not have yet, plus a gloci @PG line. The @SQ
/// lists of all files must be the same, in the same order: unlike
/// samtools, cat refuses files with different lists, because their
/// records would point at the wrong references. Every record is copied
/// as it is: only the block holding a file's header end is compressed
/// again, and the order of the records is not checked, so give sorted
/// files in sorted order. The input files may be named - (stdin, once).
/// BAM only: CRAM and SAM are refused (gloci view x.sam -o x.bam makes
/// BAM; samtools cat joins CRAM). The output is BAM, to -o or stdout.
/// --write-index also writes the .bai, for files that are already in
/// coordinate order.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci cat lane1.bam lane2.bam -o all.bam
  gloci cat part-*.bam > all.bam")]
pub struct CatArgs {
    /// BAM files to join, in output order (`-` for stdin).
    #[arg(value_name = "FILES", required = true)]
    pub files: Vec<String>,
    /// Write the BAM to FILE instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Also write a .bai index next to the output (.csi when a reference is longer than 2^29).
    /// The joined file must be coordinate-sorted.
    #[arg(long)]
    pub write_index: bool,
    #[command(flatten)]
    pub threads: ThreadArgs,
}

// P3c Task 9
const FAIDX_REGION_HELP: &str = "Another region, as the positional REGION (repeatable).";
const FAIDX_REGIONS_FILE_HELP: &str = "Read regions from a BED FILE (0-based, half-open); each prints as chrom:START-END with a 1-based start.";

/// Write regions of an indexed FASTA (samtools faidx).
///
/// REF.fa is read at random through REF.fa.fai (and REF.fa.gz.gzi for
/// BGZF); build them with gloci index. gloci never builds them itself
/// (samtools faidx does): without them, or with no region, this exits 2
/// naming gloci index REF. Each REGION is chr, chr:START, chr:START-,
/// chr:START-END, chr:-END (1-based, closed, commas allowed) or
/// {name}:START-END for a name containing a colon, read as htslib reads
/// it; --region adds more, and --regions-file reads BED lines (0-based,
/// half-open; samtools -r reads chr:START-END text instead). Each region
/// prints ">" and its text as typed, then its bases wrapped at
/// --line-width (default: that sequence's line width in the .fai; 0 is
/// one line), the case of the FASTA kept. A BED region is named
/// chrom:START-END with a 1-based start, and a zero-length one is
/// skipped with one warning that counts them. An end past the sequence
/// is cut there with a warning; a start past it writes the header alone.
/// A sequence the index lacks is an error before anything is written,
/// unless --skip-invalid skips that region entirely (samtools -c still
/// writes its header). --reverse-complement writes the reverse
/// complement, with samtools' table (IUPAC-aware, case kept), and
/// --mark-strand chooses the header suffix: rc (/rc on reverse-complemented
/// regions, the default), no, sign ((+) or (-)) or custom,POS,NEG.
/// Memory does not grow with the region length. -o FILE.fa.gz writes
/// BGZF, and --write-index also writes FILE.fai (and FILE.gzi).
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci faidx ref.fa chr1:1,000-2,000 chrM
  gloci faidx ref.fa chr2:500-900 --reverse-complement --mark-strand sign -o rc.fa
  gloci faidx ref.fa.gz --regions-file genes.bed --line-width 0 -o genes.fa.gz --write-index")]
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

// P3c Task 10
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

// P3c Task 1
/// `--column-names`: a `#`-prefixed line naming the columns comes first (samtools `-H` of
/// `depth` and `bedcov`, `-h` of `samples`).
#[derive(Args, Debug, Clone)]
pub struct ColumnNamesArg {
    /// Write a line of column names, starting with #, before the first row.
    #[arg(long)]
    pub column_names: bool,
}

/// Convert SAM FLAGS between numbers and names (samtools flags).
///
/// Each FLAGS argument is a decimal number, a 0x hexadecimal number, or
/// comma-separated names (case-insensitive), up to 0xFFFF. It prints
/// "0x<hex>", the decimal value and the names of the bits set, tab
/// separated, in bit order: PAIRED PROPER_PAIR UNMAP MUNMAP REVERSE
/// MREVERSE READ1 READ2 SECONDARY QCFAIL DUP SUPPLEMENTARY. A value with
/// none of them has an empty third field. With no argument it prints the
/// twelve flags with their descriptions. All arguments are checked before
/// anything is printed. A decimal with a leading zero (010) is an error,
/// because samtools reads it as octal 8 (00 to 07 mean the same either
/// way and are read), and so is trailing text (12abc), which samtools
/// reads as 12.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci flags 0x704 1796 secondary,dup
  gloci flags")]
pub struct FlagsArgs {
    /// A number (decimal or 0x hex) or comma-separated flag names. With none, the table of
    /// flags is printed.
    #[arg(value_name = "FLAGS")]
    pub flags: Vec<String>,
}

// P3c Task 2
/// Count the reads of each reference of a BAM, CRAM or SAM (samtools idxstats).
///
/// Writes one line per @SQ reference, in header order: name, length, mapped reads, unmapped
/// reads (tab separated), then "*", 0, 0 and the reads with no reference. A read that is
/// unmapped but placed on a reference counts as unmapped there. A BAM with a .bai or .csi
/// is answered from the index, without reading any record. Anything else (a BAM without
/// an index, a SAM, or a CRAM, whose .crai holds no counts) gets one warning and every
/// record is counted, in any order of the references; samtools stops at the first
/// reference that comes back. An index that lacks the counts of a reference also gets a
/// warning and a scan. From a CRAM only the fields counted are decoded, as samtools
/// does; a CRAM that needs a reference sequence still needs --reference, though
/// samtools reads such a file without one. Mapped and unmapped counts follow
/// FLAGS bit 0x4 only; no read filters apply. The output is the same through the index
/// and by a scan, at any number of threads.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci idxstats reads.bam
  gloci idxstats reads.cram --reference ref.fa
  gloci idxstats reads.sam -t 1")]
pub struct IdxstatsArgs {
    /// The BAM, CRAM or SAM file (- is standard input, which is always scanned).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The FASTA a CRAM input was made with, with its .fai (and .gzi when BGZF). Not needed
    /// for a CRAM written without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the counts to FILE instead of stdout; .gz compresses them.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}

// P3c Task 3
/// Count the reads of a BAM, CRAM or SAM by FLAGS category (samtools flagstat).
///
/// Writes the counts of total, primary, secondary and supplementary reads, duplicates,
/// mapped reads, paired reads, read 1 and read 2, properly paired reads, singletons and
/// mates on another reference, as "QC-passed + QC-failed" (FLAGS bit 0x200), with
/// percentages where samtools has them. -O picks the layout: text (the default), tsv or
/// json; with no -O, a -o ending in .json or .tsv picks it, and a -o ending in .gz is
/// BGZF-compressed (stats.json.gz is compressed JSON). A percentage of nothing is
/// N/A (null in JSON). Percentages are computed as samtools computes them, in 32-bit
/// floating point, so 153 of 160 is 95.63% and not 95.62%. Every record counts, unmapped
/// and unplaced ones too, and no read filters apply. From a CRAM only the fields counted
/// are decoded, as samtools does; a CRAM that needs a reference sequence still needs
/// --reference, though samtools reads such a file without one. The output is the same at
/// any number of threads.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci flagstat reads.bam
  gloci flagstat reads.bam -O json -o stats.json
  gloci flagstat reads.cram --reference ref.fa -O tsv")]
pub struct FlagstatArgs {
    /// The BAM, CRAM or SAM file (- is standard input).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Write the counts to FILE instead of stdout; .json and .tsv choose that layout, and .gz
    /// compresses it.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Layout of the counts: text, tsv or json.
    #[arg(short = 'O', long, value_name = "FMT")]
    pub output_format: Option<String>,
    /// The FASTA a CRAM input was made with, with its .fai (and .gzi when BGZF). Not needed
    /// for a CRAM written without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    #[command(flatten)]
    pub threads: ThreadArgs,
}

// P3c Task 5
/// List the samples of the read groups of BAM, CRAM or SAM files (samtools samples).
///
/// Reads only the headers. For each file it prints one line per distinct value of the --tag
/// (default SM) over the @RG lines, in the order the header lists them, then the path; a
/// file with no such value gets "." as its sample. --show-index adds Y or N for whether the
/// file has an index. --match-reference FASTA (repeatable) adds the first FASTA whose .fai
/// has the same names and lengths in the same order as the @SQ lines, or ".", trying the
/// FASTA given last first; every FASTA needs its .fai (gloci index REF). --column-names
/// writes a line such as "#SM PATH INDEX REFERENCE" first. samtools prints the samples in
/// hash order, so a file with several samples lists them in a different order.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci samples cohort.bam
  gloci samples a.bam b.cram --show-index --column-names
  gloci samples reads.bam --match-reference hg38.fa --match-reference hs37d5.fa --tag LB")]
pub struct SamplesArgs {
    /// BAM, CRAM or SAM files.
    #[arg(value_name = "FILE", required = true)]
    pub files: Vec<String>,
    /// The @RG tag to list the values of: two characters.
    #[arg(long, value_name = "TAG", default_value = "SM")]
    pub tag: String,
    /// Add a column: Y when the file has an index, N when it has none.
    #[arg(long)]
    pub show_index: bool,
    /// A FASTA (with its .fai) whose sequences are to be compared with each file's @SQ lines;
    /// adds a column with the matching FASTA. Repeatable.
    #[arg(long, value_name = "FASTA")]
    pub match_reference: Vec<PathBuf>,
    #[command(flatten)]
    pub column_names: ColumnNamesArg,
    /// Write the lines to FILE instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}

// P3c Task 6
/// Depth of coverage of alignments at each position (samtools depth).
///
/// Reads BAM, CRAM or SAM files (CRAM needs --reference) in coordinate order and prints
/// "chrom, 1-based position, depth" for each position covered by a passing read in any
/// input, with one depth column per input in the order given, references in @SQ order.
/// Inputs must have identical @SQ lists (exit 1, naming the first difference; samtools
/// compares reference ids unchecked). A read counts when none of --exclude-flags is set
/// (default UNMAP,SECONDARY,QCFAIL,DUP: this command's default, where other commands
/// exclude nothing; --exclude-flags 0 counts everything), all of --require-flags are, its
/// MAPQ is at least --min-mapq and its length (SEQ length less the soft clips at its two
/// ends; with SEQ *, the bases of M, I, = and X) is at least --min-read-length. Each M, =
/// and X base counts when its quality is at least --min-baseq (a missing quality passes,
/// and so does every base of a read with SEQ *). D bases count only with --count-deletions,
/// each gated by the quality of the next base in the read; N never counts. Positions a
/// read spans but none of its bases count (deletions, introns, filtered bases) print 0.
/// --all-positions adds every position of each reference that has a passing read;
/// --all-references adds every position of every reference. --region (repeatable) and
/// --regions-file (BED) limit the rows to their union, each position once, read through the
/// index when there is one (else the input is scanned, with a warning); with regions, the
/// --all-* options print every position in them. --column-names writes "#CHROM POS" and
/// the paths first. A read past the end of its reference is cut there with one warning,
/// where samtools prints rows beyond it. --skip-mate-overlap counts the bases two mates of
/// a pair share once. An input that goes back in (reference, position)
/// is an error telling you to run gloci sort. The inputs are read as one stream at every
/// -t, with the same output; threads decode BGZF and CRAM, shared among the inputs. Memory
/// per input follows the longest read, not the reference.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci depth sample.bam --skip-mate-overlap -o depth.tsv
  gloci depth a.bam b.bam --column-names --min-mapq 20 --min-baseq 20
  gloci depth sample.bam --region chr1:100,000-200,000 --all-positions
  gloci depth sample.bam --regions-file targets.bed --count-deletions")]
pub struct DepthArgs {
    /// BAM, CRAM or SAM files in coordinate order (- is standard input), with identical
    /// @SQ lists; one depth column each.
    #[arg(value_name = "FILE", required = true)]
    pub files: Vec<String>,
    /// Also print every position of each reference that has a passing read (samtools -a).
    #[arg(long)]
    pub all_positions: bool,
    /// Also print every position of every reference, reads or not (samtools -aa).
    #[arg(long)]
    pub all_references: bool,
    /// Count deleted bases (D), each gated by the quality of the next base in the read.
    #[arg(long)]
    pub count_deletions: bool,
    /// Count a base only when its quality is at least N. A missing quality (QUAL *) passes, and
    /// so does every base of a read with SEQ * (samtools reads the bytes after the CIGAR as
    /// qualities there).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_baseq: u8,
    /// Skip reads shorter than N bases (SEQ length less soft clips at the two ends).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_read_length: u64,
    /// Skip reads with mapping quality below N.
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_mapq: u8,
    /// Keep only reads that have all of these flags: a number (0x2) or names such as
    /// paired,proper_pair.
    #[arg(long, value_name = "FLAGS")]
    pub require_flags: Option<String>,
    /// Drop reads that have any of these flags: a number or names. As for the other commands,
    /// but this command's default is UNMAP,SECONDARY,QCFAIL,DUP (samtools depth's), where
    /// theirs is none; 0 counts every read.
    #[arg(
        long,
        value_name = "FLAGS",
        default_value = "UNMAP,SECONDARY,QCFAIL,DUP"
    )]
    pub exclude_flags: String,
    /// The FASTA a CRAM input was made with, with its .fai. Not needed for a CRAM written
    /// without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    /// Count each base of an overlapping mate pair once, per input: a read whose mate was
    /// already seen starts counting where the mate ends (samtools -s). Reads with an unknown mate
    /// position are never held, so unlike samtools a later alignment with the same name is not
    /// clipped by a stale entry.
    #[arg(long)]
    pub skip_mate_overlap: bool,
    #[command(flatten)]
    pub region: RegionArgs,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub columns: ColumnNamesArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the rows to FILE instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}

// P3c Task 8
/// Depth of coverage summed over each BED region (samtools bedcov).
///
/// For each line of BED, in any order, writes the line as it is and then, for each FILE (BAM
/// or CRAM, each with a .bai, .csi or .crai; exit 2 without one), the sum of the depth over
/// the region's bases: one column per file, files in the order given. --bases-at-depth N
/// adds one column per file of the bases whose depth is at least N, and --read-count one of
/// the reads that overlap the region (a zero-length region has no bases and counts no reads,
/// as samtools 1.21 does). Depth counts M, = and X bases always, D bases unless --split-deletions and
/// N (intron) bases unless --split, so samtools -j is --split --split-deletions. A read counts
/// when none of --exclude-flags is set (default UNMAP,SECONDARY,QCFAIL,DUP, as for depth;
/// --exclude-flags 0 counts everything but unmapped reads, which never count), all of
/// --require-flags are, and its MAPQ is at least --min-mapq. Comment, track, browser and blank
/// lines are skipped. A line on a chromosome the files lack is an error naming the line,
/// or skipped with one counted warning under --skip-invalid; samtools prints an error,
/// carries on and exits 2. Files must have identical @SQ lists. --column-names writes
/// "#chrom start end" and one column name per file and kind first. A read past the end of
/// its reference counts up to the end. Regions are read in batches of 256 lines on up to -t
/// threads, written in BED order. Memory per thread follows the longest read, not the region.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci bedcov targets.bed sample.bam
  gloci bedcov targets.bed a.bam b.bam --bases-at-depth 10 --read-count --column-names
  gloci bedcov exons.bed sample.cram --reference ref.fa --split --split-deletions -t 8")]
pub struct BedcovArgs {
    /// BED regions, in any order (- is standard input).
    #[arg(value_name = "BED")]
    pub bed: String,
    /// BAM or CRAM files with an index (.bai, .csi or .crai) and identical @SQ lists; one
    /// column each.
    #[arg(value_name = "FILE", required = true)]
    pub files: Vec<String>,
    /// Do not count N (intron) bases (samtools -j, with --split-deletions).
    #[arg(long)]
    pub split: bool,
    /// Do not count D (deleted) bases (samtools -j, with --split).
    #[arg(long)]
    pub split_deletions: bool,
    /// Also write, per file, the number of bases of the region whose depth is at least N
    /// (samtools -d; samtools then also leaves D and N bases out of both columns).
    #[arg(long, value_name = "N")]
    pub bases_at_depth: Option<u32>,
    /// Also write, per file, the number of reads that overlap the region (samtools -c).
    #[arg(long)]
    pub read_count: bool,
    /// Skip reads with mapping quality below N.
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_mapq: u8,
    /// Keep only reads that have all of these flags: a number (0x2) or names such as
    /// paired,proper_pair.
    #[arg(long, value_name = "FLAGS")]
    pub require_flags: Option<String>,
    /// Drop reads that have any of these flags: a number or names. As for the other commands,
    /// but this command's default is UNMAP,SECONDARY,QCFAIL,DUP (samtools bedcov's), where
    /// theirs is none; 0 counts every mapped read.
    #[arg(
        long,
        value_name = "FLAGS",
        default_value = "UNMAP,SECONDARY,QCFAIL,DUP"
    )]
    pub exclude_flags: String,
    /// The FASTA a CRAM input was made with, with its .fai. Not needed for a CRAM written
    /// without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub columns: ColumnNamesArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the lines to FILE instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}
// P3c Task 11
/// Replace the header of a BAM or CRAM file without decoding its records (samtools reheader).
///
/// HEADER is SAM header text, or a SAM, BAM or CRAM file whose header is
/// taken (as samtools does); - reads standard input. IN must be a regular
/// BAM or CRAM file: standard input is refused (gloci view - -o x.bam makes
/// a file) and so is SAM. The result goes to -o or stdout, in IN's format;
/// -o may name IN, which is then replaced atomically, and a failure leaves
/// it as it was. A gloci @PG line is added to the new header.
///
/// Records point at their reference by its position in the @SQ list, so the
/// new list must keep every old reference in the same place with the same
/// LN. Names (SN:) may change, and @SQ lines may be added after the last
/// one. For CRAM an M5: that both lists have must be equal. Anything else
/// is an error naming both lines, and nothing is written (samtools
/// checks nothing and would attach reads to the wrong contigs; use gloci
/// view --chrom-map to rename references safely).
///
/// BAM: the new header is compressed, and the rest of IN is copied block
/// by block: only the block that held the old header's end is compressed
/// again. CRAM (3.0 and 3.1; 2.x is refused): the file definition, a new
/// header container, then every byte of IN from its first data container. The old and new headers and
/// one BGZF block are all that is held in memory.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci reheader new-header.sam x.bam -o y.bam
  gloci reheader other.bam x.bam -o x.bam
  gloci view --header-only x.bam -O sam | sed 's/SN:chr/SN:/' | gloci reheader - x.bam -o y.bam")]
pub struct ReheaderArgs {
    /// The new header: SAM header text, or a SAM, BAM or CRAM file whose header is taken (- is
    /// standard input).
    #[arg(value_name = "HEADER")]
    pub header: String,
    /// The BAM or CRAM file to give the header to (a regular file, not standard input).
    #[arg(value_name = "IN")]
    pub input: PathBuf,
    /// Write the result to FILE instead of stdout (FILE may be IN).
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// BGZF compression level of a BAM's new header, 0-12.
    #[arg(long, value_name = "N", default_value = "6", value_parser = parse_level)]
    pub compression_level: u8,
    #[command(flatten)]
    pub threads: ThreadArgs,
}
