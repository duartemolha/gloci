use super::shared::{ChromMapArg, InvalidArg, OutputArgs, RegionArgs, ThreadArgs, parse_level};
use clap::Args;
use std::path::PathBuf;

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
#[command(about = "Join BAM files by copying their blocks (samtools cat)")]
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
#[command(about = "Count the reads of each reference (samtools idxstats)")]
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
#[command(about = "Count reads by FLAGS category (samtools flagstat)")]
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
#[command(about = "List the samples of BAM, CRAM or SAM files (samtools samples)")]
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
/// Reads BAM, CRAM or SAM files (CRAM needs --reference) in coordinate order and prints "chrom,
/// 1-based position, depth" for each position covered by a passing read in any input, with one depth
/// column per input in the order given, references in @SQ order. Inputs must have identical @SQ lists
/// (exit 1, naming the first difference; samtools compares reference ids unchecked).
#[derive(Args)]
#[command(about = "Depth of alignments at each position (samtools depth)")]
#[command(after_long_help = "\
Examples:
  gloci depth sample.bam --skip-mate-overlap -o depth.tsv
  gloci depth a.bam b.bam --column-names --min-mapq 20 --min-baseq 20
  gloci depth sample.bam --region chr1:100,000-200,000 --all-positions
  gloci depth sample.bam --regions-file targets.bed --count-deletions

Counting:
A read counts when none of --exclude-flags is set (default UNMAP,SECONDARY,QCFAIL,DUP: \
this command's default, where other commands exclude nothing; --exclude-flags 0 counts \
everything), all of --require-flags are, its MAPQ is at least --min-mapq and its length \
(SEQ length less the soft clips at its two ends; with SEQ *, the bases of M, I, = and X) \
is at least --min-read-length. Each M, = and X base counts when its quality is at least \
--min-baseq (a missing quality passes, and so does every base of a read with SEQ *). D \
bases count only with --count-deletions, each gated by the quality of the next base in the \
read; N never counts. Positions a read spans but none of its bases count (deletions, \
introns, filtered bases) print 0.

Positions and regions:
--all-positions adds every position of each reference that has a passing read; \
--all-references adds every position of every reference. --region (repeatable) and \
--regions-file (BED) limit the rows to their union, each position once, read through the \
index when there is one (else the input is scanned, with a warning); with regions, the \
--all-* options print every position in them. --column-names writes \"#CHROM POS\" and the \
paths first.

Reads and mates:
A read past the end of its reference is cut there with one warning, where samtools prints \
rows beyond it. --skip-mate-overlap counts the bases two mates of a pair share once.

Order, memory and threads:
An input that goes back in (reference, position) is an error telling you to run gloci \
sort. The inputs are read as one stream at every -t, with the same output; threads decode \
BGZF and CRAM, shared among the inputs. Memory per input follows the longest read, not the \
reference.")]
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
/// For each line of BED, in any order, writes the line as it is and then, for each FILE (BAM or CRAM,
/// each with a .bai, .csi or .crai; exit 2 without one), the sum of the depth over the region's
/// bases: one column per file, files in the order given. --bases-at-depth N adds one column per file
/// of the bases whose depth is at least N, and --read-count one of the reads that overlap the region
/// (a zero-length region has no bases and counts no reads, as samtools 1.21 does).
#[derive(Args)]
#[command(about = "Depth summed over each BED region (samtools bedcov)")]
#[command(after_long_help = "\
Examples:
  gloci bedcov targets.bed sample.bam
  gloci bedcov targets.bed a.bam b.bam --bases-at-depth 10 --read-count --column-names
  gloci bedcov exons.bed sample.cram --reference ref.fa --split --split-deletions -t 8

Depth and reads:
Depth counts M, = and X bases always, D bases unless --split-deletions and N (intron) \
bases unless --split, so samtools -j is --split --split-deletions. A read counts when none \
of --exclude-flags is set (default UNMAP,SECONDARY,QCFAIL,DUP, as for depth; \
--exclude-flags 0 counts everything but unmapped reads, which never count), all of \
--require-flags are, and its MAPQ is at least --min-mapq. A read past the end of its \
reference counts up to the end.

Regions:
Comment, track, browser and blank lines are skipped. A line on a chromosome the files lack \
is an error naming the line, or skipped with one counted warning under --skip-invalid; \
samtools prints an error, carries on and exits 2. Files must have identical @SQ lists. \
--column-names writes \"#chrom start end\" and one column name per file and kind first.

Threads and memory:
Regions are read in batches of 256 lines on up to -t threads, written in BED order. Memory \
per thread follows the longest read, not the region.")]
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
/// HEADER is SAM header text, or a SAM, BAM or CRAM file whose header is taken (as samtools does); -
/// reads standard input. IN must be a regular BAM or CRAM file: standard input is refused (gloci view -
/// -o x.bam makes a file) and so is SAM. The result goes to -o or stdout, in IN's format; -o may name
/// IN, which is then replaced atomically, and a failure leaves it as it was. A gloci @PG line is added
/// to the new header.
#[derive(Args)]
#[command(about = "Replace the header of a BAM or CRAM file (samtools reheader)")]
#[command(after_long_help = "\
Examples:
  gloci reheader new-header.sam x.bam -o y.bam
  gloci reheader other.bam x.bam -o x.bam
  gloci view --header-only x.bam -O sam | sed 's/SN:chr/SN:/' | gloci reheader - x.bam -o y.bam

Header rules:
Records point at their reference by its position in the @SQ list, so the new list must \
keep every old reference in the same place with the same LN. Names (SN:) may change, and \
@SQ lines may be added after the last one. For CRAM an M5: that both lists have must be \
equal. Anything else is an error naming both lines, and nothing is written (samtools \
checks nothing and would attach reads to the wrong contigs; use gloci view --chrom-map to \
rename references safely).

How the file is written:
BAM: the new header is compressed, and the rest of IN is copied block by block: only the \
block that held the old header's end is compressed again. CRAM (3.0 and 3.1; 2.x is \
refused): the file definition, a new header container, then every byte of IN from its \
first data container. The old and new headers and one BGZF block are all that is held in \
memory.")]
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

/// `--fraction` of `tag`: a fraction as the other commands take it (`0.5`, `50%`), held as the
/// 32-bit float bedtools compares with.
fn parse_tag_fraction(s: &str) -> Result<f32, String> {
    gloci_core::Fraction::parse(s).map_err(|e| e.to_string())?;
    let (body, scale) = match s.strip_suffix('%') {
        Some(b) => (b, 100.0),
        None => (s, 1.0),
    };
    let v: f64 = body
        .parse()
        .map_err(|_| format!("invalid fraction {s:?}"))?;
    Ok((v / scale) as f32)
}

/// Tag each alignment with the features it overlaps, as bedtools tag does.
///
/// Writes every record of ALN (BAM, CRAM or SAM) in input order, unchanged but for one Z tag
/// (--tag, default YB) added to each mapped read that overlaps features of the -b files (BED or
/// GFF, held in memory). For each -b file with a hit, in -b order, the tag holds the file's
/// label (--labels), or with --names, --scores or --intervals the name, score or whole interval
/// of each feature hit, joined by commas; the parts of the files are joined by semicolons. A read
/// overlaps a feature when they share bases of the read's span, from POS over the reference bases
/// of its CIGAR (introns and deletions included, clips not), and the bases shared are at least
/// --fraction of that span. The output is in the input's format (a CRAM gives BAM); -o and -O
/// choose BAM or SAM.
#[derive(Args)]
#[command(about = "Tag alignments with the features they overlap (bedtools tag)")]
#[command(after_long_help = "\
Examples:
  gloci tag reads.bam -b genes.bed --labels genes -o tagged.bam
  gloci tag reads.bam -b exons.bed introns.bed --labels exon,intron -s -O sam
  gloci tag reads.cram --reference ref.fa -b genes.gff3 --names --tag XG -o tagged.bam

Reads and overlap:
Unmapped reads (FLAGS 0x4) are never tagged; secondary and supplementary alignments are. The \
span counts M, D, N, = and X bases, so a feature inside an intron or a deletion overlaps the \
read and one under a soft clip does not; a read with no such base spans nothing and is never \
tagged. The shared bases are divided by the span and compared with --fraction (default: one \
base) in 32-bit floating point, as bedtools compares them. -s keeps the features on the read's \
strand and -S those not on it, so under -S a feature of unknown strand (. or no strand column) \
counts, as in bedtools tag; intersect's -S never matches one.

Features:
BED (and bedGraph, read as BED) gives the name, score and strand of columns 4, 5 and 6, and a \
zero-length record overlaps as [start-1, end+1), which --intervals shows. GFF and GTF give the \
feature type (column 3), score and strand, from start-1 to end. A column a file lacks gives an \
empty name or score. Comment, track and browser lines are skipped. Every feature is held in \
memory, about 110 bytes each plus the text the tag takes from it.

The tag:
The hits of one file come in bedtools' bin order: features that fit in one 16 kb bin first, by \
bin and then in file order, then those of each larger bin size. Duplicates are kept. An empty \
part still counts, so a value can be empty or end in a semicolon. The tag is a string added \
after the record's other tags, also when the record already has one of that name (bedtools adds \
a second one); nothing else of the record changes.")]
pub struct TagArgs {
    /// The BAM, CRAM or SAM file to tag (- is standard input).
    #[arg(value_name = "ALN")]
    pub input: String,
    /// Feature files (BED or GFF), held in memory; their parts of the tag come in this order
    /// (-b x y or -b x -b y).
    #[arg(
        short = 'b',
        long = "b",
        value_name = "FILE",
        required = true,
        num_args = 1..,
        action = clap::ArgAction::Append
    )]
    pub b: Vec<String>,
    /// One label per -b file, comma-separated: the tag's value (bedtools -labels), and the start
    /// of each --intervals entry.
    #[arg(
        long,
        value_name = "L,...",
        value_delimiter = ',',
        required_unless_present_any = ["names", "scores"],
        conflicts_with_all = ["names", "scores"]
    )]
    pub labels: Option<Vec<String>>,
    /// The tag to add: two characters, a letter and then a letter or a digit (bedtools -tag).
    #[arg(long, value_name = "XX", default_value = "YB")]
    pub tag: String,
    /// The tag holds the names of the features hit (bedtools -names).
    #[arg(long, conflicts_with_all = ["scores", "intervals"])]
    pub names: bool,
    /// The tag holds the scores of the features hit, as written (bedtools -scores).
    #[arg(long, conflicts_with = "intervals")]
    pub scores: bool,
    /// The tag holds label:chrom:start-end,name,score,strand for each feature hit (bedtools
    /// -intervals); needs --labels.
    #[arg(long)]
    pub intervals: bool,
    /// Only features on the read's strand.
    #[arg(short = 's', long = "same-strand", conflicts_with = "opposite_strand")]
    pub same_strand: bool,
    /// Only features not on the read's strand; one of unknown strand counts.
    #[arg(short = 'S', long = "opposite-strand")]
    pub opposite_strand: bool,
    /// The least share of the read's span a feature must cover, e.g. 0.5 or 50% (default: one
    /// base; bedtools -f).
    #[arg(short = 'f', long, value_name = "FRAC", value_parser = parse_tag_fraction)]
    pub fraction: Option<f32>,
    /// The FASTA a CRAM input was made with, with its .fai. Not needed for a CRAM written
    /// without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
    #[command(flatten)]
    pub out: OutputArgs,
    #[command(flatten)]
    pub threads: ThreadArgs,
}
