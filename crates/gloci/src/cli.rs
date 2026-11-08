//! clap definitions: global flags, shared option groups and the subcommands.
use clap::builder::TypedValueParser as _;
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// Parse a size such as `512M` or `2G`: an integer with a required K/M/G/T
/// suffix in binary units (powers of 1024). `10` and `1.5G` are errors.
pub fn parse_size(s: &str) -> Result<u64, String> {
    let bad = || {
        format!("invalid size {s:?}: use an integer with a K, M, G or T suffix, e.g. 512M or 2G")
    };
    let last = s.chars().last().ok_or_else(bad)?;
    let shift = match last.to_ascii_uppercase() {
        'K' => 10,
        'M' => 20,
        'G' => 30,
        'T' => 40,
        _ => return Err(bad()),
    };
    let digits = &s[..s.len() - 1];
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad());
    }
    let n: u64 = digits.parse().map_err(|_| bad())?;
    n.checked_mul(1 << shift)
        .ok_or_else(|| format!("size {s:?} is too large"))
}

fn parse_threads(s: &str) -> Result<NonZeroUsize, String> {
    s.parse::<NonZeroUsize>()
        .map_err(|_| format!("invalid thread count {s:?}: use an integer of at least 1"))
}

fn parse_level(s: &str) -> Result<u8, String> {
    match s.parse::<u8>() {
        Ok(n) if n <= 12 => Ok(n),
        _ => Err(format!("invalid compression level {s:?}: use 0-12")),
    }
}

fn parse_comment(s: &str) -> Result<u8, String> {
    match s.as_bytes() {
        [b] if b.is_ascii() => Ok(*b),
        _ => Err(format!(
            "invalid comment character {s:?}: use one ASCII character"
        )),
    }
}

#[derive(Parser)]
#[command(
    name = "gloci",
    version,
    about = "Genomic interval toolkit",
    arg_required_else_help = true
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
    /// List the built-in genomes that -g accepts, or print one.
    ///
    /// Without a subcommand, prints one line per built-in genome with the
    /// columns id, species, assembly, naming (ucsc or ensembl), sequences
    /// and total_bp, separated by tabs. The lengths come from UCSC
    /// chrom.sizes files and from the Ensembl core databases; the sources
    /// are recorded in data/genomes/SOURCES.toml in the source tree.
    #[command(after_long_help = "\
Examples:
  gloci genomes
  gloci genomes show hg38
  gloci sort peaks.bed -g hg38 > peaks.sorted.bed")]
    Genomes(GenomesArgs),
    /// Build a .tbi/.csi index for a sorted BGZF file.
    ///
    /// BED, GFF/GTF, VCF and BCF files are recognised by content. Any other
    /// tab-separated file is indexed like tabix does when you name its
    /// columns with -s (sequence), -b (begin) and -e (end, optional). The
    /// columns are numbered from 1. Without --zero-based the begin column is
    /// 1-based and a feature covers [begin-1, end]; with it, [begin, end).
    /// Without -e a feature covers one base. Lines starting with --comment
    /// (default #) and the first --skip-lines lines are not indexed; both
    /// are recorded in the index. Chromosomes must be contiguous and sorted
    /// by begin. The index is a .tbi, or a .csi when an end passes 2^29.
    #[command(after_long_help = "\
Examples:
  gloci index calls.vcf.gz
  gloci index peaks.bed.gz -t 8
  gloci index scores.tsv.gz -s 1 -b 2 -e 3 --zero-based
  gloci index sites.txt.gz -s 1 -b 2 --skip-lines 1 --comment %")]
    Index(IndexArgs),
    /// Sort BED/bedGraph/GFF/VCF/BCF files with bounded memory.
    ///
    /// Records are ordered by chromosome, then start, then end. The
    /// chromosome order comes from -g, else from the ##contig lines of a
    /// VCF/BCF header, else from --chrom-order (karyotype by default).
    /// Records that tie keep their input order unless --full-key is given.
    /// BED header lines are kept only with --header; the leading ## lines of
    /// a GFF and the header of a VCF/BCF are always kept. VCF/BCF records are
    /// ordered by POS and REF end whatever --vcf-span says, so the output is
    /// tabix-compatible. The output format is the input format: -o must
    /// match it. Use -o with the same path as an input to sort a file in
    /// place.
    ///
    /// Sorting spills to --tmp-dir when the input is larger than --max-mem.
    /// That directory must exist and be writable; it is checked before any
    /// input is read.
    #[command(after_long_help = "\
Examples:
  gloci sort peaks.bed > peaks.sorted.bed
  gloci sort calls.vcf.gz -o calls.sorted.vcf.gz --write-index
  gloci sort huge.bed --max-mem 8G --tmp-dir /scratch -t 16 -o huge.sorted.bed.gz
  gloci sort --check peaks.sorted.bed")]
    Sort(SortArgs),
    /// Print the records in regions, stream or recompress a file, convert
    /// VCF to BCF and back, or print only the header.
    ///
    /// With no regions the records are written as read (comment and track
    /// lines after the first record and blank lines are dropped, and CRLF
    /// line ends become LF), so view can recompress a file (-o FILE.gz is
    /// BGZF) or index it (--write-index).
    /// VCF and BCF convert to each other with -O or the -o extension;
    /// the header is kept as it is, plus a ##glociCommand line. BED and GFF
    /// header lines are kept by default. Conversions between other formats are not
    /// available yet (they will be gloci convert).
    ///
    /// With regions, view prints the records whose footprint (per
    /// --vcf-span) overlaps a region, each once, in file order. A point
    /// (a zero-length BED record, an insertion) overlaps a region only when
    /// strictly inside it, unless --inclusive-points is given. Regions are
    /// read through a .tbi/.csi index when the file has one, and otherwise
    /// by scanning the whole file, with a warning. Either way the order of
    /// every record read is checked, so an unsorted file is an error, never
    /// a wrong answer. With --assume-sorted, an uncompressed BED file
    /// without an index is binary-searched instead (like bedextract); that
    /// trusts the file to be sorted. Index the file (gloci sort FILE -o
    /// FILE.gz --write-index) for speed and checking both.
    #[command(after_long_help = "\
Examples:
  gloci view calls.vcf.gz chr1:1,000,000-2,000,000 chrX
  gloci view peaks.bed --regions-file targets.bed --inclusive-points
  gloci view calls.bcf -O vcf | less
  gloci view calls.vcf.gz -o calls.bcf --write-index")]
    View(ViewArgs),
    /// Validate files and print a one-screen summary of each.
    ///
    /// check reads every file to the end, keeps going past problems and
    /// reports the record count, the chromosomes, whether the records are
    /// sorted, and the first 10 errors and warnings with their line numbers.
    /// Errors: records that do not parse, records out of order, an end past
    /// the chromosome length or a chromosome missing from -g, bad BED12
    /// blocks, a bedGraph value that is not a number, a VCF line with the
    /// wrong number of columns. Warnings: exact duplicate lines at the same
    /// position, and blank lines. Exit status 1 when any file has errors.
    /// The order is checked against -g, else the ##contig lines of a
    /// VCF/BCF, else --chrom-order (karyotype by default).
    #[command(after_long_help = "\
Examples:
  gloci check peaks.bed
  gloci check calls.vcf.gz -g hg38.fa.fai
  gloci check a.bed b.bed --chrom-order lex")]
    Check(CheckArgs),
    /// Report the records of A that overlap, or do not overlap, records of B.
    ///
    /// The default mode is any: A's records that overlap at least one B
    /// record are written as read, in A's format. --mode picks another
    /// report: none (A's records with no overlap), pairs (one TSV row per
    /// overlapping A and B pair; --overlap-bases adds the shared bases and
    /// --keep-unmatched adds a row with placeholder B columns for each A
    /// without one), clip (A cut to the shared span, once per overlap) and
    /// count (A with its number of overlaps, summed over all -b files).
    /// With several -b files, pairs rows carry a label column (--names,
    /// else the file basename) before B's columns; unmatched rows get `.`.
    ///
    /// Overlap follows spec section 4.3: records overlap when they share a
    /// base, a point overlaps an interval when strictly inside it (or on
    /// its edge with --inclusive-points), and every threshold (-f, -F, -r,
    /// -e, --min-bases) is met when either side is a point. Thresholds are
    /// measured against B's own span; they cannot be combined with -w,
    /// --left or --right. A VCF/BCF record covers the bases its ALT
    /// alleles change (--vcf-span). --match alleles also requires some ALT
    /// of A to equal some ALT of B after trimming (VCF/BCF inputs only).
    ///
    /// A and every B must be sorted in one chromosome order, which is
    /// checked as they are read; a violation is an error with a fix. With
    /// --unsorted, any order works, A is still streamed and reported in its
    /// own order, and every -b file is loaded into memory, so memory grows
    /// with the size of B. Pairs and count write TSV (a VCF A gets a
    /// #CHROM line naming the columns); clip does not accept VCF/BCF.
    #[command(after_long_help = "\
Examples:
  gloci intersect -a peaks.bed -b genes.bed
  gloci intersect -a calls.vcf.gz -b exons.bed --mode none -o off_target.vcf.gz
  gloci intersect -a a.bed -b x.bed y.bed --mode pairs --overlap-bases --names X,Y
  gloci intersect -a reads.bed -b targets.bed -f 50% -s --mode count --unsorted")]
    Intersect(IntersectArgs),
}

/// What `intersect` reports.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum IntersectMode {
    /// A's records with an overlap (default).
    Any,
    /// A's records without one.
    #[value(name = "none")]
    NoHit,
    /// One TSV row per overlapping pair.
    Pairs,
    /// A cut to the shared span, once per overlap.
    Clip,
    /// A with its number of overlaps.
    Count,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum IntersectMatch {
    Overlap,
    Alleles,
}

#[derive(Args)]
pub struct IntersectArgs {
    #[command(flatten)]
    pub a: AArg,
    #[command(flatten)]
    pub b: BArgs,
    /// What to report: any, none, pairs, clip or count.
    #[arg(long, value_enum, default_value = "any", value_name = "MODE")]
    pub mode: IntersectMode,
    /// With --mode pairs: add the number of shared bases to each row.
    #[arg(long)]
    pub overlap_bases: bool,
    /// With --mode pairs: also write each A without a hit, with placeholder B columns.
    #[arg(long)]
    pub keep_unmatched: bool,
    /// What makes a hit count: overlap, or overlap plus a shared ALT allele
    /// (alleles; VCF/BCF inputs only).
    #[arg(
        long = "match",
        value_enum,
        default_value = "overlap",
        value_name = "HOW"
    )]
    pub matching: IntersectMatch,
    #[command(flatten)]
    pub o: OverlapArgs,
    #[command(flatten)]
    pub w: WindowArgs,
    #[command(flatten)]
    pub u: UnsortedArg,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

#[derive(Args)]
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
    pub format: FormatArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
}

#[derive(Args)]
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
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

#[derive(Args)]
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
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

#[derive(Args)]
#[command(args_conflicts_with_subcommands = true)]
pub struct GenomesArgs {
    #[command(subcommand)]
    pub action: Option<GenomesAction>,
}

#[derive(Subcommand)]
pub enum GenomesAction {
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
        /// Genome id (case-insensitive), as listed by gloci genomes.
        name: String,
    },
}

#[derive(Args)]
pub struct IndexArgs {
    /// Sorted, BGZF-compressed BED/GFF/VCF/BCF file to index.
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
    #[command(flatten)]
    pub threads: ThreadArgs,
}

/// Output destination and encoding.
#[derive(Args, Debug, Clone)]
pub struct OutputArgs {
    /// Write output to FILE (format from its extension) instead of stdout. A
    /// device or pipe such as /dev/null or >(cmd) is written to directly.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Output format when writing to stdout.
    #[arg(short = 'O', long, value_name = "FMT")]
    pub output_format: Option<String>,
    /// Also write a .tbi/.csi index next to the output.
    #[arg(long)]
    pub write_index: bool,
    /// BGZF compression level, 0-12.
    #[arg(long, value_name = "N", default_value = "6", value_parser = parse_level)]
    pub compression_level: u8,
}

/// Thread count.
#[derive(Args, Debug, Clone)]
pub struct ThreadArgs {
    /// Number of threads (default: all available cores).
    #[arg(short = 't', long, value_name = "N", env = "GLOCI_THREADS", value_parser = parse_threads)]
    pub threads: Option<NonZeroUsize>,
}

#[derive(Args, Debug, Clone)]
pub struct GenomeArg {
    /// Chromosome order and sizes: a built-in genome such as hg38 (see
    /// gloci genomes), or a .fai, .genome, .dict, VCF/BCF or BAM header file.
    #[arg(short, long, value_name = "NAME|FILE")]
    pub genome: Option<String>,
}

#[derive(Args, Debug, Clone)]
pub struct ChromOrderArg {
    /// Chromosome order when neither -g nor a VCF/BCF header sets one (spec
    /// section 4.5): karyotype (chr1..chr22, chrX, chrY, chrM) or lex (byte
    /// order, as sort -k1,1, sort-bed and bedtools sort).
    #[arg(
        long,
        value_name = "ORDER",
        default_value = "karyotype",
        value_parser = clap::builder::PossibleValuesParser::new(["karyotype", "lex"])
            .try_map(|s| s.parse::<gloci_core::ChromOrder>())
    )]
    pub chrom_order: gloci_core::ChromOrder,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum VcfSpan {
    Allele,
    Ref,
}

#[derive(Args, Debug, Clone)]
pub struct VcfSpanArg {
    /// Which footprint a VCF record covers: the allele span or the REF span.
    #[arg(long, value_enum, default_value = "allele", value_name = "MODE")]
    pub vcf_span: VcfSpan,
}

#[derive(Args, Debug, Clone)]
pub struct InvalidArg {
    /// Skip invalid records with a warning instead of failing.
    #[arg(long)]
    pub skip_invalid: bool,
}

#[derive(Args, Debug, Clone)]
pub struct TmpDirArg {
    /// Directory for temporary spill files; it must exist and be writable
    /// (default: $GLOCI_TMPDIR if set and non-empty, else the system temp dir).
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct RegionArgs {
    /// Restrict to a region such as chr1:1,000-2,000 (repeatable).
    #[arg(long, value_name = "REGION")]
    pub region: Vec<String>,
    /// Read regions from FILE (BED: 0-based, half-open).
    #[arg(long, value_name = "FILE")]
    pub regions_file: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct HeaderArg {
    /// Keep header lines in the output.
    #[arg(long)]
    pub header: bool,
}

#[derive(Args, Debug, Clone)]
pub struct FormatArg {
    /// Input format, overriding detection (bed, bedgraph, gff, vcf, bcf).
    #[arg(long, value_name = "FMT")]
    pub format: Option<String>,
}

/// Parse a fraction such as `0.5` or `50%` (spec 4.3).
fn parse_fraction(s: &str) -> Result<gloci_core::Fraction, String> {
    gloci_core::Fraction::parse(s).map_err(|e| e.to_string())
}

/// The A input of a relation command.
#[derive(Args, Debug, Clone)]
pub struct AArg {
    /// A: the file whose records are reported (`-` for stdin).
    #[arg(short = 'a', long = "a", value_name = "FILE", required = true)]
    pub a: String,
    /// A's format, overriding detection (bed, bedgraph, gff, vcf, bcf).
    #[arg(long, value_name = "FMT")]
    pub a_format: Option<String>,
}

/// The B inputs of a relation command.
#[derive(Args, Debug, Clone)]
pub struct BArgs {
    /// B: one or more files to relate A to, in the order given (-b x y or -b x -b y).
    #[arg(
        short = 'b',
        long = "b",
        value_name = "FILE",
        required = true,
        num_args = 1..,
        action = clap::ArgAction::Append
    )]
    pub b: Vec<String>,
    /// Labels for the B files, comma-separated, one per -b (default: the file
    /// basenames, spec 5).
    #[arg(long, value_name = "A,B,...", value_delimiter = ',')]
    pub names: Option<Vec<String>>,
    /// B's format, overriding detection (bed, bedgraph, gff, vcf, bcf).
    #[arg(long, value_name = "FMT")]
    pub b_format: Option<String>,
}

impl BArgs {
    /// One label per B file: `--names`, or the basenames. Usage error (exit
    /// 2) when `--names` does not give exactly one label per `-b`.
    pub fn labels(&self) -> Result<Vec<String>, String> {
        match &self.names {
            Some(n) if n.len() != self.b.len() => Err(format!(
                "--names gives {} label(s) but there are {} -b file(s); give one label per -b",
                n.len(),
                self.b.len()
            )),
            Some(n) => Ok(n.clone()),
            None => Ok(self
                .b
                .iter()
                .map(|p| {
                    std::path::Path::new(p)
                        .file_name()
                        .map_or_else(|| p.clone(), |f| f.to_string_lossy().into_owned())
                })
                .collect()),
        }
    }
}

/// Overlap thresholds, strand rule and block mode (spec 4.3).
#[derive(Args, Debug, Clone)]
pub struct OverlapArgs {
    /// Minimum fraction of A covered, e.g. 0.5 or 50% (4.3: fractions).
    #[arg(short = 'f', long, value_name = "FRAC", value_parser = parse_fraction)]
    pub min_frac_a: Option<gloci_core::Fraction>,
    /// Minimum fraction of B covered, e.g. 0.5 or 50% (4.3: fractions).
    #[arg(short = 'F', long, value_name = "FRAC", value_parser = parse_fraction)]
    pub min_frac_b: Option<gloci_core::Fraction>,
    /// Require both fractions to be met (4.3: -r).
    #[arg(short = 'r', long)]
    pub reciprocal: bool,
    /// Require either fraction to be met (4.3: -e).
    #[arg(short = 'e', long)]
    pub either: bool,
    /// Minimum number of overlapping bases (4.3: --min-bases).
    #[arg(long, value_name = "N")]
    pub min_bases: Option<u64>,
    /// Only overlaps on the same strand; unknown strands never match (4.3: strand).
    #[arg(short = 's', long = "same-strand", conflicts_with = "opposite_strand")]
    pub same_strand: bool,
    /// Only overlaps on opposite strands; unknown strands never match (4.3: strand).
    #[arg(short = 'S', long = "opposite-strand")]
    pub opposite_strand: bool,
    /// Compare BED12 blocks (exons) instead of whole spans (4.3: BED12 blocks).
    #[arg(long)]
    pub split: bool,
    /// A point on the edge of an interval overlaps it (4.3: points).
    #[arg(long)]
    pub inclusive_points: bool,
}

impl OverlapArgs {
    pub fn criteria(&self) -> gloci_core::OverlapCriteria {
        use gloci_core::StrandRule;
        gloci_core::OverlapCriteria {
            min_frac_a: self.min_frac_a,
            min_frac_b: self.min_frac_b,
            reciprocal: self.reciprocal,
            either: self.either,
            min_bases: self.min_bases,
            strand: match (self.same_strand, self.opposite_strand) {
                (true, _) => StrandRule::Same,
                (_, true) => StrandRule::Opposite,
                _ => StrandRule::Any,
            },
            inclusive_points: self.inclusive_points,
        }
    }
}

/// B padding (spec 4.3: windows). Cannot be combined with thresholds (D37).
#[derive(Args, Debug, Clone)]
#[command(group = clap::ArgGroup::new("pad").multiple(true).args(["window", "left", "right"]))]
pub struct WindowArgs {
    /// Pad B by N on both sides before testing (4.3: windows).
    #[arg(
        short = 'w',
        long = "window",
        value_name = "N",
        conflicts_with_all = WINDOW_CONFLICTS
    )]
    pub window: Option<u64>,
    /// Pad B by N on the left (4.3: windows).
    #[arg(long, value_name = "N", conflicts_with_all = WINDOW_CONFLICTS)]
    pub left: Option<u64>,
    /// Pad B by N on the right (4.3: windows).
    #[arg(long, value_name = "N", conflicts_with_all = WINDOW_CONFLICTS)]
    pub right: Option<u64>,
    /// Left and right follow A's strand (4.3: windows).
    #[arg(long, requires = "pad", conflicts_with_all = WINDOW_CONFLICTS)]
    pub strand_window: bool,
}

/// The threshold flags a window cannot be combined with (D37).
const WINDOW_CONFLICTS: [&str; 5] = [
    "min_frac_a",
    "min_frac_b",
    "reciprocal",
    "either",
    "min_bases",
];

impl WindowArgs {
    /// The window, when any window flag is given. `--left`/`--right`
    /// override `-w` for their side.
    pub fn window(&self) -> Option<gloci_core::Window> {
        if self.window.is_none() && self.left.is_none() && self.right.is_none() {
            return None;
        }
        let w = self.window.unwrap_or(0);
        Some(gloci_core::Window {
            left: self.left.unwrap_or(w),
            right: self.right.unwrap_or(w),
            strand_aware: self.strand_window,
        })
    }
}

/// Sort-free mode for B (spec 8.2).
#[derive(Args, Debug, Clone)]
pub struct UnsortedArg {
    /// Accept inputs in any order; B is held in memory (8.2: unsorted mode).
    #[arg(long)]
    pub unsorted: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_core::ChromOrder;

    #[derive(Parser)]
    struct T {
        #[command(flatten)]
        o: ChromOrderArg,
    }

    #[test]
    fn chrom_order_arg() {
        let order = |a: &[&str]| T::try_parse_from(a).map(|t| t.o.chrom_order);
        assert_eq!(order(&["t"]).unwrap(), ChromOrder::Karyotype);
        assert_eq!(
            order(&["t", "--chrom-order", "lex"]).unwrap(),
            ChromOrder::Lex
        );
        assert!(order(&["t", "--chrom-order", "natural"]).is_err());
    }

    #[derive(Parser)]
    struct R {
        #[command(flatten)]
        a: AArg,
        #[command(flatten)]
        b: BArgs,
        #[command(flatten)]
        o: OverlapArgs,
        #[command(flatten)]
        w: WindowArgs,
        #[command(flatten)]
        u: UnsortedArg,
    }

    fn parse(args: &[&str]) -> Result<R, clap::Error> {
        let mut v = vec!["t", "-a", "a.bed", "-b", "b.bed"];
        v.extend_from_slice(args);
        R::try_parse_from(v)
    }

    #[test]
    fn relation_groups_parse() {
        let r = R::try_parse_from([
            "t",
            "-a",
            "a.bed",
            "--a-format",
            "bed",
            "-b",
            "x.bed",
            "y.bed",
            "-b",
            "z.bed",
            "--names",
            "x,y,z",
            "--b-format",
            "bed",
            "-f",
            "50%",
            "-F",
            "0.25",
            "-r",
            "-e",
            "--min-bases",
            "3",
            "-s",
            "--split",
            "--inclusive-points",
            "--unsorted",
        ])
        .unwrap();
        assert_eq!(r.a.a, "a.bed");
        assert_eq!(r.b.b, ["x.bed", "y.bed", "z.bed"]);
        assert_eq!(r.b.labels().unwrap(), ["x", "y", "z"]);
        let c = r.o.criteria();
        assert_eq!(
            c.min_frac_a,
            Some(gloci_core::Fraction::parse("0.5").unwrap())
        );
        assert_eq!(c.min_bases, Some(3));
        assert_eq!(c.strand, gloci_core::StrandRule::Same);
        assert!(r.o.split && r.u.unsorted && r.w.window().is_none());
        assert!(parse(&["-f", "1.5"]).is_err());
        assert!(
            R::try_parse_from(["t", "-a", "a.bed"]).is_err(),
            "-b is required"
        );
    }

    #[test]
    fn names_must_match_b_count() {
        assert!(parse(&["--names", "x,y"]).unwrap().b.labels().is_err());
        let r = parse(&[]).unwrap();
        assert_eq!(r.b.labels().unwrap(), ["b.bed"]);
    }

    #[test]
    fn window_conflicts_with_thresholds() {
        use clap::error::ErrorKind::ArgumentConflict;
        for w in [["-w", "5"], ["--left", "5"], ["--right", "5"]] {
            for t in [
                &["-f", "0.5"][..],
                &["-F", "0.5"],
                &["-r"],
                &["-e"],
                &["--min-bases", "1"],
            ] {
                let args: Vec<&str> = w.iter().chain(t).copied().collect();
                let e = parse(&args)
                    .err()
                    .unwrap_or_else(|| panic!("{args:?} parsed"));
                assert_eq!(e.kind(), ArgumentConflict, "{args:?}");
                assert_eq!(e.exit_code(), 2);
            }
        }
        assert!(parse(&["--strand-window", "-r"]).is_err());
        assert!(parse(&["--strand-window"]).is_err(), "needs a padding flag");
        let r = parse(&["-w", "10", "--left", "2", "--strand-window", "-s"]).unwrap();
        let w = r.w.window().unwrap();
        assert_eq!((w.left, w.right, w.strand_aware), (2, 10, true));
    }

    #[test]
    fn same_and_opposite_strand_conflict() {
        assert!(parse(&["-s", "-S"]).is_err());
        assert!(parse(&["-S"]).is_ok());
    }
}
