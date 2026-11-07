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
    #[command(after_long_help = "\
Examples:
  gloci index calls.vcf.gz
  gloci index peaks.bed.gz -t 8")]
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
    /// With no regions the records are written exactly as read, so view can
    /// recompress a file (-o FILE.gz is BGZF) or index it (--write-index).
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
    #[command(flatten)]
    pub threads: ThreadArgs,
}

/// Output destination and encoding.
#[derive(Args, Debug, Clone)]
pub struct OutputArgs {
    /// Write output to FILE (format from its extension) instead of stdout.
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
}
