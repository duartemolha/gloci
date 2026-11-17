use super::shared::{
    AlignArgs, ChromMapArg, ChromOrderArg, FormatArg, GenomeArg, HeaderArg, InvalidArg, OutputArgs,
    ThreadArgs, VcfSpanArg, parse_size,
};
use clap::Args;
use std::path::PathBuf;

/// Aggregate columns over runs of lines with equal keys (bedtools groupby).
///
/// Reads any tab-separated text, BED included, and groups consecutive lines whose -g key columns are
/// equal; the input does not have to be sorted, and a key that comes back later starts a new group
/// (as in bedtools groupby; sort first to get one group per key). Each group writes its key columns,
/// exactly as in its first line, followed by one column per -c COL:OP, with the ops of map (count,
/// count_distinct, sum, mean, median, min, max, stdev, var, stdev_pop, var_pop, collapse, distinct,
/// distinct_sorted, first, last, indicator, cv, mad, quantile:Q, tmean:LO-HI, mode, antimode); the
/// ops that need a second input (wmean, records, min_record, max_record, frac_*, bases_uniq, bases)
/// are for map.
#[derive(Args)]
#[command(about = "Aggregate columns over runs of equal keys (bedtools groupby)")]
#[command(after_long_help = "\
Examples:
  gloci groupby scores.bed -g 1,2 -c 5:sum -c 5:mean
  gloci groupby calls.tsv -g 1 -c 4:collapse -c 5:count_distinct
  gloci groupby table.tsv --header -g gene -c expr:median -o medians.tsv

Input lines:
Blank lines, `#` lines and track/browser lines are skipped. A line with fewer columns than \
-g needs is an error naming the line.

Header and other input:
With --header the first line names the columns: -g and -c may use the names, and the \
output starts with a line holding the key names and the -c labels (bedtools -header copies \
the input header unchanged). VCF and BCF input is refused; extract the fields first with \
gloci vcf extract. --threads is accepted for uniformity with the other commands; groupby \
runs on a single thread.")]
pub struct GroupbyArgs {
    /// Tab-separated input (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Key columns: 1-based numbers (or names with --header), comma-separated.
    #[arg(short = 'g', long = "group", value_name = "COLS", required = true)]
    pub group: String,
    /// Aggregate a column per group, such as 5:sum or 4:collapse (repeatable;
    /// commas also separate).
    #[arg(short = 'c', long = "column", value_name = "COL:OP", required = true)]
    pub columns: Vec<String>,
    /// The first line names the columns (allowing names in -g and -c) and
    /// the output starts with a line naming its columns.
    #[arg(long)]
    pub header: bool,
    /// Value for an aggregation with no values.
    #[arg(long, value_name = "S", default_value = ".")]
    pub empty: String,
    /// Separator for collapse and distinct lists.
    #[arg(long, value_name = "S", default_value = ",")]
    pub delim: String,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Merge overlapping or nearby records into clusters (bedtools merge, bedops -m).
///
/// Records are grouped by their footprints: a record joins the open cluster when its start is at most
/// the cluster end plus -d/--max-gap. The default 0 joins overlapping and bookended features; a
/// negative gap such as -d -2 needs at least that many shared bases; a point joins a cluster it lies
/// inside or touches. The output is the merged regions as BED (chrom, start, end) whatever the input
/// format (VCF, GFF, BAM, ...), or TSV for -o x.tsv; the input records themselves are not written.
#[derive(Args)]
#[command(about = "Merge overlapping or near records (bedtools merge, bedops -m)")]
#[command(after_long_help = "\
Examples:
  gloci merge peaks.bed
  gloci merge a.bed b.bed -d 100 --count
  gloci merge genes.bed -s -c 5:max -c 4:distinct
  # variants in, BED regions out (-d -1: only footprints sharing a base join)
  gloci merge calls.vcf.gz -d -1 -o regions.bed.gz --write-index

Counting and columns:
--count adds the number of records per cluster and -c COL:OP adds one aggregation column \
per option (count, sum, mean, median, min, max, stdev, var, stdev_pop, var_pop, collapse, \
distinct, distinct_sorted, first, last, count_distinct, indicator, cv, mad, quantile:Q, \
tmean:LO-HI, mode, antimode; the ops that need a second input, like wmean, are for map). \
-s clusters each strand separately (unknown strand apart from + and -) without adding a \
strand column.

Sorted input:
Several files are merged together as one stream. Every input must be sorted in one \
chromosome order, which is checked as it is read; an unsorted input is an error whose hint \
is to run gloci sort or to pass --sort, which sorts the inputs first with an external sort \
bounded by --max-mem (spilling to --tmp-dir); the result equals gloci sort FILES | gloci \
merge -. Memory stays flat apart from the clusters being built.")]
pub struct MergeArgs {
    /// Input files (`-` for stdin).
    #[arg(value_name = "FILES", required = true)]
    pub files: Vec<String>,
    /// Join records up to N bases apart; negative values require that much
    /// overlap.
    #[arg(
        short = 'd',
        long = "max-gap",
        value_name = "N",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    pub max_gap: i64,
    /// Cluster each strand separately.
    #[arg(short = 's', long)]
    pub strand: bool,
    /// Add a column with the number of records in each cluster.
    #[arg(long)]
    pub count: bool,
    /// Aggregate a column per cluster, such as 5:max or 4:distinct
    /// (repeatable; commas also separate).
    #[arg(short = 'c', long = "column", value_name = "COL:OP")]
    pub columns: Vec<String>,
    /// Value for an aggregation with no values.
    #[arg(long, value_name = "S", default_value = ".")]
    pub empty: String,
    /// Separator for collapse and distinct lists.
    #[arg(long, value_name = "S", default_value = ",")]
    pub delim: String,
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
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub sort: SortArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Give every record the id of its cluster (bedtools cluster).
///
/// Writes each input line unchanged followed by a tab and a cluster id. Clusters are those of gloci
/// merge, with the same -d/--max-gap (default 0: overlapping and bookended records share a cluster; a
/// negative gap needs that many shared bases) and -s (each strand clustered separately, unknown
/// strand apart from + and -). Ids start at 1 and increase through the file, across chromosomes. With
/// -s, a cluster takes the next id when its first record is read, and the records stay in input order
/// (bedtools cluster -s regroups them by strand within each chromosome).
#[derive(Args)]
#[command(after_long_help = "Examples:
  gloci cluster peaks.bed
  gloci cluster peaks.bed -d 100 -s
  gloci cluster calls.vcf.gz --vcf-span ref -d 50

Output:
The output is TSV, not a copy of the input format: a VCF input (clustered by its \
footprints, see --vcf-span) gets a #CHROM column line ending in \"cluster\" and then its \
records with the id appended. --header keeps the header lines of a BED or GFF input.

Sorted input:
The input must be sorted in one chromosome order, checked as it is read; pass --sort to \
sort it first.")]
pub struct ClusterArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Join records up to N bases apart; negative values require that much
    /// overlap.
    #[arg(
        short = 'd',
        long = "max-gap",
        value_name = "N",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    pub max_gap: i64,
    /// Cluster each strand separately.
    #[arg(short = 's', long)]
    pub strand: bool,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub sort: SortArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Give every record its distance to the previous record (bedtools spacing).
///
/// Writes each input line unchanged followed by a tab and the number of
/// bases between the record and the one before it in the file on the
/// same chromosome: . for the first record of a chromosome, 0 when
/// they are bookended, -1 when their footprints overlap (including a
/// record inside the previous one; only the previous record counts).
/// Unlike bedtools, a zero-length record is a point on its base
/// rather than a base wide on each side. The output is TSV, not a copy
/// of the input format; for VCF input see gloci cluster. The input
/// must be sorted in one chromosome order; pass --sort to sort it first.
#[derive(Args)]
#[command(about = "Distance to the previous record (bedtools spacing)")]
#[command(after_long_help = "Examples:
  gloci spacing peaks.bed
  gloci spacing calls.vcf.gz --vcf-span ref")]
pub struct SpacingArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub sort: SortArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Split the inputs at every boundary into segments and say who covers each (bedtools multiinter and
/// unionbedg, bedops -p, -i and -s).
///
/// The footprints of all files are cut at every start and end, giving disjoint segments in sorted
/// order. Each row is chrom, start, end, the number of files covering the segment and their labels
/// joined with commas (--names, default the file names without directories). A file counts once
/// however many of its features overlap there, and points cover no bases, so they add neither
/// segments nor boundaries. Only segments covered by at least one file are written.
#[derive(Args)]
#[command(about = "Segment inputs by coverage (bedtools multiinter, bedops -p)")]
#[command(after_long_help = "\
Examples:
  gloci partition a.bed b.bed c.bed --merge-adjacent
  gloci partition a.bed b.bed --names tumour,normal --min-files 2
  gloci partition a.bed b.bed --all --merge-adjacent --no-membership
  gloci partition a.bg b.bg --bedgraph --empty 0 -g hg38

Selecting segments:
--min-files K and --max-files K keep segments covered by at least or at most K files; \
--all keeps those covered by every file. --empty implies --min-files 0 (an explicit \
--min-files overrides it); --min-files 0 needs --empty. By default a boundary is kept even \
where the membership does not change (as in bedops -p); --merge-adjacent joins \
neighbouring segments with the same output columns (as multiinter does). --no-membership \
writes coordinates only, so with --merge-adjacent every run of touching kept segments is \
joined.

bedGraph input:
With --bedgraph the inputs are bedGraph files and each row has one value per file instead \
of the membership; a file that does not cover the segment shows --empty (default 0). \
Giving --empty also writes the segments no file covers, with that value in every column; \
with -g (needs lengths) it covers whole chromosomes, including those with no features, \
otherwise only the gaps between features. The intervals of one bedGraph must not overlap.

Sorted input:
Every input must be sorted in one chromosome order, which is checked as it is read; an \
unsorted input is an error whose hint names --sort, which sorts each file first with an \
external sort bounded by --max-mem. Memory is bounded by the intervals open at one \
position.")]
pub struct PartitionArgs {
    /// Input files (`-` for stdin).
    #[arg(value_name = "FILES", required = true)]
    pub files: Vec<String>,
    /// Labels of the files, comma-separated, one per file (default: the
    /// file names without directories).
    #[arg(long, value_name = "A,B,...", value_delimiter = ',')]
    pub names: Vec<String>,
    /// Keep segments covered by at least K files (default 1).
    #[arg(long, value_name = "K", conflicts_with = "all")]
    pub min_files: Option<usize>,
    /// Keep segments covered by at most K files.
    #[arg(long, value_name = "K", conflicts_with = "all")]
    pub max_files: Option<usize>,
    /// Keep segments covered by every file.
    #[arg(long)]
    pub all: bool,
    /// Join neighbouring segments whose output columns are equal.
    #[arg(long)]
    pub merge_adjacent: bool,
    /// Write coordinates only (BED3).
    #[arg(long, conflicts_with = "bedgraph")]
    pub no_membership: bool,
    /// Inputs are bedGraphs: write one value column per file.
    #[arg(long)]
    pub bedgraph: bool,
    /// With --bedgraph: the value of a file that does not cover a segment
    /// (default 0); also writes the segments no file covers.
    #[arg(long, value_name = "VALUE", requires = "bedgraph")]
    pub empty: Option<String>,
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
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub sort: SortArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Depth of coverage of a genome by the features of a BED file or the reads of a BAM (bedtools
/// genomecov -i and -ibam).
///
/// Every feature adds one to the depth over its bases. -g is required (unless the input is BAM, CRAM
/// or SAM): the output covers every chromosome of the genome from base 0 to its length, in genome
/// order, chromosomes without features included (a built-in genome, .genome, .fai or .dict file with
/// lengths). The default --hist writes, per chromosome and then for the whole genome (chrom
/// `genome`), one row per depth present: chrom, depth, bases at that depth, chromosome length,
/// fraction (six significant digits, as bedtools).
#[derive(Args)]
#[command(about = "Coverage depth of a genome (bedtools genomecov)")]
#[command(after_long_help = "\
Examples:
  gloci genomecov peaks.bed -g hg38
  gloci genomecov reads.bed -g hg38 --bedgraph --scale 0.000001 -o cov.bedgraph
  gloci genomecov exons12.bed -g hg38 --split --per-base
  gloci genomecov reads.bam --bga --split

Output:
--bedgraph writes runs of equal depth above zero, --bga also the zero-depth runs, \
--per-base one row per base: chrom, 1-based position, depth. --scale X multiplies the \
depths of the last three (six significant digits, as bedtools); --max-depth N puts every \
depth of N or more into the bin N of --hist.

Selecting features:
--strand + or - counts only features on that strand (BED needs column 6; a feature with \
unknown strand is not counted). --split counts the BED12 blocks of each feature. --5p and \
--3p count one base per feature, its 5' or 3' end by strand (unknown counts as plus). A \
zero-length feature covers no bases and is ignored. A feature past the end of its \
chromosome is an error.

Alignments:
BAM, CRAM and SAM input (a CRAM needs --reference): the chromosomes and sizes come from \
the header (header order, empty chromosomes included) unless -g is given. A read counts \
its whole reference span, deletions included; --split counts the blocks between N \
operations, and --split-deletions also cuts at D. --split-deletions is bedtools -ibam \
-split, so use it for parity; plain -ibam (cuts at D, not N) has no gloci counterpart, see \
docs/differences.md#bam-genomecov-deletions. Unmapped reads are left out, and --min-mapq, \
--require-flags and --exclude-flags filter reads first. --strand follows the read's \
direction. An indexed BAM or CRAM at -t above 1 runs per chromosome, with the same output \
as -t 1. At most 8 runs are in flight; memory is then about 64 MiB plus 20 MiB per run in \
flight: the output waiting to be written shares 64 MiB and the rest spills to --tmp-dir.

Sorted input:
The input must be sorted in the genome order, which is checked as it is read; an unsorted \
input is an error whose hint names --sort. Memory is bounded by the features open at one \
position, whatever the output mode.")]
pub struct GenomecovArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Chromosome order and sizes (required unless the input is BAM, CRAM or SAM,
    /// which gives its header's): a built-in genome such as hg38 (see gloci
    /// genomes list), or a .fai, .genome or .dict file.
    #[arg(short, long, value_name = "NAME|FILE")]
    pub genome: Option<String>,
    /// Write the histogram of depths (the default).
    #[arg(long, conflicts_with_all = ["bedgraph", "bga", "per_base"])]
    pub hist: bool,
    /// Write runs of equal depth above zero as bedGraph.
    #[arg(long, conflicts_with_all = ["bga", "per_base"])]
    pub bedgraph: bool,
    /// Write runs of equal depth as bedGraph, zero-depth runs included.
    #[arg(long, conflicts_with = "per_base")]
    pub bga: bool,
    /// Write the depth of every base (1-based positions).
    #[arg(long)]
    pub per_base: bool,
    /// Multiply the reported depths by X (not with --hist).
    #[arg(long, value_name = "X", allow_negative_numbers = true)]
    pub scale: Option<f64>,
    /// Count only features on this strand (+ or -).
    #[arg(long, value_name = "+|-", allow_hyphen_values = true,
          value_parser = clap::builder::PossibleValuesParser::new(["+", "-"]))]
    pub strand: Option<String>,
    /// Count the BED12 blocks of each feature, or the blocks of each read
    /// between N operations (BAM, CRAM and SAM), instead of the whole feature.
    #[arg(long, conflicts_with_all = ["five_prime", "three_prime"])]
    pub split: bool,
    /// BAM, CRAM and SAM: also cut the blocks of a read at deletions (D), with or
    /// without --split; this is bedtools genomecov -ibam -split. Plain -ibam
    /// (cuts at D, not N) has no gloci counterpart, see
    /// docs/differences.md#bam-genomecov-deletions.
    #[arg(long, conflicts_with_all = ["five_prime", "three_prime"])]
    pub split_deletions: bool,
    /// Count only the 5' base of each feature (strand-aware).
    #[arg(long = "5p", conflicts_with = "three_prime")]
    pub five_prime: bool,
    /// Count only the 3' base of each feature (strand-aware).
    #[arg(long = "3p")]
    pub three_prime: bool,
    /// With --hist: put every depth of N or more into the bin N.
    #[arg(long, value_name = "N")]
    pub max_depth: Option<u64>,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub sort: SortArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Write the bases no feature covers (bedtools complement, bedops -c).
///
/// Features of all files are united, and the gaps between them are written as BED (chrom, start,
/// end). With -g, the gaps are taken over each chromosome of the genome, in genome order: a
/// chromosome with no features becomes one row over its whole length, and the gaps before the first
/// and after the last feature are included. The genome needs lengths (a built-in genome, .genome,
/// .fai or .dict file); a feature that reaches past the end of its chromosome is an error. Without -g
/// only the gaps between features of the same chromosome are written.
#[derive(Args)]
#[command(about = "Bases no feature covers (bedtools complement, bedops -c)")]
#[command(after_long_help = "\
Examples:
  gloci complement peaks.bed
  gloci complement peaks.bed -g hg38 -o uncovered.bed
  gloci complement a.bed b.bed -g sacCer3
  gloci complement shuffled.bed -g genome.txt --sort --max-mem 512M

Points and bookends:
Points cover no bases and are ignored. Bookended features leave no gap.

Sorted input:
Every input must be sorted in one chromosome order, which is checked as it is read; with \
-g, use the genome's own order (gloci sort -g). An unsorted input is an error whose hint \
names --sort, which sorts the inputs first with an external sort bounded by --max-mem; the \
result equals gloci sort FILES | gloci complement -.")]
pub struct ComplementArgs {
    /// Input files (`-` for stdin).
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
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub sort: SortArg,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// The amounts shared by `slop` and `flank`.
#[derive(Args, Debug, Clone)]
pub struct PadArgs {
    /// Pad both sides by N.
    #[arg(short = 'b', long, value_name = "N", conflicts_with_all = ["left", "right"],
          required_unless_present_any = ["left", "right"])]
    pub both: Option<String>,
    /// Pad the left side (smaller coordinates) by N.
    #[arg(long, value_name = "N")]
    pub left: Option<String>,
    /// Pad the right side by N.
    #[arg(long, value_name = "N")]
    pub right: Option<String>,
    /// Swap --left and --right for minus-strand records.
    #[arg(short = 's')]
    pub strand: bool,
    /// Read each amount as a fraction of the feature length (0.5 or 50%,
    /// floored) instead of bases.
    #[arg(long)]
    pub pct: bool,
}

/// Grow each feature on its sides (bedtools slop, bedops --range).
///
/// Every record is written with its coordinates changed to
/// [start - L, end + R), clipped at 0 and, with -g, at the end of its
/// chromosome (a built-in genome, .genome, .fai or .dict file with
/// lengths; a chromosome missing from it is an error). Amounts are bases,
/// or with --pct a decimal fraction of each feature's length (0.5 or
/// 50%, floored). -b pads both sides; --left and --right pad one each.
/// -s swaps them for minus-strand records (an unknown strand counts as
/// plus). BED and GFF/GTF only; VCF/BCF input is a usage error. The
/// input is streamed in its own order and needs no sorting; the output
/// keeps the input format, and a warning is printed once when -s or
/// --pct makes it go out of order.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci slop genes.bed -b 1000 -g hg38
  gloci slop genes.bed --left 2000 --right 500 -s -g hg38 -o promoters.bed
  gloci slop peaks.bed --pct -b 0.5 -g genome.txt
  gloci slop annotation.gff3 -b 100 -g hg38")]
pub struct SlopArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    #[command(flatten)]
    pub pad: PadArgs,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Move each feature along its chromosome (bedtools shift).
///
/// --by N moves every record by N bases (negative moves toward 0);
/// --pct reads N as a fraction of the feature's own length (0.5, -25%;
/// the product is floored using exact decimal arithmetic). -s negates the shift on
/// minus-strand records (an unknown strand counts as +), so a positive N
/// moves every feature downstream. The result is clipped like slop: at 0
/// and, with -g, at the chromosome end, so a feature pushed over an edge
/// shrinks; one pushed wholly over an edge becomes the 1-base interval at
/// that edge, as in bedtools (a zero-length feature stays zero-length
/// there). BED and GFF/GTF only. The output
/// is in input order and may not be sorted; a warning is printed once
/// when the shift itself breaks the order.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci shift genes.bed --by 100 -g hg38
  gloci shift genes.bed --by -50 -g hg38
  gloci shift peaks.bed --by 1000 -s -g hg38 | gloci sort -g hg38 -
  gloci shift peaks.bed --pct --by 0.5 -g genome.txt")]
pub struct ShiftArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Bases to move each feature by; negative moves toward the start of the
    /// chromosome (with --pct, a signed fraction of the feature length).
    #[arg(long, value_name = "N", required = true, allow_hyphen_values = true)]
    pub by: String,
    /// Negate the shift for minus-strand records (unknown strand counts as +).
    #[arg(short = 's')]
    pub strand: bool,
    /// Read --by as a fraction of the feature length (0.5, -25%, floored)
    /// instead of bases.
    #[arg(long)]
    pub pct: bool,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Write the strips beside each feature (bedtools flank).
///
/// For every record, a left strip [start - L, start) and then a right
/// strip [end, end + R) are written as records with the input's other
/// columns, clipped at 0 and, with -g, at the chromosome end. A strip
/// clipped to nothing is not written. Options as for slop: -b, --left,
/// --right, -s, --pct, -g. The output is in input order and may not be
/// sorted; a warning is printed once when it goes out of order.
#[derive(Args)]
#[command(after_long_help = "\
Examples:
  gloci flank genes.bed -b 100 -g hg38
  gloci flank genes.bed --left 2000 -s -g hg38
  gloci flank peaks.bed --pct -b 0.25 -g genome.txt | gloci sort -g genome.txt -")]
pub struct FlankArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    #[command(flatten)]
    pub pad: PadArgs,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Cut chromosomes or intervals into windows (bedtools makewindows,
/// bedops --chop).
///
/// -g cuts every chromosome of a genome (built-in, .genome, .fai or .dict
/// with lengths) in genome order; -i cuts the footprint of every record
/// of a file in input order (no sorting needed; points give no windows;
/// -i - reads stdin). -w SIZE gives windows of SIZE bases every --step
/// bases (default SIZE; a smaller step overlaps, a larger one leaves
/// gaps), the last one truncated at the interval end. -n COUNT gives
/// COUNT windows of floor(length / COUNT) bases, the last one taking
/// the remainder; an interval shorter than COUNT is skipped with a
/// warning. Output is BED3; --id-column adds the 1-based window number
/// of each source interval, counted from the end with --reverse.
#[derive(Args)]
#[command(about = "Window chromosomes or intervals (bedtools makewindows)")]
#[command(after_long_help = "\
Examples:
  gloci windows -g hg38 -w 1000000 -o bins.bed
  gloci windows -g genome.txt -w 1000 --step 500
  gloci windows -i genes.bed -n 10 --id-column --reverse
  gloci merge peaks.bed | gloci windows -i - -w 100")]
#[command(group(clap::ArgGroup::new("source").required(true).multiple(false).args(["genome", "intervals"])))]
#[command(group(clap::ArgGroup::new("split").required(true).multiple(false).args(["size", "count"])))]
pub struct WindowsArgs {
    /// Cut every chromosome of a built-in genome such as hg38 (see gloci
    /// genomes list), or a .fai, .genome or .dict file with lengths.
    #[arg(short, long, value_name = "NAME|FILE")]
    pub genome: Option<String>,
    /// Cut the footprint of every record of FILE (`-` for stdin).
    #[arg(short = 'i', long, value_name = "FILE")]
    pub intervals: Option<String>,
    /// Window size in bases.
    #[arg(short = 'w', long, value_name = "SIZE")]
    pub size: Option<u64>,
    /// Distance between window starts (default: the window size).
    #[arg(long, value_name = "S", requires = "size", conflicts_with = "count")]
    pub step: Option<u64>,
    /// Split every interval into COUNT windows.
    #[arg(short = 'n', long, value_name = "COUNT")]
    pub count: Option<u64>,
    /// Number the windows from the end of the interval.
    #[arg(long)]
    pub reverse: bool,
    /// Append the 1-based window number of each source interval.
    #[arg(long)]
    pub id_column: bool,
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

/// `--sort` and the options that bound it (D39).
#[derive(Args, Debug, Clone)]
pub struct SortArg {
    /// Sort the inputs first, for inputs that are unsorted or too large to
    /// sort in memory (external sort bounded by --max-mem). The sorted
    /// records are written to a temporary file in --tmp-dir and read back,
    /// so besides the spill runs (about the size of the input, compressed)
    /// the directory needs one uncompressed sorted copy: allow roughly twice
    /// the input size (BAM and CRAM inputs are written as BAM; a SAM input as plain SAM,
    /// several times larger). Everything is removed afterwards.
    #[arg(long)]
    pub sort: bool,
    /// Memory budget for --sort, in binary units (K, M, G, T); at least 16M.
    #[arg(long, value_name = "SIZE", default_value = "2G", value_parser = parse_size, requires = "sort")]
    pub max_mem: u64,
    /// Directory for --sort's temporary files, ignored without --sort; it
    /// must exist and be writable (default: $GLOCI_TMPDIR if set and
    /// non-empty, else the system temp dir).
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<PathBuf>,
}

/// Write an IGV batch script that takes a snapshot of each feature (bedtools igv).
///
/// Each BED, GFF or VCF record becomes a goto line for IGV's batch mode and a snapshot line
/// that names the image chrom_start_end.png, with the record's name and the slop added to the
/// name on request. The script goes to stdout or -o; run it from IGV (Tools, Run Batch
/// Script) after loading the genome and tracks, or name a session with --session. Records are
/// read in order and need no sorting. As in bedtools, the goto line shows the 0-based start,
/// a zero-length BED record is shown one base wider on each side, and every BED or GFF record
/// must have the same number of columns.
#[derive(Args)]
#[command(about = "IGV snapshot script for each feature (bedtools igv)")]
#[command(after_long_help = "\
Examples:
  gloci igv peaks.bed -o snapshots.txt
  gloci igv peaks.bed --path /data/shots --session run.xml --sort position --collapse
  gloci igv peaks.bed --name --slop 500 --image svg

Names:
With --name the image is named chrom_start_end_NAME. NAME is the BED name column, the GFF \
feature type, or for a VCF record REF/ALT followed by _ID when the ID is not a dot. A record \
with no name is an error.

Differences from bedtools:
--image takes png, eps or svg (bedtools writes any extension it is given), and --slop must \
be a whole number from 0 to 2147483647 (bedtools reads a bad value as 0 and a negative one \
as an error).")]
pub struct IgvArgs {
    /// Input BED, GFF or VCF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The directory IGV writes the images to, as IGV is to read it.
    #[arg(long, value_name = "DIR", default_value = "./")]
    pub path: String,
    /// The path of an IGV session file to load before the snapshots, as IGV is to read it.
    #[arg(long, value_name = "PATH")]
    pub session: Option<String>,
    /// Sort the reads of each image: base, position, strand, quality, sample or readGroup.
    #[arg(
        long,
        value_name = "KEY",
        value_parser = clap::builder::PossibleValuesParser::new(gloci_ops::IGV_SORTS)
    )]
    pub sort: Option<String>,
    /// Collapse the reads before each snapshot.
    #[arg(long)]
    pub collapse: bool,
    /// Put the name of each record in its image's file name.
    #[arg(long)]
    pub name: bool,
    /// Bases to show on each side of the feature (the image file name ends _slopN).
    #[arg(long, value_name = "N", default_value_t = 0, value_parser = parse_slop)]
    pub slop: i64,
    /// The image type.
    #[arg(
        long,
        value_name = "TYPE",
        default_value = "png",
        value_parser = ["png", "eps", "svg"]
    )]
    pub image: String,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the script to FILE instead of stdout; .gz compresses it. A file is written only
    /// when the command succeeds.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}

fn parse_slop(s: &str) -> Result<i64, String> {
    match s.parse::<i64>() {
        Ok(n) if (0..=i64::from(i32::MAX)).contains(&n) => Ok(n),
        _ => Err(format!(
            "\"{s}\" is not a whole number of bases from 0 to {}",
            i32::MAX
        )),
    }
}

/// Write an HTML page with a link to the UCSC Genome Browser for each feature (bedtools links).
///
/// The page has a table with one row per BED, GFF or VCF record: a link to the record's
/// region, and the record's name, score and strand when it has them. Open it in a browser. The
/// links point to the UCSC browser for human hg18, as bedtools does; --base, --org and --db
/// point them at another build or a mirror. Records are read in order and need no sorting. As
/// in bedtools, a zero-length BED record is shown one base wider on each side and every BED or
/// GFF record must have the same number of columns.
#[derive(Args)]
#[command(about = "HTML genome browser links per feature (bedtools links)")]
#[command(after_long_help = "\
Examples:
  gloci links peaks.bed -o peaks.html
  gloci links peaks.bed --org human --db hg38 -o peaks.html
  gloci links peaks.bed --base http://mirror.example.edu --org mouse --db mm9

Columns:
A BED record of 4, 5, 6, 9 or 12 columns gets the name, then the score, then the strand, as \
far as it has them. A BED record of any other length gets the link alone. A GFF record gets \
its feature type, score and strand. A VCF record gets the link alone, or for 9 or 12 columns \
REF/ALT_ID, an empty score and the strand +. Nothing is escaped, as in bedtools.

Title:
The page title is the input's name as given on the command line, or stdin for -.")]
pub struct LinksArgs {
    /// Input BED, GFF or VCF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The browser's address; the links add /cgi-bin/hgTracks.
    #[arg(long, value_name = "URL", default_value = "http://genome.ucsc.edu")]
    pub base: String,
    /// The organism.
    #[arg(long, value_name = "ORG", default_value = "human")]
    pub org: String,
    /// The genome build.
    #[arg(long, value_name = "DB", default_value = "hg18")]
    pub db: String,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// Write the page to FILE instead of stdout; .gz compresses it. A file is written only
    /// when the command succeeds.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}

/// Move intervals to another assembly with a chain (CrossMap)
///
/// Lifts BED, bedGraph, GFF3 and GTF records through a chain file by the rules of CrossMap (its `bed`,
/// `gff` and `region` modes), not those of UCSC `liftOver`, whose rules it does not follow. A
/// record that crosses a gap in the chain is split into one record per aligned piece. Records
/// that do not lift are written as read to --unmapped, or counted in a warning. The output is
/// in input order, so pipe it through `gloci sort` when the target needs sorting. CHAIN is a file
/// (plain, gzip or BGZF text) or a built-in id: GRCh37-to-GRCh38, GRCh38-to-GRCh37,
/// GRCh38-to-CHM13v2, CHM13v2-to-GRCh38, GRCh37-to-CHM13v2 or CHM13v2-to-GRCh37. Ensembl's and the
/// T2T consortium's chains are free to use; a UCSC chain file brings UCSC's licence
/// (non-commercial use only) with it.
#[derive(Args)]
#[command(about = "Move intervals to another assembly with a chain (CrossMap)")]
#[command(after_long_help = "\
Examples:
  gloci liftover peaks.bed --chain GRCh37-to-GRCh38 -o peaks.b38.bed --unmapped peaks.unmapped.bed
  gloci liftover genes.gtf --chain GRCh37_to_GRCh38.chain.gz | gloci sort - -o genes.b38.gtf
  gloci liftover cnv.bed --chain CHM13v2-to-GRCh38 --chrom-style long -o cnv.hg38.bed")]
pub struct LiftoverArgs {
    /// Input BED, bedGraph, GFF3 or GTF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The chain to lift through: a file (plain, gzip or BGZF; `-` for stdin) or a built-in id, such
    /// as GRCh37-to-GRCh38.
    #[arg(long, value_name = "CHAIN", required = true)]
    pub chain: String,
    /// Write the records that did not lift to FILE, as read, in the input's format whatever the
    /// file's name (CrossMap `--unmap-file`; `.gz` compresses it; `-` is stdout when -o names a
    /// file).
    #[arg(long, value_name = "FILE")]
    pub unmapped: Option<PathBuf>,
    /// How the target chromosomes are named (CrossMap `--chromid a|s|l|n`): as-input follows
    /// the input's `chr` prefix, short removes `chr`, long adds it, chain keeps the chain's name.
    #[arg(
        long,
        value_name = "STYLE",
        default_value = "as-input",
        value_parser = ["as-input", "short", "long", "chain"]
    )]
    pub chrom_style: String,
    /// Lift a BED12 record by columns 1 to 3 only, as any other BED record, leaving blocks and
    /// thick range as read (CrossMap `--naive-bed-parsing`). Without it, each block is lifted
    /// and must fall in one piece. Not for GFF or GTF.
    #[arg(long)]
    pub ignore_blocks: bool,
    /// Write a record that crosses chain gaps as one record from its first piece's start to its
    /// last piece's end, when enough of it maps (CrossMap `region`). Pieces on different
    /// chromosomes or strands leave it unmapped. Not for GFF or GTF, nor BED12 without
    /// --ignore-blocks.
    #[arg(long)]
    pub merge_segments: bool,
    /// The fraction of a record that must map for --merge-segments to merge it, from 0 to 1
    /// (CrossMap `region -r`).
    #[arg(
        long,
        value_name = "FRAC",
        default_value_t = 0.85,
        requires = "merge_segments",
        allow_negative_numbers = true
    )]
    pub min_ratio: f64,
    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}
