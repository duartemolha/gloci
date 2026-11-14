//! clap definitions: global flags, shared option groups and the subcommands.
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

mod files;
mod intervals;
mod relate;
mod sampling;
mod shared;

pub use files::*;
pub use intervals::*;
pub use relate::*;
pub use sampling::*;
pub use shared::*;

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

/// VCF and BCF operations: filter, extract, stats, concat, split, norm, annotate, merge and compare.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf filter cohort.bcf --pass-only --min-qual 30 -o filtered.bcf
  gloci vcf filter input.vcf --samples s1,s2 --min-dp 10 --recalc-ac -o out.vcf
  gloci vcf filter variants.bcf --types indel --biallelic")]
pub struct VcfArgs {
    #[command(subcommand)]
    pub command: VcfCommand,
}

#[derive(Subcommand, Debug, Clone)]
pub enum VcfCommand {
    #[command(before_help = EXPERIMENTAL)]
    Filter(Box<VcfFilterArgs>),
    #[command(before_help = EXPERIMENTAL)]
    Extract(VcfExtractArgs),
    #[command(before_help = EXPERIMENTAL)]
    Stats(VcfStatsArgs),
    #[command(before_help = EXPERIMENTAL)]
    Concat(VcfConcatArgs),
    #[command(before_help = EXPERIMENTAL)]
    Split(VcfSplitArgs),
    // P3b Task 21
    #[command(before_help = EXPERIMENTAL)]
    Norm(VcfNormArgs),
    // P3b Task 22
    #[command(before_help = EXPERIMENTAL)]
    Annotate(VcfAnnotateArgs),
    // P3b Task 23
    #[command(before_help = EXPERIMENTAL)]
    Merge(VcfMergeArgs),
    #[command(before_help = EXPERIMENTAL)]
    Compare(Box<VcfCompareArgs>),
}

/// `vcf compare --match`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum CompareMatchArg {
    Alleles,
    Position,
}

/// Compare the genotypes of two VCF/BCF files: per site, per sample and
/// as a matrix (vcftools --diff-site-discordance, --diff-indv-discordance,
/// --diff-discordance-matrix).
///
/// Both files are read as sorted streams (the chromosome order of -g,
/// else A's ##contig lines, else --chrom-order; an out-of-order record is
/// exit 1 and names the file), so no index is needed. Sites are paired by
/// --match: alleles (the default) pairs records at one CHROM and POS that
/// share an ALT allele after trimming, position pairs as vcftools does (the
/// same POS with the same REF, where a REF of N or . takes the other's).
/// Discordance is counted over the sites in both files and the samples in
/// both files (the union by name, after --sample-map renames B's samples).
/// --thin needs both files to declare the same ##contig lines and no -g
/// (exit 2 otherwise). A and B are always read on threads of their own, plus the thread that
/// walks them, so -t 1 and -t 2 run three threads; more -t adds selection
/// workers up to -t in all.
/// A genotype is its unordered alleles, so 0/1 equals 1/0 and 0|1, and
/// 0/1 of A>C,G equals 0/2 of A>G,C; any missing allele makes it missing
/// (vcftools counts ./1 as called). --samples and --exclude-samples take
/// names present in either file and apply to each file that has them; a
/// name in neither is exit 2. With none of --site-discordance,
/// --sample-discordance and --discordance-matrix, all three files are
/// written, atomically and only when the whole run succeeds:
/// PREFIX.diff.sites (CHROM POS FILES MATCHING_ALLELES N_COMMON_CALLED
/// N_DISCORD DISCORDANCE, one row per site in both files, and per site in
/// one file with zeros and nan), PREFIX.diff.indv (every sample of the
/// union; one in a single file is 0, 0, nan) and
/// PREFIX.diff.discordance_matrix (biallelic pairs of the same site:
/// with --match alleles the same trimmed allele, so AT>A and ATG>AG
/// count; with --match position the same REF and ALT text; haploid
/// genotypes are not counted). Not available: vcftools
/// --diff-site and --diff-indv (use gloci intersect and subtract --match
/// alleles for sites, and the #CHROM lines for samples).
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf compare truth.vcf.gz calls.vcf.gz --out-prefix cmp
  gloci vcf compare a.bcf b.bcf --out-prefix cmp --sample-discordance --match position
  gloci vcf compare a.vcf.gz b.vcf.gz --out-prefix cmp --sample-map b-to-a.txt --samples s1,s2 --min-dp 10")]
pub struct VcfCompareArgs {
    /// The first VCF or BCF file (`-` for stdin), sorted.
    #[arg(value_name = "A")]
    pub a: String,
    /// The second VCF or BCF file, sorted in the same chromosome order (`-`
    /// is allowed for one of the two).
    #[arg(value_name = "B")]
    pub b: String,

    /// Applies to both files: samples, genotype masks, site filters,
    /// expressions and regions.
    #[command(flatten)]
    pub select: VcfSelectArgs,

    /// Prefix of the output files, PREFIX.diff.sites, PREFIX.diff.indv and
    /// PREFIX.diff.discordance_matrix. An existing file of the same name is
    /// replaced.
    #[arg(long, value_name = "PREFIX", required = true)]
    pub out_prefix: PathBuf,

    /// How sites of the two files are paired: alleles (the same CHROM and
    /// POS and a shared ALT allele; each record of A takes the first unmatched
    /// match in B) or position (vcftools: the same POS and REF).
    #[arg(
        long = "match",
        value_enum,
        default_value = "alleles",
        value_name = "RULE"
    )]
    pub match_rule: CompareMatchArg,

    /// Rename samples of B to their names in A: two columns, the name in B
    /// then the name in A, whitespace-separated, `#` comments allowed
    /// (vcftools --diff-indv-map). A sample of B the file does not name keeps
    /// its name.
    #[arg(long, visible_alias = "diff-indv-map", value_name = "FILE")]
    pub sample_map: Option<PathBuf>,

    /// Write PREFIX.diff.sites: per site, the common called genotypes and
    /// the discordant ones (vcftools --diff-site-discordance).
    #[arg(long, visible_alias = "diff-site-discordance")]
    pub site_discordance: bool,
    /// Write PREFIX.diff.indv: per sample, the same over the shared sites
    /// (vcftools --diff-indv-discordance).
    #[arg(long, visible_alias = "diff-indv-discordance")]
    pub sample_discordance: bool,
    /// Write PREFIX.diff.discordance_matrix: genotype of B (rows) by genotype
    /// of A (columns), 0/0, 0/1, 1/1 and ./. (vcftools
    /// --diff-discordance-matrix).
    #[arg(long, visible_alias = "diff-discordance-matrix")]
    pub discordance_matrix: bool,

    /// Which footprint a VCF record covers, which only matters for --region;
    /// the default ref accepts symbolic ALTs without INFO/END.
    #[arg(long, value_enum, default_value = "ref", value_name = "MODE")]
    pub vcf_span: VcfSpan,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
}

/// Write VCF/BCF files one after another (bcftools concat, vcf-concat).
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf concat chr1.vcf.gz chr2.vcf.gz -o all.vcf.gz
  gloci vcf concat chr*.bcf -o all.bcf
  gloci vcf concat --naive chr1.vcf.gz chr2.vcf.gz -o all.vcf.gz")]
pub struct VcfConcatArgs {
    /// Input VCF or BCF files, in output order (`-` for stdin). The samples
    /// must be identical and in the same order, the INFO, FILTER, FORMAT,
    /// ALT and contig definitions must match, and the last record of each
    /// file must sort at or before the first record of the next (exit 1
    /// otherwise; gloci sort fixes the order). The output has the first
    /// file's header. No ##glociCommand line is added, because keeping the
    /// first file's header unchanged lets --naive and checked output stay
    /// identical.
    #[arg(value_name = "FILES", required = true, num_args = 1..)]
    pub files: Vec<String>,
    /// Copy the BGZF blocks of the inputs without recompressing them. Needs
    /// BGZF-compressed VCF files with byte-identical headers (exit 2
    /// otherwise; BCF inputs are not supported) and writes BGZF VCF. The
    /// order is checked across files only, from each file's first and last
    /// record. The output is the same as without --naive, after decompression
    /// (bcftools concat -n).
    #[arg(long)]
    pub naive: bool,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    /// Applies to every input: contig names in each header and the records
    /// are renamed before the headers are compared, the order is checked and
    /// the output is written. Not allowed with --naive, which copies blocks
    /// unchanged.
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// A record type for `vcf filter --types` / `--exclude-types` (spec 7.8).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum VcfTypeArg {
    Snv,
    Mnv,
    Ins,
    Del,
    Indel,
    Sv,
    Ref,
    Other,
}

impl From<VcfTypeArg> for gloci_vcf::VcfType {
    fn from(t: VcfTypeArg) -> Self {
        match t {
            VcfTypeArg::Snv => gloci_vcf::VcfType::Snv,
            VcfTypeArg::Mnv => gloci_vcf::VcfType::Mnv,
            VcfTypeArg::Ins => gloci_vcf::VcfType::Ins,
            VcfTypeArg::Del => gloci_vcf::VcfType::Del,
            VcfTypeArg::Indel => gloci_vcf::VcfType::Indel,
            VcfTypeArg::Sv => gloci_vcf::VcfType::Sv,
            VcfTypeArg::Ref => gloci_vcf::VcfType::Ref,
            VcfTypeArg::Other => gloci_vcf::VcfType::Other,
        }
    }
}

/// The selection options shared by the `vcf` commands: samples, genotype
/// masks, site filters, regions and thinning. `vcf filter` and `vcf extract`
/// apply them through one pipeline stage.
#[derive(Args, Debug, Clone)]
pub struct VcfSelectArgs {
    /// Keep only these sample names (repeatable or comma-separated). They are
    /// written in the input header's order, not the order given here
    /// (bcftools view -s writes them in the order given).
    #[arg(
        long,
        value_delimiter = ',',
        value_name = "SAMPLES",
        conflicts_with = "samples_file"
    )]
    pub samples: Option<Vec<String>>,
    /// File containing sample names to keep, one per line (written in the
    /// input header's order, as with --samples).
    #[arg(long, value_name = "FILE")]
    pub samples_file: Option<PathBuf>,
    /// Remove these sample names (repeatable or comma-separated).
    #[arg(
        long,
        value_delimiter = ',',
        value_name = "SAMPLES",
        conflicts_with = "exclude_samples_file"
    )]
    pub exclude_samples: Option<Vec<String>>,
    /// File containing sample names to remove, one per line.
    #[arg(long, value_name = "FILE")]
    pub exclude_samples_file: Option<PathBuf>,

    /// Mask genotypes with FORMAT/GQ below N.
    #[arg(long, value_name = "N")]
    pub min_gq: Option<u32>,
    /// Mask genotypes with FORMAT/DP below N.
    #[arg(long, value_name = "N")]
    pub min_dp: Option<u32>,
    /// Mask genotypes with FORMAT/DP above N.
    #[arg(long, value_name = "N")]
    pub max_dp: Option<u32>,
    /// Mask genotypes whose FORMAT/FT is not PASS, . or absent.
    #[arg(long)]
    pub mask_filtered_genotypes: bool,

    /// Restrict to regions.
    #[command(flatten)]
    pub region: RegionArgs,
    /// Exclude records whose footprint overlaps a region in FILE (BED).
    #[arg(long, value_name = "FILE")]
    pub exclude_regions_file: Option<PathBuf>,
    /// Keep records on these chromosomes (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',', value_name = "CHROM")]
    pub chrom: Vec<String>,
    /// Remove records on these chromosomes (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',', value_name = "CHROM")]
    pub exclude_chrom: Vec<String>,

    /// Keep records matching (chrom, POS) in FILE.
    #[arg(long, value_name = "FILE")]
    pub positions: Option<PathBuf>,
    /// Remove records matching (chrom, POS) in FILE.
    #[arg(long, value_name = "FILE")]
    pub exclude_positions: Option<PathBuf>,
    /// Keep records whose ID matches a name in FILE.
    #[arg(long, value_name = "FILE")]
    pub ids: Option<PathBuf>,
    /// Remove records whose ID matches a name in FILE.
    #[arg(long, value_name = "FILE")]
    pub exclude_ids: Option<PathBuf>,

    /// Keep records of these variant types (snv, mnv, ins, del, indel, sv, ref, other).
    #[arg(long, value_enum, value_delimiter = ',', value_name = "LIST")]
    pub types: Vec<VcfTypeArg>,
    /// Remove records of these variant types.
    #[arg(long, value_enum, value_delimiter = ',', value_name = "LIST")]
    pub exclude_types: Vec<VcfTypeArg>,
    /// Minimum number of listed alleles (REF + ALTs).
    #[arg(long, value_name = "N")]
    pub min_alleles: Option<usize>,
    /// Maximum number of listed alleles (REF + ALTs).
    #[arg(long, value_name = "N")]
    pub max_alleles: Option<usize>,
    /// Restrict to biallelic sites (equivalent to --min-alleles 2 --max-alleles 2).
    #[arg(long, conflicts_with_all = ["min_alleles", "max_alleles"])]
    pub biallelic: bool,

    /// Minimum QUAL score.
    #[arg(long, value_name = "X")]
    pub min_qual: Option<f64>,
    /// Keep only records whose FILTER is PASS or .
    #[arg(long)]
    pub pass_only: bool,
    /// Keep records containing any of these FILTER tags.
    #[arg(long, value_delimiter = ',', value_name = "TAGS")]
    pub keep_filter: Vec<String>,
    /// Remove records containing any of these FILTER tags.
    #[arg(long, value_delimiter = ',', value_name = "TAGS")]
    pub remove_filter: Vec<String>,

    /// Keep only records that define this INFO tag.
    #[arg(long, value_name = "FLAG")]
    pub require_info: Vec<String>,
    /// Remove records that define this INFO tag.
    #[arg(long, value_name = "FLAG")]
    pub exclude_info: Vec<String>,

    /// Keep records for which the bcftools-style EXPR holds, e.g.
    /// 'QUAL>30 && INFO/DP>10 && N_PASS(GT="het")>=2'. FORMAT fields are limited
    /// to GT, DP and GQ. With --samples it sees only the selected samples
    /// (bcftools view -s -i sees all of them), and a sample masked by --min-gq
    /// and friends reads as missing in it, though the output keeps its values.
    #[arg(long, value_name = "EXPR", conflicts_with = "exclude")]
    pub include: Option<String>,
    /// Remove records for which the bcftools-style EXPR holds (same grammar as
    /// --include).
    #[arg(long, value_name = "EXPR")]
    pub exclude: Option<String>,

    /// Minimum minor allele frequency among called genotypes [0, 1].
    #[arg(long, value_name = "X")]
    pub min_maf: Option<f64>,
    /// Maximum minor allele frequency among called genotypes [0, 1].
    #[arg(long, value_name = "X")]
    pub max_maf: Option<f64>,
    /// Minimum minor allele count among called genotypes.
    #[arg(long, value_name = "N")]
    pub min_mac: Option<u64>,
    /// Maximum minor allele count among called genotypes.
    #[arg(long, value_name = "N")]
    pub max_mac: Option<u64>,

    /// Maximum fraction of missing genotypes [0, 1].
    #[arg(long, value_name = "FRAC")]
    pub max_missing: Option<f64>,
    /// Maximum count of missing genotypes across samples.
    #[arg(long, value_name = "N")]
    pub max_missing_count: Option<u64>,

    /// Minimum mean read depth across samples with depth.
    #[arg(long, value_name = "X")]
    pub min_mean_dp: Option<f64>,
    /// Maximum mean read depth across samples with depth.
    #[arg(long, value_name = "X")]
    pub max_mean_dp: Option<f64>,

    /// Hardy-Weinberg p-value threshold (exact test on biallelic diploid sites).
    #[arg(long, value_name = "P")]
    pub hwe: Option<f64>,
    /// Minimum distance between kept sites in bp (first wins; requires sorted input).
    #[arg(long, value_name = "N")]
    pub thin: Option<u64>,
    /// Keep only sites where every genotype of ploidy 2 or more is written
    /// phased (`|`), called or not: `.|.` counts as phased, `./.` does not
    /// (as vcftools --phased).
    #[arg(long)]
    pub phased_only: bool,

    /// A point on a region's edge overlaps it.
    #[arg(long)]
    pub inclusive_points: bool,
}

/// The matrix `vcf extract --genotype-matrix` writes.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum GenotypeMatrix {
    /// Non-reference allele counts, vcftools' --012 files.
    #[value(name = "012")]
    Zero12,
}

/// How `vcf extract --format` lays out the values.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum ExtractLayout {
    /// One row per site, `SAMPLE:KEY` columns.
    Wide,
    /// One row per site and sample.
    Long,
}

/// Extract INFO or FORMAT values as a table, or the 012 genotype matrix.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf extract cohort.vcf.gz --info DP,AF -o info.tsv
  gloci vcf extract cohort.bcf --format GT,DP --layout long --samples s1,s2
  gloci vcf extract cohort.vcf.gz --genotype-matrix 012 --out-prefix cohort --biallelic")]
#[command(group(
    clap::ArgGroup::new("what").required(true).args(["info", "format", "genotype_matrix"])
))]
pub struct VcfExtractArgs {
    /// Input VCF or BCF file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,

    #[command(flatten)]
    pub select: VcfSelectArgs,

    /// Write these INFO keys as a TSV, one row per site: CHROM POS ID REF ALT
    /// then a column per key (comma-separated). A missing value is `.`, a
    /// flag is 1 when set and 0 when not.
    #[arg(long, value_delimiter = ',', value_name = "KEYS")]
    pub info: Option<Vec<String>>,
    /// Write these FORMAT keys (any key, as text) as a TSV; see --layout.
    /// Keys are comma-separated; a missing value is `.`.
    #[arg(long, value_delimiter = ',', value_name = "KEYS")]
    pub format: Option<Vec<String>>,
    /// Layout of --format: wide is one row per site (CHROM POS ID REF ALT,
    /// then SAMPLE:KEY columns); long is one row per site and sample (CHROM
    /// POS ID REF ALT SAMPLE, then a column per key). Default: wide.
    #[arg(long, value_enum, value_name = "LAYOUT")]
    pub layout: Option<ExtractLayout>,
    /// Write the genotype matrix to --out-prefix.012 (a row per sample: its
    /// index, then the count of non-reference alleles at each site, -1 when
    /// missing), --out-prefix.012.indv (sample names) and --out-prefix.012.pos
    /// (CHROM and POS of each site). Multi-allelic sites are left out, as in
    /// vcftools --012. Memory stays within --max-mem: the codes spill to a
    /// temporary file in --tmp-dir beyond it and the rows are written in
    /// blocks of samples. The budget cannot bound a single row: one site's
    /// codes (a sample count) or one sample's row (a site count) wider than
    /// half the budget is held whole, so memory then exceeds the budget by
    /// that row.
    #[arg(long, value_enum, value_name = "KIND")]
    pub genotype_matrix: Option<GenotypeMatrix>,
    /// Prefix of the --genotype-matrix files.
    #[arg(long, value_name = "PREFIX")]
    pub out_prefix: Option<PathBuf>,

    /// Write the table to FILE (.tsv or .tsv.gz) instead of stdout.
    #[arg(short, long, value_name = "FILE", conflicts_with = "genotype_matrix")]
    pub output: Option<PathBuf>,
    /// BGZF compression level for a .gz output, 0-12.
    #[arg(long, value_name = "N", default_value = "6", value_parser = parse_level)]
    pub compression_level: u8,
    /// Memory budget of --genotype-matrix, in binary units (K, M, G, T); at
    /// least 16M.
    #[arg(long, value_name = "SIZE", default_value = "1G", value_parser = parse_size)]
    pub max_mem: u64,
    #[command(flatten)]
    pub tmp_dir: TmpDirArg,

    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
}

// The "at least one statistic" rule is `check_request`'s, not a clap group:
// its message lists every flag.
/// Compute statistics in one pass: per-site, per-sample and whole-file,
/// written as PREFIX.<ext> in vcftools' formats, plus a summary.
///
/// FORMAT fields are validated lazily: only the fields a requested
/// statistic or a mask rule reads (GT, and DP, GQ or FT when a flag uses
/// them) are decoded, so a malformed DP, GQ or FT value that nothing
/// reads is not an error.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf stats cohort.vcf.gz --out-prefix qc --freq --missing-indv --het --summary
  gloci vcf stats cohort.bcf --out-prefix qc --depth --site-mean-depth --min-dp 8 --samples s1,s2
  gloci vcf stats cohort.vcf.gz --out-prefix qc --tstv-summary --json --pass-only")]
pub struct VcfStatsArgs {
    /// Input VCF or BCF file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,

    #[command(flatten)]
    pub select: VcfSelectArgs,

    /// Prefix of the output files: each statistic writes PREFIX.<extension>
    /// (the extensions vcftools uses: .frq, .idepth, ...). Files are
    /// written atomically and only when the whole run succeeds; an existing
    /// file of the same name is replaced.
    #[arg(long, value_name = "PREFIX", required = true)]
    pub out_prefix: PathBuf,

    /// Allele frequencies per site: PREFIX.frq (vcftools --freq).
    #[arg(long)]
    pub freq: bool,
    /// Allele counts per site: PREFIX.frq.count (vcftools --counts).
    #[arg(long)]
    pub counts: bool,
    /// Mean depth per sample: PREFIX.idepth (vcftools --depth).
    #[arg(long)]
    pub depth: bool,
    /// Total depth per site: PREFIX.ldepth (vcftools --site-depth).
    #[arg(long)]
    pub site_depth: bool,
    /// Mean depth per site: PREFIX.ldepth.mean (vcftools --site-mean-depth).
    #[arg(long)]
    pub site_mean_depth: bool,
    /// Depth of every genotype: PREFIX.gdepth (vcftools --geno-depth).
    #[arg(long)]
    pub geno_depth: bool,
    /// QUAL per site: PREFIX.lqual (vcftools --site-quality).
    #[arg(long)]
    pub site_quality: bool,
    /// Missing data per site: PREFIX.lmiss (vcftools --missing-site).
    #[arg(long)]
    pub missing_site: bool,
    /// Missing data per sample: PREFIX.imiss (vcftools --missing-indv).
    #[arg(long)]
    pub missing_indv: bool,
    /// Inbreeding coefficient per sample: PREFIX.het (vcftools --het).
    #[arg(long)]
    pub het: bool,
    /// Hardy-Weinberg test per site: PREFIX.hwe (vcftools --hardy).
    #[arg(long)]
    pub hardy: bool,
    /// Singleton and doubleton alleles: PREFIX.singletons (vcftools
    /// --singletons).
    #[arg(long)]
    pub singletons: bool,
    /// Records per FILTER value: PREFIX.FILTER.summary (vcftools
    /// --FILTER-summary).
    #[arg(long, visible_alias = "FILTER-summary")]
    pub filter_summary: bool,
    /// Substitution types and Ts/Tv: PREFIX.TsTv.summary (vcftools
    /// --TsTv-summary).
    #[arg(long, visible_alias = "TsTv-summary")]
    pub tstv_summary: bool,
    /// Histogram of indel lengths: PREFIX.indel.hist (vcftools
    /// --hist-indel-len).
    #[arg(long, visible_alias = "hist-indel-len")]
    pub indel_hist: bool,
    /// SNPs per BIN bases: PREFIX.snpden (vcftools --SNPdensity). Needs
    /// input sorted by chromosome and position.
    #[arg(long, value_name = "BIN", visible_alias = "SNPdensity")]
    pub snp_density: Option<u64>,
    /// A short overview in PREFIX.summary.txt: records by type, Ts/Tv,
    /// samples, mean depth and missingness. A record's type is that of its
    /// ALT alleles when they agree: snp, mnp or indel (a multi-allelic record
    /// of SNPs only is snp). Mixed ALTs, such as a SNP with an indel, are
    /// other, and other also takes structural variants (sv is not a type of
    /// its own). A `<*>` ALT is ignored beside a real one; a record whose
    /// only ALT is `<*>` or `.` is ref. Ts/Tv counts biallelic single-base
    /// SNPs only.
    #[arg(long)]
    pub summary: bool,
    /// The same overview as JSON in PREFIX.summary.json. Keys are sorted;
    /// "schema_version" is 1. Layout: counts_by_type {indel, mnp, other, ref,
    /// snp}, depth {mean, values}, genotypes {filtered, missing,
    /// missing_fraction, total}, multiallelic_sites, records, samples,
    /// schema_version, ts_tv {ratio, transitions, transversions}. A ratio
    /// that has no denominator is null.
    #[arg(long)]
    pub json: bool,

    /// Which footprint a VCF record covers: the allele span or the REF span.
    /// The statistics do not depend on it, so the default is ref, which
    /// accepts symbolic ALTs without INFO/END; it only matters for --region.
    #[arg(long, value_enum, default_value = "ref", value_name = "MODE")]
    pub vcf_span: VcfSpan,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
}

/// Filter VCF/BCF files with the selection pipeline.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf filter cohort.bcf --pass-only --min-qual 30 --min-maf 0.01 --max-missing 0.1 -o filtered.bcf
  gloci vcf filter input.vcf --samples s1,s2 --min-dp 10 --recalc-ac -o out.vcf
  gloci vcf filter variants.bcf --types indel --biallelic")]
pub struct VcfFilterArgs {
    /// Input VCF or BCF file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,

    #[command(flatten)]
    pub select: VcfSelectArgs,

    /// Keep only these INFO keys in output records and header.
    #[arg(long, value_delimiter = ',', value_name = "KEYS")]
    pub keep_info: Option<Vec<String>>,
    /// Remove these INFO keys from output records and header.
    #[arg(long, value_delimiter = ',', value_name = "KEYS")]
    pub drop_info: Vec<String>,
    /// Recalculate INFO/AC, AN and AF after sample selection and masking.
    #[arg(long)]
    pub recalc_ac: bool,

    #[command(flatten)]
    pub format: FormatArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub out: OutputArgs,
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

// P3b Task 20

/// Which multi-allelic records `vcf split` splits.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum SplitTypesArg {
    Snps,
    Indels,
    Both,
    Any,
}

/// Split multi-allelic records into biallelic ones (bcftools norm -m-).
///
/// A record with N ALT alleles becomes N records, one per ALT, in ALT
/// order; ID, QUAL and FILTER are copied. In each, the kept ALT is
/// allele 1 in GT and every other allele is 0 (missing stays missing;
/// phase and ploidy are kept). INFO and FORMAT fields follow their
/// declared Number: A keeps that ALT's value, R keeps REF's and that
/// ALT's, G keeps that ALT's genotype values (haploid and diploid);
/// other fields are copied. INFO and FORMAT keys the header does not
/// declare are dropped from split records. Alleles are not trimmed or
/// left-aligned. A record with one ALT is written as it is. A field
/// with the wrong number of values for its Number is an error (exit 1).
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf split calls.vcf.gz -o calls.split.vcf.gz
  gloci vcf split calls.bcf --types snps -o snps-split.bcf
  gloci vcf split calls.vcf.gz | gloci subtract -a - -b known.vcf.gz --match alleles")]
pub struct VcfSplitArgs {
    /// Input VCF or BCF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Split only records with an ALT of this type: snps (a record with a
    /// SNP ALT), indels (an insertion or deletion ALT), both or any (every
    /// multi-allelic record; bcftools splits both and any alike). A split
    /// record is split into all its ALTs; other records are written as they
    /// are.
    #[arg(long, value_enum, default_value = "any", value_name = "TYPES")]
    pub types: SplitTypesArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

// P3b Task 21

/// `--check-ref` modes (bcftools `-c`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum CheckRefArg {
    /// Fail on the first disagreement (the default).
    #[value(name = "e")]
    E,
    /// Warn and write the record as it is.
    #[value(name = "w")]
    W,
    /// Drop the record.
    #[value(name = "x")]
    X,
    /// Set REF from the FASTA (an ALT that equals it is swapped into REF).
    #[value(name = "s")]
    S,
}

/// `--rm-dup` types (bcftools `-d`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum RmDupArg {
    Snps,
    Indels,
    Both,
    All,
    Exact,
}

/// Left-align and normalise records against a reference (bcftools norm -f).
///
/// Every indel, and every record whose REF and ALT share a first or last
/// base, is trimmed and moved to its leftmost position, as `bcftools
/// norm -f REF` does without -a: the bases all alleles end with are
/// removed (the alleles grow to the left from the reference when one
/// would be empty), then the bases they all start with. Records that
/// move are re-sorted within --window bases and written in position
/// order; INFO/END follows the new REF. REF is checked against the FASTA
/// first (--check-ref). The FASTA needs a .fai (and a .gzi for BGZF):
/// run `gloci index REF` once; gloci never builds one itself. Only
/// windows of the reference are read, never whole chromosomes. A record
/// with any symbolic ALT (`<DEL>`, `<DUP>`, ...) is not trimmed or moved.
/// --skip-invalid also skips records that fail --check-ref e. `bcftools
/// norm -a`, `-m+`, `--old-rec-tag` and `-N` are not available.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf norm calls.vcf.gz -f ref.fa -o calls.norm.vcf.gz
  gloci vcf norm calls.bcf --fasta ref.fa --split -o calls.norm.bcf
  gloci vcf norm calls.vcf.gz -f ref.fa --check-ref x --rm-dup exact | gloci vcf stats -
  gloci vcf norm calls.vcf.gz -f ref.fa --check-ref s -o fixed.vcf.gz")]
pub struct VcfNormArgs {
    /// Input VCF or BCF (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Reference FASTA, with its .fai (and .gzi if BGZF-compressed); create
    /// them with `gloci index REF`.
    #[arg(short = 'f', long, visible_alias = "fasta-ref", value_name = "REF")]
    pub fasta: PathBuf,
    /// Split multi-allelic records first, as `bcftools norm -m-TYPES -f`:
    /// snps, indels, both or any (the default; bcftools splits both and any
    /// alike). The same split as `gloci vcf split`.
    #[arg(
        long,
        value_enum,
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "any",
        value_name = "TYPES"
    )]
    pub split: Option<SplitTypesArg>,
    /// What to do with a record whose REF is not the FASTA's (or has a base
    /// other than ACGTN, or an ALT equal to REF): e (fail, the default), w
    /// (warn, write the record as it is), x (drop it) or s (set REF from the
    /// FASTA; an ALT that equals the FASTA is swapped into REF and GT is
    /// recoded).
    #[arg(
        short = 'c',
        long,
        value_enum,
        default_value = "e",
        value_name = "MODE"
    )]
    pub check_ref: CheckRefArg,
    /// Remove duplicate records after normalising, keeping the first: snps
    /// (a SNP or MNP record at the position of a kept one), indels, both
    /// (each type on its own), all (any record at a kept record's position)
    /// or exact (the same REF and ALT alleles).
    #[arg(short = 'd', long, value_enum, value_name = "TYPES")]
    pub rm_dup: Option<RmDupArg>,
    /// Records that move left are re-sorted within this many bases (bcftools
    /// --site-win). A window smaller than the largest shift can order records
    /// differently from bcftools; the default 1000 matches it.
    #[arg(short = 'w', long, default_value = "1000", value_name = "BASES", value_parser = clap::value_parser!(u64).range(1..))]
    pub window: u64,
    // bcftools norm options gloci does not offer: each is a usage error that
    // names the option.
    #[arg(short = 'a', long = "atomize", hide = true)]
    pub atomize: bool,
    #[arg(
        short = 'm',
        long = "multiallelics",
        hide = true,
        value_name = "TYPE",
        allow_hyphen_values = true
    )]
    pub multiallelics: Option<String>,
    #[arg(long = "old-rec-tag", hide = true, value_name = "STR")]
    pub old_rec_tag: Option<String>,
    #[arg(short = 'N', long = "do-not-normalize", hide = true)]
    pub do_not_normalize: bool,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
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

// P3b Task 22

/// Transfer, remove and add fields of a VCF/BCF (bcftools annotate).
///
/// --annotations SRC names a VCF/BCF, a tab-delimited file (plain or
/// BGZF) or a BED file; --columns says what to take from it. FILE and SRC
/// must both be sorted in one chromosome order (FILE's header contigs,
/// else --chrom-order); both are read as a stream, no index is needed
/// and neither is held in memory. A record of a VCF SRC annotates the
/// FILE record at the same CHROM and POS that shares an allele with it
/// (bcftools' default --pair-logic); a tab file with REF and ALT columns
/// matches the same way, one with CHROM,POS matches the records whose
/// REF span covers the position, and a BED file or CHROM,FROM,TO columns
/// match overlapping records (--min-overlap restricts them). --columns
/// takes ID, QUAL, FILTER, INFO/TAG (or TAG), NEW:=OLD, INFO (every tag
/// of a VCF SRC) and ^INFO/TAG lists; +TAG fills only missing values;
/// `.TAG` carries a tab file's missing values over; `-` skips a tab
/// column. A tag the FILE header lacks is declared from a VCF SRC, or
/// must come from --header-lines. --remove drops ID, QUAL, FILTER, INFO,
/// FORMAT (GT stays), INFO/TAG, FORMAT/TAG or FILTER/NAME, and their
/// header lines; ^ keeps only the listed tags of its kind. --set-id
/// writes IDs from %CHROM, %POS, %REF, %FIRST_ALT, %ID and %INFO/TAG
/// (a + prefix fills only missing IDs). --mark-sites +TAG flags the
/// records that have a match, -TAG those that do not. INFO values are
/// copied as written, not reformatted; Number=A and R values follow
/// their allele, not their position. Not available: -s/-S, -i/-e,
/// --rename-chrs, --merge-logic, --single-overlaps, --rename-annots,
/// FORMAT transfers.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf annotate calls.vcf.gz -a dbsnp.vcf.gz -c ID,INFO/RS -o calls.annot.vcf.gz
  gloci vcf annotate calls.vcf.gz -a genes.bed.gz -c CHROM,FROM,TO,GENE --header-lines genes.hdr
  gloci vcf annotate calls.vcf.gz -x INFO,FORMAT/PL --mark-sites +CLEAN
  gloci vcf annotate calls.vcf.gz -I '+%CHROM:%POS:%REF:%FIRST_ALT' -o ids.vcf.gz")]
pub struct VcfAnnotateArgs {
    /// Input VCF or BCF (`-` for stdin), sorted.
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The annotation source: a VCF/BCF, a tab-delimited file (plain or
    /// BGZF) or a BED file, sorted in FILE's chromosome order. A tab file
    /// needs --columns; a BED file is a tab file named *.bed, *.bed.gz or
    /// *.bed.bgz (0-based starts).
    #[arg(short = 'a', long, value_name = "SRC")]
    pub annotations: Option<String>,
    /// What to take from --annotations, comma-separated: ID, QUAL, FILTER,
    /// INFO/TAG or TAG, NEW:=OLD (a VCF source's OLD stored as NEW), INFO
    /// (every tag of a VCF source), ^INFO/TAG (INFO without it), +TAG (only
    /// fill missing values), .TAG (carry a tab file's missing values) and
    /// `-` (skip a tab column). A tab file's list also names its columns:
    /// CHROM, POS (or FROM and TO) and optionally REF and ALT, in file order,
    /// e.g. CHROM,POS,REF,ALT,ID,INFO/AF.
    #[arg(short = 'c', long, value_name = "LIST", allow_hyphen_values = true)]
    pub columns: Option<String>,
    /// Remove, comma-separated: ID, QUAL, FILTER, INFO, FORMAT (GT is kept),
    /// INFO/TAG, FORMAT/TAG, FILTER/NAME, or ##key for a header line. A
    /// leading ^ (^INFO/TAG, ^FORMAT/TAG, ^FILTER/NAME) keeps only the
    /// listed tags of that kind and removes the rest. Their header lines go
    /// too. Removals run before anything is transferred.
    #[arg(short = 'x', long, value_name = "LIST")]
    pub remove: Option<String>,
    /// Append the header lines of FILE (## lines such as ##INFO=<...>) to the
    /// header; they can declare a tag --columns then uses.
    #[arg(long, value_name = "FILE")]
    pub header_lines: Option<PathBuf>,
    /// Write IDs from a format of %CHROM, %POS, %REF, %FIRST_ALT, %ID and
    /// %INFO/TAG; a leading + sets only missing IDs, e.g.
    /// '+%CHROM:%POS:%REF:%FIRST_ALT'. A tag name is letters, digits, _ and
    /// . only: write %POS:%REF, not %POS_%REF.
    #[arg(
        short = 'I',
        long,
        value_name = "[+]FORMAT",
        allow_hyphen_values = true
    )]
    pub set_id: Option<String>,
    /// Add INFO/TAG as a flag to the records that have a match in
    /// --annotations (+TAG, or TAG) or that have none (-TAG); without
    /// --annotations every record gets it.
    #[arg(short = 'm', long, value_name = "[+-]TAG", allow_hyphen_values = true)]
    pub mark_sites: Option<String>,
    /// For a BED or FROM,TO annotation file: the fraction of the annotation
    /// (ANN), of the record (:VCF) or of both (ANN:VCF) that must overlap.
    #[arg(long, value_name = "ANN[:VCF]")]
    pub min_overlap: Option<String>,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    // bcftools annotate options gloci does not offer: each is a usage error
    // that names the option.
    #[arg(long = "rename-chrs", hide = true, value_name = "FILE")]
    pub rename_chrs: Option<String>,
    #[arg(short = 's', long = "samples", hide = true, value_name = "LIST")]
    pub samples: Option<String>,
    #[arg(short = 'S', long = "samples-file", hide = true, value_name = "FILE")]
    pub samples_file: Option<String>,
    #[arg(short = 'i', long = "include", hide = true, value_name = "EXPR")]
    pub include: Option<String>,
    #[arg(short = 'e', long = "exclude", hide = true, value_name = "EXPR")]
    pub exclude: Option<String>,
    #[arg(short = 'l', long = "merge-logic", hide = true, value_name = "LOGIC")]
    pub merge_logic: Option<String>,
    #[arg(long = "single-overlaps", hide = true)]
    pub single_overlaps: bool,
    #[arg(long = "rename-annots", hide = true, value_name = "FILE")]
    pub rename_annots: Option<String>,
    #[arg(long = "pair-logic", hide = true, value_name = "LOGIC")]
    pub pair_logic: Option<String>,
    #[arg(short = 'C', long = "columns-file", hide = true, value_name = "FILE")]
    pub columns_file: Option<String>,
    #[arg(short = 'H', long = "header-line", hide = true, value_name = "STR")]
    pub header_line: Option<String>,
    #[arg(short = 'k', long = "keep-sites", hide = true)]
    pub keep_sites: bool,
    #[arg(short = 'r', long = "regions", hide = true, value_name = "REGION")]
    pub regions: Option<String>,
    #[arg(short = 'R', long = "regions-file", hide = true, value_name = "FILE")]
    pub regions_file: Option<String>,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

// P3b Task 23

/// `bcftools merge -m`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum MergeModeArg {
    None,
    Snps,
    Indels,
    Both,
    All,
    Id,
}

/// `bcftools merge -F`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum FilterLogicArg {
    /// Records are merged with all their filters (`+`).
    #[value(name = "+")]
    All,
    /// A PASS in any input removes the filters (`x`).
    #[value(name = "x")]
    Pass,
}

/// Combine the samples of several VCF/BCF files (bcftools merge).
///
/// The inputs must be sorted in one chromosome order (the first input's
/// ##contig lines, else --chrom-order); they are read as a stream, so no
/// index is needed and only the records of one position are held. Records
/// at the same CHROM and POS are merged as --merge says: ALT alleles are
/// united in input order (a shorter REF is extended to the longest), and
/// GT and Number=A, R and G FORMAT values are re-indexed to the merged
/// alleles. A sample absent from a record is ./. (or 0/0 with -0) in the
/// record's ploidy. QUAL is the largest, FILTER the union without PASS,
/// ID the distinct IDs joined with ;. Sample names must be distinct
/// (--force-samples renames a repeat to N:NAME). The header is the first
/// input's with the other inputs' new lines added. The inputs are open
/// together: merging more files than `ulimit -n` allows fails with the
/// operating system's error; raise the limit (ulimit -n 4096). bcftools
/// merge -g, -L, -r, -R, -T, --print-header and --use-header are not
/// available.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf merge a.vcf.gz b.vcf.gz c.vcf.gz -o cohort.vcf.gz
  gloci vcf merge --file-list shards.txt -m none --force-samples -o cohort.bcf
  gloci vcf merge a.bcf b.bcf -i DP:sum,AF:avg -0 | gloci vcf stats - --out-prefix qc
  gloci vcf merge a.vcf.gz b.vcf.gz -f PASS,. -m id -o by-id.vcf.gz")]
pub struct VcfMergeArgs {
    /// Input VCF or BCF files (two or more; `-` is not allowed), each
    /// sorted in the same chromosome order. Or use --file-list.
    #[arg(value_name = "FILE")]
    pub files: Vec<String>,
    /// Read the input file names from FILE, one per line (bcftools -l).
    #[arg(short = 'l', long, value_name = "FILE")]
    pub file_list: Option<PathBuf>,
    /// Which records at one position are merged into one (bcftools -m):
    /// none (only records with the same alleles), snps, indels, both (SNPs
    /// together and indels together; the default), all (any records) or id
    /// (records with the same ID).
    #[arg(
        short = 'm',
        long,
        value_enum,
        default_value = "both",
        value_name = "MODE"
    )]
    pub merge: MergeModeArg,
    /// Rename a repeated sample name to N:NAME (N is the number of its
    /// file) instead of failing.
    #[arg(long)]
    pub force_samples: bool,
    /// How INFO tags are merged, TAG:METHOD with METHOD one of sum, avg,
    /// min, max, join, comma-separated; `-` for none. The default is
    /// DP:sum,DP4:sum (and AC:sum,AN:sum when there are no samples; with
    /// samples AC and AN are recomputed from the merged genotypes). A
    /// tag with no rule takes the first record's value; Number=A and R tags
    /// are re-indexed to the merged alleles.
    #[arg(short = 'i', long, value_name = "RULES", allow_hyphen_values = true)]
    pub info_rules: Option<String>,
    /// Fill a sample absent from a record with 0/0 instead of ./. (bcftools -0).
    #[arg(short = '0', long)]
    pub missing_to_ref: bool,
    /// Skip, before merging, the records whose FILTER has none of these
    /// names (comma-separated; `.` is a missing FILTER, PASS is PASS).
    #[arg(short = 'f', long, value_name = "LIST")]
    pub apply_filters: Option<String>,
    /// What to do with the FILTER of merged records: + (all filters kept,
    /// PASS dropped when another filter is set; the default) or x (PASS in
    /// any input makes the result PASS).
    #[arg(
        short = 'F',
        long,
        value_enum,
        default_value = "+",
        value_name = "LOGIC"
    )]
    pub filter_logic: FilterLogicArg,
    /// Accepted and ignored: gloci merges sorted inputs as a stream and
    /// never needs an index (bcftools merge needs one unless --no-index).
    #[arg(long, alias = "force-no-index", hide = true)]
    pub no_index: bool,
    /// Accepted: merge a single input (bcftools merge --force-single).
    #[arg(long, hide = true)]
    pub force_single: bool,
    // bcftools merge options gloci does not offer: each is a usage error
    // that names the option.
    #[arg(short = 'g', long = "gvcf", hide = true, value_name = "REF")]
    pub gvcf: Option<String>,
    #[arg(short = 'L', long = "local-alleles", hide = true, value_name = "N")]
    pub local_alleles: Option<String>,
    #[arg(long = "print-header", hide = true)]
    pub print_header: bool,
    #[arg(long = "use-header", hide = true, value_name = "FILE")]
    pub use_header: Option<String>,
    #[arg(short = 'r', long = "regions", hide = true, value_name = "REGION")]
    pub regions: Option<String>,
    #[arg(short = 'R', long = "regions-file", hide = true, value_name = "FILE")]
    pub regions_file: Option<String>,
    #[arg(short = 'T', long = "targets-file", hide = true, value_name = "FILE")]
    pub targets_file: Option<String>,
    #[arg(long = "merge-logic", hide = true, value_name = "LOGIC")]
    pub merge_logic: Option<String>,
    #[arg(short = 'M', long = "missing-rules", hide = true, value_name = "RULES")]
    pub missing_rules: Option<String>,
    #[arg(long = "regions-overlap", hide = true, value_name = "0|1|2")]
    pub regions_overlap: Option<String>,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
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
