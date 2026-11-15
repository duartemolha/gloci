use super::EXPERIMENTAL;
use super::shared::{
    ChromMapArg, ChromOrderArg, FormatArg, GenomeArg, InvalidArg, OutputArgs, RegionArgs,
    ThreadArgs, TmpDirArg, VcfSpan, VcfSpanArg, parse_level, parse_size,
};
use super::vcf_edit::{VcfAnnotateArgs, VcfConcatArgs, VcfMergeArgs, VcfNormArgs, VcfSplitArgs};
use clap::{Args, Subcommand, ValueEnum};
use std::path::PathBuf;

/// VCF and BCF operations: filter, extract, stats, concat, split, norm, annotate, merge, compare
/// and popgen (bcftools and vcftools equivalents).
#[derive(Args, Debug, Clone)]
#[command(
    about = "VCF and BCF operations (bcftools and vcftools equivalents)",
    long_about
)]
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
    #[command(before_help = EXPERIMENTAL)]
    Popgen(Box<VcfPopgenArgs>),
}

/// `vcf compare --match`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum CompareMatchArg {
    Alleles,
    Position,
}

/// Compare the genotypes of two VCF/BCF files: per site, per sample and as a matrix (vcftools
/// --diff-site-discordance, --diff-indv-discordance, --diff-discordance-matrix).
///
/// With none of --site-discordance, --sample-discordance and --discordance-matrix, all three files
/// are written, atomically and only when the whole run succeeds: PREFIX.diff.sites (CHROM POS FILES
/// MATCHING_ALLELES N_COMMON_CALLED N_DISCORD DISCORDANCE, one row per site in both files, and per
/// site in one file with zeros and nan), PREFIX.diff.indv (every sample of the union; one in a single
/// file is 0, 0, nan) and PREFIX.diff.discordance_matrix (biallelic pairs of the same site: with
/// --match alleles the same trimmed allele, so AT>A and ATG>AG count; with --match position the same
/// REF and ALT text; haploid genotypes are not counted).
#[derive(Args, Debug, Clone)]
#[command(about = "Genotype discordance (vcftools --diff-*-discordance)")]
#[command(after_long_help = "\
Examples:
  gloci vcf compare truth.vcf.gz calls.vcf.gz --out-prefix cmp
  gloci vcf compare a.bcf b.bcf --out-prefix cmp --sample-discordance --match position
  gloci vcf compare a.vcf.gz b.vcf.gz --out-prefix cmp --sample-map b-to-a.txt --samples s1,s2 --min-dp 10

Pairing and counting:
Both files are read as sorted streams (the chromosome order of -g, else A's ##contig \
lines, else --chrom-order; an out-of-order record is exit 1 and names the file), so no \
index is needed. Sites are paired by --match: alleles (the default) pairs records at one \
CHROM and POS that share an ALT allele after trimming, position pairs as vcftools does \
(the same POS with the same REF, where a REF of N or . takes the other's). Discordance is \
counted over the sites in both files and the samples in both files (the union by name, \
after --sample-map renames B's samples). --thin needs both files to declare the same \
##contig lines and no -g (exit 2 otherwise).

Threads:
A and B are always read on threads of their own, plus the thread that walks them, so -t 1 \
and -t 2 run three threads; more -t adds selection workers up to -t in all.

Genotypes and samples:
A genotype is its unordered alleles, so 0/1 equals 1/0 and 0|1, and 0/1 of A>C,G equals \
0/2 of A>G,C; any missing allele makes it missing (vcftools counts ./1 as called). \
--samples and --exclude-samples take names present in either file and apply to each file \
that has them; a name in neither is exit 2.

Not available:
vcftools --diff-site and --diff-indv (use gloci intersect and subtract \
--match alleles for sites, and the #CHROM lines for samples).")]
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

/// Extract INFO or FORMAT values as a table, or the 012 genotype matrix (vcftools --012).
#[derive(Args, Debug, Clone)]
#[command(
    about = "Tabulate INFO/FORMAT values or genotypes (vcftools --012)",
    long_about
)]
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
#[command(about = "Per-site, per-sample and file statistics (vcftools)")]
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

// The "at least one statistic" rule is `check_popgen_request`'s, not a clap
// group: its message lists every flag.
/// Pi, Tajima's D, Fst (vcftools --*-pi, --TajimaD, --weir-fst-pop)
///
/// Nucleotide diversity per site (--site-pi) and in windows (--window-pi), Tajima's D in bins
/// (--tajima-d) and Weir and Cockerham's Fst between populations (--fst-pop, per site, in
/// windows and over the run), as vcftools 0.1.16 computes them over the selected samples. Every
/// requested statistic is computed in one pass over the records that pass the selection flags
/// (samples, genotype masks, site filters, regions), and each writes PREFIX.<extension> with
/// vcftools' file name and columns. The files are written atomically and only when the whole
/// run succeeds.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci vcf popgen cohort.vcf.gz --out-prefix pg --site-pi --window-pi 10000 --tajima-d 10000
  gloci vcf popgen cohort.vcf.gz --out-prefix pg --window-pi 50000 --window-pi-step 10000 --min-dp 8
  gloci vcf popgen cohort.bcf --out-prefix fst --fst-pop yri.txt --fst-pop ceu.txt --fst-pop chb.txt
  gloci vcf popgen cohort.bcf --out-prefix fst --fst-pop yri --fst-pop ceu --fst-window-size 50000

Sites:
  A site counts only when every selected sample that is not masked is diploid (for Fst, every
  population member). A half-called genotype (./1) counts its called allele in pi and Tajima's
  D, and is missing for Fst; a masked genotype is missing. Site and window pi and Fst take every
  allele of a multi-allelic site; Tajima's D takes biallelic sites only.

Window pi:
  Windows are [k*STEP+1, k*STEP+SIZE] on each chromosome. A polymorphic site adds its pairs of
  called alleles, and the pairs that differ, to every window that holds it; every other position
  of the window counts as monomorphic, with all 2N(2N-1) pairs of the N selected samples. pi is
  the differing pairs over all pairs. N_VARIANTS counts the polymorphic sites with a called
  allele other than REF; a window without a polymorphic site has no row.

Tajima's D:
  POS falls in bin floor(POS * (1/BIN)), the product taken in floating point as vcftools does
  (POS 49 with BIN 49 is in bin 0); BIN_START is bin * BIN. The sample size is 2 x the selected
  samples whatever is missing at a site. A chromosome's rows run from its first bin with a
  segregating site to its last bin with a biallelic site; a bin without one is 0 nan.

Fst:
  Each --fst-pop file is a population: the first word of each line names a sample (blank lines
  are skipped). Names are matched against the selected samples, and the populations select no
  samples for the other statistics or the site filters. Fewer than two files, a file that names
  no selected sample, or a sample in two files is exit 2; names that are not selected samples
  are left out, with one warning that counts them. Per site, each population's fully called
  genotypes give vcftools' a, b and c per allele (Weir and Cockerham 1984), and Fst is
  sum(a) / sum(a+b+c) over the alleles; it is nan when a population has no call or nothing
  varies. PREFIX.weir.fst has a row per site; PREFIX.weir.fst.summary gives N_SITES, MEAN_FST and
  WEIGHTED_FST (sum(a) / sum(a+b+c)) over the sites that are not nan, which vcftools only logs.
  With --fst-window-size, PREFIX.windowed.weir.fst sums them over the windows of window pi; a
  window without such a site, or whose sum(a+b+c) is 0, has no row.

Sorted input:
  --window-pi, --tajima-d and --fst-window-size write their rows as the input passes them, so
  they need input sorted by chromosome and position: a chromosome that comes back or a position
  that goes backwards is exit 1, and gloci sort fixes it. --site-pi and --fst-pop without a
  window take records in any order.")]
pub struct VcfPopgenArgs {
    /// Input VCF or BCF file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,

    #[command(flatten)]
    pub select: VcfSelectArgs,

    /// Prefix of the output files: PREFIX.sites.pi, PREFIX.windowed.pi,
    /// PREFIX.Tajima.D, PREFIX.weir.fst, PREFIX.weir.fst.summary and
    /// PREFIX.windowed.weir.fst. An existing file of the same name is
    /// replaced.
    #[arg(long, value_name = "PREFIX", required = true)]
    pub out_prefix: PathBuf,

    /// Nucleotide diversity per site: PREFIX.sites.pi (vcftools --site-pi).
    #[arg(long)]
    pub site_pi: bool,
    /// Nucleotide diversity in windows of SIZE bases: PREFIX.windowed.pi
    /// (vcftools --window-pi). Needs sorted input.
    #[arg(long, value_name = "SIZE")]
    pub window_pi: Option<u64>,
    /// Start a --window-pi window every STEP bases, at most SIZE (vcftools
    /// --window-pi-step). Default: SIZE, windows that do not overlap.
    #[arg(long, value_name = "STEP", requires = "window_pi")]
    pub window_pi_step: Option<u64>,
    /// Tajima's D in bins of BIN bases: PREFIX.Tajima.D (vcftools --TajimaD).
    /// Needs sorted input.
    #[arg(long, value_name = "BIN")]
    pub tajima_d: Option<u64>,
    /// A population for Weir and Cockerham's Fst: a file of sample names,
    /// one per line. Give it once per population, at least twice:
    /// PREFIX.weir.fst and PREFIX.weir.fst.summary (vcftools --weir-fst-pop).
    #[arg(long, value_name = "FILE")]
    pub fst_pop: Vec<PathBuf>,
    /// Fst in windows of BP bases as well: PREFIX.windowed.weir.fst
    /// (vcftools --fst-window-size). Needs sorted input.
    #[arg(long, value_name = "BP", requires = "fst_pop")]
    pub fst_window_size: Option<u64>,
    /// Start an Fst window every BP bases, at most the window size (vcftools
    /// --fst-window-step). Default: the window size, windows that do not
    /// overlap.
    #[arg(long, value_name = "BP", requires = "fst_window_size")]
    pub fst_window_step: Option<u64>,

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

/// Filter VCF/BCF files with the selection pipeline (bcftools view, vcftools).
#[derive(Args, Debug, Clone)]
#[command(about = "Filter VCF/BCF files (bcftools view, vcftools)", long_about)]
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
