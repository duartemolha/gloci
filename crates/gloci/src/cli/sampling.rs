use super::EXPERIMENTAL;
use super::shared::{
    AArg, AlignArgs, ChromMapArg, ChromOrderArg, FormatArg, GenomeArg, HeaderArg, InvalidArg,
    OutputArgs, OverlapArgs, RelationSortArg, ThreadArgs, UnsortedArg, VcfSpanArg, WindowArgs,
    parse_fraction,
};
use clap::{Args, Subcommand, ValueEnum};

/// `--seed` for the commands that draw random numbers (D143).
#[derive(Args, Debug, Clone, Copy)]
pub struct SeedArg {
    /// Seed for the random stream, so a run can be reproduced (default: a fresh seed, which is
    /// printed as a warning).
    #[arg(long, value_name = "N")]
    pub seed: Option<u64>,
}

/// Statistics of interval sets (bedtools jaccard, fisher and reldist).
#[derive(Args, Debug, Clone)]
#[command(
    about = "Interval set statistics (bedtools jaccard, fisher, reldist)",
    long_about
)]
#[command(after_long_help = "\
Examples:
  gloci stats jaccard -a peaks.bed -b genes.bed
  gloci stats jaccard -a rep1.bed -b rep2.bed -f 0.5 -r
  gloci stats fisher -a peaks.bed -b genes.bed -g hg38
  gloci stats reldist -a peaks.bed -b genes.bed --detail")]
pub struct StatsArgs {
    #[command(subcommand)]
    pub command: StatsCommand,
}

#[derive(Subcommand, Debug, Clone)]
pub enum StatsCommand {
    #[command(before_help = EXPERIMENTAL)]
    Jaccard(JaccardArgs),
    #[command(before_help = EXPERIMENTAL)]
    Fisher(FisherArgs),
    #[command(before_help = EXPERIMENTAL)]
    Reldist(ReldistArgs),
}

/// The Jaccard similarity of two interval sets (bedtools jaccard).
///
/// Each input is first merged as `gloci merge` does (overlapping and bookended intervals join; per
/// strand with -s or -S). The intersection is the sum of the overlap bases of every pair (merged A,
/// merged B) that passes -f, -F, -r, -e, --min-bases, -s, -S and --inclusive-points (measured on the
/// merged intervals), the union is the bases covered by either input, and jaccard = intersection /
/// union divided in single precision and printed as %g (six significant digits; `nan` when both are
/// empty). n_intersections counts the accepted pairs. Output is a header line and one row:
/// intersection, union, jaccard, n_intersections.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci stats jaccard -a peaks.bed -b genes.bed
  gloci stats jaccard -a rep1.bed -b rep2.bed -f 0.5 -r
  gloci stats jaccard -a calls.vcf.gz -b regions.bed --vcf-span ref
  gloci stats jaccard -a unsorted-a.bed -b unsorted-b.bed --sort

Points and sorting:
A point (zero-length interval) covers no bases, so it adds nothing to the intersection or \
the union, but a point inside a merged interval of the other input counts as an \
intersection (bedtools widens a point to two bases; see docs/differences.md). Inputs must \
be sorted (--sort sorts them; --unsorted does not apply); a VCF is read as its footprints \
(--vcf-span); BAM, CRAM and SAM reads are read as their reference spans.")]
pub struct JaccardArgs {
    #[command(flatten)]
    pub a: AArg,
    /// B: the file to compare A with (`-` for stdin).
    #[arg(short = 'b', long = "b", value_name = "FILE", required = true)]
    pub b: String,
    /// B's format, overriding detection: bed, bedgraph, gff, vcf or bcf; bam, cram or sam
    /// (alignments are read as intervals, as `intersect` reads them).
    #[arg(long, value_name = "FMT")]
    pub b_format: Option<String>,
    /// Minimum fraction of a merged A interval covered, e.g. 0.5 or 50% (4.3: fractions).
    #[arg(short = 'f', long, value_name = "FRAC", value_parser = parse_fraction)]
    pub min_frac_a: Option<gloci_core::Fraction>,
    /// Minimum fraction of a merged B interval covered, e.g. 0.5 or 50% (4.3: fractions).
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
    /// Merge each strand apart and count only overlaps on the same strand; unknown strands
    /// never match.
    #[arg(short = 's', long = "same-strand", conflicts_with = "opposite_strand")]
    pub same_strand: bool,
    /// Merge each strand apart and count only overlaps on opposite strands; unknown strands
    /// never match.
    #[arg(short = 'S', long = "opposite-strand")]
    pub opposite_strand: bool,
    /// A point on the edge of an interval overlaps it (4.3: points).
    #[arg(long)]
    pub inclusive_points: bool,
    #[command(flatten)]
    pub sort: RelationSortArg,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Fisher's exact test of the overlap of two interval sets (bedtools fisher).
///
/// Counts the records of A and B (nA, nB), their total lengths (LA, LB)
/// and the pairs that overlap under -f, -F, -r, -e, --min-bases, -s, -S
/// and --inclusive-points (n11); with --merge-overlaps each input is
/// merged first. The genome size G is the sum of the -g lengths (so -g is
/// required and must have lengths) less the merged length of the --exclude
/// BED. As bedtools does, the table is n11, n12 = nA - n11, n21 = nB -
/// n11 and n22 = N - n11 - n12 - n21, where N, the number of possible
/// intervals, is estimated as G / (2 + LA/nA + LB/nB) (rounded down, and
/// never fewer than n11 + n12 + n21). The report is bedtools' text:
/// the table in `#` lines, then the left, right and two-tail p-values of
/// klib's exact test (%.5g) and the odds ratio (%.3f, `inf` or `-nan`).
/// Inputs must be sorted in -g's order (--sort sorts A, B and the exclude
/// file); a VCF is read as its footprints (--vcf-span).
#[derive(Args, Debug, Clone)]
#[command(about = "Fisher's exact test of overlap (bedtools fisher)")]
#[command(after_long_help = "\
Examples:
  gloci stats fisher -a peaks.bed -b genes.bed -g hg38
  gloci stats fisher -a peaks.bed -b genes.bed -g sizes.genome --merge-overlaps
  gloci stats fisher -a peaks.bed -b genes.bed -g hg38 --exclude gaps.bed -f 0.5
  gloci stats fisher -a unsorted-a.bed -b unsorted-b.bed -g hg38 --sort")]
pub struct FisherArgs {
    #[command(flatten)]
    pub a: AArg,
    /// B: the file to compare A with (`-` for stdin).
    #[arg(short = 'b', long = "b", value_name = "FILE", required = true)]
    pub b: String,
    /// B's format, overriding detection: bed, bedgraph, gff, vcf or bcf; bam, cram or sam
    /// (alignments are read as intervals, as `intersect` reads them).
    #[arg(long, value_name = "FMT")]
    pub b_format: Option<String>,
    /// Genome with chromosome lengths: a built-in genome such as hg38 (see gloci genomes list),
    /// or a .genome, .fai, .dict, VCF/BCF or BAM header file. Its lengths add up to the genome
    /// size.
    #[arg(short, long, value_name = "NAME|FILE", required = true)]
    pub genome: String,
    /// Merge overlapping and bookended intervals of each input before counting (bedtools -m).
    #[arg(short = 'm', long)]
    pub merge_overlaps: bool,
    /// Leave the merged length of this sorted BED out of the genome size (bedtools -exclude).
    #[arg(long, value_name = "BED")]
    pub exclude: Option<String>,
    /// Minimum fraction of an A interval covered, e.g. 0.5 or 50% (4.3: fractions).
    #[arg(short = 'f', long, value_name = "FRAC", value_parser = parse_fraction)]
    pub min_frac_a: Option<gloci_core::Fraction>,
    /// Minimum fraction of a B interval covered, e.g. 0.5 or 50% (4.3: fractions).
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
    /// Count only overlaps on the same strand; with --merge-overlaps each strand merges apart
    /// and unknown strands never match.
    #[arg(short = 's', long = "same-strand", conflicts_with = "opposite_strand")]
    pub same_strand: bool,
    /// Count only overlaps on opposite strands (bedtools -S needs a + or - here; gloci's is a
    /// flag, as in the other commands).
    #[arg(short = 'S', long = "opposite-strand")]
    pub opposite_strand: bool,
    /// A point on the edge of an interval overlaps it (4.3: points).
    #[arg(long)]
    pub inclusive_points: bool,
    #[command(flatten)]
    pub sort: RelationSortArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// What `pair-intersect` reports (bedtools `-type`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum PairTypeArg {
    /// A with each B hitting either end (default).
    Either,
    /// A alone when neither end hits.
    Neither,
    /// A with each B, when both ends hit.
    Both,
    /// A with each B, when exactly one end hits.
    Xor,
    /// A with each B when one end hits, A alone (padded) when none does.
    Notboth,
    /// A with each B hitting the gap between the ends.
    Ispan,
    /// A with each B hitting the span from the first start to the last end.
    Ospan,
    /// A alone when nothing hits the gap between the ends.
    Notispan,
    /// A alone when nothing hits the outer span.
    Notospan,
}

/// BEDPE pairs by how their ends overlap B (bedtools pairtobed).
///
/// A is a BEDPE file (chrom1 start1 end1 chrom2 start2 end2, then optionally name, score, strand1,
/// strand2 and more; it is never detected, so a BED file of pairs is read as BEDPE only here). B is
/// BED, GFF or VCF (a VCF as its footprints, --vcf-span) and is loaded into memory; A streams in any
/// order. Each present end of a pair is tested against B as an interval under -f, -F, -r, -e,
/// --min-bases, -s, -S and --inclusive-points (section 4.3: a point is inside only strictly, and
/// unknown strands never match -s or -S); -s and -S use strand1 for end 1 and strand2 for end 2. An
/// end whose chromosome is `.` overlaps nothing. --type chooses what is written: A with the B its ends
/// hit, or A alone, by one of the rules either (default), both, xor, neither, notboth, ispan, ospan,
/// notispan or notospan.
#[derive(Args, Debug, Clone)]
#[command(after_long_help = "\
Examples:
  gloci pair-intersect -a sv.bedpe -b genes.bed
  gloci pair-intersect -a sv.bedpe -b exons.bed --type both -f 0.5
  gloci pair-intersect -a sv.bedpe -b repeats.bed --type neither -o clean.bedpe
  gloci pair-intersect -a sv.bedpe -b genes.bed --type xor -s -t 4

Output:
Rows with a B are A's line, a tab and B's line (TSV); the \
A-alone types write A's BEDPE line. either (default): A with each B hitting end 1, then \
each hitting end 2 (a B hitting both ends is written twice). both: the same, when both \
ends hit. xor: the same, when exactly one end hits. neither: A alone, when no end hits. \
notboth: A with each B, when exactly one end hits, and A alone padded with a placeholder B \
(`.`, `-1`, `-1` ...) when none does. ispan and ospan: A with each B overlapping the gap \
between the ends, or the span from the first start to the last end; notispan and notospan: \
A alone when none does. The span types need both ends on one chromosome (a pair with one \
end absent or the ends on two chromosomes is left out, counted in one warning), an empty \
gap (touching or overlapping ends) hits nothing, and -s/-S with them is an error. Hits of \
one A come in B's input order; output follows A's order whatever the thread count. Output \
is TSV or BEDPE (-o x.bedpe, x.tsv, .gz for BGZF); header lines of A are not written.")]
#[command(mut_arg("split", |a| a.hide(true)))]
pub struct PairIntersectArgs {
    /// A: the BEDPE file whose pairs are reported (`-` for stdin).
    #[arg(short = 'a', long = "a", value_name = "FILE", required = true)]
    pub a: String,
    /// B: the BED, GFF or VCF file the ends are tested against (`-` for stdin).
    #[arg(short = 'b', long = "b", value_name = "FILE", required = true)]
    pub b: String,
    /// B's format, overriding detection: bed, bedgraph, gff, vcf or bcf (B does not read
    /// alignments; convert them to BED first).
    #[arg(long, value_name = "FMT")]
    pub b_format: Option<String>,
    /// What to report: either, neither, both, xor, notboth, ispan, ospan, notispan or notospan.
    #[arg(
        long = "type",
        value_enum,
        default_value = "either",
        value_name = "TYPE"
    )]
    pub kind: PairTypeArg,
    #[command(flatten)]
    pub o: OverlapArgs,
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

/// What `pair-pair` reports (bedtools `-type`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum PairPairTypeArg {
    /// A alone when no end of A overlaps an end of B.
    Neither,
    /// A with each B having an end that overlaps an end of A.
    Either,
    /// A with each B overlapping on both ends, in either orientation (default).
    Both,
    /// A alone when no B overlaps on both ends.
    Notboth,
}

/// BEDPE pairs by how their ends overlap the ends of the pairs of B (bedtools pairtopair).
///
/// A and B are both BEDPE files (never detected). B is loaded into memory; A streams in any order.
/// Orientation 1 tests end 1 of A against end 1 of B and end 2 against end 2; orientation 2 tests end
/// 1 of A against end 2 of B and end 2 against end 1. --type both (default): A with each B that
/// overlaps on both ends in either orientation. either: A with each B with any end overlapping an end
/// of A. neither: A alone when no B end overlaps. notboth: A alone when no B overlaps on both ends.
#[derive(Args, Debug, Clone)]
#[command(about = "BEDPE pairs against BEDPE pairs (bedtools pairtopair)")]
#[command(after_long_help = "\
Examples:
  gloci pair-pair -a sv1.bedpe -b sv2.bedpe
  gloci pair-pair -a sv1.bedpe -b sv2.bedpe --type either -w 100
  gloci pair-pair -a sv1.bedpe -b sv2.bedpe --type notboth -o unique.bedpe
  gloci pair-pair -a sv.bedpe -b sv.bedpe --require-different-names -s -t 4

Overlap tests:
Ends are tested under -f, -F, -r, -e, --min-bases and --inclusive-points (section 4.3); \
strand is ignored unless -s or -S, which compare strand1 with strand1 and strand2 with \
strand2 of the tested ends (unknown strands never match). -w, --left, --right and \
--strand-window pad B's ends before testing (not with the thresholds). An end whose \
chromosome is `.` overlaps nothing. --require-different-names drops a B named like A \
(column 7); a pair without a name is then an error.

Output:
A B is written once per A, as its line reads (TSV: A's line, a tab, B's line), in B's \
input order; the A-alone types write A's BEDPE line. Output follows A's order whatever the \
thread count; -o x.bedpe, x.tsv, .gz for BGZF.")]
#[command(mut_arg("split", |a| a.hide(true)))]
pub struct PairPairArgs {
    /// A: the BEDPE file whose pairs are reported (`-` for stdin).
    #[arg(short = 'a', long = "a", value_name = "FILE", required = true)]
    pub a: String,
    /// B: the BEDPE file the pairs of A are tested against (`-` for stdin).
    #[arg(short = 'b', long = "b", value_name = "FILE", required = true)]
    pub b: String,
    /// What to report: neither, either, both or notboth.
    #[arg(long = "type", value_enum, default_value = "both", value_name = "TYPE")]
    pub kind: PairPairTypeArg,
    #[command(flatten)]
    pub o: OverlapArgs,
    #[command(flatten)]
    pub w: WindowArgs,
    /// Drop a B named like A (column 7), as bedtools -rdn; every pair then needs a name.
    #[arg(long)]
    pub require_different_names: bool,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// The relative distance of A midpoints to B midpoints (bedtools reldist).
///
/// For each A interval, finds the two B midpoints that enclose A's
/// midpoint (midpoint = floor of (start + end) / 2) and takes the distance
/// to the nearer one over the distance between the two, from 0 (A sits on
/// a B midpoint) to 0.5 (halfway). The summary bins the values by
/// floor(rel * 100) (so 0.29 lands in bin 0.28, as in bedtools) and prints
/// reldist, count, total and fraction. An A midpoint before the first or
/// after the last B midpoint of its chromosome, or on a chromosome with
/// fewer than two B midpoints, is left out, as bedtools does. With
/// --detail each qualifying A record is printed with its value (%.3f)
/// appended. A VCF is read as its footprints (--vcf-span). Inputs must be
/// sorted (--sort sorts them; --unsorted loads every B midpoint instead).
#[derive(Args, Debug, Clone)]
#[command(about = "Relative distance of A to B midpoints (bedtools reldist)")]
#[command(after_long_help = "\
Examples:
  gloci stats reldist -a peaks.bed -b genes.bed
  gloci stats reldist -a peaks.bed -b genes.bed --detail
  gloci stats reldist -a unsorted-a.bed -b unsorted-b.bed --sort
  gloci stats reldist -a peaks.bed -b genes.bed --unsorted")]
pub struct ReldistArgs {
    #[command(flatten)]
    pub a: AArg,
    /// B: the file to compare A with (`-` for stdin).
    #[arg(short = 'b', long = "b", value_name = "FILE", required = true)]
    pub b: String,
    /// B's format, overriding detection: bed, bedgraph, gff, vcf or bcf; bam, cram or sam
    /// (alignments are read as intervals, as `intersect` reads them).
    #[arg(long, value_name = "FMT")]
    pub b_format: Option<String>,
    /// Print each A record with its relative distance (%.3f) appended, instead of the summary.
    #[arg(long)]
    pub detail: bool,
    #[command(flatten)]
    pub sort: RelationSortArg,
    #[command(flatten)]
    pub u: UnsortedArg,
    #[command(flatten)]
    pub genome: GenomeArg,
    #[command(flatten)]
    pub chrom_order: ChromOrderArg,
    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

impl FisherArgs {
    pub fn criteria(&self) -> gloci_core::OverlapCriteria {
        OverlapArgs {
            min_frac_a: self.min_frac_a,
            min_frac_b: self.min_frac_b,
            reciprocal: self.reciprocal,
            either: self.either,
            min_bases: self.min_bases,
            same_strand: self.same_strand,
            opposite_strand: self.opposite_strand,
            split: false,
            inclusive_points: self.inclusive_points,
        }
        .criteria()
    }
}

impl JaccardArgs {
    pub fn criteria(&self) -> gloci_core::OverlapCriteria {
        OverlapArgs {
            min_frac_a: self.min_frac_a,
            min_frac_b: self.min_frac_b,
            reciprocal: self.reciprocal,
            either: self.either,
            min_bases: self.min_bases,
            same_strand: self.same_strand,
            opposite_strand: self.opposite_strand,
            split: false,
            inclusive_points: self.inclusive_points,
        }
        .criteria()
    }
}

/// Random intervals placed uniformly over a genome (bedtools random).
///
/// Every valid start of every chromosome is equally likely: a chromosome
/// of L bases has L - l + 1 starts for intervals of length l, and one
/// shorter than l has none. Output is BED6 (chromosome, start, end, id,
/// length, strand) in id order, not sorted (pipe to gloci sort to sort).
/// -l longer than every chromosome is an error naming the longest one.
/// The output depends on --seed only, never on --threads; without --seed
/// a fresh seed is chosen and printed as a warning, so rerunning with it
/// reproduces the output. It is not bedtools' output for the same seed
/// (see docs/differences.md#rng-stream).
#[derive(Args, Debug, Clone)]
#[command(about = "Random intervals over a genome (bedtools random)")]
#[command(after_long_help = "\
Examples:
  gloci random -g hg38 -l 500 -n 100000 --seed 42 -o random.bed
  gloci random -g genome.txt -l 1000 -n 1000 --seed 7 | gloci sort -g genome.txt - > sorted.bed
  gloci random -g hg38 -n 10000 --seed 1 -o a.bed && gloci shift a.bed --by 5 -g hg38")]
pub struct RandomArgs {
    /// Chromosomes and lengths: a built-in genome such as hg38 (see gloci genomes list), or a
    /// .genome, .fai or .dict file.
    #[arg(short, long, value_name = "NAME|FILE", required = true)]
    pub genome: String,
    /// Length of each interval in bases.
    #[arg(short = 'l', long, value_name = "BASES", default_value_t = 100)]
    pub length: u64,
    /// How many intervals to write.
    #[arg(short = 'n', long, value_name = "COUNT", default_value_t = 1_000_000)]
    pub count: u64,
    #[command(flatten)]
    pub seed: SeedArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Move every record to a random place in the genome (bedtools shuffle).
///
/// Each record keeps its length, its other columns and its strand and goes to a uniformly random
/// place that keeps all of it inside the allowed regions: --include regions (default: the whole
/// genome) less --exclude regions. (bedtools only puts the start inside an --include region; see
/// docs/differences.md#shuffle-placement.) --same-chrom keeps records on their chromosome;
/// --chrom-first picks a chromosome at random among those with room, then a place on it, instead of
/// weighting chromosomes by their free length; --no-overlapping redraws a place that overlaps an
/// earlier output record, up to --max-tries times.
#[derive(Args, Debug, Clone)]
#[command(about = "Move every record to a random place (bedtools shuffle)")]
#[command(after_long_help = "\
Examples:
  gloci shuffle peaks.bed -g hg38 --exclude gaps.bed --seed 42 -o null.bed
  gloci shuffle peaks.bed -g genome.txt --include mappable.bed --same-chrom --seed 7 | gloci sort -g genome.txt -
  gloci shuffle loops.bedpe --format bedpe -g hg38 --exclude blacklist.bed --seed 1 -o null.bedpe

Placement:
A record with no place is an error naming its line, length and the fix (nothing is written \
to a file); --skip-unplaceable drops such records and warns with a count. With --format \
bedpe each pair is moved: a pair on one chromosome keeps its end lengths and gap, a pair \
on two places its ends independently, and only --exclude, --seed, --max-tries and \
--skip-unplaceable apply.

Output:
Output is in input order and not sorted (pipe to gloci sort). WIG is read as intervals and \
written as bedGraph; VCF and alignment input are refused (gloci convert). \
--include/--exclude read the names of -g, after --chrom-map (which renames them like the \
input).

Random numbers:
The output depends on --seed only, never on --threads; without --seed a fresh seed is \
chosen and printed as a warning. It is not bedtools' output for the same seed (see \
docs/differences.md#rng-stream).")]
pub struct ShuffleArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Chromosomes and lengths: a built-in genome such as hg38 (see gloci genomes list), or a
    /// .genome, .fai or .dict file.
    #[arg(short, long, value_name = "NAME|FILE", required = true)]
    pub genome: String,
    /// Keep records inside these regions (BED, GFF or VCF; default: the whole genome).
    #[arg(long, value_name = "BED")]
    pub include: Option<String>,
    /// Keep records off these regions (BED, GFF or VCF).
    #[arg(long, value_name = "BED")]
    pub exclude: Option<String>,
    /// Keep each record on its own chromosome (bedtools -chrom).
    #[arg(long, conflicts_with = "chrom_first")]
    pub same_chrom: bool,
    /// Pick the chromosome uniformly among those with room, then a place on it
    /// (bedtools -chromFirst).
    #[arg(long)]
    pub chrom_first: bool,
    /// Place records so that no two output records overlap (bedtools
    /// -noOverlapping).
    #[arg(long)]
    pub no_overlapping: bool,
    /// Draws per record before it counts as unplaceable (--no-overlapping, and
    /// BEDPE --exclude).
    #[arg(long, value_name = "N", default_value_t = 1000,
          value_parser = clap::value_parser!(u32).range(1..))]
    pub max_tries: u32,
    /// Drop records that have no place, with a warning, instead of failing.
    #[arg(long)]
    pub skip_unplaceable: bool,
    #[command(flatten)]
    pub seed: SeedArg,
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

/// A random subset of the records of a file (bedtools sample; for
/// alignments see samtools view -s, gloci samples is a different command).
///
/// Keeps exactly N records, chosen uniformly at random, and writes them
/// in input order (bedtools writes them in reservoir order; see
/// docs/differences.md#sample-order). With --strand + or -, only records
/// on that strand are sampled (a record of unknown strand never is).
/// Fewer than N records is an error that gives the count and writes
/// nothing. BED, bedGraph, GFF/GTF and VCF/BCF are read; --header keeps
/// the header lines of text formats and a VCF always keeps its header.
/// The kept records are held in memory. The output depends on --seed
/// only, never on --threads; without --seed a fresh seed is chosen and
/// printed as a warning. It is not bedtools' output for the same seed
/// (see docs/differences.md#rng-stream).
#[derive(Args, Debug, Clone)]
#[command(about = "A random subset of the records (bedtools sample)")]
#[command(after_long_help = "\
Examples:
  gloci sample peaks.bed -n 1000 --seed 42 -o subset.bed
  gloci sample peaks.bed -n 500 --strand + --seed 1
  gloci sample calls.vcf.gz -n 10000 --seed 3 -O vcf.gz -o subset.vcf.gz")]
pub struct SampleArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// How many records to keep.
    #[arg(short = 'n', long, value_name = "COUNT", required = true)]
    pub count: u64,
    /// Sample only the records on this strand (+ or -).
    #[arg(long, value_name = "+|-", allow_hyphen_values = true,
          value_parser = clap::builder::PossibleValuesParser::new(["+", "-"]))]
    pub strand: Option<String>,
    #[command(flatten)]
    pub seed: SeedArg,
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
