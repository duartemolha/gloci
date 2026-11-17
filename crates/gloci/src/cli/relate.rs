use super::shared::{
    AArg, AlignArgs, ChromMapArg, ChromOrderArg, GenomeArg, HeaderArg, InvalidArg, OutputArgs,
    OverlapArgs, RegionArgs, RelationSortArg, ThreadArgs, UnsortedArg, VcfSpanArg, WindowArgs,
};
use clap::{Args, ValueEnum};

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

/// Report the records of A that overlap, or do not overlap, records of B (bedtools intersect).
///
/// The default mode is any: A's records that overlap at least one B record are written as read, in
/// A's format. --mode picks another report: none (A's records with no overlap), pairs (one TSV row
/// per overlapping A and B pair; --overlap-bases adds the shared bases and --keep-unmatched adds a
/// row with placeholder B columns for each A without one), clip (A cut to the shared span, once per
/// overlap) and count (A with its number of overlaps, summed over all -b files). With several -b
/// files, pairs rows carry a label column (--names, else the file basename) before B's columns;
/// unmatched rows get `.`.
#[derive(Args)]
#[command(about = "Report A records that overlap B (bedtools intersect)")]
#[command(after_long_help = "\
Examples:
  gloci intersect -a peaks.bed -b genes.bed
  gloci intersect -a calls.vcf.gz -b exons.bed --mode none -o off_target.vcf.gz
  gloci intersect -a a.bed -b x.bed y.bed --mode pairs --overlap-bases --names X,Y
  gloci intersect -a reads.bam -b targets.bed -o on_target.bam

Overlap rule:
Records overlap when they share a base, a point overlaps an interval when strictly inside \
it (or on its edge with --inclusive-points), and every threshold (-f, -F, -r, -e, \
--min-bases) is met when either side is a point. Thresholds are measured against B's own \
span; they cannot be combined with -w, --left or --right. A VCF/BCF record covers the \
bases its ALT alleles change (--vcf-span). --match alleles also requires some ALT of A to \
equal some ALT of B after trimming (VCF/BCF inputs only).

Sorted input:
A and every B must be sorted in one chromosome order, which is checked as they are read; a \
violation is an error with a fix. With --unsorted, any order works, A is still streamed \
and reported in its own order, and every -b file is loaded into memory, so memory grows \
with the size of B.

Two ways to deal with inputs that are not sorted. Use --sort when the inputs are too big \
for memory or you want sorted output: A and every -b file are sorted first with an \
external sort bounded by --max-mem (spilling to --tmp-dir), then the normal sweep runs, so \
the result equals gloci sort on each input followed by this command. It needs disk, about \
twice the input size in --tmp-dir, and removes it afterwards. Use --unsorted when you want \
A's order kept and every -b file fits in memory. The two cannot be combined.

Output:
Pairs and count write TSV (a VCF A gets a #CHROM line naming the columns); clip does not \
accept VCF/BCF.

Alignments:
BAM, CRAM and SAM (A or B) work as intervals (a CRAM needs --reference, and an indexed BAM \
or CRAM is read per chromosome at -t above 1): a read is the reference span of its M, D, \
N, = and X operations (clips never widen it), on the strand of its flag. Unmapped reads \
are skipped; --min-mapq, --require-flags and --exclude-flags filter reads first. --split \
uses the segments between N operations, --split-deletions also cuts at D. With an \
alignment A, --mode any and none write the reads themselves as BAM (SAM for a SAM A, or as \
-o/-O say); pairs and count write A's reads as BED6 rows, as gloci convert does. --mode \
clip of reads is refused. Alignments stream: --sort sorts them into a temporary BAM first, \
bounded by --max-mem.

Parallel runs:
With more than one thread, when A and every -b file are BGZF files with a .tbi or .csi \
index (a BAM with a .bai or .csi) or plain (uncompressed) BED, bedGraph or VCF files, the \
chromosomes are worked on in parallel and written in the chromosome-table order (lex \
by default; -g, a VCF header or --chrom-order set it). The output is the same as at -t 1.

A parallel run holds about 64 MiB plus 20 MiB for each run of chromosomes in flight, at \
most 8 (more with many -b files; threads beyond 8 still speed up BGZF). Output waiting to \
be written beyond 64 MiB spills to --tmp-dir.

What depends on the thread count: after an error, how much output was already written to \
standard output (-o FILE is written atomically); which input the error names when more \
than one is malformed (the exit code is the same); and, with -v, the counts of reads left \
out of an indexed -b or of anything under --region, which are of the reads read. Plain \
inputs are validated to their end at every thread count; for indexed inputs, records \
outside the queried chromosomes are validated only at -t 1. An input that is not in \
chromosome-table order is refused with the same error at every thread count.")]
#[command(mut_arg("region", |a| a.help(RELATION_REGION_HELP)))]
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
    /// BAM, CRAM and SAM: also cut the blocks of a read at deletions (D), with or without
    /// --split (bedtools bamtobed -splitD).
    #[arg(long)]
    pub split_deletions: bool,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub region: RegionArgs,
    #[command(flatten)]
    pub u: UnsortedArg,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Remove from A the bases covered by B (bedtools subtract, bedops -d).
///
/// Each A record is cut by the union of the B records that overlap it, and what is left is written as
/// A's own record with new coordinates: one record per piece, in ascending order (BED columns 2-3,
/// GFF columns 4-5). An A that nothing touches is written as read, and an A that is covered
/// completely disappears. A point has no bases to cut: it is dropped when a B contains it, and a B
/// point removes nothing from a cut; under --whole, or --sum-overlaps without -f, a B point strictly
/// inside A (or on its edge with --inclusive-points) still removes A (the point rule of intersect).
/// The criteria (-f, -F, -r, -e, --min-bases, -s, -S) choose which B count and are measured as in
/// intersect. -w, --left and --right pad B before it cuts.
#[derive(Args)]
#[command(about = "Remove the bases of B from A (bedtools subtract, bedops -d)")]
#[command(after_long_help = "\
Examples:
  gloci subtract -a peaks.bed -b blacklist.bed
  gloci subtract -a reads.bed -b exons.bed --sum-overlaps -f 0.8 -s
  gloci subtract -a calls.vcf.gz -b blacklist.bed -o filtered.vcf.gz
  gloci subtract -a calls.vcf.gz -b known.vcf.gz --match alleles -o novel.vcf.gz

Other modes:
--whole drops an A that has any counting B instead of cutting it (bedtools -A). \
--sum-overlaps drops an A when the bases shared with all the B that pass the strand rule \
add up to at least -f of A (bedtools -N); it takes only -f, and without -f it acts like \
--whole. --split cuts by the BED12 blocks of B, not its span. A VCF/BCF A is never cut: a \
record is dropped whole, in default or --whole mode when intersect --mode none would drop \
it (any accepted hit), or with --sum-overlaps -f when the summed overlaps reach -f of its \
footprint, and the output is VCF with its header.

Sorted input:
A and every B must be sorted in one chromosome order, which is checked as they are read. \
With --unsorted, any order works, A is still streamed and written in its own order, and \
every -b file is loaded into memory, so memory grows with the size of B.

Two ways to deal with inputs that are not sorted. Use --sort when the inputs are too big \
for memory or you want sorted output: A and every -b file are sorted first with an \
external sort bounded by --max-mem (spilling to --tmp-dir), then the normal sweep runs, so \
the result equals gloci sort on each input followed by this command. It needs disk, about \
twice the input size in --tmp-dir, and removes it afterwards. Use --unsorted when you want \
A's order kept and every -b file fits in memory. The two cannot be combined.

Parallel runs:
With more than one thread, when A and every -b file are BGZF files with a .tbi or .csi \
index (a BAM with a .bai or .csi) or plain (uncompressed) BED, bedGraph or VCF files, the \
chromosomes are worked on in parallel and written in the chromosome-table order (lex \
by default; -g, a VCF header or --chrom-order set it). The output is the same as at -t 1.

A parallel run holds about 64 MiB plus 20 MiB for each run of chromosomes in flight, at \
most 8 (more with many -b files; threads beyond 8 still speed up BGZF). Output waiting to \
be written beyond 64 MiB spills to --tmp-dir.

What depends on the thread count: after an error, how much output was already written to \
standard output (-o FILE is written atomically); which input the error names when more \
than one is malformed (the exit code is the same); and, with -v, the counts of reads left \
out of an indexed -b or of anything under --region, which are of the reads read. Plain \
inputs are validated to their end at every thread count; for indexed inputs, records \
outside the queried chromosomes are validated only at -t 1. An input that is not in \
chromosome-table order is refused with the same error at every thread count.")]
#[command(mut_arg("region", |a| a.help(RELATION_REGION_HELP)))]
pub struct SubtractArgs {
    #[command(flatten)]
    pub a: AArg,
    #[command(flatten)]
    pub b: BArgs,
    /// Drop an A that has any counting B instead of cutting it.
    #[arg(long)]
    pub whole: bool,
    /// Drop an A when the bases shared with all B (strand rule only) reach
    /// -f of A; takes only -f.
    #[arg(long)]
    pub sum_overlaps: bool,
    /// What makes a hit count: overlap, or overlap plus a shared ALT allele
    /// (alleles; VCF/BCF inputs only). With alleles, an A is dropped only
    /// for a B with the same position, REF and ALT; one shared ALT drops a
    /// multi-allelic A whole.
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
    pub region: RegionArgs,
    #[command(flatten)]
    pub u: UnsortedArg,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// The frame of signed distances (`--signed`).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum SignedFrame {
    /// Upstream is lower coordinates.
    Ref,
    /// Upstream follows A's strand.
    A,
    /// Upstream follows B's strand.
    B,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum ClosestTies {
    All,
    First,
    Last,
}

/// Find the nearest B to each A (bedtools closest, bedops closest-features).
///
/// For each A, reports the nearest B as TSV rows: A's record, a label column when there are several
/// -b files, B's record, and with --distance the distance. An overlapping B has distance 0, bookended
/// features 1, and otherwise the gap in bases plus 1. An A with no candidate gets one row with a
/// placeholder B (`. -1 -1 ...`) and distance -1. A VCF/BCF A gets a column-name line.
#[derive(Args)]
#[command(about = "Nearest B per A (bedtools closest, bedops closest-features)")]
#[command(after_long_help = "\
Examples:
  gloci closest -a genes.bed -b peaks.bed --distance
  gloci closest -a genes.bed -b peaks.bed -k 3 --signed a -s
  gloci closest -a snps.vcf.gz -b genes.bed --max-distance 10000 --ignore-overlaps
  gloci closest -a genes.bed -b peaks.bed --each-side --distance

Ranks and ties:
-k N reports the N nearest, and every B tied with the N-th. With --ties first or last, \
only as many of the tied B as N still needs are kept: the first or the last by (file, \
input order). Rows are ordered by distance, then by file and input order.

Distance and strand:
--signed ref|a|b prints signed distances (and implies --distance): negative is upstream, \
where ref means lower coordinates, a follows A's strand and b follows B's strand; unknown \
strand counts as +. --ignore-upstream and --ignore-downstream drop B by that sign in the \
same frame (ref when --signed is absent), --ignore-overlaps drops overlapping B, and \
--max-distance N drops B farther than N. -s and -S keep only B on the same or the opposite \
strand; B and A of unknown strand never match.

Each side:
--each-side writes one row per A with the nearest upstream B and then the nearest \
downstream B; with --distance each distance is signed (as closest-features does). An \
overlapping B counts on both sides. -k applies per side, and extra ranks and ties add \
rows. With several -b files, the nearest are chosen across all of them, or with --per-file \
in each file, whose rows are written in file order; a file with no candidate for an A \
writes no row (as bedtools -mdb each), and an A with no candidate in any file gets one \
placeholder row. --each-side with --ignore-upstream or --ignore-downstream leaves that \
side empty except for an overlapping B, which still shows there at distance 0 unless \
--ignore-overlaps is given.

Memory and a pathology:
The sorted sweep looks ahead in B for the nearest B to the right of A. Under -s/-S, or \
with --signed b and --ignore-upstream, --ignore-downstream or --each-side, a chromosome \
where one strand has no B left makes the look-ahead read and buffer all the remaining B of \
that chromosome. --max-distance bounds this. With --unsorted every -b file is loaded into \
memory instead, so memory grows with the size of B, and A is streamed and written in its \
own order.

Sorted input:
Two ways to deal with inputs that are not sorted. Use --sort when the inputs are too big \
for memory or you want sorted output: A and every -b file are sorted first with an \
external sort bounded by --max-mem (spilling to --tmp-dir), then the normal sweep runs, so \
the result equals gloci sort on each input followed by this command. It needs disk, about \
twice the input size in --tmp-dir, and removes it afterwards. Use --unsorted when you want \
A's order kept and every -b file fits in memory. The two cannot be combined.

Parallel runs:
With more than one thread, when A and every -b file are BGZF files with a .tbi or .csi \
index (a BAM with a .bai or .csi) or plain (uncompressed) BED, bedGraph or VCF files, the \
chromosomes are worked on in parallel and written in the chromosome-table order (lex \
by default; -g, a VCF header or --chrom-order set it). The output is the same as at -t 1.

A parallel run holds about 64 MiB plus 20 MiB for each run of chromosomes in flight, at \
most 8 (more with many -b files; threads beyond 8 still speed up BGZF). Output waiting to \
be written beyond 64 MiB spills to --tmp-dir.

What depends on the thread count: after an error, how much output was already written to \
standard output (-o FILE is written atomically); which input the error names when more \
than one is malformed (the exit code is the same); and, with -v, the counts of reads left \
out of an indexed -b or of anything under --region, which are of the reads read. Plain \
inputs are validated to their end at every thread count; for indexed inputs, records \
outside the queried chromosomes are validated only at -t 1. An input that is not in \
chromosome-table order is refused with the same error at every thread count.")]
#[command(mut_arg("region", |a| a.help(RELATION_REGION_HELP)))]
pub struct ClosestArgs {
    #[command(flatten)]
    pub a: AArg,
    #[command(flatten)]
    pub b: BArgs,
    /// Report the N nearest B (and their ties).
    #[arg(
        short = 'k',
        long,
        value_name = "N",
        default_value_t = 1,
        value_parser = clap::value_parser!(u64).range(1..)
    )]
    pub k: u64,
    /// Among B tied at the N-th distance: keep all, the first or the last.
    #[arg(long, value_enum, default_value = "all", value_name = "WHICH")]
    pub ties: ClosestTies,
    /// Add the distance as the last column (-1 with no candidate).
    #[arg(long)]
    pub distance: bool,
    /// Print signed distances in this frame (implies --distance): ref, a or b.
    #[arg(long, value_enum, value_name = "FRAME")]
    pub signed: Option<SignedFrame>,
    /// Ignore B that overlap A.
    #[arg(long)]
    pub ignore_overlaps: bool,
    /// Ignore B upstream of A (in the --signed frame, default ref).
    #[arg(long)]
    pub ignore_upstream: bool,
    /// Ignore B downstream of A (in the --signed frame, default ref).
    #[arg(long)]
    pub ignore_downstream: bool,
    /// Ignore B farther than N from A.
    #[arg(long, value_name = "N")]
    pub max_distance: Option<u64>,
    /// One row per A: the nearest upstream B, then the nearest downstream B.
    #[arg(long)]
    pub each_side: bool,
    /// With several -b files: choose the nearest in each file.
    #[arg(long)]
    pub per_file: bool,
    /// Only B on the same strand; unknown strands never match.
    #[arg(short = 's', long = "same-strand", conflicts_with = "opposite_strand")]
    pub same_strand: bool,
    /// Only B on the opposite strand; unknown strands never match.
    #[arg(short = 'S', long = "opposite-strand")]
    pub opposite_strand: bool,
    /// A point on the edge of an interval overlaps it (4.3: points).
    #[arg(long)]
    pub inclusive_points: bool,
    #[command(flatten)]
    pub region: RegionArgs,
    #[command(flatten)]
    pub u: UnsortedArg,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// How much of each A the B records cover (bedtools coverage and annotate).
///
/// For each A, in A's order, counts the B records that overlap it (the criteria of intersect: -f, -F,
/// -r, -e, --min-bases, -s, -S, --split, --inclusive-points) and writes A's record followed by the
/// number of hits, the bases of A covered by at least one hit, A's length and the covered fraction.
/// --hist writes one row per depth present in A (depth, bases at that depth, A's length, fraction)
/// and, after the last A, `all` rows with the totals. --per-base writes one row per base of A
/// (1-based offset, depth). --mean writes the mean depth over A. With several -b files only the
/// default report is allowed and each file gives one covered fraction (bedtools annotate).
#[derive(Args)]
#[command(about = "How much of each A the B records cover (bedtools coverage)")]
#[command(after_long_help = "\
Examples:
  gloci coverage -a genes.bed -b reads.bed
  gloci coverage -a exons.bed -b reads.bed --per-base -s
  gloci coverage -a peaks.bed -b rep1.bed rep2.bed --names r1,r2
  gloci coverage -a targets.bed -b reads.bam --split

Output format:
Fractions and the mean depth are written as bedtools writes them: seven decimals of a \
32-bit quotient, and six decimals of a 64-bit one in the report of several -b files \
(bedtools annotate).

The output is TSV with A's columns first.

Blocks and points:
With --split, A and B are their BED12 blocks: depth counts the bases shared by an A block \
and a B block, A's length is the sum of its blocks, and --per-base offsets run along the \
blocks. A point A covers no bases: it reports its hit count, 0, 0 and `.`, --mean gives \
`.`, and --hist and --per-base write no rows for it.

Sorted input:
A and every B must be sorted in one chromosome order, which is checked as they are read. \
With --unsorted, any order works, A is still streamed in its own order, and every -b file \
is loaded into memory, so memory grows with the size of B.

Two ways to deal with inputs that are not sorted. Use --sort when the inputs are too big \
for memory or you want sorted output: A and every -b file are sorted first with an \
external sort bounded by --max-mem (spilling to --tmp-dir), then the normal sweep runs, so \
the result equals gloci sort on each input followed by this command. It needs disk, about \
twice the input size in --tmp-dir, and removes it afterwards. Use --unsorted when you want \
A's order kept and every -b file fits in memory. The two cannot be combined.

Alignments:
BAM, CRAM and SAM (A or B) work as intervals (a CRAM needs --reference, and an indexed BAM \
or CRAM is read per chromosome at -t above 1): a read is the reference span of its M, D, \
N, = and X operations (clips never widen it). Unmapped reads are skipped; --min-mapq, \
--require-flags and --exclude-flags filter reads first. With --split the blocks of a read \
are the segments between N operations; --split-deletions also cuts at D. An alignment A is \
written as its BED6 row (as gloci convert does). Alignments stream: --sort sorts them into \
a temporary BAM first, bounded by --max-mem.

Parallel runs:
With more than one thread, when A and every -b file are BGZF files with a .tbi or .csi \
index (a BAM with a .bai or .csi) or plain (uncompressed) BED, bedGraph or VCF files, the \
chromosomes are worked on in parallel and written in the chromosome-table order (lex \
by default; -g, a VCF header or --chrom-order set it). The output is the same as at -t 1.

A parallel run holds about 64 MiB plus 20 MiB for each run of chromosomes in flight, at \
most 8 (more with many -b files; threads beyond 8 still speed up BGZF). Output waiting to \
be written beyond 64 MiB spills to --tmp-dir.

What depends on the thread count: after an error, how much output was already written to \
standard output (-o FILE is written atomically); which input the error names when more \
than one is malformed (the exit code is the same); and, with -v, the counts of reads left \
out of an indexed -b or of anything under --region, which are of the reads read. Plain \
inputs are validated to their end at every thread count; for indexed inputs, records \
outside the queried chromosomes are validated only at -t 1. An input that is not in \
chromosome-table order is refused with the same error at every thread count.")]
#[command(mut_arg("region", |a| a.help(RELATION_REGION_HELP)))]
pub struct CoverageArgs {
    #[command(flatten)]
    pub a: AArg,
    #[command(flatten)]
    pub b: BArgs,
    /// One row per depth present in A, then `all` rows with the totals.
    #[arg(long, group = "report")]
    pub hist: bool,
    /// One row per base of A: its 1-based offset and depth.
    #[arg(long, group = "report")]
    pub per_base: bool,
    /// The mean depth over A.
    #[arg(long, group = "report")]
    pub mean: bool,
    #[command(flatten)]
    pub o: OverlapArgs,
    /// BAM, CRAM and SAM: also cut the blocks of a read at deletions (D), with or without
    /// --split (bedtools bamtobed -splitD).
    #[arg(long)]
    pub split_deletions: bool,
    #[command(flatten)]
    pub align: AlignArgs,
    #[command(flatten)]
    pub region: RegionArgs,
    #[command(flatten)]
    pub u: UnsortedArg,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// Aggregate columns of the B records over each A (bedtools map, bedops bedmap).
///
/// For each A, in A's order, folds the B records that overlap it (the criteria of intersect: -f, -F,
/// -r, -e, --min-bases, -s, -S, --split, or a window -w/--left/--right) and writes A's record
/// followed by one column per -c COL:OP, in B file then B input order. COL is a 1-based column of B,
/// 0 for the whole record (count, bases, indicator, bases_uniq, frac_a, frac_b and records only), or
/// for a VCF B QUAL, ID, FILTER or INFO/KEY.
#[derive(Args)]
#[command(about = "Aggregate B columns over each A (bedtools map, bedops bedmap)")]
#[command(after_long_help = "\
Examples:
  gloci map -a genes.bed -b scores.bedgraph -c 4:mean
  gloci map -a windows.bed -b peaks.bed -c 0:count -c 5:max --skip-unmatched
  gloci map -a snps.vcf.gz -b regions.bed -c 4:distinct --delim ';'
  gloci map -a genes.bed -b peaks.bed -w 100 -c 0:count -c 0:bases

Operations:
Ops: count count_distinct sum mean median min max stdev var stdev_pop var_pop collapse \
distinct distinct_sorted first last bases indicator cv mad quantile:Q tmean:LO-HI wmean \
wsum mode antimode bases_uniq frac_a frac_b min_record max_record records (COL 0 also works \
with bases_uniq, frac_a, frac_b and records; records ignores COL). cv is the sample stdev \
over the mean; mad the unscaled median absolute deviation; quantile is type 7 and tmean \
the mean of the values whose rank fraction is within LO-HI (a lone value is kept); wmean \
weights by overlapping bases; wsum is the sum of value times overlapping bases (over a bigWig, \
-c 4:wsum -c 0:bases_uniq -c 4:wmean --empty 0 gives bigWigAverageOverBed's sum, covered \
and mean); ties of mode and antimode go to the bytewise smallest value, \
as in bedtools; bases_uniq counts overlapped bases once; frac_a is bases_uniq over A's \
length and frac_b the summed overlap over the summed B length (--empty when the length is \
0); records copies whole B records, tabs shown as |, and does not escape --delim where it \
occurs in the data. count counts every overlapping B; the other ops skip `.` and empty \
values, and a value that is not a number in a numeric op is an error. sum, min and max \
print integers when every value was one; stdev and var are sample statistics (one value \
gives --empty); stdev_pop and var_pop divide by n (one value gives 0), as bedtools stdev \
does. bases is the overlap of A with each B's own span, which is 0 for a B that is within \
a window only; indicator is 1 when any B counts, else 0. An A with no B gets 0 for count \
and count_distinct and --empty (default `.`) for the rest; --skip-unmatched drops such A \
instead.

Columns and output:
Columns are resolved against the first -b file, and all -b files must be one format \
family. The output is TSV with A's columns first; a VCF/BCF A gets a column-name line \
ending with the -c labels.

Sorted input:
A and every B must be sorted in one chromosome order, which is checked as they are read. \
With --unsorted, any order works, A is still streamed in its own order, and every -b file \
is loaded into memory, so memory grows with the size of B.

Two ways to deal with inputs that are not sorted. Use --sort when the inputs are too big \
for memory or you want sorted output: A and every -b file are sorted first with an \
external sort bounded by --max-mem (spilling to --tmp-dir), then the normal sweep runs, so \
the result equals gloci sort on each input followed by this command. It needs disk, about \
twice the input size in --tmp-dir, and removes it afterwards. Use --unsorted when you want \
A's order kept and every -b file fits in memory. The two cannot be combined.

Parallel runs:
With more than one thread, when A and every -b file are BGZF files with a .tbi or .csi \
index (a BAM with a .bai or .csi) or plain (uncompressed) BED, bedGraph or VCF files, the \
chromosomes are worked on in parallel and written in the chromosome-table order (lex \
by default; -g, a VCF header or --chrom-order set it). The output is the same as at -t 1.

A parallel run holds about 64 MiB plus 20 MiB for each run of chromosomes in flight, at \
most 8 (more with many -b files; threads beyond 8 still speed up BGZF). Output waiting to \
be written beyond 64 MiB spills to --tmp-dir.

What depends on the thread count: after an error, how much output was already written to \
standard output (-o FILE is written atomically); which input the error names when more \
than one is malformed (the exit code is the same); and, with -v, the counts of reads left \
out of an indexed -b or of anything under --region, which are of the reads read. Plain \
inputs are validated to their end at every thread count; for indexed inputs, records \
outside the queried chromosomes are validated only at -t 1. An input that is not in \
chromosome-table order is refused with the same error at every thread count.")]
#[command(mut_arg("region", |a| a.help(RELATION_REGION_HELP)))]
pub struct MapArgs {
    #[command(flatten)]
    pub a: AArg,
    #[command(flatten)]
    pub b: BArgs,
    /// Aggregate a column of the overlapping B, such as 5:mean or 0:count
    /// (repeatable; commas also separate).
    #[arg(short = 'c', long = "column", value_name = "COL:OP", required = true)]
    pub columns: Vec<String>,
    /// Drop an A with no overlapping B instead of writing it.
    #[arg(long)]
    pub skip_unmatched: bool,
    /// Value for an aggregation with no values.
    #[arg(long, value_name = "S", default_value = ".")]
    pub empty: String,
    /// Separator for collapse and distinct lists.
    #[arg(long, value_name = "S", default_value = ",")]
    pub delim: String,
    #[command(flatten)]
    pub o: OverlapArgs,
    #[command(flatten)]
    pub w: WindowArgs,
    #[command(flatten)]
    pub region: RegionArgs,
    #[command(flatten)]
    pub u: UnsortedArg,
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
    pub header: HeaderArg,
    #[command(flatten)]
    pub invalid: InvalidArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

/// `--region` help in a relation command, where only A is restricted
/// (applied with `mut_arg` on each relation command's args).
const RELATION_REGION_HELP: &str = "Restrict A to a region such as chr1:1,000-2,000 (repeatable). Only A is restricted, to the records that overlap a region (a point on the edge with --inclusive-points); B is read in full, so an A record reaching out of a region still meets B there.";

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
    /// basenames).
    #[arg(long, value_name = "A,B,...", value_delimiter = ',')]
    pub names: Option<Vec<String>>,
    /// B's format, overriding detection: bed, bedgraph, wig, gff, gtf, vcf or
    /// bcf; bigwig or bigbed (from a file); bam, cram or sam where the
    /// command reads alignments.
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
                    // An in-memory input is labelled by the basename of its display name.
                    let shown = gloci_io::memio::display_name(p);
                    std::path::Path::new(&shown)
                        .file_name()
                        .map_or_else(|| shown.clone(), |f| f.to_string_lossy().into_owned())
                })
                .collect()),
        }
    }
}
