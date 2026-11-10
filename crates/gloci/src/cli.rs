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
    /// header lines are kept by default. For conversions to BED from VCF/BCF,
    /// GFF3/GTF or BED12, use gloci convert.
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
    ///
    /// Two ways to deal with inputs that are not sorted. Use --sort when the
    /// inputs are too big for memory or you want sorted output: A and every
    /// -b file are sorted first with an external sort bounded by --max-mem
    /// (spilling to --tmp-dir), then the normal sweep runs, so the result
    /// equals gloci sort on each input followed by this command. It needs
    /// disk, about twice the input size in --tmp-dir, and removes it
    /// afterwards. Use --unsorted when you want A's order kept and every
    /// -b file fits in memory. The two cannot be combined.
    ///
    /// When A and every -b file are BGZF files with a .tbi or .csi index and
    /// more than one thread is used, the work runs on runs of chromosomes in
    /// parallel, which are read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF); their
    /// output waiting to be written shares 64 MiB of memory and the rest
    /// spills to --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count.
    #[command(after_long_help = "\
Examples:
  gloci intersect -a peaks.bed -b genes.bed
  gloci intersect -a calls.vcf.gz -b exons.bed --mode none -o off_target.vcf.gz
  gloci intersect -a a.bed -b x.bed y.bed --mode pairs --overlap-bases --names X,Y
  gloci intersect -a reads.bed -b targets.bed -f 50% -s --mode count --unsorted")]
    Intersect(IntersectArgs),
    /// Merge overlapping or nearby records into clusters (bedtools merge, bedops -m).
    ///
    /// Records are grouped by their footprints: a record joins the open
    /// cluster when its start is at most the cluster end plus -d/--max-gap.
    /// The default 0 joins overlapping and bookended features; a negative
    /// gap such as -d -2 needs at least that many shared bases; a point
    /// joins a cluster it lies inside or touches. The output is BED
    /// (chrom, start, end), or TSV for -o x.tsv, and is not A's records.
    /// --count adds the number of records per cluster and -c COL:OP adds
    /// one aggregation column per option (count, sum, mean, median, min,
    /// max, stdev, var, stdev_pop, var_pop, collapse, distinct, distinct_sorted, first, last,
    /// count_distinct, indicator, cv, mad, quantile:Q, tmean:LO-HI, mode,
    /// antimode; the ops that need a second input, like wmean, are for
    /// map). -s clusters each strand separately
    /// (unknown strand apart from + and -) without adding a strand column.
    ///
    /// Several files are merged together as one stream. Every input must be
    /// sorted in one chromosome order, which is checked as it is read; an
    /// unsorted input is an error whose hint is to run gloci sort or to pass
    /// --sort, which sorts the inputs first with an external sort bounded by
    /// --max-mem (spilling to --tmp-dir); the result equals gloci sort
    /// FILES | gloci merge -. Memory stays flat apart from the clusters
    /// being built.
    #[command(after_long_help = "\
Examples:
  gloci merge peaks.bed
  gloci merge a.bed b.bed -d 100 --count
  gloci merge genes.bed -s -c 5:max -c 4:distinct
  gloci merge calls.vcf.gz -d -1 -o regions.bed.gz --write-index")]
    Merge(MergeArgs),
    /// Write the bases no feature covers (bedtools complement, bedops -c).
    ///
    /// Features of all files are united, and the gaps between them are
    /// written as BED (chrom, start, end). With -g, the gaps are taken over
    /// each chromosome of the genome, in genome order: a chromosome with no
    /// features becomes one row over its whole length, and the gaps before
    /// the first and after the last feature are included. The genome needs
    /// lengths (a built-in genome, .genome, .fai or .dict file); a feature
    /// that reaches past the end of its chromosome is an error. Without -g
    /// only the gaps between features of the same chromosome are written.
    /// Points cover no bases and are ignored. Bookended features leave no
    /// gap.
    ///
    /// Every input must be sorted in one chromosome order, which is checked
    /// as it is read; with -g, use the genome's own order (gloci sort -g).
    /// An unsorted input is an error whose hint names --sort, which sorts
    /// the inputs first with an external sort bounded by --max-mem; the
    /// result equals gloci sort FILES | gloci complement -.
    #[command(after_long_help = "\
Examples:
  gloci complement peaks.bed
  gloci complement peaks.bed -g hg38 -o uncovered.bed
  gloci complement a.bed b.bed -g sacCer3
  gloci complement shuffled.bed -g genome.txt --sort --max-mem 512M")]
    Complement(ComplementArgs),
    /// Split the inputs at every boundary into segments and say who covers
    /// each (bedtools multiinter and unionbedg, bedops -p, -i and -s).
    ///
    /// The footprints of all files are cut at every start and end, giving
    /// disjoint segments in sorted order. Each row is chrom, start, end, the
    /// number of files covering the segment and their labels joined with
    /// commas (--names, default the file names without directories). A file
    /// counts once however many of its features overlap there, and points
    /// cover no bases, so they add neither segments nor boundaries. Only
    /// segments covered by at least one file are written. --min-files K and
    /// --max-files K keep segments covered by at least or at most K files;
    /// --all keeps those covered by every file. --empty implies --min-files 0
    /// (an explicit --min-files overrides it); --min-files 0 needs --empty.
    /// By default a boundary is
    /// kept even where the membership does not change (as in bedops -p);
    /// --merge-adjacent joins neighbouring segments with the same output
    /// columns (as multiinter does). --no-membership writes coordinates
    /// only, so with --merge-adjacent every run of touching kept segments is
    /// joined.
    ///
    /// With --bedgraph the inputs are bedGraph files and each row has one
    /// value per file instead of the membership; a file that does not cover
    /// the segment shows --empty (default 0). Giving --empty also writes the
    /// segments no file covers, with that value in every column; with -g
    /// (needs lengths) it covers whole chromosomes, including those with no
    /// features, otherwise only the gaps between features. The intervals of
    /// one bedGraph must not overlap.
    ///
    /// Every input must be sorted in one chromosome order, which is checked
    /// as it is read; an unsorted input is an error whose hint names --sort,
    /// which sorts each file first with an external sort bounded by
    /// --max-mem. Memory is bounded by the intervals open at one position.
    #[command(after_long_help = "\
Examples:
  gloci partition a.bed b.bed c.bed --merge-adjacent
  gloci partition a.bed b.bed --names tumour,normal --min-files 2
  gloci partition a.bed b.bed --all --merge-adjacent --no-membership
  gloci partition a.bg b.bg --bedgraph --empty 0 -g hg38")]
    Partition(PartitionArgs),
    /// Depth of coverage of a genome by the features of a BED file
    /// (bedtools genomecov -i).
    ///
    /// Every feature adds one to the depth over its bases. -g is required:
    /// the output covers every chromosome of the genome from base 0 to its
    /// length, in genome order, chromosomes without features included
    /// (a built-in genome, .genome, .fai or .dict file with lengths). The
    /// default --hist writes, per chromosome and then for the whole genome
    /// (chrom `genome`), one row per depth present: chrom, depth, bases at
    /// that depth, chromosome length, fraction (six significant digits, as
    /// bedtools). --bedgraph writes runs of equal depth above zero, --bga
    /// also the zero-depth runs, --per-base one row per base: chrom,
    /// 1-based position, depth. --scale X multiplies the depths of the
    /// last three (six significant digits, as bedtools); --max-depth N
    /// puts every depth of N or more into the bin N of --hist.
    /// --strand + or - counts only features on that strand (BED needs
    /// column 6; a feature with unknown strand is not counted). --split
    /// counts the BED12 blocks of each feature. --5p and --3p count one
    /// base per feature, its 5' or 3' end by strand (unknown counts as
    /// plus). A zero-length feature covers no bases and is ignored. A
    /// feature past the end of its chromosome is an error. BAM input is
    /// not supported yet.
    ///
    /// The input must be sorted in the genome order, which is checked as it
    /// is read; an unsorted input is an error whose hint names --sort.
    /// Memory is bounded by the features open at one position, whatever
    /// the output mode.
    #[command(after_long_help = "\
Examples:
  gloci genomecov peaks.bed -g hg38
  gloci genomecov reads.bed -g hg38 --bedgraph --scale 0.000001 -o cov.bedgraph
  gloci genomecov genes.bed -g genome.txt --bga --strand +
  gloci genomecov exons12.bed -g hg38 --split --per-base")]
    Genomecov(GenomecovArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci slop genes.bed -b 1000 -g hg38
  gloci slop genes.bed --left 2000 --right 500 -s -g hg38 -o promoters.bed
  gloci slop peaks.bed --pct -b 0.5 -g genome.txt
  gloci slop annotation.gff3 -b 100 -g hg38")]
    Slop(SlopArgs),
    /// Write the strips beside each feature (bedtools flank).
    ///
    /// For every record, a left strip [start - L, start) and then a right
    /// strip [end, end + R) are written as records with the input's other
    /// columns, clipped at 0 and, with -g, at the chromosome end. A strip
    /// clipped to nothing is not written. Options as for slop: -b, --left,
    /// --right, -s, --pct, -g. The output is in input order and may not be
    /// sorted; a warning is printed once when it goes out of order.
    #[command(after_long_help = "\
Examples:
  gloci flank genes.bed -b 100 -g hg38
  gloci flank genes.bed --left 2000 -s -g hg38
  gloci flank peaks.bed --pct -b 0.25 -g genome.txt | gloci sort -g genome.txt -")]
    Flank(FlankArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci shift genes.bed --by 100 -g hg38
  gloci shift genes.bed --by -50 -g hg38
  gloci shift peaks.bed --by 1000 -s -g hg38 | gloci sort -g hg38 -
  gloci shift peaks.bed --pct --by 0.5 -g genome.txt")]
    Shift(ShiftArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci windows -g hg38 -w 1000000 -o bins.bed
  gloci windows -g genome.txt -w 1000 --step 500
  gloci windows -i genes.bed -n 10 --id-column --reverse
  gloci merge peaks.bed | gloci windows -i - -w 100")]
    Windows(WindowsArgs),
    /// Convert VCF/BCF, GFF3/GTF, WIG or BED12 to BED (bedops convert2bed,
    /// vcf2bed, gff2bed, gtf2bed; bedtools bed12tobed6).
    ///
    /// VCF/BCF gives one row per record: chrom, start, end, ID, REF, ALT,
    /// QUAL, FILTER, with the footprint chosen by --vcf-span.
    /// --split-alleles writes one row per ALT instead, each with that ALT's
    /// own footprint (with --vcf-span ref, the record's REF footprint on
    /// every row; `*`, `<*>` and `<NON_REF>` rows are skipped when other ALTs
    /// exist). --types keeps the records
    /// (with --split-alleles, the ALTs) of the listed kinds: snv, mnv, ins,
    /// del, sv, ref, other; a record is kept when any of its ALTs is of a
    /// listed kind. GFF3/GTF gives 10 columns (bedops gff2bed): chrom,
    /// start, end, name, score, strand, source, type, phase, attributes.
    /// The name is the --name-attr attribute (default ID for GFF3, gene_id
    /// for GTF; `.` when absent). --feature keeps the listed column-3 types.
    /// WIG (fixedStep or variableStep, with optional span) gives one bedGraph
    /// row per data line: chrom, start, end, value; equal values are not
    /// merged. fixedStep needs start= and step=; values must be finite. A
    /// span larger than the step gives overlapping rows.
    /// BED12 with --blocks gives one BED6 row per block; a line that is not
    /// valid BED12 (fewer than 12 columns, block counts that disagree, blocks
    /// out of order, a last block not ending at chromEnd) is an error, also
    /// with --skip-invalid. The input is
    /// converted in its own order and needs no sorting. A conversion whose
    /// output goes out of order (allele footprints can start after their
    /// record) prints one warning; --write-index then fails.
    #[command(after_long_help = "\
Examples:
  gloci convert calls.vcf.gz -o calls.bed
  gloci convert calls.vcf.gz --split-alleles --types del,ins
  gloci convert genes.gtf --feature exon --name-attr gene_name
  gloci convert transcripts.bed --blocks -o exons.bed")]
    Convert(ConvertArgs),
    /// Remove from A the bases covered by B (bedtools subtract, bedops -d).
    ///
    /// Each A record is cut by the union of the B records that overlap it,
    /// and what is left is written as A's own record with new coordinates:
    /// one record per piece, in ascending order (BED columns 2-3, GFF
    /// columns 4-5). An A that nothing touches is written as read, and an A
    /// that is covered completely disappears. A point has no bases to cut:
    /// it is dropped when a B contains it, and a B point removes nothing
    /// from a cut; under --whole, or --sum-overlaps without -f, a B point
    /// strictly inside A (or on its edge with --inclusive-points) still
    /// removes A (the point rule of intersect).
    /// The criteria (-f, -F, -r, -e, --min-bases, -s, -S) choose which B
    /// count and are measured as in intersect. -w, --left and --right pad B
    /// before it cuts.
    ///
    /// --whole drops an A that has any counting B instead of cutting it
    /// (bedtools -A). --sum-overlaps drops an A when the bases shared with
    /// all the B that pass the strand rule add up to at least -f of A
    /// (bedtools -N); it takes only -f, and without -f it acts like
    /// --whole. --split cuts by the BED12 blocks of B, not its span.
    /// A VCF/BCF A is never cut: a record is dropped whole, in default or
    /// --whole mode when intersect --mode none would drop it (any accepted
    /// hit), or with --sum-overlaps -f when the summed overlaps reach -f of
    /// its footprint, and the output is VCF with its header.
    ///
    /// A and every B must be sorted in one chromosome order, which is
    /// checked as they are read. With --unsorted, any order works, A is
    /// still streamed and written in its own order, and every -b file is
    /// loaded into memory, so memory grows with the size of B.
    ///
    /// Two ways to deal with inputs that are not sorted. Use --sort when the
    /// inputs are too big for memory or you want sorted output: A and every
    /// -b file are sorted first with an external sort bounded by --max-mem
    /// (spilling to --tmp-dir), then the normal sweep runs, so the result
    /// equals gloci sort on each input followed by this command. It needs
    /// disk, about twice the input size in --tmp-dir, and removes it
    /// afterwards. Use --unsorted when you want A's order kept and every
    /// -b file fits in memory. The two cannot be combined.
    ///
    /// When A and every -b file are BGZF files with a .tbi or .csi index and
    /// more than one thread is used, the work runs on runs of chromosomes in
    /// parallel, which are read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF); their
    /// output waiting to be written shares 64 MiB of memory and the rest
    /// spills to --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count.
    #[command(after_long_help = "\
Examples:
  gloci subtract -a peaks.bed -b blacklist.bed
  gloci subtract -a genes.bed -b repeats.bed -f 50% --whole -o clean.bed
  gloci subtract -a reads.bed -b exons.bed --sum-overlaps -f 0.8 -s
  gloci subtract -a calls.vcf.gz -b blacklist.bed -o filtered.vcf.gz")]
    Subtract(SubtractArgs),
    /// Find the nearest B to each A (bedtools closest, bedops closest-features).
    ///
    /// For each A, reports the nearest B as TSV rows: A's record, a label
    /// column when there are several -b files, B's record, and with
    /// --distance the distance. An overlapping B has distance 0, bookended
    /// features 1, and otherwise the gap in bases plus 1. An A with no
    /// candidate gets one row with a placeholder B (`. -1 -1 ...`) and
    /// distance -1. A VCF/BCF A gets a column-name line.
    ///
    /// -k N reports the N nearest, and every B tied with the N-th. With
    /// --ties first or last, only as many of the tied B as N still needs
    /// are kept: the first or the last by (file, input order). Rows are
    /// ordered by distance, then by file and input order.
    /// --signed ref|a|b prints signed distances (and implies --distance):
    /// negative is upstream, where ref means lower coordinates, a follows
    /// A's strand and b follows B's strand; unknown strand counts as +.
    /// --ignore-upstream and --ignore-downstream drop B by that sign in the
    /// same frame (ref when --signed is absent), --ignore-overlaps drops
    /// overlapping B, and --max-distance N drops B farther than N. -s and
    /// -S keep only B on the same or the opposite strand; B and A of
    /// unknown strand never match.
    ///
    /// --each-side writes one row per A with the nearest upstream B and
    /// then the nearest downstream B; with --distance each distance is
    /// signed (as closest-features does). An overlapping B counts on both sides. -k
    /// applies per side, and extra ranks and ties add rows. With several
    /// -b files, the nearest are chosen across all of them, or with
    /// --per-file in each file, whose rows are written in file order; a
    /// file with no candidate for an A writes no row (as bedtools -mdb
    /// each), and an A with no candidate in any file gets one placeholder
    /// row. --each-side with --ignore-upstream or --ignore-downstream
    /// leaves that side empty except for an overlapping B, which still
    /// shows there at distance 0 unless --ignore-overlaps is given.
    ///
    /// Memory and a pathology: the sorted sweep looks ahead in B for the
    /// nearest B to the right of A. Under -s/-S, or with --signed b and
    /// --ignore-upstream, --ignore-downstream or --each-side, a chromosome where one strand
    /// has no B left makes the look-ahead read and buffer all the remaining
    /// B of that chromosome. --max-distance bounds this. With --unsorted
    /// every -b file is loaded into memory instead, so memory grows with
    /// the size of B, and A is streamed and written in its own order.
    ///
    /// Two ways to deal with inputs that are not sorted. Use --sort when the
    /// inputs are too big for memory or you want sorted output: A and every
    /// -b file are sorted first with an external sort bounded by --max-mem
    /// (spilling to --tmp-dir), then the normal sweep runs, so the result
    /// equals gloci sort on each input followed by this command. It needs
    /// disk, about twice the input size in --tmp-dir, and removes it
    /// afterwards. Use --unsorted when you want A's order kept and every
    /// -b file fits in memory. The two cannot be combined.
    ///
    /// When A and every -b file are BGZF files with a .tbi or .csi index and
    /// more than one thread is used, the work runs on runs of chromosomes in
    /// parallel, which are read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF); their
    /// output waiting to be written shares 64 MiB of memory and the rest
    /// spills to --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count.
    #[command(after_long_help = "\
Examples:
  gloci closest -a genes.bed -b peaks.bed --distance
  gloci closest -a genes.bed -b peaks.bed -k 3 --signed a -s
  gloci closest -a snps.vcf.gz -b genes.bed --max-distance 10000 --ignore-overlaps
  gloci closest -a genes.bed -b peaks.bed --each-side --distance")]
    Closest(ClosestArgs),
    /// Aggregate columns of the B records over each A (bedtools map, bedops bedmap).
    ///
    /// For each A, in A's order, folds the B records that overlap it (the
    /// criteria of intersect: -f, -F, -r, -e, --min-bases, -s, -S, --split,
    /// or a window -w/--left/--right) and writes A's record followed by one
    /// column per -c COL:OP, in B file then B input order. COL is a 1-based
    /// column of B, 0 for the whole record (count, bases, indicator, bases_uniq, frac_a,
    /// frac_b and records only),
    /// or for a VCF B QUAL, ID, FILTER or INFO/KEY. Ops: count
    /// count_distinct sum mean median min max stdev var stdev_pop var_pop
    /// collapse distinct distinct_sorted first last bases indicator cv mad
    /// quantile:Q tmean:LO-HI wmean mode antimode bases_uniq frac_a frac_b
    /// min_record max_record records (COL 0 also works with bases_uniq, frac_a,
    /// frac_b and records; records ignores COL). cv is the sample
    /// stdev over the mean; mad the unscaled median absolute deviation;
    /// quantile is type 7 and tmean the mean of the values whose rank
    /// fraction is within LO-HI (a lone value is kept); wmean weights by
    /// overlapping bases; ties of mode and antimode go to the bytewise
    /// smallest value, as in bedtools; bases_uniq counts overlapped
    /// bases once; frac_a is bases_uniq over A's length and frac_b the summed
    /// overlap over the summed B length (--empty when the length is 0);
    /// records copies whole B records, tabs shown as |, and does not escape --delim
    /// where it occurs in the data. count counts every
    /// overlapping B; the other ops skip `.` and empty values, and a value
    /// that is not a number in a numeric op is an error. sum, min and max
    /// print integers when every value was one; stdev and var are sample
    /// statistics (one value gives --empty); stdev_pop and var_pop divide by n
    /// (one value gives 0), as bedtools stdev does. bases is the overlap of A with
    /// each B's own span, which is 0 for a B that is within a window only;
    /// indicator is 1 when any B counts, else 0. An A with no B gets 0 for
    /// count and count_distinct and --empty (default `.`) for the rest;
    /// --skip-unmatched drops such A instead.
    ///
    /// Columns are resolved against the first -b file, and all -b files must
    /// be one format family. The output is TSV with A's columns first; a
    /// VCF/BCF A gets a column-name line ending with the -c labels.
    ///
    /// A and every B must be sorted in one chromosome order, which is
    /// checked as they are read. With --unsorted, any order works, A is
    /// still streamed in its own order, and every -b file is loaded into
    /// memory, so memory grows with the size of B.
    ///
    /// Two ways to deal with inputs that are not sorted. Use --sort when the
    /// inputs are too big for memory or you want sorted output: A and every
    /// -b file are sorted first with an external sort bounded by --max-mem
    /// (spilling to --tmp-dir), then the normal sweep runs, so the result
    /// equals gloci sort on each input followed by this command. It needs
    /// disk, about twice the input size in --tmp-dir, and removes it
    /// afterwards. Use --unsorted when you want A's order kept and every
    /// -b file fits in memory. The two cannot be combined.
    ///
    /// When A and every -b file are BGZF files with a .tbi or .csi index and
    /// more than one thread is used, the work runs on runs of chromosomes in
    /// parallel, which are read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF); their
    /// output waiting to be written shares 64 MiB of memory and the rest
    /// spills to --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count.
    #[command(after_long_help = "\
Examples:
  gloci map -a genes.bed -b scores.bedgraph -c 4:mean
  gloci map -a windows.bed -b peaks.bed -c 0:count -c 5:max --skip-unmatched
  gloci map -a snps.vcf.gz -b regions.bed -c 4:distinct --delim ';'
  gloci map -a genes.bed -b peaks.bed -w 100 -c 0:count -c 0:bases")]
    Map(MapArgs),
    /// How much of each A the B records cover (bedtools coverage and annotate).
    ///
    /// For each A, in A's order, counts the B records that overlap it (the
    /// criteria of intersect: -f, -F, -r, -e, --min-bases, -s, -S, --split,
    /// --inclusive-points) and writes A's record followed by the number of
    /// hits, the bases of A covered by at least one hit, A's length and the
    /// covered fraction. --hist writes one row per depth present in A
    /// (depth, bases at that depth, A's length, fraction) and, after the
    /// last A, `all` rows with the totals. --per-base writes one row per
    /// base of A (1-based offset, depth). --mean writes the mean depth over
    /// A. With several -b files only the default report is allowed and each
    /// file gives one covered fraction (bedtools annotate).
    ///
    /// With --split, A and B are their BED12 blocks: depth counts the bases
    /// shared by an A block and a B block, A's length is the sum of its
    /// blocks, and --per-base offsets run along the blocks. A point A covers
    /// no bases: it reports its hit count, 0, 0 and `.`, --mean gives `.`,
    /// and --hist and --per-base write no rows for it. Fractions and means
    /// print in the shortest exact decimal form (0.3, not 0.3000000).
    ///
    /// A and every B must be sorted in one chromosome order, which is
    /// checked as they are read. With --unsorted, any order works, A is
    /// still streamed in its own order, and every -b file is loaded into
    /// memory, so memory grows with the size of B. The output is TSV with
    /// A's columns first.
    ///
    /// Two ways to deal with inputs that are not sorted. Use --sort when the
    /// inputs are too big for memory or you want sorted output: A and every
    /// -b file are sorted first with an external sort bounded by --max-mem
    /// (spilling to --tmp-dir), then the normal sweep runs, so the result
    /// equals gloci sort on each input followed by this command. It needs
    /// disk, about twice the input size in --tmp-dir, and removes it
    /// afterwards. Use --unsorted when you want A's order kept and every
    /// -b file fits in memory. The two cannot be combined.
    ///
    /// When A and every -b file are BGZF files with a .tbi or .csi index and
    /// more than one thread is used, the work runs on runs of chromosomes in
    /// parallel, which are read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF); their
    /// output waiting to be written shares 64 MiB of memory and the rest
    /// spills to --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count.
    #[command(after_long_help = "\
Examples:
  gloci coverage -a genes.bed -b reads.bed
  gloci coverage -a genes.bed -b reads.bed --hist
  gloci coverage -a exons.bed -b reads.bed --per-base -s
  gloci coverage -a peaks.bed -b rep1.bed rep2.bed --names r1,r2")]
    Coverage(CoverageArgs),
    /// Aggregate columns over runs of lines with equal keys (bedtools groupby).
    ///
    /// Reads any tab-separated text, BED included, and groups consecutive
    /// lines whose -g key columns are equal; the input does not have to be
    /// sorted, and a key that comes back later starts a new group (as in
    /// bedtools groupby; sort first to get one group per key). Each group
    /// writes its key columns, exactly as in its first line, followed by one
    /// column per -c COL:OP, with the ops of map (count, count_distinct, sum,
    /// mean, median, min, max, stdev, var, stdev_pop, var_pop, collapse,
    /// distinct, distinct_sorted, first, last, indicator, cv, mad, quantile:Q,
    /// tmean:LO-HI, mode, antimode); the ops that need a second input
    /// (wmean, records, min_record, max_record, frac_*, bases_uniq, bases)
    /// are for map. Blank lines, `#` lines and track/browser lines are
    /// skipped. A line with fewer columns than -g needs is an error naming
    /// the line.
    ///
    /// With --header the first line names the columns: -g and -c may use the
    /// names, and the output starts with a line holding the key names and
    /// the -c labels (bedtools -header copies the input header unchanged).
    /// VCF and BCF input is refused; extract the fields first with gloci vcf
    /// extract. --threads is accepted for uniformity with the other commands;
    /// groupby runs on a single thread.
    #[command(after_long_help = "\
Examples:
  gloci groupby scores.bed -g 1,2 -c 5:sum -c 5:mean
  gloci groupby calls.tsv -g 1 -c 4:collapse -c 5:count_distinct
  gloci groupby table.tsv --header -g gene -c expr:median -o medians.tsv")]
    Groupby(GroupbyArgs),
    /// VCF operations (spec §7.8).
    #[command(after_long_help = "\
Examples:
  gloci vcf filter cohort.bcf --pass-only --min-qual 30 -o filtered.bcf
  gloci vcf filter input.vcf --samples s1,s2 --min-dp 10 --recalc-ac -o out.vcf
  gloci vcf filter variants.bcf --types indel --biallelic")]
    Vcf(Box<VcfArgs>),
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args, Debug, Clone)]
pub struct VcfArgs {
    #[command(subcommand)]
    pub command: VcfCommand,
}

#[derive(Subcommand, Debug, Clone)]
pub enum VcfCommand {
    /// Filter VCF/BCF files with the selection pipeline.
    #[command(after_long_help = "\
Examples:
  gloci vcf filter cohort.bcf --pass-only --min-qual 30 --min-maf 0.01 --max-missing 0.1 -o filtered.bcf
  gloci vcf filter input.vcf --samples s1,s2 --min-dp 10 --recalc-ac -o out.vcf
  gloci vcf filter variants.bcf --types indel --biallelic")]
    Filter(VcfFilterArgs),
    /// Extract INFO or FORMAT values as a table, or the 012 genotype matrix.
    #[command(after_long_help = "\
Examples:
  gloci vcf extract cohort.vcf.gz --info DP,AF -o info.tsv
  gloci vcf extract cohort.bcf --format GT,DP --layout long --samples s1,s2
  gloci vcf extract cohort.vcf.gz --genotype-matrix 012 --out-prefix cohort --biallelic")]
    Extract(VcfExtractArgs),
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

#[derive(Args, Debug, Clone)]
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
    /// vcftools --012. The matrix is held in memory: one byte per sample and
    /// site.
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

    #[command(flatten)]
    pub vcf_span: VcfSpanArg,
    #[command(flatten)]
    pub chrom_map: ChromMapArg,
    #[command(flatten)]
    pub threads: ThreadArgs,
    #[command(flatten)]
    pub invalid: InvalidArg,
}

#[derive(Args, Debug, Clone)]
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
pub struct GenomecovArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Chromosome order and sizes (required): a built-in genome such as
    /// hg38 (see gloci genomes), or a .fai, .genome or .dict file.
    #[arg(short, long, value_name = "NAME|FILE", required = true)]
    pub genome: String,
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
    /// Count the BED12 blocks of each feature instead of the whole feature.
    #[arg(long, conflicts_with_all = ["five_prime", "three_prime"])]
    pub split: bool,
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
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

/// A variant kind for `convert --types` (spec 4.2).
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum KindArg {
    Snv,
    Mnv,
    Ins,
    Del,
    Sv,
    Ref,
    Other,
}

impl From<KindArg> for gloci_core::VariantKind {
    fn from(k: KindArg) -> Self {
        use gloci_core::VariantKind as K;
        match k {
            KindArg::Snv => K::Snv,
            KindArg::Mnv => K::Mnv,
            KindArg::Ins => K::Ins,
            KindArg::Del => K::Del,
            KindArg::Sv => K::Sv,
            KindArg::Ref => K::Ref,
            KindArg::Other => K::Other,
        }
    }
}

#[derive(Args)]
pub struct ConvertArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// VCF/BCF: one row per ALT, with that ALT's footprint.
    #[arg(long)]
    pub split_alleles: bool,
    /// VCF/BCF: keep records (with --split-alleles, ALTs) of these kinds.
    #[arg(long, value_enum, value_delimiter = ',', value_name = "LIST")]
    pub types: Vec<KindArg>,
    /// GFF3/GTF: keep this column-3 type (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',', value_name = "NAME")]
    pub feature: Vec<String>,
    /// GFF3/GTF: the attribute for the name column (default ID for GFF3,
    /// gene_id for GTF).
    #[arg(long, value_name = "KEY", value_parser = clap::builder::NonEmptyStringValueParser::new())]
    pub name_attr: Option<String>,
    /// BED12: write one BED6 row per block.
    #[arg(long)]
    pub blocks: bool,
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

#[derive(Args)]
#[command(group(clap::ArgGroup::new("source").required(true).multiple(false).args(["genome", "intervals"])))]
#[command(group(clap::ArgGroup::new("split").required(true).multiple(false).args(["size", "count"])))]
pub struct WindowsArgs {
    /// Cut every chromosome of a built-in genome such as hg38 (see gloci
    /// genomes), or a .fai, .genome or .dict file with lengths.
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
    pub chrom_map: ChromMapArg,
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
    pub chrom_map: ChromMapArg,
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
    pub chrom_map: ChromMapArg,
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

/// `--chrom-map` (D62), on every command that reads intervals or variants.
#[derive(Args, Debug, Clone)]
pub struct ChromMapArg {
    /// Rename chromosomes on input with a two-column table (old name, new
    /// name; tab-separated, `#` comments allowed). Unlisted names are kept.
    /// A file that has both an old name and its new name is refused.
    /// Order checks, -g and the output all see the new names; -g itself is
    /// never renamed.
    #[arg(long, value_name = "FILE")]
    pub chrom_map: Option<String>,
}

impl ChromMapArg {
    /// Read the table, when one was given.
    pub fn load(&self) -> Result<Option<std::sync::Arc<gloci_io::ChromMap>>, gloci_io::IoError> {
        self.chrom_map
            .as_deref()
            .map(|p| gloci_io::ChromMap::load(p).map(std::sync::Arc::new))
            .transpose()
    }
}

#[derive(Args, Debug, Clone)]
pub struct InvalidArg {
    /// Skip invalid records with a warning instead of failing.
    #[arg(long)]
    pub skip_invalid: bool,
}

/// `--sort` and the options that bound it (D39).
#[derive(Args, Debug, Clone)]
pub struct SortArg {
    /// Sort the inputs first, for inputs that are unsorted or too large to
    /// sort in memory (external sort bounded by --max-mem). The sorted
    /// records are written to a temporary file in --tmp-dir and read back,
    /// so besides the spill runs (about the size of the input, compressed)
    /// the directory needs one uncompressed sorted copy: allow roughly twice
    /// the input size. Everything is removed afterwards.
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

/// `--sort` and temporary directory options for relation commands.
#[derive(Args, Debug, Clone)]
pub struct RelationSortArg {
    /// Sort the inputs first, for inputs that are unsorted or too large to
    /// sort in memory (external sort bounded by --max-mem). The sorted
    /// records are written to a temporary file in --tmp-dir and read back,
    /// so besides the spill runs (about the size of the input, compressed)
    /// the directory needs one uncompressed sorted copy: allow roughly twice
    /// the input size. Everything is removed afterwards.
    #[arg(long)]
    pub sort: bool,
    /// Memory budget for --sort, in binary units (K, M, G, T); at least 16M.
    #[arg(long, value_name = "SIZE", default_value = "2G", value_parser = parse_size, requires = "sort")]
    pub max_mem: u64,
    /// Directory for temporary files: --sort's, and the spills of
    /// per-chromosome runs on indexed inputs; ignored when nothing spills. It
    /// must exist and be writable (default: $GLOCI_TMPDIR if set and
    /// non-empty, else the system temp dir).
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<PathBuf>,
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
    /// Restrict to a region such as chr1:1,000-2,000 (repeatable). In a
    /// relation command only A is restricted, to the records that overlap a
    /// region (a point on the edge with --inclusive-points); B is read in
    /// full, so an A record reaching out of a region still meets B there.
    #[arg(long, value_name = "REGION")]
    pub region: Vec<String>,
    /// Read regions from FILE (BED: 0-based, half-open); see --region.
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
    /// Accept inputs in any order and keep A's order in the output; every -b
    /// file is held in memory. Cannot be combined with --sort.
    #[arg(long, conflicts_with = "sort")]
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
        #[command(flatten)]
        sort: SortArg,
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
