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
    /// are recorded in crates/gloci-io/data/genomes/SOURCES.toml in the
    /// source tree.
    #[command(after_long_help = "\
Examples:
  gloci genomes
  gloci genomes show hg38
  gloci sort peaks.bed -g hg38 > peaks.sorted.bed")]
    Genomes(GenomesArgs),
    /// Build a .tbi/.csi index for a sorted BGZF file, a .bai/.csi for a BAM, a .crai for a CRAM,
    /// or a .fai for a FASTA.
    ///
    /// A coordinate-sorted BAM gets PATH.bai, or PATH.csi with --csi or when a
    /// reference is longer than 2^29; samtools and gloci read each other's
    /// indexes. SAM text has no index: write BAM with gloci view x.sam -o x.bam.
    ///
    /// A CRAM file gets PATH.crai, with the entries samtools index writes;
    /// it needs no reference, even for slices that hold several references.
    /// A FASTA file (plain or BGZF) gets PATH.fai, and PATH.gzi when it is
    /// BGZF, byte-identical to samtools faidx. Every line of a sequence but
    /// the last must have the same length. Reading a FASTA at random
    /// (getfasta, nuc, a CRAM reference) needs these files and never builds
    /// them itself; maskfasta streams its FASTA and needs no index. Plain gzip is refused: recompress with
    /// bgzip. Otherwise, BED, GFF/GTF, VCF and BCF files are recognised by content. Any other
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
  gloci index sorted.bam
  gloci index reads.cram
  gloci index scores.tsv.gz -s 1 -b 2 -e 3 --zero-based")]
    Index(IndexArgs),
    /// Sort BED/bedGraph/GFF/VCF/BCF files, or BAM/SAM/CRAM alignments, with
    /// bounded memory.
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
    ///
    /// Alignments (one BAM, SAM or CRAM file, or stdin) are sorted as
    /// samtools sort does: by reference in header order (unplaced reads
    /// last), position and strand, ties in input order; -n by read name in
    /// natural order (r2 before r10) and -N by plain byte order, each with
    /// READ1 before READ2 and primary before supplementary before secondary.
    /// Every record is kept, unmapped reads included. The header's @HD line
    /// gets SO (and SS for -n and -N) and gloci adds an @PG line. The output
    /// keeps the input's format (BAM for BAM and CRAM, SAM for SAM; -o and -O
    /// override), where samtools sort always writes BAM. --write-index
    /// works for coordinate-sorted BAM output only. -g, --full-key, --unique and
    /// --duplicates do not apply: the order is the header's. To merge sorted
    /// BAM files use sort -m (or cat for files already in order); -n reads.bam sorts by name.
    #[command(after_long_help = "\
Examples:
  gloci sort reads.bam -o reads.sorted.bam --write-index -t 8
  gloci sort peaks.bed > peaks.sorted.bed
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
    /// BAM, CRAM and SAM input is written as BAM or SAM, as -O or
    /// the -o extension say (the same as the input by default; BAM for
    /// CRAM), every record, unmapped reads included, with an @PG line added
    /// to the header. --header-only, --no-header (SAM output only),
    /// --min-mapq, --require-flags and --exclude-flags apply. Regions select
    /// the reads whose footprint overlaps them; an unmapped read placed at a
    /// position counts as one base there. It is not samtools view: there is
    /// no -b, -h, -c or -L. CRAM needs the FASTA it was made with
    /// (--reference, with its .fai), unless it was written without one; MD
    /// and NM tags come out as stored, not regenerated as samtools does.
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
  gloci view reads.bam -O sam | head
  gloci view reads.cram --reference ref.fa -o reads.bam")]
    View(ViewArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci cat lane1.bam lane2.bam -o all.bam
  gloci cat part-*.bam > all.bam")]
    Cat(CatArgs),
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
    /// BAM, CRAM and SAM (A or B) work as intervals (a CRAM needs --reference, and an
    /// indexed BAM or CRAM is read per chromosome at -t above 1): a read is the reference span of its
    /// M, D, N, = and X operations (clips never widen it), on the strand of its
    /// flag. Unmapped reads are skipped; --min-mapq, --require-flags and
    /// --exclude-flags filter reads first. --split uses the segments between
    /// N operations, --split-deletions also cuts at D. With an alignment A,
    /// --mode any and none write the reads themselves as BAM (SAM for a SAM A,
    /// or as -o/-O say); pairs and count write A's reads as BED6 rows, as
    /// gloci convert does. --mode clip of reads is refused. Alignments stream:
    /// --sort sorts them into a temporary BAM first, bounded by --max-mem.
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
    /// When A and every -b file are BGZF files with a .tbi or .csi index (a
    /// BAM with a .bai or .csi), or plain (uncompressed) BED, bedGraph or VCF
    /// files, which are scanned once
    /// for where each chromosome's records lie, and more than one thread is
    /// used, the work runs on runs of chromosomes in parallel, which are
    /// read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF). Memory is
    /// about 64 MiB plus 20 MiB per run in flight (more with many -b files):
    /// the output waiting to be written shares 64 MiB and the rest spills to
    /// --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count. After
    /// an error, how much output was already written to standard output
    /// varies with the thread count (-o FILE is written atomically). Plain
    /// inputs are validated to their end at every thread count; for indexed
    /// inputs, records outside the queried chromosomes are validated only at
    /// -t 1. With -v, the counts of reads left out of an indexed -b (or of
    /// anything under --region) are of the reads read, which depend on the
    /// thread count. When more than one input is malformed, the error reported may
    /// name a different input at different thread counts; the exit code is
    /// the same.
    #[command(after_long_help = "\
Examples:
  gloci intersect -a peaks.bed -b genes.bed
  gloci intersect -a calls.vcf.gz -b exons.bed --mode none -o off_target.vcf.gz
  gloci intersect -a a.bed -b x.bed y.bed --mode pairs --overlap-bases --names X,Y
  gloci intersect -a reads.bam -b targets.bed -o on_target.bam")]
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
    /// Give every record the id of its cluster (bedtools cluster).
    ///
    /// Writes each input line unchanged followed by a tab and a cluster id.
    /// Clusters are those of gloci merge, with the same -d/--max-gap
    /// (default 0: overlapping and bookended records share a cluster; a
    /// negative gap needs that many shared bases) and -s (each strand
    /// clustered separately, unknown strand apart from + and -). Ids start
    /// at 1 and increase through the file, across chromosomes. With -s, a
    /// cluster takes the next id when its first record is read, and the
    /// records stay in input order (bedtools cluster -s regroups them by
    /// strand within each chromosome).
    ///
    /// The output is TSV, not a copy of the input format: a VCF input
    /// (clustered by its footprints, see --vcf-span) gets a #CHROM column
    /// line ending in "cluster" and then its records with the id appended.
    /// --header keeps the header lines of a BED or GFF input. The input
    /// must be sorted in one chromosome order, checked as it is read; pass
    /// --sort to sort it first.
    #[command(after_long_help = "Examples:
  gloci cluster peaks.bed
  gloci cluster peaks.bed -d 100 -s
  gloci cluster calls.vcf.gz --vcf-span ref -d 50")]
    Cluster(ClusterArgs),
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
    #[command(after_long_help = "Examples:
  gloci spacing peaks.bed
  gloci spacing calls.vcf.gz --vcf-span ref")]
    Spacing(SpacingArgs),
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
    /// Depth of coverage of a genome by the features of a BED file or the
    /// reads of a BAM (bedtools genomecov -i and -ibam).
    ///
    /// Every feature adds one to the depth over its bases. -g is required
    /// (unless the input is BAM, CRAM or SAM):
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
    /// feature past the end of its chromosome is an error.
    ///
    /// BAM, CRAM and SAM input (a CRAM needs --reference): the chromosomes and sizes come from the header
    /// (header order, empty chromosomes included) unless -g is given. A read
    /// counts its whole reference span, deletions included; --split counts
    /// the blocks between N operations, and --split-deletions also cuts at
    /// D. --split-deletions is bedtools -ibam -split, so use it for parity;
    /// plain -ibam (cuts at D, not N) has no gloci counterpart, see
    /// docs/differences.md#bam-genomecov-deletions.
    /// Unmapped reads are left out, and --min-mapq, --require-flags and
    /// --exclude-flags filter reads first. --strand follows the read's
    /// direction. An indexed BAM or CRAM at -t above 1 runs per chromosome,
    /// with the same output as -t 1. At most 8 runs are in flight; memory is
    /// then about 64 MiB plus 20 MiB per run in flight: the output waiting to
    /// be written shares 64 MiB and the rest spills to --tmp-dir.
    ///
    /// The input must be sorted in the genome order, which is checked as it
    /// is read; an unsorted input is an error whose hint names --sort.
    /// Memory is bounded by the features open at one position, whatever
    /// the output mode.
    #[command(after_long_help = "\
Examples:
  gloci genomecov peaks.bed -g hg38
  gloci genomecov reads.bed -g hg38 --bedgraph --scale 0.000001 -o cov.bedgraph
  gloci genomecov exons12.bed -g hg38 --split --per-base
  gloci genomecov reads.bam --bga --split")]
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
    /// Convert VCF/BCF, GFF3/GTF, WIG, PSL, RepeatMasker .out, BED12 or
    /// BAM/CRAM/SAM to BED (bedops convert2bed, vcf2bed, gff2bed, gtf2bed;
    /// bedtools bed12tobed6, bamtobed).
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
    /// PSL (BLAT; the header is skipped) gives bedops psl2bed's 21 columns:
    /// tName, tStart, tEnd, qName, matches, strand, then the other PSL
    /// columns in order. RepeatMasker .out gives bedops rmsk2bed's columns:
    /// sequence, begin-1, end, repeat, score, strand (`C` is `-`), the
    /// divergence, deletion and insertion percentages, the `(left)` of the
    /// query, class/family, the repeat positions as written, ID and `*`.
    /// BAM, CRAM and SAM give one BED6 row per mapped read, in input order:
    /// chrom, start, end, name (with /1 or /2 for the first or second mate),
    /// MAPQ, strand; unmapped reads are skipped and counted at -v, and
    /// --min-mapq, --require-flags and --exclude-flags filter. The span is
    /// the reference bases the CIGAR covers (M D N = X), cut at the end of
    /// the reference. --bed12 writes one BED12 row per read with the thick
    /// span equal to the read and --color as itemRgb; --blocks writes one
    /// BED6 row per block, cut at N; --split-deletions cuts at D as well.
    /// --score-tag NM takes the score from an integer tag (NM is bedtools
    /// -ed), and --cigar appends the CIGAR to BED6 rows.
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
  gloci convert transcripts.bed --blocks -o exons.bed
  gloci convert reads.bam --bed12 -o reads.bed
  gloci convert reads.bam --blocks --score-tag NM --min-mapq 30")]
    Convert(ConvertArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci getfasta --fasta ref.fa genes.bed
  gloci getfasta --fasta ref.fa genes.bed12 --split -s --name-only -o transcripts.fa
  gloci getfasta --fasta ref.fa.gz calls.vcf.gz --vcf-span ref --tab")]
    Getfasta(GetfastaArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci nuc --fasta ref.fa genes.bed
  gloci nuc --fasta ref.fa genes.bed -s --seq -o genes.nuc.tsv
  gloci nuc --fasta ref.fa.gz calls.vcf.gz --vcf-span ref --pattern CG")]
    Nuc(NucArgs),
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
    /// When A and every -b file are BGZF files with a .tbi or .csi index (a
    /// BAM with a .bai or .csi), or plain (uncompressed) BED, bedGraph or VCF
    /// files, which are scanned once
    /// for where each chromosome's records lie, and more than one thread is
    /// used, the work runs on runs of chromosomes in parallel, which are
    /// read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF). Memory is
    /// about 64 MiB plus 20 MiB per run in flight (more with many -b files):
    /// the output waiting to be written shares 64 MiB and the rest spills to
    /// --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count. After
    /// an error, how much output was already written to standard output
    /// varies with the thread count (-o FILE is written atomically). Plain
    /// inputs are validated to their end at every thread count; for indexed
    /// inputs, records outside the queried chromosomes are validated only at
    /// -t 1. With -v, the counts of reads left out of an indexed -b (or of
    /// anything under --region) are of the reads read, which depend on the
    /// thread count. When more than one input is malformed, the error reported may
    /// name a different input at different thread counts; the exit code is
    /// the same.
    #[command(after_long_help = "\
Examples:
  gloci subtract -a peaks.bed -b blacklist.bed
  gloci subtract -a reads.bed -b exons.bed --sum-overlaps -f 0.8 -s
  gloci subtract -a calls.vcf.gz -b blacklist.bed -o filtered.vcf.gz
  gloci subtract -a calls.vcf.gz -b known.vcf.gz --match alleles -o novel.vcf.gz")]
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
    /// When A and every -b file are BGZF files with a .tbi or .csi index (a
    /// BAM with a .bai or .csi), or plain (uncompressed) BED, bedGraph or VCF
    /// files, which are scanned once
    /// for where each chromosome's records lie, and more than one thread is
    /// used, the work runs on runs of chromosomes in parallel, which are
    /// read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF). Memory is
    /// about 64 MiB plus 20 MiB per run in flight (more with many -b files):
    /// the output waiting to be written shares 64 MiB and the rest spills to
    /// --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count. After
    /// an error, how much output was already written to standard output
    /// varies with the thread count (-o FILE is written atomically). Plain
    /// inputs are validated to their end at every thread count; for indexed
    /// inputs, records outside the queried chromosomes are validated only at
    /// -t 1. With -v, the counts of reads left out of an indexed -b (or of
    /// anything under --region) are of the reads read, which depend on the
    /// thread count. When more than one input is malformed, the error reported may
    /// name a different input at different thread counts; the exit code is
    /// the same.
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
    /// When A and every -b file are BGZF files with a .tbi or .csi index (a
    /// BAM with a .bai or .csi), or plain (uncompressed) BED, bedGraph or VCF
    /// files, which are scanned once
    /// for where each chromosome's records lie, and more than one thread is
    /// used, the work runs on runs of chromosomes in parallel, which are
    /// read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF). Memory is
    /// about 64 MiB plus 20 MiB per run in flight (more with many -b files):
    /// the output waiting to be written shares 64 MiB and the rest spills to
    /// --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count. After
    /// an error, how much output was already written to standard output
    /// varies with the thread count (-o FILE is written atomically). Plain
    /// inputs are validated to their end at every thread count; for indexed
    /// inputs, records outside the queried chromosomes are validated only at
    /// -t 1. With -v, the counts of reads left out of an indexed -b (or of
    /// anything under --region) are of the reads read, which depend on the
    /// thread count. When more than one input is malformed, the error reported may
    /// name a different input at different thread counts; the exit code is
    /// the same.
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
    /// BAM, CRAM and SAM (A or B) work as intervals (a CRAM needs --reference, and an
    /// indexed BAM or CRAM is read per chromosome at -t above 1): a read is the reference span of its
    /// M, D, N, = and X operations (clips never widen it). Unmapped reads are
    /// skipped; --min-mapq, --require-flags and --exclude-flags filter reads
    /// first. With --split the blocks of a read are the segments between N
    /// operations; --split-deletions also cuts at D. An alignment A is written
    /// as its BED6 row (as gloci convert does). Alignments stream: --sort sorts
    /// them into a temporary BAM first, bounded by --max-mem.
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
    /// When A and every -b file are BGZF files with a .tbi or .csi index (a
    /// BAM with a .bai or .csi), or plain (uncompressed) BED, bedGraph or VCF
    /// files, which are scanned once
    /// for where each chromosome's records lie, and more than one thread is
    /// used, the work runs on runs of chromosomes in parallel, which are
    /// read and written in the chromosome-table order
    /// (karyotype by default; -g, a VCF header or --chrom-order set it).
    /// At most 8 runs are in flight (more threads speed up BGZF). Memory is
    /// about 64 MiB plus 20 MiB per run in flight (more with many -b files):
    /// the output waiting to be written shares 64 MiB and the rest spills to
    /// --tmp-dir. Indexed and streaming runs give identical output
    /// for inputs in chromosome-table order; an input that is not in table
    /// order is refused with the same error whatever the thread count. After
    /// an error, how much output was already written to standard output
    /// varies with the thread count (-o FILE is written atomically). Plain
    /// inputs are validated to their end at every thread count; for indexed
    /// inputs, records outside the queried chromosomes are validated only at
    /// -t 1. With -v, the counts of reads left out of an indexed -b (or of
    /// anything under --region) are of the reads read, which depend on the
    /// thread count. When more than one input is malformed, the error reported may
    /// name a different input at different thread counts; the exit code is
    /// the same.
    #[command(after_long_help = "\
Examples:
  gloci coverage -a genes.bed -b reads.bed
  gloci coverage -a exons.bed -b reads.bed --per-base -s
  gloci coverage -a peaks.bed -b rep1.bed rep2.bed --names r1,r2
  gloci coverage -a targets.bed -b reads.bam --split")]
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
    #[command(after_long_help = "\
Examples:
  gloci maskfasta --fasta hg38.fa -b repeats.bed -o hg38.masked.fa
  gloci maskfasta --fasta hg38.fa.gz -b rmsk.bed segdups.bed --soft -o soft.fa.gz --write-index
  gloci maskfasta --fasta ref.fa -b low-mappability.bed --mask-char X | head")]
    Maskfasta(MaskfastaArgs),
    /// VCF operations (spec §7.8).
    #[command(after_long_help = "\
Examples:
  gloci vcf filter cohort.bcf --pass-only --min-qual 30 -o filtered.bcf
  gloci vcf filter input.vcf --samples s1,s2 --min-dp 10 --recalc-ac -o out.vcf
  gloci vcf filter variants.bcf --types indel --biallelic")]
    Vcf(Box<VcfArgs>),
    // P3b Task 2
    /// Statistics of interval sets (bedtools jaccard).
    #[command(after_long_help = "\
Examples:
  gloci stats jaccard -a peaks.bed -b genes.bed
  gloci stats jaccard -a rep1.bed -b rep2.bed -f 0.5 -r")]
    Stats(StatsArgs),
    // P3b Task 7
    /// BEDPE pairs by how their ends overlap B (bedtools pairtobed).
    ///
    /// A is a BEDPE file (chrom1 start1 end1 chrom2 start2 end2, then optionally name, score,
    /// strand1, strand2 and more; it is never detected, so a BED file of pairs is read
    /// as BEDPE only here). B is BED, GFF or VCF (a VCF as its footprints, --vcf-span) and is
    /// loaded into memory; A streams in any order. Each present end of a pair is tested against B
    /// as an interval under -f, -F, -r, -e, --min-bases, -s, -S and --inclusive-points (section
    /// 4.3: a point is inside only strictly, and unknown strands never match -s or -S); -s and -S
    /// use strand1 for end 1 and strand2 for end 2. An end whose chromosome is `.` overlaps
    /// nothing.
    ///
    /// --type chooses what is written. Rows with a B are A's line, a tab and B's line (TSV);
    /// the A-alone types write A's BEDPE line. either (default): A with each B hitting end 1, then
    /// each hitting end 2 (a B hitting both ends is written twice). both: the same, when both ends
    /// hit. xor: the same, when exactly one end hits. neither: A alone, when no end hits. notboth:
    /// A with each B, when exactly one end hits, and A alone padded with a placeholder B (`.`,
    /// `-1`, `-1` ...) when none does. ispan and ospan: A with each B overlapping the gap between
    /// the ends, or the span from the first start to the last end; notispan and notospan: A alone
    /// when none does. The span types need both ends on one chromosome (a pair with one end absent or the
    /// ends on two chromosomes is left out, counted in one warning), an empty gap (touching or overlapping ends) hits nothing, and
    /// -s/-S with them is an error. Hits of one A come in B's input order; output follows A's
    /// order whatever the thread count. Output is TSV or BEDPE (-o x.bedpe, x.tsv, .gz for
    /// BGZF); header lines of A are not written.
    #[command(after_long_help = "\
Examples:
  gloci pair-intersect -a sv.bedpe -b genes.bed
  gloci pair-intersect -a sv.bedpe -b exons.bed --type both -f 0.5
  gloci pair-intersect -a sv.bedpe -b repeats.bed --type neither -o clean.bedpe
  gloci pair-intersect -a sv.bedpe -b genes.bed --type xor -s -t 4")]
    PairIntersect(PairIntersectArgs),
    // P3b Task 8
    /// BEDPE pairs by how their ends overlap the ends of the pairs of B (bedtools pairtopair).
    ///
    /// A and B are both BEDPE files (never detected). B is loaded into memory; A streams in any
    /// order. Orientation 1 tests end 1 of A against end 1 of B and end 2 against end 2;
    /// orientation 2 tests end 1 of A against end 2 of B and end 2 against end 1. --type both
    /// (default): A with each B that overlaps on both ends in either orientation. either: A with
    /// each B with any end overlapping an end of A. neither: A alone when no B end overlaps.
    /// notboth: A alone when no B overlaps on both ends. A B is written once per A, as its line
    /// reads (TSV: A's line, a tab, B's line), in B's input order; the A-alone types write A's
    /// BEDPE line. Ends are tested under -f, -F, -r, -e, --min-bases and --inclusive-points
    /// (section 4.3); strand is ignored unless -s or -S, which compare strand1 with strand1 and
    /// strand2 with strand2 of the tested ends (unknown strands never match). -w, --left, --right
    /// and --strand-window pad B's ends before testing (not with the thresholds). An end whose
    /// chromosome is `.` overlaps nothing. --require-different-names drops a B named like A
    /// (column 7); a pair without a name is then an error. Output follows A's order whatever the
    /// thread count; -o x.bedpe, x.tsv, .gz for BGZF.
    #[command(after_long_help = "\
Examples:
  gloci pair-pair -a sv1.bedpe -b sv2.bedpe
  gloci pair-pair -a sv1.bedpe -b sv2.bedpe --type either -w 100
  gloci pair-pair -a sv1.bedpe -b sv2.bedpe --type notboth -o unique.bedpe
  gloci pair-pair -a sv.bedpe -b sv.bedpe --require-different-names -s -t 4")]
    PairPair(PairPairArgs),
    // P3b Task 5
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
    #[command(after_long_help = "\
Examples:
  gloci random -g hg38 -l 500 -n 100000 --seed 42 -o random.bed
  gloci random -g genome.txt -l 1000 -n 1000 --seed 7 | gloci sort -g genome.txt - > sorted.bed
  gloci random -g hg38 -n 10000 --seed 1 -o a.bed && gloci shift a.bed --by 5 -g hg38")]
    Random(RandomArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci sample peaks.bed -n 1000 --seed 42 -o subset.bed
  gloci sample peaks.bed -n 500 --strand + --seed 1
  gloci sample calls.vcf.gz -n 10000 --seed 3 -O vcf.gz -o subset.vcf.gz")]
    Sample(SampleArgs),
    // P3b Task 6
    /// Move every record to a random place in the genome (bedtools shuffle).
    ///
    /// Each record keeps its length, its other columns and its strand and goes
    /// to a uniformly random place that keeps all of it inside the allowed
    /// regions: --include regions (default: the whole genome) less --exclude
    /// regions. (bedtools only puts the start inside an --include region; see
    /// docs/differences.md#shuffle-placement.) --same-chrom keeps records on
    /// their chromosome; --chrom-first picks a chromosome at random among those
    /// with room, then a place on it, instead of weighting chromosomes by
    /// their free length; --no-overlapping redraws a place that overlaps an
    /// earlier output record, up to --max-tries times. A record with no place
    /// is an error naming its line, length and the fix (nothing is written to a
    /// file); --skip-unplaceable drops such records and warns with a count.
    /// With --format bedpe each pair is moved: a pair on one chromosome keeps
    /// its end lengths and gap, a pair on two places its ends independently,
    /// and only --exclude, --seed, --max-tries and --skip-unplaceable apply.
    /// Output is in input order and not sorted (pipe to gloci sort); VCF and
    /// alignment input are refused (gloci convert). --include/--exclude read
    /// the names of -g, after --chrom-map (which renames them like the input). The output depends
    /// on --seed only, never on --threads; without --seed a fresh seed is
    /// chosen and printed as a warning. It is not bedtools' output for the same
    /// seed (see docs/differences.md#rng-stream).
    #[command(after_long_help = "\
Examples:
  gloci shuffle peaks.bed -g hg38 --exclude gaps.bed --seed 42 -o null.bed
  gloci shuffle peaks.bed -g genome.txt --include mappable.bed --same-chrom --seed 7 | gloci sort -g genome.txt -
  gloci shuffle loops.bedpe --format bedpe -g hg38 --exclude blacklist.bed --seed 1 -o null.bedpe")]
    Shuffle(ShuffleArgs),
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

#[derive(Args)]
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

#[derive(Args)]
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

#[derive(Args)]
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
    Filter(Box<VcfFilterArgs>),
    /// Extract INFO or FORMAT values as a table, or the 012 genotype matrix.
    #[command(after_long_help = "\
Examples:
  gloci vcf extract cohort.vcf.gz --info DP,AF -o info.tsv
  gloci vcf extract cohort.bcf --format GT,DP --layout long --samples s1,s2
  gloci vcf extract cohort.vcf.gz --genotype-matrix 012 --out-prefix cohort --biallelic")]
    Extract(VcfExtractArgs),
    /// Compute statistics in one pass: per-site, per-sample and whole-file,
    /// written as PREFIX.<ext> in vcftools' formats, plus a summary.
    ///
    /// FORMAT fields are validated lazily: only the fields a requested
    /// statistic or a mask rule reads (GT, and DP, GQ or FT when a flag uses
    /// them) are decoded, so a malformed DP, GQ or FT value that nothing
    /// reads is not an error.
    #[command(after_long_help = "\
Examples:
  gloci vcf stats cohort.vcf.gz --out-prefix qc --freq --missing-indv --het --summary
  gloci vcf stats cohort.bcf --out-prefix qc --depth --site-mean-depth --min-dp 8 --samples s1,s2
  gloci vcf stats cohort.vcf.gz --out-prefix qc --tstv-summary --json --pass-only")]
    Stats(VcfStatsArgs),
    /// Write VCF/BCF files one after another (bcftools concat, vcf-concat).
    #[command(after_long_help = "\
Examples:
  gloci vcf concat chr1.vcf.gz chr2.vcf.gz -o all.vcf.gz
  gloci vcf concat chr*.bcf -o all.bcf
  gloci vcf concat --naive chr1.vcf.gz chr2.vcf.gz -o all.vcf.gz")]
    Concat(VcfConcatArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci vcf split calls.vcf.gz -o calls.split.vcf.gz
  gloci vcf split calls.bcf --types snps -o snps-split.bcf
  gloci vcf split calls.vcf.gz | gloci subtract -a - -b known.vcf.gz --match alleles")]
    Split(VcfSplitArgs),
    // P3b Task 21
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
    #[command(after_long_help = "\
Examples:
  gloci vcf norm calls.vcf.gz -f ref.fa -o calls.norm.vcf.gz
  gloci vcf norm calls.bcf --fasta ref.fa --split -o calls.norm.bcf
  gloci vcf norm calls.vcf.gz -f ref.fa --check-ref x --rm-dup exact | gloci vcf stats -
  gloci vcf norm calls.vcf.gz -f ref.fa --check-ref s -o fixed.vcf.gz")]
    Norm(VcfNormArgs),
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
    #[command(after_long_help = "\
Examples:
  gloci vcf annotate calls.vcf.gz -a dbsnp.vcf.gz -c ID,INFO/RS -o calls.annot.vcf.gz
  gloci vcf annotate calls.vcf.gz -a genes.bed.gz -c CHROM,FROM,TO,GENE --header-lines genes.hdr
  gloci vcf annotate calls.vcf.gz -x INFO,FORMAT/PL --mark-sites +CLEAN
  gloci vcf annotate calls.vcf.gz -I '+%CHROM:%POS:%REF:%FIRST_ALT' -o ids.vcf.gz")]
    Annotate(VcfAnnotateArgs),
    // P3b Task 23
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
    #[command(after_long_help = "\
Examples:
  gloci vcf merge a.vcf.gz b.vcf.gz c.vcf.gz -o cohort.vcf.gz
  gloci vcf merge --file-list shards.txt -m none --force-samples -o cohort.bcf
  gloci vcf merge a.bcf b.bcf -i DP:sum,AF:avg -0 | gloci vcf stats - --out-prefix qc
  gloci vcf merge a.vcf.gz b.vcf.gz -f PASS,. -m id -o by-id.vcf.gz")]
    Merge(VcfMergeArgs),
}

#[derive(Args, Debug, Clone)]
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
#[derive(Args, Debug, Clone)]
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

#[derive(Args)]
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

#[derive(Args)]
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
    /// Chromosome order and sizes (required unless the input is BAM, CRAM or SAM,
    /// which gives its header's): a built-in genome such as hg38 (see gloci
    /// genomes), or a .fai, .genome or .dict file.
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

#[derive(Args)]
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
    /// BED12, BAM, CRAM, SAM: write one BED6 row per block (for alignments the
    /// segments between N operations, as bedtools bamtobed -split).
    #[arg(long)]
    pub blocks: bool,
    /// BAM, CRAM, SAM: write one BED12 row per read, its blocks cut at N
    /// (bedtools bamtobed -bed12).
    #[arg(long)]
    pub bed12: bool,
    /// BAM, CRAM, SAM with --bed12: the itemRgb column, as R,G,B.
    #[arg(long, value_name = "R,G,B")]
    pub color: Option<String>,
    /// BAM, CRAM, SAM: take the score column from this integer tag (NM, AS, ...)
    /// instead of the MAPQ; every read needs it (bedtools bamtobed -tag, -ed is NM).
    #[arg(long, value_name = "TAG")]
    pub score_tag: Option<String>,
    /// BAM, CRAM, SAM: append the CIGAR as a seventh column (BED6 only).
    #[arg(long)]
    pub cigar: bool,
    /// BAM, CRAM, SAM: also cut blocks at deletions (bedtools bamtobed -splitD);
    /// implies --blocks unless --bed12 is given.
    #[arg(long)]
    pub split_deletions: bool,
    #[command(flatten)]
    pub align: AlignArgs,
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
    pub align: AlignArgs,
    #[command(flatten)]
    pub out: OutputArgs,
}

#[derive(Args)]
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
    /// Alignments: sort by read name, in natural order (r2 before r10), as samtools sort -n.
    #[arg(short = 'n', long, conflicts_with = "by_name_bytes")]
    pub by_name: bool,
    /// Alignments: sort by read name as bytes (r10 before r2), as samtools sort -N.
    #[arg(short = 'N', long)]
    pub by_name_bytes: bool,
    /// The FASTA a CRAM input was made with (bedtools --cram-ref), with its .fai (and .gzi
    /// when BGZF). Not needed for a CRAM written without a reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
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
    /// Sorted, BGZF-compressed BED/GFF/VCF/BCF file to index, a FASTA file (.fai, .gzi), or
    /// a CRAM file (.crai).
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
    /// Write a .csi instead of a .bai for a BAM file. A BAM with a reference longer than
    /// 2^29 always gets one.
    #[arg(long)]
    pub csi: bool,
    #[command(flatten)]
    pub threads: ThreadArgs,
    /// The FASTA a CRAM was made with. Indexing a CRAM reads no reference bases, so it is
    /// optional; when given, its sequence lengths are checked against the CRAM's header.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
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
    /// Also write a .tbi/.csi index next to the output (.bai or .csi for BAM output).
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

/// The alignment filters (D101), flattened into each command that reads alignments.
#[derive(Args, Debug, Clone, Default)]
pub struct AlignArgs {
    /// Skip alignments with a mapping quality below N (BAM, CRAM and SAM input).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_mapq: u8,
    /// Keep only alignments that have all of these flags: a number (0x2 or 2)
    /// or names such as paired,proper_pair (as samtools view -f).
    #[arg(long, value_name = "FLAGS")]
    pub require_flags: Option<String>,
    /// Drop alignments that have any of these flags: a number or names such as
    /// unmap,secondary,dup (as samtools view -F).
    #[arg(long, value_name = "FLAGS")]
    pub exclude_flags: Option<String>,
    /// The FASTA a CRAM input was made with (bedtools --cram-ref), with its
    /// .fai (and .gzi when BGZF). Not needed for a CRAM written without a
    /// reference.
    #[arg(long, value_name = "FILE")]
    pub reference: Option<PathBuf>,
}

impl AlignArgs {
    /// The reader options; `full_records` keeps whole records (for `view` and `sort`).
    pub fn options(
        &self,
        full_records: bool,
    ) -> Result<gloci_io::AlignOptions, gloci_io::UsageError> {
        let flags = |name: &str, v: &Option<String>| match v {
            None => Ok(0),
            Some(s) => {
                gloci_io::parse_flags(s).map_err(|e| gloci_io::UsageError(format!("{name}: {e}")))
            }
        };
        Ok(gloci_io::AlignOptions {
            min_mapq: self.min_mapq,
            require_flags: flags("--require-flags", &self.require_flags)?,
            exclude_flags: flags("--exclude-flags", &self.exclude_flags)?,
            full_records,
            reference: self.reference.clone(),
            score_tag: None,
        })
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
    /// per-chromosome runs on indexed and plain files; ignored when nothing spills. It
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

/// `--region` help in a relation command, where only A is restricted
/// (applied with `mut_arg` on each relation command's args).
const RELATION_REGION_HELP: &str = "Restrict A to a region such as chr1:1,000-2,000 (repeatable). Only A is restricted, to the records that overlap a region (a point on the edge with --inclusive-points); B is read in full, so an A record reaching out of a region still meets B there.";

#[derive(Args, Debug, Clone)]
pub struct RegionArgs {
    /// Restrict to a region such as chr1:1,000-2,000 (repeatable).
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
    /// Input format, overriding detection: bed, bedgraph, wig, gff, gtf, vcf
    /// or bcf; bam, cram or sam where the command reads alignments; psl or
    /// rmsk for convert.
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
    /// A's format, overriding detection: bed, bedgraph, wig, gff, gtf, vcf or
    /// bcf; bam, cram or sam where the command reads alignments.
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
    /// B's format, overriding detection: bed, bedgraph, wig, gff, gtf, vcf or
    /// bcf; bam, cram or sam where the command reads alignments.
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
    /// Compare BED12 blocks (exons) and the aligned segments of BAM/SAM reads (cut at
    /// CIGAR N) instead of whole spans (4.3: BED12 blocks).
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

// P3b Task 1
/// `--seed` for the commands that draw random numbers (D143).
#[derive(Args, Debug, Clone, Copy)]
pub struct SeedArg {
    /// Seed for the random stream, so a run can be reproduced (default: a fresh seed, which is
    /// printed as a warning).
    #[arg(long, value_name = "N")]
    pub seed: Option<u64>,
}

// P3b Task 2
#[derive(Args, Debug, Clone)]
pub struct StatsArgs {
    #[command(subcommand)]
    pub command: StatsCommand,
}

#[derive(Subcommand, Debug, Clone)]
pub enum StatsCommand {
    /// The Jaccard similarity of two interval sets (bedtools jaccard).
    ///
    /// Each input is first merged as `gloci merge` does (overlapping and
    /// bookended intervals join; per strand with -s or -S). The intersection
    /// is the sum of the overlap bases of every pair (merged A, merged B) that
    /// passes -f, -F, -r, -e, --min-bases, -s, -S and --inclusive-points
    /// (measured on the merged intervals), the union is the bases covered by
    /// either input, and jaccard = intersection / union divided in single
    /// precision and printed as %g (six significant digits; `nan` when both
    /// are empty). n_intersections counts the accepted pairs. Output is a
    /// header line and one row: intersection, union, jaccard,
    /// n_intersections. A point (zero-length interval) covers no bases, so
    /// it adds nothing to the intersection or the union, but a point inside a
    /// merged interval of the other input counts as an intersection (bedtools
    /// widens a point to two bases; see docs/differences.md). Inputs must be
    /// sorted (--sort sorts them; --unsorted does not apply); a VCF is read
    /// as its footprints (--vcf-span); BAM, CRAM and SAM reads are read as
    /// their reference spans.
    #[command(after_long_help = "\
Examples:
  gloci stats jaccard -a peaks.bed -b genes.bed
  gloci stats jaccard -a rep1.bed -b rep2.bed -f 0.5 -r
  gloci stats jaccard -a calls.vcf.gz -b regions.bed --vcf-span ref
  gloci stats jaccard -a unsorted-a.bed -b unsorted-b.bed --sort")]
    Jaccard(JaccardArgs),
    // P3b Task 3
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
    #[command(after_long_help = "\
Examples:
  gloci stats fisher -a peaks.bed -b genes.bed -g hg38
  gloci stats fisher -a peaks.bed -b genes.bed -g sizes.genome --merge-overlaps
  gloci stats fisher -a peaks.bed -b genes.bed -g hg38 --exclude gaps.bed -f 0.5
  gloci stats fisher -a unsorted-a.bed -b unsorted-b.bed -g hg38 --sort")]
    Fisher(FisherArgs),
    // P3b Task 4
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
    #[command(after_long_help = "\
Examples:
  gloci stats reldist -a peaks.bed -b genes.bed
  gloci stats reldist -a peaks.bed -b genes.bed --detail
  gloci stats reldist -a unsorted-a.bed -b unsorted-b.bed --sort
  gloci stats reldist -a peaks.bed -b genes.bed --unsorted")]
    Reldist(ReldistArgs),
}

#[derive(Args, Debug, Clone)]
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

// P3b Task 3
#[derive(Args, Debug, Clone)]
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
    /// Genome with chromosome lengths: a built-in genome such as hg38 (see gloci genomes), or a
    /// .genome, .fai, .dict, VCF/BCF or BAM header file. Its lengths add up to the genome size.
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

// P3b Task 7
#[derive(Args, Debug, Clone)]
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

// P3b Task 8
#[derive(Args, Debug, Clone)]
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

// P3b Task 4
#[derive(Args, Debug, Clone)]
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
    fn align_args_options() {
        let a = AlignArgs {
            min_mapq: 20,
            require_flags: Some("paired".into()),
            exclude_flags: Some("0x904".into()),
            reference: None,
        };
        let o = a.options(true).unwrap();
        assert_eq!(
            (o.min_mapq, o.require_flags, o.exclude_flags),
            (20, 1, 0x904)
        );
        assert!(o.full_records);
        let bad = AlignArgs {
            exclude_flags: Some("bogus".into()),
            ..AlignArgs::default()
        };
        assert!(
            bad.options(false)
                .unwrap_err()
                .0
                .contains("--exclude-flags")
        );
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

// P3b Task 20

/// Which multi-allelic records `vcf split` splits.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum SplitTypesArg {
    Snps,
    Indels,
    Both,
    Any,
}

#[derive(Args, Debug, Clone)]
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

#[derive(Args, Debug, Clone)]
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

// P3b Task 22

#[derive(Args, Debug, Clone)]
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

#[derive(Args, Debug, Clone)]
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
// P3b Task 5
#[derive(Args, Debug, Clone)]
pub struct RandomArgs {
    /// Chromosomes and lengths: a built-in genome such as hg38 (see gloci genomes), or a
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

// P3b Task 6
#[derive(Args, Debug, Clone)]
pub struct ShuffleArgs {
    /// Input file (`-` for stdin).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Chromosomes and lengths: a built-in genome such as hg38 (see gloci genomes), or a
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

#[derive(Args, Debug, Clone)]
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
