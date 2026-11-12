use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Bgzf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Bed,
    BedGraph,
    /// WIG (D63): input only; read as bedGraph-like records.
    Wig,
    Vcf,
    Bcf,
    Gff3,
    Gtf,
    /// BAM (D100): a BGZF stream that decodes to `BAM\1`.
    Bam,
    Cram,
    Sam,
    /// FASTA: sequence input for `getfasta`, `maskfasta` and `nuc`; never an interval input.
    Fasta,
    /// PSL (D114): read by `convert` only.
    Psl,
    /// RepeatMasker `.out` (D114): read by `convert` only.
    Rmsk,
    /// BEDPE (D146): paired ends; never detected, chosen with `--format bedpe`
    /// and read by `pair-intersect`, `pair-pair` and `shuffle --format bedpe`.
    Bedpe,
}

impl Format {
    /// Alignment records (BAM, CRAM, SAM): read only by the commands that take alignments.
    pub fn is_alignment(self) -> bool {
        matches!(self, Format::Bam | Format::Cram | Format::Sam)
    }

    /// Formats that only `convert` reads (PSL, RepeatMasker).
    pub fn is_convert_only(self) -> bool {
        matches!(self, Format::Psl | Format::Rmsk)
    }

    /// The name users know the format by, for messages.
    pub fn name(self) -> &'static str {
        match self {
            Format::Bed => "BED",
            Format::BedGraph => "bedGraph",
            Format::Wig => "WIG",
            Format::Vcf => "VCF",
            Format::Bcf => "BCF",
            Format::Gff3 => "GFF3",
            Format::Gtf => "GTF",
            Format::Bam => "BAM",
            Format::Cram => "CRAM",
            Format::Sam => "SAM",
            Format::Fasta => "FASTA",
            Format::Psl => "PSL",
            Format::Rmsk => "RepeatMasker",
            Format::Bedpe => "BEDPE",
        }
    }
}

impl FromStr for Format {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "bed" => Ok(Format::Bed),
            "bedgraph" => Ok(Format::BedGraph),
            "wig" => Ok(Format::Wig),
            "vcf" => Ok(Format::Vcf),
            "bcf" => Ok(Format::Bcf),
            "gff" | "gff3" => Ok(Format::Gff3),
            "gtf" => Ok(Format::Gtf),
            "bam" => Ok(Format::Bam),
            "cram" => Ok(Format::Cram),
            "sam" => Ok(Format::Sam),
            "psl" => Ok(Format::Psl),
            "rmsk" => Ok(Format::Rmsk),
            "bedpe" => Ok(Format::Bedpe),
            _ => Err(format!(
                "unknown format '{s}'; expected bed, bedgraph, wig, vcf, bcf, gff, gff3, gtf, \
bam, cram, sam, psl, rmsk or bedpe"
            )),
        }
    }
}

/// Classify raw bytes: BGZF is gzip with a FEXTRA `BC` subfield.
pub fn detect_compression(head: &[u8]) -> Compression {
    if head.len() < 2 || head[0] != 0x1f || head[1] != 0x8b {
        return Compression::None;
    }
    if head.len() >= 12 && head[3] & 4 != 0 {
        let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
        let extra = &head[12..head.len().min(12 + xlen)];
        let mut i = 0;
        while i + 4 <= extra.len() {
            let len = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            if extra[i] == b'B' && extra[i + 1] == b'C' && len == 2 {
                return Compression::Bgzf;
            }
            i += 4 + len;
        }
    }
    Compression::Gzip
}

fn is_int(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

fn is_signed_int(s: &str) -> bool {
    is_int(s.strip_prefix('-').unwrap_or(s))
}

fn is_float(s: &str) -> bool {
    s.parse::<f64>().is_ok() && s.bytes().all(|b| b.is_ascii_digit() || b == b'.')
}

/// A SAM header line: `@HD`, `@SQ`, `@RG`, `@PG` or `@CO`, then a TAB.
fn is_sam_header(line: &str) -> bool {
    let b = line.as_bytes();
    b.len() >= 4
        && b[0] == b'@'
        && b[3] == b'\t'
        && matches!(&b[1..3], b"HD" | b"SQ" | b"RG" | b"PG" | b"CO")
}

/// `*` or a CIGAR string: one or more `<count><op>` pairs.
fn is_cigar(s: &str) -> bool {
    if s == "*" {
        return true;
    }
    let mut digits = 0;
    for b in s.bytes() {
        if b.is_ascii_digit() {
            digits += 1;
        } else if b"MIDNSHP=X".contains(&b) && digits > 0 {
            digits = 0;
        } else {
            return false;
        }
    }
    !s.is_empty() && digits == 0
}

/// `*`, `=` or a reference name (SAM spec: `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!-~]*`).
fn is_rnext(s: &str) -> bool {
    let b = s.as_bytes();
    s == "*"
        || s == "="
        || (!b.is_empty()
            && (b[0].is_ascii_alphanumeric() || b"!#$%&+./:;?@^_|~-".contains(&b[0]))
            && b.iter().all(|c| (b'!'..=b'~').contains(c)))
}

/// A SAM record: at least 11 TAB fields with a numeric FLAG (<= 65535), POS and MAPQ (<= 255),
/// a CIGAR, a valid RNEXT, a non-negative PNEXT, a signed TLEN, and a SEQ and QUAL of the
/// right alphabets. A BED12 line with a numeric name fails on its CIGAR (`+`/`-`/`.`).
fn is_sam_record(cols: &[&str]) -> bool {
    cols.len() >= 11
        && is_int(cols[1])
        && cols[1].parse::<u32>().is_ok_and(|v| v <= 65535)
        && is_int(cols[3])
        && is_int(cols[4])
        && cols[4].parse::<u32>().is_ok_and(|v| v <= 255)
        && is_cigar(cols[5])
        && is_rnext(cols[6])
        && is_int(cols[7])
        && is_signed_int(cols[8])
        && (cols[9] == "*"
            || cols[9]
                .bytes()
                .all(|b| b.is_ascii_alphabetic() || b == b'=' || b == b'.'))
        && !cols[9].is_empty()
        && !cols[10].is_empty()
        && cols[10].bytes().all(|b| (b'!'..=b'~').contains(&b))
}

/// A comma-separated list of integers that ends in a comma (`51,30,`).
fn is_block_list(s: &str) -> bool {
    s.strip_suffix(',')
        .is_some_and(|l| l.split(',').all(is_int))
}

/// A PSL row: 21 or more columns, integer counts, a strand, names, sizes that
/// hold the coordinates, and block lists that end in a comma. This is tight
/// on purpose: a BED with 21 columns and an integer chromosome must not pass.
fn is_psl_record(cols: &[&str]) -> bool {
    cols.len() >= 21
        && [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 14, 15, 16, 17]
            .iter()
            .all(|&i| is_int(cols[i]))
        && matches!(cols[8], "+" | "-" | "++" | "+-" | "-+" | "--")
        && !cols[9].is_empty()
        && !cols[13].is_empty()
        && cols[10].parse::<u64>().ok() >= cols[12].parse::<u64>().ok()
        && cols[14].parse::<u64>().ok() >= cols[16].parse::<u64>().ok()
        && cols[18..21].iter().all(|c| is_block_list(c))
}

/// `psLayout` as a whole word (a chromosome called `psLayoutX` is not a header).
fn is_psl_header(line: &str) -> bool {
    line.strip_prefix("psLayout")
        .is_some_and(|r| r.is_empty() || r.starts_with(|c: char| c.is_whitespace()))
}

/// A RepeatMasker `.out` header line (space separated).
fn is_rmsk_header(line: &str) -> bool {
    let t = line.trim_start();
    (t.starts_with("SW ") && t.contains("perc")) || (t.starts_with("score ") && t.contains("div."))
}

/// A RepeatMasker `.out` data line: score, three percentages, name, begin, end, `(left)`, `+`/`C`.
fn is_rmsk_record(line: &str) -> bool {
    let f: Vec<&str> = line.split_whitespace().collect();
    f.len() >= 14
        && is_int(f[0])
        && f[1..4].iter().all(|x| is_float(x))
        && is_int(f[5])
        && is_int(f[6])
        && f[7].starts_with('(')
        && f[7].ends_with(')')
        && matches!(f[8], "+" | "C")
}

/// Detect the format from decompressed leading bytes.
pub fn detect_format(head: &[u8]) -> Result<Format, String> {
    // CRAM is checked on the raw bytes: it is not compressed as a whole.
    // A CRAM file definition is `CRAM`, major, minor and a 20-byte file id (26 bytes); the
    // major version is a small byte, so a BED line on a chromosome named `CRAM...` is not CRAM.
    if head.len() >= 26 && head.starts_with(b"CRAM") && (2..=4).contains(&head[4]) {
        return match head[4] {
            2 | 3 => Ok(Format::Cram),
            v => Err(format!(
                "{{path}} is CRAM {v}, which gloci does not read; it reads CRAM 3.0 and 3.1"
            )),
        };
    }
    if head.starts_with(b"BAM\x01") {
        return Ok(Format::Bam);
    }
    // FASTA may start with blank lines, which samtools skips.
    let first = head.iter().position(|&b| b != b'\n' && b != b'\r');
    if first.is_some_and(|i| head[i] == b'>') {
        return Ok(Format::Fasta);
    }
    if head.starts_with(b"BCF\x02") {
        return Ok(Format::Bcf);
    }
    if head.starts_with(b"##fileformat=VCF") {
        return Ok(Format::Vcf);
    }
    if head.starts_with(b"##gff-version 3") {
        return Ok(Format::Gff3);
    }
    if head.starts_with(b"##gff-version 2") {
        return Ok(Format::Gtf);
    }
    let text = String::from_utf8_lossy(head);
    let mut bedgraph = false;
    for raw in text.split('\n') {
        let line = raw.strip_suffix('\r').unwrap_or(raw);
        if line.trim().is_empty() || line.starts_with('#') || line.starts_with("browser") {
            continue;
        }
        if is_sam_header(line) {
            return Ok(Format::Sam);
        }
        if is_psl_header(line) {
            return Ok(Format::Psl);
        }
        if is_rmsk_header(line) || is_rmsk_record(line) {
            return Ok(Format::Rmsk);
        }
        if line.starts_with("track") {
            if line.contains("type=bedGraph") {
                bedgraph = true;
            }
            continue;
        }
        if crate::wig::is_step_line(line.as_bytes()) {
            return Ok(Format::Wig);
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if is_sam_record(&cols) {
            return Ok(Format::Sam);
        }
        if is_psl_record(&cols) {
            return Ok(Format::Psl);
        }
        if cols.len() >= 9 && is_int(cols[3]) && is_int(cols[4]) && !is_int(cols[1]) {
            return Ok(if cols[8].contains('"') {
                Format::Gtf
            } else {
                Format::Gff3
            });
        }
        if cols.len() >= 3 && is_int(cols[1]) && is_int(cols[2]) {
            return Ok(if bedgraph {
                Format::BedGraph
            } else {
                Format::Bed
            });
        }
        let shown: String = line.chars().take(80).collect();
        return Err(format!(
            "could not recognise the format of {{path}} (first data line: {shown}); pass --format"
        ));
    }
    Ok(Format::Bed)
}
