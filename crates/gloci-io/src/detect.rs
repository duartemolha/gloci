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
    Vcf,
    Bcf,
    Gff3,
    Gtf,
}

impl Format {
    /// The name users know the format by, for messages.
    pub fn name(self) -> &'static str {
        match self {
            Format::Bed => "BED",
            Format::BedGraph => "bedGraph",
            Format::Vcf => "VCF",
            Format::Bcf => "BCF",
            Format::Gff3 => "GFF3",
            Format::Gtf => "GTF",
        }
    }
}

impl FromStr for Format {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "bed" => Ok(Format::Bed),
            "bedgraph" => Ok(Format::BedGraph),
            "vcf" => Ok(Format::Vcf),
            "bcf" => Ok(Format::Bcf),
            "gff" | "gff3" => Ok(Format::Gff3),
            "gtf" => Ok(Format::Gtf),
            _ => Err(format!(
                "unknown format '{s}'; expected bed, bedgraph, vcf, bcf, gff, gff3 or gtf"
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

/// Detect the format from decompressed leading bytes.
pub fn detect_format(head: &[u8]) -> Result<Format, String> {
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
        if line.starts_with("track") {
            if line.contains("type=bedGraph") {
                bedgraph = true;
            }
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
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
