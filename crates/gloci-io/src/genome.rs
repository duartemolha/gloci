//! `--genome` loading and chromosome-table resolution (spec §4.5, D7, D27).

use crate::bcf::read_bcf_header;
use crate::builtin::resolve_genome;
use crate::error::{IoError, Loc};
use crate::record::{Header, RecordReader};
use crate::vcf::VcfHeader;
use gloci_core::{ChromOrder, ChromTable, Genome, GenomeError, OrderSource};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Load a chromosome list, detecting the kind by content: SAM `.dict`,
/// VCF header (plain, gzip or BGZF), BCF, BAM, CRAM, else `.fai` / `.genome`.
pub fn load_genome(path: &Path) -> Result<Genome, IoError> {
    let name = path.display().to_string();
    let read_err = |source| IoError::Read {
        path: name.clone(),
        source,
    };
    let f = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut raw = BufReader::new(f);
    let gz = raw.fill_buf().map_err(read_err)?.starts_with(&[0x1f, 0x8b]);
    let mut src: Box<dyn BufRead> = if gz {
        Box::new(BufReader::new(flate2::bufread::MultiGzDecoder::new(raw)))
    } else {
        Box::new(raw)
    };
    let magic = read_up_to(&mut src, 4).map_err(read_err)?;
    let mut src: Box<dyn BufRead> = Box::new(std::io::Cursor::new(magic.clone()).chain(src));
    let parse_err = |loc, msg: String| IoError::Parse {
        path: name.clone(),
        loc,
        msg,
    };
    let from_genome_err = |e: GenomeError| match e {
        GenomeError::Line { line, msg } => parse_err(Loc::Line(line as u64), msg),
        e => IoError::Detect {
            path: name.clone(),
            msg: format!("{name}: {e}"),
        },
    };

    if magic.starts_with(b"BCF\x02") {
        let h = read_bcf_header(&mut src).map_err(|msg| IoError::Detect {
            path: name.clone(),
            msg: format!("{name}: {msg}"),
        })?;
        return header_genome(&name, &h);
    }
    if magic == b"BAM\x01" {
        return bam_genome(&mut *src, &name);
    }
    if magic == b"CRAM" {
        // The @SQ lines of its header; no reference is needed (D103).
        let h = crate::align::read_cram_header(&mut src, &name)?;
        return Genome::new(h.refs.into_iter().map(|(n, l)| (n, Some(l))).collect())
            .map_err(from_genome_err);
    }
    if magic.starts_with(b">") {
        // A FASTA: its chromosome sizes are in its .fai, which must exist (D108).
        return fasta_genome(path, &name);
    }
    // Text: the first non-empty line decides.
    let mut first = Vec::new();
    while first.is_empty() || first == b"\n" || first == b"\r\n" {
        first.clear();
        if src.read_until(b'\n', &mut first).map_err(read_err)? == 0 {
            break;
        }
    }
    if first.starts_with(b"##fileformat=VCF") {
        let mut raw = first;
        loop {
            if raw.rsplit(|&b| b == b'\n').nth(1).is_some_and(|l| {
                let l = l.strip_suffix(b"\r").unwrap_or(l);
                l == b"#CHROM" || l.starts_with(b"#CHROM\t")
            }) {
                break;
            }
            if src.read_until(b'\n', &mut raw).map_err(read_err)? == 0 {
                break;
            }
        }
        let h = VcfHeader::parse(raw).map_err(|msg| parse_err(Loc::Line(1), msg))?;
        return header_genome(&name, &h);
    }
    let mut text = first;
    src.read_to_end(&mut text).map_err(read_err)?;
    if text.starts_with(b"@HD") || text.starts_with(b"@SQ") {
        Genome::parse_dict(&text)
    } else {
        Genome::parse_sizes(&text)
    }
    .map_err(from_genome_err)
}

/// The chromosome list of a FASTA, from `PATH.fai`. gloci never builds the
/// index implicitly; a missing one is a usage error naming `gloci index`.
fn fasta_genome(fasta: &Path, name: &str) -> Result<Genome, IoError> {
    let (fai, stale) = crate::fasta::find_fai(fasta)?;
    if stale {
        log::warn!(
            "the index {} is older than {name}; it may be out of date (rebuild it with: gloci index {name})",
            fai.display()
        );
    }
    load_genome(&fai)
}

/// Read up to `n` bytes (fewer only at end of input).
fn read_up_to(src: &mut dyn Read, n: usize) -> std::io::Result<Vec<u8>> {
    let mut v = vec![0u8; n];
    let mut got = 0;
    while got < n {
        match src.read(&mut v[got..]) {
            Ok(0) => break,
            Ok(k) => got += k,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    v.truncate(got);
    Ok(v)
}

fn bam_genome(src: &mut dyn BufRead, name: &str) -> Result<Genome, IoError> {
    let bad = |msg: &str| IoError::Detect {
        path: name.to_string(),
        msg: format!("{name}: {msg}"),
    };
    let io = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let i32le = |src: &mut dyn BufRead| -> Result<i32, IoError> {
        let mut b = [0u8; 4];
        src.read_exact(&mut b)
            .map_err(|_| bad("BAM header is truncated"))?;
        Ok(i32::from_le_bytes(b))
    };
    let mut magic = [0u8; 4];
    src.read_exact(&mut magic).map_err(io)?;
    let l_text = i32le(src)?;
    let l_text = u64::try_from(l_text).map_err(|_| bad("BAM header has a negative text length"))?;
    std::io::copy(&mut src.take(l_text), &mut std::io::sink()).map_err(io)?;
    let n_ref = i32le(src)?;
    let n_ref = usize::try_from(n_ref).map_err(|_| bad("BAM header has a negative n_ref"))?;
    let mut entries = Vec::new();
    for _ in 0..n_ref {
        let l_name = i32le(src)?;
        let l_name =
            usize::try_from(l_name).map_err(|_| bad("BAM reference name length is negative"))?;
        if l_name == 0 || l_name > 1 << 20 {
            return Err(bad("BAM reference name length is invalid"));
        }
        let mut n = vec![0u8; l_name];
        src.read_exact(&mut n)
            .map_err(|_| bad("BAM header is truncated"))?;
        n.pop(); // trailing NUL
        let l_ref = i32le(src)?;
        let l_ref = u64::try_from(l_ref).map_err(|_| bad("BAM reference length is negative"))?;
        entries.push((n, Some(l_ref)));
    }
    Genome::new(entries).map_err(|e| bad(&e.to_string()))
}

/// The chromosome order a VCF/BCF header declares. BCF uses the dictionary
/// order (matching record order and CSI ids), text VCF the `##contig` order.
fn header_entries(h: &VcfHeader) -> Vec<(Vec<u8>, Option<u64>)> {
    match &h.bcf {
        Some(d) => d
            .contigs
            .iter()
            .filter(|n| !n.is_empty())
            .map(|n| {
                let len = h.contigs.iter().find(|(c, _)| c == n).and_then(|c| c.1);
                (n.clone(), len)
            })
            .collect(),
        None => h.contigs.clone(),
    }
}

fn header_genome(name: &str, h: &VcfHeader) -> Result<Genome, IoError> {
    let entries = header_entries(h);
    if entries.is_empty() {
        return Err(IoError::Detect {
            path: name.to_string(),
            msg: format!(
                "{name} declares no ##contig lines, so it cannot give a chromosome order; pass a .fai, .genome or .dict file"
            ),
        });
    }
    Genome::new(entries).map_err(|e| IoError::Detect {
        path: name.to_string(),
        msg: format!("{name}: {e}"),
    })
}

/// Choose the chromosome table: `--genome` (the raw `-g` value: a built-in
/// genome id or a file), else the first VCF/BCF/alignment header with contigs, else
/// `fallback` (karyotype by default, spec §4.5).
pub fn chrom_table_for(
    genome: Option<&str>,
    inputs: &[&dyn RecordReader],
    fallback: ChromOrder,
) -> Result<ChromTable, IoError> {
    chrom_table(genome, inputs, fallback, false)
}

/// Like [`chrom_table_for`], for commands that need no chromosome order
/// beyond contiguity (`view`): a record on a contig the VCF header does not
/// declare is accepted, with a warning, as htslib does. `-g` stays strict.
pub fn chrom_table_lenient(
    genome: Option<&str>,
    inputs: &[&dyn RecordReader],
    fallback: ChromOrder,
) -> Result<ChromTable, IoError> {
    chrom_table(genome, inputs, fallback, true)
}

fn chrom_table(
    genome: Option<&str>,
    inputs: &[&dyn RecordReader],
    fallback: ChromOrder,
    lenient: bool,
) -> Result<ChromTable, IoError> {
    if let Some(arg) = genome {
        let (g, source) = resolve_genome(arg)?;
        return Ok(ChromTable::closed(&g, source));
    }
    for r in inputs {
        if let Header::Alignment(h) = r.header() {
            if !h.refs.is_empty() {
                let g = Genome::new(h.refs.iter().map(|(n, l)| (n.clone(), Some(*l))).collect())
                    .map_err(|e| IoError::Detect {
                        path: r.name().to_string(),
                        msg: format!("{}: {e}", r.name()),
                    })?;
                let source = OrderSource::AlignmentHeader(r.name().to_string());
                return Ok(if lenient {
                    ChromTable::declared(&g, source)
                } else {
                    ChromTable::closed(&g, source)
                });
            }
            continue;
        }
        if let Header::Vcf(h) = r.header() {
            let entries = header_entries(h);
            if !entries.is_empty() {
                let g = header_genome(r.name(), h)?;
                let source = OrderSource::VcfHeader(r.name().to_string());
                return Ok(if lenient {
                    ChromTable::declared(&g, source)
                } else {
                    ChromTable::closed(&g, source)
                });
            }
        }
    }
    Ok(match fallback {
        ChromOrder::Karyotype => ChromTable::karyotype(),
        ChromOrder::Lex => ChromTable::lexicographic(),
    })
}
