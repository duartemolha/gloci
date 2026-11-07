//! Byte-level fixture builders shared by the test suites.

use std::io::Write;
use std::path::{Path, PathBuf};

/// BGZF-compress `data`; the output ends with the 28-byte EOF block.
pub fn bgzip(data: &[u8]) -> Vec<u8> {
    let mut w = noodles_bgzf::io::Writer::new(Vec::new());
    w.write_all(data).expect("write to Vec");
    w.finish().expect("finish bgzf")
}

/// Plain (non-BGZF) gzip.
pub fn gzip(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).expect("write to Vec");
    e.finish().expect("finish gzip")
}

/// Convert VCF text to BGZF-compressed BCF 2.2 with noodles, an implementation
/// independent of gloci's BCF code.
///
/// noodles drops `IDX=` when it re-serialises a header, although it encodes
/// records with the `IDX=` dictionaries. So the BCF header text is the
/// fixture's own header (as htslib writes it) and noodles encodes the records.
/// Panics if the dictionaries the header text declares differ from the ones
/// noodles encodes with, so a fixture can never be silently inconsistent.
pub fn vcf_to_bcf(vcf_text: &[u8]) -> Vec<u8> {
    use noodles_vcf::header::StringMaps;
    use noodles_vcf::variant::io::Write as _;

    let mut reader = noodles_vcf::io::Reader::new(vcf_text);
    let header = reader.read_header().expect("fixture VCF header");
    let text_len = header_text_len(vcf_text);
    let text = &vcf_text[..text_len];
    let declared: StringMaps = std::str::from_utf8(text)
        .expect("fixture header is UTF-8")
        .parse()
        .expect("fixture header dictionaries");
    let encoded = StringMaps::try_from(&header).expect("noodles dictionaries");
    assert_eq!(
        declared, encoded,
        "the fixture header's dictionaries differ from the ones noodles encodes with; \
         list INFO, then FILTER, then FORMAT lines or give every line an IDX="
    );

    let mut w = noodles_bcf::io::Writer::from(Vec::new());
    w.write_header(&header).expect("noodles BCF header");
    w.get_mut().clear(); // keep only the records; the header text is written below
    for rec in reader.records() {
        let rec = rec.expect("fixture VCF record");
        w.write_variant_record(&header, &rec)
            .expect("noodles BCF record");
    }
    let records = w.into_inner();

    let mut out = Vec::with_capacity(text_len + records.len() + 16);
    out.extend_from_slice(b"BCF\x02\x02");
    out.extend_from_slice(&(text_len as u32 + 1).to_le_bytes());
    out.extend_from_slice(text);
    out.push(0);
    out.extend_from_slice(&records);
    bgzip(&out)
}

/// Length of the header (every line up to and including `#CHROM`).
fn header_text_len(vcf: &[u8]) -> usize {
    let mut off = 0;
    for line in vcf.split_inclusive(|&b| b == b'\n') {
        off += line.len();
        if line.starts_with(b"#CHROM") {
            return off;
        }
    }
    panic!("fixture VCF has no #CHROM line");
}

/// Write `bytes` to `dir/name` and return the path.
pub fn write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, bytes).expect("write fixture");
    p
}

/// A BGZF-wrapped BAM header (`BAM\1`, `l_text = 0`) listing `refs` as name and length.
pub fn bam_header(refs: &[(&str, u32)]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"BAM\x01");
    out.extend_from_slice(&0i32.to_le_bytes());
    out.extend_from_slice(&(refs.len() as i32).to_le_bytes());
    for (name, len) in refs {
        out.extend_from_slice(&(name.len() as i32 + 1).to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.push(0);
        out.extend_from_slice(&(*len as i32).to_le_bytes());
    }
    bgzip(&out)
}
