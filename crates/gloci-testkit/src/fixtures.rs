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

/// htslib's `hts_reg2bin` for the 0-based region `[beg, end)`.
fn reg2bin(beg: u64, end: u64, min_shift: u32, depth: u32) -> usize {
    let end = end.max(beg + 1) - 1;
    let (mut s, mut t) = (min_shift, ((1u64 << (3 * depth)) - 1) / 7);
    let mut l = depth;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as usize;
        }
        l -= 1;
        s += 3;
        t -= 1 << (3 * l);
    }
    0
}

/// A CSI index (min_shift 14, depth 5) for a BGZF BCF, built the way
/// htslib's `bcftools index` builds one: it accepts contigs in any order as
/// long as each one's records are contiguous, so reference ids need not
/// rise through the file. noodles' indexer (and gloci's) refuse such files,
/// so the index is assembled from its parts. Record spans are the REF span.
pub fn bcf_csi_any_order(bcf: &[u8]) -> Vec<u8> {
    use noodles_csi::binning_index::index::ReferenceSequence;
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};
    use noodles_vcf::variant::record::ReferenceBases as _;
    use std::collections::BTreeMap;

    const MIN_SHIFT: u8 = 14;
    const DEPTH: u8 = 5;
    let mut reader = noodles_bcf::io::Reader::new(bcf);
    let header = reader.read_header().expect("fixture BCF header");
    let n_refs = header.contigs().len();
    // rid -> (bin -> chunks, record count)
    type Refs = BTreeMap<usize, (BTreeMap<usize, Vec<Chunk>>, u64)>;
    let mut refs: Refs = BTreeMap::new();
    let mut spans = Vec::new();
    let mut rec = noodles_bcf::Record::default();
    loop {
        let start = reader.get_ref().virtual_position();
        if reader.read_record(&mut rec).expect("fixture BCF record") == 0 {
            break;
        }
        let end = reader.get_ref().virtual_position();
        let rid = rec.reference_sequence_id().expect("CHROM");
        let pos = usize::from(rec.variant_start().expect("POS").expect("POS")) as u64 - 1;
        let ref_len = rec.reference_bases().len() as u64;
        let bin = reg2bin(pos, pos + ref_len.max(1), MIN_SHIFT.into(), DEPTH.into());
        let entry = refs.entry(rid).or_default();
        let chunks = entry.0.entry(bin).or_default();
        match chunks.last_mut() {
            Some(c) if c.end() == start => *c = Chunk::new(c.start(), end),
            _ => chunks.push(Chunk::new(start, end)),
        }
        entry.1 += 1;
        spans.push((rid, start, end));
    }
    let reference_sequences = (0..n_refs)
        .map(|rid| {
            let Some((bins, count)) = refs.remove(&rid) else {
                return ReferenceSequence::new(Default::default(), Default::default(), None);
            };
            let first = spans.iter().find(|s| s.0 == rid).expect("a record").1;
            let last = spans.iter().rev().find(|s| s.0 == rid).expect("a record").2;
            let loffsets = bins
                .iter()
                .map(|(&id, chunks)| (id, chunks[0].start()))
                .collect();
            let bins = bins
                .into_iter()
                .map(|(id, chunks)| (id, Bin::new(chunks)))
                .collect();
            ReferenceSequence::new(bins, loffsets, Some(Metadata::new(first, last, count, 0)))
        })
        .collect();
    let index = noodles_csi::Index::builder()
        .set_min_shift(MIN_SHIFT)
        .set_depth(DEPTH)
        .set_reference_sequences(reference_sequences)
        .build();
    let mut w = noodles_csi::io::Writer::new(Vec::new());
    w.write_index(&index).expect("write CSI");
    w.into_inner().finish().expect("finish CSI")
}
