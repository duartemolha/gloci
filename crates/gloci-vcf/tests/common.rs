//! Shared helpers for the integration tests.
#![allow(dead_code)]
use gloci_io::BcfDicts;
use std::io::Read;

pub const HEADER: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"d\">\n\
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"n\">\n\
##FILTER=<ID=LowQual,Description=\"l\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n\
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"f\">\n";

pub fn vcf(samples: &[&str], records: &[&str]) -> String {
    let mut s = String::from(HEADER);
    s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for n in samples {
        s.push('\t');
        s.push_str(n);
    }
    s.push('\n');
    for r in records {
        s.push_str(r);
        s.push('\n');
    }
    s
}

/// Dictionaries and framed records of BGZF BCF bytes.
pub fn bcf_records(bcf: &[u8]) -> (BcfDicts, Vec<Vec<u8>>) {
    let mut r = noodles_bgzf::io::Reader::new(bcf);
    let h = gloci_io::read_bcf_header(&mut r).expect("header");
    let mut rest = Vec::new();
    r.read_to_end(&mut rest).unwrap();
    let mut recs = Vec::new();
    let mut cur = &rest[..];
    while !cur.is_empty() {
        let ls = u32::from_le_bytes(cur[0..4].try_into().unwrap()) as usize;
        let li = u32::from_le_bytes(cur[4..8].try_into().unwrap()) as usize;
        recs.push(cur[..8 + ls + li].to_vec());
        cur = &cur[8 + ls + li..];
    }
    (h.bcf.expect("dicts"), recs)
}
