//! Stable stand-in for the cargo-fuzz targets: arbitrary and structured bytes
//! (raw and BGZF-wrapped) through the same bodies must never panic.
//!
//! The first byte of every input is the config byte of `fuzz_bodies`.

use gloci_io::Format;
use gloci_io::fuzz_bodies::{drain_input, read_all_both, region_parse};
use gloci_testkit::fixtures::{bgzip, sam_to_bam, vcf_to_bcf};
use proptest::prelude::*;

pub const VCF_HEADER: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=248956422>\n\
##contig=<ID=chr2,length=242193529>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Len\">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n";

fn bytes() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(any::<u8>(), 0..512)
}

/// Text-ish input, so the parsers get past the first few bytes.
fn texty() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop::sample::select(b"chr1\t2 345-:#.,;=_\nACGT<>*/|\rHOMVCF".to_vec()),
        0..400,
    )
}

fn cfg() -> impl Strategy<Value = u8> {
    any::<u8>()
}

fn with_cfg(c: u8, mut d: Vec<u8>) -> Vec<u8> {
    d.insert(0, c);
    d
}

fn num() -> impl Strategy<Value = String> {
    prop_oneof![
        3 => (0u64..2_000_000).prop_map(|n| n.to_string()),
        1 => Just("0".to_string()),
        1 => Just("-5".to_string()),
        1 => Just("-2147483648".to_string()),
        1 => Just("18446744073709551615".to_string()),
        1 => Just("99999999999999999999999".to_string()),
        1 => Just("4294967296".to_string()),
        1 => Just("1e5".to_string()),
        1 => Just("".to_string()),
        1 => Just(".".to_string()),
        1 => Just("x".to_string()),
    ]
}

fn allele() -> impl Strategy<Value = String> {
    prop_oneof![
        3 => "[ACGTN]{1,6}",
        1 => Just("<DEL>".to_string()),
        1 => Just("<INS>".to_string()),
        1 => Just("<DUP:TANDEM>".to_string()),
        1 => Just("<*>".to_string()),
        1 => Just("*".to_string()),
        1 => Just(".".to_string()),
        1 => Just("".to_string()),
        1 => Just("A]chr1:100]".to_string()),
        1 => Just("]chr2:5]T".to_string()),
        1 => Just("G[chr1:7[".to_string()),
        1 => Just("C.".to_string()),
        1 => Just(".T".to_string()),
        1 => "[ -~]{0,8}",
    ]
}

fn info() -> impl Strategy<Value = String> {
    let item = prop_oneof![
        num().prop_map(|n| format!("END={n}")),
        num().prop_map(|n| format!("SVLEN={n}")),
        (num(), num()).prop_map(|(a, b)| format!("SVLEN={a},{b}")),
        Just("SVTYPE=DEL".to_string()),
        Just("END".to_string()),
        Just("SVLEN=".to_string()),
        "[A-Za-z0-9=,;.-]{0,10}",
    ];
    prop::collection::vec(item, 0..4).prop_map(|v| {
        if v.is_empty() {
            ".".into()
        } else {
            v.join(";")
        }
    })
}

fn vcf_line(samples: usize) -> impl Strategy<Value = String> {
    (
        prop::sample::select(vec!["chr1", "chr2", "chrUn", "1", ""]),
        num(),
        allele(),
        allele(),
        prop::collection::vec(allele(), 0..3),
        info(),
        prop::collection::vec(
            prop::sample::select(vec!["0/1", "1|1", ".", "2", "x"]),
            samples,
        ),
        0u8..12,
    )
        .prop_map(move |(chrom, pos, r, a, extra, info, gts, drop)| {
            let alts = std::iter::once(a)
                .chain(extra)
                .collect::<Vec<_>>()
                .join(",");
            let mut cols = vec![
                chrom.to_string(),
                pos,
                ".".to_string(),
                r,
                alts,
                ".".to_string(),
                "PASS".to_string(),
                info,
            ];
            if samples > 0 {
                cols.push("GT".to_string());
                cols.extend(gts.into_iter().map(String::from));
            }
            // Sometimes drop trailing columns.
            let keep = cols
                .len()
                .saturating_sub(usize::from(drop.saturating_sub(8)));
            cols.truncate(keep.max(1));
            cols.join("\t")
        })
}

fn vcf_text() -> impl Strategy<Value = Vec<u8>> {
    (0usize..=3).prop_flat_map(|n| {
        prop::collection::vec(vcf_line(n), 0..12).prop_map(move |lines| {
            let mut s = String::from(VCF_HEADER);
            s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
            if n > 0 {
                s.push_str("\tFORMAT");
                for i in 0..n {
                    s.push_str(&format!("\tS{i}"));
                }
            }
            s.push('\n');
            for l in lines {
                s.push_str(&l);
                s.push('\n');
            }
            s.into_bytes()
        })
    })
}

const BCF_SEED_VCF: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=248956422>\n\
##contig=<ID=chr2,length=242193529>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Len\">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tACGT\tA\t.\tPASS\t.\n\
chr1\t200\t.\tA\t<DEL>\t.\tPASS\tEND=300;SVLEN=-100;SVTYPE=DEL\n\
chr2\t5\t.\tG\tGT,*\t.\tPASS\t.\n\
chr2\t9\t.\tC\tC]chr1:5]\t.\tPASS\tSVTYPE=BND\n";

/// A valid BCF (BGZF-compressed) and its decompressed bytes' record offset.
fn bcf_seed() -> Vec<u8> {
    vcf_to_bcf(BCF_SEED_VCF.as_bytes())
}

fn bcf_plain() -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(&bcf_seed()[..])
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Byte-level mutations of a valid BCF (flip, truncate, insert), re-wrapped.
fn bcf_mutated() -> impl Strategy<Value = Vec<u8>> {
    let base = bcf_plain();
    let len = base.len();
    let header_end = 9 + u32::from_le_bytes(base[5..9].try_into().unwrap()) as usize;
    (
        prop::collection::vec((header_end..len, any::<u8>()), 0..6),
        prop::option::of(header_end..len),
        prop::option::of((header_end..len, bytes())),
        prop::collection::vec((0..len, any::<u8>()), 0..2),
    )
        .prop_map(move |(flips, trunc, ins, hdr_flips)| {
            let mut d = base.clone();
            for (i, b) in flips.into_iter().chain(hdr_flips) {
                d[i] = b;
            }
            if let Some((at, extra)) = ins {
                d.splice(at..at, extra);
            }
            if let Some(t) = trunc {
                d.truncate(t.min(d.len()));
            }
            bgzip(&d)
        })
}

/// Valid BCF header prefix followed by arbitrary record bytes.
fn bcf_header_then_junk() -> impl Strategy<Value = Vec<u8>> {
    let base = bcf_plain();
    let header_end = 9 + u32::from_le_bytes(base[5..9].try_into().unwrap()) as usize;
    bytes().prop_map(move |junk| {
        let mut d = base[..header_end].to_vec();
        d.extend(junk);
        bgzip(&d)
    })
}

const BAM_SEED_SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n\
a\t0\tchr1\t11\t60\t5S10M2D5M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\tNM:i:2\n\
b\t16\tchr1\t21\t30\t8M100N8M\t*\t0\t0\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\tNM:i:0\tXX:Z:hi\n\
c\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
d\t0\tchr2\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";

fn bam_plain() -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(&sam_to_bam(BAM_SEED_SAM)[..])
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Byte-level mutations of a valid BAM (flip, truncate, insert), re-wrapped.
fn bam_mutated() -> impl Strategy<Value = Vec<u8>> {
    let base = bam_plain();
    let len = base.len();
    (
        prop::collection::vec((0..len, any::<u8>()), 0..6),
        prop::option::of(0..len),
        prop::option::of((0..len, bytes())),
    )
        .prop_map(move |(flips, trunc, ins)| {
            let mut d = base.clone();
            for (i, b) in flips {
                d[i] = b;
            }
            if let Some((at, extra)) = ins {
                d.splice(at..at, extra);
            }
            if let Some(t) = trunc {
                d.truncate(t.min(d.len()));
            }
            bgzip(&d)
        })
}

fn bed_lines() -> impl Strategy<Value = Vec<u8>> {
    let field = prop_oneof![
        3 => num(),
        1 => Just("+".to_string()),
        1 => Just("-".to_string()),
        1 => Just("0,0,0".to_string()),
        1 => Just("2".to_string()),
        1 => Just("10,20".to_string()),
        1 => Just("0,50".to_string()),
        1 => "[ -~]{0,6}",
    ];
    let line = (
        prop::sample::select(vec!["chr1", "chr2", "track", "#c", "browser x", ""]),
        prop::collection::vec(field, 0..11),
    )
        .prop_map(|(c, f)| {
            let mut v = vec![c.to_string()];
            v.extend(f);
            v.join("\t")
        });
    prop::collection::vec(line, 0..10).prop_map(|l| l.join("\n").into_bytes())
}

fn gff_lines() -> impl Strategy<Value = Vec<u8>> {
    let attrs = prop_oneof![
        Just(".".to_string()),
        Just("ID=g1;Name=x".to_string()),
        Just("gene_id \"g1\"; transcript_id \"t\";".to_string()),
        "[ -~]{0,16}",
    ];
    let line = (
        prop::sample::select(vec!["chr1", "chr2", "##gff-version 3", "#c", ""]),
        prop::sample::select(vec![".", "src"]),
        prop::sample::select(vec!["gene", "exon", ".", "CDS"]),
        num(),
        num(),
        prop::sample::select(vec![".", "1.5", "x"]),
        prop::sample::select(vec!["+", "-", ".", "?"]),
        prop::sample::select(vec![".", "0", "1", "3"]),
        attrs,
        0usize..10,
    )
        .prop_map(|(a, b, c, d, e, f, g, h, i, keep)| {
            let mut v = vec![
                a.to_string(),
                b.into(),
                c.into(),
                d,
                e,
                f.into(),
                g.into(),
                h.into(),
                i,
            ];
            v.truncate(keep.max(1));
            v.join("\t")
        });
    prop::collection::vec(line, 0..10).prop_map(|l| l.join("\n").into_bytes())
}

fn all_formats(d: &[u8], c: u8) {
    for f in [
        Format::Bed,
        Format::Vcf,
        Format::Gff3,
        Format::Gtf,
        Format::Bcf,
    ] {
        read_all_both(&with_cfg(c, d.to_vec()), Some(f));
    }
    read_all_both(&with_cfg(c, d.to_vec()), None);
    read_all_both(&with_cfg(c, bgzip(d)), None);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn readers_never_panic(d in prop_oneof![bytes(), texty()], c in cfg()) {
        all_formats(&d, c);
        read_all_both(&with_cfg(c, bgzip(&d)), Some(Format::Bcf));
        read_all_both(&with_cfg(c, bgzip(&d)), Some(Format::Vcf));
    }

    #[test]
    fn structured_vcf(d in vcf_text(), c in cfg()) {
        read_all_both(&with_cfg(c, d.clone()), Some(Format::Vcf));
        read_all_both(&with_cfg(c, bgzip(&d)), Some(Format::Vcf));
        read_all_both(&with_cfg(c, d), None);
    }

    #[test]
    fn structured_vcf_converted_to_bcf_never_panics(d in vcf_text(), c in cfg()) {
        // Arbitrary record text mostly fails noodles' own parse; that is fine.
        let _ = std::panic::catch_unwind(|| vcf_to_bcf(&d)).map(|b| {
            read_all_both(&with_cfg(c, b), Some(Format::Bcf));
        });
    }

    #[test]
    fn mutated_bcf(d in bcf_mutated(), c in cfg()) {
        read_all_both(&with_cfg(c, d), Some(Format::Bcf));
    }

    #[test]
    fn mutated_bam(d in bam_mutated(), c in cfg()) {
        read_all_both(&with_cfg(c, d), Some(Format::Bam));
    }

    #[test]
    fn arbitrary_bytes_as_bam(d in prop_oneof![bytes(), texty()], c in cfg()) {
        read_all_both(&with_cfg(c, bgzip(&d)), Some(Format::Bam));
    }

    #[test]
    fn bcf_header_junk(d in bcf_header_then_junk(), c in cfg()) {
        read_all_both(&with_cfg(c, d), Some(Format::Bcf));
    }

    #[test]
    fn structured_bed(d in bed_lines(), c in cfg()) {
        read_all_both(&with_cfg(c, d.clone()), Some(Format::Bed));
        read_all_both(&with_cfg(c, d.clone()), Some(Format::BedGraph));
        read_all_both(&with_cfg(c, bgzip(&d)), None);
    }

    #[test]
    fn structured_gff(d in gff_lines(), c in cfg()) {
        read_all_both(&with_cfg(c, d.clone()), Some(Format::Gff3));
        read_all_both(&with_cfg(c, d), Some(Format::Gtf));
    }

    #[test]
    fn input_never_panics(d in prop_oneof![bytes(), texty()], c in cfg()) {
        drain_input(&with_cfg(c, d.clone()));
        drain_input(&with_cfg(c, bgzip(&d)));
    }

    #[test]
    fn region_parse_never_panics(d in prop_oneof![bytes(), texty()]) {
        region_parse(&d);
    }
}

/// Writes the seed corpus under `fuzz/seeds/`. Run with
/// `cargo test -p gloci-io --test fuzz_smoke -- --ignored regen_seeds`.
#[test]
#[ignore]
fn regen_seeds() {
    use std::fs::{create_dir_all, write};
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fuzz/seeds");
    let put = |target: &str, name: &str, cfg: u8, d: &[u8]| {
        let dir = root.join(target);
        create_dir_all(&dir).unwrap();
        let mut v = vec![cfg];
        v.extend_from_slice(d);
        write(dir.join(name), v).unwrap();
    };
    let vcf = format!(
        "{VCF_HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n\
chr1\t100\t.\tACGT\tA\t.\tPASS\t.\tGT\t0/1\n\
chr1\t200\t.\tA\t<DEL>\t.\tPASS\tEND=300;SVLEN=-100;SVTYPE=DEL\tGT\t1|1\n\
chr2\t5\t.\tG\tGT,*\t.\tPASS\t.\tGT\t.\n\
chr2\t9\t.\tC\tC]chr1:5]\t.\tPASS\tSVTYPE=BND\tGT\t0/1\n"
    );
    let bed = "chr1\t10\t20\n#c\nchr2\t5\t6\tname\t7\t+\n";
    let bed12 = "chr1\t100\t400\tg\t0\t+\t100\t400\t0,0,0\t2\t50,60\t0,240\n";
    let gff = "##gff-version 3\nchr1\tsrc\tgene\t101\t400\t.\t+\t.\tID=g1;Name=x\nchr1\tsrc\texon\t101\t200\t.\t+\t.\tParent=g1\n";
    let gtf = "chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n";
    put("vcf_reader", "vcf", 0, vcf.as_bytes());
    put("vcf_reader", "vcf_threads", 0b00110, vcf.as_bytes());
    put("bed_reader", "bed", 0, bed.as_bytes());
    put("bed_reader", "bed12", 0, bed12.as_bytes());
    put("gff_reader", "gff3", 0, gff.as_bytes());
    put("gff_reader", "gtf", 0, gtf.as_bytes());
    // Bit 0 clear: the target does not wrap, so the seed is already BGZF.
    put("bcf_reader", "bcf", 0, &bcf_seed());
    put("bgzf_input", "bgzf_vcf", 0, &bgzip(vcf.as_bytes()));
    put("bgzf_input", "plain_bed", 0, bed.as_bytes());
    for (n, r) in [
        "chr1",
        "chr1:100-200",
        "chr1:5-",
        "{chr1:2}:3-4",
        "chr1:1,000-2,000",
    ]
    .iter()
    .enumerate()
    {
        let dir = root.join("region_parse");
        create_dir_all(&dir).unwrap();
        write(dir.join(format!("r{n}")), r).unwrap();
    }
}
