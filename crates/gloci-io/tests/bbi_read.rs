//! bigWig and bigBed detection (D171).

use gloci_io::{BIGBED_MAGIC, BIGWIG_MAGIC, Format, bbi_format};

#[test]
fn bbi_magic_both_byte_orders() {
    assert_eq!(bbi_format(&[0x26, 0xFC, 0x8F, 0x88]), Some(Format::BigWig));
    assert_eq!(bbi_format(&[0x88, 0x8F, 0xFC, 0x26]), Some(Format::BigWig));
    assert_eq!(bbi_format(&[0xEB, 0xF2, 0x89, 0x87]), Some(Format::BigBed));
    assert_eq!(bbi_format(&[0x87, 0x89, 0xF2, 0xEB]), Some(Format::BigBed));
    assert_eq!(bbi_format(b"chr1\t0\t1\n"), None);
    assert_eq!(bbi_format(&[0x26, 0xFC, 0x8F]), None);
    assert_eq!(BIGWIG_MAGIC, 0x888F_FC26);
    assert_eq!(BIGBED_MAGIC, 0x8789_F2EB);
    assert!(Format::BigWig.is_bbi() && Format::BigBed.is_bbi() && !Format::Bed.is_bbi());
    assert_eq!(Format::BigWig.name(), "bigWig");
    assert_eq!("bigbed".parse::<Format>(), Ok(Format::BigBed));
}

#[test]
fn chrom_sizes_in_the_files_stored_order() {
    let bw = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/bigwig/pybigwig.bw");
    let sizes = gloci_io::bbi_chrom_sizes(bw.to_str().unwrap()).unwrap();
    let shown: Vec<(String, u64)> = sizes
        .into_iter()
        .map(|(n, l)| (String::from_utf8(n).unwrap(), l))
        .collect();
    assert_eq!(
        shown,
        [
            ("chr1".to_string(), 5000),
            ("chr2".to_string(), 3000),
            ("chr10".to_string(), 2000),
            ("chrM".to_string(), 100),
        ]
    );
}
