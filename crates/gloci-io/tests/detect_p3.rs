//! Format detection (D100) and `parse_flags`.
use gloci_io::{Format, detect_format, parse_flags};
use gloci_testkit::fixtures;

const SAM_LINE: &[u8] = b"r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n";
const PSL_LINE: &[u8] =
    b"50\t1\t0\t0\t0\t0\t0\t0\t+\tq1\t51\t0\t51\tchr1\t1000\t100\t151\t1\t51,\t0,\t100,\n";
const RMSK_LINE: &[u8] = b"  1306  15.6  6.2  0.0  chr1  10001  10468 (248945954) + (TAATCC)n Simple_repeat 1 463 (0) 1\n";

fn detect(b: &[u8]) -> Result<Format, String> {
    detect_format(b)
}

#[test]
fn detect_p3_formats() {
    assert_eq!(detect(b"BAM\x01\0\0\0\0"), Ok(Format::Bam));
    // Through the reader the BGZF layer is removed first.
    let bam = fixtures::bgzip(b"BAM\x01\0\0\0\0");
    let mut dec = Vec::new();
    std::io::Read::read_to_end(&mut noodles_bgzf::io::Reader::new(&bam[..]), &mut dec).unwrap();
    assert_eq!(detect(&dec), Ok(Format::Bam));

    let cram = |major: u8, minor: u8| {
        let mut v = b"CRAM".to_vec();
        v.extend([major, minor]);
        v.extend([b'f'; 20]);
        v
    };
    assert_eq!(detect(&cram(3, 0)), Ok(Format::Cram));
    assert_eq!(detect(&cram(2, 1)), Ok(Format::Cram));
    assert!(detect(&cram(4, 0)).unwrap_err().contains("CRAM 4"));
    // A BED line on a chromosome that starts with CRAM is BED; a short prefix is not CRAM.
    assert_eq!(detect(b"CRAM\t1\t2\n"), Ok(Format::Bed));
    assert_eq!(detect(b"CRAMx\t1\t2\n"), Ok(Format::Bed));
    assert!(detect(b"CRAM\x03\x00").is_err());

    assert_eq!(detect(b"@HD\tVN:1.6\n"), Ok(Format::Sam));
    assert_eq!(detect(b"@SQ\tSN:chr1\tLN:5\n"), Ok(Format::Sam));
    assert_eq!(detect(SAM_LINE), Ok(Format::Sam));
    // A negative TLEN is still SAM.
    assert_eq!(
        detect(b"r1\t99\tchr1\t100\t60\t10M\t=\t200\t-110\tACGTACGTAC\t*\n"),
        Ok(Format::Sam)
    );

    assert_eq!(detect(b">chr1 d\nACGT\n"), Ok(Format::Fasta));
    let gz = fixtures::gzip(b">chr1 d\nACGT\n");
    let mut dec = Vec::new();
    std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(&gz[..]), &mut dec).unwrap();
    assert_eq!(detect(&dec), Ok(Format::Fasta));

    assert_eq!(
        detect(b"psLayout version 3\n\nmatch\tmis-\trep. \n"),
        Ok(Format::Psl)
    );
    assert_eq!(detect(PSL_LINE), Ok(Format::Psl));

    let header = b"   SW  perc perc perc  query      position in query    matching  repeat\n\
score  div. del.  ins.  sequence    begin     end    (left)   repeat  class/family\n\n";
    assert_eq!(detect(header), Ok(Format::Rmsk));
    assert_eq!(detect(RMSK_LINE), Ok(Format::Rmsk));

    assert_eq!(detect(b"chr1\t10\t20\n"), Ok(Format::Bed));
    assert_eq!(detect(b"##fileformat=VCFv4.2\n"), Ok(Format::Vcf));
}

#[test]
fn format_flag_names_and_helpers() {
    for (s, f) in [
        ("bam", Format::Bam),
        ("cram", Format::Cram),
        ("sam", Format::Sam),
        ("psl", Format::Psl),
        ("rmsk", Format::Rmsk),
    ] {
        assert_eq!(s.parse::<Format>(), Ok(f));
    }
    assert!(
        Format::Bam.is_alignment() && Format::Cram.is_alignment() && Format::Sam.is_alignment()
    );
    assert!(!Format::Bed.is_alignment() && !Format::Fasta.is_alignment());
    assert!(Format::Psl.is_convert_only() && Format::Rmsk.is_convert_only());
    assert!(!Format::Sam.is_convert_only());
}

#[test]
fn parse_flags_table() {
    assert_eq!(parse_flags("0x704"), Ok(0x704));
    assert_eq!(parse_flags("1796"), Ok(1796));
    assert_eq!(parse_flags("secondary,dup"), Ok(0x500));
    assert_eq!(parse_flags("UNMAP,MUNMAP"), Ok(0xC));
    assert_eq!(parse_flags("Supplementary"), Ok(0x800));
    assert!(parse_flags("").is_err());
    assert!(parse_flags("70000").is_err());
    assert!(parse_flags("bogus").unwrap_err().contains("PAIRED"));
}

#[test]
fn numeric_name_bed12_and_bed13_are_not_sam() {
    // Cols 2, 4, 5, 8 and 9 are integers here, as in a SAM record, but col 6 is a strand.
    assert_eq!(
        detect(b"chr1\t100\t200\t77\t5\t+\t100\t200\t0\t1\t100,\t0,\n"),
        Ok(Format::Bed)
    );
    assert_eq!(
        detect(b"chr1\t100\t200\t77\t5\t.\t100\t200\t0\t1\t100,\t0,\t3\n"),
        Ok(Format::Bed)
    );
    assert_eq!(
        detect(b"chr1\t100\t200\t77\t5\t-\t100\t200\t12\t2\t50,50,\t0,100,\n"),
        Ok(Format::Bed)
    );
}

#[test]
fn sam_field_checks() {
    let sam = |f: &[&str]| detect(format!("{}\n", f.join("\t")).as_bytes());
    let ok = [
        "r1",
        "0",
        "chr1",
        "100",
        "60",
        "10M",
        "*",
        "0",
        "0",
        "ACGTACGTAC",
        "*",
    ];
    assert_eq!(sam(&ok), Ok(Format::Sam));
    let with = |i: usize, v: &'static str| {
        let mut f = ok;
        f[i] = v;
        f
    };
    assert_eq!(sam(&with(5, "5S10M2I")), Ok(Format::Sam));
    assert_eq!(sam(&with(6, "=")), Ok(Format::Sam));
    assert_eq!(sam(&with(6, "chr2")), Ok(Format::Sam));
    assert_eq!(sam(&with(10, "IIIIIIIIII")), Ok(Format::Sam));
    assert_ne!(sam(&with(1, "70000")), Ok(Format::Sam));
    assert_ne!(sam(&with(4, "300")), Ok(Format::Sam));
    assert_ne!(sam(&with(5, "10")), Ok(Format::Sam));
    assert_ne!(sam(&with(5, "M10")), Ok(Format::Sam));
    assert_ne!(sam(&with(7, "-1")), Ok(Format::Sam));
    assert_ne!(sam(&with(9, "AC1T")), Ok(Format::Sam));
}
