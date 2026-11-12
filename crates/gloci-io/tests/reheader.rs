//! The byte-level parts of `reheader` (D130): the BAM header's exact length, the header as BGZF
//! blocks, the CRAM header container, and the HEADER argument's sources.
use gloci_io::bgzf_copy::{inflate, read_block};
use gloci_io::{
    AlnHeader, BGZF_EOF, IoError, bam_header_blocks, read_bam_file_header,
    read_bam_header_with_extent, read_cram_prelude, read_new_header, write_cram_header_container,
};
use std::io::Write;
use std::path::PathBuf;

fn fixture(dir: &str, name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs")
        .join(dir)
        .join(name)
}

fn inflate_all(mut data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    while let Some(b) = read_block(&mut data).unwrap() {
        out.extend(inflate(&b).unwrap());
    }
    out
}

fn sam(text: &str) -> AlnHeader {
    AlnHeader::parse_sam(text.as_bytes().to_vec()).unwrap()
}

#[test]
fn header_blocks_end_on_a_block_edge() {
    // 200,000 references: about 6 MB of text, many blocks, the last one part full.
    let mut text = String::from("@HD\tVN:1.6\n");
    for i in 0..200_000 {
        text.push_str(&format!("@SQ\tSN:c{i}\tLN:{}\n", 100 + i));
    }
    let h = sam(&text);
    let blocks = bam_header_blocks(&h, 6).unwrap();
    assert!(!blocks.ends_with(&BGZF_EOF));
    let data = inflate_all(&blocks);
    let (back, extent) = read_bam_header_with_extent(&mut data.as_slice()).unwrap();
    assert_eq!(extent as usize, data.len());
    assert_eq!(back.refs, h.refs);
    assert_eq!(back.text, h.text);
    // Many blocks, and the sizes of the stored blocks add up to the bytes.
    let (mut src, mut n, mut stored) = (blocks.as_slice(), 0, 0);
    while let Some(b) = read_block(&mut src).unwrap() {
        n += 1;
        stored += b.raw.len();
    }
    assert!(n > 90, "{n} blocks");
    assert_eq!(stored, blocks.len());
    // The level reaches the compressor.
    assert!(bam_header_blocks(&h, 0).unwrap().len() > bam_header_blocks(&h, 9).unwrap().len());
    assert!(bam_header_blocks(&h, 13).is_err());
}

#[test]
fn extent_of_a_real_bam_is_where_the_records_start() {
    let bytes = std::fs::read(fixture("samtools", "reheader-in.bam")).unwrap();
    let data = inflate_all(&bytes);
    let (h, extent) = read_bam_header_with_extent(&mut data.as_slice()).unwrap();
    assert_eq!(h.refs.len(), 3);
    // The first record follows: block_size, refID 0, pos 9 (r1 at chr1:10).
    let rest = &data[extent as usize..];
    assert_eq!(&rest[4..8], 0i32.to_le_bytes());
    assert_eq!(&rest[8..12], 9i32.to_le_bytes());
    let (h2, e2) = read_bam_file_header(&fixture("samtools", "reheader-in.bam")).unwrap();
    assert_eq!((h2.refs, e2), (h.refs, extent));
}

/// A CRAM file: definition, `container`, and then the end-of-file container of the version.
fn cram_with(text: &[u8]) -> Vec<u8> {
    const EOF_V3: [u8; 38] = [
        0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00,
        0x00, 0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00,
        0x01, 0x00, 0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
    ];
    let mut f = b"CRAM".to_vec();
    f.extend_from_slice(&[3, 0]);
    f.extend_from_slice(&[b'x'; 20]);
    write_cram_header_container(text, &mut f).unwrap();
    f.extend_from_slice(&EOF_V3);
    f
}

#[test]
fn cram_v3_header_container_reads_back() {
    let dir = tempfile::tempdir().unwrap();
    // Texts around the ITF8 width changes: the block size is the text plus 4.
    {
        for n in [
            0usize, 1, 100, 124, 125, 126, 16_000, 16_380, 16_381, 70_000, 3_000_000,
        ] {
            let mut text = b"@HD\tVN:1.6\n".to_vec();
            while text.len() < n {
                text.extend_from_slice(b"@CO\tpadding to a length\n");
            }
            text.truncate(n.max(text.len().min(n)));
            if text.last() != Some(&b'\n') && !text.is_empty() {
                text.push(b'\n');
            }
            let p = dir.path().join("x.cram");
            let bytes = cram_with(&text);
            std::fs::write(&p, &bytes).unwrap();
            let got = read_cram_prelude(&p)
                .unwrap_or_else(|e| panic!("text of {} bytes: {e}", text.len()));
            assert_eq!(got.header_text, text);
            assert_eq!(got.major, 3);
            assert_eq!(got.data_start as usize, bytes.len() - 38);
            assert!(gloci_io::cram_eof_present(&p).unwrap());
        }
    }
}

#[test]
fn cram_v3_container_checksums_are_right() {
    let mut f = Vec::new();
    write_cram_header_container(b"@HD\tVN:1.6\n", &mut f).unwrap();
    // container header: 4 + 8 bytes then its CRC32; the block (method 0, type 0, id 0, sizes)
    // ends with its CRC32 as well.
    let crc = |b: &[u8]| {
        let mut c = flate2::Crc::new();
        c.update(b);
        c.sum().to_le_bytes()
    };
    assert_eq!(f[12..16], crc(&f[..12]));
    let block = &f[16..];
    assert_eq!(block[block.len() - 4..], crc(&block[..block.len() - 4]));
}

#[test]
fn new_header_sources() {
    let dir = tempfile::tempdir().unwrap();
    let text = "@HD\tVN:1.6\n@SQ\tSN:a\tLN:10\n@CO\tx\n";
    let write = |name: &str, bytes: &[u8]| {
        let p = dir.path().join(name);
        std::fs::write(&p, bytes).unwrap();
        p.to_string_lossy().into_owned()
    };
    // Plain text, with and without a final newline, with CRLF, and with reads after it.
    for body in [
        text.to_string(),
        text.trim_end().to_string(),
        format!("{text}r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"),
    ] {
        let h = read_new_header(&write("a.sam", body.as_bytes())).unwrap();
        assert_eq!(String::from_utf8(h.text).unwrap(), text);
        assert_eq!(h.refs, [(b"a".to_vec(), 10)]);
    }
    // Gzip and BGZF text.
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(text.as_bytes()).unwrap();
    let h = read_new_header(&write("a.sam.gz", &gz.finish().unwrap())).unwrap();
    assert_eq!(h.refs.len(), 1);
    let mut bgzf = gloci_io::bgzf_copy::deflate_blocks(text.as_bytes()).unwrap();
    bgzf.extend_from_slice(&BGZF_EOF);
    let h = read_new_header(&write("b.sam.gz", &bgzf)).unwrap();
    assert_eq!(h.refs.len(), 1);
    // A BAM header file: a header with no @SQ lines in its text gets them from its list.
    let mut raw = Vec::new();
    sam("@HD\tVN:1.6\n@SQ\tSN:q\tLN:5\n")
        .write_bam(&mut raw)
        .unwrap();
    let bare = AlnHeader::new(
        b"@HD\tVN:1.6\n".to_vec(),
        vec![(b"q".to_vec(), 5)],
        gloci_io::SortOrder::Unknown,
    );
    let mut raw2 = Vec::new();
    bare.write_bam(&mut raw2).unwrap();
    for raw in [raw, raw2] {
        let mut bam = gloci_io::bgzf_copy::deflate_blocks(&raw).unwrap();
        bam.extend_from_slice(&BGZF_EOF);
        let h = read_new_header(&write("h.bam", &bam)).unwrap();
        assert_eq!(h.refs, [(b"q".to_vec(), 5)]);
        assert!(
            String::from_utf8(h.text)
                .unwrap()
                .contains("@SQ\tSN:q\tLN:5")
        );
    }
    // A CRAM file: its header text.
    let p = dir.path().join("h.cram");
    std::fs::write(&p, cram_with(text.as_bytes())).unwrap();
    let h = read_new_header(p.to_str().unwrap()).unwrap();
    assert_eq!(h.refs, [(b"a".to_vec(), 10)]);
    // Errors: no header lines, an @SQ without LN, a missing file.
    assert!(matches!(
        read_new_header(&write("n.sam", b"r1\t4\n")),
        Err(IoError::Detect { .. })
    ));
    let e = read_new_header(&write("n2.sam", b"@SQ\tSN:a\n")).unwrap_err();
    assert!(e.to_string().contains("LN"), "{e}");
    assert!(matches!(
        read_new_header(&dir.path().join("none").to_string_lossy()),
        Err(IoError::Open { .. })
    ));
}

#[test]
fn cram_2x_is_refused_by_name_and_version() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("old.cram");
    let mut f = b"CRAM\x02\x01".to_vec();
    f.extend_from_slice(&[b'x'; 20]);
    f.extend_from_slice(&[0u8; 60]);
    std::fs::write(&p, &f).unwrap();
    // One message, wherever the version is read: the record readers' (frame.rs), with the
    // reference the rewrite needs.
    let ro = gloci_io::ReaderOptions {
        accept: gloci_io::Accept {
            alignments: Some(gloci_io::AlignOptions::default()),
            ..gloci_io::Accept::default()
        },
        chrom_map: None,
        threads: std::num::NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::default(),
        invalid: gloci_io::InvalidPolicy::Fail,
    };
    let path = p.display().to_string();
    let reader = gloci_io::open_aln_source(&path, &ro)
        .err()
        .unwrap()
        .to_string();
    for e in [
        read_cram_prelude(&p).unwrap_err(),
        gloci_io::check_cram_version(&p).unwrap_err(),
    ] {
        let m = e.to_string();
        assert!(matches!(e, IoError::Detect { .. }));
        assert!(m.contains("old.cram is CRAM 2.1"), "{m}");
        assert!(
            m.contains("CRAM 3.0 and 3.1")
                && m.contains("samtools view -C -T ref.fa --output-fmt-option version=3.0"),
            "{m}"
        );
        assert_eq!(m, reader);
    }
}
