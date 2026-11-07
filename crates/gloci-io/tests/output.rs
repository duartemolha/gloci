use gloci_io::*;
use std::fs;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::{Mutex, MutexGuard};

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

// cleanup_registered() is process-global; serialise tests that touch the registry.
static SERIAL: Mutex<()> = Mutex::new(());
fn serial() -> MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

fn opts(threads: usize, write_index: bool) -> OutputOptions {
    OutputOptions {
        threads: NonZeroUsize::new(threads).unwrap(),
        level: 6,
        write_index,
    }
}

fn no_tmp_files(dir: &Path) -> bool {
    fs::read_dir(dir)
        .unwrap()
        .all(|e| !e.unwrap().file_name().to_string_lossy().contains(".tmp-"))
}

fn write_via_output(dir: &Path, name: &str, data: &[u8], write_index: bool) -> Result<(), IoError> {
    let spec = parse_output_spec(Some(&dir.join(name)), None).unwrap();
    let mut out = open_output(&spec, spec.format.unwrap(), &opts(2, write_index))?;
    out.write_all(data).unwrap();
    out.finish()
}

fn bgzf_bytes(data: &[u8], threads: usize) -> Vec<u8> {
    let dir = tempfile::tempdir().unwrap();
    write_via_output_threads(dir.path(), "o.bed.gz", data, threads);
    fs::read(dir.path().join("o.bed.gz")).unwrap()
}

fn write_via_output_threads(dir: &Path, name: &str, data: &[u8], threads: usize) {
    let spec = parse_output_spec(Some(&dir.join(name)), None).unwrap();
    let mut out = open_output(&spec, spec.format.unwrap(), &opts(threads, false)).unwrap();
    // write in odd-sized pieces to check chunking does not leak into block boundaries
    for c in data.chunks(7919) {
        out.write_all(c).unwrap();
    }
    out.finish().unwrap();
}

#[test]
fn spec_parsing() {
    let s = parse_output_spec(Some(Path::new("o.bed.gz")), None).unwrap();
    assert_eq!((s.format, s.bgzf), (Some(OutFormat::Bed), true));
    assert_eq!(
        parse_output_spec(Some(Path::new("O.VCF")), None)
            .unwrap()
            .format,
        Some(OutFormat::Vcf)
    );
    assert!(
        parse_output_spec(Some(Path::new("o.bcf")), None)
            .unwrap()
            .bgzf
    );
    assert!(
        parse_output_spec(Some(Path::new("o.txt")), None)
            .unwrap_err()
            .0
            .contains(".bed.gz")
    );
    assert!(parse_output_spec(None, Some("vcf.gz")).unwrap().bgzf);
    assert!(parse_output_spec(Some(Path::new("a.bed")), Some("vcf")).is_err());
    assert_eq!(
        parse_output_spec(Some(Path::new("a.gff3.gz")), None)
            .unwrap()
            .format,
        Some(OutFormat::Gff3)
    );
    assert_eq!(parse_output_spec(None, None).unwrap().format, None);
    assert_eq!(OutFormat::family(Format::BedGraph), OutFormat::BedGraph);
}

#[test]
fn level_and_tty() {
    assert!(check_level(12).is_ok() && check_level(13).is_err());
    assert!(refuse_binary_to_tty(true, OutFormat::Bcf, true).is_err());
    assert!(refuse_binary_to_tty(false, OutFormat::Bcf, true).is_ok());
    assert!(refuse_binary_to_tty(true, OutFormat::Bed, false).is_ok());
    assert!(refuse_binary_to_tty(true, OutFormat::Bed, true).is_err());
}

#[test]
fn atomic_file_semantics() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let dest = dir.path().join("out.txt");
    {
        let mut f = AtomicFile::create(&dest).unwrap();
        f.write_all(b"x").unwrap();
    }
    assert!(!dest.exists() && no_tmp_files(dir.path()));
    fs::write(&dest, b"old").unwrap();
    {
        let mut f = AtomicFile::create(&dest).unwrap();
        f.write_all(b"new").unwrap();
    }
    assert_eq!(fs::read(&dest).unwrap(), b"old");
    let mut f = AtomicFile::create(&dest).unwrap();
    f.write_all(b"new").unwrap();
    f.commit().unwrap();
    assert_eq!(fs::read(&dest).unwrap(), b"new");
    assert!(no_tmp_files(dir.path()));
}

#[test]
fn stale_index_removed_and_write_index_builds() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("o.bed.gz.tbi"), b"stale").unwrap();
    fs::write(dir.path().join("o.bed.gz.csi"), b"stale").unwrap();
    write_via_output(dir.path(), "o.bed.gz", b"chr1\t1\t2\n", false).unwrap();
    assert!(!dir.path().join("o.bed.gz.tbi").exists());
    assert!(!dir.path().join("o.bed.gz.csi").exists());
    assert!(dir.path().join("o.bed.gz").exists());

    // --write-index builds the index from the finished temp file (Task 12)
    fs::write(dir.path().join("p.bed.gz.csi"), b"stale").unwrap();
    write_via_output(dir.path(), "p.bed.gz", b"chr1\t1\t2\n", true).unwrap();
    let tbi = noodles_tabix::fs::read(dir.path().join("p.bed.gz.tbi")).unwrap();
    assert_eq!(tbi.reference_sequences().len(), 1);
    assert!(!dir.path().join("p.bed.gz.csi").exists());
    // the index is published after the data, so it is not older (D26)
    let info = find_index(&dir.path().join("p.bed.gz")).unwrap();
    assert_eq!(info.kind, IndexKind::Tbi);
    assert!(!info.stale);
    // replacing an indexed file with --write-index keeps exactly one fresh index
    write_via_output(dir.path(), "p.bed.gz", b"chr1\t5\t9\nchr2\t1\t2\n", true).unwrap();
    let tbi = noodles_tabix::fs::read(dir.path().join("p.bed.gz.tbi")).unwrap();
    assert_eq!(tbi.reference_sequences().len(), 2);
    assert!(!dir.path().join("p.bed.gz.csi").exists());
    assert!(no_tmp_files(dir.path()));
    // an index that cannot be built publishes nothing
    let e =
        write_via_output(dir.path(), "u.bed.gz", b"chr1\t5\t6\nchr1\t1\t2\n", true).unwrap_err();
    assert!(matches!(e, IoError::Order(_)), "{e}");
    assert!(e.to_string().contains("u.bed.gz is not sorted"), "{e}");
    assert!(!dir.path().join("u.bed.gz").exists());
    assert!(!dir.path().join("u.bed.gz.tbi").exists());
    assert!(no_tmp_files(dir.path()));
    // validation: stdout and non-BGZF
    let stdout = OutputSpec {
        path: None,
        format: Some(OutFormat::Bed),
        bgzf: true,
    };
    assert!(matches!(
        open_output(&stdout, OutFormat::Bed, &opts(1, true)),
        Err(IoError::Usage(_))
    ));
    let e = write_via_output(dir.path(), "q.bed", b"x", true).unwrap_err();
    assert!(e.to_string().contains("BGZF"));
    assert!(no_tmp_files(dir.path()));
}

#[test]
fn bgzf_is_deterministic_and_terminated() {
    let _s = serial();
    let five_mb: Vec<u8> = (0..5_000_000u32)
        .map(|i| b"ACGT\n\t0123"[(i.wrapping_mul(2654435761) >> 7) as usize % 10])
        .collect();
    let a = bgzf_bytes(&five_mb, 1);
    assert_eq!(a, bgzf_bytes(&five_mb, 8));
    assert!(a.ends_with(&BGZF_EOF));
    // round trip
    let mut got = Vec::new();
    std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(&a[..]), &mut got).unwrap();
    assert_eq!(got, five_mb);
}

#[test]
fn bcf_output_is_always_bgzf_and_empty_output_is_valid() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    write_via_output(dir.path(), "e.bed.gz", b"", false).unwrap();
    assert_eq!(fs::read(dir.path().join("e.bed.gz")).unwrap(), BGZF_EOF);
    let spec = OutputSpec {
        path: Some(dir.path().join("x.dat")),
        format: Some(OutFormat::Bcf),
        bgzf: false,
    };
    let out = open_output(&spec, OutFormat::Bcf, &opts(1, false)).unwrap();
    out.finish().unwrap();
    assert_eq!(fs::read(dir.path().join("x.dat")).unwrap(), BGZF_EOF);
}

#[test]
fn cleanup_registry() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("spill");
    fs::create_dir(&p).unwrap();
    fs::write(p.join("a"), b"1").unwrap();
    let _g = register_temp(&p);
    let f = dir.path().join("t.tmp-abc");
    fs::write(&f, b"1").unwrap();
    let _g2 = register_temp(&f);
    cleanup_registered();
    assert!(!p.exists() && !f.exists());
}
