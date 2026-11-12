use gloci_io::cancel::{RunGuard, request_cancel};
use gloci_io::memio::{discard, display_name, is_mem_path, register_input, register_output};
use gloci_io::*;
use std::io::{Cursor, Write};
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex, MutexGuard};

// The cancel flag is process-global: serialise.
static SERIAL: Mutex<()> = Mutex::new(());
fn serial() -> MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

fn ropts() -> ReaderOptions {
    ReaderOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        format: None,
        vcf_span: Default::default(),
        invalid: Default::default(),
        chrom_map: None,
        accept: Default::default(),
    }
}

fn oopts(write_index: bool) -> OutputOptions {
    OutputOptions {
        threads: NonZeroUsize::new(2).unwrap(),
        level: 6,
        write_index,
    }
}

#[derive(Clone, Default)]
struct Shared(Arc<Mutex<Vec<u8>>>);
impl Write for Shared {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

const BED: &[u8] = b"chr1\t1\t5\nchr1\t7\t9\n";

fn count(r: &mut Box<dyn RecordReader>) -> Result<usize, IoError> {
    let mut chroms = gloci_core::ChromTable::lexicographic();
    let mut batch = RecordBatch::default();
    let mut n = 0;
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 100)?;
        n += batch.len();
        if !more {
            return Ok(n);
        }
    }
}

#[test]
fn mem_input_reads_like_stdin() {
    let _s = serial();
    let name = register_input("<list>", Box::new(Cursor::new(BED.to_vec())));
    assert!(is_mem_path(&name));
    assert!(!is_regular_file(&name));
    let mut r = open_reader(&name, &ropts()).unwrap();
    assert_eq!(r.name(), "<list>");
    assert_eq!(count(&mut r).unwrap(), 2);
    let e = open_reader(&name, &ropts()).err().unwrap();
    assert!(e.to_string().contains("can be read only once"), "{e}");
    assert!(e.to_string().contains("<list>"), "{e}");
    assert_eq!(display_name(&name), "<list>");
    discard(&name);
    assert_eq!(display_name("plain.bed"), "plain.bed");
}

#[test]
fn mem_output_receives_exact_bytes() {
    let _s = serial();
    let sink = Shared::default();
    let name = register_output("<intersect output>", Box::new(sink.clone()));
    let spec = parse_output_spec(Some(std::path::Path::new(&name)), Some("bed")).unwrap();
    assert!(spec.direct);
    let mut out = open_output(&spec, OutFormat::Bed, &oopts(false)).unwrap();
    out.write_all(BED).unwrap();
    out.finish().unwrap();
    assert_eq!(*sink.0.lock().unwrap(), BED);

    let sink = Shared::default();
    let name = register_output("<out>", Box::new(sink.clone()));
    let spec = parse_output_spec(Some(std::path::Path::new(&name)), Some("bed.gz")).unwrap();
    let mut out = open_output(&spec, OutFormat::Bed, &oopts(false)).unwrap();
    out.write_all(BED).unwrap();
    out.finish().unwrap();
    let bytes = sink.0.lock().unwrap().clone();
    assert_eq!(detect_compression(&bytes), Compression::Bgzf);
    assert_eq!(&bytes[bytes.len() - 28..], &BGZF_EOF);
    let mut text = Vec::new();
    std::io::Read::read_to_end(&mut noodles_bgzf::io::Reader::new(&bytes[..]), &mut text).unwrap();
    assert_eq!(text, BED);

    let name = register_output("<out>", Box::new(Shared::default()));
    let spec = parse_output_spec(Some(std::path::Path::new(&name)), Some("bed.gz")).unwrap();
    let e = open_output(&spec, OutFormat::Bed, &oopts(true))
        .err()
        .unwrap();
    assert!(matches!(e, IoError::Usage(_)), "{e}");
    assert!(e.to_string().contains("needs -o FILE"), "{e}");
}

#[test]
fn mem_names_never_collide_with_files() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("gloci:1");
    std::fs::write(&p, BED).unwrap();
    let _ = register_input("<list>", Box::new(Cursor::new(Vec::new())));
    let mut r = open_reader(p.to_str().unwrap(), &ropts()).unwrap();
    assert_eq!(count(&mut r).unwrap(), 2);
    assert!(!is_mem_path(p.to_str().unwrap()));
}

#[test]
fn cancel_fails_reads_and_writes() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("in.bed");
    std::fs::write(&bed, BED).unwrap();
    let dest = dir.path().join("out.bed");
    {
        let _g = RunGuard::enter();
        let spec = parse_output_spec(Some(&dest), None).unwrap();
        let mut out = open_output(&spec, OutFormat::Bed, &oopts(false)).unwrap();
        out.write_all(BED).unwrap();
        request_cancel();
        assert!(gloci_io::cancel::is_cancelled());
        let e = out.write_all(BED).unwrap_err();
        assert!(gloci_io::cancel::is_cancel_error(&e));
        drop(out);
        let e = open_reader(bed.to_str().unwrap(), &ropts()).err();
        assert!(e.is_some_and(|e| e.is_cancelled()));
    }
    assert!(!gloci_io::cancel::is_cancelled());
    assert!(!dest.exists());
    assert!(
        std::fs::read_dir(dir.path()).unwrap().all(|e| !e
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains(".tmp-"))
    );
    // A reader opened before the cancel fails with Cancelled at read_batch.
    let mut r = open_reader(bed.to_str().unwrap(), &ropts()).unwrap();
    {
        let _g = RunGuard::enter();
        request_cancel();
        let e = count(&mut r).unwrap_err();
        assert!(matches!(e, IoError::Cancelled), "{e:?}");
        assert_eq!(e.to_string(), "interrupted; the partial output was removed");
    }
    assert_eq!(count(&mut r).unwrap(), 2);
}
