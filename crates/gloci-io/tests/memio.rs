use gloci_io::cancel::{CancelToken, RunGuard};
use gloci_io::memio::{
    MEM_PREFIX, discard, display_name, is_mem_path, register_input, register_output,
};
use gloci_io::*;
use std::io::{Cursor, Write};
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex, MutexGuard};

// Serialise: tests share the registry.
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
    // A registered name, and a real file named like its printable part ("gloci:<id>"),
    // holding other data: each opens its own data.
    let name = register_input("<list>", Box::new(Cursor::new(BED.to_vec())));
    let file = dir.path().join(
        name.strip_prefix(MEM_PREFIX)
            .map(|id| format!("gloci:{id}"))
            .unwrap(),
    );
    std::fs::write(&file, b"chr1\t1\t2\n").unwrap();
    assert!(is_regular_file(file.to_str().unwrap()));
    assert!(!is_mem_path(file.to_str().unwrap()));
    let mut r = open_reader(file.to_str().unwrap(), &ropts()).unwrap();
    assert_eq!(count(&mut r).unwrap(), 1);
    let mut m = open_reader(&name, &ropts()).unwrap();
    assert_eq!(m.name(), "<list>");
    assert_eq!(count(&mut m).unwrap(), 2);
}

#[test]
fn cancel_fails_reads_and_writes() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("in.bed");
    std::fs::write(&bed, BED).unwrap();
    let dest = dir.path().join("out.bed");
    let token = CancelToken::new();
    {
        let _g = RunGuard::enter(&token);
        let spec = parse_output_spec(Some(&dest), None).unwrap();
        let mut out = open_output(&spec, OutFormat::Bed, &oopts(false)).unwrap();
        out.write_all(BED).unwrap();
        token.cancel();
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
    // A reader built inside a run fails with Cancelled at read_batch, on any thread, and a
    // later run (new token) is unaffected.
    let token = CancelToken::new();
    let mut r = {
        let _g = RunGuard::enter(&token);
        open_reader(bed.to_str().unwrap(), &ropts()).unwrap()
    };
    token.cancel();
    let e = count(&mut r).unwrap_err();
    assert!(matches!(e, IoError::Cancelled), "{e:?}");
    assert_eq!(e.to_string(), "interrupted");
    let mut again = open_reader(bed.to_str().unwrap(), &ropts()).unwrap();
    assert_eq!(count(&mut again).unwrap(), 2);
}

#[test]
fn a_cancel_stays_with_its_run() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("in.bed");
    let big: Vec<u8> = (0..20000)
        .flat_map(|i| format!("chr1\t{}\t{}\n", i * 10, i * 10 + 5).into_bytes())
        .collect();
    std::fs::write(&bed, &big).unwrap();
    let (a, b) = (CancelToken::new(), CancelToken::new());
    let open = |t: &CancelToken| {
        let _g = RunGuard::enter(t);
        open_reader(bed.to_str().unwrap(), &ropts()).unwrap()
    };
    let (mut ra, mut rb) = (open(&a), open(&b));
    // Cancelling a run that has ended (nothing alive uses `c`) changes nothing for later runs.
    let c = CancelToken::new();
    c.cancel();
    let mut rd = open(&CancelToken::new());
    a.cancel();
    assert!(matches!(count(&mut ra), Err(IoError::Cancelled)));
    assert_eq!(count(&mut rb).unwrap(), 20000);
    assert_eq!(count(&mut rd).unwrap(), 20000);
    // No token installed: nothing can cancel the read.
    let mut none = open_reader(bed.to_str().unwrap(), &ropts()).unwrap();
    assert_eq!(count(&mut none).unwrap(), 20000);
}

#[test]
fn a_cancelled_memory_output_gets_no_end_of_file_block() {
    let _s = serial();
    let token = CancelToken::new();
    let sink = Shared::default();
    {
        let _g = RunGuard::enter(&token);
        let name = register_output("<out>", Box::new(sink.clone()));
        let spec = parse_output_spec(Some(std::path::Path::new(&name)), Some("bed.gz")).unwrap();
        let mut out = open_output(&spec, OutFormat::Bed, &oopts(false)).unwrap();
        out.write_all(BED).unwrap();
        token.cancel();
        assert!(out.write_all(BED).is_err());
        drop(out);
    }
    let bytes = sink.0.lock().unwrap().clone();
    assert!(
        bytes.len() < 28 || bytes[bytes.len() - 28..] != BGZF_EOF,
        "{} bytes",
        bytes.len()
    );
}
