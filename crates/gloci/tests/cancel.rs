//! Cancellation (D91) at command level: a run whose token fires while it works fails with a
//! cancellation (never a write or read error), exits 1 and leaves no output file.
use gloci_io::cancel::{CancelToken, RunGuard};
use gloci_io::memio::register_input;
use std::ffi::OsString;
use std::io::{Cursor, Read};
use std::path::Path;
use std::sync::{Mutex, MutexGuard};

// Commands install process-wide state (the rayon pool, the logger): one at a time.
static SERIAL: Mutex<()> = Mutex::new(());
fn serial() -> MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

/// A source that cancels its run when it reaches the end of its bytes, so the cancellation
/// lands after the input was read, inside the command's own loops.
struct CancelAtEof {
    inner: Cursor<Vec<u8>>,
    token: CancelToken,
}

impl Read for CancelAtEof {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n == 0 {
            self.token.cancel();
        }
        Ok(n)
    }
}

fn mem(data: Vec<u8>, token: &CancelToken) -> String {
    register_input(
        "<generator>",
        Box::new(CancelAtEof {
            inner: Cursor::new(data),
            token: token.clone(),
        }),
    )
}

fn args(a: &[&str]) -> Vec<OsString> {
    std::iter::once("gloci")
        .chain(a.iter().copied())
        .map(OsString::from)
        .collect()
}

fn bed(n: u64) -> Vec<u8> {
    let mut s = String::new();
    for i in (0..n).rev() {
        s.push_str(&format!("chr1\t{}\t{}\n", i * 10, i * 10 + 5));
    }
    s.into_bytes()
}

fn bed_asc(n: u64) -> Vec<u8> {
    let mut v = bed(n);
    v = String::from_utf8(v)
        .unwrap()
        .lines()
        .rev()
        .flat_map(|l| format!("{l}\n").into_bytes())
        .collect();
    v
}

fn vcf(n: u64) -> Vec<u8> {
    let mut s = String::from("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n");
    s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for i in 0..n {
        s.push_str(&format!("chr1\t{}\t.\tA\tC\t50\tPASS\t.\n", i + 1));
    }
    s.into_bytes()
}

fn files_in(dir: &Path) -> Vec<String> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect()
}

/// Run `build(mem names)` under a token that the first input cancels at its end.
fn cancelled_run(inputs: Vec<Vec<u8>>, build: impl Fn(&[String], &str) -> Vec<String>) {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.bed");
    let token = CancelToken::new();
    let names: Vec<String> = inputs.into_iter().map(|d| mem(d, &token)).collect();
    let argv = build(&names, out.to_str().unwrap());
    let argv: Vec<&str> = argv.iter().map(String::as_str).collect();
    let (code, msg) = {
        let _g = RunGuard::enter(&token);
        gloci::run_reporting(args(&argv))
    };
    assert_eq!(msg.as_deref(), Some("interrupted"), "{argv:?}");
    assert_eq!(
        format!("{code:?}"),
        format!("{:?}", std::process::ExitCode::from(1))
    );
    assert!(
        files_in(dir.path()).is_empty(),
        "{:?}",
        files_in(dir.path())
    );
}

#[test]
fn sort_with_spills_is_cancelled() {
    cancelled_run(vec![bed(1_200_000)], |m, out| {
        [
            "sort",
            &m[0],
            "--max-mem",
            "16M",
            "--tmp-dir",
            Path::new(out).parent().unwrap().to_str().unwrap(),
            "-o",
            out,
        ]
        .map(String::from)
        .to_vec()
    });
}

#[test]
fn sweep_is_cancelled() {
    cancelled_run(vec![bed_asc(5000), bed_asc(5000)], |m, out| {
        ["intersect", "-a", &m[0], "-b", &m[1], "-o", out]
            .map(String::from)
            .to_vec()
    });
}

#[test]
fn selection_stage_with_workers_is_cancelled() {
    cancelled_run(vec![vcf(50_000)], |m, out| {
        [
            "vcf",
            "filter",
            &m[0],
            "-t",
            "4",
            "-o",
            &out.replace(".bed", ".vcf"),
        ]
        .map(String::from)
        .to_vec()
    });
}
