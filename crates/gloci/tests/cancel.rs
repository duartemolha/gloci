//! Cancellation (D91) at command level: a run whose token fires while it works fails with a
//! cancellation (never a write or read error), exits 1 and leaves no output file.
use gloci_io::cancel::{RunGuard, RunToken};
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
    token: RunToken,
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

/// A source that cancels its run once it has given `after` bytes, so the cancellation lands
/// while the command reads it.
struct CancelAfter {
    inner: Cursor<Vec<u8>>,
    after: u64,
    token: RunToken,
}

impl Read for CancelAfter {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        if self.inner.position() >= self.after {
            self.token.cancel();
        }
        Ok(n)
    }
}

fn mem(data: Vec<u8>, token: &RunToken) -> String {
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
    cancelled_run_named(
        |token| inputs.into_iter().map(|d| mem(d, token)).collect(),
        build,
    );
}

/// Run `build(names, output)` under a token, with the inputs `register` registers.
fn cancelled_run_named(
    register: impl FnOnce(&RunToken) -> Vec<String>,
    build: impl Fn(&[String], &str) -> Vec<String>,
) {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.bed");
    let token = RunToken::new();
    let names = register(&token);
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

/// A CRAM's file definition and header container, its data containers, and its end-of-file
/// container.
fn cram_split(cram: &[u8]) -> (&[u8], &[u8], &[u8]) {
    fn int(b: &[u8], at: &mut usize, long: bool) {
        // ITF8 takes at most four bytes after the first, LTF8 eight.
        let ones = b[*at].leading_ones() as usize;
        *at += 1 + if long { ones } else { ones.min(4) };
    }
    let mut at = 26 + 4;
    for long in [false, false, false, false, true, true, false] {
        int(cram, &mut at, long);
    }
    let landmarks = cram[at] as usize; // the header container has one or none
    int(cram, &mut at, false);
    for _ in 0..landmarks {
        int(cram, &mut at, false);
    }
    at += 4;
    let len = u32::from_le_bytes([cram[26], cram[27], cram[28], cram[29]]) as usize;
    let data = at + len;
    let eof = cram.len() - 38;
    (&cram[..data], &cram[data..eof], &cram[eof..])
}

/// `flagstat` of a CRAM, whose fields it decodes alone, cancelled while it reads the file (past
/// the 64 KiB the input peeks at open): interrupted, exit 1, no output, at -t 1 and 2.
#[test]
fn cram_flagstat_is_cancelled() {
    let cram = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/aln/mates.noref.cram"),
    )
    .unwrap();
    let (head, data, eof) = cram_split(&cram);
    let big: Vec<u8> = [head, &data.repeat(30), eof].concat();
    {
        // Uncancelled, the file is read whole: 30 times mates.noref.cram's 182 + 16 reads.
        let _s = serial();
        let dir = tempfile::tempdir().unwrap();
        let (input, out) = (dir.path().join("big.cram"), dir.path().join("out.txt"));
        std::fs::write(&input, &big).unwrap();
        let argv = [
            "flagstat",
            input.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ];
        let (code, msg) = gloci::run_reporting(args(&argv));
        assert_eq!(
            format!("{code:?}"),
            format!("{:?}", std::process::ExitCode::SUCCESS),
            "{msg:?}"
        );
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.starts_with("5460 + 480 in total"), "{text}");
    }
    for t in ["1", "2"] {
        let big = big.clone();
        cancelled_run_named(
            |token| {
                vec![register_input(
                    "<generator>",
                    Box::new(CancelAfter {
                        after: big.len() as u64 / 2,
                        inner: Cursor::new(big),
                        token: token.clone(),
                    }),
                )]
            },
            |m, out| {
                let out = out.replace(".bed", ".txt");
                ["flagstat", &m[0], "-t", t, "-o", &out]
                    .map(String::from)
                    .to_vec()
            },
        );
    }
}
