//! The stage runtime (D85): bounded in-memory pipes, and several commands running at once on
//! threads, chained by those pipes, with one error reported for the set and cancellation.
use gloci::embed::{
    self, ExecOptions, Outcome, PIPE_CHUNK, PIPE_CHUNKS, PipeFail, SiteKind, StageError, StageSet,
    StageSpec, pipe,
};
use gloci_io::Loc;
use gloci_io::cancel::{self, CancelToken, RunGuard};
use gloci_io::memio::Registered;
use std::ffi::OsString;
use std::io::{Cursor, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::{Duration, Instant};

/// Every test here runs alone: one counts the process's threads, and the others start some.
static SERIAL: Mutex<()> = Mutex::new(());

fn serial() -> MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|p| p.into_inner())
}

#[test]
fn pipe_is_bounded() {
    let _s = serial();
    let (mut w, r) = pipe(PIPE_CHUNKS, PIPE_CHUNK);
    let accepted = Arc::new(AtomicUsize::new(0));
    let writer = {
        let accepted = accepted.clone();
        std::thread::spawn(move || {
            let piece = [b'x'; 1000];
            // 64 MB, far more than the pipe holds; stops at the first error.
            for _ in 0..64_000 {
                if let Err(e) = w.write_all(&piece) {
                    return Some(e.kind());
                }
                accepted.fetch_add(piece.len(), Ordering::SeqCst);
            }
            None
        })
    };
    let full = PIPE_CHUNKS * PIPE_CHUNK;
    let load = || accepted.load(Ordering::SeqCst);
    // The pipe takes at least its chunks' worth (however slow this machine is) ...
    let deadline = Instant::now() + Duration::from_secs(30);
    while load() < full - 1000 {
        assert!(!writer.is_finished(), "the writer ended early");
        assert!(Instant::now() < deadline, "only {} bytes accepted", load());
        std::thread::sleep(Duration::from_millis(10));
    }
    // ... then the writer stops: no progress over a window ...
    let held = loop {
        let before = load();
        std::thread::sleep(Duration::from_millis(200));
        if load() == before || writer.is_finished() {
            break load();
        }
        assert!(
            Instant::now() < deadline,
            "still writing at {} bytes",
            load()
        );
    };
    assert!(
        !writer.is_finished(),
        "the writer never blocked ({held} bytes accepted)"
    );
    // ... and none over a further one.
    std::thread::sleep(Duration::from_millis(300));
    assert_eq!(load(), held, "the writer went on after stopping");
    // At most the channel's chunks, plus the one partial buffer being filled, are unread.
    assert!(
        held <= full + PIPE_CHUNK,
        "{held} bytes accepted before blocking"
    );
    // The reader going away unblocks it with BrokenPipe.
    drop(r);
    assert_eq!(writer.join().unwrap(), Some(ErrorKind::BrokenPipe));
}

#[test]
fn pipe_ends_eof_fail_and_broken_pipe() {
    let _s = serial();
    // Dropping the writer is EOF, after every byte.
    let (mut w, mut r) = pipe(2, 8);
    let t = std::thread::spawn(move || {
        for i in 0..1000u32 {
            writeln!(w, "{i}").unwrap();
        }
    });
    let mut got = String::new();
    r.read_to_string(&mut got).unwrap();
    t.join().unwrap();
    assert_eq!(got, (0..1000).map(|i| format!("{i}\n")).collect::<String>());

    // flush sends a partial chunk at once.
    let (mut w, mut r) = pipe(2, 1 << 20);
    w.write_all(b"abc").unwrap();
    w.flush().unwrap();
    let mut b = [0u8; 3];
    r.read_exact(&mut b).unwrap();
    assert_eq!(&b, b"abc");
    // fail: the reader gets the message, not EOF.
    w.write_all(b"partial").unwrap();
    w.fail("bad line 3");
    let e = r.read_to_end(&mut Vec::new()).unwrap_err();
    assert_eq!(e.to_string(), "upstream step failed: bad line 3");

    // A fail handle works once the writer is gone (moved into a sink and dropped) ...
    let (w, mut r) = pipe(2, 8);
    let f = w.fail_handle();
    drop(w);
    f.clone().fail("stage 0 failed");
    let e = r.read(&mut [0u8; 4]).unwrap_err();
    assert_eq!(e.to_string(), "upstream step failed: stage 0 failed");
    drop(f);

    // ... and the reader sees EOF only once the writer and every fail handle are gone, so a
    // failure reported after the writer was dropped is never mistaken for the end.
    let (mut w, mut r) = pipe(2, 8);
    let f = w.fail_handle();
    w.write_all(b"done").unwrap();
    drop(w);
    let reader = std::thread::spawn(move || {
        let mut v = Vec::new();
        r.read_to_end(&mut v).map(|_| v)
    });
    std::thread::sleep(Duration::from_millis(200));
    assert!(
        !reader.is_finished(),
        "EOF before the fail handle was dropped"
    );
    drop(f);
    assert_eq!(reader.join().unwrap().unwrap(), b"done");

    // The reader going away is BrokenPipe for the writer.
    let (mut w, r) = pipe(2, 8);
    drop(r);
    let e = w.write_all(b"0123456789").unwrap_err();
    assert_eq!(e.kind(), ErrorKind::BrokenPipe);
}

#[test]
fn parked_pipe_ends_see_the_cancel() {
    let _s = serial();
    // A writer parked on a full pipe, and a reader parked on an empty one, each on a thread that
    // entered the run's token.
    let token = CancelToken::new();
    let (mut w, _r) = pipe(1, 4);
    let writer = {
        let token = token.clone();
        std::thread::spawn(move || {
            let _g = RunGuard::enter(&token);
            loop {
                if let Err(e) = w.write_all(b"xxxx") {
                    return e;
                }
            }
        })
    };
    let (_w, mut r) = pipe(1, 4);
    let reader = {
        let token = token.clone();
        std::thread::spawn(move || {
            let _g = RunGuard::enter(&token);
            r.read(&mut [0u8; 4]).unwrap_err()
        })
    };
    std::thread::sleep(Duration::from_millis(200));
    assert!(!writer.is_finished() && !reader.is_finished());
    token.cancel();
    let t0 = std::time::Instant::now();
    assert!(cancel::is_cancel_error(&writer.join().unwrap()));
    assert!(cancel::is_cancel_error(&reader.join().unwrap()));
    assert!(t0.elapsed() < Duration::from_secs(1), "{:?}", t0.elapsed());
}

/// An argv: `gloci`, then each argument (`&str`, `&Path`, an endpoint's name).
macro_rules! argv {
    ($($a:expr),* $(,)?) => { vec![OsString::from("gloci"), $(OsString::from($a)),*] };
}

type StageFn = Box<dyn FnOnce() -> anyhow::Result<bool> + Send>;

fn f(run: impl FnOnce() -> anyhow::Result<bool> + Send + 'static) -> StageFn {
    Box::new(run)
}

/// A sink whose bytes the test keeps.
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

impl Shared {
    fn bytes(&self) -> Vec<u8> {
        self.0.lock().unwrap().clone()
    }
}

/// A pipe between two stages: its writer registered as an output, its reader as an input, and
/// the writer's fail handle, for the writing stage's spec.
struct Link {
    out: Registered,
    inp: Registered,
    fail: PipeFail,
}

fn link(display: &str) -> Link {
    let (w, r) = pipe(PIPE_CHUNKS, PIPE_CHUNK);
    let fail = w.fail_handle();
    Link {
        out: Registered::output(display, Box::new(w)),
        inp: Registered::input(display, Box::new(r)),
        fail,
    }
}

fn stage(
    argv: Vec<OsString>,
    endpoints: Vec<Registered>,
    fail_downstream: Vec<PipeFail>,
) -> StageSpec {
    StageSpec {
        command_line: gloci_io::shell_quote_args(&argv),
        argv,
        tmp_dir: None,
        endpoints,
        fail_downstream,
    }
}

/// The set's result, failing the test (rather than hanging it) after two minutes. (The error
/// is boxed only to keep clippy's result_large_err quiet.)
fn finish(mut set: StageSet) -> Result<bool, Box<StageError>> {
    let deadline = Instant::now() + Duration::from_secs(120);
    loop {
        if let Some(r) = set.wait(Duration::from_millis(50)) {
            return r.map_err(Box::new);
        }
        assert!(Instant::now() < deadline, "the stages hung");
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/embed_stages")
        .join(name)
}

/// `gloci intersect --mode any -a A -b B EXTRA0 | gloci merge --max-gap 50 EXTRA1 -`, into a
/// sink.
fn intersect_merge(
    a: &Path,
    b: &Path,
    extra0: &[&str],
    extra1: &[&str],
) -> (Vec<StageSpec>, Shared) {
    let p = link("<intersect output>");
    let sink = Shared::default();
    let out = Registered::output("<merge output>", Box::new(sink.clone()));
    let mut argv0 = argv![
        "intersect",
        "--mode",
        "any",
        "-a",
        a,
        "-b",
        b,
        "-o",
        p.out.name()
    ];
    argv0.extend(extra0.iter().map(OsString::from));
    let mut argv1 = argv!["merge", "--max-gap", "50", p.inp.name(), "-o", out.name()];
    argv1.extend(extra1.iter().map(OsString::from));
    let specs = vec![
        stage(argv0, vec![p.out], vec![p.fail]),
        stage(argv1, vec![p.inp, out], vec![]),
    ];
    (specs, sink)
}

/// A genome file with one 10 Mb chromosome: `windows -w 1` over it writes 10M lines, far more
/// than a pipe holds.
fn big_genome(dir: &Path) -> PathBuf {
    let g = dir.join("genome.txt");
    std::fs::write(&g, "chr1\t10000000\n").unwrap();
    g
}

/// Sorted BED with `n` records on each of chr1..chr3, from a fixed generator.
fn generated_bed(path: &Path, n: u64, seed: u64, max_gap: u64, max_len: u64, named: bool) {
    let mut x = seed;
    let mut next = |m: u64| {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (x >> 33) % m
    };
    let mut s = String::new();
    for c in ["chr1", "chr2", "chr3"] {
        let mut pos = 0;
        for i in 0..n {
            pos += next(max_gap);
            let end = pos + 1 + next(max_len);
            if named {
                let strand = if next(2) == 0 { '+' } else { '-' };
                s.push_str(&format!(
                    "{c}\t{pos}\t{end}\t{c}_{i}\t{}\t{strand}\n",
                    next(1000)
                ));
            } else {
                s.push_str(&format!("{c}\t{pos}\t{end}\n"));
            }
        }
    }
    std::fs::write(path, s).unwrap();
}

/// The threads of this process that gloci started: stage threads (`gloci-stage-N`), the runs'
/// pool workers (`gloci`), and the threads those start, which inherit the name. Counted by name
/// so that the test harness's own threads, which come and go, never count.
#[cfg(target_os = "linux")]
fn gloci_threads() -> usize {
    std::fs::read_dir("/proc/self/task")
        .unwrap()
        .filter_map(|t| std::fs::read_to_string(t.ok()?.path().join("comm")).ok())
        .filter(|c| c.starts_with("gloci"))
        .count()
}

#[cfg(target_os = "linux")]
fn all_threads() -> usize {
    let s = std::fs::read_to_string("/proc/self/status").unwrap();
    let l = s.lines().find_map(|l| l.strip_prefix("Threads:")).unwrap();
    l.trim().parse().unwrap()
}

/// The message the CLI prints for `gloci ARGS FILE`, which fails with exit status 1.
fn cli_error(args: &[&str], file: &Path) -> String {
    let cli = assert_cmd::Command::cargo_bin("gloci")
        .unwrap()
        .args(args)
        .arg(file)
        .output()
        .unwrap();
    assert_eq!(cli.status.code(), Some(1));
    let stderr = String::from_utf8(cli.stderr).unwrap();
    stderr
        .strip_prefix("error: ")
        .unwrap()
        .trim_end()
        .to_string()
}

#[cfg(unix)]
#[test]
fn chained_stages_equal_shell_pipeline() {
    let _s = serial();
    let (a, b) = (fixture("a.bed"), fixture("b.bed"));
    let (specs, sink) = intersect_merge(&a, &b, &[], &[]);
    assert!(finish(StageSet::spawn(specs)).unwrap());
    let shell = std::process::Command::new("sh")
        .arg("-c")
        .arg(r#""$0" intersect --mode any -a "$1" -b "$2" | "$0" merge --max-gap 50 -"#)
        .arg(assert_cmd::cargo::cargo_bin("gloci"))
        .arg(&a)
        .arg(&b)
        .output()
        .unwrap();
    assert!(
        shell.status.success(),
        "{}",
        String::from_utf8_lossy(&shell.stderr)
    );
    // The merge joins some of the intersected records, so both steps did work.
    let lines = shell.stdout.iter().filter(|&&c| c == b'\n').count();
    assert!(lines > 100 && lines < 740, "{lines} lines");
    assert_eq!(
        String::from_utf8_lossy(&sink.bytes()),
        String::from_utf8_lossy(&shell.stdout)
    );
}

#[test]
fn upstream_data_error_wins() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    // 200,000 good records, more than the pipe holds, then a bad one: the next stage is
    // mid-stream when this one fails.
    let mut text = String::new();
    for i in 0..200_000u64 {
        text.push_str(&format!("chr1\t{}\t{}\n", i * 100, i * 100 + 10));
    }
    text.push_str("chr1\tx\t20\n");
    let bad = dir.path().join("bad.bed");
    std::fs::write(&bad, &text).unwrap();
    let p = link("<merge output>");
    let sink = Shared::default();
    let out = Registered::output("<out>", Box::new(sink.clone()));
    let argv0 = argv!["merge", &bad, "-o", p.out.name()];
    let argv1 = argv!["merge", "--max-gap", "50", p.inp.name(), "-o", out.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![p.out], vec![p.fail]),
        stage(argv1, vec![p.inp, out], vec![]),
    ]);
    let err = finish(set).unwrap_err();
    assert_eq!(err.stage, 0, "{:#}", err.error);
    assert_eq!(err.outcome, Outcome::Failure(cli_error(&["merge"], &bad)));
    let site = err.site.expect("a parse error has a site");
    assert_eq!(site.kind, SiteKind::Parse);
    assert_eq!(site.path, bad.to_str().unwrap());
    assert_eq!(site.loc, Some(Loc::Line(200_001)));

    // The reading stage's error is derived (it read `upstream step failed`), not merely later:
    // with the reader listed first, the writer's error is still the one reported. The writer
    // fails before writing anything, so without the failure the reader (vcf filter) would see
    // an empty input, a usage error of its own.
    let first_bad = dir.path().join("first_bad.bed");
    std::fs::write(&first_bad, "chr1\tx\t20\n").unwrap();
    let p = link("<merge output>");
    let out = Registered::output("<out>", Box::new(Shared::default()));
    let argv0 = argv!["vcf", "filter", p.inp.name(), "-o", out.name()];
    let argv1 = argv!["merge", &first_bad, "-o", p.out.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![p.inp, out], vec![]),
        stage(argv1, vec![p.out], vec![p.fail]),
    ]);
    let err = finish(set).unwrap_err();
    assert_eq!(err.stage, 1, "{:#}", err.error);
    assert_eq!(
        err.outcome,
        Outcome::Failure(cli_error(&["merge"], &first_bad))
    );
}

#[test]
fn downstream_usage_error_beats_upstream_broken_pipe() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let g = big_genome(dir.path());
    let p = link("<windows output>");
    let out = Registered::output("<out>", Box::new(Shared::default()));
    let argv0 = argv!["windows", "-g", &g, "-w", "1", "-o", p.out.name()];
    let argv1 = argv!["merge", "--bogus", p.inp.name(), "-o", out.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![p.out], vec![p.fail]),
        stage(argv1, vec![p.inp, out], vec![]),
    ]);
    // Stage 1 fails at parse and drops its endpoints, so stage 0's writes fail with a broken
    // pipe instead of blocking for ever.
    let err = finish(set).unwrap_err();
    assert_eq!(err.stage, 1, "{:#}", err.error);
    let Outcome::Usage(msg) = &err.outcome else {
        panic!("{:?}", err.outcome)
    };
    assert!(msg.contains("'--bogus'"), "{msg}");
}

#[test]
fn early_usage_error_does_not_hang() {
    let _s = serial();
    // Stage 0 fails at parse, before it takes its output: the next stage, waiting on that pipe,
    // is told, and the set reports stage 0.
    let p = link("<merge output>");
    let out = Registered::output("<out>", Box::new(Shared::default()));
    let argv0 = argv![
        "merge",
        "--bogus",
        fixture("a.bed").as_path(),
        "-o",
        p.out.name()
    ];
    let argv1 = argv!["merge", p.inp.name(), "-o", out.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![p.out], vec![p.fail]),
        stage(argv1, vec![p.inp, out], vec![]),
    ]);
    let err = finish(set).unwrap_err();
    assert_eq!(err.stage, 0, "{:#}", err.error);
    assert!(
        matches!(&err.outcome, Outcome::Usage(m) if m.contains("'--bogus'")),
        "{:?}",
        err.outcome
    );
}

#[test]
fn downstream_early_close_is_quiet() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let g = big_genome(dir.path());
    // Start from no gloci thread at all: an earlier test's pool workers may still be exiting.
    #[cfg(target_os = "linux")]
    let before = {
        let t0 = Instant::now();
        while gloci_threads() > 0 && t0.elapsed() < Duration::from_secs(1) {
            std::thread::sleep(Duration::from_millis(10));
        }
        assert_eq!(
            gloci_threads(),
            0,
            "gloci threads left from an earlier test"
        );
        all_threads()
    };
    let p0 = link("<windows output>");
    let (w1, mut r1) = pipe(PIPE_CHUNKS, PIPE_CHUNK);
    let f1 = w1.fail_handle();
    let o1 = Registered::output("<view output>", Box::new(w1));
    let argv0 = argv!["windows", "-g", &g, "-w", "1", "-o", p0.out.name()];
    let argv1 = argv!["view", p0.inp.name(), "-o", o1.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![p0.out], vec![p0.fail]),
        stage(argv1, vec![p0.inp, o1], vec![f1]),
    ]);
    // The consumer reads the first records, then goes away.
    let mut head = [0u8; 1000];
    r1.read_exact(&mut head).unwrap();
    assert!(head.starts_with(b"chr1\t0\t1\nchr1\t1\t2\n"));
    #[cfg(target_os = "linux")]
    assert!(gloci_threads() > 0, "no stage thread seen");
    drop(r1);
    // Each stage ends with a broken pipe, which is quiet: the set succeeds.
    assert!(finish(set).unwrap());
    #[cfg(target_os = "linux")]
    {
        // Back to none within 1 s: the stage threads are joined, and the runs' pool workers
        // exit on their own once their pools are dropped.
        let t0 = Instant::now();
        while gloci_threads() > 0 && t0.elapsed() < Duration::from_secs(1) {
            std::thread::sleep(Duration::from_millis(10));
        }
        assert_eq!(
            gloci_threads(),
            0,
            "gloci threads left after 1 s (Threads: {} now, {before} before)",
            all_threads()
        );
    }
}

#[test]
fn deterministic_across_threads() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let (a, b) = (dir.path().join("a.bed"), dir.path().join("b.bed"));
    generated_bed(&a, 100_000, 1, 100, 80, true);
    generated_bed(&b, 30_000, 2, 300, 200, false);
    let mut outs = Vec::new();
    for t in ["1", "2", "4"] {
        // collapse keeps every A name, in order, so the merge output shows the whole stream.
        let (specs, sink) = intersect_merge(&a, &b, &["-t", t], &["-t", t, "-c", "4:collapse"]);
        assert!(finish(StageSet::spawn(specs)).unwrap(), "-t {t}");
        outs.push(sink.bytes());
    }
    // More than the pipe holds went through it: intersect's output, the merge's input.
    let piped = assert_cmd::Command::cargo_bin("gloci")
        .unwrap()
        .args(["intersect", "--mode", "any", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .unwrap()
        .stdout
        .len();
    assert!(piped > (PIPE_CHUNKS + 1) * PIPE_CHUNK, "{piped} bytes");
    assert!(!outs[0].is_empty());
    assert!(outs[1] == outs[0], "-t 2 differs from -t 1");
    assert!(outs[2] == outs[0], "-t 4 differs from -t 1");
}

#[test]
fn panic_becomes_internal_error() {
    let _s = serial();
    // A panic in a rayon pool the stage runs (rayon raises it again on the stage's thread).
    let set = StageSet::spawn_fns(vec![
        f(|| Ok(true)),
        f(|| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(2).build()?;
            pool.install(|| {
                rayon::join(
                    || (),
                    || assert!(std::hint::black_box(false), "boom in the pool"),
                )
            });
            Ok(true)
        }),
    ]);
    let err = finish(set).unwrap_err();
    assert_eq!(err.stage, 1);
    let Outcome::Failure(msg) = &err.outcome else {
        panic!("{:?}", err.outcome)
    };
    assert!(msg.contains("boom in the pool"), "{msg}");
    assert!(msg.contains("this is a bug"), "{msg}");
    // And on the stage's own thread.
    let err = finish(StageSet::spawn_fns(vec![f(|| panic!("boom"))])).unwrap_err();
    assert!(
        matches!(&err.outcome, Outcome::Failure(m) if m.contains("boom") && m.contains("this is a bug")),
        "{:?}",
        err.outcome
    );
}

#[test]
fn cancel_unblocks_parked_stages() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let g = big_genome(dir.path());
    let Link { out, inp, fail } = link("<windows output>");
    drop(fail);
    let argv0 = argv!["windows", "-g", &g, "-w", "1", "-o", out.name()];
    let mut set = StageSet::spawn_fns(vec![
        // A real command, parked on the full pipe.
        f(move || {
            let _out = out;
            let cli = embed::parse(&argv0)?;
            let o = ExecOptions {
                command_line: gloci_io::shell_quote_args(&argv0),
                tmp_dir: None,
            };
            embed::execute(&cli, &o)
        }),
        // Holds the pipe's reader and never reads it.
        f(move || {
            let _inp = inp;
            while !cancel::is_cancelled() {
                std::thread::sleep(Duration::from_millis(5));
            }
            Err(gloci_io::IoError::Cancelled.into())
        }),
    ]);
    std::thread::sleep(Duration::from_millis(500));
    assert!(
        set.wait(Duration::ZERO).is_none(),
        "the stages ended before the cancel"
    );
    set.cancel();
    let t0 = Instant::now();
    let err = finish(set).unwrap_err();
    assert!(t0.elapsed() < Duration::from_secs(1), "{:?}", t0.elapsed());
    assert_eq!(err.outcome, Outcome::Cancelled);
    assert_eq!(err.stage, 0);
}

#[test]
fn dropping_a_running_set_cancels_and_joins_it() {
    let _s = serial();
    let dir = tempfile::tempdir().unwrap();
    let g = big_genome(dir.path());
    let p0 = link("<windows output>");
    let (w1, mut r1) = pipe(PIPE_CHUNKS, PIPE_CHUNK);
    let f1 = w1.fail_handle();
    let o1 = Registered::output("<view output>", Box::new(w1));
    let argv0 = argv!["windows", "-g", &g, "-w", "1", "-o", p0.out.name()];
    let argv1 = argv!["view", p0.inp.name(), "-o", o1.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![p0.out], vec![p0.fail]),
        stage(argv1, vec![p0.inp, o1], vec![f1]),
    ]);
    // Nobody reads the last pipe, so both stages park.
    std::thread::sleep(Duration::from_millis(500));
    let t0 = Instant::now();
    drop(set);
    assert!(t0.elapsed() < Duration::from_secs(1), "{:?}", t0.elapsed());
    // The last stage's reader is told the run stopped, not given a clean end.
    let mut rest = Vec::new();
    let e = r1.read_to_end(&mut rest).unwrap_err();
    assert_eq!(e.to_string(), "upstream step failed: interrupted");
}

/// A stage thread has a stack of its own size: with `RUST_MIN_STACK=65536`, a thread spawned
/// without one gets 64 KiB, on which parsing the argv overflows. Runs a test of real commands
/// in a child process of this binary with that setting.
#[test]
fn stage_threads_do_not_take_the_default_stack() {
    let _s = serial();
    let out = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "downstream_early_close_is_quiet",
            "--test-threads",
            "1",
        ])
        .env("RUST_MIN_STACK", "65536")
        .output()
        .unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success() && stdout.contains("1 passed"),
        "{:?}\n{stdout}\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn wait_with_an_unbounded_timeout_blocks_until_the_end() {
    let _s = serial();
    // Duration::MAX does not fit after Instant::now(): wait blocks with no deadline instead.
    let mut set = StageSet::spawn_fns(vec![f(|| {
        std::thread::sleep(Duration::from_millis(100));
        Ok(true)
    })]);
    assert!(set.wait(Duration::MAX).expect("the result").unwrap());
    assert!(
        set.wait(Duration::MAX).is_none(),
        "the result is returned once"
    );
}

/// Whether `msg` shows an in-memory endpoint only by a display name: no raw `\0gloci:` name and
/// no `<memory>` (what a name that was discarded too early shows as).
fn no_raw_names(msg: &str) -> bool {
    !msg.contains('\0') && !msg.contains("<memory>")
}

#[test]
fn errors_show_display_names() {
    let _s = serial();
    // R3: the stage judges its error while its endpoints are still registered.
    let src = Registered::input(
        "<my list>",
        Box::new(Cursor::new(b"chr1\t0\t10\nchr1\tx\t20\n".to_vec())),
    );
    let out = Registered::output("<my sink>", Box::new(Shared::default()));
    let argv = argv!["merge", src.name(), "-o", out.name()];
    let err = finish(StageSet::spawn(vec![stage(argv, vec![src, out], vec![])])).unwrap_err();
    let Outcome::Failure(msg) = &err.outcome else {
        panic!("{:?}", err.outcome)
    };
    assert!(msg.contains("<my list>") && no_raw_names(msg), "{msg:?}");
    let site = err.site.as_ref().expect("a parse error has a site");
    assert_eq!(site.path, "<my list>");

    // Fed through a pipe: the unsorted error of the stage reading it names the pipe.
    let src = Registered::input(
        "<unsorted list>",
        Box::new(Cursor::new(b"chr1\t50\t60\nchr1\t0\t10\n".to_vec())),
    );
    let p = link("<view output>");
    let out = Registered::output("<my sink>", Box::new(Shared::default()));
    let argv0 = argv!["view", src.name(), "-o", p.out.name()];
    let argv1 = argv!["merge", p.inp.name(), "-o", out.name()];
    let set = StageSet::spawn(vec![
        stage(argv0, vec![src, p.out], vec![p.fail]),
        stage(argv1, vec![p.inp, out], vec![]),
    ]);
    let err = finish(set).unwrap_err();
    assert_eq!(err.stage, 1, "{:#}", err.error);
    let Outcome::Failure(msg) = &err.outcome else {
        panic!("{:?}", err.outcome)
    };
    assert!(
        msg.contains("<view output>") && no_raw_names(msg),
        "{msg:?}"
    );
    let site = err.site.as_ref().expect("an unsorted error has a site");
    assert_eq!(site.kind, SiteKind::Unsorted);
    assert_eq!(site.path, "<view output>");
    let order = site.order.as_ref().unwrap();
    assert!(no_raw_names(&order.hint), "{:?}", order.hint);

    // The two errors above already carry display names when they are made. This one holds the
    // raw name until it is judged (-g takes a path or a genome name, so an in-memory name is an
    // unknown genome): judged after the endpoints were discarded, it would say `<memory>`.
    let genome = Registered::input(
        "<my genome>",
        Box::new(Cursor::new(b"chr1\t100\n".to_vec())),
    );
    let out = Registered::output("<my sink>", Box::new(Shared::default()));
    let argv = argv!["windows", "-g", genome.name(), "-w", "10", "-o", out.name()];
    let err = finish(StageSet::spawn(vec![stage(
        argv,
        vec![genome, out],
        vec![],
    )]))
    .unwrap_err();
    assert!(
        format!("{:#}", err.error).contains('\0'),
        "the raw error no longer holds the raw name; find another case for R3"
    );
    let Outcome::Usage(msg) = &err.outcome else {
        panic!("{:?}", err.outcome)
    };
    assert!(
        msg.starts_with("unknown genome \"<my genome>\"") && no_raw_names(msg),
        "{msg:?}"
    );
}

/// A host source that yields one record, then blocks until released (a slow host iterator),
/// or for 5 s at most.
struct Blocking {
    first: bool,
    gate: std::sync::mpsc::Receiver<()>,
    unblocked: Arc<AtomicBool>,
}

impl Read for Blocking {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if std::mem::take(&mut self.first) {
            let l = b"chr1\t0\t10\n";
            out[..l.len()].copy_from_slice(l);
            return Ok(l.len());
        }
        let _ = self.gate.recv_timeout(Duration::from_secs(5));
        self.unblocked.store(true, Ordering::SeqCst);
        Ok(0)
    }
}

#[test]
fn a_cancel_wakes_a_host_reading_the_last_pipe() {
    let _s = serial();
    // Stage 0 is stuck in a host source; the host thread, which has no token, is parked
    // reading the last pipe. The set's token, bound to the pipe through its fail handle, wakes
    // that read: the thread's own token cannot.
    let (release, gate) = std::sync::mpsc::channel();
    let unblocked = Arc::new(AtomicBool::new(false));
    let src = Registered::input(
        "<slow>",
        Box::new(Blocking {
            first: true,
            gate,
            unblocked: unblocked.clone(),
        }),
    );
    let (w, mut r) = pipe(PIPE_CHUNKS, PIPE_CHUNK);
    let fail = w.fail_handle();
    let out = Registered::output("<merge output>", Box::new(w));
    let argv = argv!["merge", src.name(), "-o", out.name()];
    let set = StageSet::spawn(vec![stage(argv, vec![src, out], vec![fail])]);
    assert!(cancel::current().is_none(), "the test thread has no token");
    let (read, waited) = std::thread::scope(|s| {
        s.spawn(|| {
            std::thread::sleep(Duration::from_millis(300));
            set.cancel();
        });
        let t0 = Instant::now();
        let read = r.read(&mut [0u8; 64]);
        (read, t0.elapsed())
    });
    assert!(
        !unblocked.load(Ordering::SeqCst),
        "the read waited for the source"
    );
    let e = read.expect_err("the read ends with the cancel");
    assert!(cancel::is_cancel_error(&e), "{e}");
    assert!(waited < Duration::from_millis(1300), "{waited:?}");
    release.send(()).unwrap();
    let err = finish(set).unwrap_err();
    assert_eq!(err.outcome, Outcome::Cancelled);
}
