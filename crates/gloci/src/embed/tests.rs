//! Unit tests of the error classification and the run's thread pool.
use super::*;
use gloci_io::memio;
use std::io::{Error, ErrorKind};

fn pipe() -> Error {
    Error::from(ErrorKind::BrokenPipe)
}

/// An error whose `source()` is a `UsageError`.
#[derive(Debug)]
struct HasUsage(UsageError);
impl std::fmt::Display for HasUsage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "has usage")
    }
}
impl std::error::Error for HasUsage {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.0)
    }
}

#[test]
fn usage_error_beats_broken_pipe_in_same_chain() {
    // chain: io::Error(BrokenPipe) -> UsageError (a custom io::Error
    // forwards `source()` of its payload), so both are in one chain.
    let io = Error::new(
        ErrorKind::BrokenPipe,
        HasUsage(UsageError("bad flag".into())),
    );
    let e = anyhow::Error::new(io);
    let kinds: Vec<_> = e
        .chain()
        .map(|c| c.downcast_ref::<UsageError>().is_some())
        .collect();
    assert!(kinds.contains(&true), "fixture must hold both: {kinds:?}");
    assert!(matches!(classify(&e), Outcome::Usage(m) if m == "bad flag"));
    let wrapped = anyhow::Error::new(IoError::Usage(UsageError("u".into()))).context("ctx");
    assert!(matches!(classify(&wrapped), Outcome::Usage(_)));
}

#[test]
fn a_cancelled_write_or_open_is_cancelled_not_a_write_error() {
    let w = anyhow::Error::new(IoError::Write {
        path: "out".into(),
        source: gloci_io::cancel::io_error(),
    });
    assert!(matches!(classify(&w), Outcome::Cancelled));
    let bare = anyhow::Error::new(gloci_io::cancel::io_error()).context("while sorting");
    assert!(matches!(classify(&bare), Outcome::Cancelled));
    // Under a cancelled token any failure is the cancellation.
    let token = gloci_io::cancel::RunToken::new();
    let _g = gloci_io::cancel::RunGuard::enter(&token);
    token.cancel();
    let other = anyhow::Error::new(Error::from(ErrorKind::NotFound));
    assert!(matches!(classify(&other), Outcome::Cancelled));
}

#[test]
fn write_broken_pipe_is_quiet_success() {
    let e = anyhow::Error::new(IoError::Write {
        path: "out".into(),
        source: pipe(),
    });
    assert!(matches!(classify(&e), Outcome::BrokenPipe));
}

#[test]
fn bare_broken_pipe_is_quiet_success() {
    assert!(matches!(
        classify(&anyhow::Error::new(pipe())),
        Outcome::BrokenPipe
    ));
    let other = anyhow::Error::new(Error::from(ErrorKind::NotFound));
    assert!(matches!(classify(&other), Outcome::Failure(_)));
}

/// About `n` x 256 KiB of stack.
fn deep(n: u32) -> u64 {
    let buf = [n as u8; 256 << 10];
    let x = u64::from(std::hint::black_box(&buf)[buf.len() - 1]);
    if n == 0 { x } else { x + deep(n - 1) }
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

#[test]
fn the_command_runs_on_the_worker_stack() {
    // 8 MiB of frames: more than a test thread's 2 MiB, well within the 16 MiB.
    let none = Watch::default();
    assert_eq!(
        in_run_pool(threads(1), &none, || deep(32)).unwrap(),
        (0..=32).sum()
    );
    assert_eq!(
        in_run_pool(threads(2), &none, || std::thread::current()
            .name()
            .map(String::from))
        .unwrap(),
        Some("gloci".into())
    );
    // Every worker has it, not only the one running the command.
    let all = in_run_pool(threads(3), &none, || rayon::broadcast(|_| deep(32))).unwrap();
    assert_eq!(all, vec![(0..=32).sum::<u64>(); 3]);
}

#[test]
fn every_worker_of_the_run_pool_carries_the_run_token() {
    let token = RunToken::new();
    let watch = {
        let _g = RunGuard::enter(&token);
        Watch::current()
    };
    // Not on this thread: the guard is gone.
    assert!(cancel::current().is_none());
    let before = in_run_pool(threads(4), &watch, || {
        rayon::broadcast(|_| cancel::current().is_some() && !cancel::is_cancelled())
    })
    .unwrap();
    assert_eq!(before, vec![true; 4]);
    token.cancel();
    let after = in_run_pool(threads(4), &watch, || {
        // rayon work the command starts sees the cancellation on whichever worker runs it
        let mut v = vec![false; 64];
        rayon::scope(|s| {
            for slot in v.iter_mut() {
                s.spawn(move |_| *slot = cancel::is_cancelled());
            }
        });
        v
    })
    .unwrap();
    assert!(after.iter().all(|&c| c));
    // Another run's pool is untouched.
    let other = in_run_pool(threads(2), &Watch::default(), || {
        rayon::broadcast(|_| cancel::is_cancelled())
    })
    .unwrap();
    assert_eq!(other, vec![false; 2]);
}

#[test]
fn classify_under_judges_by_the_given_token_not_the_threads() {
    let token = RunToken::new();
    let other = anyhow::Error::new(Error::from(ErrorKind::NotFound));
    assert!(matches!(
        classify_under(&other, Some(&token)),
        Outcome::Failure(_)
    ));
    token.cancel();
    // No guard on this thread: only the token passed in says it was cancelled.
    assert!(matches!(classify(&other), Outcome::Failure(_)));
    assert_eq!(classify_under(&other, Some(&token)), Outcome::Cancelled);
    assert!(matches!(classify_under(&other, None), Outcome::Failure(_)));
}

#[test]
fn error_sites_show_display_names() {
    let name = memio::register_input("<list>", Box::new(std::io::empty()));
    let parse = anyhow::Error::new(IoError::Parse {
        path: name.clone(),
        loc: Loc::Line(3),
        msg: "bad".into(),
    })
    .context("while intersecting");
    let site = error_site(&parse).unwrap();
    assert_eq!(site.kind, SiteKind::Parse);
    assert_eq!(site.path, "<list>");
    assert_eq!(site.loc, Some(Loc::Line(3)));
    assert_eq!(site.order, None);
    let order = anyhow::Error::new(IoError::from(OrderError {
        path: name.clone(),
        prev_loc: Loc::Line(2),
        prev: "chr2\t1\t2".into(),
        loc: Loc::Line(3),
        rec: "chr1\t1\t2".into(),
        hint: format!("sort it: gloci sort {name}"),
    }));
    let site = error_site(&order).unwrap();
    assert_eq!(site.kind, SiteKind::Unsorted);
    assert_eq!(site.path, "<list>");
    let o = site.order.unwrap();
    assert_eq!((o.prev_loc, o.loc), (Loc::Line(2), Loc::Line(3)));
    assert_eq!(o.hint, "sort it: gloci sort <list>");
    assert!(!format!("{o:?}").contains('\0'));
    let write = anyhow::Error::new(IoError::Write {
        path: name.clone(),
        source: Error::from(ErrorKind::Other),
    });
    assert_eq!(error_site(&write).unwrap().kind, SiteKind::Write);
    assert_eq!(error_site(&anyhow::anyhow!("plain")), None);
    let usage = anyhow::Error::new(IoError::Usage(UsageError("u".into())));
    assert_eq!(error_site(&usage), None);
    memio::discard(&name);
}

#[test]
fn messages_show_display_names() {
    let name = memio::register_input("<generator>", Box::new(std::io::empty()));
    let e = anyhow::Error::new(IoError::Truncated { path: name.clone() });
    let Outcome::Failure(msg) = classify(&e) else {
        panic!("not a failure")
    };
    assert!(msg.starts_with("<generator>: BGZF"), "{msg:?}");
    let u = anyhow::Error::new(UsageError(format!("{name} is in memory")));
    assert_eq!(
        classify(&u),
        Outcome::Usage("<generator> is in memory".into())
    );
    memio::discard(&name);
}

#[test]
fn parse_errors_are_plain_text() {
    let argv: Vec<OsString> = ["gloci", "merge", "--max-gap", "x", "a.bed"]
        .into_iter()
        .map(OsString::from)
        .collect();
    let e = parse(&argv).err().unwrap();
    assert!(
        e.0.starts_with("invalid value 'x' for '--max-gap <N>'"),
        "{e}"
    );
    assert!(
        !e.0.contains("For more information") && !e.0.contains('\x1b'),
        "{e}"
    );
    assert!(!e.0.ends_with('\n'));
}

fn words(w: &[&str]) -> Vec<OsString> {
    w.iter().map(OsString::from).collect()
}

#[test]
fn a_needless_separator_is_dropped_from_the_shell_form() {
    // Plain: step_argv's `--` before the input, after an option given its one value.
    let merge = binding(&["merge"]).unwrap();
    let kw = [("max_gap".to_string(), KwValue::Value("50".into()))];
    let argv = step_argv(merge, &words(&["-"]), &kw, &[]).unwrap();
    assert_eq!(
        argv,
        words(&["gloci", "merge", "--max-gap", "50", "--", "-"])
    );
    assert_eq!(
        drop_needless_separator(&argv),
        words(&["gloci", "merge", "--max-gap", "50", "-"])
    );
    // A value attached because it starts with `-`, a switch, and a two-word command.
    for (with, without) in [
        (
            &["gloci", "merge", "--max-gap=-2", "--", "a.bed"][..],
            &["gloci", "merge", "--max-gap=-2", "a.bed"][..],
        ),
        (
            &["gloci", "merge", "-s", "--", "a.bed", "b.bed"],
            &["gloci", "merge", "-s", "a.bed", "b.bed"],
        ),
        (
            &["gloci", "vcf", "filter", "--include", "QUAL>10", "--", "-"],
            &["gloci", "vcf", "filter", "--include", "QUAL>10", "-"],
        ),
    ] {
        assert_eq!(drop_needless_separator(&words(with)), words(without));
    }
}

#[test]
fn a_needed_separator_is_kept() {
    for argv in [
        // A positional that looks like a flag.
        &["gloci", "merge", "--", "a.bed", "-weird.bed"][..],
        // The option before it could take more values (`-b x y`), and would swallow them.
        &[
            "gloci",
            "maskfasta",
            "--fasta",
            "ref.fa",
            "-b",
            "m.bed",
            "--",
            "x.bed",
        ],
    ] {
        assert_eq!(drop_needless_separator(&words(argv)), words(argv));
    }
    // No separator: as it is.
    let plain = words(&["gloci", "merge", "a.bed"]);
    assert_eq!(drop_needless_separator(&plain), plain);
}
