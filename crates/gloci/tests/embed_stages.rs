//! The stage runtime (D85): bounded in-memory pipes, and several commands running at once on
//! threads, chained by those pipes, with one error reported for the set and cancellation.
use gloci::embed::{PIPE_CHUNK, PIPE_CHUNKS, pipe};
use gloci_io::cancel::{self, CancelToken, RunGuard};
use std::io::{ErrorKind, Read, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::Duration;

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
    // Wait until the writer stops making progress.
    let mut last = usize::MAX;
    loop {
        std::thread::sleep(Duration::from_millis(200));
        let now = accepted.load(Ordering::SeqCst);
        if now == last || writer.is_finished() {
            break;
        }
        last = now;
    }
    let held = accepted.load(Ordering::SeqCst);
    let full = PIPE_CHUNKS * PIPE_CHUNK;
    assert!(
        !writer.is_finished(),
        "the writer never blocked ({held} bytes accepted)"
    );
    // At most the channel's chunks, plus the one partial buffer being filled, are unread.
    assert!(
        held >= full - 1000 && held <= full + PIPE_CHUNK,
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
