//! Atomic file output and the registry of temp paths removed on signals/panics.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard, Once};

/// Temp paths to remove on termination.
struct Registry {
    paths: Vec<PathBuf>,
}

static REGISTRY: Mutex<Registry> = Mutex::new(Registry { paths: Vec::new() });
static COUNTER: AtomicU64 = AtomicU64::new(0);
static INSTALL: Once = Once::new();

fn registry() -> MutexGuard<'static, Registry> {
    REGISTRY.lock().unwrap_or_else(|e| e.into_inner())
}

// INVARIANT (protected by review): every filesystem creation of a temp path
// happens in exactly one critical section that holds REGISTRY, makes one
// syscall and, for registered paths, pushes the path before unlocking. The
// termination path takes the same lock and holds it until `exit`. Keep all such
// sections in this file, never run caller-supplied code under the lock, and
// never call another registry function (or panic) while holding it.

/// Creates `path` (which must not exist) for writing and registers it for
/// removal on termination, as one step.
pub fn create_temp_file(path: &Path) -> io::Result<File> {
    let mut r = registry();
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    r.paths.push(path.to_path_buf());
    Ok(file)
}

/// Creates a uniquely named directory `PREFIXxxxxxx` in `parent` and registers
/// it for removal on termination, as one step. The guard unregisters it.
pub fn create_temp_dir_in(
    parent: &Path,
    prefix: &str,
) -> io::Result<(tempfile::TempDir, TempGuard)> {
    let mut r = registry();
    let dir = tempfile::Builder::new().prefix(prefix).tempdir_in(parent)?;
    let path = dir.path().to_path_buf();
    r.paths.push(path.clone());
    Ok((dir, TempGuard { path }))
}

/// Creates a file inside a registered temp directory. The file needs no
/// registration of its own, but creating it under the lock guarantees none
/// appears while the termination path is removing the directory.
pub fn create_file_in_temp_dir(path: &Path) -> io::Result<File> {
    let _r = registry();
    File::create(path)
}

/// TEST-ONLY hook: runs `f` with the registry lock held, passing the registered
/// paths. Never call it from production code. Calling any registry function,
/// any `create_*` function or `cleanup_registered` inside `f` self-deadlocks.
/// Not part of the API.
#[doc(hidden)]
pub fn __test_with_registry_locked<T>(f: impl FnOnce(&mut Vec<PathBuf>) -> T) -> T {
    f(&mut registry().paths)
}

fn unregister(path: &Path) {
    let mut r = registry();
    if let Some(i) = r.paths.iter().position(|p| p == path) {
        r.paths.swap_remove(i);
    }
}

/// Unregisters its path on drop. Used for spill directories.
#[must_use = "the path is unregistered when the guard drops"]
pub struct TempGuard {
    path: PathBuf,
}

impl Drop for TempGuard {
    fn drop(&mut self) {
        unregister(&self.path);
    }
}

fn remove_paths(paths: Vec<PathBuf>) {
    for p in paths {
        if p.is_dir() {
            let _ = fs::remove_dir_all(&p);
        } else {
            let _ = fs::remove_file(&p);
        }
    }
}

/// Removes every registered file and directory (best effort), holding the
/// registry lock throughout so nothing is created mid-removal.
pub fn cleanup_registered() {
    let mut r = registry();
    remove_paths(std::mem::take(&mut r.paths));
}

/// Termination path: removes everything and exits with `code` while still
/// holding the lock, so nothing can be created (and left behind) between the
/// removal and the exit. Lock holders make a single syscall, so the lock is
/// waited for at most ~2 s; if it is still unavailable (a panic raised while
/// this thread holds it, which is a bug) exit without cleanup rather than hang.
/// Trade-off: if another thread holds the lock longer than that (for example a
/// create on slow network storage), gloci exits without cleaning up.
fn cleanup_and_exit(code: i32) -> ! {
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
    loop {
        match REGISTRY.try_lock() {
            Ok(mut r) => {
                remove_paths(std::mem::take(&mut r.paths));
                std::process::exit(code)
            }
            Err(std::sync::TryLockError::Poisoned(e)) => {
                remove_paths(std::mem::take(&mut e.into_inner().paths));
                std::process::exit(code)
            }
            Err(std::sync::TryLockError::WouldBlock) => {
                if std::time::Instant::now() >= deadline {
                    std::process::exit(code);
                }
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
        }
    }
}

/// The line printed when the program panics, which is always a bug.
pub fn internal_error_line(message: &str) -> String {
    format!("error: internal error: {message}; this is a bug, please report it")
}

thread_local! {
    /// What a panic on this thread means for the user, when it is not a bug of gloci's own
    /// (a library failing on malformed input it was handed): see [`panic_note`].
    static PANIC_NOTE: std::cell::RefCell<Option<String>> = const { std::cell::RefCell::new(None) };
}

/// While the returned guard lives, a panic on this thread is reported as `note` (followed by
/// the panic's message) instead of an internal error; the exit status stays 1.
pub(crate) fn panic_note(note: String) -> PanicNote {
    let prev = PANIC_NOTE.with(|n| n.try_borrow_mut().ok().and_then(|mut n| n.replace(note)));
    PanicNote { prev }
}

/// Restores the note before it on drop.
pub(crate) struct PanicNote {
    prev: Option<String>,
}

impl Drop for PanicNote {
    fn drop(&mut self) {
        let prev = self.prev.take();
        PANIC_NOTE.with(|n| {
            if let Ok(mut n) = n.try_borrow_mut() {
                *n = prev;
            }
        });
    }
}

/// The line the panic hook prints for `message`: the thread's note, or an internal error.
fn panic_line(message: &str) -> String {
    let note = PANIC_NOTE
        .try_with(|n| n.try_borrow().ok().and_then(|n| n.clone()))
        .ok()
        .flatten();
    match note {
        Some(note) => format!("error: {note} ({message})"),
        None => internal_error_line(message),
    }
}

fn panic_message(info: &std::panic::PanicHookInfo<'_>) -> String {
    let payload = info.payload();
    payload
        .downcast_ref::<&str>()
        .map(|s| (*s).to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "a panic with no message".to_string())
}

/// Installs SIGINT/SIGTERM/SIGHUP handling (cleanup, exit 130) and a panic hook
/// that cleans up, reports an internal error and exits with status 1 (the
/// release profile aborts on panic, which would otherwise end in SIGABRT).
/// Safe to call more than once.
pub fn install_cleanup_handlers() {
    INSTALL.call_once(|| {
        let _ = ctrlc::set_handler(|| {
            cleanup_and_exit(130);
        });
        std::panic::set_hook(Box::new(|info| {
            let _ = writeln!(io::stderr(), "{}", panic_line(&panic_message(info)));
            cleanup_and_exit(1);
        }));
    });
}

fn temp_suffix() -> String {
    const ALPHABET: &[u8; 36] = b"abcdefghijklmnopqrstuvwxyz0123456789";
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let mut x =
        (u64::from(std::process::id()) << 32) ^ nanos ^ n.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    // splitmix64 finalizer
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    (0..6)
        .map(|_| {
            let c = ALPHABET[(x % 36) as usize] as char;
            x /= 36;
            c
        })
        .collect()
}

/// A file written under `DEST.tmp-XXXXXX` and renamed onto `DEST` by [`AtomicFile::commit`].
/// Dropping without committing removes the temp file. It captures the run's token when created,
/// and its commit records `DEST` as an output of the run ([`crate::cancel::CancelToken::recording`]).
pub struct AtomicFile {
    file: Option<File>,
    tmp: PathBuf,
    dest: PathBuf,
    run: crate::cancel::Watch,
}

impl AtomicFile {
    pub fn create(dest: &Path) -> io::Result<AtomicFile> {
        let name = dest.file_name().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "output path has no file name")
        })?;
        let mut last_err = None;
        for _ in 0..100 {
            let mut tname: OsString = name.to_owned();
            tname.push(format!(".tmp-{}", temp_suffix()));
            let tmp = dest.with_file_name(tname);
            match create_temp_file(&tmp) {
                Ok(file) => {
                    return Ok(AtomicFile {
                        file: Some(file),
                        tmp,
                        dest: dest.to_path_buf(),
                        run: crate::cancel::Watch::current(),
                    });
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => last_err = Some(e),
                Err(e) => return Err(e),
            }
        }
        Err(last_err.unwrap_or_else(|| io::Error::other("cannot create temp file")))
    }

    /// Path of the temporary file.
    pub fn temp_path(&self) -> &Path {
        &self.tmp
    }

    /// Flushes and renames the temp file onto the destination, and records the destination as
    /// an output of the run.
    pub fn commit(self) -> io::Result<()> {
        let (run, dest) = self.commit_unrecorded()?;
        run.record_output(&dest);
        Ok(())
    }

    /// [`AtomicFile::commit`] without the record: for [`crate::StagedSet::publish`], which
    /// records its files once all of them are in place. Returns what to record it with.
    pub(crate) fn commit_unrecorded(mut self) -> io::Result<(crate::cancel::Watch, PathBuf)> {
        let file = self.file.take().expect("file present until commit");
        file.sync_all()?;
        drop(file);
        fs::rename(&self.tmp, &self.dest)?;
        unregister(&self.tmp);
        // Disarm Drop: tmp no longer exists.
        self.tmp = PathBuf::new();
        Ok((self.run.clone(), std::mem::take(&mut self.dest)))
    }
}

impl Write for AtomicFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.file.as_mut().expect("file open").write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.file.as_mut().expect("file open").flush()
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        if !self.tmp.as_os_str().is_empty() {
            self.file.take();
            let _ = fs::remove_file(&self.tmp);
            unregister(&self.tmp);
        }
    }
}
