//! Atomic file output and the registry of temp paths removed on signals/panics.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard, Once};

static REGISTRY: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());
static COUNTER: AtomicU64 = AtomicU64::new(0);
static INSTALL: Once = Once::new();

fn registry() -> MutexGuard<'static, Vec<PathBuf>> {
    REGISTRY.lock().unwrap_or_else(|e| e.into_inner())
}

fn register(path: &Path) {
    registry().push(path.to_path_buf());
}

fn unregister(path: &Path) {
    let mut r = registry();
    if let Some(i) = r.iter().position(|p| p == path) {
        r.swap_remove(i);
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

/// Registers `path` (file or directory) for removal by [`cleanup_registered`].
pub fn register_temp(path: &Path) -> TempGuard {
    register(path);
    TempGuard {
        path: path.to_path_buf(),
    }
}

/// Removes every registered file and directory (best effort).
pub fn cleanup_registered() {
    let paths: Vec<PathBuf> = std::mem::take(&mut *registry());
    for p in paths {
        if p.is_dir() {
            let _ = fs::remove_dir_all(&p);
        } else {
            let _ = fs::remove_file(&p);
        }
    }
}

/// The line printed when the program panics, which is always a bug.
pub fn internal_error_line(message: &str) -> String {
    format!("error: internal error: {message}; this is a bug, please report it")
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
            cleanup_registered();
            std::process::exit(130);
        });
        std::panic::set_hook(Box::new(|info| {
            cleanup_registered();
            let _ = writeln!(
                io::stderr(),
                "{}",
                internal_error_line(&panic_message(info))
            );
            std::process::exit(1);
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
/// Dropping without committing removes the temp file.
pub struct AtomicFile {
    file: Option<File>,
    tmp: PathBuf,
    dest: PathBuf,
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
            match OpenOptions::new().write(true).create_new(true).open(&tmp) {
                Ok(file) => {
                    register(&tmp);
                    return Ok(AtomicFile {
                        file: Some(file),
                        tmp,
                        dest: dest.to_path_buf(),
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

    /// Flushes and renames the temp file onto the destination.
    pub fn commit(mut self) -> io::Result<()> {
        let file = self.file.take().expect("file present until commit");
        file.sync_all()?;
        drop(file);
        fs::rename(&self.tmp, &self.dest)?;
        unregister(&self.tmp);
        // Disarm Drop: tmp no longer exists.
        self.tmp = PathBuf::new();
        Ok(())
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
