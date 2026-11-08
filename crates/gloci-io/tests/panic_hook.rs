//! The panic hook cleans up, reports an internal error and exits with 1.
use gloci_io::{AtomicFile, install_cleanup_handlers, internal_error_line};
use std::io::Write;
use std::process::Command;

#[test]
fn internal_error_line_format() {
    assert_eq!(
        internal_error_line("boom"),
        "error: internal error: boom; this is a bug, please report it"
    );
}

/// Runs only inside the child process started by `panic_exits_with_status_one`.
#[test]
fn panic_hook_child() {
    let Ok(dir) = std::env::var("GLOCI_PANIC_HOOK_DIR") else {
        return;
    };
    install_cleanup_handlers();
    let mut f = AtomicFile::create(&std::path::Path::new(&dir).join("out.bed")).unwrap();
    f.write_all(b"partial").unwrap();
    panic!("something impossible happened");
}

#[test]
fn panic_exits_with_status_one_after_cleanup() {
    let dir = tempfile::tempdir().unwrap();
    let out = Command::new(std::env::current_exe().unwrap())
        .args(["panic_hook_child", "--exact", "--nocapture", "--test-threads=1"])
        .env("GLOCI_PANIC_HOOK_DIR", dir.path())
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1), "{out:?}");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains(
            "error: internal error: something impossible happened; this is a bug, please report it"
        ),
        "{err}"
    );
    assert!(!err.contains("panicked at"), "default panic output leaked: {err}");
    // the temp file is gone and nothing was published
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn internal_io_error_wording() {
    let e = gloci_io::IoError::Internal("an empty sort chunk rejected a record".into());
    assert_eq!(
        e.to_string(),
        "internal error: an empty sort chunk rejected a record; this is a bug, please report it"
    );
}
