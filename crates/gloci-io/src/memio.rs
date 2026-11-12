//! In-memory endpoints (D84). A host registers a reader or writer and gets a name that starts
//! with [`MEM_PREFIX`]; the name goes wherever a path does. No path contains NUL, so a name
//! never collides with a file. An endpoint is read or written once, like stdin, and is never
//! indexed. Messages call it by its display name (`<list>`, `<pyarrow.Table>`, ...).

use crate::error::{IoError, UsageError};
use std::collections::HashMap;
use std::io::{Read, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard};

pub const MEM_PREFIX: &str = "\0gloci:";

enum Endpoint {
    In(Option<Box<dyn Read + Send>>),
    Out(Option<Box<dyn Write + Send>>),
}

struct Entry {
    display: String,
    ep: Endpoint,
}

static NEXT: AtomicU64 = AtomicU64::new(1);
static REGISTRY: Mutex<Option<HashMap<u64, Entry>>> = Mutex::new(None);

fn lock() -> MutexGuard<'static, Option<HashMap<u64, Entry>>> {
    REGISTRY.lock().unwrap_or_else(|p| p.into_inner())
}

fn register(display: &str, ep: Endpoint) -> String {
    let id = NEXT.fetch_add(1, Ordering::Relaxed);
    lock().get_or_insert_with(HashMap::new).insert(
        id,
        Entry {
            display: display.to_string(),
            ep,
        },
    );
    format!("{MEM_PREFIX}{id}")
}

/// Register a source; the returned name opens it wherever a path is read.
pub fn register_input(display: &str, src: Box<dyn Read + Send>) -> String {
    register(display, Endpoint::In(Some(src)))
}

/// Register a sink; the returned name receives what is written to an output path.
pub fn register_output(display: &str, sink: Box<dyn Write + Send>) -> String {
    register(display, Endpoint::Out(Some(sink)))
}

pub fn is_mem_path(path: &str) -> bool {
    path.starts_with(MEM_PREFIX)
}

fn id_of(path: &str) -> Option<u64> {
    path.strip_prefix(MEM_PREFIX)?.parse().ok()
}

/// The name messages use: the registered display name, or `path` itself when it is a file.
pub fn display_name(path: &str) -> String {
    if !is_mem_path(path) {
        return path.to_string();
    }
    id_of(path)
        .and_then(|id| lock().as_ref()?.get(&id).map(|e| e.display.clone()))
        .unwrap_or_else(|| "<memory>".to_string())
}

/// Forget an endpoint (and drop its source or sink, if it was never opened).
pub fn discard(path: &str) {
    if let (Some(id), Some(m)) = (id_of(path), lock().as_mut()) {
        m.remove(&id);
    }
}

fn once(what: &str, display: String) -> IoError {
    UsageError(format!(
        "{display} {what}; an in-memory endpoint is used once, so register a fresh one"
    ))
    .into()
}

/// Take the source named `path` (once). Returns its display name with it.
pub fn take_input(path: &str) -> Result<(String, Box<dyn Read + Send>), IoError> {
    let mut g = lock();
    let Some(e) = id_of(path).and_then(|id| g.as_mut()?.get_mut(&id)) else {
        return Err(UsageError(format!(
            "unknown in-memory input; it was never registered or was discarded ({})",
            path.escape_debug()
        ))
        .into());
    };
    match &mut e.ep {
        Endpoint::In(src) => src
            .take()
            .map(|s| (e.display.clone(), s))
            .ok_or_else(|| once("can be read only once", e.display.clone())),
        Endpoint::Out(_) => {
            Err(UsageError(format!("{} is an output; it cannot be read", e.display)).into())
        }
    }
}

/// Take the sink named `path` (once).
pub fn take_output(path: &str) -> Result<(String, Box<dyn Write + Send>), IoError> {
    let mut g = lock();
    let Some(e) = id_of(path).and_then(|id| g.as_mut()?.get_mut(&id)) else {
        return Err(UsageError(format!(
            "unknown in-memory output; it was never registered or was discarded ({})",
            path.escape_debug()
        ))
        .into());
    };
    match &mut e.ep {
        Endpoint::Out(w) => w
            .take()
            .map(|s| (e.display.clone(), s))
            .ok_or_else(|| once("can be written only once", e.display.clone())),
        Endpoint::In(_) => {
            Err(UsageError(format!("{} is an input; it cannot be written", e.display)).into())
        }
    }
}

/// `text` with every in-memory name in it (`\0gloci:<id>`) replaced by its display name, so a
/// raw name never reaches a message, a log line or a command line.
pub fn scrub(text: &str) -> String {
    if !text.contains(MEM_PREFIX) {
        return text.to_string();
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(i) = rest.find(MEM_PREFIX) {
        out.push_str(&rest[..i]);
        let after = &rest[i + MEM_PREFIX.len()..];
        let digits = after.bytes().take_while(u8::is_ascii_digit).count();
        let name = &rest[i..i + MEM_PREFIX.len() + digits];
        out.push_str(&display_name(name));
        rest = &after[digits..];
    }
    out.push_str(rest);
    out
}

/// The label an input gets in output data (`--names` defaults, partition): the file name
/// without directories, or the display name of an in-memory input (D84).
pub fn label(path: &str) -> String {
    if path == "-" {
        return "stdin".into();
    }
    let shown = display_name(path);
    if is_mem_path(path) {
        return shown;
    }
    std::path::Path::new(&shown)
        .file_name()
        .map_or_else(|| shown.clone(), |f| f.to_string_lossy().into_owned())
}
