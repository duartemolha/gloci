//! Minimal stderr logger: `warning: ...` / `info: ...` / `debug: ...`.
use log::{Level, LevelFilter, Log, Metadata, Record};

struct Logger;

impl Log for Logger {
    fn enabled(&self, _: &Metadata) -> bool {
        true
    }
    fn log(&self, r: &Record) {
        if !self.enabled(r.metadata()) {
            return;
        }
        let label = match r.level() {
            Level::Error => "error",
            Level::Warn => "warning",
            Level::Info => "info",
            Level::Debug | Level::Trace => "debug",
        };
        eprintln!("{label}: {}", gloci_io::memio::scrub(&r.args().to_string()));
    }
    fn flush(&self) {}
}

/// Level for the flags: `-q` shows errors only, default warnings, `-v` info, `-vv` debug.
pub fn level_for(quiet: bool, verbose: u8) -> LevelFilter {
    if quiet {
        LevelFilter::Error
    } else {
        match verbose {
            0 => LevelFilter::Warn,
            1 => LevelFilter::Info,
            _ => LevelFilter::Debug,
        }
    }
}

pub fn init(quiet: bool, verbose: u8) {
    // Ignore "already set" so repeated `run` calls (tests) are fine.
    let _ = log::set_logger(&Logger);
    log::set_max_level(level_for(quiet, verbose));
}
