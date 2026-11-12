//! A raw in-memory name (`\0gloci:N`, D84) never reaches output data, messages or command
//! lines: it is shown by its display name everywhere.
use gloci_io::memio::{display_name, register_input, register_output, scrub};
use std::ffi::OsString;
use std::io::{Cursor, Write};
use std::sync::{Arc, Mutex, MutexGuard};

static SERIAL: Mutex<()> = Mutex::new(());
fn serial() -> MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

fn args(a: &[&str]) -> Vec<OsString> {
    std::iter::once("gloci")
        .chain(a.iter().copied())
        .map(OsString::from)
        .collect()
}

fn mem(display: &str, data: &str) -> String {
    register_input(display, Box::new(Cursor::new(data.as_bytes().to_vec())))
}

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

const BED: &str = "chr1\t1\t10\n";

#[test]
fn intersect_labels_use_the_display_basename() {
    let _s = serial();
    let a = mem("<a>", BED);
    let b1 = mem("/some/dir/first.bed", BED);
    let b2 = mem("<second>", BED);
    let sink = Shared::default();
    let out = register_output("<out>", Box::new(sink.clone()));
    let (code, msg) = gloci::run_reporting(args(&[
        "intersect",
        "-a",
        &a,
        "-b",
        &b1,
        "-b",
        &b2,
        "--mode",
        "pairs",
        "-o",
        &out,
        "-O",
        "tsv",
    ]));
    assert_eq!(msg, None);
    assert_eq!(
        format!("{code:?}"),
        format!("{:?}", std::process::ExitCode::SUCCESS)
    );
    let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
    assert!(!text.contains('\0'), "{text:?}");
    assert!(
        text.contains("\tfirst.bed\t") && text.contains("\t<second>\t"),
        "{text}"
    );
}

#[test]
fn partition_labels_use_the_display_basename() {
    let _s = serial();
    let a = mem("/x/y/one.bed", BED);
    let b = mem("<two>", BED);
    let sink = Shared::default();
    let out = register_output("<out>", Box::new(sink.clone()));
    let (_, msg) = gloci::run_reporting(args(&["partition", &a, &b, "-o", &out, "-O", "tsv"]));
    assert_eq!(msg, None);
    let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
    assert!(
        !text.contains('\0') && text.contains("one.bed,<two>"),
        "{text:?}"
    );
}

#[test]
fn messages_and_headers_show_display_names() {
    let _s = serial();
    // index of a mem name is a usage error naming the display name.
    let f = mem("<list>", BED);
    let (code, msg) = gloci::run_reporting(args(&["index", &f]));
    assert_eq!(
        format!("{code:?}"),
        format!("{:?}", std::process::ExitCode::from(2))
    );
    assert_eq!(
        msg.as_deref(),
        Some("<list> is in memory; index needs a file")
    );
    // sort -n of a non-alignment input says its display name.
    let f = mem("<generator>", BED);
    let (_, msg) = gloci::run_reporting(args(&["sort", "-n", &f]));
    let msg = msg.unwrap();
    assert!(
        msg.contains("<generator>") && !msg.contains('\0'),
        "{msg:?}"
    );
    // The command line of a header holds no NUL.
    let f = mem(
        "<list>",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t5\t.\t.\n",
    );
    let sink = Shared::default();
    let out = register_output("<out>", Box::new(sink.clone()));
    let (_, msg) = gloci::run_reporting(args(&["vcf", "filter", &f, "-o", &out, "-O", "vcf"]));
    assert_eq!(msg, None);
    let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
    assert!(!text.contains('\0'), "{text:?}");
    assert!(
        text.contains("##glociCommand") && text.contains("<list>"),
        "{text}"
    );
}

#[test]
fn scrub_replaces_every_name() {
    let _s = serial();
    let a = mem("<a>", "");
    let b = mem("<b>", "");
    assert_eq!(display_name(&a), "<a>");
    assert_eq!(scrub(&format!("x {a} y {b} z")), "x <a> y <b> z");
    assert_eq!(scrub("plain"), "plain");
}
