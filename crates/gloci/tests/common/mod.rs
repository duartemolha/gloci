#![allow(dead_code)]
use std::path::PathBuf;

pub fn gloci() -> assert_cmd::Command {
    assert_cmd::Command::cargo_bin("gloci").unwrap()
}

pub fn tempdir() -> tempfile::TempDir {
    tempfile::tempdir().unwrap()
}

pub fn bgzf_file(dir: &std::path::Path, name: &str, text: &str) -> PathBuf {
    gloci_testkit::fixtures::write(dir, name, &gloci_testkit::fixtures::bgzip(text.as_bytes()))
}

pub fn plain_file(dir: &std::path::Path, name: &str, text: &str) -> PathBuf {
    gloci_testkit::fixtures::write(dir, name, text.as_bytes())
}

/// Run `gloci ARGS` through bash, so `<(…)` process substitutions work,
/// with `stdin` on its standard input.
#[cfg(unix)]
pub fn bash_gloci(args: &str, stdin: &str) -> std::process::Output {
    use std::io::Write;
    let bin = assert_cmd::cargo::cargo_bin("gloci");
    let mut child = std::process::Command::new("bash")
        .arg("-c")
        .arg(format!("'{}' {args}", bin.display()))
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let mut sin = child.stdin.take().unwrap();
    let data = stdin.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = sin.write_all(&data);
    });
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    out
}

/// Copy `data` and its `index` (a `.bai`/`.csi` named `<data name><ext>`) into `dir`, and make the
/// index newer than the data by an hour, so that a checkout's file times never decide whether
/// the index looks stale. Returns the copy of `data`.
pub fn copy_with_fresh_index(
    dir: &std::path::Path,
    data: &std::path::Path,
    index: &std::path::Path,
) -> PathBuf {
    let name = data.file_name().unwrap().to_string_lossy().into_owned();
    let ext = index
        .file_name()
        .unwrap()
        .to_string_lossy()
        .strip_prefix(&name)
        .expect("the index is named after the data file")
        .to_string();
    let d = dir.join(&name);
    std::fs::copy(data, &d).unwrap();
    let i = dir.join(format!("{name}{ext}"));
    std::fs::copy(index, &i).unwrap();
    std::fs::File::options()
        .write(true)
        .open(&i)
        .unwrap()
        .set_modified(std::time::SystemTime::now() + std::time::Duration::from_secs(3600))
        .unwrap();
    d
}
