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
