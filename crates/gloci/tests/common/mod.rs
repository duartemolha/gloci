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
