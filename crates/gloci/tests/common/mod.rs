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

/// `xtask/diff/inputs`, the fixtures the differential cases use.
pub fn diff_inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

/// A FASTA in `dir` with the names and lengths of `fa/ref.fa` (and its `.fai`) but every base
/// an `A`: a CRAM made with ref.fa decodes against it only where no base is read.
pub fn wrong_bases_fasta(dir: &std::path::Path) -> PathBuf {
    let fa = std::fs::read_to_string(diff_inputs().join("fa/ref.fa")).unwrap();
    let other: String = fa
        .lines()
        .map(|l| {
            if l.starts_with('>') {
                format!("{l}\n")
            } else {
                format!("{}\n", "A".repeat(l.len()))
            }
        })
        .collect();
    let wrong = dir.join("other.fa");
    std::fs::write(&wrong, other).unwrap();
    std::fs::copy(
        diff_inputs().join("fa/ref.fa.fai"),
        dir.join("other.fa.fai"),
    )
    .unwrap();
    wrong
}

/// Every CRAM fixture that `gloci view` decodes in full with `fa/ref.fa`, with the SAM text of
/// that decode: its header lines and its records. The fixtures that are meant to fail (a
/// missing end-of-file container, a read past the overhang limit, another reference) are left
/// out, so what is returned is what the full decode says each file holds.
pub fn decodable_crams() -> Vec<(PathBuf, Vec<String>, Vec<String>)> {
    let reference = diff_inputs().join("fa/ref.fa");
    let mut paths: Vec<PathBuf> = std::fs::read_dir(diff_inputs().join("aln"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "cram"))
        .collect();
    paths.sort();
    let mut out = Vec::new();
    for p in paths {
        let o = gloci()
            .args(["view", "-O", "sam", "--reference"])
            .arg(&reference)
            .arg(&p)
            .output()
            .unwrap();
        if !o.status.success() {
            continue;
        }
        let text = String::from_utf8(o.stdout).unwrap();
        let (header, records): (Vec<String>, Vec<String>) = text
            .lines()
            .map(String::from)
            .partition(|l| l.starts_with('@'));
        out.push((p, header, records));
    }
    assert!(out.len() >= 20, "only {} CRAM fixtures decode", out.len());
    out
}
