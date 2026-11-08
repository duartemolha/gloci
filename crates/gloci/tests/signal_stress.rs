//! Stress test for the signal-cleanup race (task 26a). Sends SIGTERM to a
//! spilling `gloci sort` at pseudo-random 0-12 ms delays and checks that no
//! temp output or spill directory survives. Before the registry-lock fix this
//! left a file in about 3 of 400 runs under CPU load.
//!
//! Run: `cargo test -p gloci --test signal_stress -- --ignored --nocapture`
//! (`GLOCI_STRESS_RUNS=N` sets the run count, default 400; add CPU load, such
//! as a few busy loops, to widen the race window).
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;

fn leftovers(root: &Path, keep: &[&str]) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            let name = p.file_name().unwrap().to_string_lossy().into_owned();
            if keep.contains(&name.as_str()) {
                continue;
            }
            if p.is_dir() && name != "spill" {
                stack.push(p.clone());
            }
            out.push(p);
        }
    }
    out
}

#[test]
#[ignore = "stress loop; run explicitly"]
fn sigterm_at_random_moments_leaves_nothing() {
    let runs: usize = std::env::var("GLOCI_STRESS_RUNS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(400);
    let d = tempfile::tempdir().unwrap();
    let big = d.path().join("big.bed");
    let mut text = String::new();
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in 0..3_000_000u64 {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let start = (x >> 33) % 5_000_000;
        writeln!(
            text,
            "chr{}\t{start}\t{}\tr{i}",
            1 + (x >> 20) % 12,
            start + 100
        )
        .unwrap();
    }
    std::fs::write(&big, text).unwrap();
    let bin = assert_cmd::cargo::cargo_bin("gloci");
    let mut bad = Vec::new();
    for run in 0..runs {
        let spill = d.path().join("spill");
        let _ = std::fs::remove_dir_all(&spill);
        std::fs::create_dir(&spill).unwrap();
        let out = d.path().join("out.bed");
        let mut c = Command::new(&bin)
            .arg("sort")
            .arg(&big)
            .args(["--max-mem", "16M", "--tmp-dir"])
            .arg(&spill)
            .arg("-o")
            .arg(&out)
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
        std::thread::sleep(Duration::from_micros((x >> 40) % 12_000)); // the stimulus
        unsafe { libc::kill(c.id() as i32, libc::SIGTERM) };
        let st = c.wait().unwrap();
        // Killed before the handler was installed (143) or a clean 130.
        let left = leftovers(d.path(), &["big.bed", "spill"]);
        if !left.is_empty() || !matches!(st.code(), Some(130) | None) || out.exists() {
            bad.push(format!("run {run}: status {st:?}, leftovers {left:?}"));
        }
    }
    assert!(bad.is_empty(), "{} of {runs} runs: {bad:#?}", bad.len());
}
