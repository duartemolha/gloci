//! The Python stub (`__init__.pyi`) is generated from the Rust specs, so it cannot drift from
//! what the runtime accepts.
use std::fs;

fn stub_path() -> std::path::PathBuf {
    xtask::diff::workspace_root().join("crates/gloci-py/python/gloci/__init__.pyi")
}

#[test]
fn pyi_is_current() {
    let d = tempfile::tempdir().unwrap();
    let out = d.path().join("__init__.pyi");
    xtask::pystub::generate_pyi(&out).unwrap();
    let generated = fs::read_to_string(&out).unwrap();
    let committed = fs::read_to_string(stub_path()).unwrap_or_default();
    assert_eq!(
        generated, committed,
        "__init__.pyi is out of date; run: cargo xtask pystub"
    );
}

#[test]
fn every_step_and_function_is_in_the_stub() {
    let d = tempfile::tempdir().unwrap();
    let out = d.path().join("__init__.pyi");
    xtask::pystub::generate_pyi(&out).unwrap();
    let text = fs::read_to_string(&out).unwrap();
    for b in gloci::embed::BINDINGS {
        assert!(
            text.contains(&format!("    def {}(", b.method)),
            "no step method {}",
            b.method
        );
    }
    for (_, f) in gloci::embed::MODULE_FUNCTIONS {
        assert!(text.contains(&format!("\ndef {f}(")), "no function {f}");
    }
    // The reserved keyword arguments never appear on a step.
    let step = text
        .split("    def merge(")
        .nth(1)
        .and_then(|s| s.split("\n    def ").next())
        .expect("merge in the stub");
    for reserved in [
        "output:",
        "threads:",
        "quiet:",
        "verbose:",
        "output_format:",
    ] {
        assert!(!step.contains(reserved), "merge() has {reserved}");
    }
    assert!(text.contains("Source: TypeAlias ="), "the Source alias");
    assert!(text.contains("-> pyarrow.Table"), "to_arrow's type");
}
