#[test]
fn man_and_completions_generate() {
    let d = tempfile::tempdir().unwrap();
    let man = d.path().join("man");
    let comp = d.path().join("completions");
    xtask::docs::generate_man(&man).unwrap();
    xtask::docs::generate_completions(&comp).unwrap();
    for f in ["gloci.1", "gloci-index.1"] {
        assert!(man.join(f).exists(), "{f}");
    }
    for f in ["gloci.bash", "_gloci", "gloci.fish"] {
        assert!(comp.join(f).exists(), "{f}");
    }
}
