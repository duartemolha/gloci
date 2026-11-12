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

#[test]
fn skill_md_is_current() {
    let d = tempfile::tempdir().unwrap();
    let out = d.path().join("SKILL.md");
    xtask::docs::generate_skill(&out).unwrap();
    let generated = std::fs::read_to_string(&out).unwrap();
    let root = xtask::diff::workspace_root();
    let committed_path = root.join("SKILL.md");
    let committed = std::fs::read_to_string(&committed_path).unwrap_or_default();
    assert_eq!(
        generated, committed,
        "SKILL.md is out of date; run: cargo xtask docs --skill SKILL.md"
    );

    // Verify required structure of SKILL.md.
    assert!(generated.starts_with("---\nname: gloci\ndescription: "));
    let (_, after_sem) = generated.split_once("## Core semantics\n").unwrap();
    let (sem, _) = after_sem.split_once("\n## Command reference\n").unwrap();
    let sem_lines = sem.trim().lines().count();
    assert!(
        sem_lines < 40,
        "Core semantics section has {sem_lines} lines (must be < 40)"
    );
    for cmd in gloci::command().get_subcommands() {
        if cmd.get_name() != "help" {
            let hdr = format!("### `gloci {}`", cmd.get_name());
            assert!(generated.contains(&hdr), "missing {hdr}");
        }
    }
    for req in [
        "## Migration tables",
        "## Pitfalls",
        "--max-missing",
        "--vcf-span ref",
        "--chrom-order lex",
        "--mode",
    ] {
        assert!(generated.contains(req), "SKILL.md missing {req:?}");
    }

    // Verify every {#anchor} in docs/differences.md is linked in docs/migration.md.
    let diffs = std::fs::read_to_string(root.join("docs/differences.md")).unwrap();
    let migration = std::fs::read_to_string(root.join("docs/migration.md")).unwrap();
    for line in diffs.lines() {
        if let Some((_, tail)) = line.split_once("{#")
            && let Some((anchor, _)) = tail.split_once('}')
        {
            let link = format!("differences.md#{anchor}");
            assert!(
                migration.contains(&link),
                "docs/migration.md is missing link to {link}"
            );
        }
    }
}

#[test]
fn man_and_completions_tasks_take_out() {
    let d = tempfile::tempdir().unwrap();
    for (task, file) in [("man", "gloci.1"), ("completions", "gloci.bash")] {
        let out = d.path().join(task);
        let st = std::process::Command::new(env!("CARGO_BIN_EXE_xtask"))
            .args([task, "--out"])
            .arg(&out)
            .status()
            .unwrap();
        assert!(st.success(), "{task} --out");
        assert!(out.join(file).exists(), "{task}/{file}");
    }
}
