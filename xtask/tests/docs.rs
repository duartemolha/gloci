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

/// What the generated bash completion offers for `words`, the last of which is being completed
/// (empty for a fresh word); `None` when bash is not available.
fn bash_completes(script: &std::path::Path, words: &[&str]) -> Option<Vec<String>> {
    let quoted: Vec<String> = words.iter().map(|w| format!("'{w}'")).collect();
    let cword = words.len() - 1;
    let prev = if cword > 0 { words[cword - 1] } else { "" };
    let body = format!(
        "source \"$1\"; COMP_WORDS=({}); COMP_CWORD={cword}; COMP_LINE='{}'; \
         COMP_POINT=${{#COMP_LINE}}; _gloci gloci '{}' '{prev}'; printf '%s\\n' \"${{COMPREPLY[@]}}\"",
        quoted.join(" "),
        words.join(" "),
        words[cword],
    );
    let out = std::process::Command::new("bash")
        .args(["-c", &body, "_"])
        .arg(script)
        .output()
        .ok()?;
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    Some(
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect(),
    )
}

/// After a command that takes files (`gloci merge `), bash must fall back to listing files (an
/// empty reply under `-o default`), not offer only the options, whose common prefix `-` it would
/// insert. Commands with subcommands still offer them.
#[test]
fn bash_completion_offers_files_for_a_commands_first_argument() {
    let d = tempfile::tempdir().unwrap();
    xtask::docs::generate_completions(d.path()).unwrap();
    let script = d.path().join("gloci.bash");
    let Some(merge) = bash_completes(&script, &["gloci", "merge", ""]) else {
        eprintln!("bash not found; skipping");
        return;
    };
    assert!(
        merge.is_empty(),
        "gloci merge <TAB> offered {merge:?}, not files"
    );
    for words in [
        &["gloci", "merge", "a.bed", ""][..],
        &["gloci", "vcf", "filter", ""],
        &["gloci", "intersect", "-a", "a.bed", ""],
    ] {
        let got = bash_completes(&script, words).unwrap();
        assert!(got.is_empty(), "{words:?} offered {got:?}, not files");
    }
    // Every command with no subcommands of its own, not only merge.
    fn leaves(cmd: &clap::Command, path: Vec<String>, out: &mut Vec<Vec<String>>) {
        let subs: Vec<_> = cmd
            .get_subcommands()
            .filter(|s| s.get_name() != "help")
            .collect();
        if subs.is_empty() {
            out.push(path.clone());
        }
        for s in subs {
            let mut p = path.clone();
            p.push(s.get_name().to_string());
            leaves(s, p, out);
        }
    }
    let mut all = Vec::new();
    leaves(&gloci::command(), vec!["gloci".to_string()], &mut all);
    assert!(all.len() > 50, "only {} leaf commands found", all.len());
    for path in &all {
        let mut words: Vec<&str> = path.iter().map(String::as_str).collect();
        words.push("");
        let got = bash_completes(&script, &words).unwrap();
        assert!(
            got.is_empty(),
            "{} <TAB> offered {got:?}, not files",
            path.join(" ")
        );
    }
    let opts = bash_completes(&script, &["gloci", "merge", "--max"]).unwrap();
    assert!(opts.contains(&"--max-gap".to_string()), "{opts:?}");
    let top = bash_completes(&script, &["gloci", ""]).unwrap();
    assert!(
        top.contains(&"merge".to_string()) && top.contains(&"vcf".to_string()),
        "{top:?}"
    );
    let vcf = bash_completes(&script, &["gloci", "vcf", ""]).unwrap();
    assert!(vcf.contains(&"filter".to_string()), "{vcf:?}");
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
    let (sem, _) = after_sem.split_once("\n## Commands\n").unwrap();
    let sem_lines = sem.trim().lines().count();
    assert!(
        sem_lines < 40,
        "Core semantics section has {sem_lines} lines (must be < 40)"
    );
    for cmd in gloci::command().get_subcommands() {
        if cmd.get_name() != "help" {
            let entry = format!("\n- `gloci {}`", cmd.get_name());
            assert!(generated.contains(&entry), "missing {entry}");
        }
    }
    for req in [
        "## Commands",
        "## Where the details are",
        "## Pitfalls",
        "(experimental)",
        "--max-missing",
        "--vcf-span ref",
        "--chrom-order karyotype",
        "--mode",
    ] {
        assert!(generated.contains(req), "SKILL.md missing {req:?}");
    }
    // An agent loads the whole guide, so it stays an index: the rules of a command are in its
    // --help and the old-tool tables in docs/migration.md (it was 173 KB with both inlined).
    assert!(
        generated.len() < 48 * 1024,
        "SKILL.md is {} bytes; keep it under 48 KiB",
        generated.len()
    );

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
