//! The paper's real-data benchmark cases (D96): they parse, name inputs that
//! exist, use valid gloci command lines, and the peak-RSS plumbing works.
use std::collections::BTreeMap;

use clap::Parser;
use xtask::bench::{
    Case, Catalog, gloci_words, parse_cases, parse_rss_json, parse_rss_kb, rss_json,
    substitute_extra, substitute_pair,
};
use xtask::diff::workspace_root;
use xtask::paper::bench::{load_paper_cases, parse_paper_cases};

fn catalog() -> Catalog {
    let text = std::fs::read_to_string(workspace_root().join("xtask/bench/datasets.toml")).unwrap();
    Catalog::parse(&text).unwrap()
}

#[test]
fn paper_cases_parse_and_reference_known_inputs() {
    let cat = catalog();
    let cases = load_paper_cases().unwrap();
    let own: Vec<&Case> = cases.iter().filter(|c| !c.inputs.is_empty()).collect();
    assert!(own.len() >= 19, "only {} real-input cases", own.len());
    for c in &own {
        assert!(c.inputs.contains_key("a"), "{}: no input a", c.name);
        for (slot, name) in &c.inputs {
            assert!(
                ["a", "b", "c", "genome"].contains(&slot.as_str()),
                "{}: slot {slot}",
                c.name
            );
            let known = cat.datasets.iter().any(|d| &d.name == name)
                || cat.derived.iter().any(|d| &d.name == name);
            assert!(known, "{}: input `{name}` is not in datasets.toml", c.name);
        }
        assert!(
            c.table.is_some() && c.label.is_some(),
            "{}: label/table",
            c.name
        );
        assert!(c.threads.is_some(), "{}: threads", c.name);
        // Every placeholder a command uses has an input behind it.
        for (tool, cmd) in &c.commands {
            for (ph, slot) in [("{in_b}", "b"), ("{in_c}", "c"), ("{genome}", "genome")] {
                assert!(
                    !cmd.contains(ph) || c.inputs.contains_key(slot),
                    "{}: `{tool}` uses {ph} but the case has no input {slot}",
                    c.name
                );
            }
        }
        // Targets name commands of the case.
        for t in &c.targets {
            for w in t
                .split_whitespace()
                .filter(|w| w.chars().next().is_some_and(|c| c.is_ascii_alphabetic()))
            {
                assert!(
                    c.commands.contains_key(w),
                    "{}: target `{t}`: no `{w}`",
                    c.name
                );
            }
        }
    }
    // The brief's cases are all there.
    for n in [
        "sort-exons",
        "intersect-ccre-exons",
        "intersect-ccre-exons-bgzf",
        "subtract-exons-ccre",
        "merge-exons",
        "closest-peaks-exons",
        "map-ccre-signal",
        "coverage-ccre-exons",
        "genomecov-exons",
        "partition-3",
        "intersect-1kg-exons",
        "vcf-filter-gnomad",
        "vcf-filter-1kg-bcf",
        "depth-ctcf-bam",
        "flagstat-ctcf-cram",
        "bedcov-ctcf-peaks",
        "vcf-compare-1kg",
    ] {
        assert!(cases.iter().any(|c| c.name == n), "missing case {n}");
    }
}

#[test]
fn every_paper_gloci_command_parses() {
    let mut checked = 0;
    for c in load_paper_cases().unwrap() {
        for (_, cmd) in c.commands.iter().filter(|(k, _)| k.starts_with("gloci")) {
            let args: Vec<String> = gloci_words(cmd)
                .iter()
                .flat_map(|w| {
                    let w = substitute_extra(w, "c.bed", "g.genome", "python3");
                    let w = substitute_pair(&w, "gloci", "a.bed", "b.bed", "out.bed")
                        .replace("{ref}", "ref.fa");
                    if w == "{ins}" {
                        vec!["x.bam".into()]
                    } else {
                        vec![w]
                    }
                })
                .collect();
            if let Err(e) = gloci::cli::Cli::try_parse_from(&args) {
                panic!(
                    "case `{}`: `{}` does not parse:\n{e}",
                    c.name,
                    args.join(" ")
                );
            }
            checked += 1;
        }
    }
    assert!(checked >= 20, "only {checked} gloci command lines found");
}

#[test]
fn from_nightly_names_exist() {
    let root = workspace_root();
    let nightly =
        parse_cases(&std::fs::read_to_string(root.join("xtask/bench/cases.toml")).unwrap())
            .unwrap();
    let cases = load_paper_cases().unwrap();
    for n in [
        "vcf-filter-1kg-chr20",
        "vcf-stats-1kg-chr20",
        "partition-plain-10M",
        "genomecov-plain-10M",
    ] {
        assert!(
            cases.iter().any(|c| c.name == n),
            "{n} not in the paper run"
        );
        assert!(nightly.iter().any(|c| c.name == n), "{n} not in cases.toml");
    }
    // An unknown name is an error, not a silent skip.
    let bad = "from_nightly = [\"no-such-case\"]\n";
    let err = parse_paper_cases(bad, &nightly).unwrap_err().to_string();
    assert!(err.contains("no-such-case"), "{err}");
}

#[test]
fn rss_file_parses() {
    assert_eq!(parse_rss_kb("123456\n").unwrap(), 123_456);
    // A failing command: GNU time prints a status line first.
    assert_eq!(
        parse_rss_kb("Command exited with non-zero status 1\n2048\n").unwrap(),
        2048
    );
    assert!(parse_rss_kb("").is_err());
    assert!(
        parse_rss_kb("0\n").is_err(),
        "0 is a broken time, not a reading"
    );
    assert!(parse_rss_kb("time: cannot run x\n").is_err());

    let case = load_paper_cases()
        .unwrap()
        .into_iter()
        .find(|c| c.name == "merge-exons")
        .unwrap();
    let rows = vec![
        ("gloci".to_string(), Some(4096)),
        ("bedops".to_string(), None),
    ];
    let json = serde_json::to_string(&rss_json(&case, &rows)).unwrap();
    assert_eq!(parse_rss_json(&json).unwrap(), rows);
}

#[test]
fn placeholders_substitute() {
    let cmd = "{gloci} partition {in_a} {in_b} {in_c} -g {genome} > {out}; {python} x";
    let s = substitute_pair(
        &substitute_extra(cmd, "/d/c", "/d/g", "py"),
        "gloci",
        "/d/a",
        "/d/b",
        "/o",
    );
    assert_eq!(s, "gloci partition /d/a /d/b /d/c -g /d/g > /o; py x");
    assert!(!s.contains('{'));
    let _: BTreeMap<String, String> = BTreeMap::new();
}
