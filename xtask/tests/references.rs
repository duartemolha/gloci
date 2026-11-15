//! Every reference the tracked tree makes to the design documents resolves, and no comment names a
//! plan task: `D157` is defined in `docs/decisions.md`, `spec §4.3` is a numbered heading of
//! `docs/design.md`, and a source comment does not say which task of an internal plan wrote it.
use std::collections::BTreeSet;
use std::path::Path;
use std::process::Command;

/// Specifications other than gloci's own, whose section numbers a `§` may follow.
const OTHER_SPECS: &[&str] = &[
    "BCF2.2", "BCF2", "BCF", "CRAM 3", "CRAM", "SAM", "SAMv1", "VCF", "VCFv4.3", "VCFv4.4",
    "tabix", "CSI", "BAM", "htslib",
];

/// A reference that does not resolve, or a comment that must go.
struct Finding {
    file: String,
    line: usize,
    what: String,
}

fn root() -> std::path::PathBuf {
    xtask::diff::workspace_root()
}

/// The tracked files whose text can cite a decision or a section: not the internal plans, not
/// `Cargo.lock`, not the decision log itself (it defines the numbers), not (for citations) fixtures or golden
/// outputs (sample names, CIGAR strings and qualities that look like citations), not this file.
fn scanned_files() -> Vec<String> {
    let out = Command::new("git")
        .args(["ls-files"])
        .current_dir(root())
        .output()
        .expect("git ls-files");
    assert!(out.status.success(), "git ls-files failed");
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|f| {
            let ext = f.rsplit('.').next().unwrap_or("");
            let text_kind = matches!(
                ext,
                "md" | "rs" | "py" | "pyi" | "toml" | "sh" | "yml" | "yaml" | "txt"
            ) || f.contains("Dockerfile");
            text_kind
                && !f.starts_with("docs/superpowers/")
                && *f != "Cargo.lock"
                && *f != "xtask/tests/references.rs"
        })
        .map(str::to_string)
        .collect()
}

fn word_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

/// The numbers of the decisions a line cites: a `D` that does not continue a word, then digits
/// that do not continue one (so `10M2D5M` and `0xD9` are not citations).
fn decisions_cited(line: &str) -> Vec<u32> {
    let b = line.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'D' && (i == 0 || !word_byte(b[i - 1])) {
            let mut j = i + 1;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            if j > i + 1
                && (j == b.len() || !word_byte(b[j]))
                && let Ok(n) = line[i + 1..j].parse()
            {
                out.push(n);
            }
            i = j;
        } else {
            i += 1;
        }
    }
    out
}

/// The text of a section number starting at `s` (`4.3`, `6.2.1`, `9`), without a trailing dot.
fn number_at(s: &str) -> Option<&str> {
    let end = s
        .bytes()
        .position(|c| !(c.is_ascii_digit() || c == b'.'))
        .unwrap_or(s.len());
    let n = s[..end].trim_end_matches('.');
    (n.as_bytes().first().is_some_and(u8::is_ascii_digit)).then_some(n)
}

/// Does the text before a `§` name another specification (`BCF2.2 §6.3.3`, also the second of
/// `BCF2.2 §6.2.1, §6.2.2`)?
fn follows_other_spec(before: &str) -> bool {
    let mut rest = before.trim_end();
    loop {
        // Step back over a list of sections: `§6.2.1, `, `§4 and `.
        let trimmed = rest
            .trim_end_matches([',', '/', ' '])
            .trim_end_matches("and")
            .trim_end();
        let digits = trimmed.trim_end_matches(|c: char| c.is_ascii_digit() || c == '.');
        if digits.len() != trimmed.len() && digits.ends_with('§') {
            rest = digits.trim_end_matches('§').trim_end();
            continue;
        }
        return OTHER_SPECS.iter().any(|n| rest.ends_with(n));
    }
}

/// The sections of gloci's own specification a line cites: `spec §4.3`, `design §8.3`, a bare
/// `§4.3` not following another specification's name, `spec 4.3` and `spec section 4.3`.
fn sections_cited(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    for (at, _) in line.match_indices('§') {
        let after = &line[at + '§'.len_utf8()..];
        if let Some(n) = number_at(after)
            && !follows_other_spec(&line[..at])
        {
            out.push(n.to_string());
        }
    }
    let lower = line.to_ascii_lowercase();
    for (at, _) in lower.match_indices("spec") {
        if at > 0 && word_byte(lower.as_bytes()[at - 1]) {
            continue;
        }
        let mut after = &line[at + 4..];
        after = after.strip_prefix("ification").unwrap_or(after);
        let Some(mut after) = after.strip_prefix(' ') else {
            continue;
        };
        if let Some(a) = after.strip_prefix("section ") {
            after = a;
        }
        if let Some(n) = number_at(after)
            && n.contains('.')
        {
            out.push(n.to_string());
        }
    }
    out
}

/// The numbered headings of `docs/design.md` (`## 4.`, `### 4.3 Overlap`), as `4` and `4.3`.
fn design_sections(design: &str) -> BTreeSet<String> {
    design
        .lines()
        .filter(|l| l.starts_with("## ") || l.starts_with("### ") || l.starts_with("#### "))
        .filter_map(|l| number_at(l.trim_start_matches('#').trim_start()).map(str::to_string))
        .collect()
}

/// The line of the decision log without the anchor and the label that define its number.
fn without_definition(line: &str) -> &str {
    let Some(rest) = line.strip_prefix("<a id=\"d") else {
        return line;
    };
    let Some(end) = rest.find("</a>") else {
        return line;
    };
    let rest = &rest[end + 4..];
    rest.strip_prefix("**D")
        .and_then(|r| r.find("**").map(|i| &r[i + 2..]))
        .unwrap_or(rest)
}

/// The part of a line a label could sit in: the text after `//` in Rust, after `#` in scripts,
/// TOML, YAML and Python, the whole line in Markdown. Trailing comments are included.
fn comment_part<'a>(file: &str, line: &'a str) -> &'a str {
    let ext = file.rsplit('.').next().unwrap_or("");
    let marker = match ext {
        "md" => return line,
        "rs" => "//",
        "py" | "pyi" | "toml" | "sh" | "yml" | "yaml" | "txt" => "#",
        _ if file.contains("Dockerfile") => "#",
        _ => return "",
    };
    line.find(marker).map_or("", |i| &line[i..])
}

/// Does a comment (or Markdown text) name a task of an internal plan: `Task` and a number, as
/// in `Task 12`, `(Task 9, R3)` or a phase and task? The words "task" and `asyncio.Task` have
/// no number after them. The design specification is history and keeps its plan mentions.
fn labelled_with_task(file: &str, line: &str) -> bool {
    if file == "docs/design.md" {
        return false;
    }
    let text = comment_part(file, line);
    text.match_indices("Task ").any(|(at, _)| {
        let before_ok = at == 0 || !word_byte(text.as_bytes()[at - 1]);
        before_ok && text[at + 5..].starts_with(|c: char| c.is_ascii_digit())
    })
}

fn findings() -> Vec<Finding> {
    let defined: BTreeSet<u32> = std::fs::read_to_string(root().join("docs/decisions.md"))
        .unwrap()
        .lines()
        .filter_map(|l| {
            let rest = l.strip_prefix("<a id=\"d")?;
            rest[..rest.find('"')?].parse().ok()
        })
        .collect();
    let sections =
        design_sections(&std::fs::read_to_string(root().join("docs/design.md")).unwrap());
    let mut out = Vec::new();
    for f in scanned_files() {
        let Ok(text) = std::fs::read_to_string(Path::new(&root()).join(&f)) else {
            continue;
        };
        for (i, line) in text.lines().enumerate() {
            // The decision log defines its numbers: skip the anchor and label that define one.
            let line = if f == "docs/decisions.md" {
                without_definition(line)
            } else {
                line
            };
            let fixture = f.contains("/inputs/") || f.contains("/golden/");
            let mut push = |what: String| {
                out.push(Finding {
                    file: f.clone(),
                    line: i + 1,
                    what,
                })
            };
            for n in decisions_cited(line).into_iter().filter(|_| !fixture) {
                if !defined.contains(&n) {
                    push(format!(
                        "cites D{n}, which docs/decisions.md does not define: add `<a id=\"d{n}\"></a>**D{n}** ...` there, or fix the number"
                    ));
                }
            }
            for s in sections_cited(line).into_iter().filter(|_| !fixture) {
                if !sections.contains(&s) {
                    push(format!(
                        "cites section {s}, which is not a numbered heading of docs/design.md: fix the section number, or name the other specification before it (`BCF2.2 §6.3.3`)"
                    ));
                }
            }
            if f != "docs/decisions.md" && labelled_with_task(&f, line) {
                push(
                    "names a plan task (`Task N`): drop the label and keep any explanation"
                        .to_string(),
                );
            }
        }
    }
    out
}

#[test]
fn every_cited_decision_and_section_resolves_and_no_comment_names_a_plan_task() {
    let found = findings();
    let report: Vec<String> = found
        .iter()
        .map(|f| format!("{}:{}: {}", f.file, f.line, f.what))
        .collect();
    assert!(
        found.is_empty(),
        "{} unresolved reference(s):\n{}",
        found.len(),
        report.join("\n")
    );
}

#[test]
fn the_scanner_tells_citations_from_look_alikes() {
    assert_eq!(
        decisions_cited("(D37) and D8-D11, D157"),
        vec![37, 8, 11, 157]
    );
    assert!(decisions_cited("CIGAR 10M2D5M, 0xD9, ID9 and D").is_empty());
    assert_eq!(decisions_cited("D1000x D12"), vec![12]);

    assert_eq!(sections_cited("spec §4.3 and §7.8"), vec!["4.3", "7.8"]);
    assert_eq!(
        sections_cited("(spec 8.5), spec section 4.5 rule 3"),
        vec!["8.5", "4.5"]
    );
    assert_eq!(sections_cited("see design §8.3."), vec!["8.3"]);
    assert!(sections_cited("BCF2.2 §6.3.3, CRAM 3 §8, SAM §1.4").is_empty());
    assert!(sections_cited("BCF2.2 §6.2.1, §6.2.2 and §6.2.3").is_empty());
    assert_eq!(
        sections_cited("the specification is long"),
        Vec::<String>::new()
    );

    let design = design_sections("## 4. Core\n### 4.3 Overlap\n### Goals\n## 9. Correctness\n");
    assert!(design.contains("4") && design.contains("4.3") && design.contains("9"));
    assert!(!design.contains("4.4"));

    let label = format!("Task {}", 7);
    assert!(labelled_with_task("a.rs", &format!("// P3b {label}")));
    assert!(labelled_with_task(
        "a.rs",
        &format!("let x = 1; // {label}: pairs are read here")
    ));
    assert!(labelled_with_task(
        "a.sh",
        &format!("# {label}: the inputs")
    ));
    assert!(labelled_with_task(
        "a.toml",
        &format!("x = 1 # (P2b {label}, R3)")
    ));
    assert!(labelled_with_task(
        "a.md",
        &format!("**Pinned ({label}):** every input")
    ));
    assert!(!labelled_with_task(
        "a.rs",
        "// the task of the run; asyncio.Task is not one"
    ));
    assert!(!labelled_with_task(
        "a.rs",
        &format!("let s = \"{label}\";")
    ));
    assert!(!labelled_with_task(
        "docs/design.md",
        &format!("P2b {label} only")
    ));
    assert_eq!(
        without_definition("<a id=\"d12\"></a>**D12** see D13"),
        " see D13"
    );
}
