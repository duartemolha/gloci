//! The mdBook site (D94): generated reference, includes, links, and the pages that must not drift.
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Once;

fn root() -> PathBuf {
    xtask::diff::workspace_root()
}

fn src() -> PathBuf {
    root().join("docs/book/src")
}

/// Generate the book once per test binary (the tests share the gitignored output).
fn generated() -> PathBuf {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| xtask::book::generate_book(&src()).unwrap());
    src()
}

/// Every leaf command's path, in clap order.
fn leaves() -> Vec<Vec<String>> {
    fn walk(cmd: &clap::Command, prefix: &[String], out: &mut Vec<Vec<String>>) {
        for sub in cmd.get_subcommands() {
            if sub.get_name() == "help" || sub.is_hide_set() {
                continue;
            }
            let mut path = prefix.to_vec();
            path.push(sub.get_name().to_string());
            if sub.get_subcommands().any(|s| s.get_name() != "help") {
                walk(sub, &path, out);
            } else {
                out.push(path);
            }
        }
    }
    let mut out = Vec::new();
    walk(&gloci::command(), &[], &mut out);
    out
}

fn md_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for e in fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            md_files(&p, out);
        } else if p.extension().is_some_and(|x| x == "md") {
            out.push(p);
        }
    }
}

#[test]
fn book_generates_every_command() {
    let src = generated();
    let leaves = leaves();
    assert!(leaves.len() > 40, "only {} leaf commands?", leaves.len());
    let summary = fs::read_to_string(src.join("SUMMARY.md")).unwrap();
    assert!(
        !summary.contains("{{reference}}"),
        "placeholder left in SUMMARY.md"
    );
    let mut last = 0;
    for path in &leaves {
        let file = format!("reference/{}.md", path.join("-"));
        let page = fs::read_to_string(src.join(&file))
            .unwrap_or_else(|e| panic!("{file} was not generated: {e}"));
        assert!(
            page.starts_with(&format!("# gloci {}\n", path.join(" "))),
            "{file} title"
        );
        assert!(page.contains("\n## Usage\n"), "{file} has no Usage");
        assert!(
            page.contains(&format!("gloci {}", path.join(" "))),
            "{file}: usage does not name the command"
        );
        let at = summary
            .find(&format!("]({file})"))
            .unwrap_or_else(|| panic!("SUMMARY.md does not list {file}"));
        assert!(at > last, "SUMMARY.md lists {file} out of clap order");
        last = at;
    }
    // A page with examples renders them.
    let sort = fs::read_to_string(src.join("reference/sort.md")).unwrap();
    assert!(
        sort.contains("\n## Examples\n"),
        "sort has examples in its help"
    );
    assert!(sort.contains("```sh"));
}

/// The `## Examples` block of a reference page holds the examples only: every line is a gloci
/// command or a comment, never a heading or a sentence of the long help's sections.
#[test]
fn reference_examples_blocks_hold_only_commands() {
    let src = generated();
    for path in leaves() {
        let file = format!("reference/{}.md", path.join("-"));
        let page = fs::read_to_string(src.join(&file)).unwrap();
        let Some(block) = page.split("\n## Examples\n\n```sh\n").nth(1) else {
            continue;
        };
        let block = block.split("\n```").next().unwrap();
        assert!(!block.trim().is_empty(), "{file}: empty Examples block");
        for line in block.lines() {
            assert!(
                line.starts_with("gloci ") || line.starts_with("# "),
                "{file}: not an example line: {line}"
            );
        }
    }
}

#[test]
fn every_page_is_in_the_summary() {
    let src = generated();
    let summary = fs::read_to_string(src.join("SUMMARY.md")).unwrap();
    let mut files = Vec::new();
    md_files(&src, &mut files);
    for f in files {
        let rel = f
            .strip_prefix(&src)
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/");
        if rel != "SUMMARY.md" {
            assert!(
                summary.contains(&format!("]({rel})")),
                "{rel} is not in SUMMARY.md"
            );
        }
    }
}

#[test]
fn book_includes_resolve() {
    let src = generated();
    let mut files = Vec::new();
    md_files(&src, &mut files);
    let mut seen = 0;
    for f in &files {
        let text = fs::read_to_string(f).unwrap();
        for target in xtask::book::includes(&text) {
            let p = f.parent().unwrap().join(&target);
            assert!(
                p.is_file(),
                "{}: include {target} does not exist",
                f.display()
            );
            seen += 1;
        }
    }
    assert!(
        seen >= 5,
        "expected the include-only pages, saw {seen} includes"
    );
    for (page, wants) in [
        ("migration.md", "docs/migration.md"),
        ("differences.md", "docs/differences.md"),
        ("changelog.md", "CHANGELOG.md"),
        ("python.md", "crates/gloci-py/README.md"),
        ("python-api.md", "crates/gloci-py/python/gloci/__init__.pyi"),
    ] {
        let text = fs::read_to_string(src.join(page)).unwrap();
        let hit = xtask::book::includes(&text).into_iter().any(|t| {
            src.join(&t).canonicalize().unwrap() == root().join(wants).canonicalize().unwrap()
        });
        assert!(hit, "{page} must include {wants}");
    }
}

#[test]
fn book_links_resolve() {
    let problems = xtask::book::check_links(&generated()).unwrap();
    assert!(
        problems.is_empty(),
        "broken links:\n{}",
        problems.join("\n")
    );
}

#[test]
fn link_checker_finds_broken_links() {
    let d = tempfile::tempdir().unwrap();
    fs::write(
        d.path().join("a.md"),
        "# Title {#top}\n\n## Some `code` Heading!\n\n\
         [ok](b.md) [ok2](b.md#explicit-id) [ok3](#some-code-heading) [ok4](#top)\n\
         [web](https://example.org/x#y)\n\
         `[inline](nope.md)`\n\n```\n[fenced](nope.md)\n```\n\n\
         [gone](missing.md) [bad anchor](b.md#nowhere) [bad local](#nope)\n",
    )
    .unwrap();
    fs::write(d.path().join("b.md"), "# B\n\n## Thing {#explicit-id}\n").unwrap();
    let problems = xtask::book::check_links(d.path()).unwrap();
    assert_eq!(problems.len(), 3, "{problems:#?}");
    assert!(problems.iter().any(|p| p.contains("missing.md")));
    assert!(problems.iter().any(|p| p.contains("#nowhere")));
    assert!(problems.iter().any(|p| p.contains("#nope")));
}

/// The lines of the first run of table rows that follows `heading` (and, when given,
/// `after_marker` within it), as one block.
fn table_after(text: &str, heading: &str, after_marker: Option<&str>) -> String {
    let mut s = text
        .split(heading)
        .nth(1)
        .unwrap_or_else(|| panic!("no {heading}"));
    if let Some(m) = after_marker {
        s = s.split(m).nth(1).unwrap_or_else(|| panic!("no {m}"));
    }
    let rows: Vec<&str> = s
        .lines()
        .skip_while(|l| !l.starts_with('|'))
        .take_while(|l| l.starts_with('|'))
        .collect();
    assert!(rows.len() > 5, "table under {heading} not found");
    rows.join("\n")
}

#[test]
fn semantics_tables_match_spec() {
    let spec = fs::read_to_string(root().join("docs/design.md")).unwrap();
    let page = fs::read_to_string(src().join("semantics.md")).unwrap();
    let footprints = table_after(&spec, "### 4.2 Variant footprints", None);
    assert!(
        page.contains(&footprints),
        "semantics.md lacks the spec §4.2 table verbatim"
    );
    let overlap = table_after(&spec, "### 4.3 Overlap", Some("Worked examples"));
    assert!(
        page.contains(&overlap),
        "semantics.md lacks the spec §4.3 worked examples verbatim"
    );
    let coords = table_after(&spec, "### 4.1 Coordinates", None);
    assert!(
        page.contains(&coords),
        "semantics.md lacks the spec §4.1 table verbatim"
    );
}

#[test]
fn quickstart_has_the_spec_examples() {
    let spec = fs::read_to_string(root().join("docs/design.md")).unwrap();
    let block = spec
        .split("## 3. Shape of the tool")
        .nth(1)
        .unwrap()
        .split("```sh\n")
        .nth(1)
        .unwrap()
        .split("```")
        .next()
        .unwrap();
    let page = fs::read_to_string(src().join("quickstart.md")).unwrap();
    for line in block.lines().filter(|l| !l.is_empty()) {
        assert!(
            page.contains(line),
            "quickstart.md lacks the spec §3 example: {line}"
        );
    }
}

#[test]
fn install_page_names_every_channel() {
    let page = fs::read_to_string(src().join("install.md")).unwrap();
    for needle in [
        "github.com/duartemolha/gloci/releases",
        "pip install https://github.com/duartemolha/gloci/releases/download/",
        ".whl",
        "cargo install --git https://github.com/duartemolha/gloci",
        "SHA256SUMS",
    ] {
        assert!(page.contains(needle), "install.md lacks {needle:?}");
    }
    // The channels that are not live come after a heading that says so.
    let (live, later) = page
        .split_once("## Once published")
        .expect("install.md needs a `## Once published` section");
    for needle in [
        "cargo install gloci",
        "pip install gloci",
        "conda install -c conda-forge -c bioconda gloci",
        "py-gloci",
        "ghcr.io/",
    ] {
        assert!(
            later.contains(needle),
            "the once-published section lacks {needle:?}"
        );
        assert!(
            !live.contains(needle),
            "{needle:?} is outside the once-published section"
        );
    }
    assert!(
        later.contains("not published"),
        "the section must say nothing is published yet"
    );
}

#[test]
fn releasing_page_lists_the_owner_steps() {
    let page = fs::read_to_string(src().join("maintainers/releasing.md")).unwrap();
    for needle in [
        "trusted publisher",
        "`pypi` environment",
        "CARGO_REGISTRY_TOKEN",
        "GLOCI_PUBLISH_PYPI",
        "GLOCI_PUBLISH_CRATES",
        "GLOCI_PUBLISH_GHCR",
        "dry_run",
        "conda-recipe",
        "bioconda-recipes",
        "Zenodo",
        "Pages",
        "8-core",
        "quiet",
        "0.2.0",
    ] {
        assert!(page.contains(needle), "releasing.md lacks {needle:?}");
    }
}

/// The words of a shell line up to its first unquoted operator (`|`, `>`, `<`, `&&`, `;`, `#`).
fn shell_words(line: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut cur: Option<String> = None;
    let mut quote: Option<char> = None;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match (quote, c) {
            (Some(q), c) if c == q => quote = None,
            (Some(_), c) => cur.get_or_insert_with(String::new).push(c),
            (None, '\'' | '"') => {
                quote = Some(c);
                cur.get_or_insert_with(String::new);
            }
            (None, '|' | '>' | '<' | ';' | '#') => break,
            (None, '&') if chars.peek() == Some(&'&') => break,
            (None, c) if c.is_whitespace() => words.extend(cur.take()),
            (None, c) => cur.get_or_insert_with(String::new).push(c),
        }
    }
    words.extend(cur.take());
    words
}

/// Every `gloci ...` command line of the shell blocks of a Markdown file (a trailing `\` joins lines).
fn gloci_command_lines(md: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut in_sh = false;
    let mut pending = String::new();
    for line in md.lines() {
        let t = line.trim_start();
        if let Some(rest) = t.strip_prefix("```") {
            in_sh = !in_sh && matches!(rest.trim(), "sh" | "bash" | "shell" | "console");
            pending.clear();
            continue;
        }
        if !in_sh {
            continue;
        }
        let t = t.strip_prefix("$ ").unwrap_or(t);
        if pending.is_empty() && !t.starts_with("gloci ") {
            continue;
        }
        if let Some(head) = t.strip_suffix('\\') {
            pending.push_str(head);
            continue;
        }
        pending.push_str(t);
        out.push(std::mem::take(&mut pending));
    }
    out
}

#[test]
fn every_documented_gloci_command_line_parses() {
    let files = [
        src().join("quickstart.md"),
        root().join("crates/gloci-py/README.md"),
        root().join("README.md"),
    ];
    let mut checked = 0;
    let mut failures = Vec::new();
    for f in &files {
        let md = fs::read_to_string(f).unwrap();
        for line in gloci_command_lines(&md) {
            let argv: Vec<std::ffi::OsString> =
                shell_words(&line).into_iter().map(Into::into).collect();
            checked += 1;
            if let Err(e) = gloci::embed::parse(&argv) {
                failures.push(format!("{}: `{line}`: {}", f.display(), e.0));
            }
        }
    }
    assert!(
        checked >= 8,
        "only {checked} command lines found; the extractor is broken"
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
