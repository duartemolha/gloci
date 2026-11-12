//! The mdBook site (D94): the generated command reference, `SUMMARY.md`, and the link checker.
//!
//! `docs/book/src` holds the hand-written pages and the include-only pages; `generate_book`
//! adds `reference/<command-path>.md` for every leaf command and `SUMMARY.md` (from
//! `docs/book/SUMMARY.in.md`). Both outputs are gitignored.
use anyhow::{Context, Result, bail};
use std::collections::{BTreeSet, HashMap};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

const INCLUDE: &str = "{{#include ";

/// The paths named by `{{#include PATH}}` directives, in order.
pub fn includes(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(at) = rest.find(INCLUDE) {
        let after = &rest[at + INCLUDE.len()..];
        let Some(end) = after.find("}}") else { break };
        out.push(after[..end].trim().to_string());
        rest = &after[end + 2..];
    }
    out
}

/// `text` with every include replaced by the included file (recursively), as mdBook does.
fn expand(text: &str, dir: &Path, depth: u32) -> Result<String> {
    if depth > 8 {
        bail!("includes nested too deeply in {}", dir.display());
    }
    let mut out = String::new();
    let mut rest = text;
    while let Some(at) = rest.find(INCLUDE) {
        out.push_str(&rest[..at]);
        let after = &rest[at + INCLUDE.len()..];
        let Some(end) = after.find("}}") else {
            bail!("unterminated include in {}", dir.display());
        };
        let path = dir.join(after[..end].trim());
        let body = std::fs::read_to_string(&path)
            .with_context(|| format!("include {} (from {})", path.display(), dir.display()))?;
        out.push_str(&expand(&body, path.parent().unwrap_or(dir), depth + 1)?);
        rest = &after[end + 2..];
    }
    out.push_str(rest);
    Ok(out)
}

fn is_fence(line: &str) -> bool {
    let t = line.trim_start();
    t.starts_with("```") || t.starts_with("~~~")
}

/// The ids mdBook gives a page's headings: the explicit `{#id}`, or the slug of the text.
fn heading_ids(text: &str) -> BTreeSet<String> {
    let mut ids = BTreeSet::new();
    let mut seen: HashMap<String, u32> = HashMap::new();
    let mut fenced = false;
    for line in text.lines() {
        if is_fence(line) {
            fenced = !fenced;
            continue;
        }
        if fenced {
            continue;
        }
        let hashes = line.bytes().take_while(|&b| b == b'#').count();
        if !(1..=6).contains(&hashes) || !line[hashes..].starts_with(' ') {
            continue;
        }
        let mut title = line[hashes..].trim().trim_end_matches('#').trim_end();
        let mut explicit = None;
        if let Some(open) = title.rfind("{#")
            && title.ends_with('}')
        {
            explicit = Some(title[open + 2..title.len() - 1].to_string());
            title = title[..open].trim_end();
        }
        let id = explicit.unwrap_or_else(|| slug(title));
        let n = seen.entry(id.clone()).or_insert(0);
        let unique = if *n == 0 {
            id.clone()
        } else {
            format!("{id}-{n}")
        };
        *n += 1;
        ids.insert(unique);
    }
    ids
}

/// mdBook's `normalize_id` applied to the text of a heading.
fn slug(title: &str) -> String {
    let mut text = String::new();
    let mut chars = title.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '`' | '*' => {}
            // [text](url) contributes its text only.
            '[' => {}
            ']' if chars.peek() == Some(&'(') => {
                for d in chars.by_ref() {
                    if d == ')' {
                        break;
                    }
                }
            }
            ']' => {}
            c => text.push(c),
        }
    }
    text.chars()
        .filter_map(|c| {
            if c.is_alphanumeric() || c == '_' || c == '-' {
                Some(c.to_ascii_lowercase())
            } else if c.is_whitespace() {
                Some('-')
            } else {
                None
            }
        })
        .collect()
}

/// `(line number, target)` of every inline link outside code.
fn links(text: &str) -> Vec<(usize, String)> {
    let mut out = Vec::new();
    let mut fenced = false;
    for (i, line) in text.lines().enumerate() {
        if is_fence(line) {
            fenced = !fenced;
            continue;
        }
        if fenced {
            continue;
        }
        // Drop inline code spans.
        let mut plain = String::new();
        let mut in_code = false;
        for c in line.chars() {
            if c == '`' {
                in_code = !in_code;
            } else if !in_code {
                plain.push(c);
            }
        }
        let mut rest = plain.as_str();
        while let Some(at) = rest.find("](") {
            let after = &rest[at + 2..];
            let Some(end) = after.find(')') else { break };
            let target = after[..end].trim();
            let target = target.split_once(' ').map_or(target, |(t, _)| t);
            out.push((i + 1, target.trim_matches(['<', '>']).to_string()));
            rest = &after[end + 1..];
        }
    }
    out
}

fn md_files(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    let mut entries: Vec<_> = std::fs::read_dir(dir)
        .with_context(|| format!("reading {}", dir.display()))?
        .collect::<std::io::Result<_>>()?;
    entries.sort_by_key(|e| e.path());
    for e in entries {
        let p = e.path();
        if p.is_dir() {
            md_files(&p, out)?;
        } else if p.extension().is_some_and(|x| x == "md") {
            out.push(p);
        }
    }
    Ok(())
}

/// Check every relative link and `#anchor` of every page under `src` (includes expanded, so
/// links resolve from the including page as they do in the built book). Returns one line per
/// problem; an empty list means everything resolves.
pub fn check_links(src: &Path) -> Result<Vec<String>> {
    let mut files = Vec::new();
    md_files(src, &mut files)?;
    let mut pages: HashMap<PathBuf, String> = HashMap::new();
    for f in &files {
        let raw = std::fs::read_to_string(f).with_context(|| format!("reading {}", f.display()))?;
        let text = expand(&raw, f.parent().unwrap_or(src), 0)?;
        pages.insert(f.clone(), text);
    }
    let ids: HashMap<&PathBuf, BTreeSet<String>> =
        pages.iter().map(|(p, t)| (p, heading_ids(t))).collect();

    let mut problems = Vec::new();
    for f in &files {
        let shown = f.strip_prefix(src).unwrap_or(f).display();
        for (line, target) in links(&pages[f]) {
            if target.is_empty() || target.contains("://") || target.starts_with("mailto:") {
                continue;
            }
            let (path, anchor) = match target.split_once('#') {
                Some((p, a)) => (p, Some(a)),
                None => (target.as_str(), None),
            };
            let dest = if path.is_empty() {
                f.clone()
            } else {
                let joined = f.parent().unwrap_or(src).join(path);
                if !joined.exists() {
                    problems.push(format!(
                        "{shown}:{line}: link `{target}`: {path} does not exist"
                    ));
                    continue;
                }
                joined
            };
            let Some(anchor) = anchor.filter(|a| !a.is_empty()) else {
                continue;
            };
            // Anchors are checked in pages we know; a normalised path finds them.
            let key = dest.canonicalize().unwrap_or(dest.clone());
            let page_ids = ids
                .iter()
                .find(|(p, _)| p.canonicalize().ok().as_ref() == Some(&key))
                .map(|(_, v)| v);
            match page_ids {
                Some(v) if v.contains(anchor) => {}
                Some(_) => problems.push(format!(
                    "{shown}:{line}: link `{target}`: no heading or {{#id}} `#{anchor}` in {}",
                    if path.is_empty() {
                        shown.to_string()
                    } else {
                        path.to_string()
                    }
                )),
                None => {} // an anchor into a non-Markdown file
            }
        }
    }
    Ok(problems)
}

fn leaf_name(path: &[String]) -> String {
    path.join("-")
}

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn reference_page(path: &[String], cmd: &clap::Command) -> String {
    let full = format!("gloci {}", path.join(" "));
    let about = cmd
        .get_long_about()
        .or_else(|| cmd.get_about())
        .map(|s| s.to_string())
        .unwrap_or_default();
    let mut page = format!("# {full}\n\n");
    if !about.trim().is_empty() {
        let _ = writeln!(page, "{}\n", escape_html(about.trim()));
    }
    let help = cmd.clone().render_long_help().to_string();
    let _ = writeln!(page, "## Usage\n\n```text\n{}\n```\n", help.trim_end());
    let examples = crate::docs::examples_of(cmd);
    if !examples.is_empty() {
        let _ = writeln!(page, "## Examples\n\n```sh");
        for e in examples {
            let _ = writeln!(page, "{e}");
        }
        let _ = writeln!(page, "```");
    }
    page
}

/// Write `reference/<path-with-dashes>.md` per leaf command and `SUMMARY.md` into `src`.
/// `SUMMARY.md` is `../SUMMARY.in.md` with `{{reference}}` replaced by the command list, in
/// clap order (a command with subcommands is a draft chapter that nests its leaves).
pub fn generate_book(src: &Path) -> Result<()> {
    let template_path = src.parent().unwrap_or(src).join("SUMMARY.in.md");
    let template = std::fs::read_to_string(&template_path)
        .with_context(|| format!("reading {}", template_path.display()))?;
    if !template.contains("{{reference}}") {
        bail!(
            "{} has no {{{{reference}}}} placeholder",
            template_path.display()
        );
    }
    let ref_dir = src.join("reference");
    if ref_dir.exists() {
        std::fs::remove_dir_all(&ref_dir)?;
    }
    std::fs::create_dir_all(&ref_dir)?;

    let mut root = gloci::command();
    root.build();
    let mut list = String::new();
    fn walk(
        cmd: &clap::Command,
        prefix: &[String],
        depth: usize,
        ref_dir: &Path,
        list: &mut String,
    ) -> Result<()> {
        for sub in cmd.get_subcommands() {
            if sub.get_name() == "help" || sub.is_hide_set() {
                continue;
            }
            let mut path = prefix.to_vec();
            path.push(sub.get_name().to_string());
            let indent = "    ".repeat(depth);
            let title = format!("gloci {}", path.join(" "));
            if sub.get_subcommands().any(|s| s.get_name() != "help") {
                let _ = writeln!(list, "{indent}- [{title}]()");
                walk(sub, &path, depth + 1, ref_dir, list)?;
            } else {
                let file = format!("{}.md", leaf_name(&path));
                std::fs::write(ref_dir.join(&file), reference_page(&path, sub))
                    .with_context(|| format!("writing reference/{file}"))?;
                let _ = writeln!(list, "{indent}- [{title}](reference/{file})");
            }
        }
        Ok(())
    }
    walk(&root, &[], 0, &ref_dir, &mut list)?;
    let summary = template.replace("{{reference}}", list.trim_end());
    std::fs::write(src.join("SUMMARY.md"), summary).context("writing SUMMARY.md")?;
    Ok(())
}
