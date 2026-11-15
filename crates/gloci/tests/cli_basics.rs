use assert_cmd::Command;

#[test]
fn version_prints_name_and_version() {
    Command::cargo_bin("gloci")
        .unwrap()
        .arg("--version")
        .assert()
        .success()
        .stdout("gloci 0.1.0\n");
}

#[test]
fn no_args_is_usage_error() {
    Command::cargo_bin("gloci").unwrap().assert().code(2);
}

mod common;
use common::gloci as gloci_cmd;

#[test]
fn every_command_has_2_to_4_examples() {
    for sub in gloci::command().get_subcommands_mut() {
        let help = sub.render_long_help().to_string();
        let n = help
            .split("Examples:")
            .nth(1)
            .unwrap_or_else(|| panic!("{}", sub.get_name()))
            .lines()
            .skip(1)
            // The examples end at the first blank line; the headed sections after it may
            // mention other commands.
            .take_while(|l| !l.trim().is_empty())
            .filter(|l| l.trim_start().starts_with("gloci "))
            .count();
        assert!((2..=4).contains(&n), "{}: {n} examples", sub.get_name());
    }
}

/// Stability tiers: a command is in `gloci::CORE_COMMANDS` or its help (and the help of each
/// of its subcommands) starts with the experimental marker, never both; the top-level help
/// names the core.
#[test]
fn every_command_is_core_or_marked_experimental() {
    let mut root = gloci::command();
    let top = root.render_help().to_string();
    for name in gloci::CORE_COMMANDS {
        assert!(
            root.find_subcommand(name).is_some(),
            "CORE_COMMANDS names {name}, which is not a command"
        );
        assert!(
            top.contains(name),
            "the top-level help does not name {name}"
        );
    }
    assert!(top.contains("Stable core:"), "{top}");
    for sub in root.get_subcommands_mut() {
        let name = sub.get_name().to_string();
        if name == "help" {
            continue;
        }
        // The help is wrapped at the terminal width, so compare the words.
        let marked = |c: &mut clap::Command| {
            let help = c.render_help().to_string();
            let words: Vec<&str> = help.split_whitespace().collect();
            let marker: Vec<&str> = gloci::cli::EXPERIMENTAL.split_whitespace().collect();
            words.starts_with(&marker)
        };
        assert_eq!(
            marked(sub),
            !gloci::is_core(&name),
            "{name}: core commands carry no marker, every other command does"
        );
        for nested in sub.get_subcommands_mut() {
            if nested.get_name() == "help" {
                continue;
            }
            assert_eq!(
                marked(nested),
                !gloci::is_core(&name),
                "{name} {}",
                nested.get_name()
            );
        }
    }
}

/// Every command path (`["vcf", "filter"]`, ...) below the root, `help` left out.
fn command_paths(cmd: &clap::Command, prefix: &[String], out: &mut Vec<Vec<String>>) {
    for sub in cmd.get_subcommands() {
        if sub.get_name() == "help" {
            continue;
        }
        let mut path = prefix.to_vec();
        path.push(sub.get_name().to_string());
        command_paths(sub, &path, out);
        out.push(path);
    }
}

fn all_paths() -> Vec<Vec<String>> {
    let mut paths = Vec::new();
    command_paths(&gloci::command(), &[], &mut paths);
    paths
}

/// The command at `path`, built for a 100-column terminal.
fn at_path(path: &[String]) -> clap::Command {
    let mut root = gloci::command().term_width(100).max_term_width(100);
    root.build();
    let mut cmd = &mut root;
    for name in path {
        cmd = cmd.find_subcommand_mut(name).expect("a command");
    }
    cmd.clone()
}

/// The top-level list groups the commands: the stable group is `CORE_COMMANDS`, every other
/// command is in exactly one group of the table in `cli/mod.rs`, and the list prints them in
/// the groups' order.
#[test]
fn every_command_is_in_one_group() {
    let root = gloci::command();
    let all: Vec<String> = root
        .get_subcommands()
        .map(|c| c.get_name().to_string())
        .filter(|n| n != "help")
        .collect();
    let groups = gloci::cli::command_groups();
    assert_eq!(groups[0].0, "Stable core");
    assert_eq!(
        groups[0].1,
        gloci::CORE_COMMANDS,
        "the stable group is CORE_COMMANDS"
    );
    for (title, _) in &groups[1..] {
        assert!(title.ends_with("(experimental)"), "{title}");
    }
    for name in &all {
        let n = groups
            .iter()
            .filter(|(_, names)| names.contains(&name.as_str()))
            .count();
        assert_eq!(n, 1, "{name} is in {n} groups, not exactly one");
    }
    for (title, names) in &groups {
        for name in names {
            assert!(
                all.contains(&name.to_string()),
                "{title} names {name}, not a command"
            );
        }
    }
    let listed: Vec<&str> = groups.iter().flat_map(|(_, n)| n.iter().copied()).collect();
    let top = at_path(&[]).render_help().to_string();
    let mut at = 0;
    for (title, _) in &groups {
        let found = top[at..]
            .find(&format!("{title}:"))
            .unwrap_or_else(|| panic!("no heading {title} in\n{top}"));
        at += found;
    }
    let order: Vec<&str> = top
        .lines()
        .filter(|l| l.starts_with("  ") && !l.starts_with("  -"))
        .filter_map(|l| l.split_whitespace().next())
        .collect();
    assert_eq!(order, listed, "{top}");
    assert!(
        top.contains("Experimental:"),
        "the tiers line is missing:\n{top}"
    );
    // Every line of the list fits an 80-column terminal.
    for l in top.lines() {
        assert!(
            l.chars().count() <= 80,
            "{} characters: {l}",
            l.chars().count()
        );
    }
}

/// The summary of every command is one short line that still names the tool it replaces.
#[test]
fn every_summary_is_at_most_64_characters() {
    for path in all_paths() {
        let cmd = at_path(&path);
        let about = cmd.get_about().map(|a| a.to_string()).unwrap_or_default();
        assert!(!about.is_empty(), "{path:?} has no summary");
        assert!(
            about.chars().count() <= 64,
            "{path:?}: {} characters: {about}",
            about.chars().count()
        );
    }
}

/// The tools a command's summary names (the parenthesised part) are also named in the first
/// paragraph of its long help, which SKILL.md and the Python stub use.
#[test]
fn long_summary_names_every_tool_of_the_short_one() {
    const TOOLS: &[&str] = &[
        "bedtools",
        "bedops",
        "samtools",
        "bcftools",
        "vcftools",
        "tabix",
        "convert2bed",
        "vcf-concat",
    ];
    let mut missing = Vec::new();
    for path in all_paths() {
        let cmd = at_path(&path);
        let about = cmd.get_about().map(|a| a.to_string()).unwrap_or_default();
        let Some((_, tools)) = about.split_once('(') else {
            continue;
        };
        let long = cmd
            .get_long_about()
            .or_else(|| cmd.get_about())
            .map(|a| a.to_string())
            .unwrap_or_default();
        let first = long.trim().split("\n\n").next().unwrap_or("");
        for piece in tools.trim_end_matches(')').split(',') {
            let tool = piece.split_whitespace().next().unwrap_or("");
            if TOOLS.contains(&tool) && !first.contains(tool) {
                missing.push(format!("{path:?} names {tool}: {first}"));
            }
        }
    }
    assert!(
        missing.is_empty(),
        "the long summary must name the tools of the short one:\n{}",
        missing.join("\n")
    );
}

/// The summary paragraph before the usage line is at most 12 lines at 100 columns (the
/// experimental marker above it does not count).
#[test]
fn every_summary_paragraph_is_at_most_12_lines() {
    let marker: Vec<&str> = gloci::cli::EXPERIMENTAL.split_whitespace().collect();
    for path in all_paths() {
        let mut cmd = at_path(&path);
        let help = cmd.render_long_help().to_string();
        let before: Vec<&str> = help
            .lines()
            .take_while(|l| !l.starts_with("Usage:"))
            .collect();
        assert!(
            before.len() < help.lines().count(),
            "{path:?}: no usage line"
        );
        let text = before.join("\n");
        let words: Vec<&str> = text.split_whitespace().collect();
        let body = if words.starts_with(&marker) {
            // Drop the marker's own lines (its first paragraph).
            text.split_once("\n\n").map_or("", |(_, rest)| rest)
        } else {
            &text
        };
        let n = body.trim().lines().count();
        assert!(
            n <= 12,
            "{path:?}: the summary paragraph is {n} lines:\n{body}"
        );
    }
}

#[test]
fn bad_threads_is_usage_error() {
    gloci_cmd()
        .args(["index", "x.gz", "-t", "0"])
        .assert()
        .code(2);
    gloci_cmd()
        .env("GLOCI_THREADS", "0")
        .args(["index", "x.gz"])
        .assert()
        .code(2);
}

#[test]
fn parse_sizes() {
    assert_eq!(gloci::cli::parse_size("2G"), Ok(2 << 30));
    assert_eq!(gloci::cli::parse_size("512M"), Ok(512 << 20));
    assert_eq!(gloci::cli::parse_size("4k"), Ok(4 << 10));
    assert!(gloci::cli::parse_size("10").is_err());
    assert!(gloci::cli::parse_size("1.5G").is_err());
}

#[test]
fn format_flags_list_the_alignment_formats() {
    // Every --format, --a-format and --b-format help names the formats the flag takes.
    let mut cmd = gloci::command();
    for sub in cmd.get_subcommands_mut() {
        let name = sub.get_name().to_string();
        for arg in sub.get_arguments() {
            let long = arg.get_long().unwrap_or_default();
            if !matches!(long, "format" | "a-format" | "b-format") {
                continue;
            }
            let help = arg.get_help().map(|h| h.to_string()).unwrap_or_default();
            // pair-intersect's B does not read alignments (the gate refuses them), so its help
            // names the interval formats and says so instead.
            let formats: &[&str] = if name == "pair-intersect" && long == "b-format" {
                assert!(
                    help.contains("does not read alignments"),
                    "{name} --{long}: {help}"
                );
                &["bed", "vcf"]
            } else {
                &["bed", "vcf", "bam", "cram", "sam"]
            };
            for f in formats {
                assert!(help.contains(f), "{name} --{long}: {help}");
            }
        }
    }
}

#[test]
fn index_help_does_not_claim_maskfasta_needs_an_index() {
    let help = gloci::command()
        .find_subcommand_mut("index")
        .unwrap()
        .render_long_help()
        .to_string()
        .replace('\n', " ");
    assert!(help.contains("maskfasta streams"), "{help}");
    assert!(!help.contains("nuc, maskfasta"), "{help}");
}

/// merge's help says its output is BED regions whatever the input, and its VCF example says so.
#[test]
fn merge_help_says_the_output_is_regions_whatever_the_input() {
    let help = gloci::command()
        .find_subcommand_mut("merge")
        .unwrap()
        .render_long_help()
        .to_string();
    let flat = help.split_whitespace().collect::<Vec<_>>().join(" ");
    assert!(!flat.contains("A's records"), "{help}");
    assert!(flat.contains("whatever the input format"), "{help}");
    let examples = help.split("Examples:").nth(1).unwrap();
    let vcf = examples
        .lines()
        .position(|l| l.contains("calls.vcf.gz"))
        .unwrap();
    let above = examples.lines().nth(vcf - 1).unwrap().trim_start();
    assert!(above.starts_with('#') && above.contains("BED"), "{above}");
}
