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
