//! Standard input is read for at most one input of a run. A command given `-` for two of its
//! inputs refuses the run with a usage error (exit 2) that names both, before it reads or writes
//! anything. The commands and their inputs come from the command tree, so a command added later
//! is tested without editing this file.

use crate::common::*;
use crate::helpers::*;
use gloci::embed::{ArgSpec, CommandSpec, Outcome, command_specs};
use std::ffi::OsString;
use std::path::Path;

/// What a command needs besides its inputs for clap to accept it, when the command tree cannot
/// tell: one of several arguments is required, an argument requires another, or a required value
/// has no fixed list. They are added only when clap refuses the argv without them. A command
/// clap refuses for any other reason than two inputs that conflict fails the test, which says
/// what to add here.
const NEEDS: &[(&str, &[&str])] = &[
    ("convert", &["--bam"]),
    ("flank", &["--both", "1"]),
    ("map", &["-c", "5:sum"]),
    ("random", &["-l", "10", "-n", "1"]),
    ("sample", &["--count", "1"]),
    ("shift", &["--by", "1"]),
    ("slop", &["--both", "1"]),
    ("tag", &["--names"]),
    ("vcf compare", &["--out-prefix", "o"]),
    ("vcf extract", &["--info", "DP"]),
    ("vcf popgen", &["--out-prefix", "o"]),
    ("vcf query", &["-f", "%POS"]),
    ("vcf stats", &["--out-prefix", "o", "--freq"]),
    ("windows", &["-w", "10"]),
];

/// The name the message gives an argument: its short flag, else its long one, else (a
/// positional) its value name.
fn name(a: &ArgSpec) -> String {
    match (a.short, &a.long) {
        (Some(c), _) => format!("-{c}"),
        (None, Some(l)) => format!("--{l}"),
        (None, None) => a
            .value_names
            .first()
            .cloned()
            .unwrap_or_else(|| a.id.to_uppercase()),
    }
}

/// The arguments whose values are files to read: the inputs, which `-` makes standard input.
fn inputs(c: &CommandSpec) -> Vec<&ArgSpec> {
    c.args.iter().filter(|a| a.takes_file).collect()
}

/// The pairs of inputs a command is given `-` for, each as the candidates to try in order (the
/// first clap accepts is run): each input after the first with the first one, or with the next
/// input when the two conflict; and each input that takes several files with itself. So every
/// input is given `-` at least once.
fn pairs(c: &CommandSpec) -> Vec<Vec<(&ArgSpec, &ArgSpec)>> {
    let ins = inputs(c);
    let mut out: Vec<Vec<_>> = Vec::new();
    for (i, b) in ins.iter().enumerate().skip(1) {
        let others = ins.iter().enumerate().filter(|(j, _)| *j != i);
        out.push(others.map(|(_, a)| (*a, *b)).collect());
    }
    out.extend(ins.iter().filter(|a| a.multiple).map(|a| vec![(*a, *a)]));
    out
}

/// The argv that gives `-` to both inputs of `pair`, a missing file to every other required
/// input and a value to every other required argument, with the positionals first; and the two
/// names the message gives, in the order of the argv.
fn argv(c: &CommandSpec, pair: (&ArgSpec, &ArgSpec)) -> (Vec<String>, [String; 2]) {
    let twice = std::ptr::eq(pair.0, pair.1);
    let given = |a: &ArgSpec| std::ptr::eq(a, pair.0) || std::ptr::eq(a, pair.1);
    let mut argv: Vec<String> = c.path.clone();
    let mut names = Vec::new();
    let mut dash = |argv: &mut Vec<String>, a: &ArgSpec, k: usize| {
        argv.push("-".into());
        names.push(if twice {
            format!("{} (value {k})", name(a))
        } else {
            name(a)
        });
    };
    // A required value with no fixed list comes from NEEDS.
    let filler = |a: &ArgSpec| -> Option<String> {
        if a.takes_file {
            Some("missing.file".into())
        } else {
            a.possible_values.first().cloned()
        }
    };
    for a in c.args.iter().filter(|a| a.positional) {
        if given(a) {
            dash(&mut argv, a, 1);
            if twice {
                dash(&mut argv, a, 2);
            }
        } else if a.required {
            argv.extend(filler(a));
        }
    }
    for a in c.args.iter().filter(|a| !a.positional) {
        let flag = match (a.short, &a.long) {
            (Some(c), _) => format!("-{c}"),
            (_, l) => format!("--{}", l.as_ref().expect("an option has a flag")),
        };
        if given(a) {
            argv.push(flag.clone());
            dash(&mut argv, a, 1);
            if twice {
                if a.max_values == Some(1) {
                    argv.push(flag);
                }
                dash(&mut argv, a, 2);
            }
        } else if a.required
            && a.takes_value
            && let Some(v) = filler(a)
        {
            argv.push(flag);
            argv.push(v);
        }
    }
    let names: [String; 2] = names.try_into().expect("two inputs given -");
    (argv, names)
}

fn message(names: &[String; 2]) -> String {
    format!(
        "standard input is given for both {} and {}; give one of them as a file",
        names[0], names[1]
    )
}

/// `argv` as clap accepts it: as it is, else with what the command NEEDS. `Ok(None)` when the
/// two inputs conflict, so that they cannot both be given.
fn accepted(c: &CommandSpec, mut argv: Vec<String>) -> Result<Option<Vec<String>>, String> {
    if let Ok(true) = parses(&argv) {
        return Ok(Some(argv));
    }
    if let Some((_, more)) = NEEDS.iter().find(|(n, _)| *n == c.name()) {
        argv.extend(more.iter().map(|s| s.to_string()));
    }
    match parses(&argv) {
        Ok(true) => Ok(Some(argv)),
        Ok(false) => Ok(None),
        Err(e) => Err(format!(
            "clap refuses `gloci {}`; add what {} needs to NEEDS:\n{e}",
            argv.join(" "),
            c.name()
        )),
    }
}

/// Whether clap accepts `argv` (after `gloci`): `Ok(false)` when two arguments conflict.
fn parses(argv: &[String]) -> Result<bool, String> {
    let full: Vec<OsString> = std::iter::once("gloci".into())
        .chain(argv.iter().map(OsString::from))
        .collect();
    match gloci::command().try_get_matches_from(full) {
        Ok(_) => Ok(true),
        Err(e) if e.kind() == clap::error::ErrorKind::ArgumentConflict => Ok(false),
        Err(e) => Err(e.render().to_string()),
    }
}

/// Every command that takes two inputs, given `-` for both, exits 2 naming both, writes nothing
/// to stdout and creates no file.
#[test]
fn every_command_refuses_standard_input_for_two_inputs() {
    let mut ran = Vec::new();
    let mut failures = Vec::new();
    for c in command_specs() {
        for candidates in pairs(&c) {
            let mut chosen = None;
            for pair in candidates {
                let (args, names) = argv(&c, pair);
                match accepted(&c, args) {
                    Ok(Some(args)) => chosen = Some((args, names)),
                    Ok(None) => continue,
                    Err(e) => failures.push(e),
                }
                break;
            }
            let Some((args, names)) = chosen else {
                continue;
            };
            let d = tempdir();
            let out = gloci()
                .args(&args)
                .current_dir(d.path())
                .write_stdin("chr1\t1\t2\n")
                .output()
                .unwrap();
            let want = format!("error: {}\n", message(&names));
            let stderr = String::from_utf8_lossy(&out.stderr);
            let left = leftovers(d.path(), &[]);
            if out.status.code() != Some(2) || !out.stdout.is_empty() || stderr != want {
                failures.push(format!(
                    "gloci {}: exit {:?}, {} bytes on stdout, stderr {stderr:?}",
                    args.join(" "),
                    out.status.code(),
                    out.stdout.len()
                ));
            } else if !left.is_empty() {
                failures.push(format!("gloci {} created {left:?}", args.join(" ")));
            }
            ran.push(args.join(" "));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    // The commands every user reaches first are among them.
    for want in [
        "intersect -a - -b -",
        "subtract -a - -b -",
        "closest -a - -b -",
        "coverage -a - -b -",
        "vcf compare - - --out-prefix o",
        "tag - -b - --names",
        "merge - -",
        "windows -i - --chrom-map - -w 10",
    ] {
        assert!(
            ran.iter().any(|r| r == want),
            "{want} was not run: {ran:#?}"
        );
    }
}

/// An argument whose help speaks of standard input is an input of the command tree, so the
/// check sees `-` given to it.
#[test]
fn every_argument_that_reads_standard_input_is_an_input() {
    for c in command_specs() {
        for a in &c.args {
            let help = a.help.to_lowercase();
            if help.contains("stdin") || help.contains("standard input") {
                assert!(a.takes_file, "{} {}: {}", c.name(), name(a), a.help);
            }
        }
    }
}

/// The check holds for a run, not for the process: in-process runs (the embedding API, the
/// Python library) never see what an earlier run took.
#[test]
fn one_run_never_holds_standard_input_for_the_next() {
    let d = tempdir();
    let out = d.path().join("x.bed");
    let run = |args: &[&str]| {
        let argv: Vec<OsString> = std::iter::once("gloci")
            .chain(args.iter().copied())
            .map(OsString::from)
            .collect();
        let cli = gloci::embed::parse(&argv).unwrap();
        let err = gloci::embed::execute(&cli, &gloci::embed::ExecOptions::default())
            .expect_err("the run fails before it reads standard input");
        gloci::embed::classify(&err)
    };
    let twice = Outcome::Usage(message(&[
        "FILES (value 1)".into(),
        "FILES (value 2)".into(),
    ]));
    assert_eq!(run(&["sort", "-", "-"]), twice);
    // One `-` is accepted (the run then fails on its options, before it opens its input).
    let once = run(&["sort", "--check", "-o", &s_of(&out), "-"]);
    assert_eq!(
        once,
        Outcome::Usage("--check writes nothing, so -o and --write-index do not apply".into())
    );
    assert_eq!(run(&["sort", "-", "-"]), twice);
    assert!(!out.exists());
}

/// Every reader of standard input opens it through `gloci_io::open_stdin`, which holds the
/// check: no other source file calls `stdin()`.
#[test]
fn standard_input_is_opened_in_one_place() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let mut found = Vec::new();
    let mut dirs: Vec<_> = std::fs::read_dir(&crates)
        .unwrap()
        .map(|e| e.unwrap().path().join("src"))
        .filter(|p| p.is_dir())
        .collect();
    while let Some(dir) = dirs.pop() {
        for e in std::fs::read_dir(&dir).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                dirs.push(p);
            } else if p.extension().is_some_and(|x| x == "rs") {
                let text = std::fs::read_to_string(&p).unwrap();
                for (i, line) in text.lines().enumerate() {
                    let call = line.match_indices("stdin()").any(|(i, _)| {
                        !line[..i].ends_with(|c: char| c.is_alphanumeric() || c == '_')
                    });
                    if call && !line.trim_start().starts_with("//") {
                        found.push(format!(
                            "{}:{}",
                            p.strip_prefix(&crates).unwrap().display(),
                            i + 1
                        ));
                    }
                }
            }
        }
    }
    let one = found.len() == 1 && found[0].starts_with("gloci-io/src/input.rs:");
    assert!(one, "stdin() is called outside open_stdin: {found:#?}");
}
