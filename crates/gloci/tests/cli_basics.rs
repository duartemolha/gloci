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
