//! The real-data concordance cases (D97): they parse into valid gloci command
//! lines on inputs `datasets.toml` pins, their anchors are documented and
//! pinned by a golden, every comparable command of docs/migration.md has a
//! case (or a stated reason not to), and the verdicts follow the rules.
use clap::Parser;
use gloci_testkit::differential::MatchStats;
use xtask::bench::Catalog;
use xtask::diff::workspace_root;
use xtask::paper::concord::{
    self, CSV_HEADER, ConcordCase, Expect, NO_CONCORDANCE_CASE, Verdict, csv_row, gloci_command,
    gloci_invocations, load_cases, merge_csv, migration_commands, parse_cases, verdict,
};

fn catalog() -> Catalog {
    let text = std::fs::read_to_string(workspace_root().join("xtask/bench/datasets.toml")).unwrap();
    Catalog::parse(&text).unwrap()
}

/// `cmd` with placeholders filled by stand-in paths.
fn filled(cmd: &str) -> String {
    concord::fill(
        cmd,
        &concord::Slots {
            gloci: "gloci".into(),
            in_a: "a.bed".into(),
            in_b: "b.bed".into(),
            in_c: "c.bed".into(),
            genome: "g.genome".into(),
            out: "out".into(),
        },
    )
}

#[test]
fn concordance_cases_parse() {
    let cat = catalog();
    let cases = load_cases().unwrap();
    assert!(cases.len() >= 50, "only {} concordance cases", cases.len());
    // The rows of the brief's table, by command and reference tool.
    for (cmd, tool) in [
        ("sort", "sort-bed"),
        ("sort", "bedtools"),
        ("intersect", "bedtools"),
        ("intersect", "bedops"),
        ("subtract", "bedtools"),
        ("subtract", "bedops"),
        ("closest", "bedtools"),
        ("closest", "bedops"),
        ("map", "bedtools"),
        ("map", "bedops"),
        ("coverage", "bedtools"),
        ("merge", "bedtools"),
        ("merge", "bedops"),
        ("complement", "bedtools"),
        ("slop", "bedtools"),
        ("flank", "bedtools"),
        ("shift", "bedtools"),
        ("windows", "bedtools"),
        ("genomecov", "bedtools"),
        ("partition", "bedops"),
        ("partition", "bedtools"),
        ("cluster", "bedtools"),
        ("spacing", "bedtools"),
        ("groupby", "bedtools"),
        ("convert", "bedops"),
        ("view", "tabix"),
        ("vcf filter", "vcftools"),
        ("vcf filter", "bcftools"),
        ("vcf stats", "vcftools"),
        ("vcf extract", "vcftools"),
        ("vcf concat", "bcftools"),
    ] {
        assert!(
            cases
                .iter()
                .any(|c| c.command == cmd && c.reference_tool == tool),
            "no concordance case for `gloci {cmd}` against {tool}"
        );
    }
    // The brief's compat runs.
    for (cmd, flags) in [
        ("sort", "--full-key"),
        ("intersect", "--vcf-span ref"),
        ("view", "--vcf-span ref --inclusive-points"),
    ] {
        assert!(
            cases
                .iter()
                .any(|c| c.command == cmd && c.compat.as_deref() == Some(flags)),
            "no `{cmd}` case with compat `{flags}`"
        );
    }
    let vcf_stats: Vec<&ConcordCase> = cases.iter().filter(|c| c.command == "vcf stats").collect();
    assert!(
        vcf_stats.len() >= 10,
        "vcf stats has {} cases",
        vcf_stats.len()
    );

    for c in &cases {
        for name in c.inputs.values() {
            let known = cat.datasets.iter().any(|d| &d.name == name)
                || cat.derived.iter().any(|d| &d.name == name);
            assert!(known, "{}: input `{name}` is not in datasets.toml", c.name);
        }
        // Every gloci invocation of both runs is a valid command line.
        let runs = if c.compat.is_some() {
            vec![false, true]
        } else {
            vec![false]
        };
        for compat in runs {
            let cmd = filled(&gloci_command(c, compat));
            let inv = gloci_invocations(&cmd, "gloci");
            assert!(!inv.is_empty(), "{}: no gloci invocation", c.name);
            for args in inv {
                if let Err(e) = gloci::cli::Cli::try_parse_from(&args) {
                    panic!("{}: `{}` does not parse:\n{e}", c.name, args.join(" "));
                }
            }
        }
    }
}

#[test]
fn default_runs_add_chrom_order_lex_and_compat_flags() {
    let case = |gloci: &str, compat: Option<&str>| ConcordCase {
        name: "x".into(),
        command: "sort".into(),
        reference_tool: "sort-bed".into(),
        inputs: [("a".to_string(), "exons.bed".to_string())].into(),
        gloci: gloci.into(),
        reference: "sort-bed {in_a} > {out}".into(),
        compare: "exact".into(),
        compat: compat.map(String::from),
        ..ConcordCase::default()
    };
    let c = case("{gloci} sort {in_a} > {out}", Some("--full-key"));
    assert_eq!(
        gloci_command(&c, false),
        "{gloci} sort --chrom-order lex {in_a} > {out}"
    );
    assert_eq!(
        gloci_command(&c, true),
        "{gloci} sort --chrom-order lex --full-key {in_a} > {out}"
    );
    // Each invocation of a pipeline that takes the flag gets it; compat goes to the first.
    let c = case(
        "{gloci} merge {in_a} | {gloci} subtract -a - -b {in_b} > {out}",
        Some("-d 5"),
    );
    assert_eq!(
        gloci_command(&c, true),
        "{gloci} merge --chrom-order lex -d 5 {in_a} | {gloci} subtract --chrom-order lex -a - -b {in_b} > {out}"
    );
    // A subcommand without --chrom-order is left alone, and so is one that sets it.
    let c = case("{gloci} vcf filter --pass-only {in_a} > {out}", None);
    assert_eq!(gloci_command(&c, false), c.gloci);
    let c = case("{gloci} sort --chrom-order natural {in_a} > {out}", None);
    assert_eq!(gloci_command(&c, false), c.gloci);
}

#[test]
fn concordance_case_file_errors() {
    let ok = r#"
[[case]]
name = "a"
command = "sort"
reference_tool = "sort-bed"
inputs = { a = "exons.bed" }
gloci = "{gloci} sort {in_a} > {out}"
reference = "sort-bed {in_a} > {out}"
"#;
    assert_eq!(parse_cases(ok).unwrap()[0].compare, "exact");
    let bad = |from: &str, to: &str| parse_cases(&ok.replace(from, to)).unwrap_err().to_string();
    assert!(bad("sort-bed\"\n", "sort-bed\"\ncompare = \"bam\"\n").contains("compare"));
    assert!(bad("{ a = ", "{ z = ").contains("slot"));
    assert!(bad("sort {in_a} > {out}\"", "sort {in_a}\"").contains("{out}"));
    assert!(bad("sort {in_a} > {out}\"", "sort {in_b} > {out}\"").contains("{in_b}"));
    assert!(bad("command = \"sort\"", "command = \"nosuch\"").contains("nosuch"));
    assert!(bad("sort-bed\"\n", "sort-bed\"\ncompat_difference = \"#x\"\n").contains("compat"));
    assert!(
        bad("sort {in_a} > {out}\"", "sort <(cat {in_a}) > {out}\"")
            .contains("process substitution")
    );
    let with =
        |extra: &str| parse_cases(&ok.replace("sort-bed\"\n", &format!("sort-bed\"\n{extra}\n")));
    assert!(
        with("min_rate = 0.5")
            .unwrap_err()
            .to_string()
            .contains("need")
    );
    assert!(with("expected_difference = \"#sort-tie-order\"\nmin_rate = 0.5").is_ok());
    assert!(
        with("expected_difference = \"#sort-tie-order\"\nmin_rate = 1.0")
            .unwrap_err()
            .to_string()
            .contains("below 1")
    );
    assert!(
        with("expected_difference = \"#sort-tie-order\"\norder_only = true\nmin_rate = 0.5")
            .unwrap_err()
            .to_string()
            .contains("order_only")
    );
    assert!(with("compat_order_only = true").is_err());
    assert!(with("allow_empty = \"no peaks on chrY\"").is_ok());
    assert!(
        with("allow_empty = \" \"")
            .unwrap_err()
            .to_string()
            .contains("allow_empty")
    );
    let twice = format!("{ok}{ok}");
    assert!(
        parse_cases(&twice)
            .unwrap_err()
            .to_string()
            .contains("duplicate")
    );
}

#[test]
fn concordance_anchors_exist_and_are_pinned() {
    let root = workspace_root();
    let doc = std::fs::read_to_string(root.join("docs/differences.md")).unwrap();
    // Every expected_difference a golden carries.
    let mut pinned = std::collections::BTreeSet::new();
    for e in std::fs::read_dir(root.join("crates/gloci/tests/golden")).unwrap() {
        let p = e.unwrap().path().join("expected_difference");
        if let Ok(t) = std::fs::read_to_string(&p) {
            pinned.insert(t.trim().trim_start_matches('#').to_string());
        }
    }
    assert!(pinned.len() >= 20, "only {} pinned anchors", pinned.len());
    let mut n = 0;
    for c in load_cases().unwrap() {
        for a in [&c.expected_difference, &c.compat_difference]
            .into_iter()
            .flatten()
        {
            n += 1;
            gloci_testkit::differential::check_anchor(&doc, a)
                .unwrap_or_else(|e| panic!("{}: {e}", c.name));
            assert!(
                pinned.contains(a.trim_start_matches('#')),
                "{}: anchor {a} is pinned by no golden expected_difference",
                c.name
            );
        }
    }
    assert!(n >= 10, "only {n} anchors in the concordance cases");
}

#[test]
fn every_comparable_command_has_a_concordance_case() {
    let md = std::fs::read_to_string(workspace_root().join("docs/migration.md")).unwrap();
    let commands = migration_commands(&md);
    for want in [
        "sort",
        "intersect",
        "vcf filter",
        "vcf stats",
        "bedcov",
        "stats jaccard",
    ] {
        assert!(
            commands.iter().any(|c| c == want),
            "migration_commands missed `{want}`: {commands:?}"
        );
    }
    let cases = load_cases().unwrap();
    for c in &commands {
        let covered = cases.iter().any(|k| &k.command == c);
        let exempt = NO_CONCORDANCE_CASE.iter().any(|(n, _)| n == c);
        assert!(
            covered || exempt,
            "`gloci {c}` is available (docs/migration.md) but has no concordance case \
             and no entry in NO_CONCORDANCE_CASE"
        );
        assert!(
            !(covered && exempt),
            "`gloci {c}` has a concordance case: drop it from NO_CONCORDANCE_CASE"
        );
    }
    for (n, why) in NO_CONCORDANCE_CASE {
        assert!(
            commands.iter().any(|c| c == n),
            "NO_CONCORDANCE_CASE names `{n}`, which docs/migration.md does not map"
        );
        assert!(
            why.len() >= 12 && !why.contains('\n'),
            "`{n}`: the reason must be one line"
        );
    }
    // A case names a command docs/migration.md maps.
    for k in &cases {
        assert!(
            commands.iter().any(|c| c == &k.command),
            "{}: `{}` is not a command of docs/migration.md",
            k.name,
            k.command
        );
    }
}

#[test]
fn verdict_rules() {
    let s = |n_gloci, n_reference, n_matching, same_order| MatchStats {
        n_gloci,
        n_reference,
        n_matching,
        same_order,
    };
    let none = Expect::default();
    let rate = |min_rate| Expect {
        anchor: Some("#closest-ties-k"),
        min_rate,
        ..Expect::default()
    };
    let order = Expect {
        anchor: Some("#sort-tie-order"),
        order_only: true,
        ..Expect::default()
    };
    // Nothing documented.
    assert_eq!(verdict(Some(&s(10, 10, 10, true)), &none), Verdict::Match);
    // An empty reference output is a failed run, unless the case allows it.
    assert_eq!(verdict(Some(&s(0, 0, 0, true)), &none), Verdict::Failed);
    assert_eq!(verdict(Some(&s(5, 0, 0, false)), &none), Verdict::Failed);
    let empty_ok = Expect {
        allow_empty: true,
        ..Expect::default()
    };
    assert_eq!(verdict(Some(&s(0, 0, 0, true)), &empty_ok), Verdict::Match);
    assert_eq!(
        verdict(Some(&s(5, 0, 0, false)), &empty_ok),
        Verdict::Mismatch
    );
    assert_eq!(
        verdict(Some(&s(10, 10, 10, false)), &none),
        Verdict::Reordered
    );
    assert_eq!(
        verdict(Some(&s(10, 10, 9, false)), &none),
        Verdict::Mismatch
    );
    assert_eq!(
        verdict(Some(&s(10, 11, 10, false)), &none),
        Verdict::Mismatch
    );
    assert_eq!(verdict(None, &none), Verdict::Failed);
    // A documented difference that does not show.
    assert_eq!(
        verdict(Some(&s(10, 10, 10, true)), &rate(Some(0.9))),
        Verdict::AnchorUnused
    );
    assert_eq!(
        verdict(Some(&s(10, 10, 10, true)), &order),
        Verdict::AnchorUnused
    );
    // A rate difference passes down to its min_rate, not below, not without one.
    assert_eq!(
        verdict(Some(&s(10, 10, 9, false)), &rate(Some(0.9))),
        Verdict::Documented
    );
    assert_eq!(
        verdict(Some(&s(10, 10, 8, false)), &rate(Some(0.9))),
        Verdict::Mismatch
    );
    assert_eq!(
        verdict(Some(&s(10, 10, 9, false)), &rate(None)),
        Verdict::Mismatch
    );
    assert_eq!(
        verdict(Some(&s(1000, 1000, 0, false)), &rate(Some(0.0))),
        Verdict::Documented
    );
    // ...and a run that only reorders is not what it documents.
    assert_eq!(
        verdict(Some(&s(10, 10, 10, false)), &rate(Some(0.9))),
        Verdict::Reordered
    );
    // An order-only difference passes with every line matched, nothing else.
    assert_eq!(
        verdict(Some(&s(10, 10, 10, false)), &order),
        Verdict::Documented
    );
    assert_eq!(
        verdict(Some(&s(10, 10, 9, false)), &order),
        Verdict::Mismatch
    );
    assert_eq!(
        verdict(Some(&s(10, 11, 10, false)), &order),
        Verdict::Mismatch
    );
    assert_eq!(verdict(None, &order), Verdict::Failed);
    for v in [Verdict::Match, Verdict::AnchorUnused, Verdict::Documented] {
        assert!(v.passes(), "{v:?}");
    }
    for v in [Verdict::Reordered, Verdict::Mismatch, Verdict::Failed] {
        assert!(!v.passes(), "{v:?}");
    }
    assert_eq!(Verdict::AnchorUnused.as_str(), "anchor-unused");
    assert_eq!(
        CSV_HEADER,
        "case,command,reference_tool,dataset,mode,n_gloci,n_reference,n_matching,match_rate,verdict,anchors"
    );
}

#[test]
fn a_failing_step_fails_the_run() {
    let d = tempfile::tempdir().unwrap();
    let ok = |cmd: &str| concord::host_command_succeeds(cmd, d.path()).unwrap();
    // The form reshaped references use: an intermediate file in the same
    // command, so a failing step stops it.
    assert!(ok("printf 'a\\n' > out.x && cat out.x > out && rm out.x"));
    assert!(!ok(
        "awk 'BEGIN { exit 3 }' > out.x && cat out.x > out && rm out.x"
    ));
    // pipefail: a failing first stage fails the pipe.
    assert!(!ok("false | cat > out"));
    // A process substitution would hide it, so the cases may not use one.
    assert!(ok("cat <(false) > out"));
    // A failed command is a failed run, whatever is documented.
    assert_eq!(verdict(None, &Expect::default()), Verdict::Failed);
}

#[test]
fn concordance_csv_row_format() {
    let c = ConcordCase {
        name: "sort-x".into(),
        command: "sort".into(),
        reference_tool: "sort-bed".into(),
        inputs: [
            ("a".to_string(), "exons.bed".to_string()),
            ("genome".to_string(), "hg38.lex.genome".to_string()),
        ]
        .into(),
        gloci: "{gloci} sort {in_a} > {out}".into(),
        reference: "sort-bed {in_a} > {out}".into(),
        compare: "exact".into(),
        ..ConcordCase::default()
    };
    let s = MatchStats {
        n_gloci: 4,
        n_reference: 4,
        n_matching: 3,
        same_order: false,
    };
    // Anchors are plain names, match_rate a decimal from 0 to 1.
    assert_eq!(
        csv_row(
            &c,
            "default",
            Some(&s),
            Verdict::Documented,
            Some("#sort-tie-order")
        ),
        "sort-x,sort,sort-bed,exons.bed+hg38.lex.genome,default,4,4,3,0.750000,documented,sort-tie-order"
    );
    // Truncated, not rounded: a run that misses one line in ten million is
    // not 1.000000.
    let s = MatchStats {
        n_gloci: 10_000_000,
        n_reference: 10_000_000,
        n_matching: 9_999_999,
        same_order: false,
    };
    assert!(
        csv_row(&c, "default", Some(&s), Verdict::Mismatch, None).contains(",0.999999,mismatch,")
    );
    let s = MatchStats {
        n_gloci: 10,
        n_reference: 10,
        n_matching: 3,
        same_order: false,
    };
    assert!(csv_row(&c, "default", Some(&s), Verdict::Mismatch, None).contains(",0.300000,"));
    assert_eq!(
        csv_row(&c, "compat", None, Verdict::Failed, None),
        "sort-x,sort,sort-bed,exons.bed+hg38.lex.genome,compat,,,,,failed,"
    );
}

#[test]
fn concordance_csv_keeps_rows_of_cases_not_run() {
    let order = ["a".to_string(), "b".to_string(), "c".to_string()];
    let old = format!(
        "{CSV_HEADER}\nb,merge,bedtools,exons.bed,default,1,1,1,1.000000,match,\nc,x,y,z,default,1,1,0,0.000000,mismatch,\ngone,x,y,z,default,1,1,1,1.000000,match,\n"
    );
    let new = vec![
        "c,x,y,z,compat,1,1,1,1.000000,match,".to_string(),
        "c,x,y,z,default,1,1,1,1.000000,match,".to_string(),
        "a,sort,sort-bed,exons.bed,default,1,1,1,1.000000,match,".to_string(),
    ];
    let merged = merge_csv(
        Some(&old),
        &new,
        &["c".to_string(), "a".to_string()],
        &order,
    );
    assert_eq!(
        merged,
        format!(
            "{CSV_HEADER}\na,sort,sort-bed,exons.bed,default,1,1,1,1.000000,match,\nb,merge,bedtools,exons.bed,default,1,1,1,1.000000,match,\nc,x,y,z,default,1,1,1,1.000000,match,\nc,x,y,z,compat,1,1,1,1.000000,match,\n"
        )
    );
}

#[test]
fn reference_cache_key_names_the_image_and_input_digests() {
    let cat = catalog();
    // Pinned digests: a derived input's sha256, a dataset's own checksum.
    assert_eq!(
        concord::pinned_digest(&cat, "exons.bed").as_deref(),
        Some("sha256:e4f59e3acf34138e0c8bb04e5b89af0577e53bd1af83f1a7c489bcb6976e9a58")
    );
    assert_eq!(
        concord::pinned_digest(&cat, "1kg-30x-chr20").as_deref(),
        Some("md5:34a9a5dbfb285e7eb9d1f2881953014b")
    );
    assert_eq!(concord::pinned_digest(&cat, "no-such-input"), None);
    let inputs = vec![
        (
            "a".to_string(),
            "exons.bed".to_string(),
            "sha256:aa".to_string(),
        ),
        (
            "b".to_string(),
            "ccre.bed".to_string(),
            "sha256:bb".to_string(),
        ),
    ];
    let stamp = |tools: &str, inputs: &[(String, String, String)], cmd: &str| {
        concord::reference_stamp(tools, inputs, cmd)
    };
    let base = stamp("image gloci-bench:p1 sha256:1111", &inputs, "bedtools x");
    assert!(
        base.contains("sha256:1111") && base.contains("sha256:aa") && base.contains("bedtools x")
    );
    // Another image build, another input or another command is another key.
    assert_ne!(
        base,
        stamp("image gloci-bench:p1 sha256:2222", &inputs, "bedtools x")
    );
    let mut other = inputs.clone();
    other[1].2 = "sha256:bc".into();
    assert_ne!(
        base,
        stamp("image gloci-bench:p1 sha256:1111", &other, "bedtools x")
    );
    assert_ne!(
        base,
        stamp("image gloci-bench:p1 sha256:1111", &inputs, "bedtools y")
    );
    assert_eq!(
        base,
        stamp("image gloci-bench:p1 sha256:1111", &inputs, "bedtools x")
    );
}

/// A concordance-only results directory gets the run records the paper
/// tables read, so `cargo xtask paper` renders T4 from it.
#[test]
fn a_concordance_results_dir_renders_t4() {
    let dir = tempfile::tempdir().unwrap();
    let results = dir.path().join("9.9.9");
    std::fs::create_dir_all(&results).unwrap();
    let csv = format!(
        "{CSV_HEADER}\nsort-x,sort,sort-bed,exons.bed,default,4,4,3,0.750000,documented,sort-tie-order\nmap-x,map,bedtools,ccre.bed,default,,,,,failed,\n"
    );
    std::fs::write(results.join("concordance.csv"), csv).unwrap();
    let args = concord::ConcordArgs {
        only: None,
        gloci: None,
        no_docker: true,
        image: "unused".into(),
        data_dir: dir.path().join("no-data"),
        out: dir.path().to_path_buf(),
    };
    let tables = dir.path().join("tables");
    // Without the records the tables cannot be made.
    assert!(xtask::paper::tables::tables(&results, &tables).is_err());
    let cases = load_cases().unwrap();
    concord::write_run_records(
        &args,
        &results,
        std::path::Path::new("/bin/true"),
        &cases,
        "2026-01-01T00:00:00Z",
    )
    .unwrap();
    assert!(results.join("environment.json").is_file());
    assert!(results.join("repo-snapshot.json").is_file());
    xtask::paper::tables::tables(&results, &tables).unwrap();
    let t4 = std::fs::read_to_string(tables.join("tables/T4-concordance.csv"))
        .or_else(|_| std::fs::read_to_string(tables.join("T4-concordance.csv")))
        .unwrap();
    assert!(t4.contains("sort-x,sort,sort-bed,exons.bed,default,4,4,3,0.750000,sort-tie-order"));
    assert!(t4.contains("map-x,map,bedtools,ccre.bed,default,,,,,"));
    // Records a benchmark run wrote are kept.
    std::fs::write(
        results.join("environment.json"),
        "{\"gloci_version\": \"keep\"}\n",
    )
    .unwrap();
    concord::write_run_records(
        &args,
        &results,
        std::path::Path::new("/bin/true"),
        &cases,
        "2026-01-01T00:00:00Z",
    )
    .unwrap();
    assert!(
        std::fs::read_to_string(results.join("environment.json"))
            .unwrap()
            .contains("keep")
    );
}
