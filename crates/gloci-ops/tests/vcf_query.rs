//! `vcf_query` through the shared selection stage: the sample-level rule, thread
//! counts, VCF against BCF, sample order and errors (decision D172).
use gloci_io::{IoError, Loc, ReaderOptions};
use gloci_ops::{
    ExtractMode, ExtractOptions, FilterExpr, QueryOptions, VcfSelect, extract_table, open_selected,
    vcf_query,
};
use gloci_vcf::QueryFormat;
use std::num::NonZeroUsize;
use std::path::Path;

/// The shared expression fixture, read at run time: a crate's `include_str!` may not leave it.
fn main_text() -> String {
    std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../gloci-vcf/tests/golden/expr/main.vcf"
    ))
    .unwrap()
}

fn ro() -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: gloci_io::InvalidPolicy::Fail,
        chrom_map: None,
    }
}

fn names(v: &[&str]) -> Option<Vec<String>> {
    Some(v.iter().map(|s| s.to_string()).collect())
}

fn expr(src: &str, exclude: bool) -> Option<FilterExpr> {
    Some(FilterExpr {
        src: src.to_string(),
        exclude,
    })
}

fn query(
    path: &Path,
    select: VcfSelect,
    fmt: &str,
    print_header: bool,
    threads: Option<usize>,
) -> Result<String, IoError> {
    let input = open_selected(&path.to_string_lossy(), &ro(), &select)?;
    let o = QueryOptions {
        select,
        threads: threads.map(|t| NonZeroUsize::new(t).unwrap()),
        format: QueryFormat::parse(fmt).unwrap(),
        print_header,
    };
    let mut out = Vec::new();
    vcf_query(input, &o, &mut out)?;
    Ok(String::from_utf8(out).unwrap())
}

fn main_vcf(dir: &Path) -> std::path::PathBuf {
    gloci_testkit::fixtures::write(dir, "main.vcf", main_text().as_bytes())
}

#[test]
fn query_rule_crosses_the_workers() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = main_vcf(dir.path());
    for threads in [None, Some(4)] {
        let select = VcfSelect {
            include_samples: names(&["s2", "s3"]),
            expr: expr(r#"GT="het""#, true),
            ..VcfSelect::default()
        };
        let got = query(&vcf, select, "%POS[,%SAMPLE]\n", false, threads).unwrap();
        assert_eq!(got, "10,s2,s3\n20,s3\n30,s2,s3\n40,s2,s3\n", "{threads:?}");
    }
}

/// A small deterministic generator (xorshift64*).
struct Rng(u64);
impl Rng {
    fn next(&mut self, n: u64) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 33) % n
    }
}

fn cohort() -> String {
    let n = 120;
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##contig=<ID=2,length=100000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for i in 0..n {
        text.push_str(&format!("\tp{i}"));
    }
    text.push('\n');
    let gts = ["0/0", "0/1", "1/1", "./.", "1", "0|1"];
    for rec in 0..2000u64 {
        let (chrom, pos) = if rec < 1000 {
            (1, rec)
        } else {
            (2, rec - 1000)
        };
        text.push_str(&format!(
            "{chrom}\t{}\t.\tA\tC\t50\tPASS\t.\tGT:DP",
            100 + pos * 7
        ));
        for _ in 0..n {
            let gt = gts[rng.next(gts.len() as u64) as usize];
            text.push_str(&format!("\t{gt}:{}", 1 + rng.next(40)));
        }
        text.push('\n');
    }
    text
}

#[test]
fn query_output_is_the_same_at_every_thread_count() {
    let dir = tempfile::tempdir().unwrap();
    let text = cohort();
    let vcf = gloci_testkit::fixtures::write(dir.path(), "c.vcf", text.as_bytes());
    let bcf = gloci_testkit::fixtures::write(
        dir.path(),
        "c.bcf",
        &gloci_testkit::fixtures::vcf_to_bcf(text.as_bytes()),
    );
    let evens: Vec<String> = (0..120).step_by(2).map(|i| format!("p{i}")).collect();
    let selections = [
        VcfSelect {
            include_samples: Some(evens),
            mask: gloci_vcf::MaskRules {
                min_dp: Some(5),
                ..Default::default()
            },
            expr: expr(r#"GT="het" && FMT/DP>10"#, false),
            ..VcfSelect::default()
        },
        VcfSelect {
            expr: expr(r#"GT="alt""#, true),
            ..VcfSelect::default()
        },
    ];
    let fmt = r"%CHROM\t%POS\t%TYPE[\t%SAMPLE=%GT:%DP]\n";
    for select in selections {
        let run = |path: &Path, t: usize| query(path, select.clone(), fmt, false, Some(t)).unwrap();
        let base = run(&vcf, 1);
        assert!(base.lines().count() > 200, "{} lines", base.lines().count());
        for t in [2, 8] {
            assert_eq!(run(&vcf, t), base, "vcf at {t} threads");
        }
        for t in [1, 2, 8] {
            assert_eq!(run(&bcf, t), base, "bcf at {t} threads");
        }
    }
}

#[test]
fn query_prints_samples_in_the_order_given() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = main_vcf(dir.path());
    for threads in [None, Some(1), Some(4)] {
        let by_include = VcfSelect {
            include_samples: names(&["s3", "s1"]),
            ..VcfSelect::default()
        };
        let got = query(&vcf, by_include, "%POS[,%SAMPLE]\n", false, threads).unwrap();
        assert_eq!(got.lines().count(), 4);
        assert!(got.lines().all(|l| l.ends_with(",s3,s1")), "{got}");
        let by_include = VcfSelect {
            include_samples: names(&["s3", "s1"]),
            ..VcfSelect::default()
        };
        let got = query(&vcf, by_include, "%POS[,%SAMPLE]\n", true, threads).unwrap();
        assert_eq!(
            got.lines().next().unwrap(),
            "#[1]POS,[2]s3:SAMPLE,[3]s1:SAMPLE"
        );
        let by_exclude = VcfSelect {
            exclude_samples: vec!["s2".into()],
            ..VcfSelect::default()
        };
        let got = query(&vcf, by_exclude, "%POS[,%SAMPLE]\n", false, threads).unwrap();
        assert_eq!(got.lines().count(), 4);
        assert!(got.lines().all(|l| l.ends_with(",s1,s3")), "{got}");
    }
}

#[test]
fn query_header_names_samples_unless_the_filter_reads_them() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = main_vcf(dir.path());
    let fmt = r"%POS[\t%GT]\n";
    let with = |src: &str| VcfSelect {
        expr: expr(src, false),
        ..VcfSelect::default()
    };
    let got = query(&vcf, with("QUAL>5"), fmt, true, None).unwrap();
    assert_eq!(
        got.lines().next().unwrap(),
        "#[1]POS\t[2]s1:GT\t[3]s2:GT\t[4]s3:GT"
    );
    let got = query(&vcf, with(r#"GT="alt""#), fmt, true, None).unwrap();
    assert_eq!(got.lines().next().unwrap(), "#[1]POS\t[2]GT\t[3]GT\t[4]GT");
    let got = query(&vcf, with("QUAL>1000"), fmt, true, None).unwrap();
    assert_eq!(got, "#[1]POS\t[2]s1:GT\t[3]s2:GT\t[4]s3:GT\n");
}

#[test]
fn extract_keeps_the_site_rule() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = main_vcf(dir.path());
    let select = VcfSelect {
        expr: expr(r#"GT="het""#, true),
        ..VcfSelect::default()
    };
    let input = open_selected(&vcf.to_string_lossy(), &ro(), &select).unwrap();
    let o = ExtractOptions {
        select,
        threads: None,
        mode: ExtractMode::Info(vec!["DP".into()]),
        skip_invalid: false,
        max_mem: 1 << 30,
        tmp_dir: dir.path().to_path_buf(),
    };
    let mut out = Vec::new();
    extract_table(input, &o, &mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    let rows: Vec<&str> = text.lines().skip(1).collect();
    assert_eq!(rows, ["1\t30\tx1\tAT\tA\t3"], "{text}");
}

#[test]
fn query_error_names_the_record() {
    let dir = tempfile::tempdir().unwrap();
    let text = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tC\t.\t.\tDP=5\n1\t20\t.\tA\tC\t.\t.\tDP=abc\n";
    let vcf = gloci_testkit::fixtures::write(dir.path(), "bad.vcf", text.as_bytes());
    let err = query(&vcf, VcfSelect::default(), "%DP\n", false, None).unwrap_err();
    match err {
        IoError::Parse { loc, msg, .. } => {
            assert_eq!(loc, Loc::Record(2));
            assert!(msg.contains("INFO/DP"), "{msg}");
        }
        e => panic!("{e:?}"),
    }
}

#[test]
fn query_keeps_the_given_order_with_exclusions_and_repeats() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = main_vcf(dir.path());
    let run = |inc: &[&str], exc: &[&str]| {
        let select = VcfSelect {
            include_samples: names(inc),
            exclude_samples: exc.iter().map(|s| s.to_string()).collect(),
            ..VcfSelect::default()
        };
        query(&vcf, select, "%POS[,%SAMPLE]\n", false, Some(4)).unwrap()
    };
    assert!(
        run(&["s3", "s1"], &["s1"])
            .lines()
            .all(|l| l.ends_with(",s3") && !l.contains("s1"))
    );
    assert!(
        run(&["s3", "s2", "s1"], &["s2"])
            .lines()
            .all(|l| l.ends_with(",s3,s1"))
    );
    assert!(
        run(&["s3", "s1", "s3"], &[])
            .lines()
            .all(|l| l.ends_with(",s3,s1,s3"))
    );
}

#[test]
fn a_sample_subscript_with_a_sample_option_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = main_vcf(dir.path());
    for exclude in [false, true] {
        let select = VcfSelect {
            include_samples: names(&["s2", "s3"]),
            expr: expr(r#"GT[0]="het""#, exclude),
            ..VcfSelect::default()
        };
        let err = query(&vcf, select, "%POS\n", false, None).unwrap_err();
        assert!(matches!(err, IoError::Usage(_)), "{err}");
        assert!(
            err.to_string().contains("GT[0] in --include/--exclude"),
            "{err}"
        );
    }
    // Without a sample option the numbering agrees with bcftools.
    let select = VcfSelect {
        expr: expr(r#"GT[0]="het""#, false),
        ..VcfSelect::default()
    };
    query(&vcf, select, "%POS\n", false, None).unwrap();
}
