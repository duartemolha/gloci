use gloci_io::{Format, Header, VcfHeader};
use gloci_ops::agg::{Accumulators, parse_aggs, resolve};

fn s(x: &str) -> String {
    x.to_string()
}

/// Aggregate `lines` (BED-ish, value in column 5) with `spec`.
fn run_with(spec: &str, lines: &[&str], empty: &str, delim: &str) -> String {
    let specs = parse_aggs(&[s(spec)]).unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s(empty), s(delim));
    for (i, l) in lines.iter().enumerate() {
        acc.add(
            format!("chr1\t0\t10\tn\t{l}").as_bytes(),
            5,
            "b.bed",
            i as u64 + 1,
        )
        .unwrap();
    }
    let mut out = Vec::new();
    acc.finish_into(&mut out);
    String::from_utf8(out)
        .unwrap()
        .strip_prefix('\t')
        .unwrap()
        .to_string()
}

fn run(op: &str, lines: &[&str]) -> String {
    run_with(&format!("5:{op}"), lines, ".", ",")
}

#[test]
fn parse() {
    let v = parse_aggs(&[s("5:mean,5:max"), s("4:distinct")]).unwrap();
    let labels: Vec<_> = v.iter().map(|a| a.label.as_str()).collect();
    assert_eq!(labels, ["5:mean", "5:max", "4:distinct"]);
}

#[test]
fn parse_errors() {
    assert!(parse_aggs(&[s("5")]).is_err());
    let e = parse_aggs(&[s("5:avg")]).unwrap_err().0;
    assert!(e.contains("unknown op"), "{e}");
    assert!(parse_aggs(&[s("0:sum")]).is_err());
    assert!(parse_aggs(&[s("0:count")]).is_ok());
    for op in [
        "stdev_pop",
        "var_pop",
        "cv",
        "mad",
        "quantile:0.5",
        "tmean:0.1-0.9",
        "wmean",
        "mode",
        "antimode",
        "bases_uniq",
        "frac_a",
        "frac_b",
        "min_record",
        "max_record",
        "records",
    ] {
        let e = parse_aggs(&[format!("5:{op}")]).unwrap_err().0;
        let name = op.split(':').next().unwrap();
        assert!(e.contains(&format!("{name} is not available yet")), "{e}");
        assert!(!e.contains("unknown op"), "{e}");
    }
}

#[test]
fn values() {
    let v = ["1", "2", "3", "4"];
    for (op, want) in [
        ("sum", "10"),
        ("mean", "2.5"),
        ("median", "2.5"),
        ("min", "1"),
        ("max", "4"),
        ("stdev", "1.2909944487358056"),
        ("var", "1.6666666666666667"),
        ("count", "4"),
        ("count_distinct", "4"),
        ("first", "1"),
        ("last", "4"),
    ] {
        assert_eq!(run(op, &v), want, "{op}");
    }
    assert_eq!(run("bases", &v), "20");
    assert_eq!(run("indicator", &v), "1");
    assert_eq!(run("indicator", &[]), "0");
    assert_eq!(run("median", &["5", "1", "3"]), "3");
}

#[test]
fn missing_and_empty() {
    assert_eq!(run("sum", &["1.5", "2"]), "3.5");
    assert_eq!(run("mean", &[".", "3"]), "3");
    assert_eq!(run("count", &[]), "0");
    assert_eq!(run("count_distinct", &[]), "0");
    assert_eq!(run("sum", &[]), ".");
    assert_eq!(run_with("5:sum", &[], "NA", ","), "NA");
    assert_eq!(run("stdev", &["7"]), ".");
    assert_eq!(run("var", &["7"]), ".");
    assert_eq!(run("max", &["1.5", "-2"]), "1.5");
    // Integer sums are exact beyond f64 precision.
    assert_eq!(run("sum", &["9007199254740993", "1"]), "9007199254740994");
}

#[test]
fn reuse_resets() {
    let ops = [
        "count",
        "count_distinct",
        "sum",
        "mean",
        "median",
        "min",
        "max",
        "stdev",
        "var",
        "collapse",
        "distinct",
        "distinct_sorted",
        "first",
        "last",
        "bases",
        "indicator",
    ];
    let specs = parse_aggs(&[ops
        .iter()
        .map(|o| format!("5:{o}"))
        .collect::<Vec<_>>()
        .join(",")])
    .unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(","));
    let mut round = |vals: &[&str]| {
        for (i, v) in vals.iter().enumerate() {
            acc.add(
                format!("c\t0\t1\tn\t{v}").as_bytes(),
                3,
                "b.bed",
                i as u64 + 1,
            )
            .unwrap();
        }
        let mut out = Vec::new();
        acc.finish_into(&mut out);
        String::from_utf8(out).unwrap()
    };
    assert_eq!(
        round(&["1", "2", "2"]),
        "\t3\t2\t5\t1.6666666666666667\t2\t1\t2\t0.5773502691896257\t0.3333333333333333\t1,2,2\t1,2\t1,2\t1\t2\t9\t1"
    );
    assert_eq!(
        round(&["10", "30"]),
        "\t2\t2\t40\t20\t20\t10\t30\t14.142135623730951\t200\t10,30\t10,30\t10,30\t10\t30\t6\t1"
    );
    assert_eq!(
        round(&[]),
        "\t0\t0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0\t0"
    );
}

#[test]
fn mixed_int_and_float() {
    assert_eq!(run("min", &["1", "2.5"]), "1");
    assert_eq!(run("max", &["1", "2.5"]), "2.5");
    assert_eq!(run("mean", &["1", "2.5"]), "1.75");
    assert_eq!(run("sum", &["1", "2.5"]), "3.5");
}

#[test]
fn text_ops() {
    let v = ["b", "a", "b"];
    assert_eq!(run("distinct", &v), "b,a");
    assert_eq!(run("distinct_sorted", &v), "a,b");
    assert_eq!(run("collapse", &v), "b,a,b");
    assert_eq!(run("count_distinct", &v), "2");
    assert_eq!(run_with("5:collapse", &v, ".", ";"), "b;a;b");
}

#[test]
fn non_numeric() {
    let specs = parse_aggs(&[s("5:sum")]).unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(","));
    let e = acc
        .add(b"c\t0\t1\tn\tabc", 1, "b.bed", 7)
        .unwrap_err()
        .to_string();
    assert!(
        e.contains("b.bed line 7: column 5 value \"abc\" is not a number (5:sum)"),
        "{e}"
    );
    let e = acc
        .add(b"c\t0\t1\tn\tnan", 1, "b.bed", 12_345)
        .unwrap_err()
        .to_string();
    assert!(e.contains("line 12,345"), "{e}");
}

#[test]
fn header_names() {
    let h = Header::Lines(b"#chrom\tstart\tend\tscore\n".to_vec());
    let specs = parse_aggs(&[s("score:max")]).unwrap();
    let r = resolve(&specs, Format::Bed, &h).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(","));
    acc.add(b"c\t0\t1\t9", 1, "b.bed", 1).unwrap();
    let mut out = Vec::new();
    acc.finish_into(&mut out);
    assert_eq!(out, b"\t9");
    let specs = parse_aggs(&[s("nope:max")]).unwrap();
    let e = resolve(&specs, Format::Bed, &h).unwrap_err().0;
    assert!(e.contains("chrom, start, end, score"), "{e}");
}

#[test]
fn vcf_names() {
    let header = Header::Vcf(
        VcfHeader::parse(
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec(),
        )
        .unwrap(),
    );
    let recs: [&[u8]; 3] = [
        b"1\t5\trs1\tA\tC\t30\tPASS\tDP=10;DB",
        b"1\t6\t.\tA\tC\t50\tq10\tDP=20;AF=0.5",
        b"1\t7\trs3\tA\tC\t.\t.\t.",
    ];
    let go = |spec: &str| {
        let specs = parse_aggs(&[s(spec)]).unwrap();
        let r = resolve(&specs, Format::Vcf, &header).unwrap();
        let mut acc = Accumulators::new(r, s("."), s(","));
        for (i, l) in recs.iter().enumerate() {
            acc.add(l, 1, "a.vcf", i as u64 + 1).unwrap();
        }
        let mut out = Vec::new();
        acc.finish_into(&mut out);
        String::from_utf8(out).unwrap()
    };
    assert_eq!(go("INFO/DP:mean"), "\t15");
    assert_eq!(go("QUAL:max"), "\t50");
    assert_eq!(go("INFO/DB:sum"), "\t1");
    assert_eq!(go("ID:collapse"), "\trs1,rs3");
    assert_eq!(go("FILTER:distinct_sorted"), "\tPASS,q10");
    assert_eq!(go("INFO/D:sum"), "\t.");
    let specs = parse_aggs(&[s("FOO:sum")]).unwrap();
    assert!(resolve(&specs, Format::Vcf, &header).is_err());
}
