use gloci_core::Interval;
use gloci_io::{Format, Header, VcfHeader};
use gloci_ops::agg::{Accumulators, HitInfo, parse_aggs, refuse_two_input_ops, resolve};

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
    // arguments
    assert!(parse_aggs(&[s("5:quantile:0.25")]).is_ok());
    assert!(parse_aggs(&[s("5:tmean:0.2-0.8")]).is_ok());
    assert!(parse_aggs(&[s("0:records,0:bases_uniq,0:frac_a,0:frac_b")]).is_ok());
    for (spec, want) in [
        (
            "5:quantile:1.5",
            "quantile `1.5` must be a number from 0 to 1",
        ),
        ("5:quantile:-0.1", "must be a number from 0 to 1"),
        ("5:quantile:x", "must be a number from 0 to 1"),
        ("5:quantile", "quantile needs a quantile"),
        ("5:tmean:0.8-0.2", "LO below HI"),
        ("5:tmean:0.5-0.5", "LO below HI"),
        ("5:tmean:0.1-1.5", "must be a number from 0 to 1"),
        ("5:tmean:0.1", "tmean needs a range LO-HI"),
        ("5:tmean", "tmean needs a range LO-HI"),
        ("5:mad:3", "mad takes no argument"),
        ("0:cv", "column 0"),
        ("0:min_record", "column 0"),
    ] {
        let e = parse_aggs(&[s(spec)]).unwrap_err().0;
        assert!(e.contains(want), "{spec}: {e}");
    }
}

#[test]
fn two_input_ops_are_refused_without_a_second_input() {
    for op in [
        "wmean",
        "records",
        "min_record",
        "max_record",
        "frac_a",
        "frac_b",
        "bases_uniq",
        "bases",
    ] {
        let specs = parse_aggs(&[format!("5:{op}")]).unwrap();
        let e = refuse_two_input_ops(&specs, "merge").unwrap_err().0;
        assert!(
            e.contains(&format!("-c 5:{op}"))
                && e.contains("merge has none")
                && e.contains("use map"),
            "{e}"
        );
    }
    let ok = parse_aggs(&[s("5:mean,5:mode,5:quantile:0.5,5:cv,4:collapse")]).unwrap();
    assert!(refuse_two_input_ops(&ok, "merge").is_ok());
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
        ("stdev_pop", "1.118033988749895"),
        ("var_pop", "1.25"),
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
    assert_eq!(run("stdev_pop", &["7"]), "0");
    assert_eq!(run("var_pop", &["7"]), "0");
    assert_eq!(run("var_pop", &[]), ".");
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
        "stdev_pop",
        "var_pop",
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
        "\t3\t2\t5\t1.6666666666666667\t2\t1\t2\t0.5773502691896257\t0.3333333333333333\t0.4714045207910317\t0.2222222222222222\t1,2,2\t1,2\t1,2\t1\t2\t9\t1"
    );
    assert_eq!(
        round(&["10", "30"]),
        "\t2\t2\t40\t20\t20\t10\t30\t14.142135623730951\t200\t10\t100\t10,30\t10,30\t10,30\t10\t30\t6\t1"
    );
    assert_eq!(
        round(&[]),
        "\t0\t0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t0\t0"
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

// ---- the D64 ops ----

/// Independent references, written from the D64 formulas.
fn sample_stdev(v: &[f64]) -> f64 {
    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
}

fn cv_ref(v: &[f64]) -> String {
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    format!("{}", sample_stdev(v) / mean)
}

/// Mean of the sorted values whose rank fraction i/(n-1) lies in [lo, hi].
fn tmean_ref(v: &[f64], lo: f64, hi: f64) -> String {
    let mut v = v.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    let kept: Vec<f64> = (0..n)
        .filter(|&i| {
            let f = i as f64 / (n - 1) as f64;
            lo <= f && f <= hi
        })
        .map(|i| v[i])
        .collect();
    format!("{}", kept.iter().sum::<f64>() / kept.len() as f64)
}

#[test]
fn spread_and_quantiles() {
    let v = ["1", "2", "3", "4", "10"];
    let f = [1.0, 2.0, 3.0, 4.0, 10.0];
    assert_eq!(run("mad", &v), "1");
    assert_eq!(run("quantile:0.25", &v), "2");
    assert_eq!(run("quantile:0", &v), "1");
    assert_eq!(run("quantile:1", &v), "10");
    assert_eq!(run("quantile:0.5", &v), "3");
    // h = 4 * 0.9 sits between the 4th and 5th values (4 and 10)
    let h = 4.0 * 0.9_f64;
    assert_eq!(
        run("quantile:0.9", &v),
        format!("{}", 4.0 + (h - h.floor()) * 6.0)
    );
    assert_eq!(run("cv", &v), cv_ref(&f));
    assert_eq!(run("tmean:0.2-0.8", &v), tmean_ref(&f, 0.2, 0.8));
    assert_eq!(run("tmean:0.2-0.8", &v), "3");
    assert_eq!(run("tmean:0-1", &v), "4");
    // unsorted input and floats
    let g = ["2.5", "0.5", "1.5", "9.5"];
    assert_eq!(run("cv", &g), cv_ref(&[2.5, 0.5, 1.5, 9.5]));
    assert_eq!(
        run("tmean:0.2-0.8", &g),
        tmean_ref(&[2.5, 0.5, 1.5, 9.5], 0.2, 0.8)
    );
    assert_eq!(run("mad", &["1", "2"]), "0.5");
}

#[test]
fn spread_edge_cases() {
    assert_eq!(run("cv", &["5"]), ".");
    assert_eq!(run("cv", &["-1", "1"]), ".");
    assert_eq!(run("cv", &[]), ".");
    assert_eq!(run("mad", &["7"]), "0");
    assert_eq!(run("mad", &[]), ".");
    assert_eq!(run("quantile:0.3", &["7"]), "7");
    assert_eq!(run("quantile:0.3", &[]), ".");
    assert_eq!(run("tmean:0.1-0.9", &["7"]), "7");
    // no value has a rank fraction inside [0.1, 0.2] among five
    assert_eq!(run("tmean:0.1-0.2", &["1", "2", "3", "4", "5"]), ".");
    assert_eq!(run_with("5:cv", &["5"], "NA", ","), "NA");
}

#[test]
fn mode_and_antimode() {
    let v = ["b", "a", "b", "a", "c"];
    // a and b tie at two: the first seen wins
    assert_eq!(run("mode", &v), "b");
    // c is the only one seen once
    assert_eq!(run("antimode", &v), "c");
    // all tie: first seen for both
    assert_eq!(run("mode", &["z", "y"]), "z");
    assert_eq!(run("antimode", &["z", "y"]), "z");
    assert_eq!(run("mode", &["3", "1", "1", ".", "3", "3"]), "3");
    assert_eq!(run("mode", &[]), ".");
    assert_eq!(run("antimode", &[]), ".");
}

fn hit<'a>(bases: u64, a_len: u64, b_len: u64, shared: &'a [Interval]) -> HitInfo<'a> {
    HitInfo {
        bases,
        a_len,
        b_len,
        shared,
    }
}

fn iv(start: u64, end: u64) -> Interval {
    Interval { start, end }
}

/// Fold `(value, hit)` pairs (value in column 5) and return the cell.
fn run_hits(spec: &str, hits: &[(&str, HitInfo<'_>)]) -> String {
    let specs = parse_aggs(&[s(spec)]).unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(","));
    for (i, (v, h)) in hits.iter().enumerate() {
        acc.add_hit(
            format!("chr1\t0\t10\tn\t{v}").as_bytes(),
            h,
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

#[test]
fn wmean_weights_by_overlap() {
    let (s1, s2) = ([iv(0, 10)], [iv(10, 40)]);
    let hits = [("10", hit(10, 100, 10, &s1)), ("20", hit(30, 100, 30, &s2))];
    // (10*10 + 20*30) / 40
    assert_eq!(run_hits("5:wmean", &hits), "17.5");
    // a hit with no shared bases weighs nothing
    let none = [("10", hit(0, 100, 5, &[])), ("20", hit(0, 100, 5, &[]))];
    assert_eq!(run_hits("5:wmean", &none), ".");
    assert_eq!(run_hits("5:wmean", &[]), ".");
    let one = [("10", hit(0, 100, 5, &[])), ("30", hit(4, 100, 5, &s1))];
    assert_eq!(run_hits("5:wmean", &one), "30");
}

#[test]
fn bases_uniq_is_the_union_of_overlaps() {
    // overlapping B: [0,10) and [5,15) share 5 bases of A
    let (x, y) = ([iv(0, 10)], [iv(5, 15)]);
    let hits = [("1", hit(10, 20, 10, &x)), ("1", hit(10, 20, 10, &y))];
    assert_eq!(run_hits("5:bases", &hits), "20");
    assert_eq!(run_hits("5:bases_uniq", &hits), "15");
    assert_eq!(run_hits("5:frac_a", &hits), "0.75");
    assert_eq!(run_hits("5:frac_b", &hits), "1");
    // contained and out-of-order shared intervals
    let (c1, c2) = ([iv(8, 9)], [iv(0, 20)]);
    let hits = [("1", hit(1, 20, 1, &c1)), ("1", hit(20, 20, 20, &c2))];
    assert_eq!(run_hits("5:bases_uniq", &hits), "20");
    // no hits: 0 bases, no fraction
    assert_eq!(run_hits("5:bases_uniq", &[]), "0");
    assert_eq!(run_hits("5:frac_a", &[]), ".");
    assert_eq!(run_hits("5:frac_b", &[]), ".");
}

#[test]
fn frac_a_and_frac_b_edge_cases() {
    // a point A has no length to take a fraction of
    let p = [("1", hit(0, 0, 10, &[]))];
    assert_eq!(run_hits("5:frac_a", &p), ".");
    assert_eq!(run_hits("5:bases_uniq", &p), "0");
    // frac_b is summed overlap over summed B length
    let (x, y) = ([iv(0, 5)], [iv(20, 22)]);
    let hits = [("1", hit(5, 30, 10, &x)), ("1", hit(2, 30, 30, &y))];
    assert_eq!(run_hits("5:frac_b", &hits), format!("{}", 7.0 / 40.0));
    // a point B has no length
    assert_eq!(run_hits("5:frac_b", &[("1", hit(0, 30, 0, &[]))]), ".");
    // a window-only hit overlaps nothing but is still a hit
    assert_eq!(run_hits("5:frac_a", &[("1", hit(0, 30, 10, &[]))]), "0");
}

#[test]
fn records_keep_whole_b_records() {
    let specs = parse_aggs(&[s("0:records,5:min_record,5:max_record")]).unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(";"));
    for (i, l) in [
        "chr1\t0\t10\tx\t5",
        "chr1\t5\t15\ty\t2",
        "chr1\t6\t9\tz\t.",
        "chr1\t7\t8\tw\t5",
    ]
    .iter()
    .enumerate()
    {
        acc.add(l.as_bytes(), 1, "b.bed", i as u64 + 1).unwrap();
    }
    let mut out = Vec::new();
    acc.finish_into(&mut out);
    // tabs become `|`; ties keep the first record; a missing value is skipped
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "\tchr1|0|10|x|5;chr1|5|15|y|2;chr1|6|9|z|.;chr1|7|8|w|5\tchr1|5|15|y|2\tchr1|0|10|x|5"
    );
    let mut out = Vec::new();
    acc.finish_into(&mut out);
    assert_eq!(String::from_utf8(out).unwrap(), "\t.\t.\t.");
}

#[test]
fn record_ops_need_numbers() {
    let specs = parse_aggs(&[s("5:max_record")]).unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(","));
    let e = acc
        .add(b"c\t0\t1\tn\tabc", 1, "b.bed", 4)
        .unwrap_err()
        .to_string();
    assert!(
        e.contains("b.bed line 4: column 5 value \"abc\" is not a number"),
        "{e}"
    );
}

#[test]
fn new_ops_reset_between_a_records() {
    let spec = "5:cv,5:mad,5:quantile:0.5,5:tmean:0-1,5:wmean,5:mode,5:antimode,5:bases_uniq,5:frac_a,5:frac_b,5:records,5:min_record,5:max_record";
    let specs = parse_aggs(&[s(spec)]).unwrap();
    let r = resolve(&specs, Format::Bed, &Header::None).unwrap();
    let mut acc = Accumulators::new(r, s("."), s(","));
    let sh = [iv(0, 4)];
    let mut cell = |rows: &[&str]| {
        for (i, v) in rows.iter().enumerate() {
            acc.add_hit(
                format!("c\t0\t1\tn\t{v}").as_bytes(),
                &hit(4, 8, 8, &sh),
                "b.bed",
                i as u64 + 1,
            )
            .unwrap();
        }
        let mut out = Vec::new();
        acc.finish_into(&mut out);
        String::from_utf8(out).unwrap()
    };
    let first = cell(&["1", "3"]);
    assert_eq!(
        first,
        format!(
            "\t{}\t1\t2\t2\t2\t1\t1\t4\t0.5\t0.5\tc|0|1|n|1,c|0|1|n|3\tc|0|1|n|1\tc|0|1|n|3",
            cv_ref(&[1.0, 3.0])
        )
    );
    assert_eq!(cell(&[]), "\t.\t.\t.\t.\t.\t.\t.\t0\t.\t.\t.\t.\t.");
}
