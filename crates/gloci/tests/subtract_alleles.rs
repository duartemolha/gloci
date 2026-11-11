//! `gloci subtract --match alleles`: a VCF/BCF A record is dropped only by
//! a B record that overlaps it and shares an allele; equal to `intersect
//! --mode none --match alleles`; usage errors.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::Path;

const HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"length\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn vcf(recs: &[&str]) -> String {
    let mut t = HEAD.to_string();
    for r in recs {
        // a symbolic ALT needs a span
        let info = if r.contains('<') { "SVLEN=-10" } else { "." };
        t.push_str(&format!("chr1\t{r}\t.\t.\t{info}\n"));
    }
    t
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// The record lines of VCF text, the header dropped.
fn records(text: &str) -> String {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect()
}

fn to_bcf(dir: &Path, name: &str, vcf_path: &str) -> String {
    let o = s(&dir.join(name));
    out(&["view", vcf_path, "-o", &o]);
    o
}

/// The records of a BCF or VCF as VCF text.
fn text_of(path: &str) -> String {
    records(&out(&["view", path, "-O", "vcf"]))
}

#[test]
fn subtract_match_alleles_table() {
    let dir = tempdir();
    let a = vcf(&[
        "1\t.\tN\tNA",
        "100\t.\tA\tG",
        "200\t.\tC\tT",
        "300\t.\tG\tGA",
        "400\t.\tTAC\tT",
        "500\t.\tA\tC,G",
        "600\t.\tN\t<DEL>",
    ]);
    let b = vcf(&[
        "1\t.\tN\tNA",
        "100\t.\tA\tG",
        "200\t.\tC\tA",
        "300\t.\tG\tGT",
        "400\t.\tTAC\tTA",
        "500\t.\tA\tG",
        "600\t.\tN\t<DEL>",
    ]);
    let (af, bf) = (
        s(&plain_file(dir.path(), "a.vcf", &a)),
        s(&plain_file(dir.path(), "b.vcf", &b)),
    );
    let kept = out(&["subtract", "-a", &af, "-b", &bf, "--match", "alleles"]);
    let pos: Vec<&str> = kept
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(pos, ["200", "300", "400", "600"], "{kept}");
    // the default match is by overlap, and drops every record with a B on it
    let plain = out(&["subtract", "-a", &af, "-b", &bf]);
    assert_eq!(records(&plain), "", "{plain}");
    // --whole changes nothing for a VCF A
    let whole = out(&[
        "subtract", "-a", &af, "-b", &bf, "--match", "alleles", "--whole",
    ]);
    assert_eq!(records(&whole), records(&kept));
    // BCF inputs give the same records
    let (abcf, bbcf) = (
        to_bcf(dir.path(), "a.bcf", &af),
        to_bcf(dir.path(), "b.bcf", &bf),
    );
    let o = s(&dir.path().join("o.bcf"));
    out(&[
        "subtract", "-a", &abcf, "-b", &bbcf, "--match", "alleles", "-o", &o,
    ]);
    assert_eq!(text_of(&o), records(&kept));
}

#[test]
fn subtract_alleles_equals_intersect_none() {
    let dir = tempdir();
    let sets: &[(&[&str], &[&[&str]])] = &[
        // The same insertion written with two anchors: A>AG at 100 and C>GC at 101 both trim
        // to an insertion of G after 100. A deletion that trims elsewhere (GAA>GA trims to a
        // deletion of A after 101, GA>G after 100) is not the same variant.
        (
            &["100\t.\tA\tAG", "200\t.\tGA\tG", "300\t.\tC\tT"],
            &[&["101\t.\tC\tGC", "200\t.\tGAA\tGA", "300\t.\tC\tG"]],
        ),
        // several B files, a multi-allelic A, lower case, symbolic
        (
            &[
                "10\t.\tA\tC,G",
                "20\t.\tc\tt",
                "30\t.\tN\t<DEL>",
                "40\t.\tA\t*",
                "50\t.\tA\tATT",
            ],
            &[
                &["10\t.\tA\tG", "30\t.\tN\t<DEL>"],
                &["20\t.\tC\tT"],
                &["50\t.\tA\tATT"],
            ],
        ),
        (
            &["1\t.\tN\tNA", "5\t.\tA\tC"],
            &[&["1\t.\tN\tNA", "6\t.\tA\tC"]],
        ),
    ];
    // The positions each set's subtract drops (the rest are kept): a symbolic ALT (30)
    // never matches, and `*` (40) is not an allele.
    let dropped: [&[&str]; 3] = [&["100"], &["10", "20", "50"], &["1"]];
    for (i, (a, bs)) in sets.iter().enumerate() {
        let af = s(&plain_file(dir.path(), &format!("a{i}.vcf"), &vcf(a)));
        let bfs: Vec<String> = bs
            .iter()
            .enumerate()
            .map(|(j, b)| s(&plain_file(dir.path(), &format!("b{i}_{j}.vcf"), &vcf(b))))
            .collect();
        let abcf = to_bcf(dir.path(), &format!("a{i}.bcf"), &af);
        let bbcfs: Vec<String> = bfs
            .iter()
            .enumerate()
            .map(|(j, b)| to_bcf(dir.path(), &format!("b{i}_{j}.bcf"), b))
            .collect();
        for (aa, bb) in [(&af, &bfs), (&abcf, &bbcfs), (&af, &bbcfs), (&abcf, &bfs)] {
            for extra in [&[][..], &["--unsorted"][..]] {
                for t in ["1", "2", "8"] {
                    let run = |cmd: &str| {
                        let mut v = vec![cmd, "-a", aa.as_str(), "--match", "alleles", "-t", t];
                        for b in bb {
                            v.push("-b");
                            v.push(b);
                        }
                        if cmd == "intersect" {
                            v.extend(["--mode", "none"]);
                        }
                        v.extend_from_slice(extra);
                        let o = s(&dir.path().join(format!("o_{cmd}.out")));
                        let ext = if aa.ends_with(".bcf") { "bcf" } else { "vcf" };
                        let o = format!("{o}.{ext}");
                        v.extend(["-o", &o]);
                        out(&v);
                        text_of(&o)
                    };
                    let kept = run("subtract");
                    assert_eq!(
                        kept,
                        run("intersect"),
                        "set {i} {aa} {bb:?} {extra:?} -t {t}"
                    );
                    let pos = |l: &str| l.split('\t').nth(1).unwrap().to_string();
                    let want: Vec<String> = a
                        .iter()
                        .map(|r| pos(&format!("chr1\t{r}")))
                        .filter(|p| !dropped[i].contains(&p.as_str()))
                        .collect();
                    let got: Vec<String> = kept.lines().map(pos).collect();
                    assert_eq!(got, want, "set {i}: kept positions");
                }
            }
        }
    }
}

#[test]
fn subtract_alleles_usage_errors() {
    let dir = tempdir();
    let v = s(&plain_file(dir.path(), "a.vcf", &vcf(&["100\t.\tA\tG"])));
    let bed = s(&plain_file(dir.path(), "x.bed", "chr1\t0\t1000\n"));
    let sub = |a: &str, b: &str, extra: &[&str]| {
        let mut args = vec!["subtract", "-a", a, "-b", b, "--match", "alleles"];
        args.extend_from_slice(extra);
        gloci().args(args).assert().code(2)
    };
    sub(&bed, &v, &[]).stderr(contains("must be VCF or BCF").and(contains("--match overlap")));
    sub(&v, &bed, &[]).stderr(contains("must be VCF or BCF"));
    let refused: &[(&[&str], &str)] = &[
        (&["-f", "0.5"], "-f"),
        (&["-F", "0.5"], "-F"),
        (&["-r"], "-r"),
        (&["-e"], "-e"),
        (&["--min-bases", "1"], "--min-bases"),
        (&["-s"], "-s/-S"),
        (&["-S"], "-s/-S"),
        (&["-w", "5"], "-w/--left/--right"),
        (&["--left", "5"], "-w/--left/--right"),
        (&["--right", "5"], "-w/--left/--right"),
        (&["--sum-overlaps"], "--sum-overlaps"),
        (&["--sum-overlaps", "--whole"], "--sum-overlaps"),
        (&["--split"], "--split"),
        (&["--inclusive-points"], "--inclusive-points"),
    ];
    for (extra, flag) in refused {
        sub(&v, &v, extra).stderr(
            contains(*flag)
                .and(contains("--match alleles"))
                .and(contains("--match overlap")),
        );
    }
}
