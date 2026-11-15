//! Tests for `gloci vcf compare` (decision D163): the discordance files match
//! the vcftools 0.1.16 files committed in `gloci-vcf/tests/golden/compare`
//! (made by its `regen.sh`, compared by position), the differences are pinned,
//! the sample union and the sorted-input rules hold, and the output is the same
//! at every thread count.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gloci-vcf/tests/golden/compare")
}

fn g(name: &str) -> String {
    golden_dir().join(name).to_string_lossy().into_owned()
}

fn golden(name: &str) -> String {
    // vcftools prints 0/0 as nan or -nan by CPU; gloci prints nan.
    std::fs::read_to_string(golden_dir().join(name))
        .unwrap()
        .replace("-nan", "nan")
}

const EXTS: [&str; 3] = ["diff.sites", "diff.indv", "diff.discordance_matrix"];

/// Run `vcf compare A B` into a fresh directory; the prefix is returned with it.
fn compare(a: &str, b: &str, rest: &[&str]) -> (tempfile::TempDir, String) {
    let d = tempdir();
    let prefix = d.path().join("o").to_string_lossy().into_owned();
    gloci()
        .args(["vcf", "compare", a, b, "--out-prefix", &prefix])
        .args(rest)
        .assert()
        .success();
    (d, prefix)
}

fn read(prefix: &str, ext: &str) -> String {
    std::fs::read_to_string(format!("{prefix}.{ext}")).unwrap()
}

fn vcf(samples: &[&str], rows: &[&str]) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for n in samples {
        s.push('\t');
        s.push_str(n);
    }
    s.push('\n');
    for r in rows {
        s.push_str(r);
        s.push('\n');
    }
    s
}

/// `CHROM POS REF ALT GT...` as a record line.
fn rec(chrom: &str, pos: u64, r: &str, alt: &str, gts: &[&str]) -> String {
    format!(
        "{chrom}\t{pos}\t.\t{r}\t{alt}\t50\tPASS\t.\tGT\t{}",
        gts.join("\t")
    )
}

fn write(d: &Path, name: &str, text: &str) -> String {
    plain_file(d, name, text).to_string_lossy().into_owned()
}

#[test]
fn compare_matches_vcftools_by_position() {
    for (set, a, b, extra) in [
        ("main", "a.vcf", "b.vcf", vec![]),
        ("ds", "diffsamples-a.vcf", "diffsamples-b.vcf", vec![]),
        (
            "dsmap",
            "diffsamples-a.vcf",
            "diffsamples-b.vcf",
            vec!["--sample-map", "map.txt"],
        ),
    ] {
        let mut rest = vec!["--match", "position"];
        let map = g("map.txt");
        for e in &extra {
            rest.push(if *e == "map.txt" { &map } else { e });
        }
        let (_d, prefix) = compare(&g(a), &g(b), &rest);
        for ext in EXTS {
            assert_eq!(
                read(&prefix, ext),
                golden(&format!("{set}.{ext}")),
                "{set}.{ext}"
            );
        }
    }
}

#[test]
fn compare_vcftools_spellings_and_single_files() {
    let (d, prefix) = compare(
        &g("a.vcf"),
        &g("b.vcf"),
        &["--match", "position", "--diff-indv-discordance"],
    );
    assert_eq!(read(&prefix, "diff.indv"), golden("main.diff.indv"));
    let names: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(names, ["o.diff.indv"]);
    let (d, prefix) = compare(
        &g("a.vcf"),
        &g("b.vcf"),
        &[
            "--match",
            "position",
            "--site-discordance",
            "--discordance-matrix",
        ],
    );
    assert_eq!(read(&prefix, "diff.sites"), golden("main.diff.sites"));
    assert_eq!(
        read(&prefix, "diff.discordance_matrix"),
        golden("main.diff.discordance_matrix")
    );
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 2);
}

#[test]
fn compare_allele_matching_differs_from_vcftools() {
    // --match alleles is the default: a record pairs with one that shares an
    // allele (chr1:60, AT>A against ATG>AG), and a different ALT, or an N
    // REF, at one POS is not the same site (chr1:50, chr1:75), where vcftools
    // pairs them by position.
    let (_d, prefix) = compare(&g("a.vcf"), &g("b.vcf"), &[]);
    let alleles = read(&prefix, "diff.sites");
    assert_eq!(
        alleles,
        "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         chr1\t10\tB\t1\t3\t1\t0.333333\n\
         chr1\t20\tB\t1\t2\t0\t0\n\
         chr1\t30\tB\t0\t3\t1\t0.333333\n\
         chr1\t40\tB\t1\t2\t0\t0\n\
         chr1\t60\tB\t0\t3\t0\t0\n\
         chr1\t70\tB\t1\t3\t1\t0.333333\n\
         chr1\t80\tB\t1\t3\t0\t0\n\
         chr1\t90\t1\t0\t0\t0\tnan\n\
         chr1\t95\t2\t0\t0\t0\tnan\n\
         chr1\t100\tB\t1\t3\t2\t0.666667\n\
         chr1\t110\t1\t0\t0\t0\tnan\n\
         chr1\t111\t2\t0\t0\t0\tnan\n\
         chr2\t5\tB\t1\t3\t1\t0.333333\n\
         chr2\t15\t1\t0\t0\t0\tnan\n\
         chr2\t25\t2\t0\t0\t0\tnan\n"
    );
    let want = golden("main.diff.sites");
    assert_ne!(alleles, want);
    // vcftools has the chr1:50 and chr1:75 rows, gloci does not
    assert!(want.contains("chr1\t50\tB\t0\t3\t1\t"));
    assert!(want.contains("chr1\t75\tB\t1\t3\t1\t"));
    assert!(!alleles.contains("chr1\t50\t") && !alleles.contains("chr1\t75\t"));
    // and chr1:60 (a different REF) is a row in gloci, not in vcftools
    assert!(!want.contains("chr1\t60\t") && alleles.contains("chr1\t60\t"));
    // the matrix counts the pairs of the sites that paired
    let m = read(&prefix, "diff.discordance_matrix");
    assert_ne!(m, golden("main.diff.discordance_matrix"));
}

#[test]
fn compare_half_missing_genotypes_are_missing() {
    // vcftools counts ./1 as called (and the matrix stops with `Unhandled
    // case`, so there is no vcftools matrix); gloci reads any missing allele
    // as a missing genotype.
    let (_d, prefix) = compare(&g("half-a.vcf"), &g("half-b.vcf"), &["--match", "position"]);
    assert_eq!(
        read(&prefix, "diff.sites"),
        "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         chr1\t10\tB\t1\t0\t0\tnan\n\
         chr1\t20\tB\t1\t0\t0\tnan\n\
         chr1\t30\tB\t1\t0\t0\tnan\n"
    );
    assert_eq!(
        read(&prefix, "diff.indv"),
        "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         s1\t0\t0\tnan\ns2\t0\t0\tnan\ns3\t0\t0\tnan\ns4\t0\t0\tnan\n"
    );
    assert_eq!(
        read(&prefix, "diff.discordance_matrix"),
        "-\tN_0/0_file1\tN_0/1_file1\tN_1/1_file1\tN_./._file1\n\
         N_0/0_file2\t0\t0\t0\t2\n\
         N_0/1_file2\t0\t0\t0\t3\n\
         N_1/1_file2\t0\t0\t0\t1\n\
         N_./._file2\t1\t1\t0\t4\n"
    );
    let theirs = golden("half.diff.sites");
    assert!(theirs.contains("chr1\t10\tB\t1\t4\t3\t0.75"), "{theirs}");
    assert!(!golden_dir().join("half.diff.discordance_matrix").exists());
}

#[test]
fn compare_sample_union() {
    // The union by name in byte order; a sample in one file is 0, 0, nan; the
    // samples of both files count over the shared sites.
    let (_d, prefix) = compare(
        &g("diffsamples-a.vcf"),
        &g("diffsamples-b.vcf"),
        &["--sample-discordance"],
    );
    assert_eq!(
        read(&prefix, "diff.indv"),
        "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         alpha\t0\t0\tnan\nbeta\t3\t1\t0.333333\ndel2\t0\t0\tnan\n\
         delta\t0\t0\tnan\nomega\t0\t0\tnan\nzed\t3\t2\t0.666667\n"
    );
    // --sample-map renames B's samples: the pairs join, the names are A's
    let (_d, prefix) = compare(
        &g("diffsamples-a.vcf"),
        &g("diffsamples-b.vcf"),
        &["--sample-discordance", "--sample-map", &g("map.txt")],
    );
    assert_eq!(read(&prefix, "diff.indv"), golden("dsmap.diff.indv"));
    // --samples names samples of either file and applies to each that has them
    let (_d, prefix) = compare(
        &g("diffsamples-a.vcf"),
        &g("diffsamples-b.vcf"),
        &["--sample-discordance", "--samples", "zed,omega,del2"],
    );
    assert_eq!(
        read(&prefix, "diff.indv"),
        "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         del2\t0\t0\tnan\nomega\t0\t0\tnan\nzed\t3\t2\t0.666667\n"
    );
    // by a renamed sample's name in A, with the map
    let (_d, prefix) = compare(
        &g("diffsamples-a.vcf"),
        &g("diffsamples-b.vcf"),
        &[
            "--sample-discordance",
            "--sample-map",
            &g("map.txt"),
            "--exclude-samples",
            "alpha,beta",
        ],
    );
    assert_eq!(
        read(&prefix, "diff.indv"),
        "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         delta\t3\t0\t0\nzed\t3\t2\t0.666667\n"
    );
    // a name in neither file is a usage error
    let d = tempdir();
    let p = d.path().join("o").to_string_lossy().into_owned();
    gloci()
        .args([
            "vcf",
            "compare",
            &g("diffsamples-a.vcf"),
            &g("diffsamples-b.vcf"),
        ])
        .args(["--out-prefix", &p, "--samples", "zed,x"])
        .assert()
        .code(2)
        .stderr(contains("--samples names x"));
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
    // so is a sample-map that makes a name twice
    let map = write(d.path(), "dup.txt", "beta\tzed\n");
    gloci()
        .args([
            "vcf",
            "compare",
            &g("diffsamples-a.vcf"),
            &g("diffsamples-b.vcf"),
        ])
        .args(["--out-prefix", &p, "--sample-map", &map])
        .assert()
        .code(2)
        .stderr(contains("twice"));
}

#[test]
fn compare_reads_only_the_selected_samples() {
    // A has a sample whose cell is not a genotype; --samples leaves it out of
    // the selection, so nothing reads it
    let d = tempdir();
    let a = write(
        d.path(),
        "a.vcf",
        &vcf(
            &["s1", "only_a"],
            &[&format!(
                "{}\t0/1\tNOT-A-GENOTYPE",
                rec("chr1", 5, "A", "C", &[]).trim_end()
            )],
        ),
    );
    let b = write(
        d.path(),
        "b.vcf",
        &vcf(&["s1"], &[&rec("chr1", 5, "A", "C", &["0/1"])]),
    );
    let p = d.path().join("o").to_string_lossy().into_owned();
    gloci()
        .args(["vcf", "compare", &a, &b, "--out-prefix", &p])
        .assert()
        .code(1)
        .stderr(contains("NOT-A-GENOTYPE"));
    gloci()
        .args([
            "vcf",
            "compare",
            &a,
            &b,
            "--out-prefix",
            &p,
            "--samples",
            "s1",
        ])
        .assert()
        .success();
    assert_eq!(
        read(&p, "diff.sites"),
        "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         chr1\t5\tB\t1\t1\t0\t0\n"
    );
}

#[test]
fn compare_masked_genotypes_are_missing() {
    let d = tempdir();
    let row = |gt: &str, gq: &str| format!("chr1\t5\t.\tA\tC\t50\tPASS\t.\tGT:GQ\t{gt}:{gq}");
    let a = write(d.path(), "a.vcf", &vcf(&["s1"], &[&row("0/1", "5")]));
    let b = write(d.path(), "b.vcf", &vcf(&["s1"], &[&row("1/1", "50")]));
    let (_d, p) = compare(&a, &b, &["--site-discordance"]);
    assert!(read(&p, "diff.sites").ends_with("chr1\t5\tB\t1\t1\t1\t1\n"));
    let (_d, p) = compare(&a, &b, &["--site-discordance", "--min-gq", "20"]);
    assert!(read(&p, "diff.sites").ends_with("chr1\t5\tB\t1\t0\t0\tnan\n"));
}

#[test]
fn compare_sorted_input_is_checked() {
    let d = tempdir();
    let good = write(
        d.path(),
        "good.vcf",
        &vcf(
            &["s1"],
            &[
                &rec("chr1", 5, "A", "C", &["0/1"]),
                &rec("chr2", 5, "A", "C", &["0/1"]),
            ],
        ),
    );
    let bad = write(
        d.path(),
        "bad.vcf",
        &vcf(
            &["s1"],
            &[
                &rec("chr1", 9, "A", "C", &["0/1"]),
                &rec("chr1", 5, "A", "C", &["0/1"]),
            ],
        ),
    );
    let p = d.path().join("o").to_string_lossy().into_owned();
    for (a, b, name) in [(&good, &bad, "bad.vcf"), (&bad, &good, "bad.vcf")] {
        gloci()
            .args(["vcf", "compare", a, b, "--out-prefix", &p])
            .assert()
            .code(1)
            .stderr(contains(format!("{name} is not sorted")));
    }
    let backwards = write(
        d.path(),
        "back.vcf",
        &vcf(
            &["s1"],
            &[
                &rec("chr2", 5, "A", "C", &["0/1"]),
                &rec("chr1", 5, "A", "C", &["0/1"]),
            ],
        ),
    );
    gloci()
        .args(["vcf", "compare", &good, &backwards, "--out-prefix", &p])
        .assert()
        .code(1)
        .stderr(contains("back.vcf is not sorted"));
    // an error publishes nothing
    let names: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.starts_with("o."))
        .collect();
    assert!(names.is_empty(), "{names:?}");
}

#[test]
fn compare_chrom_tables() {
    // header order chr1, chr2, chr10 in both files (karyotype and lexical
    // order would put chr10 elsewhere)
    let d = tempdir();
    let header = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n##contig=<ID=chr10>\n\
                  ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
    let rows_a = [("chr1", 5), ("chr2", 5), ("chr10", 5)];
    let rows_b = [("chr2", 5), ("chr10", 5), ("chr10", 7)];
    let text = |rows: &[(&str, u64)]| {
        let mut s = header.to_string();
        for (c, p) in rows {
            s.push_str(&rec(c, *p, "A", "C", &["0/1"]));
            s.push('\n');
        }
        s
    };
    let a = write(d.path(), "a.vcf", &text(&rows_a));
    let b = write(d.path(), "b.vcf", &text(&rows_b));
    let (_o, p) = compare(&a, &b, &["--site-discordance"]);
    assert_eq!(
        read(&p, "diff.sites"),
        "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         chr1\t5\t1\t0\t0\t0\tnan\n\
         chr2\t5\tB\t1\t1\t0\t0\n\
         chr10\t5\tB\t1\t1\t0\t0\n\
         chr10\t7\t2\t0\t0\t0\tnan\n"
    );
    // -g names another order: chr10 before chr2 makes these files unsorted
    let genome = write(d.path(), "g.genome", "chr1\t100\nchr10\t100\nchr2\t100\n");
    let o = d.path().join("g").to_string_lossy().into_owned();
    gloci()
        .args(["vcf", "compare", &a, &b, "--out-prefix", &o, "-g", &genome])
        .assert()
        .code(1)
        .stderr(contains("is not sorted"));
    // --chrom-map renames both files' chromosomes before the order is read
    let map = write(d.path(), "map.tsv", "chr10\tchrX\n");
    let (_o, p) = compare(&a, &b, &["--site-discordance", "--chrom-map", &map]);
    assert!(read(&p, "diff.sites").contains("chrX\t5\tB"));
    // a chromosome the -g genome lacks is an error that names it
    let small = write(d.path(), "s.genome", "chr1\t100\nchr2\t100\n");
    gloci()
        .args(["vcf", "compare", &a, &b, "--out-prefix", &o, "-g", &small])
        .assert()
        .failure()
        .stderr(contains("chr10"));
}

#[test]
fn compare_huge_contig_count() {
    // 200,000 declared contigs, a record on each of a few
    let d = tempdir();
    let mut header = String::from("##fileformat=VCFv4.2\n");
    for i in 0..200_000 {
        header.push_str(&format!("##contig=<ID=c{i}>\n"));
    }
    header.push_str(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n",
    );
    let text = |chroms: &[u32]| {
        let mut s = header.clone();
        for c in chroms {
            s.push_str(&rec(&format!("c{c}"), 5, "A", "C", &["0/1"]));
            s.push('\n');
        }
        s
    };
    let a = write(d.path(), "a.vcf", &text(&[0, 3, 199_999]));
    let b = write(d.path(), "b.vcf", &text(&[3, 100_000, 199_999]));
    let (_o, p) = compare(&a, &b, &["--site-discordance"]);
    assert_eq!(
        read(&p, "diff.sites"),
        "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         c0\t5\t1\t0\t0\t0\tnan\n\
         c3\t5\tB\t1\t1\t0\t0\n\
         c100000\t5\t2\t0\t0\t0\tnan\n\
         c199999\t5\tB\t1\t1\t0\t0\n"
    );
}

/// A pseudo-random pair of files: `sites` records for `samples` samples.
fn random_pair(d: &Path, sites: u64, samples: usize) -> (String, String) {
    let mut x = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let names: Vec<String> = (0..samples).map(|i| format!("s{i:03}")).collect();
    let refs: Vec<&str> = names.iter().map(String::as_str).collect();
    let gt = |r: u64| match r % 9 {
        0..=2 => "0/0",
        3 | 4 => "0/1",
        5 => "1/1",
        6 => "./.",
        7 => "1",
        _ => "1|0",
    };
    let mut rows = [Vec::new(), Vec::new()];
    for i in 0..sites {
        let pos = 100 + i * 3;
        for (f, r) in rows.iter_mut().enumerate() {
            // each file skips some sites, and sometimes has another ALT
            if next() % 11 == 0 {
                continue;
            }
            let alt = if next() % 13 == 0 && f == 1 { "G" } else { "C" };
            let gts: Vec<&str> = (0..samples).map(|_| gt(next())).collect();
            r.push(rec("chr1", pos, "A", alt, &gts));
        }
    }
    let [ra, rb] = rows;
    let a = write(
        d,
        "ra.vcf",
        &vcf(&refs, &ra.iter().map(String::as_str).collect::<Vec<_>>()),
    );
    let b = write(
        d,
        "rb.vcf",
        &vcf(&refs, &rb.iter().map(String::as_str).collect::<Vec<_>>()),
    );
    (a, b)
}

#[test]
fn compare_output_is_the_same_at_every_thread_count() {
    let d = tempdir();
    // 3,000 sites of 40 samples: several batches per thread
    let (a, b) = random_pair(d.path(), 3000, 40);
    for rule in ["alleles", "position"] {
        let mut outs: Vec<Vec<String>> = Vec::new();
        for t in ["1", "2", "8"] {
            let (_o, p) = compare(&a, &b, &["--match", rule, "-t", t]);
            outs.push(EXTS.iter().map(|e| read(&p, e)).collect());
        }
        assert_eq!(outs[0], outs[1], "{rule} -t 2");
        assert_eq!(outs[0], outs[2], "{rule} -t 8");
        assert!(outs[0][0].lines().count() > 2500);
    }
}

#[test]
fn compare_reads_bcf_and_gzip() {
    let d = tempdir();
    let a = g("a.vcf");
    let b = g("b.vcf");
    let abcf = d.path().join("a.bcf").to_string_lossy().into_owned();
    gloci()
        .args(["view", &a, "-O", "bcf", "-o", &abcf])
        .assert()
        .success();
    let bgz = bgzf_file(d.path(), "b.vcf.gz", &std::fs::read_to_string(&b).unwrap());
    let (_o, p) = compare(&abcf, &bgz.to_string_lossy(), &["--match", "position"]);
    for ext in EXTS {
        assert_eq!(read(&p, ext), golden(&format!("main.{ext}")), "{ext}");
    }
}

#[test]
fn compare_regions_and_site_filters_apply_to_both_files() {
    let (_o, p) = compare(
        &g("a.vcf"),
        &g("b.vcf"),
        &[
            "--match",
            "position",
            "--site-discordance",
            "--chrom",
            "chr2",
        ],
    );
    assert_eq!(
        read(&p, "diff.sites"),
        "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\n\
         chr2\t5\tB\t1\t3\t1\t0.333333\n\
         chr2\t15\t1\t0\t0\t0\tnan\n\
         chr2\t25\t2\t0\t0\t0\tnan\n"
    );
}

#[test]
fn compare_usage_errors() {
    let d = tempdir();
    let p = d.path().join("o").to_string_lossy().into_owned();
    gloci()
        .args(["vcf", "compare", "-", "-", "--out-prefix", &p])
        .assert()
        .code(2)
        .stderr(contains("standard input"));
    let bed = write(d.path(), "x.bed", "chr1\t1\t2\n");
    gloci()
        .args(["vcf", "compare", &bed, &g("a.vcf"), "--out-prefix", &p])
        .assert()
        .code(2)
        .stderr(contains("vcf compare works on VCF and BCF"));
}

#[test]
fn compare_matrix_scope_follows_the_match_rule() {
    // AT>A against ATG>AG: one trimmed allele, so alleles mode counts it in
    // the matrix; position mode sees different REFs and does not.
    let d = tempdir();
    let a = write(
        d.path(),
        "a.vcf",
        &vcf(
            &["s1", "s2"],
            &[&rec("chr1", 60, "AT", "A", &["0/1", "1/1"])],
        ),
    );
    let b = write(
        d.path(),
        "b.vcf",
        &vcf(
            &["s1", "s2"],
            &[&rec("chr1", 60, "ATG", "AG", &["0/1", "1/1"])],
        ),
    );
    let head = "-\tN_0/0_file1\tN_0/1_file1\tN_1/1_file1\tN_./._file1\n";
    let (_o, p) = compare(&a, &b, &["--discordance-matrix"]);
    assert_eq!(
        read(&p, "diff.discordance_matrix"),
        format!(
            "{head}N_0/0_file2\t0\t0\t0\t0\nN_0/1_file2\t0\t1\t0\t0\nN_1/1_file2\t0\t0\t1\t0\nN_./._file2\t0\t0\t0\t0\n"
        )
    );
    let (_o, p) = compare(&a, &b, &["--discordance-matrix", "--match", "position"]);
    assert_eq!(
        read(&p, "diff.discordance_matrix"),
        format!(
            "{head}N_0/0_file2\t0\t0\t0\t0\nN_0/1_file2\t0\t0\t0\t0\nN_1/1_file2\t0\t0\t0\t0\nN_./._file2\t0\t0\t0\t0\n"
        )
    );
}

#[test]
fn compare_error_in_both_files_is_the_same_at_every_thread_count() {
    let d = tempdir();
    let bad = |pos: u64| format!("chr1\t{pos}\t.\tA\tC\t50\tPASS\t.\tGT\tNOPE");
    let ok = |pos: u64| rec("chr1", pos, "A", "C", &["0/1"]);
    let mut ra: Vec<String> = (1..3000).map(|p| ok(p * 2)).collect();
    ra.insert(2500, bad(7000));
    let mut rb: Vec<String> = (1..3000).map(|p| ok(p * 2)).collect();
    rb.insert(1200, bad(3000));
    let text = |rows: &[String]| {
        vcf(
            &["s1"],
            &rows.iter().map(String::as_str).collect::<Vec<_>>(),
        )
    };
    let a = write(d.path(), "a.vcf", &text(&ra));
    let b = write(d.path(), "b.vcf", &text(&rb));
    let mut msgs = Vec::new();
    for t in ["1", "2", "8"] {
        let p = d
            .path()
            .join(format!("o{t}"))
            .to_string_lossy()
            .into_owned();
        let out = gloci()
            .args(["vcf", "compare", &a, &b, "--out-prefix", &p, "-t", t])
            .assert()
            .code(1)
            .get_output()
            .stderr
            .clone();
        msgs.push(String::from_utf8(out).unwrap());
    }
    assert_eq!(msgs[0], msgs[1]);
    assert_eq!(msgs[0], msgs[2]);
    // B's error is met first in the walk (position 3000 before 7000)
    assert!(msgs[0].contains("b.vcf record 1,201"), "{}", msgs[0]);
}

#[test]
fn compare_order_error_names_the_previous_selected_record() {
    // record 2 is filtered out (not PASS); record 3 breaks the order after
    // record 1, which is the record the error names
    let d = tempdir();
    let filtered = "chr1\t50\t.\tA\tC\t50\tq10\t.\tGT\t0/1";
    let a = write(
        d.path(),
        "a.vcf",
        &vcf(
            &["s1"],
            &[
                &rec("chr1", 90, "A", "C", &["0/1"]),
                filtered,
                &rec("chr1", 5, "A", "C", &["0/1"]),
            ],
        ),
    );
    let b = write(
        d.path(),
        "b.vcf",
        &vcf(&["s1"], &[&rec("chr1", 1, "A", "C", &["0/1"])]),
    );
    let p = d.path().join("o").to_string_lossy().into_owned();
    gloci()
        .args(["vcf", "compare", &a, &b, "--out-prefix", &p, "--pass-only"])
        .assert()
        .code(1)
        .stderr(contains("record 3 (chr1:5) comes after record 1 (chr1:90)"));
}

/// `--thin` thins each file on its own, whatever the two headers list and
/// with `-g`: it checks only that each chromosome is one block.
#[test]
fn compare_thin_works_with_different_contig_lists_and_with_g() {
    let d = tempdir();
    let rows = |ps: &[u64]| -> Vec<String> {
        ps.iter()
            .map(|p| rec("chr1", *p, "A", "C", &["0/1"]))
            .collect()
    };
    let (ra, rb) = (rows(&[1, 5, 20]), rows(&[1, 12, 15, 20]));
    let ra: Vec<&str> = ra.iter().map(String::as_str).collect();
    let rb: Vec<&str> = rb.iter().map(String::as_str).collect();
    let a = write(d.path(), "a.vcf", &vcf(&["s1"], &ra));
    let b = write(
        d.path(),
        "b.vcf",
        &vcf(&["s1"], &rb).replace("##contig=<ID=chr2,length=1000>\n", ""),
    );
    let genome = write(d.path(), "g.genome", "chr1\t1000\nchr2\t1000\n");
    let sites = |extra: &[&str]| {
        let p = d.path().join("o").to_string_lossy().into_owned();
        gloci()
            .args([
                "vcf",
                "compare",
                &a,
                &b,
                "--out-prefix",
                &p,
                "--site-discordance",
            ])
            .args(["--thin", "10"])
            .args(extra)
            .assert()
            .success();
        read(&p, "diff.sites")
    };
    // A keeps 1 and 20; B keeps 1 and 12 (15 and 20 are under 10 from 12)
    for out in [sites(&[]), sites(&["-g", &genome])] {
        let rows: Vec<(&str, &str)> = out
            .lines()
            .skip(1)
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[1], f[2])
            })
            .collect();
        assert_eq!(rows, [("1", "B"), ("12", "2"), ("20", "1")], "{out}");
    }
}
