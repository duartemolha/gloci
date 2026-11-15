//! `gloci igv` and `gloci links` (D169, D170): the rules beyond the bedtools goldens
//! (`igv-*-vs-bedtools`, `links-*-vs-bedtools`).
mod common;
use common::*;
use predicates::str::contains;

fn run(args: &[&str]) -> (i32, String, String) {
    let o = gloci().args(args).output().unwrap();
    (
        o.status.code().unwrap(),
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
    )
}

fn file(dir: &tempfile::TempDir, name: &str, text: &str) -> String {
    plain_file(dir.path(), name, text)
        .to_string_lossy()
        .into_owned()
}

#[test]
fn igv_script_lines_in_input_order() {
    let d = tempdir();
    // Unsorted on purpose: the script follows the file.
    let f = file(&d, "a.bed", "chr2\t10\t20\tb\nchr1\t100\t200\tfoo/bar x\n");
    let (code, out, _) = run(&["igv", &f]);
    assert_eq!(code, 0);
    assert_eq!(
        out,
        "snapshotDirectory ./\ngoto chr2:10-20\nsnapshot chr2_10_20.png\n\
goto chr1:100-200\nsnapshot chr1_100_200.png\n"
    );
    let (_, out, _) = run(&[
        "igv",
        &f,
        "--path",
        "/x",
        "--session",
        "s.xml",
        "--sort",
        "readGroup",
        "--collapse",
        "--name",
        "--slop",
        "30",
        "--image",
        "eps",
    ]);
    assert_eq!(
        out,
        "snapshotDirectory /x\nload s.xml\n\
goto chr2:-20-50\nsort readGroup\ncollapse\nsnapshot chr2_10_20_b_slop30.eps\n\
goto chr1:70-230\nsort readGroup\ncollapse\nsnapshot chr1_100_200_foo/bar x_slop30.eps\n"
    );
}

/// A zero-length BED record is shown one base wider on each side, GFF is not, and a VCF record
/// is its REF bases whatever the ALT or END says.
#[test]
fn igv_positions_follow_the_format() {
    let d = tempdir();
    let bed = file(&d, "p.bed", "chr1\t0\t0\ta\nchr1\t5\t5\tb\nchr1\t7\t8\tc\n");
    let (_, out, _) = run(&["igv", &bed, "--name"]);
    assert_eq!(
        out,
        "snapshotDirectory ./\ngoto chr1:-1-1\nsnapshot chr1_-1_1_a.png\n\
goto chr1:4-6\nsnapshot chr1_4_6_b.png\ngoto chr1:7-8\nsnapshot chr1_7_8_c.png\n"
    );
    let gff = file(
        &d,
        "p.gff",
        "chr1\ts\tgene\t101\t100\t.\t+\t.\tID=x\nchr1\ts\texon\t1\t50\t.\t-\t.\tID=y\n",
    );
    let (_, out, _) = run(&["igv", &gff, "--name"]);
    assert_eq!(
        out,
        "snapshotDirectory ./\ngoto chr1:100-100\nsnapshot chr1_100_100_gene.png\n\
goto chr1:0-50\nsnapshot chr1_0_50_exon.png\n"
    );
    let vcf = file(
        &d,
        "p.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t200\t.\tAT\tA\t.\t.\t.\nchr1\t300\trs3\tG\t<DEL>\t.\t.\tEND=500\n",
    );
    let (_, out, _) = run(&["igv", &vcf, "--name"]);
    assert_eq!(
        out,
        "snapshotDirectory ./\ngoto chr1:199-201\nsnapshot chr1_199_201_AT/A.png\n\
goto chr1:299-300\nsnapshot chr1_299_300_G/<DEL>_rs3.png\n"
    );
}

#[test]
fn igv_name_needs_a_name() {
    let d = tempdir();
    let f = file(&d, "a.bed", "chr1\t1\t2\nchr1\t3\t4\n");
    let (code, out, err) = run(&["igv", &f, "--name"]);
    assert_eq!(code, 1);
    // What came before the first bad record stays, as with bedtools.
    assert_eq!(out, "snapshotDirectory ./\n");
    assert!(err.contains("line 1") && err.contains("--name"), "{err}");
    // A trailing tab does not make a fourth column.
    let f = file(&d, "t.bed", "chr1\t1\t2\t\n");
    assert_eq!(run(&["igv", &f, "--name"]).0, 1);
    assert_eq!(run(&["igv", &f]).0, 0);
}

#[test]
fn records_must_have_the_same_number_of_columns() {
    let d = tempdir();
    let f = file(&d, "r.bed", "#c\n\nchr1\t1\t2\ta\tb\nchr1\t3\t4\tc\n");
    for cmd in ["igv", "links"] {
        let (code, _, err) = run(&[cmd, &f]);
        assert_eq!(code, 1, "{cmd}");
        assert!(
            err.contains("line 4") && err.contains("4 columns") && err.contains("has 5"),
            "{cmd}: {err}"
        );
    }
    // A VCF is not checked.
    let v = file(
        &d,
        "r.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\tG\t.\t.\t.\nchr1\t2\t.\tA\tG\t.\t.\t.\tGT\n",
    );
    assert_eq!(run(&["igv", &v]).0, 0);
}

#[test]
fn unrecognised_input_writes_nothing() {
    let d = tempdir();
    let f = file(&d, "s.bed", "chr1 100 200 n\n");
    for cmd in ["igv", "links"] {
        let (code, out, err) = run(&[cmd, &f]);
        assert_eq!((code, out.as_str()), (1, ""), "{cmd}");
        assert!(err.contains("could not recognise"), "{err}");
    }
}

#[test]
fn igv_option_values_are_checked() {
    let d = tempdir();
    let f = file(&d, "a.bed", "chr1\t1\t2\ta\n");
    for args in [
        ["--image", "jpg"],
        ["--slop", "abc"],
        ["--slop", "-5"],
        ["--slop", "2147483648"],
        ["--sort", "Base"],
        ["--sort", ""],
    ] {
        let (code, out, _) = run(&["igv", &f, args[0], args[1]]);
        assert_eq!((code, out.as_str()), (2, ""), "{args:?}");
    }
    let (code, out, _) = run(&["igv", &f, "--slop", "2147483647"]);
    assert_eq!(code, 0);
    assert!(out.contains("goto chr1:-2147483646-2147483649\n"), "{out}");
}

#[test]
fn links_page_cells_by_column_count() {
    let d = tempdir();
    let cells = |text: &str| {
        let f = file(&d, "c.bed", text);
        let (code, out, _) = run(&["links", &f]);
        assert_eq!(code, 0, "{text}");
        out.matches("\t<td>\n").count() - 1 // the link's cell
    };
    assert_eq!(cells("chr1\t1\t2\n"), 0);
    // One tab at the end of a line is not a column.
    assert_eq!(cells("chr1\t1\t2\t\n"), 0);
    assert_eq!(cells("chr1\t1\t2\ta\t\n"), 1);
    assert_eq!(cells("chr1\t1\t2\ta\t\t\n"), 2);
    assert_eq!(cells("chr1\t1\t2\ta\n"), 1);
    assert_eq!(cells("chr1\t1\t2\ta\t5\n"), 2);
    assert_eq!(cells("chr1\t1\t2\ta\t5\t+\n"), 3);
    assert_eq!(cells("chr1\t1\t2\ta\t5\t+\tx\n"), 0);
    assert_eq!(cells("chr1\t1\t2\ta\t5\t+\tx\ty\n"), 0);
    assert_eq!(cells("chr1\t1\t2\ta\t5\t+\tx\ty\tz\n"), 3);
}

#[test]
fn links_options_and_title() {
    let d = tempdir();
    let f = file(&d, "a.bed", "chr1\t100\t200\n");
    let (_, out, _) = run(&[
        "links",
        &f,
        "--base",
        "http://m.x",
        "--org",
        "mouse",
        "--db",
        "mm9",
    ]);
    assert!(out.contains(&format!("<title>{f}</title>\n")), "{out}");
    assert!(
        out.contains(
            "\t\t<a href=http://m.x/cgi-bin/hgTracks?org=mouse&db=mm9&position=chr1:100-200>\
chr1:101-200</a>\n"
        ),
        "{out}"
    );
    let (_, out, _) = run(&["links", &f]);
    assert!(out.contains("href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg18&"));
    // Standard input is titled stdin.
    let o = gloci()
        .args(["links", "-"])
        .write_stdin("chr1\t1\t2\n")
        .output()
        .unwrap();
    assert!(
        String::from_utf8(o.stdout)
            .unwrap()
            .contains("<title>stdin</title>\n")
    );
}

#[test]
fn output_is_the_same_at_every_thread_count() {
    let d = tempdir();
    let mut text = String::new();
    for i in 0..20_000 {
        text.push_str(&format!(
            "chr{}\t{}\t{}\tn{i}\t{}\t+\n",
            i % 5 + 1,
            i * 7,
            i * 7 + 5,
            i % 9
        ));
    }
    let f = file(&d, "big.bed", &text);
    let gz = bgzf_file(d.path(), "big.bed.gz", &text)
        .to_string_lossy()
        .into_owned();
    for (cmd, extra) in [("igv", vec!["--name", "--slop", "5"]), ("links", vec![])] {
        for input in [&f, &gz] {
            let mut outs = Vec::new();
            for t in ["1", "2", "8"] {
                let mut a = vec![cmd, input.as_str(), "-t", t];
                a.extend(&extra);
                let (code, out, _) = run(&a);
                assert_eq!(code, 0);
                outs.push(out);
            }
            assert!(outs[0].len() > 100_000, "{cmd}");
            assert!(outs.iter().all(|o| *o == outs[0]), "{cmd} {input}");
        }
    }
}

#[test]
fn output_file_is_written_only_on_success() {
    let d = tempdir();
    let ok = file(&d, "a.bed", "chr1\t1\t2\ta\n");
    let out = d.path().join("s.txt");
    gloci()
        .args(["igv", &ok, "-o", out.to_str().unwrap()])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&out)
            .unwrap()
            .ends_with("snapshot chr1_1_2.png\n")
    );
    let bad = file(&d, "b.bed", "chr1\t1\t2\ta\nchr1\t1\t2\n");
    let out2 = d.path().join("t.html");
    gloci()
        .args(["links", &bad, "-o", out2.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(contains("line 2"));
    assert!(!out2.exists());
}

#[test]
fn other_formats_are_refused() {
    let d = tempdir();
    let f = file(&d, "a.bed", "chr1\t1\t2\ta\n");
    gloci()
        .args(["links", &f, "--format", "wig"])
        .assert()
        .failure();
}
