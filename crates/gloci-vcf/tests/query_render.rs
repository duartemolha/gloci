//! `vcf query` renders records and the -H line as bcftools 1.21 does (D172).

use gloci_io::VcfHeader;
use gloci_vcf::{QueryFormat, SampleCols};

const INPUTS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../xtask/diff/inputs");

fn read(path: &str) -> String {
    if path == "main.vcf" {
        let p = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/expr/main.vcf");
        return std::fs::read_to_string(p).unwrap();
    }
    std::fs::read_to_string(format!("{INPUTS}/{path}")).unwrap()
}

fn split(text: &str) -> (VcfHeader, Vec<&str>) {
    let end: usize = text
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(|l| l.len() + 1)
        .sum();
    let header = VcfHeader::parse(text.as_bytes()[..end].to_vec()).unwrap();
    (header, text[end..].lines().collect())
}

fn render_text(text: &str, fmt: &str) -> String {
    let (header, lines) = split(text);
    let q = QueryFormat::parse(fmt)
        .unwrap()
        .resolve(&header, false)
        .unwrap();
    let kept: Vec<usize> = (0..header.samples.len()).collect();
    let samples = SampleCols {
        kept: &kept,
        order: &kept,
        names: &header.samples,
        pass: None,
    };
    let mut out = Vec::new();
    for l in lines {
        q.write_record(l.as_bytes(), &samples, &mut out).unwrap();
    }
    String::from_utf8(out).unwrap()
}

fn run(path: &str, fmt: &str) -> String {
    render_text(&read(path), fmt)
}

fn run_err(text: &str, fmt: &str) -> String {
    let (header, lines) = split(text);
    let q = QueryFormat::parse(fmt)
        .unwrap()
        .resolve(&header, false)
        .unwrap();
    let kept: Vec<usize> = (0..header.samples.len()).collect();
    let samples = SampleCols {
        kept: &kept,
        order: &kept,
        names: &header.samples,
        pass: None,
    };
    let mut out = Vec::new();
    q.write_record(lines[0].as_bytes(), &samples, &mut out)
        .unwrap_err()
}

fn head(path: &str, fmt: &str, reads: bool) -> String {
    let (header, _) = split(&read(path));
    let q = QueryFormat::parse(fmt)
        .unwrap()
        .resolve(&header, true)
        .unwrap();
    let names: Vec<&[u8]> = header.samples.iter().map(|n| n.as_slice()).collect();
    let mut out = Vec::new();
    q.write_header(&names, reads, &mut out);
    String::from_utf8(out).unwrap()
}

#[test]
fn fixed_fields_match_bcftools() {
    let want = "chr1\t100\t99\t100\t99\trs1\tC\tA,G\tA\t50\tPASS\tSNP\n\
chr1\t200\t199\t202\t201\t.\tCAA\tCA,C,CAAA\tCA\t.\tq10\tINDEL\n\
chr1\t300\t299\t300\t299\t.\tA\tG,*\tG\t.\t.\tSNP,OVERLAP\n\
chr1\t400\t399\t400\t399\t.\tA\t<DEL>,T\t<DEL>\t.\t.\tSNP,OTHER\n\
chr1\t500\t499\t501\t500\t.\tAC\tGT,AT,GC\tGT\t.\t.\tSNP,MNP\n\
chr1\t600\t599\t600\t599\t.\tA\tT\tT\t.\t.\tSNP\n\
chr1\t700\t699\t700\t699\t.\tA\t.\t.\t.\t.\tREF\n\
chr1\t800\t799\t801\t800\t.\tCT\tCTA,C,CA\tCTA\t.\t.\tSNP,INDEL\n\
chr1\t900\t899\t902\t901\t.\tCAT\tCT,CAG,C\tCT\t.\t.\tSNP,INDEL\n\
chr1\t1000\t999\t1000\t999\t.\tG\tA,C,T\tA\t99\tPASS\tSNP\n\
chr1\t1100\t1099\t1100\t1099\t.\tA\tC,G\tC\t.\t.\tSNP\n";
    let got = run(
        "vcf/split/multi.vcf",
        r"%CHROM\t%POS\t%POS0\t%END\t%END0\t%ID\t%REF\t%ALT\t%FIRST_ALT\t%QUAL\t%FILTER\t%TYPE\n",
    );
    assert_eq!(got, want);
    assert_eq!(
        run(
            "vcf/query/types.vcf",
            r"%POS\t%TYPE\t%FIRST_ALT\t%ALT{1}\t%END\n"
        ),
        "5\tBND\tA[1:10[\t.\t5\n6\tSNP,BND\t]1:10]A\tC\t6\n7\tREF\t<NON_REF>\t.\t7\n8\tSNP,MNP,INDEL\tA\tACC\t9\n9\tREF\tX\t.\t9\n"
    );
}

#[test]
fn end_uses_info_end_or_ref_length() {
    let text = "##fileformat=VCFv4.2\n##contig=<ID=1>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t40\t.\tA\t<DEL>\t.\t.\tEND=50\n1\t60\t.\tAT\tA\t.\t.\tEND=10\n1\t70\t.\tA\t<DEL>\t.\t.\tEND=70\n";
    assert_eq!(
        render_text(text, r"%POS\t%END\t%END0\n"),
        "40\t50\t49\n60\t61\t60\n70\t70\t69\n"
    );
}

#[test]
fn info_tags_flags_and_subscripts() {
    let m = "vcf/split/multi.vcf";
    let got = run(m, r"%INFO/DP\t%DB\t%AC\t%AF\t%RC\t%GC\t%XS\n");
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines[0], "10\t1\t3,4\t0.25,0.5\t5,6,7\t1,2,3,4,5,6\ta,b");
    assert_eq!(lines[1], ".\t.\t1,2,3\t.\t.\t.\t.");
    assert_eq!(lines[5], ".\t.\t1\t.\t.\t.\t.");
    assert_eq!(
        lines[9],
        "20\t.\t1,2,3\t.\t1,2,3,4\t1,2,3,4,5,6,7,8,9,10\t."
    );

    let got = run(
        m,
        r"%RC{0}\t%RC{3}\t%XS{1}\t%XS{5}\t%DB{0}\t%DB{1}\t%GC{9}\t%AF{1}\t%ALT{1}\t%ALT{5}\n",
    );
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines[0], "5\t.\tb\t.\t1\t1\t.\t0.5\tG\t.");
    assert_eq!(lines[2], ".\t.\t.\t.\t.\t.\t.\t.\t*\t.");
    assert_eq!(lines[9], "1\t4\t.\t.\t.\t.\t10\t.\tC\t.");

    // gloci's rule: a subscript picks the comma-separated element of the value.
    let got = run(m, r"%POS\t%DP{1}\t%AC{1}\n");
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines[0], "100\t.\t4");
    assert_eq!(lines[5], "600\t.\t.");
}

#[test]
fn numbers_print_as_bcftools_prints_them() {
    assert_eq!(
        run(
            "vcf/query/numbers.vcf",
            r"%POS\t%QUAL\t%I\t%F[\t%GL:%FI:%FF]\n"
        ),
        "5\t3.14159\t7\t2.5,-0,nan\t-0.1,-2,-30.5:12:0.1\n6\t0\t-5\t.\t.,.,.:.:.\n"
    );
    assert_eq!(
        run(
            "filter-float-text.vcf",
            r"%CHROM\t%POS\t%QUAL\t%AF\t%AF_afr\n"
        ),
        "chr22\t10736362\t.\t0.194444\t0\nchr22\t10736399\t50\t0.490022\t0.6\nchr22\t10736400\t.\t0\t0\nchr22\t10736414\t.\t0.000719424\t0.015\n"
    );
}

#[test]
fn format_vectors_and_subscripts() {
    let m = "vcf/split/multi.vcf";
    let got = run(m, r"%POS[\t%SAMPLE=%GT:%AD:%PL:%GL:%AQ:%DP]\n");
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(
        lines[1],
        "200\ts1=0/1:1,2,3,4:.:.:.:.\ts2=2/3:.:.:.:.:.\ts3=1/3:1,2,3,4:.:.:.:.\ts4=./1:.:.:.:.:.\ts5=0/0:1,2,3,4:.:.:.:.\ts6=1:9,9,9,9:.:.:.:."
    );
    assert!(lines[9].ends_with("\ts5=.:.:.:.:.:.\ts6=1/1:.,.,.,.:.,.,.,.,.,.,.,.,.,.:.:.:3"));

    let got = run(m, r"%POS[\t%AD{0}:%AD{1}:%PL{2}]\n");
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines[0], "100\t5:6:3\t1:2:13\t1:2:23\t.:.:.\t5:6:.\t1:2:33");
    assert_eq!(lines[1], "200\t1:2:.\t.:.:.\t1:2:.\t.:.:.\t1:2:.\t9:9:.");
}

#[test]
fn gt_and_tgt() {
    let got = run("vcf/split/multi.vcf", r"%POS[\t%TGT]\n");
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines[0], "100\tC/A\tA/G\tG|G\t./.\tA\tC|G");
    assert_eq!(lines[1], "200\tCAA/CA\tC/CAAA\tCA/CAAA\t./CA\tCAA/CAA\tCA");
    assert_eq!(lines[2], "300\tG/*\tA/G\t*/*\t./.\tA\t*");
    assert_eq!(lines[3], "400\t<DEL>/T\tA/<DEL>\tT/T\t./.\tA\tT");
    assert_eq!(
        run("mixed-ploidy-missing.vcf", r"%POS[\t%GT\t%TGT]\n"),
        "10\t0/1\tA/G\t.\t.\n20\t0/1\tA/G\t1/1\tG/G\n30\t1\tG\t0|1\tA|G\n40\t.\t.\t0/1/2\tA/G/C\n50\t./.\t./.\t.\t.\n"
    );
    assert_eq!(
        run("extract-haploid.vcf", r"%POS[\t%TGT]\n"),
        "100\tA\tC\tA/C\tC/C\n200\tG\tG\tG/.\t./T\n300\t.\tG\t./.\tC/C\n"
    );
    let text = "##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n1\t5\t.\tA\tC\t.\t.\t.\tGT\t0/5\n";
    let e = run_err(text, r"[%TGT]\n");
    assert!(e.contains("allele 5") && e.contains("s1"), "{e}");
}

#[test]
fn site_fields_repeat_inside_a_block() {
    assert_eq!(
        run("expr-format.vcf", r"%POS[\t%SAMPLE:%QUAL:%DP:%INFO/DP]\n"),
        "10\ts1:50:25:20\ts2:50:10:20\ts3:50:30:20\n20\ts1:.:.:.\ts2:.:5:.\ts3:.:.:.\n30\ts1:10:.:3\ts2:10:.:3\ts3:10:.:3\n40\ts1:5:.:.\ts2:5:.:.\ts3:5:.:.\n"
    );
}

#[test]
fn escapes_and_the_automatic_newline() {
    let got = run("expr-format.vcf", r"a\\b\t%POS\n");
    assert_eq!(got.lines().next().unwrap(), "a\\b\t10");
    assert_eq!(run("expr-format.vcf", "%POS"), "10\n20\n30\n40\n");
    assert_eq!(
        run("expr-format.vcf", "[%SAMPLE]"),
        "s1s2s3\ns1s2s3\ns1s2s3\ns1s2s3\n"
    );
    assert_eq!(
        run("expr-format.vcf", r"%POS\n%ID"),
        "10\nrs120\n.30\nx140\n."
    );
}

fn first_data_line(text: &str) -> String {
    text.lines().find(|l| !l.starts_with('#')).unwrap().into()
}

#[test]
fn samples_follow_kept_and_pass() {
    let text = read("main.vcf");
    let (header, _) = split(&text);
    let full = first_data_line(&text);
    let mut cols: Vec<&str> = full.split('\t').collect();
    cols.remove(9);
    let line = cols.join("\t");
    let q = QueryFormat::parse(r"%POS[\t%SAMPLE=%GT]\n")
        .unwrap()
        .resolve(&header, false)
        .unwrap();
    let render = |pass: Option<&[bool]>| {
        let mut out = Vec::new();
        let samples = SampleCols {
            kept: &[1, 2],
            order: &[0, 1],
            names: &header.samples,
            pass,
        };
        q.write_record(line.as_bytes(), &samples, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    };
    assert_eq!(render(Some(&[true, false, true])), "10\ts3=1/1\n");
    assert_eq!(render(None), "10\ts2=0/0\ts3=1/1\n");
}

#[test]
fn order_decides_the_printing_and_naming_order() {
    let text = read("main.vcf");
    let (header, _) = split(&text);
    let line = first_data_line(&text);
    let q = QueryFormat::parse(r"%POS[\t%SAMPLE=%GT]\n")
        .unwrap()
        .resolve(&header, true)
        .unwrap();
    let samples = SampleCols {
        kept: &[0, 1, 2],
        order: &[2, 0],
        names: &header.samples,
        pass: None,
    };
    let mut out = Vec::new();
    q.write_record(line.as_bytes(), &samples, &mut out).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), "10\ts3=1/1\ts1=0/1\n");
    // The -H line names the samples in the order it is given them.
    let names: Vec<&[u8]> = vec![b"s3", b"s1"];
    let mut out = Vec::new();
    q.write_header(&names, false, &mut out);
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "#[1]POS\t[2]s3:SAMPLE=[3]s3:GT\t[4]s1:SAMPLE=[5]s1:GT\n"
    );
}

#[test]
fn info_and_line_are_written_as_read() {
    assert_eq!(
        run("vcf/query/numbers.vcf", r"%INFO\n"),
        "I=007;F=2.50,-0.0,nan\nI=-5;F=.\n"
    );
    let text = read("vcf/query/numbers.vcf");
    let want: String = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(run("vcf/query/numbers.vcf", "%LINE"), want);
}

#[test]
fn bad_typed_value_is_an_error() {
    let text = "##fileformat=VCFv4.2\n##contig=<ID=1>\n##INFO=<ID=I,Number=1,Type=Integer,Description=\"x\">\n##FORMAT=<ID=FI,Number=1,Type=Integer,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n1\t5\t.\tA\tC\t.\t.\tI=abc\tFI\tx\n";
    let e = run_err(text, r"%I\n");
    assert!(e.contains("INFO/I") && e.contains("abc"), "{e}");
    let e = run_err(text, r"[%FI]\n");
    assert!(
        e.contains("FORMAT/FI") && e.contains('x') && e.contains("s1"),
        "{e}"
    );
}

#[test]
fn header_line_matches_bcftools() {
    let f = "expr-format.vcf";
    assert_eq!(
        head(f, r"%CHROM %POS %INFO/DP %DB{0} [ %GT %DP{0}]\n", false),
        "#[1]CHROM [2]POS [3]DP [4]DB  [5]s1:GT [6]s1:DP [7]s2:GT [8]s2:DP [9]s3:GT [10]s3:DP\n"
    );
    assert_eq!(
        head(f, r"%CHROM %POS %INFO/DP %DB{0} [ %GT %DP{0}]\n", true),
        "#[1]CHROM [2]POS [3]DP [4]DB  [5]GT [6]DP [7]GT [8]DP [9]GT [10]DP\n"
    );
    assert_eq!(
        head(f, r"%POS\t[%SAMPLE\t%GT\n]", false),
        "#[1]POS\t[2]SAMPLE\t[3]GT\n"
    );
    assert_eq!(
        head(f, r"[%SAMPLE\t%GT]", false),
        "#[1]s1:SAMPLE\t[2]s1:GT[3]s2:SAMPLE\t[4]s2:GT[5]s3:SAMPLE\t[6]s3:GT\n"
    );
    assert_eq!(
        head(f, r"[%SAMPLE\t%GT\n]x\n", false),
        "#[1]SAMPLE\t[2]GTx\n\n"
    );
    assert_eq!(
        head(f, r"[%SAMPLE\n]x[\t%GT]\n", false),
        "#[1]SAMPLEx\t[2]GT\n\n"
    );
    assert_eq!(
        head(f, r"%POS[\t%SAMPLE]x[\t%GT\n]\n", false),
        "#[1]POS\t[2]SAMPLE\t[3]SAMPLE\t[4]SAMPLEx\t[5]GT\n\n"
    );
    assert_eq!(
        head(f, r"%ALT{1}\t%INFO/AF\t%AF{0}[\t%TGT]\n", false),
        "#[1]ALT\t[2]AF\t[3]AF\t[4]s1:GT\t[5]s2:GT\t[6]s3:GT\n"
    );
    assert_eq!(head(f, r"%POS\t%INFO\n", false), "#[1]POS\t[2]INFO\n");
}

#[test]
fn block_with_no_samples() {
    let text = "##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tC\t.\t.\t.\n";
    assert_eq!(render_text(text, r"%POS[\t%GT]\n"), "100\n");
    let (header, _) = split(text);
    let q = QueryFormat::parse(r"%POS[\t%GT]\n")
        .unwrap()
        .resolve(&header, true)
        .unwrap();
    let mut out = Vec::new();
    q.write_header(&[], false, &mut out);
    assert_eq!(String::from_utf8(out).unwrap(), "#[1]POS\n");
}
