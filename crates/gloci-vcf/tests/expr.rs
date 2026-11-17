//! Filter expressions (D71): parser, evaluation tables, bcftools parity from
//! committed golden results, BCF/text agreement, and the pinned differences.
mod common;
use gloci_core::Interval;
use gloci_io::VcfHeader;
use gloci_vcf::{
    BcfSite, EvalScratch, Expr, GenotypeSource, MaskRules, MaskedGenotypes, SiteStats, SiteView,
    bcf_genotypes, parse_expr, site_stats, text_genotypes,
};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

struct Counting;
thread_local! { static ALLOCS: Cell<usize> = const { Cell::new(0) }; }
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let _ = ALLOCS.try_with(|c| c.set(c.get() + 1));
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        let _ = ALLOCS.try_with(|c| c.set(c.get() + 1));
        unsafe { System.realloc(p, l, n) }
    }
}
#[global_allocator]
static A: Counting = Counting;

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/expr");

fn header_of(vcf: &str) -> VcfHeader {
    let end = vcf
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(|l| l.len() + 1)
        .sum::<usize>();
    VcfHeader::parse(vcf.as_bytes()[..end].to_vec()).unwrap()
}

fn records(vcf: &str) -> Vec<&str> {
    vcf.lines().filter(|l| !l.starts_with('#')).collect()
}

fn span_of(line: &str) -> Interval {
    let f: Vec<&str> = line.split('\t').collect();
    let pos: u64 = f[1].parse().unwrap();
    Interval {
        start: pos - 1,
        end: pos - 1 + f[3].len() as u64,
    }
}

fn n_alleles(alts: &[u8]) -> usize {
    if alts == b"." {
        1
    } else {
        2 + alts.iter().filter(|&&b| b == b',').count()
    }
}

/// POS values of the text records for which `e` holds.
fn run_text(e: &Expr, vcf: &str) -> Vec<u64> {
    let mut scratch = EvalScratch::default();
    let mut masked = Vec::new();
    let mut out = Vec::new();
    for line in records(vcf) {
        let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
        let g = text_genotypes(line.as_bytes()).unwrap();
        let stats = site_stats(
            &g,
            None,
            &MaskRules::default(),
            n_alleles(site.alts),
            &mut masked,
        )
        .unwrap();
        if e.eval(&site, Some(&g), Some(&stats), &mut scratch).unwrap() {
            out.push(site.pos);
        }
    }
    out
}

/// The framed BCF records of `vcf`'s data lines (native encoder, so mixed
/// ploidy with a missing haploid GT round-trips) and the dictionaries.
fn to_bcf(vcf: &str) -> (gloci_io::BcfDicts, Vec<Vec<u8>>) {
    let header = header_of(vcf);
    let dicts = gloci_io::BcfDicts::from_header(&header.raw).unwrap();
    let mut enc = gloci_io::TextToBcf::new("t.vcf", &header).unwrap();
    let raws = records(vcf)
        .iter()
        .enumerate()
        .map(|(i, l)| {
            let mut rec = Vec::new();
            enc.encode(l.as_bytes(), i as u64 + 1, &mut rec).unwrap();
            rec
        })
        .collect();
    (dicts, raws)
}

fn run_bcf(e: &Expr, vcf: &str) -> Vec<u64> {
    let (dicts, raws) = to_bcf(vcf);
    let mut scratch = EvalScratch::default();
    let mut masked = Vec::new();
    let mut bs = BcfSite::with_info_values();
    let mut out = Vec::new();
    for (raw, line) in raws.iter().zip(records(vcf)) {
        let g = bcf_genotypes(raw, &dicts).unwrap();
        let site = bs.view(raw, &dicts, span_of(line)).unwrap();
        let stats = site_stats(
            &g,
            None,
            &MaskRules::default(),
            n_alleles(site.alts),
            &mut masked,
        )
        .unwrap();
        if e.eval(&site, Some(&g), Some(&stats), &mut scratch).unwrap() {
            out.push(site.pos);
        }
    }
    out
}

fn fixture(name: &str) -> String {
    std::fs::read_to_string(format!("{FIXTURES}/{name}.vcf")).unwrap()
}

fn pos_list(s: &str) -> Vec<u64> {
    s.split(',')
        .filter(|p| !p.is_empty())
        .map(|p| p.parse().unwrap())
        .collect()
}

#[test]
fn golden_results_match_bcftools_on_text_and_bcf() {
    let golden = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/expr_bcftools.tsv"
    ))
    .unwrap();
    let mut n = 0;
    for row in golden.lines() {
        let [fx, expr, want] = row.split('\t').collect::<Vec<_>>()[..] else {
            panic!("bad golden row {row:?}");
        };
        let vcf = fixture(fx);
        let e = parse_expr(expr, &header_of(&vcf)).unwrap_or_else(|e| panic!("{fx} {expr}: {e}"));
        assert_eq!(run_text(&e, &vcf), pos_list(want), "text: {fx} {expr}");
        assert_eq!(run_bcf(&e, &vcf), pos_list(want), "bcf: {fx} {expr}");
        n += 1;
    }
    assert!(n > 90, "{n} golden rows");
}

fn main_vcf() -> (String, VcfHeader) {
    let v = fixture("main");
    let h = header_of(&v);
    (v, h)
}

fn hits(expr: &str) -> Vec<u64> {
    let (v, h) = main_vcf();
    let e = parse_expr(expr, &h).unwrap_or_else(|e| panic!("{expr}: {e}"));
    let t = run_text(&e, &v);
    assert_eq!(run_bcf(&e, &v), t, "bcf vs text: {expr}");
    t
}

fn err(expr: &str) -> String {
    let (_, h) = main_vcf();
    parse_expr(expr, &h).unwrap_err().0
}

#[test]
fn precedence_and_parentheses() {
    // && binds tighter than ||.
    assert_eq!(hits("POS=10 || POS=20 && POS=30"), vec![10]);
    assert_eq!(hits("(POS=10 || POS=20) && POS=30"), Vec::<u64>::new());
    assert_eq!(hits("(POS=10 || POS=20) && QUAL>1"), vec![10]);
    // ! binds tighter than &&.
    assert_eq!(hits("!POS=10 && POS<30"), vec![20]);
    assert_eq!(hits("!(POS=10 || POS=20) && POS<40"), vec![30]);
    assert_eq!(hits("!!(POS=10)"), vec![10]);
    // `=` and `==` agree.
    assert_eq!(hits("POS==10"), hits("POS=10"));
}

#[test]
fn spec_examples() {
    assert_eq!(hits("QUAL>30"), vec![10]); // missing QUAL is false
    assert_eq!(hits("!(QUAL<20 || INFO/DP<5)"), vec![10, 20]);
    assert_eq!(hits("FILTER=\"PASS\""), vec![10, 40]);
    assert_eq!(hits("ID~\"^rs\""), vec![10]);
    assert_eq!(hits("TYPE=\"indel\""), vec![30]);
    assert_eq!(hits("N_PASS(FMT/GQ>30)>=2"), vec![10, 30]);
    assert_eq!(hits("INFO/AF[1]>0.1"), vec![20]);
    assert_eq!(hits("MAF<0.3"), vec![30, 40]);
}

#[test]
fn deep_nesting_is_an_error_not_a_stack_overflow() {
    let parens = format!("{}QUAL>1{}", "(".repeat(5000), ")".repeat(5000));
    let e = err(&parens);
    assert!(e.contains("expression nested too deeply"), "{e}");
    assert!(e.contains("column 129"), "{e}");
    assert!(
        e.len() < 600,
        "the expression is cut in the message: {}",
        e.len()
    );
    let nots = format!("{}QUAL>1", "!".repeat(20000));
    let e = err(&nots);
    assert!(e.contains("expression nested too deeply"), "{e}");
    assert!(e.contains("column 128"), "{e}");
    // Nesting within the limit still parses.
    let ok = format!("{}QUAL>1{}", "(".repeat(100), ")".repeat(100));
    let (_, h) = main_vcf();
    assert!(parse_expr(&ok, &h).is_ok());
    let ok = format!("{}QUAL>1", "!".repeat(100));
    assert!(parse_expr(&ok, &h).is_ok());
}

#[test]
fn parse_errors_name_the_column() {
    let e = err("QUAL>30 && ?");
    assert!(e.contains("column 12"), "{e}");
    assert!(e.contains("\"QUAL>30 && ?\""), "{e}");
    let e = err("QUAL>30 &&");
    assert!(e.contains("column 11") && e.contains("ends here"), "{e}");
    let e = err("(QUAL>30");
    assert!(e.contains("column 9") && e.contains("`)`"), "{e}");
    let e = err("QUAL > 30 30");
    assert!(e.contains("column 11"), "{e}");
    let e = err("ID~\"(\"");
    assert!(
        e.contains("column 4") && e.contains("regular expression"),
        "{e}"
    );
    let e = err("ID=\"abc");
    assert!(e.contains("column 4") && e.contains("closing"), "{e}");
    let e = err("QUAL");
    assert!(
        e.contains("column 1") && e.contains("without a comparison"),
        "{e}"
    );
}

#[test]
fn unsupported_constructs_are_usage_errors_with_hints() {
    let e = err("QUAL>30 & POS>1");
    assert!(
        e.contains("`&` is not supported; use && / ||") && e.contains("column 9"),
        "{e}"
    );
    let e = err("QUAL>30 | POS>1");
    assert!(e.contains("`|` is not supported; use && / ||"), "{e}");
    let e = err("ABS(QUAL)>1");
    assert!(
        e.contains("ABS is not supported") && e.contains("N_PASS") && e.contains("column 1"),
        "{e}"
    );
    let e = err("QUAL+1>10");
    assert!(e.contains("arithmetic is not supported"), "{e}");
    let e = err("QUAL>POS");
    assert!(e.contains("two fields") && e.contains("column 5"), "{e}");
    let e = err("FMT/AD>1");
    assert!(e.contains("not defined in the header"), "{e}");
    let e = err("N_PASS(QUAL>1)>1");
    assert!(e.contains("per-sample"), "{e}");
    let e = err("TYPE=\"weird\"");
    assert!(e.contains("snp, mnp, indel, other, ref"), "{e}");
    let e = err("QUAL=\"x\"");
    assert!(e.contains("numeric") || e.contains("number"), "{e}");
    let e = err("FILTER=\"Nope\"");
    assert!(e.contains("not defined in the header"), "{e}");
}

#[test]
fn unknown_keys_and_ambiguity_are_usage_errors() {
    let e = err("INFO/NOPE>1");
    assert!(e.contains("INFO/NOPE is not defined"), "{e}");
    let e = err("NOPE>1");
    assert!(e.contains("\"NOPE\"") && e.contains("not defined"), "{e}");
    // DP is both INFO and FORMAT in the fixture, as in bcftools.
    let e = err("DP>10");
    assert!(
        e.contains("ambiguous") && e.contains("INFO/DP") && e.contains("FMT/DP"),
        "{e}"
    );
    let e = err("DP>10");
    assert!(e.contains("write INFO/DP or FMT/DP"), "{e}");
    let e = err("FOO/DP>1");
    assert!(e.contains("unknown prefix"), "{e}");
    // A key defined once resolves bare.
    assert_eq!(hits("AF[1]>0.25"), vec![20]);
    assert_eq!(hits("GQ>45"), vec![10]);
}

#[test]
fn needs_flags() {
    let (_, h) = main_vcf();
    let n = |s: &str| {
        let e = parse_expr(s, &h).unwrap();
        (e.needs_genotypes(), e.needs_stats(), e.needs_info_values())
    };
    assert_eq!(n("QUAL>1 && ID~\"rs\""), (false, false, false));
    assert_eq!(n("INFO/DP>1"), (false, false, true));
    assert_eq!(n("FMT/DP>1"), (true, false, false));
    assert_eq!(n("GT=\"het\""), (true, false, false));
    assert_eq!(n("N_PASS(GT=\"het\")>1"), (true, false, false));
    assert_eq!(n("MAF<0.1"), (true, true, false));
    assert_eq!(n("N_SAMPLES>1"), (false, false, false));
}

#[test]
fn no_genotypes_supplied_is_an_error_not_a_panic() {
    let (v, h) = main_vcf();
    let e = parse_expr("MAF<0.1", &h).unwrap();
    let line = records(&v)[0];
    let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
    let msg = e
        .eval(&site, None, None, &mut EvalScratch::default())
        .unwrap_err();
    assert!(msg.contains("MAF"), "{msg}");
}

#[test]
fn masked_samples_read_as_missing() {
    let (v, h) = main_vcf();
    // Mask samples with GQ < 45: record 10 keeps only s3 (GQ 50) with 1/1.
    let rules = MaskRules {
        min_gq: Some(45),
        ..Default::default()
    };
    let line = records(&v)[0];
    let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
    let g = text_genotypes(line.as_bytes()).unwrap();
    let mut masked = Vec::new();
    let stats = site_stats(&g, None, &rules, 2, &mut masked).unwrap();
    let mg = MaskedGenotypes {
        inner: &g,
        masked: &masked,
    };
    let mut sc = EvalScratch::default();
    let mut t = |s: &str| {
        let e = parse_expr(s, &h).unwrap();
        e.eval(&site, Some(&mg), Some(&stats), &mut sc).unwrap()
    };
    assert!(!t("GT=\"het\"")); // s1's 0/1 is masked
    assert!(t("GT=\"hom\""));
    assert!(t("N_PASS(GT=\"mis\")=2"));
    assert!(t("N_MISSING=2"));
    assert!(t("AN=2"));
    assert!(t("FMT/GQ>45") && !t("FMT/GQ>50"));
    assert!(t("FMT/DP=30") && !t("FMT/DP=25"));
}

#[test]
fn sample_selection_limits_fmt_and_npass() {
    let (v, h) = main_vcf();
    let line = records(&v)[0]; // s1 0/1, s2 0/0, s3 1/1
    let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
    let g = text_genotypes(line.as_bytes()).unwrap();
    let mut sc = EvalScratch::default();
    sc.set_selection(Some(&[false, true, false]));
    let mut t = |s: &str| {
        parse_expr(s, &h)
            .unwrap()
            .eval(&site, Some(&g), None, &mut sc)
            .unwrap()
    };
    assert!(!t("GT=\"het\"") && !t("GT=\"alt\"") && t("GT=\"ref\""));
    assert!(t("N_PASS(GT=\"ref\")=1") && t("N_PASS(FMT/DP>0)=1"));
    sc.set_selection(None);
    let e = parse_expr("GT=\"het\"", &h).unwrap();
    assert!(e.eval(&site, Some(&g), None, &mut sc).unwrap());
}

#[test]
fn eval_does_not_allocate_after_warmup() {
    let (v, h) = main_vcf();
    let e = parse_expr(
        "(QUAL>1 || INFO/AF[1]>0.1 || ID~\"^rs\") && N_PASS(FMT/GQ>30)>=0 && GT!=\"mis\" \
         && FMT/DP>0 && TYPE!=\"indel\" && INFO/S!=\"x\" && MAF>=0 && AC>=0",
        &h,
    )
    .unwrap();
    let mut sc = EvalScratch::default();
    let mut masked = Vec::new();
    let mut prepared = Vec::new();
    for line in records(&v) {
        let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
        let g = text_genotypes(line.as_bytes()).unwrap();
        let stats = site_stats(
            &g,
            None,
            &MaskRules::default(),
            n_alleles(site.alts),
            &mut masked,
        )
        .unwrap();
        prepared.push((site, g, stats));
    }
    for (s, g, st) in &prepared {
        e.eval(s, Some(g), Some(st), &mut sc).unwrap(); // warm up
    }
    let before = ALLOCS.with(Cell::get);
    for _ in 0..50 {
        for (s, g, st) in &prepared {
            e.eval(s, Some(g), Some(st), &mut sc).unwrap();
        }
    }
    assert_eq!(ALLOCS.with(Cell::get), before, "eval allocated");
}

// ---- pinned differences from bcftools (docs/differences.md) ---------------

fn bcftools_pinned(fx: &str, expr: &str) -> Vec<u64> {
    let golden = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/expr_differences.tsv"
    ))
    .unwrap();
    let row = golden
        .lines()
        .find(|r| r.starts_with(&format!("{fx}\t{expr}\t")))
        .unwrap_or_else(|| panic!("no pinned row for {fx} {expr}"));
    pos_list(row.rsplit('\t').next().unwrap())
}

fn gloci_on(fx: &str, expr: &str) -> Vec<u64> {
    let v = fixture(fx);
    let e = parse_expr(expr, &header_of(&v)).unwrap();
    let t = run_text(&e, &v);
    assert_eq!(run_bcf(&e, &v), t, "bcf vs text: {expr}");
    t
}

#[test]
fn difference_not_match_is_a_real_negation() {
    // bcftools 1.21 evaluates `!~` like `~`.
    let exp = "ID!~\"^rs\"";
    assert_eq!(bcftools_pinned("main", exp), vec![10]);
    assert_eq!(gloci_on("main", exp), vec![20, 30, 40]);
}

#[test]
fn difference_index_past_the_end_is_missing() {
    // bcftools 1.21 matches site 10 (one AF value) for AF[5]>0.
    let exp = "INFO/AF[5]>0";
    assert_eq!(bcftools_pinned("main", exp), vec![10]);
    assert_eq!(gloci_on("main", exp), Vec::<u64>::new());
}

#[test]
fn difference_bare_flag_is_true_when_present() {
    // bcftools 1.21 matches no record for a bare flag.
    let exp = "INFO/DB";
    assert_eq!(bcftools_pinned("main", exp), Vec::<u64>::new());
    assert_eq!(gloci_on("main", exp), vec![10]);
}

#[test]
fn difference_negated_group_parses() {
    // bcftools 1.21 cannot parse `!(`.
    let exp = "!(QUAL<20)";
    assert_eq!(bcftools_pinned("main", exp), Vec::<u64>::new());
    assert_eq!(gloci_on("main", exp), vec![10, 20]);
}

#[test]
fn difference_missing_literal_compares() {
    // bcftools 1.21 rejects `QUAL=.` (the tag "." is not defined).
    assert_eq!(bcftools_pinned("main", "QUAL=."), Vec::<u64>::new());
    assert_eq!(gloci_on("main", "QUAL=."), vec![20]);
    assert_eq!(gloci_on("main", "QUAL!=."), vec![10, 30, 40]);
}

#[test]
fn difference_field_vs_field_is_rejected() {
    // bcftools supports `INFO/DP>FMT/DP`; gloci's subset does not.
    let (_, h) = main_vcf();
    assert!(parse_expr("INFO/DP>FMT/DP", &h).is_err());
}

#[test]
fn difference_flipped_operands_are_accepted() {
    assert_eq!(gloci_on("main", "30<QUAL"), vec![10]);
}

#[test]
fn difference_maf_is_first_alt_not_second_most_common() {
    // `MAF` follows bcftools (first ALT, min(f, 1-f)); `SiteStats::maf`
    // (--min-maf) is the second most common allele and differs on m.vcf
    // site 2 (counts 2,6,2 of 10): 0.4 against 0.2.
    let v = fixture("m");
    let line = records(&v)[1];
    let g = text_genotypes(line.as_bytes()).unwrap();
    let mut masked = Vec::new();
    let st: SiteStats = site_stats(&g, None, &MaskRules::default(), 3, &mut masked).unwrap();
    assert_eq!(st.expr_maf(), Some(0.4));
    assert_eq!(st.maf(), Some(0.2));
    assert_eq!(st.ac(), &[6, 2]);
    assert_eq!(st.an(), 10);
}

#[test]
fn bcf_info_values_render_like_text() {
    let (v, h) = main_vcf();
    let (dicts, raws) = to_bcf(&v);
    let mut bs = BcfSite::with_info_values();
    let want = ["DP=20;AF=0.1;DB;S=abc", "AF=0.2,0.3", "DP=3", "."];
    for ((raw, line), w) in raws.iter().zip(records(&v)).zip(want) {
        let site = bs.view(raw, &dicts, span_of(line)).unwrap();
        assert_eq!(String::from_utf8_lossy(site.info), w);
    }
    let _ = h;
    let _: Option<&dyn GenotypeSource> = None;
}

// ---- later fixes -----------------------------------------------------------

#[test]
fn bare_ac_an_maf_read_header_info_and_compute_otherwise() {
    // infoac.vcf defines INFO AC, AN and MAF (goldens pin the bcftools rows).
    let v = fixture("infoac");
    let h = header_of(&v);
    let t = |s: &str| {
        let e = parse_expr(s, &h).unwrap();
        (e.needs_stats(), e.needs_info_values(), run_text(&e, &v))
    };
    assert_eq!(t("AN=99"), (false, true, vec![10]));
    assert_eq!(t("MAF<0.02"), (false, true, vec![10]));
    assert_eq!(t("AC=2"), (false, true, vec![10])); // any value of AC=3,2
    // An INFO-resolved name is missing where the record lacks the key.
    assert_eq!(t("AN!=99").2, vec![20, 30, 40]);
    // F_MISSING / N_MISSING stay computed even when the header defines them
    // (fmiss.vcf defines both; the goldens pin the bcftools rows).
    assert!(t("N_MISSING=0").0);
    let v = fixture("fmiss");
    let h = header_of(&v);
    let run = |s: &str| run_text(&parse_expr(s, &h).unwrap(), &v);
    assert_eq!(run("N_MISSING=1"), vec![10]); // not the INFO value 7
    assert_eq!(run("N_MISSING=7"), Vec::<u64>::new());
    assert_eq!(run("F_MISSING>0.5"), vec![20]);
    assert_eq!(run("INFO/N_MISSING=7"), vec![10]); // the INFO key stays reachable
    // Without the header key the computed value is used.
    let (v, h) = main_vcf();
    let e = parse_expr("AN=6", &h).unwrap();
    assert!(e.needs_stats() && !e.needs_info_values());
    assert_eq!(run_text(&e, &v), vec![10]);
}

#[test]
fn sample_subscripts_are_sample_numbers_over_the_whole_header() {
    // main.vcf: s1 s2 s3; goldens pin the bcftools rows.
    assert_eq!(hits("FMT/GQ[1]>35"), vec![30]);
    assert_eq!(hits("FMT/GQ[0]>35"), vec![10]);
    assert_eq!(hits("GT[1]=\"het\""), vec![20]);
    assert_eq!(hits("N_PASS(FMT/GQ[1]>30)>=1"), vec![30]);
    // bcftools -s numbers over the full header: the subscript still means the
    // original column, and a sample outside the selection can be subscripted.
    let (v, h) = main_vcf();
    let mut sc = EvalScratch::default();
    sc.set_selection(Some(&[false, true, true]));
    let line = records(&v)[0]; // GQ 40, 20, 50
    let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
    let g = text_genotypes(line.as_bytes()).unwrap();
    let mut t = |s: &str| {
        parse_expr(s, &h)
            .unwrap()
            .eval(&site, Some(&g), None, &mut sc)
            .unwrap()
    };
    assert!(t("FMT/GQ[0]>35")); // s1 is not selected but is sample 0
    assert!(!t("FMT/GQ[1]>35"));
    assert!(t("FMT/GQ>35")); // unsubscripted: any selected sample (s3)
}

#[test]
fn sample_subscript_past_the_header_is_a_parse_error() {
    let e = err("FMT/GQ[7]>1");
    assert!(
        e.contains("sample index 7") && e.contains("3 samples") && e.contains("column 8"),
        "{e}"
    );
    let e = err("GT[3]=\"het\"");
    assert!(e.contains("sample index 3"), "{e}");
    assert!(err("GT[x]=\"het\"").contains("whole number"));
}

#[test]
fn bcf_info_edge_values_render_and_evaluate_like_text() {
    let v = fixture("infotypes");
    let (dicts, raws) = to_bcf(&v);
    let mut bs = BcfSite::with_info_values();
    let want = [
        "DP=300;AF=.,0.1",
        "DP=70000;AF=0.2,.",
        "DP=7;AF=0.5",
        "AD=1,2,300,.;S=x,y",
        "AD=.,.",
    ];
    for ((raw, line), w) in raws.iter().zip(records(&v)).zip(want) {
        let site = bs.view(raw, &dicts, span_of(line)).unwrap();
        assert_eq!(String::from_utf8_lossy(site.info), w);
    }
}

#[test]
fn eval_refuses_a_keys_only_bcf_view() {
    let v = fixture("main");
    let h = header_of(&v);
    let (dicts, raws) = to_bcf(&v);
    let mut keys_only = BcfSite::default();
    let g = bcf_genotypes(&raws[0], &dicts).unwrap();
    let site = keys_only
        .view(&raws[0], &dicts, span_of(records(&v)[0]))
        .unwrap();
    assert!(!site.info_has_values);
    let mut sc = EvalScratch::default();
    let needs = parse_expr("INFO/DP>10", &h).unwrap();
    let msg = needs.eval(&site, Some(&g), None, &mut sc).unwrap_err();
    assert!(msg.contains("with_info_values"), "{msg}");
    // An expression that does not read INFO is fine on the same view.
    let fine = parse_expr("QUAL>1", &h).unwrap();
    assert!(fine.eval(&site, Some(&g), None, &mut sc).unwrap());
}

#[test]
fn unsupported_format_key_is_a_usage_error() {
    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=1,length=10>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
    let h = header_of(vcf);
    for ex in ["FMT/AD>1", "AD>1"] {
        let e = parse_expr(ex, &h).unwrap_err().0;
        assert!(
            e.contains("FMT/AD is not supported; supported FORMAT fields: GT, DP, GQ")
                && e.contains("column 1"),
            "{e}"
        );
    }
}

#[test]
fn difference_negated_bare_comparison_parses() {
    // bcftools 1.21: `!POS=10` is "Could not parse the expression".
    assert_eq!(hits("!POS=10 && POS<30"), vec![20]);
}

// ---- vcf query's per-sample rule (D172) -----------------------------------

use gloci_vcf::QueryPass;

/// `eval_query` over `vcf`, as the TSV writes it: `POS:SAMPLES` per printed
/// record, space separated (`All` is every selected sample).
fn run_query(
    e: &Expr,
    exclude: bool,
    vcf: &str,
    keep: Option<&[bool]>,
    rules: MaskRules,
    bcf: bool,
) -> String {
    let h = header_of(vcf);
    let names: Vec<String> = h
        .samples
        .iter()
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect();
    let (dicts, raws) = to_bcf(vcf);
    let mut scratch = EvalScratch::default();
    scratch.set_selection(keep);
    let mut masked = Vec::new();
    let mut bs = BcfSite::with_info_values();
    let mut pass = Vec::new();
    let mut out = Vec::new();
    for (n, line) in records(vcf).into_iter().enumerate() {
        let (site, g): (SiteView, Box<dyn GenotypeSource>) = if bcf {
            (
                bs.view(&raws[n], &dicts, span_of(line)).unwrap(),
                Box::new(bcf_genotypes(&raws[n], &dicts).unwrap()),
            )
        } else {
            (
                SiteView::from_text(line.as_bytes(), span_of(line)).unwrap(),
                Box::new(text_genotypes(line.as_bytes()).unwrap()),
            )
        };
        let stats = site_stats(&*g, keep, &rules, n_alleles(site.alts), &mut masked).unwrap();
        let mg = MaskedGenotypes {
            inner: &*g,
            masked: &masked,
        };
        let r = e
            .eval_query(
                exclude,
                &site,
                Some(&mg),
                Some(&stats),
                &mut scratch,
                &mut pass,
            )
            .unwrap();
        let sel = |i: usize| keep.is_none_or(|k| k[i]);
        let who: Vec<&str> = match r {
            QueryPass::Drop => continue,
            QueryPass::All => (0..names.len())
                .filter(|&i| sel(i))
                .map(|i| &*names[i])
                .collect(),
            QueryPass::Samples => {
                assert_eq!(pass.len(), names.len());
                (0..names.len())
                    .filter(|&i| pass[i])
                    .map(|i| &*names[i])
                    .collect()
            }
        };
        out.push(format!("{}:{}", site.pos, who.join(",")));
    }
    out.join(" ")
}

fn q(expr: &str, exclude: bool, keep: Option<&[bool]>, rules: MaskRules) -> String {
    let (v, h) = main_vcf();
    run_query(
        &parse_expr(expr, &h).unwrap(),
        exclude,
        &v,
        keep,
        rules,
        false,
    )
}

#[test]
fn query_rule_matches_bcftools_on_text_and_bcf() {
    let (v, h) = main_vcf();
    let golden = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/expr_query_bcftools.tsv"
    ))
    .unwrap();
    let mut rows = 0;
    for row in golden.lines() {
        let f: Vec<&str> = row.split('\t').collect();
        let e = parse_expr(f[0], &h).unwrap();
        for bcf in [false, true] {
            let got = run_query(&e, f[1] == "e", &v, None, MaskRules::default(), bcf);
            assert_eq!(got, f[2], "{} -{} (bcf={bcf})", f[0], f[1]);
        }
        rows += 1;
    }
    assert!(rows >= 40, "{rows} rows");
}

#[test]
fn not_flips_each_selected_sample() {
    let d = MaskRules::default;
    assert_eq!(
        q("!(GT=\"het\")", false, None, d()),
        "10:s2,s3 20:s1,s3 30:s1,s2,s3 40:s2,s3"
    );
    assert_eq!(q("!(QUAL>8)", false, None, d()), "20:s1,s2,s3 40:s1,s2,s3");
}

#[test]
fn only_selected_samples_pass() {
    let keep = [false, true, true];
    let k = Some(&keep[..]);
    let d = MaskRules::default;
    assert_eq!(q("GT=\"het\"", false, k, d()), "20:s2");
    assert_eq!(q("GT[0]=\"het\"", false, k, d()), "");
    assert_eq!(
        q("GT=\"het\" || QUAL>8", false, k, d()),
        "10:s2,s3 20:s2 30:s2,s3"
    );
    assert_eq!(
        q("GT=\"het\"", true, k, d()),
        "10:s2,s3 20:s3 30:s2,s3 40:s2,s3"
    );
}

#[test]
fn masked_samples_read_missing() {
    let rules = MaskRules {
        min_dp: Some(20),
        ..Default::default()
    };
    assert_eq!(q("GT=\"het\"", false, None, rules), "10:s1 40:s1");
}

#[test]
fn eval_query_does_not_allocate_once_warm() {
    let (v, h) = main_vcf();
    let e = parse_expr("GT=\"het\" && FMT/DP>5 || (QUAL<8 && GT=\"alt\")", &h).unwrap();
    let mut sc = EvalScratch::default();
    let mut masked = Vec::new();
    let mut pass = Vec::new();
    let mut prepared = Vec::new();
    for line in records(&v) {
        let site = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
        let g = text_genotypes(line.as_bytes()).unwrap();
        let stats = site_stats(
            &g,
            None,
            &MaskRules::default(),
            n_alleles(site.alts),
            &mut masked,
        )
        .unwrap();
        prepared.push((site, g, stats));
    }
    for exclude in [false, true] {
        for (s, g, st) in &prepared {
            e.eval_query(exclude, s, Some(g), Some(st), &mut sc, &mut pass)
                .unwrap();
        }
        let before = ALLOCS.with(Cell::get);
        for _ in 0..50 {
            for (s, g, st) in &prepared {
                e.eval_query(exclude, s, Some(g), Some(st), &mut sc, &mut pass)
                    .unwrap();
            }
        }
        assert_eq!(
            ALLOCS.with(Cell::get),
            before,
            "eval_query allocated (exclude={exclude})"
        );
    }
}
