//! `Splitter`: the Number=A/R/G rules, GT recoding and the types gate.
use gloci_vcf::{SplitOutcome, SplitTypes, Splitter};

const H: &str = "##INFO=<ID=A1,Number=A,Type=Integer,Description=\"a\">\n\
##INFO=<ID=R1,Number=R,Type=Integer,Description=\"r\">\n\
##INFO=<ID=G1,Number=G,Type=Integer,Description=\"g\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"r\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"g\">\n\
##FORMAT=<ID=AQ,Number=A,Type=Integer,Description=\"a\">\n";

fn split(line: &str, t: SplitTypes, ns: usize) -> Result<Vec<String>, String> {
    let s = Splitter::new(H.as_bytes(), t, ns);
    let (mut out, mut ends) = (Vec::new(), Vec::new());
    Ok(match s.split(line.as_bytes(), &mut out, &mut ends)? {
        SplitOutcome::Unchanged => vec![line.to_string()],
        SplitOutcome::Split(_) => {
            let mut p = 0;
            ends.iter()
                .map(|&e| {
                    let r = String::from_utf8(out[p..e].to_vec()).unwrap();
                    p = e;
                    r
                })
                .collect()
        }
    })
}

fn fields(rec: &str) -> Vec<&str> {
    rec.split('\t').collect()
}

#[test]
fn split_number_rules() {
    // Diploid, haploid, phased and missing GT; A/R/G in INFO and FORMAT.
    let line = "c\t5\tid\tC\tA,G\t9\tq\tA1=1,2;R1=5,6,7;G1=1,2,3,4,5,6;N=1\tGT:AD:PL:AQ\t\
1/2:1,2,3:1,2,3,4,5,6:8,9\t2|2:1,2,3:1,2,3,4,5,6:8,9\t1:4,5,6:7,8,9:8,9\t./.:.:.:.\t0/1:1,2,3:.:.";
    let r = split(line, SplitTypes::Any, 5).unwrap();
    assert_eq!(r.len(), 2);
    let a = fields(&r[0]);
    assert_eq!(
        &a[..8],
        ["c", "5", "id", "C", "A", "9", "q", "A1=1;R1=5,6;G1=1,2,3"]
    );
    assert_eq!(a[8], "GT:AD:PL:AQ");
    assert_eq!(a[9], "1/0:1,2:1,2,3:8");
    assert_eq!(a[10], "0|0:1,2:1,2,3:8");
    // haploid G (3 values for 3 alleles): alleles 0 and 1
    assert_eq!(a[11], "1:4,5:7,8:8");
    assert_eq!(a[12], "./.:.,.:.:.");
    assert_eq!(a[13], "0/1:1,2:.:.");
    let g = fields(&r[1]);
    assert_eq!(g[7], "A1=2;R1=5,7;G1=1,4,6");
    assert_eq!(g[9], "0/1:1,3:1,4,6:9");
    assert_eq!(g[10], "1|1:1,3:1,4,6:9");
    assert_eq!(g[11], "0:4,6:7,9:9");
    assert_eq!(g[13], "0/0:1,3:.:.");
}

#[test]
fn single_alt_and_dot_are_unchanged() {
    for alt in ["T", "."] {
        let l = format!("c\t5\t.\tA\t{alt}\t.\t.\tA1=1,2\tGT\t0/1");
        assert_eq!(split(&l, SplitTypes::Any, 1).unwrap(), vec![l.clone()]);
    }
}

#[test]
fn undeclared_keys_are_dropped_and_alleles_are_not_trimmed() {
    let l = "c\t5\t.\tCAA\tCA,CAAA\t.\t.\tNOPE=1;A1=3,4\tGT:XX\t1/2:9";
    let r = split(l, SplitTypes::Any, 1).unwrap();
    assert_eq!(r[0], "c\t5\t.\tCAA\tCA\t.\t.\tA1=3\tGT\t1/0");
    assert_eq!(r[1], "c\t5\t.\tCAA\tCAAA\t.\t.\tA1=4\tGT\t0/1");
}

#[test]
fn wrong_counts_are_errors() {
    let e = split("c\t5\t.\tA\tC,G\t.\t.\tA1=1\tGT\t0/1", SplitTypes::Any, 1).unwrap_err();
    assert!(e.contains("INFO/A1") && e.contains("needs 2"), "{e}");
    let e = split(
        "c\t5\t.\tA\tC,G\t.\t.\t.\tGT:AD\t0/1:1,2",
        SplitTypes::Any,
        1,
    )
    .unwrap_err();
    assert!(e.contains("FORMAT/AD") && e.contains("needs 3"), "{e}");
    let e = split("c\t5\t.\tA\tC,G\t.\t.\t.\tGT\t0/1", SplitTypes::Any, 2).unwrap_err();
    assert!(e.contains("sample columns"), "{e}");
}

#[test]
fn split_types() {
    let n = |ref_: &str, alts: &str, t| {
        let l = format!("c\t5\t.\t{ref_}\t{alts}\t.\t.\t.");
        split(&l, t, 0).unwrap().len()
    };
    // SNP + SNP, indel + SNP, MNP + MNP, complex + symbolic, trimmed SNP, indels
    let cases = [
        ("A", "C,G", [2, 1, 2, 2]),
        ("AT", "A,C", [1, 2, 2, 2]),
        ("AT", "CG,GA", [1, 1, 2, 2]),
        ("ACG", "TT,<X>", [1, 1, 2, 2]),
        ("TA", "GA,<X>", [2, 1, 2, 2]),
        ("AC", "A,ACC", [1, 2, 2, 2]),
    ];
    for (r, a, want) in cases {
        let got = [
            n(r, a, SplitTypes::Snps),
            n(r, a, SplitTypes::Indels),
            n(r, a, SplitTypes::Both),
            n(r, a, SplitTypes::Any),
        ];
        let want_n = a.split(',').count();
        let want: Vec<usize> = want
            .iter()
            .map(|&k| if k == 1 { 1 } else { want_n })
            .collect();
        assert_eq!(got.to_vec(), want, "{r}>{a}");
    }
}
