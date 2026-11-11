//! Genotype concordance of paired records (decision D163).

use gloci_core::ChromId;
use gloci_vcf::compare::{
    AlleleMap, DiscordanceMatrix, Gt, GtCompare, MatchRule, MatrixCell, Site,
};

fn site(pos: u64, r: &str, alts: &str) -> Site {
    Site::bare(ChromId(0), pos, r, alts)
}

fn dip(a: u32, b: u32) -> Gt {
    Gt::Dip(a, b)
}

fn cmp(rule: MatchRule, a: &Site, b: &Site, ga: &Gt, gb: &Gt) -> GtCompare {
    let mut m = AlleleMap::default();
    m.build(a, b, rule);
    m.compare(ga, gb)
}

fn same(a: &Site, b: &Site, ga: &Gt, gb: &Gt) -> GtCompare {
    cmp(MatchRule::Position, a, b, ga, gb)
}

#[test]
fn genotype_compare_table() {
    let s = site(10, "A", "C");
    // 0/1 and 1/0 are the same genotype
    assert_eq!(same(&s, &s, &dip(0, 1), &dip(1, 0)), GtCompare::Concordant);
    assert_eq!(same(&s, &s, &dip(0, 1), &dip(0, 0)), GtCompare::Discordant);
    // alleles are what they are, not their index: A>C,G against A>G,C
    let (a, b) = (site(10, "A", "C,G"), site(10, "A", "G,C"));
    assert_eq!(same(&a, &b, &dip(0, 1), &dip(0, 2)), GtCompare::Concordant);
    assert_eq!(same(&a, &b, &dip(0, 1), &dip(0, 1)), GtCompare::Discordant);
    assert_eq!(same(&a, &b, &dip(1, 2), &dip(2, 1)), GtCompare::Concordant);
    assert_eq!(same(&a, &b, &dip(2, 2), &dip(1, 1)), GtCompare::Concordant);
    assert!(!{
        let mut m = AlleleMap::default();
        m.build(&a, &b, MatchRule::Position);
        m.matching
    });
    // missing on either side, or on both
    assert_eq!(
        same(&s, &s, &dip(0, 1), &Gt::Missing),
        GtCompare::SecondMissing
    );
    assert_eq!(
        same(&s, &s, &Gt::Missing, &dip(0, 1)),
        GtCompare::FirstMissing
    );
    assert_eq!(
        same(&s, &s, &Gt::Missing, &Gt::Missing),
        GtCompare::BothMissing
    );
    // half-missing alleles make a genotype missing
    assert_eq!(Gt::from_alleles(&[None, Some(1)]), Gt::Missing);
    assert_eq!(Gt::from_alleles(&[Some(0), None]), Gt::Missing);
    assert_eq!(Gt::from_alleles(&[None]), Gt::Missing);
    assert_eq!(Gt::from_alleles(&[]), Gt::Missing);
    assert_eq!(Gt::from_alleles(&[Some(1), Some(0)]), dip(1, 0));
    // haploid: equal alleles are concordant; haploid against diploid is not
    assert_eq!(
        same(&s, &s, &Gt::Hap(1), &Gt::Hap(1)),
        GtCompare::Concordant
    );
    assert_eq!(
        same(&s, &s, &Gt::Hap(1), &Gt::Hap(0)),
        GtCompare::Discordant
    );
    assert_eq!(same(&s, &s, &Gt::Hap(1), &dip(1, 1)), GtCompare::Discordant);
    // polyploid compare as multisets
    let p = |v: &[u32]| Gt::Poly(v.to_vec().into());
    assert_eq!(
        same(&s, &s, &p(&[0, 0, 1]), &p(&[1, 0, 0])),
        GtCompare::Concordant
    );
    assert_eq!(
        same(&s, &s, &p(&[0, 0, 1]), &p(&[0, 1, 1])),
        GtCompare::Discordant
    );
    // a REF of N takes the other's REF in position mode
    let (n, r) = (site(10, "N", "C"), site(10, "ACG", "C"));
    assert_eq!(same(&n, &r, &dip(0, 0), &dip(0, 0)), GtCompare::Concordant);
    assert_eq!(same(&n, &r, &dip(0, 1), &dip(0, 1)), GtCompare::Concordant);
    // trimmed alleles are the same allele in alleles mode
    let (x, y) = (site(10, "AT", "A"), site(10, "ATG", "AG"));
    assert_eq!(
        cmp(MatchRule::Alleles, &x, &y, &dip(0, 1), &dip(1, 0)),
        GtCompare::Concordant
    );
    // different ALT with the same index: discordant
    let (c, g) = (site(10, "A", "C"), site(10, "A", "G"));
    assert_eq!(same(&c, &g, &dip(1, 1), &dip(1, 1)), GtCompare::Discordant);
    assert_eq!(same(&c, &g, &dip(0, 0), &dip(0, 0)), GtCompare::Concordant);
}

#[test]
fn discordance_matrix_rules() {
    let bi = |a: &Site, b: &Site| {
        let mut m = AlleleMap::default();
        m.build(a, b, MatchRule::Position);
        m.biallelic_same
    };
    // only biallelic pairs with the same REF and ALT
    assert!(bi(&site(1, "A", "C"), &site(1, "A", "C")));
    assert!(bi(&site(1, "N", "C"), &site(1, "A", "C")));
    assert!(!bi(&site(1, "A", "C"), &site(1, "A", "G")));
    assert!(!bi(&site(1, "A", "C,G"), &site(1, "A", "C,G")));
    assert!(!bi(&site(1, "A", "."), &site(1, "A", ".")));
    // the cells: 0/0, 0/1, 1/1, ./.; haploid and polyploid are skipped
    assert_eq!(AlleleMap::cell(&dip(0, 0)), MatrixCell::Cell(0));
    assert_eq!(AlleleMap::cell(&dip(1, 0)), MatrixCell::Cell(1));
    assert_eq!(AlleleMap::cell(&dip(1, 1)), MatrixCell::Cell(2));
    assert_eq!(AlleleMap::cell(&Gt::Missing), MatrixCell::Cell(3));
    assert_eq!(AlleleMap::cell(&Gt::Hap(1)), MatrixCell::Skip);
    assert_eq!(
        AlleleMap::cell(&Gt::Poly(vec![0, 0, 1].into())),
        MatrixCell::Skip
    );
    // file 2 is the row, file 1 the column
    let mut m = DiscordanceMatrix::default();
    m.add(1, 0);
    assert_eq!(m.cells[0][1], 1);
}
