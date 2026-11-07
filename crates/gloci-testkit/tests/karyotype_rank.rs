use gloci_core::karyotype_cmp;
use gloci_testkit::naive::karyotype_rank;

#[test]
fn naive_rank_agrees_with_core_comparator() {
    let names = [
        "chrM",
        "chr10",
        "chr1_KI270706v1_random",
        "chr2",
        "chrY",
        "chr1",
        "chrX",
        "1",
        "chr01",
        "chrUn_GL000220v1",
        "chr22",
        "chrEBV",
        "HLA-A*01:01:01:01",
        "MT",
        "chr2L",
        "chr2A",
        "chr99999999999999999999999",
        "chr100000000000000000000000",
        "scaf9",
        "scaf10",
        "scaf09",
    ];
    let rank = karyotype_rank(&names);
    for a in names {
        for b in names {
            assert_eq!(
                rank(a).cmp(&rank(b)),
                karyotype_cmp(a.as_bytes(), b.as_bytes()),
                "{a} vs {b}"
            );
        }
    }
}
