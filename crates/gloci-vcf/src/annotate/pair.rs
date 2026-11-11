//! Pairing of the records of FILE and the annotation source at one position
//! (`--pair-logic some`, bcftools annotate's default): htslib's `bcf_sr_sort`
//! in [`crate::srsort`].

pub use crate::srsort::*;

/// The rows of a position that both files have records at: each row has the
/// index of a record of `file`, of `src`, or both, in bcftools' order.
pub fn pair_rows(file: &[PairRec], src: &[PairRec]) -> Vec<(Option<usize>, Option<usize>)> {
    sort_rows(&[file, src], Pair::Some)
        .into_iter()
        .map(|row| {
            let mut r = (None, None);
            for (reader, rec) in row {
                if reader == 0 {
                    r.0 = Some(rec);
                } else {
                    r.1 = Some(rec);
                }
            }
            r
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(r: &[u8], alts: &[&[u8]]) -> PairRec {
        pair_rec(r, alts, None)
    }

    #[test]
    fn exact_match_goes_first_then_partial() {
        // FILE: C>T,G and C>G; SRC: C>G,T and C>G: the exact pair is a row of
        // its own and comes first, then the partial pair.
        let file = [rec(b"C", &[b"T", b"G"]), rec(b"C", &[b"G"])];
        let src = [rec(b"C", &[b"G", b"T"]), rec(b"C", &[b"G"])];
        assert_eq!(
            pair_rows(&file, &src),
            vec![(Some(1), Some(1)), (Some(0), Some(0))]
        );
    }

    #[test]
    fn unrelated_records_do_not_pair() {
        let file = [rec(b"A", &[b"G"])];
        let src = [rec(b"A", &[b"T"])];
        assert_eq!(
            pair_rows(&file, &src),
            vec![(Some(0), None), (None, Some(0))]
        );
    }

    #[test]
    fn types() {
        assert_eq!(alt_type(b"A", b"G"), VCF_SNP);
        assert_eq!(alt_type(b"A", b"AT"), VCF_INDEL | VCF_INS);
        assert_eq!(alt_type(b"AT", b"A"), VCF_INDEL | VCF_DEL);
        assert_eq!(alt_type(b"AC", b"GT"), VCF_MNP);
        assert_eq!(alt_type(b"AC", b"GTT"), VCF_OTHER);
        assert_eq!(alt_type(b"A", b"<DEL>"), VCF_OTHER);
        assert_eq!(alt_type(b"A", b"*"), VCF_OVERLAP);
        assert_eq!(alt_type(b"A", b"."), 0);
        assert_eq!(alt_type(b"A", b"A]1:5]"), VCF_BND);
    }
}
