//! `Normalizer` and `Dedup` on a small in-memory reference.
use gloci_vcf::{CheckRef, Dedup, Normalizer, Reference, RmDup, Verdict};

struct Mem(Vec<u8>);

impl Reference for Mem {
    fn get(
        &mut self,
        chrom: &[u8],
        start0: u64,
        len: u64,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        if chrom != b"c" {
            return Err("no such contig".into());
        }
        let a = (start0 as usize).min(self.0.len());
        let b = ((start0 + len) as usize).min(self.0.len());
        out.extend_from_slice(&self.0[a..b]);
        Ok(())
    }
}

fn norm(check: CheckRef, refseq: &str, line: &str) -> Result<(Verdict, String), String> {
    let mut r = Mem(refseq.as_bytes().to_vec());
    let mut out = Vec::new();
    let v = Normalizer::new(check).process(line.as_bytes(), &mut r, &mut out)?;
    Ok((v, String::from_utf8(out).unwrap()))
}

fn rec(pos: u64, r: &str, a: &str, info: &str) -> String {
    format!("c\t{pos}\t.\t{r}\t{a}\t.\t.\t{info}\tGT\t0/1\t1/1")
}

const REF: &str = "GATTTTTCACACAG";

#[test]
fn indel_in_a_homopolymer_moves_left() {
    // GATTTTTCACACAG: T at 3..7. Deleting the last T (pos 7, anchor 6).
    let (v, l) = norm(CheckRef::Exit, REF, &rec(6, "TT", "T", "DP=1")).unwrap();
    assert_eq!(
        v,
        Verdict::Keep {
            pos: 2,
            changed: true
        }
    );
    assert_eq!(l, rec(2, "AT", "A", "DP=1"));
    // An insertion of the same base.
    let (v, l) = norm(CheckRef::Exit, REF, &rec(7, "T", "TT", ".")).unwrap();
    assert_eq!(
        v,
        Verdict::Keep {
            pos: 2,
            changed: true
        }
    );
    assert_eq!(l, rec(2, "A", "AT", "."));
}

#[test]
fn trimming_and_untouched_records() {
    // Shared suffix and prefix.
    let (_, l) = norm(CheckRef::Exit, REF, &rec(8, "CAC", "CTC", ".")).unwrap();
    assert_eq!(l, rec(9, "A", "T", "."));
    // A SNP, an MNP and a record at the start of the contig stay as they are.
    for (p, r, a) in [(8, "C", "T"), (8, "CA", "GT"), (1, "GA", "G")] {
        let line = rec(p, r, a, ".");
        let (v, l) = norm(CheckRef::Exit, REF, &line).unwrap();
        assert_eq!(
            v,
            Verdict::Keep {
                pos: p,
                changed: false
            }
        );
        assert_eq!(l, line);
    }
}

#[test]
fn end_follows_the_new_ref_and_symbolic_alts_stop_a_record() {
    let (_, l) = norm(CheckRef::Exit, REF, &rec(6, "TT", "T", "END=7;DP=2")).unwrap();
    assert_eq!(l, rec(2, "AT", "A", "END=3;DP=2"));
    let line = rec(6, "TT", "T,<INS>", ".");
    assert_eq!(norm(CheckRef::Exit, REF, &line).unwrap().1, line);
    // A lone <DEL> is not moved (bcftools moves it to POS 1).
    let line = rec(6, "T", "<DEL>", ".");
    assert_eq!(norm(CheckRef::Exit, REF, &line).unwrap().1, line);
}

#[test]
fn rewritten_alleles_are_upper_case_and_reference_bases_keep_their_case() {
    let (_, l) = norm(CheckRef::Exit, "GAttttCA", &rec(5, "tt", "t", ".")).unwrap();
    assert_eq!(l, rec(2, "At", "A", "."));
    let (_, l) = norm(CheckRef::Exit, "GAcgtt", &rec(4, "gt", "cT", ".")).unwrap();
    assert_eq!(l, rec(4, "G", "C", "."));
}

#[test]
fn check_ref_modes() {
    let line = rec(8, "G", "T", ".");
    assert!(
        norm(CheckRef::Exit, REF, &line)
            .unwrap_err()
            .contains("c:8")
    );
    assert_eq!(norm(CheckRef::Warn, REF, &line).unwrap().1, line);
    assert_eq!(
        norm(CheckRef::Exclude, REF, &line).unwrap().0,
        Verdict::Drop
    );
    assert_eq!(
        norm(CheckRef::Set, REF, &line).unwrap().1,
        rec(8, "C", "T", ".")
    );
    // An ALT that is the reference base is swapped into REF; GT is recoded.
    let (_, l) = norm(CheckRef::Set, REF, &rec(8, "G", "C", ".")).unwrap();
    assert_eq!(l, "c\t8\t.\tC\tG\t.\t.\t.\tGT\t1/0\t0/0");
    // Non-ACGTN bases.
    assert!(
        norm(CheckRef::Exit, "GAN", &rec(2, "R", "A", "."))
            .unwrap_err()
            .contains("ACGTN")
    );
    assert!(
        norm(CheckRef::Exit, REF, &rec(8, "C", "R", "."))
            .unwrap_err()
            .contains("ACGTN")
    );
    // ALT equal to REF.
    assert!(
        norm(CheckRef::Exit, REF, &rec(8, "C", "c", "."))
            .unwrap_err()
            .contains("same as REF")
    );
    // Past the end of the contig.
    assert!(norm(CheckRef::Exit, REF, &rec(40, "A", "C", ".")).is_err());
    assert!(
        norm(CheckRef::Set, REF, &rec(40, "A", "C", "."))
            .unwrap_err()
            .contains("past the end")
    );
}

#[test]
fn bad_records_and_contigs_are_errors() {
    assert!(
        norm(CheckRef::Exit, REF, "c\t5\t.\tA")
            .unwrap_err()
            .contains("columns")
    );
    assert!(
        norm(CheckRef::Exit, REF, "c\tx\t.\tA\tC\t.\t.\t.")
            .unwrap_err()
            .contains("POS")
    );
    assert!(
        norm(
            CheckRef::Exit,
            REF,
            &rec(1, "G", "A", ".").replace('c', "z")
        )
        .is_err()
    );
}

#[test]
fn rm_dup_rules() {
    let keep = |mode, lines: &[(u64, &str, &str)]| -> Vec<bool> {
        let mut d = Dedup::new(mode);
        lines
            .iter()
            .map(|&(p, r, a)| d.keep(rec(p, r, a, ".").as_bytes()))
            .collect()
    };
    let set = [
        (5, "A", "C"),
        (5, "A", "G"),
        (5, "AT", "A"),
        (5, "AC", "A"),
        (5, "A", "a"),
        (6, "A", "C"),
    ];
    assert_eq!(
        keep(RmDup::Snps, &set),
        [true, false, true, true, true, true]
    );
    assert_eq!(
        keep(RmDup::Indels, &set),
        [true, true, true, false, true, true]
    );
    assert_eq!(
        keep(RmDup::Both, &set),
        [true, false, true, false, true, true]
    );
    assert_eq!(
        keep(RmDup::All, &set),
        [true, false, false, false, false, true]
    );
    assert_eq!(
        keep(
            RmDup::Exact,
            &[
                (5, "A", "C,G"),
                (5, "A", "G,C"),
                (5, "A", "G"),
                (5, "G", "A")
            ]
        ),
        [true, false, true, true]
    );
}

#[test]
fn any_symbolic_alt_leaves_the_record_alone() {
    for alt in ["T,<DEL>", "<DEL>", "<DEL:ME>,T"] {
        let line = rec(6, "TT", alt, ".");
        assert_eq!(norm(CheckRef::Exit, REF, &line).unwrap().1, line, "{alt}");
    }
}

#[test]
fn set_ref_on_a_non_acgtn_ref_does_not_swap() {
    // The FASTA (cleaned) has CNT where the VCF has CRT.
    let (_, l) = norm(CheckRef::Set, "GCNT", &rec(2, "CRT", "C", ".")).unwrap();
    assert_eq!(l, rec(2, "CNT", "C", "."));
}

#[test]
fn swap_recounts_ac_of_the_swapped_alt_only() {
    let (_, l) = norm(CheckRef::Set, REF, &rec(8, "G", "C,A", "AC=9,4")).unwrap();
    // GT 0/1 1/1 -> 1/0 0/0 with C now REF: ALT 1 (G) is counted in the new GT.
    assert_eq!(l, "c\t8\t.\tC\tG,A\t.\t.\tAC=1,4\tGT\t1/0\t0/0");
}

#[test]
fn rm_dup_other_length_changes_are_neither_snp_nor_indel() {
    let mut d = Dedup::new(RmDup::Both);
    for (r, a) in [("TAT", "CG"), ("CC", "ATGG"), ("T", "CGG"), ("G", "CC")] {
        assert!(d.keep(rec(5, r, a, ".").as_bytes()), "{r}>{a}");
    }
    assert!(d.keep(rec(5, "TA", "T", ".").as_bytes()));
    assert!(!d.keep(rec(5, "TA", "T", ".").as_bytes()));
}
