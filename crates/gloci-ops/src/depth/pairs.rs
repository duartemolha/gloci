//! mosdepth 0.3.10's rule for the two reads of a proper pair (D173): the bases of the second
//! read that the first read's aligned segments already cover are not counted again. The table
//! holds reads whose mate is still ahead, as [`MateTable`](super::MateTable) does for D139.

use gloci_core::Interval;
use gloci_io::AlnRef;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// The aligned segments of `r` from `pos`: maximal runs of M, = and X positions (I, S, H and P do
/// not end a run, D and N do), cut at `end`. Written to `out`, which is cleared first.
pub fn aligned_segments(r: &AlnRef<'_>, pos: u64, end: u64, out: &mut Vec<Interval>) {
    out.clear();
    let mut rpos = pos;
    // Where the open run started, if the last operation that moved `rpos` was aligned.
    let mut open: Option<u64> = None;
    let close = |out: &mut Vec<Interval>, start: u64, stop: u64| {
        let stop = stop.min(end);
        if start < stop {
            out.push(Interval { start, end: stop });
        }
    };
    r.for_each_op(&mut |op, n| {
        let n = u64::from(n);
        match op {
            b'M' | b'=' | b'X' => {
                open.get_or_insert(rpos);
                rpos = rpos.saturating_add(n);
            }
            b'D' | b'N' => {
                if let Some(start) = open.take() {
                    close(out, start, rpos);
                }
                rpos = rpos.saturating_add(n);
            }
            _ => {}
        }
    });
    if let Some(start) = open {
        close(out, start, rpos);
    }
}

/// mosdepth 0.3.10's pairs (D173): reads waiting for their mate, with their aligned segments.
/// A min-heap on the mate position drops entries nobody can claim any more; an entry taken or
/// replaced leaves a stale heap item, skipped when it surfaces.
#[derive(Default)]
pub struct PairTable {
    pending: HashMap<Box<[u8]>, (u64, Vec<Interval>)>,
    by_mpos: BinaryHeap<Reverse<(u64, Box<[u8]>)>>,
    peak: usize,
}

impl PairTable {
    pub fn new() -> PairTable {
        PairTable::default()
    }

    /// The segments `r` must not count again, written to `out` (cleared first; empty when none).
    /// Only a read with PROPER_PAIR (0x2) set and SUPPLEMENTARY (0x800) clear takes part. It is
    /// held (replacing an entry of its name) when its mate is on its reference, `end > mpos`, and
    /// `pos < mpos` or (`pos == mpos` and no entry of its name is held). Otherwise it takes the
    /// entry of its name, if any, and `out` gets that entry's segments.
    pub fn shared_with_mate(
        &mut self,
        r: &AlnRef<'_>,
        pos: u64,
        end: u64,
        out: &mut Vec<Interval>,
    ) {
        out.clear();
        let flag = r.flag();
        if flag & FLAG_PROPER_PAIR == 0 || flag & FLAG_SUPPLEMENTARY != 0 {
            return;
        }
        let mpos = match (r.mate_ref_id(), r.mate_pos()) {
            (Some(mref), Some(mpos)) if Some(mref) == r.ref_id() => Some(mpos),
            _ => None,
        };
        if let Some(mpos) = mpos
            && end > mpos
            && (pos < mpos || (pos == mpos && !self.pending.contains_key(r.name())))
        {
            let mut segments = Vec::new();
            aligned_segments(r, pos, end, &mut segments);
            let name: Box<[u8]> = r.name().into();
            self.by_mpos.push(Reverse((mpos, name.clone())));
            self.pending.insert(name, (mpos, segments));
            self.peak = self.peak.max(self.pending.len());
        } else if let Some((_, segments)) = self.pending.remove(r.name()) {
            *out = segments;
        }
    }

    /// Drop the entries whose mate position is below `pos`.
    pub fn advance(&mut self, pos: u64) {
        while let Some(Reverse((mpos, _))) = self.by_mpos.peek() {
            if *mpos >= pos {
                break;
            }
            let Some(Reverse((mpos, name))) = self.by_mpos.pop() else {
                break;
            };
            // A later read of the same name may have been held with its own mate position.
            if self.pending.get(&name).is_some_and(|&(m, _)| m == mpos) {
                self.pending.remove(&name);
            }
        }
        if self.pending.is_empty() {
            self.by_mpos.clear();
        }
    }

    /// Forget everything: a new reference starts.
    pub fn reset(&mut self) {
        self.pending.clear();
        self.by_mpos.clear();
    }

    /// Reads waiting for a mate.
    pub fn len(&self) -> usize {
        self.pending.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }

    /// The most reads that waited at once.
    pub fn peak(&self) -> usize {
        self.peak
    }
}

#[cfg(test)]
mod tests {
    use super::super::rows::{DepthSink, ReadFilter, ReadLength};
    use super::super::test_sink::{Recorder, sam_source};
    use super::super::{AllPositions, CountRules, DepthOptions, DepthRing, ModeRules, run_sources};
    use super::{Interval, aligned_segments};
    use gloci_core::VcfSpan;
    use gloci_io::{Accept, AlignOptions, AlnRef, InvalidPolicy, ReaderOptions, open_aln_source};
    use std::io;
    use std::num::NonZeroUsize;

    fn ro() -> ReaderOptions {
        ReaderOptions {
            accept: Accept {
                alignments: Some(AlignOptions {
                    exclude_flags: 0,
                    ..AlignOptions::default()
                }),
                ..Accept::default()
            },
            chrom_map: None,
            threads: NonZeroUsize::MIN,
            format: None,
            vcf_span: VcfSpan::Ref,
            invalid: InvalidPolicy::Fail,
        }
    }

    fn opts(all: AllPositions) -> DepthOptions {
        DepthOptions {
            rules: CountRules::default(),
            min_read_length: 0,
            all,
            regions: None,
            skip_mate_overlap: false,
            column_names: false,
        }
    }

    fn mosdepth_mode(exclude_flags: u16) -> ModeRules {
        ModeRules {
            read_length: ReadLength::SeqLessClips,
            filter: Some(ReadFilter {
                exclude_flags,
                min_mapq: 0,
            }),
            shared_bases: true,
        }
    }

    fn iv(start: u64, end: u64) -> Interval {
        Interval { start, end }
    }

    /// A sink that runs `f` on every read it is given.
    struct PerRead<F: FnMut(&AlnRef<'_>, Option<Interval>)>(F);

    impl<F: FnMut(&AlnRef<'_>, Option<Interval>)> DepthSink for PerRead<F> {
        fn position(&mut self, _: usize, _: u64, _: &[u32]) -> io::Result<()> {
            Ok(())
        }
        fn zeros(&mut self, _: usize, _: u64, _: u64) -> io::Result<()> {
            Ok(())
        }
        fn read(&mut self, _: usize, _: usize, a: &AlnRef<'_>, span: Option<Interval>) {
            (self.0)(a, span);
        }
    }

    /// Give every read of the one-reference SAM made of `reads` (`name flag pos cigar`) to `f`.
    fn each_read(reads: &[&str], f: impl FnMut(&AlnRef<'_>, Option<Interval>)) {
        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n");
        for r in reads {
            let p: Vec<&str> = r.split(' ').collect();
            text.push_str(&format!(
                "{}\t{}\tc1\t{}\t60\t{}\t*\t0\t0\t*\t*\n",
                p[0], p[1], p[2], p[3]
            ));
        }
        let src = sam_source(&text, &ro());
        run_sources(
            vec![src],
            &["test.sam".to_string()],
            &opts(AllPositions::Covered),
            ModeRules::depth(),
            &ro(),
            |_| Ok(PerRead(f)),
        )
        .unwrap();
    }

    #[test]
    fn aligned_segments_merge_across_insertions() {
        let cases: [(&str, u64, u64, Vec<Interval>); 4] = [
            ("5S20M2I20M", 500, 1000, vec![iv(500, 540)]),
            ("20M10D20M", 300, 1000, vec![iv(300, 320), iv(330, 350)]),
            ("10M50N10M", 400, 1000, vec![iv(400, 410), iv(460, 470)]),
            ("50M", 0, 30, vec![iv(0, 30)]),
        ];
        for (cigar, pos, end, expect) in cases {
            let read = format!("r 0 {} {cigar}", pos + 1);
            let mut got = vec![iv(9, 9)];
            each_read(&[&read], |a, _| aligned_segments(a, pos, end, &mut got));
            assert_eq!(got, expect, "{cigar} at {pos} cut at {end}");
        }
    }

    #[test]
    fn add_except_skips_the_mate_bases() {
        let mut depths = Vec::new();
        let mut ring = DepthRing::new();
        let except = [iv(300, 320), iv(330, 350)];
        each_read(&["r 0 321 40M"], |a, span| {
            ring.add_except(a, span.unwrap(), &CountRules::default(), &except);
            depths = (300..360).map(|p| ring.take(p)).collect();
        });
        let run = |from: usize, to: usize| depths[from - 300..to - 300].to_vec();
        assert_eq!(run(320, 330), [1; 10]);
        assert_eq!(run(330, 350), [0; 20]);
        assert_eq!(run(350, 360), [1; 10]);
    }

    /// The runs of the reference `rid` as `start-end:depth`.
    fn runs_of(rec: &Recorder, rid: usize) -> String {
        rec.runs
            .iter()
            .filter(|r| r.0 == rid)
            .map(|r| format!("{}-{}:{}", r.1, r.2, r.3))
            .collect::<Vec<_>>()
            .join(" ")
    }

    fn fixture_runs(file: &str, exclude_flags: u16) -> Recorder {
        let path = format!(
            "{}/../../xtask/diff/inputs/samtools/{file}",
            env!("CARGO_MANIFEST_DIR")
        );
        let src = open_aln_source(&path, &ro()).unwrap();
        let (_, rec) = run_sources(
            vec![src],
            &[path],
            &opts(AllPositions::AllRefs),
            mosdepth_mode(exclude_flags),
            &ro(),
            |_| Ok(Recorder::default()),
        )
        .unwrap();
        rec
    }

    #[test]
    fn pair_rule_matches_mosdepth_per_base() {
        let rec = fixture_runs("mosdepth.bam", 1796);
        // mosdepth 0.3.10's per-base output for the same file and filter.
        assert_eq!(
            runs_of(&rec, 0),
            "0-100:0 100-109:1 109-139:2 139-180:1 180-200:0 200-260:1 260-300:0 300-360:1 \
360-400:0 400-470:1 470-500:0 500-550:1 550-600:0 600-610:1 610-630:2 630-640:1 640-700:0 \
700-730:1 730-750:2 750-770:1 770-800:0 800-840:1 840-850:0 850-910:1 910-1000:0"
        );
        assert_eq!(runs_of(&rec, 1), "0-300:0");
        assert_eq!(runs_of(&rec, 2), "0-200:0");
        assert_eq!(
            runs_of(&rec, 3),
            "0-5000:0 5000-5001:320 5001-6000:0 6000-6010:1 6010-6020:2 6020-6030:1 6030-20000:0"
        );
    }

    #[test]
    fn pair_rule_mates_bam() {
        // mosdepth's, except where a mate lies inside the other with one-operation CIGARs
        // (template c): mosdepth subtracts up to the first mate's end, gloci the segments.
        let rec = fixture_runs("mates.bam", 1796);
        assert_eq!(
            runs_of(&rec, 0),
            "0-100:0 100-250:1 250-300:0 300-380:1 380-400:0 400-600:1 600-700:0 700-750:1 \
750-900:0 900-950:1 950-1000:0 1000-1040:1 1040-1100:0 1100-1130:1 1130-1140:0 1140-1145:1 \
1145-1200:0 1200-1250:1 1250-1400:0 1400-1500:1 1500-1519:2 1519-1530:3 1530-1539:2 \
1539-1550:1 1550-1600:0 1600-1670:1 1670-1800:0 1800-1840:1 1840-1900:0 1900-1990:1 \
1990-2000:0"
        );
        assert_eq!(
            runs_of(&rec, 1),
            "0-4:0 4-44:1 44-100:0 100-170:1 170-300:0 300-350:1 350-500:0"
        );
    }

    #[test]
    fn stale_pair_entry_is_dropped() {
        // mosdepth keeps the entry of `s` to the chromosome's end, so the secondary alignment
        // at 301 takes it and loses [150,300) of depth; gloci dropped it at `o`.
        let rec = fixture_runs("mosdepth-stale.bam", 1540);
        assert_eq!(
            runs_of(&rec, 0),
            "0-100:0 100-130:1 130-140:2 140-150:1 150-300:0 300-350:1 350-1000:0"
        );
    }

    #[test]
    fn pair_table_bounded() {
        // 100,000 proper pairs with 50 bases between the mates, then 1,000 whose second mate
        // is a duplicate that the filter drops: no entry outlives the next read.
        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:30000000\n");
        for i in 0..101_000u64 {
            let (first, second) = (200 * i + 1, 200 * i + 51);
            let flag2 = if i < 100_000 { 147 } else { 1171 };
            text.push_str(&format!(
                "p{i}\t99\tc1\t{first}\t60\t150M\t=\t{second}\t200\t*\t*\n\
p{i}\t{flag2}\tc1\t{second}\t60\t150M\t=\t{first}\t-200\t*\t*\n"
            ));
        }
        let src = sam_source(&text, &ro());
        let (summary, _) = run_sources(
            vec![src],
            &["test.sam".to_string()],
            &opts(AllPositions::Covered),
            mosdepth_mode(1796),
            &ro(),
            |_| Ok(Recorder::default()),
        )
        .unwrap();
        assert!(summary.max_pending_mates >= 1);
        assert!(
            summary.max_pending_mates <= 2,
            "peak {}",
            summary.max_pending_mates
        );
    }
}
