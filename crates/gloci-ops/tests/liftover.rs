//! Liftover of BED and bedGraph by CrossMap 0.7.4's rules (D174, R2, R4 to R6).

use gloci_core::{ChromOrder, VcfSpan};
use gloci_io::{ChainIndex, InvalidPolicy, IoError, ReaderOptions, RecordWriter, chrom_table_for};
use gloci_ops::liftover::{
    ChromStyle, LiftoverCounts, LiftoverOptions, Unmapped, liftover, target_name,
};
use gloci_testkit::mem::mem_reader;
use std::io;
use std::num::NonZeroUsize;
use std::path::PathBuf;

#[derive(Default)]
struct Sink(Vec<String>);

impl RecordWriter for Sink {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        self.0.push(String::from_utf8(raw.to_vec()).unwrap());
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }
}

fn fixture(name: &str) -> String {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/liftover")
        .join(name);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn chain() -> ChainIndex {
    ChainIndex::parse("toy.chain", fixture("toy.chain").as_bytes()).unwrap()
}

fn run_with(
    chain: &ChainIndex,
    name: &str,
    text: &str,
    o: &LiftoverOptions,
    want_unmapped: bool,
) -> Result<(Vec<String>, Vec<String>, LiftoverCounts), IoError> {
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let reader = mem_reader(name, text, &ro);
    let mut chroms = chrom_table_for(None, &[reader.as_ref()], ChromOrder::Karyotype).unwrap();
    let (mut out, mut un) = (Sink::default(), Sink::default());
    let counts = liftover(
        reader,
        &mut chroms,
        chain,
        o,
        &mut out,
        want_unmapped.then_some(&mut un as &mut dyn RecordWriter),
    )?;
    Ok((out.0, un.0, counts))
}

fn run(name: &str, text: &str, o: &LiftoverOptions) -> (Vec<String>, Vec<String>, LiftoverCounts) {
    run_with(&chain(), name, text, o, true).unwrap()
}

fn lines(s: &str) -> Vec<String> {
    s.lines().map(str::to_string).collect()
}

const TOY_BED_OUT: &str = "chr1\t1120\t1180\tr1\t0\t+
chr1\t1150\t1200\tr2\t0\t+
chr1\t1250\t1300\tr2\t0\t+
chr1\t1330\t1350\tr3\t0\t+
chr1\t1350\t1380\tr3\t0\t+
chr1\t1100\t1150\tr5\t0\t+
chr2\t3940\t3980\tr6\t0\t-
chr2\t3920\t3940\tr7\t0\t+
chr2\t3880\t3900\tr7\t0\t+
chr1\t1420\t1440\tr8\t0\t+
chr3\t20\t40\tr8\t0\t+
chr4\t30\t50\tr8\t0\t+
chr1\t1150\t1150\tr9\t0\t+
chr5\t510\t520\tr11\t0\t+
5\t530\t540\tr12\t0\t-
1\t1120\t1180\tr13
chr1\t1120\t1130
chr2\t3940\t3980\tr16\t0\t.
chr1\t1550\t1600\tr17\t0\t+
chr2\t3950\t4000\tr17\t0\t-
chr1\t1100\t1101\tr19\t0\t+
chr1\t1100\t1200\tr21\t0\t+
chr1\t1250\t1350\tr21\t0\t+
chr1\t1350\t1600\tr21\t0\t+
chr3\t0\t100\tr21\t0\t+
chr4\t10\t60\tr21\t0\t+
chr2\t3920\t4000\tr21\t0\t-
chr2\t3800\t3900\tr21\t0\t-
chrU\t50\t100\tr22\t0\t+
chrU\t850\t900\tr22\t0\t-
chr1\t1160\t1170\tr23\t0\t+\tx\t1
chr1\t1420\t1420\tr24\t0\t+
chr3\t20\t20\tr24\t0\t+
chr4\t30\t30\tr24\t0\t+
chr1\t1180\t1200\tr25\t0\t+
";

fn toy_counts() -> LiftoverCounts {
    LiftoverCounts {
        records: 27,
        written: 35,
        split: 8,
        shortened: 6,
        invalid: 0,
        unmapped: [(Unmapped::OffChain, 7)].into(),
    }
}

#[test]
fn bed_lifts_as_crossmap() {
    let bed = fixture("toy.bed");
    let (out, un, counts) = run("toy.bed", &bed, &LiftoverOptions::default());
    assert_eq!(out, lines(TOY_BED_OUT));
    let want: Vec<String> = bed
        .lines()
        .filter(|l| {
            ["r4", "r10", "r15", "r18", "r20", "r26", "r27"]
                .contains(&l.split('\t').nth(3).unwrap_or(""))
        })
        .map(str::to_string)
        .collect();
    assert_eq!(want.len(), 7);
    assert_eq!(un, want);
    assert_eq!(counts, toy_counts());
    assert_eq!(counts.unmapped_total(), 7);
}

#[test]
fn unmapped_without_writer_is_counted() {
    let bed = fixture("toy.bed");
    let (out, un, counts) = run_with(
        &chain(),
        "toy.bed",
        &bed,
        &LiftoverOptions::default(),
        false,
    )
    .unwrap();
    assert_eq!(out, lines(TOY_BED_OUT));
    assert!(un.is_empty());
    assert_eq!(counts, toy_counts());
}

#[test]
fn bedgraph_lifts_as_crossmap() {
    let (out, un, _) = run(
        "toy.bedGraph",
        &fixture("toy.bedGraph"),
        &LiftoverOptions::default(),
    );
    assert_eq!(
        out,
        [
            "chr1\t1120\t1180\t1.5",
            "chr1\t1150\t1200\t2",
            "chr1\t1250\t1300\t2",
            "chr1\t1420\t1440\t3",
            "chr3\t20\t40\t3",
            "chr4\t30\t50\t3",
            "chr2\t3940\t3980\t4",
        ]
    );
    assert_eq!(un, ["chr1\t310\t340\t9"]);
}

#[test]
fn chrom_names_follow_crossmap() {
    let names = fixture("toy-names.bed");
    let want = [
        (ChromStyle::AsInput, ["chr5", "5", "chr1", "1", "chrT", "T"]),
        (ChromStyle::Short, ["5", "5", "1", "1", "T", "T"]),
        (
            ChromStyle::Long,
            ["chr5", "chr5", "chr1", "chr1", "chrT", "chrT"],
        ),
        (
            ChromStyle::Chain,
            ["5", "5", "chr1", "chr1", "chrT", "chrT"],
        ),
    ];
    for (style, chroms) in want {
        let o = LiftoverOptions {
            chrom_style: style,
            ..Default::default()
        };
        let (out, un, _) = run("toy-names.bed", &names, &o);
        assert!(un.is_empty(), "{style:?}: {un:?}");
        let got: Vec<(&str, &str)> = out
            .iter()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0], f[3])
            })
            .collect();
        let want: Vec<(&str, &str)> = chroms
            .iter()
            .enumerate()
            .map(|(i, c)| (*c, ["n1", "n2", "n3", "n4", "n5", "n6"][i]))
            .collect();
        assert_eq!(got, want, "{style:?}");
    }
    assert_eq!(target_name(ChromStyle::AsInput, "1", "chrUn_chr7"), "Un_7");
    assert_eq!(target_name(ChromStyle::AsInput, "chrX", "MT"), "chrMT");
    assert_eq!(target_name(ChromStyle::Short, "x", "7"), "7");
}

#[test]
fn chrm_is_not_mt() {
    let c = ChainIndex::parse(
        "mt.chain",
        b"chain 100 MT 1000 + 0 100 MT 1000 + 0 100 1\n100\n",
    )
    .unwrap();
    let (out, un, counts) = run_with(
        &c,
        "m.bed",
        "chrM\t10\t20\n",
        &LiftoverOptions::default(),
        true,
    )
    .unwrap();
    assert!(out.is_empty());
    assert_eq!(un, ["chrM\t10\t20"]);
    assert_eq!(counts.unmapped[&Unmapped::OffChain], 1);
}

#[test]
fn strand_column_is_column_six() {
    let (out, _, _) = run(
        "toy-strand.bed",
        &fixture("toy-strand.bed"),
        &LiftoverOptions::default(),
    );
    assert_eq!(out, ["chr2\t3940\t3980\tn\t0\t+\t-", "chr2\t3940\t3980\t-"]);
}

#[test]
fn fields_split_on_tabs() {
    let (out, _, _) = run(
        "toy-spaces.bed",
        &fixture("toy-spaces.bed"),
        &LiftoverOptions::default(),
    );
    assert_eq!(out, ["chr1\t1120\t1180\tname with space\t0\t+"]);
}

#[test]
fn other_formats_are_refused() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t120\t.\tA\tC\t.\t.\t.\n";
    match run_with(&chain(), "x.vcf", vcf, &LiftoverOptions::default(), false) {
        Err(IoError::Usage(e)) => assert!(e.0.contains("VCF liftover"), "{}", e.0),
        other => panic!("{:?}", other.map(|_| ())),
    }
}
