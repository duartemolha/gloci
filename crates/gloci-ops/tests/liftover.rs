//! Liftover of BED and bedGraph by CrossMap 0.7.4's rules (D174, R2, R4 to R6).

use gloci_core::{ChromOrder, VcfSpan};
use gloci_io::{
    ChainIndex, Format, InvalidPolicy, IoError, ReaderOptions, RecordWriter, chrom_table_for,
};
use gloci_ops::liftover::{
    ChromStyle, LiftoverCounts, LiftoverOptions, Unmapped, liftover, source_name, target_name,
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
    run_mapped(chain, name, text, o, want_unmapped, None, None)
}

fn run_mapped(
    chain: &ChainIndex,
    name: &str,
    text: &str,
    o: &LiftoverOptions,
    want_unmapped: bool,
    chrom_map: Option<std::sync::Arc<gloci_io::ChromMap>>,
    format: Option<Format>,
) -> Result<(Vec<String>, Vec<String>, LiftoverCounts), IoError> {
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map,
        threads: NonZeroUsize::MIN,
        format,
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
        first_invalid: None,
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

#[test]
fn chrom_map_renames_before_the_lookup() {
    let map = gloci_io::ChromMap::parse("map.tsv", b"contig7\tchr1\n").unwrap();
    let (out, un, _) = run_mapped(
        &chain(),
        "m.bed",
        "contig7\t120\t180\ta\ncontig8\t120\t180\tb\n",
        &LiftoverOptions::default(),
        true,
        Some(std::sync::Arc::new(map)),
        None,
    )
    .unwrap();
    // The renamed input (`chr1`) styles the target as-input: `chr1`, not `1`.
    assert_eq!(out, ["chr1\t1120\t1180\ta"]);
    assert_eq!(un, ["contig8\t120\t180\tb"]);
}

#[test]
fn source_lookup_removes_every_chr() {
    let c = ChainIndex::parse(
        "un.chain",
        b"chain 100 Un_7 1000 + 0 100 chrUn_chr7 1000 + 0 100 1\n100\n",
    )
    .unwrap();
    assert!(source_name(&c, "chrUn_chr7").is_some());
    assert!(source_name(&c, "Un_7").is_some());
    assert!(source_name(&c, "chrUn_7").is_some());
    assert!(source_name(&c, "Un_chr7").is_some());
    assert!(source_name(&c, "chrUn_8").is_none());
}

#[test]
fn describe_reads_after_a_count() {
    use Unmapped::*;
    let want = [
        (OffChain, "on no chain"),
        (BlockSplit, "with a block not in exactly one piece"),
        (
            BlocksApart,
            "with blocks on different chromosomes or strands",
        ),
        (
            BadBlocks,
            "with blocks or a thick range that do not fit after lifting",
        ),
        (SeveralPieces, "in several pieces"),
        (Partial, "only partly on a chain (also written, shortened)"),
        (Strand, "with a strand other than +, - or ."),
        (LowRatio, "with a mapped fraction below --min-ratio"),
        (CrossChroms, "with pieces on different chromosomes"),
        (MixedStrands, "with pieces on different strands"),
    ];
    for (u, text) in want {
        assert_eq!(u.describe(), text);
    }
}

fn refused(name: &str, text: &str, format: Option<Format>) -> String {
    let o = LiftoverOptions::default();
    match run_mapped(&chain(), name, text, &o, false, None, format) {
        Err(IoError::Usage(e)) => e.0,
        other => panic!("{name}: {:?}", other.map(|_| ())),
    }
}

#[test]
fn wig_and_gff_are_refused() {
    let wig = refused("x.wig", "variableStep chrom=chr1\n121 1\n", None);
    assert!(
        wig.contains("WIG liftover") && wig.contains("gloci convert"),
        "{wig}"
    );
    let gff = refused("x.gff", "chr1\tsrc\tgene\t121\t180\t.\t+\t.\tID=g\n", None);
    assert!(gff.contains("GFF/GTF liftover"), "{gff}");
}

fn bed12(
    name: &str,
    text: &str,
    o: &LiftoverOptions,
) -> (Vec<String>, Vec<String>, LiftoverCounts) {
    run(name, text, o)
}

fn named<'a>(recs: &'a [String], name: &str) -> Vec<&'a String> {
    recs.iter()
        .filter(|l| l.split('\t').nth(3) == Some(name))
        .collect()
}

#[test]
fn bed12_lifts_as_crossmap() {
    let bed = fixture("toy-bed12.bed");
    let (out, un, counts) = bed12("toy-bed12.bed", &bed, &LiftoverOptions::default());
    assert_eq!(
        out,
        lines(
            "chr1\t1120\t1580\tt1\t0\t+\t1130\t1560\t0\t3\t60,60,60\t0,150,400
chr2\t3950\t4000\tt3\t0\t+\t3960\t3990\t0\t1\t50\t0
chr1\t1120\t1180\tt7\t0\t+\t1120\t1180\t0\t1\t60\t0
chr1\t1150\t1150\tt8\t0\t+\t1150\t1150\t0\t1\t0\t0
chr1\t1120\t1580\tt10\t0\t+\t1580\t1580\t0\t3\t60,60,60\t0,150,400
chrT\t1010\t1550\tt11\t0\t+\t1020\t1540\t0\t2\t40,40\t0,500"
        )
    );
    let want: Vec<String> = bed
        .lines()
        .filter(|l| ["t2", "t4", "t5", "t6", "t9"].contains(&l.split('\t').nth(3).unwrap()))
        .map(str::to_string)
        .collect();
    assert_eq!(un, want);
    assert_eq!(
        counts.unmapped,
        [
            (Unmapped::BadBlocks, 1),
            (Unmapped::BlockSplit, 3),
            (Unmapped::BlocksApart, 1)
        ]
        .into()
    );
}

#[test]
fn bed12_thick_ends_follow_their_blocks() {
    let bed = fixture("toy-bed12-thick.bed");
    let (out, _, _) = bed12("toy-bed12-thick.bed", &bed, &LiftoverOptions::default());
    assert_eq!(
        out,
        lines(
            "chr1\t1120\t1580\tk1\t0\t+\t1280\t1560\t0\t3\t60,60,60\t0,150,400
chr2\t3950\t4000\tk2\t0\t-\t3960\t3995\t0\t1\t50\t0
chr1\t1120\t1330\tk3\t0\t+\t1120\t1120\t0\t2\t60,60\t0,150"
        )
    );
}

#[test]
fn bed12_thick_in_an_intron_keeps_crossmap_offsets() {
    let text = "chr1\t120\t580\ti1\t0\t+\t200\t560\t0\t3\t60,60,60,\t0,100,400,\n\
chr1\t120\t580\tf1\t0\t+\t200\t250\t0\t3\t60,60,60,\t0,100,400,\n";
    let (out, _, _) = bed12("i.bed", text, &LiftoverOptions::default());
    let thick = |l: &String| {
        let f: Vec<&str> = l.split('\t').collect();
        (f[6].to_string(), f[7].to_string())
    };
    // thickStart 200 lies in an intron, so both ends keep CrossMap's offsets.
    assert_eq!(thick(&out[0]), ("1200".into(), "1560".into()), "{}", out[0]);
    assert_eq!(thick(&out[1]), ("1200".into(), "1250".into()), "{}", out[1]);
}

#[test]
fn bed12_blocks_out_of_order_do_not_lift() {
    let bed = fixture("toy-bed12-order.bed");
    let (out, un, counts) = bed12("o.bed", &bed, &LiftoverOptions::default());
    assert!(out.is_empty(), "{out:?}");
    assert_eq!(un, lines(bed.trim_end()));
    assert_eq!(counts.unmapped, [(Unmapped::BadBlocks, 1)].into());
}

#[test]
fn bed12_unknown_strand_lifts_as_plus() {
    let bed = fixture("toy-bed12-dot.bed");
    let (out, _, _) = bed12("u.bed", &bed, &LiftoverOptions::default());
    assert_eq!(
        out,
        lines("chr1\t1120\t1580\tu1\t0\t.\t1130\t1560\t0\t3\t60,60,60\t0,150,400")
    );
}

#[test]
fn bed12_extra_columns_are_kept() {
    let bed = fixture("toy-bed-extra.bed");
    let (out, _, _) = bed12("e.bed", &bed, &LiftoverOptions::default());
    assert_eq!(
        out,
        lines(
            "chr1\t1120\t1580\te1\t0\t+\t1130\t1560\t0\t3\t60,60,60\t0,150,400\tx1\tx2\tx3\tx4\tx5\tx6\tx7\tx8
chr1\t1120\t1580\te2\t0\t+\t1130\t1560\t0\t3\t60,60,60\t0,150,400\textra"
        )
    );
}

#[test]
fn ignore_blocks_lifts_columns_1_to_3_only() {
    let bed = fixture("toy-bed12.bed");
    let o = LiftoverOptions {
        ignore_blocks: true,
        ..Default::default()
    };
    let (out, _, _) = bed12("toy-bed12.bed", &bed, &o);
    assert_eq!(out.len(), 31);
    assert_eq!(
        out[0],
        "chr1\t1120\t1200\tt1\t0\t+\t130\t560\t0\t3\t60,60,60,\t0,100,400,"
    );
    assert_eq!(
        named(&out, "t3"),
        ["chr2\t3950\t4000\tt3\t0\t+\t710\t740\t0\t1\t50,\t0,"]
    );
    assert_eq!(
        out.last().unwrap(),
        "chrT\t1500\t1550\tt11\t0\t+\t20\t240\t0\t2\t40,40\t0,200"
    );
}

#[test]
fn bed12_on_an_unknown_chromosome_is_off_chain() {
    let text = "chrZ\t120\t580\tz\t0\t+\t130\t560\t0\t3\t60,60,60,\t0,100,400,\n";
    let (out, un, counts) = bed12("z.bed", text, &LiftoverOptions::default());
    assert!(out.is_empty());
    assert_eq!(un, lines(text.trim_end()));
    assert_eq!(counts.unmapped, [(Unmapped::OffChain, 1)].into());
}

#[test]
fn bed12_with_no_block_on_a_chain_is_off_chain() {
    // chr1 has no chain block in 0-100 or 600-690.
    let text = "chr1\t10\t690\tg\t0\t+\t10\t690\t0\t2\t20,20,\t0,660,\n";
    let (out, _, counts) = bed12("g.bed", text, &LiftoverOptions::default());
    assert!(out.is_empty());
    assert_eq!(counts.unmapped, [(Unmapped::OffChain, 1)].into());
}

#[test]
fn bed12_block_count_is_written_as_a_number() {
    let text = "chr1\t120\t580\tc\t0\t+\t130\t560\t0\t03\t60,60,60,\t0,100,400,\n";
    let (out, _, _) = bed12("c.bed", text, &LiftoverOptions::default());
    assert_eq!(out[0].split('\t').nth(9), Some("3"), "{}", out[0]);
}

#[test]
fn bed12_bad_thick_values_are_invalid_lines_not_panics() {
    for thick in [
        "+5",
        "-5",
        "99999999999999999999",
        "9223372036854775808",
        "581",
    ] {
        let text = format!("chr1\t120\t580\tb\t0\t+\t{thick}\t560\t0\t3\t60,60,60,\t0,100,400,\n");
        let end = format!("chr1\t120\t580\tb\t0\t+\t130\t{thick}\t0\t3\t60,60,60,\t0,100,400,\n");
        for t in [text, end] {
            match run_with(&chain(), "b.bed", &t, &LiftoverOptions::default(), false) {
                Err(IoError::Parse { msg, loc, .. }) => {
                    assert_eq!(loc, gloci_io::Loc::Line(1));
                    assert!(msg.contains("thick"), "{thick}: {msg}");
                }
                other => panic!("{thick}: {:?}", other.map(|r| r.0)),
            }
        }
    }
}
