//! The `vcf query` format language: syntax and names (D172).

use gloci_io::VcfHeader;
use gloci_vcf::QueryFormat;
use gloci_vcf::query::{Field, Piece, Tag, TagType};

fn header() -> VcfHeader {
    let vcf = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../xtask/diff/inputs/vcf/split/multi.vcf"
    ))
    .unwrap();
    let end: usize = vcf
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(|l| l.len() + 1)
        .sum();
    VcfHeader::parse(vcf.as_bytes()[..end].to_vec()).unwrap()
}

fn q(fmt: &str) -> Vec<Piece> {
    QueryFormat::parse(fmt)
        .unwrap()
        .resolve(&header(), false)
        .unwrap()
        .pieces
}

fn syntax_err(fmt: &str) -> String {
    QueryFormat::parse(fmt).unwrap_err().0
}

fn name_err(fmt: &str) -> String {
    QueryFormat::parse(fmt)
        .unwrap()
        .resolve(&header(), false)
        .unwrap_err()
        .0
}

fn t(s: &str) -> Piece {
    Piece::Text(s.into())
}

fn tag(n: &str, ty: TagType, sub: Option<usize>) -> Tag {
    Tag {
        name: n.into(),
        ty,
        sub,
    }
}

#[test]
fn parses_fields_text_blocks_and_subscripts() {
    assert_eq!(
        q(r"%CHROM:%POS0\t%INFO/DP %AC{1}[ %SAMPLE=%GT/%TGT %AD{0}]\n"),
        vec![
            Piece::Field(Field::Chrom),
            t(":"),
            Piece::Field(Field::Pos0),
            t("\t"),
            Piece::Field(Field::InfoTag(tag("DP", TagType::Integer, None))),
            t(" "),
            Piece::Field(Field::InfoTag(tag("AC", TagType::Integer, Some(1)))),
            Piece::Block(vec![
                t(" "),
                Piece::Field(Field::Sample),
                t("="),
                Piece::Field(Field::Gt),
                t("/"),
                Piece::Field(Field::Tgt),
                t(" "),
                Piece::Field(Field::FormatTag(tag("AD", TagType::Integer, Some(0)))),
            ]),
            t("\n"),
        ]
    );
    assert_eq!(
        q(r"%DP[%DP]\n"),
        vec![
            Piece::Field(Field::InfoTag(tag("DP", TagType::Integer, None))),
            Piece::Block(vec![Piece::Field(Field::FormatTag(tag(
                "DP",
                TagType::Integer,
                None
            )))]),
            t("\n"),
        ]
    );
    assert_eq!(
        q(r"[%XS %DB %AF %QUAL %INFO/DP]\n"),
        vec![
            Piece::Block(vec![
                Piece::Field(Field::InfoTag(tag("XS", TagType::String, None))),
                t(" "),
                Piece::Field(Field::InfoTag(tag("DB", TagType::Flag, None))),
                t(" "),
                Piece::Field(Field::InfoTag(tag("AF", TagType::Float, None))),
                t(" "),
                Piece::Field(Field::Qual),
                t(" "),
                Piece::Field(Field::InfoTag(tag("DP", TagType::Integer, None))),
            ]),
            t("\n"),
        ]
    );
    assert_eq!(
        q(r"%ALT{2}%FIRST_ALT%END0%TYPE%INFO%LINE"),
        vec![
            Piece::Field(Field::Alt(Some(2))),
            Piece::Field(Field::FirstAlt),
            Piece::Field(Field::End0),
            Piece::Field(Field::Type),
            Piece::Field(Field::Info),
            Piece::Field(Field::Line),
            t("\n"),
        ]
    );
    assert_eq!(q(r"a\\b%POS\n")[0], t("a\\b"));
}

#[test]
fn appends_a_newline_when_the_format_has_none() {
    assert_eq!(q("%POS"), vec![Piece::Field(Field::Pos), t("\n")]);
    assert_eq!(q("[%SAMPLE]").last(), Some(&t("\n")));
    assert_eq!(q(r"%POS\n%ID").len(), 3);
    assert_eq!(q(r"[%SAMPLE\n]").len(), 1);
}

#[test]
fn syntax_errors_name_the_construct_and_column() {
    let cases: &[(&str, &[&str])] = &[
        ("%POS %", &["column 6", "followed by a field name"]),
        ("%%", &["column 1"]),
        (r"%POS\x", &["\\x", "\\t, \\n or \\\\", "column 5"]),
        (r"%POS\", &["column 5", "backslash"]),
        ("%POS]", &["] without [", "column 5"]),
        ("[[%GT]]", &["nested", "column 2"]),
        ("[%GT", &["missing ]", "column 1"]),
        ("%AF{1", &["{N}", "column 4"]),
        ("%AC{-1}", &["{N}", "column 4"]),
        ("%AC{99999999999999999999999}", &["too large", "column 4"]),
    ];
    for (fmt, wants) in cases {
        let e = syntax_err(fmt);
        assert!(e.contains("invalid query format"), "{fmt}: {e}");
        for w in *wants {
            assert!(e.contains(w), "{fmt}: {e} lacks {w}");
        }
    }
}

#[test]
fn misplaced_and_unknown_fields_are_usage_errors() {
    let cases: &[(&str, &[&str])] = &[
        ("%CHROM %SAMPLE", &["%SAMPLE", "inside [ ]", "column 8"]),
        ("%GT", &["inside [ ]"]),
        ("%TGT", &["inside [ ]"]),
        ("%AD", &["FORMAT tag", "inside [ ]"]),
        ("[%LINE]", &["%LINE", "outside [ ]"]),
        ("%FMT/DP", &["[%DP]"]),
        ("[%FORMAT/DP]", &["[%DP]"]),
        ("%XX", &["XX", "not declared"]),
        ("%INFO/XX", &["INFO/XX", "not declared"]),
        ("[%IUPACGT]", &["IUPACGT"]),
        ("%POS{1}", &["%POS takes no subscript"]),
        ("[%GT{0}]", &["%GT takes no subscript"]),
        ("%INFO/", &["%INFO/ needs a tag name"]),
    ];
    for (fmt, wants) in cases {
        let e = name_err(fmt);
        assert!(e.contains("invalid query format"), "{fmt}: {e}");
        assert!(e.ends_with(')'), "{fmt}: {e}");
        for w in *wants {
            assert!(e.contains(w), "{fmt}: {e} lacks {w}");
        }
    }
}

#[test]
fn header_flag_refuses_line() {
    let f = QueryFormat::parse("%LINE").unwrap();
    let e = f.resolve(&header(), true).unwrap_err().0;
    assert!(e.contains("-H") && e.contains("%LINE"), "{e}");
    assert!(f.resolve(&header(), false).is_ok());
}
