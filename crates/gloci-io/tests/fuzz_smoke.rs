//! Stable stand-in for the cargo-fuzz targets: arbitrary bytes (raw and
//! BGZF-wrapped) through the same bodies must never panic.

use gloci_io::Format;
use gloci_io::fuzz_bodies::{drain_input, read_all_both, region_parse};
use gloci_testkit::fixtures::bgzip;
use proptest::prelude::*;

fn bytes() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(any::<u8>(), 0..512)
}

/// Text-ish input, so the parsers get past the first few bytes.
fn texty() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop::sample::select(b"chr1\t2 345-:#.,;=_\nACGT<>*/|\r".to_vec()),
        0..400,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn readers_never_panic(d in prop_oneof![bytes(), texty()]) {
        for f in [Format::Bed, Format::Vcf, Format::Gff3, Format::Bcf] {
            read_all_both(&d, Some(f));
            read_all_both(&bgzip(&d), Some(f));
        }
        read_all_both(&d, None);
        read_all_both(&bgzip(&d), None);
    }

    #[test]
    fn input_never_panics(d in prop_oneof![bytes(), texty()]) {
        drain_input(&d);
        drain_input(&bgzip(&d));
    }

    #[test]
    fn region_parse_never_panics(d in prop_oneof![bytes(), texty()]) {
        region_parse(&d);
    }
}
