#![no_main]
use gloci_io::Format;
use gloci_io::fuzz_bodies::{read_all_both, read_aln_source};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // data[0] is the config byte (see gloci_io::fuzz_bodies); SAM is text, so there is no
    // wrapping to choose.
    read_all_both(data, Some(Format::Sam));
    read_aln_source(data, Format::Sam);
});
