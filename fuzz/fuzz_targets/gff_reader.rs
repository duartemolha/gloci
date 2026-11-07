#![no_main]
use gloci_io::Format;
use gloci_io::fuzz_bodies::read_all_both;

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    read_all_both(data, Some(Format::Gff3));
});
