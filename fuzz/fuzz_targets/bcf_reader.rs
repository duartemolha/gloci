#![no_main]
use gloci_io::Format;
use gloci_io::fuzz_bodies::read_all_both;
use gloci_testkit::fixtures::bgzip;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Half the inputs get a BGZF wrapper so the fuzzer gets past that layer.
    if data.first().is_some_and(|b| b & 1 == 0) {
        read_all_both(&bgzip(data), Some(Format::Bcf));
    } else {
        read_all_both(data, Some(Format::Bcf));
    }
});
