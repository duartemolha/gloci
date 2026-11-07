#![no_main]
use gloci_io::fuzz_bodies::region_parse;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| region_parse(data));
