#![no_main]
use gloci_io::fuzz_bodies::read_bedpe;

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    read_bedpe(data);
});
