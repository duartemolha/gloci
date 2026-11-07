#![no_main]
use gloci_io::fuzz_bodies::drain_input;
use gloci_testkit::fixtures::bgzip;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if data.first().is_some_and(|b| b & 1 == 0) {
        drain_input(&bgzip(data));
    } else {
        drain_input(data);
    }
});
