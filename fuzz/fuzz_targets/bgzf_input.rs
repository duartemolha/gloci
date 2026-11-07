#![no_main]
use gloci_io::fuzz_bodies::drain_input;
use gloci_testkit::fixtures::bgzip;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // data[0] is the config byte; bit 0 set wraps the rest in BGZF.
    match data.split_first() {
        Some((&c, rest)) if c & 1 == 1 => {
            let mut v = vec![c];
            v.extend(bgzip(rest));
            drain_input(&v);
        }
        _ => drain_input(data),
    }
});
