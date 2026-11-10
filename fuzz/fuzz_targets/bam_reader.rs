#![no_main]
use gloci_io::fuzz_bodies::read_bam;
use gloci_testkit::fixtures::bgzip;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // data[0] is the config byte (see gloci_io::fuzz_bodies). Bit 0 set wraps the
    // rest in BGZF so the fuzzer gets past that layer; the seeds are pre-wrapped.
    match data.split_first() {
        Some((&c, rest)) if c & 1 == 1 => {
            let mut v = vec![c];
            v.extend(bgzip(rest));
            read_bam(&v);
        }
        _ => read_bam(data),
    }
});
