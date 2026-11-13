#![no_main]
use gloci_io::fuzz_bodies::read_cram;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // data[0] is the config byte (see gloci_io::fuzz_bodies). CRAM carries its own framing and
    // compression; under `--cfg fuzzing` its checksums are not enforced, so mutated bytes reach
    // the slice checks and the decoders.
    read_cram(data);
});
