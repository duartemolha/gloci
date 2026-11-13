#![no_main]
use gloci_io::fuzz_bodies::fasta_index;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // data[0] is the config byte (unused here, kept so that every seed has the same layout).
    fasta_index(data);
});
