//! `depth` shares its decode threads among its inputs (N1): each input gets `-t / inputs`, at
//! least 1, so many inputs at a high `-t` do not start inputs × `-t` threads. Its own test
//! binary, so no other test's threads are counted.

use gloci_core::VcfSpan;
use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions};
use gloci_ops::depth::{AllPositions, CountRules, DepthOptions, depth, open_depth_inputs};
use std::num::NonZeroUsize;
use std::path::PathBuf;

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/samtools")
        .join(name)
        .display()
        .to_string()
}

fn ro() -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions::default()),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
    }
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

/// The threads of this process.
#[cfg(target_os = "linux")]
fn tasks() -> usize {
    std::fs::read_dir("/proc/self/task").unwrap().count()
}

#[cfg(target_os = "linux")]
#[test]
fn depth_decode_threads_do_not_grow_with_the_inputs() {
    let bam = fixture("depth.bam");
    let baseline = tasks();
    for (inputs, t) in [(60, 8), (2, 8), (3, 4), (1, 4)] {
        // The threads of the inputs dropped last time end on their own time.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        while tasks() > baseline && std::time::Instant::now() < deadline {
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        let paths = vec![bam.clone(); inputs];
        let before = tasks();
        let opened = open_depth_inputs(&paths, &ro(), threads(t)).unwrap();
        let started = tasks().saturating_sub(before);
        // Each input with a share of 2 or more decodes with that many workers and a reader
        // thread; a share of 1 decodes on the calling thread.
        let bound = t + inputs.min(t);
        assert!(
            started <= bound,
            "{inputs} inputs at -t {t} started {started} threads (at most {bound})"
        );
        drop(opened);
    }
    depth_output_is_the_same_whatever_the_share();
}

/// Called by the test above, not a test of its own: the thread counts must not see the
/// threads of a test running beside them.
fn depth_output_is_the_same_whatever_the_share() {
    let bam = fixture("depth.bam");
    let o = DepthOptions {
        rules: CountRules::default(),
        min_read_length: 0,
        all: AllPositions::RefsWithReads,
        regions: None,
        skip_mate_overlap: false,
        column_names: true,
    };
    let run = |inputs: usize, t: usize| {
        let mut out = Vec::new();
        depth(&vec![bam.clone(); inputs], &o, &ro(), threads(t), &mut out).unwrap();
        out
    };
    for inputs in [1, 3, 12] {
        let one = run(inputs, 1);
        assert!(!one.is_empty());
        for t in [2, 4, 8] {
            assert!(run(inputs, t) == one, "{inputs} inputs at -t {t}");
        }
    }
}
