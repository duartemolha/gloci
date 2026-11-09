use gloci_engine::chunk::ChunkSink;
use gloci_engine::parallel::{all_indexed, run_chunks};
use gloci_io::{IoError, RecordWriter};
use std::num::NonZeroUsize;
use tempfile::tempdir;

struct VecWriter {
    records: Vec<Vec<u8>>,
}

impl VecWriter {
    fn new() -> Self {
        Self {
            records: Vec::new(),
        }
    }
}

impl RecordWriter for VecWriter {
    fn write_record(&mut self, raw: &[u8]) -> std::io::Result<()> {
        self.records.push(raw.to_vec());
        Ok(())
    }

    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }
}

#[test]
fn chunk_sink_in_memory_and_spill() {
    let dir = tempdir().unwrap();
    let mut sink = ChunkSink::with_mem_limit(dir.path(), 1024);

    let rec1 = b"chr1\t100\t200\trec1";
    sink.write_record(rec1).unwrap();

    let large_rec = vec![b'A'; 2000];
    sink.write_record(&large_rec).unwrap();

    let mut out = VecWriter::new();
    sink.drain_to(&mut out).unwrap();

    assert_eq!(out.records.len(), 2);
    assert_eq!(out.records[0], rec1);
    assert_eq!(out.records[1], large_rec);
}

#[test]
fn chunk_sink_drop_cleans_up() {
    let dir = tempdir().unwrap();
    {
        let mut sink = ChunkSink::with_mem_limit(dir.path(), 64);
        let large_rec = vec![b'X'; 200];
        sink.write_record(&large_rec).unwrap();
    }
    let entries: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
    assert!(
        entries.is_empty(),
        "spill dir not cleaned up on drop: {:?}",
        entries
    );
}

#[test]
fn run_chunks_order_and_determinism() {
    let dir = tempdir().unwrap();
    let chroms = vec![
        b"chr1".to_vec(),
        b"chr2".to_vec(),
        b"chr3".to_vec(),
        b"chr4".to_vec(),
        b"chr5".to_vec(),
        b"chr6".to_vec(),
    ];

    let mut out = VecWriter::new();
    let threads = NonZeroUsize::new(4).unwrap();

    let res = run_chunks(
        chroms.clone(),
        threads,
        dir.path(),
        &|chrom: &Vec<u8>, sink: &mut ChunkSink| {
            let row = format!("{}\t10\t20", String::from_utf8_lossy(chrom));
            sink.write_record(row.as_bytes()).unwrap();
            Ok(chrom.to_vec())
        },
        &mut out,
    )
    .unwrap();

    assert_eq!(res, chroms);
    assert_eq!(out.records.len(), 6);
    for (i, chrom) in chroms.iter().enumerate() {
        let expected = format!("{}\t10\t20", String::from_utf8_lossy(chrom));
        assert_eq!(out.records[i], expected.as_bytes());
    }
}

#[test]
fn run_chunks_error_stops_and_cleans_up() {
    let dir = tempdir().unwrap();
    let chroms = vec![
        b"chr1".to_vec(),
        b"chr2".to_vec(),
        b"chr3".to_vec(),
        b"chr4".to_vec(),
    ];

    let mut out = VecWriter::new();
    let threads = NonZeroUsize::new(2).unwrap();

    let err = run_chunks(
        chroms,
        threads,
        dir.path(),
        &|chrom: &Vec<u8>, sink: &mut ChunkSink| {
            if chrom == b"chr2" {
                return Err(IoError::Parse {
                    path: "test".to_string(),
                    loc: gloci_io::Loc::Line(2),
                    msg: "malformed record".to_string(),
                });
            }
            sink.write_record(b"dummy").unwrap();
            Ok(())
        },
        &mut out,
    );

    assert!(err.is_err());
    let entries: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
    assert!(
        entries.is_empty(),
        "tmp dir not empty after error: {:?}",
        entries
    );
}

#[test]
fn test_all_indexed() {
    assert!(!all_indexed(&[]));
    assert!(!all_indexed(&["-"]));
    assert!(!all_indexed(&["nonexistent_file"]));
}

/// M3: when several chunks fail, the error of the first in job order is
/// returned, whatever finishes first; a panicking worker is an error, not a
/// hang.
#[test]
fn run_chunks_reports_the_first_error_in_order_and_panics() {
    let dir = tempdir().unwrap();
    let jobs: Vec<u32> = (0..8).collect();
    for threads in [2, 4, 8] {
        let mut out = VecWriter::new();
        let err = run_chunks(
            jobs.clone(),
            NonZeroUsize::new(threads).unwrap(),
            dir.path(),
            &|&j: &u32, _sink: &mut ChunkSink| {
                // Job 1 fails late, job 5 fails at once.
                match j {
                    1 => {
                        std::thread::sleep(std::time::Duration::from_millis(100));
                        Err(IoError::Internal("job 1".into()))
                    }
                    5 => Err(IoError::Internal("job 5".into())),
                    _ => Ok(()),
                }
            },
            &mut out,
        )
        .unwrap_err();
        assert!(err.to_string().contains("job 1"), "-t {threads}: {err}");
    }
    let mut out = VecWriter::new();
    let err = run_chunks(
        jobs,
        NonZeroUsize::new(4).unwrap(),
        dir.path(),
        &|&j: &u32, _sink: &mut ChunkSink| {
            if j == 3 {
                panic!("boom in job 3");
            }
            Ok(())
        },
        &mut out,
    )
    .unwrap_err();
    assert!(err.to_string().contains("boom in job 3"), "{err}");
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
}
