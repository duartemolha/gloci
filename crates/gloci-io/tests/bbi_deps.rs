//! Dependency smoke test for bigWig/bigBed support (D171): pins the numbers bigtools gives for
//! a small bigBed with one zoom level, which exercises its writer's use of `index_list`.

use std::collections::HashMap;
use std::fs::File;

use bigtools::bed::bedparser::BedValueError;
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BedEntry, BigBedRead, BigBedWrite, InputSortType};

#[test]
fn bigbed_zoom_counts_overlaps() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bb");

    let sizes = HashMap::from([("chr1".to_string(), 1000u32)]);
    let mut writer = BigBedWrite::new(File::create(&path).unwrap(), sizes);
    writer.options.input_sort_type = InputSortType::START;
    writer.options.manual_zoom_sizes = Some(vec![1000]);
    let entry = |start, end, rest: &str| {
        Ok::<_, BedValueError>((
            "chr1".to_string(),
            BedEntry {
                start,
                end,
                rest: rest.to_string(),
            },
        ))
    };
    let items = vec![
        entry(0, 100, "a"),
        entry(50, 150, "b"),
        entry(120, 200, "c"),
    ];
    let runtime = tokio::runtime::Builder::new_current_thread()
        .build()
        .unwrap();
    writer
        .write(
            BedParserStreamingIterator::wrap_iter(items.into_iter(), true),
            runtime,
        )
        .unwrap();

    let mut reader = BigBedRead::open(File::open(&path).unwrap()).unwrap();
    let zooms: Vec<_> = reader
        .get_zoom_interval("chr1", 0, 1000, 1000)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    let covered: u64 = zooms.iter().map(|z| z.summary.bases_covered).sum();
    let sum: f64 = zooms.iter().map(|z| z.summary.sum).sum();
    let sum_squares: f64 = zooms.iter().map(|z| z.summary.sum_squares).sum();
    let max = zooms
        .iter()
        .map(|z| z.summary.max_val)
        .fold(f64::MIN, f64::max);
    let min = zooms
        .iter()
        .map(|z| z.summary.min_val)
        .fold(f64::MAX, f64::min);
    assert_eq!(covered, 200);
    assert_eq!(sum, 280.0);
    assert_eq!(sum_squares, 440.0);
    assert_eq!(max, 2.0);
    assert_eq!(min, 1.0);

    let entries: Vec<_> = reader
        .get_interval("chr1", 0, u32::MAX)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    let rests: Vec<&str> = entries.iter().map(|e| e.rest.as_str()).collect();
    assert_eq!(rests, ["a", "b", "c"]);
}
