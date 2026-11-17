//! `convert` to bigWig and bigBed (D171): every record of the input goes to a [`BbiWriter`],
//! which checks it, and a refused record is named by its place in the input.

use gloci_core::ChromTable;
use gloci_io::{
    BbiKind, BbiWriter, Format, IoError, Loc, OrderGuard, PushError, RecordBatch, RecordReader,
    UsageError, next_batch_size,
};

/// The input formats each output takes: bigWig from BedGraph, Bed (4 numeric columns), Wig,
/// BigWig; bigBed from Bed, BigBed. Anything else: a usage error naming them.
pub fn check_track_input(format: Format, kind: BbiKind) -> Result<(), UsageError> {
    let ok = match kind {
        BbiKind::BigWig => matches!(
            format,
            Format::BedGraph | Format::Bed | Format::Wig | Format::BigWig
        ),
        BbiKind::BigBed => matches!(format, Format::Bed | Format::BigBed),
    };
    if ok {
        return Ok(());
    }
    let (out, takes) = match kind {
        BbiKind::BigWig => (
            "bigWig",
            "bedGraph, a 4-column BED, WIG or bigWig (bedGraphToBigWig, wigToBigWig)",
        ),
        BbiKind::BigBed => ("bigBed", "BED or bigBed (bedToBigBed)"),
    };
    let alt = if format == Format::BedGraph && kind == BbiKind::BigBed {
        "; a bedGraph is signal, so write a bigWig (-o x.bw), or give a BED file for a bigBed"
    } else {
        ""
    };
    Err(UsageError(format!(
        "convert writes a {out} from {takes}, not from {} input{alt}",
        format.name()
    )))
}

/// Every record of `input` to `w`, in input order: each batch is checked for order (every
/// chromosome's records together and ascending), then pushed. A refused record is an
/// [`IoError::Parse`] naming the input and the record's line (text input) or record number
/// (binary input).
pub fn convert_to_track(
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    w: &mut BbiWriter,
) -> Result<(), IoError> {
    let name = input.name().to_string();
    let binary = input.format().is_binary();
    let mut guard = OrderGuard::contiguous(&name, input.format());
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, max)?;
        // Each record is checked for order and then pushed, so the first bad line is the one
        // reported, whether the order or the writer refuses it.
        for (i, m) in batch.recs.iter().enumerate() {
            guard.check_range(chroms, &batch, i..i + 1)?;
            match w.push(batch.raw(i)) {
                Ok(()) => {}
                Err(PushError::Writer(e)) => return Err(e),
                Err(PushError::Record(msg)) => {
                    return Err(IoError::Parse {
                        path: name,
                        loc: if binary {
                            Loc::Record(m.loc)
                        } else {
                            Loc::Line(m.loc)
                        },
                        msg,
                    });
                }
            }
        }
        guard.settle(chroms, &batch, batch.recs.len());
        max = next_batch_size(&batch, max);
        if !more {
            crate::skipped::log_reader(input.as_ref());
            return Ok(());
        }
    }
}
