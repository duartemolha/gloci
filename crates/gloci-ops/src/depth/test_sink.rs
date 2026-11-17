//! Test doubles for the depth core's unit tests: a sink that records what it is given, and a
//! SAM text opened as an input.

use super::rows::DepthSink;
use gloci_core::Interval;
use gloci_io::{AlnRef, AlnSource, ReaderOptions, open_aln_source_from};
use std::io::{self, Cursor};

/// What a [`Recorder`] saw besides positions, in the order it saw it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Event {
    /// `reference_done(rid)`.
    Done(usize),
    /// `read(input, rid, ..)`, with the record's name and the span it was given.
    Read {
        input: usize,
        rid: usize,
        name: String,
        span: Option<Interval>,
    },
    /// `placed(rid)`.
    Placed(usize),
}

/// A [`DepthSink`] that keeps positions and zeros as merged runs.
#[derive(Default)]
pub(crate) struct Recorder {
    /// `(rid, start, end, depth)`: adjacent positions of one reference with one depth (summed
    /// over the inputs) are one run.
    pub runs: Vec<(usize, u64, u64, u32)>,
    pub events: Vec<Event>,
    /// For each `reference_done`, how many runs had been given by then.
    pub runs_at_done: Vec<usize>,
}

impl Recorder {
    fn push(&mut self, rid: usize, a: u64, b: u64, depth: u32) {
        if let Some(last) = self.runs.last_mut()
            && last.0 == rid
            && last.2 == a
            && last.3 == depth
        {
            last.2 = b;
        } else {
            self.runs.push((rid, a, b, depth));
        }
    }

    /// The depth at `pos` of `rid`, 0 where no run covers it.
    pub(crate) fn depth_at(&self, rid: usize, pos: u64) -> u32 {
        self.runs
            .iter()
            .find(|r| r.0 == rid && r.1 <= pos && pos < r.2)
            .map_or(0, |r| r.3)
    }
}

impl DepthSink for Recorder {
    fn position(&mut self, rid: usize, pos: u64, depths: &[u32]) -> io::Result<()> {
        self.push(rid, pos, pos + 1, depths.iter().sum());
        Ok(())
    }

    fn zeros(&mut self, rid: usize, a: u64, b: u64) -> io::Result<()> {
        self.push(rid, a, b, 0);
        Ok(())
    }

    fn reference_done(&mut self, rid: usize) -> io::Result<()> {
        self.events.push(Event::Done(rid));
        self.runs_at_done.push(self.runs.len());
        Ok(())
    }

    fn read(&mut self, input: usize, rid: usize, a: &AlnRef<'_>, span: Option<Interval>) {
        self.events.push(Event::Read {
            input,
            rid,
            name: String::from_utf8_lossy(a.name()).into_owned(),
            span,
        });
    }

    fn placed(&mut self, rid: usize) {
        self.events.push(Event::Placed(rid));
    }
}

/// `text` (SAM) as an input called `test.sam`.
pub(crate) fn sam_source(text: &str, ro: &ReaderOptions) -> Box<dyn AlnSource> {
    open_aln_source_from(
        "test.sam",
        Box::new(Cursor::new(text.as_bytes().to_vec())),
        ro,
    )
    .expect("the test SAM opens")
}
