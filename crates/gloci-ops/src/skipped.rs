//! Reporting records skipped under `--skip-invalid`, also for inputs an
//! engine consumed: each reader's diagnostics are kept when it is dropped.

use gloci_core::ChromTable;
use gloci_io::{Compression, Diagnostics, Format, Header, IoError, RecordBatch, RecordReader};
use std::sync::{Arc, Mutex};

/// Passes records through and, when dropped, saves the reader's
/// diagnostics for [`Skipped::take`].
struct Captured {
    inner: Box<dyn RecordReader>,
    slot: Arc<Mutex<Diagnostics>>,
}

impl Drop for Captured {
    fn drop(&mut self) {
        let d = self.inner.diagnostics().clone();
        *self.slot.lock().unwrap_or_else(|e| e.into_inner()) = d;
    }
}

impl RecordReader for Captured {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn compression(&self) -> Compression {
        self.inner.compression()
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        self.inner.read_batch(chroms, batch, max)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
    fn consume_to_end(&self) -> bool {
        self.inner.consume_to_end()
    }
}

/// The diagnostics of a set of inputs, collected after the readers are gone.
#[derive(Default)]
pub(crate) struct Skipped {
    slots: Vec<(String, Arc<Mutex<Diagnostics>>)>,
    /// The slots (by index) that hold alignment inputs.
    alignments: Vec<usize>,
}

impl Skipped {
    /// Watch `r`, whose diagnostics are reported under `name`.
    pub(crate) fn watch_as(
        &mut self,
        name: &str,
        r: Box<dyn RecordReader>,
    ) -> Box<dyn RecordReader> {
        let slot = Arc::new(Mutex::new(Diagnostics::default()));
        if r.format().is_alignment() {
            self.alignments.push(self.slots.len());
        }
        self.slots.push((name.to_string(), slot.clone()));
        Box::new(Captured { inner: r, slot })
    }

    /// Watch `r` under its own name.
    pub(crate) fn watch(&mut self, r: Box<dyn RecordReader>) -> Box<dyn RecordReader> {
        let name = r.name().to_string();
        self.watch_as(&name, r)
    }

    /// Each watched input's name and diagnostics, in the order watched.
    /// Readers still alive report what they had when last dropped (nothing).
    pub(crate) fn take(&self) -> Vec<(String, Diagnostics)> {
        self.slots
            .iter()
            .map(|(n, d)| {
                let d = d.lock().unwrap_or_else(|e| e.into_inner()).clone();
                (n.clone(), d)
            })
            .collect()
    }
}

impl Skipped {
    /// Say at `-v` what each alignment input left out (D101, D102), once its reader is gone.
    pub(crate) fn log_alignment_counts(&self) {
        for &i in &self.alignments {
            let (name, d) = &self.slots[i];
            let d = d.lock().unwrap_or_else(|e| e.into_inner()).clone();
            d.log_alignment_counts(name);
        }
    }
}

/// Warn about every input that had invalid records skipped: the count and
/// the first 10 examples.
pub(crate) fn log_skipped(inputs: &[(String, Diagnostics)]) {
    for (name, d) in inputs {
        if let Some(s) = d.summary(name) {
            log::warn!("{s}");
        }
    }
}

/// Warn about a reader still held: see [`log_skipped`].
pub(crate) fn log_reader(r: &dyn RecordReader) {
    if let Some(s) = r.diagnostics().summary(r.name()) {
        log::warn!("{s}");
    }
}
