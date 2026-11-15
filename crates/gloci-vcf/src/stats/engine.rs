//! The `vcf stats` engine: every requested statistic in one pass.

use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::{SampleTally, SiteStats, TallyWants};
use std::io::{self, BufWriter, Write};

/// One statistic's output file. [`StatWrite::commit`] publishes it; dropping
/// it without committing discards it. `gloci-io`'s `AtomicFile` (a `Write`
/// with `commit`) is the intended implementation.
pub trait StatWrite: Write {
    fn commit(self: Box<Self>) -> io::Result<()>;
}

/// Where statistic outputs go. A statistic asks for its file by extension
/// (`frq`, `FILTER.summary`, ...) when the run starts.
pub trait StatSink {
    fn create(&mut self, ext: &str) -> io::Result<Box<dyn StatWrite>>;
}

/// A statistic's buffered output plus the row buffer it formats into. Per
/// row, a statistic clears `buf`, formats into it, and calls `row`.
pub struct Out {
    w: Option<BufWriter<Box<dyn StatWrite>>>,
    pub buf: Vec<u8>,
}

impl Out {
    pub fn closed() -> Out {
        Out {
            w: None,
            buf: Vec::with_capacity(256),
        }
    }

    pub fn open(&mut self, sink: &mut dyn StatSink, ext: &str) -> io::Result<()> {
        self.w = Some(BufWriter::with_capacity(1 << 16, sink.create(ext)?));
        Ok(())
    }

    /// Write `buf` (a finished row) to the file and clear it.
    pub fn row(&mut self) -> Result<(), String> {
        let w = self.w.as_mut().ok_or("statistic output is not open")?;
        let r = w.write_all(&self.buf);
        self.buf.clear();
        r.map_err(|e| format!("cannot write the statistic output: {e}"))
    }

    /// Flush and commit the file.
    pub fn commit(&mut self) -> io::Result<()> {
        match self.w.take() {
            Some(w) => w.into_inner().map_err(|e| e.into_error())?.commit(),
            None => Ok(()),
        }
    }
}

/// A statistic. `site` runs once per surviving record in input order;
/// per-site statistics stream their rows from it, whole-file statistics
/// accumulate in bounded maps and write in `finish`.
pub trait Stat {
    /// Open the output file. Called once before the first site; this
    /// replaces passing the sink to `finish` because per-site rows are
    /// written as they are produced.
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()>;
    /// `gt` is `Some` when any statistic needs genotypes (masked samples read
    /// as missing); `st` is the site's [`SiteStats`] under the same masking.
    fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String>;
    /// Write what remains and commit the file.
    fn finish(&mut self) -> io::Result<()>;
    fn needs_genotypes(&self) -> bool;
    /// True when `site` reads samples from `gt`. (A statistic that only needs
    /// the [`SiteStats`] does not, and a tallied one does not.)
    fn reads_genotypes(&self) -> bool {
        false
    }
    /// The per-sample counters selection workers could keep for this
    /// statistic instead (see [`SampleTally`]).
    fn tally_wants(&self) -> TallyWants {
        TallyWants::default()
    }
    /// Switch to workers' counters: `site` stops reading samples, and
    /// [`Stat::absorb`] takes the totals.
    fn use_tally(&mut self) {}
    /// A tallied record's bitset of fully called samples, after `site`
    /// (only for a site that contributes to `--het`).
    fn called(&mut self, _called: &[u64]) {}
    /// The workers' totals, before `finish`.
    fn absorb(&mut self, _t: &SampleTally) {}
    /// vcftools' file extension, without the output prefix or leading dot.
    fn ext(&self) -> &'static str;
    /// The input as a command line names it (`-` for standard input), for the
    /// `gloci sort` hint of a statistic that needs sorted input.
    fn set_input(&mut self, _name: &str) {}
}

/// Runs every requested statistic in one pass.
pub struct StatsEngine {
    stats: Vec<Box<dyn Stat>>,
}

impl StatsEngine {
    pub fn new(stats: Vec<Box<dyn Stat>>) -> StatsEngine {
        StatsEngine { stats }
    }

    /// True when the caller must supply genotypes and [`SiteStats`].
    pub fn needs_genotypes(&self) -> bool {
        self.stats.iter().any(|s| s.needs_genotypes())
    }

    /// True when `site` must be given the samples.
    pub fn reads_genotypes(&self) -> bool {
        self.stats.iter().any(|s| s.reads_genotypes())
    }

    /// The counters workers could keep for this request.
    pub fn tally_wants(&self) -> TallyWants {
        self.stats
            .iter()
            .fold(TallyWants::default(), |w, s| w.or(s.tally_wants()))
    }

    /// Take the workers' counters (see [`StatsEngine::tally_wants`]); `site`
    /// then needs no samples for those statistics.
    pub fn use_tally(&mut self) {
        for s in &mut self.stats {
            s.use_tally();
        }
    }

    /// Name the input (as a command line names it) in the messages of the
    /// statistics.
    pub fn set_input(&mut self, name: &str) {
        for s in &mut self.stats {
            s.set_input(name);
        }
    }

    /// Give each tallied statistic the bitset of fully called samples of the
    /// record just fed to `site`.
    pub fn called(&mut self, called: &[u64]) {
        for s in &mut self.stats {
            s.called(called);
        }
    }

    /// Add the workers' totals.
    pub fn absorb(&mut self, t: &SampleTally) {
        for s in &mut self.stats {
            s.absorb(t);
        }
    }

    pub fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        for s in &mut self.stats {
            s.start(sink)?;
        }
        Ok(())
    }

    /// Feed one record to every statistic. An `Err` names the record.
    pub fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        for stat in &mut self.stats {
            stat.site(s, gt, st)
                .map_err(|e| format!("{}:{}: {e}", String::from_utf8_lossy(s.chrom), s.pos))?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> io::Result<()> {
        for s in &mut self.stats {
            s.finish()?;
        }
        Ok(())
    }
}
