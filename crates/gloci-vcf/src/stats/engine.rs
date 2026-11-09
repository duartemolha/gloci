//! The `vcf stats` engine: every requested statistic in one pass.

use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::SiteStats;
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
    /// vcftools' file extension, without the output prefix or leading dot.
    fn ext(&self) -> &'static str;
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
