//! `gloci vcf stats`: every requested statistic in one pass over the records
//! that pass the shared selection stage (decision D72).
//!
//! Each statistic writes its own temporary file next to its destination.
//! Only when the whole input has been read and every statistic has finished
//! are the files renamed into place, so an error publishes nothing; the
//! outputs of a previous run that are replaced are put back if publishing
//! fails part-way.

use crate::vcf_select::{SelectStage, SelectedInput, VcfSelect, run_selection};
use gloci_io::{AtomicFile, BcfDicts, IoError, Loc, UsageError};
use gloci_vcf::{
    BcfSite, FilterSummary, Freq, GenoDepth, GtFields, Hardy, Heterozygosity, IndelHist,
    MaskedGenotypes, MissingSample, MissingSite, RewritePlan, SampleDepth, Singletons, SiteDepth,
    SiteMeanDepth, SiteQuality, SiteView, SnpDensity, Stat, StatSink, StatWrite, StatsEngine,
    Summary, TsTv,
};
use std::cell::RefCell;
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
use std::rc::Rc;

const COMMAND: &str = "vcf stats";

/// Which statistics to compute; each names the file it writes, as vcftools
/// does (`PREFIX.frq`, `PREFIX.idepth`, ...).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StatsRequest {
    pub freq: bool,
    pub counts: bool,
    pub depth: bool,
    pub site_depth: bool,
    pub site_mean_depth: bool,
    pub geno_depth: bool,
    pub site_quality: bool,
    pub missing_site: bool,
    pub missing_indv: bool,
    pub het: bool,
    pub hardy: bool,
    pub singletons: bool,
    pub filter_summary: bool,
    pub tstv_summary: bool,
    pub indel_hist: bool,
    /// Bin size of `--snp-density`.
    pub snp_density: Option<u64>,
    /// `PREFIX.summary.txt`.
    pub summary: bool,
    /// `PREFIX.summary.json`.
    pub json: bool,
}

impl StatsRequest {
    /// True when no statistic is requested.
    pub fn is_empty(&self) -> bool {
        *self == StatsRequest::default()
    }

    /// True when a requested statistic reads genotypes. Only quality, the
    /// FILTER summary, Ts/Tv, the indel histogram and SNP density do not.
    pub fn needs_genotypes(&self) -> bool {
        let site_only = StatsRequest {
            site_quality: self.site_quality,
            filter_summary: self.filter_summary,
            tstv_summary: self.tstv_summary,
            indel_hist: self.indel_hist,
            snp_density: self.snp_density,
            ..StatsRequest::default()
        };
        *self != site_only
    }
}

/// Options for `vcf stats`.
#[derive(Clone, Debug)]
pub struct StatsOptions {
    pub select: VcfSelect,
    pub request: StatsRequest,
    /// Output files are `PREFIX.<extension>`.
    pub out_prefix: PathBuf,
}

/// What a `vcf stats` run did.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct StatsSummary {
    pub records_in: u64,
    pub records_selected: u64,
    /// The files written, in publication order.
    pub files: Vec<PathBuf>,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn write_err(path: &str, source: io::Error) -> IoError {
    IoError::Write {
        path: path.to_string(),
        source,
    }
}

/// Check the request before any input is read.
pub fn check_request(o: &StatsOptions) -> Result<(), IoError> {
    o.select.validate()?;
    if o.request.is_empty() {
        return Err(usage(
            "no statistic requested; add at least one of --freq, --counts, --depth, \
             --site-depth, --site-mean-depth, --geno-depth, --site-quality, --missing-site, \
             --missing-indv, --het, --hardy, --singletons, --filter-summary, --tstv-summary, \
             --indel-hist, --snp-density BIN, --summary or --json",
        ));
    }
    if o.request.snp_density == Some(0) {
        return Err(usage(
            "--snp-density needs a bin size above 0; give the bin width in bases, such as \
             --snp-density 10000",
        ));
    }
    if o.out_prefix.as_os_str().is_empty() {
        return Err(usage(
            "--out-prefix is empty; give a prefix such as --out-prefix results/cohort",
        ));
    }
    Ok(())
}

/// The statistics a request asks for, over the selected samples
/// (`samples`: index in the record and name).
fn build_stats(r: &StatsRequest, samples: &[(usize, String)]) -> Vec<Box<dyn Stat>> {
    let mut v: Vec<Box<dyn Stat>> = Vec::new();
    if r.freq {
        v.push(Box::new(Freq::frequencies()));
    }
    if r.counts {
        v.push(Box::new(Freq::counts()));
    }
    if r.depth {
        v.push(Box::new(SampleDepth::new(samples.to_vec())));
    }
    if r.site_depth {
        v.push(Box::new(SiteDepth::new()));
    }
    if r.site_mean_depth {
        v.push(Box::new(SiteMeanDepth::new()));
    }
    if r.geno_depth {
        v.push(Box::new(GenoDepth::new(samples.to_vec())));
    }
    if r.site_quality {
        v.push(Box::new(SiteQuality::new()));
    }
    if r.missing_site {
        v.push(Box::new(MissingSite::new()));
    }
    if r.missing_indv {
        v.push(Box::new(MissingSample::new(samples.to_vec())));
    }
    if r.het {
        v.push(Box::new(Heterozygosity::new(samples.to_vec())));
    }
    if r.hardy {
        v.push(Box::new(Hardy::new()));
    }
    if r.singletons {
        v.push(Box::new(Singletons::new(samples.to_vec())));
    }
    if r.filter_summary {
        v.push(Box::new(FilterSummary::new()));
    }
    if r.tstv_summary {
        v.push(Box::new(TsTv::new()));
    }
    if r.indel_hist {
        v.push(Box::new(IndelHist::new()));
    }
    if let Some(bin) = r.snp_density {
        v.push(Box::new(SnpDensity::new(bin)));
    }
    if r.summary || r.json {
        v.push(Box::new(Summary::new(
            samples.len() as u64,
            r.summary,
            r.json,
        )));
    }
    v
}

/// Files that have been written and flushed but not yet renamed into place.
/// A statistic "commits" by handing its file here; [`Staged::publish`]
/// renames them all once every statistic has finished.
#[derive(Default)]
struct Staged {
    files: Vec<(AtomicFile, PathBuf)>,
}

impl Staged {
    /// Rename every staged file onto its destination. A destination that
    /// already exists (the previous run's output) is first moved aside to a
    /// registered temp name in the same directory. When every file is in
    /// place the backups are deleted; if a step fails, the files published in
    /// this run are removed and the backups put back, so a failed run leaves
    /// the previous outputs as they were and no partial set of new ones.
    fn publish(self) -> Result<Vec<PathBuf>, IoError> {
        // A directory on a destination cannot be replaced; say so before
        // anything is moved.
        for (_, dest) in &self.files {
            if std::fs::symlink_metadata(dest).is_ok_and(|m| m.is_dir()) {
                let d = dest.display().to_string();
                return Err(IoError::Detect {
                    path: d.clone(),
                    msg: format!(
                        "cannot write {d}: a directory with that name exists; \
remove or rename it, or choose another --out-prefix"
                    ),
                });
            }
        }
        let mut done: Vec<PathBuf> = Vec::new();
        // (destination, where its previous content waits)
        let mut backups: Vec<(PathBuf, PathBuf)> = Vec::new();
        for (file, dest) in self.files {
            let step = (|| -> io::Result<()> {
                if std::fs::symlink_metadata(&dest).is_ok() {
                    let backup = back_up(
                        &dest,
                        |a, b| std::fs::hard_link(a, b),
                        |a, b| std::fs::rename(a, b),
                    )?;
                    backups.push((dest.clone(), backup));
                }
                file.commit()
            })();
            if let Err(e) = step {
                let kept = roll_back(&done, &backups);
                let err = write_err(&dest.display().to_string(), e);
                return Err(if kept.is_empty() {
                    err
                } else {
                    IoError::Detect {
                        path: dest.display().to_string(),
                        msg: format!(
                            "{err}; the previous outputs could not all be put back: {kept}"
                        ),
                    }
                });
            }
            done.push(dest);
        }
        for (_, backup) in &backups {
            let _ = std::fs::remove_file(backup);
        }
        Ok(done)
    }
}

/// Keep the previous content of `dest` under a registered name next to it.
/// `dest` stays where it is, as a hard link, so the rename that publishes the
/// new file replaces it atomically and a signal never finds it missing. When
/// the link cannot be made (a file system without hard links) the file is
/// moved aside instead; if that fails too, the reserved name is removed.
fn back_up(
    dest: &std::path::Path,
    link: impl Fn(&std::path::Path, &std::path::Path) -> io::Result<()>,
    rename: impl Fn(&std::path::Path, &std::path::Path) -> io::Result<()>,
) -> io::Result<PathBuf> {
    let backup = reserve_backup_name(dest)?;
    // The name is registered, so it stays reserved for this run; a link needs
    // it free.
    std::fs::remove_file(&backup)?;
    if link(dest, &backup).is_ok() {
        return Ok(backup);
    }
    rename(dest, &backup).inspect_err(|_| {
        let _ = std::fs::remove_file(&backup);
    })?;
    Ok(backup)
}

/// A new, registered, empty file next to `dest` to move its previous content
/// onto (the move replaces it).
fn reserve_backup_name(dest: &std::path::Path) -> io::Result<PathBuf> {
    let name = dest.file_name().unwrap_or_default();
    for n in 0..1000u32 {
        let mut b: OsString = name.to_owned();
        b.push(format!(".bak-{}-{n}", std::process::id()));
        let path = dest.with_file_name(b);
        match gloci_io::create_temp_file(&path) {
            Ok(_) => return Ok(path),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    Err(io::Error::other("cannot create a backup file name"))
}

/// Remove the files published in this run and move the backups back. Returns
/// a description of the previous outputs that could not be restored (empty
/// when all were).
fn roll_back(done: &[PathBuf], backups: &[(PathBuf, PathBuf)]) -> String {
    for p in done {
        let _ = std::fs::remove_file(p);
    }
    let mut kept = Vec::new();
    for (dest, backup) in backups {
        if std::fs::rename(backup, dest).is_err() {
            kept.push(format!(
                "{} is kept as {}",
                dest.display(),
                backup.display()
            ));
        } else {
            // A hard-linked backup of a dest that was never replaced is the
            // same file: the rename keeps both names.
            let _ = std::fs::remove_file(backup);
        }
    }
    kept.join(", ")
}

struct PrefixSink {
    prefix: PathBuf,
    staged: Rc<RefCell<Staged>>,
}

struct StagedFile {
    file: AtomicFile,
    dest: PathBuf,
    staged: Rc<RefCell<Staged>>,
}

impl Write for StagedFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.file.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

impl StatWrite for StagedFile {
    fn commit(mut self: Box<Self>) -> io::Result<()> {
        self.file.flush()?;
        let s = *self;
        s.staged.borrow_mut().files.push((s.file, s.dest));
        Ok(())
    }
}

impl StatSink for PrefixSink {
    fn create(&mut self, ext: &str) -> io::Result<Box<dyn StatWrite>> {
        let mut name: OsString = self.prefix.clone().into_os_string();
        name.push(".");
        name.push(ext);
        let dest = PathBuf::from(name);
        let file = AtomicFile::create(&dest).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!(
                    "cannot create {}: {e}; check that the directory exists and is writable",
                    dest.display()
                ),
            )
        })?;
        Ok(Box::new(StagedFile {
            file,
            dest,
            staged: Rc::clone(&self.staged),
        }))
    }
}

/// The optional FORMAT fields the selection and the statistics read: DP for a
/// depth statistic, a summary, a DP mask or a mean-depth site filter; GQ for a
/// GQ mask; FT for an FT mask. An expression may read any, so it takes all.
fn decode_fields(o: &StatsOptions) -> GtFields {
    if o.select.expr.is_some() {
        return GtFields::ALL;
    }
    let (r, m, s) = (&o.request, &o.select.mask, &o.select.sites);
    GtFields {
        dp: r.depth
            || r.site_depth
            || r.site_mean_depth
            || r.geno_depth
            || r.summary
            || r.json
            || m.min_dp.is_some()
            || m.max_dp.is_some()
            || s.min_mean_dp.is_some()
            || s.max_mean_dp.is_some(),
        gq: m.min_gq.is_some(),
        ft: m.mask_filtered,
    }
}

/// The selected samples as (index in the record, name).
fn selected_samples(stage: &SelectStage<'_>) -> Vec<(usize, String)> {
    let all = &stage.header().samples;
    let keep = stage.plan().selection.as_ref().map(|s| &s.keep);
    all.iter()
        .enumerate()
        .filter(|(i, _)| keep.is_none_or(|k| k[*i]))
        .map(|(i, s)| (i, String::from_utf8_lossy(s).into_owned()))
        .collect()
}

/// Compute the requested statistics of the records of `input` and write
/// `PREFIX.<extension>` for each.
pub fn vcf_stats(input: SelectedInput, o: &StatsOptions) -> Result<StatsSummary, IoError> {
    check_request(o)?;
    let (reader, mut chroms) = input.into_parts();
    let base = RewritePlan::default();
    // The selection computes the site statistics (under the masking rules)
    // whenever a statistic needs them.
    let need_genotypes = o.request.needs_genotypes();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), &o.select, &base, need_genotypes)?;
    stage.share_decoded(decode_fields(o));
    let name = stage.name().to_string();
    let samples = selected_samples(&stage);
    let mut engine = StatsEngine::new(build_stats(&o.request, &samples));
    debug_assert_eq!(engine.needs_genotypes(), need_genotypes);

    let staged = Rc::new(RefCell::new(Staged::default()));
    let mut sink = PrefixSink {
        prefix: o.out_prefix.clone(),
        staged: Rc::clone(&staged),
    };
    let dest = o.out_prefix.display().to_string();
    engine.start(&mut sink).map_err(|e| write_err(&dest, e))?;

    let header = stage.header().clone();
    let dicts = if stage.input_is_bcf() {
        match &header.bcf {
            Some(d) => d.clone(),
            None => BcfDicts::from_header(&header.raw).map_err(|e| IoError::Detect {
                path: name.clone(),
                msg: e,
            })?,
        }
    } else {
        BcfDicts::default()
    };
    let is_bcf = stage.input_is_bcf();
    let mut bcf_site = BcfSite::default();
    let mut selected = 0u64;
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        selected += 1;
        let parse_err = |msg: String| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(hit.rec_no),
            msg,
        };
        let view: SiteView = if is_bcf {
            bcf_site.view(hit.raw, &dicts, hit.meta.span)
        } else {
            SiteView::from_text(hit.raw, hit.meta.span)
        }
        .map_err(parse_err)?;
        if !need_genotypes {
            return engine.site(&view, None, None).map_err(parse_err);
        }
        // Every selected sample was decoded by the selection stage.
        let shared = hit.decoded.reader();
        let masked = MaskedGenotypes {
            inner: &shared,
            masked: hit.masked,
        };
        engine
            .site(&view, Some(&masked), hit.stats)
            .map_err(parse_err)
    })?;
    engine.finish().map_err(|e| write_err(&dest, e))?;
    drop(engine);
    drop(sink);
    let staged = Rc::try_unwrap(staged)
        .map_err(|_| IoError::Internal("the staged outputs are still shared".into()))?
        .into_inner();
    let files = staged.publish()?;
    Ok(StatsSummary {
        records_in,
        records_selected: selected,
        files,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn staged(dir: &std::path::Path, name: &str, body: &str) -> (AtomicFile, PathBuf) {
        let dest = dir.join(name);
        let mut f = AtomicFile::create(&dest).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        (f, dest)
    }

    fn names(dir: &std::path::Path) -> Vec<String> {
        let mut v: Vec<_> = std::fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        v.sort();
        v
    }

    #[test]
    fn publish_puts_the_previous_outputs_back_when_a_later_rename_fails() {
        let d = tempfile::tempdir().unwrap();
        std::fs::write(d.path().join("o.frq"), "previous frq").unwrap();
        std::fs::write(d.path().join("o.hwe"), "previous hwe").unwrap();
        let first = staged(d.path(), "o.frq", "new frq");
        let second = staged(d.path(), "o.hwe", "new hwe");
        // The second temp file vanishes, so its rename fails after the first
        // file was published.
        std::fs::remove_file(second.0.temp_path()).unwrap();
        let err = Staged {
            files: vec![first, second],
        }
        .publish()
        .unwrap_err();
        assert!(err.to_string().contains("o.hwe"), "{err}");
        assert_eq!(names(d.path()), vec!["o.frq", "o.hwe"]);
        assert_eq!(
            std::fs::read_to_string(d.path().join("o.frq")).unwrap(),
            "previous frq"
        );
        assert_eq!(
            std::fs::read_to_string(d.path().join("o.hwe")).unwrap(),
            "previous hwe"
        );
    }

    #[test]
    fn back_up_hard_links_and_leaves_the_destination_in_place() {
        let d = tempfile::tempdir().unwrap();
        let dest = d.path().join("o.frq");
        std::fs::write(&dest, "previous").unwrap();
        let backup = back_up(
            &dest,
            |a, b| std::fs::hard_link(a, b),
            |a, b| std::fs::rename(a, b),
        )
        .unwrap();
        assert_eq!(std::fs::read_to_string(&dest).unwrap(), "previous");
        assert_eq!(std::fs::read_to_string(&backup).unwrap(), "previous");
    }

    #[test]
    fn back_up_falls_back_to_a_move_and_cleans_up_when_that_fails() {
        let d = tempfile::tempdir().unwrap();
        let dest = d.path().join("o.frq");
        std::fs::write(&dest, "previous").unwrap();
        let no_link = |_: &std::path::Path, _: &std::path::Path| -> io::Result<()> {
            Err(io::Error::other("no links"))
        };
        let backup = back_up(&dest, no_link, |a, b| std::fs::rename(a, b)).unwrap();
        assert!(!dest.exists());
        assert_eq!(std::fs::read_to_string(&backup).unwrap(), "previous");
        std::fs::rename(&backup, &dest).unwrap();
        let no_move = |_: &std::path::Path, _: &std::path::Path| -> io::Result<()> {
            Err(io::Error::other("no move"))
        };
        back_up(&dest, no_link, no_move).unwrap_err();
        assert_eq!(names(d.path()), vec!["o.frq"]);
    }

    #[test]
    fn publish_removes_what_it_published_when_there_was_no_previous_output() {
        let d = tempfile::tempdir().unwrap();
        let first = staged(d.path(), "o.frq", "new frq");
        let second = staged(d.path(), "o.hwe", "new hwe");
        std::fs::remove_file(second.0.temp_path()).unwrap();
        Staged {
            files: vec![first, second],
        }
        .publish()
        .unwrap_err();
        assert!(names(d.path()).is_empty(), "{:?}", names(d.path()));
    }
}
