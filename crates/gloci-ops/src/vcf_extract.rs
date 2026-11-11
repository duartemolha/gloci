//! `gloci vcf extract`: INFO and FORMAT tables and the 012 genotype matrix
//! (decision D73). The records pass the shared selection stage first.

use crate::vcf_select::{SelectStage, SelectedInput, VcfSelect, run_selection};
use gloci_io::{BcfText, IoError, Loc, UsageError, create_file_in_temp_dir, create_temp_dir_in};
use gloci_vcf::{
    InfoKey, RewritePlan, declared_keys, format_cell, format_declared, format_index, genotype_012,
    info_cell, info_kind, rewrite_text, split_columns,
};
use std::io::{BufReader, BufWriter, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

const COMMAND: &str = "vcf extract";

/// What `vcf extract` writes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExtractMode {
    /// One row per site: `CHROM POS ID REF ALT` and one column per INFO key.
    Info(Vec<String>),
    /// FORMAT keys, one row per site (`wide`) or per site and sample (`long`).
    Format { keys: Vec<String>, long: bool },
    /// The 012 matrix, `.indv` and `.pos` files.
    Matrix012,
}

/// Options for `vcf extract`.
#[derive(Clone, Debug)]
pub struct ExtractOptions {
    pub select: VcfSelect,
    /// Threads the selection runs on (`None` = 1); see `SelectStage::set_threads`.
    pub threads: Option<NonZeroUsize>,
    pub mode: ExtractMode,
    /// `--skip-invalid`: a genotype naming an allele the ALT column lacks is
    /// written as -1 (counted in a warning) instead of failing.
    pub skip_invalid: bool,
    /// Memory budget of the 012 matrix, in bytes.
    pub max_mem: u64,
    /// Where the matrix spills when it exceeds the budget.
    pub tmp_dir: PathBuf,
}

/// Counts from an extract run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct ExtractSummary {
    pub records_in: u64,
    pub sites_out: u64,
    /// Multi-allelic sites the 012 matrix leaves out.
    pub multiallelic_skipped: u64,
    /// The 012 matrix spilled its codes to the temp directory.
    pub spilled: bool,
    /// Passes over the codes that wrote the matrix rows (blocks of samples).
    pub row_blocks: u64,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".to_string(),
        source,
    }
}

/// Reject an empty or duplicated key list before any input is read.
pub fn check_keys(flag: &str, keys: &[String]) -> Result<(), IoError> {
    if keys.is_empty() || keys.iter().any(String::is_empty) {
        return Err(usage(format!(
            "{flag} needs at least one key, such as {flag} DP; list keys separated by commas"
        )));
    }
    for (i, k) in keys.iter().enumerate() {
        if keys[..i].contains(k) {
            return Err(usage(format!("{flag} names {k} twice; list each key once")));
        }
    }
    Ok(())
}

/// The records of a stage as VCF text with the selection's sample choice and
/// masking applied: a BCF record is rendered by the converter `view -O vcf`
/// uses, so VCF and BCF input give identical lines.
struct Lines {
    bcf: Option<BcfText>,
    plan: RewritePlan,
    rewritten: Vec<u8>,
}

impl Lines {
    fn for_stage(stage: &SelectStage<'_>) -> Result<Lines, IoError> {
        Ok(Lines {
            bcf: if stage.input_is_bcf() {
                Some(BcfText::for_header(stage.name(), stage.header())?)
            } else {
                None
            },
            plan: stage.plan().clone(),
            rewritten: Vec::new(),
        })
    }

    fn line<'b>(
        &'b mut self,
        name: &str,
        hit: &crate::vcf_select::Hit<'b>,
    ) -> Result<&'b [u8], IoError> {
        let text: &[u8] = match self.bcf.as_mut() {
            Some(b) => b.render(hit.raw, hit.rec_no)?,
            None => hit.raw,
        };
        if self.plan.is_passthrough() {
            return Ok(text);
        }
        self.rewritten.clear();
        rewrite_text(text, &self.plan, hit.masked, hit.stats, &mut self.rewritten).map_err(
            |msg| IoError::Parse {
                path: name.to_string(),
                loc: Loc::Record(hit.rec_no),
                msg,
            },
        )?;
        Ok(&self.rewritten)
    }
}

fn kept_samples(stage: &SelectStage<'_>) -> Vec<Vec<u8>> {
    let all = &stage.header().samples;
    match &stage.plan().selection {
        Some(sel) => all
            .iter()
            .zip(&sel.keep)
            .filter(|(_, k)| **k)
            .map(|(s, _)| s.clone())
            .collect(),
        None => all.clone(),
    }
}

/// Every key must be declared in the header, else the column would be
/// silently empty.
fn require_declared(
    flag: &str,
    section: &'static str,
    keys: &[String],
    header: &[u8],
    declared: impl Fn(&str) -> bool,
) -> Result<(), IoError> {
    for k in keys {
        if !declared(k) {
            let list = declared_keys(header, section.as_bytes());
            return Err(usage(format!(
                "{flag} {k}: the header does not declare {section}/{k}; declared {section} keys: {}",
                if list.is_empty() {
                    "none".to_string()
                } else {
                    list.join(", ")
                }
            )));
        }
    }
    Ok(())
}

/// Write the first five columns of a site: `CHROM POS ID REF ALT`.
fn site_prefix(cols: &[&[u8]], out: &mut Vec<u8>) {
    for (i, c) in [cols[0], cols[1], cols[2], cols[3], cols[4]]
        .iter()
        .enumerate()
    {
        if i > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(c);
    }
}

/// Extract INFO or FORMAT values of the records of `input` as a TSV to `w`.
pub fn extract_table(
    input: SelectedInput,
    o: &ExtractOptions,
    w: &mut dyn Write,
) -> Result<ExtractSummary, IoError> {
    o.select.validate()?;
    let (info_keys, format_keys, long): (&[String], &[String], bool) = match &o.mode {
        ExtractMode::Info(k) => (k, &[], false),
        ExtractMode::Format { keys, long } => (&[], keys, *long),
        ExtractMode::Matrix012 => {
            return Err(usage("extract_table does not write the 012 matrix"));
        }
    };
    let base = RewritePlan::default();
    let (reader, mut chroms) = input.into_parts();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), &o.select, &base, false)?;
    stage.set_threads(o.threads.unwrap_or(NonZeroUsize::MIN));
    let name = stage.name().to_string();
    let header_raw = stage.header().raw.clone();
    let samples = kept_samples(&stage);
    require_declared("--info", "INFO", info_keys, &header_raw, |k| {
        info_kind(&header_raw, k.as_bytes()).is_some()
    })?;
    require_declared("--format", "FORMAT", format_keys, &header_raw, |k| {
        format_declared(&header_raw, k.as_bytes())
    })?;
    let infos: Vec<InfoKey> = info_keys
        .iter()
        .map(|k| InfoKey {
            name: k.as_bytes().to_vec(),
            kind: info_kind(&header_raw, k.as_bytes()).expect("declared"),
        })
        .collect();

    let mut head = b"CHROM\tPOS\tID\tREF\tALT".to_vec();
    if !infos.is_empty() {
        for k in &infos {
            head.push(b'\t');
            head.extend_from_slice(&k.name);
        }
    } else if long {
        head.extend_from_slice(b"\tSAMPLE");
        for k in format_keys {
            head.push(b'\t');
            head.extend_from_slice(k.as_bytes());
        }
    } else {
        for s in &samples {
            for k in format_keys {
                head.push(b'\t');
                head.extend_from_slice(s);
                head.push(b':');
                head.extend_from_slice(k.as_bytes());
            }
        }
    }
    head.push(b'\n');
    w.write_all(&head).map_err(write_err)?;

    let mut lines = Lines::for_stage(&stage)?;
    let mut row = Vec::new();
    let mut sites_out = 0u64;
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        let line = lines.line(&name, &hit)?;
        let parse_err = |msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(hit.rec_no),
            msg,
        };
        let mut cols: Vec<&[u8]> = Vec::new();
        split_columns(line, &mut cols).map_err(parse_err)?;
        let indexes: Vec<Option<usize>> = if format_keys.is_empty() {
            Vec::new()
        } else {
            format_keys
                .iter()
                .map(|k| format_index(cols.get(8).copied().unwrap_or(b"."), k.as_bytes()))
                .collect()
        };
        let sample_col = |i: usize| cols.get(9 + i).copied().unwrap_or(b".");
        if long && !format_keys.is_empty() {
            for (i, sample) in samples.iter().enumerate() {
                row.clear();
                site_prefix(&cols, &mut row);
                row.push(b'\t');
                row.extend_from_slice(sample);
                for idx in &indexes {
                    row.push(b'\t');
                    row.extend_from_slice(format_cell(sample_col(i), *idx));
                }
                row.push(b'\n');
                w.write_all(&row).map_err(write_err)?;
            }
        } else {
            row.clear();
            site_prefix(&cols, &mut row);
            if format_keys.is_empty() {
                for k in &infos {
                    row.push(b'\t');
                    row.extend_from_slice(info_cell(cols[7], k));
                }
            } else {
                for i in 0..samples.len() {
                    for idx in &indexes {
                        row.push(b'\t');
                        row.extend_from_slice(format_cell(sample_col(i), *idx));
                    }
                }
            }
            row.push(b'\n');
            w.write_all(&row).map_err(write_err)?;
        }
        sites_out += 1;
        Ok(())
    })?;
    Ok(ExtractSummary {
        records_in,
        sites_out,
        multiallelic_skipped: 0,
        spilled: false,
        row_blocks: 0,
    })
}

/// Fail now, naming the directory, when `dir` cannot hold the spill of
/// `extract_matrix` (it must exist and be writable), so a bad `--tmp-dir`
/// is reported before any input is read even if the matrix fits in memory.
pub fn check_tmp_dir(dir: &Path) -> Result<(), IoError> {
    create_temp_dir_in(dir, "gloci-extract-")
        .map(drop)
        .map_err(|source| IoError::Create {
            path: dir.display().to_string(),
            source,
        })
}

/// The codes of the 012 matrix, site-major, two bits each (0, 1, 2 and 3 for
/// -1). They stay in memory while they fit half the budget, then spill to a
/// registered temp file, so the matrix never needs more than the budget.
struct CodeStore {
    row_bytes: usize,
    sites: u64,
    mem: Vec<u8>,
    mem_limit: usize,
    tmp_dir: PathBuf,
    spill: Option<Spill>,
}

struct Spill {
    // Dropped (and unregistered) with the store; removes the directory.
    _dir: tempfile::TempDir,
    _guard: gloci_io::TempGuard,
    path: PathBuf,
    w: BufWriter<std::fs::File>,
}

impl CodeStore {
    fn new(n_samples: usize, max_mem: u64, tmp_dir: PathBuf) -> CodeStore {
        CodeStore {
            row_bytes: n_samples.div_ceil(4),
            sites: 0,
            mem: Vec::new(),
            mem_limit: usize::try_from(max_mem / 2).unwrap_or(usize::MAX),
            tmp_dir,
            spill: None,
        }
    }

    /// Append one site: `codes` has one entry per sample.
    fn push_site(&mut self, codes: &[i8]) -> Result<(), IoError> {
        let start = self.mem.len();
        self.mem.resize(start + self.row_bytes, 0);
        for (i, &c) in codes.iter().enumerate() {
            let v = if c < 0 { 3 } else { c as u8 };
            self.mem[start + i / 4] |= v << (2 * (i % 4));
        }
        self.sites += 1;
        if self.spill.is_some() || self.mem.len() > self.mem_limit {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), IoError> {
        if self.spill.is_none() {
            let tmp_err = |source| IoError::Create {
                path: self.tmp_dir.display().to_string(),
                source,
            };
            let (dir, guard) =
                create_temp_dir_in(&self.tmp_dir, "gloci-extract-").map_err(tmp_err)?;
            let path = dir.path().join("codes");
            let f = create_file_in_temp_dir(&path).map_err(|source| IoError::Create {
                path: path.display().to_string(),
                source,
            })?;
            self.spill = Some(Spill {
                _dir: dir,
                _guard: guard,
                path,
                w: BufWriter::new(f),
            });
        }
        let sp = self.spill.as_mut().expect("just created");
        sp.w.write_all(&self.mem).map_err(|source| IoError::Write {
            path: sp.path.display().to_string(),
            source,
        })?;
        self.mem.clear();
        Ok(())
    }

    /// Write the rows of the matrix: `block` samples at a time, each pass
    /// reading every site once and gathering the block's codes.
    fn write_rows(
        mut self,
        n_samples: usize,
        max_mem: u64,
        out: &mut dyn Write,
    ) -> Result<(bool, u64), IoError> {
        let spilled = self.spill.is_some();
        let mut blocks = 0u64;
        let sites = usize::try_from(self.sites).unwrap_or(usize::MAX);
        let budget = usize::try_from(max_mem / 2).unwrap_or(usize::MAX);
        let block = (budget / sites.max(1)).clamp(1, n_samples.max(1));
        enum Src {
            Mem(Vec<u8>),
            File(PathBuf, BufReader<std::fs::File>),
        }
        let mut src = match self.spill.take() {
            None => Src::Mem(std::mem::take(&mut self.mem)),
            Some(mut sp) => {
                // Hand the rest of the codes to the file.
                sp.w.write_all(&self.mem).map_err(|source| IoError::Write {
                    path: sp.path.display().to_string(),
                    source,
                })?;
                sp.w.flush().map_err(|source| IoError::Write {
                    path: sp.path.display().to_string(),
                    source,
                })?;
                let f = std::fs::File::open(&sp.path).map_err(|source| IoError::Create {
                    path: sp.path.display().to_string(),
                    source,
                })?;
                self.spill = Some(sp);
                Src::File(
                    self.spill.as_ref().unwrap().path.clone(),
                    BufReader::with_capacity(1 << 20, f),
                )
            }
        };
        let mut cells = vec![0u8; block * sites];
        let mut row_buf = vec![0u8; self.row_bytes];
        let mut text = Vec::new();
        let mut first = 0;
        while first < n_samples {
            let k = block.min(n_samples - first);
            match &mut src {
                Src::Mem(v) => {
                    for site in 0..sites {
                        gather(
                            &v[site * self.row_bytes..(site + 1) * self.row_bytes],
                            first,
                            k,
                            site,
                            sites,
                            &mut cells,
                        );
                    }
                }
                Src::File(path, r) => {
                    use std::io::Seek;
                    r.rewind().map_err(|source| IoError::Write {
                        path: path.display().to_string(),
                        source,
                    })?;
                    for site in 0..sites {
                        r.read_exact(&mut row_buf)
                            .map_err(|source| IoError::Write {
                                path: path.display().to_string(),
                                source,
                            })?;
                        gather(&row_buf, first, k, site, sites, &mut cells);
                    }
                }
            }
            for j in 0..k {
                text.clear();
                text.extend_from_slice((first + j).to_string().as_bytes());
                for &v in &cells[j * sites..(j + 1) * sites] {
                    text.extend_from_slice(match v {
                        0 => b"\t0",
                        1 => b"\t1",
                        2 => b"\t2",
                        _ => b"\t-1",
                    });
                }
                text.push(b'\n');
                out.write_all(&text).map_err(write_err)?;
            }
            first += k;
            blocks += 1;
        }
        Ok((spilled, blocks))
    }
}

/// Copy the codes of samples `first..first + k` from one site's packed row
/// into column `site` of each of the block's rows.
fn gather(row: &[u8], first: usize, k: usize, site: usize, sites: usize, cells: &mut [u8]) {
    for j in 0..k {
        let i = first + j;
        cells[j * sites + site] = (row[i / 4] >> (2 * (i % 4))) & 3;
    }
}

/// Write vcftools' `--012` files for the records of `input`: `matrix` (one
/// row per kept sample: its 0-based index, then one code per site), `indv`
/// (the kept sample names) and `pos` (`CHROM POS` of each site written).
/// Sites with more than one ALT allele are left out, as vcftools does. The
/// codes are held within `o.max_mem` (spilling to `o.tmp_dir` beyond it).
pub fn extract_matrix(
    input: SelectedInput,
    o: &ExtractOptions,
    matrix: &mut dyn Write,
    indv: &mut dyn Write,
    pos: &mut dyn Write,
) -> Result<ExtractSummary, IoError> {
    o.select.validate()?;
    let base = RewritePlan::default();
    let (reader, mut chroms) = input.into_parts();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), &o.select, &base, false)?;
    stage.set_threads(o.threads.unwrap_or(NonZeroUsize::MIN));
    let name = stage.name().to_string();
    let samples = kept_samples(&stage);
    for s in &samples {
        indv.write_all(s).map_err(write_err)?;
        indv.write_all(b"\n").map_err(write_err)?;
    }
    let mut store = CodeStore::new(samples.len(), o.max_mem, o.tmp_dir.clone());
    let mut lines = Lines::for_stage(&stage)?;
    let mut sites_out = 0u64;
    let mut multiallelic_skipped = 0u64;
    let mut invalid_alleles = 0u64;
    let mut site_codes = vec![0i8; samples.len()];
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        let line = lines.line(&name, &hit)?;
        let parse_err = |msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(hit.rec_no),
            msg,
        };
        let mut cols: Vec<&[u8]> = Vec::new();
        split_columns(line, &mut cols).map_err(parse_err)?;
        let n_alt = if cols[4] == b"." || cols[4].is_empty() {
            0
        } else {
            cols[4].split(|&b| b == b',').count()
        };
        if n_alt > 1 {
            multiallelic_skipped += 1;
            return Ok(());
        }
        let gt_idx = format_index(cols.get(8).copied().unwrap_or(b"."), b"GT");
        for (i, slot) in site_codes.iter_mut().enumerate() {
            let gt = format_cell(cols.get(9 + i).copied().unwrap_or(b"."), gt_idx);
            *slot = match genotype_012(gt, n_alt) {
                Ok(c) => c,
                Err(e) if e.allele_out_of_range && o.skip_invalid => {
                    invalid_alleles += 1;
                    -1
                }
                Err(e) => {
                    return Err(IoError::Parse {
                        path: name.clone(),
                        loc: Loc::Record(hit.rec_no),
                        msg: format!(
                            "sample {} at {}:{}: {}",
                            String::from_utf8_lossy(&samples[i]),
                            String::from_utf8_lossy(cols[0]),
                            String::from_utf8_lossy(cols[1]),
                            e.msg
                        ),
                    });
                }
            };
        }
        store.push_site(&site_codes)?;
        pos.write_all(cols[0]).map_err(write_err)?;
        pos.write_all(b"\t").map_err(write_err)?;
        pos.write_all(cols[1]).map_err(write_err)?;
        pos.write_all(b"\n").map_err(write_err)?;
        sites_out += 1;
        Ok(())
    })?;
    if multiallelic_skipped > 0 {
        log::warn!(
            "{name}: the 012 matrix leaves out {multiallelic_skipped} multi-allelic site(s), as \
vcftools does; split them first (bcftools norm -m-) to include their alleles"
        );
    }
    if invalid_alleles > 0 {
        log::warn!(
            "{name}: {invalid_alleles} genotype(s) named an allele the site's ALT column lacks \
and were written as -1 (--skip-invalid)"
        );
    }
    let (spilled, row_blocks) = store.write_rows(samples.len(), o.max_mem, matrix)?;
    Ok(ExtractSummary {
        records_in,
        sites_out,
        multiallelic_skipped,
        spilled,
        row_blocks,
    })
}
