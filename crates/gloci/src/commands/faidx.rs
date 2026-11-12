//! `gloci faidx`: regions of an indexed FASTA (D131).
use crate::Ctx;
use crate::cli::FaidxArgs;
use gloci_core::{Interval, Region};
use gloci_io::{
    Accept, IndexedFasta, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions,
    UsageError, open_output, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::faidx::{FaidxOptions, FaidxRegion, FaidxWriter, MarkStrand, RegionsFile};
use std::io::IsTerminal;

/// The most unknown names a warning lists.
const SHOWN: usize = 5;

pub fn run(ctx: &Ctx, a: &FaidxArgs) -> anyhow::Result<()> {
    let fasta_name = a.fasta.display().to_string();
    let typed: Vec<&String> = a.regions.iter().chain(&a.region.region).collect();
    if typed.is_empty() && a.region.regions_file.is_none() {
        return Err(UsageError(format!(
            "no region given; name one or more regions (for example chr1:1,000-2,000), or pass \
--regions-file BED. To index {fasta_name} instead, run: gloci index {fasta_name} \
(samtools faidx builds the index when no region is given; gloci does not)"
        ))
        .into());
    }
    let mark: MarkStrand = a.mark_strand.parse().map_err(UsageError)?;
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    if let Some(f) = spec.format
        && f != OutFormat::Fasta
    {
        return Err(UsageError(format!(
            "faidx writes FASTA but the output format is {}; use -o FILE.fa (or .fa.gz), or -O fa \
for standard output",
            f.name()
        ))
        .into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Fasta,
        spec.bgzf,
    )?;
    let mut fasta = IndexedFasta::open(&a.fasta)?;

    // Every region is parsed and checked before the output is created.
    let mut regions: Vec<FaidxRegion> = Vec::with_capacity(typed.len());
    for t in typed {
        let r = Region::parse(t, |n| fasta.len_of(n).is_some())
            .map_err(|e| UsageError(format!("{e}; regions look like chr1:1,000-2,000")))?;
        regions.push(FaidxRegion {
            label: t.clone(),
            chrom: r.chrom,
            iv: Interval {
                start: r.start,
                end: r.end.unwrap_or(u64::MAX).max(r.start),
            },
        });
    }
    // A regions file is scanned for unknown sequences and zero-length lines
    // without keeping its regions, then streamed into the output below.
    let ro = ReaderOptions {
        accept: Accept::default(),
        chrom_map: None,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::default(),
        invalid: InvalidPolicy::Fail,
    };
    let bed = match &a.region.regions_file {
        Some(f) => Some((f, RegionsFile::scan(&f.display().to_string(), &ro, &fasta)?)),
        None => None,
    };
    let missing: Vec<&FaidxRegion> = regions
        .iter()
        .filter(|r| fasta.len_of(&r.chrom).is_none())
        .collect();
    if let Some(first) = missing.first() {
        let shown = String::from_utf8_lossy(&first.chrom);
        if !a.invalid.skip_invalid {
            return Err(IoError::Detect {
                path: fasta_name.clone(),
                msg: format!(
                    "{fasta_name}: sequence {shown} (region {}) is not in the index {fasta_name}.fai{}; \
check the sequence name, or pass --skip-invalid to skip regions on missing sequences",
                    first.label,
                    if missing.len() > 1 {
                        format!(" ({} regions name missing sequences)", missing.len())
                    } else {
                        String::new()
                    }
                ),
            }
            .into());
        }
        let mut names: Vec<String> = Vec::new();
        for r in &missing {
            let n = String::from_utf8_lossy(&r.chrom).into_owned();
            if names.len() < SHOWN && !names.contains(&n) {
                names.push(n);
            }
        }
        log::warn!(
            "skipped {} region{} on sequences not in {fasta_name}.fai (first: {}); check the \
sequence names",
            missing.len(),
            if missing.len() == 1 { "" } else { "s" },
            names.join(", ")
        );
        regions.retain(|r| fasta.len_of(&r.chrom).is_some());
    }

    if let Some((f, rf)) = &bed {
        if rf.zero_length > 0 {
            let n = rf.zero_length;
            log::warn!(
                "{n} zero-length region{} in {} {} no bases and {} skipped",
                if n == 1 { "" } else { "s" },
                f.display(),
                if n == 1 { "has" } else { "have" },
                if n == 1 { "was" } else { "were" },
            );
        }
        if let Some(first) = &rf.first_missing {
            if !a.invalid.skip_invalid {
                return Err(IoError::Detect {
                    path: fasta_name.clone(),
                    msg: format!(
                        "{fasta_name}: region {first} in {} is on a sequence not in the index \
{fasta_name}.fai ({} regions of the file); check the sequence names, or pass --skip-invalid to \
skip regions on missing sequences",
                        f.display(),
                        rf.missing
                    ),
                }
                .into());
            }
            log::warn!(
                "skipped {} region{} of {} on sequences not in {fasta_name}.fai (first: {}); \
check the sequence names",
                rf.missing,
                if rf.missing == 1 { "" } else { "s" },
                f.display(),
                rf.missing_names.join(", ")
            );
        }
    }

    let o = FaidxOptions {
        line_width: a.line_width,
        reverse_complement: a.reverse_complement,
        mark,
    };
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Fasta, &oo)?;
    let mut w = FaidxWriter::new(&o);
    let mut write_all = || -> Result<(), IoError> {
        for r in &regions {
            w.write(&mut fasta, r, &mut out)?;
        }
        if let Some((_, rf)) = &bed {
            rf.write(&mut fasta, &ro, a.invalid.skip_invalid, &mut w, &mut out)?;
        }
        Ok(())
    };
    write_all().map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
