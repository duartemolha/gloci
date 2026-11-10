//! `gloci faidx`: regions of an indexed FASTA (D131).
use crate::Ctx;
use crate::cli::FaidxArgs;
use gloci_core::{Interval, Region};
use gloci_io::{
    Accept, IndexedFasta, InvalidPolicy, IoError, OutFormat, OutputOptions, ReaderOptions,
    UsageError, open_output, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::faidx::{FaidxOptions, FaidxRegion, MarkStrand, faidx};
use gloci_ops::read_regions_file;
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
    if let Some(f) = &a.region.regions_file {
        let ro = ReaderOptions {
            accept: Accept::default(),
            chrom_map: None,
            threads: ctx.threads,
            format: None,
            vcf_span: gloci_core::VcfSpan::default(),
            invalid: InvalidPolicy::Fail,
        };
        let mut points = 0u64;
        for (chrom, iv) in read_regions_file(&f.display().to_string(), &ro)? {
            if iv.start == iv.end {
                points += 1;
                continue;
            }
            let shown = String::from_utf8_lossy(&chrom);
            regions.push(FaidxRegion {
                label: format!("{shown}:{}-{}", iv.start + 1, iv.end),
                chrom,
                iv,
            });
        }
        if points > 0 {
            log::warn!(
                "{points} zero-length region{} in {} {} no bases and {} skipped",
                if points == 1 { "" } else { "s" },
                f.display(),
                if points == 1 { "has" } else { "have" },
                if points == 1 { "was" } else { "were" },
            );
        }
    }
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
    faidx(&mut fasta, &regions, &o, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
