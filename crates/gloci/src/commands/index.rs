use crate::Ctx;
use crate::cli::IndexArgs;
use gloci_io::{
    Format, InputOptions, IoError, TsvColumns, UsageError, open_input, warn_unused_reference,
    write_alignment_index, write_crai, write_fasta_index, write_index, write_tsv_index,
};

pub fn run(ctx: &Ctx, args: &IndexArgs) -> anyhow::Result<()> {
    if args.file.as_os_str() == "-" {
        return Err(UsageError(
            "index needs a file and cannot read standard input; write the data to a BGZF file \
first (gloci view - -o x.bed.gz), then run: gloci index x.bed.gz"
                .into(),
        )
        .into());
    }
    let csi_only = |what: &str| {
        UsageError(format!(
            "--csi is for BAM files, but {what}; {}",
            "tabix files get a .csi on their own when a position is beyond 2^29"
        ))
    };
    if let Some(cols) = tsv_columns(args)? {
        if args.csi {
            return Err(csi_only("this is a TSV index").into());
        }
        warn_unused_reference(args.reference.as_deref(), false);
        if std::fs::metadata(&args.file).is_ok_and(|m| m.len() == 0) {
            return Err(IoError::empty_input(&args.file.to_string_lossy()).into());
        }
        // write_tsv_index rejects non-BGZF input with the same recompress hint.
        write_tsv_index(&args.file, &args.file, &cols, ctx.threads)?;
        return Ok(());
    }
    let mut input = open_input(
        &args.file.to_string_lossy(),
        &InputOptions {
            threads: ctx.threads,
            format: None,
        },
    )?;
    // An empty file is detected as BED (an empty interval file is valid), which would then
    // be refused as not BGZF-compressed.
    if input.is_empty()? {
        return Err(IoError::empty_input(&input.name).into());
    }
    let format = input.format;
    drop(input);
    if format == Format::Cram {
        if args.csi {
            return Err(UsageError(format!(
                "--csi is for BAM files, but {} is CRAM, whose index is a .crai; drop --csi",
                args.file.display()
            ))
            .into());
        }
        write_crai(&args.file, args.reference.as_deref(), ctx.threads)?;
        return Ok(());
    }
    warn_unused_reference(args.reference.as_deref(), false);
    if format == Format::Bam {
        write_alignment_index(&args.file, &args.file, args.csi, ctx.threads)?;
        return Ok(());
    }
    if args.csi {
        return Err(csi_only(&format!("{} is {}", args.file.display(), format.name())).into());
    }
    if format == Format::Fasta {
        write_fasta_index(&args.file, ctx.threads)?;
        return Ok(());
    }
    // write_index rejects non-BGZF input with the "recompress with: gloci view" usage error.
    write_index(&args.file, &args.file, format, ctx.threads)?;
    Ok(())
}

/// The TSV columns when any of -s/-b/-e is given; a usage error for a
/// half-given set or for TSV-only flags without them.
fn tsv_columns(args: &IndexArgs) -> Result<Option<TsvColumns>, UsageError> {
    match (args.seq_col, args.begin_col, args.end_col) {
        (None, None, None) => {
            for (given, flag) in [
                (args.zero_based, "--zero-based"),
                (args.skip_lines.is_some(), "--skip-lines"),
                (args.comment.is_some(), "--comment"),
            ] {
                if given {
                    return Err(UsageError(format!(
                        "{flag} needs the columns of a TSV file: add -s SEQ -b BEGIN [-e END]"
                    )));
                }
            }
            Ok(None)
        }
        (Some(seq), Some(begin), end) => Ok(Some(TsvColumns {
            seq,
            begin,
            end,
            zero_based: args.zero_based,
            skip_lines: args.skip_lines.unwrap_or(0),
            comment: args.comment.unwrap_or(b'#'),
        })),
        _ => Err(UsageError(
            "-s (sequence column) and -b (begin column) must be given together, for example: \
gloci index FILE -s 1 -b 2 -e 3"
                .into(),
        )),
    }
}
