use crate::Ctx;
use crate::cli::IndexArgs;
use gloci_io::{
    Format, InputOptions, TsvColumns, UsageError, open_input, write_fasta_index, write_index,
    write_tsv_index,
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
    if let Some(cols) = tsv_columns(args)? {
        // write_tsv_index rejects non-BGZF input with the same recompress hint.
        write_tsv_index(&args.file, &args.file, &cols, ctx.threads)?;
        return Ok(());
    }
    let input = open_input(
        &args.file.to_string_lossy(),
        &InputOptions {
            threads: ctx.threads,
            format: None,
        },
    )?;
    let format = input.format;
    drop(input);
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
