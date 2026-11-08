use crate::Ctx;
use crate::cli::IndexArgs;
use gloci_io::{InputOptions, UsageError, open_input, write_index};

pub fn run(ctx: &Ctx, args: &IndexArgs) -> anyhow::Result<()> {
    if args.file.as_os_str() == "-" {
        return Err(UsageError(
            "index needs a file and cannot read standard input; write the data to a BGZF file \
first (gloci view - -o x.bed.gz), then run: gloci index x.bed.gz"
                .into(),
        )
        .into());
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
    // write_index rejects non-BGZF input with the "recompress with: gloci view" usage error.
    write_index(&args.file, &args.file, format, ctx.threads)?;
    Ok(())
}
