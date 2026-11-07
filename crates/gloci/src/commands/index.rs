use crate::Ctx;
use crate::cli::IndexArgs;
use gloci_io::{InputOptions, open_input, write_index};

pub fn run(ctx: &Ctx, args: &IndexArgs) -> anyhow::Result<()> {
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
