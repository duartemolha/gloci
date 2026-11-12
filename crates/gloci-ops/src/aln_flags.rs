//! `gloci flags` (D136): SAM FLAGS as numbers and names, as `samtools flags` writes them.

use gloci_io::{SAM_FLAGS, UsageError, flag_names, parse_flags};
use std::fmt::Write;

/// The line for one FLAGS argument: `0x<hex> TAB <decimal> TAB <names>`, with a newline.
pub fn flags_line(arg: &str) -> Result<String, UsageError> {
    let flag = parse_flags(arg).map_err(|e| UsageError(format!("flags: {e}")))?;
    Ok(format!("0x{flag:x}\t{flag}\t{}\n", flag_names(flag)))
}

/// samtools' table of the twelve flags, one line each: `%#6x %5d  %-15s%s`.
pub fn flags_table() -> String {
    let mut out = String::new();
    for (bit, name, desc) in SAM_FLAGS {
        // `{:#6x}` pads after the 0x prefix as C's `%#6x` does: the whole width is 6.
        let _ = writeln!(out, "{:>6} {bit:>5}  {name:<15}{desc}", format!("{bit:#x}"));
    }
    out
}
