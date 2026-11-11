//! SAM FLAGS as numbers and names (D136): the table of the twelve flag bits, `flag_names`, and
//! the parser every FLAGS option shares (`samtools view -f` and `samtools flags` read the same
//! syntax, with the two differences D136 lists).

/// The twelve SAM flag bits: the bit, its name and samtools' description of it.
pub const SAM_FLAGS: [(u16, &str, &str); 12] = [
    (0x1, "PAIRED", "paired-end / multiple-segment sequencing technology"),
    (0x2, "PROPER_PAIR", "each segment properly aligned according to aligner"),
    (0x4, "UNMAP", "segment unmapped"),
    (0x8, "MUNMAP", "next segment in the template unmapped"),
    (0x10, "REVERSE", "SEQ is reverse complemented"),
    (0x20, "MREVERSE", "SEQ of next segment in template is rev.complemented"),
    (0x40, "READ1", "the first segment in the template"),
    (0x80, "READ2", "the last segment in the template"),
    (0x100, "SECONDARY", "secondary alignment"),
    (0x200, "QCFAIL", "not passing quality controls or other filters"),
    (0x400, "DUP", "PCR or optical duplicate"),
    (0x800, "SUPPLEMENTARY", "supplementary alignment"),
];

/// The names of the bits set in `flag`, comma-joined in bit order (empty when none is set;
/// bits above 0x800 have no name).
pub fn flag_names(flag: u16) -> String {
    SAM_FLAGS
        .iter()
        .filter(|(bit, _, _)| flag & bit != 0)
        .map(|(_, name, _)| *name)
        .collect::<Vec<_>>()
        .join(",")
}

/// Parse a SAM flag set: a decimal or `0x` hexadecimal number up to 65535, or comma-separated
/// names (`PAIRED`, `UNMAP`, `DUP`, ...), case-insensitive.
///
/// Where samtools differs, gloci refuses (D136): a decimal with a leading zero (`010`, which
/// samtools reads as octal 8) and trailing text (`12abc`, which samtools reads as 12).
pub fn parse_flags(s: &str) -> Result<u16, String> {
    let t = s.trim();
    let names = || {
        SAM_FLAGS
            .iter()
            .map(|(_, n, _)| *n)
            .collect::<Vec<_>>()
            .join(", ")
    };
    if t.is_empty() {
        return Err(format!(
            "empty flag value; give a number such as 0x704 or names such as secondary,dup \
(names: {})",
            names()
        ));
    }
    let hex = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X"));
    let number = if let Some(h) = hex {
        Some(if h.bytes().all(|b| b.is_ascii_hexdigit()) && !h.is_empty() {
            u32::from_str_radix(h, 16).ok()
        } else {
            None
        })
    } else if t.bytes().all(|b| b.is_ascii_digit()) {
        if t.len() > 1 && t.starts_with('0') {
            // samtools reads this as octal; the same digits would mean something else in gloci.
            let octal = u32::from_str_radix(t, 8).ok();
            let hint = match octal {
                Some(o) => format!(
                    "samtools reads {t} as octal {o}; write {o} (or 0x{o:x}) for that, or {} \
for decimal {t}",
                    t.trim_start_matches('0')
                ),
                None => format!(
                    "write the number without leading zeros ({}) or as 0x hex",
                    t.trim_start_matches('0')
                ),
            };
            return Err(format!(
                "flag value '{t}' has a leading zero, which samtools reads as octal and gloci \
does not guess at; {hint}"
            ));
        }
        Some(t.parse::<u32>().ok())
    } else if t.starts_with(|c: char| c.is_ascii_digit() || c == '-') {
        // Trailing text (`12abc`) or a sign (`-1`): not a number.
        Some(None)
    } else {
        None
    };
    if let Some(n) = number {
        return n
            .and_then(|n| u16::try_from(n).ok())
            .ok_or_else(|| format!("flag value '{t}' is not a number from 0 to 65535"));
    }
    let mut out = 0u16;
    for part in t.split(',') {
        let p = part.trim();
        match SAM_FLAGS.iter().find(|(_, n, _)| n.eq_ignore_ascii_case(p)) {
            Some((bit, _, _)) => out |= bit,
            None => {
                return Err(format!(
                    "unknown flag name '{p}'; use a number or names from: {}",
                    names()
                ));
            }
        }
    }
    Ok(out)
}
