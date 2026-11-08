//! Numeric parsing and formatting for aggregation values.

/// A parsed number: integer literals stay exact.
#[derive(Debug, Clone, Copy)]
pub enum Num {
    Int(i128),
    Float(f64),
}

impl Num {
    pub fn as_f64(self) -> f64 {
        match self {
            Num::Int(i) => i as f64,
            Num::Float(f) => f,
        }
    }
}

/// Parse a number: an integer-literal fast path, else a finite `f64`.
pub fn parse_num(s: &[u8]) -> Option<Num> {
    let digits = match s.first()? {
        b'-' | b'+' => &s[1..],
        _ => s,
    };
    if !digits.is_empty() && digits.len() <= 36 && digits.iter().all(u8::is_ascii_digit) {
        let mut v: i128 = 0;
        for &d in digits {
            v = v * 10 + i128::from(d - b'0');
        }
        return Some(Num::Int(if s[0] == b'-' { -v } else { v }));
    }
    let f: f64 = std::str::from_utf8(s).ok()?.parse().ok()?;
    f.is_finite().then_some(Num::Float(f))
}

/// Rust's shortest round-trip form.
pub fn fmt_f64(f: f64, out: &mut Vec<u8>) {
    use std::io::Write;
    let _ = write!(out, "{f}");
}
