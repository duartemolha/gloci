//! Allocation-free C-style number formatting for the `vcf stats` tables.
//!
//! vcftools prints through C `printf`/iostreams: `%g` (6 digits) for most
//! values, `%e` for the HWE columns, `%.2f` for expected counts, and `nan` /
//! `inf` in lower case (never `-nan`).

use std::fmt::Write as _;
use std::io::Write as _;

/// A stack buffer that implements `fmt::Write`; overflow is an error, which
/// no f64 formatted here can reach.
struct StackBuf {
    b: [u8; 48],
    n: usize,
}

impl StackBuf {
    fn new() -> StackBuf {
        StackBuf { b: [0; 48], n: 0 }
    }
    fn bytes(&self) -> &[u8] {
        &self.b[..self.n]
    }
}

impl std::fmt::Write for StackBuf {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        let e = self.n + s.len();
        if e > self.b.len() {
            return Err(std::fmt::Error);
        }
        self.b[self.n..e].copy_from_slice(s.as_bytes());
        self.n = e;
        Ok(())
    }
}

fn special(out: &mut Vec<u8>, x: f64) -> bool {
    if x.is_nan() {
        out.extend_from_slice(b"nan");
    } else if x.is_infinite() {
        out.extend_from_slice(if x < 0.0 { b"-inf" } else { b"inf" });
    } else {
        return false;
    }
    true
}

/// Split `{:.Ne}` output into (negative, digits without the point, exponent).
fn sci(buf: &StackBuf) -> (bool, [u8; 16], usize, i32) {
    let s = buf.bytes();
    let (neg, s) = match s.split_first() {
        Some((b'-', r)) => (true, r),
        _ => (false, s),
    };
    let e = s.iter().position(|&c| c == b'e').expect("exponent");
    let mut digits = [0u8; 16];
    let mut n = 0;
    for &c in &s[..e] {
        if c != b'.' {
            digits[n] = c;
            n += 1;
        }
    }
    let exp: i32 = std::str::from_utf8(&s[e + 1..])
        .expect("ascii")
        .parse()
        .expect("exponent");
    (neg, digits, n, exp)
}

fn push_exp(out: &mut Vec<u8>, exp: i32) {
    out.push(b'e');
    out.push(if exp < 0 { b'-' } else { b'+' });
    let a = exp.unsigned_abs();
    if a < 10 {
        out.push(b'0');
    }
    let mut ib = itoa::Buffer::new();
    out.extend_from_slice(ib.format(a).as_bytes());
}

/// `printf("%g", x)`: 6 significant digits, trailing zeros trimmed, an
/// exponent when it is below -4 or at least 6.
pub fn write_g6(out: &mut Vec<u8>, x: f64) {
    if special(out, x) {
        return;
    }
    if x == 0.0 {
        out.extend_from_slice(if x.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    let mut buf = StackBuf::new();
    write!(buf, "{x:.5e}").expect("fits");
    let (neg, digits, n, exp) = sci(&buf);
    let digits = &digits[..n];
    if neg {
        out.push(b'-');
    }
    let trimmed = {
        let mut t = digits.len();
        while t > 1 && digits[t - 1] == b'0' {
            t -= 1;
        }
        &digits[..t]
    };
    if (-4..6).contains(&exp) {
        if exp >= 0 {
            let int_len = exp as usize + 1;
            out.extend_from_slice(&digits[..int_len]);
            let frac = &digits[int_len..];
            let ft = frac.iter().rposition(|&c| c != b'0').map_or(0, |p| p + 1);
            if ft > 0 {
                out.push(b'.');
                out.extend_from_slice(&frac[..ft]);
            }
        } else {
            out.extend_from_slice(b"0.");
            for _ in 0..(-exp - 1) {
                out.push(b'0');
            }
            out.extend_from_slice(trimmed);
        }
    } else {
        out.push(digits[0]);
        if trimmed.len() > 1 {
            out.push(b'.');
            out.extend_from_slice(&trimmed[1..]);
        }
        push_exp(out, exp);
    }
}

/// `printf("%e", x)`: six decimals and a two-digit exponent.
pub fn write_e6(out: &mut Vec<u8>, x: f64) {
    if special(out, x) {
        return;
    }
    let mut buf = StackBuf::new();
    write!(buf, "{x:.6e}").expect("fits");
    let (neg, digits, n, exp) = sci(&buf);
    if neg {
        out.push(b'-');
    }
    out.push(digits[0]);
    out.push(b'.');
    out.extend_from_slice(&digits[1..n]);
    push_exp(out, exp);
}

/// `printf("%.2f", x)`.
pub fn write_f2(out: &mut Vec<u8>, x: f64) {
    if special(out, x) {
        return;
    }
    write!(out, "{x:.2}").expect("vec write");
}

/// `printf("%.1f", x)`.
pub fn write_f1(out: &mut Vec<u8>, x: f64) {
    if special(out, x) {
        return;
    }
    write!(out, "{x:.1}").expect("vec write");
}

/// `printf("%.5f", x)`.
pub fn write_f5(out: &mut Vec<u8>, x: f64) {
    if special(out, x) {
        return;
    }
    write!(out, "{x:.5}").expect("vec write");
}

pub fn write_u64(out: &mut Vec<u8>, v: u64) {
    let mut ib = itoa::Buffer::new();
    out.extend_from_slice(ib.format(v).as_bytes());
}

pub fn write_i64(out: &mut Vec<u8>, v: i64) {
    let mut ib = itoa::Buffer::new();
    out.extend_from_slice(ib.format(v).as_bytes());
}
