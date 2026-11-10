//! Allocation-free C-style number formatting (D145).
//!
//! The reference tools print through C `printf`/iostreams: `%g` (6 digits)
//! for most values, `%e` for the HWE columns, `%.Nf` for expected counts, and
//! `nan` / `inf` in lower case (never `-nan`). `write_fixed` relies on Rust
//! rounding exact ties to even, as glibc does.

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
fn sci(buf: &StackBuf) -> (bool, [u8; 24], usize, i32) {
    let s = buf.bytes();
    let (neg, s) = match s.split_first() {
        Some((b'-', r)) => (true, r),
        _ => (false, s),
    };
    let e = s.iter().position(|&c| c == b'e').expect("exponent");
    let mut digits = [0u8; 24];
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
    write_g(out, x, 6);
}

/// `printf("%.{prec}g", x)`: `prec` significant digits (1 to 17; `%.0g` is
/// `%.1g`), trailing zeros trimmed, an exponent when it is below -4 or at
/// least `prec`.
pub fn write_g(out: &mut Vec<u8>, x: f64, prec: usize) {
    let prec = prec.clamp(1, 17);
    if special(out, x) {
        return;
    }
    if x == 0.0 {
        out.extend_from_slice(if x.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    let mut buf = StackBuf::new();
    write!(buf, "{x:.*e}", prec - 1).expect("fits");
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
    if (-4..prec as i32).contains(&exp) {
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

/// `printf("%.{prec}f", x)`.
pub fn write_fixed(out: &mut Vec<u8>, x: f64, prec: usize) {
    if special(out, x) {
        return;
    }
    write!(out, "{x:.prec$}").expect("vec write");
}

pub fn write_u64(out: &mut Vec<u8>, v: u64) {
    let mut ib = itoa::Buffer::new();
    out.extend_from_slice(ib.format(v).as_bytes());
}

pub fn write_i64(out: &mut Vec<u8>, v: i64) {
    let mut ib = itoa::Buffer::new();
    out.extend_from_slice(ib.format(v).as_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(x: f64, p: usize) -> String {
        let mut v = Vec::new();
        write_g(&mut v, x, p);
        String::from_utf8(v).unwrap()
    }

    fn fixed(x: f64, p: usize) -> String {
        let mut v = Vec::new();
        write_fixed(&mut v, x, p);
        String::from_utf8(v).unwrap()
    }

    #[test]
    fn cfmt_table() {
        assert_eq!(g(1.234567e-4, 5), "0.00012346");
        assert_eq!(g(123456789.0, 5), "1.2346e+08");
        assert_eq!(g(1.0, 5), "1");
        assert_eq!(g(2.5e-5, 6), "2.5e-05");
        assert_eq!(g(f64::NAN, 6), "nan");
        assert_eq!(g(f64::NEG_INFINITY, 6), "-inf");
        assert_eq!(g(0.0, 6), "0");
        assert_eq!(g(100000.0, 6), "100000");
        assert_eq!(g(1000000.0, 6), "1e+06");
        assert_eq!(g(-0.5, 0), "-0.5");
        assert_eq!(g(0.1 + 0.2, 17), "0.30000000000000004");
        assert_eq!(fixed(0.125, 2), "0.12");
        assert_eq!(fixed(0.375, 2), "0.38");
        assert_eq!(fixed(f64::INFINITY, 3), "inf");
        assert_eq!(fixed(2.0, 0), "2");
    }

    #[test]
    fn g6_matches_g_with_six() {
        for x in [0.0, 1.5, 123456.7, 1.0e-5, 99999.95, -3.25e10] {
            let mut a = Vec::new();
            write_g6(&mut a, x);
            assert_eq!(String::from_utf8(a).unwrap(), g(x, 6));
        }
    }
}
