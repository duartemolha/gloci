//! A structural check of BAI and CSI bytes before noodles reads them: noodles sizes its
//! tables from the counts a file states, so a corrupt count would ask for gigabytes (and abort).
//! Every count is checked against the bytes that remain, so what noodles reads is no larger
//! than the file, and the shifts a query uses stay within what a 64-bit position can hold.

use std::io::{self, Read};

/// A cursor over the bytes of an index.
struct Cur<'a> {
    b: &'a [u8],
    at: usize,
}

impl Cur<'_> {
    fn left(&self) -> usize {
        self.b.len() - self.at
    }

    fn skip(&mut self, n: u64, what: &str) -> Result<(), String> {
        match usize::try_from(n) {
            Ok(n) if n <= self.left() => {
                self.at += n;
                Ok(())
            }
            _ => Err(format!("it ends inside {what}")),
        }
    }

    fn word<const N: usize>(&mut self, what: &str) -> Result<[u8; N], String> {
        if self.left() < N {
            return Err(format!("it ends inside {what}"));
        }
        let mut out = [0u8; N];
        out.copy_from_slice(&self.b[self.at..self.at + N]);
        self.skip(N as u64, what)?;
        Ok(out)
    }

    fn i32(&mut self, what: &str) -> Result<i32, String> {
        if self.left() < 4 {
            return Err(format!("it ends inside {what}"));
        }
        Ok(i32::from_le_bytes(self.word(what)?))
    }

    /// A count that must fit in the bytes that remain at `each` bytes per item.
    fn count(&mut self, what: &str, each: u64) -> Result<u64, String> {
        let n = u64::try_from(self.i32(what)?)
            .map_err(|_| format!("{what} is negative; the index is corrupt"))?;
        if n.saturating_mul(each) > self.left() as u64 {
            return Err(format!(
                "{what} is {n}, more than the rest of the file can hold; the index is corrupt"
            ));
        }
        Ok(n)
    }
}

/// The bins of one reference: each bin id at most `max_bin`, and every chunk offset (except
/// those of the pseudo-bin 37450 holding counts) within the `data_len` bytes of the data file.
fn check_bins(
    c: &mut Cur,
    n_bin: u64,
    max_bin: u64,
    csi: bool,
    data_len: u64,
) -> Result<(), String> {
    for _ in 0..n_bin {
        let bin = u64::from(u32::from_le_bytes(c.word("a bin")?));
        if bin > max_bin {
            return Err(format!(
                "a bin is numbered {bin}, above the largest, {max_bin}; the index is corrupt"
            ));
        }
        if csi {
            c.skip(8, "a bin")?;
        }
        let n_chunk = c.count("a chunk count", 16)?;
        for _ in 0..n_chunk {
            let (b, e) = (
                u64::from_le_bytes(c.word("a chunk")?),
                u64::from_le_bytes(c.word("a chunk")?),
            );
            if bin != 37450 && ((b >> 16) > data_len || (e >> 16) > data_len) {
                return Err(format!(
                    "a chunk lies beyond the end of the data file ({data_len} bytes); the index \
belongs to another file or is corrupt"
                ));
            }
        }
    }
    Ok(())
}

fn check_linear(c: &mut Cur, data_len: u64) -> Result<(), String> {
    let n_intv = c.count("an interval count", 8)?;
    for _ in 0..n_intv {
        if (u64::from_le_bytes(c.word("a linear index")?) >> 16) > data_len {
            return Err(format!(
                "the linear index points beyond the end of the data file ({data_len} bytes); the \
index belongs to another file or is corrupt"
            ));
        }
    }
    Ok(())
}

/// Check a BAI (`.bai`, uncompressed) of a data file of `data_len` bytes with `n_ref_header`
/// references.
pub(super) fn check_bai(
    bytes: &[u8],
    data_len: u64,
    n_ref_header: Option<usize>,
) -> Result<(), String> {
    let mut c = Cur { b: bytes, at: 0 };
    if bytes.get(..4) != Some(b"BAI\x01") {
        return Err("it is not a BAM index (no BAI magic bytes)".into());
    }
    c.at = 4;
    // each reference has at least its bin count and its interval count
    let n_ref = c.count("the reference count", 8)?;
    if let Some(h) = n_ref_header
        && n_ref != h as u64
    {
        return Err(format!(
            "it lists {n_ref} references but the BAM header has {h}; it belongs to another file"
        ));
    }
    for _ in 0..n_ref {
        let n_bin = c.count("a bin count", 8)?;
        check_bins(&mut c, n_bin, 37450, false, data_len)?;
        check_linear(&mut c, data_len)?;
    }
    Ok(())
}

/// Check a CSI (`.csi`, decompressed). `n_ref_header` is the number of references the data
/// file's header lists, when it has a list (BAM); a BCF's or tabix file's index may list fewer.
pub(super) fn check_csi(
    bytes: &[u8],
    data_len: u64,
    n_ref_header: Option<usize>,
) -> Result<(), String> {
    let mut c = Cur { b: bytes, at: 0 };
    if bytes.get(..4) != Some(b"CSI\x01") {
        return Err("it is not a CSI index (no CSI magic bytes)".into());
    }
    c.at = 4;
    let min_shift = c.i32("the minimum shift")?;
    let depth = c.i32("the depth")?;
    if !(1..=24).contains(&min_shift) || !(0..=9).contains(&depth) || min_shift + 3 * depth > 62 {
        return Err(format!(
            "its bins (minimum shift {min_shift}, depth {depth}) are not ones a position can \
use; the index is corrupt"
        ));
    }
    let l_aux = c.count("the length of the auxiliary data", 1)?;
    c.skip(l_aux, "the auxiliary data")?;
    let n_ref = c.count("the reference count", 4)?;
    if let Some(h) = n_ref_header
        && n_ref != h as u64
    {
        return Err(format!(
            "it lists {n_ref} references but the BAM header has {h}; it belongs to another file"
        ));
    }
    let max_bin = ((1u64 << (3 * (depth + 1))) - 1) / 7 + 1;
    for _ in 0..n_ref {
        let n_bin = c.count("a bin count", 16)?;
        check_bins(&mut c, n_bin, max_bin, true, data_len)?;
    }
    Ok(())
}

/// Check a TBI (`.tbi`, decompressed): the layout tabix writes, with every count held to the
/// bytes that remain.
pub(super) fn check_tbi(bytes: &[u8], data_len: u64) -> Result<(), String> {
    let mut c = Cur { b: bytes, at: 0 };
    if bytes.get(..4) != Some(b"TBI\x01") {
        return Err("it is not a tabix index (no TBI magic bytes)".into());
    }
    c.at = 4;
    // each reference has at least its bin count and its interval count
    let n_ref = c.count("the reference count", 8)?;
    // format, three columns, meta char, skip lines
    c.skip(6 * 4, "the header")?;
    let l_nm = c.count("the length of the names", 1)?;
    c.skip(l_nm, "the sequence names")?;
    for _ in 0..n_ref {
        let n_bin = c.count("a bin count", 8)?;
        check_bins(&mut c, n_bin, 37450, false, data_len)?;
        check_linear(&mut c, data_len)?;
    }
    Ok(())
}

/// Decompress a BGZF file into memory, at most `cap` bytes.
pub(super) fn read_bgzf_capped(r: impl Read, cap: u64) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(r)
        .take(cap)
        .read_to_end(&mut out)?;
    Ok(out)
}
