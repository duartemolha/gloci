//! Lazy access to the tab-separated columns of a raw record line.

/// The 0-based tab-separated field `n` of `raw`, if the line has one.
pub fn field(raw: &[u8], n: usize) -> Option<&[u8]> {
    split_fields(raw).nth(n)
}

/// The number of tab-separated fields (an empty line has one).
pub fn column_count(raw: &[u8]) -> usize {
    memchr::memchr_iter(b'\t', raw).count() + 1
}

/// The tab-separated fields of `raw`, in order.
pub fn split_fields(raw: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut rest = Some(raw);
    std::iter::from_fn(move || {
        let cur = rest?;
        match memchr::memchr(b'\t', cur) {
            Some(i) => {
                rest = Some(&cur[i + 1..]);
                Some(&cur[..i])
            }
            None => {
                rest = None;
                Some(cur)
            }
        }
    })
}
