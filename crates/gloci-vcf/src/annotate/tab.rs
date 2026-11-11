//! Rows of a tab-delimited or BED annotation file.

use std::ops::Range;

/// Where the fields are in the annotation file's columns (0-based indexes).
#[derive(Clone, Debug)]
pub struct TabLayout {
    pub chr: usize,
    /// The start column (`POS`, `FROM` or `BEG`).
    pub beg: usize,
    /// The end column (`TO`/`END`): an interval file, matched by overlap
    /// (otherwise a row is the single position `beg`).
    pub end: Option<usize>,
    pub ref_: Option<usize>,
    pub alt: Option<usize>,
    /// A `.bed` file: 0-based starts, exclusive ends.
    pub bed: bool,
}

/// One annotation row with its footprint.
#[derive(Clone, Debug)]
pub struct TabRow {
    pub line: Vec<u8>,
    cols: Vec<Range<usize>>,
    chrom: Range<usize>,
    /// 0-based inclusive footprint.
    pub start: i64,
    pub end: i64,
}

fn parse_decimal(s: &[u8]) -> Option<i64> {
    if s.is_empty() || !s.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(s).ok()?.parse().ok()
}

impl TabRow {
    /// Parse `line`; `Ok(None)` for a comment or blank line.
    pub fn parse(line: &[u8], layout: &TabLayout) -> Result<Option<TabRow>, String> {
        if line.is_empty() || line[0] == b'#' {
            return Ok(None);
        }
        let mut cols = Vec::new();
        let mut s = 0;
        for (i, &b) in line.iter().enumerate() {
            if b == b'\t' {
                cols.push(s..i);
                s = i + 1;
            }
        }
        cols.push(s..line.len());
        let need = [
            Some(layout.chr),
            Some(layout.beg),
            layout.end,
            layout.ref_,
            layout.alt,
        ]
        .into_iter()
        .flatten()
        .max()
        .unwrap_or(0)
            + 1;
        if cols.len() < need {
            return Err(format!(
                "the line has {} tab-separated columns but --columns needs {need}",
                cols.len()
            ));
        }
        let chrom = cols[layout.chr].clone();
        let from_s = &line[cols[layout.beg].clone()];
        let from = parse_decimal(from_s).ok_or_else(|| {
            format!(
                "column {} is `{}`, expected a position",
                layout.beg + 1,
                String::from_utf8_lossy(from_s)
            )
        })?;
        let (start, end) = match layout.end {
            None => {
                if from == 0 {
                    return Err(format!(
                        "column {} is 0, expected a 1-based position",
                        layout.beg + 1
                    ));
                }
                (from - 1, from - 1)
            }
            Some(e) => {
                let to_s = &line[cols[e].clone()];
                let to = parse_decimal(to_s);
                if layout.bed {
                    let to = to.ok_or_else(|| {
                        format!(
                            "column {} is `{}`, expected an end",
                            e + 1,
                            String::from_utf8_lossy(to_s)
                        )
                    })?;
                    (from, to - 1)
                } else {
                    if from == 0 || to == Some(0) {
                        return Err("positions are 1-based; a 0 start or end is not valid (use a .bed file for 0-based starts)".to_string());
                    }
                    (from - 1, to.map_or(from - 1, |t| t - 1))
                }
            }
        };
        Ok(Some(TabRow {
            line: line.to_vec(),
            cols,
            chrom,
            start,
            end,
        }))
    }

    pub fn chrom(&self) -> &[u8] {
        &self.line[self.chrom.clone()]
    }

    pub fn col(&self, i: usize) -> Option<&[u8]> {
        self.cols.get(i).map(|r| &self.line[r.clone()])
    }

    /// REF then the comma-split ALT column, with the number of alleles.
    pub fn alleles(&self, layout: &TabLayout) -> Option<(Vec<Vec<u8>>, usize)> {
        let r = self.col(layout.ref_?)?;
        let a = self.col(layout.alt?)?;
        let mut v = vec![r.to_vec()];
        v.extend(a.split(|&b| b == b',').map(<[u8]>::to_vec));
        let n = v.len();
        Some((v, n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn layout(end: Option<usize>, bed: bool) -> TabLayout {
        TabLayout {
            chr: 0,
            beg: 1,
            end,
            ref_: None,
            alt: None,
            bed,
        }
    }

    #[test]
    fn footprints() {
        let r = TabRow::parse(b"1\t100\tx", &layout(None, false))
            .unwrap()
            .unwrap();
        assert_eq!((r.start, r.end, r.chrom()), (99, 99, &b"1"[..]));
        let r = TabRow::parse(b"1\t100\t200\tx", &layout(Some(2), false))
            .unwrap()
            .unwrap();
        assert_eq!((r.start, r.end), (99, 199));
        let r = TabRow::parse(b"1\t100\t200\tx", &layout(Some(2), true))
            .unwrap()
            .unwrap();
        assert_eq!((r.start, r.end), (100, 199));
        assert!(
            TabRow::parse(b"#c", &layout(None, false))
                .unwrap()
                .is_none()
        );
        assert!(TabRow::parse(b"1\tx\ty", &layout(None, false)).is_err());
        assert!(TabRow::parse(b"1", &layout(None, false)).is_err());
    }
}
