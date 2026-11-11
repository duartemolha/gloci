//! Allele comparison across records whose REF differs in length, a port of
//! bcftools' `vcmp.c`: one REF is a prefix of the other and the extra bases
//! are appended to the alleles of the shorter one before comparing.

#[derive(Default, Debug)]
pub struct Vcmp {
    /// Positive when REF 1 is longer, negative when REF 2 is.
    ndref: isize,
    dref: Vec<u8>,
}

fn up(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

impl Vcmp {
    /// `vcmp_set_ref`: `None` when the REFs are not compatible.
    pub fn set_ref(&mut self, ref1: &[u8], ref2: &[u8]) -> Option<()> {
        self.ndref = 0;
        let n = ref1
            .iter()
            .zip(ref2)
            .take_while(|(a, b)| up(**a) == up(**b))
            .count();
        let (r1, r2) = (&ref1[n..], &ref2[n..]);
        match (r1.is_empty(), r2.is_empty()) {
            (true, true) => Some(()),
            (false, false) => None,
            (false, true) => {
                self.dref = r1.iter().map(|&b| up(b)).collect();
                self.ndref = r1.len() as isize;
                Some(())
            }
            (true, false) => {
                self.dref = r2.iter().map(|&b| up(b)).collect();
                self.ndref = -(r2.len() as isize);
                Some(())
            }
        }
    }

    /// `vcmp_find_allele`: the index of `al2` among `als1`.
    pub fn find_allele<A: AsRef<[u8]>>(&self, als1: &[A], al2: &[u8]) -> Option<usize> {
        for (i, a) in als1.iter().enumerate() {
            let a = a.as_ref();
            let n = a
                .iter()
                .zip(al2)
                .take_while(|(x, y)| up(**x) == up(**y))
                .count();
            let (ra, rb) = (&a[n..], &al2[n..]);
            if !ra.is_empty() && !rb.is_empty() {
                continue;
            }
            if self.ndref == 0 {
                if ra.is_empty() && rb.is_empty() {
                    return Some(i);
                }
                continue;
            }
            if !ra.is_empty() {
                if self.ndref < 0 {
                    continue;
                }
                if ra.len() == self.ndref as usize
                    && ra.iter().zip(&self.dref).all(|(x, y)| up(*x) == *y)
                {
                    return Some(i);
                }
                continue;
            }
            if self.ndref > 0 {
                continue;
            }
            if rb.len() == (-self.ndref) as usize
                && rb.iter().zip(&self.dref).all(|(x, y)| up(*x) == *y)
            {
                return Some(i);
            }
        }
        None
    }

    /// `vcmp_map_ARvalues`: for each of the `n` values of the destination
    /// record, the index of the source value of the same allele. `als1` are
    /// the source alleles (REF first), `als2` the destination's.
    pub fn map_ar_values(
        &mut self,
        n: usize,
        als1: &[&[u8]],
        als2: &[&[u8]],
    ) -> Option<Vec<Option<usize>>> {
        self.set_ref(als1[0], als2[0])?;
        let ifrom = if n == als2.len() { 0 } else { 1 };
        Some(
            (ifrom..als2.len())
                .map(|i| self.find_allele(&als1[ifrom..], als2[i]))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_ref_and_longer_ref() {
        let mut v = Vcmp::default();
        assert!(v.set_ref(b"C", b"c").is_some());
        assert_eq!(v.find_allele(&[&b"T"[..], b"G"], b"g"), Some(1));
        assert_eq!(v.find_allele(&[&b"T"[..]], b"G"), None);
        // REF 1 = CA is longer than REF 2 = C: ALT T of record 2 is TA there.
        assert!(v.set_ref(b"CA", b"C").is_some());
        assert_eq!(v.find_allele(&[&b"TA"[..], b"G"], b"T"), Some(0));
        assert!(v.set_ref(b"CA", b"GA").is_none());
    }

    #[test]
    fn ar_mapping() {
        let mut v = Vcmp::default();
        // Source C>G,T; destination C>T,G: A values map (T,G) -> (1,0).
        let m = v
            .map_ar_values(2, &[b"C", b"G", b"T"], &[b"C", b"T", b"G"])
            .unwrap();
        assert_eq!(m, vec![Some(1), Some(0)]);
        // R values include REF.
        let m = v
            .map_ar_values(3, &[b"C", b"G", b"T"], &[b"C", b"T", b"G"])
            .unwrap();
        assert_eq!(m, vec![Some(0), Some(2), Some(1)]);
    }
}
