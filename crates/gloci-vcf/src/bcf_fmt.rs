//! BCF genotype decoder: walks the FORMAT (individual) section of a framed
//! record (BCF2.2 §6.3.4).

use crate::gt::{GenotypeSource, GtFields, SampleFields, ft_passes};
use gloci_io::BcfDicts;
use gloci_io::bcf::typed::{self, Int, Ty};

/// One FORMAT field: per-sample element count and the data of all samples.
#[derive(Clone, Copy)]
struct Col<'a> {
    ty: Ty,
    n: usize,
    data: &'a [u8],
}

impl<'a> Col<'a> {
    /// The typed vector of sample `i`.
    fn of(&self, i: usize) -> typed::Value<'a> {
        let w = self.n * self.ty.width();
        typed::Value {
            ty: self.ty,
            n: self.n,
            data: &self.data[i * w..(i + 1) * w],
        }
    }
}

pub struct BcfGenotypes<'a> {
    n_samples: usize,
    gt: Option<Col<'a>>,
    dp: Option<Col<'a>>,
    gq: Option<Col<'a>>,
    ft: Option<Col<'a>>,
}

/// Index the FORMAT fields of one framed BCF record (`l_shared`, `l_indiv`,
/// shared, individual), as `RecordBatch::raw` returns it.
pub fn bcf_genotypes<'a>(raw: &'a [u8], dicts: &'a BcfDicts) -> Result<BcfGenotypes<'a>, String> {
    bcf_genotypes_with(raw, dicts, GtFields::ALL)
}

/// [`bcf_genotypes`] reading only GT and the fields `fields` names.
pub fn bcf_genotypes_with<'a>(
    raw: &'a [u8],
    dicts: &'a BcfDicts,
    fields: GtFields,
) -> Result<BcfGenotypes<'a>, String> {
    const TRUNC: &str = "the record is shorter than its header says";
    if raw.len() < 8 {
        return Err(TRUNC.into());
    }
    let l_shared = u32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]) as usize;
    let l_indiv = u32::from_le_bytes([raw[4], raw[5], raw[6], raw[7]]) as usize;
    let shared = raw.get(8..8 + l_shared).ok_or(TRUNC)?;
    let mut indiv = raw.get(8 + l_shared..8 + l_shared + l_indiv).ok_or(TRUNC)?;
    if shared.len() < 24 {
        return Err(format!(
            "the shared section is {} bytes; a BCF record needs at least 24",
            shared.len()
        ));
    }
    let nfs = u32::from_le_bytes([shared[20], shared[21], shared[22], shared[23]]);
    let (n_fmt, n_samples) = ((nfs >> 24) as usize, (nfs & 0x00ff_ffff) as usize);
    let mut out = BcfGenotypes {
        n_samples,
        gt: None,
        dp: None,
        gq: None,
        ft: None,
    };
    for _ in 0..n_fmt {
        let k = typed::read(&mut indiv)?;
        let key = match k.int(0) {
            Some(Int::Value(key)) if k.n == 1 && key >= 0 => key as usize,
            _ => return Err("a FORMAT key is not a single non-negative integer".into()),
        };
        let (ty, n) = typed::read_descriptor(&mut indiv)?;
        let len = n
            .checked_mul(ty.width())
            .and_then(|l| l.checked_mul(n_samples))
            .filter(|&l| l <= indiv.len())
            .ok_or("a FORMAT field runs past the end of the record")?;
        let (data, rest) = indiv.split_at(len);
        indiv = rest;
        let slot = match dicts.strings.get(key).map(Vec::as_slice) {
            Some(b"GT") => &mut out.gt,
            Some(b"DP") if fields.dp => &mut out.dp,
            Some(b"GQ") if fields.gq => &mut out.gq,
            Some(b"FT") if fields.ft => &mut out.ft,
            Some(_) => continue,
            None => return Err(format!("FORMAT key {key} is not in the header dictionary")),
        };
        slot.get_or_insert(Col { ty, n, data });
    }
    Ok(out)
}

fn int_field(col: Option<Col<'_>>, i: usize, what: &str) -> Result<Option<u32>, String> {
    let Some(col) = col else { return Ok(None) };
    match col.ty {
        Ty::Missing => Ok(None),
        Ty::Int8 | Ty::Int16 | Ty::Int32 => {
            if col.n == 0 {
                return Ok(None);
            }
            Ok(match col.of(i).int(0) {
                Some(Int::Value(v)) if v >= 0 => Some(v as u32),
                _ => None,
            })
        }
        _ => Err(format!("FORMAT/{what} is not an integer field")),
    }
}

impl GenotypeSource for BcfGenotypes<'_> {
    fn n_samples(&self) -> usize {
        self.n_samples
    }

    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        out.reset();
        out.clear_gt();
        out.dp = int_field(self.dp, i, "DP")?;
        out.gq = int_field(self.gq, i, "GQ")?;
        if let Some(ft) = self.ft {
            out.ft_pass = match ft.ty {
                Ty::Char => {
                    let d = ft.of(i).data;
                    let end = d.iter().rposition(|&b| b != 0).map_or(0, |p| p + 1);
                    ft_passes(&d[..end])
                }
                _ => true,
            };
        }
        let Some(gt) = self.gt else { return Ok(()) };
        if matches!(gt.ty, Ty::Missing) {
            return Ok(());
        }
        if !matches!(gt.ty, Ty::Int8 | Ty::Int16 | Ty::Int32) {
            return Err("FORMAT/GT is not an integer field".into());
        }
        let v = gt.of(i);
        let mut phased_seps = 0usize;
        let mut ploidy = 0usize;
        for j in 0..gt.n {
            let allele = match v.int(j) {
                Some(Int::End) | None => break,
                Some(Int::Missing) => None,
                Some(Int::Value(x)) if x < 0 => return Err(format!("GT value {x} is negative")),
                Some(Int::Value(x)) => {
                    if j > 0 && x & 1 == 1 {
                        phased_seps += 1;
                    }
                    u32::try_from(x >> 1).ok().and_then(|a| a.checked_sub(1))
                }
            };
            if ploidy == 0 {
                out.gt_buf();
            }
            ploidy += 1;
            out.gt
                .as_mut()
                .expect("buffer just created")
                .alleles
                .push(allele);
        }
        if let Some(g) = out.gt.as_mut() {
            g.phased = ploidy >= 2 && phased_seps == ploidy - 1;
        }
        Ok(())
    }
}
