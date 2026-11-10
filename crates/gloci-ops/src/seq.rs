//! Sequence helpers shared by `getfasta` and `nuc`/`maskfasta`: the
//! complement table and the summary a sequence command returns.

/// What a sequence command did with its intervals.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeqSummary {
    /// Intervals whose sequence was written.
    pub written: u64,
    /// Zero-length intervals, which have no sequence and are skipped.
    pub points: u64,
    /// Intervals outside the FASTA skipped under `--skip-invalid`.
    pub skipped: u64,
}

const fn table(rna: bool) -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    // (upper-case base, its complement); the lower-case pair is derived.
    // U complements to A, but A to U only for RNA; S, W and N are their own.
    let pairs: [(u8, u8); 13] = [
        (b'A', if rna { b'U' } else { b'T' }),
        (b'T', b'A'),
        (b'U', b'A'),
        (b'C', b'G'),
        (b'G', b'C'),
        (b'R', b'Y'),
        (b'Y', b'R'),
        (b'K', b'M'),
        (b'M', b'K'),
        (b'B', b'V'),
        (b'V', b'B'),
        (b'D', b'H'),
        (b'H', b'D'),
    ];
    let mut p = 0;
    while p < pairs.len() {
        let (b, c) = pairs[p];
        t[b as usize] = c;
        t[(b + 32) as usize] = c + 32;
        p += 1;
    }
    t
}

static DNA: [u8; 256] = table(false);
static RNA: [u8; 256] = table(true);

/// The complement of `b` as bedtools writes it: IUPAC-aware (`R`/`Y`,
/// `K`/`M`, `B`/`V`, `D`/`H`), case kept, `U` to `A`; `S`, `W`, `N` and every
/// other byte unchanged. With `rna`, `A` gives `U` instead of `T`.
pub fn complement(b: u8, rna: bool) -> u8 {
    if rna {
        RNA[b as usize]
    } else {
        DNA[b as usize]
    }
}

/// Reverse-complement `seq` in place.
pub fn revcomp_in_place(seq: &mut [u8], rna: bool) {
    let t = if rna { &RNA } else { &DNA };
    seq.reverse();
    for b in seq.iter_mut() {
        *b = t[*b as usize];
    }
}
