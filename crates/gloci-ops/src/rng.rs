//! The seeded random stream (D143).
//!
//! `rand_chacha`'s `ChaCha8Rng` is used only through `next_u64`, `set_stream`
//! and `set_word_pos`; every draw is gloci's own so the output does not
//! depend on a `rand` distribution implementation. The key is the seed in
//! little-endian bytes followed by 24 zero bytes; the stream number is the
//! unit of work (D144), so a result does not depend on thread scheduling.

use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};
use std::hash::{BuildHasher, Hasher};

/// One deterministic stream of random words.
pub struct Stream {
    rng: ChaCha8Rng,
}

impl Stream {
    /// The stream numbered `stream` of the key made from `seed`, at its start.
    pub fn new(seed: u64, stream: u64) -> Stream {
        let mut key = [0u8; 32];
        key[..8].copy_from_slice(&seed.to_le_bytes());
        let mut rng = ChaCha8Rng::from_seed(key);
        rng.set_stream(stream);
        Stream { rng }
    }

    /// Like [`Stream::new`], positioned `word` 32-bit words into the stream,
    /// so a unit of work can start at its own offset.
    pub fn at_word(seed: u64, stream: u64, word: u128) -> Stream {
        let mut s = Stream::new(seed, stream);
        s.rng.set_word_pos(word);
        s
    }

    pub fn next_u64(&mut self) -> u64 {
        self.rng.next_u64()
    }

    /// A uniform integer in `0..n` (`n > 0`; `n == 0` returns 0): Lemire's
    /// multiply-shift with rejection.
    pub fn below(&mut self, n: u64) -> u64 {
        if n == 0 {
            return 0;
        }
        let mut m = u128::from(self.next_u64()) * u128::from(n);
        if (m as u64) < n {
            let t = n.wrapping_neg() % n;
            while (m as u64) < t {
                m = u128::from(self.next_u64()) * u128::from(n);
            }
        }
        (m >> 64) as u64
    }

    /// A uniform float in `[0, 1)` with 53 random bits.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// A fair coin: the top bit of a word.
    pub fn coin(&mut self) -> bool {
        self.next_u64() >> 63 == 1
    }
}

/// The first four words of block 0 of stream `stream` (what `Stream::new`
/// yields as its first two `next_u64` draws), from one scalar ChaCha8 block.
/// `rand_chacha` fills four blocks per refill; a unit of work that draws two
/// words per stream (`random`) does not need the other three.
pub fn first_words(seed: u64, stream: u64) -> [u64; 2] {
    let mut init = [0u32; 16];
    init[..4].copy_from_slice(&[0x6170_7865, 0x3320_646e, 0x7962_2d32, 0x6b20_6574]);
    init[4] = seed as u32;
    init[5] = (seed >> 32) as u32;
    // words 6..12 (the rest of the key) and the block counter 12, 13 are 0.
    init[14] = stream as u32;
    init[15] = (stream >> 32) as u32;
    let mut x = init;
    macro_rules! qr {
        ($a:expr, $b:expr, $c:expr, $d:expr) => {
            x[$a] = x[$a].wrapping_add(x[$b]);
            x[$d] = (x[$d] ^ x[$a]).rotate_left(16);
            x[$c] = x[$c].wrapping_add(x[$d]);
            x[$b] = (x[$b] ^ x[$c]).rotate_left(12);
            x[$a] = x[$a].wrapping_add(x[$b]);
            x[$d] = (x[$d] ^ x[$a]).rotate_left(8);
            x[$c] = x[$c].wrapping_add(x[$d]);
            x[$b] = (x[$b] ^ x[$c]).rotate_left(7);
        };
    }
    for _ in 0..4 {
        qr!(0, 4, 8, 12);
        qr!(1, 5, 9, 13);
        qr!(2, 6, 10, 14);
        qr!(3, 7, 11, 15);
        qr!(0, 5, 10, 15);
        qr!(1, 6, 11, 12);
        qr!(2, 7, 8, 13);
        qr!(3, 4, 9, 14);
    }
    let w = |i: usize| u64::from(x[i].wrapping_add(init[i]));
    [w(0) | (w(1) << 32), w(2) | (w(3) << 32)]
}

/// `seed` itself, or, when none was given, a seed from the platform's
/// `RandomState` with a warning that says how to reproduce the run.
pub fn seed_or_entropy(seed: Option<u64>) -> u64 {
    if let Some(s) = seed {
        return s;
    }
    let s = std::collections::hash_map::RandomState::new()
        .build_hasher()
        .finish();
    log::warn!("no --seed given; using --seed {s} (pass it to reproduce this output)");
    s
}

#[cfg(test)]
mod first_words_tests {
    use super::*;

    #[test]
    fn first_words_match_the_stream() {
        for seed in [0u64, 1, 42, u64::MAX, 0x1234_5678_9abc_def0] {
            for stream in [0u64, 1, 2, 4095, 9_999_999, u64::MAX, 1 << 33] {
                let mut s = Stream::new(seed, stream);
                let want = [s.next_u64(), s.next_u64()];
                assert_eq!(first_words(seed, stream), want, "{seed} {stream}");
            }
        }
    }
}
