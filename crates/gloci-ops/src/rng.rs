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
