//! The windows of the windowed statistics (decision D159): `[k·STEP + 1,
//! k·STEP + SIZE]` on each chromosome, k = 0, 1, …, as vcftools lays them
//! out for `--window-pi` and `--fst-window-size`.

use std::collections::VecDeque;

/// The open windows of one chromosome, each holding a statistic's sums `W`.
///
/// A window opens when a site first adds to it. With input sorted by
/// position, a window below the first one holding a site can get nothing
/// more, so the statistic takes it out ([`OpenWindows::pop_before`]) and
/// writes it; at most `⌈SIZE/STEP⌉ + 1` windows are open at once.
pub(super) struct OpenWindows<W> {
    size: u64,
    step: u64,
    /// The index `k` of the first open window.
    base: u64,
    open: VecDeque<W>,
    /// The most windows open at once.
    peak: usize,
}

impl<W: Clone + Default> OpenWindows<W> {
    /// `size` and `step` must be positive, `step` at most `size`.
    pub(super) fn new(size: u64, step: u64) -> OpenWindows<W> {
        assert!(
            size > 0 && step > 0 && step <= size,
            "a window needs 0 < step <= size"
        );
        OpenWindows {
            size,
            step,
            base: 0,
            open: VecDeque::new(),
            peak: 0,
        }
    }

    /// The index of the first window holding `pos`: the least k ≥ 0 with
    /// `pos <= k·STEP + SIZE`.
    pub(super) fn first_holding(&self, pos: u64) -> u64 {
        if pos > self.size {
            (pos - self.size).div_ceil(self.step)
        } else {
            0
        }
    }

    /// The first and last position of window `k`.
    pub(super) fn bounds(&self, k: u64) -> (u64, u64) {
        (k * self.step + 1, k * self.step + self.size)
    }

    /// Take out the first open window, with its index, if that index is
    /// below `end` (`u64::MAX` takes every window in turn).
    pub(super) fn pop_before(&mut self, end: u64) -> Option<(u64, W)> {
        if self.base >= end {
            return None;
        }
        let w = self.open.pop_front()?;
        let k = self.base;
        self.base += 1;
        Some((k, w))
    }

    /// Apply `f` to every window holding `pos`, opening the ones not open
    /// yet. The windows below [`OpenWindows::first_holding`] must have been
    /// taken out first.
    pub(super) fn add(&mut self, pos: u64, mut f: impl FnMut(&mut W)) {
        let first = self.first_holding(pos);
        // k·STEP + 1 <= pos holds for k < ⌈pos/STEP⌉ (none for POS 0).
        let last = pos.div_ceil(self.step);
        if self.open.is_empty() {
            self.base = first;
        }
        debug_assert_eq!(self.base, first, "windows below the site are still open");
        let have = self.base + self.open.len() as u64;
        if last > have {
            self.open
                .resize(self.open.len() + (last - have) as usize, W::default());
            self.peak = self.peak.max(self.open.len());
        }
        let from = (first - self.base) as usize;
        let to = (last - self.base) as usize;
        self.open.range_mut(from..to).for_each(&mut f);
    }

    /// The most windows held open at once so far.
    pub(super) fn peak(&self) -> usize {
        self.peak
    }
}

#[cfg(test)]
mod tests {
    use super::OpenWindows;

    /// The windows each position falls in, by taking out every window below
    /// the position before adding, as a statistic does.
    fn placed(size: u64, step: u64, positions: &[u64]) -> Vec<(u64, u64, u32)> {
        let mut w: OpenWindows<u32> = OpenWindows::new(size, step);
        let mut out = Vec::new();
        for &pos in positions {
            let first = w.first_holding(pos);
            while let Some((k, n)) = w.pop_before(first) {
                let (a, b) = w.bounds(k);
                out.push((a, b, n));
            }
            w.add(pos, |n| *n += 1);
        }
        while let Some((k, n)) = w.pop_before(u64::MAX) {
            let (a, b) = w.bounds(k);
            out.push((a, b, n));
        }
        out
    }

    #[test]
    fn windows_hold_every_position_once_per_window() {
        assert_eq!(
            placed(10, 5, &[1, 10, 11, 30]),
            vec![
                (1, 10, 2),
                (6, 15, 2),
                (11, 20, 1),
                // 16-25 holds no position, so it never opens
                (21, 30, 1),
                (26, 35, 1),
            ]
        );
        // POS 0 is in no window; a step that does not divide the size.
        assert_eq!(placed(10, 3, &[0, 4]), vec![(1, 10, 1), (4, 13, 1)]);
    }

    #[test]
    fn at_most_size_over_step_plus_one_windows_are_open() {
        for (size, step) in [(10, 10), (10, 5), (10, 3), (100, 7)] {
            let mut w: OpenWindows<u32> = OpenWindows::new(size, step);
            for pos in (1..2000).step_by(3) {
                let first = w.first_holding(pos);
                while w.pop_before(first).is_some() {}
                w.add(pos, |n| *n += 1);
            }
            assert!(w.peak() as u64 <= size.div_ceil(step) + 1, "{size}/{step}");
        }
    }
}
