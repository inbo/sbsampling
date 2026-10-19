//! Self-contained, value-stable randomness for the samplers.
//!
//! Output must never change between builds, because sort keys and Owen
//! scrambles are part of the sampling design (reproducible master samples).
//! Everything here is plain integer arithmetic; the only external piece is the
//! PCG-64 MCG generator, whose output sequence is fixed by its published
//! algorithm.

use rand::RngCore;
use rand_pcg::Pcg64Mcg;

/// SplitMix64 step (Steele, Lea & Flood 2014).
#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// One-shot 64-bit mixer, for combining a seed with an identifier.
#[inline]
pub fn mix(x: u64) -> u64 {
    let mut s = x;
    splitmix64(&mut s)
}

/// Seed a PCG-64 MCG from a u64 via two SplitMix64 outputs
/// (`Pcg64Mcg::new` forces the state odd).
pub fn seeded_rng(seed: u64) -> Pcg64Mcg {
    let mut s = seed;
    let hi = splitmix64(&mut s);
    let lo = splitmix64(&mut s);
    Pcg64Mcg::new(((hi as u128) << 64) | lo as u128)
}

/// Unbiased integer in `[0, n)`, `n > 0` (Lemire 2019, multiply-shift with
/// rejection).
#[inline]
fn bounded(rng: &mut Pcg64Mcg, n: u64) -> u64 {
    let mut m = (rng.next_u64() as u128) * (n as u128);
    let mut l = m as u64;
    if l < n {
        let t = n.wrapping_neg() % n;
        while l < t {
            m = (rng.next_u64() as u128) * (n as u128);
            l = m as u64;
        }
    }
    (m >> 64) as u64
}

/// Fisher–Yates shuffle (descending index variant).
pub fn shuffle<T>(rng: &mut Pcg64Mcg, xs: &mut [T]) {
    for i in (1..xs.len()).rev() {
        let j = bounded(rng, (i + 1) as u64) as usize;
        xs.swap(i, j);
    }
}

/// High 32 bits of the next 64-bit output.
pub fn next_u32_hi(rng: &mut Pcg64Mcg) -> u32 {
    (rng.next_u64() >> 32) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    // Golden values: if any of these change, every sample changes.
    #[test]
    fn pcg_stream_is_pinned() {
        assert_eq!(seeded_rng(42).next_u64(), 0x83ae_6410_88e2_559a);
    }

    #[test]
    fn mix_is_pinned() {
        assert_eq!(mix(42), 0xbdd7_3226_2feb_6e95);
    }

    #[test]
    fn shuffle_is_pinned() {
        let mut a = [0u64, 1, 2, 3, 4, 5, 6];
        shuffle(&mut seeded_rng(0), &mut a);
        assert_eq!(a, [6, 5, 0, 1, 3, 2, 4]);

        let mut b = [0u64, 1, 2, 3, 4, 5, 6];
        shuffle(&mut seeded_rng(42), &mut b);
        assert_eq!(b, [5, 6, 1, 4, 2, 0, 3]);

        let mut c: Vec<u64> = (0..122).collect();
        shuffle(&mut seeded_rng(42), &mut c);
        assert_eq!(&c[..8], &[23, 39, 107, 26, 70, 108, 43, 17]);
    }

    #[test]
    fn shuffle_is_a_permutation() {
        let mut v: Vec<u64> = (0..122).collect();
        shuffle(&mut seeded_rng(7), &mut v);
        let mut s = v.clone();
        s.sort_unstable();
        assert_eq!(s, (0..122).collect::<Vec<u64>>());
    }
}
